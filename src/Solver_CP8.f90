!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE SOLVER8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module manages the solvers for applying matrix-vector products

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE INPUTCP
      USE RESTART
      USE BLOCKPOWER
      USE BLOCKUTILS
      USE FEAST8
      USE ALSDRVR
      USE CPr8
      USE ALS8DRVR

      implicit none
      real(kind=8), allocatable, private :: itn_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_Solver_Module_CP8()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(itn_time(mpinodes))
      itn_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_Solver_Module_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_Solver_Module_CP8()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_Solver_Module_CP8()
      call Get_MPI_Timings('* Solver CP8 module (iterations)',itn_time)
      MODULE_SETUP = .FALSE.
      deallocate(itn_time)

      end subroutine Dispose_Solver_Module_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      integer function GetSolverType_CP8(cpp,Q) result(styp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CP8), ALLOCATABLE, INTENT(INOUT) :: Q(:)

!     Determine styp
      IF ((cpp%solver.seq.'powr') .or. (cpp%solver.seq.'pALS')) THEN
         styp=1
      ELSEIF (cpp%solver .seq. 'iitn') THEN
         styp=2
      ELSEIF (cpp%solver .seq. 'iitf') THEN
         styp=3
      ELSE
         call &
         AbortWithError('GetSolverType_CP8(): Solver not recognized')
      ENDIF

      end function GetSolverType_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function useDiag(nev,nbas) result(diag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines if diagonalization should be used

      implicit none
      integer, intent(in) :: nev
      integer, intent(in) :: nbas(:)
      logical :: diag
      integer :: j,ndof,prod

      ndof=SIZE(nbas)

!     Diagonalize H directly if ALL states are requested
      diag=.FALSE.
      prod=1
      DO j=1,ndof
         prod=prod*nbas(j)
         IF (prod.gt.nev) EXIT
      ENDDO
      IF (prod.eq.nev) diag=.TRUE.

      end function useDiag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WrapSolver_CP8(eigv,delta,bounds,ML,cpp,Q,H,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for computing the eigenfunctions and 
! eigenvalues using the solver of choice

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CP), INTENT(IN)  :: H,W
      TYPE (CP8), ALLOCATABLE :: Q8(:)
      TYPE (CP8) :: H8,W8
      real(kind=8), allocatable, intent(inout) :: eigv(:),delta(:)
      real(kind=8), intent(inout) :: bounds(2)
      integer :: i,nbloc

      select case(cpp%algo)
      case(0)
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(X,A/)') &
         'Entering Solver: CPU algorithm selected...'
      case(1)
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(X,A/)') &
         'Entering Solver: GPU algorithm selected...'
      case default
        call AbortWithError('WrapSolver_CP8(): unrecognized algo')
      end select

      nbloc=SIZE(Q)
      allocate(Q8(nbloc))
      call H8%fromCP(H)
      if (cpp%algo.eq.1) call H8%copyintodevice()
      if (W%R().gt.0) call W8%fromCP(W)
      do i=1,nbloc
         call Q8(i)%fromCP(Q(i))
         call FlushCP(Q(i))
      enddo

      call set_als_settings_ALS8(cpp%alspenalty,cpp%als_linsys_alg)
      call SolveHPsi_CP8(eigv,delta,bounds,ML,cpp,Q8,H8,W8)

      do i=1,nbloc
         call Q8(i)%toCP(Q(i))
         call Q8(i)%flush()
      enddo
      deallocate(Q8)

      call H8%flush()
      if (W%R().gt.0) call W8%flush()

      end subroutine WrapSolver_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveHPsi_CP8(eigv,delta,bounds,ML,cpp,Q,H,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for computing the eigenfunctions and 
! eigenvalues using the solver of choice

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CP8), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CP8), INTENT(IN)  :: H,W
      real(kind=8), allocatable, intent(inout) :: eigv(:)
      real(kind=8), allocatable, intent(inout) :: delta(:)
      real(kind=8), allocatable  :: eigtmp(:)
      real(kind=8), intent(inout) :: bounds(2)
      integer :: i,j,nev,nup,ndown,nsame,nloc,ist,styp,nsubm,dummy
      logical :: conv,diag,tiledH,readsuccess,calcbounds
      real(kind=8)  :: rmsdelta,oldrms,maxdelta,sumdelta
      real(kind=8), parameter :: redtol=1.d-12

!     Initializations
      nev=SIZE(eigv)
      nsubm=H%D()
      ist=0
      oldrms=1.d99
      diag=useDiag(nev,Q(1)%nbas)
      styp=GetSolverType_CP8(cpp,Q)

!     Determine if H should be tiled over MPI ranks
      tiledH=(cpp%algo.ge.0 .and. cpp%lowmem.gt.2 .and. nsubm.ge.2 & 
              .and. .not.(nsubm.eq.2 .and. (cpp%red2D.seq.'SVD')) &
              .and. mpinodes.gt.1)

!     Read psi file if this is a restart
      call ReadPsi_CP8(ist,bounds,eigv,delta,Q,ML,readsuccess)

!     Calculate the spectral range of H
      IF (cpp%npow.gt.0 .and. cpp%ncycle.gt.0 .and. (.not.diag) &
          .and. ((readsuccess .and. (cpp%calcbounds.seq.'recalc')) &
          .or. (.not.readsuccess))) THEN
         calcbounds=(.not.(cpp%calcbounds.seq.'guess'))
         call GetSpectralRange_CP8(min(50,cpp%psirank),10,5,Q,H,bounds,&
                               cpp%algo,calcbounds,cpp%padbounds,tiledH)
      ENDIF

!     Extend guess vectors to target rank
      DO i=1,nev
         call Q(i)%extendrand(cpp%psirank)
         if (cpp%algo.eq.1) call Q(i)%copyintodevice()
      ENDDO

!     Initial guess and pre-diagonalization
      IF (readsuccess) THEN
         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,*)
            write(*,*) 'Eigenvalues read: ',ist,(eigv(j),j=1,nev)
         ENDIF
      ELSE
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,*) 'Initial guess   : ',0,(eigv(j),j=1,nev)
         IF (cpp%ncycle.gt.0 .and. (cpp%diag.or.diag)) THEN
            call Diagonalize_CP8(Q,H,eigv,cpp%ovrlpenalty,0.d0,dummy,&
                                 cpp%psinals,cpp%algo,cpp%reduceHQ,&
                                 cpp%update,tiledH)
            IF (mpirank.eq.mpi_prnt_rank) THEN
               write(*,*)
               write(*,*) 'Diagonalization : ',0,(eigv(j),j=1,nev)
            ENDIF
         ENDIF
         call SavePsi_CP8(0,bounds,eigv,delta,Q,ML)
      ENDIF

!     If all the eigenvalues were requested, exit here since the 
!     pre-diagonalization already gives the exact answer
      IF (diag) RETURN

      ALLOCATE(eigtmp(nev))
      IF (ist.eq.0) delta=1.d99
      nloc=0

!     Loop over solver cycles
      DO i=1,cpp%ncycle

         IF (i.le.ist) CYCLE

!        Copy eigenvalues to temp array
         eigtmp=eigv

!        Run iterations using the solver of choice
!!! UNDER CONSTRUCTION: call FEAST instead of Iterate when implemented
!         call FEAST_iterate_CP8(Q,H,W,eigtmp,eigv,cpp,nloc,i,tiledH)
!!! END CONSTRUCTION
         call Iterate_CP8(Q,H,W,eigtmp,eigv,cpp,bounds,nloc,i,styp,&
                          tiledH)

!        Check for convergence (rms change < tol for all eigenvalues  
!        excluding the top) and exit if achieved
         conv=.TRUE.
         rmsdelta=0.d0
         maxdelta=0.d0
         sumdelta=0.d0
         nup=0
         ndown=0
         nsame=0
         DO j=1,nev
            delta(j)=(eigtmp(j)-eigv(j))/abs(eigv(j))
            rmsdelta=rmsdelta+delta(j)**2
            sumdelta=sumdelta+delta(j)
            IF (abs(delta(j)).lt.cpp%solvtol) THEN
               nsame=nsame+1
            ELSEIF (eigtmp(j).gt.eigv(j)) THEN
               ndown=ndown+1
            ELSE
               nup=nup+1
            ENDIF
            IF (abs(delta(j)).gt.maxdelta) maxdelta=abs(delta(j))
         ENDDO
         rmsdelta=sqrt(rmsdelta)/nev
         IF (rmsdelta.gt.cpp%solvtol) conv=.FALSE.

         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(X,A,ES11.4,A,I5)') 'relative rms delta = ',&
                 rmsdelta,' ; nup   = ',nup
            write(*,'(X,A,ES11.4,A,I5)') 'relative max delta = ',&
                 maxdelta,' ; ndown = ',ndown
            write(*,'(X,A,ES11.4,2(A,I5))') 'relative sum delta = ',&
                 sumdelta,' ; nsame = ',nsame,' ; nlock = ',nloc
         ENDIF

!        Save the wavefunction from the current cycle
         call SavePsi_CP8(i,bounds,eigv,delta,Q,ML)

!        Exit if energies are converged
         IF (conv) THEN
            IF (mpirank.eq.mpi_prnt_rank) &
               write(*,'(/X,2A,ES11.4)') 'Relative rms error of all ',&
               'states converged to within: ',cpp%solvtol
            EXIT
         ENDIF
         oldrms=rmsdelta
      ENDDO  ! Loop over cycles

      DEALLOCATE(eigtmp)

      end subroutine SolveHPsi_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Iterate_CP8(Q,H,W,eigvo,eigv,cpp,bounds,nconv,i,styp,&
                              tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Performs iterations of various types

      implicit none
      TYPE (CPpar), INTENT(IN)  :: cpp
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8), INTENT(IN)    :: H,W
      TYPE (CP8) :: Qdummy
      integer, intent(inout) :: nconv
      integer, intent(in)    :: i,styp
      logical, intent(in)    :: tiledH
      real(kind=8), intent(inout) :: eigv(:)
      real(kind=8), intent(in)    :: eigvo(:),bounds(2)
      real(kind=8), parameter     :: tol=1.d-15
      character(len=18) :: tag
      real(kind=8)  :: Eshift,rq
      integer :: j,k,nbloc,sz,os,szmx
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_Solver_Module_CP8()

!     Easy exit for zero iterations
      IF (cpp%npow.lt.1) RETURN

      call CPU_TIME(ti1)

!     Set parameters
      nbloc=SIZE(Q)
      call GetBlockShift(eigv,bounds,Eshift)

      IF (styp.eq.1) THEN
         tag='ALS8-Power cycle: '
      ELSEIF (styp.eq.2) THEN
         tag='Inv.it.N-8 cycle: '
      ELSEIF (styp.eq.3) THEN
         tag='Inv.it.F-8 cycle: '
      ELSE
         call AbortWithError('Iterate_CP8(): invalid solver type')
      ENDIF

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc-nconv,sz,os,szmx)
      os=os+nconv

!     Run power iterations on each vector in the block
      call nvtx_start('Solver Iterations')
      DO j=1,sz
         k=os+j
         IF (styp.eq.1) THEN
            rq=ALS_pow_alg_CP8(H,1,Eshift,Q(k),cpp%npow,cpp%algo,tiledH)
         ELSEIF (styp.eq.2) THEN
            call ALS_invpow_normal_CP8(H,1,eigv(k),Q(k),cpp%psinals,&
                                     cpp%npow,cpp%algo)
         ELSEIF (styp.eq.3) THEN
            call ALS_invpow_fast_CP8(H,1,eigv(k),Q(k),cpp%psinals,&
                                     cpp%npow,cpp%algo,tiledH)
         ENDIF
         call Q(k)%updatefromdevice()
      ENDDO

!     If using the H-tiled algorithm, every rank must participate in
!     cycling portions of H. However, if the number of vectors to 
!     iterate doesn't evenly divide the number of MPI ranks, then some
!     ranks will make fewer calls to the cycling routine. In this case,
!     include a additional vector to preserve the number of MPI calls.
      IF (sz.lt.szmx .and. tiledH) THEN
         call Qdummy%copyfrom(Q(1))
         IF (styp.eq.1) THEN
            rq=ALS_pow_alg_CP8(H,1,Eshift,Qdummy,cpp%npow,cpp%algo,tiledH)
         ELSEIF (styp.eq.3) THEN
                 write(*,*) '(Eshift = ',Eshift,')'
            call ALS_invpow_fast_CP8(H,1,Eshift,Qdummy,cpp%psinals,&
                                     cpp%npow,cpp%algo,tiledH)
         ELSE
            write(*,*) 'MPI rank: ',mpirank,&
            ': tiled H only implemented for intertwined solver'
            call AbortWithError("Error in Iterate_CP8()")
         ENDIF
         call Qdummy%flush()
      ENDIF
      call nvtx_stop()

!     Sync vectors before Gram-Schmidt, updates
      DO j=1,nbloc
         call Q(j)%deletefromdevice()
      ENDDO
      call MPI_Sync_block_CP8(Q)
      if (cpp%algo.eq.1) then
         DO j=1,nbloc
            call Q(j)%copyintodevice()
         ENDDO
      endif

      call CPU_TIME(ti2)
      itn_time=itn_time+ti2-ti1

!     Orthogonalization and update/vector sort
      call Orthogonalize_CP8(Q,nconv,cpp%psinals,cpp%algo,cpp%orthogalg)
      IF (cpp%diag) call &
         Diagonalize_CP8(Q,H,eigv,cpp%ovrlpenalty,bounds(1),nconv,&
                         cpp%psinals,cpp%algo,cpp%reduceHQ,cpp%update,&
                         tiledH)

!     Test convergence on eigenvalues and "lock" converged vectors
      nconv=0
      DO j=1,nbloc
         IF (abs((eigv(j)-eigvo(j))/eigv(j)).gt.tol) EXIT
         nconv=j
      ENDDO

!     Print the eigenvalues
      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,*)
         write(*,*) tag,i,(eigv(j),j=1,nbloc)
      ENDIF

      end subroutine Iterate_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetSpectralRange_CP8(rk,npow,ncyc,Q,H,bounds,algo,&
                                      calcbounds,padbounds,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Estimates the spectral range of the Hamiltonian using power method

      implicit none
      TYPE (CP8), INTENT(IN) :: H
      TYPE (CP8) :: Q(:)
      TYPE (CP8), allocatable :: Qt(:)
      integer, intent(in) :: rk,npow,ncyc,algo
      logical, intent(in) :: calcbounds,tiledH
      real(kind=8), intent(in) :: padbounds
      real(kind=8), intent(inout) :: bounds(2)
      integer :: i,j,ndof
      integer :: ishs(2)
      real(kind=8) :: btmp(2)
      real(kind=8), allocatable :: bhist(:,:)
      real(kind=8), parameter :: smallnr = 1.d-15
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_Solver_Module_CP8()
      call CPU_TIME(ti1)
      call nvtx_start('GetSpectralRange')

      ndof=Q(1)%D()
      btmp=0.d0
      ishs=(/1,0/)

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(X,A/)') 'Power spectral range estimation'
         write(*,'(A,X,2(9X,A))') '  i','boundsl','boundsu'
         write(*,'(3X,2(X,f16.6),X,A)') bounds(1),bounds(2),&
                                         '(uncoupled guess)'
      ENDIF

      IF (calcbounds) THEN

         ALLOCATE(Qt(2))
         ALLOCATE(bhist(ncyc,2))

!        Guesses for highest and lowest eigenvectors
         DO j=1,2
            call Qt(j)%clone0(Q(1))
            Qt(j)%coef(1)=1.d0
            Qt(j)%base(:)=smallnr
         ENDDO

         DO i=1,ndof
            Qt(1)%base(Qt(1)%BS(1,i))=1.d0
            Qt(2)%base(Qt(2)%BF(1,i))=1.d0
         ENDDO

         DO j=1,2
            call Qt(j)%extendrand(rk)
            if (algo.eq.1) call Qt(j)%copyintodevice()
         ENDDO

!        Run the power method to improve the bounds
         DO j=2,1,-1
            DO i=1,ncyc
               bhist(i,j)=ALS_pow_alg_CP8(H,ishs(j),btmp(j),Qt(j),npow,&
                                          algo,tiledH)
               if (j.eq.1) btmp(1)=0.5*(bhist(i,1)+bhist(ncyc,2))
            ENDDO
            if (j.eq.2) btmp(1)=bhist(ncyc,2)
         ENDDO

         bounds(1)=minval(bhist(ncyc,:))
         bounds(2)=maxval(bhist(ncyc,:))

         IF (mpirank.eq.mpi_prnt_rank) THEN
            DO i=1,ncyc
               write(*,'(i3,2(X,f16.6))') i,bhist(i,1),bhist(i,2)
            ENDDO
         ENDIF

!        Clean up
         DO j=1,2
            call Qt(j)%flush()
         ENDDO
         DEALLOCATE(Qt,bhist)
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(/X,A,2(f15.6,A))') 'Spectral range of H = [',&
                                      bounds(1),',',bounds(2),']'
      ENDIF

      btmp(2)=padbounds*(bounds(2)-bounds(1))
      bounds(1)=bounds(1)-btmp(2)
      bounds(2)=bounds(2)+btmp(2)

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A,ES9.2,A,2(f15.6,A)/)') 'Padded by ',padbounds,&
                                    ' = [',bounds(1),',',bounds(2),']'

      call nvtx_stop()
      call CPU_TIME(ti2)
      itn_time=itn_time+ti2-ti1

      end subroutine GetSpectralRange_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE SOLVER8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
