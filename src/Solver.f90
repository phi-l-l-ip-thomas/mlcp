!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE SOLVER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module manages the solvers for applying matrix-vector products

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE INPUTCP
      USE SEPDREPN
      USE RESTART
      USE BLOCKPOWER
      USE ALSPOW
      USE ALSUTILS
      USE ALSDRVR

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      integer function GetSolverType(cpp,Q) result(styp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)

!     Determine styp
      IF (cpp%solver .seq. 'powr') THEN
         styp=1
      ELSEIF (cpp%solver .seq. 'pALS') THEN
!        Use ALS-guided power method unless 2D with SVD reduction
         IF (SIZE(Q(1)%nbas).eq.2 .and. cpp%red2D.eq.'SVD') THEN
            IF (mpirank.eq.mpi_prnt_rank) call ShowWarning(&
            'SVD reduction selected; using ordinary power iterations')
            styp=1
         ELSE
            styp=-1
         ENDIF
      ELSE
         call AbortWithError('GetSolverType(): Solver not recognized')
      ENDIF

      end function GetSolverType

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowPsiMem(eigv,cpp,Q,H,styp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for computing the eigenfunctions and 
! eigenvalues using the solver of choice:
! styp=0 -> power iteration
! styp=3 -> ALS-guided power method ("intertwining")

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (CP), INTENT(IN) :: H
      integer, intent (in) :: styp
      real*8, intent(in)   :: eigv(:)
      integer :: htrm,ndof,rF,rG,nev,ncpu,npara,veclen,nmax
      real*8  :: GB,mvGB,GSGB,QHQGB,upGB,redGB,Blen,mvlen
      logical :: useSVD

!     Initializations
      htrm=SIZE(H%coef)+1 ! Plus 1 due to E-shift
      ndof=SIZE(Q(1)%nbas)
      nmax=maxval(Q(1)%nbas)
      rF=cpp%psirank
      nev=SIZE(eigv)
      ncpu=cpp%ncpu
      npara=min(ncpu,nev)
      GB=1.d0/(2**27)
      mvGB=0.d0  ! matrix-vector product memory
      upGB=0.d0  ! vector update memory
      GSGB=0.d0  ! Gram-Schmidt memory
      redGB=0.d0 ! reduction memory

      veclen=SIZE(Q(1)%base,1)+1
      Blen=nev*REAL(veclen)*rF
      mvlen=htrm*REAL(veclen)*rF

!     Different memory requirements for reduction if SVD is used
      useSVD=.FALSE.
      IF (ndof.eq.2 .and. cpp%red2D.eq.'SVD') useSVD=.TRUE.

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

!     Print parameters affecting the memory usage
      write(*,'(3X,A)') '*** Solver memory usage ***'
      write(*,'(7X,A,30X,A,X,I11)') 'Number of CPUs,','(N): ',ncpu
      write(*,'(7X,A,31X,A,X,I11)') 'Vector length,','(V): ',veclen
      write(*,'(7X,A,34X,A,X,I11)') 'Block size,','(B): ',nev
      write(*,'(7X,A,14X,A,X,I11)') 'Terms in H (including E-shift),',&
            '(T): ',htrm
      write(*,'(7X,A,40X,A,X,I11)') 'Rank,','(R): ',rF
      write(*,'(7X,A,39X,A,X,I11)') 'Max n,','(X): ',nmax
      write(*,'(7X,A,31X,A,X,I11)') 'Number of DOF,','(D): ',ndof
      write(*,'(7X,A,4X,A,X,I11)') 'Number of parallelized vectors,',&
           'min(N,B), (P): ',npara
!     Storage required for vectors
      write(*,'(7X,A,21X,A,X,f12.6,A)') 'Single-vector storage,',&
           '(V*R):',veclen*rF*GB,' GB'
      write(*,'(7X,A,20X,A,X,f12.6,A)') 'Block vector storage,',&
              '(B*V*R):',Blen*GB,' GB'
      write(*,'(7X,A,19X,A,X,f12.6,A)') 'Matrix-vector product,',&
           '(T*V*R):',mvlen*GB,' GB'

!     Memory cost for the solver of choice (parallelized over vectors)
      IF (cpp%npow.gt.0) THEN

!        Memory for the iteration
         IF (abs(styp).eq.1) THEN
            rG=htrm*rF
            IF (styp.eq.-1) THEN  ! power or intertwining
               IF (cpp%lowmem.eq.2 .or. cpp%lowmem.eq.3) THEN
                  mvGB=(Blen+npara*nmax)*GB               
                  write(*,'(7X,A,21X,A,X,f12.6,A)') &
                      'ALS-Power iteration,','(B+P*X):',mvGB,' GB'
               ELSE  ! lowmem=1,2 or invalid choice defaults here
                  mvGB=(Blen+npara*REAL(rG)*veclen)*GB
                  write(*,'(7X,A,15X,A,X,f12.6,A)') &
                      'ALS-Power iteration,','([B+P*T]*V*R):',mvGB,' GB'
               ENDIF
               redGB=getALSPOWmem(rG,rF,nmax,ndof,npara,cpp%lowmem)*GB
               write(*,'(7X,A,16X,A,X,f12.6,A)') &
                    'Extra memory for rank-reductions',':',redGB,' GB'
               mvGB=mvGB+redGB
               write(*,'(7X,A,23X,A,X,f12.6,A)') &
                    'ALS-Power iteration TOTAL',':',mvGB,' GB'

            ELSE
               mvGB=(Blen+npara*REAL(rG)*veclen)*GB
               write(*,'(7X,A,14X,A,X,f12.6,A)') &
                    'BlockPower iteration,','([B+P*T]*V*R):',mvGB,' GB'
               redGB=getRednMem(rG,rF,Q(1)%nbas,npara,useSVD)*GB
               write(*,'(7X,A,16X,A,X,f12.6,A)') &
                    'Extra memory for rank-reductions',':',redGB,' GB'
               mvGB=mvGB+redGB
               write(*,'(7X,A,22X,A,X,f12.6,A)') &
                    'BlockPower iteration TOTAL',':',mvGB,' GB'
            ENDIF
         ELSEIF (styp.eq.-2) THEN
            write(*,*) 'Davidson mem coming soon!!!'
         ELSEIF (styp.eq.-3) THEN
            write(*,*) 'MSBII mem coming soon!!!'
         ELSEIF (styp.eq.-4) THEN
            write(*,*) 'Inverse itn mem coming soon!!!'
         ELSE
            call AbortWithError('ShowPsiMem(): Solver not recognized')
         ENDIF
                 
!        Memory for Gram-Schmidt
         rG=nev*rF
         IF (styp.ge.0) THEN
            GSGB=2*Blen*GB
            write(*,'(7X,A,19X,A,X,f12.6,A)') &
                 'Gram-Schmidt storage,','(2B*V*R):',GSGB,' GB'
            redGB=getRednMem(rG,rF,Q(1)%nbas,1,useSVD)*GB
         ELSE  !!! Add styp=-2
            GSGB=Blen*GB
            write(*,'(7X,A,20X,A,X,f12.6,A)') &
                 'Gram-Schmidt storage,','(B*V*R):',GSGB,' GB'
            redGB=getALSOrthoMem(rG,rF,nmax,cpp%lowmem)*GB
         ENDIF
         write(*,'(7X,A,16X,A,X,f12.6,A)') &
              'Extra memory for rank-reductions',':',redGB,' GB'
         GSGB=GSGB+redGB
         write(*,'(7X,A,30X,A,X,f12.6,A)') &
              'Gram-Schmidt TOTAL',':',GSGB,' GB'

!        Memory for vector updates
         IF (cpp%update) THEN
            rG=nev*rF
            IF (styp.ge.0) THEN
               QHQGB=(Blen+npara*(mvlen+REAL(veclen)*rF))*GB
               write(*,'(7X,A,7X,A,X,f12.6,A)') &
                    'QHQ calculation storage,',&
                    '((B+P*(T+1))*V*R):',QHQGB,' GB'
               redGB=getRednMem(htrm*rF,rF,Q(1)%nbas,npara,useSVD)*GB
               write(*,'(7X,A,16X,A,X,f12.6,A)') &
                    'Extra memory for rank-reductions',':',redGB,' GB'
               QHQGB=QHQGB+redGB
               write(*,'(7X,A,27X,A,X,f12.6,A)') &
                    'QHQ calculation TOTAL',':',QHQGB,' GB'
               upGB=(2+npara)*Blen*GB
               write(*,'(7X,A,13X,A,X,f12.6,A)') &
                    'Vector update storage,',&
                    '([2+P]*B*V*R):',upGB,' GB'
               redGB=getRednMem(rG,rF,Q(1)%nbas,npara,useSVD)*GB
            ELSE !!! Add styp=-2,-3
               QHQGB=(Blen+npara*REAL(veclen)*rF)*GB
               write(*,'(7X,A,13X,A,X,f12.6,A)') &
                    'QHQ calculation storage,',&
                    '((B+P)*V*R):',QHQGB,' GB'
               redGB=getALSPOWmem(htrm*rF,rF,nmax,ndof,npara,&
                                  MAX(cpp%lowmem,2))*GB
               write(*,'(7X,A,16X,A,X,f12.6,A)') &
                    'Extra memory for rank-reductions',':',redGB,' GB'
               QHQGB=QHQGB+redGB
               write(*,'(7X,A,27X,A,X,f12.6,A)') &
                    'QHQ calculation TOTAL',':',QHQGB,' GB'
               upGB=2*Blen*GB
               write(*,'(7X,A,17X,A,X,f12.6,A)') &
                    'Vector update storage,',&
                    '(2*B*V*R):',upGB,' GB'
               redGB=getALSUpdateMem(rG,rF,nmax,npara,cpp%lowmem)*GB
            ENDIF
            write(*,'(7X,A,16X,A,X,f12.6,A)') &
                 'Extra memory for rank-reductions',':',redGB,' GB'
            upGB=upGB+redGB
            write(*,'(7X,A,29X,A,X,f12.6,A)') &
                 'Vector update TOTAL',':',upGB,' GB'
         ELSE
            QHQGB=(Blen+npara*REAL(veclen))*GB
            write(*,'(7X,A,13X,A,X,f12.6,A)') &
                    'QHQ calculation storage,',&
                    '(B*V*R+P*V):',QHQGB,' GB'
            upGB=(Blen+npara*REAL(veclen)*rF)*GB
            write(*,'(7X,A,14X,A,X,f12.6,A)') &
                 'Vector sorting storage,',&
                 '([B+P]*V*R):',upGB,' GB'
         ENDIF
      ENDIF

      write(*,'(7X,A,7X,A,X,f12.6,A/)') 'PSI MEMORY REQUIRED,',&
           '(MAX OF ABOVE TOTALS):',max(mvGB,GSGB,QHQGB,upGB),' GB'

      ENDIF rank0 

      end subroutine ShowPsiMem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DetermineDiag(Q,diag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines if diagonalization should be used

      implicit none
      TYPE (CP), INTENT(IN) :: Q(:)
      logical, intent(out) :: diag
      integer :: j,nev,ndof,prod

      nev=SIZE(Q)
      ndof=SIZE(Q(1)%nbas)

!     Diagonalize H directly if ALL states are requested
      diag=.FALSE.
      prod=1
      DO j=1,ndof
         prod=prod*Q(1)%nbas(j)
         IF (prod.gt.nev) EXIT
      ENDDO
      IF (prod.eq.nev) diag=.TRUE.

      end subroutine DetermineDiag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveHPsi(eigv,delta,cpp,Q,H,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for computing the eigenfunctions and 
! eigenvalues using the solver of choice

      implicit none
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CP), INTENT(IN)  :: H,W
      real*8, allocatable, intent(inout) :: eigv(:)
      real*8, allocatable, intent(inout) :: delta(:)
      real*8, allocatable  :: eigtmp(:),ccoef(:)
      real*8  :: bounds(2)
      integer :: i,j,nev,nup,ndown,nsame,nloc,ist,styp
      logical :: conv,showFmG,diag,readsuccess
      real*8  :: rmsdelta,oldrms,maxdelta,sumdelta
      real*8, parameter :: redtol=1.d-12

!     Initializations
      nev=SIZE(eigv)
      ist=0
      showFmG=.FALSE.
      oldrms=1.d99
      bounds=0.d0
      styp=GetSolverType(cpp,Q)

      call DetermineDiag(Q,diag)
      call ShowPsiMem(eigv,cpp,Q,H,styp)
      call ReadPsi(ist,bounds,eigv,delta,Q,cpp,readsuccess)

!     Calculate the spectral range of H
      IF (cpp%npow.gt.0 .and. cpp%ncycle.gt.0 &
          .and. (.not.diag) .and. (.not.readsuccess)) THEN
         call SetReductionParameters(min(50,cpp%psirank),cpp%psinals,&
                                     redtol,showFmG,'SVD','ALS')
         call GetSpectralRange(min(50,cpp%psirank),10,5,Q,H,bounds,&
                               cpp%lowmem)
      ENDIF

      call SetReductionParameters(cpp%psirank,cpp%psinals,redtol,&
                                  showFmG,cpp%red2D,cpp%redND)

!     Initial guess and pre-diagonalization
      IF (readsuccess) THEN
         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,*)
            write(*,*) 'Eigenvalues read: ',ist,(eigv(j),j=1,nev)
         ENDIF
      ELSE
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,*) 'Initial guess   : ',0,(eigv(j),j=1,nev)
         IF (cpp%ncycle.gt.0 .and. ((cpp%update.and.(styp.ne.-2)) &
             .or.diag)) THEN
            call Diagonalize(Q,H,eigv,.FALSE.,cpp%psinals,cpp%lowmem)
            IF (mpirank.eq.mpi_prnt_rank) THEN
               write(*,*)
               write(*,*) 'Diagonalization : ',0,(eigv(j),j=1,nev)
            ENDIF
         ENDIF
         call SavePsi(0,bounds,eigv,delta,Q,cpp)
      ENDIF

!     If all the eigenvalues were requested, exit here since the 
!     pre-diagonalization already gives the exact answer
      IF (diag) RETURN

!     If intertwining is used and the number of vectors in the block
!     ('nev') is smaller than psirank, the vectors will have a maximum 
!     rank of nev. Since intertwining does not change the rank of the 
!     input vectors, the step below ensures that we apply intertwining
!     to vectors with the correct rank
      IF (styp.lt.0 .and. styp.ne.-3) call AugmentQWithRandom(Q,cpp%psirank)

      ALLOCATE(eigtmp(nev))
      IF (ist.eq.0) delta=1.d99
      nloc=0

!     Loop over solver cycles
      DO i=1,cpp%ncycle

         IF (i.le.ist) CYCLE

!        Copy eigenvalues to temp array
         eigtmp=eigv

!        Run iterations using the solver of choice
         call Iterate(Q,H,W,eigtmp,eigv,cpp,bounds,nloc,i,styp)

         IF (styp.eq.-3) EXIT

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
         call SavePsi(i,bounds,eigv,delta,Q,cpp)

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

      end subroutine SolveHPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Iterate(Q,H,W,eigvo,eigv,cpp,bounds,nconv,i,styp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Performs iterations of various types, depending on the value of styp:

      implicit none
      TYPE (CPpar), INTENT(IN)    :: cpp
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), INTENT(IN)    :: H,W
      integer, intent(inout) :: nconv
      integer, intent(in)    :: i,styp
      real*8, intent(inout)  :: eigv(:)
      real*8, intent(in)     :: eigvo(:),bounds(2)
      real*8, parameter      :: tol=1.d-15
      character(len=18)      :: tag
      real*8  :: Eshift
      integer :: j,nbloc,sz,os

!     Easy exit for zero iterations
      IF (cpp%npow.lt.1) RETURN

!     Set parameters
      nbloc=SIZE(Q)

!     Solver specific procedures
      IF (styp.eq.1) THEN
         tag='BlockPower cycle: '
!        Power method: get estimate of optimal E-shift
         call GetBlockShift(eigv,bounds,Eshift)
      ELSEIF (styp.eq.-1) THEN
         call GetBlockShift(eigv,bounds,Eshift)
         tag='ALS-Power cycle : '
      ELSE
         call AbortWithError('Iterate(): invalid solver type')
      ENDIF

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc-nconv,sz,os)
      os=os+nconv

!     Run power iterations on each vector in the block
      DO j=1,sz
         IF (styp.eq.1) THEN
            call PowrRecurse(Q(os+j),H,cpp%npow,Eshift)
         ELSEIF (styp.eq.-1) THEN
            call ALS_POW_alg(H,Q(os+j),cpp%npow,1,Eshift,cpp%lowmem)
         ENDIF
      ENDDO

!     Sync vectors before Gram-Schmidt, updates
      call MPI_Sync_CP_block(Q)

!     Orthogonalization and update/vector sort
      IF (cpp%update) THEN
         IF (styp.eq.-1) THEN
            call ALS_ORTHO_alg(Q,cpp%psinals,cpp%lowmem)
            call Diagonalize(Q,H,eigv,.TRUE.,cpp%psinals,cpp%lowmem)
         ELSE
            call GRAMORTHO(Q)
            call Diagonalize(Q,H,eigv,.FALSE.,cpp%psinals,cpp%lowmem)
         ENDIF
      ELSE
         IF (styp.eq.-1) THEN
            call ALS_ORTHO_alg(Q,cpp%psinals,cpp%lowmem)
         ELSE
            call pGRAMORTHO(Q,nconv)
         ENDIF
         eigv(nconv+1:nbloc)=0.d0
         call GetQHQdiag(Q,H,eigv,nconv)
         call SortVecs(Q,eigv,nconv,bounds(1))
      ENDIF

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

      end subroutine Iterate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetSpectralRange(rk,npow,ncyc,Q,H,bounds,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Estimates the spectral range of the Hamiltonian using power method

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP) :: Q(:),T
      TYPE (CP), allocatable :: Qt(:)
      integer, intent(in) :: rk,npow,ncyc,lowmem
      real*8, intent(out) :: bounds(2)
      integer :: i,j,gst,ndof
      real*8  :: btmp(2)

      ndof=SIZE(Q(1)%nbas)
      bounds=0.d0

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(X,A/)') 'Power spectral range estimation'
         write(*,'(A,X,2(9X,A))') '  i','boundsl','boundsu'
      ENDIF

      ALLOCATE(Qt(2))

!     Guesses for upper and lower eigenvectors
      Qt(1)=ZeroCPvec(Q(1)%nbas)
      Qt(2)=ZeroCPvec(Q(1)%nbas)
      Qt(1)%coef(1)=1.d0
      Qt(2)%coef(1)=1.d0
      Qt(1)%base=1.d-15
      Qt(2)%base=1.d-15
      gst=0
      DO i=1,ndof
         Qt(1)%base(gst+1,1)=1.d0
         gst=gst+Qt(1)%nbas(i)
         Qt(2)%base(gst,1)=1.d0
      ENDDO

!     Since the ALS-guided power method (used for > 2 DOF) does not
!     change the rank, each vector must be initiated with rank rk, so
!     fill Qt(1) and Qt(2) up to rank rk with random terms
      IF (SIZE(Q(1)%nbas).gt.2 .and. rk.gt.1) THEN
         T=RandomCP(Q(1),rk-1)
         call NORMBASE(T)
         T%coef=1.d-7
         call SUMVECVEC(Qt(1),1.d0,T,1.d0)
         call SUMVECVEC(Qt(2),1.d0,T,1.d0)
         call FlushCP(T)
      ENDIF

!     Initial bounds
!$omp parallel
!$omp do private(j) schedule(static)
      DO j=1,2
         bounds(j)=RayleighQuotient2(Qt(j),H)
      ENDDO
!$omp end do
!$omp end parallel

!     btmp is the shift for the power method
      btmp(1)=0.5*(bounds(1)+bounds(2))
      btmp(2)=0.d0

      IF (mpirank.eq.mpi_prnt_rank) &
         write(*,'(i3,2(X,2f16.6))') 0,bounds(1),bounds(2)

!     Run the power method to improve the bounds
      DO i=1,ncyc
!$omp parallel
!$omp do private(j) schedule(static)
         DO j=1,2
            IF (SIZE(Q(1)%nbas).eq.2) THEN  ! 2 DOF: use power itn + SVD
               call PowrRecurse(Qt(j),H,npow,btmp(j))
            ELSE  ! >2 DOF: use ALS-guided power method
               call ALS_POW_alg(H,Qt(j),npow,1,btmp(j),lowmem)
            ENDIF
            bounds(j)=RayleighQuotient2(Qt(j),H)
         ENDDO
!$omp end do
!$omp end parallel
!        Update the shift for the ground state
         btmp(1)=0.5*(bounds(1)+bounds(2))
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,'(i3,2(X,2f16.6))') i,bounds(1),bounds(2)
      ENDDO

      DEALLOCATE(Qt)

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A,2(f15.6,A))') 'Spectral range of H = [',&
                                   bounds(1),',',bounds(2),']'
      btmp(1)=0.001*(bounds(2)-bounds(1))
      bounds(1)=bounds(1)-btmp(1)
      bounds(2)=bounds(2)+btmp(1)

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A,2(f15.6,A)/)') 'Range padded by .1% = [',&
                                   bounds(1),',',bounds(2),']'

      end subroutine GetSpectralRange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE SOLVER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
