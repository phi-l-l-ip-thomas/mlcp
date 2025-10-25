!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ANALYZER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MODECOMB
      USE SEPDREPN
      USE HAMILSETUP
      USE REDUCTION
      USE MODVECVEC
      USE MUNKRES
      USE CPCONFIG

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_Analyzer_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_Analyzer_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_Analyzer_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_Analyzer_Module()
      call Get_MPI_Timings('Analyzer module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_Analyzer_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzePsi(im,eigv,delta,Q,Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in) :: im
      real(kind=8), intent(in)  :: eigv(:),delta(:)
      real(kind=8) :: ti1,ti2
      integer, allocatable :: qns(:,:)

      IF (.NOT. MODULE_SETUP) call Init_Analyzer_Module()

      call CPU_TIME(ti1)

      call AnalyzeConfigs(Q,im,eigv,Ham)
!!! EXPERIMENTAL: use with caution
!      call AssignConfigs(Q,qns)
!      call AssignConfigsPlus(Q,qns,im,Ham,eigv)
!!! END EXPERIMENTAL
      call AnalyzeRank1(Q,qns)
      call SetEigenbasis(Ham%nt,im,qns,eigv,delta)
      deallocate(qns)

      IF (Ham%nt(im)%nHterm().gt.0) THEN
         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(/X,2A/)') 'Eigenvectors, assignments based on ',&
                                'rank-1 approximation :'
            call Ham%nt(im)%showeigen()
         ENDIF
      ENDIF

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine AnalyzePsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzeRank1(Q,qns)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes leading configuration of the rank-1 approximation

      implicit none
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP) :: v
      integer, allocatable, intent(out) :: qns(:,:)
      character(len=64), parameter :: solver='LU'
      real*8, parameter   :: redtol=1.d-12
      integer :: i,j,k,ndof,nev,maxind,gst
      real*8  :: maxcoef

!     Set parameters
      nev=SIZE(Q)
      ndof=SIZE(Q(1)%nbas)

!     Extract the assignment from a rank-1 approximation of each
!     eigenfunction (if > 2 DOFs, use SR1 with many steps since
!     SR1 exits if the coef converges; note that als_penalty and
!     als_solver are irrelevant since SR1 does not use them.)
      call SetReductionParameters(1,100,redtol,.FALSE.,'SVD','SR1',&
                                  1.d-10,solver)

!     Loop over eigenstates and reduce each vector to rank-1. Then
!     extract the index of the most important coefficient for 
!     each sub-mode basis function
      ALLOCATE(qns(nev,ndof))
      DO i=1,nev
         v=CopyCP(Q(i))
         call reduc(v)
         call NORMALIZE(v)

         gst=0
         DO j=1,ndof
            maxind=1
            maxcoef=abs(v%base(gst+1,1))
            DO k=2,Q(i)%nbas(j)
               IF (abs(v%base(gst+k,1)).gt.maxcoef) THEN
                  maxind=k
                  maxcoef=abs(v%base(gst+k,1))
               ENDIF
            ENDDO
            qns(i,j)=maxind
            gst=gst+v%nbas(j)
         ENDDO
         call FlushCP(v)
      ENDDO

      end subroutine AnalyzeRank1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzeConfigs(Q,im,eigv,Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes/prints dominant product configurations of the wavefunction

      implicit none
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs) :: v
      real(kind=8), intent(in)   :: eigv(:)
      integer, intent(in)  :: im
      integer, allocatable :: qns(:)
      integer :: i,j,k,nev,nsubm,ndof,ncoef
      character*64 :: frmt
      integer, parameter :: ncoefmax=16
      real(kind=8), parameter :: printtol=5.d-2

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)
      ndof=Ham%nt(im)%ndof()

      IF (nsubm.lt.2) RETURN

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(/X,A,ES11.4,A,I0,A/)') &
         'Eigenvectors: configurations with |c|^2 larger than : ',&
         printtol,' x largest coef, or largest ',ncoefmax,' coefs'
         write(frmt,*) '(A,X,',ndof,'(I2,X))'
         write(*,frmt) 'Mode:',(Ham%nt(im)%dofs(j),j=1,ndof)
      ENDIF

      DO i=1,nev
!        Get the list of dominant configurations for each eigenvalue
         call GetConfigList(Q(i),100,v)

!        Get and print the full assignment of the largest coefficient
         call ConstructAssignment(Ham%nt,im,v%qns(1,:),qns)
         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(frmt,*) '(I4,A,',3*ndof+1,'X,2(f19.12,X))'
            write(*,frmt) i,')',eigv(i),eigv(i)-eigv(1)
            write(frmt,'(A,I0,A)') '(6X,',ndof,'(I2,X),4X,f6.3)'
            write(*,frmt) (qns(k),k=1,ndof),v%coef(1)**2
         ENDIF
         DEALLOCATE(qns)

!        Print other configurations if the coefficients are large enough
         ncoef=1
         DO j=2,SIZE(v%coef)
            IF (v%coef(j)**2.gt.printtol*v%coef(1)**2) THEN
               ncoef=ncoef+1
               call ConstructAssignment(Ham%nt,im,v%qns(j,:),qns)
               IF (mpirank.eq.mpi_prnt_rank) &
                  write(*,frmt) (qns(k),k=1,ndof),v%coef(j)**2
               DEALLOCATE(qns)
               IF (ncoef.eq.ncoefmax) EXIT
            ENDIF
         ENDDO
         call FlushConfigs(v)
      ENDDO

      end subroutine AnalyzeConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AssignConfigs(Q,qns)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assign states via dominant configurations, using Munkres algorithm

      implicit none
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (Configs) :: v,w
      integer, allocatable, intent(out) :: qns(:,:)
      integer, allocatable :: avec(:)
      real*8, allocatable  :: weights(:,:),wtmp(:,:)
      integer :: i,j,k,nev,nsubm,nfound,assgn
      integer, parameter :: cmax=100
      logical :: found,success

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

!!!  TO DO: change below to avoid unnecessarily large arrays
      ALLOCATE(weights(nev,nev*cmax))
      weights(:,:)=0.d0
      call NewConfigs(w,Q(1)%nbas,nev*cmax)

      nfound=0
      DO i=1,nev
!        Get the list of dominant configurations for this eigenvalue
         call GetConfigList(Q(i),cmax,v)
!!! NORMALIZATION OF v (might be needed if some configs are large) ???

         DO j=1,SIZE(v%coef)

            found=.FALSE.
            DO k=1,nfound
!              Config found in master list
               IF (ALL(v%qns(j,:).eq.w%qns(k,:))) THEN
                  weights(i,k)=v%coef(j)**2
                  found=.TRUE.
                  EXIT
               ENDIF
            ENDDO

!           Config not found: extend the master list
            IF (.not.found) THEN
               nfound=nfound+1
               w%qns(nfound,:)=v%qns(j,:)
               w%coef(nfound)=1.d0 ! Set coef to 1 for trimming later
               weights(i,nfound)=v%coef(j)**2
            ENDIF

         ENDDO
         call FlushConfigs(v)
      ENDDO

!     Resize configuration list and weights matrix and assign states
      call ResizeConfigList(w,nfound)
      ALLOCATE(wtmp(nev,nfound))
      wtmp(:,:)=weights(:,:nfound)
      DEALLOCATE(weights)
      weights=AssignMatrix(wtmp) ! <--Munkres called here
      DEALLOCATE(wtmp)

!     Extract assignments and store in qns array
      ALLOCATE(qns(nev,nsubm))
      avec=GetMunkresAssignVec(weights)
      DO i=1,nev
         IF (avec(i).eq.0) &
            call AbortWithError('AssignConfigs(): j not found')
         qns(i,:)=w%qns(avec(i),:)
      ENDDO

      DEALLOCATE(weights,avec)
      call FlushConfigs(w)

      end subroutine AssignConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AssignConfigsPlus(Q,qns,im,Ham,eigv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assign states via dominant configurations, using Munkres algorithm

      implicit none
      TYPE (CP), INTENT(IN)     :: Q(:)
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs) :: v,w
      real*8, intent(in)   :: eigv(:)
      integer, intent(in)  :: im
      integer, allocatable, intent(out) :: qns(:,:)
      integer, allocatable :: avec(:),avec2(:),qtmp(:),ntmp(:)
      real*8, allocatable  :: weights(:,:),wtmp(:,:)
      integer :: i,j,k,nev,nsubm,nfound,assgn,ndof
      integer, parameter :: cmax=100
      logical :: found,success
      real*8  :: t1,t2
      character*72 :: frmt

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

!!!  TO DO: change below to avoid unnecessarily large arrays
      ALLOCATE(weights(nev,nev*cmax))
      weights(:,:)=0.d0
      call NewConfigs(w,Q(1)%nbas,nev*cmax)

      nfound=0
      DO i=1,nev
!        Get the list of dominant configurations for this eigenvalue
         call GetConfigList(Q(i),cmax,v)
!!! NORMALIZATION OF v (might be needed if some configs are large) ???

         DO j=1,SIZE(v%coef)

            found=.FALSE.
            DO k=1,nfound
!              Config found in master list
               IF (ALL(v%qns(j,:).eq.w%qns(k,:))) THEN
                  weights(i,k)=v%coef(j)**2
                  found=.TRUE.
                  EXIT
               ENDIF
            ENDDO

!           Config not found: extend the master list
            IF (.not.found) THEN
               nfound=nfound+1
               w%qns(nfound,:)=v%qns(j,:)
               w%coef(nfound)=1.d0 ! Set coef to 1 for trimming later
               weights(i,nfound)=v%coef(j)**2
            ENDIF

         ENDDO
         call FlushConfigs(v)
      ENDDO

!     Resize configuration list and weights matrix and assign states
      call ResizeConfigList(w,nfound)
      ALLOCATE(wtmp(nev,nfound))
      wtmp(:,:)=weights(:,:nfound)
      DEALLOCATE(weights)
      weights=AssignMatrix(wtmp) ! <--Munkres called here

!     Extract assignments and store in qns array
      ALLOCATE(qns(nev,nsubm))
      avec=GetMunkresAssignVec(weights)

!     Refine assignments by energy
      IF (nsubm.gt.1) THEN
              
!        Instead of w, need configs corresponding to 1D functions from
!        ConstructAssignment(). Call using the 1st one to get correct
!        width, then fill the rest
         ndof=Ham%nt(im)%ndof()
         call ConstructAssignment(Ham%nt,im,w%qns(1,:),qtmp)
         ALLOCATE(ntmp(SIZE(qtmp)))
         ntmp(:)=16384
         call NewConfigs(v,ntmp,nfound)
         v%coef(:)=1.d0
         v%qns(1,:)=qtmp(:)
         DEALLOCATE(ntmp,qtmp)
         DO i=2,nfound
            call ConstructAssignment(Ham%nt,im,w%qns(i,:),qtmp)
            v%qns(i,:)=qtmp(:)
            DEALLOCATE(qtmp)
         ENDDO

!        Refine assignments using energies
         call RefineByEnergy(v,eigv,wtmp,weights,avec2,success)

         IF (success) THEN
            IF (.not.ALL(avec(:).eq.avec2(:))) THEN
               write(*,'(/X,A/)') 'States reassigned using energies:'
               write(frmt,*) '(A,X,',ndof,'(I2,X),3X,',ndof,&
                             '(I2,X),5X,A)'
               write(*,frmt) 'Mode:',(Ham%nt(im)%dofs(j),j=1,ndof),&
                             (Ham%nt(im)%dofs(j),j=1,ndof),'Energy'
               write(frmt,*) '(I4,A,X,',ndof,'(I2,X),A,X,',&
                             ndof,'(I2,X),f19.12)'
               DO i=1,nev
                  IF (avec2(i).ne.avec(i)) &
                     write(*,frmt) i,')',&
                     (v%qns(avec(i),j),j=1,ndof),'->',&
                     (v%qns(avec2(i),j),j=1,ndof),eigv(i)
               ENDDO
               avec(:)=avec2(:)
            ENDIF
         ENDIF

         DEALLOCATE(avec2)
         call FlushConfigs(v)

      ENDIF ! Assign-by-energy

      DO i=1,nev
         IF (avec(i).eq.0) THEN
            write(*,*) 'No assignment found for state ',i
            call AbortWithError('AssignConfigs(): config not found')
         ENDIF
         qns(i,:)=w%qns(avec(i),:)
      ENDDO

      DEALLOCATE(weights,wtmp,avec)
      call FlushConfigs(w)

      end subroutine AssignConfigsPlus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RefineByEnergy(v,eigv,wts,Mwts,avec,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Refine Munkres assignments using sums of energies

      implicit none
      TYPE (Configs), intent(in) :: v
      TYPE (Configs) :: w
      real*8, intent(in) :: eigv(:)
      real*8, allocatable :: eigc(:)
      real*8, allocatable, intent(inout) :: wts(:,:),Mwts(:,:)
      logical, intent(out) :: success
      logical :: found
      integer, allocatable, intent(out) :: avec(:)
      integer, allocatable :: sumcfg(:),rowidx(:),colidx(:)
      integer :: i,j,k,l,nev,nconfig,maxsumv,nsubm,nassng,tmp
      real*8  :: tmpE
      character*72 :: frmt

      nev=SIZE(eigv)
      nconfig=SIZE(v%coef)
      nsubm=SIZE(v%nbas)
      success=.FALSE.

!     Create copies of energies, configs, also indices (for permuting)
      call CopyConfigsWtoV(w,v)
      ALLOCATE(avec(nev),eigc(nev),rowidx(nev),colidx(nconfig))
      avec(:)=0
      eigc(:)=eigv(:)
      DO i=1,nev
         rowidx(i)=i
      ENDDO
      DO i=1,nconfig
         colidx(i)=i
      ENDDO

!     Get array of qn sums and max-sum-of-qns
      ALLOCATE(sumcfg(nconfig))
      maxsumv=0
      DO i=1,nconfig
         sumcfg(i)=-nsubm !!! so that sumcfg has min value of zero
         DO j=1,nsubm
            sumcfg(i)=sumcfg(i)+w%qns(i,j)
         ENDDO
         IF (sumcfg(i).gt.maxsumv) maxsumv=sumcfg(i)
      ENDDO

      nassng=0

!     Assign the ground state using input Munkres result
      found=.FALSE.
      DO j=1,nconfig
         IF (sumcfg(j).eq.0) THEN
!           Config found: make sure it is assigned to a state
            DO i=1,nev
               IF (Mwts(i,j).eq.1.d0) THEN
!                 Permute the energy and config lists, the qn sum array,
!                 and the weight/assignment matrices to place this
!                 config first
                  nassng=nassng+1
                  call PermuteArrays(nassng,i,j,w,eigc,wts,Mwts,&
                                     sumcfg,rowidx,colidx)
                  found=.TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (found) EXIT
         ENDIF
      ENDDO

!     If ground state is not found, exit without using the energy
!     prediction
      IF (.not.found) THEN
         DEALLOCATE(eigc,sumcfg,rowidx,colidx)
         call FlushConfigs(w)
         RETURN
      ENDIF

!     Assign singly excited states using input Munkres result    
      DO j=nassng+1,nconfig
         IF (sumcfg(j).eq.1) THEN
!           Config found: make sure it is assigned to a state
            DO i=nassng+1,nev
               IF (Mwts(i,j).eq.1.d0) THEN
!                 Permute arrays to place singly excited states after
!                 the ground state
                  nassng=nassng+1
                  call PermuteArrays(nassng,i,j,w,eigc,wts,Mwts,&
                                     sumcfg,rowidx,colidx)
                  EXIT
               ENDIF
               IF (nassng.eq.nev) EXIT
            ENDDO
            IF (nassng.eq.nev) EXIT
         ENDIF
      ENDDO


!     Main loop over sum v_i
      DO k=2,maxsumv

!        Generate guess weight matrix from earlier assignments
!        and use Munkres to get assignments for this k
!        Also replace existing block of Mwts with new one
         call EnergyBasedAssignments(nassng,w,eigc,wts,Mwts,sumcfg)

!        Assign the next group of found states    
         DO j=nassng+1,nconfig
            IF (sumcfg(j).le.k) THEN
!              Config found: make sure it is assigned to a state
               DO i=nassng+1,nev
                  IF (Mwts(i,j).eq.1.d0) THEN
!                    Permute arrays to place up-to-k-excited states next
                     nassng=nassng+1
                     call PermuteArrays(nassng,i,j,w,eigc,wts,Mwts,&
                                        sumcfg,rowidx,colidx)
                     EXIT
                  ENDIF
                  IF (nassng.eq.nev) EXIT
               ENDDO ! i loop over energies
               IF (nassng.eq.nev) EXIT
            ENDIF
         ENDDO ! j loop over configs
         IF (nassng.eq.nev) EXIT
      ENDDO ! k loop over sum v_i

!     Fill assignment vector
      IF (nassng.eq.nev) THEN
         DO i=1,nev
            avec(rowidx(i))=colidx(i)
         ENDDO
         success=.TRUE.
      ENDIF

      DEALLOCATE(eigc,sumcfg,rowidx,colidx)
      call FlushConfigs(w)

      end subroutine RefineByEnergy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine EnergyBasedAssignments(n,w,eig,wts,Mwts,sumcfg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates assignments for higher excited states from guesses based on
! lower excited states and the ground state

      implicit none
      TYPE (Configs), intent(inout) :: w
      integer, intent(in) :: sumcfg(:)
      real*8, intent(in) :: eig(:)
      real*8, intent(in) :: wts(:,:)
      real*8, intent(inout) :: Mwts(:,:)
      real*8, allocatable :: wtstmp(:,:),Mwtstmp(:,:)
      integer, intent(in) :: n
      integer :: i,j,r,c,d,redr,redc
      real*8, parameter :: fac=2.d0 ! Gaussian dropoff from Eguess
      real*8, parameter :: mix=0.25d0 ! portion of energy weight
      real*8  :: gap
      character*72 :: frmt

      r=SIZE(eig)
      c=SIZE(w%qns,1)
      d=SIZE(w%qns,2)
      redr=r-n
      redc=c-n

      ALLOCATE(wtstmp(redr,redc))
      wtstmp(:,:)=0.d0

      call GuessEnergies(n,w)

!     Compute the guess-derived weights
      DO i=n+1,r
         DO j=n+1,c
!           Gap estimate: (E_state - E_gs)/(quanta_in_config)
            gap=abs(eig(i)-eig(1))/REAL(sumcfg(j))
            wtstmp(i-n,j-n)=exp(-fac*abs(eig(i)-w%coef(j))/gap)
         ENDDO
      ENDDO

!     Add wavefunction-derived weights to energy-derived ones
      wtstmp(:,:)=mix*wtstmp(:,:)+(1.d0-mix)*wts(n+1:,n+1:)

!     Assign states for this block and overwrite previous assignments
      Mwtstmp=AssignMatrix(wtstmp)
      Mwts(n+1:,n+1:)=Mwtstmp(:,:)

      DEALLOCATE(wtstmp,Mwtstmp)

      end subroutine EnergyBasedAssignments

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GuessEnergies(n,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Guess energies for configs in w after the n-th, where the coefs of w
! hold the energies used for making the guess

      implicit none
      TYPE (Configs), intent(inout) :: w
      integer, intent(in)  :: n
      integer, allocatable :: cfg1(:),cfg2(:)
      integer :: i,j,k,c,d,nguess,ic1,ic2
      real*8  :: guess

      c=SIZE(w%qns,1)
      d=SIZE(w%qns,2)

      ALLOCATE(cfg1(d),cfg2(d))

!     Loop over configs to have their energies guessed
      DO i=n+1,c

         nguess=0
         guess=0.d0
         w%coef(i)=1.d99

!        First (preferred) guess method, average the extrapolated
!        guesses along each mode
         cfg2(:)=w%qns(i,:)
         cfg1(:)=w%qns(i,:)
         DO j=1,d
!           Second-previous configuration too low: skip this iteration
            IF (cfg2(j)-2.lt.1) CYCLE

!           Generate decremented configurations for guess
            cfg2(j)=cfg2(j)-2
            cfg1(j)=cfg1(j)-1

!           If both decremented configs are found, make the guess
            ic2=0
            DO k=1,n
               IF (ALL(w%qns(k,:).eq.cfg2(:))) THEN
                  ic2=k
                  EXIT
               ENDIF
            ENDDO
            IF (ic2.eq.0) THEN ! Config not found
               cfg2(j)=cfg2(j)+2
               cfg1(j)=cfg1(j)+1
               CYCLE
            ENDIF

            ic1=0
            DO k=1,n
               IF (ALL(w%qns(k,:).eq.cfg1(:))) THEN
                  ic1=k
                  EXIT
               ENDIF
            ENDDO
            IF (ic1.eq.0) THEN ! Config not found
               cfg2(j)=cfg2(j)+2
               cfg1(j)=cfg1(j)+1
               CYCLE
            ENDIF

!           Make the guess if extrapolation gives an energy increase
            IF (w%coef(ic1).gt.w%coef(ic2)) THEN
               guess=guess+2*w%coef(ic1)-w%coef(ic2)
               nguess=nguess+1
            ENDIF

!           Restore decremented configurations to full
            cfg2(j)=cfg2(j)+2
            cfg1(j)=cfg1(j)+1
         ENDDO

         IF (nguess.gt.0) THEN ! First method succeeded
            w%coef(i)=guess/REAL(nguess)
         ELSE ! Try second method
!           Second guess method: average sums of 1-mode, (d-1)-mode
!           energies.
            DO j=1,d
!              Generate decremented configurations for guess
               cfg2(:)=w%qns(i,:)
               cfg2(j)=1
               cfg1(:)=1
               cfg1(j)=w%qns(i,j)

!              If both decremented configs are found, make the guess
               ic2=0
               DO k=1,n
                  IF (ALL(w%qns(k,:).eq.cfg2(:))) THEN
                     ic2=k
                     EXIT
                  ENDIF
               ENDDO
               IF (ic2.eq.0) CYCLE

               ic1=0
               DO k=1,n
                  IF (ALL(w%qns(k,:).eq.cfg1(:))) THEN
                     ic1=k
                     EXIT
                  ENDIF
               ENDDO
               IF (ic1.eq.0) CYCLE

!              Make the guess
               guess=guess+w%coef(ic1)+w%coef(ic2)-w%coef(1)
               nguess=nguess+1
            ENDDO
            IF (nguess.gt.0) THEN ! Second method succeeded
               w%coef(i)=guess/REAL(nguess)
            ENDIF
         ENDIF
      ENDDO

      DEALLOCATE(cfg1,cfg2)

      end subroutine GuessEnergies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PermuteArrays(n,i,j,w,eig,M1,M2,scfg,ridx,cidx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Helper routine for rearranging arrays

      implicit none
      TYPE (Configs), intent(inout) :: w
      real*8, intent(inout)  :: eig(:)
      real*8, intent(inout)  :: M1(:,:),M2(:,:)
      integer, intent(inout) :: scfg(:),ridx(:),cidx(:)
      integer, intent(in)  :: n,i,j
      integer, allocatable :: tmpcfg(:)
      real*8, allocatable  :: rowperm(:),colperm(:)
      integer :: rows,cols,ndof,tmp,k
      real*8  :: tmpE

      rows=SIZE(eig)
      cols=SIZE(w%qns,1)
      ndof=SIZE(w%qns,2)

      ALLOCATE(rowperm(cols),colperm(rows),tmpcfg(ndof))
      rowperm(:)=0.d0

      tmpE=eig(n)
      eig(n)=eig(i)
      eig(i)=tmpE

      tmpcfg(:)=w%qns(n,:)
      w%qns(n,:)=w%qns(j,:)
      w%qns(j,:)=tmpcfg(:)

!     Copy the guess energy to coef of w as it is used later
      w%coef(n)=eig(n)

      tmp=scfg(n)
      scfg(n)=scfg(j)
      scfg(j)=tmp

      tmp=ridx(n)
      ridx(n)=ridx(i)
      ridx(i)=tmp

      tmp=cidx(n)
      cidx(n)=cidx(j)
      cidx(j)=tmp

      rowperm(:)=M1(n,:)
      M1(n,:)=M1(i,:)
      M1(i,:)=rowperm(:)

      colperm(:)=M1(:,n)
      M1(:,n)=M1(:,j)
      M1(:,j)=colperm(:)

      rowperm(:)=M2(n,:)
      M2(n,:)=M2(i,:)
      M2(i,:)=rowperm(:)

      colperm(:)=M2(:,n)
      M2(:,n)=M2(:,j)
      M2(:,j)=colperm(:)

      DEALLOCATE(rowperm,colperm,tmpcfg)

      end subroutine PermuteArrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ANALYZER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
