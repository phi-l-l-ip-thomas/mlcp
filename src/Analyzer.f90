!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ANALYZER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE MODECOMB
      USE SEPDREPN
      USE HAMILSETUP
      USE REDUCTION
      USE MODVECVEC
      USE MUNKRES
      USE CPCONFIG

      implicit none
      real*8, private  :: anal_time=0.d0
      logical, private :: ANAL_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeAnalModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      anal_time = 0.d0
      ANAL_SETUP = .TRUE.

      end subroutine InitializeAnalModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeAnalModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. ANAL_SETUP) call InitializeAnalModule()

      ANAL_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total wave-function analysis time (s)',&
                             anal_time

      end subroutine DisposeAnalModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzePsi(il,im,eigv,delta,Q,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in) :: il,im
      real*8, intent(in)  :: eigv(:),delta(:)

      IF (.NOT. ANAL_SETUP) call InitializeAnalModule()

      call AnalyzeConfigs(Q,il,im,eigv,Ham,ML)
!!! EXPERIMENTAL: use with caution
!      call AssignConfigs(Q,Ham%eig(il,im)%assgn)
!      call AssignConfigsPlus(Q,Ham%eig(il,im)%assgn,il,im,Ham,eigv,ML)
!!! END EXPERIMENTAL
      call AnalyzeRank1(Q,Ham%eig(il,im)%assgn)
      call PrintAssignments(il,im,eigv,delta,Ham,ML)

      end subroutine AnalyzePsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzeRank1(Q,qns)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes leading configuration of the rank-1 approximation

      implicit none
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP) :: v
      integer, allocatable, intent(out) :: qns(:,:)
      real*8, parameter   :: redtol=1.d-12
      integer :: i,j,k,ndof,nev,maxind,gst
      real*8  :: maxcoef,t1,t2

      call CPU_TIME(t1)

!     Set parameters
      nev=SIZE(Q)
      ndof=SIZE(Q(1)%nbas)

!     Extract the assignment from a rank-1 approximation of each
!     eigenfunction (if > 2 DOFs, use SR1 with many steps since
!     SR1 exits if the coef converges)
      call SetReductionParameters(1,100,redtol,.FALSE.,'SVD','SR1')

!     Assignment matrix
      ALLOCATE(qns(nev,ndof))

!     Loop over eigenstates and reduce each vector to rank-1. Then
!     extract the index of the most important coefficient for 
!     each sub-mode basis function

      DO i=1,nev
         v=CopyCP(Q(i))
         call CPU_TIME(t2)
         anal_time=anal_time+t2-t1
         call reduc(v)
         call NORMALIZE(v)
         call CPU_TIME(t1)

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

      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1

      end subroutine AnalyzeRank1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AnalyzeConfigs(Q,il,im,eigv,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes/prints dominant product configurations of the wavefunction

      implicit none
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs) :: v
      real*8, intent(in)   :: eigv(:)
      integer, intent(in)  :: il,im
      integer, allocatable :: qns(:)
      integer :: i,j,k,nev,nsubm,nagn,mst,mfi,ncoef
      character*64 :: frmt
      integer, parameter :: ncoefmax=16
      real*8, parameter  :: printtol=5.d-2
      real*8  :: t1,t2

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

      IF (nsubm.lt.2) RETURN

      call CPU_TIME(t1)

      write(*,'(/X,A,ES11.4,A,I0,A/)') &
            'Eigenvectors: configurations with |c|^2 larger than : ',&
             printtol,' x largest coef, or largest ',ncoefmax,' coefs'

!     Print mode numbers
      mst=firstmode(il,1,im,ML)
      mfi=lastmode(il,1,im,ML)
      nagn=mfi-mst+1
      write(frmt,*) '(A,X,',nagn,'(I2,X))'
      write(*,frmt) 'Mode:',(ML%resort(j),j=mst,mfi)

      DO i=1,nev
!        Get the list of dominant configurations for each eigenvalue
         call GetConfigList(Q(i),100,v)

!        Get and print the full assignment of the largest coefficient
         call GetFullAssignment(il,im,Ham,ML,v%qns(1,:),qns)
         write(frmt,*) '(I4,A,',3*nagn+1,'X,2(f19.12,X))'
         write(*,frmt) i,')',eigv(i),eigv(i)-eigv(1)
         write(frmt,'(A,I0,A)') '(6X,',nagn,'(I2,X),4X,f6.3)'
         write(*,frmt) (qns(k)-1,k=1,nagn),v%coef(1)**2
         DEALLOCATE(qns)

!        Print other configurations if the coefficients are large enough
         ncoef=1
         DO j=2,SIZE(v%coef)
            IF (v%coef(j)**2.gt.printtol*v%coef(1)**2) THEN
               ncoef=ncoef+1
               call GetFullAssignment(il,im,Ham,ML,v%qns(j,:),qns)
               write(*,frmt) (qns(k)-1,k=1,nagn),v%coef(j)**2
               IF (ncoef.eq.ncoefmax) EXIT
            ENDIF
         ENDDO
         call FlushConfigs(v)
      ENDDO

      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1

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
      real*8  :: t1,t2

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

      call CPU_TIME(t1)

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
      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1
      weights=AssignMatrix(wtmp) ! <--Munkres called here
      call CPU_TIME(t1)
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

      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1

      end subroutine AssignConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AssignConfigsPlus(Q,qns,il,im,Ham,eigv,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assign states via dominant configurations, using Munkres algorithm

      implicit none
      TYPE (CP), INTENT(IN)     :: Q(:)
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs) :: v,w
      real*8, intent(in)   :: eigv(:)
      integer, intent(in)  :: il,im
      integer, allocatable, intent(out) :: qns(:,:)
      integer, allocatable :: avec(:),avec2(:),qtmp(:),ntmp(:)
      real*8, allocatable  :: weights(:,:),wtmp(:,:)
      integer :: i,j,k,nev,nsubm,vsubm,nfound,assgn,mst,mfi
      integer, parameter :: cmax=100
      logical :: found,success
      real*8  :: t1,t2
      character*72 :: frmt

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

      call CPU_TIME(t1)

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
      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1
      weights=AssignMatrix(wtmp) ! <--Munkres called here
      call CPU_TIME(t1)

!     Extract assignments and store in qns array
      ALLOCATE(qns(nev,nsubm))
      avec=GetMunkresAssignVec(weights)

!     Refine assignments by energy
      IF (nsubm.gt.1) THEN
              
!        Instead of w, need configs corresponding to 1D functions from
!        GetFullAssignment(). Call using the 1st one to get correct
!        width, then fill the rest
         call GetFullAssignment(il,im,Ham,ML,w%qns(1,:),qtmp)
         ALLOCATE(ntmp(SIZE(qtmp)))
         ntmp(:)=16384
         call NewConfigs(v,ntmp,nfound)
         v%coef(:)=1.d0
         v%qns(1,:)=qtmp(:)
         DEALLOCATE(ntmp,qtmp)
         DO i=2,nfound
            call GetFullAssignment(il,im,Ham,ML,w%qns(i,:),qtmp)
            v%qns(i,:)=qtmp(:)
            DEALLOCATE(qtmp)
         ENDDO

!        Refine assignments using energies
         call RefineByEnergy(v,eigv,wtmp,weights,avec2,success)

         IF (success) THEN
            IF (.not.ALL(avec(:).eq.avec2(:))) THEN
               write(*,'(/X,A/)') 'States reassigned using energies:'
               mst=firstmode(il,1,im,ML)
               mfi=lastmode(il,1,im,ML)
               vsubm=SIZE(v%qns,2)

               write(frmt,*) '(A,X,',vsubm,'(I2,X),3X,',vsubm,&
                             '(I2,X),5X,A)'
               write(*,frmt) 'Mode:',(ML%resort(j),j=mst,mfi),&
                             (ML%resort(j),j=mst,mfi),'Energy'

               write(frmt,*) '(I4,A,X,',vsubm,'(I2,X),A,X,',&
                             vsubm,'(I2,X),f19.12)'
               DO i=1,nev
                  IF (avec2(i).ne.avec(i)) &
                     write(*,frmt) i,')',&
                     (v%qns(avec(i),j)-1,j=1,vsubm),'->',&
                     (v%qns(avec2(i),j)-1,j=1,vsubm),eigv(i)
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

      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1

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

      subroutine PrintAssignments(il,im,eigv,delta,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in)  :: il,im
      real*8, intent(in)   :: eigv(:),delta(:)
      integer, allocatable :: qns(:)
      integer :: i,j,ndof,nev,nagn,sp
      integer :: iexc,jexc,nexc,mst,mfi
      character*72 :: frmt
      character*8  :: labl
      real*8 :: t1,t2

!     Set parameters
      nev=SIZE(Ham%eig(il,im)%assgn,1)
      ndof=SIZE(Ham%eig(il,im)%assgn,2)

!     Print output only if > 1 submode
      IF (ndof.lt.2) RETURN

      call CPU_TIME(t1)

      write(*,'(/X,2A/)') 'Eigenvectors, assignments based on ',&
                          'rank-1 approximation :'

!     Print mode numbers
      mst=firstmode(il,1,im,ML)
      mfi=lastmode(il,1,im,ML)
      nagn=mfi-mst+1
      write(frmt,*) '(A,X,',nagn,'(I2,X),5X,A,14X,A,11X,A,5X,A)'
      write(*,frmt) 'Mode:',(ML%resort(j),j=mst,mfi),'Energy',&
                    'E-E0','Assignment','delta'

      DO i=1,nev
!        Get the full assignment
         call GetFullAssignment(il,im,Ham,ML,Ham%eig(il,im)%assgn(i,:),&
                                qns)

!        Label vibration if only 1 DOF is excited
         nexc=0
         iexc=0
         jexc=0
         DO j=1,nagn
            IF (qns(j).gt.1) THEN
               nexc=nexc+1
               iexc=qns(j)
               jexc=j
               IF (nexc.eq.2) EXIT
            ENDIF
         ENDDO
         IF (nexc.eq.0) THEN
            labl='ZPVE'
         ELSEIF (nexc.eq.1 .and. iexc.eq.2) THEN
            labl='FUND nu_'
         ELSEIF (nexc.eq.1 .and. iexc.gt.2) THEN
            labl=' OT '
         ELSE
            labl='    '
         ENDIF

         IF (nexc.ne.1 .or. iexc.ne.2) THEN
            write(frmt,*) '(I4,A,X,',nagn,&
                          '(I2,X),2(f19.12,X),A,5X,ES11.3)'
            write(*,frmt) i,')',(qns(j)-1,j=1,nagn),eigv(i),&
                          eigv(i)-eigv(1),labl,delta(i)
         ELSE
            sp=3-int(log10(REAL(ML%resort(mst+jexc-1))))+1
            write(frmt,*) '(I4,A,X,',nagn,&
                  '(I2,X),2(f19.12,X),A,I0,',sp,'X,ES11.3)'
            write(*,frmt) i,')',(qns(j)-1,j=1,nagn),eigv(i),&
                    eigv(i)-eigv(1),labl,ML%resort(mst+jexc-1),delta(i)
         ENDIF
         DEALLOCATE(qns)
      ENDDO

      call CPU_TIME(t2)
      anal_time=anal_time+t2-t1

      end subroutine PrintAssignments

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ANALYZER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
