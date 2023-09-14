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
!!!
      real*8, allocatable :: BB(:,:)
!!!

      nev=SIZE(Q)
      nsubm=SIZE(Q(1)%nbas)

      IF (nsubm.lt.2) RETURN

      call CPU_TIME(t1)

!     Eigenvector printing
!      write(*,*)
!      DO i=1,nev
!         write(*,*) 'Eigenvector: ',i,'; Condition Nr: ',COND(Q(i))
!         call Q(i)%printvec()
!!! TEST
!         ALLOCATE(BB(Q(i)%R(),Q(i)%R()))
!         DO j=1,nsubm
!            write(*,*) 'Eigenvector: ',i,'; BBk for mode: ',j
!            call CONSTPk(Q(i),j,BB)
!            call PrintMatrix(BB)
!            write(*,*)
!            write(*,*) 'Factor matrix SVD:'
!            call PrintSVD(Q(i)%base(Q(i)%ibas(j):Q(i)%fbas(j),:))
!         ENDDO
!         DEALLOCATE(BB)
!!!
!      ENDDO

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

      subroutine GetFullAssignment(il,im,Ham,ML,qns,qnfull)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in)    :: qns(:)
      integer, allocatable, intent(out) :: qnfull(:)
      integer, intent(in)    :: il,im
      integer, allocatable   :: modind(:),qntmp(:)
      integer :: i,j,eigind,imn,nsubm,mst,nagn

!     Get the full assignment in terms of primitive DOFs
      ALLOCATE(modind(il),qntmp(ML%nmode(1)))

      qntmp=0
      nagn=0
      modind=1
      DO
         imn=im

         IF (modind(1).gt.1) EXIT

!        Trace each assignment to the bottom layer
         DO j=il,2,-1
            mst=ML%modstart(j,imn)

!           Take number of sub-modes from input array on first pass,
!           then extract the number from stored assignment array
            IF (j.eq.il) THEN
               nsubm=SIZE(qns)
            ELSE
               nsubm=SIZE(Ham%eig(j,imn)%assgn,2)
            ENDIF

            IF (modind(il-j+2).gt.nsubm) THEN
               modind(il-j+1)=modind(il-j+1)+1
               modind(il-j+2:)=1
               EXIT
            ENDIF

!           Use input array to get 'eigind' on first pass here, too
            IF (j.eq.il) THEN
               eigind=qns(modind(2))
            ELSE
               eigind=Ham%eig(j,imn)%assgn(eigind,modind(il-j+2))
            ENDIF

            imn=mst+modind(il-j+2)-1

            IF (j.eq.2) THEN
               nagn=nagn+1
               qntmp(nagn)=eigind
               modind(il-j+2)=modind(il-j+2)+1
            ENDIF
         ENDDO
      ENDDO

!     Copy qntmp to qnfull
      ALLOCATE(qnfull(nagn))
      qnfull(:)=qntmp(:nagn)
      DEALLOCATE(modind,qntmp)

      end subroutine GetFullAssignment

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ANALYZER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
