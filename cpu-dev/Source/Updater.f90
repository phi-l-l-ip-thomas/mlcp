!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE UPDATER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module updates the Hamiltonian

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE SEPDREPN
      USE HAMILSETUP
      USE MODECOMB
      USE MODVECVEC

      implicit none
      real*8, private  :: update_time=0.d0
      logical, private :: UPDATE_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeUpdateModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      update_time = 0.d0
      UPDATE_SETUP = .TRUE.

      end subroutine InitializeUpdateModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeUpdateModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. UPDATE_SETUP) call InitializeUpdateModule()

      UPDATE_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total operator update time        (s)',&
                             update_time

      end subroutine DisposeUpdateModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateH(il,im,eigv,Q,H,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Updates H by building the operator matrices for the upper layers 
! (il > 1) and transforming the multi-mode operators into the eigenbasis

      implicit none
      TYPE (CPpar)       :: cpp
      TYPE (MLtree)      :: ML
      TYPE (Hamiltonian) :: H
      TYPE (CP), INTENT(IN) :: Q(:)
      real*8, allocatable, intent(in) :: eigv(:)
      integer, intent(in) :: il,im
      integer :: i,j,imode,nev,mpl,mstart,trm
      real*8  :: t1,t2

      IF (.NOT. UPDATE_SETUP) call InitializeUpdateModule()

      call CPU_TIME(t1)

      nev=SIZE(Q)
      mstart=ML%modstart(il,im)
      trm=GetModeHNr(il,im,H)  ! mode term

!     Store the eigenvalues for mode 'im' that were computed in
!     the previous call to the solver
      ALLOCATE(H%eig(il,im)%evals(SIZE(eigv)))
      H%eig(il,im)%evals(:)=eigv(:)

!     No operator update necessary for the last layer
      IF (il.eq.ML%nlayr) RETURN

!     Update operators if mode is not pre-solved or if mode is 
!     presolved and the basis is truncated in the current layer

      IF ((H%ndof(trm,il).gt.1 .or. H%nop(trm,il).gt.1 .or. &
          SIZE(H%eig(il-1,mstart)%evals).gt.nev) .and. &
          cpp%ncycle.gt.0) THEN

!        Transform primitive operators with mode 'im' into the eigenbasis
         DO i=1,SIZE(H%pops)
            IF (uppermodenr(il,1,H%pops(i)%dof,ML).eq.im) THEN
!              Determine imode
               imode=uppermodenr(il-1,1,H%pops(i)%dof,ML)-mstart+1
!              Update operator
               call UpdateOperMat(imode,Q,H%pops(i))
            ENDIF
         ENDDO
      ENDIF

      call CPU_TIME(t2)
      update_time=update_time+t2-t1

      end subroutine UpdateH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateOperMat(imode,Q,X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T X Q, where X is primitive operator matrix applying to
! imode of Q

      implicit none
      TYPE (OperMat), INTENT(INOUT) :: X
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (CP) :: XQ
      real*8, allocatable :: QXQ(:,:)
      integer, intent(in) :: imode
      integer :: i,j,nbloc

      nbloc=SIZE(Q)

      ALLOCATE(QXQ(nbloc,nbloc))
      QXQ=0.d0

      DO i=1,nbloc
!        XQ=X*Q(i)
         call PRODXV(imode,Q(i),XQ,X)
!$omp parallel
!$omp do private(j)
         DO j=i,nbloc
!           QXQ(i,j)=<Q(j),X(Q(i))>
            QXQ(i,j)=PRODVV(Q(j),XQ)
         ENDDO
!$omp end do
!$omp end parallel 
         call FlushCP(XQ)
      ENDDO

!     Replace the operator matrix
      DEALLOCATE(X%mat)
      call SymPackMat2Vec(X%mat,QXQ)
      DEALLOCATE(QXQ)

      end subroutine UpdateOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODXV(imode,F,G,X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product X*F = G, where X is an individual term
! in the Hamiltonian and F and G are CP-format vectors.
! imode = the mode in F to which X applies

      implicit none
      TYPE (OperMat), INTENT(IN) :: X
      TYPE (CP), INTENT(IN)   :: F
      TYPE (CP), INTENT(OUT)  :: G
      integer, intent(in) :: imode
      integer :: i,rF,gdim,gst,gi,gf

!     Set parameters
      rF=SIZE(F%coef)  ! rank of F
      gdim=F%nbas(imode)
      gst=0
      IF (imode.gt.1) THEN
         DO i=2,imode
            gst=gst+F%nbas(i-1)
         ENDDO
      ENDIF
      gi=gst+1
      gf=gst+gdim

      G=CopyCP(F)

!     Operation X*V
      DO i=1,rF
         call dspmv('U',gdim,1.d0,X%mat,F%base(gi:gf,i),1,0.d0,&
                    G%base(gi:gf,i),1)
      ENDDO

      end subroutine PRODXV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE UPDATER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
