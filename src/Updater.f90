!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE UPDATER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module updates the Hamiltonian

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE INPUTCP
      USE SEPDREPN
      USE HAMILSETUP
      USE MODECOMB
      USE MODVECVEC
      USE BLOCKUTILS

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_Updater_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_Updater_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_Updater_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_Updater_Module()
      call Get_MPI_Timings('Updater module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_Updater_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateH(im,Q,Ham,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Updates H by building the operator matrices for the upper layers 
! (il > 1) and transforming the multi-mode operators into the eigenbasis

      implicit none
      TYPE (CPpar)       :: cpp
      TYPE (Hamiltonian) :: Ham
      TYPE (CP), INTENT(IN) :: Q(:)
      integer, intent(in) :: im
      integer :: i,j,k,l,nev,nop,nsubm,sm,nsubdof,subdof,noptyp
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_Updater_Module()

      nev=SIZE(Q)

!     Set the node's done status 
      Ham%nt(im)%done=.TRUE.

!     No operator update necessary for the following cases:
!     - memory check run
!     - node has only one sub-node (without truncation)
      IF (cpp%ncycle.eq.0 .or. &
          (Ham%nt(im)%nHterm().eq.0 .and. &  ! Presolved
           Ham%nt(im)%nsubm() .eq.1 .and. &
           Ham%nt(Ham%nt(im)%subm(1))%nbas().eq.nev)) RETURN

      call CPU_TIME(ti1)

      nsubm=max(1,Ham%nt(im)%nsubm())  ! Nr of sub-nodes in this node
      nop=SIZE(Ham%ops,2)
      noptyp=SIZE(Ham%ops,3)

!     Transform primitive operators with mode 'im' into the eigenbasis
      DO i=1,nsubm
         sm=Ham%nt(im)%subm(i)  ! Sub-node index
         nsubdof=Ham%nt(sm)%ndof() ! Nr of DOFs in this sub-node
         DO j=1,nsubdof
            subdof=Ham%nt(sm)%dofs(j) ! DOF index
            DO l=1,noptyp
               DO k=1,nop          ! Nr of potential ops for this DOF
                  IF (Ham%optable(subdof,k,l)) THEN
!                    Transform primitive operators into basis in Q
                     IF (cpp%algo.lt.0) then
                        call UpdateOperMat(i,Q,Ham%ops(subdof,k,l)%mat)
                     ELSE
                        call UpdateOperMat_CP8(i,Q,Ham%ops(subdof,k,l)%mat,cpp%algo)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine UpdateH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateOperMat(imode,Q,mat)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T X Q, where X is primitive operator matrix applying to
! imode of Q

      implicit none
      real(kind=8), allocatable, intent(inout) :: mat(:)
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
         call PRODXV(imode,Q(i),XQ,mat)
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

      DEALLOCATE(mat)
      call SymPackMat2Vec(mat,QXQ)
      DEALLOCATE(QXQ)

      end subroutine UpdateOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateOperMat_CP8(imode,Q,mat,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T X Q, where X is primitive operator matrix applying to
! imode of Q

      implicit none
      real(kind=8), allocatable, intent(inout) :: mat(:)
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (CP8), allocatable :: Qmode(:)
      TYPE (CP8) :: X8
      real(kind=8), allocatable :: QXQ(:,:),tmat(:,:),tvec(:)
      integer, allocatable :: rows(:),cols(:)
      integer, intent(in)  :: imode,algo
      integer :: i,nbloc

      nbloc=SIZE(Q)

      ALLOCATE(QXQ(nbloc,nbloc),Qmode(nbloc))
      QXQ=0.d0

!     Put the operator matrix into CP8 structure
      call Vec2SymPackMat(mat,tmat)
      call Mat2Vec(tvec,tmat,.FALSE.)
      call X8%identity(Q(1)%rows,Q(1)%rows)
      call PutModeTerm_CP8(X8,1,imode,tvec)
      DEALLOCATE(mat,tmat,tvec)

!     Put Q into CP8 structure
      do i=1,nbloc
         call Qmode(i)%fromCP(Q(i))
      enddo

      if (algo.eq.1) then
         call X8%copyintodevice()
         do i=1,nbloc
            call Qmode(i)%copyintodevice()
         enddo
      endif

!     Use the block MPI code to construct QXQ
      call GetQXQ_CP8(Qmode,X8,QXQ,algo,imode)

!     Replace the operator matrix
      call SymPackMat2Vec(mat,QXQ)
      call X8%flush()
      do i=1,nbloc
         call Qmode(i)%flush()
      enddo
      DEALLOCATE(QXQ,Qmode)

      end subroutine UpdateOperMat_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODXV(imode,F,G,mat)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product X*F = G, where X is an individual term
! in the Hamiltonian and F and G are CP-format vectors.
! imode = the mode in F to which X applies

      implicit none
      TYPE (CP), INTENT(IN)    :: F
      TYPE (CP), INTENT(OUT)   :: G
      real(kind=8), intent(in) :: mat(:)
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
         call dspmv('U',gdim,1.d0,mat,F%base(gi:gf,i),1,0.d0,&
                    G%base(gi:gf,i),1)
      ENDDO

      end subroutine PRODXV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE UPDATER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
