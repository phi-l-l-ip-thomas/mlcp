!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE PRECONDITION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module builds the CP-format mode Hamiltonian

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE LINALG
      USE SEPDREPN
      USE INPUTCP
      USE MODECOMB
      USE MODVECVEC
      USE REDUCTION
      USE ALSDRVR

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1Inverse(H,Hi,ish,Esh,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes approximate H^-1 by reducing H to rank-1 and inverting each
! coordinate Hamiltonian. On input, H should be given as the full matrix 
! rep'n, not the upper triangle, (sym=.FALSE.)

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(OUT) :: Hi
      TYPE (CP) :: I
      integer, intent(in) :: ish
      real(kind=8), intent(in) :: Esh
      real(kind=8), allocatable :: tmat(:,:),tmat1(:,:)
      real(kind=8), allocatable :: tvec(:)
      character(len=64), parameter :: solver='LU'
      integer, allocatable :: nbas(:)
      logical, allocatable :: sym(:)
      integer :: j,ndof,gi,gf

      ndof=SIZE(H%nbas)
      ALLOCATE(nbas(ndof),sym(ndof))

      DO j=1,ndof
         nbas(j)=NINT(sqrt(REAL(H%nbas(j))))
      ENDDO
      sym(:)=.FALSE.

!     Hshifted = H - E*I for ish>0
!           or = E*I - H for ish<0
!     store in I (not ideal since creates full rank copy of H)
      IF (ish.ne.0) THEN 
         I=IdentityCPMatrix(nbas,nbas,sym)
         call SUMVECVEC(I,-Esh,H,1.d0)
         IF (ish.lt.0) call VecSignChange(I,1,I%R())
      ELSE
         I=CopyCP(H)
      ENDIF


!     Initial guess: reduce Hi <- H, with Hi rank-1
      call SetReductionParameters(1,cpp%hnals,1.d-12,.FALSE.,'SVD','SR1',&
                                  cpp%alspenalty,cpp%als_linsys_alg)
      call reduc(Hi,I)
      call FlushCP(I)

!     Invert each little-h in Hi
      gi=1
      DO j=1,ndof
         gf=gi+Hi%nbas(j)-1
         call Vec2Mat(Hi%base(gi:gf,1),tmat,nbas(j),nbas(j))
         call MatrixPseudoinverse(tmat,tmat1)
!        Symmetrize to correct small numerical errors from ALS
         call SymmetrizeMat(tmat1)
         call Mat2Vec(tvec,tmat1,.FALSE.)
         Hi%base(gi:gf,1)=tvec(1:Hi%nbas(j))
         DEALLOCATE(tmat,tmat1,tvec)
         gi=gf+1
      ENDDO
      Hi%coef(1)=1.d0/Hi%coef(1)

      DEALLOCATE(nbas,sym)

      end subroutine GetRank1Inverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE PRECONDITION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
