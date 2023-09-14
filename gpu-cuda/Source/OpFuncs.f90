!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE OPFUNCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module performs the Hamiltonian operations on the primitive basis

      USE ERRORTRAP
      USE UTILS

      TYPE OperMat
         integer :: dof
         character(7) :: label
         real*8, allocatable :: mat(:)
      END TYPE OperMat

      INTERFACE SumOperMats
                MODULE PROCEDURE SumOperMats_BA
                MODULE PROCEDURE SumOperMats_ABA
      END INTERFACE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetPrimitiveOperMat(dof,N,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills type OperMat with a primitive operator matrix (p_i^2 or q_i^P)
! in harmonic oscillator basis
! dof = degree-of-freedom that primitive operator acts upon
! N = number of basis fxns for that DOF

      implicit none
      TYPE (OperMat)      :: GetPrimitiveOperMat
      integer, intent(in) :: dof,N,P
      real*8, allocatable :: Md(:,:),Mt(:,:)
      character(7) :: lab

!     Operator defs and calls
      IF (P.ge.0) THEN      ! Harmonic oscillator q_i^P
         write(lab,'(A,I0,A)') '(q^',P,')_'
         call xnop(N,P,Md)
      ELSEIF (P.eq.-2) THEN ! Harmonic oscillator KEO
         lab='(p^2)_'
         call top(N,Md)
      ELSE
         write(*,*) 'Operator type not recognized'
         call AbortWithError('Error in GetPrimitiveOperMat()')
      ENDIF

      call MatDiag2Utriang(Md,Mt,MOD(abs(P),2),.FALSE.)
      call SymPackMat2Vec(GetPrimitiveOperMat%mat,Mt)
      DEALLOCATE(Md,Mt)

      GetPrimitiveOperMat%dof=dof
      GetPrimitiveOperMat%label=lab

      end function GetPrimitiveOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetEigenOperMat(dof,eigvals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills type OperMat with diagonal operator matrix
! (i.e. list of eigenvalues) provided as input
! dof = degree-of-freedom/mode that operator acts upon

      implicit none
      TYPE (OperMat)      :: GetEigenOperMat
      integer, intent(in) :: dof
      real*8, intent(in)  :: eigvals(:)
      integer :: i,s,n

      n=SIZE(eigvals)

      GetEigenOperMat%dof=dof
      GetEigenOperMat%label='eigen'

!     Fill the operator matrix with eigenvalues
      ALLOCATE(GetEigenOperMat%mat(n*(n+1)/2))
      GetEigenOperMat%mat=0.d0
      s=0
      DO i=1,n
         s=s+i
         GetEigenOperMat%mat(s)=eigvals(i)
      ENDDO

      end function GetEigenOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetIdentityOperMat(dof,n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills type OperMat with n x n identity operator matrix
! dof = degree-of-freedom/mode that operator acts upon

      implicit none
      TYPE (OperMat)      :: GetIdentityOperMat
      integer, intent(in) :: dof,n
      integer :: i,s

      GetIdentityOperMat%dof=dof
      GetIdentityOperMat%label='I'

!     Fill the operator matrix with 1s and 0s
      ALLOCATE(GetIdentityOperMat%mat(n*(n+1)/2))
      GetIdentityOperMat%mat=0.d0
      s=0
      DO i=1,n
         s=s+i
         GetIdentityOperMat%mat(s)=1.d0
      ENDDO

      end function GetIdentityOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SumOperMats_ABA(OM1,fac1,OM2,fac2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums operator matrices OM1 and OM2, replacing OM1
! Matrices OM1 and OM2 are multiplied by fac1 and fac2, respectively

      implicit none
      TYPE (OperMat), intent(inout) :: OM1
      TYPE (OperMat), intent(in)    :: OM2
      real*8, intent(in) :: fac1,fac2
      integer :: n1,n2

      n1=SIZE(OM1%mat)
      n2=SIZE(OM2%mat)

!     Error checking
      IF (n1.ne.n2) THEN
         write(*,*) 'Operator array size mismatch'
         CALL AbortWithError('Error in SumOperMatsABA()')
      ENDIF
      IF (OM1%dof.ne.OM2%dof) THEN
         write(*,*) 'Summed operators correspond to different DOFs'
         CALL AbortWithError('Error in SumOperMatsABA()')
      ENDIF

!     Sum the matrices
      OM1%mat=fac1*OM1%mat+fac2*OM2%mat

!     New parameter for sum
      OM1%label='summed'

      end subroutine SumOperMats_ABA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SumOperMats_BA(OM1,OM2,fac2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies fac2*OM2 into OM1

      implicit none
      TYPE (OperMat), intent(out) :: OM1
      TYPE (OperMat), intent(in)  :: OM2
      real*8, intent(in) :: fac2
      integer :: n

      n=SIZE(OM2%mat)

      ALLOCATE(OM1%mat(n))
      OM1%mat=fac2*OM2%mat
      OM1%dof=OM2%dof
      OM1%label=OM2%label

      end subroutine SumOperMats_BA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushOperMat(OM)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes OperMat type

      implicit none
      TYPE (OperMat) :: OM

      IF (ALLOCATED(OM%mat)) DEALLOCATE(OM%mat)
      OM%dof=-1

      end subroutine FlushOperMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      subroutine xnop(N,P,Hmat)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates the Hamiltonian matrix for <phi|x^P|phi> in harmonic
! oscillator basis. Returns an array whose columns correspond to 
! <i|x^P|i>, <i,x^P,i+-2>, <i,x^P,i+-4>, etc, for P even and
! <i,x^P,i+-1>, <i,x^P,i+-3>, etc, for P odd
! N is the basis set size (HO basis: 0-->N-1)

      implicit none
      real*8, allocatable, intent(out) :: Hmat(:,:)
      integer, intent(in) :: N,P
      integer :: i,imod,j,l,j2,tp,mp,wp,cols,ll,facp,facm
      real*8  :: fac1,fac2,sm

      tp=P/2
      mp=MOD(P,2)
      IF (N-mp.lt.1) &
         CALL AbortWithError('Error in xnop(): matrix of size zero')

      fac1=0.5d0**(0.5d0*P) * FACRL(P)
      cols=min(tp+1,(N-mp+1)/2)
      ALLOCATE(Hmat(N-mp,cols))

      DO i=1,cols
         imod=mp+2*(i-1)
         IF (i.gt.1) fac1=fac1*sqrt(imod*(imod-1.d0))
         fac2=fac1
         DO j=1,N-mp
            j2=j-1+imod
            IF (j2.gt.N-1) THEN
               Hmat(j,i)=0.d0
            ELSE
               wp=(P-j2-j+1)/2
               l=MAX(0,-wp)
               IF (j.gt.1) fac2=fac2*SQRT((j-1.d0)*j2)
               IF (l.gt.1) fac2=fac2/l
               facp=2**(wp+l)*FACRL(wp+l)
               facm=FACRL(j-1-l)*FACRL(j2-l)
               sm=0.d0
               DO
                  IF (l.ge.j) EXIT
                  sm=sm+fac2/facp/facm
                  l=l+1
                  facp=facp*2*l*(wp+l)
                  IF (j.gt.l) facm=facm/(j-l)
                  IF (j2.ge.l) facm=facm/(j2+1-l)
               ENDDO
               Hmat(j,i)=sm !*fac2
            ENDIF
         ENDDO
      ENDDO

      end subroutine xnop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      subroutine top(N,Hmat)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates the operator matrix for <phi|T|phi> in harmonic
! oscillator basis, using a call to xnop()
! N is the basis set size (HO basis: 0-->N-1)

      implicit none
      real*8, dimension(:,:), allocatable, intent(out) :: Hmat
      integer, intent(in) :: N

!     Get <phi|x^2|phi> matrix
      call xnop(N,2,Hmat)

!     Change sign of diagonal +- 2 elements
      IF (SIZE(Hmat,2).gt.1) Hmat(:,2)=-Hmat(:,2)

      end subroutine top

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE OPFUNCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
