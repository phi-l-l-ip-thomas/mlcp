!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE OPFUNCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module performs the Hamiltonian operations on the primitive basis

      USE ERRORTRAP
      USE UTILS
      USE CHEBLIB

      TYPE OperMat
         integer :: dof
         character(64) :: label
         real*8, allocatable :: mat(:)
      END TYPE OperMat

      INTERFACE SumOperMats
                MODULE PROCEDURE SumOperMats_BA
                MODULE PROCEDURE SumOperMats_ABA
      END INTERFACE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetPrimitiveOperMat(dof,N,P,typ,alpha)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills type OperMat with a primitive operator matrix (p_i^2 or q_i^P)
! in harmonic oscillator basis
! dof = degree-of-freedom that primitive operator acts upon
! N = number of basis fxns for that DOF

      implicit none
      TYPE (OperMat)      :: GetPrimitiveOperMat
      integer, intent(in) :: dof,N,P,typ
      real*8, intent(in)  :: alpha
      real*8, allocatable :: Md(:,:),Mt(:,:)
      integer, parameter  :: ip=0
      real*8, parameter   :: rp=0.d0
      character(64) :: lab

!     Operator defs and calls
      IF (P.ge.0) THEN      ! Harmonic oscillator q_i^P
         IF (typ.eq.0) THEN
            write(lab,'(A,I0,A)') '(q^',P,')_'
!           Power of q
            call xnop(N,P,Md)
            call MatDiag2Utriang(Md,Mt,MOD(abs(P),2),.FALSE.)
            DEALLOCATE(Md)

         ELSE
!           Arbitrary PES type
!           For now only one integer parameter and one real parameter
!           are used; the rest are called with dummy arguments
            write(lab,'(3A,I0,A)') '(',&
                  TRIM(ADJUSTL(GetFunctionLabel(typ))),'^',P,')_'
            call GenFunctionHOmat(N,Mt,typ,P,ip,alpha,rp,rp,rp,rp)
         ENDIF
      ELSEIF (P.eq.-2) THEN ! Harmonic oscillator KEO
         lab='(p^2)_'
         call top(N,Md)
         call MatDiag2Utriang(Md,Mt,MOD(abs(P),2),.FALSE.)
         DEALLOCATE(Md)

      ELSE
         write(*,*) 'Operator type not recognized'
         call AbortWithError('Error in GetPrimitiveOperMat()')
      ENDIF

      call SymPackMat2Vec(GetPrimitiveOperMat%mat,Mt)
      GetPrimitiveOperMat%dof=dof
      GetPrimitiveOperMat%label=lab

      DEALLOCATE(Mt)

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

      subroutine xnop(N,P,Hmat)

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

      subroutine top(N,Hmat)

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

      subroutine GenFunctionHOmat(N,Hmat,typ,ip1,ip2,&
                                  rp1,rp2,rp3,rp4,rp5)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates the operator matrix for <phi'|function|phi> for a general
! function in harmonic oscillator basis, using numerical integration
! N is the basis set size (HO basis: 0-->N-1)

      implicit none
      TYPE (ChebObj) :: chb
      real*8, allocatable, intent(out) :: Hmat(:,:)
      integer, intent(in)  :: N,typ,ip1,ip2
      real*8, intent(in)   :: rp1,rp2,rp3,rp4,rp5
      real*8, allocatable  :: v(:),phi(:),dummy(:)
      real*8, allocatable  :: phim(:,:)
      integer, parameter   :: ncheb=257
      real*8, parameter    :: tol=1.d-15
      integer :: i,j,k
      real*8  :: xlim

!     Get the limits of integration. These must exceed the classical
!     turning point of the highest-order HO function in the basis.
      xlim=sqrt(2.d0*(N-1)+1)
      DO
!        Step away from the origin until the amplitudes of the lowest,
!        highest HO functions in the basis are negligible
         IF (ABS(HObasisfxn(0,xlim)).lt.tol .and. &
             ABS(HObasisfxn(N-1,xlim)).lt.tol) EXIT
         xlim=xlim+1.d0
      ENDDO

!     Generate a new Chebyshev object for doing numerical integration
!     via Clenshaw-Curtis quadrature
      call NewChebObject(chb,ncheb-1,-xlim,xlim,.TRUE.,'fin')

!     Evaluate the PES and the HO basis at the quadrature points
      ALLOCATE(v(ncheb),phim(N,ncheb))
!$omp parallel
!$omp do private(k) schedule(static)
      DO k=1,ncheb
         v(k)=my1Dfunction(typ,chb%mpt(k),ip1,ip2,&
                           rp1,rp2,rp3,rp4,rp5)
         call HObasisseries(chb%mpt(k),phim(:,k))
      ENDDO
!$omp end do
!$omp end parallel

!     Get the upper triangle of the <phi|function|phi> matrix
      ALLOCATE(phi(ncheb),Hmat(N,N))
      Hmat=0.d0
      DO i=1,N
!        Construct v|phi_i>
         DO k=1,ncheb
             phi(k)=v(k)*phim(i,k)
         ENDDO
         DO j=i,N
!           Construct <phi_j|v|phi_i>
!$omp parallel
!$omp do private(k) schedule(static)
            DO k=1,ncheb
               chb%fpt(k)=phim(j,k)*phi(k)
            ENDDO
!$omp end do
!$omp end parallel
!           Compute the matrix element
            call ChebCalculus(chb,0,dummy)
            Hmat(i,j)=chb%defint
         ENDDO
      ENDDO

      call FlushChebObject(chb)
      DEALLOCATE(v,phim,phi)

      end subroutine GenFunctionHOmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetFunctionLabel(id)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer, intent(in) :: id
      character(len=32)   :: GetFunctionLabel

      IF (id.eq.0) THEN
         GetFunctionLabel='q'
      ELSEIF (id.eq.1) THEN
         GetFunctionLabel='tanh(a*q)'
      ELSEIF (id.eq.2) THEN
         GetFunctionLabel='1-exp(-a*q)'
      ELSEIF (id.eq.3) THEN
         GetFunctionLabel='sqrt(1-exp(-a*q^2))'
      ELSE
         write(*,*) 'function id = ',id
         call AbortWithError("GetFunctionLabel(): unknown function id")
      ENDIF

      end function GetFunctionLabel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE OPFUNCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
