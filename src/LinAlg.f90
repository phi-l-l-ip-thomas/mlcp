!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE LINALG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains wrappers for LAPACK subroutines

      USE ERRORTRAP
      USE UTILS

      implicit none
      real*8, private  :: eigen_time=0.d0
      logical, private :: EIGEN_SETUP = .FALSE.

      INTERFACE VecVecProd
         MODULE PROCEDURE VecVecProd_ABA,VecVecProd_ABC
      END INTERFACE VecVecProd

      INTERFACE MatVecProd
         MODULE PROCEDURE MatVecProd_ABB,MatVecProd_ABC
      END INTERFACE MatVecProd

      INTERFACE MatrixMult
         MODULE PROCEDURE MatrixMult_ABA,MatrixMult_ABC
      END INTERFACE MatrixMult

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeEigen()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      eigen_time = 0.d0
      EIGEN_SETUP = .TRUE.

      end subroutine InitializeEigen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeEigen()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. EIGEN_SETUP) call InitializeEigen()

      EIGEN_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total eigenvalue calculation time (s)',&
                            eigen_time

      end subroutine DisposeEigen


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixMult_ABA(M1,t1,M2,t2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M1 <- M1 x M2. Set tx (x={1,2}) to .TRUE. to do transpose

      implicit none
      real*8, allocatable, intent(inout) :: M1(:,:)
      real*8, intent(in)  :: M2(:,:)
      real*8, allocatable :: Mt(:,:)
      logical, intent(in) :: t1,t2
      integer :: r1,r2,c1,c2,ld1,ld2
      character(1) :: tN1,tN2

      ld1=SIZE(M1,1)
      ld2=SIZE(M2,1)

      r1=SIZE(M1,1)
      c1=SIZE(M1,2)
      tN1='N'
      IF (t1) THEN
         r1=SIZE(M1,2)
         c1=SIZE(M1,1)
         tN1='T'
      ENDIF
      r2=SIZE(M2,1)
      c2=SIZE(M2,2)
      tN2='N'
      IF (t2) THEN
         r2=SIZE(M2,2)
         c2=SIZE(M2,1)
         tN2='T'
      ENDIF

      IF (c1.ne.r2) THEN
         write(*,*) 'cols(M1) = ',c1,'; rows(M2) = ',r2
         write(*,*) 'Error: mismatch in dimensions between M1 and M2'
         call AbortWithError('Error in MatrixMult_ABA()')
      ENDIF

      ALLOCATE(Mt(r1,c2))
      call DGEMM(tN1,tN2,r1,c2,c1,1.d0,M1,ld1,M2,ld2,0.d0,Mt,r1)

      DEALLOCATE(M1)
      ALLOCATE(M1(r1,c2))
      M1=Mt
      DEALLOCATE(Mt)

      end subroutine MatrixMult_ABA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixMult_ABC(M1,t1,M2,t2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M3 <- M1 x M2. Set tx (x={1,2}) to .TRUE. to do transpose

      implicit none
      real*8, intent(in)  :: M1(:,:)
      real*8, intent(in)  :: M2(:,:)
      real*8, allocatable, intent(out) :: M3(:,:)
      logical, intent(in) :: t1,t2
      integer :: r1,r2,c1,c2,ld1,ld2
      character(1) :: tN1,tN2

      ld1=SIZE(M1,1)
      ld2=SIZE(M2,1)

      r1=SIZE(M1,1)
      c1=SIZE(M1,2)
      tN1='N'
      IF (t1) THEN
         r1=SIZE(M1,2)
         c1=SIZE(M1,1)
         tN1='T'
      ENDIF
      r2=SIZE(M2,1)
      c2=SIZE(M2,2)
      tN2='N'
      IF (t2) THEN
         r2=SIZE(M2,2)
         c2=SIZE(M2,1)
         tN2='T'
      ENDIF

      IF (c1.ne.r2) THEN
         write(*,*) 'cols(M1) = ',c1,'; rows(M2) = ',r2
         write(*,*) 'Error: mismatch in dimensions between M1 and M2'
         call AbortWithError('Error in MatrixMult_ABC()')
      ENDIF

      ALLOCATE(M3(r1,c2))
      call DGEMM(tN1,tN2,r1,c2,c1,1.d0,M1,ld1,M2,ld2,0.d0,M3,r1)

      end subroutine MatrixMult_ABC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatVecProd_ABB(M,t,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M*v, replacing v. Set 't' to .TRUE. to use transpose of M

      implicit none
      real*8, intent(in)  :: M(:,:)
      real*8, allocatable, intent(inout) :: v(:)
      real*8, allocatable :: vtmp(:)
      logical, intent(in) :: t
      integer :: mr,mc,i,j

      IF (t) THEN
         mr=SIZE(M,2)
         mc=SIZE(M,1)
      ELSE
         mr=SIZE(M,1)
         mc=SIZE(M,2)
      ENDIF

      IF (mc.ne.SIZE(v)) THEN
         write(*,*) 'cols(M1) = ',mc,'; rows(v) = ',SIZE(v)
         write(*,*) 'Error: mismatch in dimensions between M and v'
         call AbortWithError('Error in MatVecProd_ABB()')
      ENDIF

      ALLOCATE(vtmp(mc))
      vtmp=v
      DEALLOCATE(v)
      ALLOCATE(v(mr))
      v=0.d0

      DO i=1,mr
         DO j=1,mc
            IF (t) THEN
               v(i)=v(i)+vtmp(j)*M(j,i)
            ELSE
               v(i)=v(i)+vtmp(j)*M(i,j)
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(vtmp)

      end subroutine MatVecProd_ABB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatVecProd_ABC(M,t,v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M*v = w. Set 't' to .TRUE. to use transpose of M

      implicit none
      real*8, intent(in)  :: M(:,:)
      real*8, intent(in)  :: v(:)
      real*8, intent(out) :: w(:)
      logical, intent(in) :: t
      integer :: mr,mc,i,j

      IF (t) THEN
         mr=SIZE(M,2)
         mc=SIZE(M,1)
      ELSE
         mr=SIZE(M,1)
         mc=SIZE(M,2)
      ENDIF

      IF (mc.ne.SIZE(v)) THEN
         write(*,*) 'cols(M1) = ',mc,'; rows(v) = ',SIZE(v)
         write(*,*) 'Error: mismatch in dimensions between M and v'
         call AbortWithError('Error in MatVecProd_ABC()')
      ENDIF
      IF (mr.ne.SIZE(w)) THEN
         write(*,*) 'rows(M1) = ',mr,'; rows(w) = ',SIZE(w)
         write(*,*) 'Error: mismatch in dimensions between M and w'
         call AbortWithError('Error in MatVecProd_ABC()')
      ENDIF

      w=0.d0

      DO i=1,mr
         DO j=1,mc
            IF (t) THEN
               w(i)=w(i)+v(j)*M(j,i)
            ELSE
               w(i)=w(i)+v(j)*M(i,j)
            ENDIF
         ENDDO
      ENDDO

      end subroutine MatVecProd_ABC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WrappedMatVecProd(M,v,w,t)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M*v = w. Set 't' to .TRUE. to use transpose of M
! Here M is stored as a vector with dimensions (SIZE(v)*SIZE(w))

      implicit none
      real*8, intent(in)    :: M(:)
      real*8, intent(in)    :: v(:)
      real*8, intent(inout) :: w(:)
      logical, intent(in)   :: t
      integer :: mr,mc,i,j,k,mci

      mr=SIZE(w)
      mc=SIZE(v)

      IF (mr*mc.ne.SIZE(M)) THEN
         write(*,'(3(A,I0),A)') 'v(',mc,') x w(',mr,&
                                ') != M(',SIZE(M),')'
         call AbortWithError('WrappedMatVecProd(): size mismatch')
      ENDIF

      w(:)=0.d0

      IF (t) THEN
         DO i=1,mr
            mci=mc*(i-1)
            DO j=1,mc
               k=j+mci
               w(i)=w(i)+M(k)*v(j)
            ENDDO
         ENDDO
      ELSE
         DO i=1,mr
            DO j=1,mc
               k=i+mr*(j-1)
               w(i)=w(i)+M(k)*v(j)
            ENDDO
         ENDDO
      ENDIF

      end subroutine WrappedMatVecProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecVecProd_ABA(v1,v2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies v1=v1*v2, entrywise multiplication

      implicit none
      real*8, intent(inout) :: v1(:)
      real*8, intent(in)    :: v2(:)
      integer :: i,sv

      sv=SIZE(v1)
      IF (sv.ne.SIZE(v2)) THEN
         write(*,*) 'rows(v1) = ',sv,'; rows(v2) = ',SIZE(v2)
         write(*,*) 'Error: mismatch in dimensions between v1 and v2'
         call AbortWithError('Error in VecVecProd_ABA()')
      ENDIF

      DO i=1,sv
         v1(i)=v1(i)*v2(i)
      ENDDO

      end subroutine VecVecProd_ABA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecVecProd_ABC(v1,v2,v3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies v3=v1*v2, entrywise multiplication

      implicit none
      real*8, intent(inout) :: v3(:)
      real*8, intent(in)    :: v1(:),v2(:)
      integer :: i,sv

      sv=SIZE(v1)
      IF (sv.ne.SIZE(v2)) THEN
         write(*,*) 'rows(v1) = ',sv,'; rows(v2) = ',SIZE(v2)
         write(*,*) 'Error: mismatch in dimensions between v1 and v2'
         call AbortWithError('Error in VecVecProd_ABC()')
      ENDIF
      IF (sv.ne.SIZE(v3)) THEN
         write(*,*) 'rows(v1) = ',sv,'; rows(v3) = ',SIZE(v3)
         write(*,*) 'Error: mismatch in dimensions between v1 and v3'
         call AbortWithError('Error in VecVecProd_ABC()')
      ENDIF

      DO i=1,sv
         v3(i)=v1(i)*v2(i)
      ENDDO

      end subroutine VecVecProd_ABC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveTridiag(avec,bvec,QHQ,COMPZ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given Lanczos a and b vectors (and possibly QHQ), solve the eigenvalue
! problem with a call to the LAPACK subroutine DSTEQR
! COMPZ is the problem type:
! 'N' -> just eigvals
! 'I' ->  " + eigvecs of tridiagonal
! 'V' ->  " + eigvecs of original matrix (give QHQ for this)

      implicit none
      real*8, intent(inout) :: QHQ(:,:)
      real*8, intent(inout) :: avec(:),bvec(:)
      real*8, allocatable   :: WORK(:)
      integer     :: n,INFO,LWORK
      real*8      :: t1,t2
      character*1 :: COMPZ

      IF (.NOT. EIGEN_SETUP) call InitializeEigen()
      call CPU_TIME(t1)

      n=SIZE(avec)

      LWORK=max(1,2*n-2)
      ALLOCATE(WORK(LWORK))
      call DSTEQR(COMPZ,n,avec,bvec(2:n),QHQ,n,WORK,INFO)
      DEALLOCATE(WORK)
      if (INFO.ne.0) then
         write(*,*)'error in the diagonalization, info=',INFO
         call AbortWithError('error in DSTEQR')
      endif

      call CPU_TIME(t2)
      eigen_time=eigen_time+t2-t1

      end subroutine SolveTridiag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveGenEigval(avec,S,QHQ,COMPZ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given QHQ and S, solve the generalized eigenvalue problem using a call
! to the LAPACK subroutine DSYGV (eigvals returned in avec)
! COMPZ is the problem type:
! 'N' -> eigvals
! 'V' ->   " + eigvecs

      implicit none
      real*8, intent(inout) :: QHQ(:,:),S(:,:)
      real*8, intent(inout) :: avec(:)
      real*8, allocatable   :: WORK(:)
      integer     :: n,INFO,LWORK
      real*8      :: t1,t2
      character*1 :: COMPZ

      IF (.NOT. EIGEN_SETUP) call InitializeEigen()
      call CPU_TIME(t1)

      n=SIZE(avec)

!     Diagonalize QHQ, accounting for overlaps
      LWORK=3*n-1
      allocate(WORK(LWORK))
      call DSYGV(1,COMPZ,'U',n,QHQ,n,S,n,avec,WORK,LWORK,INFO)
      deallocate(WORK)

      if (INFO.ne.0) then
         write(*,*)'error in the diagonalization, info=',INFO
         call AbortWithError('error in DSYGV')
      endif

      call CPU_TIME(t2)
      eigen_time=eigen_time+t2-t1

      end subroutine SolveGenEigval

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveWithSVD(svals,A,U,VT)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SVD of m x n matrix, using call to LAPACK DGESVD

      implicit none
      real*8, allocatable, intent(inout) :: A(:,:)
      real*8, allocatable, intent(out)   :: U(:,:),VT(:,:)
      real*8, allocatable, intent(out)   :: svals(:)
      real*8, allocatable :: WORK(:)
      integer :: n,m,mmn,INFO,LWORK
      real*8  :: optdim(1)

      m=SIZE(A,1)
      n=SIZE(A,2)
      mmn=min(m,n)

      ALLOCATE(U(m,mmn),VT(mmn,n),svals(mmn))

!     Calculate optimal LWORK value for DGESVD
      CALL DGESVD('S','S',m,n,A,m,svals,U,m,VT,mmn,optdim,-1,INFO)
      LWORK=INT(optdim(1))
      ALLOCATE(WORK(LWORK))

!     Now calculate SVD...
      CALL DGESVD('S','S',m,n,A,m,svals,U,m,VT,mmn,WORK,LWORK,INFO)
      DEALLOCATE(WORK)

!     Error checking
      if (INFO.ne.0) then
         write(*,*) 'error in the DGESVD, info=',INFO
         call AbortWithError('error in DGESVD')
      endif

      end subroutine SolveWithSVD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine QRdecomp(Q,R)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets QR-decomposition of matrix. Matrix is given as R on input

      implicit none
      real*8, allocatable, intent(inout) :: R(:,:)
      real*8, allocatable, intent(out)   :: Q(:,:)
      real*8, allocatable :: RR(:,:),WORK(:),TAU(:)
      integer :: i,m,n,k,INFO,LWORK
      real*8  :: optdim(1)

      m=SIZE(R,1)
      n=SIZE(R,2)
      k=min(m,n)

      ALLOCATE(RR(m,n),TAU(k))
      RR=R
      DEALLOCATE(R)

!     Calculate optimal LWORK value for DGEQRF
      CALL DGEQRF(m,n,RR,m,TAU,optdim,-1,INFO)
      LWORK=INT(optdim(1))
      ALLOCATE(WORK(LWORK))

!     Get the QR-factorization as product of elementary reflectors
      CALL DGEQRF(m,n,RR,m,TAU,WORK,LWORK,INFO)
      if (INFO.ne.0) then
         write(*,*) 'error in the DGEQRF, info=',INFO
         call AbortWithError('error in DGEQRF')
      endif

      ALLOCATE(Q(m,k),R(k,n))
      R=0.d0
      DO i=1,k
         R(i,i:n)=RR(i,i:n)
      ENDDO

      Q(:,:)=RR(1:m,1:k)
      DEALLOCATE(RR)

!     Construct the Q matrix
      CALL DORGQR(m,k,k,Q,m,TAU,WORK,LWORK,INFO)
      if (INFO.ne.0) then
         write(*,*) 'error in the DORGQR, info=',INFO
         call AbortWithError('error in DORGQR')
      endif

      DEALLOCATE(TAU,WORK)

      end subroutine QRdecomp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine QRdecompALT(Q,R)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets QR-decomposition of matrix. Matrix is given as R on input

      implicit none
      real*8, intent(inout) :: R(:,:)
      real*8, allocatable, intent(out)   :: Q(:,:)
      real*8, allocatable :: WORK(:),TAU(:)
      integer :: i,j,m,n,k,INFO,LWORK
      real*8  :: optdim(1)

      m=SIZE(R,1)
      n=SIZE(R,2)
      k=min(m,n)

      ALLOCATE(TAU(k))

!     Calculate optimal LWORK value for DGEQRF
      CALL DGEQRF(m,n,R,m,TAU,optdim,-1,INFO)
      LWORK=INT(optdim(1))
      ALLOCATE(WORK(LWORK))

!     Get the QR-factorization as product of elementary reflectors
      CALL DGEQRF(m,n,R,m,TAU,WORK,LWORK,INFO)
      if (INFO.ne.0) then
         write(*,*) 'error in the DGEQRF, info=',INFO
         call AbortWithError('error in DGEQRF')
      endif

      ALLOCATE(Q(m,k))
      Q(:,:)=R(1:m,1:k)

!     Set the appropriate entries of R to zero
      DO i=2,k
         R(i,1:i-1)=0.d0
      ENDDO

!     Construct the Q matrix
      CALL DORGQR(m,k,k,Q,m,TAU,WORK,LWORK,INFO)
      if (INFO.ne.0) then
         write(*,*) 'error in the DORGQR, info=',INFO
         call AbortWithError('error in DORGQR')
      endif

      DEALLOCATE(TAU,WORK)

      end subroutine QRdecompALT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UnitaryTFM(U,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Performs unitary transformation U^T*M*U, returning the transformed M

      real*8, allocatable, intent(inout) :: M(:,:)
      real*8, intent(in)  :: U(:,:)
      real*8, allocatable :: Mtmp(:,:)
      integer :: ur,uc,mr,mc,i,j,k

      mr=SIZE(M,1)
      mc=SIZE(M,2)
      ur=SIZE(U,1)
      uc=SIZE(U,2)

      IF (mr.ne.mc) THEN
         write(*,*) 'rows(M) = ',mr,'; cols(M) = ',mc
         write(*,*) 'Error: matrix M must be square'
         call AbortWithError('Error in UnitaryTFM()')
      ENDIF
      IF (ur.ne.mc) THEN
         write(*,*) 'rows(U) = ',ur,'; cols(M) = ',mc
         write(*,*) 'Error: mismatch in dimensions between U and M'
         call AbortWithError('Error in UnitaryTFM()')
      ENDIF

      ALLOCATE(Mtmp(mr,uc))

!     Mtmp=M*U
      Mtmp=0.d0
      DO i=1,mr
         DO j=1,uc
            DO k=1,mc
               Mtmp(i,j)=Mtmp(i,j)+M(i,k)*U(k,j)
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(M)
      ALLOCATE(M(uc,uc))

!     (new) M=U^T*Mtmp
      M=0.d0
      DO i=1,uc
         DO j=1,uc
            DO k=1,mr
               M(i,j)=M(i,j)+U(k,i)*Mtmp(k,j)
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(Mtmp)

      end subroutine UnitaryTFM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixPseudoinverse(M,Mi)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes matrix pseudoinverse using SVD.

      implicit none
      real*8, intent(in) :: M(:,:)
      real*8, allocatable, intent(out) :: Mi(:,:)
      real*8, allocatable :: SIG(:,:),U(:,:),VT(:,:),svals(:)
      integer :: i,rm,cm,ns

      rm=SIZE(M,1)
      cm=SIZE(M,2)

!     Copy M to Mi since LAPACK destroys M
      ALLOCATE(Mi(rm,cm))
      Mi=M

!     SVD of M (as Mi)
      call SolveWithSVD(svals,Mi,U,VT)

!     Form SIGMA^(-1) matrix
      ns=SIZE(svals)
      ALLOCATE(SIG(ns,ns))
      SIG=0.d0
      DO i=1,ns
         SIG(i,i)=1.d0/svals(i)
      ENDDO

!     M^(-1) = V * SIGMA^(-1) * U^T
      call MatrixMult(SIG,.FALSE.,U,.TRUE.)
      call MatrixMult(VT,.TRUE.,SIG,.FALSE.,Mi)

      DEALLOCATE(U,VT,SIG,svals)

      end subroutine MatrixPseudoinverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE LINALG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
