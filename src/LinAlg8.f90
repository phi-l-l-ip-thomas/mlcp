!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE LINALG8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains wrappers for LAPACK subroutines

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MYACC
#if ACC_ENABLED
      USE OPENACC
      USE CUDAFOR
      USE CUBLAS_V2
      USE CUSOLVERDN
#endif

      implicit none
      real(kind=8), allocatable, private :: linsolv_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_LinAlg8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(linsolv_time(mpinodes))
      linsolv_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_LinAlg8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_LinAlg8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: i,ierr

      IF (.NOT. MODULE_SETUP) call Init_LinAlg8_Module()
      call Get_MPI_Timings('LinAlg8: solve LU-decomp',linsolv_time)
      MODULE_SETUP = .FALSE.
      deallocate(linsolv_time)

      end subroutine Dispose_LinAlg8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixMult8(M1,row1,col1,t1,M2,row2,col2,t2,M3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies M3 <- M1 x M2. Set tx (x={1,2}) to .TRUE. to do transpose

      implicit none
      real(kind=8), intent(in)    :: M1(:),M2(:)
      real(kind=8), intent(inout) :: M3(:)
      real(kind=8), parameter :: alpha=1.d0, beta=0.d0
      integer, intent(in) :: row1,col1,row2,col2,algo
      logical, intent(in) :: t1,t2
      integer :: r1,r2,c1,c2,ld1,ld2,tr1,tr2
      integer(4)   :: stat
      character(1) :: tN1,tN2

      ld1=row1
      ld2=row2

      r1=row1
      c1=col1
      tN1='N'
      tr1=0
      IF (t1) THEN
         r1=col1
         c1=row1
         tN1='T'
         tr1=1
      ENDIF
      r2=row2
      c2=col2
      tN2='N'
      tr2=0
      IF (t2) THEN
         r2=col2
         c2=row2
         tN2='T'
         tr2=1
      ENDIF

      IF (row1*col1 .ne. SIZE(M1)) THEN
         write(*,*) 'Bad length M1: size is ',SIZE(M1),&
         ' but dims passed are (',row1,' x ',col1,')'
         call AbortWithError('Error in MatrixMult8()')
      ENDIF

      IF (row2*col2 .ne. SIZE(M2)) THEN
         write(*,*) 'Bad length M2: size is ',SIZE(M2),&
         ' but dims passed are (',row2,' x ',col2,')'
         call AbortWithError('Error in MatrixMult8()')
      ENDIF

      IF (r1*c2 .ne. SIZE(M3)) THEN
         write(*,*) 'Bad length M3: size is ',SIZE(M3),&
         ' but required dims are (',r1,' x ',c2,')'
         call AbortWithError('Error in MatrixMult8()')
      ENDIF

      IF (c1.ne.r2) THEN
         write(*,*) 'cols(M1) = ',c1,'; rows(M2) = ',r2
         write(*,*) 'Error: mismatch in dimensions between M1 and M2'
         call AbortWithError('Error in MatrixMult8()')
      ENDIF

      select case(algo)
         case(0)
            call DGEMM(tN1,tN2,r1,c2,c1,alpha,M1,ld1,M2,ld2,beta,M3,r1)
         case(1)
#if ACC_ENABLED
            if (.not.acc_is_present(M1)) call &
               AbortWithError('MatrixMult8(): M1 not on device')
            if (.not.acc_is_present(M2)) call &
               AbortWithError('MatrixMult8(): M2 not on device')
            if (.not.acc_is_present(M3)) call &
               AbortWithError('MatrixMult8(): M3 not on device')

            !$acc data present(M1,M2,M3)
            !$acc host_data use_device(M1,M2,M3)
            stat=cublasDGEMM_v2(cublas_handle,tr1,tr2,r1,c2,c1,alpha,&
                                M1,ld1,M2,ld2,beta,M3,r1)
            if (stat.ne.0) then
               write(*,*) 'cublasDGEMM_v2 returned status ',stat
               call AbortWithError("MatrixMult8(): fail in cuBLAS call")
            endif
            !$acc end host_data
            !$acc end data
            !$acc wait
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("MatrixMult8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("MatrixMult8(): wrong algorithm")
      end select

      end subroutine MatrixMult8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CalcSVD8(A,m,n,U,svals,VT,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SVD of m x n matrix, using call to LAPACK DGESVD

      implicit none
      integer, intent(in) :: m,n,algo
      real(kind=8), intent(inout) :: A(:),U(:),VT(:)
      real(kind=8), intent(inout) :: svals(:)
      real(kind=8), allocatable   :: work(:),rwork(:)
      integer :: mmn,INFO,LWORK,err
      real(kind=8) :: optdim(1)
#if ACC_ENABLED
      integer, device :: cinfo
#endif

      mmn=min(m,n)

      IF (m*n .ne. SIZE(A)) THEN
         write(*,*) 'Bad length A: size is ',SIZE(A),&
         ' but dims passed are (',m,' x ',n,')'
         call AbortWithError('Error in CalcSVD8()')
      ENDIF

      IF (m*mmn .ne. SIZE(U)) THEN
         write(*,*) 'Bad length U: size is ',SIZE(U),&
         ' but dims must be (',m,' x ',mmn,')'
         call AbortWithError('Error in CalcSVD88()')
      ENDIF
      
      IF (mmn*n .ne. SIZE(VT)) THEN
         write(*,*) 'Bad length VT: size is ',SIZE(VT),&
         ' but dims must be (',mmn,' x ',n,')'
         call AbortWithError('Error in CalcSVD8()')
      ENDIF

      IF (mmn .ne. SIZE(svals)) THEN
         write(*,*) 'Bad length svals: size is ',SIZE(svals),&
         ' but dim must be (',mmn,')'
         call AbortWithError('Error in CalcSVD8()')
      ENDIF

      select case(algo)
         case(0)
!           Calculate optimal LWORK value for DGESVD
            CALL DGESVD('S','S',m,n,A,m,svals,U,m,VT,mmn,optdim,-1,INFO)
            LWORK=INT(optdim(1))
            ALLOCATE(WORK(LWORK))

!           Now calculate SVD...
            CALL DGESVD('S','S',m,n,A,m,svals,U,m,VT,mmn,WORK,LWORK,INFO)
            DEALLOCATE(WORK)

!           Error checking
            if (INFO.ne.0) then
               write(*,*) 'error in the DGESVD, info=',INFO
               call AbortWithError('CalcSVD8: fail in DGESVD call')
            endif

         case(1)
#if ACC_ENABLED
            if (.not.acc_is_present(A)) call &
               AbortWithError('CalcSVD8(): A not on device')
            if (.not.acc_is_present(U)) call &
               AbortWithError('CalcSVD8(): U not on device')
            if (.not.acc_is_present(svals)) call &
               AbortWithError('CalcSVD8(): svals not on device')
            if (.not.acc_is_present(VT)) call &
               AbortWithError('CalcSVD8(): VT not on device')

!           Calculate buffer size for SVD
            err=cusolverDnDgesvd_buffersize(cusolver_handle,m,n,LWORK)
            if (err.ne.0) then
               write(*,*) &
               'Error in cusolverDnDgesvd_buffersize; err = ',err
               call AbortWithError('Error in CalcSVD8()')
            endif

            allocate(work(LWORK),rwork(LWORK))
            !$acc data create(work,rwork) present(A,U,svals,VT)
            !$acc host_data use_device(A,U,svals,VT,work,rwork)
!           Calculate SVD
            err=cusolverDnDgesvd(cusolver_handle,&
               'S','S',m,n,A,m,svals,U,m,VT,mmn,work,LWORK,rwork,cinfo)
            if (err.ne.0) then
                write(*,*) &
               'Error in cusolverDnDgesvd; err = ',err
               call AbortWithError('Error in CalcSVD8_cusolve()')
            endif

            !$acc end host_data
            !$acc end data
            !$acc wait
            deallocate(work,rwork)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("CalcSVD8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("CalcSVD8(): wrong algorithm")
      end select

      end subroutine CalcSVD8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveWeightedLS8(LHS,RHS,mn,solver,valpen,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real(kind=8), intent(inout) :: LHS(:),RHS(:)
      integer, intent(in) :: mn,algo
      real(kind=8), intent(in) :: valpen
      character(len=*), intent(in) :: solver

      select case(algo)
         case(0)
            call SolveWeightedLS8_alg_cpu(LHS,RHS,mn,solver,valpen)
         case(1)
#if ACC_ENABLED
            call SolveWeightedLS8_alg_acc(LHS,RHS,mn,solver,valpen)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("SolveWeightedLS8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("SolveWeightedLS8(): wrong algorithm")
      end select

      end subroutine SolveWeightedLS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveWeightedLS8_alg_cpu(LHS,RHS,mn,solver,valpen)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves N rhs linear system (ALS case with weights) 

      implicit none
      real(kind=8), intent(inout) :: LHS(:),RHS(:)
      real(kind=8), allocatable :: rhs1(:)
      integer, intent(in) :: mn
      real(kind=8), intent(in)  :: valpen
      character(len=*), intent(in) :: solver
      integer :: i,j,k,rk,rk2,m,mrk,n,lst,lsf,r,r1

!     LHS is rk x rk * M (M=rows of F);
!     RHS is rk x M * N
      rk=SIZE(RHS)/mn
      rk2=rk*rk
      m=SIZE(LHS)/rk2
      n=SIZE(RHS)/rk/m
      mrk=m*rk

      ALLOCATE(rhs1(rk*n))

      lst=1
      lsf=rk2
      DO i=1,M
!        Extract the rhs for this value of i
         DO j=1,N
            DO k=1,rk
               r=(j-1)*Mrk+(i-1)*rk+k
               r1=(j-1)*rk+k
               rhs1(r1)=RHS(r)
            ENDDO
         ENDDO

         IF (TRIM(ADJUSTL(solver)).seq.'lu') THEN
            call SolveLinSysLU8(LHS(lst:lsf),rhs1,rk,n,valpen,0)
         ELSEIF (TRIM(ADJUSTL(solver)).seq.'svd') THEN
            call SolveLinSysSVD8(LHS(lst:lsf),rhs1,rk,n,valpen,0)
         ELSE
            write(*,*) "Linear solver algorithm '",TRIM(ADJUSTL(solver)),&
            "' not recognized; valid choices are 'LU' and 'SVD'"
            call AbortWithError("Error in SolveWeightedLS8_alg_cpu()")
         ENDIF

!        Copy the result back to RHS
         DO j=1,n
            DO k=1,rk
               r=(j-1)*mrk+(i-1)*rk+k
               r1=(j-1)*rk+k
               RHS(r)=rhs1(r1)
            ENDDO
         ENDDO
         lst=lst+rk2
         lsf=lsf+rk2
      ENDDO
      DEALLOCATE(rhs1)

      end subroutine SolveWeightedLS8_alg_cpu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveWeightedLS8_alg_acc(LHS,RHS,mn,solver,valpen)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves N rhs linear system (ALS case with weights) 

      implicit none
      real(kind=8), intent(inout) :: LHS(:),RHS(:)
      real(kind=8), allocatable :: lhs1(:),rhs1(:)
      integer, intent(in) :: mn
      real(kind=8), intent(in) :: valpen
      character(len=*), intent(in) :: solver
      integer :: i,j,k,rk,rk2,m,mrk,n,lst,lsf,r,r1

#if ACC_ENABLED
      if (.not.acc_is_present(LHS)) call AbortWithError(&
          'SolveWeightedLS8_alg_acc(): LHS not on device')
      if (.not.acc_is_present(RHS)) call AbortWithError(&
          'SolveWeightedLS8_alg_acc(): RHS not on device')

!     LHS is rk x rk * M (M=rows of F);
!     RHS is rk x M * N
      rk=SIZE(RHS)/mn
      rk2=rk*rk
      m=SIZE(LHS)/rk2
      n=SIZE(RHS)/rk/m
      mrk=m*rk

      ALLOCATE(lhs1(rk2),rhs1(rk*n))
      !$acc enter data create(lhs1,rhs1)
      lst=1
      lsf=rk2
      DO i=1,M
!        Extract the lhs for this value of i
         !$acc parallel loop gang vector private(r) async
         DO k=1,rk2
            r=lst+k-1
            lhs1(k)=LHS(r)
         ENDDO

!        Extract the rhs for this value of i
         !$acc parallel loop gang vector private(r,r1) collapse(2) async
         DO j=1,N
            DO k=1,rk
               r=(j-1)*Mrk+(i-1)*rk+k
               r1=(j-1)*rk+k
               rhs1(r1)=RHS(r)
            ENDDO
         ENDDO

         IF (TRIM(ADJUSTL(solver)).seq.'lu') THEN
            call SolveLinSysLU8(lhs1,rhs1,rk,n,valpen,1)
         ELSEIF (TRIM(ADJUSTL(solver)).seq.'svd') THEN
            call SolveLinSysSVD8(lhs1,rhs1,rk,n,valpen,1)
         ELSE
            write(*,*) "Linear solver algorithm '",TRIM(ADJUSTL(solver)),&
            "' not recognized; valid choices are 'LU' and 'SVD'"
            call AbortWithError("Error in SolveWeightedLS8_alg_acc()")
         ENDIF

!        Copy the result back to RHS
         !$acc parallel loop gang vector private(r,r1) collapse(2) async
         DO j=1,n
            DO k=1,rk
               r=(j-1)*mrk+(i-1)*rk+k
               r1=(j-1)*rk+k
               RHS(r)=rhs1(r1)
            ENDDO
         ENDDO

         lst=lst+rk2
         lsf=lsf+rk2
      ENDDO
      !$acc exit data delete(lhs1,rhs1)
      !$acc wait
      DEALLOCATE(lhs1,rhs1)

#endif

      end subroutine SolveWeightedLS8_alg_acc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSys8(A,B,m,n,solver,reg,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via DGETRF for LU-decomp + DGETRS, where
! A is an (m x m) matrix, and B,X are a group of n length-m vectors.
! Solution X replaces B, and A is destroyed by this routine
! This is the wrapper to call the {CPU,GPU} routines

      implicit none
      integer, intent(in) :: m,n,algo
      real(kind=8), intent(inout) :: A(:)
      real(kind=8), intent(inout) :: B(:)
      real(kind=8), intent(in)    :: reg
      character(len=*), intent(in) :: solver

      IF (TRIM(ADJUSTL(solver)).seq.'lu') THEN
         call SolveLinSysLU8(A,B,m,n,reg,algo)
      ELSEIF (TRIM(ADJUSTL(solver)).seq.'svd') THEN
         call SolveLinSysSVD8(A,B,m,n,reg,algo)
      ELSE
         write(*,*) "Linear solver algorithm '",TRIM(ADJUSTL(solver)),&
         "' not recognized; valid choices are 'LU' and 'SVD'"
         call AbortWithError("Error in SolveLinSys8()")
      ENDIF

      end subroutine SolveLinSys8 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysSVD8(A,B,m,n,reg,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via DGETRF for LU-decomp + DGETRS, where
! A is an (m x m) matrix, and B,X are a group of n length-m vectors.
! Solution X replaces B, and A is destroyed by this routine
! This is the wrapper to call the {CPU,GPU} routines

      implicit none
      integer, intent(in) :: m,n,algo
      real(kind=8), intent(inout) :: A(:)
      real(kind=8), intent(inout) :: B(:)
      real(kind=8), intent(in)    :: reg
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_LinAlg8_Module()

      IF (reg.lt.0.d0 .or. reg.gt.1.d0) THEN
         write(*,*) 'reg ',reg,' must be in range [0 <= reg <= 1]'
         call AbortWithError('SolveLinSysLU():reg < 0')
      ENDIF

      IF (SIZE(A).ne.m*m) THEN
         write(*,*) 'A has size = ',SIZE(A),' but must be [',m,&
                    ' x ',m,'] = ',m*m
         call AbortWithError('SolveLinSysLU8(): wrong size A')
      ENDIF

      IF (SIZE(B).ne.m*n) THEN
         write(*,*) 'B has size = ',SIZE(B),' but must be [',m,&
                    ' x ',n,'] = ',m*n
         call AbortWithError('SolveLinSysLU8(): wrong size B')
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('SolveLinSysSVD8')

      select case(algo)
         case(0)
            call SolveLinSysSVD8_cpu(A,B,m,n,reg)
         case(1)
#if ACC_ENABLED
            call SolveLinSysSVD8_cusolve(A,B,m,n,reg)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("SolveLinSysSVD8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("SolveLinSysSVD8(): wrong algorithm")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      linsolv_time(mpirank+1)=linsolv_time(mpirank+1)+ti2-ti1

      end subroutine SolveLinSysSVD8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysSVD8_cpu(A,B,m,n,reg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via SVD-calculated Moore-Penrose
! pseudoinverse, where A is an (m x m) matrix, and B,X are a group of n 
! length-m vectors. Solution X replaces B; A is destroyed by this routine

      implicit none
      integer, parameter :: algo=0
      real(kind=8), intent(in) :: reg
      real(kind=8), intent(inout) :: A(:),B(:)
      real(kind=8), allocatable :: U(:),S(:),VT(:),W(:),svals(:)
      integer :: i,j,m,n
      real(kind=8)  :: thresh

      ALLOCATE(U(m*m),S(m*m),VT(m*m),W(m*n),svals(m))

!     SVD of A (destroys A)
      call CalcSVD8(A,m,m,U,svals,VT,algo)

!     Form SIGMA^(-1) matrix, keeping only values exceeding threshold
      thresh=reg*svals(1)
      S(:)=0.d0
      DO i=1,m
         if (svals(i).lt.thresh) exit
         j=(i-1)*m+i
         S(j)=1.d0/svals(i)
      ENDDO

!     Calculate  M^(-1)*b = [V * SIGMA^(-1) * U^T] * b

!     A <- V * SIGMA^(-1)
      call MatrixMult8(VT,m,m,.TRUE.,S,m,m,.FALSE.,A,algo)
!     VT <- U^T * B
      call MatrixMult8(U,m,m,.TRUE.,B,m,n,.FALSE.,W,algo)
!     B <- [V * SIGMA^(-1)] * [U^T * b]
      call MatrixMult8(A,m,m,.FALSE.,W,m,n,.FALSE.,B,algo)

      DEALLOCATE(U,S,VT,W,svals)

      end subroutine SolveLinSysSVD8_cpu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysSVD8_cusolve(A,B,m,n,reg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via SVD-calculated Moore-Penrose
! pseudoinverse, where A is an (m x m) matrix, and B,X are a group of n 
! length-m vectors. Solution X replaces B; A is destroyed by this routine

      implicit none
      integer, parameter :: algo=1
      real(kind=8), intent(in) :: reg
      real(kind=8), intent(inout) :: A(:),B(:)
      real(kind=8), allocatable :: U(:),S(:),VT(:),W(:),svals(:)
      integer :: i,j,m,n,m2

#if ACC_ENABLED
      m2=m*m
      ALLOCATE(U(m2),S(m2),VT(m2),W(m*n),svals(m))

      !$acc data present(A,B) create(U,S,VT,W,svals)

!     SVD of A (destroys A)
      call CalcSVD8(A,m,m,U,svals,VT,algo)

!     Form SIGMA^(-1) matrix, keeping only values exceeding threshold

      !$acc parallel loop gang vector async
      do i=1,m2
         S(i)=0.d0
      enddo

      !$acc parallel loop gang vector private(j) async
      do i=1,m
         if (svals(i).gt.(reg*svals(1))) then
            j=(i-1)*m+i
            S(j)=1.d0/svals(i)
         endif
      enddo

!     Calculate  M^(-1)*b = [V * SIGMA^(-1) * U^T] * b

!     A <- V * SIGMA^(-1)
      call MatrixMult8(VT,m,m,.TRUE.,S,m,m,.FALSE.,A,algo)
!     VT <- U^T * B
      call MatrixMult8(U,m,m,.TRUE.,B,m,n,.FALSE.,W,algo)
!     B <- [V * SIGMA^(-1)] * [U^T * b]
      call MatrixMult8(A,m,m,.FALSE.,W,m,n,.FALSE.,B,algo)

      !$acc end data
      !$acc wait
      DEALLOCATE(U,S,VT,W,svals)
#endif

      end subroutine SolveLinSysSVD8_cusolve

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysLU8(A,B,m,n,reg,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via DGETRF for LU-decomp + DGETRS, where
! A is an (m x m) matrix, and B,X are a group of n length-m vectors.
! Solution X replaces B, and A is destroyed by this routine
! This is the wrapper to call the {CPU,GPU} routines

      implicit none
      integer, intent(in) :: m,n,algo
      real(kind=8), intent(inout) :: A(:)
      real(kind=8), intent(inout) :: B(:)
      real(kind=8), intent(in)    :: reg
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_LinAlg8_Module()

      IF (reg.lt.0.d0) THEN
         write(*,*) 'reg = ',reg,' must be nonnegative'
         call AbortWithError('SolveLinSysLU():reg < 0')
      ENDIF

      IF (SIZE(A).ne.m*m) THEN
         write(*,*) 'A has size = ',SIZE(A),' but must be [',m,&
                    ' x ',m,'] = ',m*m
         call AbortWithError('SolveLinSysLU8(): wrong size A')
      ENDIF

      IF (SIZE(B).ne.m*n) THEN
         write(*,*) 'B has size = ',SIZE(B),' but must be [',m,&
                    ' x ',n,'] = ',m*n
         call AbortWithError('SolveLinSysLU8(): wrong size B')
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('SolveLinSysLU8')

      select case(algo)
         case(0)
            call SolveLinSysLU8_cpu(A,B,m,n,reg)
         case(1)
#if ACC_ENABLED
            call SolveLinSysLU8_cusolve(A,B,m,n,reg)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("SolveLinSysLU8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("SolveLinSysLU8(): wrong algorithm")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      linsolv_time(mpirank+1)=linsolv_time(mpirank+1)+ti2-ti1

      end subroutine SolveLinSysLU8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysLU8_cpu(A,B,m,n,reg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via DGETRF for LU-decomp + DGETRS, where
! A is an (m x m) matrix, and B,X are a group of n length-m vectors.
! Solution X replaces B, and A is destroyed (CPU-LAPACK version)

      implicit none
      real(kind=8), intent(inout) :: A(:)
      real(kind=8), intent(inout) :: B(:)
      real(kind=8), intent(in)    :: reg
      integer, intent(in)  :: m,n
      integer, allocatable :: IPV(:)
      integer :: i,j,info

!     Add regularization penalty to avoid ill-conditioning
      DO i=1,m
         j=(i-1)*m + i
         A(j)=A(j)+reg
      ENDDO

      allocate(IPV(m))
      call dgetrf(m,m,A,m,IPV,info)
      call dgetrs('N',m,n,A,m,IPV,B,m,info)
      deallocate(IPV)

      end subroutine SolveLinSysLU8_cpu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLinSysLU8_cusolve(A,B,m,n,reg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system A*X = B, via DGETRF for LU-decomp + DGETRS, where
! A is an (m x m) matrix, and B,X are a group of n length-m vectors.
! Solution X replaces B, and A is destroyed (GPU-CUSOLVER version)

      implicit none
      real(kind=8), intent(inout) :: A(:)
      real(kind=8), intent(inout) :: B(:)
      real(kind=8), intent(in)    :: reg
      integer, intent(in)  :: m,n
      integer :: i,j,lwork,err
      integer, allocatable :: ipv(:)
      integer, parameter   :: trans=0
      real(kind=8), allocatable :: work(:)
#if ACC_ENABLED
      integer, device :: info
#endif

#if ACC_ENABLED
      if (.not.acc_is_present(A)) call &
         AbortWithError('SolveLinSysLU8_cusolve(): A not on device')
      if (.not.acc_is_present(B)) call &
         AbortWithError('SolveLinSysLU8_cusolve(): B not on device')

      !$acc data present(A,B)
!     Add regularization penalty to avoid ill-conditioning
      !$acc parallel loop gang vector private(j) async
      DO i=1,m
         j=(i-1)*m + i
         A(j)=A(j)+reg
      ENDDO

      !$acc host_data use_device(A)
      err = cusolverDnDgetrf_buffersize(cusolver_handle,&
              m,m,A,m,lwork)
      !$acc end host_data
      if (err.ne.0) then
          write(*,*) &
         'Error in cusolverDnDgetrf_buffersize; err = ',err
         call AbortWithError('Error in SolveLinSysLU8_cusolve()')
      endif

      allocate(work(lwork),ipv(m))
      !$acc data create(work,IPV)
      !$acc host_data use_device(A,B,work,IPV)
      err = cusolverDnDgetrf(cusolver_handle,&
              m,m,A,m,work,ipv,info)
      if (err.ne.0) then
          write(*,*) &
         'Error in cusolverDnDgetrf; err = ',err
         call AbortWithError('Error in SolveLinSysLU8_cusolve()')
      endif

      err = cusolverDnDgetrs(cusolver_handle,&
              trans,m,n,A,m,ipv,B,m,info)
      if (err.ne.0) then
          write(*,*) & 
         'Error in cusolverDnDgetrs; err = ',err
         call AbortWithError('Error in SolveLinSysLU8_cusolve()')
      endif
      !$acc end host_data
      !$acc end data
      deallocate(work,ipv)
      !$acc end data
      !$acc wait
#endif

      end subroutine SolveLinSysLU8_cusolve

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE LINALG8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
