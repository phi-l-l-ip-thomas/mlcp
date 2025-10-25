!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPMM8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MYACC
      USE CPr8
#if ACC_ENABLED
      USE OPENACC
      USE CUDAFOR
      USE CUBLAS_V2
      USE CUTENSOR_v2
#endif

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE MM8
         TYPE (CP8) :: I1,I2
         integer, allocatable, dimension(:) :: r1,r2,c1,c2,ld1,ld2,sd1,sd2
         integer :: ndof,rk1,rks1,rk2,rks2,rk3,ish1,ish2,algo
         real(kind=8) :: Esh1,Esh2,alpha1,alpha2
         logical :: t1,t2,t3
         logical :: setup = .FALSE.
         character(1) :: tN1,tN2
#if ACC_ENABLED
         type(cutensorDataType) :: dataType
         type(cutensorComputeDescriptor) :: descCompute
         type(cutensortensordescriptor), allocatable :: descM1(:),descM2(:),descM3(:) ! cutensor 2.x
         type(cutensortensordescriptor) :: descC1,descC2,descC3 ! cutensor 2.x
         type(cutensoroperationdescriptor), allocatable :: op_desc(:) ! cutensor 2.x
         type(cutensoroperationdescriptor) :: op_descC ! cutensor 2.x
         type(cutensorplanpreference) :: plan_pref ! cutensor 2.x
         type(cutensorplan), allocatable :: mm_plan(:) ! cutensor 2.x
         type(cutensorplan) :: mm_planC ! cutensor 2.x
         integer(8), allocatable :: work_size(:) ! cutensor 1.x and 2.x
         integer(8) :: work_sizeC ! cutensor 1.x and 2.x
#endif
         CONTAINS
            PROCEDURE :: new => New_MM8
            PROCEDURE :: flush => Flush_MM8
      END TYPE MM8

      INTERFACE CPMM_CP8
         MODULE PROCEDURE CPMM_CP8_simple_wrapper 
         MODULE PROCEDURE CPMM_CP8_noshifts_wrapper
         MODULE PROCEDURE CPMM_CP8_allmodes_wrapper
         MODULE PROCEDURE CPMM_CP8_1mode_wrapper
         MODULE PROCEDURE CPMM_CP8_general_wrapper
         MODULE PROCEDURE CPMM_CP8_allmodes
         MODULE PROCEDURE CPMM_CP8_1mode
         MODULE PROCEDURE CPMM_CP8_general
      END INTERFACE CPMM_CP8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_CPMM8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_CPMM8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_CPMM8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_CPMM8_Module()
      call Get_MPI_Timings('CPMM8: matrix multiplication',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_CPMM8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine New_MM8(U,M1,ish1,Esh1,t1,&
                           M2,ish2,Esh2,t2,&
                           M3,t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes MM8 type

      implicit none
      CLASS (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      integer, intent(in)       :: ish1,ish2,algo
      real(kind=8), intent(in)  :: Esh1,Esh2
      logical, intent(in)  :: t1,t2,t3
      integer :: d
     
!     Set parameters: # modes, ranks
      U%algo=algo
      U%t3=t3
      U%ndof=M1%D()

      U%ish1=ish1
      U%Esh1=Esh1
      U%t1=t1
      U%rk1=M1%R()
      U%rks1=M1%R()
      if (ish1.ne.0) U%rks1=U%rks1+1
      U%alpha1=1.0
      if (ish1.lt.0) U%alpha1=-U%alpha1
      U%tN1='N'
      if (t1) U%tN1='T'

      U%ish2=ish2
      U%Esh2=Esh2
      U%t2=t2
      U%rk2=M2%R()
      U%rks2=M2%R()
      if (ish2.ne.0) U%rks2=U%rks2+1
      U%alpha2=1.0
      if (ish2.lt.0) U%alpha2=-U%alpha2
      U%rk3=U%rks1*U%rks2
      U%tN2='N'
      if (t2) U%tN2='T'

!     Error checking
      IF (U%setup) &
         call AbortWithError('New_MM8(): U already set up')


      IF (M1%D().lt.1 .or. M2%D().lt.1) THEN
         write(*,*) 'M1, M2 have ',M1%D(),', ',M2%D(),&
                    ' modes, respectively; must both be > 0'
         call AbortWithError('New_MM8(): # modes must be > 0')
      ENDIF

      IF (M2%D().ne.M1%D()) THEN
         write(*,*) 'M1 has ',M1%D(),&
         ' modes but M2 has ',M1%D(),' modes'
         call AbortWithError('New_MM8(): M1, M2 mode mismatch')
      ENDIF

      ALLOCATE(U%r1(U%ndof),U%c1(U%ndof),U%ld1(U%ndof),U%sd1(U%ndof))
      ALLOCATE(U%r2(U%ndof),U%c2(U%ndof),U%ld2(U%ndof),U%sd2(U%ndof))

!     Fill the BLAS dimensions
      DO d=1,U%ndof
         U%r1(d)=M1%rows(d)
         U%c1(d)=M1%cols(d)
         U%ld1(d)=M1%rows(d)
         U%sd1(d)=M1%cols(d)
         U%r2(d)=M2%rows(d)
         U%c2(d)=M2%cols(d)
         U%ld2(d)=M2%rows(d)
         U%sd2(d)=M2%cols(d)
         if (t1) call swap(U%r1(d),U%c1(d))
         if (t2) call swap(U%r2(d),U%c2(d))

         IF (U%c1(d).ne.U%r2(d)) THEN
            write(*,*) 'Mode ',d,': M1 cols (',U%c1(d),&
                       '), M2 rows (',U%r2(d),') must match'
            call AbortWithError('New_MM8(): bad matrix sizes')
         ENDIF
      ENDDO

      IF (M3%R().gt.0) THEN
         IF (M3%R().ne.U%rk3) THEN
            write(*,*) '(existing) M3 has rank ',M3%R(),&
            ', but must be ',U%rk3
            call AbortWithError('New_MM8(): wrong M3 rank')
         ENDIF

         IF (M3%D().ne.U%ndof) THEN
            write(*,*) '(existing) M3 has ',M3%D(),&
            ' modes but must have ',U%ndof,' modes'
            call AbortWithError('New_MM8(): wrong M3 ndof')
         ENDIF

         DO d=1,U%ndof
            IF (M3%rows(d).ne.U%r1(d) .or. M3%cols(d).ne.U%c2(d)) THEN
               write(*,*) 'Mode ',d,': M3 has dimensions (',&
                          M3%rows(d),' x ',M3%cols(d),&
                         ') but must be (',U%r1(d),' x ',U%c2(d),')'
               call AbortWithError('NewMM8(): wrong M3 dimensions')
            ENDIF
         ENDDO
      ENDIF

      if (ish1.ne.0) then
         call U%I1%cloneid(M1)
         call U%I1%mult(-U%alpha1*Esh1)
#if ACC_ENABLED
         if (algo.eq.1) call U%I1%copyintodevice()
#endif
      endif

      if (ish2.ne.0) then
         call U%I2%cloneid(M2)
         call U%I2%mult(-U%alpha2*Esh2)
#if ACC_ENABLED
         if (algo.eq.1) call U%I2%copyintodevice()
#endif
      endif

      if (algo.eq.1) call Init_cutensor_v2_MM8(U)

      U%setup=.TRUE.

      end subroutine New_MM8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_cutensor_v2_MM8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes MM8 type

      implicit none
      CLASS (MM8) :: U
      integer(8) :: m1_ext(3),m2_ext(3),m3_ext(4)
      integer(8) :: c1_ext(1),c2_ext(1),c3_ext(2)
      integer(8) :: m1_str(3),m2_str(3),m3_str(4)
      integer(8) :: c1_str(1),c2_str(1),c3_str(2)
      integer    :: m1_ord(3),m2_ord(3),m3_ord(4)
      integer    :: c1_ord(1),c2_ord(1),c3_ord(2)
      integer    :: d,ndof

#if ACC_ENABLED
      ndof=U%ndof

      U%dataType = CUTENSOR_R_64F
      U%descCompute = CUTENSOR_COMPUTE_DESC_64F 

      ! Create tensor descriptors for base of M1,M2,M3
       ALLOCATE(U%descM1(ndof),U%descM2(ndof),U%descM3(ndof))
       DO d=1,ndof
          m1_ext=(/U%ld1(d),U%sd1(d),U%rks1/)
          m1_str=(/1,U%ld1(d),U%ld1(d)*U%sd1(d)/)
          m2_ext=(/U%ld2(d),U%sd2(d),U%rks2/)
          m2_str=(/1,U%ld2(d),U%ld2(d)*U%sd2(d)/)
          if (U%t3) then
             m3_ext=(/U%r1(d),U%c2(d),U%rks2,U%rks1/)
             m3_str=(/1,U%r1(d),U%r1(d)*U%c2(d),U%r1(d)*U%c2(d)*U%rks2/)
          else
             m3_ext=(/U%r1(d),U%c2(d),U%rks1,U%rks2/)
             m3_str=(/1,U%r1(d),U%r1(d)*U%c2(d),U%r1(d)*U%c2(d)*U%rks1/)
          endif
          cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
               U%descM1(d),3,m1_ext,m1_str,U%dataType,cutensor_align)
          if (cutensor_status%stat .ne.0) then
             write(*,*) 'cutensorCreateTensorDescriptor, M1 base mode',d,&
             '; exit = ',cutensor_status%stat
             call AbortWithError('Error in Init_cutensor_v2_MM8()')
          endif
          cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
               U%descM2(d),3,m2_ext,m2_str,U%dataType,cutensor_align)
          if (cutensor_status%stat .ne.0) then
             write(*,*) 'cutensorCreateTensorDescriptor, M2 base mode',d,&
             '; exit = ',cutensor_status%stat
             call AbortWithError('Error in Init_cutensor_v2_MM8()')
          endif
          cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
               U%descM3(d),4,m3_ext,m3_str,U%dataType,cutensor_align)
          if (cutensor_status%stat .ne.0) then
             write(*,*) 'cutensorCreateTensorDescriptor, M3 base mode',d,&
             '; exit = ',cutensor_status%stat
             call AbortWithError('Error in Init_cutensor_v2_MM8()')
          endif
      ENDDO

!     Create tensor descriptors for coef of M1,M2,M3
      c1_ext=(/U%rks1/)
      c1_str=(/1/)
      c2_ext=(/U%rks2/)
      c2_str=(/1/)
      if (U%t3) then
         c3_ext=(/U%rks2,U%rks1/)
         c3_str=(/1,U%rks2/)
      else
         c3_ext=(/U%rks1,U%rks2/)
         c3_str=(/1,U%rks1/)
      endif
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC1,1,c1_ext,c1_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M1 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC2,1,c2_ext,c2_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M2 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC3,2,c3_ext,c3_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M3 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif

      ! Create operator (contraction) descriptor for base of M1,M2,M3
      ALLOCATE(U%op_desc(ndof))
      m1_ord=(/1,2,3/)
      if (U%t1) call swap(m1_ord(1),m1_ord(2))
      m2_ord=(/2,4,5/)
      if (U%t2) call swap(m2_ord(1),m2_ord(2))
      if (U%t3) then
         m3_ord=(/1,4,5,3/)
      else
         m3_ord=(/1,4,3,5/)
      endif
      DO d=1,ndof
         cutensor_status=cutensorcreatecontraction(cutensor_handle,U%op_desc(d),&
                         U%descM1(d),m1_ord,CUTENSOR_OP_IDENTITY,&
                         U%descM2(d),m2_ord,CUTENSOR_OP_IDENTITY,&
                         U%descM3(d),m3_ord,CUTENSOR_OP_IDENTITY,&
                         U%descM3(d),m3_ord,U%descCompute)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreatecontraction, base mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_MM8()')
         endif
      ENDDO

      ! Create operator (contraction) descriptor for coef of M1,M2,M3
      c1_ord=(/1/)
      c2_ord=(/2/)
      if (U%t3) then
         c3_ord=(/2,1/)
      else
         c3_ord=(/1,2/)
      endif
      cutensor_status=cutensorcreatecontraction(cutensor_handle,U%op_descC,&
                      U%descC1,c1_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC2,c2_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC3,c3_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC3,c3_ord,U%descCompute)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreatecontraction, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif

      ! Create plan preference
      cutensor_status=cutensorcreateplanpreference(cutensor_handle,U%plan_pref,&
                      CUTENSOR_ALGO_DEFAULT,CUTENSOR_JIT_MODE_NONE)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreateplanpreference; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif

      ! Estimate workspace size for base contraction
      ALLOCATE(U%work_size(ndof))
      do d=1,ndof
         cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%op_desc(d),&
                         U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%work_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorestimateworkspacesize, base mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_MM8()')
         endif
      enddo

      ! Estimate workspace size for coef contraction
      cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%op_descC,&
                      U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%work_sizeC)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorestimateworkspacesize, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif

      ! Create the plan for base contraction
      ALLOCATE(U%mm_plan(ndof))
      do d=1,ndof
        cutensor_status=cutensorcreateplan(cutensor_handle,U%mm_plan(d),&
                        U%op_desc(d),U%plan_pref,U%work_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreateplan, base mode',d,&
                       '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
         endif
      enddo

      ! Create the plan for coef contraction
       cutensor_status=cutensorcreateplan(cutensor_handle,U%mm_planC,&
                       U%op_descC,U%plan_pref,U%work_sizeC)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreateplan, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_MM8()')
      endif
#endif

      end subroutine Init_cutensor_v2_MM8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_MM8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (MM8) :: U
      integer :: d,ndof

      ndof=U%ndof

      IF (.not.U%setup) &
         call AbortWithError('Flush_MM8(): U not setup')

      IF (ALLOCATED(U%r1))  DEALLOCATE(U%r1)
      IF (ALLOCATED(U%c1))  DEALLOCATE(U%c1)
      IF (ALLOCATED(U%ld1)) DEALLOCATE(U%ld1)
      IF (ALLOCATED(U%sd1)) DEALLOCATE(U%sd1)
      IF (ALLOCATED(U%r2))  DEALLOCATE(U%r2)
      IF (ALLOCATED(U%c2))  DEALLOCATE(U%c2)
      IF (ALLOCATED(U%ld2)) DEALLOCATE(U%ld2)
      IF (ALLOCATED(U%sd2)) DEALLOCATE(U%sd2)

      IF (U%ish1.ne.0) call U%I1%flush()
      IF (U%ish2.ne.0) call U%I2%flush()

#if ACC_ENABLED
      if (U%algo.eq.1) then

         ! Destroy plan for coef and base
         cutensor_status = cutensorDestroyPlan(U%mm_planC)
         DO d=1,ndof
            cutensor_status = cutensorDestroyPlan(U%mm_plan(d))
         ENDDO

         ! Destroy plan preference
         cutensor_status = cutensorDestroyPlanPreference(U%plan_pref)

         ! Destroy operator (contraction) descriptor for coef and base
         cutensor_status = cutensorDestroyOperationDescriptor(U%op_descC)
         DO d=1,ndof
            cutensor_status = cutensorDestroyOperationDescriptor(U%op_desc(d))
         ENDDO

         ! Destroy tensor descriptors for coef and base of M1,M2,M3
         cutensor_status = cutensorDestroyTensorDescriptor(U%descC1)
         cutensor_status = cutensorDestroyTensorDescriptor(U%descC2)
         cutensor_status = cutensorDestroyTensorDescriptor(U%descC3)
         DO d=1,ndof
            cutensor_status = cutensorDestroyTensorDescriptor(U%descM1(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descM2(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descM3(d))
         ENDDO

         IF (ALLOCATED(U%mm_plan)) DEALLOCATE(U%mm_plan)
         IF (ALLOCATED(U%work_size)) DEALLOCATE(U%work_size)
         IF (ALLOCATED(U%op_desc)) DEALLOCATE(U%op_desc)
         IF (ALLOCATED(U%descM1)) DEALLOCATE(U%descM1)
         IF (ALLOCATED(U%descM2)) DEALLOCATE(U%descM2)
         IF (ALLOCATED(U%descM3)) DEALLOCATE(U%descM3)
      endif
#endif

      U%setup=.FALSE.

      end subroutine Flush_MM8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Precopy_base_CP8(V,M,I,d,ish)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copy base data from (CP8) M into temporary array V before contraction

      implicit none
      TYPE (CP8), INTENT(IN) :: M,I
      real(kind=8), intent(inout) :: V(:)
      integer, intent(in) :: d,ish
      integer :: r,shm,shi,nm,rk,rknm

#if ACC_ENABLED
      shm=M%look(1,d)
      nm=M%M(d)*M%N(d)
      rk=M%R()
      rknm=rk*nm

      !$acc data present(V,M,M%base)
      if (d.eq.1 .and. ish.lt.0) then
         !$acc parallel loop gang vector async
         do r=1,rknm
            V(r)=-M%base(shm+r)
         enddo
      else
         !$acc parallel loop gang vector async
         do r=1,rknm
            V(r)=M%base(shm+r)
         enddo
      endif
      !$acc end data

      if (ish.ne.0) then
         shi=I%look(1,d)
         !$acc data present(V,I,I%base)
         !$acc parallel loop gang vector async
         do r=1,nm
            V(rknm+r)=I%base(shi+r)
         enddo
         !$acc end data
      endif
#endif
   
      end subroutine Precopy_base_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Precopy_coef_CP8(V,M,ish,Esh)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copy coef data from (CP8) M into temporary array V before contraction

      implicit none
      TYPE (CP8), INTENT(IN) :: M
      real(kind=8), intent(inout) :: V(:)
      integer, intent(in) :: ish
      real(kind=8), intent(in) :: Esh
      integer :: r,rk

#if ACC_ENABLED
      rk=M%R()

      !$acc data present(V,M,M%coef)
      !$acc parallel loop gang vector async
      do r=1,rk
         V(r)=M%coef(r)
      enddo

      if (ish.ne.0) then
         !$acc parallel loop gang vector async
         do r=1,1
            V(rk+r)=abs(Esh)
         enddo
      endif
      !$acc end data
#endif

      end subroutine Precopy_coef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_simple_wrapper(M1,M2,M3,t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version automates the process of
! creating/using/flushing the MM8 type, all-modes-no-shifts-no-transpose

      implicit none
      TYPE (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      logical, intent(in) :: t3
      integer, intent(in) :: algo

      call U%new(M1,0,0.0,.FALSE.,M2,0,0.0,.FALSE.,M3,t3,algo)
      call CPMM_CP8(U,M1,M2,M3)
      call U%flush()

      end subroutine CPMM_CP8_simple_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_noshifts_wrapper(M1,t1,M2,t2,M3,t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version automates the process of
! creating/using/flushing the MM8 type, all-modes-no-shifts

      implicit none
      TYPE (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      integer, intent(in)  :: algo
      logical, INTENT(in)  :: t1,t2,t3

      call U%new(M1,0,0.0,t1,M2,0,0.0,t2,M3,t3,algo)
      call CPMM_CP8(U,M1,M2,M3)
      call U%flush()

      end subroutine CPMM_CP8_noshifts_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_allmodes_wrapper(M1,ish1,Esh1,t1,&
                                           M2,ish2,Esh2,t2,M3,&
                                           t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version automates the process of
! creating/using/flushing the MM8 type, for all modes

      implicit none
      TYPE (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      INTEGER, INTENT(IN)       :: ish1,ish2,algo
      REAL(KIND=8), INTENT(IN)  :: Esh1,Esh2
      LOGICAL, INTENT(IN)  :: t1,t2,t3

      call U%new(M1,ish1,Esh1,t1,M2,ish2,Esh2,t2,M3,t3,algo)
      call CPMM_CP8(U,M1,M2,M3)
      call U%flush()

      end subroutine CPMM_CP8_allmodes_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_1mode_wrapper(M1,ish1,Esh1,t1,&
                                M2,ish2,Esh2,t2,M3,imode,t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version automates the process of
! creating/using/flushing the MM8 type, for selected mode

      implicit none
      TYPE (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      INTEGER, INTENT(IN)       :: ish1,ish2,imode,algo
      REAL(KIND=8), INTENT(IN)  :: Esh1,Esh2
      LOGICAL, INTENT(IN)  :: t1,t2,t3

      call U%new(M1,ish1,Esh1,t1,M2,ish2,Esh2,t2,M3,t3,algo)
      call CPMM_CP8(U,M1,M2,M3,imode)
      call U%flush()

      end subroutine CPMM_CP8_1mode_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_general_wrapper(M1,ish1,Esh1,t1,&
                                          M2,ish2,Esh2,t2,&
                                          M3,modes,t3,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version automates the process of
! creating/using/flushing the MM8 type, general routine

      implicit none
      TYPE (MM8) :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      integer, intent(in)       :: ish1,ish2,algo
      real(kind=8), intent(in)  :: Esh1,Esh2
      logical, intent(in)  :: t1,t2,t3,modes(:)

      call U%new(M1,ish1,Esh1,t1,M2,ish2,Esh2,t2,M3,t3,algo)
      call CPMM_CP8(U,M1,M2,M3,modes)
      call U%flush()

      end subroutine CPMM_CP8_general_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_allmodes(U,M1,M2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version uses existing MM8 type to do the
! operation, for all modes

      implicit none
      TYPE (MM8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      logical, allocatable      :: modes(:)

      allocate(modes(U%ndof))
      modes(:)=.TRUE.
      call CPMM_CP8(U,M1,M2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_CP8_allmodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_1mode(U,M1,M2,M3,imode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version uses existing MM8 type to do the
! operation, for selected mode

      implicit none
      TYPE (MM8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      integer, intent(in)       :: imode
      logical, allocatable      :: modes(:)

!     Create and fill matrix-multiply object
      if (imode.lt.1 .or. imode.gt.U%ndof) then
         write(*,*) 'Mode ',imode,' out of range: [',1,',',U%ndof,']'
         call AbortWithError('CPMM_CP8_1mode(): wrong imode')
      endif

      allocate(modes(U%ndof))
      modes(:)=.FALSE.
      modes(imode)=.TRUE.
      call CPMM_CP8(U,M1,M2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_CP8_1mode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_CP8_general(U,M1,M2,M3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix multiply routine; this version uses existing MM8 type to do the
! operation, general version

      implicit none
      TYPE (MM8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      logical, intent(in)       :: modes(:)
      real(kind=8) :: ti1,ti2
      integer :: d

      IF (.NOT. MODULE_SETUP) call Init_CPMM8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('CPMM_CP8(): U not setup')

      IF (SIZE(modes).ne.U%ndof) THEN
         write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',U%ndof
         call AbortWithError('CPMM_CP8(): bad size for modes array')
      ENDIF

      IF (M1%D().ne.U%ndof) THEN
         write(*,*) 'M1 has ',M1%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("CPMM_CP8(): M1, U, # modes mismatch")
      ENDIF

      IF (M2%D().ne.U%ndof) THEN
         write(*,*) 'M2 has ',M2%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("CPMM_CP8(): M2, U, # modes mismatch")
      ENDIF

      IF (M1%R().ne.U%rk1) THEN
         write(*,*) 'M1 has ',M1%R(),' terms but U has ',U%rk1,' terms'
         call AbortWithError("CPMM_CP8(): M1, U, rank mismatch")
      ENDIF

      IF (M2%R().ne.U%rk2) THEN
         write(*,*) 'M2 has ',M2%R(),' terms but U has ',U%rk2,' terms'
         call AbortWithError("CPMM_CP8(): M2, U, rank mismatch")
      ENDIF

      if (M3%R().eq.0) call M3%new(U%rk3,U%r1,U%c2)

      IF (M3%D().ne.U%ndof) THEN
         write(*,*) 'M3 has ',M3%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("CPMM_CP8(): M3, U, # modes mismatch")
      ENDIF

      IF (M3%R().ne.U%rk3) THEN
         write(*,*) 'M3 has ',M3%R(),' terms but U has ',U%rk3,' terms'
         call AbortWithError("CPMM_CP8(): M3, U, rank mismatch")
      ENDIF

      DO d=1,U%ndof
         IF (M1%M(d).ne.U%ld1(d)) THEN
            write(*,*) 'M1: ',M1%M(d),' rows; U: ',U%ld1(d),' rows'
            call AbortWithError("CPMM_CP8(): M1, U, rows mismatch")
         ENDIF
         IF (M1%N(d).ne.U%sd1(d)) THEN
            write(*,*) 'M1: ',M1%N(d),' cols; U: ',U%sd1(d),' cols'
            call AbortWithError("CPMM_CP8(): M1, U, cols mismatch")
         ENDIF

         IF (M2%M(d).ne.U%ld2(d)) THEN
            write(*,*) 'M2: ',M2%M(d),' rows; U: ',U%ld2(d),' rows'
            call AbortWithError("CPMM_CP8(): M2, U, rows mismatch")
         ENDIF
         IF (M2%N(d).ne.U%sd2(d)) THEN
            write(*,*) 'M2: ',M2%N(d),' cols; U: ',U%sd2(d),' cols'
            call AbortWithError("CPMM_CP8(): M2, U, cols mismatch")
         ENDIF

         IF (M3%M(d).ne.U%r1(d)) THEN
            write(*,*) 'M3 has ',M3%M(d),' rows; U has ',U%r1(d),' rows'
            call AbortWithError("CPMM_CP8(): M3, U, rows mismatch")
         ENDIF
         IF (M3%N(d).ne.U%c2(d)) THEN
            write(*,*) 'M3 has ',M3%N(d),' cols; U has ',U%c2(d),' cols'
            call AbortWithError("CPMM_CP8(): M3, U, cols mismatch")
         ENDIF
      ENDDO

      call CPU_TIME(ti1)
      call nvtx_start('CPMM8')

      select case(U%algo)
         case(0)
            call CPMM_alg_cpu_CP8(U,M1,M2,M3,modes)
         case(1)
#if ACC_ENABLED
            if (.not.M3%live_on_device) call M3%createondevice()
            call CPMM_alg_tensor_CP8(U,M1,M2,M3,modes)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("CPMM_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("CPMM_CP8(): wrong algorithm")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      module_time(mpirank+1)=module_time(mpirank+1)+ti2-ti1

      end subroutine CPMM_CP8_general

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_alg_cpu_CP8(U,M1,M2,M3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MM8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      logical, intent(in) :: modes(:)
      logical :: docoef
      integer :: d,i,k,l
      real(kind=8), parameter :: b=0.0
      real(kind=8) :: afac

!     docoef ensures that coefs are multiplied for first included mode
!!! future fix: make sure that coefs are done even if no modes are selected !!!
      docoef=.TRUE.

      !$omp parallel
      do d=1,U%ndof
         if (modes(d)) then
            afac=1.d0

!           Unshifted terms
            if (d.eq.1) afac=U%alpha1*U%alpha2
            !$omp do private(i,k,l) collapse(2)
            do k=1,U%rk1
               do i=1,U%rk2
                  if (U%t3) then
                     l=(k-1)*U%rks2+i ! R-rank-major
                  else
                     l=(i-1)*U%rks1+k ! L-rank-major
                  endif
                  if (docoef) M3%coef(l)=M1%coef(k)*M2%coef(i)
                  call DGEMM(U%tN1,U%tN2,U%r1(d),U%c2(d),U%c1(d),afac,&
                             M1%base(M1%BS(k,d):M1%BF(k,d)),U%ld1(d),&
                             M2%base(M2%BS(i,d):M2%BF(i,d)),U%ld2(d),b,&
                             M3%base(M3%BS(l,d):M3%BF(l,d)),U%r1(d))
               enddo
            enddo
            !$omp enddo

!           Terms from shift on M1 only
            if (U%ish1.ne.0) then
               if (d.eq.1) afac=U%alpha2
               !$omp do private(i,k,l) collapse(2)
               do k=U%rk1+1,U%rks1
                  do i=1,U%rk2
                     if (U%t3) then
                        l=(k-1)*U%rks2+i ! R-rank-major
                     else
                        l=(i-1)*U%rks1+k ! L-rank-major
                     endif
                     if (docoef) M3%coef(l)=U%I1%coef(1)*M2%coef(i)
                     call DGEMM(U%tN1,U%tN2,U%r1(d),U%c2(d),U%c1(d),afac,&
                            U%I1%base(U%I1%BS(1,d):U%I1%BF(1,d)),U%ld1(d),&
                                M2%base(M2%BS(i,d):M2%BF(i,d)),U%ld2(d),b,&
                                M3%base(M3%BS(l,d):M3%BF(l,d)),U%r1(d))
                  enddo
               enddo
               !$omp enddo
            endif

!           Terms from shift on M2 only
            if (U%ish2.ne.0) then
               if (d.eq.1) afac=U%alpha1
               !$omp do private(i,k,l) collapse(2)
               do k=1,U%rk1
                  do i=U%rk2+1,U%rks2
                     if (U%t3) then
                        l=(k-1)*U%rks2+i ! R-rank-major
                     else
                        l=(i-1)*U%rks1+k ! L-rank-major
                     endif
                     if (docoef) M3%coef(l)=M1%coef(k)*U%I2%coef(1)
                     call DGEMM(U%tN1,U%tN2,U%r1(d),U%c2(d),U%c1(d),afac,&
                                M1%base(M1%BS(k,d):M1%BF(k,d)),U%ld1(d),&
                          U%I2%base(U%I2%BS(1,d):U%I2%BF(1,d)),U%ld2(d),b,&
                                M3%base(M3%BS(l,d):M3%BF(l,d)),U%r1(d))
                  enddo
               enddo
               !$omp enddo
            endif

!           Terms from shift on M1 and M2
            if (U%ish1.ne.0.and.U%ish2.ne.0) then
               if (d.eq.1) afac=1.0
               !$omp do private(i,k,l) collapse(2)
               do k=U%rk1+1,U%rks1
                  do i=U%rk2+1,U%rks2
                     if (U%t3) then
                        l=(k-1)*U%rks2+i ! R-rank-major
                     else
                        l=(i-1)*U%rks1+k ! L-rank-major
                     endif
                     if (docoef) M3%coef(l)=U%I1%coef(1)*U%I2%coef(1)
                     call DGEMM(U%tN1,U%tN2,U%r1(d),U%c2(d),U%c1(d),afac,&
                            U%I1%base(U%I1%BS(1,d):U%I1%BF(1,d)),U%ld1(d),&
                          U%I2%base(U%I2%BS(1,d):U%I2%BF(1,d)),U%ld2(d),b,&
                                M3%base(M3%BS(l,d):M3%BF(l,d)),U%r1(d))
                  enddo
               enddo
               !$omp enddo
            endif
            docoef=.FALSE.
         endif
      enddo
      !$omp end parallel

      end subroutine CPMM_alg_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_alg_tensor_CP8(U,M1,M2,M3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MM8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: M1,M2
      TYPE (CP8), INTENT(INOUT) :: M3
      logical, intent(in) :: modes(:)
      integer :: d,m,n,k,ms
      integer(4) :: stat
      real(kind=8), parameter :: alpha=1.0, beta=0.0
      real(kind=8) ,allocatable :: v1(:),v2(:),work(:)

#if ACC_ENABLED
!     Make sure arrays are on device
      if (.not.M1%live_on_device) call &
         AbortWithError('CPMM_alg_tensor_CP8(): M1 not on device')
      if (.not.M2%live_on_device) call &
         AbortWithError('CPMM_alg_tensor_CP8(): M2 not on device')
      if (.not.M3%live_on_device) call &
         AbortWithError('CPMM_alg_tensor_CP8(): M3 not on device')

!     Multiply bases
      do d=1,U%ndof
         if (modes(d)) then
            m=U%r1(d)*U%rks1
            n=U%c2(d)*U%rks2
            k=U%r2(d)
            ms=M3%MS(d)
            allocate(v1(m*k),v2(k*n),work(U%work_size(d)))
            !$acc data create(v1,v2,work)
            call Precopy_base_CP8(v1,M1,U%I1,d,U%ish1)
            call Precopy_base_CP8(v2,M2,U%I2,d,U%ish2)

            !$acc host_data use_device(v1,v2,M3%base,work)
            cutensor_status = cutensorcontract(cutensor_handle,&
                              U%mm_plan(d),alpha,v1,v2,beta,&
                              M3%base(ms),M3%base(ms),&
                              work,U%work_size(d),acc_stream_g)
            !$acc end host_data

            if (cutensor_status%stat.ne.0) then
               write(*,*) 'cutensorcontract, base mode ',d,&
                       '; exit = ',cutensor_status%stat
               call AbortWithError('Error in CPMM_alg_tensor_CP8()')
            endif
            !$acc wait
            !$acc end data
            deallocate(v1,v2,work)
         endif
      enddo

!     Multiply coefs
      allocate(v1(U%rks1),v2(U%rks2),work(U%work_sizeC))
      !$acc data create(v1,v2,work)
      call Precopy_coef_CP8(v1,M1,U%ish1,U%Esh1)
      call Precopy_coef_CP8(v2,M2,U%ish2,U%Esh2)

      !$acc host_data use_device(v1,v2,M3%coef,work)
      cutensor_status = cutensorcontract(cutensor_handle,&
                        U%mm_planC,alpha,v1,v2,beta,&
                        M3%coef(1),M3%coef(1),&
                        work,U%work_sizeC,acc_stream_g)
      !$acc end host_data

      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcontract, coef',&
                 '; exit = ',cutensor_status%stat
         call AbortWithError('Error in CPMM_alg_tensor_CP8()')
      endif

      !$acc wait
      !$acc end data
      deallocate(v1,v2,work)
#endif

      end subroutine CPMM_alg_tensor_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPMM8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
