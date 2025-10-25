!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPVV8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MYACC
      USE CPr8
#if ACC_ENABLED
      USE OPENACC
      USE CUDAFOR
      USE CUTENSOR_v2
#endif

      implicit none
      real(kind=8), allocatable, private :: pvv_time(:),norm_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE VV8
         integer, allocatable :: nbas(:)
         integer :: ndof,rk1,rk2,algo
         logical :: setup = .FALSE.
#if ACC_ENABLED
         type(cutensorDataType) :: dataType
         type(cutensorComputeDescriptor) :: descCompute
         type(cutensortensordescriptor), allocatable :: descV1(:),descV2(:),descV3(:)
         type(cutensortensordescriptor) :: descC1,descC2,descC3
         type(cutensoroperationdescriptor), allocatable :: op_desc(:)
         type(cutensoroperationdescriptor) :: op_descC
         type(cutensorplanpreference) :: plan_pref
         type(cutensorplan), allocatable :: vv_plan(:)
         type(cutensorplan) :: vv_planC
         real(kind=8), allocatable :: P1(:)
         integer(8), allocatable   :: work_size(:)
         integer(8) :: work_sizeC
#endif
         CONTAINS
            PROCEDURE :: new => New_VV8
            PROCEDURE :: flush => Flush_VV8
      END TYPE VV8

      INTERFACE CONSTPVV_CP8
        MODULE PROCEDURE CONST_PVV_CP8_wrapper
        MODULE PROCEDURE CONST_PVV_CP8
      END INTERFACE CONSTPVV_CP8

      INTERFACE CONSTPT_CP8
         MODULE PROCEDURE CONST_PT_CP8_wrapper
         MODULE PROCEDURE CONST_PT_CP8
      END INTERFACE CONSTPT_CP8

      INTERFACE CONSTPk_CP8
        MODULE PROCEDURE CONST_Pk_CP8_wrapper
        MODULE PROCEDURE CONST_Pk_CP8
      END INTERFACE CONSTPk_CP8

      INTERFACE UPDATEP_CP8
        MODULE PROCEDURE UPDATE_P_CP8_wrapper
        MODULE PROCEDURE UPDATE_P_CP8
      END INTERFACE UPDATEP_CP8

      INTERFACE UPDATEPCoef_CP8
        MODULE PROCEDURE UPDATE_P_Coef_CP8_wrapper
        MODULE PROCEDURE UPDATE_P_Coef_CP8
      END INTERFACE UPDATEPCoef_CP8

      INTERFACE BuildPVV_CP8
         MODULE PROCEDURE BuildPVV_CP8_general_wrapper
         MODULE PROCEDURE BuildPVV_CP8_general
      END INTERFACE BuildPVV_CP8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_CPVV8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(pvv_time(mpinodes),norm_time(mpinodes))
      pvv_time(:) = 0.d0
      norm_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_CPVV8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_CPVV8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: i,ierr

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()
      call Get_MPI_Timings('CPVV8: dot products ',pvv_time)
      call Get_MPI_Timings('CPVV8: normalization',norm_time)
      MODULE_SETUP = .FALSE.
      deallocate(pvv_time,norm_time)

      end subroutine Dispose_CPVV8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine New_VV8(U,V1,V2,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes VV8 type

      implicit none
      CLASS (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      integer, intent(in) :: algo
      integer :: d

!     Set parameters: # modes, ranks
      U%algo=algo
      U%ndof=V1%D()
      U%rk1=V1%R()
      U%rk2=V2%R()

!     Error checking
      IF (U%setup) &
         call AbortWithError('New_VV8(): U already set up')

      IF (V1%D().lt.1 .or. V2%D().lt.1) THEN
         write(*,*) 'V1, V2 have ',V1%D(),', ',V2%D(),&
                    ' modes, respectively; must both be > 0'
         call AbortWithError('New_VV8(): # modes must be > 0')
      ENDIF

      IF (V2%D().ne.V1%D()) THEN
         write(*,*) 'V1 has ',V1%D(),&
         ' modes but V2 has ',V1%D(),' modes'
         call AbortWithError('New_VV8(): V1, V2 mode mismatch')
      ENDIF

!     Vector sizes
      ALLOCATE(U%nbas(U%ndof))
      DO d=1,U%ndof
         U%nbas(d)=V1%MN(d)
         IF (V1%MN(d).ne.V2%MN(d)) THEN
            write(*,*) 'Mode ',d,': V1 length (',V1%MN(d),&
                       '), V2 length (',V2%MN(d),') must match'
            call AbortWithError('New_VV8(): bad vector lengths')
         ENDIF
      ENDDO

      if (algo.eq.1) then
         call Init_cutensor_v2_VV8(U)
      endif

      U%setup=.TRUE.

      end subroutine New_VV8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_cutensor_v2_VV8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes VV8 type

      implicit none
      CLASS (VV8) :: U
      integer(8) :: m1_ext(2),m2_ext(2),m3_ext(2)
      integer(8) :: c1_ext(1),c2_ext(1),c3_ext(2)
      integer(8) :: m1_str(2),m2_str(2),m3_str(2)
      integer(8) :: c1_str(1),c2_str(1),c3_str(2)
      integer    :: m1_ord(2),m2_ord(2),m3_ord(2)
      integer    :: c1_ord(1),c2_ord(1),c3_ord(2)
      integer    :: d,ndof

#if ACC_ENABLED
      ndof=U%ndof

      U%dataType = CUTENSOR_R_64F
      U%descCompute = CUTENSOR_COMPUTE_DESC_64F

      ! Create tensor descriptors for base of M1,M2,M3
      ALLOCATE(U%descV1(ndof),U%descV2(ndof),U%descV3(ndof))
      DO d=1,ndof
         m1_ext=(/U%nbas(d),U%rk1/)
         m1_str=(/1,U%nbas(d)/)
         m2_ext=(/U%nbas(d),U%rk2/)
         m2_str=(/1,U%nbas(d)/)
         m3_ext=(/U%rk1,U%rk2/)
         m3_str=(/1,U%rk1/)

         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
              U%descV1(d),2,m1_ext,m1_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, M1 base mode',d,&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
              U%descV2(d),2,m2_ext,m2_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, M2 base mode',d,&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
              U%descV3(d),2,m3_ext,m3_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, M3 base mode',d,&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
      ENDDO

!     Create tensor descriptors for coef of M1,M2,M3
      c1_ext=(/U%rk1/)
      c1_str=(/1/)
      c2_ext=(/U%rk2/)
      c2_str=(/1/)
      c3_ext=(/U%rk1,U%rk2/)
      c3_str=(/1,U%rk1/)
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC1,1,c1_ext,c1_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M1 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC2,1,c2_ext,c2_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M2 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif
      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descC3,2,c3_ext,c3_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, M3 coef',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif

      ! Create operator (contraction) descriptor for base of M1,M2,M3
      ALLOCATE(U%op_desc(ndof))
      m1_ord=(/1,2/)
      m2_ord=(/1,3/)
      m3_ord=(/2,3/)

      DO d=1,ndof
         cutensor_status=cutensorcreatecontraction(cutensor_handle,U%op_desc(d),&
                         U%descV1(d),m1_ord,CUTENSOR_OP_IDENTITY,&
                         U%descV2(d),m2_ord,CUTENSOR_OP_IDENTITY,&
                         U%descV3(d),m3_ord,CUTENSOR_OP_IDENTITY,&
                         U%descV3(d),m3_ord,U%descCompute)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreatecontraction, base mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
      ENDDO

      ! Create operator (contraction) descriptor for coef of M1,M2,M3
      c1_ord=(/1/)
      c2_ord=(/2/)
      c3_ord=(/1,2/)
      cutensor_status=cutensorcreatecontraction(cutensor_handle,U%op_descC,&
                      U%descC1,c1_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC2,c2_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC3,c3_ord,CUTENSOR_OP_IDENTITY,&
                      U%descC3,c3_ord,U%descCompute)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreatecontraction, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif

      ! Create plan preference
      cutensor_status=cutensorcreateplanpreference(cutensor_handle,U%plan_pref,&
                      CUTENSOR_ALGO_DEFAULT,CUTENSOR_JIT_MODE_NONE)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreateplanpreference; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif

      ! Estimate workspace size for base contraction
      ALLOCATE(U%work_size(ndof))
      do d=1,ndof
         cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%op_desc(d),&
                         U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%work_size(d))
!         write(*,*) 'mode ',d,'; work size = ',U%work_size(d)
         if (cutensor_status%stat .ne.0) then
            write(*,*) 'cutensorestimateworkspacesize, base mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
      enddo

      ! Estimate workspace size for coef contraction
      cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%op_descC,&
                      U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%work_sizeC)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorestimateworkspacesize, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif

      ! Create the plan for base contraction
      ALLOCATE(U%vv_plan(ndof))
      do d=1,ndof
        cutensor_status=cutensorcreateplan(cutensor_handle,U%vv_plan(d),&
                        U%op_desc(d),U%plan_pref,U%work_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreateplan, base mode',d,&
                       '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
         endif
      enddo

      ! Create the plan for coef contraction
       cutensor_status=cutensorcreateplan(cutensor_handle,U%vv_planC,&
                       U%op_descC,U%plan_pref,U%work_sizeC)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreateplan, coef',&
                    '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_v2_VV8()')
      endif

      ALLOCATE(U%P1(U%rk1*U%rk2))
      !$acc enter data create(U%P1)
#endif

      end subroutine Init_cutensor_v2_VV8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_VV8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates VV8 type

      implicit none
      CLASS (VV8) :: U
      integer :: d,ndof

      ndof=U%ndof

      IF (.not.U%setup) &
         call AbortWithError('Flush_VV8(): U not setup')

      IF (ALLOCATED(U%nbas))  DEALLOCATE(U%nbas)

#if ACC_ENABLED
      if (U%algo.eq.1) then

         ! Destroy plan for coef and base
         cutensor_status = cutensorDestroyPlan(U%vv_planC)
         DO d=1,ndof
            cutensor_status = cutensorDestroyPlan(U%vv_plan(d))
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
            cutensor_status = cutensorDestroyTensorDescriptor(U%descV1(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descV2(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descV3(d))
         ENDDO

         IF (ALLOCATED(U%vv_plan)) DEALLOCATE(U%vv_plan)
         IF (ALLOCATED(U%work_size)) DEALLOCATE(U%work_size)
         IF (ALLOCATED(U%op_desc)) DEALLOCATE(U%op_desc)
         IF (ALLOCATED(U%descV1)) DEALLOCATE(U%descV1)
         IF (ALLOCATED(U%descV2)) DEALLOCATE(U%descV2)
         IF (ALLOCATED(U%descV3)) DEALLOCATE(U%descV3)
         !$acc exit data delete(U%P1)
         IF (ALLOCATED(U%P1))  DEALLOCATE(U%P1)
      endif
#endif

      U%setup=.FALSE.

      end subroutine Flush_VV8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Normalize_CP8(v,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes CP-vector

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: algo
      real(kind=8) :: norm

!     Normalize base of P
      call NormBase_CP8(v,algo)

!     Multiply coefs by 1/sqrt(norm)
      norm=1.d0/sqrt(PRODVV_CP8(v,v,algo))
      call ScaleCoefs_CP8(v,norm,algo)

      end subroutine Normalize_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODVV_CP8(v,w,algo) result(pvv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to do
! the operation, for all modes with coefs.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), intent(in) :: v,w
      integer, intent(in) :: algo
      real(kind=8), allocatable :: P(:)
      real(kind=8) :: pvv

!     Build P for all modes and coefs
      allocate(P(v%R()*w%R()))
      if (algo.eq.1) then
#if ACC_ENABLED
         !$acc enter data create(P)
#else
         write(*,*) 'OpenACC algorithm not enabled'
         call AbortWithError("PRODVV_CP8(): wrong algorithm")
#endif
      endif

      call U%new(v,w,algo)
      call CONSTPVV_CP8(U,v,w,P)
      call U%flush()

!     Reduce entries of P -> pvv
      pvv=ReduceP_CP8(P,algo)

      if (algo.eq.1) then
#if ACC_ENABLED
         !$acc exit data delete(P)
#endif
      endif
      deallocate(P)

      end function PRODVV_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_PVV_CP8_wrapper(V1,V2,P,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to do
! the operation, for all modes with coefs.
! P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: algo

      call U%new(V1,V2,algo) 
      call CONSTPVV_CP8(U,V1,V2,P)
      call U%flush()

      end subroutine CONST_PVV_CP8_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_PVV_CP8(U,V1,V2,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to do
! the operation, for all modes with coefs. 
! P must be allocated before calling.

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      logical, allocatable :: modes(:)
      logical :: docoef,divide,ovr

      docoef=.TRUE.
      ovr=.TRUE.
      divide=.FALSE.
      allocate(modes(U%ndof))
      modes(:)=.TRUE.
      call BuildPVV_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
      deallocate(modes)

      end subroutine CONST_PVV_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_PT_CP8_wrapper(V1,V2,k,P,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to do
! the operation, for all modes (skipping k).
! P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: k,algo

      call U%new(V1,V2,algo)       
      call CONSTPT_CP8(U,V1,V2,k,P)
      call U%flush()

      end subroutine CONST_PT_CP8_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_PT_CP8(U,V1,V2,k,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to do
! the operation, for all modes (skipping k). P must be allocated before
! calling.

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in)  :: k
      logical, allocatable :: modes(:)
      logical :: docoef,divide,ovr

      IF (k.lt.0 .or. k.gt.U%ndof) THEN
          write(*,*) 'mode: ',k,' out of range [0,',U%ndof,']'
          call AbortWithError('Error in CONST_PT_CP8()')
      ENDIF

      docoef=.FALSE.
      ovr=.TRUE.
      divide=.FALSE.
      allocate(modes(U%ndof))
      modes(:)=.TRUE.
      if (k.gt.0) modes(k)=.FALSE.
      call BuildPVV_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
      deallocate(modes)

      end subroutine CONST_PT_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_Pk_CP8_wrapper(V1,V2,k,P,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to do
! the operation, for mode k.
! P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: k,algo

      call U%new(V1,V2,algo)
      call CONSTPk_CP8(U,V1,V2,k,P)
      call U%flush()

      end subroutine CONST_Pk_CP8_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONST_Pk_CP8(U,V1,V2,k,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to do
! the operation, for mode k. P must be allocated before calling.
! Set k=0 to initialize P(:)=1

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in)  :: k
      logical, allocatable :: modes(:)
      logical :: docoef,divide,ovr

      IF (k.lt.0 .or. k.gt.U%ndof) THEN
          write(*,*) 'mode: ',k,' out of range [0,',U%ndof,']'
          call AbortWithError('Error in CONST_Pk_CP8()')
      ENDIF

      docoef=.FALSE.
      ovr=.TRUE.
      divide=.FALSE.
      allocate(modes(U%ndof))
      modes(:)=.FALSE.
      if (k.gt.0) modes(k)=.TRUE.
      call BuildPVV_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
      deallocate(modes)

      end subroutine CONST_Pk_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UPDATE_P_CP8_wrapper(V1,V2,k,P,divide,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to
! update/downdate P, for mode k. P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: k,algo
      logical, intent(in) :: divide

      call U%new(V1,V2,algo)
      call UPDATEP_CP8(U,V1,V2,k,P,divide)
      call U%flush()

      end subroutine UPDATE_P_CP8_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UPDATE_P_CP8(U,V1,V2,k,P,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to
! update/downdate P, for mode k. P must be allocated before calling.

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in)  :: k
      logical, intent(in)  :: divide
      logical, allocatable :: modes(:)
      logical :: docoef,ovr

      IF (k.lt.1 .or. k.gt.U%ndof) THEN
          write(*,*) 'mode: ',k,' out of range [1,',U%ndof,']'
          call AbortWithError('Error in UPDATE_P_CP8()')
      ENDIF

      docoef=.FALSE.
      ovr=.FALSE.
      allocate(modes(U%ndof))
      modes(:)=.FALSE.
      if (k.gt.0) modes(k)=.TRUE.
      call BuildPVV_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
      deallocate(modes)

      end subroutine UPDATE_P_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UPDATE_P_coef_CP8_wrapper(V1,V2,P,divide,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to
! update/downdate P, using coefs. P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: algo
      logical, intent(in) :: divide

      call U%new(V1,V2,algo)
      call UPDATEPCoef_CP8(U,V1,V2,P,divide)
      call U%flush()

      end subroutine UPDATE_P_coef_CP8_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UPDATE_P_coef_CP8(U,V1,V2,P,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to
! update/downdate P, using coefs. P must be allocated before calling.

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      logical, intent(in)  :: divide
      logical, allocatable :: modes(:)
      logical :: docoef,ovr

      docoef=.TRUE.
      ovr=.FALSE.
      allocate(modes(U%ndof))
      modes(:)=.FALSE.
      call BuildPVV_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
      deallocate(modes)

      end subroutine UPDATE_P_Coef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildPVV_CP8_general_wrapper(V1,V2,P,modes,docoef,ovr,&
                                              divide,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version creates a VV8 type to
! update/downdate P, using coefs. P must be allocated before calling.

      implicit none
      TYPE (VV8) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      integer, intent(in) :: algo
      logical, intent(in) :: modes(:)
      logical, intent(in) :: docoef,divide,ovr

      call U%new(V1,V2,algo)
      call BuildPVV_CP8_general(U,V1,V2,P,modes,docoef,ovr,divide)
      call U%flush()

      end subroutine BuildPVV_CP8_general_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildPVV_CP8_general(U,V1,V2,P,modes,docoef,ovr,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine; this version uses existing VV8 type to do
! the operation, general version. P must be allocated before calling.

      implicit none
      TYPE (VV8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      logical, intent(in) :: modes(:)
      logical, intent(in) :: docoef,divide,ovr
      real(kind=8) :: ti1,ti2
      integer :: d

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('BuildPVV_CP8(): U not setup')

      IF (SIZE(modes).ne.U%ndof) THEN
         write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',U%ndof
         call AbortWithError('BuildPVV_CP8(): bad modes array size')
      ENDIF

      IF (V1%D().ne.U%ndof) THEN
         write(*,*) 'V1 has ',V1%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("BuildPVV_CP8(): V1, U, # modes mismatch")
      ENDIF

      IF (V2%D().ne.U%ndof) THEN
         write(*,*) 'V2 has ',V2%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("BuildPVV_CP8(): V2, U, # modes mismatch")
      ENDIF

      IF (V1%R().ne.U%rk1) THEN
         write(*,*) 'V1 has ',V1%R(),' terms but U has ',U%rk1,' terms'
         call AbortWithError("BuildPVV_CP8(): V1, U, rank mismatch")
      ENDIF

      IF (V2%R().ne.U%rk2) THEN
         write(*,*) 'V2 has ',V2%R(),' terms but U has ',U%rk2,' terms'
         call AbortWithError("BuildPVV_CP8(): V2, U, rank mismatch")
      ENDIF

      DO d=1,U%ndof
         IF (V1%MN(d).ne.U%nbas(d)) THEN
            write(*,*) 'V1: ',V1%MN(d),' len; U: ',U%nbas(d),' len'
            call AbortWithError("BuildPVV_CP8(): V1, U, len mismatch")
         ENDIF

         IF (V2%MN(d).ne.U%nbas(d)) THEN
            write(*,*) 'V2: ',V2%MN(d),' len; U: ',U%nbas(d),' len'
            call AbortWithError("BuildPVV_CP8(): V2, U, len mismatch")
         ENDIF
      ENDDO

      IF (SIZE(P).ne.U%rk1*U%rk2) THEN
         write(*,*) 'P size: ',SIZE(P),' must be ',U%rk1*U%rk2
         call AbortWithError("BuildPVV_CP8(): wrong size P")
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('CPVV8: BuildPVV')

      select case(U%algo)
         case(0)
            call PVV_alg_cpu_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
         case(1)
#if ACC_ENABLED
            call PVV_alg_tensor_CP8(U,V1,V2,P,modes,docoef,ovr,divide)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("Build_PVV_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("Build_PVV_CP8(): wrong algorithm")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      pvv_time(mpirank+1)=pvv_time(mpirank+1)+ti2-ti1

      end subroutine BuildPVV_CP8_general

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PVV_alg_cpu_CP8(U,V1,V2,P,modes,docoef,ovr,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine, CPU version

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), INTENT(INOUT) :: P(:)
      real(kind=8), parameter :: one=1.d0
      real(kind=8) :: dp
      logical, intent(in) :: modes(:)
      logical, intent(in) :: docoef,divide,ovr
      integer :: d,r1,r2,bs1,bf1,bs2,bf2,ind

      !$omp parallel
!     Set ovr=.TRUE. to compute inner product matrix from scratch
      if (ovr) then
         !$omp do private(r1)
         do r1=1,U%rk1*U%rk2
            P(r1)=one
         enddo
         !$omp enddo
      endif

!     Inner product of coefficients
      if (docoef) then
         !$omp do private(r1,r2,ind,dp) collapse(2)
         do r2=1,U%rk2
            do r1=1,U%rk1
               ind=(r2-1)*U%rk1+r1
               dp=V1%coef(r1)*V2%coef(r2)
               if (divide) dp=one/dp
               P(ind)=P(ind)*dp
            enddo
         enddo
         !$omp enddo
      endif

!     Inner product of bases
      do d=1,U%ndof
         if (modes(d)) then
            !$omp do private(r1,r2,ind,bs1,bf1,bs2,bf2,dp) collapse(2)
            do r2=1,U%rk2
               do r1=1,U%rk1
                  bs1=V1%BS(r1,d)
                  bf1=V1%BF(r1,d)
                  bs2=V2%BS(r2,d)
                  bf2=V2%BF(r2,d)
                  ind=(r2-1)*U%rk1+r1
                  dp=dot_product(V1%base(bs1:bf1),V2%base(bs2:bf2))
                  if (divide) dp=one/dp 
                  P(ind)=P(ind)*dp
               enddo
            enddo
            !$omp enddo
         endif
      enddo
      !$omp end parallel

      end subroutine PVV_alg_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PVV_alg_tensor_CP8(U,V1,V2,P,modes,docoef,ovr,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP inner product routine, GPU (CUTENSOR) version

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: V1,V2
      real(kind=8), intent(inout) :: P(:)
      logical, intent(in) :: modes(:)
      logical, intent(in) :: docoef,divide,ovr
      integer :: d,r1,r2,ms1,mf1,ms2,mf2

#if ACC_ENABLED
!     Make sure arrays are on device
      if (.not.V1%live_on_device) call &
         AbortWithError('PVV_alg_tensor_CP8(): V1 not on device')
      if (.not.V2%live_on_device) call &
         AbortWithError('PVV_alg_tensor_CP8(): V2 not on device')
      if (.not.acc_is_present(P)) call &
         AbortWithError('PVV_alg_tensor_CP8(): P not on device')
      if (.not.acc_is_present(U%P1)) call &
         AbortWithError('PVV_alg_tensor_CP8(): U%P1 not on device')

!     Set ovr=.TRUE. to compute inner product matrix from scratch
      if (ovr) call PVV_setP_CP8(P)

!     cutensor mult: coefs
      if (docoef) then
         call PVV_contraction_wrapper_CP8(U,V1%coef,V2%coef,U%P1,0)
         call PVV_updateP_CP8(P,U%P1,divide)
      endif

!     cutensor mult: base
      do d=1,U%ndof
         if (modes(d)) then
            ms1=V1%MS(d)
            mf1=V1%MF(d)
            ms2=V2%MS(d)
            mf2=V2%MF(d)
            call PVV_contraction_wrapper_CP8(U,V1%base(ms1:mf1),V2%base(ms2:mf2),U%P1,d)
            call PVV_updateP_CP8(P,U%P1,divide)
         endif
      enddo
#endif
      end subroutine PVV_alg_tensor_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PVV_contraction_wrapper_CP8(U,v1,v2,P,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for calling cuTensor to compute inner products

      implicit none
      TYPE (VV8), INTENT(IN) :: U
      integer, intent(in) :: d
      real(kind=8), intent(in) :: v1(:),v2(:)
      real(kind=8), intent(inout) :: P(:)
      real(kind=8), allocatable :: work(:)
      real(kind=8), parameter   :: alpha=1.0, beta=0.0

#if ACC_ENABLED
      if (d.gt.0) then
         allocate(work(U%work_size(d))) ! base mode d
      else
         allocate(work(U%work_sizeC))   ! coefs
      endif
      !$acc data create(work)
      !$acc host_data use_device(v1,v2,P,work)
      if (d.gt.0) then ! base mode d
         cutensor_status = cutensorcontract(cutensor_handle,&
                           U%vv_plan(d),alpha,v1,v2,beta,P,P,&
                           work,U%work_size(d),acc_stream_g)
         if (cutensor_status%stat.ne.0) write(*,*) &
            'cutensorcontract, base mode ',d,'; exit = ',cutensor_status%stat
      else             ! coefs
         cutensor_status = cutensorcontract(cutensor_handle,&
                           U%vv_planC,alpha,v1,v2,beta,P,P,&
                           work,U%work_sizeC,acc_stream_g)
         if (cutensor_status%stat.ne.0) write(*,*) &
            'cutensorcontract, coef; exit = ',cutensor_status%stat
      endif
      if (cutensor_status%stat.ne.0) &
         call AbortWithError('Error in PVV_contraction_wrapper_CP8()')
      !$acc end host_data
      !$acc wait
      !$acc end data
      deallocate(work)
#endif
      end subroutine PVV_contraction_wrapper_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PVV_setP_CP8(P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! GPU kernel for initializing inner products

      implicit none
      real(kind=8), intent(inout) :: P(:)
      integer :: i,l
      real(kind=8), parameter :: one=1.d0

#if ACC_ENABLED
      l=SIZE(P)
      !$acc data present(P)
      !$acc parallel loop gang vector async
      do i=1,l
         P(i)=one
      enddo
      !$acc end data
      !$acc wait
#endif 

      end subroutine PVV_setP_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PVV_updateP_CP8(P,P1,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! GPU kernel to update or downdate inner products

      implicit none
      real(kind=8), intent(inout) :: P(:)
      real(kind=8), intent(in) :: P1(:)
      logical, intent(in) :: divide
      integer :: i,l

#if ACC_ENABLED
      l=SIZE(P)
      !$acc data present(P,P1)
      if (divide) then
         !$acc parallel loop gang vector async
         do i=1,l
            P(i)=P(i)/P1(i)
         enddo
      else
         !$acc parallel loop gang vector async
         do i=1,l
            P(i)=P(i)*P1(i)
         enddo
      endif
      !$acc end data
      !$acc wait
#endif 

      end subroutine PVV_updateP_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NormBase_CP8(v,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes base of all modes

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: algo
      integer :: d
      logical :: ovr

      do d=1,v%D()
         call NormBaseD_CP8(v,d,.FALSE.,algo)
      enddo

      end subroutine NormBase_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NormBaseD_CP8(v,d,ovr,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to normalize base of single mode. Set ovr=.TRUE. to overwrite 
! the coef, and ovr=.FALSE. to multiply the coef by the norm

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: d,algo
      logical, intent(in) :: ovr
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()

      IF (d.lt.1 .or. d.gt.v%D()) THEN
          write(*,*) 'mode: ',d,' out of range [1,',v%D(),']'
          call AbortWithError('Error in NormBaseD_CP8()')
      ENDIF

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            call NormBaseD_cpu_CP8(v,d,ovr)
         case(1)
#if ACC_ENABLED
            call NormBaseD_acc_CP8(v,d,ovr)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("NormBaseD(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("NormBaseD(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      norm_time(mpirank+1)=norm_time(mpirank+1)+ti2-ti1

      end subroutine NormBaseD_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NormBaseD_cpu_CP8(v,d,ovr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to normalize base of single mode, CPU version.

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: d
      logical, intent(in) :: ovr
      integer :: r,rk,bs,bf
      real(kind=8) :: prod

      rk=v%R()

      do r=1,rk
         bs=v%BS(r,d)
         bf=v%BF(r,d)
         prod=sqrt(abs(dot_product(v%base(bs:bf),v%base(bs:bf))))
         v%base(bs:bf)=v%base(bs:bf)/prod
         if (.not.ovr) prod=v%coef(r)*prod
         v%coef(r)=prod
      enddo

      end subroutine NormBaseD_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NormBaseD_acc_CP8(v,d,ovr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to normalize base of single mode, OpenACC version.

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: d
      logical, intent(in) :: ovr
      integer :: r,rk,i,n,st
      real(kind=8) :: prod

#if ACC_ENABLED
!     Make sure array is on device
      if (.not.v%live_on_device) call &
         AbortWithError('NormBaseD_acc_CP8(): v not on device')

      rk=v%R()

      !$acc data present(v%base,v%coef,v%look,v%nbas)
      !$acc parallel loop gang private(st,n,prod) async
      do r=1,rk
         st=v%look(r,d)
         n=v%nbas(d)
         prod=0.d0
         !$acc loop vector reduction(+:prod)
         do i=1,n
            prod=prod+v%base(st+i)**2
         enddo

         prod=1.d0/sqrt(prod)

         !$acc loop vector
         do i=1,n
            v%base(st+i)=v%base(st+i)*prod
         enddo

         prod=1.d0/prod
         if (.not.ovr) prod=prod*v%coef(r)
         v%coef(r)=prod
      enddo
      !$acc end data
      !$acc wait
#endif

      end subroutine NormBaseD_acc_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ReduceP_CP8(P,algo) result(prod)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to sum elements of P, for normalization

      implicit none
      real(kind=8), intent(in) :: P(:)
      integer, intent(in) :: algo
      real(kind=8) :: prod
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            prod=ReduceP_cpu_CP8(P)
         case(1)
#if ACC_ENABLED
            prod=ReduceP_acc_CP8(P)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("ReduceP_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("ReduceP_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      norm_time(mpirank+1)=norm_time(mpirank+1)+ti2-ti1

      end function ReduceP_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ReduceP_cpu_CP8(P) result(prod)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to sum elements of P, for normalization, cpu version

      implicit none
      real(kind=8), intent(in) :: P(:)
      real(kind=8) :: prod
      integer :: i,n

      n=SIZE(P)
      prod=0.d0
      do i=1,n
         prod=prod+P(i)
      enddo

      end function ReduceP_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ReduceP_acc_CP8(P) result(prod)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to sum elements of P, for normalization, OpenACC version

      implicit none
      real(kind=8), intent(in) :: P(:)
      real(kind=8) :: ptmp,prod
      integer :: i,iout,n
      
#if ACC_ENABLED
      if (.not.acc_is_present(P)) call &
         AbortWithError('ReduceP_acc_CP8(): P not on device')

      n=SIZE(P)
      prod=0.d0
      !$acc data present(P)
      !$acc parallel loop gang private(ptmp) vector_length(256) async
      do iout=1,n,256
         ptmp=0.d0
         !$acc loop vector reduction(+:ptmp) shortloop
         do i=iout,min(iout+255,n)
            ptmp=ptmp+P(i)
         enddo
         !$acc atomic
         prod=prod+ptmp
      enddo
      !$acc end data
      !$acc wait
#endif

      end function ReduceP_acc_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCNorm_CP8(v,algo) result(prod)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to compute coefficient norm of v

      implicit none
      TYPE (CP8), intent(in) :: v
      integer, intent(in) :: algo
      real(kind=8) :: prod
      integer :: i,rk
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            prod=dot_product(v%coef(:),v%coef(:))
         case(1)
#if ACC_ENABLED
            rk=v%R()
            prod=0.d0
            !$acc data present(v%coef)
            !$acc parallel loop vector reduction(+:prod) async
            do i=1,rk
               prod=prod+v%coef(i)*v%coef(i)
            enddo
            !$acc end data
            !$acc wait
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("ReduceP_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("GetCNorm_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      norm_time(mpirank+1)=norm_time(mpirank+1)+ti2-ti1

      end function GetCNorm_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ScaleCoefs_CP8(v,fac,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to scale coefficients of v, for normalization

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: algo
      real(kind=8), intent(in) :: fac
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPVV8_Module()

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            call v%mult(fac)
         case(1)
#if ACC_ENABLED
            call ScaleCoefs_acc_CP8(v,fac)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("ScaleCoefs_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("ScaleCoefs_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      norm_time(mpirank+1)=norm_time(mpirank+1)+ti2-ti1

      end subroutine ScaleCoefs_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ScaleCoefs_acc_CP8(v,fac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Kernel to scale coefficients of v, for normalization

      implicit none
      TYPE (CP8) :: v
      real(kind=8), intent(in) :: fac
      integer :: i,rk

#if ACC_ENABLED
!     Make sure array is on device
      if (.not.v%live_on_device) call &
         AbortWithError('ScaleCoefs_acc_CP8(): v not on device')

      rk=v%R()

      !$acc data present(v%coef)
      !$acc parallel loop gang vector async
      do i=1,rk
         v%coef(i)=v%coef(i)*fac
      enddo
      !$acc end data
      !$acc wait
#endif

      end subroutine ScaleCoefs_acc_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPVV8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
