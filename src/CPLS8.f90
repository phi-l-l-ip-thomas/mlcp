!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPLS8

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
      real(kind=8), allocatable, private :: lhs_time(:),putsoln_time(:)
      real(kind=8), allocatable, private :: rhs_time(:),prepg_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE LS8
         integer, allocatable :: nbas(:),m(:),n(:)
         integer :: algo,ndof,rksAp,rksA,rkF,rkG,rkW,ish
         real(kind=8) :: Esh
         logical :: useA,useW
         logical :: setup = .FALSE.
         logical :: normal
#if ACC_ENABLED
         type(cutensorDataType) :: dataType
         type(cutensorComputeDescriptor) :: descCompute
         type(cutensortensordescriptor), allocatable :: descWb(:),descLHS(:)
         type(cutensortensordescriptor), allocatable :: descGb(:),descRHS(:)
         type(cutensortensordescriptor) :: descB,descP
         type(cutensoroperationdescriptor), allocatable :: Lop_desc(:),Rop_desc(:)
         type(cutensorplanpreference) :: plan_pref
         type(cutensorplan), allocatable :: Lplan(:),Rplan(:)
         integer(8), allocatable :: Lwork_size(:),Rwork_size(:)
#endif
         CONTAINS
            PROCEDURE :: new => New_LS8
            PROCEDURE :: flush => Flush_LS8
      END TYPE LS8

      INTERFACE GetLHS_CP8
         MODULE PROCEDURE GetLHS_LinSys_CP8
         MODULE PROCEDURE GetLHS_ALS_CP8
      END INTERFACE
      
      INTERFACE PrepG_CP8
         MODULE PROCEDURE PrepG_allmodes_CP8
         MODULE PROCEDURE PrepG_1mode_CP8
         MODULE PROCEDURE PrepG_general_CP8
      END INTERFACE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_CPLS8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(lhs_time(mpinodes),rhs_time(mpinodes))
      allocate(putsoln_time(mpinodes),prepg_time(mpinodes))
      lhs_time(:) = 0.d0
      rhs_time(:) = 0.d0
      putsoln_time(:)=0.d0
      prepg_time(:)=0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_CPLS8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_CPLS8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: i,ierr

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()
      call Get_MPI_Timings('CPLS8:  left-hand sides',lhs_time)
      call Get_MPI_Timings('CPLS8: right-hand sides',rhs_time)
      call Get_MPI_Timings('CPLS8: copy LS solution',putsoln_time)
      call Get_MPI_Timings('CPLS8: prepare G vector',prepg_time)

      MODULE_SETUP = .FALSE.
      deallocate(lhs_time,rhs_time,putsoln_time,prepg_time)

      end subroutine Dispose_CPLS8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine New_LS8(U,A,F,G,W,ish,Esh,algo,usenormal)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes LS8 type for ALS or linear solver
! Call with dummy arguments A,W for ALS

      implicit none
      CLASS (LS8) :: U
      TYPE (CP8), INTENT(IN)   :: A,F,G,W
      logical, intent(in), optional :: usenormal
      integer, intent(in)      :: ish,algo
      real(kind=8), intent(in) :: Esh
      integer :: d
      logical, parameter :: normalform=.TRUE.

!     Set parameters: # modes, ranks
      IF (present(usenormal)) THEN
         U%normal=usenormal
      ELSE
         U%normal=.TRUE.
      ENDIF
      U%algo=algo
      U%ndof=G%D()
      U%useA=.TRUE.
      U%useW=.TRUE.

      if (A%R().eq.0) then
         U%useA=.FALSE.
      endif

      if (W%R().eq.0) then
         U%useW=.FALSE.
      endif

      U%rksA=A%R()
      if (ish.ne.0) U%rksA=U%rksA+1
      U%rksA=max(U%rksA,1)
      U%rksAp=1
      if (U%normal) U%rksAp=U%rksA

      U%ish=ish
      U%Esh=Esh
      U%rkF=F%R()
      U%rkG=G%R()
      U%rkW=max(W%R(),1)

!     Error checking
      IF (U%setup) &
         call AbortWithError('New_LS8(): U already set up')

      IF (G%D().lt.1) THEN
         write(*,*) 'G has ',G%D(),' modes; must be > 0'
         call AbortWithError('New_LS8(): # modes must be > 0')
      ENDIF

      IF (F%D().ne.G%D()) THEN
         write(*,*) 'F has ',F%D(),&
         ' modes but G has ',G%D(),' modes'
         call AbortWithError('New_LS8(): F,G mode mismatch')
      ENDIF

      IF (A%D().ne.G%D() .and. A%D().ne.0) THEN
         write(*,*) 'A has ',A%D(),&
         ' modes but G has ',G%D(),' modes'
         call AbortWithError('New_LS8(): A,G mode mismatch')
      ENDIF

      IF (W%D().ne.G%D() .and. W%D().ne.0) THEN
         write(*,*) 'W has ',W%D(),&
         ' modes but G has ',G%D(),' modes'
         call AbortWithError('New_LS8(): W,G mode mismatch')
      ENDIF

!     Vector sizes
      ALLOCATE(U%nbas(U%ndof),U%m(U%ndof),U%n(U%ndof))
      DO d=1,U%ndof
         U%nbas(d)=G%MN(d)
         U%m(d)=F%m(d)
         if (U%useA) then
            U%n(d)=F%m(d)
         else
            U%n(d)=1
         endif
      ENDDO

      if (algo.eq.1) then
         call Init_cutensor_RHS_LS8(U) ! Must be done before LHS
         if (U%useA.or.U%useW) then
            ! Ordinary ALS just copies B; wALS and LS uses cutensor
            call Init_cutensor_LHS_LS8(U)
         endif
      endif

      U%setup=.TRUE.

      end subroutine New_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_cutensor_LHS_LS8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes LS8 type

      implicit none
      CLASS (LS8) :: U
      integer(8) :: b_ext(4),wb_ext(4),lhs_ext(4)
      integer(8) :: b_str(4),wb_str(4),lhs_str(4)
      integer    :: b_ord(4),wb_ord(4),lhs_ord(4)
      integer    :: d,ndof

#if ACC_ENABLED
      ndof=U%ndof
      U%dataType = CUTENSOR_R_64F
      U%descCompute = CUTENSOR_COMPUTE_DESC_64F

      ! Create tensor descriptors for A^TA or A^TWA base and l.h.s.
      ALLOCATE(U%descWb(ndof),U%descLHS(ndof))
      DO d=1,ndof
         wb_ext=(/U%rksA*U%rkW,U%rksAp,U%n(d),U%m(d)/)
         if (U%useA) then ! LS
            wb_str=(/1,U%rksA*U%rkW,U%rksA*U%rkW*U%rksAp,U%rksA*U%rkW*U%rksAp*U%n(d)/)
         else             ! wALS
            wb_str=(/1,U%rksA*U%rkW,U%rksA*U%rkW*U%rksAp,U%rksA*U%rkW*U%rksAp*(U%m(d)+1)/)
         endif

         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
              U%descWb(d),SIZE(wb_ord),wb_ext,wb_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, W base mode',d,&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_LHS_LS8()')
         endif

         lhs_ext=(/U%rkF,U%n(d),U%rkF,U%m(d)/)
         lhs_str=(/1,U%rkF,U%n(d)*U%rkF,U%n(d)*U%rkF**2/)
         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
           U%descLHS(d),SIZE(lhs_ord),lhs_ext,lhs_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, LHS',&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_LHS_LS8()')
         endif
      ENDDO

      ! Create tensor descriptor for P matrix <WAF,AF>
      b_ext=(/U%rkW*U%rksA,U%rkF,U%rksAp,U%rkF/)
      b_str=(/1,U%rkW*U%rksA,U%rkW*U%rksA*U%rkF,U%rkW*U%rksA*U%rkF*U%rksAp/)

      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descB,SIZE(b_ord),b_ext,b_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, B matrix',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_LHS_LS8()')
      endif

      ! Create operator descriptors for B*A^TWA mult
      ALLOCATE(U%Lop_desc(ndof))
      wb_ord =(/5,6,1,2/) ! rWrA,rA',M,M'
      b_ord  =(/5,4,6,3/) ! rWrA,rF,rA',rF'
      lhs_ord=(/3,1,4,2/) ! rF',M,rF,M'
      
      DO d=1,ndof
         cutensor_status=cutensorcreatecontraction(cutensor_handle,U%Lop_desc(d),&
                         U%descB,b_ord,CUTENSOR_OP_IDENTITY,&
                         U%descWb(d),wb_ord,CUTENSOR_OP_IDENTITY,&
                         U%descLHS(d),lhs_ord,CUTENSOR_OP_IDENTITY,&
                         U%descLHS(d),lhs_ord,U%descCompute)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreatecontraction, B*A^TWA, mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_LHS_LS8()')
         endif
      ENDDO

      ! Estimate workspace sizes for B*AT^WA mult
      ALLOCATE(U%Lwork_size(ndof))
      do d=1,ndof
         cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%Lop_desc(d),&
                         U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%Lwork_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorestimateworkspacesize, B*A^TWA, mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_LHS_LS8()')
         endif
      enddo

      ! Create the plans for B*AT^WA mult
      ALLOCATE(U%Lplan(ndof))
      do d=1,ndof
         cutensor_status=cutensorcreateplan(cutensor_handle,U%Lplan(d),&
                        U%Lop_desc(d),U%plan_pref,U%Lwork_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreateplan, B*A^TWA, mode',d,&
                       '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_LHS_LS8()')
         endif
      enddo
#endif

      end subroutine Init_cutensor_LHS_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_cutensor_RHS_LS8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes LS8 type

      implicit none
      CLASS (LS8) :: U
      integer(8) :: p_ext(2),gb_ext(2),rhs_ext(2)
      integer(8) :: p_str(2),gb_str(2),rhs_str(2)
      integer    :: p_ord(2),gb_ord(2),rhs_ord(2)
      integer    :: d,ndof

#if ACC_ENABLED
      ndof=U%ndof
      U%dataType = CUTENSOR_R_64F
      U%descCompute = CUTENSOR_COMPUTE_DESC_64F

      ! Create tensor descriptors for G base and r.h.s.
      ALLOCATE(U%descGb(ndof),U%descRHS(ndof))
      DO d=1,ndof
         gb_ext=(/U%rkG*U%rkW*U%rksAp,U%nbas(d)/)
         gb_str=(/1,U%rkG*U%rkW*U%rksAp/)

         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
              U%descGb(d),2,gb_ext,gb_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, G base mode',d,&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_RHS_LS8()')
         endif

         rhs_ext=(/U%rkF,U%nbas(d)/)
         rhs_str=(/1,U%rkF/)
         cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
           U%descRHS(d),2,rhs_ext,rhs_str,U%dataType,cutensor_align)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorCreateTensorDescriptor, RHS',&
            '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_RHS_LS8()')
         endif
      ENDDO

!     Create tensor descriptor for P matrix <AF,G>
      p_ext=(/U%rkG*U%rkW*U%rksAp,U%rkF/)
      p_str=(/1,U%rkG*U%rkW*U%rksAp/)

      cutensor_status=cutensorCreateTensorDescriptor(cutensor_handle,&
        U%descP,2,p_ext,p_str,U%dataType,cutensor_align)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorCreateTensorDescriptor, P matrix',&
         '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_RHS_LS8()')
      endif

      ! Create operator descriptors for P*ATG mult
      ALLOCATE(U%Rop_desc(ndof))
      p_ord=(/2,1/)   ! rGrWrA',rF 
      gb_ord=(/2,3/)  ! rGrWrA',N*N
      rhs_ord=(/1,3/) ! rF,N*N

      DO d=1,ndof
         cutensor_status=cutensorcreatecontraction(cutensor_handle,U%Rop_desc(d),&
                         U%descP,p_ord,CUTENSOR_OP_IDENTITY,&
                         U%descGb(d),gb_ord,CUTENSOR_OP_IDENTITY,&
                         U%descRHS(d),rhs_ord,CUTENSOR_OP_IDENTITY,&
                         U%descRHS(d),rhs_ord,U%descCompute)
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorcreatecontraction, P*G mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_RHS_LS8()')
         endif
      ENDDO

      ! Create plan preference
      cutensor_status=cutensorcreateplanpreference(cutensor_handle,U%plan_pref,&
                      CUTENSOR_ALGO_DEFAULT,CUTENSOR_JIT_MODE_NONE)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensorcreateplanpreference; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_RHS_LS8()')
      endif

      ! Estimate workspace sizes for P*ATG mult
      ALLOCATE(U%Rwork_size(ndof))
      do d=1,ndof
         cutensor_status=cutensorestimateworkspacesize(cutensor_handle,U%Rop_desc(d),&
                         U%plan_pref,CUTENSOR_WORKSPACE_DEFAULT,U%Rwork_size(d))
         if (cutensor_status%stat.ne.0) then
            write(*,*) 'cutensorestimateworkspacesize, P*G, mode',d,&
                       '; exit = ',cutensor_status%stat
            call AbortWithError('Error in Init_cutensor_RHS_LS8()')
         endif
      enddo

      ! Create the plans for P*ATG mult
      ALLOCATE(U%Rplan(ndof))
      do d=1,ndof
         cutensor_status=cutensorcreateplan(cutensor_handle,U%Rplan(d),&
                        U%Rop_desc(d),U%plan_pref,U%Rwork_size(d))
         if (cutensor_status%stat .ne.0) then
            write(*,*) 'cutensorcreateplan, P*G, mode',d,&
                       '; exit = ',cutensor_status%stat
         call AbortWithError('Error in Init_cutensor_RHS_LS8()')
         endif
      enddo
#endif

      end subroutine Init_cutensor_RHS_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_LS8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (LS8) :: U

      IF (.not.U%setup) &
         call AbortWithError('Flush_LS8(): U not setup')

      IF (ALLOCATED(U%nbas))  DEALLOCATE(U%nbas)
      IF (ALLOCATED(U%m))     DEALLOCATE(U%m)
      IF (ALLOCATED(U%n))     DEALLOCATE(U%n)

      call Flush_cutensor_LHS_LS8(U)
      call Flush_cutensor_RHS_LS8(U)
      
      U%setup=.FALSE.

      end subroutine Flush_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_cutensor_LHS_LS8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes cutensor items for LHS in LS8 type

      implicit none
      TYPE (LS8) :: U
      integer :: d,ndof

      ndof=U%ndof

#if ACC_ENABLED
      ! No LHS cutensor items to flush if unweighted ALS is used
      if (U%algo.eq.1 .and. (U%useA.or.U%useW)) then
         ! Destroy plans
         DO d=1,ndof
            cutensor_status = cutensorDestroyPlan(U%Lplan(d))
         ENDDO

         ! Destroy operator (contraction) descriptors
         DO d=1,ndof
            cutensor_status = cutensorDestroyOperationDescriptor(U%Lop_desc(d))
         ENDDO

         ! Destroy tensor descriptors
         cutensor_status = cutensorDestroyTensorDescriptor(U%descB)
         DO d=1,ndof
            cutensor_status = cutensorDestroyTensorDescriptor(U%descWb(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descLHS(d))
         ENDDO

         IF (ALLOCATED(U%Lplan))      DEALLOCATE(U%Lplan)
         IF (ALLOCATED(U%Lwork_size)) DEALLOCATE(U%Lwork_size)
         IF (ALLOCATED(U%Lop_desc))   DEALLOCATE(U%Lop_desc)
         IF (ALLOCATED(U%descWb))     DEALLOCATE(U%descWb)
         IF (ALLOCATED(U%descLHS))    DEALLOCATE(U%descLHS)
      endif
#endif

      end subroutine Flush_cutensor_LHS_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_cutensor_RHS_LS8(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes cutensor items for RHS in LS8 type

      implicit none
      TYPE (LS8) :: U
      integer :: d,ndof

      ndof=U%ndof

#if ACC_ENABLED
      if (U%algo.eq.1) then
         ! Destroy plans
         DO d=1,ndof
            cutensor_status = cutensorDestroyPlan(U%Rplan(d))
         ENDDO

         ! Destroy plan preferences
         cutensor_status = cutensorDestroyPlanPreference(U%plan_pref)

         ! Destroy operator (contraction) descriptors
         DO d=1,ndof
            cutensor_status = cutensorDestroyOperationDescriptor(U%Rop_desc(d))
         ENDDO

         ! Destroy tensor descriptors
         cutensor_status = cutensorDestroyTensorDescriptor(U%descP)
         DO d=1,ndof
            cutensor_status = cutensorDestroyTensorDescriptor(U%descGb(d))
            cutensor_status = cutensorDestroyTensorDescriptor(U%descRHS(d))
         ENDDO

         IF (ALLOCATED(U%Rplan))      DEALLOCATE(U%Rplan)
         IF (ALLOCATED(U%Rwork_size)) DEALLOCATE(U%Rwork_size)
         IF (ALLOCATED(U%Rop_desc))   DEALLOCATE(U%Rop_desc)
         IF (ALLOCATED(U%descGb))     DEALLOCATE(U%descGb)
         IF (ALLOCATED(U%descRHS))    DEALLOCATE(U%descRHS)

      endif
#endif

      end subroutine Flush_cutensor_RHS_LS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetLHS_LinSys_CP8(U,B,ATWA,LHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets LHS for linear solver, general case
! wALS: pass W for 'ATWA'; B is <WF,F>
!   LS: pass A^T*A for 'ATWA'; B is <AF,AF>
!  wLS: pass A^T*W*A for 'ATWA'; B is <WAF,AF>

      implicit none
      TYPE (LS8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: ATWA
      real(kind=8), intent(in)  :: B(:)
      real(kind=8), allocatable, intent(out) :: LHS(:)
      integer, intent(in) :: d

      call InitLHS_LinSys_CP8(U,LHS,d)
      call AccumulateLHS_LinSys_CP8(U,B,ATWA,LHS,d)

      end subroutine GetLHS_LinSys_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitLHS_LinSys_CP8(U,LHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates + initializes LHS for linear solver

      implicit none
      TYPE (LS8), INTENT(IN)    :: U
      real(kind=8), allocatable, intent(out) :: LHS(:)
      integer, intent(in) :: d
      real(kind=8), parameter :: zero=0.d0
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('InitLHS_LinSys_CP8(): U not setup')

      IF (d.lt.1 .or. d.gt.U%ndof) THEN
          write(*,*) 'mode: ',d,' out of range [1,',U%ndof,']'
          call AbortWithError('Error in InitLHS_LinSys_CP8()')
      ENDIF

      ALLOCATE(LHS(U%rkF**2*U%m(d)*U%n(d)))

      call CPU_TIME(ti1)

      select case(U%algo)
         case(0)
            LHS(:)=zero
         case(1)
#if ACC_ENABLED
            !$acc enter data create(LHS)
            call SetVal_acc_r8(LHS,zero)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("InitLHS_LinSys_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("InitLHS_LinSys_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      lhs_time(mpirank+1)=lhs_time(mpirank+1)+ti2-ti1

      end subroutine InitLHS_LinSys_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_HS_CP8(HS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates LHS or RHS for linsys

      implicit none
      real(kind=8), allocatable, intent(inout) :: HS(:)

      IF (allocated(HS)) THEN
#if ACC_ENABLED
         IF (acc_is_present(HS)) THEN
            !$acc exit data delete(HS)
         ENDIF
#endif
         deallocate(HS)
      ENDIF

      end subroutine Flush_HS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateLHS_LinSys_CP8(U,B,ATWA,LHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets LHS for linear solver, general case
! wALS: pass W for 'ATWA'; B is <WF,F>
!   LS: pass A^T*A for 'ATWA'; B is <AF,AF>
!  wLS: pass A^T*W*A for 'ATWA'; B is <WAF,AF>

      implicit none
      TYPE (LS8), INTENT(IN)    :: U
      TYPE (CP8), INTENT(IN)    :: ATWA
      real(kind=8), intent(in)  :: B(:)
      real(kind=8) :: LHS(:)
      integer, intent(in) :: d
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('AccumulateLHS_ALS_CP8(): U not setup')

      IF (d.lt.1 .or. d.gt.U%ndof) THEN
          write(*,*) 'mode: ',d,' out of range [1,',U%ndof,']'
          call AbortWithError('Error in AccumulateLHS_LinSys_CP8()')
      ENDIF

      IF (ATWA%R().ne.U%rksA*U%rkW*U%rksAp) THEN
         write(*,*) 'ATWA has ',ATWA%R(),' terms but U requires ',&
         U%rksA*U%rkW*U%rksAp,' terms'
         call AbortWithError(&
         "AccumulateLHS_LinSys_CP8(): ATWA, U, rank mismatch")
      ENDIF

      IF (ATWA%D().ne.U%ndof) THEN
         write(*,*) 'ATWA has ',ATWA%D(),' modes but U has ',&
         U%ndof,' modes'
         call AbortWithError(&
         "AccumulateLHS_LinSys_CP8(): ATWA, U, # modes mismatch")
      ENDIF

      IF (ATWA%M(d).ne.U%m(d) .and. ATWA%N(d).ne.U%m(d)) THEN
         write(*,*) 'mode ',d,': ATWA has basis len [',ATWA%M(d),&
         ' x ',ATWA%N(d),']; must be [',U%m(d),' x ',U%m(d),']'
         call AbortWithError(&
         "AccumulateLHS_LinSys_CP8(): ATWA, U, len mismatch")
      ENDIF

      IF (SIZE(B).ne.U%rkW*U%rksA*U%rkF*U%rksAp*U%rkF) THEN
         write(*,*) 'SIZE(B) = ',SIZE(B),'; must be ',&
                    U%rkW*U%rksA*U%rkF*U%rksAp*U%rkF
         call AbortWithError("AccumulateLHS_LinSys_CP8(): wrong size B")
      ENDIF

      IF (SIZE(LHS).ne.U%rkF**2*U%m(d)*U%n(d)) THEN
         write(*,*) 'SIZE(LHS) = ',SIZE(LHS),'; must be ',&
         U%rkF**2*U%m(d)*U%n(d)
         call AbortWithError(&
         "AccumulateLHS_LinSys_CP8(): wrong size LHS")
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('CPLS8: AccumulateLHS')

      select case(U%algo)
         case(0)
            call AccumulateLHS_alg_cpu_CP8(U,B,ATWA,LHS,d)
         case(1)
#if ACC_ENABLED
            call AccumulateLHS_alg_tensor_CP8(U,B,ATWA,LHS,d)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("Error in AccumulateLHS_LinSys_CP8()")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("Error in AccumulateLHS_LinSys_CP8()")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      lhs_time(mpirank+1)=lhs_time(mpirank+1)+ti2-ti1

      end subroutine AccumulateLHS_LinSys_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetLHS_ALS_CP8(U,B,LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets LHS for linear solver, special case for unweighted ALS

      implicit none
      TYPE (LS8), INTENT(IN)    :: U
      real(kind=8), intent(in)  :: B(:)
      real(kind=8), allocatable, intent(out) :: LHS(:)
      integer :: r,rk
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('GetLHS_ALS_CP8(): U not setup')
     
      rk=U%rkF**2

      IF (SIZE(B).ne.rk) THEN
         write(*,*) 'B size = ',SIZE(B),'; U rkF**2 = ',rk,' must match'
         call AbortWithError("Error in GetLHS_ALS_CP8()")
      ENDIF

      ALLOCATE(LHS(rk))

      call CPU_TIME(ti1)
      call nvtx_start('CPLS8: GetLHS_ALS')

      select case(U%algo)
         case(0)
            LHS(:)=B(:)
         case(1)
#if ACC_ENABLED
            if (.not.acc_is_present(B)) call &
               AbortWithError('GetLHS_ALS_CP8(): B not on device')

            !$acc enter data create(LHS)
            !$acc data present(B,LHS)
            !$acc parallel loop gang vector async
            do r=1,rk
               LHS(r)=B(r)
            enddo
            !$acc end data
            !$acc wait
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("Error in GetLHS_CP8()")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("Error in GetLHS_CP8()")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      lhs_time(mpirank+1)=lhs_time(mpirank+1)+ti2-ti1

      end subroutine GetLHS_ALS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRHS_CP8(U,P,G,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets RHS for linear system used by ALS or linear solver, for mode d

      implicit none
      TYPE (LS8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: G
      real(kind=8), intent(in) :: P(:)
      real(kind=8), allocatable, intent(out) :: RHS(:)
      integer, intent(in) :: d

      call InitRHS_CP8(U,RHS,d)
      call AccumulateRHS_CP8(U,P,G,RHS,d)

      end subroutine GetRHS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitRHS_CP8(U,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates + initializes RHS for linear solver

      implicit none
      TYPE (LS8), INTENT(IN) :: U
      real(kind=8), allocatable, intent(out) :: RHS(:)
      integer, intent(in) :: d
      real(kind=8), parameter :: zero=0.d0
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('InitRHS_CP8(): U not setup')

      IF (d.lt.1 .or. d.gt.U%ndof) THEN
          write(*,*) 'mode: ',d,' out of range [1,',U%ndof,']'
          call AbortWithError('Error in InitRHS_CP8()')
      ENDIF

      ALLOCATE(RHS(U%rkF*U%nbas(d)))

      call CPU_TIME(ti1)

      select case(U%algo)
         case(0)
            RHS(:)=zero
         case(1)
#if ACC_ENABLED
            !$acc enter data create(RHS)
            call SetVal_acc_r8(RHS,zero)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("InitRHS_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("InitRHS_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      rhs_time(mpirank+1)=rhs_time(mpirank+1)+ti2-ti1

      end subroutine InitRHS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateRHS_CP8(U,P,G,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets RHS for linear system used by ALS or linear solver, for mode d
! ALS:  call with P(rG,rF),           G(rG,nbas(d))
! wALS: call with P(rG,rW ,rF),      WG(rG,rW,nbas(d))
! LS:   call with P(rG,rA',rF),    A^TG(rG,rA',nbas(d))
! wLS:  call with P(rG,rW,rA',rF) A^TWG(rG,rW,rA',nbas(d))

      implicit none
      TYPE (LS8), INTENT(IN) :: U
      TYPE (CP8), INTENT(IN) :: G
      real(kind=8), intent(in) :: P(:)
      real(kind=8) :: RHS(:)
      integer, intent(in) :: d
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

!     Error checking: make sure that other inputs match data in U
      IF (.not.U%setup) &
         call AbortWithError('AccumulateRHS_CP8(): U not setup')

      IF (d.lt.1 .or. d.gt.U%ndof) THEN
          write(*,*) 'mode: ',d,' out of range [1,',U%ndof,']'
          call AbortWithError('Error in AccumulateRHS_CP8()')
      ENDIF

      IF (G%R().ne.U%rkG*U%rkW*U%rksAp) THEN
         write(*,*) 'G has ',G%R(),' terms but U requires ',&
         U%rkG*U%rkW*U%rksAp,' terms'
         call AbortWithError("AccumulateRHS_CP8(): G, U, rank mismatch")
      ENDIF

      IF (G%D().ne.U%ndof) THEN
         write(*,*) 'G has ',G%D(),' modes but U has ',U%ndof,' modes'
         call AbortWithError("Error in AccumulateRHS_CP8()")
      ENDIF

      IF (G%MN(d).ne.U%nbas(d)) THEN
         write(*,*) 'mode ',d,': G: ',G%MN(d),' len; U: ',U%nbas(d),&
         ' len'
         call AbortWithError("AccumulateRHS_CP8(): G, U, len mismatch")
      ENDIF

      IF (SIZE(P).ne.U%rkF*U%rksAp*U%rkW*U%rkG) THEN
         write(*,*) 'SIZE(P) = ',SIZE(P),'; must be ',&
                    U%rkF*U%rksAp*U%rkW*U%rkG
         call AbortWithError("AccumulateRHS_CP8(): wrong size P")
      ENDIF

      IF (SIZE(RHS).ne.U%rkF*U%nbas(d)) THEN
         write(*,*) 'SIZE(RHS) = ',SIZE(RHS),'; must be ',&
         U%rkF*U%nbas(d)
         call AbortWithError("AccumulateRHS_CP8(): wrong size RHS")
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('CPLS8: AccumulateRHS')

      select case(U%algo)
         case(0)
            call AccumulateRHS_alg_cpu_CP8(U,P,G,RHS,d)
         case(1)
#if ACC_ENABLED
            call AccumulateRHS_alg_tensor_CP8(U,P,G,RHS,d)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("AccumulateRHS_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',U%algo
            call AbortWithError("AccumulateRHS_CP8(): wrong algorithm")
      end select

      call nvtx_stop()
      call CPU_TIME(ti2)
      rhs_time(mpirank+1)=rhs_time(mpirank+1)+ti2-ti1

      end subroutine AccumulateRHS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateLHS_alg_cpu_CP8(U,B,W,LHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! LHS routine, CPU version (alternate version)
! This version assumes that:
! B is <WAF,AF>, with indices (fastest) <- (rA,rW,rF,rA',rF) -> (slowest)
! W (or A^TA or A^TWA) is rearranged and coef-premultiplied, with indices 
!                             (fastest) <- (rA,rW,rA'(^T),N,N) -> (slowest)
! outputs LHS, with indices   (fastest) <- (rF,M,rF,N) -> (slowest)

      implicit none
      TYPE (LS8), INTENT(IN)      :: U
      TYPE (CP8), INTENT(IN)      :: W
      real(kind=8), intent(in)    :: B(:)
      real(kind=8), intent(inout) :: LHS(:)
      real(kind=8), parameter :: zero=0.d0
      integer, intent(in) :: d
      integer :: i,j,k,l,s,v,t,M,N,rF,rA,rAp,rW,ibas
      integer :: ll,il,ip,ib,i1,s1,stride
      integer :: NrF,rFNrF,rFrWrA,rArWrA,rWrA,rFrArWrA

      M=U%m(d)
      N=U%n(d)
      rF=U%rkF
      rA=U%rksA
      rAp=U%rksAp
      rW=U%rkW

      ll=W%look(1,d)
      NrF=N*rF
      rFNrF=rF*NrF
      rWrA=rW*rA
      rArWrA=rAp*rWrA
      rFrWrA=rF*rWrA
      rFrArWrA=rF*rArWrA

!     B is <WAF,AF>, or (fastest) <- rW,rA,rF,rA',rF' -> (slowest)
!     W is A^T*W*A, or  (fastest) <- rW,rA,rA',N,N'   -> (slowest)
!     LHS  has          (fastest) <-    rF',N,rF,N'   -> (slowest)
      !$omp parallel
      !$omp do private(i,j,k,l,il,ibas,ip,ib,i1,s1,stride) collapse(4)
      do l=1,M
         do i=1,rF
            do k=1,N
               do j=1,rF
                  if (U%useA) then ! LS
                     stride=k
                  else             ! wALS
                     stride=l
                  endif
                  il=(l-1)*rFNrF+(i-1)*NrF + (k-1)*rF + j ! Entry in LHS
                  i1=(j-1)*rFrArWrA + (i-1)*rWrA
                  ibas=(l-1)*M+stride   ! Basis entry
                  ib=ll+(ibas-1)*rArWrA
                  s1=0
                  do s=1,rAp
                     ip=s1+i1
                     do t=1,rA
                        do v=1,rW
                           ip=ip+1
                           ib=ib+1
                           LHS(il)=LHS(il)+B(ip)*W%base(ib)
                        enddo
                     enddo
                     s1=s1+rFrWrA
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp enddo
      !$omp end parallel

      end subroutine AccumulateLHS_alg_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateRHS_alg_cpu_CP8(U,P,G,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! RHS routine, CPU version (alternate format)
! This version assumes that:
! P is <WG,AF>, with indices (fastest) <- (rG,rW,rA,rF) -> (slowest)
! G (or WG or A^TWG) is rearranged and coef-premultiplied, with indices 
!                            (fastest) <- (rG,rW,rA,nbas) -> (slowest)

      implicit none
      TYPE (LS8), INTENT(IN)      :: U
      TYPE (CP8), INTENT(IN)      :: G
      real(kind=8), intent(in)    :: P(:)
      real(kind=8), intent(inout) :: RHS(:)
      real(kind=8), parameter :: zero=0.d0
      integer, intent(in) :: d
      integer :: i,j,k,ir,ip,ib,rkATWG,nbas

      rkATWG=U%rksAp*U%rkW*U%rkG
      nbas=U%nbas(d)

      !$omp parallel
      !$omp do private(i,j,k,ir,ip,ib) collapse(2)
      do i=1,nbas
        do j=1,U%rkF
           ir=(i-1)*U%rkF+j
           ip=(j-1)*rkATWG
           ib=G%look(1,d)+(i-1)*rkATWG
           do k=1,rkATWG
              ip=ip+1
              ib=ib+1
              RHS(ir)=RHS(ir)+P(ip)*G%base(ib)
           enddo
        enddo
      enddo
      !$omp enddo
      !$omp end parallel

      end subroutine AccumulateRHS_alg_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateLHS_alg_tensor_CP8(U,B,ATWA,LHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! LHS routine for Linear Systems, CPU version

      implicit none
      TYPE (LS8), INTENT(IN)      :: U
      TYPE (CP8), INTENT(IN)      :: ATWA
      real(kind=8), intent(in)    :: B(:)
      real(kind=8), intent(inout) :: LHS(:)
      real(kind=8), allocatable   :: work(:)
      real(kind=8), parameter     :: alpha=1.0, beta=1.0
      integer, intent(in) :: d
      integer :: ms

#if ACC_ENABLED
!     Make sure arrays are on device
      if (.not.ATWA%live_on_device) call AbortWithError(&
          'AccumulateLHS_alg_tensor_CP8(): ATWA not on device')
      if (.not.acc_is_present(B)) call AbortWithError(&
          'AccumulateLHS_alg_tensor_CP8(): B not on device')
      if (.not.acc_is_present(LHS)) call AbortWithError(&
          'AccumulateLHS_alg_tensor_CP8(): LHS not on device')

      ms=ATWA%MS(d)
      allocate(work(U%Lwork_size(d)))
      !$acc data create(work)
      !$acc host_data use_device(B,ATWA%base,work,LHS)
      cutensor_status = cutensorcontract(cutensor_handle,&
                           U%Lplan(d),alpha,B,ATWA%base(ms),beta,LHS,LHS,&
                           work,U%Lwork_size(d),acc_stream_g)
      !$acc end host_data
      if (cutensor_status%stat.ne.0) then
          write(*,*) &
         'cutensorcontract, base mode ',d,'; exit = ',cutensor_status%stat
         call AbortWithError(&
             'Error in AccumulateLHS_LinSys_alg_tensor_CP8()')
      endif
      !$acc wait
      !$acc end data
      deallocate(work)
#endif

      end subroutine AccumulateLHS_alg_tensor_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateRHS_alg_tensor_CP8(U,P,G,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! RHS routine, GPU (CUTENSOR) version

      implicit none
      TYPE (LS8), INTENT(IN)      :: U
      TYPE (CP8), INTENT(IN)      :: G
      real(kind=8), intent(in)    :: P(:)
      real(kind=8), intent(inout) :: RHS(:)
      real(kind=8), allocatable :: work(:)
      real(kind=8), parameter   :: alpha=1.0, beta=1.0
      integer, intent(in) :: d
      integer :: ms
       
#if ACC_ENABLED
!     Make sure arrays are on device
      if (.not.G%live_on_device) call AbortWithError(&
         'AccumulateRHS_alg_tensor_CP8(): G not on device')
      if (.not.acc_is_present(P)) call AbortWithError(&
         'AccumulateRHS_alg_tensor_CP8(): P not on device')
      if (.not.acc_is_present(RHS)) call AbortWithError(&
         'AccumulateRHS_alg_tensor_CP8(): RHS not on device')

      ms=G%MS(d)
      allocate(work(U%Rwork_size(d)))
      !$acc data create(work)
      !$acc host_data use_device(P,G%base,work,RHS)
      cutensor_status = cutensorcontract(cutensor_handle,&
                           U%Rplan(d),alpha,P,G%base(ms),beta,RHS,RHS,&
                           work,U%Rwork_size(d),acc_stream_g)
      !$acc end host_data
      if (cutensor_status%stat.ne.0) then
          write(*,*) &
         'cutensorcontract, base mode ',d,'; exit = ',cutensor_status%stat
         call AbortWithError('Error in AccumulateRHS_alg_tensor_CP8()')
      endif
      !$acc wait
      !$acc end data
      deallocate(work)
#endif

      end subroutine AccumulateRHS_alg_tensor_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSolninF_CP8(F,RHS,d,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies entries from RHS after linear solve into F

      implicit none
      TYPE (CP8) :: F
      real(kind=8), intent(in) :: RHS(:)
      integer, intent(in) :: d,algo
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

      IF (d.lt.1 .or. d.gt.F%D()) THEN
          write(*,*) 'mode: ',d,' out of range [1,',F%D(),']'
          call AbortWithError('Error in PutSolninF_CP8()')
      ENDIF

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            call PutSolninF_cpu_CP8(F,RHS,d)
         case(1)
#if ACC_ENABLED
            call PutSolninF_acc_CP8(F,RHS,d)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("PutSolninF_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("PutSolninF_CP8_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      putsoln_time(mpirank+1)=putsoln_time(mpirank+1)+ti2-ti1

      end subroutine PutSolninF_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSolninF_cpu_CP8(F,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies entries from RHS after linear solve into F, CPU algorithm

      implicit none
      TYPE (CP8) :: F
      real(kind=8), intent(in) :: RHS(:)
      integer, intent(in) :: d
      integer :: r,rk,i,j,l,mn

      rk=F%R()
      mn=F%MN(d)

      do r=1,rk
         l=F%look(r,d)
         j=r
         do i=1,mn
            F%base(l+i)=RHS(j)
            j=j+rk
         enddo
      enddo

      end subroutine PutSolninF_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSolninF_acc_CP8(F,RHS,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies entries from RHS after linear solve into F, OpenACC algorithm

      implicit none
      TYPE (CP8) :: F
      real(kind=8), intent(in) :: RHS(:)
      integer, intent(in) :: d
      integer :: r,rk,i,j,l,mn

#if ACC_ENABLED
!     Make sure arrays are on device
      if (.not.F%live_on_device) call &
         AbortWithError('PutSolninF_acc_CP8(): F not on device')
      if (.not.acc_is_present(RHS)) call &
         AbortWithError('PutSolninF_acc_CP8(): RHS not on device')

      rk=F%R()
      mn=F%MN(d)

      !$acc data present(F%base,F%look,RHS)
      !$acc parallel loop gang vector private(l,j) async
      do r=1,rk
         do i=1,mn
            l=F%look(r,d)
            j=(i-1)*rk + r
            F%base(l+i)=RHS(j)
         enddo
      enddo
      !$acc wait
      !$acc end data
#endif

      end subroutine PutSolninF_acc_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrepG_allmodes_CP8(G,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Premultiplies coef*base(:) for all modes and rank-transposes base
! This wrapper performs operation for all ndof modes

      implicit none
      TYPE (CP8), INTENT(INOUT) :: G
      integer, intent(in) :: algo
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=G%D()
      allocate(modes(ndof))
      modes(:)=.TRUE.
      call PrepG_CP8(G,modes,algo)
      deallocate(modes)

      end subroutine PrepG_allmodes_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrepG_1mode_CP8(G,d,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Premultiplies coef*base(:) for all modes and rank-transposes base
! This wrapper performs operation for selected mode d

      implicit none
      TYPE (CP8), INTENT(INOUT) :: G
      integer, intent(in) :: d,algo
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=G%D()

      IF (d.lt.1 .or. d.gt.ndof) THEN
         write(*,*) 'Mode ',d,' out of range: [',1,',',ndof,']'
         call AbortWithError('PrepG_1mode_CP8(): wrong mode d')
      ENDIF

      allocate(modes(ndof))
      modes(:)=.FALSE.
      modes(d)=.TRUE.
      call PrepG_CP8(G,modes,algo)
      deallocate(modes)

      end subroutine PrepG_1mode_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrepG_general_CP8(G,modes,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Premultiplies coef*base(:) for all modes and rank-transposes base

      implicit none
      TYPE (CP8), INTENT(INOUT) :: G
      logical, intent(in) :: modes(:)
      integer, intent(in) :: algo
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPLS8_Module()

      IF (SIZE(modes).ne.G%D()) THEN
         write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',G%D()
         call AbortWithError('PrepG_CP8(): bad size for modes array')
      ENDIF

      call CPU_TIME(ti1)

      select case(algo)
         case(0)
            call PrepG_cpu_CP8(G,modes)
         case(1)
#if ACC_ENABLED
            call PrepG_acc_CP8(G,modes)
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("PrepG_CP8(): wrong algorithm")
#endif
         case default
            write(*,*) 'Unrecognized algo choice:',algo
            call AbortWithError("PrepG_CP8(): wrong algorithm")
      end select

      call CPU_TIME(ti2)
      prepg_time(mpirank+1)=prepg_time(mpirank+1)+ti2-ti1

      end subroutine PrepG_general_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrepG_cpu_CP8(G,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Premultiplies coef*base(:) for all modes and rank-transposes base

      implicit none
      TYPE (CP8), INTENT(INOUT) :: G
      logical, intent(in) :: modes(:)
      integer :: d,ndof,r,rk,i,bs,bf,ms,mf,mn,ll
      real(kind=8), allocatable :: base(:)

      ndof=G%D()
      rk=G%R()

!     Scale the base for all modes by the coef
      do d=1,ndof
         if (modes(d)) then
            do r=1,rk
               bs=G%BS(r,d)
               bf=G%BF(r,d)
               G%base(bs:bf)=G%coef(r)*G%base(bs:bf)
            enddo
         endif
      enddo

!     Rank-transpose the base
      do d=1,ndof
         if (modes(d)) then
            ms=G%MS(d)
            mf=G%MF(d)
            mn=G%MN(d)
            allocate(base(mn*rk))
            do r=1,rk
               ll=G%look(r,d)
               do i=1,mn
                  base((i-1)*rk+r)=G%base(ll+i)
               enddo
            enddo
            G%base(ms:mf)=base(1:mn)
            deallocate(base)
         endif
      enddo

      end subroutine PrepG_cpu_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrepG_acc_CP8(G,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Premultiplies coef*base(:) for all modes and rank-transposes base
! OpenACC (GPU) version

      implicit none
      TYPE (CP8), INTENT(INOUT) :: G
      logical, intent(in) :: modes(:)
      integer :: d,ndof,r,rk,i,ib,ll,ld,mn,mxmn,rkmn
      real(kind=8), allocatable :: base(:)

#if ACC_ENABLED
!     Make sure array is on device
      if (.not.G%live_on_device) call &
         AbortWithError('PrepG_acc_CP8(): G not on device')

      ndof=G%D()
      rk=G%R()
      mxmn=MAXVAL(G%nbas)
      allocate(base(mxmn*rk))

      !$acc data present(G%base,G%coef,G%look) create(base)
      do d=1,ndof
         if (modes(d)) then
            mn=G%nbas(d)
            rkmn=rk*mn
            ld=G%look(1,d)

!           Copy base + rank transpose into temp array, scaling by coefs
            !$acc parallel loop gang vector async collapse(2)
            do r=1,rk
               do i=1,mn
                  ll=G%look(r,d)+i
                  ib=(i-1)*rk+r
                  base(ib)=G%base(ll)*G%coef(r)
               enddo
            enddo

!           Copy temp array back into G
            !$acc parallel loop gang vector async
            do i=1,rkmn
               G%base(ld+i)=base(i)
            enddo
         endif
      enddo
      !$acc end data
      !$acc wait
      deallocate(base)
#endif

      end subroutine PrepG_acc_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetVal_acc_r8(v,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets entries of v to val, OpenACC kernel

      implicit none
      real(kind=8), intent(inout) :: v(:)
      real(kind=8), intent(in) :: val
      integer :: i,n

#if ACC_ENABLED
      if (.not.acc_is_present(v)) call &
         AbortWithError('SetVal_acc_r8(): v not on device')

      n=SIZE(v)

      !$acc data present(v)
      !$acc parallel loop gang vector async
      do i=1,n
         v(i)=val
      enddo
      !$acc end data
      !$acc wait
#endif

      end subroutine SetVal_acc_r8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPLS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
