!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MYACC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MPI wrapper functions

      USE ERRORTRAP
      USE MYMPI
      USE iso_c_binding
#if ACC_ENABLED
      USE OPENACC
      USE CUDAFOR
      USE CUBLAS_V2
      USE CUTENSOR_V2
      USE CUSOLVERDN
      USE NVTX
#endif
      implicit none
#if ACC_ENABLED
      integer(kind=cuda_stream_kind) :: acc_stream_g
      integer(4), parameter :: cutensor_align=256
      type(cublasHandle)     :: cublas_handle
      type(cusolverDnHandle) :: cusolver_handle
      type(cutensorstatus)   :: cutensor_status
      type(cutensorHandle)   :: cutensor_handle
      logical :: device_setup=.FALSE.

#endif

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine init_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up device (GPU)

      implicit none
      integer :: thegpu,ndev,err

#if ACC_ENABLED
      ndev=acc_get_num_devices(acc_device_nvidia)
      IF (mpirank.eq.mpi_prnt_rank) then
         write(*,'(X,A,X,I0,X,A,X,I0))') 'OpenACC enabled;',ndev,&
         'devices visible to rank',mpi_prnt_rank
      ENDIF

      thegpu=MOD(mpirank,ndev)
      call acc_set_device_num(thegpu,acc_device_nvidia)
!      write(*,'(X,3(A,I0),A)') 'MPI rank ',mpirank,' on device (',thegpu+1,'/',ndev,')'

!     Cuda stream where work takes place
      err = cudaStreamCreateWithFlags(acc_stream_g, cudaStreamNonBlocking)
      if (err.ne.0) call &
         AbortWithError('init_device(): error creating cuda stream')
      call acc_set_cuda_stream(0, acc_stream_g)

!     cuBLAS handle
      err = cublasCreate(cublas_handle)
      if (err.ne.0) call &
         AbortWithError('init_device(): error creating cublas handle')
      err = cublasSetStream(cublas_handle,acc_stream_g)
      if (err.ne.0) call &
         AbortWithError('init_device(): error-map cublas to stream')

!     cuSolver handle
      err = cusolverDnCreate(cusolver_handle)
      if (err.ne.0) call &
         AbortWithError('init_device(): error creating cusolver handle')
      err = cusolverDnSetStream(cusolver_handle,acc_stream_g)
      if (err.ne.0) call &
         AbortWithError('init_device(): error-map cusolver to stream')

!     cuTensor handle
      cutensor_status = cutensorCreate(cutensor_handle)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensor_status error:',cutensor_status%stat
         call &
         AbortWithError('init_device(): error creating cutensor handle')
      endif

      device_setup = .TRUE.
#endif

      end subroutine init_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine flush_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up device (GPU)

      implicit none
      integer :: err

#if ACC_ENABLED
!     cuBLAS handle
      err = cublasDestroy(cublas_handle)
      if (err.ne.0) call &
         AbortWithError('flush_device(): error flushing cublas_handle')

!     cuSolver handle
      err = cusolverDnDestroy(cusolver_handle)
      if (err.ne.0) call &
         AbortWithError('init_device(): error flushing cusolver handle')

!     cuTensor handle
      cutensor_status = cutensorDestroy(cutensor_handle)
      if (cutensor_status%stat.ne.0) then
         write(*,*) 'cutensor_status error:',cutensor_status%stat
         call &
         AbortWithError('flush_device(): error flushing cutensor handle')
      endif

      device_setup = .FALSE.
#endif

      end subroutine flush_device

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine nvtx_start(tag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper to initiate an nvtx region

      implicit none
      character(len=*), intent(in) :: tag

#if ACC_ENABLED
      call nvtxStartRange(tag)
#endif

      end subroutine nvtx_start

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine nvtx_stop()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper to initiate an nvtx region

      implicit none

#if ACC_ENABLED
      call nvtxEndRange()
#endif

      end subroutine nvtx_stop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MYACC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
