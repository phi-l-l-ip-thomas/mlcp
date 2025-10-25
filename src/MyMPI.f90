!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MYMPI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MPI wrapper functions

      USE ERRORTRAP
      USE MPI

      implicit none 

      INTERFACE bcast
        MODULE PROCEDURE bcast_string_ch
        MODULE PROCEDURE bcast_l_0d,bcast_l_1d
        MODULE PROCEDURE bcast_i4_0d,bcast_i4_1d,bcast_i4_2d,bcast_i4_3d
        MODULE PROCEDURE bcast_r4_0d,bcast_r4_1d,bcast_r4_2d,bcast_r4_3d
        MODULE PROCEDURE bcast_r8_0d,bcast_r8_1d,bcast_r8_2d,bcast_r8_3d
      END INTERFACE

      integer, parameter :: mpi_comm_wd=MPI_COMM_WORLD
      integer, parameter :: mpi_i4=MPI_INTEGER
      integer, parameter :: mpi_i8=MPI_INTEGER8
      integer, parameter :: mpi_r4=MPI_REAL4
      integer, parameter :: mpi_r8=MPI_DOUBLE_PRECISION
      integer, parameter :: mpi_ch=MPI_CHARACTER
      integer, parameter :: mpi_lg=MPI_LOGICAL
      integer, parameter :: mpi_sm=MPI_SUM
      integer, parameter :: mpi_prnt_rank = 0
      integer, parameter :: mpi_io_rank = 0
      integer :: mpirank, mpinodes, mpierr, nomp_threads

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine prepare_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
#if OMP_ENABLED
      integer, external :: omp_get_max_threads
#endif

      call mpi_init(mpierr)
      if (mpierr.ne.0) &
         call AbortWithError('prepare_mpi(): mpi_init failed')

      call mpi_comm_rank(mpi_comm_wd, mpirank,  mpierr)
      if (mpierr.ne.0) &
         call AbortWithError('prepare_mpi(): mpi_comm_rank failed')

      call mpi_comm_size(mpi_comm_wd, mpinodes, mpierr)
      if (mpierr.ne.0) &
         call AbortWithError('prepare_mpi(): mpi_comm_size failed')

      if (mpirank.eq.mpi_prnt_rank) then
         write(*,'(X,A,I0,A)') 'MLCP runs on ',mpinodes,' MPI processes'
#if OMP_ENABLED
         nomp_threads=omp_get_max_threads()
         write(*,'(X,A,X,I0,X,A)') 'with',nomp_threads,&
          'OpenMP threads/process'
#else
         nomp_threads=1
#endif
      endif

      end subroutine prepare_mpi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sync_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      call mpi_barrier(mpi_comm_wd, mpierr)
      if (mpierr.ne.0) &
         call AbortWithError('sync_mpi(): mpi_barrier failed')

      end subroutine sync_mpi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine finalize_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      call mpi_finalize(mpierr)
      if (mpierr.ne.0) &
         call AbortWithError('finalize_mpi(): mpi_finzlize failed')

      end subroutine finalize_mpi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_string_ch(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      character(len=*) :: p
      integer :: n
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n=len(p)
      call bcast_ch(p,n,iin)

      end subroutine bcast_string_ch

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_ch(p,n,iin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n,iin,ierr
      character(len=n) :: p

      call sync_mpi
      call mpi_bcast(p,n,mpi_ch,iin,mpi_comm_wd,ierr)
      if (ierr.ne.0) &
         call AbortWithError('bcast_ch(): mpi_bcast failed')
      call sync_mpi

      end subroutine bcast_ch

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_l_0d(m,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      logical, intent(inout) :: m
      logical :: ma(1)
      integer :: n
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      ma=.false.
      if (mpirank.eq.0) ma=m
      n=1
      call bcast_l(ma,n,iin)
      m=ma(1)

      end subroutine bcast_l_0d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_l_1d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n
      logical :: p(:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n=SIZE(p)
      call bcast_l(p,n,iin)

      end subroutine bcast_l_1d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_l(p,n,iin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n,iin,ierr
      logical :: p(n)

      call sync_mpi
      call mpi_bcast(p,n,mpi_lg,iin,mpi_comm_wd,ierr)
      if (ierr.ne.0) &
         call AbortWithError('bcast_l(): mpi_bcast failed')
      call sync_mpi

      end subroutine bcast_l

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_i4_0d(m,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      integer, intent(inout) :: m
      integer :: ma(1)
      integer :: n
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      ma=0
      if (mpirank.eq.0) ma=m
      n=1
      call bcast_i4(ma,n,iin)
      m=ma(1)

      end subroutine bcast_i4_0d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_i4_1d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n
      integer :: p(:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n=SIZE(p)
      call bcast_i4(p,n,iin)

      end subroutine bcast_i4_1d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_i4_2d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(2)
      integer :: p(:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      call bcast_i4(p,PRODUCT(n),iin)

      end subroutine bcast_i4_2d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_i4_3d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(3)
      integer :: p(:,:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      n(3)=SIZE(p,3)
      call bcast_i4(p,PRODUCT(n),iin)

      end subroutine bcast_i4_3d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_i4(p,n,iin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n,iin,ierr,p(n)

      call sync_mpi
      call mpi_bcast(p,n,mpi_i4,iin,mpi_comm_wd,ierr)
      if (ierr.ne.0) &
         call AbortWithError('bcast_i(): mpi_bcast failed')
      call sync_mpi

      end subroutine bcast_i4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r4_0d(m,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real*4, intent(inout) :: m
      real*4  :: ma(1)
      integer :: n
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      ma=0.d0
      if (mpirank.eq.0) ma=m
      n=1
      call bcast_r4(ma,n,iin)
      m=ma(1)

      end subroutine bcast_r4_0d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r4_1d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n
      real*4  :: p(:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n=SIZE(p)
      call bcast_r4(p,n,iin)

      end subroutine bcast_r4_1d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r4_2d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(2)
      real*4  :: p(:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      call bcast_r4(p,PRODUCT(n),iin)

      end subroutine bcast_r4_2d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r4_3d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(3)
      real*4  :: p(:,:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      n(3)=SIZE(p,3)
      call bcast_r4(p,PRODUCT(n),iin)

      end subroutine bcast_r4_3d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r4(p,n,iin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n,iin,ierr
      real*4  :: p(n)

      call sync_mpi
      call mpi_bcast(p,n,mpi_r4,iin,mpi_comm_wd,ierr)
      if (ierr.ne.0) &
         call AbortWithError('bcast_r4(): mpi_bcast failed')
      call sync_mpi

      end subroutine bcast_r4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r8_0d(m,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      real*8, intent(inout) :: m
      real*8  :: ma(1)
      integer :: n
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      ma=0.d0
      if (mpirank.eq.0) ma=m
      n=1
      call bcast_r8(ma,n,iin)
      m=ma(1)

      end subroutine bcast_r8_0d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r8_1d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n
      real*8  :: p(:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n=SIZE(p)
      call bcast_r8(p,n,iin)

      end subroutine bcast_r8_1d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r8_2d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(2)
      real*8  :: p(:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      call bcast_r8(p,PRODUCT(n),iin)

      end subroutine bcast_r8_2d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r8_3d(p,iproc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n(3)
      real*8  :: p(:,:,:)
      integer, optional :: iproc
      integer :: iin

      iin=0
      if (present(iproc)) iin=iproc

      n(1)=SIZE(p,1)
      n(2)=SIZE(p,2)
      n(3)=SIZE(p,3)
      call bcast_r8(p,PRODUCT(n),iin)

      end subroutine bcast_r8_3d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine bcast_r8(p,n,iin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer :: n,iin,ierr
      real*8  :: p(n)

      call sync_mpi
      call mpi_bcast(p,n,mpi_r8,iin,mpi_comm_wd,ierr)
      if (ierr.ne.0) &
         call AbortWithError('bcast_r8(): mpi_bcast failed')
      call sync_mpi

      end subroutine bcast_r8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine calc_mpi_partition(b,sz,os,szmx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Helper function for computing size and offset of b objects partitioned
! over MPI ranks

      implicit none
      integer, intent(in)  :: b
      integer, intent(out) :: sz,os,szmx
      integer :: bp,mp

      bp=b/mpinodes
      mp=mod(b,mpinodes)

      szmx=bp
      if (mp.gt.0) szmx=szmx+1

      sz=bp
      if (mpirank.lt.mp) sz=sz+1

      os=mpirank*bp+min(mpirank,mp)

      end subroutine calc_mpi_partition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Get_MPI_Timings(tag,timings)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      character(*), intent(in) :: tag
      character(LEN=36) :: ftag
      real(kind=8) :: timings(:)
      real(kind=8) :: maxtime,mintime,avetime
      integer, allocatable :: mstarts(:), mwidths(:)
      integer :: i,ierr

      allocate(mstarts(mpinodes),mwidths(mpinodes))
      do i=1,mpinodes
         mstarts(i)=i-1
         mwidths(i)=1
      enddo

!     Gather the timings from the MPI ranks
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,timings,&
                          mwidths,mstarts,mpi_r8,mpi_comm_wd,ierr)

      IF (mpirank.eq.mpi_prnt_rank) THEN
         mintime=MINVAL(timings)
         maxtime=MAXVAL(timings)
         avetime=SUM(timings)/mpinodes
         write(ftag,*) TRIM(ADJUSTL(tag))
         write(*,'(X,A,A,3(X,f14.3))') &
         ftag,':',mintime,maxtime,avetime
      ENDIF
      deallocate(mstarts,mwidths)

      end subroutine Get_MPI_Timings

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MYMPI

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
