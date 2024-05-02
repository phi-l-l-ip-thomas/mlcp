!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE RANDOM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the random number generator

      USE ERRORTRAP
      USE UTILS
      USE MYMPI

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitRandom(t1,d,t,rsinp,rs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes random seed

      implicit none
      real*8, intent(in)  :: t1
      integer, intent(in) :: d(3), t(3)
      integer, intent(in) :: rsinp(:)
      integer, allocatable, intent(out) :: rs(:)
      integer :: i,j,m,n,tp
      real*8  :: rt
      character(len=64) :: frmt

      m=SIZE(rsinp)

      call random_seed(SIZE=n)
      allocate(rs(n))
      rs(:)=1

      IF (ANY(rsinp.ne.0)) THEN
!        Copy seed from input file
         do i=1,min(m,n)
            rs(i)=rsinp(i)
         enddo

         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,'(/X,A/)') 'Random seed used from input file...'

      ELSE
!        Fill seed using values derived from date and time
         j=0
         do i=1,3
            if (j.ge.n) exit
            j=j+1
            rs(j)=t(i)
         enddo
         do i=1,3
            if (j.ge.n) exit
            j=j+1
            rs(j)=d(i)
         enddo
!        Convert t1 to a large int and fill remainder of seed
         rt=t1*10**(8-int(log10(t1)))
         j=INT(rt)
         tp=1
         do i=1,n-6
            rs(i+6)=mod(j,tp+1)+1
            tp=2*tp
         enddo

         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(frmt,'(A,I0,A)') '(/X,A,',n,'(X,I0)/)'
            write(*,frmt) 'Random seed generated: ',(rs(i),i=1,n)
         ENDIF

      ENDIF

!     Send seed to all MPI ranks and set
      call bcast(rs)
      call random_seed(PUT=rs)

      end subroutine InitRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE RANDOM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
