!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKPOWER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TISE solver using the block power method

      use ERRORTRAP
      use UTILS
      use MODHVEC
      use MODVECVEC
      use SEPDREPN
      use HAMILSETUP
      use REDUCTION

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetBlockShift(avec,bounds,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Calculates approximate optimal E-shift for block power method

      implicit none
      real*8, intent(in)  :: avec(:),bounds(2)
      real*8, intent(out) :: Eshift
      integer :: nev
      real*8  :: a

      nev=SIZE(avec)

!     Estimate the first eigenvalue outside the block
      IF (nev.gt.1) THEN
         a=2*avec(nev)-avec(nev-1)
      ELSE
         a=2*avec(1)-bounds(1)
      ENDIF

!     Shift is average of estimate above and upper bound of spectrum
      Eshift=0.5*(a+bounds(2))

!      write(*,'(/X,A,f16.8/)') 'Eshift = ',Eshift

      end subroutine GetBlockShift

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PowrRecurse(v,H,npow,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies Hamiltonian (H-E)^npow*v to vector v
! The ishift parameter controls the shifting

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: v
      TYPE (CPvec) :: w
      integer, intent(in) :: npow
      real*8, intent(in)  :: Eshift
      integer :: i

      do i=1,npow
!        w <-- H*v; then v <-- w
         call PRODHV(v,w,H,1,Eshift)

!        Reduce and normalize v
         call reduc(v,w)
         call FlushCPvec(w)
         call NORMCOEF(v)
      enddo

      end subroutine PowrRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      END MODULE BLOCKPOWER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
