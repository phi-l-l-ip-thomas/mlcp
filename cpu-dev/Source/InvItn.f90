!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE INVITN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains inverse iteration

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE ALSUTILS
      USE ALSDRVR
      USE LINSOLVER

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InvItn_wrapper(H,F,shift,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Inverse iteration solver

      implicit none
      TYPE (CP), INTENT(IN)    :: H
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: G
      integer, intent(in)   :: nals
      real*8, intent(inout) :: shift


!     First version: ordinary solver
!      G=CopyCP(F)
!      G=RandomCP(F,F%R())
!      call NORMALIZE(G)
!      call LinSolver_1(H,F,G,nals,1,shift,.TRUE.) ! Beylkin vsn
!      call LinSolver_2(H,F,G,nals,1,shift,.TRUE.) ! Espig vsn
!      call FlushCP(G)

!     Second version: intertwined solver (shift update optional)
      call LintertwinedInvItn_1(H,F,nals,shift)
!      call LintertwinedInvItn_1(H,G,nals,shift)
!      call FlushCP(G)

      write(*,*)

      end subroutine InvItn_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InvItn_shift_test(H,F,thisshift,nextshift,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Inverse iteration solver test: make shifts closer to target eigenvalue
! to see how close the shift is to true value before solver fails
! PT 02/20/2023

      implicit none
      TYPE (CP), INTENT(IN)    :: H
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: G
      integer, intent(in)   :: nals
      real*8, intent(in) :: thisshift,nextshift
      real*8 :: shift,fac
      integer :: i, nshift

      nshift=8
      fac=2.d0
     
      DO i=1,nshift
         shift=thisshift + (nextshift-thisshift)/fac 
         write(*,*) "SHIFT: this = ",thisshift,"; next = ",nextshift,&
                    "; shift = ",shift
!        Use a version which does not update shift since we are 
!        modulating it here


!     First version: ordinary solver
       G=CopyCP(F)
!      G=RandomCP(F,F%R())
!      call NORMALIZE(G)
!      call LinSolver_1(H,F,G,nals,1,shift,.TRUE.) ! Beylkin vsn
!      call LinSolver_2(H,F,G,nals,1,shift,.TRUE.) ! Espig vsn
!      call FlushCP(G)

!     Second version: intertwined solver (shift update optional)
!      call LintertwinedInvItn_1(H,F,nals,shift)
      call LintertwinedInvItn_1(H,G,nals,shift)
      call FlushCP(G)

        write(*,*)
        fac=2.d0*fac

      ENDDO

      end subroutine InvItn_shift_test


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INVITN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
