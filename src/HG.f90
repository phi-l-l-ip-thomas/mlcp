!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module HG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains the routines for CP-matrix diagonalization

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE RESTART
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE BLOCKUTILS
      USE LINSOLVER
      USE ALSOO
      USE ALSDRVR
      USE CPMATH
      USE HGORTHO

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function HGGuess(H,rk,cols) result (U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Guesses a starting matrix for CP diagonalization
! First term is identity; small random terms follow

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP) :: U
      integer, intent(in),optional :: cols(:)
      integer, intent(in)  :: rk
      logical, allocatable :: sym(:)

      allocate(sym(SIZE(H%sym)))
      sym(:)=.FALSE.

      IF (present(cols)) THEN
!         U=IdentityCPMatrix(H%rows,cols,sym)
         U=RandomCP(rk,H%rows,cols,H%sym)
      ELSE
!         U=IdentityCPMatrix(H%rows,H%cols,sym)
         U=RandomCP(H,rk)
      ENDIF
      call AugmentVWithRandom(U,rk) 

      call U%print()

      deallocate(sym)

      end function HGGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HGPowrRecurse(H,U,npow,nals,ishift,Eshift,Wt)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies Hamiltonian (H-E)^npow*U to matrix U
! Gotta start somewhere...

      implicit none
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN), OPTIONAL :: Wt
      TYPE (CP) :: V
      integer, intent(in) :: npow,nals,ishift
      real*8, intent(in)  :: Eshift
      character(len=64)   :: tag
      integer :: i,rk,redstat

      rk=U%R()

      write(*,*)
      write(*,*) '*** HGPowrRecurse() ***'
      write(*,*)

      write(*,*) 'Iterations:'

      do i=1,npow
!        Matrix-matrix product: V <-- (H-E1)*U; reduce rank U <-- V
         call CPMM(H,1,Eshift,.FALSE.,U,0,0.d0,.FALSE.,V)
         write(tag,'(A,I0)') 'HGPowrRecurse(): U<-W, i = ',i

         IF (present(Wt)) THEN
!           Transpose U,V; do column-weighted reduction; transpose back
            call CPMatrixTranspose(V)
            call CPMatrixTranspose(U)
            redstat=ALS_reduce(U,V,Wt,nals,tag)
            call CPMatrixTranspose(U)
            call CPMatrixTranspose(V)
         ELSE
            redstat=ALS_reduce(U,V,nals,tag)
         ENDIF

         IF (present(Wt)) THEN
            call OrthogCPmat(U,Wt,10)
         ELSE
            call OrthogCPmat(U,10)
         ENDIF

!!!      Check the QR decomp
         write(*,*) 'Here is QTQ, reduced to rank-1'
         call ShowRank1UTU(U)
         write(*,*) 'Here is R, reduced to rank-1'
         call ShowRank1R(V,U)
         call FlushCP(V)
         call ShowRQ(H,U)
!!!
      enddo

      call AbortWithError("Done (HGPowrRecurse)")

      end subroutine HGPowrRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module HG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
