!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module CHEB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MODHVEC
      USE MODVECVEC
      USE SEPDREPN
      USE HAMILSETUP
      USE REDUCTION

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetChebParams(avec,bounds,av,df)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts the parameters for the Chebyshev iteration from a 
! previously-computed list of eigenvalues

      implicit none
      real*8, intent(in)  :: avec(:),bounds(2)
      real*8, intent(out) :: av(:),df(:)
      integer :: i,nev
      real*8  :: a,eps

      eps=1.d-14
      nev=SIZE(avec)

      DO i=1,nev
         IF (i.lt.nev) THEN 
!           lower bound just below the 1st unwanted eigenvalue
            a=avec(i+1)-eps
         ELSEIF (nev.gt.1) THEN
!           highest eigenvalue in the block: guess a lower bound
            a=2*avec(nev)-avec(nev-1)
         ELSE
!           only one eigenvalue: lower bound just above this eigenvalue
            a=avec(1)+eps
         ENDIF
         av(i)=(a+bounds(2))/2
         df(i)=(bounds(2)-a)/2
      ENDDO

!      write(*,'(4X,A,2(f10.4,A)/)') 'Chebyshev bounds: [',a,',',b,']'

      end subroutine GetChebParams

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebRecurse(u,H,av,df,E0,npow)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Performs the Chebyshev recursion on a vector in CP-format

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: u
      TYPE (CPvec) :: v,w
      integer, intent(in) :: npow
      real*8,  intent(in) :: av,df,E0
      integer :: i
      real*8  :: sig,sigold,ts1,tc

      IF (npow.eq.0) RETURN

!     Initialize recurrence by computing sig and v=(sig/df)*(H-av)*u
      tc=2.d0/df
      sig=df/(av-E0)
      ts1=2.d0/sig
      call PRODHV(u,v,H,1,av)
      v%coef=v%coef*(sig/df)
      call reduc(v)

!     Chebyshev recurrence
      DO i=2,npow
         sigold=sig
         sig=1.d0/(ts1-sigold)
!        w = 2*sigma_n+1/c*(H-d*I)*v - sigma_n*sigma_n-1*u
         call PRODHV(v,w,H,1,av)
         call SUMVECVEC(w,tc*sig,u,-sig*sigold)

!        Recursion: u <- v, reduce v <- w
         call FlushCPvec(u)
         call CopyWtoV(u,v)
         call reduc(v,w)
         call FlushCPvec(w)
      ENDDO

!     Copy u <- v after final step since u is returned; normalize
      call ReplaceVwithW(u,v)
      call NORMCOEF(u)

      end subroutine ChebRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module CHEB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
