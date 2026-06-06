!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKPOWER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TISE solver using the block power method

      USE ERRORTRAP
      USE UTILS
      USE CHEBLIB
      USE LINALG
      USE CPMMM
      USE MODVECVEC
      USE SEPDREPN
      USE REDUCTION
      USE LINSOLVER
      USE BLOCKUTILS

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetBlockShift(avec,bounds,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates approximate optimal E-shift for block power method

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

      end subroutine GetBlockShift

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InverseRecurse(v,H,npow,nals,Eshift,which)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves (H-E1)v=w for vector v

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow,nals,which
      real*8, intent(in)  :: Eshift
      integer :: i

      write(*,*) '*----------------*'
      write(*,*) 'Eshift = ',Eshift
      write(*,*) '*----------------*'

      do i=1,npow
         write(*,'(A,I0,A,I0,A)') 'ipow = (',i,'/',npow,')'
         w=CopyCP(v)
         call LinSolver_alg(H,w,v,nals,1,Eshift,which,.TRUE.)
         call ReplaceVwithW(v,w)
         call NORMALIZE(v)
      enddo

      end subroutine InverseRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RayleighQuotientItn(v,H,npow,nals,Eshift,which)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves (H-E1)v=w for vector v

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow,nals,which
      real(kind=8), intent(in) :: Eshift
      real(kind=8) :: l,m
      integer :: i

      m=Eshift

      do i=1,npow
!         write(*,'(A,I0,A,I0,A,f15.8)') 'ipow = (',i,'/',npow,&
!               ') mu = ',m
         w=CopyCP(v)
         call LinSolver_alg(H,w,v,nals,1,m,which,.FALSE.)
         l=PRODVV(v,w)
         m=m+1.d0/l
         call ReplaceVwithW(v,w)
         call NORMALIZE(v)
      enddo

      end subroutine RayleighQuotientItn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RayleighQuotientIntertwining(v,H,npow,nals,Eshift,which)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves (H-E1)v=w for vector v

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      integer, intent(in) :: npow,nals,which
      real(kind=8), intent(in) :: Eshift

      if (which.eq.1) then
         call LintertwinedInvItn_1(H,v,nals,1,Eshift,.TRUE.)
      elseif (which.eq.2) then
         call LintertwinedInvItn_2(H,v,nals,1,Eshift,.TRUE.)
      else
         write(*,*) 'Must choose {1,2} for RQ intertwining, not ',which
         call &
         AbortWithError('RayleighQuotientIntertwining(): bad choice')
      endif

      end subroutine RayleighQuotientIntertwining

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FoldedRecurse(v,H,npow,Eshift,bounds)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies Hamiltonian (H-E)^npow*v to vector v
! The ishift parameter controls the shifting

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow
      real(kind=8), intent(in) :: Eshift
      real(kind=8), intent(in) :: bounds(2)
      real(kind=8) :: sh(2),s1,s2,sig
      integer :: i

      s1=(Eshift-bounds(1))**2
      s2=(bounds(2)-Eshift)**2
      sig=sqrt(0.5*(s1+s2))
      sh(1)=Eshift+sig
      sh(2)=Eshift-sig

!      write(*,*) 'Etarget = ',Eshift,'; sh(1) = ',sh(1),'; sh(2) = ',sh(2)

      do i=1,npow
!        w <- (H-sigma(1)*I)*v; then v <- w
         call CPMM(H,1,sh(1),.FALSE.,v,0,0.d0,.FALSE.,w)
         call reduc(v,w)
         call FlushCP(w)
!         call ReplaceVwithW(v,w)
!        w <- (H-sigma(2)*I)*v; then v <- w
         call CPMM(H,1,sh(2),.FALSE.,v,0,0.d0,.FALSE.,w)
         call reduc(v,w)
         call FlushCP(w)
!        Normalize v only after matrix-vector product pair
         call NORMALIZE(v)
      enddo

      end subroutine FoldedRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PowrRecurse(v,H,npow,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies Hamiltonian (H-E)^npow*v to vector v
! The ishift parameter controls the shifting

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow
      real*8, intent(in)  :: Eshift
      integer :: i

      do i=1,npow
!        w <-- H*v; then v <-- w
         call CPMM(H,1,Eshift,.FALSE.,v,0,0.d0,.FALSE.,w)

!        Reduce and normalize v
         call reduc(v,w)
         call FlushCP(w)
         call NORMALIZE(v)
      enddo

      end subroutine PowrRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinfinityNorm(F,npow,nals,rown,coln,norm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compute largest element of F via updated power iteration. Returns
! value of largest element and row/col indices

      implicit none
      TYPE (CP), intent(in) :: F
      TYPE (CP) :: Ft,Fr,G
      integer, allocatable, intent(out) :: rown(:),coln(:)
      logical, parameter :: show=.TRUE.
      real*8, parameter :: tol=1.d-4
      real*8, intent(out) :: norm
      integer, intent(in) :: npow,nals
      integer :: i,d,ndof
      real*8 :: val,conv

      ndof=SIZE(F%nbas)

      Ft=CopyCP(F)
      Fr=RandomCP(F,1)
      call NORMALIZE(Ft)

      DO i=1,npow
!        Power iteration Ft*Ft=G, reduction Ft<-ALS-G, normalization
         call CPMM_vec(Ft,Ft,G)
         call reduc_ALS(G,Ft,nals)
         call NORMCOEF(Ft)

!        Reduce Ft to a normalized rank-1 CP and get dominant index
         call reduc_SR1(Ft,Fr,nals)
         Fr%coef(1)=1.d0  ! normalize coef of Fr
         call GetRank1DominantEntry(Fr,rown,coln,val)

!        Check convergence
         conv=abs(abs(val)-1.d0)
         IF (show) &
            write(*,'(X,A,ES14.8)') 'Distance from L_infty = ',conv
         IF (conv.lt.tol) EXIT
         deallocate(rown,coln)
      ENDDO
      call FlushCP(G)
      call FlushCP(Ft)
      call FlushCP(Fr)

!     Extract the norm from the original vector using the indices
      norm=ExtractCPmatrixElement(F,rown,coln)

      end subroutine LinfinityNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      END MODULE BLOCKPOWER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
