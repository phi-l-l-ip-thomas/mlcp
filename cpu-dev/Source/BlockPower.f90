!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKPOWER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TISE solver using the block power method

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE CPMMM
      USE MODVECVEC
      USE SEPDREPN
      USE REDUCTION
!!!
      USE LINSOLVER
!!!
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

!      write(*,'(/X,A,f16.8/)') 'Eshift = ',Eshift

      end subroutine GetBlockShift

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InverseRecurse(v,H,npow,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves (H-E1)v=w for vector v

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow
      real*8, intent(in)  :: Eshift
      integer :: i

      w=CopyCP(v)
      do i=1,npow
!         w=RandomCP(v)
         call LinSolver_alg(H,w,v,3,0,Eshift,1)
         call FlushCP(v)
         v=CopyCP(w)
         call NORMALIZE(v)
!         call FlushCP(w)
      enddo
!      call FlushCP(w)

      end subroutine InverseRecurse

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

!!! TEST
         if (i.gt.npow) then !!! .eq.npow to do test
            call NORMALIZE(w)
!            call CPifylittleg(w,(/SIZE(v%coef),SIZE(H%coef)+1/))
         endif

         if (i.eq.0) then !!! .eq.1 to do test
            write(*,*) 'ipow :',i
            write(*,*)
            call TestPenalty(v,w)
            write(*,*)
         endif
!!! END TEST

!        Reduce and normalize v
         call reduc(v,w)
         call FlushCP(w)
         call NORMALIZE(v)
      enddo

      end subroutine PowrRecurse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TestPenalty(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test various regularization values

      implicit none
      TYPE (CP), INTENT(IN) :: v,w
      TYPE (CP)  :: F
      TYPE (ALS) :: X
      character(len=64) :: nm
      integer :: i,nals,conv
      real*8  :: pen

      nals=100
      pen=1.d-15

      do i=1,15
         F=CopyCP(v)
         call NewALS(X,v,w)

         write(nm,'(A,ES8.1)') 'w -> v redn, pen = ',pen
         call X%setname(nm)
         call X%setvalue('penalty',pen)
         conv=ALS_reduce(X,F,w,nals,nm)

         pen=10.d0*pen
         call FlushCP(F)
         call X%flush()
!         write(*,*)
      enddo

      end subroutine TestPenalty

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PowrRecurseA(v,H,npow,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies Hamiltonian (H-E)^npow*v to vector v
! The ishift parameter controls the shifting
! No rank-reduction until the end !!!

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: w
      integer, intent(in) :: npow
      real*8, intent(in)  :: Eshift
      integer, allocatable :: strides(:)
      integer :: i,rkv,rkH

      rkH=SIZE(H%coef)+1
      rkv=SIZE(v%coef)

      ALLOCATE(strides(npow+1))
      strides(:)=rkH
      strides(1)=rkv

      do i=1,npow
!        w <-- H*v; then v <-- w
         call CPMM(H,1,Eshift,.FALSE.,v,0,0.d0,.FALSE.,w)
         call NORMALIZE(w)
         call ReplaceVwithW(v,w)
      enddo
!!! TEST
!      call CPifylittleg(v,strides)
!!! END TEST

      DEALLOCATE(strides)

!     Reduce and normalize v
      w=RandomCP(v,rkv)
      call reduc(w,v)
      call ReplaceVwithW(v,w)
      call NORMALIZE(v)

      end subroutine PowrRecurseA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPifylittleg(G,rnx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Check the rank of each row of R vals in G

      implicit none
      TYPE (CP), INTENT(IN) :: G
      TYPE (CP), ALLOCATABLE :: B(:)
      TYPE (CP) :: F
      integer, intent(in) :: rnx(:)
      integer, allocatable :: cols(:)
      real*8, allocatable :: rfacs(:,:)
      integer :: d,i,j,ndof,bsize,srank,conv
      character(len=64) :: nm

      srank=10
      ndof=SIZE(rnx)
      bsize=SIZE(G%base,1)+1
      ALLOCATE(B(bsize))

      write(*,*)
      write(*,*) 'CPifylittleg...'
      write(*,*)

      if (ndof.lt.3) then
         write(*,*) 'Here is G:'
         call G%printvec()
      else
         write(*,*) 'Rank of G = ',G%R()
         call PrintCPmat_gen(G,.FALSE.)
      endif

      ALLOCATE(rfacs(G%R(),1),cols(SIZE(rnx)))
      cols(:)=1

!     Convert coefs of G into a CP-vec
      rfacs(:,1)=G%coef(:)
      F=Matrix2CP(rfacs,rnx,cols)

      if (ndof.lt.3) then
         write(*,*) 'Coefs unCPified'
         call PrintMatrix(rfacs)
         write(*,*) 'Coefs CP-ified'
         write(*,*)
         call F%print()
         write(*,*) 'RANK (spacer)'
      else
         write(nm,'(A)') 'Coefs'

!        Reduce rank of coefs
!        Coefs go in last entry of block
         B(bsize)=RandomCP(F,srank)
         conv=ALS_reduce(B(bsize),F,50,nm)
      endif

      call FlushCP(F)

      DO d=1,G%D()
         DO i=1,G%rows(d)
            j=G%ibas(d)-1+i

!           Convert this row of base of G into a CP-vec
            rfacs(:,1)=G%base(j,:)
            F=Matrix2CP(rfacs,rnx,cols)

            if (ndof.lt.3) then
               write(*,*) 'Base unCPified'
               call PrintMatrix(rfacs)
               write(*,*) 'Base CP-ified, d = ',d,'; i = ',i
               write(*,*)
               call F%print()
               write(*,*) 'RANK (spacer)'
            else
               write(nm,'(2(A,I0))') 'Base, d=',d,'; i=',i

!              Reduce rank of this row of base
!              Coefs go in j-th entry of block
               B(j)=RandomCP(F,srank)
               conv=ALS_reduce(B(j),F,50,nm)

            endif
            call FlushCP(F)
!            write(*,*)
         ENDDO
      ENDDO

      DEALLOCATE(rfacs)

!     Now reconstitute G from CP-vecs in B
      F=NewCP(G,G%R())
      rfacs=CP2Matrix(B(bsize))
      F%coef(:)=rfacs(:,1)
      DEALLOCATE(rfacs)
      DO d=1,G%D()
         DO i=1,G%rows(d)
            j=G%ibas(d)-1+i
            rfacs=CP2Matrix(B(j))
            F%base(j,:)=rfacs(:,1)
            DEALLOCATE(rfacs)
         ENDDO
      ENDDO

      write(*,*)
      write(*,*) 'CPR ||F-G|| = ',calc_FmG(G,F)
      write(*,*)

      call FlushCP(F)
      DEALLOCATE(B)

!      call AbortWithError('done with little-g grilch')

      end subroutine CPifylittleg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PowSpectralRange(npow,ncyc,Q,H,bounds)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Estimates the spectral range of the Hamiltonian using power method

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP) :: Q(:)
      TYPE (CP), allocatable :: Qt(:)
      integer, intent(in) :: npow,ncyc
      real*8, intent(out) :: bounds(2)
      integer :: i,j,gst,ndof
      real*8  :: btmp(2)

      ndof=SIZE(Q(1)%nbas)
      bounds=0.d0

      write(*,'(X,A/)') 'Power spectral range estimation'
      write(*,'(A,X,2(9X,A))') '  i','boundsl','boundsu'

      ALLOCATE(Qt(2))

!     Guesses for upper and lower eigenvectors
      Qt(1)=ZeroCPvec(Q(1)%nbas)
      Qt(2)=ZeroCPvec(Q(1)%nbas)
      Qt(1)%coef(1)=1.d0
      Qt(2)%coef(1)=1.d0
      gst=0
      DO i=1,ndof
         Qt(1)%base(gst+1,1)=1.d0
         gst=gst+Qt(1)%nbas(i)
         Qt(2)%base(gst,1)=1.d0
      ENDDO

!     Initial bounds
!$omp parallel
!$omp do private(j) schedule(static)
      DO j=1,2
         bounds(j)=RayleighQuotient2(Qt(j),H)
      ENDDO
!$omp end do
!$omp end parallel

!     btmp is the shift for the power method
      btmp(1)=0.5*(bounds(1)+bounds(2))
      btmp(2)=0.d0
      write(*,'(i3,2(X,2f16.6))') 0,bounds(1),bounds(2)

!     Run the power method to improve the bounds
      DO i=1,ncyc
!$omp parallel
!$omp do private(j) schedule(static)
         DO j=1,2
            call PowrRecurse(Qt(j),H,npow,btmp(j))
            bounds(j)=RayleighQuotient2(Qt(j),H)
         ENDDO
!$omp end do
!$omp end parallel
!        Update the shift for the ground state
         btmp(1)=0.5*(bounds(1)+bounds(2))
         write(*,'(i3,2(X,2f16.6))') i,bounds(1),bounds(2)
      ENDDO

      DEALLOCATE(Qt)

      write(*,'(/X,A,2(f15.6,A))') 'Spectral range of H = [',&
                                   bounds(1),',',bounds(2),']'
      btmp(1)=0.001*(bounds(2)-bounds(1))
      bounds(1)=bounds(1)-btmp(1)
      bounds(2)=bounds(2)+btmp(1)
      write(*,'(X,A,2(f15.6,A)/)') 'Range padded by .1% = [',&
                                   bounds(1),',',bounds(2),']'

      end subroutine PowSpectralRange

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
