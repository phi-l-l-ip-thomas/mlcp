!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines which operate on CP-blocks

      USE ERRORTRAP
      USE UTILS
      USE MODHVEC
      USE MODVECVEC
      USE SEPDREPN
      USE HAMILSETUP
      USE REDUCTION
      USE ALSPOW
      USE ALSUTILS

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GRAMORTHO(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gram-Schmidt orthogonalizes a block of CP-format vectors

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      TYPE (CPvec)  :: v
      real*8, allocatable :: vivj(:)
      integer :: i,j,nbloc

      nbloc=size(Q)
      ALLOCATE(vivj(nbloc))

!     Normalize the first vector
      call NORMCOEF(Q(1))

!     Loop through the remaining vectors
!     for each, project out all preceding vectors
      do i=2,nbloc

!        Normalize i-th vector
         call NORMCOEF(Q(i))

!        Get weights
         vivj(i)=1.d0

!$omp parallel
!$omp do private(j) schedule(static)
         do j=1,i-1
!           Q(i) <- Q(i) - <Q(i),Q(j)>*Q(j)
            vivj(j)=-PRODVV(Q(i),Q(j))
         enddo
!$omp end do
!$omp end parallel

!        Compute the orthogonalized vector
         call SUMLCVEC(v,Q(1:i),vivj(1:i))

!        Reduce v into Q(i), normalize
         call reduc(Q(i),v)
         call NORMCOEF(Q(i))
      enddo

      DEALLOCATE(vivj)

      end subroutine GRAMORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine pGRAMORTHO(Q,k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Orthogonalizes a block of CP-format vectors via Gram-Schmidt, assuming
! the first k vectors are already orthonormal

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      TYPE (CPvec)  :: v
      integer, intent(in) :: k
      real*8, allocatable :: vivj(:)
      integer :: i,j,nbloc

      nbloc=size(Q)

      IF (k.eq.nbloc) RETURN
      IF (k.lt.0 .or. k.gt.nbloc) &
         call AbortWithError('pGRAMORTHO(): k is out of range')

      ALLOCATE(vivj(nbloc))

!     Loop through the remaining vectors
!     for each, project out all preceding vectors
      do i=k+1,nbloc

!        Normalize i-th vector
         call NORMCOEF(Q(i))

!        Get weights
         vivj(i)=1.d0

!$omp parallel
!$omp do private(j) schedule(static)
         do j=1,i-1
!           Q(i) <- Q(i) - <Q(i),Q(j)>*Q(j)
            vivj(j)=-PRODVV(Q(i),Q(j))
         enddo
!$omp end do
!$omp end parallel

!        Compute the orthogonalized vector
         call SUMLCVEC(v,Q(1:i),vivj(1:i))

!        Reduce v into Q(i), normalize
         call reduc(Q(i),v)
         call NORMCOEF(Q(i))
         call FlushCPvec(v)
      enddo

      DEALLOCATE(vivj)

      end subroutine pGRAMORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetEigenFxn(v,Q,QHQ,k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given the eigenvectors (QHQ), the CP-format vectors (in Q), get the
! k-th eigenfunction
! No reduction or normalization is done here.

      implicit none
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      TYPE (CPvec), INTENT(OUT) :: v
      real*8, intent(in)  :: QHQ(:,:)
      real*8, allocatable :: weights(:)
      integer, intent(in) :: k
      integer :: n

      n=SIZE(QHQ,1)
      ALLOCATE(weights(n))
      weights=0.d0
      weights(k)=1.d0
      call GetLCEigenFxn(v,Q,QHQ,weights)
      DEALLOCATE(weights)

      end subroutine GetEigenFxn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetLCEigenFxn(v,Q,QHQ,weights)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given the eigenvectors (QHQ), the CP-format vectors (in Q), build a
! linear combination of eigenfunctions weighed by 'weights'. 
! No reduction or normalization is done here.

      implicit none
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      TYPE (CPvec), INTENT(OUT) :: v
      real*8, intent(in)  :: QHQ(:,:)
      real*8, intent(in)  :: weights(:)
      real*8, allocatable :: facs(:)
      integer :: i,j,n

      n=SIZE(QHQ,1)
      ALLOCATE(facs(n))
      facs=0.d0
      DO j=1,n
         DO i=1,n
            facs(j)=facs(j)+weights(i)*QHQ(j,i)
         ENDDO
      ENDDO

      call SUMLCVEC(v,Q,facs)
      DEALLOCATE(facs)

      end subroutine GetLCEigenFxn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Diagonalize(Q,H,avec,intw,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! One-step function to solve the generalized eigenvalue problem and 
! update vectors. Setting intw=.T. uses intertwining for the update

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      integer, intent(in)   :: nitn,lm
      logical, intent(in)   :: intw
      real*8, intent(inout) :: avec(:)
      real*8, allocatable :: QHQ(:,:),S(:,:)
      integer :: nbloc,i

!     Set parameters
      nbloc=SIZE(Q)

      allocate(QHQ(nbloc,nbloc),S(nbloc,nbloc))

!     Calculate QHQ and S matrices
      IF (intw) THEN
          call GetQHQ3(Q,H,QHQ,nitn,lm)  ! Reduces H*Q, then calcs Q^THQ
      ELSE
         call GetQHQ(Q,H,QHQ)
      ENDIF
      call GetOverlaps(Q,S)

!     Diagonalize QHQ, accounting for overlaps
      call SolveGenEigval(avec,S,QHQ,'V')

!     Update block vectors after diagonalization
!     q^n_{new} <- sum_{i=1} ^ m U_{im} q^i_{old}
      IF (intw) THEN
         call UpdateVecs2(Q,QHQ,nitn,lm) ! Avoids long vectors
      ELSE
         call UpdateVecs(Q,QHQ)
      ENDIF

      deallocate(QHQ,S)

      end subroutine Diagonalize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ(Q,H,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      TYPE (CPvec), ALLOCATABLE :: HQ(:)
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,nbloc

      nbloc=SIZE(Q)
      ALLOCATE(HQ(nbloc))
      QHQ=0.d0

!$omp parallel
!$omp do private(i,j)
      do i=1,nbloc
!        HQ=H*Q(i)
         call PRODHV(Q(i),HQ(i),H,0,0.d0)
!        Reduce the rank of HQ as this greatly accelerates computing the
!        inner products below
         call reduc(HQ(i))
         do j=i,nbloc
!           QHQ(i,j)=<Q(j),H(Q(i))>
            QHQ(i,j)=PRODVV(Q(j),HQ(i))
         enddo
         call FlushCPvec(HQ(i))
      enddo
!$omp enddo
!$omp end parallel

      DEALLOCATE(HQ)

      end subroutine GetQHQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ2(Q,H,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. Only one term of H*Q is used
! at a time to save memory

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: Q(:)
      TYPE (CPvec) :: HQ
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,k,l,rF,rH,nbloc

      nbloc=SIZE(Q)
      rH=SIZE(H%opcoef,1)

      QHQ=0.d0

      do i=1,nbloc
!        HQ=H*Q(i), term by term
         rF=SIZE(Q(i)%coef)
         do l=1,rF
            do k=1,rH
               call PRODHVR1(Q(i),HQ,H,l,k)
!$omp parallel
!$omp do private(j)
               do j=i,nbloc
!                 QHQ(i,j)=<Q(j),H(Q(i))>
                  QHQ(i,j)=QHQ(i,j)+PRODVV(Q(j),HQ)
               enddo
!$omp enddo
!$omp end parallel
               call FlushCPvec(HQ)
            enddo
         enddo
      enddo

      end subroutine GetQHQ2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ3(Q,H,QHQ,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. The matrix-vector product H*Q
! is computed and reduced via intertwining to save memory and time

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      TYPE (CPvec), ALLOCATABLE :: HQ(:)
      integer, intent(in)   :: nitn,lm
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,nbloc

      nbloc=SIZE(Q)
      ALLOCATE(HQ(nbloc))
      QHQ=0.d0

!$omp parallel
!$omp do private(i,j)
      do i=1,nbloc
!        HQ=H*Q(i)
         call PRODHV_ALS_alg(Q(i),HQ(i),H,0,0.d0,nitn,lm)
         do j=i,nbloc
!           QHQ(i,j)=<Q(j),H(Q(i))>
            QHQ(i,j)=PRODVV(Q(j),HQ(i))
         enddo
         call FlushCPvec(HQ(i))
      enddo
!$omp enddo
!$omp end parallel

      DEALLOCATE(HQ)

      end subroutine GetQHQ3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetOverlaps(Q,S)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T Q for a block of vectors

      implicit none
      TYPE (CPvec), INTENT(IN) :: Q(:)
      real*8, intent(inout) :: S(:,:)
      integer, allocatable  :: ijc(:,:)
      integer :: i,j,ij,nij,nbloc

      nbloc=SIZE(Q)
      nij=nbloc*(nbloc+1)/2

      ALLOCATE(ijc(nij,2))
      ij=0
      DO i=1,nbloc
         DO j=i,nbloc
            ij=ij+1
            ijc(ij,:)=(/i,j/)
         ENDDO
      ENDDO

      S=0.d0
!$omp parallel
!$omp do private(ij)
      DO ij=1,nij
         S(ijc(ij,1),ijc(ij,2))=PRODVV(Q(ijc(ij,2)),Q(ijc(ij,1)))
      ENDDO
!$omp enddo
!$omp end parallel

      DEALLOCATE(ijc)

      end subroutine GetOverlaps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQdiag(Q,H,QHQd,nconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes diagonal elements of Q^T H Q for a block of vectors
! The first nconv vectors are skipped

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: Q(:)
      integer, intent(in)   :: nconv
      real*8, intent(inout) :: QHQd(:)
      integer :: i,nbloc

      nbloc=SIZE(Q)

!$omp parallel
!$omp do private(i)
      do i=nconv+1,nbloc
         QHQd(i)=RayleighQuotient2(Q(i),H)   
      enddo
!$omp enddo
!$omp end parallel

      end subroutine GetQHQdiag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetOverlapQ1Q2(Q1,Q2,S)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T Q for a block of vectors

      implicit none
      TYPE (CPvec), INTENT(IN) :: Q1(:),Q2(:)
      real*8, intent(inout) :: S(:,:)
      integer :: i,j,n1,n2

      n1=SIZE(Q1)
      n2=SIZE(Q2)

      S=0.d0
      do i=1,n1
         do j=1,n2
!           S(i,j)=<Q1(i),Q2(j)>
            S(i,j)=PRODVV(Q1(i),Q2(j))
         enddo
      enddo

      end subroutine GetOverlapQ1Q2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateVecs(Q,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces a block of vectors Q with the eigenvectors whose coefficients
! are contained in QHQ

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      TYPE (CPvec), ALLOCATABLE   :: Qold(:)
      real*8, intent(in) :: QHQ(:,:)
      integer :: i,nbloc

      nbloc=SIZE(Q)

!     Copy Q to Qold
      ALLOCATE(Qold(nbloc))

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call ReplaceVwithW(Qold(i),Q(i))
      ENDDO
!$omp enddo
!$omp end parallel

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call GetEigenFxn(Q(i),Qold,QHQ,i)
!        reduce and normalize Q(i)
         call reduc(Q(i))
         call NORMCOEF(Q(i))
      ENDDO
!$omp enddo
!$omp end parallel

      DEALLOCATE(Qold)

      end subroutine UpdateVecs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateVecs2(Q,QHQ,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces a block of vectors Q with the eigenvectors whose coefficients
! are stored in QHQ. This version uses intertwining to reduce the rank.

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      TYPE (CPvec), ALLOCATABLE   :: Qold(:)
      real*8, intent(in)  :: QHQ(:,:)
      integer, intent(in) :: nitn,lm
      integer :: i,nbloc

      nbloc=SIZE(Q)

!     Copy Q to Qold
      ALLOCATE(Qold(nbloc))

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call CopyWtoV(Qold(i),Q(i))
      ENDDO
!$omp enddo
!$omp end parallel

!     Build and reduce the eigenfunction
!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call ALS_SUMLCVEC_alg(Q(i),Qold,QHQ(:,i),nitn,lm)
      ENDDO
!$omp enddo
!$omp end parallel

      DEALLOCATE(Qold)

      end subroutine UpdateVecs2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortVecs(Q,eigv,nconv,Eref)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts a block of vectors Q in terms of increasing distance from Eref
! The first nconv vectors are assumed to already be in order

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      TYPE (CPvec), ALLOCATABLE   :: Qold(:)
      integer, intent(inout) :: nconv
      real*8, intent(inout)  :: eigv(:)
      real*8, intent(in)     :: Eref
      real*8, allocatable    :: ix(:),edif(:)
      integer :: i,nbloc

      nbloc=SIZE(Q)
      IF (size(eigv).ne.nbloc) &
         call AbortWithError('SortVecs() size(eigv) != nbloc')

!     List of ordered indices
      allocate(ix(nbloc),edif(nbloc))
      do i=1,nbloc
         ix(i)=i
         edif(i)=abs(eigv(i)-Eref)
      enddo

!     Sort the eigenvalues
      call dsort(edif,ix,nbloc,2)
      edif(:)=eigv(:)
      DO i=1,nbloc
         IF (i.ne.int(ix(i))) eigv(i)=edif(int(ix(i)))
      ENDDO

!     Copy the eigenvectors into the temp array
      ALLOCATE(Qold(nbloc))
!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         IF (i.ne.int(ix(i))) call ReplaceVwithW(Qold(i),Q(i))
      ENDDO
!$omp enddo
!$omp end parallel

!     Eigenvector resort
!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         IF (i.ne.int(ix(i))) call ReplaceVwithW(Q(i),Qold(int(ix(i))))
      ENDDO
!$omp enddo
!$omp end parallel

!     If a 'converged' vector swaps places with another vector
!     then it is not really converged, so change nconv accordingly!
      DO i=1,nconv
         IF (i.ne.int(ix(i))) THEN
            nconv=i-1
            EXIT
         ENDIF
      ENDDO

      DEALLOCATE(Qold,ix,edif)

      end subroutine SortVecs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AugmentQWithRandom(Q,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments each vector in a block to rank nrk by adding random terms.

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      integer, intent(in) :: nrk
      integer :: i,nbloc

      nbloc=SIZE(Q)

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call AugmentVWithRandom(Q(i),nrk)
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine AugmentQWithRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NoisifyQ(Q,mag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments each vector in a block to rank nrk by adding random terms.

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      real*8, intent(in) :: mag
      integer :: i,nbloc

      nbloc=SIZE(Q)

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         call Noisify(Q(i),mag)
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine NoisifyQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RayleighQuotient(v,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient for vector v and Hamiltonian H

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: v
      TYPE (CPvec) :: w
      real*8 :: RayleighQuotient

!     w = H*v       
      call PRODHV(v,w,H,0,0.d0)

!     Rayleigh quotient = <v|w>/<v|v>
      RayleighQuotient=PRODVV(v,w)/PRODVV(v)
!
      call FlushCPvec(w)

      end function RayleighQuotient

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RayleighQuotient2(v,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient for vector v and Hamiltonian H

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: v
      TYPE (CPvec) :: w
      real*8  :: RayleighQuotient2
      integer :: i,k,rF,rH

      rH=SIZE(H%opcoef,1)
      rF=SIZE(v%coef)

      RayleighQuotient2=0.d0

!     Calc w=H*v, term by term
      do i=1,rF
         do k=1,rH
            call PRODHVR1(v,w,H,i,k)
!           Build <w|v>
            RayleighQuotient2=RayleighQuotient2+PRODVV(v,w)
            call FlushCPvec(w)
         enddo
      enddo

!     Rayleigh quotient = <v|w>/<v|v>
      RayleighQuotient2=RayleighQuotient2/PRODVV(v)

      end function RayleighQuotient2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RayleighResidual(v,H,eig)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient residual for vector v and Hamiltonian H
! v should be normalized using NORMCOEF() before calling

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: v
      TYPE (CPvec) :: w
      real*8, intent(in) :: eig
      real*8 :: RayleighResidual

!     w = H*v - eig*v
      call PRODHV(v,w,H,0,0.d0)
      call SUMVECVEC(w,1.d0,v,-eig)

!     If the rank of w is large, then the call to PRODVV below will be
!     expensive. If lower accuracy is tolerable, one can reduce and then
!     call PRODVV with a smaller vector
!      call reduc(w)

!     Calculate ||w||
      RayleighResidual=sqrt(abs(PRODVV(w)))

      call FlushCPvec(w)

      end function RayleighResidual

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
