!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines which operate on CP-blocks

      USE ERRORTRAP
      USE UTILS
      USE MUNKRES
      USE CPMMM
      USE MODVECVEC
      USE SEPDREPN
      USE REDUCTION
      USE ALSPOW
      USE ALSUTILS
      USE ALSDRVR

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GRAMORTHO(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gram-Schmidt orthogonalizes a block of CP-format vectors

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP)  :: v
      real*8, allocatable :: vivj(:)
      integer :: i,j,nbloc

      nbloc=size(Q)
      ALLOCATE(vivj(nbloc))

!     Normalize the first vector
      call NORMALIZE(Q(1))

!     Loop through the remaining vectors
!     for each, project out all preceding vectors
      do i=2,nbloc

!        Normalize i-th vector
         call NORMALIZE(Q(i))

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
         call NORMALIZE(Q(i))
      enddo

      DEALLOCATE(vivj)

      end subroutine GRAMORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine pGRAMORTHO(Q,k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Orthogonalizes a block of CP-format vectors via Gram-Schmidt, assuming
! the first k vectors are already orthonormal

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP)  :: v
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
         call NORMALIZE(Q(i))

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
         call NORMALIZE(Q(i))
         call FlushCP(v)
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
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), INTENT(OUT) :: v
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
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), INTENT(OUT) :: v
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), INTENT(IN) :: H
      integer, intent(in)   :: nitn,lm
      logical, intent(in)   :: intw
      real*8, intent(inout) :: avec(:)
      real*8, allocatable :: QHQ(:,:),S(:,:)
      integer :: nbloc

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

!!!   TEST
!      call TestAssignQHQ(H,QHQ,avec,Q(1)%rows,Q(1)%rows)
!!!

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

      subroutine TestAssignQHQ(H,M,eig,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assigns QHQ by applying Munkres algorithm on squared eigenvectors

      implicit none
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP) :: U,F
      integer, intent(in) :: rows(:),cols(:)
      real*8, intent(in)  :: M(:,:),eig(:)
      real*8, allocatable :: QHQ(:,:),QHQ2(:,:),AGN(:,:),ET(:,:)
      integer :: nbloc,i,j,npars,ngrp

      nbloc=SIZE(M,1)

      ALLOCATE(QHQ(nbloc,nbloc),QHQ2(nbloc,nbloc))
      QHQ(:,:)=M(:,:)
      QHQ2(:,:)=0.d0
      do i=1,nbloc
         do j=1,nbloc
            QHQ2(i,j)=QHQ(i,j)**2
         enddo
      enddo

      write(*,*)
      write(*,*) 'Eigenvector magns:'
      npars=10
      ngrp=(nbloc+npars-1)/npars
      do i=1,ngrp
         call PrintMatrix(QHQ2(:,((i-1)*npars+1):min(i*npars,nbloc)))
      enddo

      write(*,*)
      write(*,*) 'Assignment matrix:'
      AGN=AssignMatrix(QHQ2)
      do i=1,ngrp
         call PrintMatrix(AGN(:,((i-1)*npars+1):min(i*npars,nbloc)))
      enddo

      write(*,*)
      write(*,*) 'Assigned eigenvector magns:'
      call MatrixMult(QHQ2,.FALSE.,AGN,.TRUE.)
      do i=1,ngrp
         call PrintMatrix(QHQ2(:,((i-1)*npars+1):min(i*npars,nbloc)))
      enddo
      DEALLOCATE(QHQ2)

      write(*,*)
      write(*,*) 'Assigned eigenvector coefs:'
      call MatrixMult(QHQ,.FALSE.,AGN,.TRUE.)
!     Sign align the eigenvectors
      DO i=1,nbloc
         IF (QHQ(i,i).lt.0.d0) QHQ(:,i)=-QHQ(:,i)
      ENDDO
      do i=1,ngrp
         call PrintMatrix(QHQ(:,((i-1)*npars+1):min(i*npars,nbloc)))
      enddo

      write(*,*)
      write(*,*) 'Assigned eigenvalue matrix:'
      ALLOCATE(ET(1,nbloc))
      ET(1,:)=eig(:)
      call MatrixMult(ET,.FALSE.,AGN,.TRUE.)
      do i=1,ngrp
         call PrintMatrix(ET(:,((i-1)*npars+1):min(i*npars,nbloc)))
      enddo
      ALLOCATE(QHQ2(nbloc,nbloc))
      QHQ2(:,:)=0.d0
      DO i=1,nbloc
         QHQ2(i,i)=ET(1,i)
      ENDDO

      write(*,*)
      write(*,*) 'CP-ified eigenvectors:'
      U=Matrix2CP(QHQ,rows,cols)
      call U%print()

      call ShowVecRanks(H,U,30)
!      call TestRankSuccessive(U)
      call FlushCP(U)

      write(*,*)
      write(*,*) 'CP-ified eigenvalues:'
      F=Matrix2CP(QHQ2,rows,cols)
      call F%print()
      call TestRankSuccessive(F)
      call FlushCP(F)

      DEALLOCATE(QHQ,QHQ2,AGN,ET)

      end subroutine TestAssignQHQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TestRankSuccessive(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduce rank of F from 1 to rkF-1, and check error

      implicit none
      TYPE (CP), INTENT(IN) :: F
      TYPE (CP) :: Fr
      integer :: i,redstat
      character(len=64) :: tag

      DO i=1,F%R()-1
         write(tag,'(A,I0)') 'Test rank = ',i
!         Fr=IdentityCPMatrix(F%rows,F%cols,F%sym)
!         call AugmentVWithRandom(Fr,i)
         Fr=RandomCP(F,i)
         redstat=ALS_reduce(Fr,F,100,tag)
         call FlushCP(Fr)
!         write(*,*)
         IF (redstat.eq.1) EXIT ! Convergence reached
      ENDDO

      end subroutine TestRankSuccessive

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowVecRanks(H,U,maxrank)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints normalization values for all the vectors in U

      implicit none
      TYPE (CP), INTENT(IN) :: H,U
      TYPE (CP) :: F,G,Z
      integer, intent(in)  :: maxrank
      integer, allocatable :: indx(:)
      integer :: i,j,k,l,ndof,redstat
      real*8  :: norm,rq
      character(len=64) :: frmt,tag

      ndof=SIZE(U%cols)

      allocate(indx(ndof))
      indx(:)=U%cols(:)

      write(*,*)
      DO i=1,PRODUCT(U%cols)
         call NextIndex(indx,U%cols)
         G=ExtractCPvec(U,indx,.TRUE.)
         norm=sqrt(abs(PRODVV(G)))
         call NORMALIZE(G)
         rq=RayleighQuotient(G,H)
         write(frmt,'(A,I0,A)') '(A,I8,A,',ndof,'(I2,X),2(A,f16.8))'
         write(*,frmt) "Vec: ",i," (",(indx(k),k=1,ndof),&
         "); RQ = ",rq,'; ||G|| = ',norm
!         F=ALS_reduce_adaptive(mxrk,G,1.d-7,100,nm)
         DO j=1,maxrank
            write(frmt,'(A,I0,A)') '(A,I8,A,',ndof,'(I2,X),A,I0)'
            write(tag,frmt) "Vec: ",i," (",(indx(k),k=1,ndof),&
            "); rk = ",j
            F=RandomCP(G,j)
            redstat=ALS_reduce(F,G,100,tag)
            IF (redstat.eq.1) THEN
               DO l=1,SIZE(G%coef)
                  Z=NewCP(G,1)
                  call GenCopyWtoV(Z,G,1,1,l,l)
                  write(*,*) 'Term: ',l,'<F,G(l)> = ',PRODVV(F,Z)
               ENDDO
               call FlushCP(F)
               EXIT ! Convergence reached
            ENDIF
            call FlushCP(F)
         ENDDO
         write(*,*)
         call FlushCP(G)
      ENDDO
      deallocate(indx)

      end subroutine ShowVecRanks

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ(Q,H,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), ALLOCATABLE :: HQ(:)
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,nbloc

      nbloc=SIZE(Q)
      ALLOCATE(HQ(nbloc))
      QHQ=0.d0

!$omp parallel
!$omp do private(i,j)
      do i=1,nbloc
!        HQ=H*Q(i)
         call CPMM(H,.FALSE.,Q(i),.FALSE.,HQ(i))
!        Reduce the rank of HQ as this greatly accelerates computing the
!        inner products below
         call reduc(HQ(i))
         do j=i,nbloc
!           QHQ(i,j)=<Q(j),H(Q(i))>
            QHQ(i,j)=PRODVV(Q(j),HQ(i))
         enddo
         call FlushCP(HQ(i))
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
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (CP) :: HQ
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,k,l,rF,rH,nbloc

      nbloc=SIZE(Q)
      rH=SIZE(H%coef)

      QHQ=0.d0

      do i=1,nbloc
!        HQ=H*Q(i), term by term
         rF=SIZE(Q(i)%coef)
         do l=1,rF
            do k=1,rH
               call CPMM(H,k,k,.FALSE.,Q(i),l,l,.FALSE.,HQ)
!$omp parallel
!$omp do private(j)
               do j=i,nbloc
!                 QHQ(i,j)=<Q(j),H(Q(i))>
                  QHQ(i,j)=QHQ(i,j)+PRODVV(Q(j),HQ)
               enddo
!$omp enddo
!$omp end parallel
               call FlushCP(HQ)
            enddo
         enddo
      enddo

      end subroutine GetQHQ2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ2a(Q,H,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. Only one term of H is used
! at a time to save memory

      implicit none
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN) :: Q(:)
      TYPE (CP) :: HQ
      real*8, intent(inout) :: QHQ(:,:)
      integer :: i,j,k,rF,rH,nbloc

      nbloc=SIZE(Q)
      rH=SIZE(H%coef)

      QHQ=0.d0

      do i=1,nbloc
         rF=SIZE(Q(i)%coef)
         do k=1,rH
            call CPMM(H,k,k,.FALSE.,Q(i),1,rF,.FALSE.,HQ)
!$omp parallel
!$omp do private(j)
            do j=i,nbloc
!              QHQ(i,j)=<Q(j),H(Q(i))>
               QHQ(i,j)=QHQ(i,j)+PRODVV(Q(j),HQ)
            enddo
!$omp enddo
!$omp end parallel
            call FlushCP(HQ)
         enddo
      enddo

      end subroutine GetQHQ2a

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ3(Q,H,QHQ,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. The matrix-vector product H*Q
! is computed and reduced via intertwining to save memory and time

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), ALLOCATABLE :: HQ(:)
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
         call FlushCP(HQ(i))
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
      TYPE (CP), INTENT(IN) :: Q(:)
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
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
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
      TYPE (CP), INTENT(IN) :: Q1(:),Q2(:)
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE   :: Qold(:)
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
         call NORMALIZE(Q(i))
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE   :: Qold(:)
      real*8, intent(in)  :: QHQ(:,:)
      integer, intent(in) :: nitn,lm
      integer :: i,nbloc

      nbloc=SIZE(Q)

!     Copy Q to Qold
      ALLOCATE(Qold(nbloc))

!$omp parallel
!$omp do private(i)
      DO i=1,nbloc
         Qold(i)=CopyCP(Q(i))
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE   :: Qold(:)
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
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

      function RayleighQuotient(v,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient for vector v and Hamiltonian H

      implicit none
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN) :: v
      TYPE (CP) :: w
      real*8 :: RayleighQuotient

!     w = H*v       
      call CPMM(H,.FALSE.,v,.FALSE.,w)

!     Rayleigh quotient = <v|w>/<v|v>
      RayleighQuotient=PRODVV(v,w)/PRODVV(v)
!
      call FlushCP(w)

      end function RayleighQuotient

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RayleighQuotient2(v,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient for vector v and Hamiltonian H

      implicit none
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN) :: v
      TYPE (CP) :: w
      real*8  :: RayleighQuotient2
      integer :: i,k,rF,rH

      rH=SIZE(H%coef)
      rF=SIZE(v%coef)

      RayleighQuotient2=0.d0

!     Calc w=H*v, term by term
      do i=1,rF
         do k=1,rH
            call CPMM(H,k,k,.FALSE.,v,i,i,.FALSE.,w)
!           Build <w|v>
            RayleighQuotient2=RayleighQuotient2+PRODVV(v,w)
            call FlushCP(w)
         enddo
      enddo

!     Rayleigh quotient = <v|w>/<v|v>
      RayleighQuotient2=RayleighQuotient2/PRODVV(v)

      end function RayleighQuotient2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RayleighResidual(v,H,eig)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the Rayleigh quotient residual for vector v and Hamiltonian H
! v should be normalized using NORMALIZE() before calling

      implicit none
      TYPE (CP), INTENT(IN) :: H
      TYPE (CP), INTENT(IN) :: v
      TYPE (CP) :: w
      real*8, intent(in) :: eig
      real*8 :: RayleighResidual

!     w = H*v - eig*v
      call CPMM(H,.FALSE.,v,.FALSE.,w)
      call SUMVECVEC(w,1.d0,v,-eig)

!     If the rank of w is large, then the call to PRODVV below will be
!     expensive. If lower accuracy is tolerable, one can reduce and then
!     call PRODVV with a smaller vector
!      call reduc(w)

!     Calculate ||w||
      RayleighResidual=sqrt(abs(PRODVV(w)))

      call FlushCP(w)

      end function RayleighResidual

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
