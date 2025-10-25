!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines which operate on CP-blocks

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MUNKRES
      USE CPMMM
      USE MODVECVEC
      USE SEPDREPN
      USE REDUCTION
      USE ALSPOW
      USE ALSUTILS
      USE ALSDRVR

      USE CPr8
      USE ALS8DRVR

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
      integer :: nbloc,algo

!     Set parameters
      nbloc=SIZE(Q)

      allocate(QHQ(nbloc,nbloc),S(nbloc,nbloc))

!     Calculate QHQ and S matrices
      IF (intw) THEN
         call GetQHQ_intw(Q,H,QHQ,nitn,lm)  ! Reduces H*Q, then calcs Q^THQ
         call GetOverlaps(Q,S)
      ELSE
         call GetQHQ(Q,H,QHQ)
         call GetOverlaps(Q,S)
      ENDIF

!     Diagonalize QHQ, accounting for overlaps
      call SolveGenEigval(avec,S,QHQ,'V')

!     Update block vectors after diagonalization
!     q^n_{new} <- sum_{i=1} ^ m U_{im} q^i_{old}
      IF (intw) THEN
         call UpdateVecs_intw(Q,QHQ,nitn,lm) ! Avoids long vectors
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
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), ALLOCATABLE :: HQ(:)
      real*8, intent(inout)  :: QHQ(:,:)
      integer, allocatable   :: mvecs(:),moffs(:)
      integer :: i,j,nbloc,sz,os,szmx,ierr

      nbloc=SIZE(Q)
      QHQ=0.d0
      
!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)
      ALLOCATE(HQ(sz))

      do j=1,sz
!        HQ=H*Q(os+j)
         call CPMM(H,.FALSE.,Q(os+j),.FALSE.,HQ(j))
!        Reduce the rank of HQ to accelerate computing the dot products
         call reduc(HQ(j))
!$omp parallel
!$omp do
         do i=1,os+j
!           QHQ(i,j)=<Q(i),H(Q(os+j))>
            QHQ(i,os+j)=PRODVV(Q(i),HQ(j))
         enddo
!$omp enddo
!$omp end parallel
         call FlushCP(HQ(j))
      enddo

      DEALLOCATE(HQ)

      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do j=1,nbloc
         i=mod(j-1,mpinodes)+1
         mvecs(i)=mvecs(i)+nbloc
      enddo
      moffs(1)=0
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

!     Gather the columns from all processors
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
           QHQ,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)
 
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

      subroutine GetQHQ_intw(Q,H,QHQ,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. The matrix-vector product H*Q
! is computed and reduced via intertwining to save memory and time

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
      TYPE (CP), ALLOCATABLE :: HQ(:)
      integer, intent(in)    :: nitn,lm
      real*8, intent(inout)  :: QHQ(:,:)
      integer, allocatable   :: mvecs(:),moffs(:)
      integer :: i,j,nbloc,sz,os,szmx,ierr

      nbloc=SIZE(Q)
      QHQ=0.d0

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

      ALLOCATE(HQ(sz))

      do j=1,sz
!        HQ=H*Q(os+j)
         call PRODHV_ALS_alg(Q(os+j),HQ(j),H,0,0.d0,nitn,lm)
!$omp parallel
!$omp do private(i)
         do i=1,os+j
!           QHQ(i,j)=<Q(i),H(Q(os+j))>
            QHQ(i,os+j)=PRODVV(Q(i),HQ(j))
         enddo
!$omp enddo
!$omp end parallel
      enddo

      DEALLOCATE(HQ)

      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do j=1,nbloc
         i=mod(j-1,mpinodes)+1
         mvecs(i)=mvecs(i)+nbloc
      enddo
      moffs(1)=0
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

!     Gather the columns from all processors
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
           QHQ,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)

      end subroutine GetQHQ_intw

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetOverlaps(Q,S)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T Q for a block of vectors

      implicit none
      TYPE (CP), INTENT(IN) :: Q(:)
      real*8, intent(inout) :: S(:,:)
      integer, allocatable   :: mvecs(:),moffs(:)
      integer :: i,j,nbloc,sz,os,szmx,ierr

      nbloc=SIZE(Q)
      S=0.d0

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Compute inner products between vectors on this rank and
!     those preceding them in the list
      do j=1,sz
!$omp parallel
!$omp do private(i)
         do i=1,os+j
!           S(i,os+j)=<Q(i),Q(os+j)>
            S(i,os+j)=PRODVV(Q(i),Q(os+j))
         enddo
!$omp enddo
!$omp end parallel
      enddo

      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do j=1,nbloc
         i=mod(j-1,mpinodes)+1
         mvecs(i)=mvecs(i)+nbloc
      enddo
      moffs(1)=0
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

!     Gather the columns from all processors
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
           S,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)

      end subroutine GetOverlaps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQdiag(Q,H,QHQd,nconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes diagonal elements of Q^T H Q for a block of vectors
! The first nconv vectors are skipped

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(IN)  :: Q(:)
      integer, intent(in)    :: nconv
      real*8, intent(inout)  :: QHQd(:)
      integer, allocatable   :: mvecs(:),moffs(:)
      integer :: b,i,nbloc,sz,os,szmx,ierr

      nbloc=SIZE(Q)

!     Shifts, number of (unconverged) vectors for MPI
      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do b=1,nbloc-nconv
         i=mod(b-1,mpinodes)+1
         mvecs(i)=mvecs(i)+1
      enddo
      moffs(1)=nconv
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc-nconv,sz,os,szmx)
      os=os+nconv

      do i=1,sz
         QHQd(os+i)=RayleighQuotient2(Q(os+i),H)   
      enddo

      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
              QHQd,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)

      end subroutine GetQHQdiag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateVecs(Q,QHQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces a block of vectors Q with the eigenvectors whose coefficients
! are contained in QHQ

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE   :: Qold(:)
      real*8, intent(in) :: QHQ(:,:)
      integer :: i,nbloc,sz,os,szmx

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

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!$omp parallel
!$omp do private(i)
      DO i=1,sz
         call GetEigenFxn(Q(i+os),Qold,QHQ,i+os)
!        reduce and normalize Q(i)
         call reduc(Q(i+os))
         call NORMALIZE(Q(i+os))
      ENDDO
!$omp enddo
!$omp end parallel

      DEALLOCATE(Qold)

      call MPI_Sync_CP_block(Q)

      end subroutine UpdateVecs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateVecs_intw(Q,QHQ,nitn,lm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces a block of vectors Q with the eigenvectors whose coefficients
! are stored in QHQ. This version uses intertwining to reduce the rank.

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE   :: Qold(:)
      real*8, intent(in)  :: QHQ(:,:)
      integer, intent(in) :: nitn,lm
      integer :: i,nbloc,sz,os,szmx

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

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Build and reduce the eigenfunction
!$omp parallel
!$omp do private(i)
      DO i=1,sz
         call ALS_SUMLCVEC_alg(Q(os+i),Qold,QHQ(:,os+i),nitn,lm)
      ENDDO
!$omp enddo
!$omp end parallel

      DEALLOCATE(Qold)

      call MPI_Sync_CP_block(Q)

      end subroutine UpdateVecs_intw

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

      subroutine pGRAMORTHO_CP8(Q,k,nals,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Orthogonalizes a block of CP-format vectors via Gram-Schmidt, assuming
! the first k vectors are already orthonormal

      implicit none
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8)  :: w
      integer, intent(in)  :: k,nals,algo
      integer, allocatable :: tab(:,:)
      real(kind=8), allocatable :: vivj(:)
      real(kind=8) :: conv
      integer :: i,j,nbloc

      nbloc=size(Q)

      IF (k.eq.nbloc) RETURN
      IF (k.lt.0 .or. k.gt.nbloc) &
         call AbortWithError('pGRAMORTHO(): k is out of range')

      call nvtx_start('pGRAMORTHO')

      ALLOCATE(vivj(nbloc))
      call GetRankOffsetTable(Q,tab)

!     Loop through the remaining vectors
!     for each, project out all preceding vectors
      do i=k+1,nbloc

!        Normalize i-th vector
         call NORMALIZE_CP8(Q(i),algo)

         IF (i.eq.1) CYCLE

!        Pack vectors into w; compute overlaps
         vivj(:)=1.d0
         call w%sumlccp(Q(1:i),vivj(1:i))
         call GetBlockOverlaps_CP8(Q(i),w,tab(1:i,:),vivj(1:i),algo)

!        Scale each vector in w by negative of overlap from previous step
!        Q(i) <- Q(i) - <Q(i),Q(j)>*Q(j)
         do j=1,i-1
            call w%mult_terms(-vivj(j),tab(j,1),tab(j,2))
         enddo

!        Reduce w into Q(i), normalize
         conv=ALS_reduce_CP8(Q(i),w,nals,algo)
         call w%flush()
         call NORMALIZE_CP8(Q(i),algo)
      enddo

      DEALLOCATE(vivj,tab)
      call nvtx_stop()

      end subroutine pGRAMORTHO_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Diagonalize_CP8(Q,H,avec,nals,algo,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! One-step function to solve the generalized eigenvalue problem and 
! update vectors.

      implicit none
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8), INTENT(IN) :: H
      integer, intent(in)   :: nals,algo
      logical, intent(in)   :: tiledH
      real(kind=8), intent(inout) :: avec(:)
      real(kind=8), allocatable :: QHQ(:,:),S(:,:)
      integer :: nbloc

!     Set parameters
      nbloc=SIZE(Q)

      allocate(QHQ(nbloc,nbloc),S(nbloc,nbloc))

!     Calculate QHQ and S matrices
      call GetQHQ_CP8(Q,H,QHQ,nals,algo,tiledH)
      call GetOverlaps_CP8(Q,S,algo)

!     Diagonalize QHQ, accounting for overlaps
      call SolveGenEigval(avec,S,QHQ,'V')

!     Update block vectors after diagonalization
!     q^n_{new} <- sum_{i=1} ^ m U_{im} q^i_{old}
      call UpdateVecs_CP8(Q,QHQ,nals,algo)

      deallocate(QHQ,S)

      end subroutine Diagonalize_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQ_CP8(Q,H,QHQ,nitn,algo,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T H Q for a block of vectors. The matrix-vector product H*Q
! is computed and reduced via intertwining to save memory and time

      implicit none
      TYPE (CP8), INTENT(IN)  :: H
      TYPE (CP8), INTENT(IN)  :: Q(:)
      TYPE (CP8) :: HQ,HQr,w,Qdummy
      integer, intent(in) :: nitn,algo
      logical, intent(in) :: tiledH
      real(kind=8), intent(inout) :: QHQ(:,:)
      real(kind=8), allocatable   :: vivj(:)
      real(kind=8) :: qhqdummy
      integer, allocatable   :: tab(:,:)
      integer :: i,j,k,nbloc,sz,os,szmx,conv
      logical :: reduceH

      reduceH=.TRUE.
      nbloc=SIZE(Q)
      QHQ=0.d0

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Rank offsets for vectors in packed array
      if (szmx.gt.0) then
         call nvtx_start('GetQHQ')
         ALLOCATE(vivj(nbloc))
         vivj(:)=1.d0
         call GetRankOffsetTable(Q,tab)

         do j=1,sz
            k=os+j
            if (tiledH) then 
               call HQ%copyfrom(Q(k))
               qhqdummy=ProdHV_ALS_tiled_CP8(H,0,0.d0,Q(k),HQ,nitn,algo)
            else
!              HQ=H*Q(k), then reduce rank
               call CPMM_CP8(H,Q(k),HQ,.FALSE.,algo)
               if (reduceH) then
                  call HQr%copyfrom(Q(k))
                  conv=ALS_reduce_CP8(HQr,HQ,nitn,algo)
                  call HQ%replace(HQr)
               endif
            endif

!           Pack vectors in Q into w
            call w%sumlccp(Q(1:k),vivj(1:k))
!           Get overlaps between HQr(k) and those in w
            call GetBlockOverlaps_CP8(HQ,w,tab(1:k,:),QHQ(1:k,k),algo)
            call w%flush()
            call HQ%flush()
         enddo

!        Dummy vector to preserve number of MPI calls for this rank when
!        using the tiled H algorithm
         if (sz.lt.szmx .and. tiledH) then
            call Qdummy%copyfrom(Q(1))
            call HQ%copyfrom(Qdummy)
            qhqdummy=ProdHV_ALS_tiled_CP8(H,0,0.d0,Qdummy,HQ,nitn,algo)
            call Qdummy%flush()
            call HQ%flush()
         endif

         DEALLOCATE(tab,vivj)
         call nvtx_stop()
      endif

      call MPIGatherBlockedMatrix(QHQ)

      end subroutine GetQHQ_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQXQ_CP8(Q,X,QXQ,algo,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T X Q for a block of vectors, where X is an operator acting
! on mode d. 

      implicit none
      TYPE (CP8), INTENT(IN) :: X
      TYPE (CP8), INTENT(IN) :: Q(:)
      TYPE (CP8) :: XQ,w
      integer, intent(in) :: algo,d
      real(kind=8), intent(inout) :: QXQ(:,:)
      real(kind=8), allocatable   :: vivj(:)
      integer, allocatable :: tab(:,:)
      integer :: i,j,k,nbloc,sz,os,szmx

      nbloc=SIZE(Q)
      QXQ=0.d0

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Rank offsets for vectors in packed array
      if (sz.gt.0) then
         call nvtx_start('GetQXQ')
         ALLOCATE(vivj(nbloc))
         vivj(:)=1.d0
         call GetRankOffsetTable(Q,tab)

         do j=1,sz
            k=os+j
!           XQ=X*Q(k)
            call XQ%copyfrom(Q(k))
            call CPMM_CP8(X,0,0.d0,.FALSE.,Q(k),0,0.d0,.FALSE.,&
                          XQ,d,.FALSE.,algo)
!           Pack vectors in Q into w
            call w%sumlccp(Q(1:k),vivj(1:k))
!           Get overlaps between HQr(k) and those in w
            call GetBlockOverlaps_CP8(XQ,w,tab(1:k,:),QXQ(1:k,k),algo)
            call w%flush()
            call XQ%flush()
         enddo

         DEALLOCATE(tab,vivj)
         call nvtx_stop()
      endif

      call MPIGatherBlockedMatrix(QXQ)

      end subroutine GetQXQ_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQHQdiag_CP8(Q,H,QHQd,nconv,algo,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes diagonal elements of Q^T H Q for a block of vectors
! The first nconv vectors are skipped

      implicit none
      TYPE (CP8), INTENT(IN)  :: H
      TYPE (CP8), INTENT(IN)  :: Q(:)
      TYPE (CP8) :: Qdummy
      integer, intent(in) :: nconv,algo
      logical, intent(in) :: tiledH
      real(kind=8), intent(inout) :: QHQd(:)
      real(kind=8) :: qhqdummy
      integer, allocatable   :: mvecs(:),moffs(:)
      integer :: b,i,k,nbloc,sz,os,szmx,ierr

      nbloc=SIZE(Q)

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc-nconv,sz,os,szmx)
      os=os+nconv

      call nvtx_start('GetQHQdiag')

      do i=1,sz
         k=os+i
         QHQd(k)=ALS_pow_alg_CP8(H,0,0.d0,Q(k),0,algo,tiledH)
      enddo

!     Dummy vector to preserve number of MPI calls for this rank when
!     using the tiled H algorithm
      if (sz.lt.szmx .and. tiledH) then
         call Qdummy%copyfrom(Q(1))
         qhqdummy=ALS_pow_alg_CP8(H,0,0.d0,Qdummy,0,algo,tiledH)
         call Qdummy%flush()
      endif

      call nvtx_stop()

!     Shifts, number of (unconverged) vectors for MPI
      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do b=1,nbloc-nconv
         i=mod(b-1,mpinodes)+1
         mvecs(i)=mvecs(i)+1
      enddo
      moffs(1)=nconv
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
              QHQd,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)

      end subroutine GetQHQdiag_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetOverlaps_CP8(Q,S,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes Q^T Q for a block of vectors.

      implicit none
      TYPE (CP8), INTENT(IN) :: Q(:)
      TYPE (CP8) :: w
      integer, intent(in) :: algo
      real(kind=8), intent(inout) :: S(:,:)
      real(kind=8), allocatable   :: vivj(:)
      integer, allocatable :: tab(:,:)
      integer :: i,j,k,nbloc,sz,os,szmx

      nbloc=SIZE(Q)
      S=0.d0

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Rank offsets for vectors in packed array
      if (sz.gt.0) then
         call nvtx_start('GetOverlaps')
         ALLOCATE(vivj(nbloc))
         vivj(:)=1.d0
         call GetRankOffsetTable(Q,tab)

         do j=1,sz
            k=os+j
!           Pack vectors in Q into w
            call w%sumlccp(Q(1:k),vivj(1:k))
!           Get overlaps between Q(k) and those in w
            call GetBlockOverlaps_CP8(Q(k),w,tab(1:k,:),S(1:k,k),algo)
            call w%flush()
         enddo

         DEALLOCATE(vivj,tab)
         call nvtx_stop()
      endif

      call MPIGatherBlockedMatrix(S)

      end subroutine GetOverlaps_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateVecs_CP8(Q,QHQ,nals,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces a block of vectors Q with the eigenvectors whose coefficients
! are stored in QHQ.

      implicit none
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8), ALLOCATABLE   :: Qold(:)
      TYPE (CP8) :: w
      real(kind=8), intent(in)  :: QHQ(:,:)
      integer, intent(in) :: nals,algo
      integer :: i,k,nbloc,sz,os,szmx,conv

      call nvtx_start('UpdateVecs')
      nbloc=SIZE(Q)

!     Copy Q to Qold
      ALLOCATE(Qold(nbloc))
      DO i=1,nbloc
         call Qold(i)%copyfrom(Q(i))
      ENDDO

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nbloc,sz,os,szmx)

!     Build and reduce each eigenfunction on this rank
      DO i=1,sz
         k=os+i
         call w%sumlccp(Qold,QHQ(:,k))
         conv=ALS_reduce_CP8(Q(k),w,nals,algo)
         call Q(k)%updatefromdevice()
         call w%flush()
      ENDDO

!     Clean up
      DO i=1,nbloc
         call Qold(i)%flush
      ENDDO
      DEALLOCATE(Qold)
      call nvtx_stop()

!     MPI update: remove vectors from device before the MPI call since
!     vectors on other MPI ranks may have had their ranks modified.
!     After MPI call, copy updated vecs back to device if applicable
      DO i=1,nbloc
         call Q(i)%deletefromdevice()
      ENDDO
      call MPI_Sync_block_CP8(Q)
      if (algo.eq.1) then
         DO i=1,nbloc
            call Q(i)%copyintodevice()
         ENDDO
      endif

      end subroutine UpdateVecs_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortVecs_CP8(Q,eigv,nconv,Eref)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts a block of vectors Q in terms of increasing distance from Eref
! The first nconv vectors are assumed to already be in order

      implicit none
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8), ALLOCATABLE   :: Qold(:)
      integer, intent(inout) :: nconv
      real(kind=8), intent(inout)  :: eigv(:)
      real(kind=8), intent(in)     :: Eref
      real(kind=8), allocatable    :: ix(:),edif(:)
      integer :: i,nbloc

      call nvtx_start('SortVecs')
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
      DO i=1,nbloc
         IF (i.ne.int(ix(i))) call Qold(i)%replace(Q(i))
      ENDDO

!     Eigenvector resort
      DO i=1,nbloc
         IF (i.ne.int(ix(i))) call Q(i)%replace(Qold(int(ix(i))))
!        Update all vectors here to ensure that the mpi_io_rank has the
!        current sorted vectors on host before SavePsi() is called
         call Q(i)%updatefromdevice()
      ENDDO

!     If a 'converged' vector swaps places with another vector
!     then it is not really converged, so change nconv accordingly!
      DO i=1,nconv
         IF (i.ne.int(ix(i))) THEN
            nconv=i-1
            EXIT
         ENDIF
      ENDDO

      DEALLOCATE(Qold,ix,edif)
      call nvtx_stop()

      end subroutine SortVecs_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetBlockOverlaps_CP8(v,w,tab,vivj,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes overlaps between v and vectors packed in w. The rank offsets
! for each vector in w are stored in the tab array

      implicit none
      TYPE (CP8), INTENT(IN) :: v,w
      TYPE (CP8) :: P
      integer, intent(in) :: tab(:,:)
      integer, intent(in) :: algo
      real(kind=8), intent(inout) :: vivj(:)
      integer :: i,nbloc
      integer :: rk,rows(1),cols(1),pst,pfi

      nbloc=SIZE(tab,1)

!     Create P matrix in CP8 format to manage data movement, indexing
      rows=(/v%R()/)
      cols=(/1/)
      rk=w%R()
      call P%new(rk,rows,cols)
      if (algo.eq.1) call P%createondevice()

!     Compute inner products on all vecs in w
      call CONSTPVV_CP8(v,w,P%base,algo)

!     Accumulate terms from each Q(i) in QHQ
      do i=1,nbloc
!        vivj=<v,w(pst:pfi)>
         pst=P%BS(tab(i,1),1)
         pfi=P%BF(tab(i,2),1)
         vivj(i)=ReduceP_CP8(P%base(pst:pfi),algo)
      enddo
      call P%flush()

      end subroutine GetBlockOverlaps_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRankOffsetTable(Q,tab)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets array of rank offsets for block Q

      implicit none
      TYPE (CP8), INTENT(IN) :: Q(:)
      integer, allocatable, intent(out) :: tab(:,:)
      integer :: i,nbloc

      nbloc=SIZE(Q)

      if (nbloc.lt.1) call &
         AbortWithError('GetRankOffsetTable(): nbloc < 1')

      ALLOCATE(tab(nbloc,2))
      tab(1,1)=1
      tab(1,2)=Q(1)%R()
      do i=2,nbloc
         tab(i,1)=tab(i-1,1)+Q(i-1)%R()
         tab(i,2)=tab(i-1,2)+Q(i)%R()
      enddo

      end subroutine GetRankOffsetTable

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MPIGatherBlockedMatrix(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gathers matrix M with cols stored on different MPI ranks

      implicit none
      real(kind=8), intent(inout) :: M(:,:)
      integer, allocatable :: mvecs(:),moffs(:),tab(:,:)
      integer :: i,j,nbloc,sz,os,ierr

      nbloc=SIZE(M,2)

      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do j=1,nbloc
         i=mod(j-1,mpinodes)+1
         mvecs(i)=mvecs(i)+nbloc
      enddo
      moffs(1)=0
      do i=2,mpinodes
         moffs(i)=moffs(i-1)+mvecs(i-1)
      enddo

!     Gather the columns from all processors
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
           M,mvecs,moffs,mpi_r8,mpi_comm_wd,ierr)

      DEALLOCATE(mvecs,moffs)

      end subroutine MPIGatherBlockedMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module BLOCKUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
