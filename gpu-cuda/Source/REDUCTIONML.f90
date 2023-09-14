!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module REDUCTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces the size of a CP-tensor using alternating least squares,
! singular value decomposition, etc.

      use ERRORTRAP
      use UTILS
      use SEPDREPN
      use MODVECVEC
      use LINALG
      use REDORTHO

      implicit none
      integer, private :: newrank,nloop
      real*8, private  :: redn_time=0.d0,eps
      logical, private :: RED_SETUP=.FALSE.
      logical, private :: printredn=.FALSE.
      character(3), private :: red2D='N/A',redND='N/A'

      INTERFACE reduc
         MODULE PROCEDURE reduc_AA,reduc_AB
      END INTERFACE reduc

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetReductionParameters(rrnk,nals,redtol,prnt,rtyp2D,&
                                        rtypND)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer, intent(in) :: rrnk,nals
      real*8, intent(in)  :: redtol
      logical, intent(in) :: prnt
      character*3, intent(in) :: rtyp2D,rtypND

      IF (.NOT. RED_SETUP) call SetupReduction()

      newrank=rrnk
      nloop=nals
      eps=redtol
      printredn=prnt
      red2D=rtyp2D
      redND=rtypND

      end subroutine SetReductionParameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupReduction()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      redn_time=0.d0
      RED_SETUP=.TRUE.

      end subroutine SetupReduction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeReduction()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. RED_SETUP) call SetupReduction()

      RED_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total reduction time (REDUCTION)  (s)',&
                            redn_time

      end subroutine DisposeReduction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getRednMem(rG,rF,nbas,ncpu,useSVD)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines memory needed for reduction

      implicit none
      integer, intent(in) :: rG,rF,ncpu
      integer, intent(in) :: nbas(:)
      logical, intent(in) :: useSVD
      integer :: n
      real*8  :: getRednMem

      IF (.NOT. RED_SETUP) call SetupReduction()

      IF (useSVD) THEN ! Reduce ranks with SVD
         n=MINVAL(nbas)
!!!      WORK ARRAY SIZE NOT YET INCLUDED BELOW
         getRednMem=REAL(ncpu)*(nbas(1)*nbas(2)+n*(nbas(1)+nbas(2)+1))
      ELSE ! Reduce ranks with ALS
         n=MAXVAL(nbas)
!        Memory allocated for PS,BB,BBmem and either PSk or bjk+IPV
!        (PSk and bjk+IPV are never allocated simultaneously)
         getRednMem=REAL(ncpu)*rF*(rG+2*rF+MAX(rG,n+1))
      ENDIF

      end function getRednMem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_AA(G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for various reduction algorithms. This version takes a single
! vector as input and generates the initial guess in the case of ALS

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: G
      TYPE (CPvec) :: F
      real*8 :: rederr

!     No reduction requested
      IF (newrank.eq.0 .or. SIZE(G%nbas).eq.1 .or. SIZE(G%coef).eq.1 &
         .or. (SIZE(G%coef).le.newrank .and. &
         .not.(SIZE(G%nbas).eq.2 .and. red2D.eq.'SVD'))) THEN
         rederr=0.d0
         IF (printredn) write(*,'(7X,A,f20.12)') '||F-G||_2 = ',rederr
      ELSE
         call reduc_AB(F,G)
         call ReplaceVwithW(G,F)
      ENDIF

      end subroutine reduc_AA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_AB(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for various reduction algorithms. This version uses an
! existing vector (F) as the initial guess if ALS is requested
! G may be normalized if a switch to RID is needed

      implicit none
      TYPE (CPvec), INTENT(IN) :: G
      TYPE (CPvec), INTENT(INOUT) :: F
      real*8  :: rederr,t1,t2
      logical :: useA,ok
      integer :: badF

      IF (.NOT. RED_SETUP) call SetupReduction()

!     No reduction requested
      IF (newrank.eq.0 .or. SIZE(G%nbas).eq.1 .or. SIZE(G%coef).eq.1) &
         THEN
         call FlushCPvec(F)
         call CopyWtoV(F,G)
         rederr=0.d0
         IF (printredn) write(*,'(7X,A,f20.12)') '||F-G|| = ',rederr
         RETURN
      ENDIF

      call CPU_TIME(t1)

!     SVD reduction for 2D case: determine which algorithm to use
!     (depends on base-size and rank of G)
      IF (SIZE(G%nbas).eq.2 .and. red2D.eq.'SVD') THEN
         call FlushCPvec(F)
         useA=whichSVDcode(G)
         IF (useA) THEN
            call reduc_SVDa(G,F,newrank)
         ELSE
            call reduc_SVDb(G,F,newrank)
         ENDIF

!     ALS or RID reduction (set NALS=0 to do RID by itself)
      ELSEIF ((SIZE(G%nbas).eq.2 .and. red2D.eq.'ALS') .or. &
              (SIZE(G%nbas).gt.2 .and. redND.eq.'ALS')) THEN

!        The initial guess is the existing vector F, but if the rank
!        of F is too small then we must generate a new guess
         IF (.not.ALLOCATED(F%coef)) THEN
            call reduc_RID(G,F,newrank)
!            call GetRandomCPvec(F,G%nbas,newrank)
         ELSEIF (SIZE(F%coef).ne.newrank .or. nloop.eq.0) THEN
            call FlushCPvec(F)
            call reduc_RID(G,F,newrank)
!            call GetRandomCPvec(F,G%nbas,newrank)
         ENDIF
         IF (SIZE(G%coef).gt.newrank) call reduc_ALS(G,F,abs(nloop))

!     Successive rank-1 approximations (uses ALS code)
      ELSEIF ((SIZE(G%nbas).eq.2 .and. red2D.eq.'SR1') .or. &
              (SIZE(G%nbas).gt.2 .and. redND.eq.'SR1')) THEN

         IF (.not.ALLOCATED(F%coef)) THEN
!            call reduc_RID(G,F,newrank)
            call GetRandomCPvec(F,G%nbas,newrank)
         ELSEIF (SIZE(F%coef).ne.newrank .or. nloop.eq.0) THEN
            call FlushCPvec(F)
!            call reduc_RID(G,F,newrank)
            call GetRandomCPvec(F,G%nbas,newrank)
         ENDIF
         IF (SIZE(G%coef).gt.newrank) call reduc_SR1(G,F,abs(nloop))

!     Orthogonal rank-1 projections plus sorting
      ELSEIF ((SIZE(G%nbas).eq.2 .and. red2D.eq.'OPS') .or. &
              (SIZE(G%nbas).gt.2 .and. redND.eq.'OPS')) THEN
         call FlushCPvec(F)
         call reduc_orthog(G,F,5*newrank)
         call reduc_bysorting(F,newrank)

!     Rank-1 reduction with orthogonal rank-1 projection
      ELSEIF ((SIZE(G%nbas).eq.2 .and. red2D.eq.'ROP') .or. &
              (SIZE(G%nbas).gt.2 .and. redND.eq.'ROP')) THEN
         call FlushCPvec(F)
         call reduc_rop(G,F,newrank,abs(nloop))

!     Error out if type is not recognized or SVD with >2D
      ELSEIF (SIZE(G%nbas).gt.2 .and. redND.eq.'SVD') THEN
         call AbortWithError('reduc(): SVD works only for D = 2')
      ELSE
         write(*,*) 'Reduction types :',red2D,',',redND
         call AbortWithError('reduc(): unrecognized reduction type')
      ENDIF

!     Calculate ||F-G|| after reduction if requested
      IF (printredn) THEN
         rederr=calc_FmG(G,F)
         write(*,'(7X,A,f20.12)') '||F-G||_2 = ',rederr
         IF (rederr.ne.rederr) THEN
            call AbortWithError("Reduction failed with NaN values")
         ENDIF
      ENDIF

      call CPU_TIME(t2)
      redn_time=redn_time+t2-t1

      end subroutine reduc_AB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function whichSVDcode(G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the operations counts of the SVD-a and SVD-b algorithms
! SVD-a is the standard algorithm; SVD-b is from Oseledets. On return,
! whichSVDcode=.TRUE. indicates that standard algorithm is to be used.

      implicit none
      TYPE (CPvec), INTENT(IN) :: G
      real*8  :: n1,n2,rnk,k,k1,k2,l,l1,l2,kk,lk,aops,bops
      logical :: whichSVDcode

      n1=REAL(G%nbas(1))
      n2=REAL(G%nbas(2))
      rnk=REAL(SIZE(G%coef))
      k=min(n1,n2)
      k1=min(n1,rnk)
      k2=min(n2,rnk)
      l=max(n1,n2)
      l1=max(n1,rnk)
      l2=max(n2,rnk)
      kk=min(k1,k2)
      lk=max(k1,k2)

!     Calculate the number of operations (aops) for the 'a' algorithm
!     a-1) Convert G to matrix form: G ->  M
      aops=k*(rnk+l*(rnk+1))
!     a-2) SVD of M: M = U * S * V
      aops=aops+4*k**2*(l-k/3)

!     Calculate the number of operations (bops) for the 'b' algorithm
!     b-1) copy arrays: G -> U,W
      bops=k*rnk
!     b-2) QR-step: U = Qu * Ru; V = Qv * Rv
      bops=bops+2*k1**2*(l1-k1/3)+2*k2**2*(l2-k2/3)
!     b-3) form P matrix: P = R_u * R_v^T
      bops=bops+k1*rnk*k2
!     b-4) SVD of P: P = X * D * Y
      bops=bops+4*kk**2*(lk-kk/3)
!     b-5) form U' and V' matrices: U' = Qu * X; V' = Qv * Y
      bops=bops+n1*k1**2+n2*k2**2

!      write(*,*) 'aops = ,',aops,'; bops = ',bops

      whichSVDcode=.TRUE.
      IF (aops.gt.bops) whichSVDcode=.FALSE.

      end function whichSVDcode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_SVDa(G,F,rF)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduction using SVD (2D system only). The 'a' algorithm generates a 
! matrix of size nbas(1)*nbas(2), independent of rank, and is faster
! when nbas(1)*nbas(2) is small (<~1600)

      implicit none
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (CPvec), INTENT(IN)  :: G
      real*8, allocatable :: A(:,:),U(:,:),VT(:,:),svals(:)
      integer, intent(in) :: rF
      integer :: nr,nc,nrk,ncomp,i,j,k
      real*8  :: cuttol

!     Convert CP-format vector to matrix form
      call CP2DtoMat(G,A)

!     Use SVD to get optimal CP-repn
      call SolveWithSVD(svals,A,U,VT)
      DEALLOCATE(A) 

!     Reconstruct the CP-vec, compressing by removing components that
!     are smaller than eps*largest_sval
      cuttol=eps*svals(1)
      nr=G%nbas(1)
      nc=G%nbas(2)
      ncomp=1
      DO i=2,SIZE(svals)
         IF (abs(svals(i)).lt.cuttol) EXIT
         ncomp=i
      ENDDO
      nrk=min(nr,nc,rF,ncomp)

      call RebuildfromSVD(F,U,svals,VT,nrk,.TRUE.)
      DEALLOCATE(U,VT,svals)

      end subroutine reduc_SVDa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_SVDb(G,F,rF)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduction using SVD (2D system only); the 'b' algorithm avoids
! multiplying out the CP-format array by using QR decomposition, and
! performs SVD on a smaller matrix. This is faster if nbas(1)*nbas(2)
! is large (>~1600)

      implicit none
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (CPvec), INTENT(IN)  :: G
      real*8, allocatable :: U(:,:),W(:,:),svals(:)
      real*8, allocatable :: Qu(:,:),Qw(:,:),X(:,:),Y(:,:)
      integer, intent(in) :: rF
      integer :: nr,nc,nrk,ncomp,i,j,k
      real*8  :: cuttol

      nr=G%nbas(1)
      nc=G%nbas(2)

!     Separate G into 2 mode-dependent matrices
      CALL CP2DtoUW(G,U,W)

!     QR decomposition of U and W
      CALL QRdecomp(Qu,U)
      CALL QRdecomp(Qw,W)

!     U <- U*W^T (product of two 'R' matrices from QR-decomposition)
      CALL MatrixMult(U,.FALSE.,W,.TRUE.)
      DEALLOCATE(W)

!     SVD of U: U=X*svals*Y^T
      CALL SolveWithSVD(svals,U,X,Y)
      DEALLOCATE(U)

!     Build new bases 
      CALL MatrixMult(Qu,.FALSE.,X,.FALSE.)
      CALL MatrixMult(Qw,.FALSE.,Y,.TRUE.)
      DEALLOCATE(X,Y)

!     Reconstruct the CP-vec, compressing by removing components that
!     are smaller than eps*largest_sval
      cuttol=eps*svals(1)
      ncomp=1
      DO i=2,SIZE(svals)
         IF (abs(svals(i)).lt.cuttol) EXIT
         ncomp=i
      ENDDO
      nrk=min(nr,nc,rF,ncomp)

      call RebuildfromSVD(F,Qu,svals,Qw,nrk,.FALSE.)
      DEALLOCATE(Qu,Qw,svals)

      end subroutine reduc_SVDb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RebuildfromSVD(F,U,svals,VT,nrk,tv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reconstruct a vector in CP-format from U*S*V^T SVD arrays
! If VT is transposed (the usual case) then set 'tv' to .TRUE.

      implicit none
      TYPE (CPvec), INTENT(OUT) :: F
      real*8, intent(in) :: U(:,:),VT(:,:),svals(:)
      integer, intent(in) :: nrk
      logical, intent(in) :: tv
      integer :: nbas(2)
      integer :: nr,nc,i,j

      IF (nrk.lt.1 .or. nrk.gt.SIZE(svals)) THEN
         write(*,'(X,A,I0,A,2(I0))') 'Requested rank of ',nrk,&
              ', is out of bounds: ',1,SIZE(svals)
         call AbortWithError("Error in RebuildfromSVD()")
      ENDIF

      nr=SIZE(U,1)
      nc=SIZE(VT,2)
      IF (.not.tv) THEN ! VT not transposed after all...
         nc=SIZE(VT,1)
      ENDIF
      nbas=(/nr,nc/)
      call NewCPvec(F,nbas,nrk)

!     Fill F
      F%coef(1:nrk)=svals(1:nrk)
      DO i=1,nrk
         DO j=1,nr
            F%base(j,i)=U(j,i)
         ENDDO
         IF (tv) THEN
            DO j=1,nc
               F%base(nr+j,i)=VT(i,j)
            ENDDO
         ELSE
            DO j=1,nc
               F%base(nr+j,i)=VT(j,i)
            ENDDO
         ENDIF
      ENDDO

      end subroutine RebuildfromSVD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_rop(G,F,nrk,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Takes vector G as input, and returns reduced vector F by reducing to 
! rank-1 and projecting onto the nrk dominant configurations

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (CPvec)   :: T
      TYPE (Configs) :: B,C
      integer, intent(in) :: nrk,nitn
      integer :: minn

      minn=MINVAL(G%nbas)

!     Get the "best" rank-1 approximation to G using SR1
      call GetRandomCPvec(T,G%nbas,1)
      call reduc_SR1(G,T,nitn)

!     Find the dominant configurations in T and store in B
      call GetConfigList(T,nrk*minn,B)

!     Compute the overlaps of B with G, store in C
      call ConfigSetOverlap(B,C,G)

!     Convert C to a CP-vector
      call Config2CP(F,C)

!     Reduce by sorting
      call reduc_bysorting(F,nrk)

      call FlushCPvec(T)
      call FlushConfigs(B)
      call FlushConfigs(C)

      end subroutine reduc_rop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_RID(G,F,rF)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Takes vector G as input, and returns reduced vector F using
! randomized interpolative decomposition, as described in Biagioni and
! Beylkin^2, J. Comput. Phys. 281 (2015) 116-134.

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (CPvec), INTENT(OUT) :: F
      integer, intent(in)   :: rF
      integer, allocatable  :: Lk(:)
      real*8, allocatable   :: Y(:,:),P(:)
      integer :: rG,rexc,rrx,i,j,jmod,n
      character*64 :: frmt

!     Set parameters
      rG=SIZE(G%coef)
      rexc=8  ! Excess rank
      rrx=min(rF+rexc,rG)
      n=rG-rF

!     If G is already small enough, no need to reduce
      IF (rG.le.rF) THEN
         call CopyWtoV(F,G)
         RETURN
      ENDIF

!     Guess a random start vector (eqn. 3.11)
      call GetRandomCPvec(F,G%nbas,rrx)

!     Normalize F
      call NORMBASE(F)
      F%coef=1.d0

!     Form the projection matrix Y (eqn. 1.15)
      ALLOCATE(Y(rrx,rG))
      call CONSTPSTOT(G,F,0,Y)

!     Multiply entries of PSTOT by the coefficients of G
!     No need to do this for F since coefs are already 1
      DO j=1,rG
         Y(:,j)=Y(:,j)*G%coef(j)
      ENDDO
      call FlushCPvec(F)

!     Calc. matrix interpolative decomposition of Y
      call MatrixID(Y,rF,Lk,P)

!     Construct the reduced-rank CP-vec, store in F
      call NewCPvec(F,G%nbas,rF)
      DO i=1,rF
         F%base(:,i)=G%base(:,Lk(i))
         F%coef(i)=1.d0
         DO j=1,n
            jmod=(j-1)*rF+i
            F%coef(i)=F%coef(i)+P(jmod)
         ENDDO
!        Make sure coef is positive
         F%coef(i)=F%coef(i)*G%coef(Lk(i))
         IF (F%coef(i).lt.0.d0) THEN
            F%coef(i)=abs(F%coef(i))
            F%base(1:F%nbas(1),i)=-F%base(1:F%nbas(1),i)
         ENDIF
      ENDDO

      DEALLOCATE(Y,Lk,P)

      end subroutine reduc_RID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixID(A,k,Lk,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for Mark Tygert's matrix interpolative decomposition code
! Computes columnwise matrix interpolative decomposition of matrix A,
! see Cheng et al, SIAM J. Sci. Comput. 26 (2005) 1389-1404.; and
! Martinsson et al, Appl. Comput. Harm. Anal. 20 (2011) 47-68.

      implicit none
      real*8, intent(in)   :: A(:,:)
      integer, intent(in)  :: k
      integer, allocatable, intent(out) :: Lk(:)
      real*8, allocatable, intent(out)  :: P(:)
      real*8, allocatable :: W(:)
      integer :: m,n,wdim

      m=SIZE(A,1)
      n=SIZE(A,2)

!     Initialize the w array
      wdim=MAX(10000000,(2*k+17)*n+27*m+100)
      ALLOCATE(W(wdim),Lk(n),P(k*(n-k)))

! The Matrix-ID code includes FORTRAN 'save' variables and is not thread
! safe, so calls to it must be OMP CRITICAL

!$OMP CRITICAL
      call iddr_aidi(m,n,k,W)
!     Columnwise matrix-ID of A
      call iddr_aid(m,n,A,k,W,Lk,P)
!$OMP END CRITICAL

      DEALLOCATE(W)

      end subroutine MatrixID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_ALS(G,F,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Receives G (represented in reduced form) and F (trial vector to be
! reduced), and approximates G by F using alternating least squares.
! The subroutine optimizes F in dimensions k=1,ndim in succession.
! Setting kortho equal to one of the DOFs will result in an orthogonal
! set of fs for that DOF provided that the rank, rF, is equal to n, the
! basis size for that DOF
 
      implicit none
      TYPE (CPvec), INTENT(IN) :: G
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in)  :: nitn
      integer, allocatable :: IPV(:)
      real*8, dimension (:,:), allocatable :: BB,PS,bjk,BBmem
      real*8  :: valpen,rnorm,gtmp
      integer :: rG,rF,ndim,i,j,ir,imod,k,l,n,info,gst,itn,kp
      logical :: update

      TYPE (CPvec) :: v

      IF (nitn.eq.0) return

!     Set parameters
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      ndim=SIZE(G%nbas)

!     Update the ALS matrices BB and PS if ndim > 3. For ndim <= 3 it
!     is faster to build BB and PS at each iteration.
!     Building also avoids possible zero division during the update.
      update=.TRUE.
      IF (ndim.le.3)  update=.FALSE.

      allocate(BB(rF,rF),BBmem(rF,rF),PS(rG,rF))

!     Penalty to avoid bad conditioning
      valpen=maxval(G%coef)*1.d-15
!      write(*,*) 'reduc_ALS(): valpen = ',valpen

76    continue

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
!     PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!     If the ALS matrices BB and PS are to be updated, initialize
!     them here with the first DOF removed

      IF (update) THEN
         call CONSTBBTOT(F,1,BBmem)
         call CONSTPSTOT(F,G,1,PS)
      ENDIF

!     Main loop over ALS iterations
      kp=0
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Update the BB and PS matrices of the linear system. This
!           also requires copying BBmem <-- BB since LAPACK destroys BB
            IF (update) THEN
               IF (kp.ne.0) then
                  call UpdateBB(F,k,BBmem,.TRUE.)
                  call UpdatePS(F,G,k,PS,.TRUE.)
               ENDIF
               BB=BBmem

!           Alternatively, build BB and PS from scratch
            ELSE 
               call CONSTBBTOT(F,k,BB)
               call CONSTPSTOT(F,G,k,PS)
            ENDIF

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4)
            allocate(bjk(rF,n))
            bjk=0.d0
            do ir=1,n
               imod=gst+ir
               do j=1,rG
                  gtmp=G%coef(j)*G%base(imod,j)
                  do i=1,rF
                     bjk(i,ir)=bjk(i,ir)+gtmp*PS(j,i)
                  enddo
               enddo
            enddo

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            ALLOCATE(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            DEALLOCATE(IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'reduc_ALS(): NaN on update; itn = ',itn
                  call FlushCPvec(F)
                  call GetRandomCPvec(F,G%nbas,rF)
                  update=.FALSE.
                  GOTO 76
               ENDIF
            ENDDO

!           Update BB and PS with new Fs
            IF (update) THEN
               IF (itn.lt.nitn .or. k.lt.ndim) THEN
                  call UpdateBB(F,k,BBmem,.FALSE.)
                  call UpdatePS(F,G,k,PS,.FALSE.)
               ENDIF
            ENDIF

            gst=gst+n
            kp=k
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      deallocate(BB,BBmem,PS)

      end subroutine reduc_ALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_ALT(G,F,nitn,kortho)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Receives G (represented in reduced form) and F (trial vector to be
! reduced), and approximates G by F using the orthogonal low-rank tensor
! approximation method of Wang et al. SIAM J Mat Anal Appl, 36 (2015) 1.
! The subroutine optimizes F in dimensions k=1,ndim in succession.
! Setting kortho equal to one of the DOFs will result in an orthogonal
! set of fs for that DOF provided that the rank, rF, is <= n, the
! basis size for that DOF
! Unlike ALS, matrices are reconstructed instead of updated since the
! orthogonality condition makes zero division more likely in the update

      implicit none
      TYPE (CPvec), INTENT(IN) :: G
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in) :: nitn,kortho
      real*8, dimension(:,:), allocatable :: PS,bjk,U,VT
      real*8, allocatable :: svals(:)
      real*8  :: rnorm,gtmp
      integer :: rG,rF,ndim,i,j,ir,imod,k,l,n,gi,gf,itn

      IF (nitn.eq.0) return

!     Set parameters
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      ndim=SIZE(G%nbas)

!     Error checking
      IF (rF.gt.G%nbas(kortho)) &
         call AbortWithError('reduc_ALT(): rF > nbas(kortho)')

      allocate(PS(rG,rF))

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension k
         gi=1
         do k=1,ndim
            n=F%nbas(k)
            gf=gi+n-1

!           Build PS matrix from scratch
            call CONSTPSTOT(F,G,k,PS)

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4)
            allocate(bjk(n,rF))
            bjk=0.d0
            do ir=1,n
               imod=gi+ir-1
               do j=1,rG
                  gtmp=G%coef(j)*G%base(imod,j)
                  do i=1,rF
                     bjk(ir,i)=bjk(ir,i)+gtmp*PS(j,i)
                  enddo
               enddo
            enddo

!           Normalize the bjk vectors
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(:,i),bjk(:,i))))
               bjk(:,i)=bjk(:,i)/rnorm
            enddo

!           Polar decomposition orthogonalizes bjk for the kortho-th DOF
            IF (k.eq.kortho) THEN
               call SolveWithSVD(svals,bjk,U,VT)
               deallocate(bjk)
               call MatrixMult(U,.FALSE.,VT,.FALSE.,bjk)
               deallocate(U,VT,svals)
            ENDIF

!           New base of F
            F%base(gi:gf,1:rF)=bjk(1:n,1:rF)
            deallocate(bjk)

            gi=gf+1
         enddo  ! loop over k
      ENDDO  ! loop over iterations

!     Compute the final coefficients of F. Most of the work required to
!     compute these is already in PS, so one only needs to multiply PS
!     by PSk for k = ndim and then scale by the coefs of G to get <F,G>
      call UpdatePS(F,G,ndim,PS,.FALSE.)

      ALLOCATE(svals(rG))
      svals(:)=G%coef(:)
      call MatVecProd(PS,.TRUE.,svals)
      F%coef(:)=svals(:)

      deallocate(PS,svals)

      end subroutine reduc_ALT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_SR1(G,F,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Receives CP-vector G and F (trial vector to be reduced), and 
! approximates G by F using successive rank-1 approximations. 
 
      implicit none
      TYPE (CPvec), INTENT(IN)    :: G
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec) :: T
      integer, intent(in) :: nitn
      real*8, dimension (:,:), allocatable :: BB,PS,&
      PSk,PSkm1,BBk,BBkm1
      real*8  :: valpen,oldcoef,Tcoef,snorm
      integer :: rF,rG,rrG,ndim,i,j,ir,k,l,n,itn,kp,gi,gf,nbad
      real*8, parameter :: tol=1.d-15
      logical :: update

      IF (nitn.eq.0) return

!     Set parameters
      nbad=0
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      ndim=SIZE(G%nbas)
      snorm=1.d0
      valpen=maxval(G%coef)*1.d-15

!     Update the ALS matrices BB and PS if ndim > 3. For ndim <= 3 it
!     is faster to build BB and PS at each iteration.
!     Building also avoids possible zero division during the update.
      update=.TRUE.
      IF (ndim.le.3)  update=.FALSE.

      allocate(BB(1,1),BBk(1,1),BBkm1(1,1))
      allocate(PS(rG+rF,1),PSk(rG+rF,1),PSkm1(rG+rF,1))

      call NewCPvec(T,G%nbas,1)

!     Loop over rF rank-1 terms to be computed
      DO i=1,rF 
         rrG=rG+i-1

!        If a restart is needed due to update failure, restart on
!        the term where the update failed
78       continue

!        Copy a term in F to T for the initial guess
         call GenCopyWtoV(T,F,1,1,i,i)

!        Construct BB and PS "matrices"
         IF (update) THEN
            call CONSTBBk(T,2,BB)
            call CONSTPSk(T,G,2,PS(1:rG,:))
            IF (i.gt.1) call CONSTPSk(T,F,2,PS(rG+1:,:))
            DO k=3,ndim
               call CONSTBBk(T,k,BBk)
               call CONSTPSk(T,G,k,PSk(1:rG,:))
               IF (i.gt.1) call CONSTPSk(T,F,k,PSk(rG+1:,:))
               BB=BB*BBk
               PS(1:rrG,:)=PS(1:rrG,:)*PSk(1:rrG,:)
            ENDDO
         ENDIF

!        Loop over ALS iterations
         kp=0
         oldcoef=1.d99
         DO itn=1,nitn

!           Loop over DOFs
            gi=1
            do k=1,ndim

               n=T%nbas(k)
               gf=gi+n-1

!              Update the BB and PS matrices of the linear system
!              Terms in F and G must be used to update PS
               IF (update) THEN
                  IF (kp.ne.0) THEN
                     call CONSTBBk(T,k,BBk)
                     call CONSTBBk(T,kp,BBkm1)
                     call CONSTPSk(T,G,k,PSk(1:rG,:))
                     call CONSTPSk(T,G,kp,PSkm1(1:rG,:))
                     IF (i.gt.1) THEN 
                        call CONSTPSk(T,F,k,PSk(rG+1:,:))
                        call CONSTPSk(T,F,kp,PSkm1(rG+1:,:))
                     ENDIF
                     BB=BB*BBkm1/BBk
                     PS(1:rrG,:)=PS(1:rrG,:)*PSkm1(1:rrG,:)/PSk(1:rrG,:)
                  ENDIF

!              For 2D and 3D, build BB and PS from scratch each time ->
!              faster than updating and avoids possible zero division
               ELSE 
                  kp=mod(k,ndim)+1
                  call CONSTBBk(T,kp,BB)
                  call CONSTPSk(T,G,kp,PS(1:rG,:))
                  IF (i.gt.1) call CONSTPSk(T,F,kp,PS(rG+1:,:))
                  DO l=2,ndim-1
                     kp=mod(k+l-1,ndim)+1
                     call CONSTBBk(T,kp,BBk)
                     call CONSTPSk(T,G,kp,PSk(1:rG,:))
                     IF (i.gt.1) call CONSTPSk(T,F,kp,PSk(rG+1:,:))
                     BB=BB*BBk
                     PS(1:rrG,:)=PS(1:rrG,:)*PSk(1:rrG,:)
                  ENDDO
               ENDIF

!              Penalty to avoid ill-conditioning (section 3.2, Beylkin)
               BB(1,1)=BB(1,1)+valpen

!              Apply the ALS update to T
               T%base(gi:gf,1)=0.d0
               T%coef(1)=0.d0
               do ir=gi,gf
                  do j=1,rG ! Sum over terms in G
                     T%base(ir,1)=T%base(ir,1) &
                                 +G%coef(j)*G%base(ir,j)*PS(j,1)
                  enddo
                  do j=1,i-1 ! Sum over 'extra' terms of G stored in F
                     T%base(ir,1)=T%base(ir,1) &
                                 +F%coef(j)*F%base(ir,j)*PS(rG+j,1)
                  enddo
                  T%base(ir,1)=T%base(ir,1)/BB(1,1)
                  T%coef(1)=T%coef(1)+T%base(ir,1)**2
               enddo
!              Normalization
               T%coef(1)=sqrt(abs(T%coef(1)))
               T%base(gi:gf,1)=T%base(gi:gf,1)/T%coef(1)

!              Remove penalty on BB after update
               BB(1,1)=BB(1,1)-valpen

               IF (T%coef(1).ne.T%coef(1)) THEN
!                  call AbortWithError('reduc_SR1(): NaN on update')
                  write(*,*) 'reduc_SR1(): NaN on update; i = ',i,&
                             'itn = ',itn
                  nbad=nbad+1
                  IF (nbad.gt.100) THEN
                      write(*,*) 'NaN reducing G on rank:',i
                      call PrintCPvec(G)
                      call AbortWithError('reduc_SR1(): NaN on update')
                  ENDIF
                  call FlushCPvec(T)
                  call NewCPvec(T,G%nbas,1)
                  update=.FALSE.
                  GOTO 78
               ENDIF

               gi=gi+n
               kp=k
            enddo  ! loop over k

!           Check for convergence on the coefficient
            IF (abs(oldcoef-T%coef(1))/snorm.lt.tol) EXIT

            oldcoef=T%coef(1)
         ENDDO  ! loop over iterations

!        After computing the first term, set the s-norm
         IF (i.eq.1) snorm=T%coef(1)

!        T must be subtracted from G before computing the next term.
!        Instead of resizing G it is cheaper to accumulate -T in F
         call GenCopyWtoV(F,T,i,i,1,1)
         call VecSignChange(F,i,i)
      ENDDO  ! loop over i

      call FlushCPvec(T)
      deallocate(BB,BBk,BBkm1,PS,PSk,PSkm1)

!     Change back sign of F
      call VecSignChange(F,1,rF)

      end subroutine reduc_SR1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module REDUCTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
