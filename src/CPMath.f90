!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module CPMATH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Math operations involving CP vectors

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE ALSDRVR
      USE LINSOLVER
      USE BLOCKUTILS

      INTERFACE CPMatNormA
         MODULE PROCEDURE CPMatNormA_noweights,CPMatNormA_weights
      END INTERFACE CPMatNormA

      INTERFACE CPsqrt
         MODULE PROCEDURE CPsqrt_noweights,CPsqrt_weights
      END INTERFACE CPsqrt

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatNorm(U,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes CP-tensor U

      implicit none
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP) :: S,R,W,W2,W2T
      integer, intent(in)  :: nals
      integer, allocatable :: cols(:)
      logical, parameter   :: reduceW=.TRUE.
      integer :: rk

      rk=SIZE(U%coef)

!     Compute the tensor of vector self-products
      W=CPMatSelfProdTot(U,U,.TRUE.)
      call CPMatrixTranspose(W)
!     (Optional) reduction of W
      IF (reduceW) THEN
         W2=RandomCP(W,rk)
         call reduc_ALS(W,W2,nals)
         call ReplaceVwithW(W,W2)
      ENDIF

      allocate(cols(SIZE(U%nbas)))
      cols(:)=1

!     Initial guess the square root, then run the square root iteration
      S=RandomCP(rk,W%rows,cols,W%sym)
      call SquareRootGuess(W,S,rk,102)
!!!   WAS                nals 100--^
!!!   WAS    nals 100--v   v--20 nals
      call CPsqrt(W,S,102,22)

!     Guess the reciprocal square root, then invert the square root
      R=RandomCP(rk,W%rows,cols,W%sym)
      call ReciprocalGuess(S,R,rk,nals)
      call CPreciprocal(S,R,nals)
      call FlushCP(S)
      call FlushCP(W)

!     Convert 1/sqrt(W) from vector repn to diagonal matrix, store as W
      W=CPVtoDiagMatrix(R,R%sym)
      call FlushCP(R)

!     Replace U <- W*U
      call CPMM(U,0,0.d0,.FALSE.,W,0,0.d0,.FALSE.,R)
      call FlushCP(W)
      call reduc(U,R)
!      call ReplaceVwithW(U,R)
      call FlushCP(R)

      deallocate(cols)

      end subroutine CPMatNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatNormA_noweights(U,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes CP-tensor U. S -> R inversion + matrix multiply with U is 
! replaced by solving for normalized U

      implicit none
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP) :: W,W2,S,Un
      integer, intent(in)  :: nals
      integer, allocatable :: cols(:)
      logical, parameter   :: reduceW=.TRUE.
      character(len=64) :: tag
      integer :: i,rk,redstat,solstat

      rk=SIZE(U%coef)

!     Compute the tensor of vector self-products
      W=CPMatSelfProdTot(U,U,.TRUE.)
      call CPMatrixTranspose(W)

!     (Optional) reduction of W
      IF (reduceW) THEN
         W2=RandomCP(W,rk)
         tag='CPMatNormA(): W2<-W'
         redstat=ALS_reduce(W2,W,nals,tag)
         call ReplaceVwithW(W,W2)
      ENDIF

      allocate(cols(SIZE(U%nbas)))
      cols(:)=1

!     Initial guess the square root, then run the square root iteration
      S=RandomCP(rk,W%rows,cols,W%sym)
      call SquareRootGuess(W,S,rk,10*nals)
      call CPsqrt(W,S,10*nals,nals)
      call FlushCP(W)

!     Convert sqrt(W) from vector repn to diagonal matrix, store as W
      W=CPVtoDiagMatrix(S,S%sym)

!     Solve W*Un^T = U^T for Un^T
      call CPMatrixTranspose(U)
      Un=CopyCP(U)
      tag='CPMatNormA(): solve W*Un^T=U^T'
      solstat=ALS_solve(W,Un,U,nals,0,0.d0,tag)

!     Transpose Un^T and replace U <- Un      
      call CPMatrixTranspose(Un)
      call ReplaceVwithW(U,Un)

      deallocate(cols)

      end subroutine CPMatNormA_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatNormA_weights(U,Wt,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes CP-tensor U. S -> R inversion + matrix multiply with U is 
! replaced by solving for normalized U

      implicit none
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP), INTENT(IN) :: Wt
      TYPE (CP) :: W,W2,S,Un
      integer, intent(in)  :: nals
      integer, allocatable :: cols(:)
      logical, parameter   :: reduceW=.TRUE.
      character(len=64) :: tag
      integer :: i,rk,redstat,solstat

      rk=SIZE(U%coef)

!     Compute the tensor of vector self-products
      W=CPMatSelfProdTot(U,U,.TRUE.)
      call CPMatrixTranspose(W)

!     (Optional) reduction of W
      IF (reduceW) THEN
         W2=RandomCP(W,rk)
         tag='CPMatNormA(): W2<-W'
!         redstat=ALS_reduce(W2,W,nals,tag)
         redstat=ALS_reduce(W2,W,Wt,nals,tag)
!!! TEST
!         call CPEntrywiseCompare(W2,W)
!         call ShowRank1Diff(W2,W)
!!!
         call ReplaceVwithW(W,W2)
      ENDIF

!      write(*,*) 'Here is W (CPMatNormA(), b4 CPsqrt())'
!      call W%print()

      allocate(cols(SIZE(U%nbas)))
      cols(:)=1

!     Initial guess the square root, then run the square root iteration
      S=RandomCP(rk,W%rows,cols,W%sym)
      call SquareRootGuess(W,S,rk,10*nals)
      call CPsqrt(W,S,Wt,10*nals,nals)
      call FlushCP(W)

!     Convert sqrt(W) from vector repn to diagonal matrix, store as W
      W=CPVtoDiagMatrix(S,S%sym)

!     Solve W*Un^T = U^T for Un^T
      call CPMatrixTranspose(U)
      Un=CopyCP(U)
      tag='CPMatNormA(): solve W*Un^T=U^T'
      solstat=ALS_solve(W,Un,U,Wt,nals,0,0.d0,tag)

!     Transpose Un^T and replace U <- Un      
      call CPMatrixTranspose(Un)
!!! TEST
!!!
      call ReplaceVwithW(U,Un)

      deallocate(cols)

      end subroutine CPMatNormA_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPsqrt_noweights(F,G,nitn,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes square root of CP-vector F, and stores in G
! G is provided as an initial guess

      implicit none
      TYPE (CP), INTENT(IN)    :: F
      TYPE (CP), INTENT(INOUT) :: G
      TYPE (CP) :: G1,Gdiag,X
      integer, intent(in) :: nitn,nals
      real*8, parameter :: thresh=8.d-8
      real*8, parameter :: delthresh=1.d-5
      integer :: itn,lowmem,oscct,redstat,solstat
      logical :: testmode=.TRUE.
      real*8  :: fmag,oldconver,conver,del
      character(len=64) :: frmt,tag

      IF (testmode) THEN
         frmt='(3X,A,I3,A,ES16.8)'
         write(*,*)
         write(*,*) "CPsqrt() iterations..."
         fmag=sqrt(abs(PRODVV(F)))
         call CPMM_vec(G,0,0.d0,G,0,0.d0,X) ! X = G**2
         oldconver=calc_FmG(X,F)/fmag
         write(*,frmt) 'itn ',0,': ||F - G**2||/||F|| = ',oldconver
         call FlushCP(X)
         IF (abs(oldconver).lt.thresh) THEN
            write(*,*) 'CPsqrt() iterations converged!'
            RETURN
         ENDIF
      ENDIF

      lowmem=1
      oscct=0

      DO itn=1,nitn
         G1=CopyCP(G)

!        Generate diagonal matrix repn of G
         Gdiag=CPVtoDiagMatrix(G,G%sym)

!        Solve for new G
!         call LinSolver_alg(Gdiag,G,F,nals,0,0.d0,lowmem)
         write(tag,'(A,I0)') 'CPsqrt (solve), itn = ',itn
         solstat=ALS_solve(Gdiag,G,F,nals,0,0.d0,tag)

!        Update: G1 <- (G1 + G)/2
         call SUMVECVEC(G1,0.5d0,G,0.5d0)

!        Reduce rank of new G
!         call reduc_ALS(G1,G,nals)
         write(tag,'(A,I0)') 'CPsqrt (redn), itn = ',itn
         redstat=ALS_reduce(G,G1,nals,tag)

!        Clean up
         call FlushCP(G1)
         call FlushCP(Gdiag)

         IF (testmode) THEN
            frmt='(3X,A,I3,A,ES16.8,A,ES16.8)'
            call CPMM_vec(G,0,0.d0,G,0,0.d0,X) ! X = G**2
            conver=calc_FmG(X,F)/fmag
            del=abs(conver-oldconver)/abs(conver)
            IF (conver.gt.oldconver) oscct=oscct+1
            write(*,frmt) 'itn ',itn,': ||F - G**2||/||F|| = ',conver,&
                          '; del = ',del
            call FlushCP(X)
            IF (abs(conver).lt.thresh) THEN 
               write(*,*) 'CPsqrt() iterations converged!'
               EXIT
            ENDIF
            IF (del.lt.delthresh) THEN 
               write(*,*) 'CPsqrt(): exit due to small residual change'
               EXIT
            ENDIF
            IF (oscct.eq.2) THEN
               write(*,*) 'CPsqrt(): exit due to oscillating residual'
               EXIT
            ENDIF
            oldconver=conver
         ENDIF
      ENDDO

      end subroutine CPsqrt_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPsqrt_weights(F,G,Wt,nitn,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes square root of CP-vector F, and stores in G, weighted version
! G is provided as an initial guess

      implicit none
      TYPE (CP), INTENT(IN)    :: Wt
      TYPE (CP), INTENT(INOUT) :: G,F
      TYPE (CP) :: G1,Gdiag,X
      integer, intent(in) :: nitn,nals
      real*8, parameter :: thresh=8.d-8
      real*8, parameter :: delthresh=1.d-5
      integer :: itn,lowmem,oscct,redstat,solstat
      logical :: testmode=.TRUE.
      real*8  :: fmag,oldconver,conver,del
      character(len=64) :: frmt,tag

      IF (testmode) THEN
         frmt='(3X,A,I3,A,ES16.8)'
         write(*,*)
         write(*,*) "CPsqrt() iterations..."
         fmag=sqrt(abs(PRODVV(F)))
         call CPMM_vec(G,0,0.d0,G,0,0.d0,X) ! X = G**2
         oldconver=calc_FmG(X,F)/fmag
         write(*,frmt) 'itn ',0,': ||F - G**2||/||F|| = ',oldconver
         call FlushCP(X)
         IF (abs(oldconver).lt.thresh) THEN
            write(*,*) 'CPsqrt() iterations converged!'
            RETURN
         ENDIF
      ENDIF

      lowmem=1
      oscct=0

      DO itn=1,nitn
         G1=CopyCP(G)

!        Generate diagonal matrix repn of G
         Gdiag=CPVtoDiagMatrix(G,G%sym)

!        Solve for new G
!         call LinSolver_alg(Gdiag,G,F,nals,0,0.d0,lowmem)
         write(tag,'(A,I0)') 'CPsqrt (solve), itn = ',itn
         solstat=ALS_solve(Gdiag,G,F,Wt,nals,0,0.d0,tag)

!        Update: G1 <- (G1 + G)/2
         call SUMVECVEC(G1,0.5d0,G,0.5d0)

!        Reduce rank of new G
!         call reduc_ALS(G1,G,nals)
          write(tag,'(A,I0)') 'CPsqrt (redn), itn = ',itn
          redstat=ALS_reduce(G,G1,Wt,nals,tag)
!!! TEST
!          call CPEntrywiseCompare(G,G1)
!          call ShowRank1Diff(G,G1)
!!!

!        Clean up
         call FlushCP(G1)
         call FlushCP(Gdiag)

         IF (testmode) THEN
            frmt='(3X,A,I3,A,ES16.8,A,ES16.8)'
            call CPMM_vec(G,0,0.d0,G,0,0.d0,X) ! X = G**2
            conver=calc_FmG(X,F)/fmag
            del=abs(conver-oldconver)/abs(conver)
            IF (conver.gt.oldconver) oscct=oscct+1
            write(*,frmt) 'itn ',itn,': ||F - G**2||/||F|| = ',conver,&
                          '; del = ',del
            call FlushCP(X)
            IF (abs(conver).lt.thresh) THEN 
               write(*,*) 'CPsqrt() iterations converged!'
               EXIT
            ENDIF
            IF (del.lt.delthresh) THEN 
               write(*,*) 'CPsqrt(): exit due to small residual change'
               EXIT
            ENDIF
            IF (oscct.eq.2) THEN
               write(*,*) 'CPsqrt(): exit due to oscillating residual'
               EXIT
            ENDIF
            oldconver=conver
         ENDIF
      ENDDO

      end subroutine CPsqrt_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPreciprocal(F,G,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes reciprocal of CP-vector F, and stores in G
! G is provided as an initial guess

      implicit none
      TYPE (CP), INTENT(IN)    :: F
      TYPE (CP), INTENT(INOUT) :: G
      TYPE (CP) :: I,Fdiag
      real*8 :: imag
      integer, intent(in) :: nals
      integer :: lowmem,j
!!! TEST
      TYPE (CP) :: Z,ZZ
      logical   :: testmode=.TRUE.
!!! END TEST

      lowmem=1

!     Generate diagonal matrix repn of F
      Fdiag=CPVtoDiagMatrix(F,F%sym)

!     Identity matrix
      I=NewCP(1,F%rows,F%cols,F%sym)
      I%base(:,:)=1.d0
      I%coef(:)=1.d0

!     Solve: diag(F)*G = I for G
      call LinSolver_alg(Fdiag,G,I,nals,0,0.d0,lowmem)
      call FlushCP(Fdiag)

!!! TEST
      IF (testmode) THEN
         call CPMM_vec(F,0,0.d0,G,0,0.d0,Z)  ! Z = F*G
         ZZ=CopyCP(I)
         call reduc_SR1(Z,ZZ,10*nals) ! ZZ <-ALS- Z (reduce to rank-1)
         call NORMBASE(ZZ)
         call DistributeCoef(ZZ)
         imag=sqrt(abs(PRODVV(I)))
         write(*,*)
         write(*,*) "CPreciprocal(): F*G -> rank-1"
         write(*,*)
         call PrintCPvec(ZZ)
         write(*,*)
         write(*,*) 'CP reciprocal(): ||(F*G)_R - I||/||I|| = ',&
         calc_FmG(Z,I)/imag
         write(*,*) 'CP reciprocal(): ||(F*G)_1 - I||/||I|| = ',&
         calc_FmG(ZZ,I)/imag

         call FlushCP(Z)
         call FlushCP(ZZ)
      ENDIF
!!! END TEST

      call FlushCP(I)

      end subroutine CPreciprocal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RecipSqrtRank1(F,G,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes rank-1 approx to G = F**(-1/2) via ALS and simple operations

      implicit none
      TYPE (CP), INTENT(IN)  :: F
      TYPE (CP), INTENT(OUT) :: G
      integer, intent(in) :: nals
      integer :: j

      G=RandomCP(F,1)

      IF (nals.lt.1) RETURN

!     Reduce F to rank-1 as G
      call reduc_SR1(F,G,nals)

!     Calc inverse square root of base, coefficient
      DO j=1,SIZE(G%base,1)
         G%base(j,1)=1.d0/sqrt(abs(G%base(j,1)))
      ENDDO
      G%coef(1)=1.d0/sqrt(abs(G%coef(1)))

      end subroutine RecipSqrtRank1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReciprocalGuessM(F,G,rk,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes guess to inverse of diagonal matrix stored as a vector.
! Reduce to rank-1. Invert. Add 'noise' to reach the proper rank. Done.
! This version receives and returns a diagonal matrix

      implicit none
      TYPE (CP), INTENT(IN)  :: F
      TYPE (CP), INTENT(OUT) :: G
      TYPE (CP) :: Fd,Gd
      integer, intent(in) :: rk,nals

      Fd=ExtractDiagfromCPMatrix(F,.FALSE.)
      call ReciprocalGuess(Fd,Gd,rk,nals)
      call FlushCP(Fd)
      G=CPVtoDiagMatrix(Gd,F%sym)
      call FlushCP(Gd)

      end subroutine ReciprocalGuessM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReciprocalGuess(F,G,rk,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes guess to inverse of diagonal matrix stored as a vector.
! Reduce to rank-1. Invert. Add 'noise' to reach the proper rank. Done.

      implicit none
      TYPE (CP), INTENT(IN)  :: F
      TYPE (CP), INTENT(OUT) :: G
      integer, intent(in) :: rk,nals
      integer :: j

      IF (rk.lt.1) &
         call AbortWithError("ReciprocalGuess(): rk must be > 0")

      G=NewCP(1,F%rows,F%cols,F%sym)
      G%base(:,:)=1.d0
      G%coef(:)=1.d0

      IF (nals.gt.0) THEN
!        Reduce to rank-1
         call reduc_SR1(F,G,nals)

!        Invert base, coefficient
         DO j=1,SIZE(G%base,1)
            G%base(j,1)=1.d0/G%base(j,1)
         ENDDO
         G%coef(1)=1.d0/G%coef(1)

      ENDIF
      call NORMBASE(G)

!     Add small "noise" terms up to rank rk
      call AugmentVWithRandom(G,rk)

      end subroutine ReciprocalGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SquareRootGuessM(F,G,rk,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes guess to square root of diagonal matrix.
! Reduce to rank-1. Square root. Add 'noise' to reach the proper rank.
! This version receives and returns a diagonal matrix

      implicit none
      TYPE (CP), INTENT(IN)  :: F
      TYPE (CP), INTENT(OUT) :: G
      TYPE (CP) :: Fd,Gd
      integer, intent(in) :: rk,nals

      Fd=ExtractDiagfromCPMatrix(F,.FALSE.)
      call SquareRootGuess(Fd,Gd,rk,nals)
      call FlushCP(Fd)
      G=CPVtoDiagMatrix(Gd,F%sym)
      call FlushCP(Gd)

      end subroutine SquareRootGuessM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SquareRootGuess(F,G,rk,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes guess to square root of diagonal matrix stored as a vector.
! Reduce to rank-1. Square root. Add 'noise' to reach the proper rank.

      implicit none
      TYPE (CP), INTENT(IN)  :: F
      TYPE (CP), INTENT(OUT) :: G
      integer, intent(in) :: rk,nals
      integer :: j

      IF (rk.lt.1) &
         call AbortWithError("SquareRootGuess(): rk must be > 0")

      G=NewCP(1,F%rows,F%cols,F%sym)
      G%base(:,:)=1.d0
      G%coef(:)=1.d0

      IF (nals.gt.0) THEN
!        Reduce to rank-1
         call reduc_SR1(F,G,nals)

!        Square root of base, coefficient
         DO j=1,SIZE(G%base,1)
            G%base(j,1)=sqrt(abs(G%base(j,1)))
         ENDDO
         G%coef(1)=sqrt(abs(G%coef(1)))

      ENDIF
      call NORMBASE(G)

!     Add small "noise" terms up to rank rk
      call AugmentVWithRandom(G,rk)

      end subroutine SquareRootGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowNorm(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints normalization values for all the vectors in U

      implicit none
      TYPE (CP), INTENT(IN)  :: U
      TYPE (CP) :: W
      integer, allocatable :: indx(:)
      integer :: i,j,ndof
      real*8  :: norm
      character(len=64) :: frmt

      ndof=SIZE(U%cols)

      write(*,*)
      allocate(indx(ndof))
      indx(:)=U%cols(:)

      write(frmt,'(A,I0,A)') '(A,I8,A,',ndof,'(I2,X),A,f16.8)'
      write(*,*)
      DO i=1,PRODUCT(U%cols)
         call NextIndex(indx,U%cols)
         W=ExtractCPvec(U,indx,.TRUE.)
         norm=sqrt(abs(PRODVV(W)))
         write(*,frmt) "Vec: ",i," (",(indx(j),j=1,ndof),&
                       "); ||F|| = ",norm
         call FlushCP(W)
      ENDDO
      deallocate(indx)

      end subroutine ShowNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowRQ(H,U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints normalization values for all the vectors in U

      implicit none
      TYPE (CP), INTENT(IN) :: H,U
      TYPE (CP) :: W
      integer, allocatable :: indx(:)
      integer :: i,j,ndof
      real*8  :: norm,rq
      character(len=64) :: frmt

      ndof=SIZE(U%cols)

      write(*,*)
      allocate(indx(ndof))
      indx(:)=U%cols(:)

      write(frmt,'(A,I0,A)') '(A,I8,A,',ndof,'(I2,X),2(A,f16.8))'
      write(*,*)

      DO i=1,PRODUCT(U%cols)
         call NextIndex(indx,U%cols)
         W=ExtractCPvec(U,indx,.TRUE.)
         norm=sqrt(abs(PRODVV(W)))
         call NORMALIZE(W)
         rq=RayleighQuotient(W,H)
         write(*,frmt) "Vec: ",i," (",(indx(j),j=1,ndof),&
         "); RQ = ",rq,'; ||F|| = ',norm
         call FlushCP(W)
      ENDDO
      deallocate(indx)

      end subroutine ShowRQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowRank1UTU(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints rank-1 approx. to UTU as a quick check of the overlap matrix

      implicit none
      TYPE (CP), INTENT(IN)  :: U
      TYPE (CP) :: W,X,I
      integer, parameter :: nals=100

      call CPMM(U,.TRUE.,U,.FALSE.,W)
      X=RandomCP(W,1)
      call reduc_SR1(W,X,nals)
      call PrintCPmat(X)
      call FlushCP(X)

      call ShowDistanceFromIdentity(W)
      call FlushCP(W)

      end subroutine ShowRank1UTU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowRank1R(U,Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints rank-1 approx. of R from Q^T*U=R

      implicit none
      TYPE (CP), INTENT(IN) :: U,Q
      TYPE (CP) :: W,X
      integer, parameter :: nals=100

      call CPMM(Q,.TRUE.,U,.FALSE.,W)
      X=RandomCP(W,1)
      call reduc_SR1(W,X,nals)
      call FlushCP(W)
      call PrintCPmat(X)
      call FlushCP(X)

      end subroutine ShowRank1R

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowRank1Diff(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints rank-1 approx. of F-G

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      TYPE (CP) :: W,X
      integer, parameter :: nals=100

      W=CopyCP(F)
      call SUMVECVEC(W,1.d0,G,-1.d0)
      X=RandomCP(W,1)
      call reduc_SR1(W,X,nals)
      call FlushCP(W)
      call PrintCPmat(X)
      call FlushCP(X)

      end subroutine ShowRank1Diff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowDistanceFromIdentity(U)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes distance from identity matrix

      implicit none
      TYPE (CP), INTENT(IN) :: U
      TYPE (CP) :: I,X
      real*8 :: imag
      character(len=64) :: frmt

      frmt='(A,ES14.8)'

      I=IdentityCPMatrix(U%rows,U%cols,U%sym)
      imag=sqrt(abs(PRODVV(I)))
      write(*,frmt) 'Distance from Identity: ||U_R - I||/||I|| = ',&
      calc_FmG(U,I)/imag

!      X=RandomCP(U,1)
!      call reduc_SR1(U,X,100)
!      write(*,*) 'Distance from Identity: ||U_1 - I||/||I|| = ',&
!      calc_FmG(X,I)/imag

!      call FlushCP(X)
      call FlushCP(I)

      end subroutine ShowDistanceFromIdentity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module CPMATH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
