!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module CP2CP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP vectors with CP-in-rank representation

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE BLOCKUTILS
      USE LINSOLVER
      USE ALSOO
      USE ALSDRVR
      USE CPMATH

      TYPE CPR
         TYPE (CP) :: coef
         TYPE (CP), ALLOCATABLE :: base(:)
         INTEGER, ALLOCATABLE :: nbas(:),ibas(:),fbas(:)
         INTEGER, ALLOCATABLE :: rows(:),cols(:)
         CONTAINS
            PROCEDURE :: D => GetCPRndof
            PROCEDURE :: sD => GetCPRsndof
            PROCEDURE :: M => GetCPRrows
            PROCEDURE :: N => GetCPRcols
            PROCEDURE :: print => PrintCPR
            PROCEDURE :: flush => FlushCPR 
      END TYPE CPR

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPRsamesrk(rows,cols,srows,scols,crk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CPR type with all base and coef s-ranks equal

      implicit none
      TYPE (CPR) :: v
      INTEGER, INTENT(IN) :: crk
      INTEGER, INTENT(IN) :: rows(:),cols(:),srows(:),scols(:)
      INTEGER, ALLOCATABLE :: srk(:)
      INTEGER :: i,ndof,nrdim

      ndof=SIZE(rows)

!     Check for correct sizes and consistency in dimensions. This is
!     also done in NewCPRgen, but must be duplicated here to avoid a
!     segfault if the rows, cols array sizes differ
      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

      IF (SIZE(cols).ne.ndof) THEN
         write(*,*) 'Error: size of "cols" array differs from ndof!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

!     Calc size of srk array (nrdim)
      nrdim=0
      DO i=1,ndof
         nrdim=nrdim+rows(i)*cols(i)
      ENDDO

!     Set base s-ranks same as coefficient rank (crk)
      ALLOCATE(srk(nrdim))
      srk(:)=crk
      v=NewCPRgen(rows,cols,srows,scols,srk,crk)
      DEALLOCATE(srk)

      end function NewCPRsamesrk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPRgen(rows,cols,srows,scols,srk,crk,ran) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CPR type

      implicit none
      TYPE (CPR) :: v
      INTEGER, INTENT(IN) :: crk
      INTEGER, INTENT(IN) :: rows(:),cols(:),srows(:),scols(:),srk(:)
      LOGICAL, OPTIONAL   :: ran
      LOGICAL :: rcp
      INTEGER :: ndof,nrdim,srdim,i

      ndof=SIZE(rows)
      srdim=SIZE(srk)
      IF (present(ran)) THEN
         rcp=ran
      ELSE
         rcp=.FALSE.
      ENDIF

!     Check for correct sizes and consistency in dimensions
      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

      IF (SIZE(cols).ne.ndof) THEN
         write(*,*) 'Error: size of "cols" array differs from ndof!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

      DO i=1,ndof
         IF (rows(i).lt.1) THEN
            write(*,*) 'Error: number of rows must be at least 1!'
            call AbortWithError('Error in NewCPR()')
         ENDIF

         IF (cols(i).lt.1) THEN
            write(*,*) 'Error: number of cols must be at least 1!'
            call AbortWithError('Error in NewCPR()')
         ENDIF
      ENDDO

!     Allocate arrays
      ALLOCATE(v%nbas(ndof),v%ibas(ndof),v%fbas(ndof))
      ALLOCATE(v%rows(ndof),v%cols(ndof))
      v%rows(:)=rows(:)
      v%cols(:)=cols(:)

      nrdim=0
      DO i=1,ndof
         v%ibas(i)=nrdim+1
         v%nbas(i)=rows(i)*cols(i)
         nrdim=nrdim+v%nbas(i)
         v%fbas(i)=nrdim
      ENDDO

!     Error check s-rank of coefs and base
      IF (crk.lt.1) THEN
         write(*,*) 'Error: crk (',crk,') must be at least 1!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

      IF (srdim.ne.nrdim) THEN
         write(*,*) 'Error: srk array length (',srdim,') must equal',&
         ' size of base (',nrdim,')!'
         call AbortWithError('Error in NewCPR()')
      ENDIF

!     Initialize coefficients (random without shift to keep positive)
      v%coef=NewCP(crk,srows,scols)
      IF (rcp) THEN
         call random_number(v%coef%coef(:))
         DO i=1,crk
            call random_number(v%coef%base(:,i))
         ENDDO
      ENDIF

!     Initialize base
      ALLOCATE(v%base(nrdim))
      DO i=1,nrdim
         IF (rcp) THEN
            v%base(i)=RandomCP(srk(i),srows,scols)
         ELSE
            v%base(i)=NewCP(srk(i),srows,scols)
         ENDIF
      ENDDO

      end function NewCPRgen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPRref(w,crk,srk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CPR-matrix v, using sizes in w (including the rank, if not
! passed as an optional argument)

      implicit none
      TYPE (CPR) :: v
      TYPE (CPR), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: crk
      INTEGER, INTENT(IN), OPTIONAL :: srk(:)
      INTEGER, ALLOCATABLE :: srks(:)
      INTEGER :: i,nrdim

!     Use the s-rank structure provided
      IF (present(crk)) THEN
         IF (present(srk)) THEN ! All ranks set to crk
            v=NewCPRsamesrk(w%rows,w%cols,w%coef%rows,w%coef%cols,crk)
         ELSE ! crk and srk set to input vals
            v=NewCPRgen(w%rows,w%cols,w%coef%rows,w%coef%cols,srk,crk)
         ENDIF
!     Use the s-rank structure in w
      ELSE
         nrdim=SIZE(w%base)
         ALLOCATE(srks(nrdim))
         DO i=1,nrdim
            srks(i)=w%base(i)%R()
         ENDDO
         v=NewCPRgen(w%rows,w%cols,w%coef%rows,w%coef%cols,&
                     srks,w%coef%R())
         DEALLOCATE(srks)
      ENDIF

      end function NewCPRref

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushCPR(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Frees memory from CPR type

      implicit none
      CLASS (CPR) :: v

      call FlushCP(v%coef)
      IF (ALLOCATED(v%base)) DEALLOCATE(v%base)
      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)
      IF (ALLOCATED(v%ibas)) DEALLOCATE(v%ibas)
      IF (ALLOCATED(v%fbas)) DEALLOCATE(v%fbas)
      IF (ALLOCATED(v%rows)) DEALLOCATE(v%rows)
      IF (ALLOCATED(v%cols)) DEALLOCATE(v%cols)

      end subroutine FlushCPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCPRndof(v) result(ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CPR), INTENT(IN) :: v
      integer :: ndof

      IF (ALLOCATED(v%nbas)) THEN
         ndof=SIZE(v%nbas)
      ELSE
         ndof=0
      ENDIF

      end function GetCPRndof

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCPRsndof(v) result(ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of s-modes of v

      implicit none
      CLASS (CPR), INTENT(IN) :: v
      integer :: ndof

      IF (ALLOCATED(v%coef%nbas)) THEN
         ndof=SIZE(v%coef%nbas)
      ELSE
         ndof=0
      ENDIF

      end function GetCPRsndof

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCPRrows(v,d) result(rows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CPR), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: rows

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) &
         'GetCPRrows(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%rows)) THEN
         rows=v%rows(d)
      ELSE
         rows=0
      ENDIF

      end function GetCPRrows

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCPRcols(v,d) result(cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CPR), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: cols

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) &
         'GetCPRcols(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%cols)) THEN
         cols=v%cols(d)
      ELSE
         cols=0
      ENDIF

      end function GetCPRcols

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPR(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CPR vector v

      implicit none
      CLASS (CPR) :: v
      integer :: i,d

      write(*,*) 'CPR: CP-vec of coefs:'
      write(*,*)
      call v%coef%printvec()
      DO d=1,v%D()
         DO i=v%ibas(d),v%fbas(d)
            write(*,'(2(A,I0/))') &
            'CPR: CP-vec of base: d = ',d,'; i = ',i
            call v%base(i)%printvec()
         ENDDO
      ENDDO

      end subroutine PrintCPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPR_calcPk(v,w,k) result(P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes mode k inner product matrix of two CPR-vecs (no rank red'n)

      implicit none
      TYPE (CPR), INTENT(IN) :: v,w
      TYPE (CP) :: P,T
      integer, intent(in) :: k
      integer, allocatable :: rows(:),cols(:)
      integer :: i,d,sd

      d=MAX(v%D(),w%D()) !!! d must be same for v,w
      IF (k.lt.1 .or. k.gt.d) THEN
         write(*,*) 'k = ',k,', but must be in [1,',d,']'
         call AbortWithError('CPR_calcPk(): mode k out of range')
      ENDIF

!     Create P as a zero matrix with the max number of modes of v,w,
!     where rows,cols of 'excess' modes are set to 1
      sd=MAX(v%sD(),w%sD()) !!! sd can be different for v,w
      ALLOCATE(rows(sd),cols(sd))
      rows(:)=1
      cols(:)=1
      rows(:v%sD())=v%coef%rows(:)
      cols(:w%sD())=w%coef%rows(:)
      P=NewCP(1,rows,cols)
      P%coef(1)=0.d0
      DEALLOCATE(rows,cols)

!     Calc v * w^T for each basis item and accumulate the sum
      DO i=v%ibas(k),v%fbas(k)
         call CPMM(v%base(i),.FALSE.,w%base(i),.TRUE.,T)
         call SUMVECVEC(P,1.d0,T,1.d0)
         call FlushCP(T)
      ENDDO

      end function CPR_calcPk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPR_calcPtot(v,w,kskip,mulvc,mulwc,rk) result(P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes inner product matrix of two CPR-vecs, all modes
! ALS rank reduction is applied to keep memory under control
! Skips the kskip-th DOF (set kskip=0 to include all DOF)

      implicit none
      TYPE (CPR), INTENT(IN) :: v,w
      TYPE (CP) :: P,T,T2
      integer, intent(in) :: kskip,rk
      logical, intent(in) :: mulvc,mulwc
      integer, parameter  :: nals=50
      logical, parameter  :: reduceT=.FALSE.
      integer, allocatable :: rows(:),cols(:)
      character(len=64) :: tag
      integer :: k,d,sd,redstat

      IF (.not. CPRCHECKNBAS(v,w)) THEN
         write(*,*) 'v and w must have the same # of rows and cols'
         call AbortWithError('CPR_calcPtot(): v,w size mismatch')
      ENDIF

!     Create P as a matrix of ones with the max number of modes of v,w,
!     where rows,cols of 'excess' modes are set to 1
      d=MAX(v%D(),w%D())
      sd=MAX(v%sD(),w%sD()) !!! sd can be different for v,w

      ALLOCATE(rows(sd),cols(sd))
      rows(:)=1
      cols(:)=1
      rows(:v%sD())=v%coef%rows(:)
      cols(:w%sD())=w%coef%rows(:)
      P=NewCP(1,rows,cols)
      P%coef(1)=1.d0
      P%base(:,1)=1.d0
      DEALLOCATE(rows,cols)

!     Loop over modes, generate (CP-format) P_k for each, and accumulate
!     the product over k in P
      DO k=1,d
         IF (k.eq.kskip) CYCLE

!        Compute Pk for mode i
         T=CPR_calcPk(v,w,k)

!        Optional rank-reduction T <-rename- T2 <-ALS- T
         IF (reduceT) THEN
            T2=RandomCP(T,rk)
            write(tag,'(A,I0)') 'CPR_calcPtot(): T2<-T, mode ',k
            redstat=ALS_reduce(T2,T,nals,tag)
            call ReplaceVwithW(T,T2)
         ENDIF

!        Hadamard product of P,T
         T2=CPHadamardProd(P,T)
         call FlushCP(T)

!        Reduce rank: P <-ALS- T2 (P is initial guess)
         call AugmentVWithRandom(P,rk)
         write(tag,'(A,I0)') 'CPR_calcPtot(): P<-T2, mode ',k
         redstat=ALS_reduce(P,T2,nals,tag)
         call FlushCP(T2)
      ENDDO

!     Premultiply by diagonal of w coefficients if option is set
      IF (mulvc) THEN
         T=CPVtoDiagMatrix(v%coef)
         call CPMM(T,.FALSE.,P,.FALSE.,T2)
         call FlushCP(T)
         write(tag,'(A,I0)') 'CPR_calcPtot(): P<-T2, coefs (v)'
         redstat=ALS_reduce(P,T2,nals,tag)
         call FlushCP(T2)
      ENDIF

!     Postmultiply by diagonal of w coefficients if option is set
      IF (mulwc) THEN
         T=CPVtoDiagMatrix(w%coef)
         call CPMM(P,.FALSE.,T,.FALSE.,T2)
         call FlushCP(T)
         write(tag,'(A,I0)') 'CPR_calcPtot(): P<-T2, coefs (w)'
         redstat=ALS_reduce(P,T2,nals,tag)
         call FlushCP(T2)
      ENDIF

      end function CPR_calcPtot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPR_NORMBASE1(v,k,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Normalizes base of mode k of CPR vector v. Coefs of v are replaced by
! those required to normalize mode k

      implicit none
      TYPE (CPR), INTENT(INOUT) :: v
      TYPE (CP) :: T,T2
      integer, intent(in) :: k,nals
      logical, parameter  :: reduceT=.TRUE.
      character(len=64) :: tag
      integer :: i,n,stat

      IF (k.lt.1 .or. k.gt.v%D()) THEN
         write(*,*) 'k = ',k,', but must be in [1,',v%D(),']'
         call AbortWithError('CPR_NORMBASE1(): mode k out of range')
      ENDIF

!     Compute sum-of-squares over basis elements
      T=NewCP(v%coef,1)
      T%coef(1)=0.d0
      DO i=v%ibas(k),v%fbas(k)
         T2=CPHadamardProd(v%base(i),v%base(i))
         call SUMVECVEC(T,1.d0,T2,1.d0)
         call FlushCP(T2)
      ENDDO

!     (Optional) reduction of T
      IF (reduceT) THEN
         T2=RandomCP(T,v%coef%R())
         tag='CPR_NORMBASE1(): T2<-T'
         stat=ALS_reduce(T2,T,nals,tag)
         call ReplaceVwithW(T,T2)
      ENDIF

!     Run square root itn on T to get new coefs
      call CPsqrt(T,v%coef,10*nals,nals)
      call FlushCP(T)

!     Convert coefs from vector repn to diagonal matrix, store as T
      T=CPVtoDiagMatrix(v%coef)

!     Solve the linear system to get the normalized v%base(i) for all i
      DO i=v%ibas(k),v%fbas(k)
         call CPMatrixTranspose(v%base(i))
         T2=CopyCP(v%base(i))
         tag='CPMatNormA(): solve T*(v%base(i)_norm)^T=(v%base(i))^T'
         stat=ALS_solve(T,v%base(i),T2,nals,0,0.d0,tag)
!        Transpose v%base(i) back to original repn
         call CPMatrixTranspose(v%base(i))
         call FlushCP(T2)
      ENDDO

      call FlushCP(T)

      end subroutine CPR_NORMBASE1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPR_ALS(v,w,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank v <- w , where v,w are CPR vectors

      implicit none
      TYPE (CPR), INTENT(INOUT) :: v
      TYPE (CPR), INTENT(IN) :: w
      TYPE (CP) :: B,P,rhs
      integer, intent(in) :: nals
      integer, parameter :: rk=80
      integer :: itn,d,i,stat
      character(len=64) :: tag

!     Error checking

      DO itn=1,nals

         DO d=1,v%D()

!           Compute B = PI_k <v,v>, P = PI_k <v,w>; for k != d
            B=CPR_calcPtot(v,v,d,.FALSE.,.FALSE.,rk)
            P=CPR_calcPtot(v,w,d,.FALSE.,.TRUE.,rk) ! Weighted by w coef

!!! Hey, I want to reuse the work done by linear solver as concerns B,
!!! as it is same for ALL rhs. Possible with existing code? For now just
!!! redo the work until I make it functional
            DO i=w%ibas(d),w%fbas(d)
!              Compute rhs
               call CPMM(P,.FALSE.,w%base(i),.FALSE.,rhs)

!              Solve linear system via call to CP linear solver
               write(tag,'(3(A,I0))') &
               'CPR_ALS(): solve linear system: itn = ',itn,'; d = ',d,&
               '; i = ',i
               stat=ALS_solve(B,v%base(i),rhs,nals,0,0.d0,tag)

               call FlushCP(rhs)
            ENDDO ! basis components (i)

!           Normalization
!!! 3 steps here:
!!! 1) Compute CP-vec of diagonal entries of B(d) (has N*S^2 terms) -
!!!    these (after possible rank-redn) are the normalization coefs
!!! 2) Solve the linear system for each basis component and the vec from
!!!    step 1 to get the normalized basis component CP-vecs
!!! 3) replace the old basis CP-vecs with the new ones
            call CPR_NORMBASE1(v,d,nals)

!!! Test convergence here?

         ENDDO ! modes (d)

      ENDDO ! ALS iterations (itn)

      end subroutine CPR_ALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPRCHECKNBAS(v1,v2) result (same)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CPR-format vectors to make sure they are the same
! This checks the n-basis structure, not the s-basis structure

      implicit none
      TYPE (CPR), INTENT(IN) :: v1,v2
      LOGICAL :: same
      INTEGER :: i

      same=.TRUE.

      IF (v1%D().ne.v2%D()) THEN
         same=.FALSE.
      ELSE
         DO i=1,v1%D()
            IF ((v1%rows(i).ne.v2%rows(i)) .or. &
                (v1%cols(i).ne.v2%cols(i))) THEN
               same=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      end function CPRCHECKNBAS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module CP2CP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
