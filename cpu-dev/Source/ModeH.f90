!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODEH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module builds the CP-format mode Hamiltonian

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
      USE INPUTFIELDS
      USE MODECOMB
      USE MODVECVEC
      USE REDUCTION
!!!
      USE ALSDRVR
!!!

      implicit none
      real*8, private  :: init_time=0.d0
      logical, private :: INIT_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeInitModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      init_time = 0.d0
      INIT_SETUP = .TRUE.

      end subroutine InitializeInitModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeInitModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. INIT_SETUP) call InitializeInitModule()

      INIT_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total H initialization time       (s)',&
                             init_time

      end subroutine DisposeInitModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildModeHamiltonian(il,im,H,Ham,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs Hamiltonian in CP-format for mode 'im' in layer 'il'

      implicit none
      TYPE (CPpar) :: cpp
      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CP), ALLOCATABLE :: opcp(:,:)
      TYPE (CP), INTENT(OUT) :: H
!!!   TEST
      TYPE (CP) :: Hnew
!!!
      integer, intent(in)  :: il,im
      integer, allocatable :: nbas(:),poplist(:,:),popct(:)
      integer :: i,j,k,l,trm,nsubm,mst,nbloc,plo,primop,jmode
      integer :: gst,thedof,dofind,oldrank
      real*8  :: t1,t2,Hstor
      real*8, parameter :: redtol=1.d-12
      character*64 :: frmt
      logical      :: showFmG

      IF (.NOT. INIT_SETUP) call InitializeInitModule()

      call CPU_TIME(t1)

!     Set parameters
      showFmG=.FALSE.
      trm=GetModeHNr(il,im,Ham) ! term which applies to mode 'im'
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      mst=ML%modstart(il,im)    ! start index of sub-modes of 'im'
      nbloc=ML%gdim(il,im)      ! block size (# eigfxns for 'im')
      ALLOCATE(nbas(nsubm))
      nbas(1:nsubm)=ML%gdim(il-1,mst:mst+nsubm-1)

      IF (nsubm.gt.1) write(*,'(3X,A)') 'Building mode Hamiltonian...'

!     Get the list of primitive operator IDs for the sub-modes
      call GetPrimOpList(il,im,poplist,popct,Ham,ML)

!     Build the mode Hamiltonian as sums-of-products of operators,
!     indexed by their operator IDs and stored in CP-format as "opcp"
      ALLOCATE(opcp(Ham%nop(trm,il),nsubm))
      DO k=1,Ham%nop(trm,il)
         plo=Ham%mops(trm,il)%v(k) ! previous layer operator

!        (Pre-solved) single-mode operator
         IF (Ham%ndof(plo,il-1).eq.1) THEN

            DO i=1,nsubm
               thedof=Ham%dofs(plo,1,il-1)
               dofind=thedof-mst+1
               opcp(k,i)=NewCP(1,(/popct(i)/))
               opcp(k,i)%base=0.d0
               opcp(k,i)%coef=1.d0

!              Correct mode: add the (dummy) eigen-operator to opcp
               IF (dofind.eq.i) THEN
                  DO j=1,popct(i)
                     IF (poplist(i,j).eq.-i) opcp(k,i)%base(j,1)=1.d0
                  ENDDO
               ENDIF
            ENDDO

!        Multi-mode operator
         ELSE

!           Trace the previous layer operator to the bottom layer
!           ('plo' becomes the index value in the bottom layer)
            DO j=il-1,2,-1
               plo=Ham%mops(plo,j)%v(1)
            ENDDO

            DO i=1,nsubm
               opcp(k,i)=NewCP(1,(/popct(i)/))
               opcp(k,i)%base=0.d0
               opcp(k,i)%coef=1.d0
               IF (i.eq.1) opcp(k,i)%coef=Ham%facs(plo,1)

!              If the operator applies to the i-th sub-mode, include
               DO j=1,Ham%ndof(plo,1)
                  primop=Ham%ops(plo,j,1)
                  jmode=uppermodenr(il-1,1,Ham%pops(primop)%dof,ML)&
                        -mst+1
                  IF (jmode.eq.i) THEN
                     DO l=1,popct(i)
                        IF (primop.eq.poplist(i,l)) &
                           opcp(k,i)%base(l,1)=1.d0
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

!     Combine terms to reduce the rank of opcp "by hand"
      oldrank=SIZE(opcp,1)
      call condenseH(opcp)

      IF (nsubm.gt.1) THEN
         write(*,'(/3X,A)') '*** ModeH memory usage ***'
         write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
               oldrank,' to ',SIZE(opcp,1),' by sorting'

!        Calculate memory requirement
         Hstor=0.d0
         DO i=1,nsubm
            Hstor=Hstor+nbas(i)*(2*nbas(i)+1) !!! This is in SVD repn
         ENDDO
         write(*,'(7X,A,f12.6,A)') 'H memory (sorted)   : ',&
               Hstor*SIZE(opcp,1)/2**27,' GB'
      ENDIF

      IF (cpp%ncycle.gt.0) THEN

!        Now assemble H in matrix representation from the list in opcp
         call GetHMats(il,im,H,opcp,nbas,poplist,Ham,ML)

!        Additional reduction of H if desired
         call CPU_TIME(t2)
         init_time=init_time+t2-t1
         oldrank=SIZE(H%coef)

         IF (cpp%hrank.gt.0 .and. cpp%hrank.lt.oldrank) THEN
            call NORMBASE(H)
            call ordre(H)
!           Set the Hamiltonian reduction parameters first
            call SetReductionParameters(cpp%hrank,cpp%hnals,redtol,&
                                        showFmG,cpp%red2D,cpp%redND)
            Hnew=NewCP(cpp%hrank,H%rows,H%cols,H%sym)
            call GenCopyWtoV(Hnew,H,1,cpp%hrank,1,cpp%hrank)
            call reduc(Hnew,H)
            call ReplaceVwithW(H,Hnew)
         ENDIF
!!!

         IF (nsubm.gt.1 .and. SIZE(H%coef).lt.oldrank) THEN
            write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
                  oldrank,' to ',SIZE(H%coef),' by reduc()'
            write(*,'(7X,A,f12.6,A)') 'H memory (reduced)  : ',&
                  Hstor*SIZE(H%coef)/2**27,' GB'
         ENDIF

         call CPU_TIME(t1)

      ELSE
!        Allocate a "dummy" CP-vector with the correct size
!        since the solver needs this for the memory check
         oldrank=SIZE(opcp,1)
         IF (cpp%hrank.gt.0) oldrank=MIN(oldrank,cpp%hrank)
         H=NewCP(oldrank,nbas,.FALSE.)
         IF (nsubm.gt.1 .and. SIZE(H%coef).lt.SIZE(opcp,1)) THEN
            write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
                  SIZE(opcp,1),' to ',SIZE(H%coef),' by reduc()'
            write(*,'(7X,A,f12.6,A)') 'H memory (reduced)  : ',&
                  Hstor*SIZE(H%coef)/2**27,' GB'
         ENDIF

      ENDIF

!      write(*,*)
!      DO i=1,SIZE(opcp,1)
!         DO j=1,nsubm
!            write(*,('(A)')) 'Term DOF Rnk'
!            write(*,'(3I4)') i,j,SIZE(Hc(i,j)%coef)
!            call PrintCPvec(Hc(i,j))
!         ENDDO
!      ENDDO
!      call AbortWithError('Done generating the Hamiltonian of doom!!!')

      DEALLOCATE(poplist,popct,nbas,opcp)

      call CPU_TIME(t2)
      init_time=init_time+t2-t1

      end subroutine BuildModeHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPrimOpList(il,im,poplist,popct,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets list of primitive operators which apply to each sub-mode within
! a mode

      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      INTEGER, INTENT(IN) :: il,im
      INTEGER, ALLOCATABLE, INTENT(OUT) :: poplist(:,:),popct(:)
      INTEGER :: j,k,trm,nsubm,plo,primop,jmode

      trm=GetModeHNr(il,im,Ham) ! term which applies to mode 'im'
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'

!     Initialize the poplist to include one operator by default
!     (this corresponds to the diagonal eigenvalue matrix for the mode)
      ALLOCATE(poplist(nsubm,SIZE(Ham%pops)+1),popct(nsubm))
      poplist=0
      DO k=1,nsubm
         poplist(k,1)=-k
      ENDDO
      popct=1

!     Loop over the terms in a mode operator
      DO k=1,Ham%nop(trm,il)
         plo=Ham%mops(trm,il)%v(k) ! previous layer operator

!        Multi-mode operator: see if any primitive operators are to
!        be added to the list
         IF (Ham%ndof(plo,il-1).gt.1) THEN

!           Trace the previous layer operator to the bottom layer
!           ('plo' becomes the index value in the bottom layer)
            DO j=il-1,2,-1
               plo=Ham%mops(plo,j)%v(1)
            ENDDO

            DO j=1,Ham%ndof(plo,1)
               primop=Ham%ops(plo,j,1)
               jmode=uppermodenr(il-1,1,Ham%pops(primop)%dof,ML)&
                     -ML%modstart(il,im)+1

!              If the primitive operator has not yet been encountered,
!              add to list
               IF (.NOT.(ANY(poplist(jmode,:).eq.primop))) THEN
                  popct(jmode)=popct(jmode)+1
                  poplist(jmode,popct(jmode))=primop
               ENDIF  
            ENDDO
         ENDIF
      ENDDO

      end subroutine GetPrimOpList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine condenseH(H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Does "by-hand" reduction of H by combining terms with common factors

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: H(:,:)
      TYPE (CP), ALLOCATABLE :: T(:,:)
      INTEGER :: i,j,htrm,ttrm,nsubm,iadd
      LOGICAL :: add
      INTEGER :: k,pass

!     Set parameters
      htrm=SIZE(H,1)
      nsubm=SIZE(H,2)
      pass=0

!     If there is only one term in H, no need to reduce!
      IF (htrm.eq.1) RETURN

      ALLOCATE(T(htrm,nsubm))

!     Main loop over condensing cycles
      DO
         ttrm=0
         DO i=1,htrm
            add=.FALSE.
            DO j=1,ttrm
!              Determine if terms can be summed-at-constant-rank
               call compareHT(H(i,:),T(j,:),add,iadd)
               IF (add) THEN
                  call SUMVECVEC(T(j,iadd),1.d0,H(i,iadd),1.d0)
                  EXIT
               ENDIF
            ENDDO

!           If term in H cannot be summed with any term in tmp,
!           increase the rank of tmp by adding the H-term to the end
            IF (.NOT.add) THEN
               ttrm=ttrm+1
               DO j=1,nsubm
                  T(ttrm,j)=CopyCP(H(i,j))
               ENDDO
            ENDIF
         ENDDO

!        Replace H <-- tmp
         DEALLOCATE(H)
         ALLOCATE(H(ttrm,nsubm))
         DO i=1,ttrm
            DO j=1,nsubm
               call ReplaceVwithW(H(i,j),T(i,j))
            ENDDO
         ENDDO

!        If the size of HT cannot be reduced further, exit
         IF (ttrm.ge.htrm) EXIT
         htrm=ttrm
         pass=pass+1
      ENDDO

      end subroutine condenseH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine compareHT(H,T,add,iadd)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares Hamiltonian terms by their primitive operator matrices
! (indexed in H%base and T%base) and coefficients, in order to determine
! if the terms can be condensed. The coefs of H and T may be modified.

      implicit none
      TYPE (CP), INTENT(INOUT) :: H(:),T(:)
      INTEGER, INTENT(OUT) :: iadd
      LOGICAL, INTENT(OUT) :: add
      LOGICAL, ALLOCATABLE :: sameb(:)
      INTEGER :: i,j,k,l,nsubm,nonm,inm(2),nrkH,nrkT,nbas
      LOGICAL :: allfound,found,bmatch,allcmatch,cmatch
      REAL*8  :: fac
      REAL*8, PARAMETER  :: tol=1.d-15

      nsubm=SIZE(H)

!     Start by assuming terms are different in base and coef
      ALLOCATE(sameb(nsubm))
      sameb=.FALSE.

!     Loop over sub-modes and compare the bases/coefs
      nonm=0
      DO i=1,nsubm
         nrkH=SIZE(H(i)%coef)
         nrkT=SIZE(T(i)%coef)
         nbas=H(i)%nbas(1)

!        If the ranks are the same, compare term-by-term
         IF (nrkH.eq.nrkT) THEN
!           See if each term in H is somewhere in T (not necessarily in
!           the same order as in H)
            allfound=.TRUE.
            allcmatch=.TRUE.
            DO j=1,nrkH
               found=.FALSE.
               cmatch=.FALSE.
!              See if a term in the H-base matches any in the T-base
               DO k=1,nrkT
                  bmatch=.TRUE.
                  DO l=1,nbas
                     IF (abs(H(i)%base(l,j)-T(i)%base(l,k)).gt.tol) THEN
                        bmatch=.FALSE.
                        EXIT
                     ENDIF
                  ENDDO
!                 If the base matches, compare the coefficients
                  IF (bmatch) THEN
                     found=.TRUE.
                     IF (abs(H(i)%coef(j)-T(i)%coef(k)).lt.tol) &
                        cmatch=.TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (.not.found) THEN
                  allfound=.FALSE.
                  allcmatch=.FALSE.
                  EXIT
               ENDIF
               IF (.not.cmatch) allcmatch=.FALSE.
            ENDDO
         ELSE
            allfound=.FALSE.
            allcmatch=.FALSE.
         ENDIF

!        If base or coef does not match, add to count. Also, record the
!        first TWO non-matches in case a reduction is made possible by a
!        coefficient swap
         IF (allfound) sameb(i)=.TRUE.
         IF (.not.(allfound .and. allcmatch)) THEN
            nonm=nonm+1
            IF (nonm.le.2) THEN
               inm(nonm)=i
            ENDIF
         ENDIF        
      ENDDO

!     Exactly one non-matching mode --> always sum
      IF (nonm.eq.1) THEN
         add=.TRUE.
         iadd=inm(1)

!     TWO non-matching modes --> sum only if a coefficient swap can
!     convert to a one-non-matching mode case

!     First non-matching mode has the same base
      ELSEIF ((nonm.eq.2) .and. sameb(inm(1)) .and. &
              (SIZE(H(inm(1))%coef).eq.1)) THEN

!           Coef swap in H
            fac=H(inm(1))%coef(1)
            H(inm(1))%coef=H(inm(1))%coef/fac
            H(inm(2))%coef=H(inm(2))%coef*fac

!           Coef swap in H produces a match: sum
            IF (abs(H(inm(1))%coef(1)-T(inm(1))%coef(1)).lt.tol) THEN
               add=.TRUE.
               iadd=inm(2)

!           Coef swap in H produces no match: do coef swap in T also
            ELSE
               fac=T(inm(1))%coef(1)
               T(inm(1))%coef=T(inm(1))%coef/fac
               T(inm(2))%coef=T(inm(2))%coef*fac

!              Coef swap in T produces a match: sum
               IF (abs(H(inm(1))%coef(1)-T(inm(1))%coef(1)).lt.tol) THEN
                  add=.TRUE.
                  iadd=inm(2)

!              Coef swap in T produces no match: do not sum
               ELSE
                  add=.FALSE.
               ENDIF
            ENDIF

!     Second non-matching mode has the same base
      ELSEIF ((nonm.eq.2) .and. sameb(inm(2)) .and. &
              (SIZE(H(inm(2))%coef).eq.1)) THEN

!           Coef swap in H
            fac=H(inm(2))%coef(1)
            H(inm(2))%coef=H(inm(2))%coef/fac
            H(inm(1))%coef=H(inm(1))%coef*fac

!           Coef swap in H produces a match: sum
            IF (abs(H(inm(2))%coef(1)-T(inm(2))%coef(1)).lt.tol) THEN
               add=.TRUE.
               iadd=inm(1)

!           Coef swap in H produces no match: do coef swap in T also
            ELSE
               fac=T(inm(2))%coef(1)
               T(inm(2))%coef=T(inm(2))%coef/fac
               T(inm(1))%coef=T(inm(1))%coef*fac

!              Coef swap in T produces a match: sum
               IF (abs(H(inm(2))%coef(1)-T(inm(2))%coef(1)).lt.tol) THEN
                  add=.TRUE.
                  iadd=inm(1)

!              Coef swap in T produces no match: do not sum
               ELSE
                  add=.FALSE.
               ENDIF
            ENDIF

!     Too many non-matching modes --> do not sum
      ELSE
         add=.FALSE.
         iadd=0
      ENDIF

      DEALLOCATE(sameb)

      end subroutine compareHT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetHMats(il,im,H,opcp,nbas,poplist,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds H in CP-matrix form by multiplying and adding operator matrices
! The operator IDs are stored in 'poplist'
! The sum-of-products scheme for the operators is stored in 'opcp'

      implicit none
      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CP), INTENT(OUT) :: H
      TYPE (CP), INTENT(IN)  :: opcp(:,:)
      INTEGER, INTENT(IN) :: nbas(:),poplist(:,:)
      INTEGER, INTENT(IN) :: il,im
      REAL*8, ALLOCATABLE :: tvec(:),tmat(:,:),tmat1(:,:),tmat2(:,:)
      INTEGER :: i,j,k,l,m,gst,ntrm,nsubm,nbasop,nrkop,primop,thedof
      REAL*8, PARAMETER :: tol=1.d-12
      TYPE (CP) :: Hi

!     Set parameters
      ntrm=SIZE(opcp,1)
      nsubm=SIZE(opcp,2)

!     Allocate the Hamiltonian array for non-symmetric matrices
      H=NewCP(ntrm,nbas,.FALSE.)
!      H=NewCP(ntrm,nbas,.TRUE.)

      H%coef=1.d0

      DO i=1,ntrm
         gst=0
         DO j=1,nsubm
            nrkop=SIZE(opcp(i,j)%coef)
            nbasop=opcp(i,j)%nbas(1)
            ALLOCATE(tmat(nbas(j),nbas(j)))
            tmat=0.d0
            
            DO k=1,nrkop
!              Start building the operator with an identity matrix
!              scaled by the operator coefficient
               call GetIdentityMatrix(tmat1,nbas(j),.FALSE.)
               DO l=1,nbas(j)
                  tmat1(l,l)=tmat1(l,l)*opcp(i,j)%coef(k)
               ENDDO

!              Multiply the intra-sub-mode product operators
               DO l=1,nbasop
!                 If base(l,k)=1, then operator is present
                  IF (abs(opcp(i,j)%base(l,k)-1.d0).lt.tol) THEN
                     primop=poplist(j,l)
!                    Pre-solved mode operator (list of eigenvalues)
                     IF (primop.lt.0) THEN
                        thedof=ML%modstart(il,im)-(1+primop)
                        call GetIdentityMatrix(tmat2,nbas(j),.FALSE.)
                        DO m=1,nbas(j)
                           tmat2(m,m)=Ham%eig(il-1,thedof)%evals(m)
                        ENDDO

!                    All other primitive operators
                     ELSE
                        call Vec2SymPackMat(Ham%pops(primop)%mat,tmat2)
                     ENDIF
!                    tmat1 = tmat1 * tmat2
                     call MatrixMult(tmat1,.FALSE.,tmat2,.FALSE.)
                     DEALLOCATE(tmat2)

                  ENDIF
               ENDDO
!              Add tmat1 to the sum
               tmat=tmat+tmat1
               DEALLOCATE(tmat1)

            ENDDO

!           Symmetrize to scrub non-Hermitian errors due to truncating
!           basis as one moves up the tree
            call SymmetrizeMat(tmat)

!           Put the sum-of-products term into H
            call Mat2Vec(tvec,tmat,.FALSE.)

            H%base(gst+1:gst+H%nbas(j),i)=tvec
            DEALLOCATE(tmat,tvec)
            gst=gst+H%nbas(j)
         ENDDO
      ENDDO

      end subroutine GetHMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1Inverse(H,Hi)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes approximate H^-1 by reducing H to rank-1 and inverting each
! coordinate Hamiltonian. On input, H should be given as the full matrix 
! rep'n, not the upper triangle, (sym=.FALSE.)

      implicit none
      TYPE (CP), INTENT(IN)  :: H
      TYPE (CP), INTENT(OUT) :: Hi
      real*8, allocatable  :: tmat(:,:),tmat1(:,:)
      real*8, allocatable  :: tvec(:)
      integer, allocatable :: nbas(:)
      integer :: j,ndof,gi,gf

      ndof=SIZE(H%nbas)
      ALLOCATE(nbas(ndof))

      DO j=1,ndof
         nbas(j)=NINT(sqrt(REAL(H%nbas(j))))
      ENDDO

!     Initial guess: reduce Hi <- H, with Hi rank-1
      call SetReductionParameters(1,30,1.d-12,.FALSE.,'SVD','SR1')
      call reduc(Hi,H)

!     Invert each little-h in Hi
      gi=1
      DO j=1,ndof
         gf=gi+Hi%nbas(j)-1
         call Vec2Mat(Hi%base(gi:gf,1),tmat,nbas(j),nbas(j))
         call MatrixPseudoinverse(tmat,tmat1)
!        Symmetrize to correct small numerical errors that might result
!        from ALS
         call SymmetrizeMat(tmat1)
         call Mat2Vec(tvec,tmat1,.FALSE.)
         Hi%base(gi:gf,1)=tvec(1:Hi%nbas(j)) 
         DEALLOCATE(tmat,tmat1,tvec)
         gi=gf+1
      ENDDO
      Hi%coef(1)=1.d0/Hi%coef(1)

      end subroutine GetRank1Inverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShiftHbyE(H,E,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Adds a term to H corresponding to a shift of energy E. If only the
! upper triangle of Hsym is stored, set sym to .TRUE.

      implicit none
      TYPE (CP), INTENT(INOUT) :: H
      TYPE (CP) :: I
      logical, intent(in)  :: sym
      real*8, intent(in)   :: E
      logical, allocatable :: symm(:)
      integer, allocatable :: nbas(:)
      real*8, allocatable  :: tmat(:,:),tvec(:)
      integer :: j,ndof,gi,gf

      ndof=SIZE(H%nbas)

!     Find the number of basis functions per DOF from H
      ALLOCATE(nbas(ndof),symm(ndof))
      DO j=1,ndof
         IF (sym) THEN
            nbas(j)=GetSymN(H%nbas(j))
         ELSE
            nbas(j)=NINT(sqrt(REAL(H%nbas(j))))
         ENDIF
      ENDDO
      symm(:)=sym

!     Get an identity matrix in CP format
      I=IdentityCPMatrix(nbas,nbas,symm)

!     Hshifted = H + E*I
      call SUMVECVEC(H,1.d0,I,E)
      call FlushCP(I)
      DEALLOCATE(nbas,symm)

      end subroutine ShiftHbyE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RepackHmats(T)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Takes H, the Hamiltonian in CP-format, and converts component matrices 
! into symmetric-packed format

      implicit none
      TYPE (CP), INTENT(INOUT) :: T
      REAL, ALLOCATABLE :: M(:,:),v(:)
      INTEGER :: i,j,j1,j2,gst,ndof,nrk

!     Set parameters
      nrk=SIZE(T%coef)
      ndof=SIZE(T%nbas)

      DO i=1,nrk
         gst=0
         DO j=1,ndof
            IF (j.gt.1) gst=gst+T%nbas(j-1)
            j1=gst+1
            j2=gst+T%nbas(j)

!           Unwrap matrix in H, which is in stored diagonal format
            call Vec2Mat(T%base(j1:j2,i),M)

!           Convert matrix to symmetric-packed format; store in H
            call SymPackMat2Vec(v,M)
            T%base(j1:j2,i)=v(:)

            DEALLOCATE(M,v)
         ENDDO
      ENDDO

      end subroutine RepackHmats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODEH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
