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
      USE GPUINTERTWINE

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

      subroutine BuildModeHamiltonian(il,im,Ham,HL,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs Hamiltonian in CP-format for mode 'im' in layer 'il'

      implicit none
      TYPE (CPpar) :: cpp
      TYPE (MLtree), INTENT(IN)   :: ML
      TYPE (HopList), INTENT(OUT) :: HL
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (CPvec), ALLOCATABLE :: opcp(:,:)
      TYPE (CPvec) :: H,Hnew
      integer, intent(in)  :: il,im
      integer, allocatable :: nbas(:),poplist(:,:),popct(:)
      integer :: i,j,k,l,trm,nsubm,mst,nbloc,plo,primop,jmode
      integer :: thedof,dofind,oldrank
      real*8  :: t1,t2
      real*8, parameter :: redtol=1.d-12
      logical :: usegpu=.TRUE.
      logical :: showFmG=.TRUE.

      IF (.NOT. INIT_SETUP) call InitializeInitModule()

      call CPU_TIME(t1)

!     Set parameters
      trm=GetModeHNr(il,im,Ham) ! term which applies to mode 'im'
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      mst=ML%modstart(il,im)    ! start index of sub-modes of 'im'
      nbloc=ML%gdim(il,im)      ! block size (# eigfxns for 'im')
      ALLOCATE(nbas(nsubm))
      nbas(1:nsubm)=ML%gdim(il-1,mst:mst+nsubm-1)

      IF (nsubm.gt.1) write(*,'(3X,A)') 'Building mode Hamiltonian...'

      usegpu=.FALSE.
      IF ((cpp%solver .seq. 'pALS') .and. cpp%lowmem.eq.4 .and. &
          (nsubm.gt.2 .or. (nsubm.eq.2 .and. &
          .not.(cpp%red2D .seq. 'SVD')))) usegpu=.TRUE.

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
               call NewCPvec(opcp(k,i),(/popct(i)/),1)
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
               call NewCPvec(opcp(k,i),(/popct(i)/),1)
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
      Ham%redterms=SIZE(opcp,1)

      IF (nsubm.gt.1) THEN
         write(*,'(/3X,A)') '*** ModeH memory usage ***'
         write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
               oldrank,' to ',SIZE(opcp,1),' by sorting'

       IF (cpp%ncycle.gt.0) THEN

!           Hamiltonian rank reduction requested
            IF (cpp%hrank.gt.0 .and. (cpp%hrank.lt.SIZE(opcp,1))) THEN

!              Build CP-matrix repn of H
               call GetHMatsCP(il,im,H,opcp,nbas,poplist,Ham,ML)

!              Reduce rank to get compressed rep'n of H:
!              GPU code compiled in single precision -> use CPU code
!              GPU code compiled in double precision -> use GPU code
!              ...since double precision is needed for THIS reduction.
               call CPU_TIME(t2)
               init_time=init_time+t2-t1
               oldrank=SIZE(H%coef)
               Ham%redterms=cpp%hrank
               call NORMBASE(H)
               call ordre(H)
               call SetReductionParameters(cpp%hrank,cpp%hnals,redtol,&
                                           showFmG,cpp%red2D,cpp%redND)
               call NewCPvec(Hnew,H%nbas,cpp%hrank)
               call GenCopyWtoV(Hnew,H,1,cpp%hrank,1,cpp%hrank)
#ifdef USE_DOUBLES
               IF (useGPU) THEN
                  call GPU_ALS(Hnew,H,cpp%hnals)
               ELSE
                  call reduc(Hnew,H)
               ENDIF
#else
               call reduc(Hnew,H)
#endif
               call ReplaceVwithW(H,Hnew)
               call CPU_TIME(t1)
               write(*,'(7X,2(A,I0),A)') &
                    'Hamiltonian rank reduced from ',oldrank,' to ',&
                    SIZE(H%coef),' by reduc()'
               write(*,'(7X,A,f12.6,A)') 'H memory (reduced)  : ',&
                    REAL(SIZE(H%base)+SIZE(H%coef))/REAL(2**27),' GB'

!              Generate/index the operator matrices inside type Ham
               call HMatsCP2List(H,Ham)
               call GetHopListReduced(Ham,HL,usegpu)
               call FlushCPvec(H)

!           No Hamiltonian reduction
            ELSE

!              Generate/index the operator matrices inside type Ham
               call GetHMatsList(il,im,opcp,nbas,poplist,Ham,ML)

!              Operator List for GPU code
               call GetHopList(Ham,HL,nsubm,usegpu)
            ENDIF

            call PrintHamMem(Ham,HL,usegpu)

         ELSE
!           Allocate a "dummy" CP-vector with the correct size
!           since the solver needs this for the memory check
            oldrank=SIZE(opcp,1)

!           Hamiltonian rank reduction requested (memory check)
            IF (cpp%hrank.gt.0 .and. cpp%hrank.lt.oldrank) THEN

               oldrank=MIN(oldrank,cpp%hrank)
               Ham%redterms=cpp%hrank
               call NewCPmat(H,nbas,oldrank,.TRUE.)
               write(*,'(7X,2(A,I0),A)') &
                     'Hamiltonian rank reduced from ',SIZE(opcp,1),&
                     ' to ',SIZE(H%coef),' by reduc()'
               write(*,'(7X,A,f12.6,A)') 'H memory (reduced)  : ',&
                     REAL(SIZE(H%base)+SIZE(H%coef))/REAL(2**27),' GB'
               
!              Generate/index the operator matrices inside type Ham
!!!            SKIP THESE NEXT 2 STEPS fOR NOW SINCE THEY REQUIRE
!!!            UPDATER TO TRANSFORM THE OPERATOR MATRICES, WHICH WE WANT
!!!            TO AVOID IN A MEMORY CHECK
!              call HMatsCP2List(H,Ham)
!              call GetHopListReduced(Ham,HL,usegpu)
               call FlushCPvec(H)

!           No Hamiltonian reduction (memory check)
            ELSE
!!!            NOTHING DONE HERE YET SINCE THIS STEP REQUIRES UPDATER TO
!!!            TRANSFORM THE OPERATOR MATRICES, WHICH WE WANT TO AVOID
!!!            IN A MEMORY CHECK

!              Generate/index the operator matrices inside type Ham
!              call GetHMatsList(il,im,opcp,nbas,poplist,Ham,ML)

!              Operator List for GPU code
!              call GetHopList(Ham,HL,nsubm,usegpu)
            ENDIF
         ENDIF
      ENDIF

      DEALLOCATE(opcp,poplist,popct,nbas)

      call CPU_TIME(t2)
      init_time=init_time+t2-t1

      end subroutine BuildModeHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintHamMem(Ham,HL,usegpu)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints memory used to store operator matrices in Hamiltonian and
! HopList structures

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (HopList), INTENT(IN)     :: HL
      logical, intent(in) :: usegpu
      real*8, parameter   :: GB=1.d0/(2**27)
      real*8  :: nelem(5),nelemtot
      integer :: i,optot

!     Count the number of elements in the operator matrices stored in
!     TYPE Hamiltonian
      nelem=0.d0
      optot=0
      IF (ALLOCATED(Ham%iops)) THEN
         DO i=1,SIZE(Ham%iops)
            nelem(1)=nelem(1)+REAL(SIZE(Ham%iops(i)%mat))
         ENDDO
         optot=optot+SIZE(Ham%iops)
      ENDIF
      IF (ALLOCATED(Ham%eops)) THEN
         DO i=1,SIZE(Ham%eops)
            nelem(2)=nelem(2)+REAL(SIZE(Ham%eops(i)%mat))
         ENDDO
         optot=optot+SIZE(Ham%eops)
      ENDIF
      DO i=1,SIZE(Ham%pops)
         nelem(3)=nelem(3)+REAL(SIZE(Ham%pops(i)%mat))
      ENDDO
      DO i=1,SIZE(Ham%cops)
         nelem(4)=nelem(4)+REAL(SIZE(Ham%cops(i)%mat))
      ENDDO

      optot=optot+SIZE(Ham%pops)+SIZE(Ham%cops)
      nelemtot=nelem(1)+nelem(2)+nelem(3)+nelem(4)

!     Print the memory requirement for operators in TYPE Hamiltonian
      IF (ALLOCATED(Ham%iops)) THEN
         write(*,'(7X,A,2X,I6,X,A,X,f12.6,X,A)') 'Cost of storing',&
              SIZE(Ham%iops),'identity  operators:',nelem(1)*GB,'GB'
      ENDIF
      IF (ALLOCATED(Ham%eops)) THEN
          write(*,'(7X,A,2X,I6,X,A,X,f12.6,X,A)') 'Cost of storing',&
              SIZE(Ham%eops),'eigen     operators:',nelem(2)*GB,'GB'
      ENDIF
      write(*,'(7X,A,2X,I6,X,A,X,f12.6,X,A)') 'Cost of storing',&
           SIZE(Ham%pops),'primitive operators:',nelem(3)*GB,'GB'
      write(*,'(7X,A,2X,I6,X,A,X,f12.6,X,A)') 'Cost of storing',&
           SIZE(Ham%cops),'compound  operators:',nelem(4)*GB,'GB'
      write(*,'(7X,A,2X,I6,X,A,X,f12.6,X,A)') 'TOTAL operators',&
                    optot,'requiring          :',nelemtot*GB,'GB'

!     Print the memory requirement for operators in TYPE HopList
      IF (useGPU) THEN
         write(*,'(7X,2(A,X),A)') '------------------',&
               '(required for GPU code)','-----------------'
         write(*,'(7X,A,X,I6,X,A,X,f12.6,X,A)') 'HopList contains',&
               SIZE(HL%mptr),'of these, adding   :',SIZE(HL%mat)*GB,'GB'
         nelemtot=nelemtot+REAL(SIZE(HL%mat))
         write(*,'(7X,A)') &
         '------------------------------------------------------------'
      ENDIF

!     Print the total cost of storing H
      write(*,'(7X,A,16X,A,X,f12.6,X,A)') &
            'HAMILTONIAN MEMORY REQUIRED',':',nelemtot*GB,'GB'

      end subroutine PrintHamMem

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

      IMPLICIT NONE
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: H(:,:)
      TYPE (CPvec), ALLOCATABLE :: T(:,:)
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
                  call CopyWtoV(T(ttrm,j),H(i,j))
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

      IMPLICIT NONE
      TYPE (CPvec), INTENT(INOUT) :: H(:),T(:)
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

      subroutine GetHMatsCP(il,im,H,opcp,nbas,poplist,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds H in CP-matrix form by multiplying and adding operator matrices
! The operator IDs are stored in 'poplist'
! The sum-of-products scheme for the operators is stored in 'opcp'

      implicit none
      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CPvec), INTENT(OUT) :: H
      TYPE (CPvec), INTENT(IN)  :: opcp(:,:)
      INTEGER, INTENT(IN) :: nbas(:),poplist(:,:)
      INTEGER, INTENT(IN) :: il,im
      REAL*8, ALLOCATABLE :: tvec(:),tmat(:,:),tmat1(:,:),tmat2(:,:)
      INTEGER :: i,j,k,l,m,gst,ntrm,nsubm,nbasop,nrkop,primop,thedof
      REAL*8, PARAMETER :: tol=1.d-12

!     Set parameters
      ntrm=SIZE(opcp,1)
      nsubm=SIZE(opcp,2)

!     Allocate the Hamiltonian array for symmetric matrices
      call NewCPmat(H,nbas,ntrm,.TRUE.)

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

!           Put the sum-of-products term into H
            call Mat2Vec(tvec,tmat,.TRUE.)

            H%base(gst+1:gst+H%nbas(j),i)=tvec
            DEALLOCATE(tmat,tvec)
            gst=gst+H%nbas(j)
         ENDDO
      ENDDO

      end subroutine GetHMatsCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HMatsCP2List(H,Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds H in CP-matrix form by multiplying and adding operator matrices
! The operator IDs are stored in 'poplist'
! The sum-of-products scheme for the operators is stored in 'opcp'

      IMPLICIT NONE
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (CPvec), INTENT(IN)  :: H
      REAL*8, ALLOCATABLE :: tmat(:,:)
      INTEGER :: i,j,gi,gf,ntrm,nsubm,npops,icop

!     Set parameters
      ntrm=SIZE(H%coef)
      nsubm=SIZE(H%nbas)
      npops=SIZE(Ham%pops)

!     Flush the old operators and allocate space for new ones
      IF (ALLOCATED(Ham%opindx)) DEALLOCATE(Ham%opindx)
      IF (ALLOCATED(Ham%opcoef)) DEALLOCATE(Ham%opcoef)
      IF (ALLOCATED(Ham%eops)) DEALLOCATE(Ham%eops)
      IF (ALLOCATED(Ham%cops)) DEALLOCATE(Ham%cops)
      IF (ALLOCATED(Ham%iops)) DEALLOCATE(Ham%iops)
      ALLOCATE(Ham%opcoef(ntrm),Ham%opindx(ntrm,nsubm))
      ALLOCATE(Ham%cops(ntrm*nsubm))

!     Put the matrices in H into the Ham%cops array. Since the
!     Hamiltonian was previously put into CP-format and reduced, ALL
!     component matrices are now compound operators
      icop=1
      DO i=1,ntrm
         Ham%opcoef(i)=H%coef(i)
         gi=1
         DO j=1,nsubm
            gf=gi+H%nbas(j)-1
!           The operator index must be offset by the number of primitive
!           operators since HVBaseProd() looks in either Ham%pops()
!           or Ham%cops() depending on the value of Ham%opindx()
            Ham%opindx(i,j)=npops+icop
            Ham%cops(icop)%dof=j

!           Convert the base of H to the right format
!           Expand from 'stored diagonal' format in H%base (from ALS)
            call Vec2SymMat(H%base(gi:gf,i),tmat)
!           Convert to 'symmetric packed' format in Ham%cops (for BLAS)
            call SymPackMat2Vec(Ham%cops(icop)%mat,tmat)

            DEALLOCATE(tmat)
            icop=icop+1
            gi=gf+1
         ENDDO
      ENDDO

      end subroutine HMatsCP2List

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetHMatsList(il,im,opcp,nbas,poplist,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds H in CP-matrix form by multiplying and adding operator matrices
! The operator IDs are stored in 'poplist'
! The sum-of-products scheme for the operators is stored in 'opcp'

      IMPLICIT NONE
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (OperMat), ALLOCATABLE :: cops(:)
      TYPE (CPvec), INTENT(IN)  :: opcp(:,:)
      INTEGER, INTENT(IN) :: nbas(:),poplist(:,:)
      INTEGER, INTENT(IN) :: il,im
      REAL*8, ALLOCATABLE :: tmat(:,:),tmat1(:,:),tmat2(:,:)
      INTEGER :: i,j,k,l,m,ntrm,nsubm,nbasop,nrkop,primop,thedof
      INTEGER :: nop,iop,npops,icop
      REAL*8, PARAMETER :: tol=1.d-12

!     Set parameters
      ntrm=SIZE(opcp,1)
      nsubm=SIZE(opcp,2)
      npops=SIZE(Ham%pops)

!     Flush the old operators and allocate space for new ones
      IF (ALLOCATED(Ham%opindx)) DEALLOCATE(Ham%opindx)
      IF (ALLOCATED(Ham%opcoef)) DEALLOCATE(Ham%opcoef)
      IF (ALLOCATED(Ham%eops)) DEALLOCATE(Ham%eops)
      IF (ALLOCATED(Ham%cops)) DEALLOCATE(Ham%cops)
      IF (ALLOCATED(Ham%iops)) DEALLOCATE(Ham%iops)
      ALLOCATE(Ham%opcoef(ntrm),Ham%opindx(ntrm,nsubm))
      ALLOCATE(Ham%eops(nsubm),Ham%iops(nsubm))
      ALLOCATE(cops(ntrm*nsubm))

      Ham%opcoef=1.d0

      icop=0
      DO i=1,ntrm
         DO j=1,nsubm
            nrkop=SIZE(opcp(i,j)%coef)
            nbasop=opcp(i,j)%nbas(1)

!           Determine if term is an identity/eigen/primitive operator
!           These cases do not require constructing an operator matrix
            IF (nrkop.eq.1) THEN
               nop=0
               DO l=1,nbasop
                  IF (nint(opcp(i,j)%base(l,1)).ne.0) THEN
                     nop=nop+1
                     iop=l
                  ENDIF
               ENDDO
!              No operator terms listed --> identity operator
               IF (nop.eq.0) THEN
                  Ham%opindx(i,j)=0
                  Ham%opcoef(i)=Ham%opcoef(i)*opcp(i,j)%coef(1)
                  CYCLE

!              One operator term --> eigen or primitive operator
               ELSEIF (nop.eq.1) THEN
                  Ham%opindx(i,j)=poplist(j,iop)
                  Ham%opcoef(i)=Ham%opcoef(i)*opcp(i,j)%coef(1)

!                 If this is an eigen-operator, get the matrix here
                  IF (poplist(j,iop).lt.0) THEN
                     thedof=ML%modstart(il,im)+abs(poplist(j,iop))-1
                     Ham%eops(abs(poplist(j,iop)))=GetEigenOperMat(&
                     abs(poplist(j,iop)),Ham%eig(il-1,thedof)%evals)
                  ENDIF
                  CYCLE
               ENDIF
            ENDIF

!           Build compound operators and store in the 'cops' array
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
!                    This should never occur here
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

!           Put the sum-of-products term into temporary cops array
            icop=icop+1
            call SymPackMat2Vec(cops(icop)%mat,tmat) 
            cops(icop)%dof=j

            Ham%opindx(i,j)=npops+icop
            DEALLOCATE(tmat)
         ENDDO
      ENDDO

!     Copy the compound operator matrices into the Hamiltonian type
      ALLOCATE(Ham%cops(icop))
      Ham%cops(:)=cops(1:icop)
      DEALLOCATE(cops)

!     Get the matrices for the identity operators
      DO j=1,nsubm
         Ham%iops(j)=GetIdentityOperMat(j,nbas(j))
      ENDDO

      end subroutine GetHMatsList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetHopList(Ham,HL,ndof,usegpu)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts all of the operator matrices stored in TYPE Hamiltonian and
! stores them in TYPE HopList, along with associated pointer arrays
! This is needed in order to assemble the operator matrices in a form 
! that can be passed to C++ since Fortran derived types with deferred 
! shape arrays are not interoperatable with C++

      IMPLICIT NONE
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (HopList), INTENT(OUT)    :: HL
      logical, intent(in)  :: usegpu
      integer, intent(in)  :: ndof
      integer, allocatable :: pptr(:,:),pind(:)
      integer :: i,j,k,l,ls,n,msize,nrk,npop,ncop,psum,opct
      logical :: found

!     Filling HopList is only needed for GPU code
      IF (.not.usegpu) RETURN

!     Initialize values
      nrk=SIZE(Ham%opcoef,1)
      npop=SIZE(Ham%pops)
      ncop=SIZE(Ham%cops)

!     Error checking
      IF (SIZE(Ham%iops).ne.ndof) THEN
         write(*,*) 'SIZE(Ham%iops): ',SIZE(Ham%iops),'; ndof: ',ndof,&
                    ' must match'
         call AbortWithError('Error in HopList()')
      ENDIF
      IF (SIZE(Ham%eops).ne.ndof) THEN
         write(*,*) 'SIZE(Ham%eops): ',SIZE(Ham%eops),'; ndof: ',ndof,&
                    ' must match'
         call AbortWithError('Error in GetHopList()')
      ENDIF

!     We need to copy the primitive operator matrices contained in the
!     mode Hamiltonian, which is a subset of the operator matrices
!     stored in Ham%pops. So here we identify the unique primitive
!     operators in this mode Hamiltonian
      ALLOCATE(pptr(nrk,ndof),pind(ndof))
      pptr(:,:)=0
      pind(:)=0
      DO k=1,nrk
         DO j=1,ndof
            IF (Ham%opindx(k,j).ge.1 .and. &
                Ham%opindx(k,j).le.npop) THEN
               found=.FALSE.
               DO i=1,pind(j)
                  IF (Ham%opindx(k,j).eq.pptr(i,j)) THEN
                     found=.TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (.not.found) THEN
                  pind(j)=pind(j)+1
                  pptr(pind(j),j)=Ham%opindx(k,j)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!     Count the total number of operators to be stored
      opct=2*ndof+ncop
      DO j=1,ndof
         opct=opct+pind(j)
      ENDDO

!     Fill HL%opid, which points to the location of each operator in the
!     CP-format Hamiltonian (as it will appear in HL)
      ALLOCATE(HL%opid(nrk*ndof))
      l=1
      DO k=1,nrk
         psum=0
         DO j=1,ndof
            IF (Ham%opindx(k,j).eq.0) THEN  ! Identity operator
               HL%opid(l)=j
            ELSEIF (Ham%opindx(k,j).lt.0) THEN  ! Eigen operator
               HL%opid(l)=ndof+abs(Ham%opindx(k,j))
            ELSEIF (Ham%opindx(k,j).gt.npop) THEN  ! Compound operator
               HL%opid(l)=2*ndof+Ham%opindx(k,j)-npop
            ELSE  ! Primitive operator
               found=.FALSE.
               DO i=1,pind(j)
                  IF (Ham%opindx(k,j).eq.pptr(i,j)) THEN
                     found=.TRUE.
                     HL%opid(l)=2*ndof+ncop+psum+i
                     EXIT
                  ENDIF
               ENDDO
               IF (.not.found) THEN
                  call AbortWithError("primitive operator not found")
               ENDIF
            ENDIF
            psum=psum+pind(j)
            l=l+1
         ENDDO
      ENDDO

!     Determine the size of the HL%mat array, which stores matrices for 
!     identity, eigen operators for all DOF, all compound operators, and
!     the unique primitive operators
      msize=0
      DO j=1,ndof  ! Identity and eigen operators
         msize=msize+SIZE(Ham%iops(j)%mat)+SIZE(Ham%eops(j)%mat)
      ENDDO
      DO j=1,SIZE(Ham%cops)  ! Compound operators
         msize=msize+SIZE(Ham%cops(j)%mat)
      ENDDO
      DO j=1,ndof  ! Primitive operators
         DO k=1,pind(j)
            msize=msize+SIZE(Ham%pops(pptr(k,j))%mat)
         ENDDO
      ENDDO

!     Build the HL%mat array and the list of pointers, HL%mptr, which
!     point to the starting index of each matrix
      ALLOCATE(HL%mat(msize),HL%mptr(opct))
      l=1
      ls=1
      DO j=1,ndof  ! Identity operators
         n=SIZE(Ham%iops(j)%mat)
         HL%mat(ls:ls+n-1)=Ham%iops(j)%mat(1:n)
         HL%mptr(l)=ls
         l=l+1
         ls=ls+n
      ENDDO
      DO j=1,ndof  ! Eigen operators
         n=SIZE(Ham%eops(j)%mat)
         HL%mat(ls:ls+n-1)=Ham%eops(j)%mat(1:n)
         HL%mptr(l)=ls
         l=l+1
         ls=ls+n
      ENDDO
      DO j=1,ncop  ! Compound operators
         n=SIZE(Ham%cops(j)%mat)
         HL%mat(ls:ls+n-1)=Ham%cops(j)%mat(1:n)
         HL%mptr(l)=ls
         l=l+1
         ls=ls+n
      ENDDO
      DO j=1,ndof  ! Primitive operators
         DO k=1,pind(j)
            n=SIZE(Ham%pops(pptr(k,j))%mat)
            HL%mat(ls:ls+n-1)=Ham%pops(pptr(k,j))%mat(1:n)
            HL%mptr(l)=ls
            l=l+1
            ls=ls+n
         ENDDO
      ENDDO

!     Finally, fill HopList with the coefficients
      ALLOCATE(HL%coef(nrk))
      HL%coef(:)=Ham%opcoef(:)

      DEALLOCATE(pptr,pind)

      end subroutine GetHopList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetHopListReduced(Ham,HL,usegpu)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts all of the operator matrices stored in TYPE Hamiltonian and
! stores them in TYPE HopList, along with associated pointer arrays
! This is needed in order to assemble the operator matrices in a form 
! that can be passed to C++ since Fortran derived types with deferred 
! shape arrays are not interoperatable with C++. This version only
! stores compound operators resulting from rank-reducing H in CP-format;
! primitive, eigen, and identity operators are not copied

      IMPLICIT NONE
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (HopList), INTENT(OUT)    :: HL
      logical, intent(in)  :: usegpu
      integer :: j,k,l,ls,n,msize,nrk,ndof,ncop,npops
      logical :: found

!     Filling HopList is only needed for GPU code
      IF (.not.usegpu) RETURN

!     Initialize values
      nrk=SIZE(Ham%opindx,1)
      ndof=SIZE(Ham%opindx,2)
      ncop=SIZE(Ham%cops)
      npops=SIZE(Ham%pops)

!     Fill HL%opid, which points to the location of each operator in the
!     CP-format Hamiltonian (as it will appear in HL)
      ALLOCATE(HL%opid(nrk*ndof))
      l=1
      DO k=1,nrk
         DO j=1,ndof
            HL%opid(l)=Ham%opindx(k,j)-npops
            l=l+1
         ENDDO
      ENDDO

!     Determine the size of the HL%mat array, which stores matrices for 
!     compound operators resulting from reducing the rank of the 
!     CP-format Hamiltonian
      msize=0
      DO j=1,SIZE(Ham%cops)
         msize=msize+SIZE(Ham%cops(j)%mat)
      ENDDO

!     Build the HL%mat array and the list of pointers, HL%mptr, which
!     point to the starting index of each matrix
      ALLOCATE(HL%mat(msize),HL%mptr(ncop))
      l=1
      ls=1
      DO j=1,ncop  ! Compound operators
         n=SIZE(Ham%cops(j)%mat)
         HL%mat(ls:ls+n-1)=Ham%cops(j)%mat(1:n)
         HL%mptr(l)=ls
         l=l+1
         ls=ls+n
      ENDDO

!     Finally, fill HopList with the coefficients
      ALLOCATE(HL%coef(nrk))
      HL%coef(:)=Ham%opcoef(:)

      end subroutine GetHopListReduced

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODEH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
