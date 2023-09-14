!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
      USE GENCONFIG

      implicit none
      real*8, private  :: guess_time=0.d0
      logical, private :: GUESS_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeGuessModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      guess_time = 0.d0
      GUESS_SETUP = .TRUE.

      end subroutine InitializeGuessModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeGuessModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. GUESS_SETUP) call InitializeGuessModule()

      GUESS_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total wave-function guess time    (s)',&
                             guess_time

      end subroutine DisposeGuessModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GuessPsi(il,im,trunctyp,evalsND,Q,H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (Hamiltonian), INTENT(IN)  :: H
      TYPE (MLtree), INTENT(IN)       :: ML
      TYPE (CPvec), ALLOCATABLE, INTENT(OUT) :: Q(:)
      integer, intent(in)  :: il,im,trunctyp
      integer, allocatable :: qns(:,:),nbas(:)
      real*8, allocatable  :: evals1D(:,:)
      real*8, allocatable, intent(out) :: evalsND(:)
      integer :: i,j,mi,nsubm,mstart,nbloc,maxbas
      real *8 :: t1,t2
      character*64 :: frmt

      IF (.NOT. GUESS_SETUP) call InitializeGuessModule()

      call CPU_TIME(t1)

!     Set parameters
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      mstart=ML%modstart(il,im) ! start index of sub-modes of 'im'
      nbloc=ML%gdim(il,im)      ! block size (# eigfxns for 'im')
      ALLOCATE(nbas(nsubm))
      nbas(1:nsubm)=ML%gdim(il-1,mstart:mstart+nsubm-1)
      maxbas=MAXVAL(nbas)    ! max # b.f. of sub-modes in mode 'im'

!     Allocate the block (Q) and eigenvalue arrays
      ALLOCATE(Q(nbloc),evalsND(nbloc))

!     Copy the eigenvalues of the sub-modes to evals1D
      ALLOCATE(evals1D(nsubm,maxbas))
      DO i=1,nsubm
         mi=mstart+i-1
         IF (nbas(i).ne.SIZE(H%eig(il-1,mi)%evals)) &
            call AbortWithError("Mismatch in block, eigval list sizes")
         evals1D(i,:nbas(i))=H%eig(il-1,mi)%evals(:)
      ENDDO

!     Get the list of sorted N-D eigenvalues and quantum numbers
      call sortDPeigvals(nbloc,evalsND,qns,evals1D,nbas)
      call TruncateByQN(il,im,trunctyp,evalsND,qns,H,ML)

      IF (nsubm.gt.1) THEN
         write(*,'(/3X,A/)') 'Initial guess product functions:'
         write(frmt,*) '(3X,',nsubm,'(I3,X),f19.12)'
         DO i=1,nbloc
            write(*,frmt) (qns(i,j)-1,j=1,nsubm),evalsND(i)
         ENDDO
         write(*,*)
      ENDIF

!     Build the N-D separable eigenfunctions
      call BuildProdFunctions(Q,nbas,qns)

      DEALLOCATE(qns,nbas,evals1D)

      call CPU_TIME(t2)
      guess_time=guess_time+t2-t1

      end subroutine GuessPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortDPeigvals(nbloc,evalsND,qns,evals1D,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts eigenvalues of a direct product separable wavefunction

      implicit none
      integer, allocatable, intent(out) :: qns(:,:)
      real*8, allocatable, intent(inout):: evalsND(:)
      integer, allocatable, intent(in)  :: nbas(:)
      real*8, allocatable, intent(in)   :: evals1D(:,:)
      integer, allocatable :: qnlist(:,:),tmpqn(:)
      real*8, allocatable  :: sums(:)
      integer, intent(in)  :: nbloc
      integer :: i,j,k,l,ndof,prodND,nguess,ii,jj
      integer :: minind(1),nallowed
      real*8  :: sumt,bignr
      logical :: allowed,add2end

      ndof=SIZE(nbas)
      bignr=1.0d99

!     Error if more states are requested than allowed by the basis
      prodND=1
      DO j=1,ndof
         prodND=prodND*nbas(j)
         IF (prodND.ge.nbloc) EXIT
      ENDDO
      IF (nbloc.gt.prodND) THEN
         write(*,*) 'Error: nbloc is too large for product basis size'
         call AbortWithError('Error in sortDPeigvals()')
      ENDIF

      ALLOCATE(qns(nbloc,ndof),tmpqn(ndof))
      ALLOCATE(qnlist(ndof*nbloc,ndof),sums(ndof*nbloc))

!     Initialize the list of q.n. list and the sums with the ground state
      nguess=1
      qnlist(1,:)=1
      sums(1)=0.d0
      DO j=1,ndof
         sums(1)=sums(1)+evals1D(j,1)
      ENDDO

!     Build the list of states
      DO i=1,nbloc

!        Get the state in the list with the lowest energy
         evalsND(i)=MINVAL(sums(1:nguess))
         minind=MINLOC(sums(1:nguess))
         qns(i,:)=qnlist(minind(1),:)

!        Find the allowed excitations out of the previously-found
!        lowest energy state and add to qnlist
         add2end=.FALSE.
         nallowed=0
         DO j=1,ndof
            sumt=evalsND(i)
            tmpqn=qns(i,:)
            call ExciteQN(j,evals1D(j,:),nbas(j),tmpqn,sumt)
            allowed=checkallowed(qnlist(1:nguess,:),tmpqn,minind(1),&
                    nbas)
            IF (allowed) THEN
               nallowed=nallowed+1
               IF (add2end) THEN ! Allowed new q.n.: add to end
                  nguess=nguess+1
                  qnlist(nguess,:)=tmpqn
                  sums(nguess)=sumt
               ELSE  ! First new allowed q.n.: replace old one
                  qnlist(minind(1),:)=tmpqn
                  sums(minind(1))=sumt
                  add2end=.TRUE.
               ENDIF
            ENDIF
         ENDDO
!        IF none of the excitations are allowed, remove from the list
         IF (nallowed.eq.0) THEN
            qnlist(minind(1),:)=nbas(:)+1
            sums(minind(1))=bignr
         ENDIF
      ENDDO

      DEALLOCATE(qnlist,sums,tmpqn)

      end subroutine sortDPeigvals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine showprodfxns(D,nmax)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates list of states to keep

      implicit none
      integer, intent(in) :: D,nmax
      integer, allocatable :: qn(:),qmax(:)
      integer :: i,j,oldsum,newsum
      logical :: success
      character*60 :: frmt

      ALLOCATE(qn(D),qmax(D))
      qn(:)=1
      qmax(:)=nmax

      write(frmt,*) '(I5,A,X,',D,'(I2,X))'
      i=1
      oldsum=SUM(qn)
      DO
         write(*,frmt) i,')',(qn(j)-1,j=1,D)
!        Generate the next function in the "canonical" order
         call nextconfig(qn,qmax,D,.FALSE.,success)
         newsum=SUM(qn)
         IF (.not.success) EXIT
         IF (newsum.gt.oldsum) write(*,*) 'nextqn'
         oldsum=newsum
         i=i+1
         IF (i.gt.99999) EXIT
      ENDDO

      call AbortWithError('Done')

      end subroutine showprodfxns

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TruncateByQN(il,im,trunctyp,evals,qns,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates list of states to keep for single-mode truncation layers,
! where states are kept by canonical quantum number order rather than by
! increasing energy

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      real*8, intent(inout)  :: evals(:)
      integer, intent(inout) :: qns(:,:)
      integer, intent(in) :: il,im,trunctyp
      integer, allocatable :: qnc(:),qnt(:),maxqn(:),qnf(:,:)
      integer, allocatable :: istate(:),pstate(:),dupd(:)
      real*8, allocatable  :: eigv(:)
      integer :: i,j,k,nmode,nsubm,mstart,nbloc,nblocprev,tmp
      integer :: idup,fdup,ndup,mst,mfi,nagn,nmmax,cstate
      logical :: nmcfg,success
      character*64 :: frmt

!     Set parameters
      mstart=ML%modstart(il,im) ! start index of sub-modes of 'im'
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      nbloc=ML%gdim(il,im)      ! block size (# eigfxns for 'im')
      nblocprev=ML%gdim(il-1,mstart)    ! # prev-layer basis functions
      mst=firstmode(il,1,im,ML) ! start index of bottom-layer modes
      mfi=lastmode(il,1,im,ML)  ! end index of bottom layer modes
      nagn=mfi-mst+1            ! number of bottom layer modes in 'im'

!     Exit if this is not a truncation layer or if truncating the bottom
!     layer functions (in which case the states are already in order) or
!     if truncating by energy (trunctyp=0)
      IF (nsubm.gt.1 .or. nbloc.eq.nblocprev .or. il.le.2 .or. &
          trunctyp.eq.0) RETURN

      write(*,'(3X,2(A,I0),A)') &
      'Truncating basis from ',nblocprev,' to ',nbloc,' functions'

!     Get energies for printing and for ordering duplicate states
      ALLOCATE(eigv(nblocprev))
      eigv(:)=Ham%eig(il-1,mstart)%evals(:)

!     Get the list of nblocprev functions computed in the previous layer
      ALLOCATE(qnf(nblocprev,nagn),istate(nblocprev),dupd(nblocprev))
      dupd(:)=0
      DO i=1,nblocprev
         istate(i)=i
         call GetFullAssignment(il-1,mstart,Ham,ML,&
              Ham%eig(il-1,mstart)%assgn(i,:),qnt)
         qnf(i,:)=qnt(:)
         DEALLOCATE(qnt)
      ENDDO

      ALLOCATE(maxqn(nagn))
      nmmax=0
      DO i=1,nagn
         maxqn(i)=MAXVAL(qnf(:,i))
         IF (maxqn(i).gt.1) nmmax=nmmax+1
      ENDDO

!     Truncation options
      IF (trunctyp.lt.0) THEN ! Truncate by sum of qns
         nmode=MIN(nagn,nmmax)
         nmcfg=.FALSE.
         write(*,'(4X,A/)') &
         'prioritizing lower sums of quantum numbers...'
      ELSE ! Truncate by # of active modes, then by sum of qns
         nmode=MIN(MIN(trunctyp,nagn),nmmax)
         nmcfg=.TRUE.
         write(*,'(4X,A,I0,A/)') &
         'prioritizing limiting nr. of excited modes to ',nmode,'...'
      ENDIF

!!! TEST
!      write(frmt,*) '(X,A,X,',nagn,'(I2,X))'
!      write(*,frmt) 'Max qns:',(maxqn(k)-1,k=1,nagn)
!      write(*,*) 'Looking for states...'
!!!

!     Generate list of functions in canonical order
      ALLOCATE(qnc(nagn),qnt(nagn))
      qnc(:)=1 ! canonical state to match
      cstate=1 ! current state index
      DO i=1,nblocprev

!        Match "canonical" state to previously computed states, keeping
!        track of duplicate matches
         ndup=0
         DO j=cstate,nblocprev
            IF (ALL(qnc(:).eq.qnf(j,:))) THEN

!!! TEST
!               write(frmt,*) '(I4,X,I4,A,X,',nagn,'(I2,X))'
!               write(*,frmt) i,j,') found ',(qnc(k)-1,k=1,nagn)
!!!
               ndup=ndup+1
               IF (j.gt.cstate) THEN ! swap cstate and j-th state
                  qnt(:)=qnf(cstate,:)
                  qnf(cstate,:)=qnf(j,:)
                  qnf(j,:)=qnt(:)
                  tmp=istate(cstate)
                  istate(cstate)=istate(j)
                  istate(j)=tmp
               ENDIF
               cstate=cstate+1
            ENDIF

         ENDDO

!!! TEST
!         IF (ndup.eq.0) THEN
!            write(frmt,*) '(I4,X,I4,A,X,',nagn,'(I2,X))'
!            write(*,frmt) i,j,') not found ',(qnc(k)-1,k=1,nagn)
!         ENDIF
!!!

!        Reorder duplicate states by increasing energy
         IF (ndup.gt.1) THEN
            idup=cstate-ndup
            fdup=cstate-1
            dupd(idup:fdup)=1
            DO j=idup,fdup-1
               DO k=j+1,fdup
                  IF (eigv(istate(j)).gt.eigv(istate(k))) THEN
                     tmp=istate(j)
                     istate(j)=istate(k)
                     istate(k)=tmp
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

!        Generate the next function in the "canonical" order
         call nextconfig(qnc,maxqn,nmode,nmcfg,success)
         IF (.not.success) EXIT
      ENDDO
      DEALLOCATE(qnc,qnt,maxqn)

!     Get priority ordering of original states. This step is necessary
!     since we still want to keep states in order of increasing energy, but
!     only retaining those with priorities not exceeding nbloc
      ALLOCATE(pstate(nblocprev))
      DO i=1,nblocprev
         pstate(istate(i))=i
      ENDDO

!     Print state info
      write(frmt,*) '(A,X,',nagn,'(I2,X),5X,A,14X,A,11X,A)'
      write(*,frmt) 'Mode:',(ML%resort(j),j=mst,mfi),'Energy',&
                    'E-E0','Priority'
      k=1
      evals(:)=0
      qns(:,:)=0
      DO i=1,nblocprev
          IF (pstate(i).le.nbloc) THEN
            evals(k)=eigv(i)
            qns(k,1)=i
            k=k+1
            write(frmt,*) '(I4,A,X,',nagn,&
                 '(I2,X),2(f19.12,X),I5)'
            write(*,frmt) i,')',(qnf(pstate(i),j)-1,j=1,nagn),&
                 eigv(i),eigv(i)-eigv(1),pstate(i)
         ELSE
            write(frmt,*) '(I4,A,X,',nagn,&
                 '(I2,X),2(f19.12,X),I5,X,A)'
            write(*,frmt) i,')',(qnf(pstate(i),j)-1,j=1,nagn),&
                 eigv(i),eigv(i)-eigv(1),pstate(i),'(discard)'
         ENDIF
      ENDDO
      write(*,*)

      DEALLOCATE(qnf,eigv,istate,pstate,dupd)

      end subroutine TruncateByQN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildProdFunctions(Q,nbas,qns)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills block with products of eigenfunctions

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      integer, intent(in) :: qns(:,:),nbas(:)
      integer :: i,j,jmod,nbloc,ndof,gstart

      nbloc=SIZE(Q)
      ndof=SIZE(nbas)

      DO i=1,nbloc
         call NewCPvec(Q(i),nbas,1)
         Q(i)%coef(1)=1.d0
         Q(i)%base(:,1)=5.d-16 ! Not exactly zero so ALS does not crash
         gstart=0
         DO j=1,ndof
            IF (j.gt.1) gstart=gstart+nbas(j-1)
            jmod=gstart+qns(i,j)
            Q(i)%base(jmod,1)=1.d0
         ENDDO
      ENDDO

      end subroutine BuildProdFunctions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RandomOrthogonalGuess(Q,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills block with random orthogonal vectors
! Warning: the orthogonality of the vectors is not checked

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      real*8, allocatable :: svals(:),M(:,:),U(:,:),VT(:,:)
      integer, intent(in) :: nbas(:)
      integer :: i,j,n,nbloc,ndof,gi,gf,rind

      nbloc=SIZE(Q)
      ndof=SIZE(nbas)

!     Initialize the columns of Q
      DO i=1,nbloc
         call NewCPvec(Q(i),nbas,1)
         Q(i)%coef(1)=1.d0
      ENDDO

!     Build a random orthogonal basis for each DOF
      gi=1
      DO j=1,ndof
         n=nbas(j)
         gf=gi+n-1

!        Make a random set of vectors
         ALLOCATE(M(n,n))
         DO i=1,n
            call random_number(M(:,i))
         ENDDO

!        Construct orthogonal basis with SVD
         call SolveWithSVD(svals,M,U,VT)

!        Put the columns of U into Q
         DO i=1,nbloc
            rind=GetRandomInteger(n)
            Q(i)%base(gi:gf,1)=U(:,rind)
         ENDDO

         DEALLOCATE(svals,M,U,VT)
         gi=gi+n
      ENDDO

      end subroutine RandomOrthogonalGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ExciteQN(dofnr,evals,nbas,qnlist,dsum)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Excites a quantum number in a multi-D list, and replaces the sum 
! resulting from the old value with the sum resulting from the new value

      implicit none
      integer, intent(in) :: dofnr,nbas
      integer, intent(inout) :: qnlist(:)
      real*8, intent(in) :: evals(:)
      real*8, intent(inout) :: dsum
      real*8 :: bignr,subtr,add

      bignr=1.0d99

!     If qnlist already exceeds the basis, do not excite
      IF (qnlist(dofnr).lt.nbas) THEN ! excite and update sum
         subtr=evals(qnlist(dofnr))
         add=evals(qnlist(dofnr)+1)
         dsum=dsum-subtr+add
         qnlist(dofnr)=qnlist(dofnr)+1
      ELSE
         dsum=bignr
         qnlist(dofnr)=nbas+1
      ENDIF

      end subroutine ExciteQN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function checkallowed(qnlist,tmpqn,iskip,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Return true if excitation of quantum number is allowed, i.e. if there
! are no states in the q.n. list with ALL q.n.s <= the ones in tmpqn

      implicit none
      integer, intent(in) :: qnlist(:,:),tmpqn(:),nbas(:)
      integer, intent(in) :: iskip
      integer :: k,l,ndof,nterms
      logical :: checkallowed

      ndof=SIZE(tmpqn)
      nterms=SIZE(qnlist,1)

      checkallowed=.TRUE.

      DO k=1,nterms
         IF (k.eq.iskip) CYCLE
         DO l=1,ndof
            IF (tmpqn(l).gt.nbas(l)) THEN
               checkallowed=.TRUE.
               EXIT
            ENDIF
            IF (tmpqn(l).lt.qnlist(k,l)) EXIT
            IF (l.eq.ndof) checkallowed=.FALSE.
         ENDDO
         IF (.not. checkallowed) EXIT
      ENDDO

      end function checkallowed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
