!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
!!!
      USE TARGETEDSTATES
!!!
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

      subroutine GuessPsi(il,im,evalsND,Q,H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (Hamiltonian), INTENT(IN)  :: H
      TYPE (MLtree), INTENT(IN)       :: ML
      TYPE (CP), ALLOCATABLE, INTENT(OUT) :: Q(:)
      integer, intent(in)  :: il,im
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
!!!
!      call GetStatesinWindow(nbloc,evalsND,qns,evals1D,nbas,0.d0)
!!!
      call sortDPeigvals(nbloc,evalsND,qns,evals1D,nbas)
!      call structuredDPeigvals(nbloc,evalsND,qns,evals1D,nbas)

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

      function GuessWeights(il,im,T,H,ML) result(W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Guesses Boltzmann weights for separable H energies at temperature T,
! for weighting rank-reduction of CP matrices with multiple states. If 
! temperature is negative then identity matrix is guessed.

      implicit none
      TYPE (Hamiltonian), INTENT(IN)  :: H
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP) :: W
      real*8, intent(in)   :: T
      real*8, parameter    :: kb=0.69503476 ! cm^-1/K
      integer, intent(in)  :: il,im
      integer, allocatable :: nbas(:)
      logical, allocatable :: sym(:)
      integer :: i,j,mi,nsubm,mstart
      real *8 :: t1,t2,val

      IF (.NOT. GUESS_SETUP) call InitializeGuessModule()

      call CPU_TIME(t1)

!     Set parameters
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      mstart=ML%modstart(il,im) ! start index of sub-modes of 'im'
      ALLOCATE(nbas(nsubm),sym(nsubm))
      nbas(1:nsubm)=ML%gdim(il-1,mstart:mstart+nsubm-1)
      sym(:)=.FALSE.

      W=IdentityCPMatrix(nbas,nbas,sym)

      IF (T.ge.0.d0) THEN
!        Use the eigenvalues of the sub-modes to compute the weights
!        If T=0, force weights of 1 for g.s. and 0 for all other states
         DO i=1,nsubm
            mi=mstart+i-1
            DO j=1,nbas(i)
               IF (T.gt.0.d0) THEN
                  val=exp((H%eig(il-1,mi)%evals(1)-&
                           H%eig(il-1,mi)%evals(j))/(kb*T))
               ELSEIF (j.eq.1) THEN
                  val=1.d0
               ELSE
                  val=0.d0
               ENDIF
               W%base(W%ibas(i)+(j-1)*(W%rows(i)+1),1)=val   
            ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(nbas,sym)

      call CPU_TIME(t2)
      guess_time=guess_time+t2-t1

      end function GuessWeights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine structuredDPeigvals(nbloc,evalsND,qns,evals1D,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates all vecs in the direct-product representation (for testing)

      implicit none
      integer, allocatable, intent(out) :: qns(:,:)
      real*8, intent(inout):: evalsND(:)
      integer, intent(in)  :: nbas(:)
      real*8, intent(in)   :: evals1D(:,:)
      integer, allocatable :: indx(:)
      integer, intent(in)  :: nbloc
      integer :: i,j,ndof

      ndof=SIZE(nbas)

      IF (nbloc.ne.PRODUCT(nbas)) THEN
         write(*,'(2(A,I0),A)') 'nbloc (',nbloc,&
         ') must equal size of DP basis (',PRODUCT(nbas),')'
      ENDIF

      ALLOCATE(qns(nbloc,ndof),indx(ndof))
      indx(:)=nbas(:)

!     Run through DP-basis via calls to NextIndex and assemble the guess
!     eigenvalues
      DO i=1,nbloc
         call NextIndex(indx,nbas)
         qns(i,:)=indx(:)
         evalsND(i)=0.d0
         DO j=1,ndof
            evalsND(i)=evalsND(i)+evals1D(j,indx(j))
         ENDDO
      ENDDO

      DEALLOCATE(indx)

      end subroutine structuredDPeigvals

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

      subroutine BuildProdFunctions(Q,nbas,qns)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills block with products of eigenfunctions

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      integer, intent(in) :: qns(:,:),nbas(:)
      integer :: i,j,jmod,nbloc,ndof,gstart

      nbloc=SIZE(Q)
      ndof=SIZE(nbas)

      DO i=1,nbloc
         Q(i)=NewCP(1,nbas)
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
      TYPE (CP), INTENT(INOUT) :: Q(:)
      real*8, allocatable :: svals(:),M(:,:),U(:,:),VT(:,:)
      integer, intent(in) :: nbas(:)
      integer :: i,j,n,nbloc,ndof,gi,gf,rind

      nbloc=SIZE(Q)
      ndof=SIZE(nbas)

!     Initialize the columns of Q
      DO i=1,nbloc
         Q(i)=NewCP(1,nbas)
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
