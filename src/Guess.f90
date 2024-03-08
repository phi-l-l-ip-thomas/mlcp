!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
      USE INPUTFIELDS
      USE TARGETEDSTATES

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

      subroutine GuessPsi(il,im,evalsND,Q,H,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (CPpar), INTENT(IN)        :: cpp
      TYPE (Hamiltonian), INTENT(IN)  :: H
      TYPE (MLtree), INTENT(IN)       :: ML
      TYPE (CP), ALLOCATABLE, INTENT(OUT) :: Q(:)
      integer, intent(in)  :: il,im
      integer, allocatable :: qns(:,:),nbas(:),nmode(:,:),nexci(:,:)
      integer, allocatable :: qnfull(:)
      real*8, allocatable  :: evals1D(:,:)
      real*8, allocatable, intent(out) :: evalsND(:)
      integer :: i,j,mi,nsubm,mstart,nbloc,maxbas,nagn
      integer :: prodND,nmsum,nesum,neibas,mst,mfi,constraint
      real *8 :: t1,t2,Etarget
      character*64 :: frmt,tag

      IF (.NOT. GUESS_SETUP) call InitializeGuessModule()

      call CPU_TIME(t1)

!     Set parameters
      nsubm=ML%modcomb(il,im)   ! number of sub-modes in mode 'im'
      mstart=ML%modstart(il,im) ! start index of sub-modes of 'im'
      nbloc=ML%gdim(il,im)      ! block size (# eigfxns for 'im')
      mst=firstmode(il,1,im,ML) ! index of 1st mode on bottom layer
      mfi=lastmode(il,1,im,ML)  ! index of last mode on bottom layer
      nagn=mfi-mst+1            ! number of bottom-layer modes
      ALLOCATE(nbas(nsubm))
      nbas(1:nsubm)=ML%gdim(il-1,mstart:mstart+nsubm-1)
      maxbas=MAXVAL(nbas)    ! max # b.f. of sub-modes in mode 'im'

!     Truncation layer: filter energies by criterion from input file
      IF (nsubm.eq.1 .and. nbloc.lt.maxbas) THEN
         constraint=cpp%truncation
         Etarget=0.d0
      ELSE
         constraint=0 ! Energy criterion
         Etarget=0.d0
      ENDIF

!     Error if more states are requested than allowed by the basis
      prodND=1
      DO j=1,nsubm
         prodND=prodND*nbas(j)
         IF (prodND.ge.nbloc) EXIT
      ENDDO
      IF (nbloc.gt.prodND) THEN
         write(*,*) 'Error: nbloc is too large for product basis size'
         call AbortWithError('Error in GuessPsi()')
      ENDIF

!     Allocate the block vectors (Q), guess energy, and guess quantum
!     number arrays
      ALLOCATE(Q(nbloc),evalsND(nbloc),qns(nbloc,nsubm))
      ALLOCATE(nmode(maxbas,nsubm),nexci(maxbas,nsubm))

!     Copy the eigenvalues of the sub-modes to evals1D
      ALLOCATE(evals1D(nsubm,maxbas))
      DO i=1,nsubm
         mi=mstart+i-1
         IF (nbas(i).ne.SIZE(H%eig(il-1,mi)%evals)) &
            call AbortWithError("Mismatch in block, eigval list sizes")
         evals1D(i,:nbas(i))=H%eig(il-1,mi)%evals(:)
      ENDDO

!     N-mode coupling and excitation data for submodes in each state
      call getsubmodedata(il,im,H,ML,nbas,nmode,nexci,nmsum,nesum)

!     Generate the guess states
      call sortDPeigvalsGen(nbloc,evalsND,qns,evals1D,nbas,&
                            nmode,nexci,nmsum,nesum,constraint,&
                            Etarget)
!      call sortDPeigvals(nbloc,evalsND,qns,evals1D,nbas)
!      call structuredDPeigvals(nbloc,evalsND,qns,evals1D,nbas)

      IF (nsubm.gt.1 .or. (nsubm.eq.1 .and. nbloc.lt.maxbas)) THEN

         IF (nsubm.gt.1) THEN
            write(*,'(/3X,A/)') 'Initial guess product functions:'
         ELSE
            select case (cpp%truncation)
            case(0)
               tag='energy'
            case(1)
               tag='n-mode coupling'
            case(2)
               tag='total excitation'
            case default
               tag='invalid truncation choice'
            end select
            write(*,'(/3X,2(A,I0),A,A/)') &
            'Truncating basis from ',maxbas,' to ',nbloc,&
            ' functions by ',trim(adjustl(tag))
         ENDIF

         write(frmt,'(2(A,I0),A)') '(',4*nsubm+26,'X,',nagn,'(I2,X))'
         write(*,frmt) (ML%resort(j),j=mst,mfi)
         write(frmt,'(2(A,I0),A)') &
         '(3X,',nsubm,'(I3,X),f19.12,X,A,X',nagn,'(I2,X))'
         DO i=1,nbloc
            call GetFullAssignment(il,im,H,ML,qns(i,:),qnfull)
            write(*,frmt) (qns(i,j)-1,j=1,nsubm),evalsND(i),&
            '->',(qnfull(j)-1,j=1,nagn)
            deallocate(qnfull)
         ENDDO
         write(*,*)
      ENDIF

!     Build the N-D separable eigenfunctions
      call BuildProdFunctions(Q,nbas,qns)

      DEALLOCATE(qns,nbas,evals1D,nmode,nexci)

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

      subroutine getsubmodedata(il,im,H,ML,nbas,nmode,nexci,nmsum,nesum)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (MLtree), INTENT(IN)      :: ML
      integer, intent(in)  :: il,im
      integer, intent(inout) :: nbas(:),nmode(:,:),nexci(:,:)
      integer, intent(out)   :: nmsum,nesum
      integer, allocatable :: subqns(:)
      integer :: i,j,k,nsubm,neibas

      nsubm=ML%modcomb(il,im)
      nmsum=0
      nesum=0
      DO j=1,nsubm
         neibas=0
         DO i=1,nbas(j)
!           Get assignment of the single-mode function from the tree
            call GetPartialAssignment(il,im,H,ML,j,i,subqns)

!           Calculate the nmode and the nexci values for each single mode fxn
            nmode(i,j)=0
            nexci(i,j)=0
            DO k=1,SIZE(subqns)
               IF (subqns(k).gt.1) nmode(i,j)=nmode(i,j)+1
               nexci(i,j)=nexci(i,j)+subqns(k)-1
            ENDDO

            if (i.eq.1) nmsum=nmsum+SIZE(subqns)
            if (nexci(i,j).gt.neibas) neibas=nexci(i,j)
            DEALLOCATE(subqns)
         ENDDO
         nesum=nesum+neibas
      ENDDO

      end subroutine getsubmodedata

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortDPeigvalsGen(nbloc,evalsND,qns,evals1D,nbas,&
                                  nmode,nexci,nmsum,nesum,constraint,&
                                  Etarget)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates all vecs in the direct-product representation (for testing)

      implicit none
      integer, intent(inout) :: qns(:,:)
      real*8, intent(inout)  :: evalsND(:)
      integer, intent(in)    :: nmode(:,:),nexci(:,:)
      integer, intent(in)    :: nbas(:)
      real*8, intent(in)     :: evals1D(:,:)
      integer, intent(in)    :: nbloc,nmsum,nesum,constraint
      real*8, intent(in)     :: Etarget
      integer, allocatable   :: qnstmp(:,:)
      real*8,  allocatable   :: evalsNDtmp(:),tabindx(:)
      integer, allocatable   :: subqns(:)
      integer :: nmtarget(2),netarget(2)
      integer :: j,k,ndof,nstate,mstate

      ndof=SIZE(nbas)

      select case(constraint)
      case(0) ! Sort only by energy
         nmtarget=(/0,nmsum/) ! Default: any number of modes coupled
         netarget=(/0,nesum/) ! Default: any excitation level
         call GetStatesinWindow(nbloc,evalsND,qns,evals1D,nbas,&
                                nmode,nexci,nmtarget,netarget,Etarget,&
                                nstate)

      case(1) ! Sort by increasing nmode, by energy for each nmode value
         netarget=(/0,nesum/)
         nstate=0
         allocate(evalsNDtmp(nbloc),qnstmp(nbloc,ndof))
         do j=0,nmsum
            nmtarget(1:2)=j
!           Get states having j modes coupled
            call GetStatesinWindow(nbloc,evalsNDtmp,qnstmp,evals1D,nbas,&
                                   nmode,nexci,nmtarget,netarget,Etarget,&
                                   mstate)

!           Add the states to the master list
            do k=1,mstate
               nstate=nstate+1
               qns(nstate,:)=qnstmp(k,:)
               evalsND(nstate)=evalsNDtmp(k)
               if (nstate.ge.nbloc) exit
            enddo
            if (nstate.ge.nbloc) exit
         enddo
         deallocate(evalsNDtmp,qnstmp)

      case(2) ! Sort by increasing nexci, by energy for each nexci value
         nmtarget=(/0,nmsum/)
         nstate=0
         allocate(evalsNDtmp(nbloc),qnstmp(nbloc,ndof))
         do j=0,nesum
            netarget(1:2)=j
!           Get states having j excitations
            call GetStatesinWindow(nbloc,evalsNDtmp,qnstmp,evals1D,nbas,&
                                   nmode,nexci,nmtarget,netarget,&
                                   Etarget,mstate)

!           Add the states to the master list
            do k=1,mstate
               nstate=nstate+1
               qns(nstate,:)=qnstmp(k,:)
               evalsND(nstate)=evalsNDtmp(k)
               if (nstate.ge.nbloc) exit
            enddo
            if (nstate.ge.nbloc) exit
         enddo
         deallocate(evalsNDtmp,qnstmp)

      case default
         call AbortWithError('sortDPeigvalsGen(): unknown constraint')
      end select

!     Make sure nbloc states were found
      if (nstate.lt.nbloc) then
         write(*,*) 'nstate (',nstate,' < ',nbloc,') nbloc'
         call AbortWithError('sortDPeigvalsGen(): found too few states')
      endif

!     Sort states by increasing energy
      allocate(qnstmp(nbloc,ndof),tabindx(nbloc))
      qnstmp(:,:)=qns(:,:)
      do k=1,nbloc
         tabindx(k)=k
      enddo
      call dsort(evalsND,tabindx,nbloc,2)
      do k=1,nbloc
         qns(k,:)=qnstmp(int(tabindx(k)),:)
      enddo

      deallocate(qnstmp,tabindx)

      end subroutine sortDPeigvalsGen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine structuredDPeigvals(nbloc,evalsND,qns,evals1D,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates all vecs in the direct-product representation (for testing)

      implicit none
      integer, intent(inout) :: qns(:,:)
      real*8, intent(inout)  :: evalsND(:)
      integer, intent(in)  :: nbas(:)
      real*8, intent(in)   :: evals1D(:,:)
      integer, allocatable :: indx(:)
      integer, intent(in)  :: nbloc
      integer :: i,j,ndof

      ndof=SIZE(nbas)
      ALLOCATE(indx(ndof))
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
      integer, intent(inout) :: qns(:,:)
      real*8, intent(inout)  :: evalsND(:)
      integer, allocatable, intent(in) :: nbas(:)
      real*8, allocatable, intent(in)  :: evals1D(:,:)
      integer, allocatable :: qnlist(:,:),tmpqn(:)
      real*8, allocatable  :: sums(:)
      integer, intent(in)  :: nbloc
      integer :: i,j,k,l,ndof,prodND,nguess,ii,jj
      integer :: minind(1),nallowed
      real*8  :: sumt,bignr
      logical :: allowed,add2end

      ndof=SIZE(nbas)
      bignr=1.0d99

      ALLOCATE(qnlist(ndof*nbloc,ndof),sums(ndof*nbloc),tmpqn(ndof))

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
