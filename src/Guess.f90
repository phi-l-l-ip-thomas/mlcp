!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
      USE INPUTCP
      USE TARGETEDSTATES

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_Guess_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_Guess_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_Guess_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_Guess_Module()
      call Get_MPI_Timings('Guess module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_Guess_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GuessPsi(im,evalsND,delta,bounds,Q,Ham,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(OUT) :: Q(:)
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CPpar), INTENT(IN)  :: cpp
      integer, intent(in)  :: im
      real(kind=8), allocatable, intent(out) :: evalsND(:),delta(:)
      real(kind=8), intent(inout) :: bounds(2)
      integer :: nsubm

      nsubm=Ham%nt(im)%nsubm()

      IF (nsubm.eq.0) THEN ! Bottom-layer nodes
         call GuessPsi_primitive(im,evalsND,delta,bounds,Q,Ham)
      ELSE                 ! All other nodes
         call GuessPsi_node(im,evalsND,delta,bounds,Q,Ham,ML,cpp)
      ENDIF

      end subroutine GuessPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GuessPsi_primitive(im,evals,delta,bounds,Q,Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(OUT) :: Q(:)
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in)  :: im
      integer :: nbas(1)
      real(kind=8), allocatable, intent(out) :: evals(:),delta(:)
      real(kind=8), intent(inout) :: bounds(2)
      integer, allocatable :: qns(:,:)
      integer :: i,nbloc,nHterm
      character*64 :: frmt

      nbloc=Ham%nt(im)%B
      nbas(1)=nbloc
      bounds(:)=0.d0

!     Use harmonic oscillator functions as guesses
      ALLOCATE(Q(nbloc),evals(nbloc),delta(nbloc),qns(nbloc,1))

      DO i=1,nbloc
         qns(i,1)=i
         evals(i)=(2*i-1)*Ham%nt(im)%Hfacs(1) ! extract omega from Ham
         delta(i)=0.d0
      ENDDO

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(3X,A/)') 'Initial guess harmonic functions:'
         write(frmt,'(A)') '(26X,I2,X)'
         write(*,frmt) Ham%nt(im)%dofs(1)
         write(frmt,'(A)') '(3X,f19.12,X,A,X,I2,X)'
         DO i=1,nbloc
            write(*,frmt) evals(i),'->',qns(i,1)-1
         ENDDO
         write(*,*)
      ENDIF

      call BuildProdFunctions(Q,nbas,qns)
      DEALLOCATE(qns)

      end subroutine GuessPsi_primitive

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GuessPsi_node(im,evalsND,delta,bounds,Q,Ham,ML,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the master routine for generating the initial guess
! wavefunction

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(OUT) :: Q(:)
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CPpar), INTENT(IN)  :: cpp
      integer, intent(in)  :: im
      integer, allocatable :: qns(:,:),nmode(:,:),nexci(:,:),nexmx(:,:)
      integer, allocatable :: qnfull(:),nbas(:)
      real(kind=8), allocatable  :: evals1D(:,:)
      real(kind=8), allocatable, intent(out) :: evalsND(:),delta(:)
      real(kind=8), intent(inout) :: bounds(2)
      integer :: i,j,nsubm,nbloc,obloc,nblocref,sm,prodND,prodNDref,mlil
      integer :: maxbas,ndof,nmsum,nesum,neibas,nemax,mlim
      real(kind=8) :: ti1,ti2
      character*64 :: frmt

      IF (.NOT. MODULE_SETUP) call Init_Guess_Module()

      call CPU_TIME(ti1)

!     Set parameters
      nsubm=Ham%nt(im)%nsubm()
      nbloc=Ham%nt(im)%B
      ndof=Ham%nt(im)%ndof()
      mlil=Ham%nt(im)%mlil
      mlim=Ham%nt(im)%mlim
      ALLOCATE(nbas(nsubm))
      DO j=1,nsubm
         sm=Ham%nt(im)%subm(j)
         nbas(j)=Ham%nt(sm)%nbas()
      ENDDO
      maxbas=MAXVAL(nbas)    ! max # b.f. of sub-nodes in node 'im'

!     Check block size against product basis and reference and
!     dynamically adjust if needed
      nblocref=ML%gdim(mlil,mlim)
      prodNDref=getboundingbaseprod(ML%gdim(mlil-1,&
                ML%modstart(mlil,mlim):&
                ML%modstart(mlil,mlim)+nsubm-1),nblocref)
      prodND=getboundingbaseprod(nbas,nbloc)

      IF (nbloc.gt.prodND) THEN
!        Constraints in the $node namelist for earlier nodes can 
!        reduce the product basis size so that it exceeds nbloc for this
!        node; if this happens then resize nbloc to fit the new size
         nbloc=prodND
         IF (mpirank.eq.mpi_prnt_rank) THEN
            IF (nsubm.eq.1) THEN
               write(*,'(X,A,X,I0,A,X,I0,X,A,X,I0,X,A/)') &
                '  * Reduced block size inherited from node',sm,':',&
               nbloc,'(overwrites',ML%gdim(mlil,mlim),'from layers.inp)'
            ELSE
               write(*,'(X,A,X,I0,X,A,X,I0,X,A/)') &
                '  * Block size reduced to fit truncated product basis:',&
               nbloc,'(overwrites',ML%gdim(mlil,mlim),'from layers.inp)'
            ENDIF
         ENDIF
         call Ham%nt(im)%setb(nbloc)

      ELSEIF (nbloc.lt.prodND .and. nblocref.eq.prodNDref) THEN
!        Contraints in the $node namelist for earlier nodes can
!        increase the product basis size so that the current nbloc no
!        longer equals than the product basis size, so resize nbloc
         prodND=1
         DO j=1,nsubm
            prodND=prodND*nbas(j)
         ENDDO
         nbloc=prodND
         IF (mpirank.eq.mpi_prnt_rank) &
         write(*,'(X,A,X,I0,X,A,X,I0,X,A/)') &
            '  * Block size increased to fit expanded product basis:',&
         nbloc,'(overwrites',ML%gdim(mlil,mlim),'from layers.inp)'
         call Ham%nt(im)%setb(nbloc)

      ENDIF

!     Allocate the arrays containing the mode information
      ALLOCATE(nmode(maxbas,nsubm),nexci(maxbas,nsubm),nexmx(maxbas,nsubm))

!     Copy the eigenvalues of the sub-modes to evals1D, and compute the
!     guess spectral range
      ALLOCATE(evals1D(nsubm,maxbas))
      bounds(:)=0.d0
      DO i=1,nsubm
         sm=Ham%nt(im)%subm(i)
         IF (nbas(i).ne.Ham%nt(sm)%nbas()) &
            call AbortWithError("Mismatch in block, eigval list sizes")
         evals1D(i,:nbas(i))=Ham%nt(sm)%eig(:)
         bounds(1)=bounds(1)+evals1D(i,1)
         bounds(2)=bounds(2)+evals1D(i,nbas(i))
      ENDDO

!     N-mode coupling and excitation data for submodes in each state
      call getsubmodedata(im,Ham,cpp,nbas,nmode,nexci,nexmx,&
                          nmsum,nesum,nemax)

!     Use pre-truncated node block size as nbloc
      obloc=Ham%nt(im)%b

      ALLOCATE(evalsND(obloc),delta(obloc),qns(obloc,nsubm))
      delta(:)=0.d0

!     Generate the guess states
      call sortDPeigvalsGen(obloc,evalsND,qns,evals1D,nbas,&
                            nmode,nexci,nexmx,nmsum,nesum,nemax,&
                            cpp%Etarget)

!     Pass number of states found in call above to tree node type
      nbloc=SIZE(evalsND)
      IF (mpirank.eq.mpi_prnt_rank .and. nbloc.ne.obloc) THEN
         write(*,'(X,A,X,I0)') '  * Block size initially set to    :',&
            obloc
         write(*,'(X,A,X,I0,X,A,X,I0,X,A/)') &
            '  * Block resized from constraints :',&
            nbloc,'(overwrites',ML%gdim(mlil,mlim),'from layers.inp)'
      ENDIF
      call Ham%nt(im)%setb(nbloc)

      IF (mpirank.eq.mpi_prnt_rank) THEN
         IF (nsubm.gt.1 .or. (nsubm.eq.1 .and. nbloc.lt.maxbas)) THEN

            IF (nsubm.gt.1) THEN
               write(*,'(3X,A/)') 'Initial guess product functions:'
            ELSE
               write(*,'(3X,2(A,I0),A/)') &
              'Truncating basis from ',maxbas,' to ',nbloc,&
               ' functions'
            ENDIF

            write(frmt,'(2(A,I0),A)') '(',4*nsubm+26,'X,',ndof,'(I2,X))'
            write(*,frmt) (Ham%nt(im)%dofs(j),j=1,ndof)
            write(frmt,'(2(A,I0),A)') &
            '(3X,',nsubm,'(I3,X),f19.12,X,A,X,',ndof,'(I2,X))'
            DO i=1,nbloc
               call ConstructAssignment(Ham%nt,im,qns(i,:),qnfull)
               write(*,frmt) (qns(i,j)-1,j=1,nsubm),evalsND(i),&
               '->',(qnfull(j),j=1,ndof)
               deallocate(qnfull)
            ENDDO
            write(*,*)
         ENDIF
      ENDIF

!     Build the N-D separable eigenfunctions
      ALLOCATE(Q(nbloc))
      call BuildProdFunctions(Q,nbas,qns)

      DEALLOCATE(qns,nbas,evals1D,nmode,nexci,nexmx)

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine GuessPsi_node

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GuessWeights(im,T,Ham) result(W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Guesses Boltzmann weights for separable H energies at temperature T,
! for weighting rank-reduction of CP matrices with multiple states. If 
! temperature is negative then identity matrix is guessed.

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CP) :: W
      real*8, intent(in)   :: T
      real*8, parameter    :: kb=0.69503476 ! cm^-1/K
      integer, intent(in)  :: im
      integer, allocatable :: nbas(:)
      logical, allocatable :: sym(:)
      integer :: i,j,sm,nsubm,msubm,mstart
      real*8 :: val
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_Guess_Module()

      call CPU_TIME(ti1)

!     Set parameters
      nsubm=Ham%nt(im)%nsubm()
      msubm=max(1,nsubm)

      ALLOCATE(nbas(msubm),sym(msubm))
      DO j=1,msubm
         sm=Ham%nt(im)%subm(j)
         if (nsubm.eq.0) then
            nbas(j)=Ham%nt(sm)%B
         else
            nbas(j)=Ham%nt(sm)%nbas()
         endif
         sym(j)=.FALSE.
      ENDDO

      W=IdentityCPMatrix(nbas,nbas,sym)

      IF (T.ge.0.d0 .or. nsubm.gt.0) THEN
!        Use the eigenvalues of the sub-modes to compute the weights
!        If T=0, force weights of 1 for g.s. and 0 for all other states
         DO i=1,nsubm
            sm=Ham%nt(im)%subm(i)
            DO j=1,nbas(i)
               IF (T.gt.0.d0) THEN
                  val=exp((Ham%nt(sm)%eig(1)-Ham%nt(sm)%eig(j))/(kb*T))
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

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end function GuessWeights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine getsubmodedata(im,Ham,cpp,nbas,nmode,nexci,nexmx,&
                                nmsum,nesum,nemax)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CPpar), INTENT(IN)  :: cpp
      integer, intent(in)  :: im
      integer, intent(inout) :: nbas(:),nmode(:,:),nexci(:,:),nexmx(:,:)
      integer, intent(out)   :: nmsum,nesum,nemax
      integer, allocatable :: subqns(:)
      integer :: i,j,k,sm,nsubm,neibas,mlil,mlim

      nsubm=Ham%nt(im)%nsubm()
      mlil=Ham%nt(im)%mlil
      mlim=Ham%nt(im)%mlim
      nmsum=0 ! Max number of coupled modes
      nesum=0 ! Max sum of quantum numbers
      nemax=0 ! Max value of individual quantum number

      DO j=1,nsubm
         sm=Ham%nt(im)%subm(j)
         neibas=0
         DO i=1,nbas(j)
!           Calculate the nmode and the nexci values for each single mode fxn
            nmode(i,j)=0
            nexci(i,j)=0
            nexmx(i,j)=MAXVAL(Ham%nt(sm)%assgn(i,:))
            nemax=MAX(nemax,MAXVAL(Ham%nt(sm)%assgn(i,:)))
            DO k=1,Ham%nt(sm)%ndof()
               IF (Ham%nt(sm)%assgn(i,k).gt.0) nmode(i,j)=nmode(i,j)+1
               nexci(i,j)=nexci(i,j)+Ham%nt(sm)%assgn(i,k)
            ENDDO
            if (i.eq.1) nmsum=nmsum+Ham%nt(sm)%ndof()
            if (nexci(i,j).gt.neibas) neibas=nexci(i,j)
         ENDDO
         nesum=nesum+neibas
      ENDDO

!     Replace node, sum, and qn limits with those from the input file
      if (cpp%max_nmode.ge.0) nmsum=min(nmsum,cpp%max_nmode)
      if (cpp%max_sum.ge.0)   nesum=min(nesum,cpp%max_sum)
      if (cpp%max_qn.ge.0)    nemax=min(nemax,cpp%max_qn)

      if (nesum.gt.nmsum*nemax) then
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(A,2(A,I0))') &
         '   * sum-max contraint exceeds (nmode-max * q.n.-max), ',&
         'so modifying sum-max: ',&
         nesum,' -> ',nmsum*nemax
         nesum=nmsum*nemax
      endif

      if (nemax.gt.nesum) then
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(A,2(A,I0))') &
         '   * q.n.-max constraint exceeds sum-max, ',&
         'so modifying q.n.-max: ',&
         nemax,' -> ',nesum
         nemax=nesum
      endif

      end subroutine getsubmodedata

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortDPeigvalsGen(nbloc,evalsND,qns,evals1D,nbas,&
                                  nmode,nexci,nexmx,nmsum,nesum,nemax,&
                                  Etarget)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates all vecs in the direct-product representation (for testing)

      implicit none
      integer, intent(in)  :: nbloc
      integer, allocatable, intent(inout) :: qns(:,:)
      real*8,  allocatable, intent(inout) :: evalsND(:)
      integer, intent(in)  :: nmode(:,:),nexci(:,:),nexmx(:,:)
      integer, intent(in)  :: nbas(:)
      real*8, intent(in)   :: evals1D(:,:)
      integer, intent(in)  :: nmsum,nesum,nemax
      real*8, intent(in)   :: Etarget
      integer, allocatable :: qnstmp(:,:)
      real*8,  allocatable :: evalsNDtmp(:),tabindx(:)
      integer, allocatable :: subqns(:)
      integer :: nmtarget(2),netarget(2),mxtarget(2)
      integer :: i,j,k,ndof,nstate,mstate

      ndof=SIZE(nbas)
      nmtarget=(/0,nmsum/)
      netarget=(/0,nesum/)
      mxtarget=(/0,nemax/)
!     Find all states allowed by limits imposed in layerfile
      call GetStatesinWindow(nbloc,evalsND,qns,evals1D,nbas,&
                             nmode,nexci,nexmx,nmtarget,netarget,&
                             mxtarget,Etarget,nstate)

!     Truncate block by nr. of states found and sort by energy
      allocate(evalsNDtmp(nstate),qnstmp(nstate,ndof),tabindx(nstate))
      evalsNDtmp(:)=evalsND(1:nstate)
      qnstmp(:,:)=qns(1:nstate,:)
      deallocate(evalsND,qns)
      do k=1,nstate
         tabindx(k)=k
      enddo
      call dsort(evalsNDtmp,tabindx,nstate,2)
      allocate(evalsND(nstate),qns(nstate,ndof))
      do k=1,nstate
         qns(k,:)=qnstmp(int(tabindx(k)),:)
      enddo
      evalsND(:)=evalsNDtmp(:)
      deallocate(evalsNDtmp,qnstmp,tabindx)

      end subroutine sortDPeigvalsGen

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

      integer function getboundingbaseprod(nbas,nbloc) result(prodND)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares product basis in 'nbas' with block size 'nbloc'
! Returns [-1,0,1] for PI_j nbas_j [.lt.,.eq.,.gt.] nbloc, respectively

      implicit none
      integer, intent(in) :: nbas(:)
      integer, intent(in) :: nbloc
      integer :: j,nsubm

      nsubm=SIZE(nbas)
      prodND=1
      DO j=1,nsubm
         prodND=prodND*nbas(j)
         IF (prodND.gt.nbloc) EXIT
      ENDDO

      end function getboundingbaseprod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE GUESS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
