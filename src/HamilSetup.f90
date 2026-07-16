!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE HAMILSETUP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs the multi-layer Hamiltonian and manages the operators

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MODECOMB
      USE NODETREE
      USE OPFUNCS
      USE LINALG
      USE CPCONFIG
      USE FFPES
      USE REDUCTION
      USE SEPDREPN

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE Hamiltonian
         TYPE (OperMat), ALLOCATABLE :: ops(:,:,:)
         TYPE (TN), ALLOCATABLE :: nt(:)
         LOGICAL, ALLOCATABLE   :: optable(:,:,:)
      END TYPE Hamiltonian

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_HamilSetup_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_HamilSetup_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_HamilSetup_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_HamilSetup_Module()
      call Get_MPI_Timings('HamilSetup module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_HamilSetup_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupHamiltonian(Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for building Hamiltonian matrix

      implicit none
      TYPE (MLtree)       :: ML
      TYPE (Hamiltonian)  :: Ham
      TYPE (Configs), ALLOCATABLE :: V(:),vtype(:)
      integer, allocatable :: opmap(:)
      real(kind=8), allocatable :: omega(:),alpha(:),beta(:)
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_HamilSetup_Module()

      call CPU_TIME(ti1)

      if (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A)') 'Hamiltonian setup...'

!     Compute PES or read from file
      call GetPotential(V,omega,alpha,beta,ML%system,ML%pes_path,&
                        ML%nmode(1),ML%dividefc,ML%pe_transform,&
                        ML%dpp%verbosity)

!     Generate the operator map and table
      call AllocHamilOp(Ham,V,ML%pe_transform,opmap,ML%dpp%verbosity)

!     PES coordinate transformation (if requested)
      call TransformPES(V,vtype,alpha,beta,ML%pe_transform,ML%pe_trans_fac,opmap,Ham%optable)

!     Print out PES info (for debugging)
      call ShowPESInfo(Ham,V,vtype,ML%pe_transform,opmap,ML%dpp%verbosity)

!     Construct node tree from ML tree
      call BuildNodeTree(Ham%nt,ML%modcomb,ML%modstart,ML%gdim,ML%nmode,ML%resort)

!     Sort Hamiltonian into multilayer format
      call FillHamilNodeTree(Ham%nt,V,vtype,omega,ML%dpp%verbosity)

!     Get the unique primitive operator matrices
      call GetPrimitiveOperators(Ham,ML,V,alpha,beta,opmap,ML%dpp%verbosity)

!     Construct bottom-layer mode operators from primitive operator
!     matrices, then solve and update primitive operators
!      call SolveandUpdateFirstLayer(Ham)

      DEALLOCATE(V,vtype,omega,alpha,beta,opmap)

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine SetupHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_Hamiltonian(Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes Hamiltonian module

      implicit none
      TYPE (Hamiltonian) :: Ham

      IF (ALLOCATED(Ham%ops)) DEALLOCATE(Ham%ops)
      IF (ALLOCATED(Ham%optable)) DEALLOCATE(Ham%optable)
      call FlushNodeTree(Ham%nt)

      end subroutine Flush_Hamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AllocHamilOp(Ham,V,trans,opmap,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates and sets arrays in Hamiltonian structure

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (Configs), INTENT(IN) :: V(:)
      character(len=*), intent(in) :: trans
      integer, allocatable, intent(out) :: opmap(:)
      integer, intent(in) :: verbosity
      integer :: i,ndof,ncoup,mxpow,ntyp

!     Set parameters
      ndof=V(1)%nbas(1)
      ncoup=SIZE(V)

!     Determine the max possible power, y, of q^y, from the PES
      mxpow=0
      DO i=1,ncoup
!        Skip if there are no terms for this i
         IF (SIZE(V(i)%coef).eq.1 .and. V(i)%coef(1).eq.0.d0) CYCLE
         mxpow=i
      ENDDO

!     Determine the opmap array from the transformation type
      IF (trans .seq. 'none') THEN
         ntyp=2
         allocate(opmap(ntyp))
         opmap=(/-1,0/)
      ELSEIF (trans .seq. 'poly-tanh') THEN
         ntyp=3
         allocate(opmap(ntyp))
         opmap=(/-1,0,1/)
      ELSEIF ((trans .seq. 'morse-tanh') .or. &
              (trans .seq. 'read-morse-tanh')) THEN
         ntyp=4
         allocate(opmap(ntyp))
         opmap=(/-1,0,1,2/)
      ELSE
         write(*,*) "Unrecognized coordinate transformation: '",trans,&
         "' ;allowed values are 'none', 'morse-tanh', and 'poly-tanh'"
         call AbortWithError("AllocHamilOp(): bad transformation type")
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.1) THEN
         write(*,'(X,A/)') 'Operator types present:'
         DO i=1,ntyp
            select case (opmap(i))
               case(-1)
                   write(*,'(X,A,I0,A)') '[',i,'] = (p_i)^n'
               case(0)
                   write(*,'(X,A,I0,A)') '[',i,'] = (q_i)^n'
               case(1)
                   write(*,'(X,A,I0,A)') &
                   '[',i,'] = (tanh(alpha_i*q_i))^n'
               case(2)
                   write(*,'(X,A,I0,A)') &
                   '[',i,'] = (1-exp(-beta_i*q_i))^n'
            end select
         ENDDO
         write(*,*)
      ENDIF

      ALLOCATE(Ham%ops(ndof,mxpow,ntyp),Ham%optable(ndof,mxpow,ntyp))
      Ham%optable(:,:,:)=.FALSE.

!     Set optable true for kinetic energy (p^2) terms
      Ham%optable(:,2,1)=.TRUE.

      end subroutine AllocHamilOp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowPESInfo(Ham,V,vtype,trans,opmap,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates and sets arrays in Hamiltonian structure

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs), INTENT(IN) :: V(:),vtype(:)
      character(len=64), intent(in) :: trans
      integer, intent(in) :: opmap(:)
      integer, intent(in) :: verbosity
      integer, allocatable :: modpowr(:,:)
      integer :: i,j,k,mjk,ndof,npow,nopt,ndf,ncoup
      integer, allocatable :: maptag(:,:)
      character(len=64) :: frmt
      character(len=64), allocatable, dimension(:) :: labl
      character(len=3), dimension(2) :: tag=(/' X ','   '/)

      IF (mpirank.eq.mpi_prnt_rank) THEN

         ndof=SIZE(Ham%optable,1)
         npow=SIZE(Ham%optable,2)
         nopt=SIZE(Ham%optable,3)
         ncoup=SIZE(V)

!        Print the operator table
         IF (verbosity.ge.1) THEN
            ALLOCATE(maptag(npow,nopt))

            write(*,'(/X,2A/)') 'Primitive operator table: ',&
                                '(X = operator is present)'
            write(frmt,'(2(A,I0),A)') &
                  '(X,A,',nopt,'(I2,',3*(npow-1)+1,'X,A))'
            write(*,frmt) 'Type|',(k,'|',k=1,nopt)
            write(frmt,'(2(A,I0),A)') &
                  '(A,',nopt,'(',npow,'(I2,X),A))'
            write(*,frmt) 'order|',((j,j=1,npow),'|',k=1,nopt)
            write(frmt,'(2(A,I0),A)') & 
                  '(X,A,',nopt,'(',npow,'A,A))'
            write(*,frmt) '-DOF|',(('---',j=1,npow),'|',k=1,nopt)

            DO i=1,ndof
               DO k=1,nopt
                  DO j=1,npow
                     mjk=(k-1)*npow+j
                     IF (Ham%optable(i,j,k)) THEN
                        maptag(j,k)=1
                     ELSE
                        maptag(j,k)=2
                     ENDIF
                  ENDDO
               ENDDO
               write(frmt,'(2(A,I0),A)') &
                    '(I4,X,A,',nopt,'(',npow,'A,A))'
               write(*,frmt) i,'|',((tag(maptag(j,k)),j=1,npow),'|',k=1,nopt)
            ENDDO
            write(*,*)
            DEALLOCATE(maptag)
         ENDIF

!        Print the transformed PES
         IF (verbosity.ge.2) THEN
!           Loop over quadratic, cubic, quartic, ... terms in the PES
            DO k=1,ncoup

!              If V(k) is a zero vector, skip
               IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

               write(*,'(X,2A,X,I0/)') 'Potential energy terms ',&
                                       '(post-processed), order:',k

               DO i=1,SIZE(V(k)%coef)
                  call DistribModePower(V(k)%qns(i,:),modpowr)
                  ndf=SIZE(modpowr,1)
                  allocate(labl(ndf))

                  do j=1,ndf
                     write(labl(j),'(3A,2(I0,A))') '(',&
                     TRIM(ADJUSTL(GetFunctionLabel(opmap(vtype(k)%qns(i,j))))),&
                     '^',modpowr(j,2),')_',modpowr(j,1)
                  enddo

                  write(frmt,'(A,I0,A)') '(ES18.10,',ndf,'(X,A,X,A))'
                  write(*,frmt) V(k)%coef(i),&
                  ('*',TRIM(ADJUSTL(labl(j))),j=1,ndf)

                  deallocate(labl,modpowr)
               ENDDO
               write(*,*)
            ENDDO

         ENDIF
      ENDIF

      end subroutine ShowPESInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FillHamilNodeTree(nt,V,vtype,omega,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills nt with terms from configuration array V

      implicit none
      TYPE (TN), INTENT(INOUT) :: nt(:)
      TYPE (Configs), INTENT(IN) :: V(:),vtype(:)
      real(kind=8), intent(in) :: omega(:)
      integer, intent(in)  :: verbosity
      integer, allocatable :: modpowr(:,:),commonnode(:,:)
      integer, allocatable :: dofresort(:),opcounts(:),opwidths(:)
      integer :: nnode,nlayr,thenode,nunassigned
      integer :: ipass,i,k,l,l2,ndof,ndf,ncoup

!     Set parameters
      ncoup=SIZE(V)
      ndof=V(1)%nbas(1)
      nnode=SIZE(nt)

!     Count layers by traversing the tree upwards from first node
      nlayr=1
      thenode=1
      do
         thenode=nt(thenode)%supern
         if (thenode.eq.0) exit
         nlayr=nlayr+1
      enddo

!     Get list of bottom-layer nodes for each dof
      allocate(dofresort(ndof))
      do i=1,ndof
         dofresort(nt(i)%dofs(1))=i
      enddo

      allocate(opcounts(nnode),opwidths(nnode))
      nunassigned=0

!     Loop over operators and assign to nodes
!     First pass: count number in each node and allocate arrays
!     Second pass: fill arrays
      DO ipass=1,2

         opcounts(:)=0 ! Number of operators in node
         opwidths(:)=0 ! Max number of product ops in node
         opcounts(1:ndof)=1 ! Number of operators in node
         opwidths(1:ndof)=1 ! Max number of product ops in node

!        Loop over potential terms from V
         DO k=1,ncoup
!           Skip if there are no terms for this k
            IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

!           Loop over terms with k multiplied operators
            DO i=1,SIZE(V(k)%coef)
               call DistribModePower(V(k)%qns(i,:),modpowr)
               ndf=SIZE(modpowr,1)
               allocate(commonnode(ndf,2))

!              Assign each operator to bottom-layer node
               do l=1,ndf
                  commonnode(l,1)=dofresort(modpowr(l,1))
                  commonnode(l,2)=1
               enddo

!              Migrate upwards through tree until reaching a node which
!              containing all operators in the term
               do l2=1,nlayr
                  if (ALL(commonnode(:,1).eq.commonnode(1,1))) exit
                  do l=1,ndf
                     thenode=commonnode(l,1)
                     commonnode(l,1)=nt(thenode)%supern
                     commonnode(l,2)=nt(thenode)%superi
                  enddo
               enddo

               thenode=commonnode(1,1)

               if (ipass.eq.1) then
                  if (thenode.gt.0) then
                     opcounts(thenode)=opcounts(thenode)+1
                     opwidths(thenode)=max(opwidths(thenode),ndf)
                  else ! Unassigned operator
                     nunassigned=nunassigned+1
                  endif
               else
                  if (thenode.gt.0) then ! fill nodetree
                     opcounts(thenode)=opcounts(thenode)+1
                     nt(thenode)%Hfacs(opcounts(thenode))=V(k)%coef(i)
                     nt(thenode)%Hnop(opcounts(thenode))=ndf
                     nt(thenode)%Hops(opcounts(thenode),1:ndf,1:2)=modpowr(1:ndf,1:2)
                     nt(thenode)%Hops(opcounts(thenode),1:ndf,3)=vtype(k)%qns(i,1:ndf)
                     nt(thenode)%Hsubm(opcounts(thenode),1:ndf)=commonnode(1:ndf,2)
                  endif
               endif

               deallocate(modpowr,commonnode)
            ENDDO
         ENDDO

         IF (ipass.eq.1) THEN
!           Allocate arrays for storing operators
            DO i=1,nnode
               IF (opcounts(i).gt.0) THEN
                  ALLOCATE(nt(i)%Hfacs(opcounts(i)))
                  ALLOCATE(nt(i)%Hnop(opcounts(i)))
                  ALLOCATE(nt(i)%Hops(opcounts(i),opwidths(i),3))
                  ALLOCATE(nt(i)%Hsubm(opcounts(i),opwidths(i)))
               ENDIF
            ENDDO
!           Include KEO terms in first pass
            DO i=1,ndof
               thenode=dofresort(i)
               nt(thenode)%Hfacs(1)=omega(i)
               nt(thenode)%Hnop(1)=1
               nt(thenode)%Hops(1,1,1)=i
               nt(thenode)%Hops(1,1,2)=2
               nt(thenode)%Hops(1,1,3)=1
               nt(thenode)%Hsubm(1,1)=1
            ENDDO
         ENDIF
      ENDDO

      IF (nunassigned.gt.0) THEN
         write(*,*) 'WARNING: ',nunassigned,&
         'operators not assigned to a node in the tree!'
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.2) THEN
         DO i=1,nnode
            call nt(i)%showstats()
            call nt(i)%showhamil()
         ENDDO
      ENDIF

      deallocate(dofresort,opcounts)

      end subroutine FillHamilNodeTree

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPrimitiveOperators(Ham,ML,V,alpha,beta,opmap,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines which operators are unique in the Hamiltonian and gets
! their operator matrices

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      TYPE (MLtree), INTENT(IN)  :: ML
      TYPE (Configs), INTENT(IN) :: V(:)
      real(kind=8), intent(in)   :: alpha(:),beta(:)
      integer, intent(in) :: opmap(:)
      integer, intent(in) :: verbosity
      integer :: i,j,k,l,m,gdim,ndof,ncoup,oppowmax,noptyp

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A)') "--> Generating primitive operators"

!     Set parameters
      ncoup=SIZE(V)
      ndof=SIZE(Ham%ops,1)
      oppowmax=SIZE(Ham%ops,2)
      noptyp=SIZE(Ham%ops,3)

!     Determine the max possible power, y, of q^y, from the PES
      oppowmax=0
      DO i=1,ncoup
!        Skip if there are no terms for this i
         IF (SIZE(V(i)%coef).eq.1 .and. V(i)%coef(1).eq.0.d0) CYCLE
         oppowmax=i
      ENDDO

!     Generate operator matrices for this mode
      DO i=1,ndof
         m=ML%resort(i)
         gdim=ML%gdim(1,i)
         DO k=1,noptyp
            l=opmap(k)
            DO j=1,oppowmax
               IF (l.lt.2) THEN
                  if (Ham%optable(m,j,k)) &
                  Ham%ops(m,j,k)=GetPrimitiveOperMat(m,gdim,j,l,alpha(m))
               ELSEIF (l.eq.2) THEN
                  if (Ham%optable(m,j,k)) &
                  Ham%ops(m,j,k)=GetPrimitiveOperMat(m,gdim,j,l,beta(m))
               ELSE
                  write(*,*) 'l has current max value of 2; l = ',l
                  call AbortWithError('GetPrimitiveOperators() bad l')
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.1) THEN
         write(*,'(/X,A/)') '--- Primitive H Operators ---'
         DO i=1,ndof
            m=ML%resort(i)
            DO k=1,noptyp
               l=opmap(k)
               DO j=1,oppowmax
                  IF (Ham%optable(m,j,k)) THEN
                     write(*,'(X,A)') Ham%ops(m,j,k)%label
                     if (verbosity.ge.2) &
                        call PrintVector(Ham%ops(m,j,k)%mat)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         if (verbosity.eq.1) write(*,*)
      ENDIF

      end subroutine GetPrimitiveOperators

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveandUpdateFirstLayer(Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles terms in the bottom-layer Hamiltonian by summing primitive
! operator matrices.

      implicit none
      TYPE (Hamiltonian) :: Ham
      TYPE (OperMat)     :: OM
      real(kind=8), allocatable :: Tmati(:,:),Tmatj(:,:),S(:,:)
      real(kind=8), allocatable :: eigv(:)
      integer, allocatable :: qns(:,:)
      integer :: i,j,k,n,dofi,nnodes,pow,maxpow,typ,maxtyp
      real(kind=8) :: fac

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A)') "--> Solving layer 1 Hamiltonian..."

      nnodes=SIZE(Ham%nt)
      maxpow=SIZE(Ham%ops,2)
      maxtyp=SIZE(Ham%ops,3)

!     Loop over bottom layer nodes
      DO i=1,nnodes

!        Exit upon encountering the first non-bottom-layer node
         IF (Ham%nt(i)%nsubm().gt.0) EXIT

         dofi=Ham%nt(i)%dofs(1)

!        Assemble the single-mode Hamiltonian
         DO j=1,Ham%nt(i)%nHterm()
            fac=Ham%nt(i)%Hfacs(j)
            pow=Ham%nt(i)%Hops(j,1,2)
            typ=Ham%nt(i)%Hops(j,1,3)
            call SumOperMats(OM,1.d0,Ham%ops(dofi,pow,typ),fac)
         ENDDO

!        Convert summed operator matrix to upper triangular form
         call Vec2SymPackMat(OM%mat,Tmati)

!        Solve generalzed eigenvalue problem
         n=SIZE(Tmati,1)
         ALLOCATE(eigv(n),S(n,n),qns(n,1))
         S=0.d0
         DO j=1,n
            S(j,j)=1.d0
            qns(j,1)=j
         ENDDO
         call SolveGenEigval(eigv,S,Tmati,'V')

!        Store eigenvalues/assignments in Hamiltonian type
         call SetEigenbasis(Ham%nt,i,qns,eigv)

!        Transform operators on this mode into eigen-basis
         DO k=1,maxtyp
            DO j=1,maxpow
               IF (Ham%optable(dofi,j,k)) THEN
                  call Vec2SymPackMat(Ham%ops(dofi,j,k)%mat,Tmatj)
                  call UnitaryTFM(Tmati,Tmatj)

!                 Replace operator matrix with transformed one
                  DEALLOCATE(Ham%ops(dofi,j,k)%mat)
                  call SymPackMat2Vec(Ham%ops(dofi,j,k)%mat,Tmatj)
                  DEALLOCATE(Tmatj)
              ENDIF
            ENDDO
          ENDDO

         DEALLOCATE(Tmati,OM%mat,eigv,S,qns)

         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(/,X,A,I0,X,A)') '--- NODE ',i,'---'
            write(*,'(/X,2A/)') 'Eigenvalues and assignments from ',&
                                'diagonalizing the 1-mode Hamiltonian:'

            call Ham%nt(i)%showeigen()
         ENDIF
      ENDDO

      end subroutine SolveandUpdateFirstLayer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine setnodeready(Ham,im)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets node readiness by making sure sub-nodes are ready

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in) :: im
      integer :: j,nsubm,sm
      logical :: isready

      nsubm=Ham%nt(im)%nsubm()

      isready=.TRUE.
      do j=1,nsubm
         sm=Ham%nt(im)%subm(j)
         isready=(isready.and.Ham%nt(sm)%done)
      enddo

      Ham%nt(im)%ready=isready

      end subroutine setnodeready

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE HAMILSETUP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
