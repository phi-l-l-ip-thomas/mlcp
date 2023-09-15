!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE HAMILSETUP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs the multi-layer Hamiltonian and manages the operators

      use ERRORTRAP
      use UTILS
      use MODECOMB
      use OPFUNCS
      use LINALG
      use CPCONFIG
      use FFPES
      use REDUCTION
      use HAMILOPT
      use SEPDREPN

      IMPLICIT NONE
      REAL*8, PRIVATE  :: Ham_time=0.d0

      TYPE EigList
         real*8, allocatable  :: evals(:)
         integer, allocatable :: assgn(:,:)
      END TYPE Eiglist

      TYPE Hamiltonian
         TYPE (OperMat), ALLOCATABLE :: pops(:),eops(:),cops(:),iops(:)
         TYPE (EigList), ALLOCATABLE :: eig(:,:)
         TYPE (IVEC), ALLOCATABLE :: mops(:,:)
         integer, allocatable :: nterms(:),nop(:,:),ops(:,:,:)
         integer, allocatable :: ndof(:,:),dofs(:,:,:),opindx(:,:)
         real*8, allocatable  :: facs(:,:),opcoef(:)
         integer :: redterms
      END TYPE Hamiltonian

      TYPE HopList
         integer, allocatable :: opid(:),mptr(:)
         real*8, allocatable  :: coef(:),mat(:)
      END TYPE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupHamiltonian(sys,opt,H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for building Hamiltonian matrix

      implicit none
      TYPE (MLtree)       :: ML
      TYPE (Hamiltonian)  :: H
      TYPE (Configs), ALLOCATABLE :: V(:)
      logical, intent(in) :: opt
      character(5), intent(in) :: sys
      real*8  :: t1,t2

      call CPU_TIME(t1)

      call GetPotential(V,sys,ML%nmode(1))

!     Rotate Hamiltonian to use "optimized" coordinates
      IF (opt) call OptimizePesdirections(V)

!     Store PES in type Hamiltonian, construct KEO
      call FillHamilType(V,H)

      DEALLOCATE(V)

!     Sort Hamiltonian into multilayer format
      call sortHamiltonian(H,ML)

!     Get the unique primitive operator matrices
      call FindUniqueOperators(H,ML)

!     Print structure of Hamiltonian
      call printHamiltonianInfo(H)

!     Construct bottom-layer mode operators from primitive operator
!     matrices, then solve and update primitive operators
      call SolveandUpdateFirstLayer(H)

      call CPU_TIME(t2)
      Ham_time=Ham_time+t2-t1

      end subroutine SetupHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushHamiltonian(H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes Hamiltonian module

      implicit none
      TYPE (Hamiltonian) :: H

      IF (ALLOCATED(H%pops)) DEALLOCATE(H%pops)
      IF (ALLOCATED(H%eops)) DEALLOCATE(H%eops)
      IF (ALLOCATED(H%cops)) DEALLOCATE(H%cops)
      IF (ALLOCATED(H%iops)) DEALLOCATE(H%iops)
      IF (ALLOCATED(H%eig)) DEALLOCATE(H%eig)
      IF (ALLOCATED(H%nop)) DEALLOCATE(H%nop)
      IF (ALLOCATED(H%ops)) DEALLOCATE(H%ops)
      IF (ALLOCATED(H%mops)) DEALLOCATE(H%mops)
      IF (ALLOCATED(H%nterms)) DEALLOCATE(H%nterms)
      IF (ALLOCATED(H%facs)) DEALLOCATE(H%facs)
      IF (ALLOCATED(H%ndof)) DEALLOCATE(H%ndof)
      IF (ALLOCATED(H%dofs)) DEALLOCATE(H%dofs)
      IF (ALLOCATED(H%opcoef)) DEALLOCATE(H%opcoef)
      IF (ALLOCATED(H%opindx)) DEALLOCATE(H%opindx)

      write(*,'(X,A,X,f20.3)') 'Total Hamiltonian generation time (s)',&
                             Ham_time

      end subroutine FlushHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushHopList(HL)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes Hamiltonian operator list type

      implicit none
      TYPE (HopList) :: HL

      IF (ALLOCATED(HL%opid)) DEALLOCATE(HL%opid)
      IF (ALLOCATED(HL%mptr)) DEALLOCATE(HL%mptr)
      IF (ALLOCATED(HL%coef)) DEALLOCATE(HL%coef)
      IF (ALLOCATED(HL%mat)) DEALLOCATE(HL%mat)

      end subroutine FlushHopList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine printHamiltonianInfo(H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes info about Hamiltonian to 'Hamiltonian.out'

      implicit none
      TYPE (Hamiltonian) :: H
      integer      :: u,il,nlayr,i,j,k,maxdof,poplen
      character*64 :: frmt

!      u = LookForFreeUnit()
!      open(u,status='unknown',file='Hamiltonian.out')

      nlayr=SIZE(H%nterms)
      write(*,'(/A/)') '********** Hamiltonian operator **********'

      poplen=SIZE(H%pops)
      write(*,*) '---- Primitive Operators ----'
      write(*,*) '  ID  Label'
      DO i=1,poplen
         write(*,'(I5,2X,A,I0)') i,&
         TRIM(H%pops(i)%label),H%pops(i)%dof
      ENDDO
      write(*,*)

      il=1
      maxdof=MAXVAL(H%ndof(:,il))
      write(*,'(A,I2,A/)') '---- Layer ',il,' ---'
      write(frmt,'(A,I2,A)') '(A,',3*(maxdof-1)+1,'X,A)'
      write(*,frmt) ' Term |DOF','|    Factor     Operators'
      DO i=1,H%nterms(il)
!        Generate the print format 
         write(frmt,'(A,I0,A,I0,A,I0,A)') &
           '(I5,X,A,',H%ndof(i,il),'(I2,X),',3*(maxdof-H%ndof(i,il))+1,&
           'X,A,f13.6,X,',H%ndof(i,il),'(I4,X))'
!        Print the line
         write(*,frmt) &
           i,'|',(H%dofs(i,j,il),j=1,H%ndof(i,il)),'|',&
           H%facs(i,1),(H%ops(i,j,1),j=1,H%ndof(i,il))
!        Loop over grouped terms (which operate on the same dofs)
         DO k=2,H%nop(i,il)
            write(frmt,'(A,I0,A,I0,A)') '(6X,A,',&
            3*maxdof+1,'X,A,f13.6,X,',H%ndof(i,il),'(I4,X))'
            write(*,frmt) '|','|',H%facs(i,k),&
            (H%ops(i,j,k),j=1,H%ndof(i,il))
         ENDDO
      ENDDO

      DO il=2,nlayr
         maxdof=MAXVAL(H%ndof(1:H%nterms(il),il))
         write(*,'(/A,I2,A)') '---- Layer ',il,' ---'
         write(frmt,'(A,I2,A)') '(A,',3*(maxdof-1)+1,'X,A)'
         write(*,frmt) ' Term |Mod','| Previous Layer Operators'
         DO i=1,H%nterms(il)
!           Generate the print format 
            write(frmt,'(A,I0,A,I0,A,I0,A,I0,A)') &
              '(I5,X,A,',H%ndof(i,il),'(I2,X),',&
              3*(maxdof-H%ndof(i,il))+1,'X,A,',&
              H%nop(i,il),'(I0,X))'
!           Print the line
            write(*,frmt) &
              i,'|',(H%dofs(i,j,il),j=1,H%ndof(i,il)),'|',&
              (H%mops(i,il)%v(k),k=1,H%nop(i,il))

         ENDDO
      ENDDO

!      write(*,'(/A/)') '--- Primitive H matrices ---'
!      DO i=1,poplen
!         write(*,'(X,2(A,I0))') 'ID # ',i,', dof = ',H%pops(i)%dof
!         call PrintVector(H%pops(i)%mat)
!      ENDDO
      write(*,'(/A/)') '******************************************'

      end subroutine printHamiltonianInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortHamiltonian(H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Organizes the Hamiltonian by sorting and combining terms for each dof

      implicit none
      TYPE (MLtree)        :: ML
      TYPE (Hamiltonian)   :: H
      real*8, allocatable  :: facs(:,:)
      integer, allocatable :: modecnt(:),ndof(:),nop(:),iop(:)
      integer, allocatable :: dofs(:,:),ops(:,:,:)
      integer :: i,j,k,l,pass,ndim,ntrm,inewtrm,nnewtrm,maxdof,maxop,ip
      integer :: il,mil,nm,nlayr,thismode
      logical :: unique

      write(*,'(/X,A)') "--> Sorting Hamiltonian terms into layers"

      nlayr=ML%nlayr
      maxdof=MAXVAL(H%ndof(:,1))  ! Max # coupled DOF in H

!     Rearrange the modes in the Hamiltonian before sorting
      call ResortHmodes(H,ML)

!     Loop over layers in multilayer wavefunction
      DO il=1,nlayr

         mil=max(il-1,1)
         ntrm=H%nterms(mil)             ! # terms before sorting
         maxdof=MAXVAL(H%ndof(1:ntrm,mil))  ! Max # coupled modes in H
         allocate(dofs(ntrm,maxdof))    ! list of coupled modes per term
         allocate(ndof(ntrm),nop(ntrm)) ! # modes, # operators per term
         nop=0

!        First pass: check terms in H to calc new array dimensions
!        Second pass: combine terms operating on a common set of modes

         DO pass=1,2

            inewtrm=0

!           Loop over terms in previous layer of Hamiltonian
            DO i=1,ntrm
               unique=.TRUE.
!              Determine the current-layer modes of the i-th term by
!              examining the previous-layer modes
               ALLOCATE(modecnt(H%ndof(i,mil)))
!              Bottom layer: no previous layer to examine, so just
!              build the mode-count list from the DOF-list
               IF (il.eq.1) THEN 
                  modecnt(:)=H%dofs(i,1:H%ndof(i,mil),il)
                  nm=H%ndof(i,mil)
!              Higher layers: find the current-layer mode containing
!              the previous-layer mode
               ELSE
                  modecnt(1)=ML%whichmod(mil,H%dofs(i,1,mil))
                  nm=1
                  DO k=2,H%ndof(i,mil)
                     thismode=ML%whichmod(mil,H%dofs(i,k,mil))
!                    If current-layer mode not duplicated, add to list
                     IF (.NOT.ANY(modecnt(1:nm).eq.thismode)) THEN
                        nm=nm+1
                        modecnt(nm)=thismode
                     ENDIF
                  ENDDO
               ENDIF

!              Loop over stored terms in current layer Hamiltonian to
!              determine if the i-th term in the previous layer 
!              Hamiltonian can be summed with one of the existing terms
!              in the current layer. Otherwise a new term is added to
!              the current layer

               DO j=1,inewtrm
                  IF (IntListsRSame(modecnt(1:nm),dofs(j,1:ndof(j)))) &
                     unique=.FALSE.

!                 For now, all multi-mode operators are unique
                  IF (ndof(j).gt.1) unique=.TRUE.

!                 Not unique: sum with a previous term in current layer
                  IF (.not.unique) THEN
                     IF (pass.eq.1) THEN  ! First pass
!                       Bottom layer: increment nop by # ops in term
                        IF (il.eq.1) THEN
                           nop(j)=nop(j)+H%nop(i,1)
!                       Higher layers: increment nop by 1 since only the
!                       pointer to the previous layer is needed
                        ELSE
                           nop(j)=nop(j)+1
                        ENDIF
                     ELSE  ! Second pass
                        IF (il.eq.1) THEN  ! add to facs and ops
                           DO k=1,H%nop(i,1)
                              facs(j,iop(j))=H%facs(i,k)
                              DO l=1,H%ndof(i,1)
                                 ops(j,l,iop(j))=H%ops(i,l,k)
                              ENDDO
                              iop(j)=iop(j)+1
                           ENDDO
                        ELSE  ! add the pointer to the previous layer
                           H%mops(j,il)%v(iop(j))=i
                           iop(j)=iop(j)+1
                        ENDIF
                     ENDIF
                     EXIT
                  ENDIF
               ENDDO  ! loop over j (prev. terms in current layer)

!              Unique combination of DOFs: new term on current layer
               IF (unique) THEN ! The term is counted in the dofs array
                  inewtrm=inewtrm+1
                  IF (pass.eq.1) THEN  ! First pass: fill dofs, ndof, nop
                     ndof(inewtrm)=nm
                     dofs(inewtrm,1:nm)=modecnt(1:nm)
!                    Bottom layer: nop holds number of operators in term
                     IF (il.eq.1) THEN
                        nop(inewtrm)=H%nop(i,1)
!                    Higher layers: nop holds number of pointers to 
!                    operators on previous layer
                     ELSE
                        nop(inewtrm)=1
                     ENDIF
                  ELSE  ! Second pass
                     IF (il.eq.1) THEN  ! Bottom layer: fill facs, ops
                        DO k=1,H%nop(i,1)
                           facs(inewtrm,iop(inewtrm))=H%facs(i,k)   
                           DO l=1,nm
                              ops(inewtrm,l,iop(inewtrm))=H%ops(i,l,k)
                           ENDDO
                           iop(inewtrm)=iop(inewtrm)+1
                        ENDDO
                     ELSE ! Higher layers: just add pointer to mops
                        H%mops(inewtrm,il)%v(iop(inewtrm))=i
                        iop(inewtrm)=iop(inewtrm)+1
                     ENDIF
                  ENDIF
               ENDIF
               DEALLOCATE(modecnt)
            ENDDO  ! loop over i (terms in H)

!           After the 1st pass, allocate arrays used in the 2nd pass
!           nnewtrm = number of terms in the 'condensed' H
!           maxop = number of summed terms for each term in H

            IF (pass.eq.1) THEN
               maxop=MAXVAL(nop)
               nnewtrm=inewtrm
               allocate(iop(nnewtrm))
               IF (il.eq.1) THEN
                  allocate(facs(nnewtrm,maxop))
                  allocate(ops(nnewtrm,maxdof,maxop))
               ELSE
                  IF (il.eq.2) THEN
                     allocate(H%mops(nnewtrm,2:nlayr))
                  ENDIF
                  DO k=1,nnewtrm
                     call H%mops(k,il)%new(nop(k))
                  ENDDO
               ENDIF
               iop=1
            ENDIF

         ENDDO  ! loop over pass

!        Allocate arrays in H and fill with contents of temp arrays
         IF (il.eq.1) THEN
            deallocate(H%ndof,H%nop,H%dofs,H%facs,H%ops,H%nterms)
            allocate(H%ndof(nnewtrm,nlayr),H%nop(nnewtrm,nlayr))
            allocate(H%dofs(nnewtrm,maxdof,nlayr),H%facs(nnewtrm,maxop))
            allocate(H%ops(nnewtrm,maxdof,maxop),H%nterms(nlayr))

!           facs and ops apply only to the bottom layer, so fill here
            DO i=1,nnewtrm
               DO k=1,nop(i)
                  H%facs(i,k)=facs(i,k)
                  DO j=1,ndof(i)
                     H%ops(i,j,k)=ops(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
            deallocate(facs,ops)
         ENDIF

!        Copy temp arrays to Hamiltonian
         DO i=1,nnewtrm
            H%ndof(i,il)=ndof(i)
            H%nop(i,il)=nop(i)
            DO j=1,ndof(i)
               H%dofs(i,j,il)=dofs(i,j)
            ENDDO
         ENDDO
         deallocate(ndof,nop,dofs,iop)
         H%nterms(il)=nnewtrm

      ENDDO  ! loop over il

!     Allocate the array for storing the eigenvalues
      ALLOCATE(H%eig(nlayr,ML%nmode(1)))

      end subroutine sortHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ResortHmodes(H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rearranges the modes on the bottom layer of the Hamiltonian
! This is done to allow arbitrary mode combination schemes

      implicit none
      TYPE (MLtree)        :: ML
      TYPE (Hamiltonian), INTENT(INOUT)   :: H
      integer, allocatable :: dummy(:,:)
      integer :: i,j,k,l,ndum
      logical :: usedummy,resortfound

      ALLOCATE(dummy(ML%nmode(1),2))
      ndum=0

      DO k=1,ML%nmode(1)
         IF (k.ne.ML%resort(k)) THEN
            resortfound=.FALSE.
            DO l=1,ndum
!              resorted index matches the dummy, so replace dummy
               IF (dummy(l,1).eq.ML%resort(k)) THEN
                  resortfound=.TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT.resortfound) THEN
!              If the k-th DOF is found, replace with a dummy value
               usedummy=.FALSE.
               DO i=1,H%nterms(1)
                  DO j=1,H%ndof(i,1)
                     IF (H%dofs(i,j,1).eq.k) THEN
                        H%dofs(i,j,1)=-k
                        IF (.NOT. usedummy) THEN
                           usedummy=.TRUE.
                           ndum=ndum+1
                           dummy(ndum,1)=k
                           dummy(ndum,2)=-k
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
!              Now replace the resorted index with the k-th one
               DO i=1,H%nterms(1)
                  DO j=1,H%ndof(i,1)
                     IF (H%dofs(i,j,1).eq.ML%resort(k)) THEN
                        H%dofs(i,j,1)=k
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO

!     Resolve indices mapped to dummies
      DO k=1,ML%nmode(1)
         DO l=1,ndum
!           resorted index matches the dummy, so replace dummy
            IF (dummy(l,1).eq.ML%resort(k)) THEN
               DO i=1,H%nterms(1)
                  DO j=1,H%ndof(i,1)
                     IF (H%dofs(i,j,1).eq.dummy(l,2)) THEN
                        H%dofs(i,j,1)=k
                     ENDIF
                  ENDDO
               ENDDO
               EXIT
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(dummy)

      end subroutine ResortHmodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FindUniqueOperators(H,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines which operators are unique in the Hamiltonian and gets
! their operator matrices

      implicit none
      TYPE (Hamiltonian)   :: H
      TYPE (MLtree)        :: ML
      integer, allocatable :: opid(:,:,:),udof(:),uop(:)
      integer :: il,pass,i,j,k,l,opct,maxops
      logical :: unique

      write(*,'(X,A)') "--> Determining unique primitive operators"

!     Find terms for bottom layer only
      il=1

!     Array for operator ID numbers
      ALLOCATE(opid(SIZE(H%ops,1),SIZE(H%ops,2),SIZE(H%ops,3)))

      maxops=0
      DO pass=1,2
         opct=0
         DO i=1,H%nterms(il)
            DO j=1,H%ndof(i,il)
               DO k=1,H%nop(i,il)
!                 First pass: count max possible # operators
                  IF (pass.eq.1) THEN
                     maxops=maxops+1
!                 Second pass: make list of unique operators
                  ELSEIF (pass.eq.2) THEN
                     unique=.TRUE.
                     DO l=1,opct
                        IF (H%dofs(i,j,il).eq.udof(l) .and. &
                            H%ops(i,j,k).eq.uop(l)) THEN
                           unique=.FALSE.
                           opid(i,j,k)=l
                        ENDIF
                     ENDDO
!                    Add the unique operator to the list
                     IF (unique) THEN
                        opct=opct+1
                        udof(opct)=H%dofs(i,j,il)
                        uop(opct)=H%ops(i,j,k)
                        opid(i,j,k)=opct
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         IF (pass.eq.1) THEN
            ALLOCATE(udof(maxops),uop(maxops))
         ELSEIF (pass.eq.2) THEN
            ALLOCATE(H%pops(opct))
         ENDIF

      ENDDO

!     Now fill the unique primitive operators array
      DO l=1,opct
         H%pops(l)=&
         GetPrimitiveOperMat(udof(l),ML%gdim(il,udof(l)),uop(l))
      ENDDO

!     Copy opid array to H
      H%ops=opid
      DEALLOCATE(udof,uop,opid)

      end subroutine FindUniqueOperators

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveandUpdateFirstLayer(H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles terms in the bottom-layer Hamiltonian by summing primitive
! operator matrices.

      implicit none
      TYPE (Hamiltonian)  :: H
      TYPE (OperMat)      :: OM
      real*8, allocatable :: Tmati(:,:),Tmatj(:,:),S(:,:)
      real*8, allocatable :: eigvals(:)
      integer :: i,j,k,colsi,symi,dofi,symj,dofj
      real*8  :: fac

      write(*,'(X,A)') "--> Solving layer 1 Hamiltonian..."

!     Loop over terms in Hamiltonian
      DO i=1,H%nterms(1)

!        If term is a single-dof operator, assemble and diagonalize
         IF (H%ndof(i,1).eq.1) THEN 

!           Sum the operators which act only on this DOF
            fac=H%facs(i,1)
            call SumOperMats(OM,H%pops(H%ops(i,1,1)),fac)
            DO k=2,H%nop(i,1)
               fac=H%facs(i,k)
               call SumOperMats(OM,1.d0,H%pops(H%ops(i,1,k)),fac)
            ENDDO

            dofi=OM%dof

!           Convert summed operator matrix to upper triangular form
            call Vec2SymPackMat(OM%mat,Tmati)

!           Solve generalzed eigenvalue problem
            ALLOCATE(eigvals(SIZE(Tmati,1)))
            ALLOCATE(S(SIZE(Tmati,1),SIZE(Tmati,1)))
            S=0.d0
            DO j=1,SIZE(Tmati,1)
               S(j,j)=1.d0
            ENDDO
            call SolveGenEigval(eigvals,S,Tmati,'V')

!           Store eigenvalues/assignments in Hamiltonian type
            ALLOCATE(H%eig(1,dofi)%evals(SIZE(Tmati,1)))
            ALLOCATE(H%eig(1,dofi)%assgn(SIZE(Tmati,1),1))
            DO j=1,SIZE(Tmati,1)
               H%eig(1,dofi)%assgn(j,1)=j
            ENDDO
            H%eig(1,dofi)%evals(:)=eigvals(:)
            DEALLOCATE(eigvals,S)

!           Use the eigenvector matrix to transform the primitive
!           operators for the same DOF/mode
            DO j=1,SIZE(H%pops)
               dofj=H%pops(j)%dof
               IF (dofj.eq.dofi) THEN  ! Transform the matrix
!                 Convert operator matrix to regular form and
!                 transform using the eigenbasis
                  call Vec2SymPackMat(H%pops(j)%mat,Tmatj)
                  call UnitaryTFM(Tmati,Tmatj)

!                 Replace operator matrix with transformed one
                  DEALLOCATE(H%pops(j)%mat)
                  call SymPackMat2Vec(H%pops(j)%mat,Tmatj)
                  DEALLOCATE(Tmatj)
               ENDIF
            ENDDO
            DEALLOCATE(Tmati)
            DEALLOCATE(OM%mat)
         ENDIF
      ENDDO

      end subroutine SolveandUpdateFirstLayer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetModeHNr(il,im,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Find the single-mode term in the Hamiltonian corresponding to mode im
! in layer il

      implicit none
      TYPE (Hamiltonian)   :: H
      integer, intent(in)  :: il,im
      integer :: i,GetModeHNr

      GetModeHNr=-1

      DO i=1,H%nterms(il)
         IF (H%ndof(i,il).eq.1 .and. H%dofs(i,1,il).eq.im) GetModeHNr=i
      ENDDO

      end function GetModeHNr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FillHamilType(V,H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills TYPE (Hamiltonian) with terms from configuration array V

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: H
      TYPE (Configs), INTENT(IN) :: V(:)
      integer, allocatable :: modpowr(:,:)
      integer :: i,k,l,ndof,ndf,ncoup,hterms

!     Set parameters
      ncoup=SIZE(V)
      ndof=V(1)%nbas(1)

!     Start with ndof terms for the KEO
      hterms=ndof

!     Count the potential terms from V
      DO k=1,ncoup
!        Skip if there are no terms for this k
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE
         hterms=hterms+SIZE(V(k)%coef)
      ENDDO

      allocate(H%nterms(1),H%facs(hterms,1))
      allocate(H%dofs(hterms,ndof,1),H%ops(hterms,ndof,1))
      allocate(H%ndof(hterms,1),H%nop(hterms,1))
      H%nterms(1)=hterms

      l=1
      DO k=1,ncoup

!        If V(k) is a zero vector, skip
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

         DO i=1,SIZE(V(k)%coef)
            call DistribModePower(V(k)%qns(i,:),modpowr)
            ndf=SIZE(modpowr,1)

!           KEO: duplicate quadratic terms, but with KEO operator flag
            IF (k.eq.2 .and. ndf.eq.1) THEN
               H%facs(l,1)=V(k)%coef(i)
               H%ndof(l,1)=1
               H%nop(l,1)=1
               H%dofs(l,1,1)=modpowr(1,1)
               H%ops(l,1,1)=-2 ! Flag for KEO
               l=l+1
            ENDIF

            H%facs(l,1)=V(k)%coef(i)
            H%ndof(l,1)=ndf
            H%nop(l,1)=1
            H%dofs(l,1:ndf,1)=modpowr(1:ndf,1)
            H%ops(l,1:ndf,1)=modpowr(1:ndf,2)
            l=l+1
            deallocate(modpowr)
         ENDDO
      ENDDO

      end subroutine FillHamilType

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetFullAssignment(il,im,Ham,ML,qns,qnfull)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets configuration of primitive basis function underlying config in
! qns by tracing assignments stored in Ham to the bottom layer.

      implicit none
      TYPE (MLtree), INTENT(IN)      :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in)    :: qns(:)
      integer, allocatable, intent(out) :: qnfull(:)
      integer, intent(in)    :: il,im
      integer, allocatable   :: modind(:),qntmp(:)
      integer :: i,j,eigind,imn,nsubm,mst,nagn

!     Get the full assignment in terms of primitive DOFs
      ALLOCATE(modind(il),qntmp(ML%nmode(1)))

      qntmp=0
      nagn=0
      modind=1
      DO
         imn=im

         IF (modind(1).gt.1) EXIT

!        Trace each assignment to the bottom layer
         DO j=il,2,-1
            mst=ML%modstart(j,imn)

!           Take number of sub-modes from input array on first pass,
!           then extract the number from stored assignment array
            IF (j.eq.il) THEN
               nsubm=SIZE(qns)
            ELSE
               nsubm=SIZE(Ham%eig(j,imn)%assgn,2)
            ENDIF

            IF (modind(il-j+2).gt.nsubm) THEN
               modind(il-j+1)=modind(il-j+1)+1
               modind(il-j+2:)=1
               EXIT
            ENDIF

!           Use input array to get 'eigind' on first pass here, too
            IF (j.eq.il) THEN
               eigind=qns(modind(2))
            ELSE
               eigind=Ham%eig(j,imn)%assgn(eigind,modind(il-j+2))
            ENDIF

            imn=mst+modind(il-j+2)-1

            IF (j.eq.2) THEN
               nagn=nagn+1
               qntmp(nagn)=eigind
               modind(il-j+2)=modind(il-j+2)+1
            ENDIF
         ENDDO
      ENDDO

!     Copy qntmp to qnfull
      ALLOCATE(qnfull(nagn))
      qnfull(:)=qntmp(:nagn)
      DEALLOCATE(modind,qntmp)

      end subroutine GetFullAssignment

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE HAMILSETUP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
