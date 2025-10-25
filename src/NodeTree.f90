!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE NODETREE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains node class and functions

      USE ERRORTRAP
      USE UTILS
      USE MYMPI

      implicit none

      TYPE TN
         integer :: nid,supern,superi,mlil,mlim,B
         integer, allocatable :: subs(:),dofs(:)
         integer, allocatable :: Hnop(:)   ! Nr. operators in each H term
         integer, allocatable :: Hops(:,:,:) ! Prim. Op IDs in each H term
         integer, allocatable :: Hsubm(:,:)! Sub-mode in each H term
         real(kind=8), allocatable :: Hfacs(:) ! Factors in each H term
         real(kind=8), allocatable :: eig(:), delta(:) ! Eigenvalue list
         integer, allocatable :: assgn(:,:)  ! Assignment in submode basis
         CONTAINS
            PROCEDURE :: new => NewTreeNode
            PROCEDURE :: flush => FlushTreeNode
            PROCEDURE :: showstats => ShowTNStats
            PROCEDURE :: showhamil => ShowTNHamil
            PROCEDURE :: showeigen => ShowTNEigen
            PROCEDURE :: subm => GetTNsubm
            PROCEDURE :: nsubm => GetTNnsubm
            PROCEDURE :: ndof => GetTNndof
            PROCEDURE :: nHterm => GetTNnHterm
            PROCEDURE :: nbas => GetTNnbas
            PROCEDURE :: setb => SetBlockSize
      END TYPE TN

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewTreeNode(n,nid,supern,superi,subs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs node with indices indicating its placement in tree

      implicit none
      CLASS (TN) :: n
      integer, intent(in) :: nid,supern,superi
      integer, allocatable, intent(in) :: subs(:)

!     List of sub-mode indices
      if (ALLOCATED(subs)) then
         ALLOCATE(n%subs(SIZE(subs)))
         n%subs(:)=subs(:)
      endif

!     Super-node and index of this node in super-node
      n%nid=nid
      n%supern=supern
      n%superi=superi
      n%B=0

      end subroutine NewTreeNode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushTreeNode(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates node type

      implicit none
      CLASS (TN) :: n

      if (ALLOCATED(n%subs)) DEALLOCATE(n%subs)
      if (ALLOCATED(n%dofs)) DEALLOCATE(n%dofs)
      if (ALLOCATED(n%Hnop)) DEALLOCATE(n%Hnop)
      if (ALLOCATED(n%Hops)) DEALLOCATE(n%Hops)
      if (ALLOCATED(n%Hsubm)) DEALLOCATE(n%Hsubm)
      if (ALLOCATED(n%Hfacs)) DEALLOCATE(n%Hfacs)
      if (ALLOCATED(n%eig)) DEALLOCATE(n%eig)
      if (ALLOCATED(n%assgn)) DEALLOCATE(n%assgn)
      n%nid=0
      n%supern=-1
      n%superi=-1
      n%B=0

      end subroutine FlushTreeNode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowTNStats(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates node type

      implicit none
      CLASS (TN) :: n
      integer :: i,nsubm,nHterm,nbas,ndof
      character*64 :: frmt


      write(*,'(X,A,I0)') 'Stats for node :  ',n%nid
      if (n%supern.ne.0) then
         write(*,'(X,2(A,I0),A)') 'Node in super  :  ',n%supern,&
                                  ' (',n%superi,')'
      else
         write(*,*) 'Node in super  :  N/A (top layer)'
      endif

      nsubm=n%nsubm()
      ndof=n%ndof()
      nHterm=n%nHterm()
      nbas=n%nbas()
      if (nsubm.gt.0) then
         write(frmt,*) '(X,A,',nsubm,'(X,I0))'
         write(*,frmt) 'Node sub-modes : ',(n%subs(i),i=1,nsubm)
      else
         write(*,*) 'Node sub-modes :  N/A (bottom layer)'
      endif
      write(frmt,*) '(X,A,',n%ndof(),'(X,I0))'
      write(*,frmt) 'Node DOFs      : ',(n%dofs(i),i=1,ndof)
      write(*,'(X,A,I0)') 'Node H terms   : ',nHterm
      write(*,'(X,A,I0,A,L)') 'Node Block size: ',n%B,' written: ',&
                              (nbas.gt.0)
      write(*,*)

      end subroutine ShowTNStats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowTNHamil(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates node type

      implicit none
      CLASS (TN) :: n
      integer :: i,j,nsubm,nHterm,nbas,ndof
      character*64 :: frmt


      nsubm=n%nsubm()
      nHterm=n%nHterm()

      write(*,'(X,A,I0,A,I0/)') 'List of H terms: (',nHterm+nsubm,&
                               ' terms) for node ',n%nid
      
      do i=1,nsubm
         write(*,'(X,A,I0,A)') '[Node ',n%subs(i),' eigen-operator]'
      enddo
      do i=1,nHterm
         write(frmt,*) '(X,ES18.10,X,',n%Hnop(i),'(3(A,I0)))'
         write(*,frmt) n%Hfacs(i),('* ([',n%Hops(i,j,3),']_',&
                       n%Hops(i,j,1),')^',n%Hops(i,j,2),j=1,n%Hnop(i))
      enddo
      write(*,*)

      end subroutine ShowTNHamil

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowTNEigen(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates node type

      implicit none
      CLASS (TN) :: n
      integer :: i,j,nsubm,nHterm,nbas,ndof
      integer, allocatable :: exc(:,:)
      integer, parameter :: ndofmax=23
      integer :: nexc,iexc,wmod
      character*72 :: frmt
      character(len=4) :: labl

      nsubm=n%nsubm()
      ndof=n%ndof()
      nbas=n%nbas()
      wmod=int(log10(REAL(MAXVAL(n%dofs))))+1

      IF (nbas.eq.0) RETURN

!     Construct header
      IF (ndof.le.ndofmax) THEN
         write(frmt,*) '(A,X,',ndof,'(I2,X),5X,A,14X,A,13X,A,5X,A)'
         write(*,frmt) 'Mode:',(n%dofs(j),j=1,ndof),'Energy',&
                       'E-E0','delta','Assignment'
      ELSE
         write(frmt,*) '(A,6X,A,14X,A,13X,A,5X,A)'
         write(*,frmt) 'Mode:','Energy','E-E0','delta','Assignment'
      ENDIF

      ALLOCATE(exc(ndof,2))

      DO i=1,nbas
!        Get list of excited modes
         nexc=0
         iexc=0
         DO j=1,ndof
            IF (n%assgn(i,j).gt.0) THEN
               nexc=nexc+1
               iexc=n%assgn(i,j)
               exc(nexc,1)=n%dofs(j)
               exc(nexc,2)=n%assgn(i,j)
            ENDIF
         ENDDO
!        Label vibration if only 1 DOF is excited
         IF (nexc.eq.0) THEN
            labl='ZPVE'
         ELSEIF (nexc.eq.1 .and. iexc.eq.1) THEN
            labl='FUND'
         ELSEIF (nexc.eq.1 .and. iexc.gt.1) THEN
            labl=' OT '
         ELSE
            labl='    '
         ENDIF

         IF (ndof.le.ndofmax) THEN ! Few modes: use tabular format
            IF (nexc.ne.1) THEN ! ZPVE,combinations
               write(frmt,*) '(I4,A,X,',ndof,&
                          '(I2,X),2(f19.12,X),ES11.3,X,A)'
               write(*,frmt) i,')',(n%assgn(i,j),j=1,ndof),n%eig(i),&
                           n%eig(i)-n%eig(1),n%delta(i),labl
            ELSEIF (iexc.eq.1) THEN ! Fundamentals
               write(frmt,*) '(I4,A,X,',ndof,&
                  '(I2,X),2(f19.12,X),ES11.3,X,A,X,A,I0)'
               write(*,frmt) i,')',(n%assgn(i,j),j=1,ndof),n%eig(i),&
                    n%eig(i)-n%eig(1),n%delta(i),labl,'v_',exc(1,1)
            ELSE ! Overtones
               write(frmt,*) '(I4,A,X,',ndof,&
                   '(I2,X),2(f19.12,X),ES11.3,X,A,X,I0,A,I0)'
               write(*,frmt) i,')',(n%assgn(i,j),j=1,ndof),n%eig(i),&
                    n%eig(i)-n%eig(1),n%delta(i),labl,&
                    exc(1,2),'v_',exc(1,1)
            ENDIF

         ELSE ! Many modes: print only excited modes
            IF (nexc.eq.0) THEN
               write(frmt,*) '(I4,A,X,2(f19.12,X),ES11.3,X,A,X)'
               write(*,frmt) i,')',n%eig(i),n%eig(i)-n%eig(1),&
                             n%delta(i),labl
            ELSE
               write(frmt,'(3(A,I0),A)') &
                             '(I4,A,X,2(f19.12,X),ES11.3,X,A,X,',&
                             nexc,'(A,I',wmod,'.',wmod,',A,I0,A,X))'
               write(*,frmt) i,')',n%eig(i),n%eig(i)-n%eig(1),n%delta(i),&
                            labl,('v_',exc(j,1),'(',exc(j,2),')',j=1,nexc)
            ENDIF
         ENDIF
      ENDDO

      DEALLOCATE(exc)

      end subroutine ShowTNEigen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetTNsubm(n,i) result (subm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets number of sub-nodes owned by node

      implicit none
      CLASS (TN) :: n
      integer, intent(in) :: i
      integer :: nsubm,subm

      nsubm=n%nsubm()

      if (i.eq.1 .and. nsubm.eq.0) then
         subm=n%nid
      elseif (i.ge.1 .and. i.le.nsubm) then
         subm=n%subs(i)
      else
         write(*,'(2(A,I0),A)') &
         'Sub-mode ',i,' must be in range [1,',nsubm,']'
         call AbortWithError('GetTNsubm(): i out of range')
      endif

      end function GetTNsubm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetTNnsubm(n) result (nsubm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets number of sub-nodes owned by node

      implicit none
      CLASS (TN) :: n
      integer :: nsubm

      if (ALLOCATED(n%subs)) then
         nsubm=SIZE(n%subs)
      else
         nsubm=0
      endif

      end function GetTNnsubm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetTNndof(n) result (ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets number of DOFs contained in mode

      implicit none
      CLASS (TN) :: n
      integer :: ndof

      if (ALLOCATED(n%dofs)) then
         ndof=SIZE(n%dofs)
      else
         ndof=0
      endif

      end function GetTNndof

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetTNnHterm(n) result (nHterm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets number of H terms assigned to node

      implicit none
      CLASS (TN) :: n
      integer :: nHterm

      if (ALLOCATED(n%Hfacs)) then
         nHterm=SIZE(n%Hfacs)
      else
         nHterm=0
      endif

      end function GetTNnHterm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetTNnbas(n) result (nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets number of basis functions assigned to node

      implicit none
      CLASS (TN) :: n
      integer :: nbas

      if (ALLOCATED(n%eig)) then
         nbas=SIZE(n%eig)
      else
         nbas=0
      endif

      end function GetTNnbas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetBlockSize(n,b)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets block size for eigenvalue, assignment arrays

      implicit none
      CLASS (TN) :: n
      integer, intent(in) :: b

      if (allocated(n%eig).or.allocated(n%assgn)) then
         call AbortWithError("SetBlockSize(): basis already set")
      endif

      n%B=b

      end subroutine SetBlockSize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetEigenbasis(nt,inode,qns,eig,delta)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocates arrays for eigenvalues and assignments, Stores eigenvalues

      implicit none
      TYPE (TN), intent(inout) :: nt(:)
      integer, intent(in) :: inode
      integer, intent(in) :: qns(:,:)
      real(kind=8), intent(in) :: eig(:)
      real(kind=8), intent(in), optional :: delta(:)
      integer :: i,j,k,sm,neig,nsubm,ndof,nsubdof

      neig=nt(inode)%B
      nsubm=nt(inode)%nsubm()
      ndof=nt(inode)%ndof()

      if (allocated(nt(inode)%eig).or.allocated(nt(inode)%assgn)) then
         call AbortWithError("SetEigenbasis(): basis already set")
      endif

      if (SIZE(eig).ne.neig) then
         write(*,*) 'Number of eigenvalues (',SIZE(eig),&
                    ') differs from expected number: ',neig
         call AbortWithError("SetEigenbasis(): wrong eig length")
      endif

      if (present(delta)) then
         if (SIZE(eig).ne.neig) then
            write(*,*) 'Number of deltas (',SIZE(delta),&
                       ') differs from expected number: ',neig
            call AbortWithError("SetEigenbasis(): wrong delta length")
         endif
      endif

      if (SIZE(qns,1).ne.neig .or. SIZE(qns,2).ne.max(nsubm,1)) then
         write(*,*) 'qns(',SIZE(qns,1),',',SIZE(qns,2),') dimensions ',&
         'differ from expected dimensions (',neig,',',nsubm,')'
         call AbortWithError("SetEigenbasis(): wrong qns length")
      endif

      allocate(nt(inode)%eig(neig),nt(inode)%delta(neig))
      allocate(nt(inode)%assgn(neig,ndof))
      nt(inode)%eig(:)=eig(:)

      if (present(delta)) then
         nt(inode)%delta(:)=delta(:)
      else
         nt(inode)%delta(:)=0.d0
      endif

      if (nsubm.eq.0) then ! Bottom layer assignments
         do k=1,neig     
            nt(inode)%assgn(k,1)=k-1
         enddo
      else
         j=0
         do i=1,nsubm
            sm=nt(inode)%subm(i)
            nsubdof=nt(sm)%ndof()
            do k=1,neig
               if (qns(k,i).gt.nt(sm)%nbas()) then
                  write(*,*) 'Eigenvector ',k,': item (',qns(k,i),'/',&
                  nt(sm)%nbas(),') of node ',sm,' exceeds basis'
                  call AbortWithError(&
                       'SetAssignment(): bad submode basis index')
               endif
               nt(inode)%assgn(k,j+1:j+nsubdof)=&
                  nt(sm)%assgn(qns(k,i),1:nsubdof)
            enddo
            j=j+nsubdof
         enddo
      endif

      end subroutine SetEigenbasis

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ConstructAssignment(nt,inode,qns,qnfull)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds assignment from qns, which point to basis functions in
! sub-nodes of inode

      implicit none
      TYPE (TN), intent(inout) :: nt(:)
      integer, intent(in) :: inode
      integer, intent(in) :: qns(:)
      integer, allocatable, intent(out) :: qnfull(:)
      integer :: i,j,sm,nsubm,ndof,nsubdof

      nsubm=nt(inode)%nsubm()
      ndof=nt(inode)%ndof()

      if (SIZE(qns).ne.max(nsubm,1)) then
         write(*,*) 'qns(',SIZE(qns),') dimensions ',&
         'differ from expected dimensions (',max(nsubm,1),')'
         call AbortWithError("SetEigenbasis(): wrong qns length")
      endif

      allocate(qnfull(ndof))

      if (nsubm.eq.0) then ! Bottom layer assignments
         qnfull(1)=qns(1)-1
      else
         j=0
         do i=1,nsubm
            sm=nt(inode)%subm(i)
            nsubdof=nt(sm)%ndof()
            if (qns(i).gt.nt(sm)%nbas()) then
               write(*,*) 'Basis item (',qns(i),'/',&
               nt(sm)%nbas(),') of node ',sm,' exceeds basis'
               call AbortWithError(&
                    'ConstructAssignment(): bad submode basis index')
            endif
            qnfull(j+1:j+nsubdof)=nt(sm)%assgn(qns(i),1:nsubdof)
            j=j+nsubdof
         enddo
      endif

      end subroutine ConstructAssignment

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildNodeTree(NT,modcomb,modstart,gdim,nmode,resort)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs node tree (from data from ML type, for now)

      implicit none
      TYPE (TN), allocatable :: nt(:)
      integer, intent(in) :: modcomb(:,:),modstart(:,:),gdim(:,:)
      integer, intent(in) :: nmode(:),resort(:)
      integer :: il,im,j,k,nnode,inode,lnode,nlayr,nsubm,ssubm,dofct,os

      nlayr=SIZE(nmode)

!     Count the total number of nodes
      nnode=0
      do il=1,nlayr
         do im=1,nmode(il)
            nnode=nnode+1
         enddo
      enddo

      ALLOCATE(nt(nnode))

!     Fill each node type
      inode=0
      lnode=0
      do il=1,nlayr
         do im=1,nmode(il)
            inode=inode+1

!           For cross-referencing ML and TN formats
            nt(inode)%nid=inode
            nt(inode)%mlil=il
            nt(inode)%mlim=im

!           Initialize to 0; overwrite when super node is processed
            nt(inode)%supern=0
            nt(inode)%superi=0

            call nt(inode)%setb(gdim(il,im))

!           Bottom layer: get the DOF from resort array
            if (il.eq.1) then
               ALLOCATE(nt(inode)%dofs(1))
               nt(inode)%dofs(1)=resort(im)

!           Upper layers
            else
               nsubm=modcomb(il,im)
               ssubm=lnode-nmode(il-1)+modstart(il,im)-1

!              Count bottom layer dofs contained in node
               dofct=0
               do k=1,nsubm
                  dofct=dofct+nt(ssubm+k)%ndof()
               enddo

!              Fill the node structure
               ALLOCATE(nt(inode)%subs(nsubm),nt(inode)%dofs(dofct))
               os=0
               do k=1,nsubm
                  dofct=nt(ssubm+k)%ndof()
                  nt(inode)%dofs(os+1:os+dofct)=nt(ssubm+k)%dofs(1:dofct)
                  nt(inode)%subs(k)=ssubm+k
                  nt(ssubm+k)%supern=inode
                  nt(ssubm+k)%superi=k
                  os=os+dofct
               enddo
            endif

         enddo
         lnode=inode ! index of node beginning this layer
      enddo

      end subroutine BuildNodeTree

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushNodeTree(nt)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs node tree (from data from ML type, for now)

      implicit none
      TYPE (TN), allocatable :: nt(:)
      integer :: i,nnodes

      IF (ALLOCATED(nt)) THEN
         nnodes=SIZE(nt)
         DO i=1,nnodes
            call nt(i)%flush()
         ENDDO
         DEALLOCATE(nt)
      ENDIF

      end subroutine FlushNodeTree

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TraverseNT_dofs(nt,ind)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prototype code to traverse sub-nodes of a node tree node

      implicit none
      CLASS (TN) :: nt(:)
      integer, intent(in) :: ind
      integer :: i,j,k,ns

      if (ind.lt.1 .or. ind.gt.SIZE(nt)) then
         write(*,'(2(A,I0),A)') 'ind = ',ind,'; must be in range [1,',&
                    SIZE(nt),']'
         call AbortWithError("TraverseNT_dofs: ind out of range")
      endif

      j=1
      i=ind
      k=0
      do
        ns=nt(i)%nsubm()
!        write(*,*) 'node',i,'has ',ns,' submodes; working on',j
        if (j.le.ns) then ! Descend in the tree
           i=nt(i)%subm(j)
           j=1
!           write(*,*) 'Migrating down to i=',i,'; j = ',j
        else ! Reached bottom; increment counter and backtrack
           if (ns.eq.0) k=k+1
           if (i.eq.ind) exit
           j=nt(i)%superi+1
           i=nt(i)%supern
!           write(*,*) 'Migrating   up to i=',i,'; j = ',j
        endif
      enddo

!      write(*,*) 'dofs found: ',k

      end subroutine TraverseNT_dofs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE NODETREE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
