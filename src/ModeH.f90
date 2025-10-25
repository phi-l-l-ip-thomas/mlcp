!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODEH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module builds the CP-format mode Hamiltonian

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE LINALG
      USE SEPDREPN
      USE HAMILSETUP
      USE INPUTCP
      USE MODECOMB
      USE MODVECVEC
      USE REDUCTION
      USE ALSDRVR

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_ModeH_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_ModeH_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_ModeH_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_ModeH_Module()
      call Get_MPI_Timings('ModeH module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_ModeH_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildModeHamiltonian(im,H,Ham,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs Hamiltonian in CP-format for mode 'im' in layer 'il'

      implicit none
      TYPE (CPpar) :: cpp
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (Configs), ALLOCATABLE :: pop(:),opcs(:,:)
      TYPE (CP), INTENT(OUT) :: H
      TYPE (CP) :: Hnew
      integer, intent(in)  :: im
      integer, allocatable :: nbas(:)
      integer :: i,j,nsubm,msubm,sm
      integer :: oldrank,tilesize
      real(kind=8) :: Hstor
      real(kind=8), parameter :: redtol=1.d-12
      logical      :: showFmG,tiledH
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeH_Module()

      call CPU_TIME(ti1)

!     Set parameters
      showFmG=.FALSE.
      tiledH=.FALSE.
      nsubm=Ham%nt(im)%nsubm()  ! Also equals nr. eigen terms
      msubm=max(nsubm,1)

      ALLOCATE(nbas(msubm))
      DO i=1,msubm
         sm=Ham%nt(im)%subm(i)
         if (nsubm.eq.0) then
            nbas(i)=Ham%nt(sm)%B
         else
            nbas(i)=Ham%nt(sm)%nbas()
         endif
      ENDDO

      IF (nsubm.ne.1) THEN
         IF (mpirank.eq.mpi_prnt_rank) &
         write(*,'(3X,A)') 'Building mode Hamiltonian...'
      ENDIF

!     Get the primitive and compound operator list
      call GetPrimOpList(Ham,im,pop,opcs,cpp%verbosity)

      IF (mpirank.eq.mpi_prnt_rank .and. cpp%verbosity.ge.2) THEN
         write(*,'(/X,A/)') 'Compound operator list, before sorting:'
         call ShowOPCS(opcs)
      ENDIF

!     Combine terms to reduce the rank of opcs "by hand"
      oldrank=SIZE(opcs,1)
      IF (cpp%h_sort_alg .seq. 'sort') THEN
         call CondenseHsort(opcs,.FALSE.)
      ELSEIF (cpp%h_sort_alg .seq. 'pack') THEN
         call CondenseHsort(opcs,.TRUE.)
      ELSEIF (cpp%h_sort_alg .seq. 'compare') THEN
         call CondenseHcompare(opcs)
      ELSE
         write(*,*) "Unrecognized Hamiltonian sort algorithm: '",&
                    TRIM(ADJUSTL(cpp%h_sort_alg)),&
                    "', must be either 'compare','sort', or 'pack'"
         call AbortWithError('BuildModeHamiltonian(): bad sort algo')
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank .and. cpp%verbosity.ge.2) THEN
         write(*,'(/X,A/)') 'Compound operator list, after sorting:'
         call ShowOPCS(opcs)
      ENDIF

!     Determine if H should be tiled over MPI ranks
      tiledH=(cpp%algo.ge.0 .and. cpp%lowmem.gt.2 .and. nsubm.ge.2 & 
              .and. .not.(nsubm.eq.2 .and. (cpp%red2D.seq.'SVD')) &
              .and. mpinodes.gt.1)

      IF (nsubm.ne.1) THEN

!        Calculate memory requirement
         Hstor=1.d0
         DO i=1,msubm
            Hstor=Hstor+nbas(i)*nbas(i)
         ENDDO

         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(/3X,A)') '*** ModeH memory usage ***'
            write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
                  oldrank,' to ',SIZE(opcs,1),' by combining like terms'
            write(*,'(7X,A,f12.6,A)') 'H memory (condensed): ',&
               Hstor*SIZE(opcs,1)/2**27,' GB'
            IF (tiledH) THEN
               tilesize=(SIZE(opcs,1)+mpinodes-1)/mpinodes
               write(*,'(7X,2(A,I0),A)') 'Hamiltonian tiled over ',&
                     mpinodes,' MPI processes, with ',tilesize,&
                     ' terms per process'
               write(*,'(7X,A,f12.6,A)') 'H memory     (tiled): ',&
                     Hstor*tilesize*3/2**27,&
                     ' GB (original + 2 copies)'
            ENDIF
            write(*,*)
         ENDIF

      ENDIF

      IF (cpp%ncycle.gt.0) THEN

!        Now assemble H in matrix representation from the list in opcp
         call GetHMats(im,H,opcs,nbas,pop,Ham,tiledH)

!        Additional reduction of H if desired
         oldrank=SIZE(H%coef)

         IF (cpp%hrank.gt.0 .and. cpp%hrank.lt.oldrank .and. &
             (.not.tiledH)) THEN
            call NORMBASE(H)
            call ordre(H)
!           Set the Hamiltonian reduction parameters first
            call SetReductionParameters(cpp%hrank,cpp%hnals,redtol,&
                 showFmG,cpp%red2D,cpp%redND,cpp%alspenalty,cpp%als_linsys_alg)
            Hnew=NewCP(cpp%hrank,H%rows,H%cols,H%sym)
            call GenCopyWtoV(Hnew,H,1,cpp%hrank,1,cpp%hrank)
            call reduc(Hnew,H)
            call ReplaceVwithW(H,Hnew)
         ENDIF

         IF (nsubm.ne.1 .and. SIZE(H%coef).lt.oldrank .and. &
            mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
                  oldrank,' to ',SIZE(H%coef),' by reduc()'
            write(*,'(7X,A,f12.6,A/)') 'H memory   (reduced): ',&
                  Hstor*SIZE(H%coef)/2**27,' GB'
         ENDIF

      ELSE
!        Allocate a "dummy" CP-vector with the correct size
!        since the solver needs this for the memory check
         oldrank=SIZE(opcs,1)
         IF (cpp%hrank.gt.0) oldrank=MIN(oldrank,cpp%hrank)
         H=NewCP(oldrank,nbas,.FALSE.)
         IF (nsubm.ne.1 .and. SIZE(H%coef).lt.SIZE(opcs,1) .and. &
            mpirank.eq.mpi_prnt_rank) THEN
            write(*,'(7X,2(A,I0),A)') 'Hamiltonian rank reduced from ',&
                  SIZE(opcs,1),' to ',SIZE(H%coef),' by reduc()'
            write(*,'(7X,A,f12.6,A/)') 'H memory (reduced)  : ',&
                  Hstor*SIZE(H%coef)/2**27,' GB'
         ENDIF

      ENDIF

      DEALLOCATE(nbas,opcs)

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine BuildModeHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPrimOpList(Ham,im,pop,opcs,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets list of primitive operators which apply to sub-nodes of 'inode'

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      INTEGER, INTENT(IN) :: im,verbosity
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: pop(:),opcs(:,:)
      TYPE (Configs) :: T,T2
      integer, allocatable :: nbas(:)
      integer :: ipass,i,j,k,l,nsubm,msubm,nHterm,nop,sm,nsubdof,rk,idx

      nsubm=Ham%nt(im)%nsubm()
      msubm=max(1,nsubm)
      nHterm=Ham%nt(im)%nHterm()

!     Generate the mode terms
      ALLOCATE(pop(msubm),opcs(nsubm+nHterm,msubm),nbas(1))
      nbas(:)=1
      DO j=1,msubm
         DO i=1,nsubm
            call NewConfigs(opcs(i,j),nbas,1)
            if (i.eq.j) opcs(i,j)%qns(1,1)=-1
            opcs(i,j)%coef(1)=1.d0
         ENDDO
      ENDDO
      DEALLOCATE(nbas)

      IF (nHterm.eq.0) RETURN

      DO ipass=1,2

!        Find the list of unique primitive operators
         DO j=1,msubm
            sm=Ham%nt(im)%subm(j)
            nsubdof=Ham%nt(sm)%ndof()
            ALLOCATE(nbas(2*nsubdof)) ! 2x for operator types
            nbas(:)=1

            IF (ipass.eq.1) THEN
               call NewConfigs(pop(j),nbas,1)
               pop(j)%coef(1)=1.d0
            ENDIF

            call NewConfigs(T,nbas,1)
            T%coef(1)=1.d0
            DEALLOCATE(nbas)

!           Loop over terms in H, find unique ones for this mode
            DO i=1,nHterm
               T%qns(1,:)=0
               nop=Ham%nt(im)%Hnop(i)
               rk=SIZE(pop(j)%coef)

!              Construct the Config representation of the operator
               DO k=1,nop
!                 Operator belongs to subnode j
                  IF (Ham%nt(im)%Hsubm(i,k).eq.j) THEN
                     DO l=1,nsubdof
                        IF (Ham%nt(im)%Hops(i,k,1).eq.Ham%nt(sm)%dofs(l)) THEN
                           T%qns(1,l)=Ham%nt(im)%Hops(i,k,2)
                           T%qns(1,l+nsubdof)=Ham%nt(im)%Hops(i,k,3)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

               idx=findconfigindex(pop(j),T%qns(1,:),(/1,rk/)) 
               IF (idx.eq.0) THEN

                  IF (ipass.eq.1) THEN
!                    First pass: increase the basis count (if a new 
!                    larger index is found), add operator to list of
!                    primitive operators, and resort list
                     DO l=1,2*nsubdof
                        IF (T%qns(1,l).gt.pop(j)%nbas(l)) &
                            pop(j)%nbas(l)=T%qns(1,l)
                     ENDDO

                     call NewConfigs(T2,pop(j)%nbas,rk+1)
                     call GenCopyConfigsWtoV(T2,pop(j),1,rk,1,rk)
                     call GenCopyConfigsWtoV(T2,T,rk+1,rk+1,1,1)
                     call SortConfigsByIndex(T2)
                     call ReplaceConfigsVwithW(pop(j),T2)
                  ELSE
!                    Second pass: error if operator not found
                     call AbortWithError(&
                          "GetPrimOpList(): operator not found")
                  ENDIF

               ELSE
!                 Second pass: record the operator in the list
                  IF (ipass.eq.2) THEN
                      call NewConfigs(opcs(i+nsubm,j),(/SIZE(pop(j)%coef)/),1)
                      if (j.eq.1) then
                         opcs(i+nsubm,j)%coef(1)=Ham%nt(im)%Hfacs(i)
                      else
                         opcs(i+nsubm,j)%coef(1)=1.d0
                      endif
                      opcs(i+nsubm,j)%qns(1,1)=idx
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.2) &
         call ShowPrimOpList(im,Ham,pop)

      end subroutine GetPrimOpList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowPrimOpList(im,Ham,pop)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank of OPCS by collecting common factors

      implicit none
      TYPE (Configs), INTENT(IN) :: pop(:)
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      INTEGER, INTENT(IN) :: im
      integer :: i,j,k,nsubm,nop,ndof,sm
      character*64 :: frmt
      character*3, dimension(:,:), allocatable :: tag

      nsubm=SIZE(pop)

      write(*,*)
      do i=1,nsubm
         sm=Ham%nt(im)%subm(i)
         nop=SIZE(pop(i)%qns,1)
!         ndof=SIZE(pop(i)%qns,2)/2
         ndof=Ham%nt(sm)%ndof()
         allocate(tag(ndof,2))
         write(*,'(X,A,X,I0/)') &
               'Primitive operators found for sub-mode:',i
         write(frmt,'(A,I0,A)') '(6X,A,',3*ndof-1,'X,A)'
         write(*,frmt) 'DOF','Type'
         write(frmt,'(A,I0,A)') '(6X,',ndof,'(I2,X))'
         write(*,frmt) (Ham%nt(sm)%dofs(k),k=1,ndof)
         write(frmt,'(A,I0,A,I0,A)') &
                    '(I4,A,X,',ndof,'A,A,X,',ndof,'A)'
         do j=1,nop
            do k=1,ndof
               if (pop(i)%qns(j,k) .gt. 0) then
                  write(tag(k,1),'(I2,X)') pop(i)%qns(j,k)
                  write(tag(k,2),'(I2,X)') pop(i)%qns(j,ndof+k)
               else
                  write(tag(k,1),'(3X)')
                  write(tag(k,2),'(3X)')
               endif
            enddo
            write(*,frmt) j,')',(tag(k,1),k=1,ndof),'|',&
                                (tag(k,2),k=1,ndof)
         enddo
         write(*,*)
         deallocate(tag)
      enddo

      end subroutine ShowPrimOpList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CondenseHsort(H,dopack)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank of OPCS by collecting common factors

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(INOUT) :: H(:,:)
      logical, intent(in) :: dopack
      integer :: rk,ipass

      ipass=0
      DO
        rk=SIZE(H,1)
        call SortOPCSmaster(H,(ipass.eq.0),dopack)
        IF (SIZE(H,1).eq.rk) EXIT
        ipass=ipass+1
        IF (mpirank.eq.mpi_prnt_rank) write(*,'(5X,3(A,I0))') &
           'sorting pass ',ipass,': ',rk,' -> ',SIZE(H,1)
      ENDDO

      end subroutine CondenseHsort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortOPCSmaster(opcs,firstpass,dopack)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank of OPCS by collecting common factors

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(INOUT) :: opcs(:,:)
      TYPE (Configs), ALLOCATABLE :: TS(:,:),TU(:,:)
      logical, intent(in) :: firstpass,dopack
      integer, allocatable :: nbas(:)
      integer :: i,i2,i3,j,k,l,mapk,ipass,ngroup,isz,nsolo,iref,radd,rtot
      integer :: nsubm,nHterm,mingrpsz,maxgrpsz
      logical :: same

      nHterm=SIZE(opcs,1)
      nsubm=SIZE(opcs,2)

      IF (nHterm.eq.1) RETURN

      allocate(TS(nHterm,nsubm))

      ngroup=0

      IF (firstpass) THEN
         mingrpsz=max(nHterm,2)
      ELSE
         mingrpsz=2
      ENDIF

!     Loop over group size thresholding
      DO

         IF (mingrpsz.lt.2) EXIT
         maxgrpsz=0

!         write(*,*) 'collecting terms with min group size of ',mingrpsz

         DO j=1,nsubm
!           Sort opcs hierarchically by terms-per-factor, operator ID, and coefs
            call SortOPCSouter(opcs,j,firstpass)

            DO ipass=1,2
!              1st ipass: count term groups exceeding size threshold
!              2nd ipass: copy: term groups above threshold -> TS;
!                               term groups below threshold -> TU
               nsolo=0
               iref=1
               isz=1

               DO i=2,SIZE(opcs,1)+1 ! +1 to close open group
               
                  same=(i.le.SIZE(opcs,1))
                  DO k=1,nsubm-1
                     if (.not.same) exit
!                    Sorted by mode j means that j is fastest, 
!                                              j-1 is next fastest, ...
!                    mapk below ensures that first non-matching index
!                    is reached as early as possible
                     mapk=mod(j-k+nsubm-1,nsubm)+1
                     same=CompareConfigs(opcs(i,mapk),opcs(iref,mapk))
                  ENDDO

                  IF (same) THEN ! Expand the existing group
                     isz=isz+1
                  ELSE           ! Start a new group
                     maxgrpsz=max(maxgrpsz,isz)

!                    Group size below thresh -> copy to TU on ipass=2
                     if (isz.lt.mingrpsz) then

                        if (ipass.eq.2) THEN
                           DO k=1,nsubm
                              DO i2=1,isz
                                 i3=iref+i2-1
                                 call CopyConfigsWtoV(TU(nsolo+i2,k),opcs(i3,k))
                                 call FlushConfigs(opcs(i3,k))
                              ENDDO
                           ENDDO
                        endif
                        nsolo=nsolo+isz

!                    Group size above thresh -> copy to TS on ipass=2
                     else

                        if (ipass.eq.2) then
                           ngroup=ngroup+1
!                          Consolidate group of configs into one entry of TS
                           DO k=1,nsubm
!                             Sum terms along mode j
                              IF (k.eq.j) THEN
                                 rtot=0
                                 DO i2=1,isz
                                    i3=iref+i2-1
                                    radd=SIZE(opcs(i3,k)%coef)
                                    rtot=rtot+radd
                                 ENDDO
                                 call NewConfigs(TS(ngroup,k),(/1/),rtot)
                                 rtot=0
                                 DO i2=1,isz
                                    i3=iref+i2-1
                                    radd=SIZE(opcs(i3,k)%coef)
                                    call GenCopyConfigsWtoV(TS(ngroup,k),&
                                         opcs(i3,k),rtot+1,rtot+radd,1,radd)
                                    call FlushConfigs(opcs(i3,k))
                                    rtot=rtot+radd
                                 ENDDO
!                             Non-j modes all share reference config
                              ELSE
                                 call CopyConfigsWtoV(TS(ngroup,k),opcs(iref,k))
                                 call FlushConfigs(opcs(iref,k))
                              ENDIF

                           ENDDO
                        endif
                        
                     endif
                     iref=i
                     isz=1
                  ENDIF
               ENDDO

               IF (ipass.eq.1) THEN
                  if (nsolo.gt.0) allocate(TU(nsolo,nsubm))
               ENDIF
            ENDDO

!           If no terms go to TU, then sorting is complete!
            IF (nsolo.eq.0) EXIT

!           Recycle terms in TU for sorting by next submode
            deallocate(opcs)
            allocate(opcs(SIZE(TU,1),nsubm))
            do k=1,nsubm
               do i=1,SIZE(TU,1)
                  call CopyConfigsWtoV(opcs(i,k),TU(i,k))
               enddo
            enddo
            deallocate(TU)

         ENDDO ! j

!        Sorting complete, so exit thresholding loop
         IF (nsolo.eq.0) EXIT

!         write(*,*) 'max group size found was ',maxgrpsz
         IF (dopack) THEN
            mingrpsz=min(mingrpsz-1,maxgrpsz)
         ELSE
            mingrpsz=min((mingrpsz+1)/2,maxgrpsz)
         ENDIF
      ENDDO ! Group size thresholding

!     Copy any leftover terms in opcs to TS
      DO k=1,nsubm
         DO i=1,nsolo
            call CopyConfigsWtoV(TS(ngroup+i,k),opcs(i,k))
            call FlushConfigs(opcs(i,k))
         ENDDO
      ENDDO
      ngroup=ngroup+nsolo
      DEALLOCATE(opcs)

!     Replace opcs with grouped terms in TS
      ALLOCATE(opcs(ngroup,nsubm))
      DO k=1,nsubm
         DO i=1,ngroup
            call CopyConfigsWtoV(opcs(i,k),TS(i,k))
         ENDDO
      ENDDO

      end subroutine SortOPCSmaster

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortOPCSouter(opcs,ifast,firstpass)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts opcs hierarchically, with mode 'ifast' iterating fastest

      implicit none
      TYPE (Configs), INTENT(INOUT) :: opcs(:,:)
      integer, intent(in) :: ifast
      logical, intent(in) :: firstpass
      integer, allocatable :: ntrm(:,:)
      TYPE (Configs), allocatable :: optmp(:)
      integer, allocatable :: key(:),kx(:),ist(:),iend(:)
      integer :: i,j,k,nsubm,nHterm
      real(kind=8) ::fac
      character*64 :: frmt

      nHterm=SIZE(opcs,1)
      nsubm=SIZE(opcs,2)

      IF (nHterm.eq.1) RETURN

!     Build the key and array with the term counts
      ALLOCATE(key(nHterm),ntrm(nHterm,nsubm+1))
      DO i=1,nHterm
         key(i)=i
         ntrm(i,nsubm+1)=0
         DO j=1,nsubm
            k=mod(ifast+j-1,nsubm)+1
            ntrm(i,j)=SIZE(opcs(i,k)%coef)
!           Move the coefficient to the fastest iterating mode
            IF (j.lt.nsubm) THEN
               fac=opcs(i,k)%coef(1)
               opcs(i,ifast)%coef(:)=opcs(i,ifast)%coef(:)*fac
               opcs(i,k)%coef(:)=opcs(i,k)%coef(:)/fac
            ENDIF
         ENDDO
      ENDDO

!     On the first pass ntrm=1 for all submodes, so just sort once
      IF (firstpass) THEN
         call SortOPCSmiddle(key,opcs,ifast,1,nHterm)

!     Perform the hierarchical sort by number of configs per term
      ELSE
         ALLOCATE(ist(nsubm+1),iend(nsubm+1))
         iend=0
         iend(1)=nHterm
         ist(1)=1
         j=1
         DO
         
            IF (j.le.nsubm) THEN
               k=mod(ifast+j-1,nsubm)+1
!              Sort configurations by number of terms for (mapped) mode j
               kx=getsortkey(ntrm(ist(j):iend(j),j))
               call sortbykey(ntrm(ist(j):iend(j),:),kx)
               call sortbykey(key(ist(j):iend(j)),kx)
               deallocate(kx)
            ELSE
!              Once j exceeds nsubm, sort block of terms by sub-terms
               call SortOPCSmiddle(key,opcs,ifast,ist(j),iend(j))
            ENDIF

!           Update DOF index j
            IF (j.lt.nsubm+1) j=j+1
            DO
               IF (j.eq.1) EXIT
               IF (iend(j).lt.iend(j-1)) EXIT
               j=j-1
            ENDDO

!           When j returns to 1, the entire vector is sorted
            IF (j.eq.1) EXIT

!           Update the sort ranges
            ist(j)=iend(j)+1
            iend(j)=ibisect(ntrm(ist(j):iend(j-1),j-1),1)+ist(j)-1
         ENDDO
         DEALLOCATE(ist,iend)

      ENDIF

!     Rearrange items in opcs in order of the key
      ALLOCATE(optmp(nHterm))
      DO j=1,nsubm
         DO i=1,nHterm
            IF (key(i).ne.i) THEN
               call CopyConfigsWtoV(optmp(i),opcs(i,j))
               call FlushConfigs(opcs(i,j))
            ENDIF
         ENDDO

         DO i=1,nHterm
            IF (key(i).ne.i) THEN
               call CopyConfigsWtoV(opcs(i,j),optmp(key(i)))
               call FlushConfigs(optmp(key(i)))
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(key,ntrm)

      end subroutine SortOPCSouter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortOPCSmiddle(key,opcs,ifast,ist,iend)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rearranges the sort key for opcs in ascending order of number of terms

      implicit none
      TYPE (Configs), INTENT(IN) :: opcs(:,:)
      integer, intent(inout) :: key(:)
      integer, intent(in) :: ist,iend,ifast
      integer, allocatable :: cfg(:)
      real(kind=8), allocatable :: arr(:,:)
      integer :: i,i2,j,j2,k,l,os,aw,nHterm,nsubm

      IF (iend-ist.eq.0) RETURN
      
      nHterm=iend-ist+1
      nsubm=SIZE(opcs,2)
     
!     Each item in the block of opcs passed to this subroutine has the
!     same number of terms, so just check the first for the count
      allocate(cfg(nsubm))
      aw=0
      DO j=1,nsubm
         k=mod(ifast+j-1,nsubm)+1
         cfg(j)=SIZE(opcs(key(ist),k)%coef)
         aw=aw+cfg(j)
      ENDDO

!     Build the table to be sorted
      allocate(arr(nHterm,2*aw))
      do i=1,nHterm
         i2=ist+i-1
         os=0
         do j=1,nsubm
            k=mod(ifast+j-1,nsubm)+1
            do l=1,cfg(j)
               j2=os+l
               arr(i,j2)=opcs(key(i2),k)%qns(l,1)
               arr(i,j2+aw)=opcs(key(i2),k)%coef(l)
            enddo
            os=os+cfg(j)
         enddo
      enddo

      call SortOPCSinner(key(ist:iend),arr)

      deallocate(cfg,arr)

      end subroutine SortOPCSmiddle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortOPCSmiddleA(key,opcs,ifast)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rearranges the sort key for opcs in ascending order of number of terms

      implicit none
      TYPE (Configs), INTENT(IN) :: opcs(:,:)
      integer, intent(inout) :: key(:)
      integer, intent(in) :: ifast
      integer, allocatable :: cfg(:)
      real(kind=8), allocatable :: arr(:,:)
      integer :: i,i2,j,j2,k,l,os,aw,nHterm,nsubm

      nHterm=SIZE(opcs,1)
      nsubm=SIZE(opcs,2)

      IF (nHterm.eq.0) RETURN
      
!     Find the max number of terms for each submode
      allocate(cfg(nsubm))
      cfg(:)=0
      aw=0
      DO j=1,nsubm
         k=mod(ifast+j-1,nsubm)+1
         DO i=1,nHterm
            cfg(j)=MAX(SIZE(opcs(i,k)%coef),cfg(j))
         ENDDO
         aw=aw+cfg(j)
      ENDDO

!     Build the table to be sorted
      allocate(arr(nHterm,2*aw))
      arr(:,:)=0.d0
      do i=1,nHterm
         os=0
         do j=1,nsubm
            os=os+cfg(j)
            k=mod(ifast+j-1,nsubm)+1
            do l=1,SIZE(opcs(i,k)%coef)
               j2=os+1-l
               arr(i,j2)=opcs(i,k)%qns(l,1)
               arr(i,j2+aw)=opcs(i,k)%coef(l)
            enddo
         enddo
      enddo

      call SortOPCSinner(key,arr)

      deallocate(cfg,arr)

      end subroutine SortOPCSmiddleA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortOPCSinner(key,arr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rearranges the sort key for opcs in ascending order of operator ids
! and coefs

      implicit none
      integer, intent(inout) :: key(:)
      real(kind=8), intent(inout) :: arr(:,:)
      integer, allocatable :: ist(:),iend(:),kx(:)
      integer :: kl,kw,j,i

      kl=SIZE(arr,1)
      kw=SIZE(arr,2)
      
      IF (kl.eq.1) RETURN

      ALLOCATE(ist(kw),iend(kw))
      iend=0
      iend(1)=kl
      ist(1)=1
      j=1
      DO
!        Sort the configurations by index j
         kx=getsortkey(arr(ist(j):iend(j),j))
         call sortbykey(arr(ist(j):iend(j),:),kx)
         call sortbykey(key(ist(j):iend(j)),kx)
         deallocate(kx)

!        Update DOF index j
         IF (j.lt.kw) j=j+1
         DO
            IF (j.eq.1) EXIT
            IF (iend(j).lt.iend(j-1)) EXIT
            j=j-1
         ENDDO

!        When j returns to 1, the entire array is sorted
         IF (j.eq.1) EXIT

!        Update the sort ranges
         ist(j)=iend(j)+1
         iend(j)=rbisectH(arr(ist(j):iend(j-1),j-1),&
                          arr(ist(j),j-1))+ist(j)-1
      ENDDO

      DEALLOCATE(ist,iend)

      end subroutine SortOPCSinner

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowOPCS(opcs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank of OPCS by collecting common factors

      implicit none
      TYPE (Configs), intent(in) :: opcs(:,:)
      integer :: i,j,nsubm,nrk

      nrk=SIZE(opcs,1)
      nsubm=SIZE(opcs,2)

      do i=1,nrk
         do j=1,nsubm
            write(*,'(2(X,A,X,I0),A)') 'Term',i,'(sub-mode',j,'):'
            call PrintConfigs(opcs(i,j))
         enddo
      enddo

      end subroutine ShowOPCS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CondenseHcompare(H)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Does "by-hand" reduction of H by combining terms with common factors

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(INOUT) :: H(:,:)
      TYPE (Configs), ALLOCATABLE :: T(:,:)
      INTEGER :: i,j,htrm,ttrm,nsubm,iadd,hrk,trk
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
                  trk=SIZE(T(j,iadd)%coef)
                  hrk=SIZE(H(i,iadd)%coef)
                  call ResizeConfigList(T(j,iadd),trk+hrk)
                  call GenCopyConfigsWtoV(T(j,iadd),H(i,iadd),trk+1,trk+hrk,1,hrk)
                  EXIT
               ENDIF
            ENDDO

!           If term in H cannot be summed with any term in tmp,
!           increase the rank of tmp by adding the H-term to the end
            IF (.NOT.add) THEN
               ttrm=ttrm+1
               DO j=1,nsubm
                  call CopyConfigsWtoV(T(ttrm,j),H(i,j))
               ENDDO
            ENDIF
         ENDDO

!        Replace H <-- tmp
         DEALLOCATE(H)
         ALLOCATE(H(ttrm,nsubm))
         DO i=1,ttrm
            DO j=1,nsubm
               call ReplaceConfigsVwithW(H(i,j),T(i,j))
            ENDDO
         ENDDO

!        If the size of HT cannot be reduced further, exit
         IF (ttrm.ge.htrm) EXIT
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(5X,3(A,I0))') &
               'comparing pass ',pass+1,': ',htrm,' -> ',ttrm
         htrm=ttrm
         pass=pass+1
      ENDDO

      end subroutine CondenseHcompare

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine compareHT(H,T,add,iadd)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares Hamiltonian terms by their primitive operator matrices
! (indexed in H%base and T%base) and coefficients, in order to determine
! if the terms can be condensed. The coefs of H and T may be modified.

      implicit none
      TYPE (Configs), INTENT(INOUT) :: H(:),T(:)
      INTEGER, INTENT(OUT) :: iadd
      LOGICAL, INTENT(OUT) :: add
      LOGICAL, ALLOCATABLE :: sameb(:)
      INTEGER :: i,j,k,l,nsubm,nonm,inm(2),nrkH,nrkT
      LOGICAL :: allfound,found,allcmatch,cmatch
      real(kind=8)  :: fac
      real(kind=8), PARAMETER  :: tol=1.d-15

      nsubm=SIZE(H)

!     Start by assuming terms are different in base and coef
      ALLOCATE(sameb(nsubm))
      sameb=.FALSE.

!     Loop over sub-modes and compare the bases/coefs
      nonm=0
      DO i=1,nsubm
         nrkH=SIZE(H(i)%coef)
         nrkT=SIZE(T(i)%coef)

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
!                 If the base matches, compare the coefficients
                  IF (H(i)%qns(j,1).eq.T(i)%qns(k,1)) THEN
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

      subroutine GetHMats(im,H,opcs,nbas,pop,Ham,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds H in CP-matrix form by multiplying and adding operator matrices
! The operator IDs are stored in 'poplist'
! The sum-of-products scheme for the operators is stored in 'opcp'

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      TYPE (CP), INTENT(OUT) :: H
      TYPE (Configs), INTENT(IN)  :: opcs(:,:),pop(:)
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: im
      LOGICAL, INTENT(IN) :: tiledH
      real(kind=8), ALLOCATABLE :: tvec(:),tmat(:,:),tmat1(:,:),tmat2(:,:)
      INTEGER :: ii,i,j,k,l,gst,ntrm,nsubm,nbasop,nrkop,primop
      INTEGER :: jdof,jpow,jopt,sm
      INTEGER :: ntiles,trmsizenode,trmthisnode,trmexcess,itrm,ftrm,irnk
      real(kind=8), PARAMETER :: smallnr=1.d-15

!     Set parameters
      ntrm=SIZE(opcs,1)
      nsubm=SIZE(opcs,2)

      IF (tiledH) THEN
         ntiles=mpinodes
         irnk=mpirank
      ELSE
         ntiles=1
         irnk=0
      ENDIF
      trmsizenode=(ntrm+ntiles-1)/ntiles
      trmthisnode=ntrm/ntiles
      trmexcess=mod(ntrm,ntiles)
      if (irnk.lt.trmexcess) trmthisnode=trmthisnode+1
      itrm=min(irnk,trmexcess)*trmsizenode+&
           max(irnk-trmexcess,0)*trmthisnode+1
      ftrm=itrm+trmthisnode-1

      H=NewCP(trmsizenode,nbas,.FALSE.)
      H%coef=1.d0
      ii=0
      DO i=itrm,ftrm
         gst=0
         ii=ii+1
         DO j=1,nsubm
            nrkop=SIZE(opcs(i,j)%coef)
            nbasop=SIZE(pop(j)%qns,2)/2
            sm=Ham%nt(im)%subm(j)
            ALLOCATE(tmat(nbas(j),nbas(j)))
            tmat=0.d0
            
            DO k=1,nrkop
               primop=opcs(i,j)%qns(k,1)

!              Pre-solved mode operator (list of eigenvalues)
               IF (primop.eq.-1) THEN
                  call GetIdentityMatrix(tmat1,nbas(j),.FALSE.)
                  DO l=1,nbas(j)
                     tmat1(l,l)=Ham%nt(sm)%eig(l)*opcs(i,j)%coef(k)
                  ENDDO

               ELSE
!                 Start building the operator with an identity matrix
!                 scaled by the operator coefficient
                  call GetIdentityMatrix(tmat1,nbas(j),.FALSE.)
                  DO l=1,nbas(j)
                     tmat1(l,l)=tmat1(l,l)*opcs(i,j)%coef(k)
                  ENDDO

!                 Construct the primitive operator product
                  IF (primop.gt.0) THEN ! primop=0 is identity

!                    Multiply the intra-sub-mode product operators
                     DO l=1,nbasop
                        jdof=Ham%nt(sm)%dofs(l)
                        jpow=pop(j)%qns(primop,l)
                        jopt=pop(j)%qns(primop,l+nbasop)
                          
                        IF (pop(j)%qns(primop,l).gt.0) THEN
                           call Vec2SymPackMat(Ham%ops(jdof,jpow,jopt)%mat,tmat2)
!                          tmat1 = tmat1 * tmat2
                           call MatrixMult(tmat1,.FALSE.,tmat2,.FALSE.)
                           DEALLOCATE(tmat2)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

!              Add tmat1 to the sum
               tmat=tmat+tmat1
               DEALLOCATE(tmat1)

            ENDDO

!           Symmetrize to scrub non-Hermitian errors due to truncating
!           basis as one moves up the tree
            call SymmetrizeMat(tmat)

!           Put the sum-of-products term into H
            call Mat2Vec(tvec,tmat,.FALSE.)

            H%base(gst+1:gst+H%nbas(j),ii)=tvec
            DEALLOCATE(tmat,tvec)
            gst=gst+H%nbas(j)
         ENDDO
      ENDDO

!     If the number of MPI ranks do not evenly divide out the number of
!     terms in H, some ranks will have one fewer terms. For these ranks,
!     add near-zero terms of opposite signs to avoid zero division in
!     the ALS solver
      IF (ii.lt.trmsizenode) THEN
         ii=ii+1
         H%coef(ii)=smallnr
         gst=0
         DO j=1,nsubm
            call GetIdentityMatrix(tmat,nbas(j),.FALSE.)
            call Mat2Vec(tvec,tmat,.FALSE.)
            if (j.eq.1 .and. mod(mpirank,2).eq.1) tvec=-tvec
            H%base(gst+1:gst+H%nbas(j),ii)=tvec
            DEALLOCATE(tmat,tvec)
            gst=gst+H%nbas(j)
         ENDDO
      ENDIF

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
      real(kind=8), allocatable :: tmat(:,:),tmat1(:,:)
      real(kind=8), allocatable :: tvec(:)
      character(len=64), parameter :: solver='LU'
      integer, allocatable :: nbas(:)
      integer :: j,ndof,gi,gf

      ndof=SIZE(H%nbas)
      ALLOCATE(nbas(ndof))

      DO j=1,ndof
         nbas(j)=NINT(sqrt(REAL(H%nbas(j))))
      ENDDO

!     Initial guess: reduce Hi <- H, with Hi rank-1
      call SetReductionParameters(1,30,1.d-12,.FALSE.,'SVD','SR1',&
                                  1.d-10,solver)
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
      real(kind=8), intent(in)   :: E
      logical, allocatable :: symm(:)
      integer, allocatable :: nbas(:)
      real(kind=8), allocatable  :: tmat(:,:),tvec(:)
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
      real(kind=8), allocatable :: M(:,:),v(:)
      integer :: i,j,j1,j2,gst,ndof,nrk

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
