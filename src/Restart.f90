!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE RESTART

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Restart a crashed calculation (provided it doesn't crash too hard!)

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MODECOMB
      USE INPUTCP
      USE HAMILSETUP
      USE SEPDREPN
      USE CPr8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RestartSetup(Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master routine for reading restart files

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer :: i
      logical :: success
      logical, allocatable :: s1(:),s2(:)

      if (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A)') 'Restart setup...'

!     Check for restart files to determine if this is a restart run
      call CheckForRestartFile(ML)

      IF (ML%dorestart) THEN
         call ValidateRestart(ML,s1,s2)

!        Read the file containing eigenvalues from finished nodes
         call ReadEigenvalues(ML,Ham,1,success)
         IF (.not.success) THEN
            call ReadEigenvalues(ML,Ham,2,success)
            IF (.not.success) call AbortWithError(&
               "Eigenvalue restart file could not be read")
         ENDIF

!        Read the file containing operator matrices
         call ReadOperMats(ML,Ham,1,success)
         IF (.not.success) THEN
            call ReadOperMats(ML,Ham,2,success)
            IF (.not.success) call AbortWithError(&
               "Operator matrix restart file could not be read")
         ENDIF

!        Make sure parameters have not changed for finished nodes
         DO i=1,SIZE(Ham%nt)
            IF (Ham%nt(i)%done .and. ((.not.s1(i)).or.(.not.s2(i)))) &
            call AbortWithError(&
            'RestartSetup(): finished node settings cannot be changed!')
         ENDDO
         DEALLOCATE(s1,s2)
      ENDIF

!     Save 'CP.inp' and 'layers.inp' unless restart file is 'none'
      IF (ML%resfile(1:4).seq.'none') RETURN

      call SaveModeDat(ML)

      end subroutine RestartSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CheckForRestartFile(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      character(len=128) :: fnm_layers
      logical :: found

!     Look for _layers.rst. If present, set dorestart=.TRUE.
      write(fnm_layers,'(2A)') TRIM(ADJUSTL(ML%resfile)),'_layers.rst'

      IF (mpirank.eq.mpi_io_rank) &
         INQUIRE(FILE=TRIM(ADJUSTL(fnm_layers)), EXIST=found)
     
      call bcast(found,mpi_io_rank)
      ML%dorestart=found

      IF (mpirank.eq.mpi_prnt_rank) THEN
         IF (found) THEN
            write(*,'(/5X,2A)') TRIM(ADJUSTL(fnm_layers)),' exists'
            write(*,'(5X,A/)') &
            'This job is a restart. Validating restart data...'
         ENDIF
      ENDIF

      end subroutine CheckForRestartFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateRestart(ML,s1,s2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares restart data with job parameters to determine if calculation
! can be successfully restarted

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (MLtree) :: MLrst
      TYPE (CPpar)  :: cpp,cpprst
      character(len=128) :: fnm
      integer :: il,im,j,nnode
      logical :: same1,same2
      logical, allocatable, intent(out) :: s1(:),s2(:)

!     Read the restart input file
      write(fnm,'(2A)') TRIM(ADJUSTL(ML%resfile)),'_layers.rst'
      call ReadModeDat(MLrst,fnm)
      call BcastModeDat(MLrst)

!     Validate input
      IF (ML%system.ne.MLrst%system) THEN
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,*) 'Old system: ',MLrst%system,&
                     '; New system: ',ML%system
         call AbortWithError(&
         'ValidateRestart(): System change not allowed in a restart!')
      ENDIF

      IF (ML%dividefc.ne.MLrst%dividefc) THEN
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,*) 'Old dividefc: ',MLrst%dividefc,&
                     '; New dividefc: ',ML%dividefc
         call AbortWithError(&
         'ValidateRestart(): dividefc must be same in restart!')
      ENDIF

      IF (ML%pe_transform.ne.MLrst%pe_transform) THEN
         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,*) 'Old pe_transform: ',MLrst%pe_transform,&
                     '; New pe_transform: ',ML%pe_transform
         call AbortWithError(&
         'ValidateRestart(): pe_transform must be same in restart!')
      ENDIF

      IF (.not.(ML%pe_transform.seq. 'none')) THEN
         IF (ML%pe_trans_fac.ne.MLrst%pe_trans_fac) THEN
            IF (mpirank.eq.mpi_prnt_rank) &
               write(*,*) 'Old pe_trans_fac: ',MLrst%pe_trans_fac,&
                        '; New pe_trans_fac: ',ML%pe_trans_fac
            call AbortWithError(&
            'ValidateRestart(): pe_trans_fac must be same in restart!')
         ENDIF
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank) THEN
         IF (.not.ALL(ML%rs.eq.MLrst%rs)) write(*,'(X,2(A,33(I0,X)))') &
            ' * rs changed from ',(ML%rs(j),j=1,33),&
                           ' to ',(MLrst%rs(j),j=1,33)         
      ENDIF

!     Validate input against layers.rst. Make sure that the ordering of
!     DOFs and the tree structure are the same. The basis sizes are
!     checked later
      IF (ML%ndof.ne.MLrst%ndof) &
         call AbortWithError("ValidateRestart(): ndof do not match")
      IF (ML%nlayr.ne.MLrst%nlayr) &
         call AbortWithError("ValidateRestart(): nlayr do not match")

      DO im=1,ML%ndof
         IF (ML%resort(im).ne.MLrst%resort(im)) call &
            AbortWithError("ValidateRestart(): $resort do not match")
      ENDDO

      DO il=1,ML%nlayr
         IF (ML%nmode(il).ne.MLrst%nmode(il)) call &
            AbortWithError("ValidateRestart(): $layers do not match")
         DO im=1,ML%nmode(il)
            IF (ML%modcomb(il,im).ne.MLrst%modcomb(il,im)) call &
               AbortWithError("ValidateRestart(): $layers do not match")
         ENDDO
      ENDDO

      DO il=2,ML%nlayr
         DO im=1,ML%nmode(il)
            IF (ML%gdim(il,im).ne.MLrst%gdim(il,im)) THEN
               IF (mpirank.eq.mpi_prnt_rank) write(*,'(X,4(A,I0),A)') &
                 'layer-mode: ',il,'-',im,': input (',ML%gdim(il,im),&
                 ') and restart (',MLrst%gdim(il,im),') basis mismatch'
               call AbortWithError(&
                    'ValidateRestart(): $basis do not match')
            ENDIF
         ENDDO
      ENDDO

!     Check if CPpar settings changed. Default settings must all match
      call compareCPP(ML%dpp,MLrst%dpp,same1,same2)
      IF (.not.same1 .or. .not.same2) THEN
         IF (mpirank.eq.mpi_prnt_rank) write(*,'(2X,A)') &
            'Changes found to parameters set in $control:'
         call showCPPcomparisons(ML%dpp,MLrst%dpp)
         call AbortWithError(&
             'ValidateRestart(): restart $control settings must match')
      ENDIF

!     CPpar settings for individual nodes *can* change if node has not
!     yet been processed, and, in certain cases, if node is in progress
      nnode=0
      DO il=1,ML%nlayr
         nnode=nnode+ML%nmode(il)
      ENDDO

      ALLOCATE(s1(nnode),s2(nnode))
      j=0
      DO il=1,ML%nlayr
         DO im=1,ML%nmode(il)
            j=j+1
            cpp=getnodecppar(ML,j)
            cpprst=getnodecppar(MLrst,j)
            call compareCPP(cpp,cpprst,s1(j),s2(j))
            ! Record which node has first instance of not being same as
            ! reference, using 'same1' and 'same2' criteria
            IF (.not.s1(j) .or. .not.s2(j)) THEN
               IF (mpirank.eq.mpi_prnt_rank) THEN
                  write(*,'(2X,2(A,I0),A)') &
                    'Changes found to parameters for layer-mode ',&
                    il,'-',im,':'
                  call showCPPcomparisons(cpp,cpprst)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      call Flush_ModeComb(MLrst)

      end subroutine ValidateRestart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveEigenInfo(ML,Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Saves eigenvalues, assignments, and operator matrices

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham

!     If restart file is 'none', exit without saving
      IF (ML%resfile(1:4).seq.'none') RETURN

!     Save data twice to prevent a potential corrupted file write
      call SaveEigenvalues(ML,Ham,1)
      call SaveEigenvalues(ML,Ham,2)
      call SaveOperMats(ML,Ham,1)
      call SaveOperMats(ML,Ham,2)

      end subroutine SaveEigenInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadOperMats(ML,Ham,nr,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of operator matrices to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in)  :: nr
      character(len=128)    :: fnm
      logical, intent(out) :: success
      logical :: successitems(4)
      integer :: u,i,j,k,l,ndof,npow,nopt,n,InpStat

      successitems(:)=.FALSE.
      write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
              '_',nr,'_oper.rst'

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Read file containing the operator matrices
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
         successitems(1)=(InpStat.eq.0)
        
         if (successitems(1)) then

!           Read the number of matrices
            read(u,*)
            read(u,*,IOSTAT=InpStat) ndof,npow,nopt
            successitems(2)=(InpStat.eq.0)
            read(u,*)

!           Make sure that the number of operator matrices in the restart
!           file matches the number allocated in the current run
            if (successitems(2)) then
               successitems(3)=((ndof.eq.SIZE(Ham%ops,1)) .and. &
                                (npow.eq.SIZE(Ham%ops,2)) .and. &
                                (nopt.eq.SIZE(Ham%ops,3)))

               if (successitems(3)) then

!                 Cycle through the primitive operator matrices
                  do l=1,nopt
                     do j=1,npow
                        do i=1,ndof
                           if (Ham%optable(i,j,l)) then
                              read(u,*,IOSTAT=InpStat) Ham%ops(i,j,l)%dof,n
                              successitems(4)=(InpStat.eq.0)
                              if (.not.successitems(4)) exit

!                             The operator matrices are allocated and written when the 1D
!                             (bottom layer) problems are solved, so these must be 
!                             deallocated and reallocated in a restarted run
                              IF (allocated(Ham%ops(i,j,l)%mat)) &
                                 DEALLOCATE(Ham%ops(i,j,l)%mat)
                              ALLOCATE(Ham%ops(i,j,l)%mat(n))

!                             Read the matrix for each operator
                              read(u,*,IOSTAT=InpStat) (Ham%ops(i,j,l)%mat(k),k=1,n)
                              successitems(4)=(InpStat.eq.0)
                              if (.not.successitems(4)) exit
                           endif
                        enddo
                     enddo
                  enddo
               endif ! successitems(3)
            endif ! successitems(2)
         endif ! successitems(1)
      
         close(u)

      ENDIF rank0

      call bcast(ndof,mpi_io_rank)
      call bcast(npow,mpi_io_rank)
      call bcast(nopt,mpi_io_rank)
      call bcast(successitems,mpi_io_rank)
      success=ALL(successitems(:))

      IF (mpirank.eq.mpi_prnt_rank) THEN

         if (.not.successitems(1)) then
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         elseif (.not.successitems(2)) then
            write(*,*) 'Number of primitive operators could not be read' 
         elseif (.not.successitems(3)) then
            write(*,*) 'Dimensions of operator array, current run : ',&
                       SIZE(Ham%ops,1),' x ',SIZE(Ham%ops,2),' x ',SIZE(Ham%ops,3)
            write(*,*) 'Dimensions of operator array, restart file: ',&
                       ndof,' x ',npow,' x ',nopt
            write(*,*) 'Wrong number of operators in restart file'
         elseif (.not.successitems(4)) then
            write(*,*) 'Primitive operator matrices could not be read'
         else
            write(*,'(5X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'
         endif

      ENDIF

!     Broadcast operator matrices if read succeeded
      if (success) call BcastOperMats(Ham)

      end subroutine ReadOperMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastOperMats(Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts Operator matrices from MPI rank used to read to other ranks
        
      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, allocatable :: ndim(:,:,:)
      integer :: i,j,k,ndof,npow,nopt,n

      ndof=SIZE(Ham%ops,1)
      npow=SIZE(Ham%ops,2)
      nopt=SIZE(Ham%ops,3)

      ALLOCATE(ndim(ndof,npow,nopt))
      DO k=1,nopt
         DO j=1,npow
            DO i=1,ndof
               ndim(i,j,k)=SIZE(Ham%ops(i,j,k)%mat)
            ENDDO
         ENDDO
      ENDDO

      call bcast(ndim,mpi_io_rank)

      DO k=1,nopt
         DO j=1,npow
            DO i=1,ndof
               IF (Ham%optable(i,j,k)) THEN
                  n=ndim(i,j,k)
                  Ham%ops(i,j,k)%dof=i

                  IF (mpirank.ne.mpi_io_rank) THEN
                     IF (allocated(Ham%ops(i,j,k)%mat)) &
                        DEALLOCATE(Ham%ops(i,j,k)%mat)
                        ALLOCATE(Ham%ops(i,j,k)%mat(n))
                  ENDIF

                  call bcast(Ham%ops(i,j,k)%mat,mpi_io_rank)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(ndim)

      end subroutine BcastOperMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveOperMats(ML,Ham,nr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of operator matrices to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in) :: nr
      character(len=128)  :: fnm,frmt
      integer :: u,i,j,k,l,n,ndof,npow,nopt

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Open output file
         write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
                 '_',nr,'_oper.rst'
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!        Record the number of matrices
         ndof=SIZE(Ham%ops,1)
         npow=SIZE(Ham%ops,2)
         nopt=SIZE(Ham%ops,3)
         write(u,'(A)') 'Number of operators:'
         write(u,'(I0,X,I0,X,I0)') ndof,npow,nopt
         write(u,'(A)') 'DOF# / n'

!        Write each operator matrix to file
         DO l=1,nopt
            DO j=1,npow
               DO i=1,ndof
                  IF (Ham%optable(i,j,l)) THEN
                     n=SIZE(Ham%ops(i,j,l)%mat)
                     write(u,'(4(I0,X))') Ham%ops(i,j,l)%dof,n
                     write(frmt,'(A,I0,A)') '(',n,'(E23.16,X))'
                     write(u,frmt) (Ham%ops(i,j,l)%mat(k),k=1,n)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         write(u,*)
         close(u)

      ENDIF rank0

      end subroutine SaveOperMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadEigenvalues(ML,Ham,nr,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads list of eigenvalues/assignments/deltas from restart file

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(in)  :: nr
      logical, intent(out) :: success
      logical :: successitems(5)
      integer :: gotvals(3),expvals(3)
      character(len=128) :: fnm, frmt
      integer :: u,i,k,l,nev,ndof,InpStat,itmp

      successitems(:)=.FALSE.
      write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
              '_',nr,'_eigv.rst'

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Read file containing the eigenvalue lists
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
         successitems(1)=(InpStat.eq.0)

         if (successitems(1)) then

            successitems(2)=.TRUE. !!! OLD read im
            read(u,*)

            if (successitems(2)) then

!              Loop over nodes
               DO i=1,SIZE(Ham%nt)

!                 Read the layer #, mode #, # eigenvalues and # sub-modes
                  read(u,*,IOSTAT=InpStat) itmp,nev,ndof
                  successitems(3)=(InpStat.eq.0)
                  if (.not.successitems(3)) exit

!                 Error checking
                  gotvals=(/itmp,nev,ndof/)
                  expvals=(/i,Ham%nt(i)%B,Ham%nt(i)%ndof()/)
                  successitems(4)=(&
                     (itmp.eq.i .and. ndof.eq.Ham%nt(i)%ndof()) .or. &
                     (itmp.eq.i .and. ndof.eq.0 .and. nev.eq.0))
                  if (.not.successitems(4)) exit

!                 Skip nodes that are not yet processed
                  IF (ndof.eq.0 .and. nev.eq.0) CYCLE

!                 Overwrite the block size with the restart file value
!                 which may have been determined from $truncate namelist
                  call Ham%nt(i)%setb(nev) 

!                 Make sure eigenvalue and assignment arrays are allocated
!                 The bottom layer should be already allocated
                  IF (.not.ALLOCATED(Ham%nt(i)%assgn)) &
                     ALLOCATE(Ham%nt(i)%assgn(nev,ndof))
                  IF (.not.ALLOCATED(Ham%nt(i)%eig)) &
                     ALLOCATE(Ham%nt(i)%eig(nev))
                  IF (.not.ALLOCATED(Ham%nt(i)%delta)) &
                     ALLOCATE(Ham%nt(i)%delta(nev))

!                 Read the eigenvalues/assignments/deltas for each node
                  DO k=1,nev
                     read(u,*,IOSTAT=InpStat) &
                     (Ham%nt(i)%assgn(k,l),l=1,ndof),&
                      Ham%nt(i)%eig(k),Ham%nt(i)%delta(k)
                     successitems(5)=(InpStat.eq.0)
                     if (.not.successitems(5)) exit
                  ENDDO
                  if (.not.ALL(successitems)) exit
                  Ham%nt(i)%done=.TRUE.
               ENDDO
            endif ! successitems(2)
         endif ! successitems(1)

         close(u)

      ENDIF rank0

      call bcast(gotvals,mpi_io_rank)
      call bcast(expvals,mpi_io_rank)
      call bcast(successitems,mpi_io_rank)
      success=ALL(successitems(:))

      IF (mpirank.eq.mpi_prnt_rank) THEN

         if (.not.successitems(1)) then
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         elseif (.not.successitems(2)) then
            write(*,*) 'Final node number could not be read'
         elseif (.not.successitems(3)) then
            write(*,*) 'node-nev-ndof designations', &
                       ' could not be read'
         elseif (.not.successitems(4)) then
            write(*,*) 'node :',gotvals(1),'; nev :',gotvals(2),&
                       '; nsubm :',gotvals(3),' read, but'
            write(*,*) 'node :',expvals(1),'; nev :',expvals(2),&
                       '; nsubm :',expvals(3),' expected'
         elseif (.not.successitems(5)) then
            write(*,*) 'Eigenvalues/assignments/delta could not be read'
         else
            write(*,'(5X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'
         endif

      ENDIF

!     Broadcast eigenvalues if read succeeded
      if (success) call BcastEigenvalues(Ham)

      end subroutine ReadEigenvalues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastEigenvalues(Ham)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts eigenvalues from MPI rank used to write to other ranks

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer :: i,nev,ndof

      DO i=1,SIZE(Ham%nt)
!        Broadcast the done state, which determines if anything else
!        needs to be broadcasted
         call bcast(Ham%nt(i)%done,mpi_io_rank)

         IF (.not.Ham%nt(i)%done) CYCLE

!        Broadcast the block size from the restart file
         call bcast(Ham%nt(i)%B,mpi_io_rank)

         nev=Ham%nt(i)%B
         ndof=Ham%nt(i)%ndof()

         IF (mpirank.ne.mpi_io_rank) THEN
            IF (.not.ALLOCATED(Ham%nt(i)%assgn)) &
               ALLOCATE(Ham%nt(i)%assgn(nev,ndof))
            IF (.not.ALLOCATED(Ham%nt(i)%eig)) &
               ALLOCATE(Ham%nt(i)%eig(nev))
            IF (.not.ALLOCATED(Ham%nt(i)%delta)) &
               ALLOCATE(Ham%nt(i)%delta(nev))
         ENDIF

         call bcast(Ham%nt(i)%assgn,mpi_io_rank)
         call bcast(Ham%nt(i)%eig,mpi_io_rank)
         call bcast(Ham%nt(i)%delta,mpi_io_rank)
      ENDDO

      end subroutine BcastEigenvalues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveEigenvalues(ML,Ham,nr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of eigenvalues/assignments to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in) :: nr
      character(len=128) :: fnm,frmt
      integer :: u,i,k,l,nm,nev,ndof

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Open output file
         write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
                 '_',nr,'_eigv.rst'
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!        Record the number of matrices
         write(u,'(A)') 'Node / eigenvalues / sub-modes'

!        Loop over nodes
         DO i=1,SIZE(Ham%nt)
!           Save the node #, # eigenvalues, # assignment dofs
            IF (Ham%nt(i)%done) THEN
               nev=SIZE(Ham%nt(i)%assgn,1)
               ndof=SIZE(Ham%nt(i)%assgn,2)
            ELSE
               nev=0
               ndof=0
            ENDIF
            write(u,'(4(I0,X))') i,nev,ndof

            IF (.not.Ham%nt(i)%done) CYCLE

!           Write the eigenvalues/assignments/deltas for each node
            write(frmt,'(A,I0,A)') '(',ndof,'(I0,X),2(E23.16,X))'
            DO k=1,nev
               write(u,frmt) &
               (Ham%nt(i)%assgn(k,l),l=1,ndof),&
                Ham%nt(i)%eig(k),Ham%nt(i)%delta(k)
            ENDDO
         ENDDO
         write(u,*)
         close(u)

      ENDIF rank0

      end subroutine SaveEigenvalues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadActivationCounter(im,ia,na,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads wavefunction from file during a restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      integer, intent(in)  :: im,na
      integer, intent(out) :: ia
      logical :: successitems(3)
      integer :: u,InpStat,iaread,naread
      character(len=128) :: fnm

      ia=1

!     Do nothing if restart file is 'none' or only 1 activation
      IF ((ML%resfile(1:4).seq.'none') .or. (na.eq.1)) RETURN

      successitems(:)=.FALSE.
      write(fnm,'(2A,I0,A)') &
      TRIM(ADJUSTL(ML%resfile)),'_act_node_',im,'.rst'

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Read file containing the eigenvalue list
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", &
              IOSTAT=InpStat)
         successitems(1)=(InpStat.eq.0)

         if (successitems(1)) then
            read(u,*)
            read(u,*,IOSTAT=InpStat) iaread,naread
            successitems(2)=(InpStat.eq.0)
            if (successitems(2)) then
               successitems(3)=(naread.eq.na .and. iaread.ge.1 &
                          .and. iaread.le.na)
               if (successitems(3)) ia=iaread
            endif
         endif

         CLOSE(u)

      ENDIF rank0

      call bcast(ia,mpi_io_rank)
      call bcast(successitems,mpi_io_rank)

      IF (mpirank.eq.mpi_prnt_rank) THEN
         if (.not. successitems(1)) then
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         elseif (.not. successitems(2)) then
            write(*,'(5X,2A,I0)') TRIM(ADJUSTL(fnm)),&
                                 ' could not be read for node ',im
         elseif (.not. successitems(3)) then
            write(*,'(5X,3(A,I0),A)') &
            'Bad activation counter data for node ',im,&
            '; starting activations at (',ia,'/',na,')'
         else
            write(*,'(5X,3A,I0,A,I0,A)') &
            TRIM(ADJUSTL(fnm)),' read successfully; ',&
            'restarting at H_coup activation (',ia,'/',na,')'
         endif
      ENDIF

      end subroutine ReadActivationCounter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveActivationCounter(im,ia,na,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads wavefunction from file during a restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      integer, intent(in) :: im,ia,na
      integer :: u
      character(len=128) :: fnm

!     Do nothing if restart file is 'none' or only 1 activation
      IF ((ML%resfile(1:4).seq.'none') .or. (na.eq.1)) RETURN

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Open output file
         write(fnm,'(2A,I0,A)') &
         TRIM(ADJUSTL(ML%resfile)),'_act_node_',im,'.rst'
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!        Record the activation state for this node
         write(u,'(A)') 'iact, nact:'
         write(u,'(I0,X,I0)') ia,na
         close(u)

      ENDIF rank0 

      end subroutine SaveActivationCounter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsi(isavi,bounds,eigv,delta,Q,ML,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads wavefunction from file during a restart

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      real*8, intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isavi
      logical, intent(out) :: success

      success=.FALSE.

      IF (.not.ML%dorestart) RETURN

!     Try to read the first psi file
      call ReadPsiFile(isavi,bounds,eigv,delta,Q,ML,1,success)

!     If that didn't work, try the second file
      IF (.not.success) THEN
         call ReadPsiFile(isavi,bounds,eigv,delta,Q,ML,2,success)
      ENDIF

!     If this job is a restart, after a read, successful or not, the
!     flag dorestart should be set to .FALSE. to prevent future reads
      ML%dorestart=.FALSE.

      end subroutine ReadPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsiFile(isvi,bounds,eigv,delta,Q,ML,nr,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads wavefunction file

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CP), ALLOCATABLE :: Qt(:)
      integer, intent(in)    :: nr
      real*8, intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isvi
      character(len=128)    :: fnm
      logical, intent(out) :: success
      logical :: successitems(12)
      integer, allocatable :: nbas(:)
      real*8, allocatable  :: eigt(:),deltt(:)
      real*8  :: boundt(2)
      integer :: gotvals(5)
      integer :: u,i,j,ndof,nev,nrk,isavti,InpStat

      successitems(:)=.FALSE.
      write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
              '_',nr,'_psi.rst'

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Read file containing the eigenvalue list
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", &
              FORM='UNFORMATTED',IOSTAT=InpStat)
         successitems(1)=(InpStat.eq.0)

         if (successitems(1)) then
!           Read the iteration, # eigenvalues
            read(u,IOSTAT=InpStat) isavti,nev
            successitems(2)=(InpStat.eq.0)
         endif

         if (successitems(2)) then     
!           Make sure # eigenvalues matches what is in the block
            successitems(3)=(nev.eq.SIZE(eigv))
            gotvals(1)=nev
         endif

         if (successitems(3)) then
!           Read spectral bounds
            read(u,IOSTAT=InpStat) boundt
            successitems(4)=(InpStat.eq.0)
         endif

         if (successitems(4)) then
!           Read the eigenvalues
            ALLOCATE(eigt(nev),deltt(nev))
            read(u,IOSTAT=InpStat) eigt
            successitems(5)=(InpStat.eq.0)
         endif

         if (successitems(5)) then
!           Read the deltas
            read(u,IOSTAT=InpStat) deltt
            successitems(6)=(InpStat.eq.0)
         endif

         if (successitems(6)) then
!           Read the wavefunction
            ALLOCATE(Qt(nev))
            DO i=1,nev
               gotvals(2)=i
               read(u,IOSTAT=InpStat) ndof,nrk
               successitems(7)=(InpStat.eq.0)
               if (.not.successitems(7)) exit

!              Make sure ndof of w.f. read matches that of Q
               successitems(8)=(ndof.eq.SIZE(Q(i)%nbas))
               gotvals(3)=ndof
               if (.not.successitems(8)) exit

               ALLOCATE(nbas(ndof))
               read(u,IOSTAT=InpStat) nbas
               successitems(9)=(InpStat.eq.0)
               if (.not.successitems(9)) exit

               successitems(10)=.TRUE.
!              Make sure nbas of w.f. to be read matches that of Q
               DO j=1,ndof
                  gotvals(4)=j
                  IF (nbas(j).ne.Q(i)%nbas(j)) THEN
                     successitems(10)=.FALSE.
                     gotvals(5)=nbas(j)
                     EXIT
                  ENDIF
               ENDDO
               if (.not.successitems(10)) exit

               Qt(i)=NewCP(nrk,nbas)
               DEALLOCATE(nbas)
               read(u,IOSTAT=InpStat) Qt(i)%coef
               successitems(11)=(InpStat.eq.0)
               if (.not.successitems(11)) exit

               read(u,IOSTAT=InpStat) Qt(i)%base
               successitems(12)=(InpStat.eq.0)
               if (.not.successitems(12)) exit

            ENDDO
         endif ! successitems(6)

         close(u)

         if (ALL(successitems(:))) then
            isvi=isavti
            eigv=eigt
            delta=deltt
            bounds=boundt
            DO i=1,nev
               call ReplaceVwithW(Q(i),Qt(i))
            ENDDO
         endif

      ENDIF rank0

      call bcast(gotvals,mpi_io_rank)
      call bcast(successitems,mpi_io_rank)
      success=ALL(successitems(:))

      IF (mpirank.eq.mpi_prnt_rank) THEN

         if (.not.successitems(1)) then
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         elseif (.not.successitems(2)) then
            write(*,*) 'Could not write iteration nr, nr eigenvalues'
         elseif (.not.successitems(3)) then
            write(*,*) 'Number of eigenvalues read:',gotvals(1)
            write(*,*) 'Number of eigenvalues expected:',SIZE(eigv)
         elseif (.not.successitems(4)) then
            write(*,*) 'Could not read spectral bounds'
         elseif (.not.successitems(5)) then
            write(*,*) 'Could not read eigenvalues'
         elseif (.not.successitems(6)) then
            write(*,*) 'Could not read deltas'
         elseif (.not.successitems(7)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read ndof, nrk'
         elseif (.not.successitems(8)) then
            write(*,'(X,A,I0,A,I0)') 'Psi(',gotvals(2),&
              ': number of degrees-of-freedom read:',gotvals(3)
            write(*,*) 'Number of degrees-of-freedom expected:',&
              SIZE(Q(i)%nbas)
         elseif (.not.successitems(9)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read nbas'
         elseif (.not.successitems(10)) then
            write(*,'(X,3(A,I0),A)') 'Psi(',gotvals(2),&
              '): nbas(',gotvals(4),') = ',gotvals(5),' read'
            write(*,'(X,3(A,I0),A)') 'Psi(',gotvals(2),&
              '): nbas(',gotvals(4),') = ',&
              Q(gotvals(2))%nbas(gotvals(4)),' expected'
         elseif (.not.successitems(11)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read coefficients'
         elseif (.not.successitems(12)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read base'
         else
            write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                              ' read successfully!'
         endif
         
         if (.not.success) &
            write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                              ' could not be read'

      ENDIF

      IF (ALLOCATED(Qt)) DEALLOCATE(Qt)
      IF (ALLOCATED(eigt)) DEALLOCATE(eigt)
      IF (ALLOCATED(deltt)) DEALLOCATE(deltt)

      if (success) call BcastPsi(isvi,bounds,eigv,delta,Q)

      end subroutine ReadPsiFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastPsi(isvi,bounds,eigv,delta,Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts eigenvalues from MPI rank 0 to other ranks

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      real*8, intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isvi
      integer :: i,nev

      call bcast(isvi,mpi_io_rank)
      call bcast(bounds,mpi_io_rank)
      call bcast(eigv,mpi_io_rank)
      call bcast(delta,mpi_io_rank)

      nev=SIZE(eigv)
      DO i=1,nev
         call BcastCP(Q(i),mpi_io_rank)
      ENDDO

      end subroutine BcastPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsi(isavi,bounds,eigv,delta,Q,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP), ALLOCATABLE, INTENT(IN) :: Q(:)
      real*8, intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi

!     If restart file is 'none', exit without saving
      IF (ML%resfile(1:4).seq.'none') RETURN

!     Save the psi file TWICE just in case job crashes during a write,
!     resulting in a corrupt psi file
      call SavePsiFile(isavi,bounds,eigv,delta,Q,ML,1)
      call SavePsiFile(isavi,bounds,eigv,delta,Q,ML,2)

      end subroutine SavePsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsiFile(isavi,bounds,eigv,delta,Q,ML,nr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP), ALLOCATABLE, INTENT(IN) :: Q(:)
      real*8, intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi,nr
      character(len=128) :: fnm,frmt
      integer :: u,i,j,k,l,nm,nev,nsubm

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Set parameters
         nev=SIZE(eigv)

!        Open output file
         write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
                 '_',nr,'_psi.rst'
         u = LookForFreeUnit()
         open(u, FILE=TRIM(ADJUSTL(fnm)),FORM='UNFORMATTED',&
              STATUS="UNKNOWN")

         write(u) isavi,nev
         write(u) bounds
         write(u) eigv
         write(u) delta
         DO i=1,nev
            write(u) SIZE(Q(i)%nbas),SIZE(Q(i)%coef)
            write(u) Q(i)%nbas
            write(u) Q(i)%coef
            write(u) Q(i)%base
         ENDDO

         close(u)

      ENDIF rank0

      end subroutine SavePsiFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsi_CP8(isavi,bounds,eigv,delta,Q,ML,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      TYPE (CP8), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      real(kind=8), intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isavi
      logical, intent(out) :: success

      success=.FALSE.

      IF (.not.ML%dorestart) RETURN

!     Try to read the first psi file
      call ReadPsiFile_CP8(isavi,bounds,eigv,delta,Q,ML,1,success)

!     If that didn't work, try the second file
      IF (.not.success) THEN
         call ReadPsiFile_CP8(isavi,bounds,eigv,delta,Q,ML,2,success)
      ENDIF

!     If this job is a restart, after a read, successful or not, the
!     flag dorestart should be set to .FALSE. to prevent future reads
      ML%dorestart=.FALSE.

      end subroutine ReadPsi_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsiFile_CP8(isvi,bounds,eigv,delta,Q,ML,nr,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP8), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CP8), ALLOCATABLE :: Qt(:)
      integer, intent(in)    :: nr
      real(kind=8), intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isvi
      character(len=128)    :: fnm
      logical, intent(out) :: success
      logical :: successitems(14)
      integer, allocatable :: rows(:),cols(:)
      real(kind=8), allocatable  :: eigt(:),deltt(:)
      real(kind=8) :: boundt(2)
      integer :: gotvals(7)
      integer :: u,i,j,ndof,nev,nrk,isavti,InpStat

      successitems(:)=.FALSE.
      write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
              '_',nr,'_psi.rst'

      rank0 : IF (mpirank.eq.0) THEN

!        Read file containing the eigenvalue list
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", &
              FORM='UNFORMATTED',IOSTAT=InpStat)
         successitems(1)=(InpStat.eq.0)

         if (successitems(1)) then
!           Read the iteration, # eigenvalues
            read(u,IOSTAT=InpStat) isavti,nev
            successitems(2)=(InpStat.eq.0)
         endif

         if (successitems(2)) then     
!           Make sure # eigenvalues matches what is in the block
            successitems(3)=(nev.eq.SIZE(eigv))
            gotvals(1)=nev
         endif

         if (successitems(3)) then
!           Read spectral bounds
            read(u,IOSTAT=InpStat) boundt
            successitems(4)=(InpStat.eq.0)
         endif

         if (successitems(4)) then
!           Read the eigenvalues
            ALLOCATE(eigt(nev),deltt(nev))
            read(u,IOSTAT=InpStat) eigt
            successitems(5)=(InpStat.eq.0)
         endif

         if (successitems(5)) then
!           Read the deltas
            read(u,IOSTAT=InpStat) deltt
            successitems(6)=(InpStat.eq.0)
         endif

         if (successitems(6)) then
!           Read the wavefunction
            ALLOCATE(Qt(nev))
            DO i=1,nev
               gotvals(2)=i
               read(u,IOSTAT=InpStat) ndof,nrk
               successitems(7)=(InpStat.eq.0)
               if (.not.successitems(7)) exit

!              Make sure ndof of w.f. read matches that of Q
               successitems(8)=(ndof.eq.SIZE(Q(i)%nbas))
               gotvals(3)=ndof
               if (.not.successitems(8)) exit

               ALLOCATE(rows(ndof))
               read(u,IOSTAT=InpStat) rows
               successitems(9)=(InpStat.eq.0)
               if (.not.successitems(9)) exit

               ALLOCATE(cols(ndof))
               read(u,IOSTAT=InpStat) cols
               successitems(10)=(InpStat.eq.0)
               if (.not.successitems(10)) exit

!              Make sure rows and cols of w.f. to be read matches that of Q
               successitems(11)=.TRUE.
               DO j=1,ndof
                  gotvals(4)=j
                  IF (rows(j).ne.Q(i)%rows(j)) THEN
                     successitems(11)=.FALSE.
                     gotvals(5)=rows(j)
                     EXIT
                  ENDIF
               ENDDO
               if (.not.successitems(11)) exit

               successitems(12)=.TRUE.
               DO j=1,ndof
                  gotvals(6)=j
                  IF (cols(j).ne.Q(i)%cols(j)) THEN
                     successitems(12)=.FALSE.
                     gotvals(7)=cols(j)
                     EXIT
                  ENDIF
               ENDDO
               if (.not.successitems(12)) exit

               call Qt(i)%new(nrk,rows,cols)
               DEALLOCATE(rows,cols)
               read(u,IOSTAT=InpStat) Qt(i)%coef
               successitems(13)=(InpStat.eq.0)
               if (.not.successitems(13)) exit

               read(u,IOSTAT=InpStat) Qt(i)%base
               successitems(14)=(InpStat.eq.0)
               if (.not.successitems(14)) exit

            ENDDO
         endif ! successitems(6)

         close(u)

         if (ALL(successitems(:))) then
            isvi=isavti
            eigv=eigt
            delta=deltt
            bounds=boundt
            DO i=1,nev
               call Q(i)%replace(Qt(i))
            ENDDO
         endif

      ENDIF rank0

      call bcast(gotvals,mpi_io_rank)
      call bcast(successitems,mpi_io_rank)
      success=ALL(successitems(:))

      IF (mpirank.eq.mpi_prnt_rank) THEN

         if (.not.successitems(1)) then
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         elseif (.not.successitems(2)) then
            write(*,*) 'Could not write iteration nr, nr eigenvalues'
         elseif (.not.successitems(3)) then
            write(*,*) 'Number of eigenvalues read:',gotvals(1)
            write(*,*) 'Number of eigenvalues expected:',SIZE(eigv)
         elseif (.not.successitems(4)) then
            write(*,*) 'Could not read spectral bounds'
         elseif (.not.successitems(5)) then
            write(*,*) 'Could not read eigenvalues'
         elseif (.not.successitems(6)) then
            write(*,*) 'Could not read deltas'
         elseif (.not.successitems(7)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read ndof, nrk'
         elseif (.not.successitems(8)) then
            write(*,'(X,A,I0,A,I0)') 'Psi(',gotvals(2),&
              ': number of degrees-of-freedom read:',gotvals(3)
            write(*,*) 'Number of degrees-of-freedom expected:',&
                       Q(i)%D()
         elseif (.not.successitems(9)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read rows'
         elseif (.not.successitems(10)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read cols'
         elseif (.not.successitems(11)) then
            write(*,'(X,3(A,I0,A))') 'Psi(',gotvals(2),&
              '): rows(',gotvals(4),') = ',gotvals(5),'read'
            write(*,'(X,3(A,I0,A))') 'Psi(',gotvals(2),&
              '): rows(',gotvals(4),') = ',&
              Q(gotvals(2))%M(gotvals(4)),'expected'
         elseif (.not.successitems(12)) then
            write(*,'(X,3(A,I0,A))') 'Psi(',gotvals(2),&
              '): cols(',gotvals(6),') = ',gotvals(7),'read'
            write(*,'(X,3(A,I0,A))') 'Psi(',gotvals(2),&
              '): cols(',gotvals(6),') = ',&
              Q(gotvals(2))%N(gotvals(6)),'expected'
         elseif (.not.successitems(13)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read coefficients'
         elseif (.not.successitems(14)) then
            write(*,'(X,A,I0,A)') 'Psi(',gotvals(2),&
              '): could not read base'
         else
            write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                              ' read successfully!'
         endif
         
         if (.not.success) &
            write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                              ' could not be read'

      ENDIF

      IF (ALLOCATED(Qt)) DEALLOCATE(Qt)
      IF (ALLOCATED(eigt)) DEALLOCATE(eigt)
      IF (ALLOCATED(deltt)) DEALLOCATE(deltt)

      if (success) call BcastPsi_CP8(isvi,bounds,eigv,delta,Q)

      end subroutine ReadPsiFile_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastPsi_CP8(isvi,bounds,eigv,delta,Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts eigenvalues from MPI rank 0 to other ranks

      implicit none
      TYPE (CP8), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      real(kind=8), intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isvi
      integer :: i,nev

      call bcast(isvi,mpi_io_rank)
      call bcast(bounds,mpi_io_rank)
      call bcast(eigv,mpi_io_rank)
      call bcast(delta,mpi_io_rank)

      nev=SIZE(eigv)
      DO i=1,nev
         call Bcast_CP8(Q(i),mpi_io_rank)
      ENDDO

      end subroutine BcastPsi_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsi_CP8(isavi,bounds,eigv,delta,Q,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP8), ALLOCATABLE, INTENT(IN) :: Q(:)
      real(kind=8), intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi

!     If restart file is 'none', exit without saving
      IF (ML%resfile(1:4).seq.'none') RETURN

!     Save the psi file TWICE just in case job crashes during a write,
!     resulting in a corrupt psi file
      call SavePsiFile_CP8(isavi,bounds,eigv,delta,Q,ML,1)
      call SavePsiFile_CP8(isavi,bounds,eigv,delta,Q,ML,2)

      end subroutine SavePsi_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsiFile_CP8(isavi,bounds,eigv,delta,Q,ML,nr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CP8), ALLOCATABLE, INTENT(IN) :: Q(:)
      real(kind=8), intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi,nr
      character(len=128) :: fnm,frmt
      integer :: u,i,j,k,l,nm,nev,nsubm

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Set parameters
         nev=SIZE(eigv)

!        Open output file
         write(fnm,'(2A,I0,A)') TRIM(ADJUSTL(ML%resfile)),&
                 '_',nr,'_psi.rst'
         u = LookForFreeUnit()
         open(u, FILE=TRIM(ADJUSTL(fnm)),FORM='UNFORMATTED',&
              STATUS="UNKNOWN")

         write(u) isavi,nev
         write(u) bounds
         write(u) eigv
         write(u) delta
         DO i=1,nev
            write(u) SIZE(Q(i)%nbas),SIZE(Q(i)%coef)
            write(u) Q(i)%rows
            write(u) Q(i)%cols
            write(u) Q(i)%coef
            write(u) Q(i)%base
         ENDDO

         close(u)

      ENDIF rank0

      end subroutine SavePsiFile_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE RESTART

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
