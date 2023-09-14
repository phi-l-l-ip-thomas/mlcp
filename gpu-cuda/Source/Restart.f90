!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE RESTART

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Restart a crashed calculation (provided it doesn't crash too hard!)

      use ERRORTRAP
      use UTILS
      use MODECOMB
      use INPUTFIELDS
      use HAMILSETUP
      use SEPDREPN

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RestartSetup(il,im,cpp,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master routine for reading restart files

      implicit none
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (MLtree), INTENT(IN)   :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: Ham
      integer, intent(out) :: il,im
      character(len=64) :: fnm
      logical :: success,found

!     Look for _CP.rst and _layers.rst. If both are present
!     then set dorestart=.TRUE.

      cpp%dorestart=.FALSE.
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_CP.rst'
      INQUIRE(FILE=TRIM(ADJUSTL(fnm)), EXIST=found)
      IF (found) THEN
         write(*,'(/X,2A)') TRIM(ADJUSTL(fnm)),' exists'
         write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_layers.rst'
         INQUIRE(FILE=TRIM(ADJUSTL(fnm)), EXIST=found)
         IF (found) THEN
            write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' exists'
            cpp%dorestart=.TRUE.
         ELSE
            write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' is missing'
         ENDIF
      ENDIF 

!     Check restart data
      IF (cpp%dorestart) THEN
         write(*,'(X,A)') &
               'This job is a restart. Validating restart data...'
         call ValidateRestart(cpp,ML)

!        Read the file containing eigenvalues from finished nodes
         write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_eigv.rst'
         call ReadEigenvalues(il,im,Ham,ML,fnm,success)
         IF (.not.success) THEN
            write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_eigv.rst'
            call ReadEigenvalues(il,im,Ham,ML,fnm,success)
            IF (.not.success) call AbortWithError(&
               "Eigenvalue restart file could not be read")
         ENDIF
         write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'

!        Read the file containing operator matrices
         write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_oper.rst'
         call ReadOperMats(Ham,fnm,success)
         IF (.not.success) THEN
            write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_oper.rst'
            call ReadOperMats(Ham,fnm,success)
            IF (.not.success) call AbortWithError(&
               "Operator matrix restart file could not be read")
         ENDIF
         write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'
      ELSE
         il=0
         im=0
      ENDIF

!     Save 'CP.inp' and 'layers.inp' unless restart file is 'none'
      IF (cpp%resfile(1:4).seq.'none') RETURN

      call SaveMLCPInputFile(cpp)
      call SaveModeDat(ML,cpp%resfile)

      end subroutine RestartSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateRestart(cpp,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares restart data with job parameters to determine if calculation
! can be successfully restarted

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp
      TYPE (CPpar)  :: cprst
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (MLtree) :: MLrst
      character(len=64) :: fnm
      integer :: il,im

!     Read the restart input files
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_CP.rst'
      call ReadMLCPInputFile(cprst,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_layers.rst'
      call ReadModeDat(MLrst,fnm)

!     Validate input against CP.rst
      IF (cpp%system.ne.cprst%system) THEN
         write(*,*) 'Old system: ',cprst%system,&
                    '; New system: ',cpp%system
         call AbortWithError(&
         'ValidateRestart(): System change not allowed in a restart!')
      ENDIF
      IF (cpp%truncation.ne.cprst%truncation) THEN
         write(*,*) 'Old choice of truncation: ',cprst%truncation,&
                    '; New choice of truncation: ',cpp%truncation
         call AbortWithError(&
         'ValidateRestart(): truncation must not change in a restart!')
      ENDIF
      IF (cpp%opt.neqv.cprst%opt) THEN
         write(*,*) 'Old choice of opt: ',cprst%opt,&
                    '; New choice of opt: ',cpp%opt
         call AbortWithError(&
         'ValidateRestart(): Coord. change not allowed in a restart!')
      ENDIF

!     Print changes to job parameters
      IF (cpp%ncpu.ne.cprst%ncpu) write(*,*) &
         ' * ncpu changed from ',cprst%ncpu,' to ',cpp%ncpu
      IF (cpp%red2D.ne.cprst%red2D) write(*,*) &
         ' * red2D changed from ',cprst%red2D,' to ',cpp%red2D
      IF (cpp%redND.ne.cprst%redND) write(*,*) &
         ' * redND changed from ',cprst%redND,' to ',cpp%redND
      IF (cpp%psirank.ne.cprst%psirank) write(*,*) &
         ' * psirank changed from ',cprst%psirank,' to ',cpp%psirank
      IF (cpp%hrank.ne.cprst%hrank) write(*,*) &
         ' * hrank changed from ',cprst%hrank,' to ',cpp%hrank
      IF (cpp%psinals.ne.cprst%psinals) write(*,*) &
         ' * psinals changed from ',cprst%psinals,' to ',cpp%psinals
      IF (cpp%hnals.ne.cprst%hnals) write(*,*) &
         ' * hnals changed from ',cprst%hnals,' to ',cpp%hnals
      IF (cpp%solver.ne.cprst%solver) write(*,*) &
         ' * solver changed from ',cprst%solver,' to ',cpp%solver
      IF (cpp%ncycle.ne.cprst%ncycle) write(*,*) &
         ' * ncycle changed from ',cprst%ncycle,' to ',cpp%ncycle
      IF (cpp%npow.ne.cprst%npow) write(*,*) &
         ' * npow changed from ',cprst%npow,' to ',cpp%npow
      IF (cpp%lowmem.ne.cprst%lowmem) write(*,*) &
         ' * lowmem changed from ',cprst%lowmem,' to ',cpp%lowmem
      IF (cpp%update.neqv.cprst%update) write(*,*) &
         ' * update changed from ',cprst%update,' to ',cpp%update
      IF (cpp%solvtol.ne.cprst%solvtol) write(*,*) &
         ' * solvtol changed from ',cprst%solvtol,' to ',cpp%solvtol
      write(*,*)

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

      call Flush_ModeComb(MLrst)

      end subroutine ValidateRestart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveMLCPInputFile(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('CP.inp') to another file for restart

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64) :: fnm
      integer :: u

      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_CP.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     System = Hamiltonian to use
      write(u,'(A)') 'System'
      write(u,'(A16)') cpp%system
!     NCPU = number of processors      
      write(u,'(A)') 'NCPU'
      write(u,'(I16)') cpp%ncpu
!     reduction type, 2-D modes
      write(u,'(A)') 'red2D'
      write(u,'(A16)') cpp%red2D
!     reduction type, >2-D modes
      write(u,'(A)') 'redND'
      write(u,'(A16)') cpp%redND
!     reduction rank for wavefunction
      write(u,'(A)') 'psirank'
      write(u,'(I16)') cpp%psirank
!     reduction rank for Hamiltonian
      write(u,'(A)') 'hrank'
      write(u,'(I16)') cpp%hrank
!     number of ALS iterations for wavefunction
      write(u,'(A)') 'psinals'
      write(u,'(I16)') cpp%psinals
!     number of ALS iterations for Hamiltonian
      write(u,'(A)') 'hnals'
      write(u,'(I16)') cpp%hnals
!     Eigensolver algorithm
      write(u,'(A)') 'solver'
      write(u,'(A16)') cpp%solver
!     number of power/Cheb iteration cycles
      write(u,'(A)') 'ncycle'
      write(u,'(I16)') cpp%ncycle
!     number of power iterations per cycle
      write(u,'(A)') 'npow'
      write(u,'(I16)') cpp%npow
!     low memory calculation type
      write(u,'(A)') 'lowmem'
      write(u,'(I16)') cpp%lowmem
!     basis truncation option
      write(u,'(A)') 'truncation'
      write(u,'(I16)') cpp%truncation
!     use vector updates
      write(u,'(A)') 'update'
      write(u,'(L16)') cpp%update
!     PES optimization by coordinate rotation
      write(u,'(A)') 'optimize PES'
      write(u,'(L16)') cpp%opt
!     solver tolerance (relative rms error of all states)
      write(u,'(A)') 'solvtol'
      write(u,'(ES16.1)') cpp%solvtol
!     restart file name
      write(u,'(A)') 'resfile'
      write(u,'(A48)') cpp%resfile

      close(u)

      end subroutine SaveMLCPInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveModeDat(ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('layers.inp') with mode combination data to
! file (i.e. a restart file)

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      character(len=48), intent(in) :: fnm
      character(len=64) :: fname,frmt
      integer :: u,il,im

      write(fname,'(2A)') TRIM(ADJUSTL(fnm)),'_layers.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fname)), STATUS="UNKNOWN")

!     Write resort section
      write(u,*)
      write(u,'(A)') '$resort'
      write(frmt,'(A,I0,A)') '(',ML%ndof,'(I0,X),A)'
      write(u,frmt) (ML%resort(im),im=1,ML%ndof),'/'
      write(u,'(A)') '$end-resort'

!     Write basis section
      write(u,*)
      write(u,'(A)') '$basis'
      DO il=1,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X),A)'
         write(u,frmt) (ML%gdim(il,im),im=1,ML%nmode(il)),'/'
      ENDDO
      write(u,'(A)') '$end-basis'

!     Write layers section
      write(u,*)
      write(u,'(A)') '$layers'
      DO il=2,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X),A)'
         write(u,frmt) (ML%modcomb(il,im),im=1,ML%nmode(il)),'/'
      ENDDO
      write(u,'(A)') '$end-layers'
      write(u,*)

      close(u)

      end subroutine SaveModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveEigenInfo(il,im,cpp,Ham,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Saves eigenvalues, assignments, and operator matrices

      implicit none
      TYPE (CPpar), INTENT(IN)  :: cpp
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: Ham
      integer, intent(in) :: il,im
      character(len=64) :: fnm

!     If restart file is 'none', exit without saving
!     Also, no need to save after the last layer is finished
      IF (il.eq.ML%nlayr .or. (cpp%resfile(1:4).seq.'none')) RETURN

!     Save data twice to prevent a potential corrupted file write
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_eigv.rst'
      call SaveEigenvalues(il,im,Ham,ML,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_eigv.rst'
      call SaveEigenvalues(il,im,Ham,ML,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_oper.rst'
      call SaveOperMats(Ham,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_oper.rst'
      call SaveOperMats(Ham,fnm)

      end subroutine SaveEigenInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadOperMats(H,fnm,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of operator matrices to file for restart

      implicit none
      TYPE (Hamiltonian), INTENT(INOUT) :: H
      character(len=64), intent(in)     :: fnm
      logical, intent(out) :: success
      integer :: u,i,k,nmat,n,InpStat

      success=.TRUE.

!     Open file containing the operator matrices
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         success=.FALSE.
         RETURN
      ENDIF

!     Read the number of matrices
      read(u,*)
      read(u,*,IOSTAT=InpStat) nmat
      IF (InpStat /= 0) success=.FALSE.
      read(u,*)

!     Make sure that the number of operator matrices in the restart
!     file matches the number allocated in the current run

      IF (nmat.ne.SIZE(H%pops)) THEN
         write(*,*) 'Number of operators, current run :',SIZE(H%pops)
         write(*,*) 'Number of operators, restart file:',nmat
         success=.FALSE.
      ENDIF

!     Cycle through the primitive operator matrices
      IF (success) THEN
         DO i=1,nmat
            read(u,*,IOSTAT=InpStat) H%pops(i)%dof,n
            IF (InpStat /= 0) THEN
               success=.FALSE.
               EXIT
            ENDIF

!           The operator matrices are allocated and written when the 1D
!           (bottom layer) problems are solved, so these must be 
!           deallocated and reallocated in a restarted run
            IF (allocated(H%pops(i)%mat)) DEALLOCATE(H%pops(i)%mat)
            ALLOCATE(H%pops(i)%mat(n))

!           Read the matrix for each operator
            read(u,*,IOSTAT=InpStat) (H%pops(i)%mat(k),k=1,n)
            IF (InpStat /= 0) THEN
               success=.FALSE.
               EXIT
            ENDIF
            IF (.not.success) EXIT
         ENDDO
      ENDIF

      close(u)

      end subroutine ReadOperMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveOperMats(H,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of operator matrices to file for restart

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      character(len=64), intent(in)  :: fnm
      character(len=64) :: frmt
      integer :: u,i,k,n

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     Record the number of matrices
      write(u,'(A)') 'Number of operators:'
      write(u,'(I0)') SIZE(H%pops)
      write(u,'(A)') 'DOF# / n'

!     Write each operator matrix to file
      DO i=1,SIZE(H%pops)
         n=SIZE(H%pops(i)%mat)
         write(u,'(4(I0,X))') H%pops(i)%dof,n
         write(frmt,'(A,I0,A)') '(',n,'(E23.16,X))'
         write(u,frmt) (H%pops(i)%mat(k),k=1,n)
      ENDDO
      write(u,*)
      close(u)

      end subroutine SaveOperMats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadEigenvalues(il,im,H,ML,fnm,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads list of eigenvalues/assignments from restart file

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(INOUT) :: H
      integer, intent(out) :: il,im
      logical, intent(out) :: success
      character(len=64), intent(in) :: fnm
      character(len=64) :: frmt
      integer :: u,i,j,k,l,nm,nev,nsubm,InpStat,itmp,jtmp

      success=.TRUE.

!     Open file containing the eigenvalue lists
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         success=.FALSE.
         RETURN
      ENDIF

!     Read the number of matrices
      read(u,*)
      read(u,*) il,im
      IF (InpStat /= 0) success=.FALSE.
      read(u,*)

!     Loop over layers and modes
      DO i=1,il
         nm=ML%nmode(i)
         IF (i.eq.il) nm=im
         DO j=1,nm
!           Read the layer #, mode #, # eigenvalues and # sub-modes
            read(u,*) itmp,jtmp,nev,nsubm
            IF (InpStat /= 0) THEN
               success=.FALSE.
               EXIT
            ENDIF

!           Error checking
            IF (itmp.ne.i .or. jtmp.ne.j .or. & 
               nsubm.ne.ML%modcomb(i,j)) THEN
               write(*,*) 'layer :',itmp,'; mode :',jtmp,&
                          '; nsubm :',nsubm,' read'
               write(*,*) 'layer :',i,'; mode :',j,&
                          '; nsubm :',ML%modcomb(i,j),' expected'
               success=.FALSE.
               EXIT
            ENDIF
            IF (nev.ne.ML%gdim(i,j)) THEN
               write(*,*) 'layer :',itmp,'; mode :',jtmp
               write(*,*) ' nev = ',nev,' read'
               write(*,*) ' nev = ',ML%gdim(i,j),' expected'
               success=.FALSE.
               EXIT
            ENDIF

!           Make sure eigenvalue and assignment arrays are allocated
!           The bottom layer should be already allocated
            IF (.not.ALLOCATED(H%eig(i,j)%assgn)) &
               ALLOCATE(H%eig(i,j)%assgn(nev,nsubm))
            IF (.not.ALLOCATED(H%eig(i,j)%evals)) &
               ALLOCATE(H%eig(i,j)%evals(nev))

!           Read the eigenvalues/assignments for each layer/mode
            DO k=1,nev
               read(u,*) &
               (H%eig(i,j)%assgn(k,l),l=1,nsubm),H%eig(i,j)%evals(k)
               IF (InpStat /= 0) THEN
                  success=.FALSE.
                  EXIT
               ENDIF
            ENDDO
            IF (.not.success) EXIT
         ENDDO
         IF (.not.success) EXIT
      ENDDO

      close(u)

      end subroutine ReadEigenvalues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveEigenvalues(il,im,H,ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes list of eigenvalues/assignments to file for restart

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (Hamiltonian), INTENT(IN) :: H
      integer, intent(in) :: il,im
      character(len=64), intent(in)  :: fnm
      character(len=64) :: frmt
      integer :: u,i,j,k,l,nm,nev,nsubm

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     Record the number of matrices
      write(u,'(A)') 'Last layer/mode to be solved:'
      write(u,'(2(I0,X))') il,im
      write(u,'(A)') 'Layer / Mode / eigenvalues / sub-modes'

!     Loop over layers and modes
      DO i=1,il
         nm=ML%nmode(i)
         IF (i.eq.il) nm=im
         DO j=1,nm
!           Save the layer #, mode #, # eigenvalues, # sub-modes
            nev=SIZE(H%eig(i,j)%assgn,1)
            nsubm=SIZE(H%eig(i,j)%assgn,2)
            write(u,'(4(I0,X))') i,j,nev,nsubm

!           Write the eigenvalues/assignments for each layer/mode
            write(frmt,'(A,I0,A)') '(',nsubm,'(I0,X),E23.16,X)'
            DO k=1,nev
               write(u,frmt) &
               (H%eig(i,j)%assgn(k,l),l=1,nsubm),H%eig(i,j)%evals(k)
            ENDDO

         ENDDO
      ENDDO
      write(u,*)
      close(u)

      end subroutine SaveEigenvalues

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsi(isavi,bounds,eigv,delta,Q,cpp,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (CPpar), INTENT(INOUT) :: cpp
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      real*8, intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isavi
      logical, intent(out) :: success
      character(len=64) :: fnm

      success=.FALSE.

      IF (.not.cpp%dorestart) RETURN

!     Try to read the first psi file
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_psi.rst'
      call ReadPsiFile(isavi,bounds,eigv,delta,Q,fnm,success)

!     If that didn't work, try the second file
      IF (.not.success) THEN
         write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_psi.rst'
         call ReadPsiFile(isavi,bounds,eigv,delta,Q,fnm,success)
      ENDIF

!     If this job is a restart, after a read, successful or not, the
!     flag dorestart should be set to .FALSE. to prevent future reads
      cpp%dorestart=.FALSE.

      end subroutine ReadPsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPsiFile(isvi,bounds,eigv,delta,Q,fnm,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: Q(:)
      TYPE (CPvec), ALLOCATABLE :: Qt(:)
      real*8, intent(inout)  :: bounds(2),eigv(:),delta(:)
      integer, intent(inout) :: isvi
      character(len=64), intent(in) :: fnm
      logical, intent(out) :: success
      integer, allocatable :: nbas(:)
      real*8, allocatable  :: eigt(:),deltt(:)
      real*8  :: boundt(2)
      integer :: u,i,j,ndof,nev,nrk,isavti,InpStat

      success=.TRUE.

!     Open file containing the eigenvalue list
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", &
           FORM='UNFORMATTED',IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,'(/X,A/)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                            ' not found'
         success=.FALSE.
         RETURN
      ENDIF

!     Try to read the iteration, # eigenvalues, and spectral bounds
      read(u,IOSTAT=InpStat) isavti,nev
      IF (Inpstat /=0) success=.FALSE.

!     Make sure # eigenvalues matches what is in the block
      IF (nev.ne.SIZE(eigv)) success=.FALSE.

      read(u,IOSTAT=InpStat) boundt
      IF (Inpstat /=0) success=.FALSE.

!     If previous reads succeeded, read the eigenvalues and deltas
      IF (success) THEN
         ALLOCATE(eigt(nev),deltt(nev))
         read(u,IOSTAT=InpStat) eigt
         IF (Inpstat /=0) success=.FALSE.
         read(u,IOSTAT=InpStat) deltt
         IF (Inpstat /=0) success=.FALSE.
      ENDIF

!     If eigenvalues and deltas read successfully, read the wavefunction
      IF (success) THEN
         ALLOCATE(Qt(nev))
         DO i=1,nev
            read(u,IOSTAT=InpStat) ndof,nrk
            IF (Inpstat /=0) THEN
               success=.FALSE.
               EXIT
            ENDIF

!           Make sure ndof of w.f. to be read matches that of Q
            IF (ndof.ne.SIZE(Q(i)%nbas)) THEN
               success=.FALSE.
               EXIT
            ENDIF

            ALLOCATE(nbas(ndof))
            read(u,IOSTAT=InpStat) nbas
            IF (Inpstat /=0) THEN
               success=.FALSE.
               EXIT
            ENDIF

!           Make sure nbas of w.f. to be read matches that of Q
            DO j=1,ndof
               IF (nbas(j).ne.Q(i)%nbas(j)) success=.FALSE.
            ENDDO
            IF (.not.success) EXIT

            call NewCPvec(Qt(i),nbas,nrk)
            DEALLOCATE(nbas)
            read(u,IOSTAT=InpStat) Qt(i)%coef
            IF (Inpstat /=0) THEN
               success=.FALSE.
               EXIT
            ENDIF
            read(u,IOSTAT=InpStat) Qt(i)%base
            IF (Inpstat /=0) THEN
               success=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      close(u)

!     If all was read successfully, replace existing data with read data
      IF (success) THEN
         isvi=isavti
         eigv=eigt
         delta=deltt
         bounds=boundt
         DO i=1,nev
            call ReplaceVwithW(Q(i),Qt(i))
         ENDDO
         write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                            ' read successfully!'
      ELSE
         write(*,'(X,3A)') 'Psi restart file ',TRIM(ADJUSTL(fnm)),&
                            ' could not be read'
      ENDIF

      IF (ALLOCATED(Qt)) DEALLOCATE(Qt)
      IF (ALLOCATED(eigt)) DEALLOCATE(eigt)
      IF (ALLOCATED(deltt)) DEALLOCATE(deltt)

      end subroutine ReadPsiFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsi(isavi,bounds,eigv,delta,Q,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp
      TYPE (CPvec), ALLOCATABLE, INTENT(IN) :: Q(:)
      real*8, intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi
      character(len=64)   :: fnm

!     If restart file is 'none', exit without saving
      IF (cpp%resfile(1:4).seq.'none') RETURN

!     Save the psi file TWICE just in case job crashes during a write,
!     resulting in a corrupt psi file
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_1_psi.rst'
      call SavePsiFile(isavi,bounds,eigv,delta,Q,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_2_psi.rst'
      call SavePsiFile(isavi,bounds,eigv,delta,Q,fnm)

      end subroutine SavePsi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePsiFile(isavi,bounds,eigv,delta,Q,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes wavefunction to file for restart

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(IN) :: Q(:)
      real*8, intent(in)  :: bounds(2),eigv(:),delta(:)
      integer, intent(in) :: isavi
      character(len=64), intent(in) :: fnm
      character(len=64) :: frmt
      integer :: u,i,j,k,l,nm,nev,nsubm

!     Set parameters
      nev=SIZE(eigv)

!     Open output file
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

      end subroutine SavePsiFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE RESTART

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
