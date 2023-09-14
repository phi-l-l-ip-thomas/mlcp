!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE RESTARTCSPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Restart a crashed CSPES calculation

      use ERRORTRAP
      use UTILS
      use MODEGRID
      use INPUTFIELDS
      use SAMPLEPES
      use CHEBLIB
      use SEPDREPN
      use CPCONFIG

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CSRestartSetup(csp,GP,CGrid,X,B,BB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master routine for reading restart files for the CSPES program

      implicit none
      TYPE (CSpar), INTENT(INOUT) :: csp
      TYPE (GridPar), INTENT(IN)  :: GP
      TYPE (ChebObj), INTENT(IN)  :: CGrid(:)
      TYPE (Configs), INTENT(OUT) :: X,B
      TYPE (CPvec), INTENT(OUT)   :: BB
      character(len=64) :: fnm
      logical :: success,found

!     Look for _CS.rst and _grids.rst. If both are present
!     then set dorestart=.TRUE.

      csp%dorestart=.FALSE.
!     Look for _CS.rst
      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_CS.rst'
      INQUIRE(FILE=TRIM(ADJUSTL(fnm)), EXIST=found)
      IF (found) THEN
         write(*,'(/X,2A)') TRIM(ADJUSTL(fnm)),' exists'
!        Look for _grids.rst
         write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_grids.rst'
         INQUIRE(FILE=TRIM(ADJUSTL(fnm)), EXIST=found)
         IF (found) THEN
            write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' exists'
!           If this is a 'read' run, no need to look for a _pes.rst file
!           since the PES data is already stored in the readfile
            IF (csp%sys(1:4).seq.'read') THEN
               csp%dorestart=.TRUE.
            ELSE
!              Look for _pes.rst
               write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_pes.rst'
               INQUIRE(FILE=TRIM(ADJUSTL(fnm)), EXIST=found)
               IF (found) THEN
                  write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' exists'
                  csp%dorestart=.TRUE.
               ELSE
                  write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' is missing'
               ENDIF
            ENDIF
         ELSE
            write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' is missing'
         ENDIF
      ENDIF

!     Check restart data
      IF (csp%dorestart) THEN
         write(*,'(/X,A/)') &
               'This job is a restart. Validating restart data...'
         call ValidateCSRestart(csp,GP)

!        Read the _pes.rst file if this is not a 'read' run
         IF (csp%sys(1:4).seq.'read') THEN
!           If this is a read run, the readfile is read later
         ELSE
            write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_pes.rst'
            call ReadPES(CGrid,B,BB,fnm,csp%nsamples,csp%usegrid)
            write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'
         ENDIF

!        Read the file containing PES coefficients
         write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_1_coef.rst'
         call ReadPESConfigFile(CGrid,X,fnm,csp%pesrank,success)
         IF (.not.success) THEN
            write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_2_coef.rst'
            call ReadPESConfigFile(CGrid,X,fnm,csp%pesrank,success)
            IF (.not.success) call AbortWithError(&
               "PES coefficient restart file could not be read")
         ENDIF
         write(*,'(X,2A)') TRIM(ADJUSTL(fnm)),' read successfully!'
      ENDIF

!     Save 'CS.inp' and 'grids.inp' unless restart file is 'none'
      IF (csp%resfile(1:4).seq.'none') RETURN

      call SaveCSPESInputFile(csp)
      call SaveGridDat(GP,csp%resfile)

      end subroutine CSRestartSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateCSRestart(csp,GP)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares restart data with job parameters to determine if calculation
! can be successfully restarted

      implicit none
      TYPE (CSpar), INTENT(IN) :: csp
      TYPE (CSpar)  :: csrst
      TYPE (GridPar), INTENT(IN) :: GP
      TYPE (GridPar) :: GPrst
      character(len=64) :: fnm
      integer :: i

!     Read the restart input files
      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_CS.rst'
      call ReadCSPESInputFile(csrst,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_grids.rst'
      call ReadGridDat(GPrst,fnm)

!     Validate input against CS.rst
      IF (csp%sys.ne.csrst%sys) THEN
         write(*,*) 'Old system: ',csrst%sys,&
                    '; New system: ',csp%sys
         call AbortWithError(&
         'ValidateCSRestart(): System change not allowed in a restart!')
      ENDIF
      IF (csp%usegrid.neqv.csrst%usegrid) THEN
         write(*,*) 'Old choice of usegrid: ',csrst%usegrid,&
                    '; New choice of usegrid: ',csp%usegrid
         call AbortWithError(&
         'ValidateCSRestart(): Sample rule must be same in a restart!')
      ENDIF
      IF (csp%nmcfg.neqv.csrst%nmcfg) THEN
         write(*,*) 'Old choice of nmcfg: ',csrst%nmcfg,&
                    '; New choice of nmcfg: ',csp%nmcfg
         call AbortWithError(&
         'ValidateCSRestart(): Guess method must be same in a restart!')
      ENDIF
      IF (csp%maxnmode.ne.csrst%maxnmode) THEN
         write(*,*) 'Old maxnmode: ',csrst%maxnmode,&
                    '; New maxnmode: ',csp%maxnmode
         call AbortWithError(&
         'ValidateCSRestart(): maxnmode must be same in a restart!')
      ENDIF

!     Print changes to job parameters
      IF (csp%ncpu.ne.csrst%ncpu) write(*,'(2(A,I0))') &
         ' * ncpu changed from ',csrst%ncpu,' to ',csp%ncpu
      IF (csp%nsamples.ne.csrst%nsamples) write(*,'(2(A,I0))') &
         ' * nsamples changed from ',csrst%nsamples,' to ',csp%nsamples
      IF (csp%pesrank.ne.csrst%pesrank) write(*,'(2(A,I0))') &
         ' * pesrank changed from ',csrst%pesrank,' to ',csp%pesrank
      IF (csp%ncycle.ne.csrst%ncycle) write(*,'(2(A,I0))') &
         ' * ncycle changed from ',csrst%ncycle,' to ',csp%ncycle
      IF (csp%sigma.ne.csrst%sigma) write(*,'(2(A,ES10.3))') &
         ' * sigma changed from ',csrst%sigma,' to ',csp%sigma
      IF (csp%ctol.ne.csrst%ctol) write(*,'(2(A,ES10.3))') &
         ' * ctol changed from ',csrst%ctol,' to ',csp%ctol
      IF (csp%gtol.ne.csrst%gtol) write(*,'(2(A,ES10.3))') &
         ' * gtol changed from ',csrst%gtol,' to ',csp%gtol
      IF (csp%dtol.ne.csrst%dtol) write(*,'(2(A,ES10.3))') &
         ' * dtol changed from ',csrst%dtol,' to ',csp%dtol
      IF (csp%sys(1:4).seq.'read') THEN
         IF (csp%readfile.ne.csrst%readfile) write(*,*) &
            ' * readfile changed from ',TRIM(ADJUSTL(csrst%readfile)),&
            ' to ',TRIM(ADJUSTL(csp%readfile))
      ENDIF

!     Validate against grids.rst. All grid parameters must be the same

!     First compare ndof by comparing array sizes. Only one array
!     comparison must be made since ReadGridDat() ensures that array
!     sizes are consistent within a given GridPar structure
      IF (SIZE(GP%npts).ne.SIZE(GPrst%npts)) &
         call AbortWithError('ValidateCSRestart(): inconsistent # DOF')

      DO i=1,SIZE(GP%npts)
        IF (GP%npts(i).ne.GPrst%npts(i)) call AbortWithError(&
            'ValidateCSRestart(): $basis do not match')
        IF (GP%bc(i).ne.GPrst%bc(i)) call AbortWithError(&
            'ValidateCSRestart(): $boundary do not match')
        IF (GP%egrid(i).neqv.GPrst%egrid(i)) call AbortWithError(&
            'ValidateCSRestart(): $gridtype do not match')
        IF (GP%bounds(i,1).ne.GPrst%bounds(i,1) .or. &
            GP%bounds(i,2).ne.GPrst%bounds(i,2)) call AbortWithError(&
            'ValidateCSRestart(): $parameters do not match')
      ENDDO

      call FlushGridDat(GPrst)

      end subroutine ValidateCSRestart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveCSPESInputFile(csp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('CS.inp') to another file for restart

      implicit none
      TYPE (CSpar) :: csp
      character(len=64) :: fnm
      integer :: u

      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_CS.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     System = PES to fit
      write(u,'(A)') 'System'
      write(u,'(A16)') csp%sys
!     Sample from predetermined grid
      write(u,'(A)') 'usegrid'
      write(u,'(L16)') csp%usegrid
!     Configuration guess rule
      write(u,'(A)') 'nmcfg'
      write(u,'(L16)') csp%nmcfg
!     number of processors
      write(u,'(A)') 'NCPU'
      write(u,'(I16)') csp%ncpu
!     number of PES sample points
      write(u,'(A)') 'NSAMPLES'
      write(u,'(I16)') csp%nsamples
!     reduction rank for PES fitting
      write(u,'(A)') 'PESRANK'
      write(u,'(I16)') csp%pesrank
!     max number of fitting cycles
      write(u,'(A)') 'NCYCLE'
      write(u,'(I16)') csp%ncycle
!     maximum n-mode coupling
      write(u,'(A)') 'MAXNMODE'
      write(u,'(I16)') csp%maxnmode
!     basis pursuit convergence criterion
      write(u,'(A)') 'sigma'
      write(u,'(ES16.1)') csp%sigma
!     coefficient tolerance
      write(u,'(A)') 'ctol'
      write(u,'(ES16.1)') csp%ctol
!     gradient tolerance
      write(u,'(A)') 'gtol'
      write(u,'(ES16.1)') csp%gtol
!     fractional residual decrease tolerance
      write(u,'(A)') 'dtol'
      write(u,'(ES16.1)') csp%dtol
!     PES data file
      write(u,'(A)') 'readfile'
      write(u,'(A48)') csp%readfile
!     restart file name
      write(u,'(A)') 'resfile'
      write(u,'(A48)') csp%resfile

      close(u)

      end subroutine SaveCSPESInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveGridDat(GP,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('grids.inp') to a restart file

      implicit none
      TYPE (GridPar), INTENT(IN) :: GP
      character(len=48), intent(in) :: fnm
      character(len=64) :: fname,frmt
      integer :: u,i,j,ndof
      integer, allocatable :: tmp(:)

      write(fname,'(2A)') TRIM(ADJUSTL(fnm)),'_grids.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fname)), STATUS="UNKNOWN")

      ndof=SIZE(GP%npts)

!     Write basis section
      write(u,*)
      write(u,'(A)') '$basis'
      write(frmt,'(A,I0,A)') '(',ndof,'(I0,X),A)'
      write(u,frmt) (GP%npts(i),i=1,ndof),'/'
      write(u,'(A)') '$end-basis'

      ALLOCATE(tmp(ndof))

!     Convert character array 'bc' to integer array for write
      DO i=1,ndof
         IF (GP%bc(i).seq.'fin') THEN
            tmp(i)=0
         ELSEIF (GP%bc(i).seq.'sem') THEN
            tmp(i)=1
         ELSEIF (GP%bc(i).seq.'inf') THEN
            tmp(i)=2
         ELSEIF (GP%bc(i).seq.'per') THEN
            tmp(i)=3
         ELSEIF (GP%bc(i).seq.'cos') THEN
            tmp(i)=4
         ELSE
            call AbortWithError(&
                 'SaveGridDat(): unrecognized boundary condition')
         ENDIF
      ENDDO

!     Write boundary section
      write(u,*)
      write(u,'(A)') '$boundary'
      write(frmt,'(A,I0,A)') '(',ndof,'(I0,X),A)'
      write(u,frmt) (tmp(i),i=1,ndof),'/'
      write(u,'(A)') '$end-boundary'

!     Convert logical array 'egrid' to integer array for write
      tmp(:)=0
      DO i=1,ndof
         IF (GP%egrid(i)) tmp(i)=1
      ENDDO

!     Write gridtype section
      write(u,*)
      write(u,'(A)') '$gridtype'
      write(frmt,'(A,I0,A)') '(',ndof,'(I0,X),A)'
      write(u,frmt) (tmp(i),i=1,ndof),'/'
      write(u,'(A)') '$end-gridtype'

      DEALLOCATE(tmp)

!     Write parameters section
      write(u,*)
      write(u,'(A)') '$parameters'
      write(frmt,'(A,I0,A)') '(',ndof,'(ES22.15,X),A)'
      DO j=1,2
         write(u,frmt) (GP%bounds(i,j),i=1,ndof),'/'
      ENDDO
      write(u,'(A)') '$end-parameters'

      close(u)

      end subroutine SaveGridDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePESsamples(CGrid,csp,B,BB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Saves PES samples used in fitting to file for restart

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (CSpar), INTENT(IN)   :: csp
      TYPE (CPvec), INTENT(IN)   :: BB
      TYPE (Configs), INTENT(IN) :: B
      character(len=64) :: fnm

!     Exit without saving if the PES is in a readfile
      IF (csp%sys(1:4).seq.'read') RETURN
!     Exit without saving if no restart file is specified
      IF (csp%resfile(1:4).seq.'none') RETURN

      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_pes.rst'
      call WritePES(CGrid,B,BB,fnm,csp%usegrid)

      end subroutine SavePESsamples

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SavePESCoef(csp,X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Saves PES configurations and coefficients from optimization (so far)
! to file for restart

      implicit none
      TYPE (CSpar), INTENT(IN)   :: csp
      TYPE (Configs), INTENT(IN) :: X
      character(len=64) :: fnm

!     Exit without saving if no restart file is specified
      IF (csp%resfile(1:4).seq.'none') RETURN

      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_1_coef.rst'
      call WritePESConfigFile(X,fnm)
      write(fnm,'(2A)') TRIM(ADJUSTL(csp%resfile)),'_2_coef.rst'
      call WritePESConfigFile(X,fnm)

      end subroutine SavePESCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE RESTARTCSPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
