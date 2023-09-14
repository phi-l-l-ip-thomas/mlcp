!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE INPUTFIELDS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use ERRORTRAP
      use UTILS

      implicit none
      TYPE CPpar
           integer :: ncycle,npow,lowmem,truncation
           integer :: ncpu,psirank,hrank,psinals,hnals
           real*8  :: solvtol
           logical :: update,dorestart,opt
           character(len=48) :: resfile
           character(len=5)  :: system
           character(len=4)  :: solver
           character(len=3)  :: red2D,redND
      END TYPE CPpar

      TYPE CSpar
           integer :: nsamples,ncycle
           integer :: ncpu,pesrank,maxnmode
           real*8  :: sigma,ctol,gtol,dtol
           logical :: usegrid,dorestart,nmcfg
           character(len=64) :: readfile
           character(len=64) :: resfile
           character(len=5)  :: sys
      END TYPE CSpar

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadMLCPInputFile(cpp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Reads input file for CP-format code for computing eigenvalues

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64), intent(in) :: fnm
      integer      :: u,InpStat

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         call AbortWithError("Oh, no! Error reading input file")
      ENDIF

!     list of parameters

!     System = Hamiltonian to use
      read(u,*)
      read(u,*) cpp%system
!     NCPU = number of processors      
      read(u,*)
      read(u,*) cpp%ncpu
!     reduction type, 2-D modes
      read(u,*)
      read(u,*) cpp%red2D
!     reduction type, >2-D modes
      read(u,*)
      read(u,*) cpp%redND
!     reduction rank for wavefunction
      read(u,*)
      read(u,*) cpp%psirank
!     reduction rank for Hamiltonian
      read(u,*)
      read(u,*) cpp%hrank
!     number of ALS iterations for wavefunction
      read(u,*)
      read(u,*) cpp%psinals
!     number of ALS iterations for Hamiltonian
      read(u,*)
      read(u,*) cpp%hnals
!     Solver type
      read(u,*)
      read(u,*) cpp%solver
!     number of power/Cheb iteration cycles
      read(u,*)
      read(u,*) cpp%ncycle
!     number of power iterations per cycle
      read(u,*)
      read(u,*) cpp%npow
!     low memory calculation type
      read(u,*)
      read(u,*) cpp%lowmem
!     basis truncation option
      read(u,*)
      read(u,*) cpp%truncation
!     do vector updates
      read(u,*)
      read(u,*) cpp%update
!     PES optimization by coordinate rotation
      read(u,*)
      read(u,*) cpp%opt
!     solver tolerance (relative rms error of all states)
      read(u,*)
      read(u,*) cpp%solvtol
!     restart file name
      read(u,*)
      read(u,*) cpp%resfile

      CLOSE(u)

      end subroutine ReadMLCPInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CP.inp

      implicit none
      TYPE (CPpar),INTENT(IN) :: cpp

      write(*,'(X,A/)') '********** Input parameters read: ***********'
      write(*,'(X,A,2X,A5)') 'The Hamiltonian will be set up for    :',&
                             cpp%system
      write(*,'(X,A,2X,I5)') 'Number of processors            (ncpu):',&
                             cpp%ncpu
      write(*,'(X,A,4X,A3)') 'Reduction type for 2-D modes          :',&
                             cpp%red2D
      write(*,'(X,A,4X,A3)') 'Reduction type for >2-D modes         :',&
                             cpp%redND
      write(*,'(X,A,2X,I5)') 'Wavefunction reduced rank    (psirank):',&
                             cpp%psirank
      write(*,'(X,A,2X,I5)') 'Hamiltonian reduced rank       (hrank):',&
                             cpp%hrank
      write(*,'(X,A,2X,I5)') 'Number of ALS iterations-w.f.(psinals):',&
                             cpp%psinals
      write(*,'(X,A,2X,I5)') 'Number of ALS iterations-H     (hnals):',&
                             cpp%hnals
      write(*,'(X,A,2X,A5)') 'Eigensolver algorithm to use  (solver):',&
                             cpp%solver
      write(*,'(X,A,2X,I5)') 'Number of solver cycles       (ncycle):',&
                             cpp%ncycle
      write(*,'(X,A,2X,I5)') 'Number of Power iteratons       (npow):',&
                             cpp%npow
      write(*,'(X,A,2X,I5)') 'Low-memory calculation type   (lowmem):',&
                             cpp%lowmem
      write(*,'(X,A,2X,I5)') 'Basis truncation option   (truncation):',&
                             cpp%truncation
      write(*,'(X,A,2X,L5)') 'Use vector updates            (update):',&
                             cpp%update
      write(*,'(X,A,2X,L5)') 'Optimize PES by coord. rotation  (opt):',&
                             cpp%opt
      write(*,'(X,A,2X,ES11.4)') &
                'Solver convergence criterion (solvtol):  ',cpp%solvtol
      write(*,'(X,A,2X,A)') 'Restart file name            (resfile):',&
                             cpp%resfile
      write(*,'(/X,A)') '*********************************************'

      end subroutine PrintMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadCSPESInputFile(csp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Reads input file for CP-format code for computing eigenvalues

      implicit none
      TYPE (CSpar) :: csp
      integer      :: u,InpStat
      character(len=64), intent(in) :: fnm

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         call AbortWithError("Oh, no! Error reading input file")
      ENDIF

!     list of parameters
!     PES to fit
      read(u,*)
      read(u,*) csp%sys
!     Sample from predetermined grid
      read(u,*)
      read(u,*) csp%usegrid
!     Generate configs by increasing nmode
      read(u,*)
      read(u,*) csp%nmcfg
!     number of processors
      read(u,*)
      read(u,*) csp%ncpu
!     number of PES sample points
      read(u,*)
      read(u,*) csp%nsamples
!     reduction rank for PES fitting
      read(u,*)
      read(u,*) csp%pesrank
!     max number of fitting cycles (per tau)
      read(u,*)
      read(u,*) csp%ncycle
!     maximum n-mode coupling
      read(u,*)
      read(u,*) csp%maxnmode
!     basis pursuit convergence criterion
      read(u,*)
      read(u,*) csp%sigma
!     coefficient tolerance
      read(u,*)
      read(u,*) csp%ctol
!     gradient tolerance
      read(u,*)
      read(u,*) csp%gtol
!     fractional residual decrease tolerance
      read(u,*)
      read(u,*) csp%dtol
!     PES data file
      read(u,*)
      read(u,*) csp%readfile
!     restart file name
      read(u,*)
      read(u,*) csp%resfile

      CLOSE(u)

      end subroutine ReadCSPESInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCSPESInputs(csp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CS.inp

      implicit none
      TYPE (CSpar),INTENT(IN) :: csp

      write(*,'(/X,A/)') '********** Input parameters read: ***********'
      write(*,'(X,A,2X,A5)') 'The PES will be fit for               :',&
                             csp%sys
      write(*,'(X,A,2X,L5)') 'Random sample from grid      (usegrid):',&
                             csp%usegrid
      write(*,'(X,A,2X,L5)') 'Guess configs by incr. nmode   (nmcfg):',&
                             csp%nmcfg
      write(*,'(X,A,2X,I5)') 'Number of processors            (ncpu):',&
                             csp%ncpu
      write(*,'(X,A,2X,I5)') 'Number of PES samples       (nsamples):',&
                             csp%nsamples
      write(*,'(X,A,2X,I5)') 'PES reduced rank             (pesrank):',&
                             csp%pesrank
      write(*,'(X,A,2X,I5)') 'Number of PES fitting cycles  (ncycle):',&
                             csp%ncycle
      write(*,'(X,A,2X,I5)') 'Maximum coupled DOFs        (maxnmode):',&
                             csp%maxnmode
      write(*,'(X,A,2X,ES11.4)') &
                'Basis pursuit criterion        (sigma):  ',csp%sigma
      write(*,'(X,A,2X,ES11.4)') &
                'Coefficient tolerance           (ctol):  ',csp%ctol
      write(*,'(X,A,2X,ES11.4)') &
                'Gradient tolerance              (gtol):  ',csp%gtol
      write(*,'(X,A,2X,ES11.4)') &
                'Fractional resid. decrease tol. (dtol):  ',csp%dtol
      write(*,'(X,A,2X,A)') 'PES data file name          (readfile):',&
                             csp%readfile
      write(*,'(X,A,2X,A)') 'Restart file name            (resfile):',&
                             csp%resfile
      write(*,'(/X,A)') '*********************************************'

      end subroutine PrintCSPESInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INPUTFIELDS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
