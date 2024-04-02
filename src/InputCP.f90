!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
       
      implicit none

      TYPE CPpar
           integer :: ncycle,npow,lowmem,truncation
           integer :: ncpu,psirank,hrank,psinals,hnals
           integer :: rs(33)
           real*8  :: solvtol
           logical :: update,dorestart,opt
           character(len=48) :: resfile
           character(len=5)  :: system
           character(len=4)  :: solver
           character(len=3)  :: red2D,redND
      END TYPE CPpar

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine StartInputCP(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main routine for reading and processing mode combination data

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64) :: inpfile

      IF (mpirank.eq.0) &
      write(*,'(/X,A/)') 'Reading input file (CP.inp)...'

      inpfile='CP.inp'
      CALL ReadMLCPInputs(cpp,inpfile)
      CALL PrintMLCPInputs(cpp)
      CALL BcastMLCPInputs(cpp)

      end subroutine StartInputCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadMLCPInputs(cpp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Reads input file for CP-format code for computing eigenvalues

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64), intent(in) :: fnm
      integer      :: i,u,InpStat

!     Read from MPI rank 0
      rank0 : IF (mpirank.eq.0) THEN

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
!     truncation layer options
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
!     random seed
      read(u,*)
      read(u,*) (cpp%rs(i),i=1,33)

      CLOSE(u)

      ENDIF rank0

      end subroutine ReadMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CP.inp

      implicit none
      TYPE (CPpar),INTENT(IN) :: cpp
      integer :: i

      rank0 : IF (mpirank.eq.0) THEN

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
      write(*,'(X,A,2X,I5)') 'Truncation criterion      (truncation):',&
                             cpp%truncation
      write(*,'(X,A,2X,L5)') 'Use vector updates            (update):',&
                             cpp%update
      write(*,'(X,A,2X,L5)') 'Optimize PES by coord. rotation  (opt):',&
                             cpp%opt
      write(*,'(X,A,2X,ES11.4)') &
                'Solver convergence criterion (solvtol):  ',cpp%solvtol
      write(*,'(X,A,2X,A)') 'Restart file name            (resfile):',&
                             cpp%resfile
      write(*,'(X,A,2X,33(I0,X))') &
                            'Random seed                       (rs):',&
                             (cpp%rs(i),i=1,33)
      write(*,'(/X,A)') '*********************************************'

      ENDIF rank0

      end subroutine PrintMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts MLtree to all MPI ranks

      implicit none
      TYPE (CPpar) :: cpp

!     Broadcast variables
      call bcast(cpp%ncycle)
      call bcast(cpp%npow)
      call bcast(cpp%lowmem)
      call bcast(cpp%truncation)
      call bcast(cpp%ncpu)
      call bcast(cpp%psirank)
      call bcast(cpp%hrank)
      call bcast(cpp%psinals)
      call bcast(cpp%hnals)
      call bcast(cpp%rs)
      call bcast(cpp%solvtol)
      call bcast(cpp%update)
      call bcast(cpp%dorestart)
      call bcast(cpp%opt)
      call bcast(cpp%resfile)
      call bcast(cpp%system)
      call bcast(cpp%solver)
      call bcast(cpp%red2D)
      call bcast(cpp%redND)

      end subroutine BcastMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
