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

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A/)') 'Reading input file (CP.inp)...'

      inpfile='CP.inp'
      CALL ReadMLCPInputs(cpp,inpfile)
      CALL BcastMLCPInputs(cpp)
      CALL PrintMLCPInputs(cpp)

      end subroutine StartInputCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadMLCPInputs(cpp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Reads input file for CP-format code for computing eigenvalues

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64), intent(in) :: fnm
      integer      :: i,u,InpStat

!     Read from mpi_io_rank 
      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

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

      subroutine SaveMLCPInputFile(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('CP.inp') to another file for restart

      implicit none
      TYPE (CPpar), intent(in) :: cpp
      character(len=64) :: fnm
      integer :: u,j

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

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
!     truncation layer criterion
      write(u,'(A)') 'truncation'
      write(u,'(I16)') cpp%truncation
!     USE vector updates
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
!     random seed
      write(u,'(A)') 'random'
      write(u,'(33(I0,X))') (cpp%rs(j),j=1,33)
      close(u)

      ENDIF rank0

      end subroutine SaveMLCPInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CP.inp

      implicit none
      TYPE (CPpar),INTENT(IN) :: cpp
      integer :: i

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

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
      call bcast(cpp%ncycle,mpi_io_rank)
      call bcast(cpp%npow,mpi_io_rank)
      call bcast(cpp%lowmem,mpi_io_rank)
      call bcast(cpp%truncation,mpi_io_rank)
      call bcast(cpp%ncpu,mpi_io_rank)
      call bcast(cpp%psirank,mpi_io_rank)
      call bcast(cpp%hrank,mpi_io_rank)
      call bcast(cpp%psinals,mpi_io_rank)
      call bcast(cpp%hnals,mpi_io_rank)
      call bcast(cpp%rs,mpi_io_rank)
      call bcast(cpp%solvtol,mpi_io_rank)
      call bcast(cpp%update,mpi_io_rank)
      call bcast(cpp%dorestart,mpi_io_rank)
      call bcast(cpp%opt,mpi_io_rank)
      call bcast(cpp%resfile,mpi_io_rank)
      call bcast(cpp%system,mpi_io_rank)
      call bcast(cpp%solver,mpi_io_rank)
      call bcast(cpp%red2D,mpi_io_rank)
      call bcast(cpp%redND,mpi_io_rank)

      end subroutine BcastMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
