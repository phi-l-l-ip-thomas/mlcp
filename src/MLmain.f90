!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE RANDOM
      USE MODECOMB
      USE SEPDREPN
      USE HAMILSETUP 
      USE INPUTCP
      USE RESTART
      USE REDUCTION
      USE ALSPOW
      USE LINSOLVER
      USE CPMMM
      USE BLOCKUTILS
      USE MODEH
      USE GUESS
      USE SOLVER
      USE UPDATER
      USE ANALYZER
      USE CHEBLIB
!!!
      USE TESTCPR8
!!!

      implicit none
      TYPE (CPpar)        :: cpp
      TYPE (MLtree)       :: ML
      TYPE (Hamiltonian)  :: Ham
      TYPE (CP), ALLOCATABLE :: Q(:)
      TYPE (CP) :: H,W
      real*8, allocatable  :: eigv(:),delta(:)
      integer, allocatable :: rs(:)
      integer :: d(3),t(3)
      integer :: il,im,j,trm,ilrst,imrst
      real*8  :: t1,t2
      character(len=64) :: frmt

      CALL prepare_mpi()

      call idate(d)
      call itime(t)
      call CPU_TIME(t1)

      IF (mpirank.eq.mpi_prnt_rank) then
         write(*,'(X,A/)') '##########################################'
         write(*,*)        '    Multi-layer CP-format TISE solver     '
         write(*,*)        '          by Phillip S. Thomas            '
         write(*,*)        '      based on the CP-format solver       '
         write(*,*)        '            of Arnaud Leclerc             '
         write(*,*)        '         Version ML2f 03-07-2024          '
         write(*,'(/X,A/)') '##########################################'
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank) THEN
         call PrintWallTime('MLCP initialized')
         write(*,'(X,A,I0,A)') 'running on ',mpinodes,' MPI processes'
      ENDIF

!     Set up the mode combination module, read input
      CALL StartModeComb(ML)     

!     Read input file, assign parameters
      CALL StartInputCP(cpp)

!     Initialize random number generator
      CALL CPU_TIME(t2)
      CALL InitRandom(t2-t1,d,t,cpp%rs,rs)

!!!   TEST
!      call maintestcpr8
!!!

!     Set up and sort operators into layers; solve bottom layer nodes
      CALL SetupHamiltonian(cpp%system,cpp%opt,Ham,ML)

!     Parallelization setup
      CALL omp_set_num_threads(cpp%ncpu)

!     Restart a previous run
      call RestartSetup(ilrst,imrst,cpp,Ham,ML)
      IF (ilrst.lt.1) call SaveEigenInfo(1,ML%nmode(1),cpp,Ham,ML)

!!!!!! --- MAIN RUN --- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (mpirank.eq.mpi_prnt_rank) &
         write(*,'(/X,A)') '***** MAIN RUN *****'

      DO il=2,ML%nlayr

         IF (il.lt.ilrst) CYCLE

         DO im=1,ML%nmode(il)

            IF (il.eq.ilrst .and. im.le.imrst) CYCLE

            IF (mpirank.eq.mpi_prnt_rank) &
               write(*,'(/,X,A,I0,A,I0,/)') 'LAYER-MODE: ',il,'-',im

!           Build the mode block (Q) and Hamiltonian matrix (H) here
            call BuildModeHamiltonian(il,im,H,Ham,ML,cpp)

!           Make the initial guess
            call GuessPsi(il,im,eigv,delta,Q,Ham,ML,cpp)
            W=GuessWeights(il,im,4000.0,Ham,ML)

!           Calculate the mode eigenfunctions with the solver of choice
!           If the mode on the current layer contains only one mode from
!           the previous layer, we already have the eigenvalues and 
!           eigenfunctions so no need to run the solver

            trm=GetModeHNr(il,im,Ham)  ! mode term
            IF (Ham%ndof(trm,il).eq.1 .and. Ham%nop(trm,il).eq.1) THEN
               IF (mpirank.eq.mpi_prnt_rank) &
                  write(*,'(3X,A)') '(Mode solved previously)'
            ELSE
               call SolveHPsi(eigv,delta,cpp,Q,H,W)

!              Print the wall time upon completion of the solver
               IF (mpirank.eq.mpi_prnt_rank) THEN
                  write(frmt,'(X,2(A,I0),A)') &
                  'Layer-Mode ',il,'-',im,': solver finished'
                  write(*,*)
               ENDIF
               IF (mpirank.eq.mpi_prnt_rank) call PrintWallTime(frmt)
            ENDIF

!           Analyze wavefunction and assign levels
            call AnalyzePsi(il,im,eigv,delta,Q,Ham,ML)

!           Transform operators to the mode eigenfunction basis
            call UpdateH(il,im,eigv,Q,Ham,ML,cpp)

!           Save eigenvalues and operator matrices for restart
            call SaveEigenInfo(il,im,cpp,Ham,ML)

            DEALLOCATE(eigv,Q)
            call FlushCP(H)
            call FlushCP(W)
         ENDDO
      ENDDO

      call sync_mpi

      IF (mpirank.eq.mpi_prnt_rank) write(*,*)

!     Free memory
      call DisposeModeComb(ML)
      call FlushHamiltonian(Ham)

!     Dispose modules and get CPU time
      call DisposeInitModule()
      call DisposeGuessModule()
      call DisposePRODHVModule()
      call DisposeMVV()
      call DisposeEigen()
      call DisposeReduction()
      call DisposeALSPow()
      call DisposeALSUtils()
      call DisposeLinSolver()
      call DisposeMunkres()
      call DisposeUpdateModule()
      call DisposeAnalModule()

      call CPU_TIME(t2)
      call idate(d)
      call itime(t)
      deallocate(rs)

      IF (cpp%ncpu.eq.1 .and. mpirank.eq.mpi_prnt_rank) &
          write(*,'(/X,A,11X,f20.3)') &
         'MLCP total CPU run time (s)',t2-t1

      IF (mpirank.eq.mpi_prnt_rank) then      
         write(*,*)
         call PrintWallTime('MLCP finished')
         write(*,*)
      ENDIF

      call finalize_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
