!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
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
      USE TESTCPR


      implicit none
      TYPE (CPpar)        :: cpp
      TYPE (MLtree)       :: ML
      TYPE (Hamiltonian)  :: Ham
      TYPE (CP), ALLOCATABLE :: Q(:)
      TYPE (CP) :: H,W
      real*8, allocatable :: eigv(:),delta(:)
      integer :: rs(33),d(3),t(3)
      integer :: il,im,j,trm,ilrst,imrst
      real*8  :: t1,t2
      character(len=64) :: frmt

      CALL prepare_mpi()

      call idate(d)
      call itime(t)
      call CPU_TIME(t1)
      rs=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
          mod(INT(t1),7),d(1),d(2),d(3),t(1),t(2),t(3),mod(INT(t1),5)/)
      call random_seed(PUT=rs)

      IF (mpirank.eq.0) then
      write(*,'(X,A/)') '############################################'
      write(*,*)        '     Multi-layer CP-format TISE solver      '
      write(*,*)        '           by Phillip S. Thomas             '
      write(*,*)        '       based on the CP-format solver        '
      write(*,*)        '             of Arnaud Leclerc              '
      write(*,*)        '          Version ML2f 03-07-2024           '
      write(*,'(/X,A/)') '############################################'
      ENDIF

      call sync_mpi()
      write(*,'(A,I0,A,I0,A)') 'running MLCP from rank (',&
                                mpirank,'/',mpinodes,')...'
      call sync_mpi()

      IF (mpirank.eq.0) call PrintWallTime('MLCP initialized')

!     Set up the mode combination module, read input
      CALL StartModeComb(ML)     

!     Read input file, assign parameters
      CALL StartInputCP(cpp)

!     Set up and sort operators into layers; solve bottom layer nodes
      CALL SetupHamiltonian(cpp%system,cpp%opt,Ham,ML)

!     Parallelization setup
      CALL omp_set_num_threads(cpp%ncpu)

!     Restart a previous run
      call RestartSetup(ilrst,imrst,cpp,Ham,ML)
      IF (ilrst.lt.1) call SaveEigenInfo(1,ML%nmode(1),cpp,Ham,ML)

      IF (mpirank.eq.0) THEN !!! TEST-RK0

      IF (ANY(cpp%rs.ne.0)) THEN
         rs(1:33)=cpp%rs(1:33)
         write(*,'(/X,A/)') 'Random seed used from input file...'     
      ELSE
         rs=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         mod(INT(t1),7),d(1),d(2),d(3),t(1),t(2),t(3),mod(INT(t1),5)/)
         write(*,'(/X,A,33(X,I0)/)') 'Random seed generated: ',&
                                    (rs(j),j=1,33)
      ENDIF
      call random_seed(PUT=rs)

!!!!!! --- MAIN RUN --- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,'(/X,A)') '***** MAIN RUN *****'

      DO il=2,ML%nlayr

         IF (il.lt.ilrst) CYCLE

         DO im=1,ML%nmode(il)

            IF (il.eq.ilrst .and. im.le.imrst) CYCLE

            write(*,'(/,X,A,I0,A,I0,/)') 'LAYER-MODE: ',il,'-',im

!           Build the mode block (Q) and Hamiltonian matrix (H) here
            call BuildModeHamiltonian(il,im,H,Ham,ML,cpp)

!           Make the initial guess
            call GuessPsi(il,im,eigv,Q,Ham,ML,cpp)
            W=GuessWeights(il,im,4000.0,Ham,ML)

!           Calculate the mode eigenfunctions with the solver of choice
!           If the mode on the current layer contains only one mode from
!           the previous layer, we already have the eigenvalues and 
!           eigenfunctions so no need to run the solver

            trm=GetModeHNr(il,im,Ham)  ! mode term
            IF (Ham%ndof(trm,il).eq.1 .and. Ham%nop(trm,il).eq.1) THEN
               write(*,'(3X,A)') '(Mode solved previously)'
            ELSE
               call SolveHPsi(eigv,delta,cpp,Q,H,W)

!              Print the wall time upon completion of the solver
               write(frmt,'(X,2(A,I0),A)') &
               'Layer-Mode ',il,'-',im,': solver finished'
               write(*,*)
               call PrintWallTime(frmt)
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
      write(*,*)

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

      IF (cpp%ncpu.eq.1) write(*,'(/X,A,11X,f20.3)') &
         'MLCP total CPU run time (s)',t2-t1

      endif !!! TEST-RK0

      IF (mpirank.eq.0) then      
         write(*,*)
         call PrintWallTime('MLCP finished')
        write(*,*)
      ENDIF

      call finalize_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
