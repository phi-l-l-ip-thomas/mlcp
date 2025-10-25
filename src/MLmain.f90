!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MYACC
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
      USE TIMINGS
!!!
      USE TESTCPR8
      USE CPMM8
      USE CPVV8
!!!

      implicit none
      TYPE (CPpar)        :: cpp
      TYPE (MLtree)       :: ML
      TYPE (Hamiltonian)  :: Ham
      TYPE (CP), ALLOCATABLE :: Q(:)
      TYPE (CP) :: H,W
      real(kind=8), allocatable :: eigv(:),delta(:)
      real(kind=8) :: bounds(2)
      integer, allocatable :: rs(:)
      integer :: d(3),t(3)
      integer :: im,j,imrst,ndof,nnode
      real(kind=8)  :: t1,t2
      character(len=64) :: frmt

      call prepare_mpi()
      call init_device()

      call idate(d)
      call itime(t)
      call CPU_TIME(t1)

      IF (mpirank.eq.mpi_prnt_rank) then
         write(*,'(X,A/)') &
                    '##################################################'
         write(*,*) '                       _                          '
         write(*,*) '                      | |                         '
         write(*,*) '             _ __ ___ | | ___  _ __               '
         write(*,*) "            | '_ ` _ `| |/ __`| '_ `              "
         write(*,*) '            | | | | | | | |__ | |_) |             '
         write(*,*) '            |_| |_| |_|_|`___/| .__/              '
         write(*,*) '                              | |                 '
         write(*,*) '                              |_|                 '
         write(*,*) ' ------------------------------------------------ '
         write(*,*) '    Multi-Layer CP-format solver for the Time-    '
         write(*,*) '   Independent Vibrational Schrodinger Equation   '
         write(*,*) ' ------------------------------------------------ '
         write(*,*) ' If using this code results in publication, please'
         write(*,*) ' cite as:'
         write(*,*)
         write(*,*) ' "Using Nested Contractions and a Hierarchical'
         write(*,*) ' Tensor Format to Compute Vibrational Spectra of'
         write(*,*) ' Molecules with Seven Atoms,” P. S. Thomas and T.'
         write(*,*) ' Carrington. J. Phys. Chem. A, 2015, 119(52),'
         write(*,*) ' 13074-13091. [DOI:10.1021/acs.jpca.5b10015]'
         write(*,*) '  (https://doi.org/10.1021/acs.jpca.5b10015)'
         write(*,'(/X,A/)') &
                    '##################################################'
      ENDIF

      IF (mpirank.eq.mpi_prnt_rank) THEN
         call PrintWallTime('MLCP initialized')
      ENDIF

!     Read input file, assign parameters
      CALL StartInputCP(cpp)

!     Read the mode combination data
      CALL StartModeComb(ML,cpp%verbosity)

!     Initialize random number generator
      CALL CPU_TIME(t2)
      CALL InitRandom(t2-t1,d,t,cpp%rs,rs)
!!!   TEST
!      call maintestcpr8
!      call testmpicycle()
!!!

!     Set up and sort operators into layers; solve bottom layer nodes
      CALL SetupHamiltonian(cpp%system,cpp%pe_transform,cpp%verbosity,&
                                       cpp%pe_trans_fac,Ham,ML)

      ndof=SIZE(Ham%ops,1)
      nnode=SIZE(Ham%nt)

!     Restart a previous run
      call RestartSetup(imrst,cpp,Ham,ML)
!      IF (imrst.lt.ndof) call SaveEigenInfo(ndof,cpp,Ham)

!!!!!! --- MAIN RUN --- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (mpirank.eq.mpi_prnt_rank) &
         write(*,'(/X,A)') '***** MAIN RUN *****'

      DO im=1,nnode !ndof+1,nnode
!      DO im=ndof+1,nnode

         IF (im.eq.nnode .and. (.not.cpp%dotopnode)) THEN
            IF (mpirank.eq.mpi_prnt_rank) write(*,'(/X,A)') &
               'Calculation *paused* at top layer node!'
            EXIT
         ENDIF

         IF (im.le.imrst) CYCLE

         IF (mpirank.eq.mpi_prnt_rank) &
            write(*,'(/,X,A,I0,A,I0,A,I0,A,/)') &
            '--- NODE ',im,' (LAYER-MODE: ',&
             Ham%nt(im)%mlil,'-',Ham%nt(im)%mlim,') ---'

!        Build the mode block (Q) and Hamiltonian matrix (H) here
         call BuildModeHamiltonian(im,H,Ham,cpp)

!        Make the initial guess
         call GuessPsi(im,eigv,delta,bounds,Q,Ham,ML,cpp)
         W=GuessWeights(im,4000.0,Ham)

!        Calculate the node eigenfunctions with the solver of choice
!        If the node on the current layer contains only one sub-node from
!        the previous layer, we already have the eigenvalues and 
!        eigenfunctions so no need to run the solver
         IF (Ham%nt(im)%nHterm().eq.0 .and. &
            Ham%nt(im)%nsubm().eq.1) THEN
            IF (mpirank.eq.mpi_prnt_rank) &
               write(*,'(3X,A)') '(Mode solved previously)'
         ELSE
            call SolveHPsi(eigv,delta,bounds,cpp,Q,H,W)
         ENDIF

!        Analyze wavefunction and assign levels
         call AnalyzePsi(im,eigv,delta,Q,Ham)

!        Transform operators to the mode eigenfunction basis
         call UpdateH(im,Q,Ham,cpp)

!        Save eigenvalues and operator matrices for restart
         call SaveEigenInfo(im,cpp,Ham)

!        Print the wall time
         IF (mpirank.eq.mpi_prnt_rank) THEN
            write(*,*)
            write(frmt,'(A,I0,A)') &
            'node ',im,': finished'
            call PrintWallTime(frmt)
         ENDIF

         DEALLOCATE(eigv,Q)
         call FlushCP(H)
         call FlushCP(W)
      ENDDO

      call sync_mpi

      IF (mpirank.eq.mpi_prnt_rank) write(*,*)

!     Free memory and print timings
      call Flush_ModeComb(ML)
      call Flush_Hamiltonian(Ham)

!      call DisposeReduction() ! to be removed

      call Show_section_times()
      call Show_module_times()

      call CPU_TIME(t2)
      call idate(d)
      call itime(t)
      deallocate(rs)

      IF (nomp_threads.eq.1 .and. mpirank.eq.mpi_prnt_rank) &
          write(*,'(/X,A,11X,f20.3)') &
         'MLCP total CPU run time (s)',t2-t1

      IF (mpirank.eq.mpi_prnt_rank) then      
         write(*,*)
         call PrintWallTime('MLCP finished')
         write(*,*)
      ENDIF

      call flush_device()
      call finalize_mpi()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
