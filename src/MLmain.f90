!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MODECOMB
      USE SEPDREPN
      USE HAMILSETUP 
      USE INPUTFIELDS
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
      character(len=64) :: inpfile,frmt

      call idate(d)
      call itime(t)
      call CPU_TIME(t1)
      rs=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
          mod(INT(t1),7),d(1),d(2),d(3),t(1),t(2),t(3),mod(INT(t1),5)/)
      call random_seed(PUT=rs)
!      write(*,*) 'Random seed = ',rs

      write(*,'(X,A/)') '############################################'
      write(*,*)        '     Multi-layer CP-format TISE solver      '
      write(*,*)        '           by Phillip S. Thomas             '
      write(*,*)        '       based on the CP-format solver        '
      write(*,*)        '             of Arnaud Leclerc              '
      write(*,*)        '          Version ML2f 04-27-2020           '
      write(*,'(/X,A/)') '############################################'

      call PrintWallTime('MLCP initialized')

!      call Test_calcPk()

!     Set up the mode combination module, read input
      write(*,'(/X,A/)') 'Setting up mode-combination module...'
      inpfile='layers.inp'
      CALL Init_ModeComb
      CALL ReadModeDat(ML,inpfile)
      CALL ValidateModeDat(ML)
      CALL PrintModeDat(ML)

!     Read input file, assign parameters
      write(*,'(/X,A/)') 'Reading input file (CP.inp)...' 
      inpfile='CP.inp'
      CALL ReadMLCPInputFile(cpp,inpfile)
      CALL PrintMLCPInputs(cpp)

      write(*,'(/X,A/)') 'Hamiltonian setup...'

      CALL SetupHamiltonian(cpp%system,cpp%opt,Ham,ML)

!     Parallelization setup
      CALL omp_set_num_threads(cpp%ncpu)

!     Print the eigenvalues for the bottom layer
      DO im=1,ML%nmode(1)
         write(*,'(/,X,A,I0,A,I0,/)') 'LAYER-MODE: ',1,'-',im
         write(*,*) 'Eigenvalues   : ',0,&
         (Ham%eig(1,im)%evals(j),j=1,SIZE(Ham%eig(1,im)%evals))
         write(*,*)
         DO j=1,SIZE(Ham%eig(1,im)%evals)
            write(*,'(I4,A,X,I2,X,f19.12)') j,')',j-1,&
                 Ham%eig(1,im)%evals(j)-Ham%eig(1,im)%evals(1)
         ENDDO
      ENDDO

      call RestartSetup(ilrst,imrst,cpp,Ham,ML)
      IF (ilrst.lt.1) call SaveEigenInfo(1,ML%nmode(1),cpp,Ham,ML)

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
            call GuessPsi(il,im,eigv,Q,Ham,ML)
            W=GuessWeights(il,im,4000.0,Ham,ML)

!           Calculate the mode eigenfunctions with the solver of choice
!           If the mode on the current layer contains only one mode from
!           the previous layer, we already have the eigenvalues and 
!           eigenfunctions so no need to run the solver

            trm=GetModeHNr(il,im,Ham)  ! mode term
            IF (Ham%ndof(trm,il).eq.1 .and. Ham%nop(trm,il).eq.1) THEN
               write(*,'(3X,A)') '(Mode solved previously)'
            ELSE
               call SolverAlg(eigv,delta,cpp,Q,H,W,il,ML%nlayr)

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
      call Flush_ModeComb(ML)
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

      write(*,*)
      call PrintWallTime('MLCP finished')
      write(*,*)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM MULTILAYERCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
