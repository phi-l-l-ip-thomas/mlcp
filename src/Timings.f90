!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE TIMINGS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets timings from parts of the code

      USE ERRORTRAP
      USE UTILS
      USE MYMPI

      USE MODECOMB
      USE HAMILSETUP
      USE MODEH
      USE GUESS
      USE SOLVER
      USE SOLVER8
      USE ANALYZER
      USE UPDATER

      USE LINALG
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE ALSPOW
      USE ALSUTILS
      USE LINSOLVER
      USE CPMM8
      USE CPVV8
      USE CPLS8
      USE LINALG8
      USE MUNKRES

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Show_section_times()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints timings of modules

      implicit none

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(/X,A,29X,3(7X,A))') &
         'Section:','Min time','Max time','Ave time'
      ENDIF

      call Dispose_ModeComb_Module()
      call Dispose_HamilSetup_Module()
      call Dispose_ModeH_Module()
      call Dispose_Guess_Module()
      call Dispose_Solver_Module()
      call Dispose_Solver_Module_CP8()
      call Dispose_BlockUtils_Module()
      call Dispose_Analyzer_Module()
      call Dispose_Updater_Module()

      end subroutine Show_section_times

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Show_module_times()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints timings of modules

      implicit none

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(/X,A,29X,3(7X,A))') &
         'Routine:','Min time','Max time','Ave time'
      ENDIF

      call Dispose_Eigen_Module()
      call Dispose_MVV_Module()
      call Dispose_CPMM_Module()
      call Dispose_Reduction_Module()
      call Dispose_ALSPOW_Module
      call Dispose_ALSUtils_Module()
      call Dispose_LinSolver_Module()

      call Dispose_CPr8_Module()
      call Dispose_CPMM8_Module()
      call Dispose_CPVV8_Module()
      call Dispose_CPLS8_Module()
      call Dispose_LinAlg8_Module()
      call Dispose_Munkres_Module()


      end subroutine Show_module_times

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TIMINGS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
