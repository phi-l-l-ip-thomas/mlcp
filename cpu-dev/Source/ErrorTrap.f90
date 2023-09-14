!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ERRORTRAP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        CONTAINS:           SetupErrorTrapModule()
!                            DisposeErrorTrapModule()
!                            ERROR(Test,Message)
!                            WARN(Test,Message)
!                            AbortWithError(Message)
!                            ShowWarning(Message)
!                            ShowError(Message)
!        AUTHOR(S):          1st: M.F.Somers, may 2002.
!***********************************************************************

      INTEGER :: ErrorTrapModuleSetupFlag=0

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupErrorTrapModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up the module

      implicit none

      ErrorTrapModuleSetupFlag=ErrorTrapModuleSetupFlag+1

      IF (ErrorTrapModuleSetupFlag>1) RETURN

      end subroutine SetupErrorTrapModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeErrorTrapModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes the module

      implicit none
    
      IF (ErrorTrapModuleSetupFlag<=0) RETURN

      ErrorTrapModuleSetupFlag=ErrorTrapModuleSetupFlag-1

      IF (ErrorTrapModuleSetupFlag>0) RETURN

      end subroutine DisposeErrorTrapModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowError(Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays an error message to stdout

      implicit none
      CHARACTER(*) :: Message

      IF (ErrorTrapModuleSetupFlag<= 0) CALL SetupErrorTrapModule()

      PRINT '(/A,A/)','ERROR: ',TRIM(Message)

      end subroutine ShowError

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowWarning(Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays a warning message to stdout

      implicit none
      CHARACTER(*) :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      PRINT '(/A,A/)','WARNING: ',TRIM(Message)

      end subroutine ShowWarning

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AbortWithError(String)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays message to stdout and stops program

      implicit none
      CHARACTER(*) :: String

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      CALL ShowError(String)

      STOP

      end subroutine AbortWithError

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ERROR(Test,Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests the given test and shows error message and stop program if true

      implicit none
      LOGICAL       :: Test
      CHARACTER(*)  :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      IF (Test) CALL AbortWithError(Message)

      end subroutine ERROR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WARN(Test,Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests the given test and prints warning message if true

      implicit none
      LOGICAL       :: Test
      CHARACTER(*)  :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      IF (Test) CALL ShowWarning(Message)

      end subroutine WARN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ERRORTRAP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
