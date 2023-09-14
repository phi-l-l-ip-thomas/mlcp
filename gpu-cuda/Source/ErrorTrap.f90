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

      SUBROUTINE SetupErrorTrapModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets up the module

      IMPLICIT NONE

      ErrorTrapModuleSetupFlag=ErrorTrapModuleSetupFlag+1

      IF (ErrorTrapModuleSetupFlag>1) RETURN

      END SUBROUTINE SetupErrorTrapModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE DisposeErrorTrapModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes the module

      IMPLICIT NONE
    
      IF (ErrorTrapModuleSetupFlag<=0) RETURN

      ErrorTrapModuleSetupFlag=ErrorTrapModuleSetupFlag-1

      IF (ErrorTrapModuleSetupFlag>0) RETURN

      END SUBROUTINE DisposeErrorTrapModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ShowError(Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays an error message to stdout

      IMPLICIT NONE
      CHARACTER(*) :: Message

      IF (ErrorTrapModuleSetupFlag<= 0) CALL SetupErrorTrapModule()

      PRINT '(/A,A/)','ERROR: ',TRIM(Message)

      END SUBROUTINE ShowError

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ShowWarning(Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays a warning message to stdout

      IMPLICIT NONE
      CHARACTER(*) :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      PRINT '(/A,A/)','WARNING: ',TRIM(Message)

      END SUBROUTINE ShowWarning

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE AbortWithError(String)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Displays message to stdout and stops program

      IMPLICIT NONE
      CHARACTER(*) :: String

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      CALL ShowError(String)

      STOP

      END SUBROUTINE AbortWithError

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ERROR(Test,Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests the given test and shows error message and stop program if true

      IMPLICIT NONE
      LOGICAL       :: Test
      CHARACTER(*)  :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      IF (Test) CALL AbortWithError(Message)

      END SUBROUTINE ERROR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE WARN(Test,Message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests the given test and prints warning message if true

      IMPLICIT NONE
      LOGICAL       :: Test
      CHARACTER(*)  :: Message

      IF (ErrorTrapModuleSetupFlag<=0) CALL SetupErrorTrapModule()

      IF (Test) CALL ShowWarning(Message)

      END SUBROUTINE WARN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ERRORTRAP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
