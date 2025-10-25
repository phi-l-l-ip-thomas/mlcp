!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ALS8DRVR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Drivers for ALS object-oriented code

      USE ERRORTRAP
      USE UTILS
      USE MYACC
      USE CPr8
      USE CPVV8
      USE CPMM8
      USE CPLS8
      USE ALSOO8

      implicit none
      real(kind=8), private :: als_penalty=-1.d0
      character(len=64), private :: als_solver='uninitialized'

      INTERFACE ALS_reduce_CP8
         MODULE PROCEDURE ALS_reduce_wrap_CP8,ALS_reduce_X_CP8
         MODULE PROCEDURE ALS_Wreduce_wrap_CP8,ALS_Wreduce_X_CP8
      END INTERFACE ALS_reduce_CP8

      INTERFACE ALS_fastsolve_CP8
         MODULE PROCEDURE ALS_fastsolve_wrap_CP8,ALS_fastsolve_X_CP8
         MODULE PROCEDURE ALS_fastWsolve_wrap_CP8,ALS_fastWsolve_X_CP8
      END INTERFACE ALS_fastsolve_CP8

      INTERFACE ALS_fastsolve_tiled_CP8
         MODULE PROCEDURE ALS_fastsolve_tiled_wrap_CP8
         MODULE PROCEDURE ALS_fastsolve_tiled_X_CP8
      END INTERFACE ALS_fastsolve_tiled_CP8

      INTERFACE ALS_solve_CP8
         MODULE PROCEDURE ALS_solve_wrap_CP8,ALS_solve_X_CP8
         MODULE PROCEDURE ALS_Wsolve_wrap_CP8,ALS_Wsolve_X_CP8
      END INTERFACE ALS_solve_CP8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine set_als_settings_ALS8(penalty,solver)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets the ALS regularization penalty for the module

      implicit none
      real(kind=8), intent(in) :: penalty
      character(len=64), intent(in) :: solver

      if (penalty.gt.1.d0 .or. penalty.lt.0.d0) then
         write(*,'(A,ES11.4,A)') 'ALS regularization penalty ',&
         penalty,' must be in range: 0 <= penalty <= 1'
         call AbortWithError('set_als_settings_ALS8(): wrong value')
      endif

      als_penalty=penalty
      als_solver=solver

      end subroutine set_als_settings_ALS8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_wrap_CP8(F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: G
      integer, intent(in) :: nals,algo
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_reduce_CP8(X,F,G,nals,algo,nm)
      call X%flush()

      end function ALS_reduce_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_X_CP8(X,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS rank reduction

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: G
      TYPE (CP8) :: dummy
      integer, intent(in) :: nals,algo
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      call nvtx_start('ALS_reduce')

      ndof=F%D()
      conv=0

      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(dummy,0,0.d0,F,G,dummy,.FALSE.,.FALSE.,.FALSE.,algo)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
      call X%setlogical('calcGG',.FALSE.) ! .F. safe for large case
!      call X%setlogical('prtconv',.TRUE.)

      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(F,F,dummy,G)
            conv=X%checkconv(F,G,G)
            IF (conv.gt.0) EXIT
            call X%downdateBP(F,F,dummy,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call X%updateBP(F,F,dummy,G,d)
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%resetconv()
      call X%accumprods(F,F,dummy,G)
      conv=X%checkconv(F,G,G)

      call nvtx_stop()

      end function ALS_reduce_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wreduce_wrap_CP8(F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: G,W
      integer, intent(in) :: nals,algo
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_reduce_CP8(X,F,G,W,nals,algo,nm)
      call X%flush()

      end function ALS_Wreduce_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wreduce_X_CP8(X,F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS rank reduction with weights

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: G,W
      TYPE (CP8) :: dummy,WF,WG
      TYPE (MM8) :: uWF
      integer, intent(in) :: nals,algo
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      call nvtx_start('ALS_Wreduce')

      ndof=F%D()
      conv=0

      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(dummy,0,0.d0,F,G,W,.FALSE.,.FALSE.,.FALSE.,algo)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
      call X%setlogical('calcGG',.FALSE.) ! .F. safe for large case
!      call X%setlogical('prtconv',.TRUE.)

!     Multiply G by weights: W*G = WG
      call CPMM_CP8(W,.FALSE.,G,.FALSE.,WG,.TRUE.,algo)
!     Multiply F by weights: W*F = WF; create object for repeated use
      call uWF%new(W,0,0.0,.FALSE.,F,0,0.0,.FALSE.,WF,.FALSE.,algo)
      call CPMM_CP8(uWF,W,F,WF)

      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(WF,F,dummy,WG)
            conv=X%checkconv(F,WG,G)
            IF (conv.gt.0) EXIT
            call X%downdateBP(WF,F,dummy,WG,d)
            call X%constls(d)
            call X%solvels(F,d)
            call CPMM_CP8(uWF,W,F,WF,d)
            call X%updateBP(WF,F,dummy,WG,d)
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%resetconv()
      call X%accumprods(WF,F,dummy,WG)
      conv=X%checkconv(F,WG,G)

      call WG%flush()
      call WF%flush()
      call uWF%flush()
      call nvtx_stop()

      end function ALS_Wreduce_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastsolve_alg_CP8(A,ish,Esh,F,G,nals,algo,tiled,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining"), H-tiled version.
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (CP8), intent(in) :: A,G
      TYPE (CP8), intent(inout) :: F
      integer, intent(in) :: ish,nals,algo
      logical, intent(in) :: tiled
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      IF (tiled) THEN
         conv=ALS_fastsolve_tiled_CP8(A,ish,Esh,F,G,nals,algo,nm)
      ELSE
         conv=ALS_fastsolve_CP8(A,ish,Esh,F,G,nals,algo,nm)
      ENDIF

      end function ALS_fastsolve_alg_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastsolve_wrap_CP8(A,ish,Esh,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_fastsolve_CP8(X,A,ish,Esh,F,G,nals,algo,nm)
      call X%flush

      end function ALS_fastsolve_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastsolve_X_CP8(X,A,ish,Esh,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: A,G
      TYPE (CP8) :: dummy,AF
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in)  :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      ndof=F%D()
      conv=0

      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(A,ish,Esh,F,G,dummy,.FALSE.,.FALSE.,.FALSE.,algo,.FALSE.)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
!      call X%setlogical('prtconv',.TRUE.)

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)

!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call X%recalcnormal(A,ish,Esh,dummy,G,d)
            call X%resetconv()
            call X%accumprods(AF,F,F,G)
            conv=X%checkconv(F,G,G)
            call X%downdateBP(AF,F,F,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call CPMM_CP8(uAF,A,F,AF,d)
            call X%updateBP(AF,F,F,G,d)
         ENDDO

         call X%resetconv()
         call X%accumprods(AF,F,F,G)
         conv=X%checkconv(AF,G,G)
      ENDDO

      call X%resetconv()
      call X%accumprods(AF,F,F,G)
      conv=X%checkconv(F,G,G)

      call AF%flush()
      call uAF%flush()

      end function ALS_fastsolve_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastsolve_tiled_wrap_CP8(A,ish,Esh,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_fastsolve_tiled_CP8(X,A,ish,Esh,F,G,nals,algo,nm)
      call X%flush

      end function ALS_fastsolve_tiled_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastsolve_tiled_X_CP8(X,A,ish,Esh,F,G,nals,algo,nm) &
      result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining"), H-tiled version.
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G
      character(len=*), intent(in), optional :: nm
      TYPE (CP8) :: A0,A1,AF,dummy
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: TEsh
      integer :: reqs(4)
      integer :: torank,fromrank,ntiles
      integer :: i,d,ndof,conv,j,ierr
      logical :: swapA

      call nvtx_start('ALS_fastsolve_tiled')

      ndof=F%D()
      ntiles=mpinodes

!     Make 2 copies of A, for double-buffering
      call A0%copyfrom(A)
      call A1%copyfrom(A)

!     Shift divided by ntiles since each matrix-vector product must
!     apply the shift in the same way
      TEsh=Esh/mpinodes

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,TEsh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
!     Need "dummy" AF here to create ALS object below
      call CPMM_CP8(uAF,A,F,AF)

!     Create the ALS object
      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(A,ish,TEsh,F,G,dummy,.TRUE.,.FALSE.,.FALSE.,algo,.FALSE.)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('calcGG',.FALSE.)
      call X%setlogical('chkconv',.TRUE.)
!      call X%setlogical('prtconv',.TRUE.)

      swapA=.FALSE.
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            DO j=1,ntiles
!              Source and destination ranks
               torank=mod(mpirank+mpinodes-1,mpinodes)
               fromrank=mod(mpirank+1,mpinodes)
               
               if (swapA) then ! Use A1, send A1, receive in A0
                  call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A1,F,AF)
                   call X%recalcnormal(A1,ish,TEsh,dummy,G,d)
               else            ! Use A0, send A0, receive in A1
                  call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A0,F,AF)
                  call X%recalcnormal(A0,ish,TEsh,dummy,G,d)
               endif
               
               call X%resetdofstate()
               call X%accumprods(AF,F,F,G)
               call X%downdateBP(AF,F,F,G,d)
               call X%constls(d)

               call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
               if (ierr.ne.0) &
                  write(*,*) 'rank ',mpirank,': ierr = ',ierr
               swapA=(.not.swapA)
            ENDDO
            conv=X%checkconv(F,G,G)
            call X%solvels(F,d)
         ENDDO
      ENDDO

!     Final convergence check
      call X%resetconv()
      do j=1,ntiles
         torank=mod(mpirank+mpinodes-1,mpinodes)
         fromrank=mod(mpirank+1,mpinodes)

         if (swapA) then ! Use A1, send A1, receive in A0
            call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A1,F,AF)
         else            ! Use A0, send A0, receive in A1
            call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A0,F,AF)
         endif

         call X%resetdofstate()
         call X%accumprods(AF,F,F,G)

         call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
         if (ierr.ne.0) &
            write(*,*) 'rank ',mpirank,': ierr = ',ierr
            swapA=(.not.swapA)
         conv=X%checkconv(F,G,G)
      enddo

      call AF%flush()
      call uAF%flush()
      call A1%flush()
      call A0%flush()
      call nvtx_stop()

      end function ALS_fastsolve_tiled_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastWsolve_wrap_CP8(A,ish,Esh,F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G,W
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_fastsolve_CP8(X,A,ish,Esh,F,G,W,nals,algo,nm)
      call X%flush

      end function ALS_fastWsolve_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_fastWsolve_X_CP8(X,A,ish,Esh,F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver with weights

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: A,G,W
      TYPE (CP8) :: AF,WAF,WG
      TYPE (MM8) :: uAF,uWAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      ndof=F%D()
      conv=0

!     A^T*A and A^T*G are computed here
      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(A,ish,Esh,F,G,W,.FALSE.,.FALSE.,.FALSE.,algo,.FALSE.)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
!      call X%setlogical('prtconv',.TRUE.)

!     Multiply G by weights: W*G = WG
      call CPMM_CP8(W,.FALSE.,G,.FALSE.,WG,.TRUE.,algo)
!     Multiply F by A: A*F = AF and
!                      W*AF = WAF; create objects for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)
      call uWAF%new(W,0,0.d0,.FALSE.,AF,0,0.d0,.FALSE.,WAF,.FALSE.,algo)
      call CPMM_CP8(uWAF,W,AF,WAF)

!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(WAF,F,F,WG)
            conv=X%checkconv(F,WG,G)
            IF (conv.gt.0) EXIT
            call X%downdateBP(WAF,F,F,WG,d)
            call X%recalcnormal(A,ish,Esh,W,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call CPMM_CP8(uAF,A,F,AF,d)
            call CPMM_CP8(uWAF,W,AF,WAF,d)
            call X%updateBP(WAF,F,F,WG,d)
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%resetconv()
      call X%accumprods(WAF,F,F,WG)
      conv=X%checkconv(F,WG,G)

      call WG%flush()
      call AF%flush()
      call WAF%flush()
      call uAF%flush()
      call uWAF%flush()

      end function ALS_fastWsolve_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_solve_wrap_CP8(A,ish,Esh,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_solve_CP8(X,A,ish,Esh,F,G,nals,algo,nm)
      call X%flush

      end function ALS_solve_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_solve_X_CP8(X,A,ish,Esh,F,G,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: A,G
      TYPE (CP8) :: dummy,AF
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in)  :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      call nvtx_start('ALS_solve')

      ndof=F%D()
      conv=0

      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(A,ish,Esh,F,G,dummy,.FALSE.,.FALSE.,.FALSE.,algo)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
!      call X%setlogical('prtconv',.TRUE.)

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)

!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(AF,AF,F,G)
            conv=X%checkconv(F,G,G)
            IF (conv.gt.0) EXIT
            call X%downdateBP(AF,AF,F,G,d)
            call X%recalcnormal(A,ish,Esh,dummy,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call CPMM_CP8(uAF,A,F,AF,d)
            call X%updateBP(AF,AF,F,G,d)
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%resetconv()
      call X%accumprods(AF,AF,F,G)
      conv=X%checkconv(F,G,G)

      call AF%flush()
      call uAF%flush()
      call nvtx_stop()

      end function ALS_solve_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wsolve_wrap_CP8(A,ish,Esh,F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8), intent(in) :: A,G,W
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: conv

      conv=ALS_solve_CP8(X,A,ish,Esh,F,G,W,nals,algo,nm)
      call X%flush

      end function ALS_Wsolve_wrap_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wsolve_X_CP8(X,A,ish,Esh,F,G,W,nals,algo,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver with weights

      implicit none
      TYPE (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      TYPE (CP8), intent(in) :: A,G,W
      TYPE (CP8) :: AF,WAF,WG
      TYPE (MM8) :: uAF,uWAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      call nvtx_start('ALS_Wsolve')
      ndof=F%D()
      conv=0

!     A^T*A and A^T*G are computed here
      IF (.not.ALLOCATED(X%dofincluded)) &
         call X%new(A,ish,Esh,F,G,W,.FALSE.,.FALSE.,.FALSE.,algo)
      IF (present(nm)) call X%setname(nm)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
!      call X%setlogical('prtconv',.TRUE.)

!     Multiply G by weights: W*G = WG
      call CPMM_CP8(W,.FALSE.,G,.FALSE.,WG,.TRUE.,algo)
!     Multiply F by A: A*F = AF and
!                      W*AF = WAF; create objects for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)
      call uWAF%new(W,0,0.d0,.FALSE.,AF,0,0.d0,.FALSE.,WAF,.FALSE.,algo)
      call CPMM_CP8(uWAF,W,AF,WAF)

!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(WAF,AF,F,WG)
            conv=X%checkconv(F,WG,G)
            IF (conv.gt.0) EXIT
            call X%downdateBP(WAF,AF,F,WG,d)
            call X%recalcnormal(A,ish,Esh,W,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call CPMM_CP8(uAF,A,F,AF,d)
            call CPMM_CP8(uWAF,W,AF,WAF,d)
            call X%updateBP(WAF,AF,F,WG,d)
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%resetconv()
      call X%accumprods(WAF,AF,F,WG)
      conv=X%checkconv(F,WG,G)

      call WG%flush()
      call AF%flush()
      call WAF%flush()
      call uAF%flush()
      call uWAF%flush()
      call nvtx_stop()

      end function ALS_Wsolve_X_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_pow_alg_CP8(A,ish,Esh,F,nals,algo,tiled) result(RQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining"), H-tiled version.
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (CP8), intent(in) :: A
      TYPE (CP8), intent(inout) :: F
      integer, intent(in) :: ish,nals,algo
      logical, intent(in) :: tiled
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: RQ

      IF (tiled) THEN
         RQ=ALS_pow_tiled_CP8(A,ish,Esh,F,nals,algo)
      ELSE
         RQ=ALS_pow_CP8(A,ish,Esh,F,nals,algo)
      ENDIF

      end function ALS_pow_alg_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_pow_CP8(A,ish,Esh,F,nals,algo) result(RQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining").
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(in) :: A
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8) :: AF,dummy
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: RQ
      integer :: i,d,ndof,conv

      call nvtx_start('ALS_pow')

      ndof=F%D()

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)

!     Create the ALS object
      call X%new(dummy,0,0.d0,F,AF,dummy,.FALSE.,.TRUE.,.FALSE.,algo)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('chkconv',.TRUE.)
      call X%setlogical('calcGG',.FALSE.)
!      call X%setlogical('prtconv',.TRUE.)

      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            call X%accumprods(F,F,dummy,AF)
            conv=X%checkconv(F,AF,AF)
            call X%downdateBP(F,F,dummy,AF,d)
            call X%recalcnormal(dummy,0,0.d0,dummy,AF,d)
            call X%constls(d)
            call X%solvels(F,d)
            IF (d.eq.ndof) call Normalize_CP8(F,algo)
            call CPMM_CP8(uAF,A,F,AF,d)
            call X%updateBP(F,F,dummy,AF,d)
         ENDDO
!         write(*,*) 'itn = ',i,'; RQ = ',(X%FG+ish*Esh*X%AFAF)/X%AFAF
      ENDDO

      call X%resetconv()
      call X%accumprods(F,F,dummy,AF)
      conv=X%checkconv(F,AF,AF)

      RQ=(X%FG+ish*Esh*X%AFAF)/X%AFAF

      call AF%flush()
      call uAF%flush()
      call X%flush
      call nvtx_stop()

      end function ALS_pow_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_pow_tiled_CP8(A,ish,Esh,F,nals,algo) result(RQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining"), H-tiled version.
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(in) :: A
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8) :: A0,A1,AF,dummy
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: RQ,TEsh
      integer :: reqs(4)
      integer :: torank,fromrank,ntiles
      integer :: i,d,ndof,conv,j,ierr
      logical :: swapA

      call nvtx_start('ALS_pow_tiled')

      ndof=F%D()
      ntiles=mpinodes

!     Make 2 copies of A, for double-buffering
      call A0%copyfrom(A)
      call A1%copyfrom(A)

!     Shift divided by ntiles since each matrix-vector product must
!     apply the shift in the same way
      TEsh=Esh/mpinodes

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,TEsh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
!     Need "dummy" AF here to create ALS object below
      call CPMM_CP8(uAF,A,F,AF)

!     Create the ALS object
      call X%new(dummy,0,0.d0,F,AF,dummy,.FALSE.,.TRUE.,.FALSE.,algo)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('calcGG',.FALSE.)
      call X%setlogical('chkconv',.TRUE.)

      swapA=.FALSE.
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            DO j=1,ntiles
!              Source and destination ranks
               torank=mod(mpirank+mpinodes-1,mpinodes)
               fromrank=mod(mpirank+1,mpinodes)

               if (swapA) then ! Use A1, send A1, receive in A0
                  call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A1,F,AF)
               else            ! Use A0, send A0, receive in A1
                  call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A0,F,AF)
               endif
               
               call X%resetdofstate()
               call X%accumprods(F,F,dummy,AF)
               call X%recalcnormal(dummy,0,0.d0,dummy,AF,d)
               call X%downdateBP(F,F,dummy,AF,d)
               call X%constls(d)

               call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
               if (ierr.ne.0) &
                  write(*,*) 'rank ',mpirank,': ierr = ',ierr
               swapA=(.not.swapA)
            ENDDO

            conv=X%checkconv(F,AF,AF)
            call X%solvels(F,d)
            IF (d.eq.ndof) call Normalize_CP8(F,algo)
         ENDDO
      ENDDO

!     Calculate the RQ of F in case subroutine is called with nals=0 and
!     since RQ is now calculated before the solve, not after
      call X%resetconv()
      do j=1,ntiles
         torank=mod(mpirank+mpinodes-1,mpinodes)
         fromrank=mod(mpirank+1,mpinodes)

         if (swapA) then ! Use A1, send A1, receive in A0
            call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A1,F,AF)
         else            ! Use A0, send A0, receive in A1
            call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A0,F,AF)
         endif

         call X%resetdofstate()
         call X%accumprods(F,F,dummy,AF)

         call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
         if (ierr.ne.0) &
            write(*,*) 'rank ',mpirank,': ierr = ',ierr
            swapA=(.not.swapA)
      enddo

      RQ=(X%FG+ish*Esh*X%AFAF)/X%AFAF
!      write(*,*) 'rank :',mpirank,'; RQ = ',RQ,'; FG = ',X%FG,&
!                 'AFAF = ',X%AFAF

      call AF%flush()
      call uAF%flush()
      call X%flush()
      call A1%flush()
      call A0%flush()
      call nvtx_stop()

      end function ALS_pow_tiled_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ProdHV_ALS_tiled_CP8(A,ish,Esh,Fo,F,nals,algo) result(RQ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining"), H-tiled version.
! Call with nals=0 to calculate Rayleigh Quotient of F

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(in) :: A,Fo
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8) :: A0,A1,AF,dummy
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,algo
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: RQ,TEsh
      integer :: reqs(4)
      integer :: torank,fromrank,ntiles
      integer :: i,d,ndof,conv,j,ierr
      logical :: swapA

      call nvtx_start('ProdHV_ALS_tiled')

      ndof=F%D()
      ntiles=mpinodes

!     Make 2 copies of A, for double-buffering
      call A0%copyfrom(A)
      call A1%copyfrom(A)

!     Shift divided by ntiles since each matrix-vector product must
!     apply the shift in the same way
      TEsh=Esh/mpinodes

!     Multiply Fo by A: A*Fo = AF; create object for repeated use
      call uAF%new(A,ish,TEsh,.FALSE.,Fo,0,0.0,.FALSE.,AF,.FALSE.,algo)
!     Need "dummy" AF here to create ALS object below
      call CPMM_CP8(uAF,A,Fo,AF)

!     Create the ALS object
      call X%new(dummy,0,0.d0,Fo,AF,dummy,.FALSE.,.TRUE.,.FALSE.,algo)
      call X%setrealval('penalty',als_penalty)
      call X%setstringval('solver',als_solver)
      call X%setlogical('calcGG',.FALSE.)
      call X%setlogical('chkconv',.TRUE.)

      swapA=.FALSE.
      DO i=1,nals
         DO d=1,ndof
            call X%resetconv()
            DO j=1,ntiles
!              Source and destination ranks
               torank=mod(mpirank+mpinodes-1,mpinodes)
               fromrank=mod(mpirank+1,mpinodes)

               if (swapA) then ! Use A1, send A1, receive in A0
                  call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A1,Fo,AF)
               else            ! Use A0, send A0, receive in A1
                  call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
                  call CPMM_CP8(uAF,A0,Fo,AF)
               endif
               
               call X%resetdofstate()
               call X%accumprods(F,F,dummy,AF)
               call X%recalcnormal(dummy,0,0.d0,dummy,AF,d)
               call X%downdateBP(F,F,dummy,AF,d)
               call X%constls(d)

               call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
               if (ierr.ne.0) &
                  write(*,*) 'rank ',mpirank,': ierr = ',ierr
               swapA=(.not.swapA)
            ENDDO

            conv=X%checkconv(F,AF,AF)
            call X%solvels(F,d)
         ENDDO
      ENDDO

!     Calculate the RQ of F in case subroutine is called with nals=0 and
!     since RQ is now calculated before the solve, not after
      call X%resetconv()
      do j=1,ntiles
         torank=mod(mpirank+mpinodes-1,mpinodes)
         fromrank=mod(mpirank+1,mpinodes)

         if (swapA) then ! Use A1, send A1, receive in A0
            call MPI_ISendRecv_CP8(A1,A0,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A1,Fo,AF)
         else            ! Use A0, send A0, receive in A1
            call MPI_ISendRecv_CP8(A0,A1,torank,fromrank,reqs,algo)
            call CPMM_CP8(uAF,A0,Fo,AF)
         endif

         call X%resetdofstate()
         call X%accumprods(F,F,dummy,AF)

         call MPI_WAITALL(4,reqs,MPI_STATUSES_IGNORE,ierr)
         if (ierr.ne.0) &
            write(*,*) 'rank ',mpirank,': ierr = ',ierr
            swapA=(.not.swapA)
      enddo

      RQ=(X%FG+ish*Esh*X%AFAF)/X%AFAF
!      write(*,*) 'rank :',mpirank,'; RQ = ',RQ,'; FG = ',X%FG,&
!                 'AFAF = ',X%AFAF

      call AF%flush()
      call uAF%flush()
      call X%flush()
      call A1%flush()
      call A0%flush()
      call nvtx_stop()

      end function ProdHV_ALS_tiled_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_invpow_normal_CP8(A,ish,Esh,F,nals,npow,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for inverse ALS power method.
! TODO: add code to ALSOO so that X only has to be initialized once
!      (also requires recomputing P matrix)

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(in) :: A
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8) :: AF,G,dummy
      TYPE (MM8) :: uAF
      integer, intent(in) :: ish,nals,npow,algo
      real(kind=8), intent(in) :: Esh
      real(kind=8) :: RQ
      integer :: i,k,d,ndof,conv

      call nvtx_start('ALS_invpow_normal')
      ndof=F%D()

      write(*,*) '*----------------*'
      write(*,*) 'Eshift = ',Esh
      write(*,*) '*----------------*'

!     Multiply F by A: A*F = AF; create object for repeated use
      call uAF%new(A,ish,Esh,.FALSE.,F,0,0.0,.FALSE.,AF,.FALSE.,algo)
      call CPMM_CP8(uAF,A,F,AF)

      DO k=1,npow
         write(*,'(A,I0,A,I0,A)') 'ipow = (',k,'/',npow,')'
!        Create the ALS object
         call G%copyfrom(F)
         call X%new(A,ish,Esh,F,G,dummy,.FALSE.,.FALSE.,.FALSE.,algo)

!         call X%initprods(F,F,AF)
         call X%setrealval('penalty',als_penalty)
         call X%setstringval('solver',als_solver)
         call X%setlogical('chkconv',.TRUE.)
         call X%setlogical('prtconv',.TRUE.)

         DO i=1,nals
            DO d=1,ndof
               call X%recalcnormal(A,ish,Esh,dummy,G,d)
               call X%resetconv()
               call X%accumprods(AF,AF,F,G)
               conv=X%checkconv(AF,G,G)
               call X%downdateBP(AF,AF,F,G,d)
               call X%constls(d)
               call X%solvels(F,d)
               call CPMM_CP8(uAF,A,F,AF,d)
               call X%updateBP(AF,AF,F,G,d)
            ENDDO
         ENDDO

         call X%resetconv()
         call X%accumprods(AF,AF,F,G)
         conv=X%checkconv(AF,G,G)

         call G%flush()
         call X%flush()
      ENDDO ! npow

      call AF%flush()
      call uAF%flush()
      call nvtx_stop()

      end subroutine ALS_invpow_normal_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_invpow_fast_CP8(A,ish,Esh,F,nals,npow,algo,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for inverse ALS power method, non-normal form (experimental).
! TODO: add code to ALSOO so that X only has to be initialized once
!      (also requires recomputing P matrix)

      implicit none
      TYPE (ALS8) :: X
      TYPE (CP8), intent(in) :: A
      TYPE (CP8), intent(inout) :: F
      TYPE (CP8) :: G
      integer, intent(in) :: ish,nals,npow,algo
      logical, intent(in) :: tiledH
      real(kind=8), intent(in) :: Esh
      integer :: k,conv

      call nvtx_start('ALS_invpow_fast')

      write(*,*) '*----------------*'
      write(*,*) 'Eshift = ',Esh
      write(*,*) '*----------------*'

      DO k=1,npow
         write(*,'(A,I0,A,I0,A)') 'ipow = (',k,'/',npow,')'
         call G%copyfrom(F)
         conv=ALS_fastsolve_alg_CP8(A,ish,Esh,F,G,nals,algo,tiledH)
         call G%flush()
         call Normalize_CP8(F,algo)
      ENDDO
      call nvtx_stop()

!      write(*,*) 'filtered <F,F> = ',PRODVV_CP8(F,F,algo)

      end subroutine ALS_invpow_fast_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ALS8DRVR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
