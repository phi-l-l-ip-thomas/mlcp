!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ALSDRVR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Drivers for ALS object-oriented code

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE LINALG
      USE MODVECVEC
      USE CPMMM
      USE ALSOO

      INTERFACE ALS_trials
         MODULE PROCEDURE ALS_reduce_trials
         MODULE PROCEDURE ALS_Wreduce_trials
      END INTERFACE ALS_trials

      INTERFACE ALS_reduce
         MODULE PROCEDURE ALS_reduce_wrap,ALS_reduce_X
         MODULE PROCEDURE ALS_Wreduce_wrap,ALS_Wreduce_X
      END INTERFACE ALS_reduce

      INTERFACE ALS_solve
         MODULE PROCEDURE ALS_solve_wrap,ALS_solve_X
         MODULE PROCEDURE ALS_Wsolve_wrap,ALS_Wsolve_X
      END INTERFACE ALS_solve

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_trials(rk,G,ntrials,nals,nm) result(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction with weights. This runs
! multiple trials with random starting points and keeps the best

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(in) :: G
      TYPE (CP) :: F,Ftmp
      integer, intent(in) :: rk,ntrials,nals
      integer :: i,conv
      real*8  :: conver
      character(len=*), intent(in), optional :: nm

      conver=1.d99
      F=RandomCP(G,rk)

      DO i=1,ntrials
         Ftmp=RandomCP(G,rk)
         call NewALS(X,Ftmp,G)
         IF (present(nm)) call X%setname(nm)
         conv=ALS_reduce(X,Ftmp,G,nals)
         if (X%conver.lt.conver) then
            call ReplaceVwithW(F,Ftmp)
            conver=X%conver
         endif
         call X%Flush
      ENDDO
      write(*,*) 'Best ||F-G||/||G|| found: ',conver

      end function ALS_reduce_trials

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_adaptive(rk,G,conver,nals,nm) result(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes F from ranks 1,2...rk until threshold

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(in) :: G
      TYPE (CP) :: F
      integer, intent(in) :: rk,nals
      real*8, intent(in)  :: conver
      integer :: i,conv,mxrk
      character(len=*), intent(in), optional :: nm

      mxrk=min(rk,G%R())
      F=RandomCP(G,1)

      DO i=1,mxrk
         call NewALS(X,F,G)
         IF (present(nm)) call X%setname(nm)
         write(*,'(A,I0,A,I0,A)') 'rank(',i,'/',G%R(),')'
         conv=ALS_reduce(X,F,G,nals)
!         call PrintCPmat_gen(F,.FALSE.)
         call X%Flush
         if (X%conver.lt.conver .or. F%R().ge.mxrk) exit
         call AugmentVWithRandom(F,i+1)
!         call FlushCP(F)
!         F=RandomCP(G,i+1)
      ENDDO

      end function ALS_reduce_adaptive

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_wrap(F,G,nals,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(inout) :: F
      TYPE (CP), intent(in) :: G
      integer, intent(in) :: nals
      character(len=*), intent(in), optional :: nm
      integer :: conv

      call NewALS(X,F,G)
      IF (present(nm)) call X%setname(nm)
      conv=ALS_reduce(X,F,G,nals)
      call X%Flush

      end function ALS_reduce_wrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_reduce_X(X,F,G,nals,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS rank reduction

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F
      TYPE (CP), intent(in) :: G
      integer, intent(in) :: nals
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      ndof=F%D()
      conv=0

      IF (.not.ALLOCATED(X%dofincluded)) call NewALS(X,F,G)
      IF (present(nm)) call X%setname(nm)
      call X%setoption('prtconv',.FALSE.)

      DO i=1,nals
         DO d=1,ndof
!           Print the last iteration only
            call ProdMatsDD(X,F,G,d)
            call SolveLS(X,F,G,d)
            call ProdMatsUD(X,F,G,d,conv)
            IF (conv.gt.0) EXIT
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call X%showconv()

      end function ALS_reduce_X

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wreduce_trials(rk,G,W,ntrials,nals,nm) result(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction with weights. This runs
! multiple trials with random starting points and keeps the best

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(in) :: G,W
      TYPE (CP) :: F,Ftmp
      integer, intent(in) :: rk,ntrials,nals
      integer :: i,conv
      real*8  :: conver
      character(len=*), intent(in), optional :: nm

      conver=1.d99
      F=RandomCP(G,rk)

      DO i=1,ntrials
         Ftmp=RandomCP(G,rk)
         call NewALS(X,Ftmp,G,W)
         IF (present(nm)) call X%setname(nm)
         conv=ALS_reduce(X,Ftmp,G,W,nals)
         if (X%conver.lt.conver) then
            call ReplaceVwithW(F,Ftmp)
            conver=X%conver
            write(*,'(A,ES14.8)') ' * TRIALS: new best: ',X%conver
         else
            write(*,'(A,ES14.8)') ' * TRIALS: dud...  : ',X%conver
         endif
         call X%Flush
      ENDDO

      end function ALS_Wreduce_trials

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wreduce_wrap(F,G,W,nals,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS rank reduction with weights

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(inout) :: F
      TYPE (CP), intent(in) :: G,W
      integer, intent(in) :: nals
      character(len=*), intent(in), optional :: nm
      integer :: conv

      call NewALS(X,F,G,W)
      IF (present(nm)) call X%setname(nm)
      conv=ALS_reduce(X,F,G,W,nals)
      call X%Flush

      end function ALS_Wreduce_wrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wreduce_X(X,F,G,W,nals,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS rank reduction with weights

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F
      TYPE (CP), intent(in) :: G,W
      TYPE (CP) :: WF,WG
      integer, intent(in) :: nals
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,conv

      ndof=F%D()
      conv=0

!     Multiply weights with F,G
      call CPMM(W,.FALSE.,F,.FALSE.,WF)
      call CPMM(W,.FALSE.,G,.FALSE.,WG)

      IF (.not.ALLOCATED(X%dofincluded)) call NewALS(X,F,G,W)
      IF (present(nm)) call X%setname(nm)

!      call X%setoption('prtconv',.FALSE.)

      DO i=1,nals
         DO d=1,ndof
            call ProdMatsDD(X,F,WF,G,WG,d)
            call SolveLS(X,F,WG,W,d)
            call CPMM(W,0,0.d0,.FALSE.,F,0,0.d0,.FALSE.,WF,d)
            call ProdMatsUD(X,F,WF,G,WG,d,conv)
            IF (conv.gt.0) EXIT
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call FlushCP(WF)
      call FlushCP(WG)

      end function ALS_Wreduce_X

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_solve_adaptive(rk,A,G,conver,nals,ishift,Eshift,nm) &
               result(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes F from ranks 1,2...rk until threshold

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(in) :: A,G
      TYPE (CP) :: F,AF
      integer, intent(in) :: rk,nals
      real*8, intent(in)  :: conver
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      integer :: i,ish,conv!,mxrk
      real*8  :: Esh
      character(len=*), intent(in), optional :: nm

      F=RandomCP(1,A%rows,G%cols,G%sym)

      ish=0
      Esh=0.d0

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
         ENDIF
      ENDIF

      DO i=1,rk
         write(*,*) '   rank = ',i
         call NewLinSolver(X,A,F,G,.FALSE.,.FALSE.,ish,Esh)
         IF (present(nm)) call X%setname(nm)
         call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF)
         conv=ALS_solve(X,A,F,AF,G,nals,ish,Esh)
         call X%Flush
         call FlushCP(AF)
         if (X%conver.lt.conver .or. F%R().ge.rk) exit
         call AugmentVWithRandom(F,i+1)
      ENDDO

      end function ALS_solve_adaptive

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_solve_wrap(A,F,G,nals,ishift,Eshift,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(inout) :: F
      TYPE (CP), intent(in) :: A,G
      TYPE (CP) :: AF
      integer, intent(in) :: nals
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      character(len=*), intent(in), optional :: nm
      integer :: ish,conv
      real*8  :: Esh

      ish=0
      Esh=0.d0

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
         ENDIF
      ENDIF

      call NewLinSolver(X,A,F,G,.FALSE.,.FALSE.,ish,Esh)
      IF (present(nm)) call X%setname(nm)
      call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF)
      conv=ALS_solve(X,A,F,AF,G,nals,ish,Esh)
      call X%flush
      call FlushCP(AF)

      end function ALS_solve_wrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_solve_X(X,A,F,AF,G,nals,ishift,Eshift,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F,AF
      TYPE (CP), intent(in) :: A,G
      integer, intent(in) :: nals
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,ish,conv
      real*8 :: Esh

      ndof=SIZE(F%nbas)
      conv=0
      ish=0
      Esh=0.d0

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
         ENDIF
      ENDIF

!     A^T*A and A^T*G are computed here
      IF (.not.ALLOCATED(X%dofincluded)) &
         call NewLinSolver(X,A,F,G,.FALSE.,.FALSE.,ish,Esh)
      IF (present(nm)) call X%setname(nm)

!     A*F = AF
      IF (.not.ALLOCATED(AF%coef)) &
         call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call RecomputeNormalEquations(X,A,ish,Esh,G,d)
            call ProdMatsDD(X,AF,G,d)
            call SolveLS(X,F,G,d)
            call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call ProdMatsUD(X,AF,G,d,conv)
!            IF (conv.gt.0) EXIT
         ENDDO
!         IF (conv.gt.0) EXIT
      ENDDO

      end function ALS_solve_X

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wsolve_wrap(A,F,G,W,nals,ishift,Eshift,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for driver for ALS linear solver with weights

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(inout) :: F
      TYPE (CP), intent(in) :: A,G,W
      TYPE (CP) :: AF
      integer, intent(in) :: nals
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      character(len=*), intent(in), optional :: nm
      integer :: ish,conv
      real*8  :: Esh

      ish=0
      Esh=0.d0

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
         ENDIF
      ENDIF

      call NewLinSolver(X,A,F,G,W,.FALSE.,.FALSE.,.FALSE.,ish,Esh)
      IF (present(nm)) call X%setname(nm)
      call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF)
      conv=ALS_solve(X,A,F,AF,G,W,nals,ish,Esh)
      call X%flush
      call FlushCP(AF)

      end function ALS_Wsolve_wrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ALS_Wsolve_X(X,A,F,AF,G,W,nals,ishift,Eshift,nm) result(conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS linear solver with weights

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F,AF
      TYPE (CP), intent(in) :: A,G,W
      TYPE (CP) :: WAF,WG
      integer, intent(in) :: nals
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      character(len=*), intent(in), optional :: nm
      integer :: i,d,ndof,ish,conv
      real*8 :: Esh

      ndof=SIZE(F%nbas)
      conv=0
      ish=0
      Esh=0.d0

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
         ENDIF
      ENDIF

!     A^T*A and A^T*G are computed here
      IF (.not.ALLOCATED(X%dofincluded)) &
         call NewLinSolver(X,A,F,G,W,.FALSE.,.FALSE.,.FALSE.,ish,Esh)
      IF (present(nm)) call X%setname(nm)

!     A*F = AF
      IF (.not.ALLOCATED(AF%coef)) &
         call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Multiply weights with A*F,G
      call CPMM(W,.FALSE.,AF,.FALSE.,WAF)
      call CPMM(W,.FALSE.,G,.FALSE.,WG)

      call X%setoption('prtconv',.FALSE.)
      
!     ALS iterations
      DO i=1,nals
         DO d=1,ndof
            call RecomputeNormalEquations(X,A,ish,Esh,WG,W,d)
            call ProdMatsDD(X,AF,WAF,G,WG,d)
            call SolveLS(X,F,WG,W,d)
            call CPMM(A,ish,Esh,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call CPMM(W,0,0.d0,.FALSE.,AF,0,0.d0,.FALSE.,WAF,d)
            call ProdMatsUD(X,AF,WAF,G,WG,d,conv)
            IF (conv.gt.0) EXIT
         ENDDO
         IF (conv.gt.0) EXIT
      ENDDO

      call FlushCP(WAF)
      call FlushCP(WG)

      end function ALS_Wsolve_X

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_pow(A,F,nals,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Driver for ALS power method ("intertwining")

      implicit none
      TYPE (ALS) :: X
      TYPE (CP), intent(in) :: A
      TYPE (CP), intent(inout) :: F
      TYPE (CP) :: AF
      integer, intent(in) :: nals,ishift
      real*8, intent(in) :: Eshift
      integer :: i,d,ndof,conv

      ndof=SIZE(F%nbas)

!     Calc A*F = AF, then create the ALS object
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)
      call NewALS(X,F,AF)
!      call X%show

      DO i=1,nals
         DO d=1,ndof
            call ProdMatsDD(X,F,AF,d)
            call SolveLS(X,F,AF,d)
            IF (d.eq.ndof) call NORMCOEF(F)
            call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call ProdMatsUD(X,F,AF,d,conv)
         ENDDO
      ENDDO

      call FlushCP(AF)
      call X%Flush

      end subroutine ALS_pow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ALSDRVR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
