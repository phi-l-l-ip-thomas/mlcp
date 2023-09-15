!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ALSOO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains object-oriented code for ALS

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE LINALG
      USE MODVECVEC
      USE CPMMM

      TYPE ALS
         TYPE (CP) :: ATA,ATG
         REAL*8, ALLOCATABLE :: B(:,:),P(:,:)
         LOGICAL, ALLOCATABLE :: dofincluded(:)
         REAL*8  :: FF,FG,GG,CN,conver,del
         CHARACTER(LEN=64) :: alsnm
         ! Options
         REAL*8  :: penalty ! penalty on LHS to prevent ill-conditioning
         REAL*8  :: thresh  ! fractional tolerance for convergence
         REAL*8  :: dthresh ! fractional derivative tolerance
         LOGICAL :: AisI    ! .T. for ALS, .F. for linear solver
         LOGICAL :: weights ! .T. weights rows by weights matrix
         LOGICAL :: update  ! .T. up/downdates B,P; .F. always builds
         LOGICAL :: useSVD  ! .T. for SVD solver, .F. for LU solver
         LOGICAL :: parall  ! .T. for OpenMP code (not yet implemented)
         LOGICAL :: chkconv ! .T. to check convergence after each update
         LOGICAL :: prtconv ! .T. to print convergence after each update
         LOGICAL :: Achange ! .T. to update A^T*A, A^T*G each call
         LOGICAL :: Gchange ! .T. to update A^T*G or G each call
         LOGICAL :: Wchange ! .T. to update A^T*W*A each call
         CONTAINS
            PROCEDURE :: setname => ALSSetName
            PROCEDURE :: setoption => ALSSetOption_logical
            PROCEDURE :: setvalue => ALSSetOption_real
            PROCEDURE :: show => ShowALSParameters
            PROCEDURE :: showconv => PrintALSConvergence 
            PROCEDURE :: flush => FlushALS
            PROCEDURE :: FlushNormalEquations
      END TYPE ALS

      INTERFACE NewLinSolver
         MODULE PROCEDURE NewLinSolver_noweights
         MODULE PROCEDURE NewLinSolver_weights
      END INTERFACE NewLinSolver

      INTERFACE RecomputeNormalEquations 
         MODULE PROCEDURE RecomputeNormalEquations_noweights
         MODULE PROCEDURE RecomputeNormalEquations_weights
      END INTERFACE RecomputeNormalEquations

      INTERFACE GetLHSforALS
         MODULE PROCEDURE GetLHSforALS_noweights
         MODULE PROCEDURE GetLHSforALS_weights
      END INTERFACE GetLHSforALS

      INTERFACE SolveLS
         MODULE PROCEDURE SolveLS_noweights
         MODULE PROCEDURE SolveLS_weights
      END INTERFACE SolveLS

      INTERFACE ProdMatsDD
         MODULE PROCEDURE ProdMatsDD_noweights
         MODULE PROCEDURE ProdMatsDD_weights
      END INTERFACE ProdMatsDD

      INTERFACE ProdMatsUD
         MODULE PROCEDURE ProdMatsUD_noweights
         MODULE PROCEDURE ProdMatsUD_weights
      END INTERFACE ProdMatsUD

      INTERFACE CheckALSConvergence
         MODULE PROCEDURE CheckALSConvergence_noweights
         MODULE PROCEDURE CheckALSConvergence_weights
      END INTERFACE CheckALSConvergence

!!!!!! TO DO:
!      Print memory usage
!      Auto-compute Ptot/Btot at start
!      Low-mem intertwining options:
!         ALS:
!            -use operator to generate G from F
!            -build rhs directly from F and operator, bypassing P
!         LinSolv:
!            -Recompute AF, use with above to generate LHS, RHS directly
!             (worth doing? LHS itself could be large and unbypassable)
!!!!!!

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewLinSolver_noweights(X,A,F,G,Achange,Gchange,&
                                        ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the ALS object

      implicit none
      TYPE (ALS), intent(out) :: X
      TYPE (CP), intent(in)   :: A,F,G
      logical, intent(in) :: Achange,Gchange
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      integer :: d,ndof,rA,rF,rG,ish
      real*8  :: Esh

      ndof=F%D()
      rA=A%R()
      rF=F%R()
      rG=G%R()
      ish=0
      Esh=0.d0

!     Error checking
      IF (A%D().ne.ndof) &
         call AbortWithError('NewLinSolver(): ndof mismatch A and F')
      IF (G%D().ne.ndof) &
         call AbortWithError('NewLinSolver(): ndof mismatch F and G')
      DO d=1,ndof
         IF (A%cols(d).ne.F%rows(d)) &
            call AbortWithError('NewLinSolver(): A cols != F rows')
         IF (A%rows(d).ne.G%rows(d)) &
            call AbortWithError('NewLinSolver(): A rows != G rows')
         IF (F%cols(d).ne.G%cols(d)) &
            call AbortWithError('NewLinSolver(): F cols != G cols')
      ENDDO

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
            rA=rA+1
         ENDIF
      ENDIF

      allocate(X%B(rA*rF,rA*rF),X%P(rG,rA*rF),X%dofincluded(ndof))
      X%dofincluded(:)=.FALSE.

!     Compute A^T*A if A does not change, A^T*G if A and G do not change
      IF (.not.Achange) THEN
         call CPMM(A,ish,Esh,.TRUE.,A,ish,Esh,.FALSE.,X%ATA)
         IF (.not.Gchange) THEN
            call CPMM(A,ish,Esh,.TRUE.,G,0,0.d0,.FALSE.,X%ATG)
         ENDIF
      ENDIF

      X%Achange=Achange
      X%Gchange=Gchange
      X%Wchange=.FALSE.
      X%AisI=.FALSE.
      X%weights=.FALSE.

      call ALSOptionDefaults(X,ndof)

      end subroutine NewLinSolver_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewLinSolver_weights(X,A,F,G,W,Achange,Gchange,&
                                      Wchange,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the ALS object

      implicit none
      TYPE (ALS), intent(out) :: X
      TYPE (CP), intent(in)   :: A,F,G,W
      TYPE (CP) :: Wmult
      logical, intent(in) :: Achange,Gchange,Wchange
      integer, intent(in), optional :: ishift
      real*8, intent(in), optional  :: Eshift
      integer :: d,ndof,rA,rF,rG,rW,ish
      real*8  :: Esh

      ndof=F%D()
      rA=A%R()
      rF=F%R()
      rG=G%R()
      rW=W%R()
      ish=0
      Esh=0.d0

!     Error checking
      IF (A%D().ne.ndof) &
         call AbortWithError('NewLinSolver(): ndof mismatch A and F')
      IF (G%D().ne.ndof) &
         call AbortWithError('NewLinSolver(): ndof mismatch F and G')
      IF (W%D().ne.ndof) &
         call AbortWithError('NewLinSolver(): ndof mismatch W and F')
      DO d=1,ndof
         IF (A%cols(d).ne.F%rows(d)) &
            call AbortWithError('NewLinSolver(): A cols != F rows')
         IF (A%rows(d).ne.G%rows(d)) &
            call AbortWithError('NewLinSolver(): A rows != G rows')
         IF (F%cols(d).ne.G%cols(d)) &
            call AbortWithError('NewLinSolver(): F cols != G cols')
         IF (W%rows(d).ne.A%rows(d)) &
            call AbortWithError('NewLinSolver(): W rows != A rows')
         IF (W%cols(d).ne.A%rows(d)) &
            call AbortWithError('NewLinSolver(): W cols != A rows')
      ENDDO

      IF (present(ishift)) THEN
         ish=ishift
         IF ((ish.ne.0) .and. present(Eshift)) THEN
            Esh=Eshift
            rA=rA+1
         ENDIF
      ENDIF

      allocate(X%B(rA*rF,rW*rA*rF),X%P(rW*rG,rA*rF),X%dofincluded(ndof))
      X%dofincluded(:)=.FALSE.

!     Calc A^T*W*A if A,W don't change, A^T*G if A,G don't change
!     Since this subroutine accepts A,G and not WA, WG, calc both here
      IF ((.not.Achange) .and. (.not.Wchange)) THEN
         call CPMM(W,0,0.d0,.FALSE.,A,ish,Esh,.FALSE.,Wmult)
         call CPMM(A,ish,Esh,.TRUE.,Wmult,0,0.d0,.FALSE.,X%ATA)
         call FlushCP(Wmult)
      ENDIF
      IF ((.not.Achange) .and. (.not.Gchange) .and. (.not.Wchange)) THEN
         call CPMM(W,0,0.d0,.FALSE.,G,0,0.d0,.FALSE.,Wmult)
         call CPMM(A,ish,Esh,.TRUE.,Wmult,0,0.d0,.FALSE.,X%ATG)
         call FlushCP(Wmult)
      ENDIF

      X%Achange=Achange
      X%Gchange=Gchange
      X%Wchange=Wchange
      X%AisI=.FALSE.
      X%weights=.TRUE.

      call ALSOptionDefaults(X,ndof)

      end subroutine NewLinSolver_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewALS(X,F,G,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the ALS object

      implicit none
      TYPE (ALS), intent(out) :: X
      TYPE (CP), intent(in)   :: F,G
      TYPE (CP), optional, intent(in) :: W
      integer :: d,ndof,rF,rG,rW

      ndof=F%D()
      rF=F%R()
      rG=G%R()
      rW=1
      IF (present(W)) rW=W%R()

!     Error checking
      IF (.NOT.CHECKNBAS(F,G)) &
         call AbortWithError('NewALS(): F,G dimension mismatch')

      allocate(X%B(rF,rW*rF),X%P(rW*rG,rF),X%dofincluded(ndof))
      X%dofincluded(:)=.FALSE.
      call ALSOptionDefaults(X,ndof)
!     Achange, Gchange, Wchange are irrelevant for ALS since A^T*W*A,
!     A^T*G are never computed
      X%Achange=.FALSE. 
      X%Gchange=.FALSE.
      X%Wchange=.FALSE.
      X%AisI=.TRUE.
      X%weights=present(W)

      end subroutine NewALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushALS(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Frees memory from ALS object

      implicit none
      CLASS (ALS), intent(inout) :: X

      X%alsnm=''
      call FlushCP(X%ATA)
      call FlushCP(X%ATG)
      IF (ALLOCATED(X%P)) DEALLOCATE(X%P)
      IF (ALLOCATED(X%B)) DEALLOCATE(X%B)
      IF (ALLOCATED(X%dofincluded)) DEALLOCATE(X%dofincluded)

      end subroutine FlushALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSOptionDefaults(X,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets default values for the parameters for TYPE ALS

      implicit none
      TYPE (ALS), intent(inout) :: X
      integer, intent(in) :: ndof
      logical :: update,useSVD,parall,chkconv

      X%alsnm=''
      X%GG=-1.d0
      X%conver=1.d99
      X%penalty=1.d-15
      X%thresh=8.d-8
      X%dthresh=1.d-5

!     Updating risks zero division, but is cheaper for ndof > 3
      update=(ndof.gt.3)
      usesvd=.TRUE.

!     For now, parallel code not implemented
      parall=.FALSE.
      chkconv=.FALSE.

      call ALSSetOptions(X,update,usesvd,parall,chkconv)

      end subroutine ALSOptionDefaults

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetName(X,nm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS), intent(inout) :: X
      character(len=*), intent(in) :: nm
      
      X%alsnm=TRIM(ADJUSTL(nm))

      end subroutine ALSSetName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOption_logical(X,opt,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS), intent(inout) :: X
      character(len=*), intent(in) :: opt
      logical, intent(in) :: val

      IF (TRIM(ADJUSTL(opt)).eq.'update') THEN
         X%update=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'usesvd') THEN
         X%useSVD=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'parall') THEN
         X%parall=val
         call WARN(val,'Parallel ALS not yet implemented')
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'chkconv') THEN
         X%chkconv=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'prtconv') THEN
         X%prtconv=val
      ELSE
         write(*,*) "Unrecognized logical option: '",opt,"'"
         write(*,*) "Choose from 'update', 'usesvd', "&
                    "'parall', 'chkconv','prtconv'"
         call AbortWithError('ALSSetOption(): unrecognized option')
      ENDIF

      end subroutine ALSSetOption_logical

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOption_real(X,opt,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS), intent(inout) :: X
      character(len=*), intent(in) :: opt
      real*8, intent(in) :: val

      IF (TRIM(ADJUSTL(opt)).eq.'penalty') THEN
         X%penalty=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'thresh') THEN
         X%thresh=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'dthresh') THEN
         X%dthresh=val
      ELSE
         write(*,*) "Unrecognized real option: '",opt,"'"
         write(*,*) "Choose from 'penalty', 'thresh', 'dthresh' "
         call AbortWithError('ALSSetOption(): unrecognized option')
      ENDIF

      end subroutine ALSSetOption_real

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOptions(X,update,usesvd,parall,chkconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      TYPE (ALS), intent(inout) :: X
      logical, intent(in) :: update,usesvd,parall,chkconv

      X%update=update
      X%useSVD=usesvd
      X%parall=parall
      X%chkconv=chkconv
      X%prtconv=chkconv ! Init same as chkconv

      call WARN(parall,'Parallel ALS not yet implemented')

      end subroutine ALSSetOptions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowALSParameters(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints ALS parameters

      implicit none
      CLASS (ALS), intent(in) :: X
      character(len=64) :: frmt
      integer :: d

      write(*,*)
      IF (LEN(TRIM(ADJUSTL(X%alsnm))).gt.0) THEN
         write(*,*) '*** ALS Object Parameters for ',&
         TRIM(ADJUSTL(X%alsnm)),' ***'
      ELSE
         write(*,*) '*** ALS Object Parameters ***'
      ENDIF
      write(*,*)
      write(*,*) 'Parameter  Value  Description'
      write(*,*) 'update  :',X%update,&
      ' (.T. downdates/updates P,B matrices, .F. builds each iteration)'
      write(*,*) 'useSVD  :',X%useSVD,&
      ' (.T. for Moore-Penrose pseudoinverse, .F. for LU decomposition)'
      write(*,*) 'parall  :',X%parall,&
      ' (.T. for OpenMP M*v, v*v products and LHS, RHS building)'
      write(*,*) 'chkconv :',X%chkconv,&
      ' (.T. to check LinSolv,ALS convergence, .F. skips)'
      write(*,*) 'prtconv :',X%prtconv,&
      ' (.T. to print LinSolv,ALS convergence, .F. skips)'
      write(*,*) 'Achange :',X%Achange,&
      ' (.T. if A changes each iteration, .F. for constant A)'
      write(*,*) 'Gchange :',X%Gchange,&
      ' (.T. if G changes each iteration, .F. for constant G)'
      if (X%weights) write(*,*) 'Wchange :',X%Wchange,&
      ' (.T. if W changes each iteration, .F. for constant W)'
      if (.not.X%useSVD) write(*,*) 'penalty :',X%penalty,&
      ' (value added to LHS diagonal to prevent ill-conditioning)'
      if (X%chkconv) then
         write(*,*) 'thresh  :',X%thresh,&
         ' (convergence threshold for ||F-G||/||G||)'
         write(*,*) 'dthresh :',X%dthresh,&
         ' (convergence threshold for abs((conver-X%conver)/conver))'
      endif
      write(*,*)
      write(*,*) 'Arrays:'
      IF (ALLOCATED(X%ATA%coef)) THEN
         write(*,*) 'A^T*A : ALLOCATED, with rank = ',X%ATA%R()
      ELSE
         write(*,*) 'A^T*A : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%ATG%coef)) THEN
         write(*,*) 'A^T*G : ALLOCATED, with rank = ',X%ATG%R()
      ELSE
         write(*,*) 'A^T*G : NOT ALLOCATED'
      ENDIF
      frmt='(X,A,I0,A,I0,A)'
      IF (ALLOCATED(X%B)) THEN
         write(*,frmt) '    B : [',SIZE(X%B,1),' x ',SIZE(X%B,2),']'
      ELSE
         write(*,*) '    B : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%P)) THEN
         write(*,frmt) '    P : [',SIZE(X%P,1),' x ',SIZE(X%P,2),']'
      ELSE
         write(*,*) '    P : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%dofincluded)) THEN
         frmt='(2X,I4,3X,L4)'
         write(*,*) ' mode   incl'
         DO d=1,SIZE(X%dofincluded)
            write(*,frmt) d,X%dofincluded(d)
         ENDDO
      ELSE
         write(*,*) 'dofinc: NOT ALLOCATED'
      ENDIF
      write(*,*) '*****************************'

      end subroutine ShowALSParameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RecomputeNormalEquations_noweights(X,A,ish,Esh,G,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Recomputes A^T*A and A^T*G when A or G changes

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: A,G
      integer, intent(in) :: ish,d
      real*8, intent(in) :: Esh

      IF (X%Achange.or.X%Gchange) THEN
!        Recompute A^T*G for mode d if either A or G changes
         call CPMM(A,ish,Esh,.TRUE.,G,0,0.d0,.FALSE.,X%ATG,d)

!        Recompute A^T*A for mode d only if A changes
         IF (X%Achange) THEN
            call CPMM(A,ish,Esh,.TRUE.,A,ish,Esh,.FALSE.,X%ATA,d)
         ENDIF
      ENDIF

      end subroutine RecomputeNormalEquations_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RecomputeNormalEquations_weights(X,A,ish,Esh,WG,W,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Recomputes A^T*W*A and A^T*(W*G), as necessary, when A,G,W changes. 
! Pass W*G instead of G (with newest W,G) to this routine

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: A,WG,W
      TYPE (CP) :: Wmult
      integer, intent(in) :: ish,d
      real*8, intent(in) :: Esh

      IF (X%Achange.or.X%Gchange.or.X%Wchange) THEN
!        Recompute A^T*W*G for mode d if either A,G,W change
         call CPMM(A,ish,Esh,.TRUE.,WG,0,0.d0,.FALSE.,X%ATG,d)
!         call CPMM(W,0,0.d0,.FALSE.,G,0,0.d0,.FALSE.,Wmult,d)
!         call CPMM(A,ish,Esh,.TRUE.,Wmult,0,0.d0,.FALSE.,X%ATG,d)
!         call FlushCP(Wmult)
      ENDIF

      IF (X%Achange.or.X%Wchange) THEN
!        Recompute A^T*W*A for mode d if either A,W change
         call CPMM(W,0,0.d0,.FALSE.,A,ish,Esh,.FALSE.,Wmult,d)
         call CPMM(A,ish,Esh,.TRUE.,Wmult,0,0.d0,.FALSE.,X%ATA,d)
         call FlushCP(Wmult)
      ENDIF

      end subroutine RecomputeNormalEquations_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushNormalEquations(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates A^T*A and A^T*G after solving linear system for case when 
! A or G changes 

      implicit none
      CLASS (ALS), intent(inout) :: X

      IF (X%Achange.or.X%Gchange) THEN
!        Deallocate A^T*G for mode d as it changes before next use
         call FlushCP(X%ATG)
      ENDIF

      IF (X%Achange.or.X%Wchange) THEN
!        Deallocate A^T*(W)*A for mode d as it changes before next use
         call FlushCP(X%ATA)
      ENDIF

      end subroutine FlushNormalEquations

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetAGmode(X,d,amode,gmode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets correct mode for A^T*(W)*A and A^T*G depending on whether A,G,W 
! change

      implicit none
      TYPE (ALS), intent(in) :: X
      integer, intent(in)  :: d
      integer, intent(out) :: amode,gmode

      IF (X%Achange) THEN
         amode=1
         gmode=1
      ELSE
         IF (X%Wchange) THEN
            amode=1
         ELSE
            amode=d
         ENDIF
         IF (X%Gchange) THEN
            gmode=1
         ELSE
            gmode=d
         ENDIF
      ENDIF

      end subroutine GetAGmode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckModeIncludedState(X,d) result(ok)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns .TRUE. if dofincluded state is .TRUE. for all modes except d,
! which is the desired state for solving the linear system

      implicit none
      TYPE (ALS), intent(in) :: X
      integer, intent(in) :: d
      logical :: ok
      integer :: j

      IF (.NOT.ALLOCATED(X%dofincluded)) THEN
         write(*,*) 'X%dofincluded() array is not allocated'
         call AbortWithError('CheckModeIncludedState(): no dofincluded')
      ELSEIF ((d.lt.1) .or. (d.gt.SIZE(X%dofincluded))) THEN
         write(*,*) 'mode d (',d,') must be in range: [1,'&
                    ,SIZE(X%dofincluded),']'
         call AbortWithError('CheckModeIncludedState(): d out of range')
      ELSE
         ok=(.NOT.X%dofincluded(d))
         DO j=1,SIZE(X%dofincluded)
            IF (j.eq.d) CYCLE
            ok=ok.and.X%dofincluded(j)
         ENDDO
      ENDIF

      end function CheckModeIncludedState

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetLHSforLinSys(X,F,d) result(LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets left hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, for general A.

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in)  :: F
      real*8, allocatable :: LHS(:,:)
      integer, intent(in) :: d
      integer :: rF,N

      rF=SIZE(F%coef)
      N=F%rows(d)

      allocate(LHS(rF*N,rF*N))
      LHS=0.d0

      call AccumulateLHSforLinSys(X,F,d,LHS)

      end function GetLHSforLinSys

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateLHSforLinSys(X,F,d,LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates left hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, for general A.

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in)  :: F
      real*8, intent(inout)  :: LHS(:,:)
      integer, intent(in) :: d
      integer :: i,j,k,l,s,t,u,rW,rA,rF,N,atark,ataind,amode,gmode
      integer :: browa,bcola,brow,bcol,lrowa,lcola,lrow,lcol
      real*8  :: tmp

!     Pick the correct mode to index A^T*W*A
      call GetAGmode(X,d,amode,gmode)

      IF (.not.allocated(X%ATA%coef)) &
         call AbortWithError('AccumulateLHSforLinSys(): no ATWA allocd')

!     Calculate ranks. B is (rA*rF) x (rW*rA*rF)
      rF=F%R()
      rA=SIZE(X%B,1)/rF
      rW=SIZE(X%B,2)/SIZE(X%B,1)
      N=F%rows(d)

!     Error checking
      if ((SIZE(LHS,1).ne.(rF*N)).or.(SIZE(LHS,2).ne.(rF*N))) then
         write(*,*) 'LHS is [',SIZE(LHS,1),' x ',SIZE(LHS,2),&
                 '] but must be [',rF,'*',N,' x ',rF,'*',N,']'
         call AbortWithError('AccumulateLHSforLinSys(): bad dimensions')
      endif

!     Compute the left-hand sides of the linear system
      do s=1,rA
         browa=(s-1)*rF ! row group of B from A^T
         do u=1,rW 
            do t=1,rA
               bcola=((u-1)*rA+t-1)*rF ! col group of B from W*A
               atark=((s-1)*rW+(u-1))*rA+t ! rank of AT*W*A
               do k=1,N
                  lrowa=(k-1)*rF ! row group of lhs from basis
                  do l=1,N
                     lcola=(l-1)*rF ! col group of lhs from basis
                     ataind=X%ATA%ibas(amode)+(l-1)*N+k-1 ! ATWA elem indx
                     tmp=X%ATA%coef(atark)*X%ATA%base(ataind,atark)
                     do i=1,rF
                        brow=browa+i ! row of B from AT,F
                        lrow=lrowa+i ! row of lhs from A,F
                        do j=1,rF
                           bcol=bcola+j ! col of B from W,A,F
                           lcol=lcola+j ! col of lhs from A,F
                           LHS(lrow,lcol)=&
                           LHS(lrow,lcol)+tmp*X%B(brow,bcol)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      end subroutine AccumulateLHSforLinSys

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetLHSforALS_noweights(X,d) result(LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates left hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, optimized for A = I

      implicit none
      TYPE (ALS), intent(in) :: X
      integer, intent(in) :: d
      real*8, allocatable :: LHS(:,:)
      integer :: rF,N

      rF=SIZE(X%B,1)
      allocate(LHS(rF,rF))

!     Just copy B as the left-hand sides of the linear system
      LHS(:,:)=X%B(:,:)

      end function GetLHSforALS_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetLHSforALS_weights(X,F,W,d) result(LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates left hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, optimized for A = I

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in) :: F,W
      integer, intent(in) :: d
      real*8, allocatable :: LHS(:,:)
      integer :: rF,M,N

      rF=SIZE(X%B,1)
      N=W%rows(d) ! Also equals cols-of-W and rows-of-F
      M=F%cols(d)
      allocate(LHS(rF,rF*M*N))

!     Accumulate left-hand sides of the linear system using weights
      LHS=0.d0
      call AccumulateLHSforALS(X,F,W,d,LHS)

      end function GetLHSforALS_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateLHSforALS(X,F,W,d,LHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates left hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, for general A.

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in)  :: F,W
      real*8, intent(inout)  :: LHS(:,:)
      integer, intent(in) :: d
      integer :: i,j,k,l,p,rW,rF,M,N,wind
      integer :: bcola,bcol,lcola,lcola2,lcol
      real*8  :: tmp

      rW=W%R()
      rF=SIZE(X%B,1) ! Also equals cols-of-W and rows-of-F
      N=W%rows(d)
      M=F%cols(d)

!     Error checking
      if ((SIZE(LHS,1).ne.rF).or.(SIZE(LHS,2).ne.(rF*M*N))) then
         write(*,*) 'LHS is [',SIZE(LHS,1),' x ',SIZE(LHS,2),&
                 '] but must be [',rF,' x ',rF,'*',N,']'
         call AbortWithError('AccumulateLHSforALS(): bad dimensions')
      endif

!     Compute the left-hand sides of the linear system
      do k=1,rW
         bcola=(k-1)*rF ! col group of B from W
         do l=1,N
            wind=W%ibas(d)+(l-1)*(N+1)
            tmp=W%coef(k)*W%base(wind,k)
            do p=1,M
               lcola=((l-1)+(p-1)*N)*rF ! col group of lhs from basis
               do i=1,rF
                  do j=1,rF
                     bcol=bcola+j ! col of B from W and F
                     lcol=lcola+j ! col of lhs from W and F
                     LHS(i,lcol)=LHS(i,lcol)+tmp*X%B(i,bcol)
                  enddo
               enddo
            enddo
         enddo
      enddo

      end subroutine AccumulateLHSforALS 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRHSforLinSys(X,d) result(RHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets right hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, for general A.
! Pass W*G to this routine instead of G for weighted solve.

      implicit none
      TYPE (ALS), intent(in) :: X
      real*8, allocatable :: RHS(:,:)
      integer, intent(in) :: d
      integer :: rA,rF,rWG,M,N,amode,gmode

!     Pick the correct mode to index A^T*A
      call GetAGmode(X,d,amode,gmode)

      IF (.not.allocated(X%ATG%coef)) &
         call AbortWithError('GetRHSforLinSys(): no ATG allocated')

!     Dimensions: ATG has rG*rA terms; P is [rG x rArF] 
      rWG=SIZE(X%P,1)
      rA=SIZE(X%ATG%coef)/rWG
      rF=SIZE(X%P,2)/rA
      M=X%ATG%cols(gmode)
      N=X%ATG%rows(gmode)

      allocate(RHS(rF*N,M))
      RHS=0.d0

      call AccumulateRHSforLinSys(X,d,RHS)

      end function GetRHSforLinSys

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateRHSforLinSys(X,d,RHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates right hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, for general A.
! Pass W*G to this routine instead of G for weighted solve.

      implicit none
      TYPE (ALS), intent(in) :: X
      real*8, intent(inout)  :: RHS(:,:)
      integer, intent(in) :: d
      integer :: j,k,l,s,t,rA,rF,rWG,M,N
      integer :: trk,pcola,rras,rraf,amode,gmode

!     Pick the correct mode to index A^T*A
      call GetAGmode(X,d,amode,gmode)

      IF (.not.allocated(X%ATG%coef)) &
         call AbortWithError('AccumulateRHSforLinSys(): no ATG allocd')

!     Dimensions: ATG has rG*rA terms; P is [rG x rArF] 
      rWG=SIZE(X%P,1)
      rA=SIZE(X%ATG%coef)/rWG
      rF=SIZE(X%P,2)/rA
      M=X%ATG%cols(gmode)
      N=X%ATG%rows(gmode)

!     Error checking
      if ((SIZE(RHS,1).ne.(rF*N)).or.(SIZE(RHS,2).ne.M)) then
         write(*,*) 'RHS is [',SIZE(RHS,1),' x ',SIZE(RHS,2),&
                 '] but must be [',rF,'*',N,' x ',M,']'
         call AbortWithError('AccumulateRHSforLinSys(): bad dimensions')
      endif

!     Compute the right-hand sides of the linear system
      do s=0,rA-1
         pcola=s*rF ! col group of P from A
         do t=1,rWG
            trk=s*rWG+t ! rank of ATG from AT and G
            l=X%ATG%ibas(gmode)
            do j=1,M
               rras=1
               rraf=rF
               do k=1,N
                  RHS(rras:rraf,j)=RHS(rras:rraf,j)+&
                  X%ATG%coef(trk)*X%ATG%base(l,trk)*&
                  X%P(t,pcola+1:pcola+rF)
                  l=l+1
                  rras=rras+rF
                  rraf=rraf+rF
               enddo
            enddo
         enddo
      enddo

      end subroutine AccumulateRHSforLinSys

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRHSforALS(X,G,d) result(RHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets right hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, optimized for A = I
! Pass W*G to this routine instead of G for weighted ALS.

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in)  :: G
      real*8, allocatable :: RHS(:,:)
      integer, intent(in) :: d
      integer :: rF,N

      rF=SIZE(X%P,2)
      N=G%nbas(d)

      allocate(RHS(rF,N))
      RHS=0.d0

      call AccumulateRHSforALS(X,G,d,RHS)

      end function GetRHSforALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateRHSforALS(X,G,d,RHS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates right hand side of CP-format linear system normal equations 
! A^T*A*F = A^T*G, optimized for A = I.
! Pass W*G to this routine instead of G for weighted ALS.

      implicit none
      TYPE (ALS), intent(in) :: X
      TYPE (CP), intent(in)  :: G
      real*8, intent(inout)  :: RHS(:,:)
      integer, intent(in) :: d
      integer :: j,k,l,rF,rG,N

      rF=SIZE(X%P,2)
      rG=SIZE(X%P,1)
      N=G%nbas(d)

!     Error checking
      if ((SIZE(RHS,1).ne.rF).or.(SIZE(RHS,2).ne.N)) then
         write(*,*) 'RHS is [',SIZE(RHS,1),' x ',SIZE(RHS,2),&
                 '] but must be [',rF,' x ',N,']'
         call AbortWithError('AccumulateRHSforALS(): bad dimensions')
      endif

!     Compute the right-hand sides of the linear system
      l=G%ibas(d)
      do k=1,N
         do j=1,rG
            RHS(1:rF,k)=RHS(1:rF,k)+&
            G%coef(j)*G%base(l,j)*X%P(j,1:rF)
         enddo
         l=l+1
      enddo

      end subroutine AccumulateRHSforALS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLS_noweights(X,F,G,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system for ALS object X, for mode d, given CP vectors F
! and G. The values for mode d of F, are replaced by the solution.

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F
      TYPE (CP), intent(in) :: G
      integer, intent(in) :: d
      real*8, allocatable :: LHS(:,:),RHS(:,:)
      integer :: i

      IF (.not.allocated(X%B)) &
         call AbortWithError('SolveLS(): B is not allocated')
      IF (.not.allocated(X%P)) &
         call AbortWithError('SolveLS(): P is not allocated')
      IF (.not.CheckModeIncludedState(X,d)) THEN
         do i=1,SIZE(X%dofincluded)
            if (i.eq.d) then
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .FALSE.)'
            else
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .TRUE.)'
            endif
         enddo
         call AbortWithError('SolveLS(): wrong dof-included state')
      ENDIF

!     Calculate LHS and RHS of linear system
      IF (X%AisI) THEN
         LHS=GetLHSforALS(X,d)
         RHS=GetRHSforALS(X,G,d)
      ELSE
         LHS=GetLHSforLinSys(X,F,d)
         RHS=GetRHSforLinSys(X,d)
         call X%FlushNormalEquations
      ENDIF

!     Solve linear system
      IF (X%useSVD) THEN
         call SolveLinSysSVD(LHS,RHS,X%penalty)
      ELSE
         call SolveLinSysLU(LHS,RHS,X%penalty)
      ENDIF
      deallocate(LHS)

!     Put solution into F
      call UpdateFfromSoln(F,RHS,d)
      deallocate(RHS)

!     Check coefs of F for NaN values resulting from zero division
      IF (.NOT. CHECKCOEFS(F)) THEN
         write(*,*) 'SolveLS(): NaN on update: mode ',d
         call AbortWithError('SolveLS() crashed')
      ENDIF

      end subroutine SolveLS_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLS_weights(X,F,WG,W,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves linear system for ALS object X, for mode d, given CP vectors F
! and G. The values for mode d of F, are replaced by the solution.
! Pass W*G to this routine instead of G for weighted ALS or solve.

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(inout)  :: F
      TYPE (CP), intent(in) :: WG,W
      integer, intent(in) :: d
      real*8, allocatable :: LHS(:,:),RHS(:,:)
      integer :: i

      IF (.not.allocated(X%B)) &
         call AbortWithError('SolveLS(): B is not allocated')
      IF (.not.allocated(X%P)) &
         call AbortWithError('SolveLS(): P is not allocated')
      IF (.not.CheckModeIncludedState(X,d)) THEN
         do i=1,SIZE(X%dofincluded)
            if (i.eq.d) then
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .FALSE.)'
            else
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .TRUE.)'
            endif
         enddo
         call AbortWithError('SolveLS(): wrong dof-included state')
      ENDIF

!     Calculate LHS and RHS of linear system
      IF (X%AisI) THEN
         LHS=GetLHSforALS(X,F,W,d)
         RHS=GetRHSforALS(X,WG,d)
         call SolveWeightedLS(LHS,RHS,X%useSVD,X%penalty)
      ELSE
         LHS=GetLHSforLinSys(X,F,d)
         RHS=GetRHSforLinSys(X,d)
         call X%FlushNormalEquations

!        Solve linear system
         IF (X%useSVD) THEN
            call SolveLinSysSVD(LHS,RHS,X%penalty)
         ELSE
            call SolveLinSysLU(LHS,RHS,X%penalty)
         ENDIF
      ENDIF

      deallocate(LHS)

!     Put solution into F
      call UpdateFfromSoln(F,RHS,d)
      deallocate(RHS)

!     Check coefs of F for NaN values resulting from zero division
      IF (.NOT. CHECKCOEFS(F)) THEN
         write(*,*) 'SolveLS(): NaN on update: mode ',d
         call AbortWithError('SolveLS() crashed')
      ENDIF

      end subroutine SolveLS_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveWeightedLS(LHS,RHS,useSVD,valpen)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Solves N rhs linear system (ALS case with weights) 

      implicit none
      real*8, intent(inout) :: LHS(:,:),RHS(:,:)
      real*8, allocatable :: lhs1(:,:),rhs1(:,:)
      logical, intent(in) :: useSVD
      real*8, intent(in)  :: valpen
      integer :: i,j,rF,N,cst,cfi

!     LHS is rF x N*rF; RHS is rF x N
      rF=SIZE(LHS,1)
      N=SIZE(LHS,2)/rF

      ALLOCATE(lhs1(rF,rF),rhs1(rF,1))

      cst=1
      cfi=rF
      DO i=1,N
!        Solve a separate linear equation for each i
         lhs1(1:rF,1:rF)=LHS(1:rF,cst:cfi)
         rhs1(1:rF,1:1)=RHS(1:rF,i:i)

         IF (useSVD) THEN
            call SolveLinSysSVD(lhs1,rhs1,valpen)
         ELSE
            call SolveLinSysLU(lhs1,rhs1,valpen)
         ENDIF

!        Copy the result back to RHS
         RHS(1:rF,i:i)=rhs1(1:rF,1:1)

         cst=cst+rF
         cfi=cfi+rF
      ENDDO
      DEALLOCATE(lhs1,rhs1)

      end subroutine SolveWeightedLS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsDD_noweights(X,F,G,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! For linear systems call using A*F instead of F

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,G
      integer, intent(in)   :: d

!     Call with all modes included: case for downdating
      if (ALL(X%dofincluded(:))) then
         call UPDATEP(F,d,X%B,.TRUE.)
         call UPDATEP(F,G,d,X%P,.TRUE.)
         X%dofincluded(d)=.FALSE.

!     Call with no modes included: case for building or first-call-only
!     for downdating
      elseif (ALL(.NOT.X%dofincluded(:))) then
         call CONSTPT(F,d,X%B)
         call CONSTPT(F,G,d,X%P)
         X%dofincluded(:)=.TRUE.
         X%dofincluded(d)=.FALSE.

!     Unexpected state (error out)
      else
         write(*,*) 'When building B,P matrices, must call this ',&
                    'routine with either ALL or NO modes included.'
         call AbortWithError('ProdMatsDD(): unexpected B,P state')
      endif

      end subroutine ProdMatsDD_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsDD_weights(X,F,WF,G,WG,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! For linear systems call using A*F instead of F
! Call with W*G instead of G and W*F (W*A*F) for ALS (linear solver)

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,WF,G,WG
      integer, intent(in)   :: d

!     Call with all modes included: case for downdating
      if (ALL(X%dofincluded(:))) then
         call UPDATEP(WF,F,d,X%B,.TRUE.)
         call UPDATEP(F,WG,d,X%P,.TRUE.)
         X%dofincluded(d)=.FALSE.

!     Call with no modes included: case for building or first-call-only
!     for downdating
      elseif (ALL(.NOT.X%dofincluded(:))) then
         call CONSTPT(WF,F,d,X%B)
         call CONSTPT(F,WG,d,X%P)
         X%dofincluded(:)=.TRUE.
         X%dofincluded(d)=.FALSE.

!     Unexpected state (error out)
      else
         write(*,*) 'When building B,P matrices, must call this ',&
                    'routine with either ALL or NO modes included.'
         call AbortWithError('ProdMatsDD(): unexpected B,P state')
      endif

      end subroutine ProdMatsDD_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsUD_noweights(X,F,G,d,conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! For linear systems call using A*F instead of F

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,G
      integer, intent(out)  :: conv
      integer, intent(in)   :: d
      integer :: i
      logical :: err

      conv=0

!     Make sure B, P are in expected state
!     (all modes except d included)
      err=.FALSE.
      do i=1,SIZE(X%dofincluded)
         if (i.eq.d) cycle
         if (.not.X%dofincluded(i)) err=.TRUE.
      enddo
      if (X%dofincluded(d)) err=.TRUE.
      if (err) then
         write(*,*) 'All modes must be included except mode d'
         call AbortWithError('ProdMatsUD(): unexpected B,P state')
      endif

!     Update or scrub the B, P matrices
      if (X%update) then
         call UPDATEP(F,d,X%B,.FALSE.)
         call UPDATEP(F,G,d,X%P,.FALSE.)
         X%dofincluded(d)=.TRUE.
         if (X%chkconv) conv=CheckALSConvergence(X,F,G)
      else
!        If convergence check is requested, update B, P before reseting
         if (X%chkconv) then
            call UPDATEP(F,d,X%B,.FALSE.)
            call UPDATEP(F,G,d,X%P,.FALSE.)
            X%dofincluded(d)=.TRUE.
            conv=CheckALSConvergence(X,F,G)
         endif
!        Reset all modes to "not included" but do not deallocate here
         X%dofincluded(:)=.FALSE.
      endif

      end subroutine ProdMatsUD_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsUD_weights(X,F,WF,G,WG,d,conv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! For linear systems call using A*F instead of F
! Call with W*G instead of G and W*F (W*A*F) for ALS (linear solver)

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,WF,G,WG
      integer, intent(out)  :: conv
      integer, intent(in)   :: d
      integer :: i
      logical :: err

      conv=0

!     Make sure B, P are in expected state
!     (all modes except d included)
      err=.FALSE.
      do i=1,SIZE(X%dofincluded)
         if (i.eq.d) cycle
         if (.not.X%dofincluded(i)) err=.TRUE.
      enddo
      if (X%dofincluded(d)) err=.TRUE.
      if (err) then
         write(*,*) 'All modes must be included except mode d'
         call AbortWithError('ProdMatsDD(): unexpected B,P state')
      endif

!     Update or scrub the B, P matrices
      if (X%update) then
         call UPDATEP(WF,F,d,X%B,.FALSE.)
         call UPDATEP(F,WG,d,X%P,.FALSE.)
         X%dofincluded(d)=.TRUE.
         if (X%chkconv) conv=CheckALSConvergence(X,F,WF,G,WG)
      else
!        If convergence check is requested, update B, P before reseting
         if (X%chkconv) then
            call UPDATEP(WF,F,d,X%B,.FALSE.)
            call UPDATEP(F,WG,d,X%P,.FALSE.)
            X%dofincluded(d)=.TRUE.
            conv=CheckALSConvergence(X,F,WF,G,WG)
         endif
!        Reset all modes to "not included" but do not deallocate here
         X%dofincluded(:)=.FALSE.
      endif

      end subroutine ProdMatsUD_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckALSConvergence_noweights(X,F,G) result(converged)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||F-G|| or ||AF-G|| as a convergence test
! Call using A*F instead of F when using linear solver instead of ALS

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,G
      integer :: i,rF,rG,converged
      real*8  :: conver

      rF=F%R()
      rG=G%R()

!     Compute <G,G> if not done already
      IF (X%GG.le.0.d0) X%GG=PRODVV(G)

!     <F,F> or <AF,AF> is in X%B
      X%FF=0.d0
      X%CN=0.d0
      DO i=1,rF
         X%FF=X%FF+F%coef(i)*dot_product(F%coef(1:rF),X%B(1:rF,i))
         X%CN=X%CN+F%coef(i)*F%coef(i)
      ENDDO

!     <F,G> or <AF,G> is in X%P
      X%FG=0.d0
      DO i=1,rF
         X%FG=X%FG+F%coef(i)*dot_product(G%coef(1:rG),X%P(1:rG,i))
      ENDDO

!     Compute result and check convergence
      converged=0
      conver=sqrt(abs((X%FF+X%GG-2*X%FG)/X%GG))
      X%del=abs((conver-X%conver)/conver)
      IF (X%del.lt.X%dthresh) THEN 
         converged=2
      ENDIF
      IF (conver.lt.X%thresh) THEN
         converged=1
      ENDIF
      X%conver=conver

      IF (X%prtconv) call PrintALSConvergence(X)

      end function CheckALSConvergence_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckALSConvergence_weights(X,F,WF,G,WG) &
               result(converged)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||F-G|| or ||AF-G|| as a convergence test
! Call using A*F instead of F when using linear solver instead of ALS
! G is actually WG in this call

      implicit none
      TYPE (ALS), intent(inout) :: X
      TYPE (CP), intent(in) :: F,WF,G,WG
      integer :: i,rF,rWF,rG,rWG,converged
      real*8  :: conver

      rF=F%R()
      rG=G%R()
      rWF=WF%R()
      rWG=WG%R()

!     Compute <WG,WG> if not done already
      IF (X%GG.le.0.d0) X%GG=PRODVV(G,WG)

!     <F,WF> or <AF,WAF> is in X%B
      X%FF=0.d0
      DO i=1,rWF
         X%FF=X%FF+WF%coef(i)*dot_product(F%coef(1:rF),X%B(1:rF,i))
      ENDDO

      X%CN=0.d0
      DO i=1,rF
         X%CN=X%CN+F%coef(i)*F%coef(i)
      ENDDO

!     <F,WG> or <AF,WG> is in X%P
      X%FG=0.d0
      DO i=1,rF
         X%FG=X%FG+F%coef(i)*dot_product(WG%coef(1:rWG),X%P(1:rWG,i))
      ENDDO

!     Compute result and check convergence. The more important
!     condition is error < thresh, so always set convergence to 1 if
!     this is achieved
      converged=0
      conver=sqrt(abs((X%FF+X%GG-2*X%FG)/X%GG))
      X%del=abs((conver-X%conver)/conver)
      IF (X%del.lt.X%dthresh) THEN
         converged=2
      ENDIF
      IF (conver.lt.X%thresh) THEN
         converged=1
      ENDIF
      X%conver=conver

      IF (X%prtconv) call PrintALSConvergence(X)

      end function CheckALSConvergence_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintALSConvergence(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out ALS or linear solver convergence information in ALS object

      implicit none
      CLASS (ALS), intent(in) :: X
      character(len=64) :: frmt,tag,ctag,wtag

      IF (.not.X%chkconv) THEN
         write(*,*) 'X%chkconv is .F.; must be .T. to print convergence'
         call AbortWithError('PrintALSConvergence(): no conv check')
      ENDIF

!     Tag for name of object
      IF (LEN(TRIM(ADJUSTL(X%alsnm))).gt.0) THEN
         write(tag,'(2A)') TRIM(ADJUSTL(X%alsnm)),':'
      ELSE
         write(tag,'(A)') ''
      ENDIF

!     Tag for convergence
      ctag=''
      IF (X%del.lt.X%dthresh) THEN
         ctag='delta < dtol'
      ENDIF
      IF (X%conver.lt.X%thresh) THEN
         ctag='conver < ctol!'
      ENDIF

!     Tag for weights
      IF (X%weights) THEN
         wtag='(wtd)'
      ELSE
         wtag=''
      ENDIF

      IF (X%AisI) THEN
         frmt='(A,X,A,3(A,ES14.8),X,A)'
         write(*,frmt) &
         TRIM(ADJUSTL(tag)),TRIM(ADJUSTL(wtag)),' ||F-G||/||G|| = ',&
         X%conver,'; ||F|| = ',sqrt(abs(X%FF)),'; CN = ',&
         sqrt(abs(X%CN/X%FF)),TRIM(ADJUSTL(ctag))
      ELSE
         frmt='(A,X,A,A,ES14.8,X,A)'
!        Linear solver passes A*F, not F, so must compute here
         write(*,frmt) TRIM(ADJUSTL(tag)),TRIM(ADJUSTL(wtag)),&
         ' ||A*F-G||/||G|| = ',X%conver,TRIM(ADJUSTL(ctag))
      ENDIF

      end subroutine PrintALSConvergence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ALSOO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
