!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE ALSOO8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains object-oriented code for ALS

      USE ERRORTRAP
      USE UTILS
      USE CPr8
      USE MYACC
      USE LinAlg8
      USE CPVV8
      USE CPMM8
      USE CPLS8

      TYPE ALS8
         TYPE (LS8) :: LS
         TYPE (CP8) :: ATWA,ATWG
         TYPE (VV8) :: VB,VP,VB1,VB2
         real(kind=8), allocatable :: B(:),P(:),B1(:),B2(:),LHS(:),RHS(:)
         LOGICAL, ALLOCATABLE :: dofincluded(:)
         real(kind=8)  :: AFAF,FAF,FF,FG,GG,CN,RQ,conver,del
         integer :: accum_HS_ct,accum_conv_ct
         CHARACTER(LEN=64) :: alsnm
         CHARACTER(LEN=64) :: solver
         ! Options
         integer      :: algo    ! 0 for CPU, 1 for OpenACC
         real(kind=8) :: penalty ! penalty on LHS to prevent ill-conditioning
         real(kind=8) :: thresh  ! fractional tolerance for convergence
         real(kind=8) :: dthresh ! fractional derivative tolerance
         LOGICAL :: normal  ! .T. uses normal form in linear solver
         LOGICAL :: AisI    ! .T. for ALS, .F. for linear solver
         LOGICAL :: WisI    ! .T. when not using weights
         LOGICAL :: update  ! .T. up/downdates B,P; .F. always builds
         LOGICAL :: chkconv ! .T. to check convergence after each update
         LOGICAL :: prtconv ! .T. to print convergence after each update
         LOGICAL :: calcGG  ! .T. to calculate <G,G> in convergence test
         LOGICAL :: Achange ! .T. to update A^T*A, A^T*G each call
         LOGICAL :: Gchange ! .T. to update A^T*G or G each call
         LOGICAL :: Wchange ! .T. to update A^T*W*A each call
         CONTAINS
            PROCEDURE :: new => NewALS_CP8
            PROCEDURE :: flush => FlushALS_CP8
            PROCEDURE :: setdefaults => ALSOptionDefaults_CP8
            PROCEDURE :: setname => ALSSetName_CP8
            PROCEDURE :: setlogical => ALSSetOption_logical_CP8
            PROCEDURE :: setrealval => ALSSetOption_real_CP8
            PROCEDURE :: setstringval => ALSSetOption_string_CP8
            PROCEDURE :: setoptions => ALSSetOptions_CP8
            PROCEDURE :: show => ShowALSParameters_CP8
            PROCEDURE :: showprodmats => ShowProductMatrices_CP8
            PROCEDURE :: showls => ShowLS_CP8
            PROCEDURE :: resetdofstate => ResetDOFIncludedState_CP8
            PROCEDURE :: accumprods => AccumulateProducts_CP8
            PROCEDURE :: recalcnormal => RecomputeNormalEquations_CP8
            PROCEDURE :: flushnormal => FlushNormalEquations_CP8
            PROCEDURE :: checkmodestate => CheckModeIncludedState_CP8
            PROCEDURE :: constls => ConstructLS_CP8
            PROCEDURE :: solvels => SolveLS_CP8
            PROCEDURE :: downdateBP => ProdMatsDD_CP8
            PROCEDURE :: updateBP => ProdMatsUD_CP8
            PROCEDURE :: checkconv => CheckALSConvergence_CP8
            PROCEDURE :: showconv => PrintALSConvergence_CP8
            PROCEDURE :: resetconv => ResetALSConvergence_CP8 
      END TYPE ALS8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewALS_CP8(X,A,ish,Esh,F,G,W,Achange,Gchange,Wchange,&
                            algo,usenormal)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the ALS object

      implicit none
      CLASS (ALS8) :: X
      TYPE (CP8), intent(in) :: A,F,G,W
      TYPE (CP8) :: Wmult
      logical, intent(in) :: Achange,Gchange,Wchange
      logical, intent(in), optional :: usenormal
      integer, intent(in) :: ish,algo
      real(kind=8), intent(in)  :: Esh
      integer :: ndof,rA,rAp,rF,rG,rW
      real(kind=8) :: sgn

      IF (present(usenormal)) THEN
         X%normal=usenormal
      ELSE
         X%normal=.TRUE.
      ENDIF

!     Initialize the derived types
      call X%LS%new(A,F,G,W,ish,Esh,algo,X%normal) ! LS type

!     Additional error checking (most is done in X%LS%new())
      IF (Achange .and. (.not.X%LS%useA)) THEN
         write(*,*) 'Achange must be .FALSE. for ALS'
         call AbortWithError('Error in NewALS_CP8()')
      ENDIF

      IF (Wchange .and. (.not.X%LS%useW)) THEN
         write(*,*) 'Wchange must be .FALSE. when not using weights'
         call AbortWithError('Error in NewALS_CP8()')
      ENDIF
      
      X%algo=algo
      X%AisI=(.not.X%LS%useA)
      X%WisI=(.not.X%LS%useW)
      X%Achange=Achange
      X%Gchange=Gchange
      X%Wchange=Wchange
      X%accum_HS_ct=0
      X%accum_conv_ct=0

      ndof=X%LS%ndof
      rA=X%LS%rksA
      rAp=X%LS%rksAp
      rF=X%LS%rkF
      rG=X%LS%rkG
      rW=X%LS%rkW

      allocate(X%B(rW*rA*rF*rAp*rF),X%P(rW*rG*rAp*rF),X%dofincluded(ndof))
      IF (.not.X%AisI) THEN
         allocate(X%B1(rAp*rF*rF))
         IF (X%normal) THEN
            allocate(X%B2(rF*rF))
         ENDIF
      ENDIF

      X%dofincluded(:)=.FALSE.

      if (algo.eq.1) then
#if ACC_ENABLED
         !$acc enter data create(X%B,X%P)
         if (.not.X%AisI) then
            !$acc enter data create(X%B1)
            if (X%normal) then
               !$acc enter data create(X%B2)
            endif
         endif
#else
         write(*,*) 'OpenACC algorithm not enabled'
         call AbortWithError("NewALS_CP8(): wrong algorithm")
#endif
      endif

!     Calc A^T*W*A if A,W do not change
      IF (X%LS%useA .and. X%LS%useW .and. &
         (.not.Achange) .and. (.not.Wchange)) THEN ! wLS
         IF (X%normal) THEN
            call CPMM_CP8(W,0,0.0,.FALSE.,A,ish,Esh,.FALSE.,Wmult,.FALSE.,algo)
            call CPMM_CP8(A,ish,Esh,.TRUE.,Wmult,0,0.0,.FALSE.,X%ATWA,.TRUE.,algo)
            call Wmult%flush()
         ELSE
            call CPMM_CP8(W,0,0.0,.FALSE.,A,ish,Esh,.FALSE.,X%ATWA,.FALSE.,algo)
         ENDIF
         call PrepG_CP8(X%ATWA,algo)
      ELSEIF (X%LS%useA .and. (.not.Achange)) THEN ! LS
         IF (X%normal) THEN
            call CPMM_CP8(A,ish,Esh,.TRUE.,A,ish,Esh,.FALSE.,X%ATWA,.TRUE.,algo)
         ELSE
             call Wmult%identity(A%cols,A%cols)
             if (algo.eq.1) call Wmult%copyintodevice()
             call CPMM_CP8(A,ish,Esh,.FALSE.,Wmult,0,0.0,.FALSE.,X%ATWA,.TRUE.,algo)
             call Wmult%flush()
         ENDIF
         call PrepG_CP8(X%ATWA,algo)
      ELSEIF (X%LS%useW .and. (.not.Wchange)) THEN ! wALS
         call X%ATWA%copyfrom(W)
         call PrepG_CP8(X%ATWA,algo)
      ENDIF

!     Calc A^T*W*G if A,G,W do not change
      IF (.not.Gchange) THEN
         IF (X%LS%useA .and. X%LS%useW) THEN ! wLS
            IF ((.not.Achange) .and. (.not.Wchange)) THEN
               IF (X%normal) THEN
                  call CPMM_CP8(W,.FALSE.,G,.FALSE.,Wmult,.TRUE.,algo)
                  call CPMM_CP8(A,ish,Esh,.TRUE.,Wmult,0,0.0,.FALSE.,X%ATWG,.TRUE.,algo)
                  call Wmult%flush()
               ELSE
                  call CPMM_CP8(W,.FALSE.,G,.FALSE.,X%ATWG,.TRUE.,algo)
               ENDIF
               call PrepG_CP8(X%ATWG,algo)
            ENDIF
         ELSEIF (X%LS%useA) THEN ! LS
            IF (.not.Achange) THEN 
               IF (X%normal) THEN
                  call CPMM_CP8(A,ish,Esh,.TRUE.,G,0,0.0,.FALSE.,X%ATWG,.TRUE.,algo)
               ELSE
                  call X%ATWG%copyfrom(G)
               ENDIF
               call PrepG_CP8(X%ATWG,algo)
            ENDIF
         ELSEIF (X%LS%useW) THEN ! wALS
            IF (.not.Wchange) THEN
               call CPMM_CP8(W,.FALSE.,G,.FALSE.,X%ATWG,.TRUE.,algo)
               call PrepG_CP8(X%ATWG,algo)
            ENDIF
         ELSE ! ALS
            call X%ATWG%copyfrom(G)
            call PrepG_CP8(X%ATWG,algo)
         ENDIF
      ENDIF

      call X%setdefaults(ndof)

      end subroutine NewALS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushALS_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Frees memory from ALS object

      implicit none
      CLASS (ALS8), intent(inout) :: X

      X%alsnm=''
      call X%ATWA%flush()
      call X%ATWG%flush()
      call X%VB%flush()
      call X%VP%flush()
      if ((.not.X%AisI).and.X%chkconv) then
         call X%VB1%flush()
         if (X%normal) then
            call X%VB2%flush()
         endif
      endif
      call X%LS%flush()
      if (X%algo.eq.1) then
#if ACC_ENABLED
         !$acc exit data delete(X%B,X%P)
         if (.not.X%AisI) then
            !$acc exit data delete(X%B1)
            if (X%normal) then
               !$acc exit data delete(X%B2)
            endif
         endif
#endif
      endif
      IF (ALLOCATED(X%P)) DEALLOCATE(X%P)
      IF (ALLOCATED(X%B)) DEALLOCATE(X%B)
      IF (ALLOCATED(X%B1)) DEALLOCATE(X%B1)
      IF (ALLOCATED(X%B2)) DEALLOCATE(X%B2)
      IF (ALLOCATED(X%dofincluded)) DEALLOCATE(X%dofincluded)

      end subroutine FlushALS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSOptionDefaults_CP8(X,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets default values for the parameters for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout) :: X
      integer, intent(in) :: ndof
      logical :: update,chkconv,calcGG

      X%alsnm=''
      X%GG=-1.d0
      X%conver=1.d99
      X%thresh=8.d-8
      X%dthresh=8.d-7

!     Solver parameters should be set by values in input file
      X%solver='uninitialized'
      X%penalty=-1.d0

!     Updating risks zero division, but is cheaper for ndof > 3
      update=(ndof.gt.3)
      chkconv=.FALSE.
      calcGG=.TRUE.

      call X%setoptions(update,chkconv,calcGG)

      end subroutine ALSOptionDefaults_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetName_CP8(X,nm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout) :: X
      character(len=*), intent(in) :: nm
      
      X%alsnm=TRIM(ADJUSTL(nm))

      end subroutine ALSSetName_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOption_logical_CP8(X,opt,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout) :: X
      character(len=*), intent(in) :: opt
      logical, intent(in) :: val

      IF (TRIM(ADJUSTL(opt)).eq.'update') THEN
         X%update=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'chkconv') THEN
         X%chkconv=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'calcGG') THEN
         X%calcGG=val
      ELSEIF (TRIM(ADJUSTL(opt)).eq.'prtconv') THEN
         X%prtconv=val
      ELSE
         write(*,*) "Unrecognized logical option: '",opt,"'"
         write(*,*) "Choose from 'update', 'usesvd', "&
                    "'chkconv','prtconv'"
         call AbortWithError('ALSSetOption(): unrecognized option')
      ENDIF

      end subroutine ALSSetOption_logical_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOption_real_CP8(X,opt,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout) :: X
      character(len=*), intent(in) :: opt
      real(kind=8), intent(in) :: val

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

      end subroutine ALSSetOption_real_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOption_string_CP8(X,opt,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout)  :: X
      character(len=*), intent(in) :: opt,val

      IF (TRIM(ADJUSTL(opt)).seq.'solver') THEN
         X%solver=val
      ELSE
         write(*,*) "Unrecognized string option: '",&
                    TRIM(ADJUSTL(opt)),"'"
         write(*,*) "Choose from 'solver' "
         call AbortWithError('ALSSetOption(): unrecognized option')
      ENDIF

      end subroutine ALSSetOption_string_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALSSetOptions_CP8(X,update,chkconv,calcGG)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets options for TYPE ALS

      implicit none
      CLASS (ALS8), intent(inout) :: X
      logical, intent(in) :: update,chkconv,calcGG

      X%update=update
      X%chkconv=chkconv
      X%calcGG=calcGG
      X%prtconv=chkconv ! Init same as chkconv

      end subroutine ALSSetOptions_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowALSParameters_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints ALS parameters

      implicit none
      CLASS (ALS8), intent(in) :: X
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
      write(*,*) 'algo    :',X%algo,&
      ' (0 for CPU, including OpenMP, 1 for GPU OpenACC algorithm)'
      write(*,*) 'normal  :',X%normal,&
      ' (.T. solves normal form of linear eqns, .F. uses faster algo)'
      write(*,*) 'update  :',X%update,&
      ' (.T. downdates/updates P,B matrices, .F. builds each iteration)'
      write(*,*) 'solver  :',TRIM(ADJUSTL(X%solver)),&
      "('SVD' ,Moore-Penrose pseudoinverse); or 'LU', LU decomposition)"
      write(*,*) 'chkconv :',X%chkconv,&
      ' (.T. to check LinSolv,ALS convergence, .F. skips)'
      write(*,*) 'calcGG :',X%calcGG,&
      ' (.T. to calculate <G,G> in convergence test, .F. skips)'
      write(*,*) 'prtconv :',X%prtconv,&
      ' (.T. to print LinSolv,ALS convergence, .F. skips)'
      write(*,*) 'Achange :',X%Achange,&
      ' (.T. if A changes each iteration, .F. for constant A)'
      write(*,*) 'Gchange :',X%Gchange,&
      ' (.T. if G changes each iteration, .F. for constant G)'
      write(*,*) 'Wchange :',X%Wchange,&
      ' (.T. if W changes each iteration, .F. for constant W)'
      write(*,*) 'penalty :',X%penalty,&
      ' (regularization value to prevent ill-conditioning)'
      if (X%chkconv) then
         write(*,*) 'thresh  :',X%thresh,&
         ' (convergence threshold for ||F-G||/||G||)'
         write(*,*) 'dthresh :',X%dthresh,&
         ' (convergence threshold for abs((conver-X%conver)/conver))'
      endif
      write(*,*)
      write(*,*) 'Arrays:'
      IF (ALLOCATED(X%ATWA%coef)) THEN
         write(*,*) 'A^T*W*A : ALLOCATED, with rank = ',X%ATWA%R()
      ELSE
         write(*,*) 'A^T*W*A : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%ATWG%coef)) THEN
         write(*,*) 'A^T*W*G : ALLOCATED, with rank = ',X%ATWG%R()
      ELSE
         write(*,*) 'A^T*W*G : NOT ALLOCATED'
      ENDIF
      frmt='(X,A,I0,A)'
      IF (ALLOCATED(X%B)) THEN
         write(*,frmt) '    B : [',SIZE(X%B),']'
      ELSE
         write(*,*) '    B : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%P)) THEN
         write(*,frmt) '    P : [',SIZE(X%P),']'
      ELSE
         write(*,*) '    P : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%B1)) THEN
         write(*,frmt) '    B1 : [',SIZE(X%B1),']'
      ELSE
         write(*,*) '    B1 : NOT ALLOCATED'
      ENDIF
      IF (ALLOCATED(X%B2)) THEN
         write(*,frmt) '    B2 : [',SIZE(X%B2),']'
      ELSE
         write(*,*) '    B2 : NOT ALLOCATED'
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

      end subroutine ShowALSParameters_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ResetDOFIncludedState_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resets DOF-included state to all .FALSE. for accumulated A version

      implicit none
      CLASS (ALS8), intent(inout) :: X

      X%dofincluded(:)=.FALSE.

      end subroutine ResetDOFIncludedState_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AccumulateProducts_CP8(X,WAF,AF,F,WG)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs B, P matrices when no modes are included (for convergence
! pre-check)
! Call using: F, F for B;  G, F for P ( ALS)
!            WF, F for B; WG, F for P (wALS)
!            AF, F for B;  G, F for P (  LS)
!           WAF, F for B; WG, F for P ( wLS)
!            AF,AF for B;  G,AF for P (  LS, normal form)
!           WAF,AF for B; WG,AF for P ( wLS, normal form)

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(in) :: WAF,AF,F,WG

      if ((.not.X%AisI) .and. X%normal .and. X%accum_conv_ct.gt.0) then
         write(*,*) 'LS accumulation not implemented for normal form'
         call AbortWithError('Error in AccumulateProducts_CP8()')
      endif

!     Initialize VV types on first call
!     Construct inner product arrays from scratch on first call or if
!     building in each ALS microiteration
      if (ALL(.NOT.X%dofincluded(:))) then

!        Arrays used in ALS:
         if (.not.X%VB%setup) call X%VB%new(WAF,AF,X%algo)
         call CONST_PT_CP8(X%VB,WAF,AF,0,X%B)
         if (.not.X%VP%setup) call X%VP%new(WG,AF,X%algo)
         call CONST_PT_CP8(X%VP,WG,AF,0,X%P)

!        Arrays used to check convergence:
         if ((.not.X%AisI).and.X%chkconv) then
            if (.not.X%VB1%setup) call X%VB1%new(AF,F,X%algo)
            call CONST_PT_CP8(X%VB1,AF,F,0,X%B1)
            if (X%normal) then
               if (.not.X%VB2%setup) call X%VB2%new(F,F,X%algo)
               call CONST_PT_CP8(X%VB2,F,F,0,X%B2)
            endif
         endif

         X%dofincluded(:)=.TRUE.

!     Already in correct state
      elseif (ALL(X%dofincluded(:))) then
         continue
      else
         write(*,*) 'When building B,P matrices, must call this ',&
                    'routine with either ALL or NO modes included.'
         call &
         AbortWithError('AccumulateProducts_CP8(): unexpected B,P state')
      endif

!     For convergence check: multiply coefs with inner product arrays
!     Accumulate inner products:
!      <WF,F> or <WAF,AF> from X%B
!      <WG,F> or  <WG,AF> from X%P
!     Divide out coefs to prepare for solver
      if (X%chkconv) then

!        ALS, wALS: only compute AFAF on 1st call
!        LS, wLS (fast): accumulate AFAF
         if (.not.X%AisI .or. X%accum_conv_ct.eq.0) then
            call UPDATEPCoef_CP8(X%VB,WAF,AF,X%B,.FALSE.)
            X%AFAF=X%AFAF+ReduceP_CP8(X%B,X%algo)
            call UPDATEPCoef_CP8(X%VB,WAF,AF,X%B,.TRUE.)
         endif

!        ALS, wALS: accumulate FG (intertwining case)
!        LS, wLS (fast): compute FG only on 1st call
         if (X%AisI .or. X%accum_conv_ct.eq.0) then
            call UPDATEPCoef_CP8(X%VP,WG,AF,X%P,.FALSE.)
            X%FG=X%FG+ReduceP_CP8(X%P,X%algo)
            call UPDATEPCoef_CP8(X%VP,WG,AF,X%P,.TRUE.)
         endif

!        LS, wLS: compute these inner products on 1st call only
         if (.not.X%AisI .and. X%accum_conv_ct.eq.0) then
            call UPDATEPCoef_CP8(X%VB1,AF,F,X%B1,.FALSE.)
            X%FAF=X%FAF+ReduceP_CP8(X%B1,X%algo)
            call UPDATEPCoef_CP8(X%VB1,AF,F,X%B1,.TRUE.)
            if (X%normal) then
               call UPDATEPCoef_CP8(X%VB2,F,F,X%B2,.FALSE.)
               X%FF=X%FF+ReduceP_CP8(X%B2,X%algo)
               call UPDATEPCoef_CP8(X%VB2,F,F,X%B2,.TRUE.)
            endif
         endif

         X%accum_conv_ct=X%accum_conv_ct+1
      endif

      end subroutine AccumulateProducts_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RecomputeNormalEquations_CP8(X,A,ish,Esh,W,G,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Recomputes A^T*W*A and A^T*(W*G), as necessary, when A,G,W changes. 

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(in) :: A,W,G
      TYPE (CP8) :: Wmult
      integer, intent(in) :: ish,d
      real(kind=8), intent(in) :: Esh

!     Recompute A^T*W*A for mode d
      IF (X%LS%useA .and. X%LS%useW) THEN
         IF (X%Achange .or. X%Wchange) THEN ! wLS
            IF (X%normal) THEN
               call CPMM_CP8(W,0,0.0,.FALSE.,A,ish,Esh,.FALSE.,Wmult,d,.FALSE.,X%algo)
               call CPMM_CP8(A,ish,Esh,.TRUE.,Wmult,0,0.0,.FALSE.,X%ATWA,d,.TRUE.,X%algo)
               call Wmult%flush()
            ELSE
               call CPMM_CP8(W,0,0.0,.FALSE.,A,ish,Esh,.FALSE.,X%ATWA,d,.FALSE.,X%algo)
            ENDIF
            call PrepG_CP8(X%ATWA,d,X%algo)
         ENDIF
      ELSEIF (X%LS%useA) THEN ! LS
         IF (X%Achange) THEN
            IF (X%normal) THEN
               call CPMM_CP8(A,ish,Esh,.TRUE.,A,ish,Esh,.FALSE.,X%ATWA,d,.TRUE.,X%algo)
            ELSE
               call Wmult%identity(A%rows,A%rows)
               if (X%algo.eq.1) call Wmult%copyintodevice()
               call CPMM_CP8(A,ish,Esh,.FALSE.,Wmult,0,0.0,.FALSE.,X%ATWA,d,.TRUE.,X%algo)
               call Wmult%flush()
            ENDIF
            call PrepG_CP8(X%ATWA,d,X%algo)
         ENDIF
      ELSEIF (X%LS%useW) THEN
         IF (X%Wchange) THEN ! wALS
            call X%ATWA%copy_1mode(W,d)
            call PrepG_CP8(X%ATWA,d,X%algo)
         ENDIF
      ENDIF

!     Recompute A^T*W*G for mode d
      IF (X%LS%useA .and. X%LS%useW) THEN ! wLS
         IF (X%Achange .or. X%Wchange .or. X%Gchange) THEN
            IF (X%normal) THEN
               call CPMM_CP8(W,0,0.0,.FALSE.,G,0,0.0,.FALSE.,Wmult,d,.TRUE.,X%algo)
               call CPMM_CP8(A,ish,Esh,.TRUE.,Wmult,0,0.0,.FALSE.,X%ATWG,d,.TRUE.,X%algo)
               call Wmult%flush()
            ELSE
               call CPMM_CP8(W,0,0.0,.FALSE.,G,0,0.0,.FALSE.,X%ATWG,d,.TRUE.,X%algo)
            ENDIF
            call PrepG_CP8(X%ATWG,d,X%algo)
         ENDIF
      ELSEIF (X%LS%useA) THEN ! LS
         IF (X%Achange .or. X%Gchange) THEN
            IF (X%normal) THEN
               call CPMM_CP8(A,ish,Esh,.TRUE.,G,0,0.0,.FALSE.,X%ATWG,d,.TRUE.,X%algo)
            ELSE
               call X%ATWG%copy_1mode(G,d)
            ENDIF
            call PrepG_CP8(X%ATWG,d,X%algo)
         ENDIF
      ELSEIF (X%LS%useW) THEN ! wALS
         IF (X%Wchange .or. X%Gchange) THEN
            call CPMM_CP8(W,0,0.0,.FALSE.,G,0,0.0,.FALSE.,X%ATWG,d,.TRUE.,X%algo)
            call PrepG_CP8(X%ATWG,d,X%algo)
         ENDIF
      ELSE ! ALS
         IF (X%Gchange) THEN
            call X%ATWG%copy_1mode(G,d)
            call PrepG_CP8(X%ATWG,d,X%algo)
         ENDIF
      ENDIF

      end subroutine RecomputeNormalEquations_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushNormalEquations_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocates A^T*A and A^T*G after solving linear system for case when 
! A or G changes 

      implicit none
      CLASS (ALS8), intent(inout) :: X

      IF (X%Achange.or.X%Gchange.or.X%Wchange) THEN
!        Deallocate A^T*W*G as it changes before next use
         call X%ATWG%flush()
      ENDIF

      IF (X%Achange.or.X%Wchange) THEN
!        Deallocate A^T*W*A as it changes before next use
         call X%ATWA%flush()
      ENDIF

      end subroutine FlushNormalEquations_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckModeIncludedState_CP8(X,d) result(ok)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns .TRUE. if dofincluded state is .TRUE. for all modes except d,
! which is the desired state for solving the linear system

      implicit none
      CLASS (ALS8), intent(in) :: X
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

      end function CheckModeIncludedState_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ConstructLS_CP8(X,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs linear system for ALS object X, for mode d, given CP 
! vector F

      implicit none
      CLASS (ALS8), intent(inout) :: X
      integer, intent(in) :: d
      integer :: i

      IF (.not.allocated(X%B)) &
         call AbortWithError('ConstructLS_CP8(): B is not allocated')
      IF (.not.allocated(X%P)) &
         call AbortWithError('ConstructLS_CP8(): P is not allocated')
      IF (.not.X%checkmodestate(d)) THEN
         do i=1,SIZE(X%dofincluded)
            if (i.eq.d) then
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .FALSE.)'
            else
               write(*,'(A,I0,A,L0,A)') 'Mode ',i,' included in B,P? ',&
               X%dofincluded(i),' (must be .TRUE.)'
            endif
         enddo
         call &
         AbortWithError('ConstructLS_CP8(): wrong dof-included state')
      ENDIF

!     Calculate LHS of linear system
      IF (X%AisI.and.X%WisI) THEN
!        unweighted ALS: only compute LHS on 1st pass
         IF (X%accum_HS_ct.eq.0) call GetLHS_CP8(X%LS,X%B,X%LHS)
      ELSEIF (X%AisI.and.(.not.X%WisI)) THEN
!        weighted ALS: different call, but only compute LHS on 1st pass
         IF (X%accum_HS_ct.eq.0) call GetLHS_CP8(X%LS,X%B,X%ATWA,X%LHS,d)
      ELSE
!        Linear Solver: build and accumulate LHS
         IF (X%accum_HS_ct.eq.0) THEN
            call GetLHS_CP8(X%LS,X%B,X%ATWA,X%LHS,d)
         ELSE
            call AccumulateLHS_LinSys_CP8(X%LS,X%B,X%ATWA,X%LHS,d)
         ENDIF
      ENDIF

!     Calculate RHS of linear system
      IF (X%accum_HS_ct.eq.0) THEN
         call GetRHS_CP8(X%LS,X%P,X%ATWG,X%RHS,d)
      ELSEIF (X%AisI) THEN
!        ALS only: RHS is accumulated (intertwining case)
         call AccumulateRHS_CP8(X%LS,X%P,X%ATWG,X%RHS,d)
      ENDIF

      call X%flushnormal()

!     Increment the LHS-or-RHS accumulator count
      X%accum_HS_ct=X%accum_HS_ct+1

      end subroutine ConstructLS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveLS_CP8(X,F,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs linear system for ALS object X, for mode d, given CP 
! vectors F

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(inout)  :: F
      integer, intent(in) :: d
      integer :: i,rF,mn

      IF (.not.allocated(X%LHS)) &
         call AbortWithError('SolveLS_CP8(): LHS is not allocated')
      IF (.not.allocated(X%RHS)) &
         call AbortWithError('ConstLS_CP8(): RHS is not allocated')

      IF (X%AisI) THEN
         rF=F%R()
         mn=F%MN(d)
      ELSE
         rF=F%R()*F%M(d)
         mn=F%N(d)
      ENDIF

!     Solve linear system
      IF (X%AisI.and. (.not.X%WisI)) THEN
         call SolveWeightedLS8(X%LHS,X%RHS,mn,X%solver,X%penalty,X%algo)
      ELSE
         call SolveLinSys8(X%LHS,X%RHS,rF,mn,X%solver,X%penalty,X%algo)
      ENDIF

!     Put solution into F
      call PutSolninF_CP8(F,X%RHS,d,X%algo)
      call NormBaseD_CP8(F,d,.TRUE.,X%algo)

!     Check coefs of F for NaN values resulting from zero division
      IF (.NOT.F%ok()) THEN
         write(*,*) 'SolveLS(): NaN on update: mode ',d
         call AbortWithError('SolveLS() crashed')
      ENDIF

!     The GetLHS, InitLHS, and InitRHS routines initialize the LHS,RHS
!     on device, so remove them here.
      call Flush_HS_CP8(X%LHS)
      call Flush_HS_CP8(X%RHS)

!     Reset the LHS,RHS accumulator count
      X%accum_HS_ct=0

      end subroutine SolveLS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowLS_CP8(X,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out LHS and RHS for debugging

      implicit none
      CLASS (ALS8), intent(in) :: X
      integer, intent(in) :: d
      integer :: rkF,m,n,nbas

      rkF=X%LS%rkF
      nbas=X%LS%nbas(d)
      m=X%LS%m(d)
      n=X%LS%n(d)
      if (X%AisI .and. X%WisI) m=1

      if (allocated(X%LHS)) then
         write(*,*) 'LHS:'
         call PrintMatrix(X%LHS,rkF*n,rkF*m)
         write(*,*)
      endif

      if (allocated(X%RHS)) then
         write(*,*) 'RHS:'
         call PrintMatrix(X%RHS,rkF,nbas)
         write(*,*)
      endif

      end subroutine ShowLS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsDD_CP8(X,WAF,AF,F,WG,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! Call using: F, F for B;  G, F for P ( ALS)
!            WF, F for B; WG, F for P (wALS)
!            AF, F for B;  G, F for P (  LS)
!           WAF, F for B; WG, F for P ( wLS)
!            AF,AF for B;  G,AF for P (  LS, normal form)
!           WAF,AF for B; WG,AF for P ( wLS, normal form)

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(in) :: WAF,AF,F,WG
      integer, intent(in)   :: d

!     Make sure B, P are in expected state (all modes included)
      if (.not.(ALL(X%dofincluded(:)))) then
         write(*,*) 'When building B,P matrices, must call this ',&
                    'routine with ALL modes included.'
         call AbortWithError('ProdMatsDD(): unexpected B,P state')
      endif

      call UPDATE_P_CP8(X%VB,WAF,AF,d,X%B,.TRUE.)
      call UPDATE_P_CP8(X%VP,WG,AF,d,X%P,.TRUE.)
      if ((.not.X%AisI).and.X%chkconv) then
         call UPDATE_P_CP8(X%VB1,AF,F,d,X%B1,.TRUE.)
         if (X%normal) then
            call UPDATE_P_CP8(X%VB2,F,F,d,X%B2,.TRUE.)
         endif
      endif

      X%dofincluded(d)=.FALSE.

      end subroutine ProdMatsDD_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProdMatsUD_CP8(X,WAF,AF,F,WG,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets B, P matrices by constructing from scratch or by downdating
! Call using: F, F for B;  G, F for P ( ALS)
!            WF, F for B; WG, F for P (wALS)
!            AF, F for B;  G, F for P (  LS)
!           WAF, F for B; WG, F for P ( wLS)
!            AF,AF for B;  G,AF for P (  LS, normal form)
!           WAF,AF for B; WG,AF for P ( wLS, normal form)

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(in) :: WAF,AF,F,WG
      integer, intent(in)    :: d

!     Make sure B, P are in expected state (all modes except d included)
      if (.not.X%checkmodestate(d)) then
         write(*,*) 'All modes must be included except mode d'
         call AbortWithError('ProdMatsUD(): unexpected B,P state')
      endif

!     Update or scrub the B, P matrices
      if (X%update) then
         call UPDATE_P_CP8(X%VB,WAF,AF,d,X%B,.FALSE.)
         call UPDATE_P_CP8(X%VP,WG,AF,d,X%P,.FALSE.)
         if ((.not.X%AisI).and.X%chkconv) then
            call UPDATE_P_CP8(X%VB1,AF,F,d,X%B1,.FALSE.)
            if (X%normal) then
               call UPDATE_P_CP8(X%VB2,F,F,d,X%B2,.FALSE.)
            endif
         endif

         X%dofincluded(d)=.TRUE.
      else
!        Reset all modes to "not included" but do not deallocate here
         X%dofincluded(:)=.FALSE.
      endif

      end subroutine ProdMatsUD_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowProductMatrices_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out LHS and RHS for debugging

      implicit none
      CLASS (ALS8), intent(in) :: X
      integer :: rA,rAp,rF,rG,rW

      rA=X%LS%rksA
      rAp=X%LS%rksAp
      rF=X%LS%rkF
      rG=X%LS%rkG
      rW=X%LS%rkW


      if (allocated(X%B)) then
         write(*,*) 'B:'
         call PrintMatrix(X%B,rW*rA*rF,rAp*rF)
         write(*,*)
      endif

      if (allocated(X%P)) then
         write(*,*) 'P:'
         call PrintMatrix(X%P,rW*rG,rAp*rF)
         write(*,*)
      endif

      if (allocated(X%B1)) then
         write(*,*) 'B1:'
         call PrintMatrix(X%B1,rAp*rF,rF)
         write(*,*)
      endif

      if (allocated(X%B2)) then
         write(*,*) 'B2:'
         call PrintMatrix(X%B2,rF,rF)
         write(*,*)
      endif

      end subroutine ShowProductMatrices_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckALSConvergence_CP8(X,F,WG,G) result(converged)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||F-G|| or ||AF-G|| as a convergence test
! Call using A*F instead of F when using linear solver instead of ALS
! G is actually WG in this call

      implicit none
      CLASS (ALS8), intent(inout) :: X
      TYPE (CP8), intent(in) :: F,WG,G
      integer :: converged,i,ish
      real(kind=8) :: conver,Esh

      converged=0

      IF (.not.X%chkconv) RETURN

      ish=X%LS%ish
      Esh=X%LS%Esh !!! Modify to account for tiling

!     Compute <G,WG> if not done already
      IF (X%calcGG .and. X%GG.le.0.d0) X%GG=PRODVV_CP8(G,WG,X%algo)

!     Calculate the coefficient norm (ALS only)
      IF (X%AisI) THEN
         X%CN=GetCNorm_CP8(F,X%algo)
!     Calculate the overlaps needed for Rayleigh Quotient (LS only)
      ELSE
         IF (X%normal) THEN
!           Normal form: <F,AF> in X%B1, <F,F> in X%B2
            X%RQ=(X%FAF+ish*Esh*X%FF)/X%FF
            if (ish.lt.0) X%RQ=-X%RQ
         ELSE
!           Fast form: <F,AF> in X%B, <F,F> in X%B1
            X%RQ=(X%AFAF+ish*Esh*X%FAF)/X%FAF
            if (ish.lt.0) X%RQ=-X%RQ
         ENDIF
      ENDIF

!     Compute result and check convergence. The more important
!     condition is error < thresh, so always set convergence to 1 if
!     this is achieved
      converged=0
      IF (X%AisI.or.X%normal) THEN
         conver=sqrt(abs((X%AFAF+X%GG-2*X%FG)/X%GG))
      ELSE
         conver=X%RQ
      ENDIF

      X%del=abs((conver-X%conver)/conver)
      IF (X%del.lt.X%dthresh) THEN
         converged=2
      ENDIF
      IF (conver.lt.X%thresh) THEN
         converged=1
      ENDIF
      X%conver=conver

      X%accum_conv_ct=0

      IF (X%prtconv) call X%showconv()

      end function CheckALSConvergence_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintALSConvergence_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out ALS or linear solver convergence information in ALS object

      implicit none
      CLASS (ALS8), intent(in) :: X
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
      IF (X%LS%useW) THEN
         wtag='(wtd)'
      ELSE
         wtag=''
      ENDIF

      IF (X%AisI) THEN
         frmt='(A,X,A,3(A,ES14.8),X,A)'
         write(*,frmt) &
         TRIM(ADJUSTL(tag)),TRIM(ADJUSTL(wtag)),' ||F-G||/||G|| = ',&
         X%conver,'; ||F|| = ',sqrt(abs(X%AFAF)),'; CN = ',&
         sqrt(abs(X%CN/X%AFAF)),TRIM(ADJUSTL(ctag))
      ELSEIF (X%normal) THEN
         frmt='(A,X,A,2(A,ES14.8,X),(A,F16.8,X),A))'
         write(*,frmt) TRIM(ADJUSTL(tag)),TRIM(ADJUSTL(wtag)),&
         ' ||A*F-G||/||G|| = ',X%conver,&
         '; <(W)AF,AF> = ',sqrt(abs(X%AFAF)),'; RQ = ',&
         X%RQ,TRIM(ADJUSTL(ctag))
      ELSE
         frmt='(A,X,A,(A,ES14.8,X),(A,F16.8,X),A)'
         if (mpirank.eq.0) &
!         write(*,*) ' <AF,F>= ',X%AFAF,'; <G,F> = ',X%FG,';<F,F> = ',X%FAF
         write(*,frmt) TRIM(ADJUSTL(tag)),TRIM(ADJUSTL(wtag)),&
         ' |<AF,F>-<G,F>|/||F|| = ',abs((X%AFAF-X%FG)/X%FAF),'; RQ = ',&
         X%RQ,TRIM(ADJUSTL(ctag))
      ENDIF

      end subroutine PrintALSConvergence_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ResetALSConvergence_CP8(X)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resets accumulated ALS convergence metrics

      implicit none
      CLASS (ALS8) :: X

      X%AFAF=0.d0
      X%FG=0.d0
      X%FAF=0.d0
      X%FF=0.d0

      end subroutine ResetALSConvergence_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE ALSOO8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
