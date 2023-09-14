!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BPSIMPLECF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE SEPDREPN
      USE CPCONFIG
      USE CSTOOLS
      USE RESTARTCSPES
      USE MODVECVEC
      USE REDORTHO

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveBPsimpleCF(C,B,A,X,R,G,csp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Minimizes ||Ax - b||_2 with respect to x. X is modified along the way
! by replacing configurations with small coefficients with new ones
! generated using nextconfig()

      implicit none
      TYPE (CSpar), INTENT(IN)      :: csp
      TYPE (Configs), INTENT(INOUT) :: C
      TYPE (Configs), INTENT(IN)    :: B    ! Data to fit
      TYPE (Configs), INTENT(OUT)   :: R    ! Residual of fit
      TYPE (CPvec), INTENT(IN)      :: A    ! Transformation matrix
      TYPE (CPvec), INTENT(OUT)     :: X,G  ! Solution, Gradient
      TYPE (Configs) :: Rnew,Xc,Gc,Gc2,Gc3
      TYPE (CPvec)   :: Ax,Xnew,Gnew,Xt
      integer :: iter,maxit,rkx,i,maxnmode,l0
      real*8  :: sigma,gtol,ctol,dtol,fold,fnew,fdelta
      real*8  :: gg,dgg,gam,l1
      character(len=72) :: messg

!     Set parameters from input file
      rkx=csp%pesrank
      maxnmode=csp%maxnmode
      maxit=csp%ncycle
      sigma=csp%sigma ! Converge coefs to within sigma
      ctol=csp%ctol   ! Configuration considered zero if below ctol
      gtol=csp%gtol   ! Gradient considered zero if below gtol
      dtol=csp%dtol   ! Tol on fractional decrease of obj. fxn.

!     Initial values
      iter=0

      IF (ALLOCATED(C%coef)) THEN
!        Restart if C already exists
         call ModifyCSGuess(C,X,rkx,maxnmode,csp%nmcfg)
         call ConfigSetOverlap(C,Xc,X)
      ELSE
!        Initial guess for C and X
         call NewCSGuess(C,X,B%nbas,rkx,maxnmode,csp%nmcfg)
      ENDIF

      write(*,'(/X,A/)') 'Generating guess configs...'
      call PrintConfigs(C)
      write(*,*)

!     Initialize residual config. list and objective function value
      call CPMatVecProd(A,X,Ax,.FALSE.)
      call GetResidualCF(B,Ax,R,fnew)
      call FlushCPvec(Ax)
      fold=fnew
      fdelta=1.d0

!     Calculate the gradient: G = -A^T*R
      call GetCFGradient(R,A,G)

!     Project gradient onto guess configurations
      call ConfigSetOverlap(C,Gc,G)

!     Save 2 copies of the gradient configs for CG minimization
      call CopyConfigsWtoV(Gc2,Gc)
      call CopyConfigsWtoV(Gc3,Gc)

      call FlushCPvec(G)
!      write(*,'(X,A/)') 'Projection of initial G onto guess configs:'
!      call PrintConfigs(Gc)
      call Config2CP(G,Gc)

      write(*,'(X,A)') 'CG optimization of potential coefficients...'
      write(*,'(/X,A,3(3X,A,3X),2X,A)') 'Iteration',&
            '||R||_2','f-delta','||X||_1','||X||_0'

      DO
         iter=iter+1

!        Exit conditions
         IF (sqrt(2*fnew).lt.sigma) THEN
            write(*,'(/X,A/)') 'Successful exit due to small residual!'
            EXIT
         ELSEIF (SIZE(G%coef).eq.1 .and. G%coef(1).eq.0.d0) THEN
            write(*,'(/X,A/)') 'EXIT due to zero projected gradient'
            EXIT
         ELSEIF (fdelta.lt.dtol) THEN
            write(*,'(/X,A/)') 'EXIT due to small decrease of residual'
            EXIT
         ELSEIF (iter.ge.maxit) THEN
            write(*,'(/X,A/)') 'EXIT due to too many iterations'
            EXIT
         ENDIF

!        Optimize X by projection onto B
         call CalcBvecOperlapCF(G,Rnew,A,B,X,Xnew,fnew)

!        Project solution onto guess configurations
         call FlushConfigs(Xc)
         call ConfigSetOverlap(C,Xc,Xnew)
         call SavePESCoef(csp,Xc)
         call FlushCPvec(Xnew)
         call Config2CP(Xnew,Xc)

         l1=onenorm(Xc)
         l0=zeronorm(Xc,ctol) 

!        Calculate the residual and objective function value
!        for the projected Xnew 
         call CPMatVecProd(A,Xnew,Ax,.FALSE.)
         call GetResidualCF(B,Ax,Rnew,fnew)
         call FlushCPvec(Ax)

         fdelta=(fold-fnew)/fnew
         write(*,'(I10,3(X,ES12.4),X,I8)') &
              iter,sqrt(2*fnew),fdelta,l1,l0
         fold=fnew

!        Get the new gradient
         call GetCFGradient(Rnew,A,Gnew)

!        Project gradient onto guess configurations
         call FlushConfigs(Gc)
         call ConfigSetOverlap(C,Gc,Gnew)
         call FlushCPvec(Gnew)

!        Modify the line search direction using the Polak-Ribiere update
         gg=0.d0
         dgg=0.d0
         DO i=1,SIZE(Gc%coef)
            gg=gg+Gc2%coef(i)**2
            dgg=dgg+Gc%coef(i)*(Gc%coef(i)+Gc2%coef(i))
         ENDDO
         gam=dgg/gg
         DO i=1,SIZE(Gc%coef)
            Gc2%coef(i)=-Gc%coef(i)
            Gc3%coef(i)=Gc2%coef(i)+gam*Gc3%coef(i)
            Gc%coef(i)=Gc3%coef(i)
         ENDDO

         call Config2CP(Gnew,Gc)

!        Replace X <-- Xnew, R <-- Rnew, G <-- Gnew
         call ReplaceVwithW(X,Xnew)
         call ReplaceVwithW(G,Gnew)
         call ReplaceConfigsVwithW(R,Rnew)
      ENDDO

      call ordre(X)
      call TrimZeros(X,sigma)

!      write(*,'(X,A/)') 'Solution vector, in terms of guess configs:'
!      call PrintConfigs(Xc)
!      write(*,'(/X,A/)') 'Residual:'
!      call PrintConfigs(R)
!      write(*,'(/X,A/)') 'Gradient, in terms of guess configs:'
!      call PrintConfigs(Gc)

      call FlushConfigs(Xc)
      call FlushConfigs(Gc)
      call FlushConfigs(Gc2)
      call FlushConfigs(Gc3)

      end subroutine SolveBPsimpleCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetCFGradient(R,A,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the gradient of the CP-format CS solution

      implicit none
      TYPE (Configs), INTENT(IN) :: R
      TYPE (CPvec), INTENT(IN)   :: A
      TYPE (CPvec), INTENT(OUT)  :: G
      TYPE (CPvec) :: Rcp

      call Config2CP(Rcp,R)

!     Convert to gradient using matrix-vector product and sign change
      call CPMatVecProd(A,Rcp,G,.TRUE.)
      call VecSignChange(G,1,SIZE(G%coef))
      call FlushCPvec(Rcp)

      end subroutine GetCFGradient

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetResidualCF(b,Ax,r,f)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the residual r = b - Ax, in configuration format, along
! with its l-2 norm.

      implicit none
      TYPE (Configs), INTENT(IN)  :: b
      TYPE (Configs), INTENT(OUT) :: r
      TYPE (CPvec), INTENT(IN)    :: Ax
      real*8, intent(out) :: f
      integer :: i,rk

      rk=SIZE(b%coef)

!     Get the overlap of Ax with b
      call ConfigSetOverlap(b,r,Ax)

!     Calculate b - Ax for the configurations in r
      DO i=1,rk
         r%coef(i)=b%coef(i)-r%coef(i)
      ENDDO

!     Calculate the l-2 norm
      f=0.5d0*twonorm(r)**2

      end subroutine GetResidualCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CalcBvecOperlapCF(G,R,A,B,X,Xnew,fnew)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the residual r = b - Ax, in CP-format, and its l-2 norm.

      implicit none
      TYPE (Configs), INTENT(IN)  :: B
      TYPE (Configs), INTENT(OUT) :: R
      TYPE (Configs) :: Bg
      TYPE (CPvec), INTENT(IN)  :: X,G,A
      TYPE (CPvec), INTENT(OUT) :: Xnew
      TYPE (CPvec) :: Ax,Ag
      real*8, intent(inout) :: fnew
      real*8, allocatable   :: solncoef(:)
      real*8  :: aax,bx,cx,step,l2,rav,gav,fac
      integer :: rkb,i

      rkb=SIZE(B%coef)

!     Get rid of the old Xnew and residual, if they are allocated
      call FlushCPvec(Xnew)
      call FlushConfigs(R)

!     Matrix vector product to put X,G in basis of B
      call CPMatVecProd(A,X,Ax,.FALSE.)
      call CPMatVecProd(A,G,Ag,.FALSE.)

!     Calculate the overlap of B with Ax and Ag
      call GetResidualCF(B,Ax,R,l2)
      call ConfigSetOverlap(B,Bg,Ag)
      Bg%coef(:)=-Bg%coef(:)

!     For a system with many DOFs, the coefs-of-Bg may be less than
!     (machine-precision)*coefs-of-R, in which case Bg will get rounded
!     to zero in the line search algorithm. To prevent this, average the
!     coefs of R and Bg, and scale Bg and step accordingly
      rav=0.d0
      gav=0.d0
      DO i=1,rkb
         rav=rav+abs(R%coef(i))
         gav=gav+abs(Bg%coef(i))
      ENDDO
      fac=rav/gav
      Bg%coef(:)=Bg%coef(:)*fac

      call FlushCPvec(Ax)
      call FlushCPvec(Ag)

!     NR line search
      ALLOCATE(solncoef(rkb))
      call PESCPbracket(R%coef,Bg%coef,solncoef,fnew,aax,bx,cx)
      call PESCPminimize(R%coef,Bg%coef,solncoef,aax,bx,cx,step,fnew)

      call FlushConfigs(Bg)

!     New coefs of R
      R%coef(:)=solncoef(:)

!     Calculate Xnew = X + step*G
      call CopyWtoV(Xnew,X)
      call SUMVECVEC(Xnew,1.d0,G,step*fac)

      DEALLOCATE(solncoef)

      end subroutine CalcBvecOperlapCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE BPSIMPLECF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
