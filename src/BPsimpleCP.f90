!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE BPSIMPLECP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE SEPDREPN
      USE CPCONFIG
      USE CSTOOLS
      USE RESTARTCSPES
      USE MODVECVEC

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SolveBPsimpleCP(C,B,A,X,R,G,csp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Minimizes ||Ax - b||_2 with respect to x. X is modified along the way
! by replacing configurations with small coefficients with new ones
! generated using nextconfig()

      implicit none
      TYPE (CSpar), INTENT(IN)      :: csp
      TYPE (Configs), INTENT(INOUT) :: C
      TYPE (CP), INTENT(IN)      :: A,B
      TYPE (CP), INTENT(OUT)     :: X,G,R
      TYPE (Configs) :: Xc,Gc,Gc2,Gc3
      TYPE (CP)   :: Ax,Xnew,Gnew,Rnew,Xt
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
      call GetResidual(B,Ax,R,fnew)
      call FlushCP(Ax)
      fold=fnew
      fdelta=1.d0

!     Calculate the gradient: G = -A^T*R
      call GetCPGradient(R,A,G)

!     Project gradient onto guess configurations
      call ConfigSetOverlap(C,Gc,G)

!     Save 2 copies of the gradient configs for CG minimization
      call CopyConfigsWtoV(Gc2,Gc)
      call CopyConfigsWtoV(Gc3,Gc)

      call FlushCP(G)
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
         call CalcBvecOperlaps(G,Rnew,A,B,X,Xnew,fnew)

!        Project solution onto guess configurations
         call FlushConfigs(Xc)
         call ConfigSetOverlap(C,Xc,Xnew)
         call SavePESCoef(csp,Xc)
         call FlushCP(Xnew)
         call Config2CP(Xnew,Xc)

         l1=onenorm(Xc)
         l0=zeronorm(Xc,ctol) 

!        Calculate the residual and objective function value
!        for the projected Xnew 
         call CPMatVecProd(A,Xnew,Ax,.FALSE.)
         call GetResidual(B,Ax,Rnew,fnew)
         call FlushCP(Ax)

         fdelta=(fold-fnew)/fnew
         write(*,'(I10,3(X,ES12.4),X,I8)') &
              iter,sqrt(2*fnew),fdelta,l1,l0
         fold=fnew

!        Get the new gradient
         call GetCPGradient(Rnew,A,Gnew)

!        Project gradient onto guess configurations
         call FlushConfigs(Gc)
         call ConfigSetOverlap(C,Gc,Gnew)
         call FlushCP(Gnew)

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
         call ReplaceVwithW(R,Rnew)
      ENDDO

      call ordre(X)
      call TrimZeros(X,sigma)

      write(*,'(X,A/)') 'Solution vector, in terms of guess configs:'
      call PrintConfigs(Xc)
!      write(*,*) 'Residual:'
!      call PrintCPvec(R)
!      write(*,'(/X,A/)') 'Gradient, in terms of guess configs:'
!      call PrintConfigs(Gc)

      call FlushConfigs(Xc)
      call FlushConfigs(Gc)
      call FlushConfigs(Gc2)
      call FlushConfigs(Gc3)

      end subroutine SolveBPsimpleCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetCPGradient(R,A,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the gradient of the CP-format CS solution

      implicit none
      TYPE (CP), INTENT(IN)   :: R,A
      TYPE (CP), INTENT(OUT)  :: G

!     Convert to gradient using matrix-vector product and sign change
      call CPMatVecProd(A,R,G,.TRUE.)
      call VecSignChange(G,1,SIZE(G%coef))

      end subroutine GetCPGradient

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetResidual(b,Ax,r,f)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the residual r = b - Ax, in CP-format, and its l-2 norm.

      implicit none
      TYPE (CP), INTENT(IN)  :: Ax,b
      TYPE (CP), INTENT(OUT) :: r
      real*8, intent(out) :: f
      integer :: i,j,rkb,rkA

      rkb=SIZE(b%coef)
      rkA=SIZE(Ax%coef)

      r=CopyCP(b)

      f=0.d0
!$omp parallel
!$omp do private(i,j) schedule(static)
      DO i=1,rkb
         DO j=1,rkA
            r%coef(i)=r%coef(i)-Ax%coef(j)*&
                      PRODND(r%base(:,i),Ax%base(:,j),r%nbas,0)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel

      DO i=1,rkb
         IF (r%coef(i).lt.0) THEN
            r%coef(i)=abs(r%coef(i))
            call VecSignChange(r,i,i)
         ENDIF
         f=f+r%coef(i)**2
      ENDDO

      f=0.5d0*f

      end subroutine GetResidual

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CalcBvecOperlaps(G,R,A,B,X,Xnew,fnew)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the residual r = b - Ax, in CP-format, and its l-2 norm.

      implicit none
      TYPE (CP), INTENT(IN)  :: X,G,A,B
      TYPE (CP), INTENT(OUT) :: Xnew,R
      TYPE (CP) :: Ax,Ag
      real*8, intent(inout) :: fnew
      real*8, allocatable   :: bxcoef(:),bgcoef(:),solncoef(:)
      real*8  :: aax,bx,cx,step,rav,gav,fac
      integer :: i,j,rkb,rkx,rkg

      rkb=SIZE(B%coef)
      rkg=SIZE(G%coef)
      rkx=SIZE(X%coef)

!     Get rid of the old Xnew and residual, if they are allocated
      call FlushCP(Xnew)
      call FlushCP(R)

      ALLOCATE(bxcoef(rkb),bgcoef(rkb),solncoef(rkb))
      bxcoef(:)=B%coef(:)
      bgcoef(:)=0.d0

      call CPMatVecProd(A,X,Ax,.FALSE.)
      call CPMatVecProd(A,G,Ag,.FALSE.)

!     Calculate the overlap of B with Ax and Ag
!$omp parallel
!$omp do private(i,j)
      DO i=1,rkb
         DO j=1,rkx
            bxcoef(i)=bxcoef(i)-Ax%coef(j)*&
                      PRODND(B%base(:,i),Ax%base(:,j),B%nbas,0)
         ENDDO
      ENDDO
!$omp end do
!$omp do private(i,j)
      DO i=1,rkb
         DO j=1,rkg
            bgcoef(i)=bgcoef(i)-Ag%coef(j)*&
                      PRODND(B%base(:,i),Ag%base(:,j),B%nbas,0)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel

      call FlushCP(Ax)
      call FlushCP(Ag)

!     For a system with many DOFs, the coefs-of-Bg may be less than
!     (machine-precision)*coefs-of-R, in which case Bg will get rounded
!     to zero in the line search algorithm. To prevent this, average the
!     coefs of R and Bg, and scale Bg and step accordingly
      rav=0.d0
      gav=0.d0
      DO i=1,rkb
         rav=rav+abs(bxcoef(i))
         gav=gav+abs(bgcoef(i))
      ENDDO
      fac=rav/gav
      bgcoef(:)=bgcoef(:)*fac

!     NR line search
      call PESCPbracket(bxcoef,bgcoef,solncoef,fnew,aax,bx,cx)
      call PESCPminimize(bxcoef,bgcoef,solncoef,aax,bx,cx,step,fnew)

!     Calculate Xnew = X + step*G
      Xnew=CopyCP(X)
      call SUMVECVEC(Xnew,1.d0,G,fac*step)

!     Calculate Rnew = B with new coefs
      R=CopyCP(B)
      R%coef(:)=solncoef(:)
!     Change sign to keep coefs positive
      DO i=1,rkb
         IF (R%coef(i).lt.0) THEN
            R%coef(i)=abs(R%coef(i))
            call VecSignChange(R,i,i)
         ENDIF
      ENDDO

      DEALLOCATE(bxcoef,bgcoef,solncoef)

      end subroutine CalcBvecOperlaps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE BPSIMPLECP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
