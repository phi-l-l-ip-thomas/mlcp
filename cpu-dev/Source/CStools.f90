!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CSTOOLS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE GENCONFIG
      USE SEPDREPN
      USE CPCONFIG

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewCSGuess(C,X,nbas,rkx,maxnmode,nmcfg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates initial guess for CSPES from scratch

      implicit none
      TYPE (Configs), INTENT(OUT) :: C
      TYPE (CP), INTENT(OUT)   :: X
      integer, intent(in)  :: nbas(:)
      integer, intent(in)  :: rkx,maxnmode
      logical, intent(in)  :: nmcfg
      integer, allocatable :: cfg(:)
      integer :: i,ndof
      logical :: success

      ndof=SIZE(nbas)

!     Guess vector
      ALLOCATE(cfg(ndof))
      cfg(:)=1
      call NewConfigs(C,nbas,rkx)
      C%coef=1.d0
      C%qns(1,:)=cfg(:)
      DO i=2,rkx
         call nextconfig(cfg,nbas,maxnmode,nmcfg,success)
         C%qns(i,:)=cfg(:)
      ENDDO

!     Initialize zero solution vector for BP-sigma
      X=ZeroCPvec(nbas)

      end subroutine NewCSGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ModifyCSGuess(C,X,rkx,maxnmode,nmcfg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Modifies CSPES guess read from file, if needed

      implicit none
      TYPE (Configs), INTENT(INOUT) :: C
      TYPE (CP), INTENT(OUT)     :: X
      integer, intent(in)  :: rkx,maxnmode
      logical, intent(in)  :: nmcfg
      integer, allocatable :: cfg(:)
      integer :: i,ndof,rkt
      logical :: success

      ndof=SIZE(C%nbas)
      rkt=SIZE(C%coef)

!     Change the size of C if needed
      IF (rkt.ne.rkx) THEN
         call ResizeConfigList(C,rkx)

!        Add new configs if the size of C is to be increased
         IF (rkt.lt.rkx) THEN
!           Send the last config of the existing C to nextconfig() since
!           this function generates new configs by recursion
            ALLOCATE(cfg(ndof))
            cfg(:)=C%qns(rkt,:)
            DO i=rkt+1,rkx
               call nextconfig(cfg,C%nbas,maxnmode,nmcfg,success)
               C%qns(i,:)=cfg(:)
               C%coef(i)=0.d0
            ENDDO
         ENDIF
      ENDIF

!     Initialize X as CP-representation of C
      call Config2CP(X,C)

!     Change coefs of C array to 1.d0
      C%coef(:)=1.d0

      end subroutine ModifyCSGuess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PESCPbracket(XC,GC,SC,fx,ax,bx,cx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates values of the rotation angle to bracket the minimum of the
! function, GetXonGline()
! See Bracketmethod() in Numerical Recipes, 3rd Ed (2007) p. 491

      implicit none
      real*8, intent(in)  :: XC(:),GC(:)
      real*8, intent(out) :: SC(:)
      real*8, intent(in)  :: fx
      real*8, intent(out) :: ax,bx,cx
      real*8, parameter :: GRAT=1.61803398875
      real*8, parameter :: GLIM=100.
      real*8, parameter :: smallnr=1.d-18
      real*8 :: q,r,u,fu,ulim,fa,fb,fc,lgrad
      integer :: fcalls

      fcalls=0

!     ax and bx are the first 2 points of the bracketing trio
      ax=0.d0
      bx=1.d-2
      fa=fx

      call GetXonGline(XC,GC,SC,bx,fb)
      fcalls=fcalls+1

!     Swap a and b if fb > fa
      IF (fb.gt.fa) THEN
         u=ax
         fu=fa
         ax=bx
         fa=fb
         bx=u
         fb=fu
      ENDIF

!     Guess for the third bracketing point
      cx=bx+GRAT*(bx-ax)
      call GetXonGline(XC,GC,SC,cx,fc)
      fcalls=fcalls+1

!     Refine the bracketing points until they enclose a minimum
      DO
         IF (fb.lt.fc) EXIT  ! Success

!        Parabolic extrapolation to guess next point
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2*SIGN(MAX(abs(q-r),smallnr),q-r))
         ulim=bx+GLIM*(cx-bx)

!        Bracketing cx is between b and current cx
         IF ((bx-u)*(u-cx).gt.0.d0) THEN
            call GetXonGline(XC,GC,SC,u,fu)
            fcalls=fcalls+1

!           Minimum between bx and cx
            IF (fu.lt.fc) THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               EXIT

!           Minimum between ax and u
            ELSEIF (fu.gt.fb) THEN
               cx=u
               fc=fu
               EXIT
            ENDIF

!           Discard parabolic guess
            u=cx+GRAT*(cx-bx)
            call GetXonGline(XC,GC,SC,u,fu)
            fcalls=fcalls+1

!        Bracketing cx between current cx and the max value
         ELSEIF ((cx-u)*(u-ulim).gt.0.d0) THEN
            call GetXonGline(XC,GC,SC,u,fu)
            fcalls=fcalls+1

            IF (fu.lt.fc) THEN
                bx=cx
                fb=fc
                cx=u
                fc=fu
                u=u+GRAT*(u-bx)
                call GetXonGline(XC,GC,SC,u,fu)
                fcalls=fcalls+1
            ENDIF

!        Use maximum permissible value of u to guess next cx
         ELSEIF ((u-ulim)*(ulim-cx).ge.0.d0) THEN
            u=ulim
            call GetXonGline(XC,GC,SC,u,fu)
            fcalls=fcalls+1

!        Use Golden Ratio to guess next cx
         ELSE
            u=cx+GRAT*(cx-bx)
            call GetXonGline(XC,GC,SC,u,fu)
            fcalls=fcalls+1
         ENDIF

         ax=bx
         fa=fb
         bx=cx
         fb=fc
         cx=u
         fc=fu
      ENDDO

      end subroutine PESCPbracket

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PESCPminimize(XC,GC,SC,ax,bx,cx,xmin,fmin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds the minimum of the function, GetXonGline(), via Brent's method
! See brent() in Numerical Recipes, 3rd Ed (2007) p. 498

      implicit none
      real*8, intent(in)  :: XC(:),GC(:)
      real*8, intent(out) :: SC(:)
      real*8, intent(in)  :: ax,bx,cx
      real*8, intent(out) :: xmin,fmin
      integer, parameter  :: maxit=100
      real*8, parameter   :: smallnr=1.d-16
      real*8, parameter   :: tol=1.0d-15
      real*8, parameter   :: CGOLD=0.38196601125 
      real*8 :: a,b,u,v,w,x,xm,fu,fv,fw,fx,dd,p,q,r
      real*8 :: tol1,tol2,ee,etmp
      integer :: k,fcalls

      fcalls=0

!     Initializations
      dd=0.d0
      ee=0.d0
      a=MIN(ax,cx)
      b=MAX(ax,cx)
      x=bx
      v=bx
      w=bx

!     Evaluate function and derivative
      call GetXonGline(XC,GC,SC,x,fx)
      fcalls=fcalls+1
      fv=fx
      fw=fx

!     Main loop over iterations
      DO k=1,maxit

         xm=0.5d0*(a+b)
         tol1=tol*abs(x)+smallnr
         tol2=2*tol1

!        If next point is close enough to this one, min is converged
         IF (abs(x-xm).le.(tol2-0.5d0*(b-a))) THEN
            fmin=fx
            xmin=x
            EXIT
         ENDIF

         IF (abs(ee).gt.tol1) THEN
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            IF (q.gt.0.d0) p=-p
            q=abs(q)
            etmp=ee
            ee=dd

            IF (abs(p).ge.abs(0.5d0*q*etmp) .or. p.le.q*(a-x) &
                                            .or. p.ge.q*(b-x)) THEN
               IF (x.ge.xm) THEN
                  ee=a-x
               ELSE
                  ee=b-x
               ENDIF
               dd=CGOLD*ee
            ELSE
               dd=p/q
               u=x+dd
               IF (u-a.le.tol2 .or. b-u.lt.tol2) dd=SIGN(tol1,xm-x)
            ENDIF
         ELSE
            IF (x.ge.xm) THEN
               ee=a-x
            ELSE
               ee=b-x
            ENDIF
            dd=CGOLD*ee
         ENDIF

!        Update u and evaluate the objective function
         IF (abs(dd).ge.tol1) THEN
            u=x+dd
         ELSE
            u=x+SIGN(tol1,dd)
         ENDIF
         call GetXonGline(XC,GC,SC,u,fu)
         fcalls=fcalls+1

!        Update x and f values
         IF (fu.le.fx) THEN
            IF (u.ge.x) THEN
               a=x
            ELSE
               b=x
            ENDIF
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu

         ELSE
            IF (u.lt.x) THEN
               a=u
            ELSE
               b=u
            ENDIF
            IF (fu.le.fw .or. w.eq.x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSEIF (fu.lt.fv .or. v.eq.x .or. v.eq.w) THEN
               v=u
               fv=fu
            ENDIF
         ENDIF
      ENDDO

      end subroutine PESCPminimize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetXonGline(XC,GC,SC,step,f)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the new value of CP-vector X, the residual, and the objective
! function value along the gradient line

      implicit none
      real*8, intent(in)  :: XC(:),GC(:)
      real*8, intent(out) :: SC(:)
      real*8, intent(in)  :: step
      real*8, intent(out) :: f

!     Compute the new coefficients and objective function value
      SC(:)=XC(:)+step*GC(:)
      f=0.5*dot_product(SC(:),SC(:))

      end subroutine GetXonGline

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CSTOOLS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
