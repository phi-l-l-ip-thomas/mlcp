!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE HAMILOPT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rotates Hamiltonian to optimal coordinate system

      use ERRORTRAP
      use UTILS
      use MODECOMB
      use OPFUNCS
      use LINALG
      use MODVECVEC
      use CPCONFIG
      use REDORTHO

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine OptimizePESdirections(W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts config list W into CP-vector V, scales coefficients,
! rotates PES, unscales coefficients, and replaces W with rotated PES

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(INOUT) :: W(:)
      TYPE (CPvec), ALLOCATABLE :: V(:)
      real*8, allocatable :: freqs(:)
      integer :: i,j,k,ncoup,ndof

!     Set parameters
      ncoup=SIZE(W)
      ndof=W(1)%nbas(1)

!     Store potential constants in CP-vec 'V'
      ALLOCATE(V(ncoup))
      DO k=1,ncoup
         call Config2CP(V(k),W(k))
      ENDDO

!     Get the sqrt(omega_j) values for scaling
      ALLOCATE(freqs(ndof))
      DO i=1,SIZE(V(2)%coef)
         IF (W(2)%qns(i,1).eq.W(2)%qns(i,2)) THEN
            freqs(W(2)%qns(i,1))=sqrt(2*W(2)%coef(i))
         ENDIF
      ENDDO

!     Scale coefs in V by the (omega_j)^(1/2) values before rotating
      DO k=1,ncoup
!        If V(k) is a zero vector, do not scale coefficients
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE
         DO i=1,SIZE(V(k)%coef)
            DO j=1,k
               V(k)%coef(i)=V(k)%coef(i)*freqs(W(k)%qns(i,j))
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE(W)

!     Optimize H by coordinate rotation
      call RotateHamil(V)

!     Transfer rotated potential coefficients in V back to W

      ALLOCATE(W(ncoup))
      DO k=1,ncoup
         call CP2Config(V(k),W(k))
      ENDDO

!     Extract (omega_j)^(1/2) values for scaling
      freqs(:)=0.d0
      DO i=1,SIZE(W(2)%coef)
         IF (W(2)%qns(i,1).eq.W(2)%qns(i,2)) THEN
            freqs(W(2)%qns(i,1))=(2.d0*W(2)%coef(i))**0.25d0
         ENDIF
      ENDDO

!     Scale coefs in V by the (omega_j)^(1/2) values before rotating
      DO k=1,ncoup

!        If W(k) is a zero vector, skip scaling/printing
         IF (SIZE(W(k)%coef).eq.1 .and. W(k)%coef(1).eq.0.d0) CYCLE

         DO i=1,SIZE(W(k)%coef)
            DO j=1,k
               W(k)%coef(i)=W(k)%coef(i)/freqs(W(k)%qns(i,j))
            ENDDO
         ENDDO

         write(*,'(/X,A,I0/)') 'Rotated potential terms, order: ',k
         call SortConfigsByCoef(W(k))
         call PrintConfigs(W(k))

      ENDDO
      DEALLOCATE(V,freqs)

      end subroutine OptimizePESdirections

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RotateHamil(V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for transforming H into new coordinate system

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: V(:)
      TYPE (CPvec), ALLOCATABLE :: W(:)
      TYPE (Configs) :: wc
      real*8, allocatable :: theta(:),dtheta(:)
      integer :: i,j,k,m,n,nsweep,ndof,ntns,npairs,ip,nnewcoef
      real*8  :: pi,diagamp,ax,bx,cx
      real*8, parameter :: tol=1.d-7
      character*64 :: frmt

!     Set parameters
      pi=4.d0*atan(1.d0)
      nsweep=10
      ndof=V(1)%nbas(1)
      ntns=SIZE(V)
      npairs=ndof*(ndof-1)/2

      ALLOCATE(W(ntns),theta(npairs),dtheta(npairs))
      theta=0.d0

      write(*,'(/X,A/)') 'Optimizing PES by coordinate rotation...'

      DO k=1,nsweep
         write(*,*) 'nsweep = ',k
         write(frmt,'(A,I0,A)') '(X,A,5X,A,5X,A,2X,',ntns,&
                                '(X,A,I0,A),4X,A)'
         write(*,frmt) 'Pair','Rotation angle','Delta angle',&
                       ('rank(',m+1,')',m=1,ntns),'f(theta)'
         ip=1
         DO i=1,ndof-1
            DO j=i+1,ndof
!              Get the points which bracket the maximum, then maximize!
               CALL PESCPbracket(V,W,i,j,diagamp,ax,bx,cx)
               CALL PESCPminimize(V,W,i,j,ax,bx,cx,dtheta(ip),diagamp)
               theta(ip)=theta(ip)+dtheta(ip)

!              Update V
               DO m=1,ntns
                  IF (SIZE(V(m)%coef).eq.1 .and. V(m)%coef(1).eq.0.d0) &
                     CYCLE
                  call ReplaceVwithW(V(m),W(m))
               ENDDO

               write(frmt,'(A,I0,A)') '(X,2(I3,X),f15.10,X,f15.10,X,',&
                    ntns,'(X,I7),X,f23.11)'
               write(*,frmt) i,j,theta(ip)*180.d0/pi,&
                    dtheta(ip)*180.d0/pi,(SIZE(V(m)%coef),m=1,ntns),&
                    diagamp
               ip=ip+1
            ENDDO
         ENDDO
!        Test convergence of sweeps
         IF (ALL(abs(dtheta(:)).lt.tol)) EXIT
      ENDDO

!     Convert V back into a list of rank-1 unit tensors
      DO m=1,ntns
         IF (SIZE(V(m)%coef).eq.1 .and. V(m)%coef(1).eq.0.d0) CYCLE
         nnewcoef=ndof**m
         call GetConfigList(V(m),nnewcoef,wc)
         call collectconfigs(wc,.FALSE.)
         call FlushCPvec(V(m))
         call Config2CP(V(m),wc)
         call FlushConfigs(wc)
      ENDDO

      DEALLOCATE(W,theta,dtheta)

      end subroutine RotateHamil

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RotateV(V,W,i,j,dtheta,nodiag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rotates CP-vec for dofs i and j by angle dtheta (in radians), and
! computes the ratios of the sum-of-squares of the diagonal to off-
! diagonal terms 

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(IN)    :: V(:)
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: W(:)
      TYPE (CPvec) :: F,G
      TYPE (Configs)       :: wc
      integer, intent(in)  :: i,j
      real*8, intent(in)   :: dtheta
      real*8, intent(out)  :: nodiag
      integer, allocatable :: modpowr(:,:)
      real*8, allocatable  :: freqs(:),diags(:,:)
      integer :: mm,k,l,m,p,ntns,ndof,dof1,nnewcoef,gst
      real*8  :: num,denom,sclcoef

!     Parameters
      ntns=SIZE(V)
      ndof=V(1)%nbas(1)

      ALLOCATE(freqs(ndof),diags(ndof,ntns))

!     Copy and rotate each term in V, then convert to configuration
!     list,combine like terms, and convert back to CP-format in W
      nodiag=0.d0
      diags=0.d0
      DO mm=1,ntns
!        Count m=2,3...ntens,1, so that omega values are extracted
!        during the first iteration of the loop
         m=mod(mm,ntns)+1

         IF (SIZE(V(m)%coef).eq.1 .and. V(m)%coef(1).eq.0.d0) CYCLE

!        Maximize the sum-of-squares of the diagonal elements
         call FlushCPvec(W(m))
         call CopyWtoV(W(m),V(m))
         call ApplyHRotation(W(m),i,j,dtheta)

!        Extract (omega_j)^(1/2) values for scaling
         IF (m.eq.2) THEN
            call diagelems(W(m),freqs)
            DO k=1,ndof
               freqs(k)=(2.d0*freqs(k))**0.25d0
            ENDDO
         ENDIF

!        Get the diagonal potential constants, then subtract them from W
         call diagelems(W(m),diags(:,m))
         call NewConfigs(wc,V(m)%nbas,ndof)
         DO k=1,ndof
            wc%qns(k,:)=k
            wc%coef(k)=diags(k,m)
         ENDDO

         call Config2CP(F,wc)
         call FlushConfigs(wc)
         call SUMVECVEC(F,-1.d0,W(m),1.d0)

         gst=0
         DO k=1,m
            DO l=1,SIZE(F%coef)
               DO p=1,ndof
                  F%base(gst+p,l)=F%base(gst+p,l)/freqs(p)
               ENDDO
            ENDDO
            gst=gst+F%nbas(k)
         ENDDO

!        Objective fxn: Frobenius norm of off-diagonal elements
         nodiag=nodiag+PRODVV(F)/REAL(ndof**m)

         call FlushCPvec(F)
         call FlushConfigs(wc)
      ENDDO

      DEALLOCATE(freqs,diags)

      end subroutine RotateV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ApplyHRotation(V,di,dj,dtheta)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies 2D rotation of dofs di and dj by angle dtheta (in radians)

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: V
      integer, intent(in) :: di,dj
      real*8, intent(in)  :: dtheta
      integer :: i,j,k,ntns,rV,otens,gst
      real*8  :: ct,st,tempi,tempj

!     Parameters
      rV=SIZE(V%coef)
      otens=SIZE(V%nbas)

!     Get the rotation matrix elements
      ct=COS(dtheta)
      st=SIN(dtheta)

!     Rotate each of the rV terms
      DO j=1,rV
         gst=0
!        Loop over tensor order
         DO k=1,otens
!           Apply rotation to each D-dof unit in the tensor
            tempi=ct*V%base(gst+di,j)-st*V%base(gst+dj,j)
            tempj=st*V%base(gst+di,j)+ct*V%base(gst+dj,j)
            V%base(gst+di,j)=tempi
            V%base(gst+dj,j)=tempj
            gst=gst+V%nbas(k)
         ENDDO
      ENDDO

      end subroutine ApplyHRotation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PESCPbracket(V,W,i,j,fx,ax,bx,cx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates values of the rotation angle to bracket the minimum of the
! function, RotateV()
! See Bracketmethod() in Numerical Recipes, 3rd Ed (2007) p. 491

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(IN)    :: V(:)
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: W(:)
      integer, intent(in)   :: i,j
      real*8, intent(in)    :: fx
      real*8, intent(out)   :: ax,bx,cx
      real*8, parameter :: GRAT=1.61803398875
      real*8, parameter :: GLIM=100.
      real*8, parameter :: smallnr=1.d-18
      real*8 :: q,r,u,fu,ulim,fa,fb,fc,lgrad
      integer :: fcalls

      fcalls=0

!     ax and bx are the first 2 points of the bracketing trio
      ax=0.d0
      bx=1.d0
      fa=fx

      call RotateV(V,W,i,j,bx,fb)
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
      call RotateV(V,W,i,j,cx,fc)
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
            call RotateV(V,W,i,j,u,fu)
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
            call RotateV(V,W,i,j,u,fu)
            fcalls=fcalls+1

!        Bracketing cx between current cx and the max value
         ELSEIF ((cx-u)*(u-ulim).gt.0.d0) THEN
            call RotateV(V,W,i,j,u,fu)
            fcalls=fcalls+1

            IF (fu.lt.fc) THEN
                bx=cx
                fb=fc
                cx=u
                fc=fu
                u=u+GRAT*(u-bx)
                call RotateV(V,W,i,j,u,fu)
                fcalls=fcalls+1
            ENDIF

!        Use maximum permissible value of u to guess next cx
         ELSEIF ((u-ulim)*(ulim-cx).ge.0.d0) THEN
            u=ulim
            call RotateV(V,W,i,j,u,fu)
            fcalls=fcalls+1

!        Use Golden Ratio to guess next cx
         ELSE
            u=cx+GRAT*(cx-bx)
            call RotateV(V,W,i,j,u,fu)
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

      subroutine PESCPminimize(F,G,i,j,ax,bx,cx,xmin,fmin)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds the minimum of the function, RotateV(), via Brent's method
! See brent() in Numerical Recipes, 3rd Ed (2007) p. 498

      implicit none
      TYPE (CPvec), ALLOCATABLE, INTENT(IN)    :: F(:)
      TYPE (CPvec), ALLOCATABLE, INTENT(INOUT) :: G(:)
      integer, intent(in)   :: i,j
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
      call RotateV(F,G,i,j,x,fx)
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
         call RotateV(F,G,i,j,u,fu)
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

      subroutine diagelems(F,diags)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes sum-of-squares of a-th and b-th diagonal elements of a tensor 
! with all dimensions equal

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      integer :: i,j,k,l,n,rF,ndof
      real*8, intent(inout) :: diags(:)
      real*8  :: prod1

!     Parameters
      ndof=SIZE(F%nbas)
      rF=SIZE(F%coef)
      n=F%nbas(1)

!     Error if dimensions are not all the same
      DO j=2,ndof
         IF (F%nbas(j).ne.n) THEN
            write(*,*) 'n, bad-DOF, n_bad-DOF = ',n,j,F%nbas(j)
            call AbortWithError('diagelems(): unequal dimensions')
         ENDIF
      ENDDO

!     Compute the multiplicitive trace
      DO k=1,n
         diags(k)=0.d0
         DO i=1,rF
            prod1=F%coef(i)
            l=k
            DO j=1,ndof
               prod1=prod1*F%base(l,i)
               l=l+n
            ENDDO
            diags(k)=diags(k)+prod1
         ENDDO
      ENDDO

      end subroutine diagelems

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetDiagonalZPVE(diags,ZPVE)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes sum-of-squares of a-th and b-th diagonal elements of a tensor
! with all dimensions equal

      implicit none
      real*8, intent(in)  :: diags(:,:)
      real*8, intent(out) :: ZPVE
      TYPE (OperMat) :: HM,OM
      real*8, allocatable :: eigv(:),Hmat(:,:),S(:,:)
      integer :: ndof,ntns,i,j,m
      integer, parameter :: N=10

      ndof=SIZE(diags,1)
      ntns=SIZE(diags,2)

      ALLOCATE(eigv(N),S(N,N))

      ZPVE=0.d0
      DO j=1,ndof
!        Beginning with the KEO, add the potential operator matrices
         HM=GetPrimitiveOperMat(j,N,-2)
         HM%mat=diags(j,1)*HM%mat
!         call SumOperMats(HM,OM,diags(j,1))
!         call FlushOperMat(OM)
         DO m=1,ntns
            OM=GetPrimitiveOperMat(j,N,m+1)
            call SumOperMats(HM,1.d0,OM,diags(j,m))
            call FlushOperMat(OM)
         ENDDO

!        Convert operator matrix to upper triangular form
         call Vec2SymPackMat(HM%mat,Hmat)

!        Solve the eigenvalue problem
         eigv=0.d0
         S=0.d0
         DO i=1,N
            S(i,i)=1.d0
         ENDDO
         call SolveGenEigval(eigv,S,Hmat,'N')

!        The smallest eigenvalue is added to the ZPVE
         ZPVE=ZPVE+eigv(1)

         DEALLOCATE(Hmat)
         call FlushOperMat(HM)
      ENDDO

      DEALLOCATE(eigv,S)

      end subroutine GetDiagonalZPVE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE HAMILOPT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
