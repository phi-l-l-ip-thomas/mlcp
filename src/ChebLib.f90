!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module CHEBLIB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP

      real, parameter :: PI=3.1415926535897932384626433832795028841971

      TYPE ChebObj
         integer :: n,intlevl
         character(LEN=3) :: bc
         real*8, allocatable, dimension(:) :: pt,acospt,mpt,fpt,coef,wt
         real*8, allocatable, dimension(:,:) :: deriv,intgl
         real*8  :: avdf(2),defint
         logical :: defintcalcd,egrid
      END TYPE ChebObj

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewChebObject(chb,n,st,fn,egrid,bc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes Chebyshev object by setting parameters and computing 
! Chebyshev points along with their mapped values

      implicit none
      TYPE (ChebObj) :: chb
      integer, intent(in) :: n
      real*8, intent(in)  :: st,fn
      logical, intent(in) :: egrid
      character(3), intent(in) :: bc
      integer :: i,nmod
      real*8, parameter :: bignr=1.d45

      chb%n=n
      chb%intlevl=0   ! Level of integration(>0) or differentiation(<0)
      chb%egrid=egrid ! .TRUE. for extrema grid, .FALSE. for nodes grid
      chb%bc=bc       ! boundary condition: 'fin'=finite, 
!                     'sem'=semi-infinite, 'inf'=infinite

!     Compute Chebyshev nodes or extrema
      IF (egrid) THEN
         nmod=n+1
         call ChebExtr(n,chb%pt)
      ELSE
         nmod=n
         call ChebNodes(n,chb%pt)
      ENDIF

!     Compute 'avdf' from the values of 'st' and 'fn'
      IF (bc.eq.'fin') THEN
         chb%avdf=averdif((/st,fn/))
      ELSE
         chb%avdf=(/st,fn/)
      ENDIF

      ALLOCATE(chb%mpt(nmod),chb%fpt(nmod))
      chb%fpt=0.d0
      DO i=1,nmod
         chb%mpt(i)=mapx(chb%pt(i),chb%avdf,bc,.FALSE.)
      ENDDO

!     Patch to avoid placing quadrature points at infinity for
!     semi-infinite or infinite b.c. cases
      IF ((bc.eq.'sem' .or. bc.eq.'inf') .and. egrid) THEN
         IF (nmod.lt.2) call AbortWithError( &
         "NewChebObject(): < 2 pts for bc='inf'/'sem', extrema grid")
         chb%mpt(1)=chb%mpt(2)+bignr
      ENDIF
      IF (bc.eq.'inf' .and. egrid) THEN
         IF (nmod.lt.3) call AbortWithError( &
         "NewChebObject(): < 3 pts for bc='inf', extrema grid")
         chb%mpt(nmod)=chb%mpt(nmod-1)-bignr
      ENDIF

      chb%defintcalcd=.FALSE.
      chb%defint=-1.d-300

      end subroutine NewChebObject

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CalcChebCoefs(chb)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Use to compute the Chebyshev coefs once the 'chb%fpts' array is filled

      implicit none
      TYPE (ChebObj) :: chb

      IF (chb%egrid) THEN
         call slowDCTe(chb%fpt,chb%coef,.TRUE.)
      ELSE
         call slowDCTn(chb%fpt,chb%coef,.TRUE.)
      ENDIF

      end subroutine CalcChebCoefs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebCalculus(chb,n,iconsts)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes derivatives, indefinite and definite integrals for Chebyshev-
! approximated function
! Set n=0 to compute the definite integral. This is computed using
! Chebyshev coefs if available; otherwise quadrature weights are 
! computed and the integral is done using Gaussian quadrature
! If n<0, the first |n| derivatives are computed;
! if n>0, the first n integrals are computed (which use the integration
! constants, if provided)

      implicit none
      TYPE (ChebObj) :: chb
      integer, intent(in) :: n
      real, allocatable, intent(inout) :: iconsts(:)
      real, allocatable :: tv(:),tv2(:)
      real*8  :: ic(1)
      integer :: i,j,nmod

      nmod=chb%n
      IF (chb%egrid) nmod=nmod+1

      IF (n.eq.0) THEN  ! Definite integral
         IF (ALLOCATED(chb%coef) .and. chb%bc.eq.'fin') THEN
            chb%defint=ChebDefIntgrl(chb%coef,chb%avdf,chb%egrid,chb%bc)
!            write(*,*)
!            write(*,*) 'Def. Int. from coefs.:', chb%defint
            chb%defintcalcd=.TRUE.
         ELSEIF (chb%egrid) THEN
            call acosExtr(chb%n,chb%acospt)
            call GetQuadWeights(chb%n,chb%acospt,chb%wt,chb%avdf,&
                                chb%egrid,chb%bc)
            chb%defint=ChebIntQuad(chb%fpt,chb%wt)
!            write(*,*) 'Def. Int. from quadr.:', chb%defint
            chb%defintcalcd=.TRUE.
         ELSE
            write(*,*) "ChebCalculus(): error for definite integral"
            write(*,*) " Implemented integration rules: use either of"
            write(*,*) " 1) integrate coefs + finite b.c. + any grid,  "
            write(*,*) " 2) Gauss quadrature + any b.c. + extrema grid "
            call AbortWithError('Error in ChebCalculus()')
         ENDIF

      ELSEIF (n.lt.0) THEN  ! derivative
         ALLOCATE(chb%deriv(nmod,abs(n)),tv(nmod))
         tv=chb%coef
         DO i=1,abs(n)
            call ChebDeriv(tv,tv2,1,chb%avdf,chb%bc)
            chb%deriv(:,i)=tv2(:)
            tv=tv2
            DEALLOCATE(tv2)
         ENDDO
         DEALLOCATE(tv)

      ELSE  ! indefinite integral
         ALLOCATE(chb%intgl(nmod+n,n))
         ic(1)=0.d0
!        Copy the integration constant if provided
         IF (ALLOCATED(iconsts)) ic(1)=iconsts(n)
         call ChebIndefIntgrl(chb%coef,tv,ic,chb%avdf,chb%bc)
         chb%intgl(1:nmod+1,1)=tv(1:nmod+1)
         DEALLOCATE(tv)
         DO i=2,n
            IF (ALLOCATED(iconsts)) ic(1)=iconsts(n-i+1)
            call ChebIndefIntgrl(chb%intgl(1:nmod+i-1,i-1),tv,ic,&
                                 chb%avdf,chb%bc)
            chb%intgl(1:nmod+i,i)=tv(1:nmod+i)
            DEALLOCATE(tv)
            DO j=1,i-1
               chb%intgl(nmod+i,j)=0.d0
            ENDDO
         ENDDO
      ENDIF

      end subroutine ChebCalculus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintChebInfo(chb)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints info about the Chebyshev object

      implicit none
      TYPE (ChebObj) :: chb
      integer :: i,j,ilen,nmod
      character*64 :: frmt

      write(*,'(/A)') '################################################'
      write(*,'(/X,A/)') 'Information about the Chebyshev object:'

      nmod=chb%n
      IF (chb%egrid) THEN
         nmod=nmod+1
         write(*,'(X,3A)') 'Boundary conditions = ',chb%bc,&
                            '; extrema grid'
      ELSE
         write(*,'(X,3A)') 'Boundary conditions = ',chb%bc,&
                            '; nodes grid'
      ENDIF

!     Print info about the points if the function has been evaluated
      IF (ALLOCATED(chb%fpt)) THEN
         write(*,'(/A)') 'Quadrature points and evaluations'
         write(*,'(3X,A,5X,A,7X,A,12X,A)') 'i','Cheb. points','x_i',&
               'f(x_i)'
         DO i=1,SIZE(chb%pt)
            write(*,'(I4,3(X,f16.8))') i-1,chb%pt(i),chb%mpt(i),&
            chb%fpt(i)
         ENDDO
      ENDIF

!     Print Chebyshev coefficients (and possibly those of the integrals
!     or derivatives if available
      IF (ALLOCATED(chb%coef)) THEN
         IF ((.not.ALLOCATED(chb%deriv)) .and. &
             (.not.ALLOCATED(chb%intgl))) THEN
            write(*,'(/A)') 'Chebyshev coefficients'
            write(*,'(3X,A,2X,A)') 'i','Coefficient'
            DO i=1,nmod
               write(*,'(I4,X,f16.8)') i-1,chb%coef(i)
            ENDDO
         ELSE
            IF (ALLOCATED(chb%deriv)) THEN
               ilen=SIZE(chb%deriv,2)
               write(*,'(/A)') 'Chebyshev function and derivative coefs'
               write(frmt,'(A,I0,A)') '(3X,A,8X,A,',ilen,'(8X,A,I0))'
               write(*,frmt) 'i','coef-f(x)',('coef-der',j,j=1,ilen)
               write(frmt,'(A,I0,A)') '(I4,',ilen+1,'(X,f16.8))'
               DO i=1,nmod
                  write(*,frmt) i-1,chb%coef(i),(chb%deriv(i,j),j=1,ilen)
               ENDDO
            ENDIF
            IF (ALLOCATED(chb%intgl)) THEN
               ilen=SIZE(chb%intgl,2)
               write(*,'(/A)') 'Chebyshev function and integral coefs'
               write(frmt,'(A,I0,A)') '(3X,A,8X,A,',ilen,'(8X,A,I0))'
               write(*,frmt) 'i','coef-f(x)',('coef-igl',j,j=1,ilen)
               write(frmt,'(A,I0,A)') '(I4,',ilen+1,'(X,f16.8))'
               DO i=1,nmod
                  write(*,frmt) i-1,chb%coef(i),(chb%intgl(i,j),j=1,ilen)
               ENDDO
!              Print extra lines for integral coefs since integral
!              coef array is longer than original array
               DO i=1,ilen
                  write(frmt,'(A,I0,A,I0,A)') '(I4,',17*i,'X,',&
                  ilen-1+1,'(X,f16.8))'
                  write(*,frmt) nmod+i-1,(chb%intgl(nmod+i,j),j=i,ilen)
               ENDDO
            ENDIF
         ENDIF
      ENDIF

      IF (ALLOCATED(chb%wt)) THEN
         write(*,'(/A)') 'Clenshaw-Curtis quadrature points and weights'
         write(*,'(3X,A,5X,A,7X,A,12X,A)') 'i',' arccos(x_i)','w_i',&
               'f(x_i)'
         DO i=1,nmod
            write(*,'(I4,3(X,f16.8))') i-1,chb%acospt(i),chb%wt(i),&
            chb%fpt(i)
         ENDDO
      ENDIF

      IF (chb%defintcalcd) THEN
         write(*,*)
         write(*,'(/A,f20.12)') 'Definite Integral:',chb%defint
      ENDIF

      write(*,'(/A/)') &
            '################################################'


      end subroutine PrintChebInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushChebObject(chb)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes Chebyshev object

      implicit none
      TYPE (ChebObj) :: chb

      IF (ALLOCATED(chb%pt)) DEALLOCATE(chb%pt)
      IF (ALLOCATED(chb%mpt)) DEALLOCATE(chb%mpt)
      IF (ALLOCATED(chb%fpt)) DEALLOCATE(chb%fpt)
      IF (ALLOCATED(chb%coef)) DEALLOCATE(chb%coef)
      IF (ALLOCATED(chb%acospt)) DEALLOCATE(chb%acospt)
      IF (ALLOCATED(chb%deriv)) DEALLOCATE(chb%deriv)
      IF (ALLOCATED(chb%intgl)) DEALLOCATE(chb%intgl)
      IF (ALLOCATED(chb%wt)) DEALLOCATE(chb%wt)
      chb%defintcalcd=.FALSE.

      end subroutine FlushChebObject

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Tcheb(n,x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Chebyshev polynomial of the first kind

      implicit none
      integer, intent(in) :: n
      real*8, intent(in)  :: x
      real*8 :: Tcheb

      IF (n.lt.0) THEN
         CALL AbortWithError("Tcheb(): n must be >= 0")
      ELSEIF (n.eq.0) THEN
         Tcheb=1.d0
      ELSEIF (n.eq.1) THEN
         Tcheb=x
      ELSEIF (abs(x).le.1.d0) THEN
         Tcheb=cos(n*acos(x))
      ELSEIF (x.lt.-1.d0) THEN
         Tcheb=cosh(n*acosh(-x))
         IF (mod(n,2).eq.1) Tcheb=-Tcheb
      ELSE
         Tcheb=cosh(n*acosh(x))
      ENDIF

      end function Tcheb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebNodes(n,cn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns Chebyshev nodes with indices 0...n-1

      implicit none
      integer, intent(in) :: n
      real*8, allocatable, intent(out) :: cn(:)
      integer :: i
      real*8  :: pN

      IF (n.lt.1) &
         CALL AbortWithError("ChebNodes(): n must be > 0")

      pN=PI/n

      ALLOCATE(cn(n))
      DO i=1,n
         cn(i)=cos(pN*(i-0.5d0))
      ENDDO

      end subroutine ChebNodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebExtr(n,ce)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns Chebyshev extrema with indices 0...n

      implicit none
      integer, intent(in) :: n
      real*8, allocatable, intent(out) :: ce(:)
      integer :: i
      real*8  :: pN

      IF (n.lt.0) &
         CALL AbortWithError("ChebExtr(): n must be >= 0")

      pN=PI/n

      ALLOCATE(ce(n+1))
      DO i=1,n+1
         ce(i)=cos(pN*(i-1.d0))
      ENDDO

      end subroutine ChebExtr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine acosNodes(n,acn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns acos of Chebyshev nodes with indices 0...n-1

      implicit none
      integer, intent(in) :: n
      real*8, allocatable, intent(out) :: acn(:)
      integer :: i
      real*8  :: pN

      IF (n.lt.1) &
         CALL AbortWithError("acosNodes(): n must be > 0")

      pN=PI/n

      ALLOCATE(acn(n))
      DO i=1,n
         acn(i)=pN*(i-0.5d0)
      ENDDO

      end subroutine acosNodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine acosExtr(n,ace)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns acos of Chebyshev extrema with indices 0...n

      implicit none
      integer, intent(in) :: n
      real*8, allocatable, intent(out) :: ace(:)
      integer :: i
      real*8  :: pN

      IF (n.lt.0) &
         CALL AbortWithError("acosExtr(): n must be >= 0")

      pN=PI/n

      ALLOCATE(ace(n+1))
      DO i=1,n+1
         ace(i)=pN*(i-1.d0)
      ENDDO

      end subroutine acosExtr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ChebCard(n,j,x,egrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns j-th (out of n) Chebyshev cardinal function, evaluated at x

      implicit none
      integer, intent(in) :: n,j
      real*8, intent(in)  :: x
      logical, intent(in) :: egrid
      real*8 :: ChebCard

      IF (egrid) THEN
         ChebCard=ChebCardEx(n,j,x)
      ELSE
         ChebCard=ChebCardNd(n,j,x)
      ENDIF

      end function ChebCard

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ChebCardNd(n,j,x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns j-th (out of n) Chebyshev node cardinal function, evaluated 
! at point x

      implicit none
      integer, intent(in) :: n,j
      real*8, intent(in)  :: x
      real*8, parameter :: tol=1.d-12
      real*8 :: xx,x1,ChebCardNd

      IF (n.lt.0) &
         CALL AbortWithError("ChebCardNd(): n must be >=0")
      IF (j.lt.0 .or. j.ge.n) &
         CALL AbortWithError("ChebCardNd(): j is out of range")

      xx=PI*(j+0.5d0)/n
      IF (abs(x).lt.1.d0) THEN
        x1=acos(x)
      ELSEIF (x.le.-1.d0) THEN
        x1=acosh(-x)
      ELSE
        x1=acosh(x)
      ENDIF

      IF (abs(x1-xx).lt.tol) THEN
         ChebCardNd=1.d0
      ELSEIF (xx.eq.0.d0) THEN
         IF (abs(x).lt.1.d0) THEN
            ChebCardNd=cos(n*x1)/(n*(cos(x1)-cos(xx)))
         ELSEIF (x.le.-1.d0) THEN
            ChebCardNd=cosh(n*x1)/(n*(-cosh(x1)-cos(xx)))
            IF (mod(n,2).eq.1) ChebCardNd=-ChebCardNd
         ELSE
            ChebCardNd=cosh(n*x1)/(n*(cosh(x1)-cos(xx)))
         ENDIF
      ELSE
         IF (abs(x).lt.1.d0) THEN
            ChebCardNd=cos(n*x1)*sin(xx)/(n*(cos(x1)-cos(xx))*sin(n*xx))
         ELSEIF (x.le.-1.d0) THEN
            ChebCardNd=cosh(n*x1)*sin(xx)/ &
            (n*(-cosh(x1)-cos(xx))*sin(n*xx))
            IF (mod(n,2).eq.1) ChebCardNd=-ChebCardNd
         ELSE
            ChebCardNd=cosh(n*x1)*sin(xx)/ &
            (n*(cosh(x1)-cos(xx))*sin(n*xx))
         ENDIF
      ENDIF

      end function ChebCardNd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ChebCardEx(n,j,x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns j-th (out of n) Chebyshev node cardinal function, evaluated 
! at point x

      implicit none
      integer, intent(in) :: n,j
      real*8, intent(in)  :: x
      real*8, parameter :: tol=1.d-12
      integer :: cj
      real*8  :: xj,xd,st,acx,ChebCardEx

      IF (n.lt.0) &
         CALL AbortWithError("ChebCardEx(): n must be >=0")
      IF (j.lt.0 .or. j.gt.n) &
         CALL AbortWithError("ChebCardEx(): j is out of range")

      xj=cos(PI*j/n)
      xd=x-xj

      IF (abs(xd).lt.tol) THEN
         ChebCardEx=1.d0
      ELSE
         IF (abs(x).lt.1.d0) THEN
            st=sqrt(1.d0-x*x)
         ELSEIF (x.le.-1.d0) THEN
            st=sqrt(x*x-1)
         ELSE
            st=-sqrt(x*x-1)
         ENDIF
         IF (x.ne.0.d0) THEN
            IF (abs(x).lt.1.d0) THEN
               acx=atan(st/x)
            ELSEIF (x.le.-1.d0) THEN
               acx=acosh(-x)
            ELSE
               acx=acosh(x)
            ENDIF
         ELSE
            acx=0.5d0*PI
         ENDIF
         cj=n
         IF (j.eq.0 .or. j.eq.n) cj=2*cj
         IF (mod(j,2).eq.0) cj=-cj
         IF (mod(n,2).eq.1 .and. x.lt.0) cj=-cj
         IF (abs(x).lt.1.d0) THEN
            ChebCardEx=st*sin(n*acx)/(cj*xd)
         ELSEIF (x.le.-1.d0) THEN
            ChebCardEx=st*sinh(n*acx)/(cj*xd)
         ELSE
            ChebCardEx=st*sinh(n*acx)/(cj*xd)
         ENDIF
      ENDIF

      end function ChebCardEx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine slowDCT(v,dct,forw,egrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Discrete cosine transform, slow version, wrapper for both nodes and
! extrema grids

      implicit none
      real*8, intent(in)  :: v(:)
      logical, intent(in) :: forw,egrid
      real*8, allocatable, intent(out) :: dct(:)

      IF (egrid) THEN
         call slowDCTe(v,dct,forw)
      ELSE
         call slowDCTn(v,dct,forw)
      ENDIF

      end subroutine slowDCT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine slowDCTe(v,dct,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Discrete cosine transform, extrema grid (slow version)

      implicit none
      real*8, intent(in)  :: v(:)
      logical, intent(in) :: forw
      real*8, allocatable, intent(out) :: dct(:)
      integer :: i,j,lv
      real*8  :: pl,plj,ve,vo,sm

      lv=SIZE(v)

      IF (lv.lt.2) &
         CALL AbortWithError('slowDCTe(): v must have length >= 2')

      pl=PI/(lv-1.d0)
      ve=0.5d0*(v(1)+v(lv))
      vo=0.5d0*(v(1)-v(lv))

      ALLOCATE(dct(lv))

      DO j=1,lv
         IF (mod(j-1,2).eq.0) THEN
            sm=ve
         ELSE
            sm=vo
         ENDIF
         plj=pl*(j-1)
         DO i=2,lv-1
            sm=sm+v(i)*cos(plj*(i-1))
         ENDDO
         dct(j)=sm
         IF (forw) dct(j)=dct(j)*2.d0/(lv-1.d0)
      ENDDO

      end subroutine slowDCTe

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine slowDCTn(v,dct,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Discrete cosine transform, nodes grid (slow version)

      implicit none
      real*8, intent(in)  :: v(:)
      logical, intent(in) :: forw
      real*8, allocatable, intent(out) :: dct(:)
      real*8, allocatable :: cnodes(:)
      logical, parameter  :: egrid=.FALSE.
      integer :: i,j,lv
      real*8  :: pl,sm,plj

      lv=SIZE(v)
      pl=PI/lv

      ALLOCATE(dct(lv))

      IF (forw) THEN
         DO j=1,lv
            sm=0.d0
            plj=pl*(j-1)
            DO i=1,lv
               sm=sm+v(i)*cos(plj*(i-0.5))
            ENDDO
            dct(j)=2*sm/lv
         ENDDO
      ELSE
         call ChebNodes(lv,cnodes)
         DO j=1,lv
            dct(j)=ClenshawRecur(v,cnodes(j),lv-1,egrid)
         ENDDO
         DEALLOCATE(cnodes)
      ENDIF

      end subroutine slowDCTn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine cheb2poly(v,avdf,forw,egrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Replaces Chebyshev coefs with scaled-and-shifted polynomial coefs
! See NR, 3rd ed. p.242-243

      implicit none
      real*8, intent(inout) :: v(:)
      real*8, intent(in)    :: avdf(2)
      logical, intent(in) :: forw,egrid
      real*8, allocatable :: w(:),w2(:)
      integer :: i,j,n
      real*8  :: sv

      n=SIZE(v)
      ALLOCATE(w(n),w2(n))
      w=0.d0
      w2=0.d0

      IF (forw) THEN

!        Convert Chebyshev coefs to polynomial coefs
         w(1)=v(n)
         IF (egrid) w(1)=0.5d0*w(1)

         DO i=n-1,2,-1
            DO j=n-i+2,2,-1
               sv=w(j)
               w(j)=2.d0*w(j-1)-w2(j)
               w2(j)=sv
            ENDDO
            sv=w(1)
            w(1)=v(i)-w2(1)
            w2(1)=sv
         ENDDO

         DO i=n,2,-1
            w(i)=w(i-1)-w2(i)
         ENDDO
         w(1)=0.5d0*v(1)-w2(1)

!        Map polynomial coefficients from [-1,1] --> [a,b]

!        Scale "width" of coefs to width of [a,b]
         sv=avdf(2)
         DO i=2,n
            w(i)=w(i)/sv
            sv=sv*avdf(2)
         ENDDO

!        Shift coefs to new position
         sv=avdf(1)
         DO i=1,n-1
            DO j=n-1,i,-1
               w(j)=w(j)-sv*w(j+1)
            ENDDO
         ENDDO

!        Replace v with w
         v=w
      ELSE
         call AbortWithError("cheb2poly() rev tfm not yet implemented")
      ENDIF

      DEALLOCATE(w,w2)

      end subroutine cheb2poly

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DCTmat(n,M,forw,egrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets the discrete cosine transform matrix, wrapper for both nodes and
! extrema grids

      implicit none
      integer, intent(in) :: n
      logical, intent(in) :: forw,egrid
      real*8, allocatable, intent(out) :: M(:,:)

      IF (egrid) THEN
          call DCTmate(n,M,forw)
      ELSE
          call DCTmatn(n,M,forw)
      ENDIF

      end subroutine DCTmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DCTmate(n,M,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix representation of discrete cosine transform, extrema grid

      implicit none
      integer, intent(in) :: n
      logical, intent(in) :: forw
      real*8, allocatable, intent(out) :: M(:,:)
      integer :: i,j
      real*8  :: pn,pkn

      IF (n.lt.1) CALL AbortWithError('DCTmate(): n < 1')

      pn=PI/n
      ALLOCATE(M(n+1,n+1))

!     Fill the DCT matrix
      DO i=1,n+1
         M(i,1)=0.5d0
         pkn=pn*(i-1)
         DO j=2,n
            M(i,j)=cos(pkn*(j-1))
         ENDDO
         M(i,n+1)=0.5d0
         IF (mod(i,2).eq.0) M(i,n+1)=-M(i,n+1)
      ENDDO

!     Forward transform only: multiply by 2/n
      IF (forw) THEN
         pn=2.d0/n
         M=pn*M
      ENDIF

      end subroutine DCTmate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DCTmatn(n,M,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix representation of discrete cosine transform, nodes grid

      implicit none
      integer, intent(in) :: n
      logical, intent(in) :: forw
      real*8, allocatable, intent(out) :: M(:,:)
      integer :: i,j
      real*8  :: pn,pkn,tn

      IF (n.lt.1) CALL AbortWithError('DCTmatn(): n < 1')

      pn=PI/n
      tn=2.d0/n
      ALLOCATE(M(n,n))

!     Forward transform
      IF (forw) THEN
         DO i=1,n
            pkn=pn*(i-1)
            DO j=1,n
               M(i,j)=tn*cos(pkn*(j-0.5d0))
            ENDDO
         ENDDO
!     Reverse transform
      ELSE
         DO i=1,n
            pkn=pn*(i-0.5d0)
            M(i,1)=0.5d0
            DO j=2,n
               M(i,j)=cos(pkn*(j-1.d0))
            ENDDO
         ENDDO
      ENDIF

      end subroutine DCTmatn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ClenshawRecur(cc,x,n,egrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Evaluates Chebyshev series (with coefficients in 'cc') up to order n
! at point x. 'ty' can be 'N' (nodes) or 'E' (extrema)

      implicit none
      real*8, intent(in)  :: cc(:)
      real*8, intent(in)  :: x
      integer, intent(in) :: n
      logical, intent(in) :: egrid
      real*8  :: tx,ov,nv,sv,ClenshawRecur
      integer :: k

      IF (n.ge.SIZE(cc)) &
         call AbortWithError("ClenshawRecur(): n > SIZE(v)")

      tx=2.d0*x
      ov=0.d0
      nv=0.d0

      DO k=n+1,2,-1
         sv=nv
         IF (k.eq.SIZE(cc) .and. egrid) THEN
            nv=0.5*cc(k)
         ELSE
            nv=tx*nv-ov+cc(k)
         ENDIF
         ov=sv
      ENDDO

      ClenshawRecur=x*nv-ov+0.5d0*cc(1)

      end function ClenshawRecur

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function averdif(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns average and difference formulas from endpoints

      implicit none
      real*8, intent(in) :: v(:)
      real*8  :: averdif(2)
      integer :: lv

      lv=SIZE(v)
      averdif(1)=0.5d0*(v(lv)+v(1))
      averdif(2)=0.5d0*(v(lv)-v(1))

      end function averdif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function mapx(x,avdf,bc,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Maps x to (forw=.True.) or from [-1,1] using parameters in avdf
! The string 'bc' specifies the boundary conditions:
! Finite interval: bc='fin': [a,b]<->[-1,1]
! Semi-infinite: bc='sem': [a,inf]<->[-1,1]
! Infinite: bc='inf' [-inf,inf]<->[-1,1]
! Periodic: bc='per' [a,b]<->[-1,1]
! Cosine:   bc='cos' [a,b]<->[-1,1]

      implicit none
      character(3), intent(in) :: bc
      real*8, intent(in)  :: x
      real*8, intent(in)  :: avdf(2)
      logical, intent(in) :: forw
      real*8  :: mapx,tmp
      integer :: i,n

!     cosine mapping (use higher order for sharper peaked distribution)
      n=2  ! cosine order

      IF (bc.ne.'fin' .and. bc.ne.'sem' .and. &
          bc.ne.'inf' .and. bc.ne.'per' .and. & 
          bc.ne.'cos') &
         call AbortWithError('mapx(): unrecognized b.c.')

      IF (forw) THEN
         IF (bc.eq.'fin') THEN
            mapx=(x-avdf(1))/avdf(2)
         ELSEIF (bc.eq.'sem') THEN
            tmp=x-avdf(1) ! avdf(1)=start of interval
            mapx=(tmp-avdf(2))/(tmp+avdf(2)) ! avdf(2)=width param.
         ELSEIF (bc.eq.'inf') THEN
            tmp=x-avdf(1) ! avdf(1)=center of interval
            mapx=tmp/sqrt(avdf(2)**2+tmp**2) ! avdf(2)=width param.
         ELSEIF (bc.eq.'per') THEN
            mapx=x
            tmp=abs(avdf(2)-avdf(1))
!           If x < a, add periodic interval until x >= a
            DO
               IF (mapx.ge.min(avdf(1),avdf(2))) EXIT
               mapx=mapx+tmp
            ENDDO
!           If x > a, subtract periodic interval until x <= b
            DO
               IF (mapx.le.max(avdf(1),avdf(2))) EXIT
               mapx=mapx-tmp
            ENDDO
            tmp=(avdf(2)+avdf(1))/2
            mapx=cos((mapx-tmp)*PI/(avdf(2)-avdf(1))+PI/2)
         ELSEIF (bc.eq.'cos') THEN
            mapx=(x-0.5d0*(avdf(1)+avdf(2)))*2.d0/(avdf(2)-avdf(1))
            DO i=1,n
               mapx=cos(0.5*PI*(mapx+1.d0))
            ENDDO
         ENDIF
      ELSE
         IF (bc.eq.'fin') THEN
            mapx=x*avdf(2)+avdf(1)
         ELSEIF (bc.eq.'sem') THEN
            mapx=avdf(2)*(1.d0+x)/(1.d0-x)+avdf(1)
         ELSEIF (bc.eq.'inf') THEN
            mapx=avdf(2)*x/sqrt(1.d0-x**2)+avdf(1)
         ELSEIF (bc.eq.'per') THEN
            tmp=(avdf(2)+avdf(1))/2
            mapx=(acos(x)-PI/2)*(avdf(2)-avdf(1))/PI+tmp
         ELSEIF (bc.eq.'cos') THEN
            mapx=x
            DO i=1,n
               mapx=(acos(mapx)/PI-0.5)*2
            ENDDO
            mapx=0.5d0*(mapx*(avdf(2)-avdf(1))+(avdf(1)+avdf(2)))
         ENDIF
      ENDIF

      end function mapx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ChebDefIntgrl(cc,avdf,egrid,bc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given Chebyshev coefficients in 'cc', this returns the definite
! integral Intgrl_start^finish f(x) dx

      implicit none
      real*8, intent(in)      :: cc(:),avdf(2)
      character*3, intent(in) :: bc
      logical, intent(in) :: egrid
      integer :: lv,j
      real*8  :: tmp,ChebDefIntgrl

      IF (bc.ne.'fin' .and. bc.ne.'sem' .and. bc.ne.'inf' .and. &
          bc.ne.'per' .and. bc.ne.'cos') &
         call AbortWithError('ChebDefIntgrl(): unrecognized b.c.')

      IF (bc.ne.'fin') THEN
         write(*,*) 'ChebDefIntgrl: only finite b.c. are implemented'
         call AbortWithError('Error in ChebDefIntgrl()')
      ENDIF

      lv=SIZE(cc)
      ChebDefIntgrl=0.5d0*cc(1)
      IF (lv.gt.2) THEN
         DO j=3,lv,2
            tmp=cc(j)/(1.d0-(j-1)**2)
!           Check this next line !!!
            IF (j.eq.lv .and. egrid) tmp=0.5*tmp
            ChebDefIntgrl=ChebDefIntgrl+tmp
         ENDDO
      ENDIF

!     Normalization (NEEDS TO BE CHANGED FOR ALL BUT 'fin' !!!)
      IF (bc.eq.'fin') THEN
         ChebDefIntgrl=ChebDefIntgrl*2*avdf(2)
      ELSEIF (bc.eq.'sem') THEN
         ChebDefIntgrl=ChebDefIntgrl*2*avdf(2)
         call AbortWithError('ChebDefIntgrl(): "sem" not implemented')
      ELSEIF (bc.eq.'inf') THEN
         ChebDefIntgrl=ChebDefIntgrl*2*avdf(2)
         call AbortWithError('ChebDefIntgrl(): "inf" not implemented')
      ELSEIF (bc.eq.'per') THEN
         ChebDefIntgrl=ChebDefIntgrl*2*avdf(2)
         call AbortWithError('ChebDefIntgrl(): "per" not implemented')
      ELSEIF (bc.eq.'cos') THEN
         ChebDefIntgrl=ChebDefIntgrl*2*avdf(2)
         call AbortWithError('ChebDefIntgrl(): "cos" not implemented')
      ENDIF

      end function ChebDefIntgrl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebIndefIntgrl(cc,cint,c0,avdf,bc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given Chebyshev coefficients in 'cc', this returns the coefficients
! of the n-th indefinite integral, Intgrl^n f(x) dx^n. The length of c0,
! which holds the constants of integration, determines the number of 
! integrations

      implicit none
      character(3), intent(in) :: bc
      real*8, allocatable, intent(out) :: cint(:)
      real*8, intent(in)  :: cc(:),c0(:),avdf(2)
      real*8, allocatable :: tv(:)
      integer :: n,ccs,i,j
      real*8  :: norm,sm,fac

      IF (bc.ne.'fin' .and. bc.ne.'sem' .and. &
          bc.ne.'inf' .and. bc.ne.'per' .and. &
          bc.ne.'cos') &
         call AbortWithError('ChebIndefIntgrl(): unrecognized b.c.')

      IF (bc.ne.'fin') THEN
         write(*,*) 'ChebIndefIntgrl: only finite b.c. are implemented'
         call AbortWithError('Error in ChebIndefIntgrl()')
      ENDIF

!     Order of integration
      n=SIZE(c0)
      norm=0.5*avdf(2)

      ALLOCATE(tv(SIZE(cc)))
      tv=cc

      DO i=1,n
         fac=1.d0
         sm=0.d0
         ccs=SIZE(cc)+i
         ALLOCATE(cint(ccs))
         DO j=2,ccs-2
            cint(j)=norm*(tv(j-1)-tv(j+1))/(j-1.d0)
!        Uncomment the following lines for the NR-version
!            sm=sm+fac*cint(j)
!            fac=-fac
         ENDDO
         cint(ccs-1)=norm*tv(ccs-2)/(ccs-2.d0)
         cint(ccs)=norm*tv(ccs-1)/(ccs-1.d0)
!        Uncomment the following lines for the NR-version
!         sm=sm+fac*(cint(ccs-1)-cint(ccs))
!         cint(1)=2*sm
         cint(1)=c0(n-i+1) ! constant of integration from input
         DEALLOCATE(tv)
         IF (i.lt.n) THEN  ! Update temp array
            ALLOCATE(tv(ccs))
            tv=cint
            DEALLOCATE(cint)
         ENDIF
      ENDDO

      end subroutine ChebIndefIntgrl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ChebDeriv(cc,cder,n,avdf,bc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given Chebyshev coefficients in 'cc', this returns the coefficients of
! the n-th derivative d^n/dx^n f(x)

      implicit none
      character(3), intent(in) :: bc
      real*8, allocatable, intent(out) :: cder(:)
      real*8, intent(in)  :: cc(:),avdf(2)
      integer, intent(in) :: n
      real*8, allocatable :: tv(:)
      integer :: lv,i,j

      IF (bc.ne.'fin' .and. bc.ne.'sem' .and. &
          bc.ne.'inf' .and. bc.ne.'per' .and. &
          bc.ne.'cos') &
         call AbortWithError('ChebDeriv(): unrecognized b.c.')

      IF (bc.ne.'fin') THEN
         write(*,*) 'ChebDeriv: only finite b.c. are implemented'
         call AbortWithError('Error in ChebDeriv()')
      ENDIF

      lv=SIZE(cc)
      ALLOCATE(tv(lv),cder(lv))
      tv=cc

      DO i=1,n
         cder(lv)=0.d0
         cder(lv-1)=2*(lv-1)*tv(lv)
         DO j=lv-2,1,-1
            cder(j)=cder(j+2)+2*j*tv(j+1)
         ENDDO
         DO j=1,lv
            cder(j)=cder(j)/avdf(2)
         ENDDO
         IF (i.lt.n) tv=cder
         lv=lv-1
      ENDDO

      DEALLOCATE(tv)

      end subroutine ChebDeriv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetQuadWeights(n,xs,wts,avdf,egrid,bc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get the Clenshaw-Curtis quadrature weights for boundary conditions:
! Finite interval: bc='fin': [a,b]<->[-1,1]
! Semi-infinite: bc='sem': [a,inf]<->[-1,1]
! Infinite: bc='inf' [-inf,inf]<->[-1,1]
! The array 'xs' is the full array generated by acosExtr(), where
! only the interior extrema are used (Fejer form)

      implicit none
      character(3), intent(in) :: bc
      logical, intent(in) :: egrid
      integer, intent(in) :: n
      real*8, intent(in)  :: avdf(2),xs(:)
      real*8, allocatable, intent(out) :: wts(:)
      integer :: i,j
      real*8  :: tn,sm

      IF (bc.ne.'fin' .and. bc.ne.'sem' .and. &
          bc.ne.'inf' .and. bc.ne.'per' .and. &
          bc.ne.'cos') &
         call AbortWithError('GetQuadWeights(): unrecognized b.c.')

      IF (.not.egrid) THEN
         write(*,*) "GetQuadWeights(): nodes grid not yet implemented"
         call AbortWithError("Error in GetQuadWeights()")
      ENDIF

      tn=2.d0*avdf(2)/n
      ALLOCATE(wts(n+1))

!     Compute the prefactor and (dummy) endpoint weights
      IF (bc.eq.'sem') THEN
         tn=2*tn
      ELSEIF (bc.eq.'inf') THEN
         tn=0.5d0*PI*tn
      ENDIF
      wts(1)=0.d0
      wts(n+1)=0.d0

      IF (bc.eq.'per') &
         call AbortWithError('GetQuadWeights(): "per" not implemented')
      IF (bc.eq.'cos') &
         call AbortWithError('GetQuadWeights(): "cos" not implemented')

!     Compute the interior weights
      DO i=2,n
         IF (bc.ne.'inf') THEN  ! finite and semi-inf. b.c. cases
            sm=0.d0
            DO j=1,n-1
               sm=sm+sin(j*xs(i))*(1.d0-cos(j*PI))/j
            ENDDO
            wts(i)=sin(xs(i))*tn*sm
            IF (bc.eq.'sem') THEN
               wts(i)=wts(i)/((1.d0-cos(xs(i)))**2)
            ENDIF
         ELSE  ! infinite b.c. case
            wts(i)=tn/((sin(xs(i)))**2)
         ENDIF
      ENDDO

      end subroutine GetQuadWeights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ChebIntQuad(pts,wts)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real*8, intent(in) :: pts(:),wts(:)
      real*8  :: ChebIntQuad
      integer :: i,n

      n=SIZE(wts)
      IF (SIZE(pts).ne.n) &
         call AbortWithError("ChebIntQuad(): size mismatch pts/wts")

      ChebIntQuad=0.d0
      DO i=1,n
         ChebIntQuad=ChebIntQuad+pts(i)*wts(i)
      ENDDO

      end function ChebIntQuad

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TestChebLib()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests the Chebyshev library

      implicit none
      TYPE (ChebObj) :: chbn,chbe
      real*8, allocatable :: xset(:),yset(:),xmap(:),ychbe(:),ychbn(:)
      real*8, allocatable :: iconsts(:)
      character(3) :: bc
      logical :: egrid
      real*8  :: avdf(2)
      integer :: i,j,npts,nord
      real*8  :: xstart,xend,p1,p2

      write(*,'(/A/)') '***********************************************'
      write(*,*) 'Chebyshev library test'

! #####################################################################
!     First, test Tcheb, ChebCardEx, and ChebCardNd
! #####################################################################

      npts=301
      xstart=-1.5
      xend=1.5
      write(*,'(/A)') 'Testing Tcheb...'
      write(*,*) 'npts = ',npts,'; xstart = ',xstart,'; xend = ',xend

      ALLOCATE(xset(npts),yset(npts))
      DO i=1,npts
         xset(i)=xstart+(i-1)*(xend-xstart)/(npts-1.d0)
      ENDDO

      write(*,*) 'Tcheb: n = 0...5: test_Tcheb.dat'

      OPEN(40, FILE='test_Tcheb.dat', STATUS="unknown")
      DO i=1,npts
         write(40,*) xset(i),(Tcheb(j-1,xset(i)) ,j=1,6)
      ENDDO
      CLOSE(40)

      write(*,'(/A)') 'Testing ChebCardNd: test_ChebCardNd.dat'
      write(*,*) 'Over range: [',xstart,',',xend,']...'

      OPEN(42, FILE='test_ChebCardNd.dat', STATUS="unknown")
      DO nord=1,6
         call ChebNodes(nord,xmap)
         DO i=1,npts
            write(42,*) xset(i),(ChebCardNd(nord,j-1,xset(i)),j=1,nord)
         ENDDO
         write(42,*) 
         DO i=1,nord
            write(42,*) xmap(i),(ChebCardNd(nord,j-1,xmap(i)),j=1,nord)
         ENDDO
         write(42,*)
         DEALLOCATE(xmap)
      ENDDO
      CLOSE(42)

      write(*,'(/A)') 'Testing ChebCardEx: test_ChebCardEx.dat'
      write(*,*) 'Over range: [',xstart,',',xend,']...'

      OPEN(44, FILE='test_ChebCardEx.dat', STATUS="unknown")
      DO nord=1,6
         call ChebExtr(nord,xmap)
         DO i=1,npts
            write(44,*) xset(i),(ChebCardEx(nord,j-1,xset(i)),j=1,nord+1)
         ENDDO
         write(44,*)
         DO i=1,nord+1
            write(44,*) xmap(i),(ChebCardEx(nord,j-1,xmap(i)),j=1,nord+1)
         ENDDO
         DEALLOCATE(xmap)
         write(44,*)
      ENDDO
      CLOSE(44)
      DEALLOCATE(xset,yset)

! #####################################################################
!     Test Chebyshev approx subroutines for finite interval
! #####################################################################

!     Set parameters
      npts=301
      xstart=-13.
      xend=13.
      nord=256
      bc='inf'
!     Parameters for boundary conditions
      IF (bc.eq.'fin') THEN
         p1=xstart
         p2=xend
      ELSEIF (bc.eq.'sem') THEN
         p1=xstart
         p2=1.d0
      ELSEIF (bc.eq.'inf') THEN
         p1=0.5d0*(xend+xstart)
         p2=1.d0
      ENDIF

      write(*,'(/A,A)') 'Testing Cheb approx: bc = ',bc
      write(*,*) 'npts = ',npts,'; xstart = ',xstart,'; xend = ',xend

      ALLOCATE(xset(npts),xmap(npts),yset(npts),ychbe(npts),ychbn(npts))
      DO i=1,npts  ! Eval fxn at series of points as a reference
         xset(i)=xstart+(i-1)*(xend-xstart)/(npts-1.d0)
         yset(i)=ftest(xset(i))
      ENDDO

!     Extrema grid
      egrid=.TRUE.
      call NewChebObject(chbe,nord,p1,p2,egrid,bc)
      DO i=1,nord+1
         chbe%fpt(i)=ftest(chbe%mpt(i))
      ENDDO
      call CalcChebCoefs(chbe)
      call ChebCalculus(chbe,0,iconsts)
      call ChebCalculus(chbe,-2,iconsts)
      call ChebCalculus(chbe,2,iconsts)
      call PrintChebInfo(chbe)

!     Nodes grid
      egrid=.FALSE.
      call NewChebObject(chbn,nord,p1,p2,egrid,bc)
      DO i=1,nord
         chbn%fpt(i)=ftest(chbn%mpt(i))
      ENDDO
      call CalcChebCoefs(chbn)
      call ChebCalculus(chbn,0,iconsts)
      call ChebCalculus(chbn,-2,iconsts)
      call ChebCalculus(chbn,2,iconsts)
      call PrintChebInfo(chbn)

      write(*,*) 'Reference function, extrema-approx, nodes-approx'
      write(*,'(5(11X,A))') '   x_i','mapx_i','   y_i','ex y_i','nd y_i'
      DO i=1,npts  ! Print reference set, node/extrema approximations
         xmap(i)=mapx(xset(i),chbe%avdf,bc,.TRUE.)
         egrid=.TRUE.
         ychbe(i)=ClenshawRecur(chbe%coef,xmap(i),nord,egrid)
         egrid=.FALSE.
         ychbn(i)=ClenshawRecur(chbn%coef,xmap(i),nord-1,egrid)
         write(*,'(5(X,f16.8))') xset(i),xmap(i),yset(i),ychbe(i),ychbn(i)
      ENDDO

      call FlushChebObject(chbe)
      call FlushChebObject(chbn)
      DEALLOCATE(xset,xmap,yset,ychbe,ychbn)

      write(*,'(/A/)') '***********************************************'

      end subroutine TestChebLib

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ftest(x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real*8, intent(in) :: x
      real*8 :: ftest
      integer :: i

!      ftest=(x-1.d0)*x**2+1.d0
!      ftest=gaussian(x,2.d0,PI,5.d0)
!      ftest=gaussian(x,-2.d0*PI*x,PI,5.d0)
!      ftest=HObasisfxn(3,x)
!      ftest=HObasisfxn(2,x)**2
!      ftest=morse1D(x)
      ftest=HObasisfxn(40,x)*x**3

!      ftest=0.d0
!      DO i=1,8
!         ftest=ftest+i*Tcheb(i-1,x)
!      ENDDO

      end function ftest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function my1Dfunction(typ,x,ip1,ip2,rp1,rp2,rp3,rp4,rp5)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1D function driver. Call with integer parameters ip1, ip2, and 
! real parameters rp1, rp2, rp3, rp4, rp5

      implicit none
      integer, intent(in) :: typ,ip1,ip2
      real*8, intent(in)  :: x,rp1,rp2,rp3,rp4,rp5
      real*8 :: my1Dfunction

      IF (typ.eq.0) THEN
!        q^ip1
         my1Dfunction=x**ip1
      ELSEIF (typ.eq.1) THEN
         my1Dfunction=tanh(rp1*x)**ip1
      ELSEIF (typ.eq.2) THEN
         my1Dfunction=morse1D(x,ip1,1.d0,rp1,0.d0)
      ELSEIF (typ.eq.3) THEN
         my1Dfunction=sqgauss(x,ip1,rp1)
      ELSE
         call AbortWithError('my1Dfunction(): unrecognized function')
      ENDIF

      end function my1Dfunction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function sqgauss(x,pow,a)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1D power of the square-root of a Gaussian, symmetrized to have odd
! symmetry for pow=1

      implicit none
      integer, intent(in) :: pow
      real*8, intent(in)  :: x,a
      real*8 :: sqgauss

      IF (mod(pow,2).ne.0) &
         call AbortWithError('sqgauss(): odd value of pow')
      sqgauss=(1.d0-exp(-a*x**2))**(pow/2)

      end function sqgauss

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function morse1D(x,pow,b,a,c)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1D Morse function

      implicit none
      integer, intent(in) :: pow
      real*8, intent(in)  :: x,b,a,c
      real*8 :: morse1D

      morse1D=b*(1.d0-exp(-a*(x-c)))**pow

      end function morse1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HObasisseries(x,HOx)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills array HOx with harmonic oscillator basis functions evaluated
! at x, orders 0...n-1, where n is the length of HOx

      implicit none
      real*8, intent(inout) :: HOx(:)
      real*8, intent(in)    :: x
      real*8  :: fac
      integer :: i,n

      n=SIZE(HOx)
      fac=gaussian(x,1.d0,0.5d0,0.d0)/PI**(0.25d0)

!     n = 0 function
      HOx(1)=fac

      IF (n.gt.1) THEN
!        n = 1 function
         HOx(2)=fac*sqrt(2.d0)*x
         IF (n.gt.2) THEN
!           Recursion for n > 1
            DO i=3,n
               HOx(i)=sqrt(2.d0/(i-1))*&
                      (x*HOx(i-1)-sqrt((i-2)/2.d0)*HOx(i-2))
            ENDDO
         ENDIF
      ENDIF

      end subroutine HObasisseries

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function HObasisfxn(n,x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Harmonic oscillator basis function

      implicit none
      integer, intent(in) :: n
      real*8, intent(in)  :: x
      real*8  :: fac,rov,ov,HObasisfxn
      integer :: i

      fac=gaussian(x,1.d0,0.5d0,0.d0)/PI**(0.25d0)

      IF (n.lt.0) THEN
         call AbortWithError("HObasisfxn(): n must be >= 0")
      ELSEIF (n.eq.0) THEN
         HObasisfxn=fac
      ELSEIF (n.eq.1) THEN
         HObasisfxn=fac*sqrt(2.d0)*x
      ELSE
         rov=fac
         ov=fac*sqrt(2.d0)*x
         DO i=1,n-1
            HObasisfxn=sqrt(2.d0/(i+1))*(x*ov-sqrt(i/2.d0)*rov)
            rov=ov
            ov=HObasisfxn
!           If recursion is guaranteed to give zero, exit early!
            IF (abs(ov).eq.0.d0 .and. abs(rov).eq.0.d0) EXIT
         ENDDO
      ENDIF

      end function HObasisfxn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gaussian(x,b,a,c)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gaussian function

      implicit none
      real*8, intent(in) :: b,a,c,x
      real*8 :: gaussian

      gaussian=b*exp(-a*(x-c)**2)

      end function gaussian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function hermite(n,x)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Hermite polynomial

      implicit none
      integer, intent(in) :: n
      real*8, intent(in)  :: x
      integer :: i
      real*8  :: tx,hermite,rov,ov

      IF (n.eq.0) THEN
         hermite=1.d0
      ELSEIF (n.eq.1) THEN
         hermite=2*x
      ELSE
         tx=2*x
         rov=1.d0
         ov=tx
         DO i=1,n-1
            hermite=tx*ov-2*i*rov
            rov=ov
            ov=hermite
         ENDDO
      ENDIF

      end function hermite

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module CHEBLIB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
