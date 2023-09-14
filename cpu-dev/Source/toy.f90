!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module TOY

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP-matrix diagonalization, toy problem

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE RESTART
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE BLOCKUTILS
      USE LINSOLVER
      USE ALSDRVR
      USE CPMATH
      USE HGORTHO
      USE HG

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine toyproblem()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Toy problem driver

      implicit none
      TYPE (CP) :: T

      write(*,*)
      write(*,*) '##############################'
      write(*,*) '###       toy problem      ###'
      write(*,*) '##############################'

!      T=gentoy1()
!      T=gentoy2()
!      T=gentoy2a()
      T=gentoy3()
!      T=gentoy4()
!      T=gentoy5()
!      T=gentoy6()

      call processtoy(T)
      call FlushCP(T)

      call AbortWithError('toyproblem(): play time is over!')

      end subroutine toyproblem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy1() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #1
! rk = 2, ndof = 2, rows = [2,2], cols = [2,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=2
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=2
      cols(:)=2
      sym(:)=.FALSE.

      F=NewCP(rk,rows,cols,sym)

      F%coef(:)=1.d0
      F%base(:,1)=(/2.d0,1.d0,1.d0,5.d0,-1.d0,-2.d0,-3.d0,-4.d0/)
      F%base(:,2)=(/2.d0,-0.5d0,0.25d0,2.5d0,-8.d0,5.d0,2.d0,-6.d0/)

      deallocate(rows,cols,sym)

      end function gentoy1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy2() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #2
! rk = 2, ndof = 2, rows = [3,2], cols = [3,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=2
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=(/3,2/)
      cols(:)=(/3,2/)
      sym(:)=.FALSE.

      F=NewCP(rk,rows,cols,sym)

      F%coef(:)=1.d0
      F%base(:,1)=(/2.0,1.0,3.5,1.0,5.0,-0.15,3.4,-2.1,8.0,&
                   -1.d0,-2.d0,-3.d0,-4.d0/)
      F%base(:,2)=(/2.0,-0.5,0.5,0.25,2.5,-5.0,-1.3,2.7,-3.3,&
                   -8.d0,5.d0,2.d0,-6.d0/)

      deallocate(rows,cols,sym)

      end function gentoy2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy2a() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #2, modes swapped
! rk = 2, ndof = 2, rows = [2,3], cols = [2,3]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=2
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=(/2,3/)
      cols(:)=(/2,3/)
      sym(:)=.FALSE.

      F=NewCP(rk,rows,cols,sym)

      F%coef(:)=1.d0
      F%base(:,1)=(/-1.d0,-2.d0,-3.d0,-4.d0,&
                    2.0,1.0,3.5,1.0,5.0,-0.15,3.4,-2.1,8.0/)
      F%base(:,2)=(/-8.d0,5.d0,2.d0,-6.d0,&
                    2.0,-0.5,0.5,0.25,2.5,-5.0,-1.3,2.7,-3.3/)

      deallocate(rows,cols,sym)

      end function gentoy2a

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy3() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #3
! rk = 2, ndof = 3, rows = [2,2,2], cols = [2,2,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=3
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=2
      cols(:)=2
      sym(:)=.FALSE.

      F=NewCP(rk,rows,cols,sym)

      F%coef(:)=1.d0
      F%base(:,1)=(/2.d0,1.d0,1.d0,5.d0,&
                   -1.d0,-2.d0,-3.d0,-4.d0,&
                   1.5,-0.8,2.2,3.4/)
      F%base(:,2)=(/2.d0,-0.5d0,0.25d0,2.5d0,&
                   -8.d0,5.d0,2.d0,-6.d0,&
                   4.0,2.1,-3.7,1.1/)

      deallocate(rows,cols,sym)

      end function gentoy3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy4() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #4 (toggles ALS update-instead-of-build)
! rk = 2, ndof = 4, rows = [2,2,2,2], cols = [2,2,2,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=4
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=2
      cols(:)=2
      sym(:)=.FALSE.

      F=NewCP(rk,rows,cols,sym)

      F%coef(:)=1.d0
      F%base(:,1)=(/2.d0,1.d0,1.d0,5.d0,&
                   -1.d0,-2.d0,-3.d0,-4.d0,&
                   1.5,-0.8,2.2,3.4,&
                   -0.15,3.4,-2.1,8.0/)
      F%base(:,2)=(/2.d0,-0.5d0,0.25d0,2.5d0,&
                   -8.d0,5.d0,2.d0,-6.d0,&
                   4.0,2.1,-3.7,1.1,&
                   -5.0,-1.3,2.7,-3.3/)

      deallocate(rows,cols,sym)

      end function gentoy4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy5() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #5, random vals
! rk = X, ndof = 3, rows = [2,2,2], cols = [2,2,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=3
      integer :: rk

!     Generate the toy problem matrix
      rk=8

      allocate(rows(n),cols(n),sym(n))
      rows(:)=2
      cols(:)=2
      sym(:)=.FALSE.

      F=RandomCP(rk,rows,cols,sym)

      deallocate(rows,cols,sym)

      end function gentoy5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function gentoy6() result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates toy problem #5, random vals
! rk = X, ndof = 3, rows = [2,2,2], cols = [2,2,2]

      implicit none
      TYPE (CP) :: F
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer, parameter   :: n=3
      integer :: rk

!     Generate the toy problem matrix
      rk=2

      allocate(rows(n),cols(n),sym(n))
      rows(:)=(/4,4,3/)
      cols(:)=(/4,4,3/)
      sym(:)=.FALSE.

      F=RandomCP(rk,rows,cols,sym)

      deallocate(rows,cols,sym)

      end function gentoy6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine processtoy(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Plays with toy problem

      implicit none
      TYPE (CP), INTENT(IN) :: F
      TYPE (CP) :: Fn,Fs,Fi,I
      real*8, allocatable  :: M(:,:),Mn(:,:),Ms(:,:),Qmn(:,:),Mi(:,:)
      integer, parameter   :: nals=14
      integer :: solstat
      real*8 :: norm

      write(*,*)
      write(*,*) 'Initial F:'
      call PrintCPmat(F)

!     Convert to matrix
      M=CP2Matrix(F)

      write(*,*)
      write(*,*) 'F, as a matrix:'
      call PrintMatrix(M)

!!!
      write(*,*)
      write(*,*) 'F-inverse via Solver:'
      I=IdentityCPMatrix(F%rows,F%cols,F%sym)
      Fi=CopyCP(F)
      call AugmentVWithRandomN(Fi,6)
      call CPMatrixTranspose(Fi)
      solstat=ALS_solve(F,Fi,I,2*nals)
      call FlushCP(I)
      call PrintCPmat(Fi)

!     Convert to matrix
      Mi=CP2Matrix(Fi)
      write(*,*)
      write(*,*) 'Fi, as a matrix:'
      call PrintMatrix(Mi)
      call MatrixMult(Mi,.FALSE.,M,.FALSE.,Ms)
      write(*,*)
      write(*,*) 'Fi*F = I?:'
      call PrintMatrix(Ms)
!      call AbortWithError('done inverting')
!!!

      Fn=CopyCP(F)
      call AugmentVWithRandom(Fn,15)
!!!
!      I=IdentityCPMatrix(F%rows,F%cols,F%sym)
!      call HGintertwine(I,Fn,10,10,0,0.d0)
!      call AbortWithError('done HG on toy()')
!!!
      call CPMatNormA(Fn,nals)
      Mn=CP2Matrix(Fn)

      write(*,*)
      write(*,*) 'normalized F (Fn):'
      call PrintCPmat(Fn)
      write(*,*)
      write(*,*) 'normalized F, as a matrix (Mn):'
      call PrintMatrix(Mn)
      write(*,*)
      write(*,*) 'norms of Fs vecs:'
      call ShowNorm(Fn)
      write(*,*)
      write(*,*) 'rank-1 overlap matrix:'
      call ShowRank1UTU(Fn)

!      write(*,*) 'rank-1 overlap matrix:'
!      call ShowRank1UTU(F)

!      Mn=CP2Matrix(F) !!!
      call QRdecomp(Qmn,Mn) ! Mn overwritten with R
      write(*,*) 'QR of Mn: Q'
      call PrintMatrix(Qmn)
      write(*,*) 'QR of Mn: R'
      call PrintMatrix(Mn)
      deallocate(Mn)

      Fs=CopyCP(Fn)
      Ms=CP2Matrix(Fs)
      write(*,*) 'Try to orthog Fn using OrthogCPmat():'
      call OrthogCPmat(Fn,6)
      write(*,*) 'Does Q^T*M = R (in CP)?'
      call ShowRank1R(Fs,Fn)
      write(*,*)
      write(*,*) 'Mn (Q) as a matrix:'
      Mn=CP2Matrix(Fn)
      call PrintMatrix(Mn)

      write(*,*) 'Does Q^T*M = R (matrix)?'
      call MatrixMult(Mn,.TRUE.,Ms,.FALSE.)
      call PrintMatrix(Mn)
      call distancefromR(Mn)

      call FlushCP(Fn)
      deallocate(M,Mn,Ms,Qmn)

      end subroutine processtoy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine distancefromR(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real*8, intent(in) :: M(:,:)
      real*8 :: lower,upper,ratio
      integer :: i,j,r,c

      r=SIZE(M,1)
      c=SIZE(M,2)

      upper=0.d0
      lower=0.d0

      DO j=1,c
         DO i=1,j
            upper=upper+M(i,j)**2
         ENDDO
         DO i=j+1,r
            lower=lower+M(i,j)**2
         ENDDO
      ENDDO
      ratio=sqrt(abs(lower/(upper+lower)))

      write(*,*) 'sqrt(l^2/(u^2+l^2)) = ',ratio

      end subroutine distancefromR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module TOY

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
