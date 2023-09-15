!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GPUINTERTWINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for the CUDA GPU intertwining code

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE SEPDREPN
      USE HAMILSETUP
      USE MODVECVEC
      USE ISO_C_BINDING

      implicit none
      real*8, private  :: gpumvp_time=0.d0,gpuals_time=0.d0
      logical, private :: GPUMVP_SETUP=.FALSE.,GPUALS_SETUP=.FALSE.

      INTERFACE GPU_ALS
         MODULE PROCEDURE GPU_ALS_AA,GPU_ALS_AB
      END INTERFACE GPU_ALS

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupGPUmvp()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      gpumvp_time=0.d0
      gpuals_time=0.d0
      GPUMVP_SETUP=.TRUE.
      GPUALS_SETUP=.TRUE.

      end subroutine SetupGPUmvp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeGPUmvp()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. GPUMVP_SETUP) call SetupGPUmvp()

!     Set up the module if it was not set up already
      GPUMVP_SETUP = .FALSE.
      GPUALS_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total mat-vec time on GPU (GPUMVP)(s)',&
                            gpumvp_time
      write(*,'(X,A,X,f20.3)') 'Total     ALS time on GPU (GPUALS)(s)',&
                            gpuals_time

      end subroutine DisposeGPUmvp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Get_GPU_Memory(Q,HL,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for interfacing CUDA intertwining algorithm

      implicit none
      TYPE (CPpar), INTENT(IN)   :: cpp
      TYPE (HopList), INTENT(IN) :: HL
      TYPE (CPvec), INTENT(IN)   :: Q(:)
      real*8  :: t1,t2

      IF (.NOT. GPUMVP_SETUP) call SetupGPUmvp()

      call CPU_TIME(t1)

      call IntwDeviceMem(cpp%psirank,SIZE(Q(1)%nbas),SIZE(Q),&
                         SIZE(HL%opid),SIZE(HL%mptr),SIZE(HL%mat),&
                         cpp%ncpu,cpp%update,Q(1)%nbas)

      call CPU_TIME(t2)
      gpumvp_time=gpumvp_time+t2-t1

      end subroutine Get_GPU_Memory

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GPU_Intertwine(Q,HL,eigv,cpp,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for interfacing CUDA intertwining algorithm

      implicit none
      TYPE (CPpar), INTENT(IN)    :: cpp
      TYPE (HopList), INTENT(IN)  :: HL
      TYPE (CPvec), INTENT(INOUT) :: Q(:)
      real*8, intent(inout) :: eigv(:)
      integer, intent(in)   :: ishift
      real*8, intent(in)    :: Eshift
      integer :: j,nbloc,rk,rst
      real*8  :: t1,t2
#ifdef USE_DOUBLES
      real*8, allocatable :: hlcoef(:),hlmat(:),eig(:),qcoef(:)
      real*8, allocatable :: qbase(:,:)
      real*8 :: Esft
      real*8, parameter  :: noise=1.d-14
      integer, parameter :: PREC=8
#else
      real*4, allocatable :: hlcoef(:),hlmat(:),eig(:),qcoef(:)
      real*4, allocatable :: qbase(:,:)
      real*4 :: Esft
      real*8, parameter  :: noise=1.d-7
      integer, parameter :: PREC=4
#endif
      integer, allocatable :: nbas(:)
      integer :: qlen

      qlen=SIZE(Q(1)%base,1)
      ALLOCATE(nbas(SIZE(Q(1)%nbas)))
      nbas(:)=Q(1)%nbas(:)

      IF (.NOT. GPUMVP_SETUP) call SetupGPUmvp()

      call CPU_TIME(t1)

!     Initialize values
      nbloc=SIZE(Q)
      rk=cpp%psirank

!     Pack the CP-block, operator, eigenvalue arrays into the correct 
!     precision arrays to pass to the C code (TYPE CPvec is not
!     interoperable with C, and for single precision runs the data
!     must be extracted from the double-precision Fortran type)
      ALLOCATE(hlcoef(SIZE(HL%coef)),hlmat(SIZE(HL%mat)))
      ALLOCATE(eig(SIZE(eigv)))
      ALLOCATE(qcoef(rk*nbloc),qbase(qlen,rk*nbloc))
      hlcoef(:)=REAL(HL%coef(:),PREC)
      hlmat(:)=REAL(HL%mat(:),PREC)
      eig(:)=REAL(eigv(:),PREC)
      rst=0
      DO j=1,nbloc
         qcoef(rst+1:rst+rk)=REAL(Q(j)%coef(:),PREC)+noise
         qbase(:,rst+1:rst+rk)=REAL(Q(j)%base(:,:),PREC)
         call FlushCPvec(Q(j))  ! Free memory
         rst=rst+rk
      ENDDO
      Esft=REAL(Eshift,PREC)

!     Make sure the ranks of all vectors are the same as required by
!     the CUDA code
      DO j=1,nbloc
         IF (SIZE(Q(j)%coef).ne.rk) THEN
            write(*,*) 'psirank = ',rk,'; for j = ',j,&
                       ', rank(Q(j)) = ',SIZE(Q(j)%coef)
            call AbortWithError('GPU_Intertwine(): bad rank')
         ENDIF
      ENDDO

!     Call the C++ intertwining module with the component arrays in
!     TYPE CP-vec and TYPE HopList Hamiltonian to the C++ module.
      call IntertwineCPP(qbase,qcoef,nbas,&
                         rk,SIZE(nbas),nbloc,&
                         HL%opid,HL%mptr,hlcoef,hlmat,&
                         SIZE(HL%opid),SIZE(HL%mptr),SIZE(hlmat),&
                         eig,cpp%npow,cpp%ncpu,cpp%psinals,cpp%update,&
                         ishift,Esft)

      DEALLOCATE(hlcoef,hlmat)
      eigv(:)=REAL(eig(:),8)
      DEALLOCATE(eig)

!     Repack the block from IntertwineCPP() into Q
      rst=0
      DO j=1,nbloc
         call NewCPvec(Q(j),nbas,rk)
         Q(j)%coef(:)=REAL(qcoef(rst+1:rst+rk),8)
         Q(j)%base(:,:)=REAL(qbase(:,rst+1:rst+rk),8)
         rst=rst+rk
      ENDDO

      DEALLOCATE(nbas,qcoef,qbase)

      call CPU_TIME(t2)
      gpumvp_time=gpumvp_time+t2-t1

      end subroutine GPU_Intertwine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GPU_ALS_AA(G,rk,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for interfacing CUDA ALS algorithm. Compute F which
! approximates G to minimize ||F-G||_2

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: G
      TYPE (CPvec) :: F
      integer, intent(in)  :: rk,nals

      IF ((rk.gt.0) .and. (rk.lt.SIZE(G%coef)) .and. (nals.gt.0)) THEN
         call GetRandomCPvec(F,G%nbas,rk)
         call GPU_ALS_AB(F,G,nals)
         call ReplaceVwithW(G,F)
      ENDIF

      end subroutine GPU_ALS_AA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GPU_ALS_AB(F,G,nals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Master subroutine for interfacing CUDA ALS algorithm. Compute F which
! approximates G to minimize ||F-G||_2

      implicit none
      TYPE (CPvec), INTENT(IN)    :: G
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in)  :: nals
      integer, allocatable :: nbas(:)
      integer :: D,nD,rF,rG,rGaug,mult
      real*8  :: t1,t2
#ifdef USE_DOUBLES
      real*8, allocatable :: coef(:)
      real*8, allocatable :: base(:,:)
      real*8, parameter   :: noise=1.d-14
      integer, parameter  :: PREC=8
#else
      real*4, allocatable :: coef(:)
      real*4, allocatable :: base(:,:)
      real*8, parameter   :: noise=1.d-7
      integer, parameter  :: PREC=4
#endif

      IF (nals.lt.1) RETURN

      IF (.NOT. GPUMVP_SETUP) call SetupGPUmvp()

      call CPU_TIME(t1)

!     Initialize values
      D=SIZE(F%nbas)
      nD=SIZE(F%base,1)
      rF=SIZE(F%coef)

!     Excess rank of G to work with GPU code
      mult=(SIZE(G%coef)+rF-1)/rF
      rGaug=mult*rF
      rG=MAX(rGaug,SIZE(G%coef))

!     Pack F and G into the correct precision arrays to pass to the C 
!     code (TYPE CPvec is not interoperable with C, and for single 
!     precision runs the data must be extracted from the double-
!     precision Fortran type). Since the GPU code only accepts G with
!     a rank that is a multiple of rF, "pad" the arrays containing G 
!     with random terms with zero coefficients
      ALLOCATE(coef(rF+rG),base(nD,rF+rG))
      ALLOCATE(nbas(D))
      nbas(:)=F%nbas(:)
      coef(:)=0.0
      coef(1:rF)=REAL(F%coef(:),PREC)
      coef(rF+1:rF+SIZE(G%coef))=REAL(G%coef(1:SIZE(G%coef)),PREC)
      base(:,1:rF)=REAL(F%base(:,:),PREC)
      base(:,rF+1:rF+SIZE(G%coef))=REAL(G%base(:,1:SIZE(G%coef)),PREC)
      IF (rG.gt.SIZE(G%coef)) &
         call random_number(base(:,rF+SIZE(G%coef)+1:rF+rG))

      call FlushCPvec(F) ! Free memory

!     Call the C++ ALS code with the component arrays from TYPE CP-vec
      call ALSCPP(coef,base,nbas,D,rF,rG,nals)

!     Repack F from fcoef and fbase from ALSCPP()
      call NewCPvec(F,nbas,rF)
      F%coef(:)=REAL(coef(1:rF),8)
      F%base(:,:)=REAL(base(:,1:rF),8)

      DEALLOCATE(nbas,coef,base)

      call CPU_TIME(t2)
      gpuals_time=gpuals_time+t2-t1

      end subroutine GPU_ALS_AB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE GPUINTERTWINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
