!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE SAMPLEPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Samples PES either by calling an individual PES module or by reading a
! data file. To implement a new PES module, make 2 changes to this file:
! 1) add a statement to GetPES() so that the subroutine that you want to
! call can be identified by a string in your input file
! 2) add a USE statement below to direct the compiler to your PES module

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE CPCONFIG
      USE FFPES
      USE MODVECVEC
      USE CHEBLIB
!!!
      USE LINALG
!!!

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetCPChebTfmMat(A,CGrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get the gridpoint <-> Chebyshev polynomial transformation matrix

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (CPvec), INTENT(OUT)  :: A
      real*8, allocatable  :: tvec(:),tmat(:,:)
      integer, allocatable :: nbas(:)
      integer :: i,gi,gf,ndof

!     Set parameters
      ndof=SIZE(CGrid)

      ALLOCATE(nbas(ndof))
      DO i=1,ndof
         nbas(i)=CGrid(i)%n
!        Add a point if the extrema grid is used
         IF (CGrid(i)%egrid) nbas(i)=nbas(i)+1
      ENDDO

!     Set up the Chebyshev grids for each DOF and build the DCT matrices
      call NewCPmat(A,nbas,1,.FALSE.)
      A%coef(1)=1.d0

      gi=1
      DO i=1,ndof
         gf=gi+A%nbas(i)-1

!        Get the DCT matrices, unwrap into vector form, and store in H
         call DCTmat(CGrid(i)%n,tmat,.FALSE.,CGrid(i)%egrid)
         call Mat2Vec(tvec,tmat,.FALSE.)
         A%base(gi:gf,1)=tvec

         write(*,*) 'Transformation matrix for DOF #',i
         call PrintMatrix(tmat)

         DEALLOCATE(tmat,tvec)
         gi=gi+A%nbas(i)
      ENDDO

!      write(*,*) 'Transformation matrix in CP-format'
!      call PrintCPvec(A)

      DEALLOCATE(nbas)

      end subroutine GetCPChebTfmMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetCPArbTfmMat(A,CGrid,CV)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get the gridpoint <-> Chebyshev polynomial transformation matrix

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (Configs), INTENT(IN) :: CV
      TYPE (CPvec), INTENT(OUT)  :: A
      real*8, allocatable  :: tvec(:),tmat(:,:),F(:,:),VT(:,:),svals(:)
      integer, allocatable :: nbas(:)
      integer :: i,gi,gf,ndof,j

!     Set parameters
      ndof=SIZE(CGrid)

      ALLOCATE(nbas(ndof))
      DO i=1,ndof
         nbas(i)=CGrid(i)%n
!        Add a point if the extrema grid is used
         IF (CGrid(i)%egrid) nbas(i)=nbas(i)+1
      ENDDO

!     Set up the Chebyshev grids for each DOF and build the DCT matrices
      call NewCPmat(A,nbas,1,.FALSE.)
      A%coef(1)=1.d0

      gi=1
      DO i=1,ndof
         gf=gi+A%nbas(i)-1

         call SampleFiberPointRandom(CGrid,F,20000,i,CV)
         call SolveWithSVD(svals,F,tmat,VT)
         write(*,*) 'svals:',(svals(j),j=1,SIZE(svals))
         DEALLOCATE(F,VT,svals)

         call Mat2Vec(tvec,tmat,.FALSE.)
         A%base(gi:gf,1)=tvec

         write(*,*) 'Transformation matrix for DOF #',i
         call PrintMatrix(tmat)

         DEALLOCATE(tmat,tvec)
         gi=gi+A%nbas(i)
      ENDDO

!      write(*,*) 'Transformation matrix in CP-format'
!      call PrintCPvec(A)

      DEALLOCATE(nbas)

      end subroutine GetCPArbTfmMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPES(CGrid,B,BB,nsamp,samplegrid,sys,fnm,A)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (ChebObj), INTENT(IN)    :: CGrid(:)
      TYPE (Configs), INTENT(INOUT) :: B
      TYPE (Configs), ALLOCATABLE   :: V(:)
      TYPE (Configs) :: CV
      TYPE (CPvec), INTENT(INOUT)   :: BB
      character(len=5), intent(in)  :: sys
      character(len=64), intent(in) :: fnm
      integer, intent(in) :: nsamp
      logical, intent(in) :: samplegrid
!!!
      TYPE (CPvec), INTENT(INOUT)  :: A
!!!

!     Read points from a data file...
      IF (sys(1:4).seq.'read') THEN
         write(*,*) 'Reading PES from data file...'
         call ReadPES(CGrid,B,BB,fnm,nsamp,samplegrid)
!     ... or compute a PES with an existing potential routine
      ELSE
         write(*,*) 'Sampling a pre-existing PES...'
         call GetPotential(V,sys,SIZE(CGrid))
         call FF2Configs(V,CV)
         call CalcPES(CGrid,B,BB,nsamp,samplegrid,CV,A)
      ENDIF

      end subroutine GetPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPES(CGrid,B,BB,fnm,nsamp,samplegrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read PES master routine

      implicit none
      TYPE (ChebObj), INTENT(IN)    :: CGrid(:)
      TYPE (Configs), INTENT(OUT)   :: B
      TYPE (CPvec), INTENT(OUT)     :: BB
      character(len=64), intent(in) :: fnm
      integer, intent(in)  :: nsamp
      logical, intent(in)  :: samplegrid
      logical :: success

!      IF (samplegrid) THEN
!         call ReadPESConfigFile(CGrid,B,fnm,nsamp,success)
!      ELSE
         call ReadPESPointsFile(CGrid,BB,fnm,nsamp,success)
!      ENDIF

      IF (.not. success) call &
         AbortWithError('ReadPES(): PES data could not be read')

      IF (samplegrid) THEN
         call CP2Config(BB,B)
         call FlushCPvec(BB)
      ENDIF

      end subroutine ReadPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WritePES(CGrid,B,BB,fnm,samplegrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Write PES master routine

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (Configs), INTENT(IN) :: B
      TYPE (CPvec), INTENT(IN)   :: BB
      TYPE (CPvec) :: Balt
      character(len=64), intent(in) :: fnm
      logical, intent(in) :: samplegrid
      real*8, allocatable :: X(:,:),V(:)

      IF (samplegrid) THEN
!         call WritePESConfigFile(B,fnm)
         call Config2CP(Balt,B)
         call CP2Points(CGrid,Balt,X,V)
         call FlushCPvec(Balt)
         call WritePESPointsFile(X,V,fnm)
      ELSE
         call CP2Points(CGrid,BB,X,V)
      ENDIF

      call WritePESPointsFile(X,V,fnm)
      DEALLOCATE(X,V)

      end subroutine WritePES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPESConfigFile(CGrid,B,fnm,nsamp,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read PES file with potential defined at grid points on CGrid

      implicit none
      TYPE (ChebObj), INTENT(IN)    :: CGrid(:)
      TYPE (Configs), INTENT(OUT)   :: B
      character(len=64), intent(in) :: fnm
      logical, intent(out) :: success
      integer, intent(in)  :: nsamp
      integer, allocatable :: nbas(:)
      integer :: u,InpStat,i,j,ndof,nsampt,ndoft

      ndof=SIZE(CGrid)
      success=.TRUE.

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         success=.FALSE.
         RETURN
      ENDIF

!     Make sure the #DOFs and #samples match what are expected
      nsampt=0
      ndoft=0
      read(u,*,IOSTAT=InpStat) nsampt,ndoft
      IF (InpStat /= 0) THEN
         write(*,*) 'ReadPESConfigFile(): error reading dimensions'
         success=.FALSE.
      ELSEIF (ndoft.ne.ndof) THEN
         write(*,*) 'ReadPESConfigFile(): input dimension mismatch'
         success=.FALSE.
      ENDIF
      IF (.not.success) RETURN

      IF (nsamp.lt.nsampt) THEN
         write(*,'(X,2(A,I0),A)') '# samples available, ',nsampt,&
             '; of these, ',nsamp,' will be used'
      ELSEIF (nsamp.gt.nsampt) THEN
         write(*,'(X,2(A,I0),A)') '# samples requested, ',nsamp,&
         ', exceeds number found in file! Only ',nsampt,' will be used'
      ENDIF
      nsampt=min(nsamp,nsampt)

!     Allocate the configurations array
      ALLOCATE(nbas(ndof))
      DO j=1,ndof
         nbas(j)=CGrid(j)%n
!        Add a point if the extrema grid is used
         IF (CGrid(j)%egrid) nbas(j)=nbas(j)+1
      ENDDO
      call NewConfigs(B,nbas,nsampt)

!     Read the file
      DO i=1,nsampt
         read(u,*,IOSTAT=InpStat) (B%qns(i,j),j=1,ndof),B%coef(i)
         IF (InpStat /= 0) THEN
            write(*,*) 'ReadPESConfigFile(): error reading configs'
            success=.FALSE.
            EXIT
         ENDIF
!        Error if point is outside of the range of nbas
         DO j=1,ndof
            IF (B%qns(i,j).lt.1 .or. B%qns(i,j).gt.nbas(j)) THEN
               write(*,*) 'ReadPESConfigFile(): config exceeds nbas'
               success=.FALSE.
               EXIT
            ENDIF
         ENDDO
         IF (.not.success) EXIT
      ENDDO

      close(u)

      DEALLOCATE(nbas)

      IF (.not.success) call FlushConfigs(B)

      end subroutine ReadPESConfigFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WritePESConfigFile(B,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Write list of PES samples as configurations

      implicit none
      TYPE (Configs), INTENT(IN)    :: B
      character(len=64), intent(in) :: fnm
      character(len=64)  :: frmt
      integer :: u,i,j,ndof,nsamp

!     Set parameters
      ndof=SIZE(B%nbas)
      nsamp=SIZE(B%coef)

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     Write the file
      write(u,'(2(I0,X))') nsamp,ndof
      write(frmt,'(A,I0,A)') '(',ndof,'(I0,X),ES22.15)'
      DO i=1,nsamp
         write(u,frmt) (B%qns(i,j),j=1,ndof),B%coef(i)
      ENDDO

      close(u)

      end subroutine WritePESConfigFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadPESPointsFile(CGrid,B,fnm,nsamp,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read PES file with potential sampled anywhere

      implicit none
      TYPE (ChebObj), INTENT(IN)    :: CGrid(:)
      TYPE (CPvec), INTENT(OUT)     :: B
      character(len=64), intent(in) :: fnm
      logical, intent(out) :: success
      integer, intent(in)  :: nsamp
      real*8, allocatable  :: X(:,:),V(:)
      integer :: u,InpStat,i,j,ndof,nsampt,ndoft

      ndof=SIZE(CGrid)
      success=.TRUE.

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         success=.FALSE.
         RETURN
      ENDIF

!     Make sure the #DOFs and #samples match what are expected
      nsampt=0
      ndoft=0
      read(u,*,IOSTAT=InpStat) nsampt,ndoft
      IF (InpStat /= 0) THEN
         write(*,*) 'ReadPESPointsFile(): error reading dimensions'
         success=.FALSE.
      ELSEIF (ndoft.ne.ndof) THEN
         write(*,*) 'ReadPESPointsFile(): input dimension mismatch'
         success=.FALSE.
      ENDIF
      IF (.not.success) RETURN

      IF (nsamp.lt.nsampt) THEN
         write(*,'(X,2(A,X,I0,X),A)') '# samples available,',nsampt,&
              '; of these,',nsamp,'will be used'
      ELSEIF (nsamp.gt.nsampt) THEN
         write(*,'(X,2(A,X,I0,X),A)') '# samples requested,',nsamp,&
              'exceeds number found in file! Only',nsampt,'will be used'
      ENDIF
      nsampt=min(nsamp,nsampt)

      ALLOCATE(X(nsampt,ndof),V(nsampt))

!     Read the file
      DO i=1,nsampt
         read(u,*,IOSTAT=InpStat) (X(i,j),j=1,ndof),V(i)
         IF (InpStat /= 0) THEN
            write(*,*) 'ReadPESPointsFile(): error reading points'
            success=.FALSE.
            EXIT
         ENDIF
      ENDDO

      close(u)

      IF (success) call Points2CP(CGrid,B,X,V)

      DEALLOCATE(X,V)

      end subroutine ReadPESPointsFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WritePESPointsFile(X,V,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Write list of PES samples

      implicit none
      real*8, intent(in) :: X(:,:),V(:)
      character(len=64), intent(in) :: fnm
      character(len=64)  :: frmt
      integer :: u,i,j,ndof,nsamp

!     Set parameters
      ndof=SIZE(X,2)
      nsamp=SIZE(X,1)

!     Error checking
      IF (SIZE(V).ne.nsamp) call &
         AbortWithError('WritePESPointsFile(): inconsistent X,V sizes')

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     Write the file
      write(u,'(2(I0,X))') nsamp,ndof
      write(frmt,'(A,I0,A)') '(',ndof+1,'(ES22.15,X))'
      DO i=1,nsamp
         write(u,frmt) (X(i,j),j=1,ndof),V(i)
      ENDDO

      close(u)

      end subroutine WritePESPointsFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CalcPES(CGrid,B,BB,nsamp,samplegrid,CV,A)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sample PES master routine

      implicit none
      TYPE (ChebObj), INTENT(IN)    :: CGrid(:)
      TYPE (Configs), INTENT(IN)    :: CV
      TYPE (Configs), INTENT(INOUT) :: B
      TYPE (CPvec), INTENT(INOUT)   :: BB
      TYPE (Configs) :: Bt
      TYPE (CPvec)   :: BBt
      integer, intent(in) :: nsamp
      logical, intent(in) :: samplegrid
      integer :: samples,bdim,tdim
!!!
      TYPE (CPvec), INTENT(INOUT)  :: A
!!!

      IF (samplegrid) THEN
         IF (.not.ALLOCATED(B%coef)) THEN
            call SamplePESPointRandom(CGrid,B,nsamp,CV)

!!!
!            call FlushCPvec(A)
!            call GetCPArbTfmMat(A,CGrid,CV)
!!!

         ELSEIF (SIZE(B%coef).lt.nsamp) THEN
            bdim=SIZE(B%coef)
            samples=nsamp-bdim
            call SamplePESPointRandom(CGrid,Bt,samples,CV)

            write(*,'(X,A,I9)') 'Samples read from file: ',bdim
            write(*,'(X,A,I9)') 'Samples added         : ',samples

!           Increase size of B and add the samples in Bt
            call ResizeConfigList(B,nsamp)
            call GenCopyConfigsWtoV(B,Bt,bdim+1,bdim+samples,1,samples)
            call FlushConfigs(Bt)
         ENDIF

!        Remove duplicate samples, if present
         call removeduplicateconfigs(B,0)
         write(*,'(X,A,I0,A)') 'After removing duplicates, there are '&
                             ,SIZE(B%coef),' samples'
!        The config list was sorted by removeduplicateconfigs() and
!        should be re-randomized here. If the samples are saved to file
!        and only a subset are read later, this guarantees that the
!        subset that gets read is a true random subset
         call randomizeconfiglist(B)
      ELSE
         IF (.not.ALLOCATED(BB%coef)) THEN
            call SamplePESFullRandom(CGrid,BB,nsamp,CV)
         ELSEIF (SIZE(BB%coef).lt.nsamp) THEN
            samples=nsamp-SIZE(BB%coef)
            call SamplePESFullRandom(CGrid,BBt,samples,CV)

            write(*,'(X,A,I9)') 'Samples read from file: ',SIZE(BB%coef)
            write(*,'(X,A,I9)') 'Samples added         : ',samples

            call SUMVECVEC(BB,1.d0,BBt,1.d0)
            call FlushCPvec(BBt)
         ENDIF
      ENDIF

      end subroutine CalcPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SamplePESPointRandom(CGrid,B,nsamp,CV)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Samples PES by evaluating at random grid points on CGrid

      implicit none
      TYPE (ChebObj), INTENT(IN)  :: CGrid(:)
      TYPE (Configs), INTENT(IN)  :: CV
      TYPE (Configs), INTENT(OUT) :: B
      real*8, allocatable  :: x(:)
      integer, allocatable :: nbas(:)
      integer, intent(in)  :: nsamp
      integer :: i,j,ndof

      write(*,'(/X,A)') 'Sampling points randomly from product grid...'

      ndof=SIZE(Cgrid)

!     Compute the number of points-per-dof in the product basis
      ALLOCATE(nbas(ndof))
      DO i=1,ndof
         nbas(i)=CGrid(i)%n
!        Add a point if the extrema grid is used
         IF (CGrid(i)%egrid) nbas(i)=nbas(i)+1
      ENDDO

      ALLOCATE(x(ndof))
      call NewConfigs(B,nbas,nsamp)

!     Sample the PES
      DO i=1,nsamp
!        Pick a random point for each DOF, then map to the desired range
         call random_number(x)
         DO j=1,ndof
            B%qns(i,j)=nint(nbas(j)*x(j)+0.5)
            x(j)=CGrid(j)%mpt(B%qns(i,j))
         ENDDO

!        Evaluate the PES at the sample point
         B%coef(i)=EvalPowerPES(x,CV)
      ENDDO

      DEALLOCATE(x)

      end subroutine SamplePESPointRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SampleFiberPointRandom(CGrid,B,nsamp,k,CV)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Samples PES by evaluating at random grid points on CGrid

      implicit none
      TYPE (ChebObj), INTENT(IN)  :: CGrid(:)
      TYPE (Configs), INTENT(IN)  :: CV
      real*8, allocatable, intent(out) :: B(:,:)
      real*8, allocatable  :: x(:)
      real, parameter :: PI=3.1415926535897932384626433832795028841971
      integer, allocatable :: nbas(:)
      integer, intent(in)  :: nsamp,k
      integer :: i,j,ndof,pt

      write(*,'(/X,A)') 'Sampling fibers randomly from product grid...'

      ndof=SIZE(Cgrid)

      IF (k.lt.1 .or. k.gt.ndof) &
         call AbortWithError('SampleFiberPointRandom(): k out of range')

!     Compute the number of points-per-dof in the product basis
      ALLOCATE(nbas(ndof))
      DO i=1,ndof
         nbas(i)=CGrid(i)%n
!        Add a point if the extrema grid is used
         IF (CGrid(i)%egrid) nbas(i)=nbas(i)+1
      ENDDO

      ALLOCATE(x(ndof),B(nbas(k),nsamp))

!     Sample the PES
      DO i=1,nsamp
!        Pick a random point for each DOF except the k-th,
!        then map to the desired range
         call random_number(x)
         DO j=1,ndof
            IF (j.eq.k) CYCLE

!           "Point random" version
!            pt=nint(nbas(j)*x(j)+0.5)
!            x(j)=CGrid(j)%mpt(pt)

!           "Full random" version
            x(j)=cos(x(j)*PI)  ! maps random point to [-1,1]
            x(j)=mapx(x(j),CGrid(j)%avdf,CGrid(j)%bc,.FALSE.)
         ENDDO

!        Sample a fiber for the k-th DOF
         DO j=1,nbas(k)
            x(k)=CGrid(k)%mpt(j)
!           Evaluate the PES at the sample point
            B(j,i)=EvalPowerPES(x,CV)
         ENDDO
      ENDDO

      DEALLOCATE(x)

      end subroutine SampleFiberPointRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SamplePESFullRandom(CGrid,B,nsamp,CV)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Samples PES by evaluating randomly. Points are mapped onto CGrid via
! the cardinal functions

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (Configs), INTENT(IN) :: CV
      TYPE (CPvec), INTENT(OUT)  :: B
      real*8, allocatable :: X(:,:),V(:)
      integer, intent(in) :: nsamp
      integer :: i,j,ndof
      real, parameter :: PI=3.1415926535897932384626433832795028841971
      character(len=64) :: frmt

      write(*,'(X,A/)') 'Sampling the PES at random points...'

      ndof=SIZE(Cgrid)

      ALLOCATE(X(nsamp,ndof),V(nsamp))

!     Sample the PES
!      write(frmt,'(A,I0,A)') '(8X,A1,',ndof,'(2X,A2,I4.4,7X),2X,A)'
!      write(*,frmt) 'i',('x_',j,j=1,ndof),'VPES'
!      write(frmt,'(A,I0,A)') '(X,I8,',ndof+1,'(X,ES14.7))'
      DO i=1,nsamp
!        Pick a random point for each DOF, then map to the desired range
         call random_number(X(i,:))
         DO j=1,ndof
!           Map to the range specified by the boundary conditions
!           Start with a random uniform distribution 0 <= x <= PI
            X(i,j)=cos(X(i,j)*PI)  ! maps random point to [-1,1]
            X(i,j)=mapx(X(i,j),CGrid(j)%avdf,CGrid(j)%bc,.FALSE.)
         ENDDO

!        Evaluate the PES at the sample point
         V(i)=EvalPowerPES(X(i,:),CV)

!         write(*,frmt) i,(X(i,j),j=1,ndof),V(i)
      ENDDO
!      write(*,*)

!     Convert points to CP-vector
      call Points2CP(CGrid,B,X,V)

      DEALLOCATE(X,V)

      end subroutine SamplePESFullRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Points2CP(CGrid,B,X,V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts data points sampled at arbitrary locations to a CP-vector
! whose basis functions are DVR points stored in CGrid

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (CPvec), INTENT(OUT)  :: B
      real*8, intent(in)   :: X(:,:),V(:)
      integer, allocatable :: nbas(:)
      integer :: i,j,k,n,gi,ndof,nsamp
      real*8  :: z

!     Set parameters
      nsamp=SIZE(X,1)
      ndof=SIZE(X,2)

!     Error checking
      IF (SIZE(V).ne.nsamp) call &
         AbortWithError('Points2CP(): inconsistent X,V sizes')

!     Fill nbas based on the structure of CGrid
      ALLOCATE(nbas(ndof))
      DO j=1,ndof
         nbas(j)=CGrid(j)%n
!        Add a point if the extrema grid is used
         IF (CGrid(j)%egrid) nbas(j)=nbas(j)+1
      ENDDO

      call NewCPvec(B,nbas,nsamp)
      DEALLOCATE(nbas)

      DO i=1,nsamp
         B%coef(i)=V(i)
!        Map the point to [-1,1]
         gi=0
         DO j=1,ndof
            n=B%nbas(j)
            z=mapx(X(i,j),CGrid(j)%avdf,CGrid(j)%bc,.TRUE.)

!           Make sure xt is inside Chebyshev range
            IF (abs(z).gt.1.d0) call &
              AbortWithError('Points2CP(): z is outside [-1,1]')
!           Convert z to a linear combination of cardinal
!           interpolating functions and store in B%base
            DO k=1,n
               B%base(gi+k,i)=ChebCard(CGrid(j)%n,k-1,z,CGrid(j)%egrid)
            ENDDO
            gi=gi+n
         ENDDO
      ENDDO

      end subroutine Points2CP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CP2Points(CGrid,B,X,V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts CP-vector of data defined on linear combinations of DVR
! points back into original list of points

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (CPvec), INTENT(IN)   :: B
      real*8, allocatable, intent(out) :: X(:,:),V(:)
      integer, allocatable :: nbas(:)
      integer :: i,j,k,n,gi,ndof,nsamp
      real*8  :: z

!     Set parameters
      nsamp=SIZE(B%coef)
      ndof=SIZE(B%nbas)

      ALLOCATE(X(nsamp,ndof),V(nsamp))

      DO i=1,nsamp
         V(i)=B%coef(i)
         gi=0
         DO j=1,ndof
            n=B%nbas(j)
!           Calc z inside of [-1,1] as a linear combination of grid points
            z=0.d0
            DO k=1,n
               z=z+B%base(gi+k,i)*CGrid(j)%pt(k)
            ENDDO

!           Map z back to original interval
            X(i,j)=mapx(z,CGrid(j)%avdf,CGrid(j)%bc,.FALSE.)

            gi=gi+n
         ENDDO
      ENDDO

      end subroutine CP2Points

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Cheb2PolyCP(F,CGrid,forw)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts CP-vector of Chebyshev coefs to coefs of powers of x.
! Normalization is done at the end.

      implicit none
      TYPE (ChebObj), INTENT(IN)  :: CGrid(:)
      TYPE (CPvec), INTENT(INOUT) :: F
      logical, intent(in) :: forw
      integer :: i,j,gi,gf,rk,ndof
      real*8  :: fac

!     Set parameters
      rk=SIZE(F%coef)
      ndof=SIZE(F%nbas)

!     Convert coefs for each 'f' in the base
      gi=1
      DO j=1,ndof
         gf=gi+F%nbas(j)-1
         DO i=1,rk
            call cheb2poly(F%base(gi:gf,i),CGrid(j)%avdf,&
                           forw,CGrid(j)%egrid)
         ENDDO
         gi=gf+1
      ENDDO

      call NORMBASE(F)

      end subroutine Cheb2PolyCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function EvalChebCPPES(x,C,CGrid)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Evaluates PES stored as a CP-vector of Chebyshev coefficients, using
! Clenshaw recursion, at point x = {x_1,x_2,...,x_d}.

      implicit none
      TYPE (ChebObj), INTENT(IN) :: CGrid(:)
      TYPE (CPvec), INTENT(IN)   :: C
      real*8, intent(in)   :: x(:)
      real*8, allocatable  :: y(:),clprods(:,:),prods(:)
      integer, allocatable :: gi(:),gf(:)
      integer :: i,j,k,n,ndof,rk,rtot
      real*8  :: EvalChebCPPES

!     Set parameters
      ndof=SIZE(C%nbas)
      rk=SIZE(C%coef)
      rtot=ndof*rk

!     Error checking
      IF (SIZE(x).ne.ndof) &
         call AbortWithError("EvalChebCPPES(): x and C size mismatch")

      ALLOCATE(clprods(ndof,rk),prods(rk),gi(ndof),gf(ndof),y(ndof))

!     Set the initial/final index ranges for each DOF, map x to [-1,1]
      gi(1)=1
      gf(1)=C%nbas(1)
      y(1)=mapx(x(1),CGrid(1)%avdf,CGrid(1)%bc,.TRUE.)
      DO j=2,ndof
         gi(j)=gi(j-1)+C%nbas(j-1)
         gf(j)=gf(j-1)+C%nbas(j)
         y(j)=mapx(x(j),CGrid(j)%avdf,CGrid(j)%bc,.TRUE.)
      ENDDO

!     Evaluate the PES via Clenshaw recurrence for each term/DOF
!$omp parallel
!$omp do private(i,j,k,n)
      DO k=1,rtot
         i=(k-1)/ndof+1
         j=mod(k-1,ndof)+1
         n=C%nbas(j)-1
         clprods(j,i)=ClenshawRecur(C%base(gi(j):gf(j),i),&
                                    y(j),n,CGrid(j)%egrid)
      ENDDO
!$omp end do
!$omp end parallel

!     Form the products for each term
!$omp parallel
!$omp do private(i,j)
      DO i=1,rk
         prods(i)=C%coef(i)
         DO j=1,ndof
            prods(i)=prods(i)*clprods(j,i)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel

!     Accumulate the sum over terms
      EvalChebCPPES=0.d0
      DO i=1,rk
         EvalChebCPPES=EvalChebCPPES+prods(i)
      ENDDO

      DEALLOCATE(clprods,prods,gi,gf,y)

      end function EvalChebCPPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function EvalPowerPES(y,C)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Evaluates PES stored as a configuration vector of powers of y, at 
! point y = {y_1,y_2,...,y_d}.

      implicit none
      TYPE (Configs), INTENT(IN) :: C
      real*8, intent(in)  :: y(:)
      real*8, allocatable :: x(:,:),prods(:)
      integer :: i,j,ndof,rk,nord
      real*8  :: EvalPowerPES

!     Set parameters
      ndof=SIZE(C%nbas)
      rk=SIZE(C%coef)
      nord=MAXVAL(C%nbas)

!     Error checking
      IF (SIZE(y).ne.ndof) &
         call AbortWithError("EvalPowerPES(): y and C size mismatch")

      ALLOCATE(x(ndof,nord),prods(rk))

!     Pre-calculate powers of y
      DO i=1,ndof
         x(i,1)=1.d0
         DO j=2,nord
            x(i,j)=x(i,j-1)*y(i)
         ENDDO
      ENDDO
    
!     Form the products for each term
!$omp parallel
!$omp do private(i,j)
      DO i=1,rk
         prods(i)=C%coef(i)
         DO j=1,ndof
            prods(i)=prods(i)*x(j,C%qns(i,j))
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel

!     Accumulate the sum over terms
      EvalPowerPES=0.d0
      DO i=1,rk
         EvalPowerPES=EvalPowerPES+prods(i)
      ENDDO

      DEALLOCATE(x,prods)

      end function EvalPowerPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE SAMPLEPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
