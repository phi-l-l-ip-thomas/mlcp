!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM CSPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE MODEGRID
      USE SEPDREPN
      USE CPCONFIG
      USE CHEBLIB
      USE SAMPLEPES
      USE RESTARTCSPES
      USE BPSIMPLECF
      USE BPSIMPLECP
      USE REDUCTION

      implicit none
      TYPE (CSpar)   :: csp
      TYPE (GridPar) :: gp
      TYPE (ChebObj), ALLOCATABLE :: CGrid(:)
      TYPE (CPvec)   :: A,X,G,BB,RR,XX
      TYPE (Configs) :: C,B,R
      integer :: rs(33),d(3),t(3),i,ndof
      real*8  :: t1,t2
      character(len=64) :: inpfile

      call idate(d)
      call itime(t)
      call CPU_TIME(t1)
      rs=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
          mod(INT(t1),7),d(1),d(2),d(3),t(1),t(2),t(3),mod(INT(t1),5)/)
      call random_seed(PUT=rs)

      write(*,'(X,A/)') '##############################################'
      write(*,*)       '    Compressed Sensing PES Fitting Program  '
      write(*,*)       '           by Phillip S. Thomas             '
      write(*,*)       '          Version ML2b 06-29-2016           '
      write(*,'(/X,A)') '##############################################'
      write(*,'(/X,4(A,I0),2(A,I2.2))') 'CSPES initialized on ',d(2),&
           '/',d(1),'/',d(3),' at ',t(1),':',t(2),':',t(3)

!     Read the CSPES input file
      inpfile='CS.inp'
      write(*,'(/X,A)') 'Reading Compressed Sensing input file...'
      call ReadCSPESInputFile(csp,inpfile)
      call PrintCSPESInputs(csp)

!     Parallelization setup
      CALL omp_set_num_threads(csp%ncpu)

!     Read the grid parameters input file
      inpfile='grids.inp'
      write(*,'(/X,A/)') 'Reading grids input file...'
      call ReadGridDat(gp,inpfile)
      ndof=SIZE(gp%npts)

!     Create the Chebyshev objects for each DOF
      ALLOCATE(CGrid(ndof))
      write(*,'(A/)') '################################################'
      DO i=1,ndof
         write(*,*) 'Cheb info for DOF #',i
         call NewChebObject(CGrid(i),gp%npts(i),gp%bounds(i,1),&
                            gp%bounds(i,2),gp%egrid(i),gp%bc(i))
         call PrintChebInfo(CGrid(i))
      ENDDO

!     Get the Discrete Cosine Transform matrices for each DOF
      write(*,'(X,A/)') 'Getting transformation matrices...'
      call GetCPChebTfmMat(A,CGrid)

!     Check for restarts and possibly read files
      write(*,'(X,A/)') 'Restart setup...'
      call CSRestartSetup(csp,GP,CGrid,C,B,BB)

!     Get the PES from a data file or by calling a fortran subroutine
      write(*,'(X,A/)') 'Getting the PES...'
      call GetPES(CGrid,B,BB,csp%nsamples,csp%usegrid,&
                  csp%sys(1:5),csp%readfile,A) !!! A added here

!     Save the PES if points were computed
      call SavePESsamples(CGrid,csp,B,BB)

      write(*,'(/X,A)') 'Solving the BP_simple problem:'
      IF (csp%usegrid) THEN
         call SolveBPsimpleCF(C,B,A,X,R,G,csp)
      ELSE
         call SolveBPsimpleCP(C,BB,A,X,RR,G,csp)
      ENDIF

!      write(*,*) 'Reducing the converged potential'
!      write(*,*) 'rank(X) = ',SIZE(X%coef),' before sorting'
!      call reduc_bysorting(X,SIZE(X%coef))
!      write(*,*) 'rank(X) = ',SIZE(X%coef),' after sorting'
!      DO i=1,SIZE(X%coef)
!         call SetReductionParameters(i,30,1.d-12,.FALSE.,'SVD','ALS')
!         call reduc(XX,X)
!         write(*,*) 'rank(X) = ',i,calc_FmG(X,XX)
!         call FlushCPvec(XX)
!      ENDDO

!     Polynomial coefficients as powers of x_i
!     (only valid for finite boundary conditions)
!      write(*,'(/X,A/)') 'Polynomial coefficients (powers of x_i)'
!      call Cheb2PolyCP(X,CGrid,.TRUE.) 
!      call FlushConfigs(C)
!      call GetConfigList(X,csp%pesrank,C)
!      call TrimZeroConfigs(C,csp%sigma)
!      call PrintConfigs(C)

!     Clean up
      call FlushCPvec(BB)
      call FlushCPvec(RR)
      call FlushCPvec(X)
      call FlushCPvec(G)
      call FlushCPvec(A)

      call FlushConfigs(B)
      call FlushConfigs(R)
      call FlushConfigs(C)
      call FlushGridDat(gp)

      call CPU_TIME(t2)
      call idate(d)
      call itime(t)

      write(*,'(/X,A,X,f20.3/)') 'CSPES total CPU run time (s):',t2-t1

      write(*,'(X,4(A,I0),2(A,I2.2))') 'CSPES finished on ',d(2),'/',&
           d(1),'/',d(3),' at ',t(1),':',t(2),':',t(3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM CSPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
