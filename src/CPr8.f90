!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains structures needed for CP representation, real(kind=8)

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE MYACC
      USE SEPDREPN
      USE, INTRINSIC :: ISO_C_BINDING
#if ACC_ENABLED
      USE OPENACC
      USE CUDAFOR
#endif

      IMPLICIT NONE
      real(kind=8), allocatable, private :: copy_time(:),mpi_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE CP8
         logical :: live_on_host = .FALSE.
         logical :: live_on_device = .FALSE.
         INTEGER, ALLOCATABLE :: look(:,:)
         INTEGER, ALLOCATABLE :: nbas(:)
         INTEGER, ALLOCATABLE :: rows(:)
         INTEGER, ALLOCATABLE :: cols(:)
         REAL(KIND=8), ALLOCATABLE :: base(:),coef(:)
         CONTAINS
            PROCEDURE :: new => NewGen_CP8
            PROCEDURE :: new0 => ZeroGen_CP8
            PROCEDURE :: newrand => RandomGen_CP8
            PROCEDURE :: newvec => NewVec_CP8
            PROCEDURE :: newsquarematrix => NewSqmat_CP8
            PROCEDURE :: identity => IdentityMatrix_CP8
            PROCEDURE :: clone => NewRef_CP8
            PROCEDURE :: clone0 => ZeroRef_CP8
            PROCEDURE :: cloneid => IdentityRef_CP8
            PROCEDURE :: clonerand => RandomRef_CP8
            PROCEDURE :: flush => Flush_CP8
            PROCEDURE :: zero => SetZero_CP8
            PROCEDURE :: show => ShowStats_CP8
            PROCEDURE :: printvec => PrintVec_CP8
            PROCEDURE :: printmat => PrintMat_CP8
            PROCEDURE :: C => GetCoef_CP8
            PROCEDURE :: R => GetRank_CP8
            PROCEDURE :: D => GetNdof_CP8
            PROCEDURE :: MN => GetNbas_CP8
            PROCEDURE :: M => GetRows_CP8
            PROCEDURE :: N => GetCols_CP8
            PROCEDURE :: BS => GetBasisStart_CP8
            PROCEDURE :: BF => GetBasisFinish_CP8
            PROCEDURE :: MS => GetModeStart_CP8
            PROCEDURE :: MF => GetModeFinish_CP8
            PROCEDURE :: Get => GetEntry_CP8
            PROCEDURE :: Put => PutEntry_CP8
            PROCEDURE :: ok => CheckCoefs_CP8
            PROCEDURE :: same => CheckNbas_CP8
            PROCEDURE :: copy_general => GenCopyWtoV_CP8
            PROCEDURE :: copy_terms => Copy_terms_CP8
            PROCEDURE :: copy_modes => Copy_modes_CP8
            PROCEDURE :: copy_1mode => Copy_1mode_CP8
            PROCEDURE :: copyfrom => Copy_all_CP8
            PROCEDURE :: replace => ReplaceVwithW_CP8
            PROCEDURE :: resize => Resize_CP8
            PROCEDURE :: extendrand => RandomExtend_CP8 
            PROCEDURE :: submatrix => ExtractSubmatrix_CP8
            PROCEDURE :: changesign => VecSignChange_CP8
            PROCEDURE :: mult => VecScalarMult_all_CP8
            PROCEDURE :: mult_terms => VecScalarMult_gen_CP8
            PROCEDURE :: sumcp => Sum_CP8
            PROCEDURE :: sumlccp => SumLinearCombination_CP8
            PROCEDURE :: diag => VectoDiagMatrix_CP8
            PROCEDURE :: transpose => MatrixTranspose_CP8
            PROCEDURE :: trim => TrimZeros_CP8
            PROCEDURE :: zerooffdiagonalall => MatrixZeroOffDiag_all_CP8
            PROCEDURE :: zerooffdiagonalmode => MatrixZeroOffDiag_one_CP8
            PROCEDURE :: zerooffdiagonalgen => MatrixZeroOffDiag_gen_CP8
            PROCEDURE :: multoutcoef => MultOutCoefSmallest_CP8
            PROCEDURE :: multoutcoefmode => MultOutCoefbyMode_CP8
            PROCEDURE :: distributecoef => DistributeCoef_CP8
            PROCEDURE :: largestentryinterm => GetRank1DominantEntry_gen_CP8
            PROCEDURE :: largestentryrank1 => GetRank1DominantEntry_1_CP8
            PROCEDURE :: extractvec => ExtractVec_CP8
            PROCEDURE :: extractdiagonal => ExtractDiagfromMatrix_CP8
            PROCEDURE :: putsubmatrix => PutSubmatrix_CP8
            PROCEDURE :: tensorelement => ExtractTensorElement_CP8
            PROCEDURE :: modejoin => ModeJoin_CP8
            PROCEDURE :: compareentries => EntrywiseCompare_CP8
            PROCEDURE :: createondevice => createondevice_CP8
            PROCEDURE :: copyintodevice => copyintodevice_CP8
            PROCEDURE :: updateondevice => updateondevice_CP8
            PROCEDURE :: updatefromdevice => updatefromdevice_CP8
            PROCEDURE :: deletefromdevice => deletefromdevice_CP8
            PROCEDURE :: copyoutfromdevice => copyoutfromdevice_CP8
            PROCEDURE :: tocp => thistoCP_CP8 
            PROCEDURE :: fromcp => thisfromCP_CP8
      END TYPE CP8

!!!   MOVE to Reduction or SVD modules
      INTERFACE CP2DtoMat
        MODULE PROCEDURE CP2DtoMat_CP8
!        MODULE PROCEDURE CP2DtoMat_CP4
      END INTERFACE CP2DtoMat

!!!   MOVE to Reduction or SVD modules
      INTERFACE CP2DtoUW
        MODULE PROCEDURE CP2DtoUW_CP8
!        MODULE PROCEDURE CP2DtoUW_CP4
      END INTERFACE CP2DtoUW

      INTERFACE Bcast_CP
        MODULE PROCEDURE Bcast_CP8
!        MODULE PROCEDURE Bcast_CP4
      END INTERFACE Bcast_CP

      INTERFACE MPI_Sync_block_CP
        MODULE PROCEDURE MPI_Sync_block_CP8
!        MODULE PROCEDURE MPI_Sync_block_CP4
      END INTERFACE MPI_Sync_block_CP

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_CPr8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(copy_time(mpinodes),mpi_time(mpinodes))
      copy_time(:) = 0.d0
      mpi_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_CPr8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_CPr8_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_CPr8_Module()
      call Get_MPI_Timings('CPr8: sum+copy vectors',copy_time)
      call Get_MPI_Timings('CPr8: CP8 MPI synchronization',mpi_time)
      MODULE_SETUP = .FALSE.
      deallocate(copy_time,mpi_time)

      end subroutine Dispose_CPr8_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewGen_CP8(v,rk,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! General subroutine for initializing outer product of CP-vectors or
! matrices. Symmetric storage is possible for square matrix factors

      implicit none
      CLASS (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:), cols(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER :: d,r,ndof,i,ndim

      ndof=SIZE(rows)

!     Check for correct sizes and consistency in dimensions
      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewGen_CP8()')
      ENDIF

      IF (rk.lt.1) THEN
         write(*,*) 'Error: rk must be at least 1!'
         call AbortWithError('Error in NewGen_CP8()')
      ENDIF

      IF (SIZE(cols).ne.ndof) THEN
         write(*,*) 'Error: size of "cols" array differs from ndof!'
         call AbortWithError('Error in NewGen_CP8()')
      ENDIF

      DO i=1,ndof
         IF (rows(i).lt.1) THEN
            write(*,*) 'Error: number of rows must be at least 1!'
            call AbortWithError('Error in NewGen_CP8()')
         ENDIF

         IF (cols(i).lt.1) THEN
            write(*,*) 'Error: number of cols must be at least 1!'
            call AbortWithError('Error in NewGen_CP8()')
         ENDIF
      ENDDO

!     Error if already allocated
      IF (ALLOCATED(v%nbas) .or. ALLOCATED(v%rows) &
      .or. ALLOCATED(v%cols) .or. ALLOCATED(v%look) &
      .or. ALLOCATED(v%base) .or. ALLOCATED(v%coef)) THEN
         write(*,*) 'Attempt to reallocate existing CP vector'
         call AbortWithError('Error in New_CP()')
      ENDIF

!     Construct lookup table for terms in base; count base entries
      ALLOCATE(v%look(rk,ndof))
      i=0
      DO d=1,ndof
         ndim=rows(d)*cols(d)
         DO r=1,rk
            v%look(r,d)=i
            i=i+ndim
         ENDDO
      ENDDO

      ALLOCATE(v%base(i),v%coef(rk))

!     Set arrays containing dimensions
      ALLOCATE(v%nbas(ndof),v%rows(ndof),v%cols(ndof))
      v%nbas(:)=rows(:)*cols(:)
      v%rows(:)=rows(:)
      v%cols(:)=cols(:)

      v%live_on_host=.TRUE.

      end subroutine NewGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewRef_CP8(v,w,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CP-matrix v, using sizes in w (including the rank, if not
! passed as an optional argument)

      implicit none
      CLASS (CP8)  :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk

      IF (v%live_on_host) &
         call AbortWithError('NewRef_CP8(): attempt to reallocate v')

      IF (present(rk)) THEN
         call v%new(rk,w%rows,w%cols)
      ELSE
         call v%new(w%R(),w%rows,w%cols)
      ENDIF

      IF (w%live_on_device) call v%createondevice()

      end subroutine NewRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes CP-format type

      implicit none
      CLASS (CP8) :: v

!     Remove from device if applicable
      if (v%live_on_device) call v%deletefromdevice

!     Deallocate arrays
      IF (ALLOCATED(v%look)) DEALLOCATE(v%look)
      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)
      IF (ALLOCATED(v%rows)) DEALLOCATE(v%rows)
      IF (ALLOCATED(v%cols)) DEALLOCATE(v%cols)
      IF (ALLOCATED(v%base)) DEALLOCATE(v%base)
      IF (ALLOCATED(v%coef)) DEALLOCATE(v%coef)

      v%live_on_host=.FALSE.

      end subroutine Flush_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetZero_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros CP-vector

      implicit none
      CLASS (CP8) :: v
      integer :: j,ndof

      IF (.not.ALLOCATED(v%base)) THEN
         call AbortWithError("SetZero_CP8(): v not allocated")
      ENDIF
      v%base=0.d0
      v%coef=0.d0

      end subroutine SetZero_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ShowStats_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Shows mode sizesm sym, rank of v

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: d
      character*64 :: frmt

      frmt='(X,I4,X,I5,X,X,I5,2X,L3)'
      write(*,*) '------- CP stats -------'
      write(*,*) 'live on host  :',v%live_on_host
      write(*,*) 'live on device:',v%live_on_device
      write(*,'(X,A,I4,A,I0)') 'ndof = ',v%D(),'; rank = ',v%R()
      write(*,*) 'mode [rows x cols]'
      do d=1,v%D()
         write(*,frmt) d,v%M(d),v%N(d)
      enddo
      write(*,*) '------------------------'

      end subroutine ShowStats_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintVec_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP8-vector in neat format

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: r,rk,d,ndof,i,k,n,m
      character*64 :: frmt

      rk=v%R()
      ndof=v%D()

      write(frmt,'(A,I0,A)') '(A,2X,',rk,'(ES17.10,X))'
      write(*,frmt) ' Vcoef =',(v%coef(r),r=1,rk)
      write(frmt,'(A,I0,A)') '(2(I4),X,',rk,'f18.10)'
      do d=1,ndof
         m=v%M(d)
         n=v%N(d)
         do k=1,n
            do i=1,m
               write(*,frmt) d,(k-1)*m+i,(v%base(v%look(r,d)+(k-1)*m+i),r=1,rk)
            enddo
         enddo
      enddo
      write(*,*)

      end subroutine PrintVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMat_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format, general routine

      implicit none
      CLASS (CP8), INTENT(IN), TARGET :: v
      integer :: r,rk,d,ndof,i,n,k,m
      character*64 :: frmt

      rk=v%R()
      ndof=v%D()

      write(*,*)
      DO r=1,rk
         write(*,'(A,I0,A,ES23.16)') 'RANK: ',r,'; coef = ',v%coef(r)
            DO d=1,ndof
               write(*,'(/A,I0)') 'dof : ',d
               m=v%M(d)
               n=v%N(d)
               write(frmt,'(A,I0,A)') '(',n,'(X,f14.6))'
               DO i=1,m
                  write(*,frmt) (v%base(v%look(r,d)+(k-1)*m+i),k=1,n)
               ENDDO
            ENDDO
            write(*,*)
      ENDDO
      write(*,*)

      end subroutine PrintMat_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCoef_CP8(v,r) result(coef)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns coefficient of term 'r' of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: r
      real(kind=8) :: coef

      IF (r.lt.1 .or. r.gt.v%R()) THEN
         write(*,*) 'GetCoef(): r (',r,') out of range: [1,',v%R(),']'
      ENDIF

      IF (ALLOCATED(v%base)) THEN
         coef=v%coef(r)
      ELSE
         write(*,*) 'Error: rk must be at least 1!'
         call AbortWithError('Error in GetCoef_CP8()')
      ENDIF

      end function GetCoef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRank_CP8(v) result(rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns rank of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: rk

      IF (ALLOCATED(v%look)) THEN
         rk=SIZE(v%look,1)
      ELSE
         rk=0
      ENDIF

      end function GetRank_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetNdof_CP8(v) result(ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: ndof

      IF (ALLOCATED(v%look)) THEN
         ndof=SIZE(v%look,2)
      ELSE
         ndof=0
      ENDIF

      end function GetNdof_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetNbas_CP8(v,d) result(nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of basis components in mode 'd' of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: nbas

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'Getnbas(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%nbas)) THEN
         nbas=v%nbas(d)
      ELSE
         nbas=0
      ENDIF

      end function GetNbas_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRows_CP8(v,d) result(rows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of basis rows in mode 'd' of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: rows

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'Getrows(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%rows)) THEN
         rows=v%rows(d)
      ELSE
         rows=0
      ENDIF

      end function GetRows_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetCols_CP8(v,d) result(cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of basis columns in mode 'd' of vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: cols

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'Getcols(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%cols)) THEN
         cols=v%cols(d)
      ELSE
         cols=0
      ENDIF

      end function GetCols_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetBasisStart_CP8(v,r,d) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns starting index of basis of term 'r', mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: r,d
      integer :: res

      IF (r.lt.1 .or. r.gt.v%R()) THEN
         write(*,*) 'GetBaseStart(): r (',r,') out of range: [1,',v%R(),']'
      ENDIF
      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'GetBaseStart(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%look)) THEN
         res=v%look(r,d)+1
      ELSE
         res=0
      ENDIF

      end function GetBasisStart_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetBasisFinish_CP8(v,r,d) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns finishing index of basis of term 'r', mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: r,d
      integer :: res

      IF (r.lt.1 .or. r.gt.v%R()) THEN
         write(*,*) 'GetBaseFinish(): r (',r,') out of range: [1,',v%R(),']'
      ENDIF
      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'GetBaseFinish(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%look)) THEN
         res=v%look(r,d)+v%MN(d)
      ELSE
         res=0
      ENDIF

      end function GetBasisFinish_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetModeStart_CP8(v,d) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns starting index of basis of first term, mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: res

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'GetBaseStart(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%look)) THEN
         res=v%look(1,d)+1
      ELSE
         res=0
      ENDIF

      end function GetModeStart_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetModeFinish_CP8(v,d) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns finishing index of basis of last term, mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: d
      integer :: rk,res

      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'GetBaseFinish(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      IF (ALLOCATED(v%look)) THEN
         res=v%look(v%R(),d)+v%MN(d)
      ELSE
         res=0
      ENDIF

      end function GetModeFinish_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetEntry_CP8(v,i,k,r,d) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets (i,k)-th matrix entry for rank 'r', mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: i,k,r,d
      integer :: m
      real(kind=8) :: res

      IF (r.lt.1 .or. r.gt.v%R()) THEN
         write(*,*) 'GetEntry(): r (',r,') out of range: [1,',v%R(),']'
      ENDIF
      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'GetEntry(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      m=v%M(d)

      IF (i.lt.1 .or. i.gt.m) THEN
         write(*,*) 'GetEntry(): i (',i,') out of range: [1,',m,']'
      ENDIF
      IF (k.lt.1 .or. k.gt.v%N(d)) THEN
         write(*,*) 'GetEntry(): k (',k,') out of range: [1,',v%N(d),']'
      ENDIF

      res=v%base(v%look(r,d)+(k-1)*m+i)

      end function GetEntry_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutEntry_CP8(v,i,k,r,d,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets (i,k)-th matrix entry for rank 'r', mode 'd' in vector 'v'

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      integer, intent(in) :: i,k,r,d
      integer :: m
      real(kind=8) :: val

      IF (r.lt.1 .or. r.gt.v%R()) THEN
         write(*,*) 'PutEntry(): r (',r,') out of range: [1,',v%R(),']'
      ENDIF
      IF (d.lt.1 .or. d.gt.v%D()) THEN
         write(*,*) 'PutEntry(): d (',d,') out of range: [1,',v%D(),']'
      ENDIF

      m=v%M(d)

      IF (i.lt.1 .or. i.gt.m) THEN
         write(*,*) 'PutEntry(): i (',i,') out of range: [1,',m,']'
      ENDIF
      IF (k.lt.1 .or. k.gt.v%N(d)) THEN
         write(*,*) 'PutEntry(): k (',k,') out of range: [1,',v%N(d),']'
      ENDIF

      v%base(v%look(r,d)+(k-1)*m+i)=val

      end subroutine PutEntry_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckCoefs_CP8(v) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks coefficients for good (non-NaN) values

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      LOGICAL :: res
      INTEGER :: i,rk

      res=.TRUE.
      rk=v%R()

      if (v%live_on_device) then
#if ACC_ENABLED
         !$acc data present(v%coef)
         !$acc parallel loop gang vector async
         do i=1,rk
            IF (v%coef(i).ne.v%coef(i)) res=.FALSE.
         enddo
         !$acc end data
#endif
      else
         DO i=1,rk
            IF (v%coef(i).ne.v%coef(i)) THEN
               res=.FALSE.
               EXIT
            ENDIF
         ENDDO
      endif

      end function CheckCoefs_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CheckNbas_CP8(v1,v2) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      implicit none
      CLASS (CP8) :: v1
      TYPE (CP8), INTENT(IN) :: v2
      LOGICAL :: res
      INTEGER :: i

      res=.TRUE.

      IF (v1%D().ne.v2%D()) THEN
         res=.FALSE.
      ELSE
         DO i=1,v1%D()
            IF ((v1%rows(i).ne.v2%rows(i)) .or. &
                (v1%cols(i).ne.v2%cols(i))) THEN
               res=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      end function CheckNbas_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GenCopyWtoV_CP8(v,w,vi,ve,wi,we,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V
! leaving W intact. v must be allocated beforehand.

      implicit none
      CLASS(CP8) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN) :: vi,ve,wi,we
      LOGICAL, INTENT(IN) :: modes(:)
      INTEGER :: i,n,rkv,rkw,d,ndof,bsv,bsw,bfv,bfw
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPr8_Module()

      ndof=w%D()
      rkw=w%R()
      rkv=v%R()

      IF (SIZE(modes).ne.w%D()) THEN
         write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',w%D()
         call AbortWithError('GenCopyWtoV(): bad size for modes array')
      ENDIF

      IF (.not.v%same(w)) THEN
         write(*,*) 'v, w must have same dimensions'
         write(*,*) 'v:'
         call v%show()
         write(*,*) 'w:'
         call w%show()
         CALL AbortWithError('Error in GenCopyWtoV()')
      ENDIF

      IF (vi.lt.1 .or. ve.gt.rkv .or. vi.gt.ve .or. &
          wi.lt.1 .or. we.gt.rkw .or. wi.gt.we .or. &
          we-wi.ne.ve-vi) THEN
          write(*,'(2A,6(X,I0))') 'Bad rank indices: ',&
          'vi,ve,rkv,wi,we,rkw =',vi,ve,rkv,wi,we,rkw
          CALL AbortWithError('Error in GenCopyWtoV()')
      ENDIF

      call CPU_TIME(ti1)
      call nvtx_start('CP8 copy')

!     Copy on device if both vectors present; otherwise copy on host
      if (v%live_on_device .and. w%live_on_device) then
#if ACC_ENABLED
         !$acc data present(v%base,v%coef,w%base,w%coef)
         do d=1,ndof
            if (modes(d)) then
               bsv=v%BS(vi,d)
               bsw=w%BS(wi,d)
               bfv=v%BF(ve,d)
               bfw=w%BF(we,d)
               n=bfv-bsv
               !$acc parallel loop gang vector async
               do i=0,n
                  v%base(bsv+i)=w%base(bsw+i)
               enddo
            endif
         enddo
         n=ve-vi
         !$acc parallel loop gang vector async
         do i=0,n
            v%coef(vi+i)=w%coef(wi+i)
         enddo
         !$acc wait
         !$acc end data
#endif
      else
         do d=1,ndof
            if (modes(d)) then
               bsv=v%BS(vi,d)
               bsw=w%BS(wi,d)
               bfv=v%BF(ve,d)
               bfw=w%BF(we,d)
               v%base(bsv:bfv)=w%base(bsw:bfw)
            endif
         enddo
         v%coef(vi:ve)=w%coef(wi:we)
      endif

      call nvtx_stop()
      call CPU_TIME(ti2)
      copy_time=copy_time+ti2-ti1

      end subroutine GenCopyWtoV_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Copy_terms_CP8(v,w,vi,ve,wi,we)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V
! leaving W intact, for all modes. v must be allocated beforehand.

      implicit none
      CLASS(CP8) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN) :: vi,ve,wi,we
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=w%D()
      allocate(modes(ndof))
      modes(:)=.TRUE.
      call v%copy_general(w,vi,ve,wi,we,modes)
      deallocate(modes)

      end subroutine Copy_terms_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Copy_modes_CP8(v,w,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies mode terms of W into V, for all ranks. W and V must have same
! ranks; if V is not yet allocated then it is allocated by this routine

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8), INTENT(IN) :: w
      LOGICAL, INTENT(IN) :: modes(:)
      integer :: wi,we

      if (v%R().eq.0) call v%clone(w)
      
      if (v%R().ne.w%R()) then
         write(*,*) 'v, w have ranks ',v%R(),', ',w%R(),&
         ' respectively; must be equal'
         call AbortWithError('Copy_modes_CP8(): ranks mismatch')
      endif

      wi=1
      we=w%R()
   
      call v%copy_general(w,wi,we,wi,we,modes)

      end subroutine Copy_modes_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Copy_1mode_CP8(v,w,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN) :: d
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=W%D()

      IF (d.lt.1 .or. d.gt.ndof) THEN
         write(*,*) 'Mode ',d,' out of range: [',1,',',ndof,']'
         call AbortWithError('Copy_1mode_CP8(): wrong mode d')
      ENDIF

      allocate(modes(ndof))
      modes(:)=.FALSE.
      modes(d)=.TRUE.
      call v%copy_modes(w,modes)
      deallocate(modes)

      end subroutine Copy_1mode_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Copy_all_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER   :: rk

      rk=w%R()
      call v%clone(w,rk)
      call v%copy_terms(w,1,rk,1,rk)

      end subroutine Copy_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReplaceVwithW_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, disposing W afterwards

      implicit none
      CLASS (CP8) :: v
      TYPE (CP8)  :: w

      call v%flush
      call v%copyfrom(w)
      call w%flush

      end subroutine ReplaceVwithW_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Resize_CP8(v,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resizes vector v, either by truncating at a smaller rank or by
! adding space for extra terms.

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8) :: w
      INTEGER, INTENT(IN) :: rk
      INTEGER :: rkv

      rkv=v%R()

      IF (rk.lt.1) &
         call AbortWithError('Error in ResizeV(): rk < 1')

      IF (rk.ne.rkv) THEN
         call w%clone(v,rk)
         call w%copy_terms(v,1,MIN(rkv,rk),1,MIN(rkv,rk))
         call v%replace(w)
      ENDIF

      end subroutine Resize_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractSubmatrix_CP8(W,irs,irf,ics,icf) result(V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-submatrix W from list of row, col starting and ending indices

      implicit none
      CLASS (CP8), INTENT(IN) :: W
      TYPE (CP8) :: V
      integer, intent(in)  :: irs(:),irf(:),ics(:),icf(:)
      integer, allocatable :: rows(:),cols(:)
      integer :: d,ndof,r,rk,i,k,m,n

      ndof=W%D()
      rk=W%R()

!     Loads of error checking
      IF (SIZE(irs).ne.ndof .or. SIZE(irf).ne.ndof .or. &
          SIZE(ics).ne.ndof .or. SIZE(icf).ne.ndof) THEN
          write(*,*) 'ndof of M = ',ndof,'; must equal ndof of ALL of',&
          ' irs,irf,ics,icf, which = ',&
          SIZE(irs),SIZE(irf),SIZE(ics),SIZE(icf)
          call AbortWithError('ExtractCPsubmatrix(): bad ranges ndof')
      ENDIF
      DO d=1,ndof
         IF (irs(d).lt.1 .or. irf(d).gt.W%M(d) .or. &
             irs(d).gt.irf(d)) THEN
            write(*,*) 'Row ranges: [irs(',d,'),irf(',d,')] = [',&
            irs(d),',',irf(d),'] must be in range [1,',W%M(d),']'
            call AbortWithError('ExtractCPsubmatrix(): bad rows ranges')
         ENDIF
         IF (ics(d).lt.1 .or. icf(d).gt.W%N(d) .or. &
             ics(d).gt.icf(d)) THEN
            write(*,*) 'Col ranges: [ics(',d,'),icf(',d,')] = [',&
            ics(d),',',icf(d),'] must be in range [1,',W%N(d),']'
            call AbortWithError('ExtractCPsubmatrix(): bad cols ranges')
         ENDIF
      ENDDO

      allocate(rows(ndof),cols(ndof))
      rows(:)=irf(:)-irs(:)+1
      cols(:)=icf(:)-ics(:)+1
      call V%new(rk,rows,cols)
      deallocate(rows,cols)

!     Copy the selected portion v <- w
      V%coef(:)=W%coef(:)
      DO d=1,ndof
         m=V%M(d)
         n=V%N(d)
         DO r=1,rk
            DO k=1,n
               DO i=1,m
                  call V%put(i,k,r,d,W%get(irs(d)+i-1,ics(d)+k-1,r,d))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      end function ExtractSubmatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecSignChange_CP8(v,ri,re)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Changes sign of CP-vec by negating the base for the first DOF. ri and 
! re are the rank indices over which over which to change the sign

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      integer, intent(in) :: ri,re
      integer :: bs,bf,i

      bs=v%BS(ri,1)
      bf=v%BF(re,1)
      if (v%live_on_device) then
#if ACC_ENABLED
         !$acc data present(v%base)
         !$acc parallel loop gang vector async
         do i=bs,bf
            v%base(i)=-v%base(i)
         enddo
         !$acc wait
         !$acc end data
#endif
      else
         v%base(bs:bf)=-v%base(bs:bf)
      endif

      end subroutine VecSignChange_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecScalarMult_all_CP8(v,fac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies CP object by scalar factor

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      real(kind=8), intent(in)  :: fac

      call v%mult_terms(fac,1,v%R())

      end subroutine VecScalarMult_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecScalarMult_gen_CP8(v,fac,ri,re)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies CP coefs of terms 'ri' through 're' by scalar factor.
! If the factor is negative then the sign is also changed.

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      integer, intent(in) :: ri,re
      real(kind=8), intent(in) :: fac
      integer :: i

      if (v%live_on_device) then
#if ACC_ENABLED
         !$acc data present(v%coef)
         !$acc parallel loop gang vector async
         do i=ri,re
            v%coef(i)=abs(fac)*v%coef(i)
         enddo
         !$acc wait
         !$acc end data
#endif
      else
         v%coef(ri:re)=abs(fac)*v%coef(ri:re)
      endif

      IF (fac.lt.0.d0) call v%changesign(ri,re)

      end subroutine VecScalarMult_gen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Sum_CP8(v,vfac,w,wfac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums vfac*v and wfac*w in CP-format. The summed vector replaces v;
! w is unchanged.

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8), INTENT(IN)     :: w
      real(kind=8), intent(in)  :: vfac,wfac
      integer :: rv,rw

      IF (.NOT.v%same(w)) THEN
         write(*,*) 'Dimensions or type of v and w do not match'
         write(*,*) 'v:'
         call v%show()
         write(*,*) 'w:'
         call w%show()
         call AbortWithError('Error in sum_cp8()')
      ENDIF

      rv=v%R()
      rw=w%R()

!     If one of the two terms is zero, it is replaced by the other
      IF (rv.eq.1 .and. v%coef(1).eq.0.d0) THEN
         call v%flush
         call v%copyfrom(w)
         call v%mult(wfac)

      ELSE IF (rw.eq.1 .and. w%coef(1).eq.0.d0) THEN
         call v%mult(vfac)

!     General case
      ELSE
         call v%resize(rv+rw)
         call v%copy_terms(w,rv+1,rv+rw,1,rw)
         call v%mult_terms(vfac,1,rv)
         call v%mult_terms(wfac,rv+1,rv+rw)
      ENDIF

      end subroutine Sum_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SumLinearCombination_CP8(v,Q,facs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums a linear combination of vectors in Q and stores in v
! Faster than Sum_CP8 for summing many vectors with coefs since only 
! one allocate is performed

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8), INTENT(IN)     :: Q(:)
      real(kind=8), intent(in)   :: facs(:)
      integer :: i,nbloc,ri,rq

      nbloc=SIZE(facs)
      IF (SIZE(Q).ne.nbloc) THEN
         write(*,*) 'Mismatch in dimension of Q, facs'
         call AbortWithError("Error in SumLinearCombination_CP8()")
      ENDIF

!     Get the rank of F (=sum of ranks in Q)
      rq=0
      DO i=1,nbloc
         rq=rq+Q(i)%R()
      ENDDO

      call v%clone(Q(1),rq)

!     Build v from Q and the factors
      ri=0
      DO i=1,nbloc
         rq=Q(i)%R()
         call v%copy_terms(Q(i),ri+1,ri+rq,1,rq)
         call v%mult_terms(facs(i),ri+1,ri+rq)
         ri=ri+rq
      ENDDO

      end subroutine SumLinearCombination_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewVec_CP8(v,rk,rows,trans)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a new vector in CP-format

      implicit none
      CLASS (CP8) :: v
      LOGICAL, INTENT(IN)  :: trans
      INTEGER, INTENT(IN)  :: rows(:)
      INTEGER, INTENT(IN)  :: rk
      INTEGER, ALLOCATABLE :: cols(:)
      INTEGER :: ndof

      ndof=SIZE(rows)
      ALLOCATE(cols(ndof))
      cols(:)=1

      IF (trans) then ! Column vector
         call v%new(rk,cols,rows)
      ELSE ! Row vector
         call v%new(rk,rows,cols)
      ENDIF

      DEALLOCATE(cols)

      end subroutine NewVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NewSqmat_CP8(v,rk,rows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a CP-format outer product of square matrices

      implicit none
      CLASS (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER   :: ndof

      call v%new(rk,rows,rows)

      end subroutine NewSqmat_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ZeroRef_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP of rank 1 from reference CP

      implicit none
      CLASS (CP8), intent(inout) :: v
      TYPE (CP8), intent(in) :: w

      call v%clone(w,1)
      call v%zero

      end subroutine ZeroRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ZeroGen_CP8(v,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP of rank 1 from row/col dims

      implicit none
      CLASS (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:), cols(:)

      call v%new(1,rows,cols)
      call v%zero()

      end subroutine ZeroGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RandomRef_CP8(v,w,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries, with dimensions of reference w

      implicit none
      CLASS (CP8) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk
      INTEGER :: rv,d,ndof,ms,mf
      REAL(kind=8)  :: fac

      IF (present(rk)) THEN
         rv=rk
      ELSE
         rv=w%R()
      ENDIF

      ndof=w%D()

!     Generate v with random entries and equal coefs for all terms
      call v%clone(w,rv)
      v%coef(:)=1.d0/sqrt(REAL(rv))

!     Shift, scale entries for each mode to make rms norm ~ unity
      DO d=1,ndof
         ms=v%ms(d)
         mf=v%mf(d)
         call random_number(v%base(ms:mf))
         fac=sqrt(12.d0/REAL(v%nbas(d)))
         v%base(ms:mf)=fac*(v%base(ms:mf)-0.5d0)
      ENDDO

      end subroutine RandomRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RandomGen_CP8(v,rk,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries. No normalization is done.

      implicit none
      CLASS (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:),cols(:)
      INTEGER, INTENT(IN) :: rk

      call v%new(rk,rows,cols)
      call random_number(v%base)
      v%base=v%base-0.5d0
      v%coef(:)=1.d0

      end subroutine RandomGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine RandomExtend_CP8(v,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments v to rank rk by adding random terms with small coefficients

      implicit none
      CLASS (CP8) :: v
      TYPE (CP8)  :: w
      integer, intent(in) :: rk
      integer :: rkv
      real(kind=8), parameter :: smallnr=1.d-12

      rkv=v%R()
      if (rkv.ge.rk) return      
      call w%clonerand(v,rk-rkv)
      call v%sumcp(1.d0,w,smallnr)
      call w%flush()

      end subroutine RandomExtend_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine IdentityRef_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP of rank 1 from reference CP

      implicit none
      CLASS (CP8), intent(inout) :: v
      TYPE (CP8), intent(in) :: w

      call v%identity(w%rows,w%cols)

      end subroutine IdentityRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine IdentityMatrix_CP8(v,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets rank-1 CP outer-product-of-identity-matrices

      implicit none
      CLASS (CP8) :: v
      integer, intent(in) :: rows(:),cols(:)
      integer :: d,i,l,ndof,m

      call v%new0(rows,cols)
      ndof=v%D()

!     Get an identity matrix for each DOF
      DO d=1,ndof
         m=min(v%M(d),v%N(d))
         l=v%look(1,d)
         DO i=1,m
            call v%put(i,i,1,d,1.d0)
         ENDDO
      ENDDO
      v%coef(1)=1.d0

      end subroutine IdentityMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function VectoDiagMatrix_CP8(v) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates square CP-matrix with elements of v on the diagonal

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      TYPE (CP8) :: w
      integer :: d,ndof,r,rk,i,m,vl,wl

      ndof=v%D()
      rk=v%R()

      call w%new(rk,v%nbas,v%nbas)
      call w%zero

      DO d=1,ndof
         m=w%M(d)
         DO r=1,rk
            vl=v%look(r,d)
            wl=w%look(r,d)
!           Copy elements of v to diagonal of w
            DO i=1,v%nbas(d)
               w%base(wl+(i-1)*m+i)=v%base(vl+i)
            ENDDO
         ENDDO
      ENDDO
      w%coef(1:rk)=v%coef(1:rk)

      end function VectoDiagMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixTranspose_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Transposes CP-matrix by reordering base

      implicit none
      CLASS (CP8), intent(inout) :: v
      integer :: d,ndof,r,rk,i,k,l,m,n,nbas,tmp
      real(kind=8), allocatable :: btmp(:)

      rk=v%R()
      ndof=v%D()

      DO d=1,ndof
!        If either dimension is 1, just swap row and col numbers
         IF ((v%rows(d).gt.1) .and. (v%cols(d).gt.1)) THEN
!           Reorder the elements from the temporary array
            nbas=v%nbas(d)
            m=v%rows(d)
            n=v%cols(d)
            allocate(btmp(nbas))
            DO r=1,rk
               l=v%look(r,d)
               btmp(1:nbas)=v%base(l+1:l+nbas)
               DO i=1,m
                  DO k=1,n
                     v%base(l+(i-1)*n+k)=btmp((k-1)*m+i)
                  ENDDO
               ENDDO
            ENDDO
            deallocate(btmp)
         ENDIF
         tmp=v%rows(d)
         v%rows(d)=v%cols(d)
         v%cols(d)=tmp
      ENDDO

      end subroutine MatrixTranspose_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixZeroOffDiag_all_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for all modes

      implicit none
      CLASS (CP8), intent(inout) :: v
      logical, allocatable :: domode(:)

      ALLOCATE(domode(v%D()))
      domode(:)=.TRUE.
      call v%zerooffdiagonalgen(domode)
      DEALLOCATE(domode)

      end subroutine MatrixZeroOffDiag_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixZeroOffDiag_one_CP8(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for mode d

      implicit none
      CLASS (CP8), intent(inout) :: v
      integer, intent(in)  :: d
      logical, allocatable :: domode(:)

      IF ((d.lt.1) .or. (d.gt.v%D())) THEN
         write(*,*) 'mode d (',d,') must be in range: [1,',v%D(),']'
         call AbortWithError('CPMatrixZeroOffDiag(): d out of range')
      ENDIF

      ALLOCATE(domode(v%D()))
      domode(:)=.FALSE.
      domode(d)=.TRUE.
      call v%zerooffdiagonalgen(domode)
      DEALLOCATE(domode)

      end subroutine MatrixZeroOffDiag_one_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixZeroOffDiag_gen_CP8(v,domode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for each mode matrix if domode
! is set. Mode matrix does not need to be square.

      implicit none
      CLASS (CP8), intent(inout) :: v
      logical, intent(in) :: domode(:)
      integer :: d,ndof,r,rk,i,k,l,m,n

      ndof=v%D()
      rk=v%R()

      IF (SIZE(domode).ne.ndof) THEN
         write(*,*) 'domode has ',SIZE(domode),' entries but ',&
                    'must have ndof (',ndof,') entries'
         call AbortWithError('CPMatrixZeroOffDiag(): ndof mismatch')
      ENDIF

      DO d=1,ndof
         IF (domode(d)) THEN
            m=v%rows(d)
            n=v%cols(d)
            DO r=1,rk
               l=v%look(r,d)
               DO k=1,n
                  DO i=1,m
                     IF (k.ne.i) v%base(l+(k-1)*m+i)=0.d0
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      end subroutine MatrixZeroOffDiag_gen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TrimZeros_CP8(v,tol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resize vector v by trimming columns with coefs smaller than tol
! If all columns are zero, then a rank-1 zero vector is returned

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      TYPE (CP8) :: w
      real(kind=8), intent(in)   :: tol
      integer, allocatable :: iok(:)
      integer :: i,nok,rk

      rk=v%R()
      ALLOCATE(iok(rk))

!     Check the cols for nonzero coef
      nok=0
      DO i=1,rk
         IF (abs(v%coef(i)).gt.abs(tol)) THEN
            nok=nok+1
            iok(nok)=i
         ENDIF
      ENDDO

      IF (nok.lt.rk) THEN
         IF (nok.gt.0) THEN
            call w%clone(v,nok)
            DO i=1,nok
               call w%copy_terms(v,i,i,iok(i),iok(i))
            ENDDO
         ELSE
            call w%clone0(v)
         ENDIF
         call v%replace(w)
      ENDIF

      DEALLOCATE(iok)

      end subroutine TrimZeros_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefSmallest_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out coefficient using base for DOF with the smallest
! basis. Coefficients are then set to 1.0

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      INTEGER :: imode(1)

!     Pick the mode with the smallest basis (fewest multiplies)
      imode=MINLOC(v%nbas)
      call v%multoutcoefmode(imode(1))

      end subroutine MultOutCoefSmallest_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefbyMode_CP8(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies base of mode 'd' by coefficients; resets coefficients to 1

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      integer, intent(in) :: d
      integer :: r,rk,ndof,bs,bf

      ndof=v%D()
      rk=v%R()

!     Error checking
      IF ((d.lt.1) .or. (d.gt.ndof)) THEN
         write(*,*) 'Mode ',d,' must be in range [1,',ndof,']'
         call AbortWithError('MultOutCoefbyMode(): d out of range')
      ENDIF

!     Multiply out coefficients
      DO r=1,rk
         bs=v%BS(r,d)
         bf=v%BF(r,d)
         v%base(bs:bf)=v%coef(r)*v%base(bs:bf)
         v%coef(r)=1.d0
      ENDDO

      end subroutine MultOutCoefbyMode_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DistributeCoef_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies coef by base for all DOFs, and resets coefs to 1.0
! Each mode's base is scaled by the mode size

      implicit none
      CLASS (CP8), INTENT(INOUT) :: v
      REAL(kind=8), ALLOCATABLE :: pows(:)
      REAL(kind=8) :: fac,div
      INTEGER :: d,ndof,r,rk,bs,bf

      ndof=v%D()
      rk=v%R()

!     Get the exponents for modes with different nbas values
      allocate(pows(ndof))
      div=0.d0
      DO d=1,ndof
         pows(d)=sqrt(REAL(v%nbas(d)))
         div=div+pows(d)
      ENDDO
      pows(:)=pows(:)/div

!     Multiply out coefficient
      DO r=1,rk
         DO d=1,ndof
            bs=v%BS(r,d)
            bf=v%BF(r,d)
            fac=v%coef(r)**pows(d)
            v%base(bs:bf)=fac*v%base(bs:bf)
         ENDDO
         v%coef(r)=1.d0
      ENDDO

      deallocate(pows)

      end subroutine DistributeCoef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_1_CP8(v,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real(kind=8), intent(out) :: val

      IF (v%R().ne.1) THEN
         write(*,*) 'Rank of v (',v%R(),') must be 1'
         call AbortWithError('Error in GetRank1DominantEntry_1()')
      ENDIF

      call v%largestentryinterm(1,rowi,coli,val)
!      call GetRank1DominantEntry_gen_CP8(v,1,rowi,coli,val)

      end subroutine GetRank1DominantEntry_1_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_gen_CP8(v,irk,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer, intent(in) :: irk
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real(kind=8), intent(out) :: val
      integer :: d,ndof,bs,bf
      integer :: imx(1)
      real(kind=8) :: vmx

      IF (irk.lt.1 .or. irk.gt.v%R()) THEN
         write(*,*) 'irk (',irk,') out of range: [1,',v%R(),']'
         call AbortWithError('Error in GetRank1DominantEntry_gen()')
      ENDIF

      ndof=v%D()
      ALLOCATE(rowi(ndof),coli(ndof))

      val=v%coef(irk)
      DO d=1,ndof
         bs=v%BS(irk,d)
         bf=v%BF(irk,d)
!        Use imx (range [1:v%nbas(d)]) to extract row,col indices
         imx=MAXLOC(ABS(v%base(bs:bf)))
         rowi(d)=mod(imx(1)-1,v%M(d))+1
         coli(d)=(imx(1)-1)/v%M(d)+1
         vmx=v%base(imx(1)-1+bs)
         val=val*vmx
      ENDDO

      end subroutine GetRank1DominantEntry_gen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ExtractVec_CP8(v,M,indx,getcol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-vector from CP-matrix M from index list, returns as v.
! Extracts a column by default, set getcol=.FALSE. to extract a row

      implicit none
      CLASS (CP8), INTENT(OUT) :: v
      TYPE (CP8), INTENT(IN)  :: M
      integer, intent(in)  :: indx(:)
      logical, intent(in)  :: getcol
      integer, allocatable :: one(:)
      integer :: r,d,ndof,rk,i,p

      ndof=M%D()
      rk=M%R()

      allocate(one(ndof))
      one(:)=1
      IF (getcol) THEN
         call v%new(rk,M%rows,one)
      ELSE
         call v%new(rk,one,M%cols)
      ENDIF
      deallocate(one)

!     Extract base from M
      DO d=1,ndof
!        Error checking
         IF ((indx(d).lt.1) .or. &
             (getcol.and.(indx(d).gt.M%N(d))) .or. &
             ((.not.getcol).and.(indx(d).gt.M%M(d)))) THEN
             IF (getcol) THEN
                write(*,*) 'indx(',d,') = ',indx(d),&
                           ' must be in [1,',M%N(d),']'
             ELSE
                write(*,*) 'indx(',d,') = ',indx(d),&
                           ' must be in [1,',M%M(d),']'
             ENDIF
             call AbortWithError('ExtractCPvec(): index out of range')
         ENDIF

         IF (getcol) THEN
            p=M%M(d)
            DO r=1,rk
               DO i=1,p
                  call v%put(i,1,r,d,M%get(i,indx(d),r,d))
               ENDDO
            ENDDO
         ELSE
            p=M%N(d)
            DO r=1,rk
               DO i=1,p
                  call v%put(1,i,r,d,M%get(indx(d),i,r,d))
               ENDDO
            ENDDO
         ENDIF
      END DO

!     Coefs copy directly
      v%coef(:)=M%coef(:)

      end subroutine ExtractVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ExtractDiagfromMatrix_CP8(w,v,trans)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts diagonal from CP-matrix (must be square); returns a column
! vector. Set trans to true to get a row vector.

      implicit none
      CLASS (CP8), INTENT(OUT) :: w
      TYPE (CP8), INTENT(IN) :: v
      logical, intent(in)  :: trans
      integer, allocatable :: one(:)
      integer :: d,ndof,r,rk,j

      rk=v%R()
      ndof=v%D()

      allocate(one(ndof))
      one(:)=1

      IF (trans) THEN
         call w%new(rk,one,v%cols)
      ELSE
         call w%new(rk,v%rows,one)
      ENDIF

!     Base is taken from diagonal elements of square matrix
      DO d=1,ndof
         IF (v%rows(d).ne.v%cols(d)) THEN
            write(*,*) 'Matrix for mode ',d,' is (',v%rows(d),' x ',&
            v%cols(d),') but must be square'
            call AbortWithError('ExtractDiagfromCPMatrix: matrix dims')
         ENDIF
         DO r=1,rk
            IF (trans) THEN
               DO j=1,v%rows(d)
                  call w%put(1,j,r,d,v%get(j,j,r,d))
               ENDDO
            ELSE
               DO j=1,v%rows(d)
                  call w%put(j,1,r,d,v%get(j,j,r,d))
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     Copy coefs directly
      w%coef(:)=v%coef(:)

      deallocate(one)

      end subroutine ExtractDiagfromMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSubmatrix_CP8(W,V,irs,ics)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Puts CP-submatrix V into W, overwriting that portion of W, from the
! row, col starting indices provided

      implicit none
      CLASS (CP8), INTENT(INOUT) :: W
      TYPE (CP8), INTENT(IN) :: V
      integer, intent(in)  :: irs(:),ics(:)
      integer, allocatable :: irf(:),icf(:)
      integer :: d,ndof,r,rk,i,m,k,n

      ndof=V%D()
      rk=V%R()

!     Error checking extravaganza!
      IF (W%D().ne.ndof) &
         call AbortWithError('PutCPsubmatrix(): M,V must have = ndof')
      IF (W%R().ne.rk) &
         call AbortWithError('PutCPsubmatrix(): M,V must have = rank')

      IF (SIZE(irs).ne.ndof .or. SIZE(ics).ne.ndof ) THEN
          write(*,*) 'ndof of V = ',ndof,'; must equal ndof of ALL of',&
          ' irs,ics, which = ',SIZE(irs),SIZE(ics)
          call AbortWithError('PutCPsubmatrix(): bad ranges ndof')
      ENDIF

      allocate(irf(ndof),icf(ndof))
      irf(:)=irs(:)+V%rows(:)-1
      icf(:)=ics(:)+V%cols(:)-1

      DO d=1,ndof
         IF (irs(d).lt.1 .or. irf(d).gt.W%M(d) .or. &
             irs(d).gt.irf(d)) THEN
            write(*,*) 'Row ranges: [irs(',d,'),irf(',d,')] = [',&
            irs(d),',',irf(d),'] must be in range [1,',W%M(d),']'
            call AbortWithError('PutCPsubmatrix(): bad rows ranges')
         ENDIF
         IF (ics(d).lt.1 .or. icf(d).gt.W%N(d) .or. &
             ics(d).gt.icf(d)) THEN
            write(*,*) 'Col ranges: [ics(',d,'),icf(',d,')] = [',&
            ics(d),',',icf(d),'] must be in range [1,',W%N(d),']'
            call AbortWithError('PutCPsubmatrix(): bad cols ranges')
         ENDIF
      ENDDO

!     Overwrite the selected portion of M
      W%coef(:)=V%coef(:)
      DO d=1,ndof
         m=V%M(d)
         n=V%N(d)
         DO r=1,rk
            DO k=1,n
               DO i=1,m
                  call W%put(irs(d)+i-1,ics(d)+k-1,r,d,V%get(i,k,r,d))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      deallocate(irf,icf)

      end subroutine PutSubmatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractTensorElement_CP8(W,ir,ic) result(res)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets tensor element of CP-tensor W from list of row, col indices

      implicit none
      CLASS (CP8), INTENT(IN) :: W
      integer, intent(in) :: ir(:),ic(:)
      integer :: d,ndof,r,rk
      real(kind=8) :: res
      real(kind=8), allocatable :: prod(:)
      real(kind=8), parameter   :: s=1/sqrt(2.d0)

      ndof=W%D()
      rk=W%R()

!     Error checking: index list sizes, indices in range
      IF (SIZE(ir).ne.ndof) THEN
         write(*,*) 'M has ',ndof,' modes, but "ir" array contains ',&
                    SIZE(ir),' indices'
         call AbortWithError('ExtractCPmatrixelement(): wrong # indices')
      ENDIF
      IF (SIZE(ic).ne.ndof) THEN
         write(*,*) 'M has ',ndof,' modes, but "ic" array contains ',&
                    SIZE(ic),' indices'
         call AbortWithError('ExtractCPmatrixelement(): wrong # indices')
      ENDIF
      DO d=1,ndof
         IF (ir(d).lt.1 .or. ir(d).gt.W%M(d)) THEN
            write(*,*) 'row index ir(',d,') = ',ir(d),&
                       ' is outside of range 1:',W%M(d)
            call AbortWithError('GetCPmatrixelement(): bad row index')
         ENDIF
         IF (ic(d).lt.1 .or. ic(d).gt.W%N(d)) THEN
            write(*,*) 'col index ic(',d,') = ',ic(d),&
                       ' is outside of range 1:',W%N(d)
            call AbortWithError('GetCPmatrixelement(): bad col index')
         ENDIF
      ENDDO

!     Compute the products for each term
      ALLOCATE(prod(rk))
      prod(:)=W%coef(:)

      DO d=1,ndof
         DO r=1,rk
            prod(r)=prod(r)*W%get(ir(d),ic(d),r,d)
         ENDDO
      ENDDO

!     Now accumulate the sum-over-term-products
      res=SUM(prod)

      DEALLOCATE(prod)

      end function ExtractTensorElement_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ModeJoin_CP8(u,v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Merges v and w by joining into a tensor of same rank with the sum of
! the dimensions

      implicit none
      CLASS(CP8), INTENT(INOUT) :: u
      TYPE (CP8), INTENT(IN)    :: v,w
      integer, allocatable :: rows(:),cols(:)
      integer :: i,j,rk,ndof,ndofv,ndofw

      IF (v%R().ne.w%R()) THEN
         write(*,*) 'Ranks of v,w (',v%R(),',',w%R(),&
                    ') must match'
         call AbortWithError('CPModeJoin(): v,w rank mismatch')
      ENDIF

!     Construct u with same rank and sum of mode sizes from v and w
      rk=v%R()
      ndofv=v%D()
      ndofw=w%D()
      ndof=ndofv+ndofw

      ALLOCATE(rows(ndof),cols(ndof))

      rows(1:ndofv)=v%rows(:)
      rows(ndofv+1:ndofv+ndofw)=w%rows(:)
      cols(1:ndofv)=v%cols(:)
      cols(ndofv+1:ndofv+ndofw)=w%cols(:)

      call u%new(rk,rows,cols)
      u%coef(:)=v%coef(:)*w%coef(:)
      DO j=1,ndofv
         u%base(1:SIZE(v%base))=v%base(:)
      ENDDO
      DO j=1,ndofw
         u%base(SIZE(v%base)+1:)=w%base(:)
      ENDDO

      DEALLOCATE(rows,cols)

      end subroutine ModeJoin_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine EntrywiseCompare_CP8(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares F and G, entry by entry.
! WARNING: runtime is exponential in ndof

      implicit none
      CLASS (CP8), INTENT(IN) :: F
      TYPE (CP8), INTENT(IN) :: G
      integer, allocatable  :: indx(:),inmx(:),ir(:),ic(:)
      real(kind=8)  :: valf,valg,vdif,mdif,norm,div
      integer :: d,ndof
      character(len=64) :: frmt

      ndof=F%D()

!     Make sure F and G have same structure
      IF (.not. F%same(G)) THEN
         write(*,*) 'F and G must have same shape'
         write(*,*) "F:"
         call F%show()
         write(*,*) "G:"
         call G%show()
         call AbortWithError("Error in EntrywiseCompare_CP8")
      ENDIF

!     Index values for NextIndex
      allocate(indx(2*ndof),inmx(2*ndof),ir(ndof),ic(ndof))
      div=1.d0
      DO d=1,ndof
         inmx(2*d-1)=F%M(d)
         inmx(2*d)=F%N(d)
         div=div*REAL(F%M(d))*REAL(F%N(d))
      ENDDO

      write(*,*)
      write(*,*) '*** Element-wise compare: F,G,delta ***'
      write(*,*)
      write(frmt,'(A,I0,A)') '(',ndof,'(I0,X,I0,2X),3(f16.8,2X))'
      norm=0.d0
      mdif=0.d0
      indx(:)=1
      DO
         DO d=1,ndof
            ir(d)=indx(2*d-1)
            ic(d)=indx(2*d)
         ENDDO
         valf=F%tensorelement(ir,ic)
         valg=G%tensorelement(ir,ic)
         vdif=valf-valg
         IF (abs(vdif).gt.abs(mdif)) mdif=vdif
         norm=norm+vdif**2
         write(*,frmt) (ir(d),ic(d),d=1,ndof),valf,valg,vdif
         call NextIndex(indx,inmx)
         IF (ALL(indx.eq.1)) EXIT
      ENDDO
      norm=sqrt(norm/div)

      write(*,*) 'RMS diff = ',norm
      write(*,*) 'MAX diff = ',mdif
      write(*,*)
      write(*,*) '****** Element-wise compare done ******'

      deallocate(indx,inmx,ir,ic)

      end subroutine EntrywiseCompare_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! OpenACC functions for GPU
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine createondevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_host) then
         !$acc enter data create(v,v%coef,v%base) &
         !$acc copyin(v%look,v%nbas,v%rows,v%cols)
         v%live_on_device=.TRUE.
      else
         call &
         AbortWithError("createondevice_CP8: v not allocated on host")
      endif
#endif

      end subroutine createondevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine copyintodevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_host) then
         !$acc enter data create(v) copyin(v%coef,v%base) &
         !$acc copyin(v%look,v%nbas,v%rows,v%cols)
         v%live_on_device=.TRUE.
      else
         call &
         AbortWithError("copyintodevice_CP8: v not allocated on host")
      endif
#endif

      end subroutine copyintodevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine updateondevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_device) then
         !$acc update device(v%coef,v%base)
!      else
!         call AbortWithError("updateondevice_CP8: v is not on device")
      endif
#endif

      end subroutine updateondevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine updatefromdevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_device) then
         !$acc update self(v%coef,v%base)
!      else
!         call AbortWithError("updatefromdevice_CP8: v is not on device")
      endif
#endif

      end subroutine updatefromdevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine deletefromdevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_device) then
         !$acc exit data delete(v,v%coef,v%base,v%look,v%nbas,v%rows,v%cols)
         v%live_on_device=.FALSE.
!      else
!         call AbortWithError("deletefromdevice_CP8: v is not on device")
      endif
#endif

      end subroutine deletefromdevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine copyoutfromdevice_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: v

#if ACC_ENABLED
      if (v%live_on_device) then
         !$acc exit data copyout(v%coef,v%base) delete(v,v%look,v%nbas,v%rows,v%cols)
         v%live_on_device=.FALSE.
!      else
!         call AbortWithError("copyoutfromdevice_CP8: v is not on device")
      endif
#endif

      end subroutine copyoutfromdevice_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CP2DtoMat_CP8(v) result(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out 2D CP-format vector 'v' to give matrix 'M'
! Here, both modes 1 and 2 are flattened to give nr1nc1 x nr2nc2 matrix
!!! Consider moving to Reduction.f90

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      REAL(kind=8), ALLOCATABLE :: M(:,:)
      INTEGER :: r,rk,j,k,nr,nc
      REAL(kind=8) :: fac

      IF (v%D().ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoMat()')
      ENDIF

      nr=v%MN(1)
      nc=v%MN(2)
      rk=v%R()

      ALLOCATE(M(nr,nc))
      M=0.d0
      DO r=1,rk
         IF (nc.le.nr) THEN
            DO k=1,nc
!               fac=v%coef(i)*v%data(2)%vec(k,i)
               fac=v%coef(r)*v%base(v%look(r,2)+k)
               DO j=1,nr
!                  M(j,k)=M(j,k)+fac*v%data(1)%vec(j,i)
                  M(j,k)=M(j,k)+fac*v%base(v%look(r,1)+j)
               ENDDO
            ENDDO
         ELSE
            DO j=1,nr
!               fac=v%coef(i)*v%data(1)%vec(j,i)
               fac=v%coef(r)*v%base(v%look(r,1)+j)
               DO k=1,nc
!                  M(j,k)=M(j,k)+fac*v%data(2)%vec(k,i)
                  M(j,k)=M(j,k)+fac*v%base(v%look(r,2)+k)
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      end function CP2DtoMat_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CP2DtoUW_CP8(v,U,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts directional component matrices U and W from 2D CP-format
! vector 'v'
! Arrays for modes 1, 2 are each flattened to give resultant arrays
! U(nr1,nc1 x rk), W(nr2,nc2 x rk)
!!! Consider moving to Reduction.f90

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      REAL(kind=8), ALLOCATABLE, INTENT(OUT) :: U(:,:),W(:,:)
      INTEGER :: r,rk,j,nu,nw,us,uf,ws,wf

      IF (v%D().ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoUW()')
      ENDIF

      nu=v%MN(1)
      nw=v%MN(2)
      rk=v%R()

      ALLOCATE(U(nu,rk),W(nw,rk))

!     Multiply the coef by the base with fewer elements
      DO r=1,rk
        us=v%BS(r,1)
        uf=v%BF(r,1)
        ws=v%BS(r,2)
        wf=v%BF(r,2)
        U(:,r)=v%base(us:uf)
        W(:,r)=v%base(ws:wf) 
        IF (nu.le.nw) THEN
           DO j=1,nu
              U(j,r)=U(j,r)*v%coef(r)
           ENDDO
        ELSE
           DO j=1,nw
              W(j,r)=W(j,r)*v%coef(r)
           ENDDO
        ENDIF
      ENDDO

      end subroutine CP2DtoUW_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Bcast_CP8(v,irank)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts CP-vec v from MPI rank irank to other ranks

      implicit none
      TYPE (CP8) :: v
      integer, intent(in)  :: irank
      integer, allocatable :: rows(:),cols(:)
      integer :: rk,ndof,d
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPr8_Module()

      call CPU_TIME(ti1)

      IF (mpirank.eq.irank) THEN
         rk=v%R()
         ndof=v%D()
      ENDIF

      call bcast(rk,irank)
      call bcast(ndof,irank)

      ALLOCATE(rows(ndof),cols(ndof))

      IF (mpirank.eq.irank) THEN
         rows(:)=v%rows(:)
         cols(:)=v%cols(:)
      ENDIF

      call bcast(rows,irank)
      call bcast(cols,irank)

      IF (mpirank.ne.irank) THEN
         call v%flush
         call v%new(rk,rows,cols)
      ENDIF

      call bcast(v%base,irank)
      call bcast(v%coef,irank)

      DEALLOCATE(rows,cols)

      call CPU_TIME(ti2)
      mpi_time=mpi_time+ti2-ti1

      end subroutine Bcast_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MPI_ISendRecv_CP8(v,w,torank,fromrank,req,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sends CP-vec 'v' to overwrite CP-vec 'w' residing on 'torank', and
! receives CP-vec 'w' from 'v' in 'fromrank', overwriting 'w', via non-
! blocking call. CP objects v, w must have the same rank and dimensions.

      implicit none
      TYPE (CP8) :: v,w
      integer, intent(in) :: torank,fromrank,algo
      integer, intent(out) :: req(4)
      integer :: ierrcs,ierrcr,ierrbs,ierrbr,rk,nb

      IF ((.not.v%same(w)) .or. (v%R().ne.w%R())) THEN
         write(*,*) 'v, w must have same dimensions'
         write(*,*) 'v:'
         call v%show()
         write(*,*) 'w:'
         call w%show()
         CALL AbortWithError('Error in MPI_ISendRecv_CP8()')
      ENDIF

      rk=v%R()
      nb=SIZE(v%base)
      req(:)=0

      select case(algo)
         case(0)
            call MPI_ISEND(v%coef(:),rk,mpi_r8,  torank,1,mpi_comm_wd,req(1),ierrcs)
            call MPI_IRECV(w%coef(:),rk,mpi_r8,fromrank,1,mpi_comm_wd,req(2),ierrcr)
            call MPI_ISEND(v%base(:),nb,mpi_r8,  torank,2,mpi_comm_wd,req(3),ierrbs)
            call MPI_IRECV(w%base(:),nb,mpi_r8,fromrank,2,mpi_comm_wd,req(4),ierrbr)
         case(1)
#if ACC_ENABLED
            if (.not.v%live_on_device) call &
               AbortWithError('MPI_ISendRecv_CP8(): v not on device')
            if (.not.w%live_on_device) call &
               AbortWithError('MPI_ISendRecv_CP8(): w not on device')

            !$acc host_data use_device(v%coef,w%coef,v%base,w%base)
            call MPI_ISEND(v%coef(:),rk,mpi_r8,  torank,1,mpi_comm_wd,req(1),ierrcs)
            call MPI_IRECV(w%coef(:),rk,mpi_r8,fromrank,1,mpi_comm_wd,req(2),ierrcr)
            call MPI_ISEND(v%base(:),nb,mpi_r8,  torank,2,mpi_comm_wd,req(3),ierrbs)
            call MPI_IRECV(w%base(:),nb,mpi_r8,fromrank,2,mpi_comm_wd,req(4),ierrbr)
            !$acc end host_data
#else
            write(*,*) 'OpenACC algorithm not enabled'
            call AbortWithError("MPI_ISendRecv_CP8(): wrong algorithm")
#endif
      end select

      if (ierrcs.ne.0) call AbortWithError(&
          'MPI_ISendRecv_NB_CP8(): coef send failed')
      if (ierrcr.ne.0) call AbortWithError(&
          'MPI_ISendRecv_NB_CP8(): coef recv failed')
      if (ierrbs.ne.0) call AbortWithError(&
          'MPI_ISendRecv_NB_CP8(): base send failed')
      if (ierrbr.ne.0) call AbortWithError(&
          'MPI_ISendRecv_NB_CP8(): base recv failed')

      end subroutine MPI_ISendRecv_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MPI_Sync_block_CP8(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     MPI CP-format block of vectors

      IMPLICIT NONE
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      integer, allocatable :: dims(:,:),ranks(:),starts(:),widths(:)
      integer, allocatable :: mvecs(:),moffs(:),mstarts(:),mwidths(:)
      real*8, allocatable  :: cpblock(:)
      integer :: nbloc,ndofs,termlen
      integer :: p,b,i,j,ierr,totlen,ist,nbas,bs,bf
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_CPr8_Module()

      call CPU_TIME(ti1)

      nbloc=SIZE(Q)

!     Number of CP-vecs on each MPI rank, and mpi offsets
      ALLOCATE(mvecs(mpinodes),moffs(mpinodes))
      mvecs(:)=0
      do b=1,nbloc
         i=mod(b-1,mpinodes)+1
         mvecs(i)=mvecs(i)+1
      enddo
      moffs(1)=0
      do p=2,mpinodes
         moffs(p)=moffs(p-1)+mvecs(p-1)
      enddo

!     Array of CP-ranks depends on the CP-ranks of vectors distributed
!     over different MPI-ranks
      ALLOCATE(ranks(nbloc))
      ranks(:)=0
      do i=1,mvecs(mpirank+1)
         b=moffs(mpirank+1)+i
         ranks(b)=Q(b)%R()
      enddo

!     Gather the list of ranks
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                          ranks,mvecs,moffs,mpi_i4,mpi_comm_wd,ierr)

!     Extract size information from the vecs
      if (mpirank.eq.0) then
         ndofs=Q(1)%D()
      endif
      call bcast(ndofs)

!     Arrays of rows, cols, are same for all vectors in the block
      ALLOCATE(dims(ndofs,2))
      if (mpirank.eq.0) then
         dims(:,1)=Q(1)%rows(:)
         dims(:,2)=Q(1)%cols(:)
      endif
      call bcast(dims)

!     termlen is the number of elements (including coefficient) of a rank-1 CP-vec
      termlen=1 ! Start at 1 for coef
      do i=1,ndofs
         termlen=termlen+dims(i,1)*dims(i,2)
      enddo

!     Compute offsets and widths for CP-vecs in the big array
      ALLOCATE(starts(nbloc),widths(nbloc))
      starts(1)=0
      widths(1)=ranks(1)*termlen
      do b=2,nbloc
         starts(b)=starts(b-1)+widths(b-1)
         widths(b)=ranks(b)*termlen
      enddo
      totlen=starts(nbloc)+widths(nbloc)

!     Compute offsets and widths for MPI blocks in the big array
      ALLOCATE(mstarts(mpinodes),mwidths(mpinodes))
      mstarts(:)=0
      do p=1,mpinodes
         mwidths(p)=0
         do i=1,mvecs(p)
            b=moffs(p)+i
            if (i.eq.1) mstarts(p)=starts(b)
            mwidths(p)=mwidths(p)+widths(b)
         enddo
      enddo

!     Pack the base and coefs into the big array
      ALLOCATE(cpblock(totlen))
      cpblock(:)=0.d0
      do p=1,mvecs(mpirank+1)
         b=moffs(mpirank+1)+p
!        Copy coefs to big array first
         ist=starts(b)
         cpblock(ist+1:ist+ranks(b))=Q(b)%coef(1:ranks(b))
         ist=ist+ranks(b)
!        Copy factor matrices to big array
         do j=1,ndofs
            nbas=dims(j,1)*dims(j,2)*ranks(b)
            bs=Q(b)%BS(1,j)
            bf=Q(b)%BF(ranks(b),j)
            cpblock(ist+1:ist+nbas)=Q(b)%base(bs:bf)
            ist=ist+nbas
         enddo
      enddo

!     Sync the big array over all MPI ranks
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,cpblock,&
                          mwidths,mstarts,mpi_r8,mpi_comm_wd,ierr)

!     Reconstruct the block of CP vectors on all MPI ranks
      do b=1,nbloc
         call Q(b)%flush
         call Q(b)%new(ranks(b),dims(:,1),dims(:,2))
!        Copy coefs from big array
         ist=starts(b)
         Q(b)%coef(1:ranks(b))=cpblock(ist+1:ist+ranks(b))
         ist=ist+ranks(b)
!        Copy factor matrices from big array
         do j=1,ndofs
            nbas=dims(j,1)*dims(j,2)*ranks(b)
            bs=Q(b)%BS(1,j)
            bf=Q(b)%BF(ranks(b),j)
            Q(b)%base(bs:bf)=cpblock(ist+1:ist+nbas)
            ist=ist+nbas
         enddo
      enddo

      deallocate(dims,cpblock)
      deallocate(mvecs,moffs,ranks,starts,widths,mstarts,mwidths)

      call CPU_TIME(ti2)
      mpi_time=mpi_time+ti2-ti1

      end subroutine MPI_Sync_block_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Helper functions for debug
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine thistoCP_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      TYPE (CP) :: w
      integer   :: rk,ndof,r,d,bs,bf,gst,n
      
      rk=v%R()
      ndof=v%D()

      w=NewCP(rk,v%rows,v%cols)

!     Copy CP8 data to CP type
      w%coef(:)=v%coef(:)

      gst=0
      do d=1,ndof
         n=w%nbas(d)
         do r=1,rk
            bs=v%BS(r,d)
            bf=v%BF(r,d)
            w%base(gst+1:gst+n,r)=v%base(bs:bf)
         enddo
         gst=gst+n
      enddo

      end subroutine thistoCP_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine thisfromCP_CP8(w,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      CLASS (CP8) :: w
      TYPE (CP), INTENT(IN) :: v
      integer   :: rk,ndof,r,d,bs,bf,gst,n

      rk=v%R()
      ndof=v%D()

      call w%new(rk,v%rows,v%cols)

!     Copy CP8 data to CP type
      w%coef(:)=v%coef(:)

      gst=0
      do d=1,ndof
         n=v%nbas(d)
         do r=1,rk
            bs=w%BS(r,d)
            bf=w%BF(r,d)
            w%base(bs:bf)=v%base(gst+1:gst+n,r)
         enddo
         gst=gst+n
      enddo

      end subroutine thisfromCP_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
