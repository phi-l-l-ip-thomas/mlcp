!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains structures needed for CP representation, real*8

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT NONE

      TYPE POINTERMAT8(prec)
         INTEGER, KIND :: prec
         REAL(kind=prec), ALLOCATABLE :: mat(:,:,:)
         REAL(kind=prec), POINTER :: vec(:,:) => null()
      END TYPE POINTERMAT8

      TYPE CP8
         TYPE (POINTERMAT8(8)), ALLOCATABLE :: data(:)
         REAL*8 , ALLOCATABLE :: coef(:)
         INTEGER, ALLOCATABLE :: dims(:,:)
         INTEGER, POINTER :: nbas(:) => null()
         INTEGER, POINTER :: ibas(:) => null(), fbas(:) => null()
         INTEGER, POINTER :: rows(:) => null(), cols(:) => null()
         CONTAINS
            PROCEDURE :: R => GetRank_CP8
            PROCEDURE :: D => GetNdof_CP8
            PROCEDURE :: M => GetRows_CP8
            PROCEDURE :: N => GetCols_CP8
            PROCEDURE :: show => ShowStats_CP8
            PROCEDURE :: printvec => PrintVec_CP8
            PROCEDURE :: ok => checkcoefs_CP8
      END TYPE CP8

      INTERFACE New_CP8
         MODULE PROCEDURE NewGen_CP8,NewRef_CP8,NewVec_CP8,NewSqmat_CP8
      END INTERFACE New_CP8

      INTERFACE Zero_CP8
         MODULE PROCEDURE ZeroGen_CP8,ZeroRef_CP8
      END INTERFACE Zero_CP8

      INTERFACE Random_CP8
         MODULE PROCEDURE RandomGen_CP8,RandomRef_CP8
      END INTERFACE Random_CP8

      INTERFACE PrintMat_CP8
         MODULE PROCEDURE PrintMat_all_CP8,PrintMat_gen_CP8
      END INTERFACE PrintMat_CP8

      INTERFACE Copy_CP8
         MODULE PROCEDURE Copy_all_CP8, ExtractSubmatrix_CP8
      END INTERFACE Copy_CP8

      INTERFACE MatrixZeroOffDiag_CP8
         MODULE PROCEDURE MatrixZeroOffDiag_all_CP8
         MODULE PROCEDURE MatrixZeroOffDiag_one_CP8
         MODULE PROCEDURE MatrixZeroOffDiag_gen_CP8
      END INTERFACE MatrixZeroOffDiag_CP8

      INTERFACE MultOutCoef_CP8
         MODULE PROCEDURE MultOutCoefSmallest_CP8,MultOutCoefbyMode_CP8
      END INTERFACE MultOutCoef_CP8

      INTERFACE GetRank1DominantEntry_CP8
         MODULE PROCEDURE GetRank1DominantEntry_gen_CP8
         MODULE PROCEDURE GetRank1DominantEntry_1_CP8
      END INTERFACE GetRank1DominantEntry_CP8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewPointerMat8(rk,row,col) result(P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Create a 2D array pointing to data in M

      implicit none
      TYPE (POINTERMAT8(8)), TARGET :: P
      integer, intent(in) :: rk,row,col
      integer :: i

      ALLOCATE(P%mat(row,col,rk))
!     'vec' points to each row x col portion as a 1D array
      do i=1,rk
         P%vec(1:row*col,1:rk) => P%mat(:,:,:)
      enddo

      end function NewPointerMat8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushPointerMat8(P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Deallocate a 2D array

      implicit none
      TYPE (POINTERMAT8(8)) :: P

      DEALLOCATE(P%mat)
      P%vec=>null()

      end subroutine FlushPointerMat8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewGen_CP8(rk,rows,cols) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! General subroutine for initializing outer product of CP-vectors or
! matrices. Symmetric storage is possible for square matrix factors

      implicit none
      TYPE (CP8), TARGET  :: v
      INTEGER, INTENT(IN) :: rows(:), cols(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER :: ndof,nrdim,i,j

      ndof=SIZE(rows)

!     Check for correct sizes and consistency in dimensions
      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCP()')
      ENDIF

      IF (rk.lt.1) THEN
         write(*,*) 'Error: rk must be at least 1!'
         call AbortWithError('Error in NewCP()')
      ENDIF

      IF (SIZE(cols).ne.ndof) THEN
         write(*,*) 'Error: size of "cols" array differs from ndof!'
         call AbortWithError('Error in NewCP()')
      ENDIF

      DO i=1,ndof
         IF (rows(i).lt.1) THEN
            write(*,*) 'Error: number of rows must be at least 1!'
            call AbortWithError('Error in NewCP()')
         ENDIF

         IF (cols(i).lt.1) THEN
            write(*,*) 'Error: number of cols must be at least 1!'
            call AbortWithError('Error in NewCP()')
         ENDIF
      ENDDO

!     Set arrays containing dimensions
      ALLOCATE(v%dims(ndof,0:2))
      v%dims(:,0)=rows(:)*cols(:)
      v%dims(:,1)=rows(:)
      v%dims(:,2)=cols(:)

!     Assign pointers
      v%nbas(1:ndof) => v%dims(1:ndof,0)
      v%rows(1:ndof) => v%dims(1:ndof,1)
      v%cols(1:ndof) => v%dims(1:ndof,2)

!     Allocate the factor matrices
      ALLOCATE(v%data(ndof),v%coef(rk))
      DO i=1,ndof
         v%data(i)=NewPointerMat8(rk,v%rows(i),v%cols(i))
      ENDDO

      end function NewGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewRef_CP8(w,rk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CP-matrix v, using sizes in w (including the rank, if not
! passed as an optional argument)

      implicit none
      TYPE (CP8) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk

      IF (present(rk)) THEN
         v=New_CP8(rk,w%rows,w%cols)
      ELSE
         v=New_CP8(w%R(),w%rows,w%cols)
      ENDIF

      end function NewRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewVec_CP8(rk,rows,trans) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a new vector in CP-format

      implicit none
      TYPE (CP8) :: v
      LOGICAL, INTENT(IN)  :: trans
      INTEGER, INTENT(IN)  :: rows(:)
      INTEGER, INTENT(IN)  :: rk
      INTEGER, ALLOCATABLE :: cols(:)
      INTEGER :: ndof

      ndof=SIZE(rows)
      ALLOCATE(cols(ndof))
      cols(:)=1

      IF (trans) then ! Column vector
         v=New_CP8(rk,cols,rows)
      ELSE ! Row vector
         v=New_CP8(rk,rows,cols)
      ENDIF

      DEALLOCATE(cols)

      end function NewVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewSqmat_CP8(rk,rows) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a CP-format outer product of square matrices

      implicit none
      TYPE (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER   :: ndof

      v=New_CP8(rk,rows,rows)

      end function NewSqmat_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes CP-format type

      implicit none
      TYPE (CP8) :: v
      integer :: i,j,rk,ndof

      rk=v%R()
      ndof=v%D()

!     Dereference pointers
      DO i=1,ndof
         call FlushPointerMat8(v%data(i))
      ENDDO

!     Deallocate arrays
      IF (ALLOCATED(v%data)) DEALLOCATE(v%data)
      IF (ALLOCATED(v%dims)) DEALLOCATE(v%dims)
      IF (ALLOCATED(v%coef)) DEALLOCATE(v%coef)

      end subroutine Flush_CP8

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
      write(*,'(X,A,I4,A,I0)') 'ndof = ',v%D(),'; rank = ',v%R()
      write(*,*) 'mode [rows x cols] sym'
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
      integer :: r,rk,j,d,i,n
      character*64 :: frmt

      rk=v%R()
      d=v%D()

      write(frmt,'(A,I0,A)') '(A,2X,',rk,'(ES17.10,X))'
      write(*,frmt) ' Vcoef =',(v%coef(r),r=1,rk)
      write(frmt,'(A,I0,A)') '(2(I4),X,',rk,'f18.10)'
      do j=1,d
         n=v%nbas(j)
         do i=1,n
            write(*,frmt) j,i,(v%data(j)%vec(i,r),r=1,rk)
         enddo
      enddo
      write(*,*)

      end subroutine PrintVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMat_all_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format, wrapper for usual case

      implicit none
      CLASS (CP8), INTENT(IN) :: v

      call PrintMat_gen_CP8(v,.TRUE.)

      end subroutine PrintMat_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMat_gen_CP8(v,showbase)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format, general routine

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      logical, intent(in) :: showbase
      integer :: r,rk,j,d,i,n,k,m
      character*64 :: frmt

      rk=v%R()
      d=v%D()

      write(*,*)
      DO r=1,rk
         write(*,'(A,I0,A,ES23.16)') 'RANK: ',r,'; coef = ',v%coef(r)
         IF (showbase) THEN
            DO j=1,d
               write(*,'(/A,I0)') 'dof : ',j
               m=v%M(j)
               n=v%N(j)
               write(frmt,'(A,I0,A)') '(',n,'(X,f14.6))'
               DO i=1,m
                  write(*,frmt) (v%data(j)%mat(i,k,r),k=1,n)
               ENDDO
            ENDDO
            write(*,*)
         ENDIF
      ENDDO
      write(*,*)

      end subroutine PrintMat_gen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRank_CP8(v) result(rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: rk

      IF (ALLOCATED(v%coef)) THEN
         rk=SIZE(v%coef)
      ELSE
         rk=0
      ENDIF

      end function GetRank_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetNdof_CP8(v) result(ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      integer :: ndof

      IF (ALLOCATED(v%data)) THEN
         ndof=SIZE(v%data,1)
      ELSE
         ndof=0
      ENDIF

      end function GetNdof_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetRows_CP8(v,d) result(rows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

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
! Returns number of modes of v

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

      function CHECKNBAS_CP8(v1,v2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      implicit none
      TYPE (CP8), INTENT(IN) :: v1,v2
      LOGICAL :: CHECKNBAS_CP8
      INTEGER :: i

      CHECKNBAS_CP8=.TRUE.

      IF (v1%D().ne.v2%D()) THEN
         CHECKNBAS_CP8=.FALSE.
      ELSE
         DO i=1,v1%D()
            IF ((v1%rows(i).ne.v2%rows(i)) .or. &
                (v1%cols(i).ne.v2%cols(i))) THEN
               CHECKNBAS_CP8=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      end function CHECKNBAS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CHECKCOEFS_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks coefficients for good (non-NaN) values

      implicit none
      CLASS (CP8), INTENT(IN) :: v
      LOGICAL :: CHECKCOEFS_CP8
      INTEGER :: i

      CHECKCOEFS_CP8=.TRUE.

      DO i=1,v%R()
         IF (v%coef(i).ne.v%coef(i)) THEN
            CHECKCOEFS_CP8=.FALSE.
            EXIT
         ENDIF
      ENDDO

      end function CHECKCOEFS_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ZeroRef_CP8(w) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP of rank 1 from reference CP

      implicit none
      TYPE (CP8), intent(in) :: w
      TYPE (CP8) :: v

      v=New_CP8(w,1)
      call SetZero_CP8(v)

      end function ZeroRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ZeroGen_CP8(rows,cols) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP of rank 1 from row/col dims

      implicit none
      TYPE (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:), cols(:)

      v=New_CP8(1,rows,cols)
      call SetZero_CP8(v)

      end function ZeroGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetZero_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros CP-vector

      implicit none
      TYPE (CP8), intent(inout) :: v
      integer :: j,ndof

      ndof=v%D()
      v%coef=0.d0
      do j=1,ndof
         v%data(j)%mat=0.d0
      enddo

      end subroutine SetZero_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RandomRef_CP8(w,rk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries, with dimensions of reference w

      implicit none
      TYPE (CP8) :: v
      TYPE (CP8), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk
      INTEGER :: rv,d,ndof
      REAL*8  :: fac

      IF (present(rk)) THEN
         rv=rk
      ELSE
         rv=w%R()
      ENDIF

      ndof=v%D()

!     Generate v with random entries and equal coefs for all terms
      v=New_CP8(w,rv)
      v%coef(:)=1.d0/sqrt(REAL(rv))

!     Shift, scale entries for each mode to make rms norm ~ unity
      DO d=1,ndof
         call random_number(v%data(d)%vec(:,:))
         fac=sqrt(12.d0/REAL(v%nbas(d)))
         v%data(d)%vec=fac*v%data(d)%vec-0.5d0
      ENDDO

      end function RandomRef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RandomGen_CP8(rk,rows,cols) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries. No normalization is done.

      implicit none
      TYPE (CP8) :: v
      INTEGER, INTENT(IN) :: rows(:),cols(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER :: d,ndof

      ndof=v%D()

      v=New_CP8(rk,rows,cols)
      v%coef(:)=1.d0

      DO d=1,ndof
         call random_number(v%data(d)%vec(:,:))
         v%data(d)%vec=v%data(d)%vec-0.5d0
      ENDDO

      end function RandomGen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function IdentityMatrix_CP8(rows) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets rank-1 CP outer-product-of-identity-matrices

      implicit none
      TYPE (CP8) :: v
      integer, intent(in) :: rows(:)
      integer :: i,j

      v=New_CP8(1,rows)

!     Get an identity matrix for each DOF
      DO j=1,v%D()
         v%data(j)%mat=0.d0
         DO i=1,v%M(j)
            v%data(j)%mat(i,i,1)=1.d0
         ENDDO
      ENDDO
      v%coef=1.d0

      end function IdentityMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function VectoDiagMatrix_CP8(v) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates square CP-matrix with elements of v on the diagonal

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      TYPE (CP8) :: w
      integer :: d,ndof,i,rk,j

      ndof=v%D()
      rk=v%R()

      w=New_CP8(rk,v%nbas,v%nbas)
      w%coef(:)=v%coef(:)

      DO d=1,ndof
         DO i=1,rk
!           Copy elements of v to diagonal of w
            w%data(d)%mat(:,:,i)=0.d0
            DO j=1,v%nbas(d)
               w%data(d)%mat(j,j,i)=v%data(d)%vec(j,i)
            ENDDO
         ENDDO
      ENDDO

      end function VectoDiagMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReplaceVwithW_CP8(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, disposing W afterwards

      implicit none
      TYPE (CP8) :: v,w

      call Flush_CP8(v)
      v=Copy_CP8(w)
      call Flush_CP8(w)

      end subroutine ReplaceVwithW_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Copy_all_CP8(w) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      implicit none
      TYPE (CP8), INTENT(IN) :: w
      TYPE (CP8) :: v
      INTEGER   :: rk

      rk=w%R()
      v=New_CP8(w,rk)
      call GenCopyWtoV_CP8(v,w,1,rk,1,rk)

      end function Copy_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GenCopyWtoV_CP8(v,w,vi,ve,wi,we)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V
! leaving W intact. v must be allocated beforehand.

      implicit none
      TYPE (CP8) :: v,w
      INTEGER, INTENT(IN) :: vi,ve,wi,we
      INTEGER :: rkv,rkw,d,ndof

      rkv=v%R()
      rkw=w%R()
      ndof=w%D()

      IF (.not.CHECKNBAS_CP8(v,w)) THEN
         write(*,*) 'v,w dimension mismatch'
         CALL AbortWithError('Error in GenCopyWtoV()')
      ENDIF

      IF (vi.lt.1 .or. ve.gt.rkv .or. vi.gt.ve .or. &
          wi.lt.1 .or. we.gt.rkw .or. wi.gt.we .or. &
          we-wi.ne.ve-vi) THEN
          write(*,'(2A,6(X,I0))') 'Bad rank indices: ',&
          'vi,ve,rkv,wi,we,rkw =',vi,ve,rkv,wi,we,rkw
          CALL AbortWithError('Error in GenCopyWtoV()')
      ENDIF

      do d=1,ndof
         v%data(d)%vec(:,vi:ve)=w%data(d)%vec(:,wi:we)
      enddo

      end subroutine GenCopyWtoV_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Resize_CP8(v,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resizes vector v, either by truncating at a smaller rank or by
! adding space for extra terms.

      implicit none
      TYPE (CP8), INTENT(INOUT) :: v
      TYPE (CP8) :: w
      INTEGER, INTENT(IN) :: rk
      INTEGER :: rkv

      rkv=v%R()

      IF (rk.lt.1) &
         call AbortWithError('Error in ResizeV(): rk < 1')

      call ReplaceVwithW_CP8(w,v)
      v=New_CP8(w,rk)
      call GenCopyWtoV_CP8(v,w,1,MIN(rkv,rk),1,MIN(rkv,rk))
      call Flush_CP8(w)

      end subroutine Resize_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixTranspose_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Transposes CP-matrix by reordering base

      implicit none
      TYPE (CP8), intent(inout) :: v
      integer :: d,i,j,k,rk,ndof,tmp
      real*8, allocatable :: btmp(:,:)

      rk=v%R()
      ndof=v%D()

      DO d=1,ndof
!        If either dimension is 1, just swap row and col numbers
         IF ((v%rows(d).gt.1) .and. (v%cols(d).gt.1)) THEN
!           Reorder the elements in temporary array
            allocate(btmp(v%rows(d),v%cols(d)))
            DO i=1,rk
               btmp(:,:)=v%data(d)%mat(1:v%rows(d),1:v%cols(d),i)
               call FlushPointerMat8(v%data(d))
               v%data(d)=NewPointerMat8(rk,v%cols(d),v%rows(d))
               DO j=1,v%rows(d)
                  DO k=1,v%cols(d)
                     v%data(d)%mat(k,j,i)=btmp(j,k)
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
      TYPE (CP8), intent(inout) :: v
      logical, allocatable :: domode(:)

      ALLOCATE(domode(v%D()))
      domode(:)=.TRUE.
      call MatrixZeroOffDiag_CP8(v,domode)
      DEALLOCATE(domode)

      end subroutine MatrixZeroOffDiag_all_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixZeroOffDiag_one_CP8(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for mode d

      implicit none
      TYPE (CP8), intent(inout) :: v
      integer, intent(in)  :: d
      logical, allocatable :: domode(:)

      IF ((d.lt.1) .or. (d.gt.v%D())) THEN
         write(*,*) 'mode d (',d,') must be in range: [1,',v%D(),']'
         call AbortWithError('CPMatrixZeroOffDiag(): d out of range')
      ENDIF

      ALLOCATE(domode(v%D()))
      domode(:)=.FALSE.
      domode(d)=.TRUE.
      call MatrixZeroOffDiag_CP8(v,domode)
      DEALLOCATE(domode)

      end subroutine MatrixZeroOffDiag_one_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MatrixZeroOffDiag_gen_CP8(v,domode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for each mode matrix if domode
! is set. Mode matrix does not need to be square.

      implicit none
      TYPE (CP8), intent(inout) :: v
      logical, intent(in) :: domode(:)
      integer :: d,ndof,i,rk,j,k

      ndof=v%D()
      rk=v%R()

      IF (SIZE(domode).ne.ndof) THEN
         write(*,*) 'domode has ',SIZE(domode),' entries but ',&
                    'must have ndof (',ndof,') entries'
         call AbortWithError('CPMatrixZeroOffDiag(): ndof mismatch')
      ENDIF

      DO d=1,ndof
         IF (domode(d)) THEN
            DO i=1,rk           
               DO k=1,v%cols(d)
                  DO j=1,v%rows(d)
                     IF (k.ne.j) v%data(d)%mat(j,k,i)=0.d0
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
      TYPE (CP8), INTENT(INOUT) :: v
      TYPE (CP8) :: w
      real*8, intent(in)   :: tol
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
            w=New_CP8(v,rk)
            DO i=1,nok
               call GenCopyWtoV_CP8(w,v,i,i,iok(i),iok(i))
            ENDDO
         ELSE
            w=Zero_CP8(v)
         ENDIF
         call ReplaceVwithW_CP8(v,w)
      ENDIF

      DEALLOCATE(iok)

      end subroutine TrimZeros_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefSmallest_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out coefficient using base for DOF with the smallest
! basis. Coefficients are then set to 1.0

      implicit none
      TYPE (CP8), INTENT(INOUT) :: v
      INTEGER :: imode(1)

!     Pick the mode with the smallest basis (fewest multiplies)
      imode=MINLOC(v%nbas)
      call MultOutCoefbyMode_CP8(v,imode(1))

      end subroutine MultOutCoefSmallest_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefbyMode_CP8(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies base of mode 'd' by coefficients; resets coefficients to 1

      implicit none
      TYPE (CP8), INTENT(INOUT) :: v
      integer, intent(in) :: d
      integer :: i,ndof,rk

      ndof=v%D()
      rk=v%R()

!     Error checking
      IF ((d.lt.1) .or. (d.gt.ndof)) THEN
         write(*,*) 'Mode ',d,' must be in range [1,',ndof,']'
         call AbortWithError('MultOutCoefbyMode(): d out of range')
      ENDIF

!     Multiply out coefficients
      DO i=1,rk
         v%data(d)%mat(:,:,i)=v%coef(i)*v%data(d)%mat(:,:,i)
         v%coef(i)=1.d0
      ENDDO

      end subroutine MultOutCoefbyMode_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DistributeCoef_CP8(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies coef by base for all DOFs, and resets coefs to 1.0
! Each mode's base is scaled by the mode size

      implicit none
      TYPE (CP8), INTENT(INOUT) :: v
      REAL*8, ALLOCATABLE :: pows(:)
      REAL*8  :: fac,div
      INTEGER :: d,i,ndof,rk

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
      DO i=1,rk
         DO d=1,ndof
            fac=v%coef(i)**pows(d)
            v%data(d)%mat(:,:,i)=fac*v%data(d)%mat(:,:,i)
         ENDDO
         v%coef(i)=1.d0
      ENDDO

      deallocate(pows)

      end subroutine DistributeCoef_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_1_CP8(v,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real*8, intent(out) :: val

      IF (v%R().ne.1) THEN
         write(*,*) 'Rank of v (',v%R(),') must be 1'
         call AbortWithError('Error in GetRank1DominantEntry_1()')
      ENDIF

      call GetRank1DominantEntry_gen_CP8(v,1,rowi,coli,val)

      end subroutine GetRank1DominantEntry_1_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_gen_CP8(v,irk,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      integer, intent(in) :: irk
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real*8, intent(out) :: val
      integer :: d,ndof
      integer :: imx(1)
      real*8  :: vmx

      IF (irk.lt.1 .or. irk.gt.v%R()) THEN
         write(*,*) 'irk (',irk,') out of range: [1,',v%R(),']'
         call AbortWithError('Error in GetRank1DominantEntry_gen()')
      ENDIF

      ndof=v%D()
      ALLOCATE(rowi(ndof),coli(ndof))

      val=v%coef(irk)
      DO d=1,ndof
!        Use imx (range [1:v%nbas(d)]) to extract row,col indices
         imx=MAXLOC(ABS(v%data(d)%vec(:,irk)))
         rowi(d)=mod(imx(1)-1,v%M(d))+1
         coli(d)=(imx(1)-1)/v%M(d)+1
         vmx=v%data(d)%mat(rowi(d),coli(d),irk)
         val=val*vmx
      ENDDO

      end subroutine GetRank1DominantEntry_gen_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractVec_CP8(M,indx,getcol) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-vector from CP-matrix M from index list, returns as v.
! Extracts a column by default, set getcol=.FALSE. to extract a row

      implicit none
      TYPE (CP8), INTENT(IN)  :: M
      TYPE (CP8)  :: v
      integer, intent(in)  :: indx(:)
      logical, intent(in)  :: getcol
      integer, allocatable :: one(:)
      integer :: i,d,ndof,rk

      ndof=M%D()
      rk=M%R()

      allocate(one(ndof))
      one(:)=1
      IF (getcol) THEN
         v=New_CP8(rk,M%rows,one)
      ELSE
         v=New_CP8(rk,one,M%cols)
      ENDIF
      deallocate(one)

!     Extract base from M      
      DO d=1,ndof
!        Error checking
         IF ((indx(d).lt.1) .or. &
             (getcol.and.(indx(d).gt.M%cols(d))) .or. &
             ((.not.getcol).and.(indx(d).gt.M%rows(d)))) THEN
             IF (getcol) THEN
                write(*,*) 'indx(',d,') = ',indx(d),&
                           ' must be in [1,',M%cols(d),']'
             ELSE
                write(*,*) 'indx(',d,') = ',indx(d),&
                           ' must be in [1,',M%rows(d),']'
             ENDIF
             call AbortWithError('ExtractCPvec(): index out of range')
         ENDIF

         IF (getcol) THEN
            DO i=1,rk
               v%data(d)%mat(:,1,i)=M%data(d)%mat(:,indx(d),i)
            ENDDO
         ELSE
            DO i=1,rk
               v%data(d)%mat(1,:,i)=M%data(d)%mat(indx(d),:,i)
            ENDDO
         ENDIF
      END DO

!     Coefs copy directly
      v%coef(:)=M%coef(:)

      end function ExtractVec_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractDiagfromMatrix_CP8(v,trans) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts diagonal from CP-matrix (must be square); returns a column
! vector. Set trans to true to get a row vector.

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      TYPE (CP8) :: w
      logical, intent(in)  :: trans
      integer, allocatable :: one(:)
      integer :: d,ndof,i,rk,j

      rk=v%R()
      ndof=v%D()

      allocate(one(ndof))
      one(:)=1

      IF (trans) THEN
         w=New_CP8(rk,one,v%cols)
      ELSE
         w=New_CP8(rk,v%rows,one)
      ENDIF

!     Base is taken from diagonal elements of square matrix
      DO d=1,ndof
         IF (v%rows(d).ne.v%cols(d)) THEN
            write(*,*) 'Matrix for mode ',d,' is (',v%rows(d),' x ',&
            v%cols(d),') but must be square'
            call AbortWithError('ExtractDiagfromCPMatrix: matrix dims')
         ENDIF
         DO i=1,rk
            IF (trans) THEN
               DO j=1,v%rows(d)
                  w%data(d)%mat(1,j,i)=v%data(d)%mat(j,j,i)
               ENDDO
            ELSE
               DO j=1,v%rows(d)
                  w%data(d)%mat(j,1,i)=v%data(d)%mat(j,j,i)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     Copy coefs directly
      w%coef(:)=v%coef(:)

      deallocate(one)

      end function ExtractDiagfromMatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractSubmatrix_CP8(M,irs,irf,ics,icf) result(V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-submatrix M from list of row, col starting and ending indices

      implicit none
      TYPE (CP8), INTENT(IN) :: M
      TYPE (CP8) :: V
      integer, intent(in)  :: irs(:),irf(:),ics(:),icf(:)
      integer, allocatable :: rows(:),cols(:)
      integer :: d,ndof,i,rk

      ndof=M%D()
      rk=M%R()

!     Loads of error checking
      IF (SIZE(irs).ne.ndof .or. SIZE(irf).ne.ndof .or. &
          SIZE(ics).ne.ndof .or. SIZE(icf).ne.ndof) THEN
          write(*,*) 'ndof of M = ',ndof,'; must equal ndof of ALL of',&
          ' irs,irf,ics,icf, which = ',&
          SIZE(irs),SIZE(irf),SIZE(ics),SIZE(icf)
          call AbortWithError('ExtractCPsubmatrix(): bad ranges ndof')
      ENDIF
      DO d=1,ndof
         IF (irs(d).lt.1 .or. irf(d).gt.M%rows(d) .or. &
             irs(d).gt.irf(d)) THEN
            write(*,*) 'Row ranges: [irs(',d,'),irf(',d,')] = [',&
            irs(d),',',irf(d),'] must be in range [1,',M%rows(d),']'
            call AbortWithError('ExtractCPsubmatrix(): bad rows ranges')
         ENDIF
         IF (ics(d).lt.1 .or. icf(d).gt.M%cols(d) .or. &
             ics(d).gt.icf(d)) THEN
            write(*,*) 'Col ranges: [ics(',d,'),icf(',d,')] = [',&
            ics(d),',',icf(d),'] must be in range [1,',M%cols(d),']'
            call AbortWithError('ExtractCPsubmatrix(): bad cols ranges')
         ENDIF
      ENDDO

      allocate(rows(ndof),cols(ndof))
      rows(:)=irf(:)-irs(:)+1
      cols(:)=icf(:)-ics(:)+1
      V=New_CP8(rk,rows,cols)
      deallocate(rows,cols)

!     Copy the selected portion v <- M
      V%coef(:)=M%coef(:)
      DO d=1,ndof
         V%data(d)%mat(1:v%rows(d),1:v%cols(d),:)=&
         M%data(d)%mat(irs(d):irf(d),ics(d):icf(d),:)
      ENDDO

      end function ExtractSubmatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSubmatrix_CP8(V,irs,ics,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Puts CP-submatrix V into M, overwriting that portion of M, from the 
! row, col starting indices provided

      implicit none
      TYPE (CP8), INTENT(IN) :: V
      TYPE (CP8), INTENT(INOUT) :: M
      integer, intent(in)  :: irs(:),ics(:)
      integer, allocatable :: irf(:),icf(:)
      integer :: d,ndof,i,rk

      ndof=SIZE(V%nbas)
      rk=SIZE(V%coef)

!     Error checking extravaganza!
      IF (M%D().ne.ndof) &
         call AbortWithError('PutCPsubmatrix(): M,V must have = ndof')
      IF (M%R().ne.rk) &
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
         IF (irs(d).lt.1 .or. irf(d).gt.M%rows(d) .or. &
             irs(d).gt.irf(d)) THEN
            write(*,*) 'Row ranges: [irs(',d,'),irf(',d,')] = [',&
            irs(d),',',irf(d),'] must be in range [1,',M%rows(d),']'
            call AbortWithError('PutCPsubmatrix(): bad rows ranges')
         ENDIF
         IF (ics(d).lt.1 .or. icf(d).gt.M%cols(d) .or. &
             ics(d).gt.icf(d)) THEN
            write(*,*) 'Col ranges: [ics(',d,'),icf(',d,')] = [',&
            ics(d),',',icf(d),'] must be in range [1,',M%cols(d),']'
            call AbortWithError('PutCPsubmatrix(): bad cols ranges')
         ENDIF
      ENDDO

!     Copy coefs of V -> M directly, base by reshaping
      M%coef(:)=V%coef(:)
      DO d=1,ndof
         M%data(d)%mat(irs(d):irf(d),ics(d):icf(d),:)=&
         V%data(d)%mat(1:v%rows(d),1:v%cols(d),:)
      ENDDO

      deallocate(irf,icf)

      end subroutine PutSubmatrix_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractMatrixElement_CP8(M,ir,ic)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets element of CP-matrix M from list of row, col indices

      implicit none
      TYPE (CP8), INTENT(IN) :: M
      integer, intent(in) :: ir(:),ic(:)
      integer :: d,ndof,i,rk
      real*8 :: ExtractMatrixElement_CP8
      real*8, allocatable :: prod(:)
      real*8, parameter   :: s=1/sqrt(2.d0)

      ndof=M%D()
      rk=M%R()

!     Error checking: index list sizes, indices in range
      IF (SIZE(ir).ne.ndof) THEN
         write(*,*) 'M has ',ndof,' modes, but "ir" array contains ',&
                    SIZE(ir),' indices'
         call AbortWithError('GetCPmatrixelement(): wrong # indices')
      ENDIF
      IF (SIZE(ic).ne.ndof) THEN
         write(*,*) 'M has ',ndof,' modes, but "ic" array contains ',&
                    SIZE(ic),' indices'
         call AbortWithError('GetCPmatrixelement(): wrong # indices')
      ENDIF
      DO d=1,ndof
         IF (ir(d).lt.1 .or. ir(d).gt.M%rows(d)) THEN
            write(*,*) 'row index ir(',d,') = ',ir(d),&
                       ' is outside of range 1:',M%rows(d)
            call AbortWithError('GetCPmatrixelement(): bad row index')
         ENDIF
         IF (ic(d).lt.1 .or. ic(d).gt.M%cols(d)) THEN
            write(*,*) 'col index ic(',d,') = ',ic(d),&
                       ' is outside of range 1:',M%cols(d)
            call AbortWithError('GetCPmatrixelement(): bad col index')
         ENDIF
      ENDDO

!     Compute the products for each term
      ALLOCATE(prod(rk))
      prod(:)=M%coef(:)

      DO d=1,ndof
         DO i=1,rk
            prod(i)=prod(i)*M%data(d)%mat(ir(d),ic(d),i)           
         ENDDO
      ENDDO

!     Now accumulate the sum-over-term-products
      ExtractMatrixElement_CP8=SUM(prod)

      DEALLOCATE(prod)

      end function ExtractMatrixElement_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ModeJoin_CP8(v,w) result(u)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Merges v and w by joining into a tensor of same rank with the sum of
! the dimensions

      implicit none
      TYPE (CP8), INTENT(IN) :: v,w
      TYPE (CP8) :: u
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

      u=New_CP8(rk,rows,cols)
      u%coef(:)=v%coef(:)*w%coef(:)
      DO j=1,ndofv
         u%data(j)%mat(:,:,:)=v%data(j)%mat(:,:,:)
      ENDDO
      DO j=1,ndofw
         u%data(ndofv+j)%mat(:,:,:)=w%data(j)%mat(:,:,:)
      ENDDO

      DEALLOCATE(rows,cols)

      end function ModeJoin_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine EntrywiseCompare_CP8(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares F and G, entry by entry.
! WARNING: runtime is exponential in ndof

      implicit none
      TYPE (CP8), INTENT(IN) :: F,G
      integer, allocatable  :: indx(:),inmx(:),ir(:),ic(:)
      real*8  :: valf,valg,vdif,mdif,norm,div
      integer :: d,ndof
      character(len=64) :: frmt

      ndof=F%D()

!     Index values for NextIndex
      allocate(indx(2*ndof),inmx(2*ndof),ir(ndof),ic(ndof))
      div=1.d0
      DO d=1,ndof
         inmx(2*d-1)=F%rows(d)
         inmx(2*d)=F%cols(d)
         div=div*REAL(F%rows(d))*REAL(F%cols(d))
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
         valf=ExtractMatrixElement_CP8(F,ir,ic)
         valg=ExtractMatrixElement_CP8(G,ir,ic)
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

      function CP2DtoMat_CP8(v) result(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out 2D CP-format vector 'v' to give matrix 'M'
! Here, both modes 1 and 2 are flattened to give nr1nc1 x nr2nc2 matrix
!!! Consider moving to Reduction.f90

      implicit none
      TYPE (CP8), INTENT(IN) :: v
      REAL*8, ALLOCATABLE :: M(:,:)
      INTEGER :: i,j,k,nr,nc,rk
      REAL*8  :: fac

      IF (v%D().ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoMat()')
      ENDIF

      nr=v%nbas(1)
      nc=v%nbas(2)
      rk=v%R()

      ALLOCATE(M(nr,nc))
      M=0.d0
      DO i=1,rk
         IF (nc.le.nr) THEN
            DO k=1,nc
               fac=v%coef(i)*v%data(2)%vec(k,i)
               DO j=1,nr
                  M(j,k)=M(j,k)+fac*v%data(1)%vec(j,i)
               ENDDO
            ENDDO
         ELSE
            DO j=1,nr
               fac=v%coef(i)*v%data(1)%vec(j,i)
               DO k=1,nc
                  M(j,k)=M(j,k)+fac*v%data(2)%vec(k,i)
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
      REAL*8, ALLOCATABLE, INTENT(OUT) :: U(:,:),W(:,:)
      INTEGER :: i,j,nu,nw,rk

      IF (v%D().ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoUW()')
      ENDIF

      nu=v%nbas(1)
      nw=v%nbas(2)
      rk=v%R()

      ALLOCATE(U(nu,rk),W(nw,rk))
      U(:,:)=v%data(1)%vec(:,:)
      W(:,:)=v%data(2)%vec(:,:)

!     Multiply the coef by the base with fewer elements
      DO i=1,rk
        IF (nu.le.nw) THEN
           DO j=1,nu
              U(j,i)=U(j,i)*v%coef(i)
           ENDDO
        ELSE
           DO j=1,nw
              W(j,i)=W(j,i)*v%coef(i)
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
         call Flush_CP8(v)
         v=New_CP8(rk,rows,cols)
      ENDIF

      DO d=1,ndof
         call bcast(v%data(d)%mat,irank)
      ENDDO
      call bcast(v%coef,irank)

      DEALLOCATE(rows,cols)

      end subroutine Bcast_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE MPI_Sync_block_CP8(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     MPI CP-format block of vectors

      IMPLICIT NONE
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      integer, allocatable :: dims(:,:),ranks(:),starts(:),widths(:)
      integer, allocatable :: mvecs(:),moffs(:),mstarts(:),mwidths(:)
      real*8, allocatable  :: cpblock(:)
      integer :: nbloc,ndofs,termlen
      integer :: p,b,i,j,ierr,totlen,ist,nbas

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
      ALLOCATE(cpblock(0:totlen-1))
      cpblock(:)=0.d0
      do p=1,mvecs(mpirank+1)
         b=moffs(mpirank+1)+p
!        Copy coefs to big array first
         ist=starts(b)
         cpblock(ist+1:ist+ranks(b))=Q(b)%coef(1:ranks(b))
         ist=ist+ranks(b)
!        Copy factor matrices to big array
         do j=1,ndofs
            nbas=dims(j,1)*dims(j,2)
            do i=1,ranks(b)
               cpblock(ist+1:ist+nbas)=Q(b)%data(j)%vec(:,i)
               ist=ist+nbas
            enddo
         enddo
      enddo

!     Sync the big array over all MPI ranks
      call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,cpblock,&
                          mwidths,mstarts,mpi_r8,mpi_comm_wd,ierr)

!     Reconstruct the block of CP vectors on all MPI ranks
      do b=1,nbloc
         call Flush_CP8(Q(b))
         Q(b)=New_CP8(ranks(b),dims(:,1),dims(:,2))
!        Copy coefs from big array
         ist=starts(b)
         Q(b)%coef(1:ranks(b))=cpblock(ist+1:ist+ranks(b))
         ist=ist+ranks(b)
!        Copy factor matrices from big array
         do j=1,ndofs
            nbas=dims(j,1)*dims(j,2)
            do i=1,ranks(b)
               Q(b)%data(j)%vec(:,i)=cpblock(ist+1:ist+nbas)
               ist=ist+nbas
            enddo
         enddo
      enddo

      deallocate(dims,cpblock)
      deallocate(mvecs,moffs,ranks,starts,widths,mstarts,mwidths)

      end subroutine MPI_Sync_block_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
