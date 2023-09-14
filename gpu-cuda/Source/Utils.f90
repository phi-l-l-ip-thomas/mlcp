!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE UTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Utility functions

      USE ERRORTRAP

!     Vector classes (hack for creating ragged multi-dimensional arrays)
      TYPE RVEC
         REAL*8, ALLOCATABLE :: v(:)
         CONTAINS
            PROCEDURE :: new => Newrvec
            PROCEDURE :: flush => Flushrvec
            PROCEDURE :: n => Getrvsize
      END TYPE RVEC

      TYPE IVEC
         INTEGER, ALLOCATABLE :: v(:)
         CONTAINS
            PROCEDURE :: new => Newivec
            PROCEDURE :: flush => Flushivec
            PROCEDURE :: n => Getivsize
      END TYPE IVEC

      INTERFACE OPERATOR( .SEQ. )
          MODULE PROCEDURE CompareStringsForEQ
      END INTERFACE

      INTERFACE Vec2Mat
         MODULE PROCEDURE Vec2Matrix,Vec2SymMat
      END INTERFACE Vec2Mat

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Newrvec(vec,n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of real vec

      implicit none
      CLASS (RVEC), INTENT(INOUT) :: vec
      integer, intent(in) :: n

      ALLOCATE(vec%v(n))

      end subroutine Newrvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flushrvec(vec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of real vec

      implicit none
      CLASS (RVEC), INTENT(INOUT) :: vec

      IF (ALLOCATED(vec%v)) DEALLOCATE(vec%v)

      end subroutine Flushrvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getrvsize(vec) result(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of real vec

      implicit none
      CLASS (RVEC), INTENT(IN) :: vec
      integer :: n

      n=0
      IF (ALLOCATED(vec%v)) n=SIZE(vec%v)

      end function Getrvsize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Newivec(vec,n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of real vec

      implicit none
      CLASS (IVEC), INTENT(INOUT) :: vec
      integer, intent(in) :: n

      ALLOCATE(vec%v(n))

      end subroutine Newivec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flushivec(vec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of real vec

      implicit none
      CLASS (IVEC), INTENT(INOUT) :: vec

      IF (ALLOCATED(vec%v)) DEALLOCATE(vec%v)

      end subroutine Flushivec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getivsize(vec) result(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns size of integer vec

      implicit none
      CLASS (IVEC), INTENT(IN) :: vec
      integer :: n

      n=0
      IF (ALLOCATED(vec%v)) n=SIZE(vec%v)

      end function Getivsize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER FUNCTION LookForFreeUnit()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set the smallest integer equal or greater than 20 that is 
! available as unit i/o. 
! returns the integer number of the smallest free unit.


      IMPLICIT NONE
      LOGICAL               :: File_Opened
      INTEGER, PARAMETER    :: UnitMax = 300

      DO LookForFreeUnit = 20, UnitMax ! Look for non-opened I/O unit
         INQUIRE( UNIT=LookForFreeUnit, OPENED=File_Opened )
         IF (.NOT. File_Opened) EXIT
      END DO

!     If an available unit has been found, use as function results 
      CALL ERROR((LookForFreeUnit.eq.UnitMax), &
                 "PrintTools: No free I/O unit available")

      END FUNCTION LookForFreeUnit

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintWallTime(message)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Print out the current wall time

      implicit none
      character(len=*), intent(in) :: message
      integer :: d(3),t(3)

      call idate(d)
      call itime(t)
      write(*,'(X,A,4(A,I0),2(A,I2.2))') TRIM(ADJUSTL(message)),&
            ' on ',d(2),'/',d(1),'/',d(3),' at ',t(1),':',t(2),':',t(3)

      end subroutine PrintWallTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMatrix(mat)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Print out a matrix

      implicit none
      real*8, intent(in) :: mat(:,:)
      integer   :: ir,ic,rows,cols

      rows=SIZE(mat,1)
      cols=SIZE(mat,2)

      write(*,*)
      do ir=1,rows
!         write(*,'(40(f15.8))') (mat(ir,ic),ic=1,cols)
         write(*,'(40(f13.6))') (mat(ir,ic),ic=1,cols)
      enddo
      write(*,*)

      end subroutine PrintMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintVector(vec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Print out a vector

      implicit none
      real*8, intent(in) :: vec(:)
      integer   :: ir,rows

      rows=SIZE(vec)

      write(*,*)
      do ir=1,rows
         write(*,'(I6,f24.14)') ir,vec(ir)
      enddo
      write(*,*)

      end subroutine PrintVector

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE MatDiag2Utriang(dmat,tmat,sym,refl)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts matrix in stored diagonal form to upper triangular form

      implicit none
      real*8, allocatable, intent(out) :: tmat(:,:)
      real*8, intent(in)  :: dmat(:,:)
      integer, intent(in) :: sym
      logical, intent(in) :: refl
      integer :: i,j,dimr,dimc,dl

      dimr=SIZE(dmat,1)
      dimc=SIZE(dmat,2)
      IF (sym.eq.1) dimr=dimr+1
      ALLOCATE(tmat(dimr,dimr))
      tmat=0.d0

      DO i=1,dimc
         IF (sym.eq.0 .or. sym.eq.1) THEN
            dl=2*(i-1)+sym
         ELSE
            dl=i-1
         ENDIF
         DO j=1,dimr-dl
            tmat(j,j+dl)=dmat(j,i)
         ENDDO
      ENDDO

!     Reflect to lower half if needed
      IF (refl) THEN
         DO i=2,dimr
            DO j=1,i-1
               tmat(i,j)=tmat(j,i)
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE MatDiag2Utriang

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE MatUtriang2Diag(dmat,tmat,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts matrix from upper triangular to stored diagonal form

      implicit none
      real*8, allocatable, intent(out) :: dmat(:,:)
      real*8, intent(in)     :: tmat(:,:)
      integer, intent(inout) :: sym
      integer :: i,j,dimr,dl

      dimr=SIZE(tmat,1)
      ALLOCATE(dmat(dimr,dimr))
      dmat=0.d0

      DO i=1,dimr
         dl=i-1
         DO j=1,dimr-dl
            dmat(j,i)=tmat(j,j+dl)
         ENDDO
      ENDDO

!     Trim the matrix based on symmetry and/or zero columns. Start with
!     nosym; let the following call determine the proper designation
      sym=2
      call TrimMatDiagForm(dmat,sym)

      END SUBROUTINE MatUtriang2Diag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE TrimMatDiagForm(dmat,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Trims matrix in stored-diagonal form by removing zero columns
! (i.e. those with all elements with absolute values smaller than 'tol')
! If the input symmetry is '2' (nosym), a check is made to see if
! symmetry can be imposed, reducing the size further

      implicit none
      real*8, allocatable, intent(inout) :: dmat(:,:)
      real*8, allocatable  :: tmpmat(:,:)
      integer, intent(inout) :: sym
      real*8  :: tol,elem
      integer :: i,j,nr,nc,finalcol,newr,newc,newsym

      tol=1.0d-14

      nr=SIZE(dmat,1)
      nc=SIZE(dmat,2)
      finalcol=nc

      DO i=nc,1,-1
!        Check column to see if any element exceeds tol
         elem=0.d0
         DO j=1,nr
            elem=abs(dmat(j,i))
            IF (elem.gt.tol) EXIT
         ENDDO
!        If no element in the column is smaller than tol, trim it;
!        otherwise, exit
         IF (elem.gt.tol) EXIT
         finalcol=finalcol-1
      ENDDO

!     Have symmetry already: just trim zero columns
      IF (sym.eq.0 .or. sym.eq.1) THEN
         IF (finalcol.eq.nc) THEN
            RETURN
         ELSE
            ALLOCATE(tmpmat(nr,finalcol))
            tmpmat(:,1:finalcol)=dmat(:,1:finalcol)
         ENDIF
!     Check to see if symmetry can be imposed
      ELSE
         elem=0.d0
         IF (mod(finalcol,2).eq.1) THEN  ! even or no sym
            newsym=0
            DO i=2,finalcol,2
               DO j=1,nr
                  elem=abs(dmat(j,i))
                  IF (elem.gt.tol) EXIT
               ENDDO
               IF (elem.gt.tol) EXIT
            ENDDO
            IF (elem.gt.tol) newsym=2
         ELSE  ! odd or no sym
            newsym=1
            DO i=1,finalcol,2  
               DO j=1,nr
                  elem=abs(dmat(j,i))
                  IF (elem.gt.tol) EXIT
               ENDDO
               IF (elem.gt.tol) EXIT
            ENDDO
            IF (elem.gt.tol) newsym=2
         ENDIF

!        If no columns are to be trimmed, simply return...
         IF (finalcol.eq.nc .and. newsym.eq.2) RETURN

         IF (newsym.eq.0 .or. newsym.eq.1) THEN
            newr=nr-newsym
            newc=(finalcol+1-newsym)/2
            ALLOCATE(tmpmat(newr,newc))
            DO i=1,newc
               j=2*i-1+newsym
               tmpmat(1:newr,i)=dmat(1:newr,j)
            ENDDO
         ELSE
            ALLOCATE(tmpmat(nr,finalcol))
            tmpmat(:,1:finalcol)=dmat(:,1:finalcol)
         ENDIF
         sym=newsym
      ENDIF

      DEALLOCATE(dmat)
      ALLOCATE(dmat(SIZE(tmpmat,1),SIZE(tmpmat,2)))
      dmat=tmpmat
      DEALLOCATE(tmpmat)

      END SUBROUTINE TrimMatDiagForm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Vec2Matrix(v,M,nr,nc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wraps a vector into matrix format. nr,nc = number of rows,columns. The
! row index changes fastest

      implicit none
      real*8, intent(in)  :: v(:)
      real*8, allocatable, intent(out) :: M(:,:)
      integer, intent(in) :: nr,nc
      integer :: i,ir,ic

      IF (SIZE(v).ne.nr*nc) THEN
         write(*,*) 'The dimensions of v do not match those given for M'
         call AbortWithError("Error in Vec2Matrix()")
      ENDIF

      ALLOCATE(M(nr,nc))

      DO i=1,SIZE(v)
         ir=mod(i-1,nr)+1
         ic=(i-1)/nr+1
         M(ir,ic)=v(i)
      ENDDO

      END SUBROUTINE Vec2Matrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Mat2Vec(v,M,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wrapper for matrix-to-vector unwrapping. Setting sym=.TRUE. stores the
! upper triangular portion with off-diagonal elements multiplied by
! sqrt(2)

      implicit none
      real*8, allocatable, intent(out):: v(:)
      real*8, intent(in)  :: M(:,:)
      logical, intent(in) :: sym

      IF (sym) THEN
         call SymMat2Vec(v,M)
      ELSE
         call Matrix2Vec(v,M)
      ENDIF

      END SUBROUTINE Mat2Vec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Matrix2Vec(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Unwraps a matrix to vector. The vector is filled column-at-a-time

      implicit none
      real*8, allocatable, intent(out):: v(:)
      real*8, intent(in) :: M(:,:)
      integer :: nr,nc,i,j,k

      nr=SIZE(M,1)
      nc=SIZE(M,2)

      ALLOCATE(v(nr*nc))
      k=1
      DO j=1,nc
         DO i=1,nr
            v(k)=M(i,j)
            k=k+1
         ENDDO
      ENDDO

      END SUBROUTINE Matrix2Vec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Vec2Dmat(v,M,nr,nc,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wraps a vector to matrix format, then converts to stored diagonal form

      implicit none
      real*8, intent(in)  :: v(:)
      real*8, allocatable, intent(out) :: M(:,:)
      real*8, allocatable :: D(:,:)
      integer, intent(in)  :: nr,nc
      integer, intent(out) :: sym

      call Vec2Mat(v,M,nr,nc) 
      call MatUtriang2Diag(D,M,sym)
      DEALLOCATE(M)
      ALLOCATE(M(SIZE(D,1),SIZE(D,2)))
      M=D
      DEALLOCATE(D)

      END SUBROUTINE Vec2Dmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE SymMat2Vec(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Unwraps a symmetric matrix to vector, keeping the upper triangle. 
! Vector v contains the diagonal elements first, one-off-diagonal, 
! two-off-diagonal, etc.

      implicit none
      real*8, allocatable, intent(out) :: v(:)
      real*8, intent(in) :: M(:,:)
      integer :: n,i,j,k,l
      real*8  :: s

      s=sqrt(2.d0)

      n=SIZE(M,1)

      IF (SIZE(M,2).ne.n) &
         call AbortWithError('SymMat2Vec(): M is not symmetric')
     
      ALLOCATE(v(n*(n+1)/2)) 

      k=1
      DO j=n,1,-1
         l=n-j
         DO i=1,j
            v(k)=M(i,i+l)
            k=k+1
         ENDDO
      ENDDO

!     Multiply by sqrt(2) so that reduc() will preserve Frobenius norm
      v(n+1:)=s*v(n+1:)

      END SUBROUTINE SymMat2Vec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Vec2SymMat(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts a vector, in the form generated by SymMat2Vec(), to a square 
! matrix

      implicit none
      real*8, intent(in) :: v(:)
      real*8, allocatable, intent(out) :: M(:,:)
      integer :: n,i,j,k,l
      real*8  :: s

      s=1/sqrt(2.d0)
      n=GetSymN(SIZE(v))
      ALLOCATE(M(n,n))

!     Fill the upper triangle
      k=1
      DO j=n,1,-1
         l=n-j
         DO i=1,j
            M(i,i+l)=v(k)
            k=k+1
         ENDDO
      ENDDO

!     Copy upper triangle to lower triangle, rescaling by 1/sqrt(2)
      DO i=2,n
         DO j=1,i-1
            M(j,i)=s*M(j,i)
            M(i,j)=M(j,i)
         ENDDO
      ENDDO

      END SUBROUTINE Vec2SymMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE SymPackMat2Vec(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Store a symmetric matrix as a vector in symmetric-packed form. The 
! entries are the lower triangle ordered by columns

      implicit none
      real*8, allocatable, intent(out) :: v(:)
      real*8, intent(in) :: M(:,:)
      integer :: n,i,j,k

      n=SIZE(M,1)

      IF (SIZE(M,2).ne.n) &
         call AbortWithError('SymPackMat2Vec(): M is not symmetric')

      ALLOCATE(v(n*(n+1)/2))

!     Lower triangle version
!      k=1
!      DO j=1,n
!         DO i=j,n
!            v(k)=M(i,j)
!            k=k+1
!         ENDDO
!      ENDDO

!     Upper triangle version
      k=1
      DO j=1,n
         DO i=1,j
            v(k)=M(i,j)
            k=k+1
         ENDDO
      ENDDO

      END SUBROUTINE SymPackMat2Vec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Vec2SymPackMat(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Unpacks the vector version of a symmetric packed matrix to matrix form

      implicit none
      real*8, intent(in) :: v(:)
      real*8, allocatable, intent(out) :: M(:,:)
      integer :: n,i,j,k

      n=GetSymN(SIZE(v))
      ALLOCATE(M(n,n))

!     Fill the lower triangle
!      k=1
!      DO j=1,n
!         DO i=j,n
!            M(i,j)=v(k)
!            k=k+1
!         ENDDO
!      ENDDO

!     Copy lower triangle to upper triangle
!      DO i=1,n-1
!         DO j=i+1,n
!            M(i,j)=M(j,i)
!         ENDDO
!      ENDDO

!!!   ALT
!     Fill the upper triangle
      k=1
      DO j=1,n
         DO i=1,j
            M(i,j)=v(k)
            k=k+1
         ENDDO
      ENDDO

!     Copy upper triangle to lower triangle
      DO i=2,n
         DO j=1,i-1
            M(i,j)=M(j,i)
         ENDDO
      ENDDO

      END SUBROUTINE Vec2SymPackMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE SymmetrizeMat(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Symmetrizes a matrix by averaging the i,j-th element with the j,i-th

      implicit none
      real*8, intent(inout) :: M(:,:)
      integer :: n,i,j

      n=SIZE(M,1)
      IF (SIZE(M,2).ne.n) &
         call AbortWithError('SymmetrizeMat(): M is not square')

!     Average the off-diagonal elements
      DO i=2,n
         DO j=1,i-1
            M(j,i)=0.5d0*(M(j,i)+M(i,j))
            M(i,j)=M(j,i)
         ENDDO
      ENDDO

      END SUBROUTINE SymmetrizeMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ExtractRowsColsMatrix(M,row,col)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts rows and cols of matrix, removing those for which row(i) and
! col(j) are less than one

      implicit none
      real*8, allocatable, intent(inout) :: M(:,:)
      integer, allocatable, intent(inout) :: row(:),col(:)
      real*8, allocatable  :: Mt(:,:)
      integer, allocatable :: rowt(:),colt(:)
      integer :: i,j,k,l,r,c,nzr,nzc

      r=SIZE(row)
      c=size(col)

!     Error checking
      IF (SIZE(M,1).ne.r) THEN
         write(*,*) 'M(r,c) has ',SIZE(M,1),' rows, but must have ',&
                    r,' rows'
         call AbortWithError('ExtractRowsColsMatrix(): bad dimensions')
      ENDIF
      IF (SIZE(M,2).ne.c) THEN
         write(*,*) 'M(r,c) has ',SIZE(M,2),' cols, but must have ',&
                    c,' rows'
         call AbortWithError('ExtractRowsColsMatrix(): bad dimensions')
      ENDIF

!     Count number of rows and cols to keep
      nzr=0
      DO i=1,r
         IF (row(i).gt.0) nzr=nzr+1
      ENDDO
      nzc=0
      DO j=1,c
         IF (col(j).gt.0) nzc=nzc+1
      ENDDO

!     Deallocate and exit if row or col has no entries to keep
      IF (nzr.eq.0 .or. nzc.eq.0) THEN
         DEALLOCATE(M)
         IF (nzr.eq.0) DEALLOCATE(row)
         IF (nzc.eq.0) DEALLOCATE(col)
         RETURN
      ENDIF

!     Copy greater-than-zero portions of row and col to tmp array
      ALLOCATE(Mt(nzr,nzc),rowt(nzr),colt(nzc))
      k=1
      DO i=1,r
         IF (row(i).gt.0) THEN
            rowt(k)=row(i)
            k=k+1
         ENDIF
      ENDDO
      l=1
      DO j=1,c
         IF (col(j).gt.0) THEN
            colt(l)=col(j)
            l=l+1
         ENDIF
      ENDDO

!     Extract the rows and cols for which BOTH row(i) and col(j) are
!     greater than zero
      l=1
      DO j=1,c
         IF (col(j).gt.0) THEN
            k=1
            DO i=1,r
               IF (row(i).gt.0) THEN
                  Mt(k,l)=M(i,j)
                  k=k+1
               ENDIF
            ENDDO
            l=l+1
         ENDIF
      ENDDO

      DEALLOCATE(M,row,col)
      ALLOCATE(M(nzr,nzc),row(nzr),col(nzc))
      M(:,:)=Mt(:,:)
      row(:)=rowt(:)
      col(:)=colt(:)

      END SUBROUTINE ExtractRowsColsMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE GetIdentityMatrix(M,n,norm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      real*8, allocatable, intent(out) :: M(:,:)
      integer, intent(in) :: n
      logical, intent(in) :: norm
      integer :: i
      real*8  :: val

      val=1.d0
      IF (norm) val=val/sqrt(REAL(n))

      ALLOCATE(M(n,n))
      M=0.d0
      DO i=1,n
         M(i,i)=val
      ENDDO

      END SUBROUTINE GetIdentityMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Cascade(v,up)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Moves all entries in a vector up or down one position

      implicit none
      real*8, intent(inout) :: v(:)
      logical, intent(in) :: up
      integer :: i,n

      n=SIZE(v)

      IF (up) THEN
         DO i=n-1,1,-1
            v(i+1)=v(i)
         ENDDO
         v(1)=0.d0
      ELSE
         DO i=2,n
            v(i-1)=v(i)
         ENDDO
         v(n)=0.d0
      ENDIF

      END SUBROUTINE Cascade

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION ibisect(v,i) RESULT(jl)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds index of last item of v with the same index as i, by bisection
! Data values must be in monotonically-ascending order

      implicit none
      integer, intent(in) :: v(:)
      integer, intent(in) :: i
      integer :: n,ju,jl,t

      n=SIZE(v)

      CALL ERROR(i<1 .or. i>n, 'ibisect(): i out of range')

      jl=i
      ju=n+1
      DO
        IF (ju-jl.le.1) EXIT
        t=(ju+jl)/2
        IF (v(t).eq.v(i)) THEN
           jl=t
        ELSE
           ju=t
        ENDIF
      ENDDO

      END FUNCTION ibisect

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ifind(v,i,jset)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds index of last item of v with the same index as i, by bisection
! Data values must be in monotonically-ascending order

      implicit none
      integer, intent(in) :: v(:)
      integer, intent(in) :: i
      integer, intent(out) :: jset(2)
      integer :: n,jf,ju,jl,t

      n=SIZE(v)

!     Early exit if i is outside the upper/lower bounds of the list
      IF (i.lt.v(1) .or. i.gt.v(n)) THEN
          jset=(/0,0/)
          RETURN
      ENDIF

!     Bisect to find the last value in the list equal to i
      jf=0
      jl=1
      ju=n+1
      DO
         IF (ju-jl.le.1) EXIT
         t=(ju+jl)/2
         IF (v(t).le.i) THEN
            jl=t
            IF (v(t).lt.i) jf=t
         ELSE
            ju=t
         ENDIF
      ENDDO

!     Exit if i is not found
      IF (v(jl).ne.i) THEN 
         jset=(/0,0/)
         RETURN
      ENDIF

      jset(2)=jl
!     Bisect to find the first value in the list equal to i
!     jf,jl from the previous bisection are used
      DO
         IF (jl-jf.le.1) EXIT
         t=(jl+jf)/2
         IF (v(t).lt.i) THEN
            jf=t
         ELSE
            jl=t
         ENDIF
      ENDDO
      jset(1)=jl

      END SUBROUTINE ifind

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER FUNCTION FACRL(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes n factorial

      implicit none
      integer, intent(in) :: n
      integer             :: j, prod

      CALL ERROR(n<0,'FACRL(nn) for nn < 0')

      prod=1
      IF (n.GE.2) THEN
         DO j=n,1,-1
            prod=prod*j
         END DO
      ENDIF

      FACRL=prod

      END FUNCTION FACRL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      LOGICAL FUNCTION IntListsRSame(l1,l2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares lists of integers to see if they contain the same values

      implicit none
      integer, intent(in)  :: l1(:),l2(:)
      integer :: len1,len2,i

      len1=SIZE(l1)
      len2=SIZE(l2)

      IF (len1.ne.len2) THEN
         IntListsRSame=.FALSE.
      ELSE
         IntListsRSame=.TRUE.
         DO i=1,len1
            IF(.NOT.ANY(l2.eq.l1(i))) IntListsRSame=.FALSE.
         ENDDO
      ENDIF

      END FUNCTION IntListsRSame

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      LOGICAL FUNCTION checkfornans(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Search a matrix for NaN values, return .TRUE. if there are any

      implicit none
      real*8, intent(in) :: M(:,:)
      integer :: nr,nc,i,j
      
      nr=SIZE(M,1)
      nc=SIZE(M,2)
      
      checkfornans=.FALSE.
      DO i=1,nr
         DO j=1,nc
            IF (M(i,j).ne.M(i,j)) THEN
               checkfornans=.TRUE.
               RETURN
            ENDIF
         ENDDO
      ENDDO
      
      END FUNCTION checkfornans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER FUNCTION GetRandomInteger(n)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Picks a random integer between 1 and n

      implicit none
      integer, intent(in) :: n
      real*8 :: x(1)

      CALL ERROR(n<1,'GetRandomIndex(): n < 1')

      call random_number(x)

      GetRandomInteger=nint(n*x(1)+0.5)

      END FUNCTION GetRandomInteger

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER FUNCTION GetSymN(s)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds n, the dimension for a symmetric matrix wrapped to vector with
! size s

      implicit none
      integer, intent(in) :: s
      integer :: i,n

!     Make sure s is consistent with a square matrix and compute n
      i=1
      n=1
      DO
        IF (i.ge.s) EXIT
        n=n+1
        i=i+n
      ENDDO
      IF (i.ne.s) &
         call AbortWithError('GetSymN(): SIZE(v) -X-> n x n matrix')

      GetSymN=n

      END FUNCTION GetSymN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION compareStringsForEQ( s1, s2 ) RESULT( equal_result )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares strings and does not differentiate between lower and upper
! case. Also, leading and/or trailing spaces are ignored.

      implicit none
      character(LEN=*), intent(in) :: s1, s2
      logical                      :: equal_result
      integer(4)                   :: is, offset1, offset2, L, distance
      integer(4)                   :: length1, length2
      character(LEN=1)             :: char1, char2
      logical                      :: equal

      equal = .TRUE.

!     Determine the length of each string, NOT counting empty spaces.
      length1 = LEN_TRIM(ADJUSTL(s1))
      length2 = LEN_TRIM(ADJUSTL(s2))

      IF ( length1 /= length2 ) then
         equal = .FALSE.
      ELSE

!        Left adjust each string.
         offset1 = VERIFY( s1, " " ) - 1
         offset2 = VERIFY( s2, " " ) - 1

!        Lexographic distance between 'a' and 'A' in ASCII list.
         L = ABS( ICHAR('a') - ICHAR('A') )

         do is = 1 , LEN_TRIM(s1)
!           Next characters in strings.
            char1 = s1(offset1+is:offset1+is)
            char2 = s2(offset2+is:offset2+is)

            distance = ABS( ICHAR(char1) - ICHAR(char2) )
            if ( distance /= 0 .AND. distance /= L ) then
               equal = .FALSE.
               EXIT
            endif
         enddo

      ENDIF

      equal_result = equal

      END FUNCTION compareStringsForEQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE UTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
