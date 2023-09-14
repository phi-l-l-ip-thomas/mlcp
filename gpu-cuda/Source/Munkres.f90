!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MUNKRES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains the Munkres algorithm for solving the assignment
! problem

      USE ERRORTRAP
      USE UTILS

      implicit none
      real*8, private  :: munkres_time=0.d0
      logical, private :: MUNKRES_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeMunkres()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      munkres_time = 0.d0
      MUNKRES_SETUP = .TRUE.

      end subroutine InitializeMunkres

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeMunkres()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MUNKRES_SETUP) call InitializeMunkres()

      MUNKRES_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total munkres assignment time     (s)',&
                            munkres_time

      end subroutine DisposeMunkres

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function AssignMatrix(M) result (A)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates assignment matrix via Hungarian (Munkres) algorithm. M, the
! cost matrix, must correspond to maximum matching case (all elements
! non-negative, largest value matched.

      implicit none
      real*8, intent(in)  :: M(:,:)
      real*8, allocatable :: A(:,:)
      real*8, parameter :: tol=1.d-15
      real*8 :: tmp
      logical, allocatable :: covrow(:),covcol(:)
      integer, allocatable :: Z(:,:)
      integer :: i,j,nr,nc,prow,pcol
      real*8  :: t1,t2

      IF (.NOT. MUNKRES_SETUP) call InitializeMunkres()
      call CPU_TIME(t1)

      nr=SIZE(M,1)
      nc=SIZE(M,2)

      IF (nc.lt.nr) THEN
         write(*,*) '#cols (',nc,') must equal or exceed #rows (',nr,')'
         call AbortWithError('AssignMatrix(): nc < nr')
      ENDIF
      IF (ANY(M(:,:).lt.0.d0)) THEN
         DO j=1,nc
            DO i=1,nr
               IF (M(i,j).lt.0.d0) write(*,*) i,j,M(i,j)
            ENDDO
         ENDDO
         call AbortWithError('AssignMatrix(): M must be non-negative')
      ENDIF

      ALLOCATE(A(nr,nc),Z(nr,nc))
      ALLOCATE(covrow(nr),covcol(nc))

!     Process A so that largest value corresponds to minimum cost
      A(:,:)=-M(:,:)

!     (Step 1) subtract minimum row element
      DO i=1,nr
         tmp=MINVAL(A(i,:))
         A(i,:)=A(i,:)-tmp
      ENDDO

!     (Step 2) locate zeros in A. Set Z(i,j)=1 for zero found in cross,
!     and Z(i,j)=-1 for zero in cross of previously found zero
      Z(:,:)=0
      covrow(:)=.FALSE.
      covcol(:)=.FALSE.
      DO i=1,nr
         DO j=1,nc
            IF (abs(A(i,j)).lt.tol) THEN ! Zero found
               IF ((.not.covrow(i)) .and. (.not.covcol(j))) THEN
                  Z(i,j)=1  ! "starred" zero
                  covrow(i)=.TRUE.
                  covcol(j)=.TRUE.
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO
!        (Step 3) count number of cols that are assigned
         covrow(:)=.FALSE.
         covcol(:)=.FALSE.
         DO j=1,nc
            IF (ANY(Z(:,j).eq.1)) covcol(j)=.TRUE.
         ENDDO

!        Exit if all rows or cols are assigned
         IF (COUNT(covcol).ge.min(nr,nc)) EXIT

!        (Steps 4,6): check solutions and rebalance cost
         call Munkres_steps4and6(A,Z,covrow,covcol,prow,pcol)

!        (Step 5): alternate assignments to find minimum cost
         call Munkres_step5(Z,prow,pcol)
      ENDDO

!     (Step 7): assignment matrix -> permutation matrix
      A(:,:)=0.d0
      DO j=1,nc
         DO i=1,nr
            IF (Z(i,j).eq.1) A(i,j)=1.d0
         ENDDO
      ENDDO

      DEALLOCATE(covrow,covcol)
      DEALLOCATE(Z)

      call CPU_TIME(t2)
      munkres_time=munkres_time+t2-t1

      end function AssignMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetMunkresAssignVec(A) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Convenience function for extracting vector with col indices 
! corresponding to assignments from Munkres assignment matrix

      implicit none
      real*8, intent(in) :: A(:,:)
      real, parameter :: tol=1.d-15
      integer, allocatable :: w(:)
      integer :: i,j,r,c

      r=SIZE(A,1)
      c=SIZE(A,2)

      ALLOCATE(w(r))
      w(:)=0

      DO i=1,r
         DO j=1,c
            IF (abs(A(i,j)-1.d0) .lt. tol) THEN
               w(i)=j
               EXIT
            ENDIF
         ENDDO
      ENDDO

      end function GetMunkresAssignVec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Munkres_steps4and6(A,Z,covrow,covcol,pathrow,pathcol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Steps 4,6 in Hungarian algorithm: find a non-covered zero (if present)
! prime it, and then look for starred zero to add to linear program path

      implicit none
      real*8, intent(inout)  :: A(:,:)
      integer, intent(inout) :: Z(:,:)
      logical, intent(inout) :: covrow(:),covcol(:)
      integer, intent(out) :: pathrow,pathcol
      logical :: doneo,donei
      integer :: r,c,thecol

! Check for uncovered zeros (step 4) and resulting solutions with zero
! covered (via step 5 elsewhere). If none, then rebalance cost (step 6).

      doneo=.FALSE.
      DO
         IF (doneo) EXIT ! Loop over composite of steps 4,6
         donei=.FALSE.

!        (Step 4): look for uncovered zeros
         DO
            IF (donei) EXIT ! Loop over step 4

!           Find a zero in a column that is not covered.
            call findazero(A,covrow,covcol,r,c)

            IF (r.eq.0) THEN
               donei=.TRUE.
            ELSE
               Z(r,c)=-1 ! Mark as "primed" zero

!              Get col of "starred" zero in this row, if present
               thecol=get_first_ind(Z(r,:),1)             
               IF (thecol.ne.0) THEN
                  covrow(r)=.TRUE.
                  covcol(thecol)=.FALSE.
               ELSE
                  donei=.TRUE.
                  doneo=.TRUE.
                  pathrow=r
                  pathcol=c
               ENDIF
            ENDIF         
         ENDDO

!        (Step 6): balance cost based on minimum uncovered value
         call recostify(A,covrow,covcol)
      ENDDO

      end subroutine Munkres_steps4and6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine findazero(A,covrow,covcol,row,col)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Find row,col of first zero encountered, skipping covered rows,cols

      implicit none
      real*8, intent(in)  :: A(:,:)
      logical, intent(in) :: covrow(:),covcol(:)
      integer, intent(out) :: row,col
      real*8, parameter :: tol=1.d-15
      integer :: r,c,nr,nc

      nr=SIZE(A,1)
      nc=SIZE(A,2)
      row=0
      col=0

      DO r=1,nr
         IF (covrow(r)) CYCLE
         DO c=1,nc
            IF (covcol(c)) CYCLE
            IF (abs(A(r,c)).lt.tol) THEN
               row=r
               col=c
               EXIT
            ENDIF
         ENDDO
         IF (row.gt.0) EXIT
      ENDDO

      end subroutine findazero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine recostify(A,covrow,covcol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Modifies A by adding/subtracting min uncovered value to rebalance cost

      implicit none
      real*8, intent(inout) :: A(:,:)
      logical, intent(in)   :: covrow(:),covcol(:)
      integer :: r,c,nr,nc
      real*8  :: minA

      nr=SIZE(A,1)
      nc=SIZE(A,2)
      minA=1.d99

!     Find minimum uncovered value
      DO c=1,nc
         IF (covcol(c)) CYCLE
         DO r=1,nr
            IF (covrow(r)) CYCLE
            IF (A(r,c).lt.minA) minA=A(r,c)
         ENDDO
      ENDDO

      IF (minA.eq.0.d0) RETURN

!     Add minimum uncovered value to covered rows
!     Subtract minimum uncovered value from uncovered cols
      DO c=1,nc
         DO r=1,nr
            IF (covrow(r)) A(r,c)=A(r,c)+minA
            IF (.not.covcol(c)) A(r,c)=A(r,c)-minA
         ENDDO
      ENDDO

      end subroutine recostify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Munkres_step5(Z,pathrow,pathcol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Steps 4,6 in Hungarian algorithm: find a non-covered zero (if present)
! prime it, and then look for starred zero to add to linear program path

      implicit none
      integer, intent(inout) :: Z(:,:)
      integer, intent(in)  :: pathrow,pathcol
      integer, allocatable :: path(:,:)
      logical :: done
      integer :: i,j,m,n,thez

      m=MIN(SIZE(Z,1),SIZE(Z,2))
      n=1

      ALLOCATE(path(m,2))
      path(n,:)=(/pathrow,pathcol/)

      DO
!        Get row of "starred" zero in column of most recent path item
         thez=get_first_ind(Z(:,path(n,2)),1)
         IF (thez.eq.0) EXIT
         n=n+1
         call append_munkres_path(path,m,n,thez,path(n-1,2))

!        Get col of "primed" zero in row of most recent path item
         thez=get_first_ind(Z(path(n,1),:),-1)
         n=n+1
         call append_munkres_path(path,m,n,path(n-1,1),thez)
      ENDDO

!     Toggle the "starred" zeros of the path in Z
      DO i=1,n
         IF (Z(path(i,1),path(i,2)).eq.1) THEN
            Z(path(i,1),path(i,2))=0
         ELSE
            Z(path(i,1),path(i,2))=1
         ENDIF
      ENDDO

!     Flush all primes in Z
      DO j=1,SIZE(Z,2)
         DO i=1,SIZE(Z,1)
            IF (Z(i,j).eq.-1) Z(i,j)=0
         ENDDO
      ENDDO

      DEALLOCATE(path)

      end subroutine Munkres_step5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine append_munkres_path(path,expandby,newn,r,c)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Add a new element pair to array 'path'. If it exceeds the array bounds
! then expand array by size 'expandby'

      implicit none
      integer, allocatable, intent(inout) :: path(:,:)
      integer, intent(in)  :: expandby,newn,r,c
      integer, allocatable :: ptmp(:,:)
      integer :: n

      n=SIZE(path,1)

!     Replace path with larger array if path is too small
      IF (newn.gt.n) THEN
         ALLOCATE(ptmp(n,2))
         ptmp(1:n,1:2)=path(1:n,1:2)
         DEALLOCATE(path)
         ALLOCATE(path(max(newn,n+expandby),2))
         path(1:n,1:2)=ptmp(1:n,1:2)
         DEALLOCATE(ptmp)
      ENDIF

      path(newn,:)=(/r,c/)

      end subroutine append_munkres_path

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function get_first_ind(v,val) result(ind)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds first index of item in v matching v. Returns 0 if not found.

      implicit none
      integer, intent(in) :: v(:)
      integer :: i,val,ind

      ind=0
      DO i=1,SIZE(v)
         IF (v(i).eq.val) THEN
            ind=i
            EXIT
         ENDIF
      ENDDO

      end function get_first_ind

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MUNKRES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
