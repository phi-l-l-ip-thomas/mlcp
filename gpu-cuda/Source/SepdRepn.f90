!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE SEPDREPN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains structures needed for CP representation

      USE ERRORTRAP
      USE UTILS

      TYPE CPvec
         INTEGER, ALLOCATABLE :: nbas(:)
         REAL*8,  ALLOCATABLE :: coef(:),base(:,:)
      END TYPE CPvec

      TYPE CPptr
         INTEGER, ALLOCATABLE :: nbas(:)
         REAL*8,  POINTER     :: coef(:),base(:,:)
      END TYPE CPptr

      TYPE CPblock
         INTEGER, ALLOCATABLE :: nbas(:)
         REAL*8, ALLOCATABLE  :: coef(:,:),base(:,:)
      END TYPE CPblock

      INTERFACE CHECKNBAS
         MODULE PROCEDURE CHECKNBASvv
         MODULE PROCEDURE CHECKNBASvQ
         MODULE PROCEDURE CHECKNBASQQ
      END INTERFACE

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NewCPvec(v,nbas,rnk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Initializes a new vector in CP-format

      IMPLICIT NONE
      TYPE (CPvec)        :: v
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rnk
      INTEGER :: nrdim,ndof,i

      IF (rnk.lt.1) THEN
         write(*,*) 'Error: rnk must be at least 1!'
         call AbortWithError('Error in NewCPvec()')
      ENDIF

!     Determine how the base is arranged from 'nbas', which holds the
!     number of basis functions for each DOF
      ndof=0
      nrdim=0
      DO i=1,SIZE(nbas)
         IF (nbas(i).lt.1) EXIT
         nrdim=nrdim+nbas(i)
         ndof=ndof+1
      ENDDO

      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPvec()')
      ENDIF

      ALLOCATE(v%base(nrdim,rnk),v%coef(rnk),v%nbas(ndof))
      v%nbas(1:ndof)=nbas(1:ndof)

      END SUBROUTINE NewCPvec


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NewCPptr(v,nbas,rnk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Initializes a new pointer in CP-format

      IMPLICIT NONE
      TYPE (CPptr)        :: v
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rnk
      INTEGER :: nrdim,ndof,i

      IF (rnk.lt.1) THEN
         write(*,*) 'Error: rnk must be at least 1!'
         call AbortWithError('Error in NewCPptr()')
      ENDIF

!     Determine how the base is arranged from 'nbas', which holds the
!     number of basis functions for each DOF
      ndof=0
      nrdim=0
      DO i=1,SIZE(nbas)
         IF (nbas(i).lt.1) EXIT
         nrdim=nrdim+nbas(i)
         ndof=ndof+1
      ENDDO

      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPptr()')
      ENDIF

!     Same as NewCPvec() except that coef/base arrays are not allocated
      ALLOCATE(v%nbas(ndof))
      v%nbas(1:ndof)=nbas(1:ndof)

      END SUBROUTINE NewCPptr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NewCPmat(v,nbas,rnk,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Initializes a CP-format outer product of square matrices

      IMPLICIT NONE
      TYPE (CPvec)        :: v
      LOGICAL, INTENT(IN) :: sym
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rnk
      INTEGER :: nrdim,ndof,i

      IF (rnk.lt.1) THEN
         write(*,*) 'Error: rnk must be at least 1!'
         call AbortWithError('Error in NewCPmat()')
      ENDIF

!     The structure of the base has strides of nbas(i)^2 elements (or
!     nbas(i)*(nbas(i)+1)/2 elements if symmetry is present)
      ndof=0
      nrdim=0
      DO i=1,SIZE(nbas)
         IF (nbas(i).lt.1) EXIT
         IF (sym) THEN
            nrdim=nrdim+nbas(i)*(nbas(i)+1)/2
         ELSE
            nrdim=nrdim+nbas(i)*nbas(i)
         ENDIF
         ndof=ndof+1
      ENDDO

      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPmat()')
      ENDIF

      ALLOCATE(v%base(nrdim,rnk),v%coef(rnk),v%nbas(ndof))

      IF (sym) THEN
         DO i=1,ndof
            v%nbas(i)=nbas(i)*(nbas(i)+1)/2
         ENDDO
      ELSE
         DO i=1,ndof
            v%nbas(i)=nbas(i)*nbas(i)
         ENDDO
      ENDIF

      END SUBROUTINE NewCPmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetIdentityCPMatrix(v,nbas,sym)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets rank-1 CP outer-product-of-identity-matrices. If only the upper
! triangle of Hsym is stored, set sym to .TRUE.

      implicit none
      TYPE (CPvec), INTENT(OUT) :: v
      integer, intent(in)  :: nbas(:)
      logical, intent(in)  :: sym
      real*8, allocatable  :: tmat(:,:),tvec(:)
      integer :: j,ndof,gi,gf

      ndof=SIZE(nbas)

      call NewCPmat(v,nbas,1,sym)
      v%coef=1.d0

!     Get an identity matrix for each DOF
      gi=1
      DO j=1,ndof
         gf=gi+v%nbas(j)-1
         call GetIdentityMatrix(tmat,nbas(j),.FALSE.)
         call Mat2Vec(tvec,tmat,sym)
         v%base(gi:gf,1)=tvec
         DEALLOCATE(tmat,tvec)
         gi=gi+v%nbas(j)
      ENDDO

      end subroutine GetIdentityCPMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPvec(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-vector in neat format

      implicit none
      TYPE (CPvec), INTENT(IN) :: v
      integer :: r,rk,j,d,i,n,inr
      character*64 :: frmt

      rk=SIZE(v%coef)
      d=SIZE(v%nbas)

      write(frmt,'(A,I0,A)') '(A,2X,',rk,'(X,ES13.6,X))'
      write(*,frmt) ' Vcoef =',(v%coef(r),r=1,rk)
      write(frmt,'(A,I0,A)') '(2(I4,X),',rk,'(X,f13.10,X))'
      inr=1
      do j=1,d
         n=v%nbas(j)
         do i=1,n
            write(*,frmt) j,i,(v%base(inr,r),r=1,rk)
            inr=inr+1
         enddo
      enddo
      write(*,*)

      end subroutine PrintCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPptr(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints pointer to CP-vector in neat format

      implicit none
      TYPE (CPptr), INTENT(IN) :: v
      integer :: r,rk,j,d,i,n,inr
      character*64 :: frmt

      rk=SIZE(v%coef)
      d=SIZE(v%nbas)

      write(frmt,'(A,I0,A)') '(A,2X,',rk,'(X,ES13.6,X))'
      write(*,frmt) ' Vcoef =',(v%coef(r),r=1,rk)
      write(frmt,'(A,I0,A)') '(2(I4,X),',rk,'(X,f13.10,X))'
      inr=1
      do j=1,d
         n=v%nbas(j)
         do i=1,n
            write(*,frmt) j,i,(v%base(inr,r),r=1,rk)
            inr=inr+1
         enddo
      enddo
      write(*,*)

      end subroutine PrintCPptr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPmat(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format

      implicit none
      TYPE (CPvec), INTENT(IN) :: v
      integer :: r,rk,j,d,i,k,n,inr
      character*64 :: frmt

      rk=SIZE(v%coef)
      d=SIZE(v%nbas)

      DO r=1,rk
         write(*,'(/A,I0,A,f14.6)') 'RANK: ',r,'; coef = ',v%coef(r)
         inr=0
         DO j=1,d
            write(*,'(/A,I0)') 'dof : ',j
            n=nint(sqrt(real(v%nbas(j))))
            write(frmt,'(A,I0,A)') '(',n,'(X,f14.6))'
            DO i=1,n
               write(*,frmt) (v%base(inr+i+(k-1)*n,r),k=1,n)
            ENDDO
            inr=inr+v%nbas(j)
         ENDDO
      ENDDO
      write(*,*)

      end subroutine PrintCPmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE FlushCPvec(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Disposes CP-format vector

      IMPLICIT NONE
      TYPE (CPvec) :: v

      IF (ALLOCATED(v%base)) DEALLOCATE(v%base)
      IF (ALLOCATED(v%coef)) DEALLOCATE(v%coef)
      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)

      END SUBROUTINE FlushCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE FlushCPptr(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Disposes CP-format pointer

      IMPLICIT NONE
      TYPE (CPptr) :: v

      IF (ASSOCIATED(v%base)) NULLIFY(v%base)
      IF (ASSOCIATED(v%coef)) NULLIFY(v%coef)
      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)

      END SUBROUTINE FlushCPptr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE GetZeroCPvec(v,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Builds a zero CP-vec of rank 1

      IMPLICIT NONE
      TYPE (CPvec), INTENT(OUT) :: v
      INTEGER, INTENT(IN) :: nbas(:)

      call NewCPvec(v,nbas,1)
      v%base=0.d0
      v%coef=0.d0

      END SUBROUTINE GetZeroCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE GetRandomCPvec(v,nbas,rank)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Builds a CP-vec with random entries. No normalization is done.

      IMPLICIT NONE
      TYPE (CPvec) :: v
      INTEGER, ALLOCATABLE, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rank
      INTEGER :: i

      call NewCPvec(v,nbas,rank)
      DO i=1,rank
         call random_number(v%base(:,i))
         v%coef(i)=1.d0
      ENDDO
      v%base=v%base-0.5d0

      end SUBROUTINE GetRandomCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE Noisify(v,mag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Adds noise to base of CP-vector, noise in range [-mag,mag]

      IMPLICIT NONE
      TYPE (CPvec) :: v
      REAL*8, INTENT(IN)  :: mag
      REAL*8, ALLOCATABLE :: facs(:)
      INTEGER :: i

      ALLOCATE(facs(SIZE(v%base,1)))

      DO i=1,SIZE(v%coef)
         call random_number(facs(:))
         facs(:)=2*mag*(facs(:)-0.5d0)
         v%base(:,i)=v%base(:,i)+facs(:)
      ENDDO

      DEALLOCATE(facs)

      END SUBROUTINE Noisify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ReplaceVwithW(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Copies W into V, disposing W afterwards

      IMPLICIT NONE
      TYPE (CPvec) :: v,w

      call FlushCPvec(v)
      call CopyWtoV(v,w)
      call FlushCPvec(w)

      END SUBROUTINE ReplaceVwithW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CopyWtoV(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      IMPLICIT NONE
      TYPE (CPvec) :: v,w
      INTEGER :: nrk

      nrk=SIZE(w%coef)
      call NewCPvec(v,w%nbas,nrk)
      call GenCopyWtoV(v,w,1,nrk,1,nrk)

      END SUBROUTINE CopyWtoV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE GenCopyWtoV(v,w,vi,ve,wi,we)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V, 
! leaving W intact. v must be allocated beforehand.

      IMPLICIT NONE
      TYPE (CPvec) :: v,w
      INTEGER, INTENT(in) :: vi,ve,wi,we
      INTEGER :: rkv,rkw

      rkv=SIZE(v%coef)
      rkw=SIZE(w%coef)

      IF (vi.lt.1 .or. ve.gt.rkv .or. vi.gt.ve .or. &
          wi.lt.1 .or. we.gt.rkw .or. wi.gt.we .or. &
          we-wi.ne.ve-vi) THEN
          write(*,'(2A,6(X,I0))') 'Bad rank indices: ',&
          'vi,ve,rkv,wi,we,rkw =',vi,ve,rkv,wi,we,rkw
          CALL AbortWithError('Error in GenCopyWtoV()')
      ENDIF

      v%base(:,vi:ve)=w%base(:,wi:we)
      v%coef(vi:ve)=w%coef(wi:we)

      END SUBROUTINE GenCopyWtoV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ResizeV(v,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resizes vector v, either by truncating at a smaller rank or by adding
! space for extra terms.

      IMPLICIT NONE
      TYPE (CPvec), INTENT(INOUT) :: v
      TYPE (CPvec) :: w
      INTEGER, INTENT(IN) :: nrk
      INTEGER :: rkv

      rkv=SIZE(v%coef)

      IF (nrk.lt.1) &
         call AbortWithError('Error in ResizeV(): nrk < 1')

      call ReplaceVwithW(w,v)
      call NewCPvec(v,w%nbas,nrk)
      call GenCopyWtoV(v,w,1,MIN(rkv,nrk),1,MIN(rkv,nrk))
      call FlushCPvec(w)

      END SUBROUTINE ResizeV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE TrimZeros(v,tol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resize vector v by trimming columns with coefficients smaller than tol
! If all columns are zero, then a rank-1 zero vector is returned

      IMPLICIT NONE
      TYPE (CPvec), INTENT(INOUT) :: v
      TYPE (CPvec) :: w
      real*8, intent(in)   :: tol
      integer, allocatable :: iok(:)
      integer :: i,nok,rk

      rk=SIZE(v%coef)
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
            call NewCPvec(w,v%nbas,nok)
            DO i=1,nok
               call GenCopyWtoV(w,v,i,i,iok(i),iok(i))
            ENDDO
         ELSE
            call GetZeroCPvec(w,v%nbas)
         ENDIF
         call ReplaceVwithW(v,w)
      ENDIF

      DEALLOCATE(iok)

      END SUBROUTINE TrimZeros

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CP2DtoMat(v,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out 2D CP-format vector 'v' to give matrix 'M' on output

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN) :: v
      REAL*8, ALLOCATABLE, INTENT(OUT) :: M(:,:)
      INTEGER :: i,j,k,nr,nc,nrk
      REAL*8  :: fac

      IF (SIZE(v%nbas).ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoMat()')
      ENDIF

      nr=v%nbas(1)
      nc=v%nbas(2)
      nrk=SIZE(v%coef)

      ALLOCATE(M(nr,nc))
      M=0.d0
      DO i=1,nrk
         IF (nc.le.nr) THEN
            DO k=1,nc
               fac=v%coef(i)*v%base(nr+k,i)
               DO j=1,nr
                  M(j,k)=M(j,k)+fac*v%base(j,i)
               ENDDO
            ENDDO
         ELSE
            DO j=1,nr
               fac=v%coef(i)*v%base(j,i)
               DO k=1,nc
                  M(j,k)=M(j,k)+fac*v%base(nr+k,i)
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      END SUBROUTINE CP2DtoMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CP2DtoUW(v,U,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts directional component matrices U and W from 2D CP-format 
! vector 'v' 

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN) :: v
      REAL*8, ALLOCATABLE, INTENT(OUT) :: U(:,:),W(:,:)
      INTEGER :: i,j,nu,nw,nrk

      IF (SIZE(v%nbas).ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoUW()')
      ENDIF

      nu=v%nbas(1)
      nw=v%nbas(2)
      nrk=SIZE(v%coef)

      ALLOCATE(U(nu,nrk),W(nw,nrk))
      U(:,:)=v%base(1:nu,:)
      W(:,:)=v%base(nu+1:nu+nw,:)

!     Multiply the coef by the base with fewer elements
      DO i=1,nrk
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

      END SUBROUTINE CP2DtoUW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE MultOutCoef(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out coefficient using base for DOF with the smallest basis.
! Coefficients are then set to 1.0

      IMPLICIT NONE
      TYPE (CPvec), INTENT(INOUT) :: v
      INTEGER :: i,gst,gi,gf,imode(1)

!     Find the mode with the smallest basis and the range of this mode
!     in the base
      imode=MINLOC(v%nbas)
      gst=0
      IF (imode(1).gt.1) THEN
         DO i=2,imode(1)
            gst=gst+v%nbas(i-1)
         ENDDO
      ENDIF
      gi=gst+1
      gf=gst+v%nbas(imode(1))

!     Multiply out coefficient
      DO i=1,SIZE(v%coef)
         v%base(gi:gf,i)=v%base(gi:gf,i)*v%coef(i)
         v%coef(i)=1.d0
      ENDDO

      END SUBROUTINE MultOutCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE DistributeCoef(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies coef by base for all DOFs, and resets coefficients to 1.0

      IMPLICIT NONE
      TYPE (CPvec), INTENT(INOUT) :: v
      REAL*8  :: fac
      INTEGER :: i,ndof

      ndof=SIZE(v%nbas)

!     Multiply out coefficient
      DO i=1,SIZE(v%coef)
         fac=v%coef(i)**(1.d0/ndof)
         v%base(:,i)=fac*v%base(:,i)
         v%coef(i)=1.d0
      ENDDO

      END SUBROUTINE DistributeCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION CHECKNBASvv(v1,v2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN) :: v1,v2
      LOGICAL :: CHECKNBASvv
      INTEGER :: ndof,i

      CHECKNBASvv=.TRUE.

      IF (SIZE(v1%nbas).ne.SIZE(v2%nbas)) THEN 
         CHECKNBASvv=.FALSE.
      ELSE
         DO i=1,SIZE(v1%nbas)
            IF (v1%nbas(i).ne.v2%nbas(i)) THEN
               CHECKNBASvv=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      END FUNCTION CHECKNBASvv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION CHECKNBASvQ(v,Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      IMPLICIT NONE
      TYPE (CPblock), INTENT(IN) :: Q
      TYPE (CPvec), INTENT(IN)   :: v
      LOGICAL :: CHECKNBASvQ
      INTEGER :: ndof,i

      CHECKNBASvQ=.TRUE.

      IF (SIZE(v%nbas).ne.SIZE(Q%nbas)) THEN
         CHECKNBASvQ=.FALSE.
      ELSE
         DO i=1,SIZE(v%nbas)
            IF (v%nbas(i).ne.Q%nbas(i)) THEN
               CHECKNBASvQ=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      END FUNCTION CHECKNBASvQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION CHECKNBASQQ(Q1,Q2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      IMPLICIT NONE
      TYPE (CPblock), INTENT(IN) :: Q1,Q2
      LOGICAL :: CHECKNBASQQ
      INTEGER :: ndof,i

      CHECKNBASQQ=.TRUE.

      IF (SIZE(Q1%nbas).ne.SIZE(Q2%nbas)) THEN
         CHECKNBASQQ=.FALSE.
      ELSE
         DO i=1,SIZE(Q1%nbas)
            IF (Q1%nbas(i).ne.Q2%nbas(i)) THEN
               CHECKNBASQQ=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      END FUNCTION CHECKNBASQQ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      FUNCTION GetVfromBlock(Q,l)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extract a CP-vector from the block. 'l' holds the list of quantum
! numbers for the vector to be extracted    

      IMPLICIT NONE
      TYPE (CPblock), INTENT(IN) :: Q
      TYPE (CPvec) :: GetVfromBlock
      INTEGER, INTENT(IN)  :: l(:)
      INTEGER, ALLOCATABLE :: nbas(:)
      INTEGER :: ndof,nrk,n,i,i2,j,k,r
      REAL*8  :: coef

      ndof=SIZE(Q%nbas)
      nrk=SIZE(Q%coef,2)

!     Error checking
      IF (size(l).ne.ndof) THEN
         write(*,*) 'Mismatch in # DOFS between l and Q'
         call AbortWithError('Error in GetVfromBlock()')
      ENDIF

      ALLOCATE(nbas(ndof))

      DO i=1,ndof
         IF (l(i).le.0 .or. l(i).gt.Q%nbas(i)) THEN
            write(*,*) 'Extracted vector exceeds basis: dof: ',&
                       i,'; vec #: ',l(i)
            call AbortWithError('Error in GetVfromBlock()')
         ENDIF
         nbas(i)=nint(sqrt(real(Q%nbas(i))))
      ENDDO

      call NewCPvec(GetVfromBlock,nbas,nrk)
      DO r=1,nrk
         coef=1.d0
         j=0
         k=0
         DO i=1,ndof
            n=nbas(i)
            coef=coef*Q%coef(j+l(i),r)
            DO i2=1,n
               GetVfromBlock%base(j+i2,r)=Q%base(k+(i2-1)*n+l(i),r)
            ENDDO
            j=j+n
            k=k+n*n
         ENDDO
         GetVfromBlock%coef(r)=coef
      ENDDO

      DEALLOCATE(nbas)

      END FUNCTION GetVfromBlock

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NewCPblock(Q,nbas,rnk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Initializes a new block vector in CP-format

      IMPLICIT NONE
      TYPE (CPblock)      :: Q
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rnk
      INTEGER :: ndof,ndim,nrdim,i

      IF (rnk.lt.1) THEN
         write(*,*) 'Error: rnk must be at least 1!'
         call AbortWithError('Error in NewCPblock()')
      ENDIF

!     The structure of the base has strides of nbas(i)^2 elements
      ndof=0
      ndim=0
      nrdim=0
      DO i=1,SIZE(nbas)
         IF (nbas(i).lt.1) EXIT
         ndim=ndim+nbas(i)
         nrdim=nrdim+nbas(i)*nbas(i)
         ndof=ndof+1
      ENDDO

      IF (ndof.lt.1) THEN
         write(*,*) 'Error: no degrees of freedom!'
         call AbortWithError('Error in NewCPblock()')
      ENDIF

      ALLOCATE(Q%base(nrdim,rnk),Q%coef(ndim,rnk),Q%nbas(ndof))
      DO i=1,ndof
         Q%nbas(i)=nbas(i)*nbas(i)
      ENDDO

      END SUBROUTINE NewCPblock

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPblock(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format

      implicit none
      TYPE (CPblock), INTENT(IN) :: Q
      integer :: r,rk,j,d,i,k,n,inr,jnr
      character*64 :: frmt

      rk=SIZE(Q%coef,2)
      d=SIZE(Q%nbas)

      DO r=1,rk
         inr=0
         jnr=0
         DO j=1,d
            write(*,'(/2(A,I0))') 'RANK: ',r,'; dof: ',j
            n=nint(sqrt(real(Q%nbas(j))))
            write(frmt,'(A,I0,A)') '(A,',n,'(X,f14.6))'
            write(*,frmt) 'Coef =',(Q%coef(jnr+k,r),k=1,n)
            write(frmt,'(A,I0,A)') '(6X',n,'(X,f14.6))'
            DO i=1,n
               write(*,frmt) (Q%base(inr+i+(k-1)*n,r),k=1,n)
            ENDDO
            inr=inr+Q%nbas(j)
            jnr=jnr+n
         ENDDO
      ENDDO
      write(*,*)

      end subroutine PrintCPblock

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE FlushCPblock(Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Disposes CP-format block of vectors

      IMPLICIT NONE
      TYPE (CPblock) :: Q

      IF (ALLOCATED(Q%base)) DEALLOCATE(Q%base)
      IF (ALLOCATED(Q%coef)) DEALLOCATE(Q%coef)
      IF (ALLOCATED(Q%nbas)) DEALLOCATE(Q%nbas)

      END SUBROUTINE FlushCPblock

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE SEPDREPN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
