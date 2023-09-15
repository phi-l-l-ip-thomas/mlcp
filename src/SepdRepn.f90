!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE SEPDREPN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains structures needed for CP representation

      USE ERRORTRAP
      USE UTILS
      USE, INTRINSIC :: ISO_C_BINDING

      TYPE CP
         REAL*8,  ALLOCATABLE :: base(:,:), coef(:)
         INTEGER, ALLOCATABLE :: nbas(:), ibas(:), fbas(:)
         INTEGER, ALLOCATABLE :: rows(:), cols(:)
         LOGICAL, ALLOCATABLE :: sym(:)
         CONTAINS
            PROCEDURE :: R => Getrank
            PROCEDURE :: D => Getndof
            PROCEDURE :: M => Getrows
            PROCEDURE :: N => Getcols
            PROCEDURE :: show => CPShowStats
            PROCEDURE :: print => PrintCPmat_all
            PROCEDURE :: printvec => PrintCPvec
      END TYPE CP

      INTERFACE NewCP
         MODULE PROCEDURE NewCPgen,NewCPref,NewCPvec,NewCPsqmat
      END INTERFACE NewCP

      INTERFACE RandomCP
         MODULE PROCEDURE RandomCPgen,RandomCPref
      END INTERFACE RandomCP

      INTERFACE PrintCPmat
         MODULE PROCEDURE PrintCPmat_all,PrintCPmat_gen
      END INTERFACE PrintCPmat

      INTERFACE CopyCP
         MODULE PROCEDURE CopyCP_all,ExtractCPsubmatrix 
      END INTERFACE CopyCP

      INTERFACE CPMatrixZeroOffDiag
         MODULE PROCEDURE CPMatrixZeroOffDiag_all
         MODULE PROCEDURE CPMatrixZeroOffDiag_one
         MODULE PROCEDURE CPMatrixZeroOffDiag_gen
      END INTERFACE CPMatrixZeroOffDiag

      INTERFACE MultOutCoef
         MODULE PROCEDURE MultOutCoefSmallest,MultOutCoefbyMode
      END INTERFACE MultOutCoef

      INTERFACE GetRank1DominantEntry
         MODULE PROCEDURE GetRank1DominantEntry_gen
         MODULE PROCEDURE GetRank1DominantEntry_1
      END INTERFACE GetRank1DominantEntry

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPgen(rk,rows,cols,sym) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! General subroutine for initializing outer product of CP-vectors or
! matrices. Symmetric storage is possible for square matrix factors

      implicit none
      TYPE (CP) :: v
      LOGICAL, INTENT(IN), OPTIONAL :: sym(:)
      INTEGER, INTENT(IN) :: rows(:), cols(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER :: ndof,nrdim,i

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

      IF (present(sym)) THEN
         IF (SIZE(sym).ne.ndof) THEN
            write(*,*) 'Error: size of "sym" array differs from ndof!'
            call AbortWithError('Error in NewCP()')
         ENDIF
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

         IF (present(sym)) THEN
            IF (sym(i) .and. (rows(i).ne.cols(i))) THEN
               write(*,*) 'Error: rows must equal cols to use symmetry!'
               call AbortWithError('Error in NewCP()')
            ENDIF
         ENDIF
      ENDDO

!     Allocate arrays
      ALLOCATE(v%nbas(ndof),v%ibas(ndof),v%fbas(ndof))
      ALLOCATE(v%rows(ndof),v%cols(ndof),v%sym(ndof))
      v%rows(:)=rows(:)
      v%cols(:)=cols(:)
      IF (present(sym)) THEN
         v%sym(:)=sym(:)
      ELSE
         v%sym(:)=.FALSE.
      ENDIF

      nrdim=0
      DO i=1,ndof
         v%ibas(i)=nrdim+1
         IF (v%sym(i)) THEN
            v%nbas(i)=rows(i)*(rows(i)+1)/2
         ELSE
            v%nbas(i)=rows(i)*cols(i)
         ENDIF
         nrdim=nrdim+v%nbas(i)
         v%fbas(i)=nrdim
      ENDDO

      ALLOCATE(v%base(nrdim,rk),v%coef(rk))

      end function NewCPgen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPref(w,rk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes CP-matrix v, using sizes in w (including the rank, if not
! passed as an optional argument)

      implicit none
      TYPE (CP) :: v
      TYPE (CP), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk

      IF (present(rk)) THEN
         v=NewCP(rk,w%rows,w%cols,w%sym)
      ELSE
         v=NewCP(SIZE(w%coef),w%rows,w%cols,w%sym)
      ENDIF

      end function NewCPref

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPvec(rk,rows) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a new vector in CP-format

      implicit none
      TYPE (CP) :: v
      INTEGER, INTENT(IN)  :: rows(:)
      INTEGER, INTENT(IN)  :: rk
      INTEGER, ALLOCATABLE :: cols(:)
      LOGICAL, ALLOCATABLE :: sym(:)
      INTEGER   :: ndof

      ndof=SIZE(rows)
      ALLOCATE(cols(ndof),sym(ndof))
      cols(:)=1
      sym(:)=.FALSE.

      v=NewCP(rk,rows,cols,sym)

      DEALLOCATE(cols,sym)

      end function NewCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function NewCPsqmat(rk,rows,symm) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes a CP-format outer product of square matrices, where
! symmetry is used on either-all-or-none of the DOFs

      implicit none
      TYPE (CP) :: v
      LOGICAL, INTENT(IN) :: symm
      INTEGER, INTENT(IN) :: rows(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER, ALLOCATABLE :: cols(:)
      LOGICAL, ALLOCATABLE :: sym(:)
      INTEGER   :: ndof

      ndof=SIZE(rows)
      ALLOCATE(cols(ndof),sym(ndof))
      cols(:)=rows(:)
      sym(:)=symm

      v=NewCP(rk,rows,cols,sym)

      DEALLOCATE(cols,sym)

      end function NewCPsqmat
              
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushCP(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes CP-format type

      implicit none
      TYPE (CP) :: v

      IF (ALLOCATED(v%base)) DEALLOCATE(v%base)
      IF (ALLOCATED(v%coef)) DEALLOCATE(v%coef)
      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)
      IF (ALLOCATED(v%ibas)) DEALLOCATE(v%ibas)
      IF (ALLOCATED(v%fbas)) DEALLOCATE(v%fbas)
      IF (ALLOCATED(v%rows)) DEALLOCATE(v%rows)
      IF (ALLOCATED(v%cols)) DEALLOCATE(v%cols)
      IF (ALLOCATED(v%sym)) DEALLOCATE(v%sym)

      end subroutine FlushCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPShowStats(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Shows mode sizesm sym, rank of v

      implicit none
      CLASS (CP), INTENT(IN) :: v
      integer :: d
      character*64 :: frmt

      frmt='(X,I4,X,I5,X,X,I5,2X,L3)'
      write(*,*) '------- CP stats -------'
      write(*,'(X,A,I4,A,I0)') 'ndof = ',v%D(),'; rank = ',v%R()
      write(*,*) 'mode [rows x cols] sym'
      do d=1,v%D()
         write(*,frmt) d,v%rows(d),v%cols(d),v%sym(d)
      enddo
      write(*,*) '------------------------'

      end subroutine CPShowStats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPvec(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-vector in neat format

      implicit none
      CLASS (CP), INTENT(IN) :: v
      integer :: r,rk,j,d,i,n,inr
      character*64 :: frmt

      rk=SIZE(v%coef)
      d=SIZE(v%nbas)

!      write(frmt,'(A,I0,A)') '(A,X,',rk,'f18.10)'
      write(frmt,'(A,I0,A)') '(A,2X,',rk,'(ES17.10,X))'
      write(*,frmt) ' Vcoef =',(v%coef(r),r=1,rk)
      write(frmt,'(A,I0,A)') '(2(I4),X,',rk,'f18.10)'
!      write(frmt,'(A,I0,A)') '(2(I4),X,',rk,'(E23.16,X))'
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

      subroutine PrintCPmat_all(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format, wrapper for usual case

      implicit none
      CLASS (CP), INTENT(IN) :: v

      call PrintCPmat_gen(v,.TRUE.)

      end subroutine PrintCPmat_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPmat_gen(v,showbase)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints CP-matrix in neat format, general routine

      implicit none
      CLASS (CP), INTENT(IN) :: v
      logical, intent(in) :: showbase
      integer :: r,rk,j,d,i,n,k,m,inr
      character*64 :: frmt

      rk=SIZE(v%coef)
      d=SIZE(v%nbas)

      write(*,*)
      DO r=1,rk
!         write(*,'(A,I0,A,F23.16)') 'RANK: ',r,'; coef = ',v%coef(r)
         write(*,'(A,I0,A,ES23.16)') 'RANK: ',r,'; coef = ',v%coef(r)
         IF (showbase) THEN
            inr=0
            DO j=1,d
               write(*,'(/A,I0)') 'dof : ',j
               n=v%rows(j)
               m=v%cols(j)
               IF (v%sym(j)) THEN
                  DO i=1,n
                     write(frmt,'(A,I0,A)') '(',i,'(X,f14.6))'
                     write(*,frmt) &
                          (v%base(inr+i+(k-1)*n-k*(k-1)/2,r),k=i,1,-1)
                  ENDDO
               ELSE
                  write(frmt,'(A,I0,A)') '(',m,'(X,f14.6))'
                  DO i=1,n
                     write(*,frmt) (v%base(inr+i+(k-1)*n,r),k=1,m)
                  ENDDO
               ENDIF
               inr=inr+v%nbas(j)
            ENDDO
            write(*,*)
         ENDIF
      ENDDO
      write(*,*)

      end subroutine PrintCPmat_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getrank(v) result(rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP), INTENT(IN) :: v
      integer :: rk

      IF (ALLOCATED(v%coef)) THEN
         rk=SIZE(v%coef)
      ELSE
         rk=0
      ENDIF

      end function Getrank

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getndof(v) result(ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP), INTENT(IN) :: v
      integer :: ndof

      IF (ALLOCATED(v%nbas)) THEN
         ndof=SIZE(v%nbas)
      ELSE
         ndof=0
      ENDIF

      end function Getndof

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getrows(v,d) result(rows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP), INTENT(IN) :: v
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

      end function Getrows

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function Getcols(v,d) result(cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns number of modes of v

      implicit none
      CLASS (CP), INTENT(IN) :: v
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

      end function Getcols

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ZeroCPvec(nbas) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a zero CP-vec of rank 1

      implicit none
      TYPE (CP) :: v
      INTEGER, INTENT(IN) :: nbas(:)

      v=NewCP(1,nbas)
      v%base=0.d0
      v%coef=0.d0

      end function ZeroCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RandomCPref(w,rk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries, with dimensions of reference w

      implicit none
      TYPE (CP) :: v
      TYPE (CP), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk
      INTEGER :: i,rv,d
      REAL*8  :: fac

      IF (present(rk)) THEN
         rv=rk
      ELSE
         rv=w%R()
      ENDIF

!     Generate v with random entries and equal coefs for all terms
      v=NewCP(rv,w%rows,w%cols,w%sym)
      v%coef(:)=1.d0/sqrt(REAL(rv))
      call random_number(v%base(:,:))

!     Shift, scale entries for each mode to make rms norm ~ unity
      DO d=1,v%D()
         fac=sqrt(12.d0/REAL(v%nbas(d)))
         v%base(v%ibas(d):v%fbas(d),:)=fac*&
         (v%base(v%ibas(d):v%fbas(d),:)-0.5d0)
      ENDDO

      end function RandomCPref

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RandomCPgen(rk,rows,cols,sym) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries. No normalization is done.

      implicit none
      TYPE (CP) :: v
      LOGICAL, INTENT(IN), OPTIONAL :: sym(:)
      INTEGER, INTENT(IN) :: rows(:),cols(:)
      INTEGER, INTENT(IN) :: rk
      INTEGER :: i

      IF (present(sym)) THEN
         v=NewCP(rk,rows,cols,sym)
      ELSE
         v=NewCP(rk,rows,cols)
      ENDIF

      DO i=1,rk
         call random_number(v%base(:,i))
         v%coef(i)=1.d0
      ENDDO
      v%base=v%base-0.5d0

      end function RandomCPgen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function IdentityCPMatrix(rows,cols,sym) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets rank-1 CP outer-product-of-identity-matrices. If only the
! upper triangle of Hsym is stored, set sym to .TRUE.

      implicit none
      TYPE (CP) :: v
      integer, intent(in)  :: rows(:),cols(:)
      logical, intent(in)  :: sym(:)
      integer :: i,j,k,ndof

      ndof=SIZE(rows)

      v=NewCP(1,rows,cols,sym)
      v%coef=1.d0
      v%base=0.d0

!     Get an identity matrix for each DOF
      DO j=1,ndof
          IF (sym(j)) THEN
             v%base(v%ibas(j):v%ibas(j)+v%rows(j)-1,1)=1.d0
          ELSE
             k=v%ibas(j)
             DO i=1,MIN(v%rows(j),v%cols(j))
                v%base(k,1)=1.d0
                k=k+v%rows(j)+1
             ENDDO
          ENDIF
      ENDDO

      end function IdentityCPMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPVtoDiagMatrix(v,thesym) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates square CP-matrix with elements of v on the diagonal

      implicit none
      TYPE (CP), INTENT(IN)  :: v
      TYPE (CP) :: w
      logical, intent(in), optional :: thesym(:)
      logical, allocatable :: sym(:)
      integer :: d,ndof,i,rk,j,nbas,vst,wst,nbo

      ndof=SIZE(v%nbas)
      rk=SIZE(v%coef)

      ALLOCATE(sym(ndof))
      IF (present(thesym)) THEN
         sym(:)=thesym(:)
      ELSE
         sym(:)=.FALSE.
      ENDIF

      w=NewCP(rk,v%nbas,v%nbas,sym)
      w%base(:,:)=0.d0
      w%coef(:)=v%coef(:)

      DO d=1,ndof
         nbas=v%nbas(d)
         nbo=v%nbas(d)+1
         DO i=1,rk
            IF (sym(d)) THEN
!              Since symmetry stores the diagonal first, the elements of
!              v map directly to the first elements of w, per coordinate
               w%base(w%ibas(d):w%ibas(d)+nbas-1,i)=&
               v%base(v%ibas(d):v%fbas(d),i)
            ELSE
!              Copy elements of v to diagonal of w
               wst=w%ibas(d)
               vst=v%ibas(d)
               DO j=0,nbas-1
                  w%base(wst,i)=v%base(vst,i)
                  wst=wst+nbo
                  vst=vst+1
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(sym)

      end function CPVtoDiagMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractDiagfromCPMatrix(v,trans) result(w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts diagonal from CP-matrix (must be square); returns a column
! vector. Set trans to true to get a row vector.

      implicit none
      TYPE (CP), INTENT(IN) :: v
      TYPE (CP) :: w
      logical, intent(in)  :: trans
      logical, allocatable :: sym(:)
      integer, allocatable :: one(:)
      integer :: d,ndof,i,rk,j,m,n

      rk=SIZE(v%coef)
      ndof=SIZE(v%nbas)

      allocate(one(ndof),sym(ndof))
      one(:)=1
      sym(:)=.FALSE.

      IF (trans) THEN
         w=NewCP(rk,one,v%cols,sym)
      ELSE
         w=NewCP(rk,v%rows,one,sym)
      ENDIF

!     Base is taken from diagonal elements of square matrix
      DO d=1,ndof
         IF (v%rows(d).ne.v%cols(d)) THEN
            write(*,*) 'Matrix for mode ',d,' is (',v%rows(d),' x ',&
            v%cols(d),') but must be square'
            call AbortWithError('ExtractDiagfromCPMatrix: matrix dims')
         ENDIF
         n=v%rows(d)
         IF (v%sym(d)) THEN
!           First n elements are the diagonal when symmetry is used
            w%base(w%ibas(d):w%fbas(d),:)=&
            v%base(v%ibas(d):v%ibas(d)+n-1,:)
         ELSE
            j=v%ibas(d)
            m=v%rows(d)+1
            DO i=0,n-1
               w%base(w%ibas(d)+i,:)=v%base(j,:)
               j=j+m
            ENDDO
         ENDIF
      ENDDO

!     Copy coefs directly
      w%coef(:)=v%coef(:)

      deallocate(one,sym)

      end function ExtractDiagfromCPMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatrixTranspose(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Transposes CP-matrix by reordering base

      implicit none
      TYPE (CP), intent(inout) :: v
      integer :: d,i,j,k,ik,tmp,tmp1
      real*8, allocatable :: btmp(:)

      DO d=1,SIZE(v%nbas)
!        Nothing needed for sym=.true.!
         IF (.not.v%sym(d)) THEN
!           If either dimension is 1, just swap row and col numbers
            IF ((v%rows(d).gt.1) .and. (v%cols(d).gt.1)) THEN
!              Reorder the elements in temporary array
               allocate(btmp(v%nbas(d)))
               DO j=1,SIZE(v%coef)
                  ik=1
                  tmp=v%ibas(d)
                  DO i=0,v%rows(d)-1
                     tmp1=tmp
                     DO k=0,v%cols(d)-1
                        btmp(ik)=v%base(tmp1,j)
                        ik=ik+1
                        tmp1=tmp1+v%rows(d)
                     ENDDO
                     tmp=tmp+1
                  ENDDO
                  v%base(v%ibas(d):v%fbas(d),j)=btmp(1:v%nbas(d))
               ENDDO
               deallocate(btmp)
            ENDIF
            tmp=v%rows(d)
            v%rows(d)=v%cols(d)
            v%cols(d)=tmp
         ENDIF
      ENDDO

      end subroutine CPMatrixTranspose

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatrixZeroOffDiag_all(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for all modes

      implicit none
      TYPE (CP), intent(inout) :: v
      logical, allocatable :: domode(:)

      ALLOCATE(domode(v%D()))
      domode(:)=.TRUE.
      call CPMatrixZeroOffDiag(v,domode)
      DEALLOCATE(domode)

      end subroutine CPMatrixZeroOffDiag_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatrixZeroOffDiag_one(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for mode d

      implicit none
      TYPE (CP), intent(inout) :: v
      integer, intent(in)  :: d
      logical, allocatable :: domode(:)

      IF ((d.lt.1) .or. (d.gt.v%D())) THEN
         write(*,*) 'mode d (',d,') must be in range: [1,',v%D(),']'
         call AbortWithError('CPMatrixZeroOffDiag(): d out of range')
      ENDIF

      ALLOCATE(domode(v%D()))
      domode(:)=.FALSE.
      domode(d)=.TRUE.
      call CPMatrixZeroOffDiag(v,domode)
      DEALLOCATE(domode)

      end subroutine CPMatrixZeroOffDiag_one

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMatrixZeroOffDiag_gen(v,domode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Zeros off-diagonal CP-matrix elements for each mode matrix if domode
! is set. Mode matrix does not need to be square.

      implicit none
      TYPE (CP), intent(inout) :: v
      logical, intent(in) :: domode(:)
      integer :: d,ndof,i,j,k

      IF (SIZE(domode).ne.SIZE(v%nbas)) THEN
         write(*,*) 'domode has ',SIZE(domode),' entries but ',&
                    'must have ndof (',SIZE(v%nbas),') entries'
         call AbortWithError('CPMatrixZeroOffDiag(): ndof mismatch')
      ENDIF

      ndof=SIZE(v%nbas)
      DO d=1,ndof
         IF (domode(d)) THEN
            IF (v%sym(d)) THEN
!              Symmetry: skip diagonal entries, which are listed first
               v%base(v%ibas(d)+v%rows(d):v%fbas(d),:)=0.d0
            ELSE
               k=v%ibas(d)
               DO i=1,v%cols(d)
                  DO j=1,v%rows(d)
                     IF (i.ne.j) v%base(k,:)=0.d0
                     k=k+1
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO

      end subroutine CPMatrixZeroOffDiag_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReplaceVwithW(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, disposing W afterwards

      implicit none
      TYPE (CP) :: v,w

      call FlushCP(v)
      v=CopyCP(w)
      call FlushCP(w)

      end subroutine ReplaceVwithW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CopyCP_all(w) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      implicit none
      TYPE (CP), INTENT(IN) :: w
      TYPE (CP) :: v
      INTEGER   :: rk

      rk=SIZE(w%coef)
      v=NewCP(rk,w%rows,w%cols,w%sym)
      call GenCopyWtoV(v,w,1,rk,1,rk)

      end function CopyCP_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GenCopyWtoV(v,w,vi,ve,wi,we)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V
! leaving W intact. v must be allocated beforehand.

      implicit none
      TYPE (CP) :: v,w
      INTEGER, INTENT(IN) :: vi,ve,wi,we
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

      end subroutine GenCopyWtoV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ResizeV(v,rk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resizes vector v, either by truncating at a smaller rank or by
! adding space for extra terms.

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP) :: w
      INTEGER, INTENT(IN) :: rk
      INTEGER :: rkv

      rkv=SIZE(v%coef)

      IF (rk.lt.1) &
         call AbortWithError('Error in ResizeV(): rk < 1')

      call ReplaceVwithW(w,v)
      v=NewCP(rk,w%rows,w%cols,w%sym)
      call GenCopyWtoV(v,w,1,MIN(rkv,rk),1,MIN(rkv,rk))
      call FlushCP(w)

      end subroutine ResizeV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TrimZeros(v,tol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resize vector v by trimming columns with coefs smaller than tol
! If all columns are zero, then a rank-1 zero vector is returned

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      TYPE (CP) :: w
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
            w=NewCP(rk,v%rows,v%cols,v%sym)
            DO i=1,nok
               call GenCopyWtoV(w,v,i,i,iok(i),iok(i))
            ENDDO
         ELSE
            w=NewCP(1,v%rows,v%cols,v%sym)
            w%base=0.d0
            w%coef=0.d0
         ENDIF
         call ReplaceVwithW(v,w)
      ENDIF

      DEALLOCATE(iok)

      end subroutine TrimZeros

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CP2DtoMat(v) result(M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out 2D CP-format vector 'v' to give matrix 'M'

      implicit none
      TYPE (CP), INTENT(IN) :: v
      REAL*8, ALLOCATABLE :: M(:,:)
      INTEGER :: i,j,k,nr,nc,rk
      REAL*8  :: fac

      IF (SIZE(v%nbas).ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoMat()')
      ENDIF

      nr=v%nbas(1)
      nc=v%nbas(2)
      rk=SIZE(v%coef)

      ALLOCATE(M(nr,nc))
      M=0.d0
      DO i=1,rk
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

      end function CP2DtoMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CP2DtoUW(v,U,W)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts directional component matrices U and W from 2D CP-format 
! vector 'v' 

      implicit none
      TYPE (CP), INTENT(IN) :: v
      REAL*8, ALLOCATABLE, INTENT(OUT) :: U(:,:),W(:,:)
      INTEGER :: i,j,nu,nw,rk

      IF (SIZE(v%nbas).ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-to-matrix transform'
         call AbortWithError('Error in CP2DtoUW()')
      ENDIF

      nu=v%nbas(1)
      nw=v%nbas(2)
      rk=SIZE(v%coef)

      ALLOCATE(U(nu,rk),W(nw,rk))
      U(:,:)=v%base(1:nu,:)
      W(:,:)=v%base(nu+1:nu+nw,:)

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

      end subroutine CP2DtoUW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefSmallest(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies out coefficient using base for DOF with the smallest
! basis. Coefficients are then set to 1.0

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      INTEGER :: imode(1)

!     Pick the mode with the smallest basis (fewest multiplies)
      imode=MINLOC(v%nbas)
      call MultOutCoefbyMode(v,imode(1))

      end subroutine MultOutCoefSmallest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MultOutCoefbyMode(v,d)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies base of mode 'd' by coefficients; resets coefficients to 1

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      integer, intent(in) :: d
      integer :: i

!     Error checking
      IF ((d.lt.1) .or. (d.gt.SIZE(v%nbas))) THEN
         write(*,*) 'Mode ',d,' must be in range [1,',SIZE(v%nbas),']'
         call AbortWithError('MultOutCoefbyMode(): d out of range')
      ENDIF

!     Multiply out coefficients
      DO i=1,SIZE(v%coef)
         v%base(v%ibas(d):v%fbas(d),i)=v%coef(i)*&
         v%base(v%ibas(d):v%fbas(d),i)
         v%coef(i)=1.d0
      ENDDO

      end subroutine MultOutCoefbyMode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DistributeCoef(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies coef by base for all DOFs, and resets coefs to 1.0
! Each mode's base is scaled by the mode size

      implicit none
      TYPE (CP), INTENT(INOUT) :: v
      REAL*8, ALLOCATABLE :: pows(:)
      REAL*8  :: fac,div
      INTEGER :: d,i,ndof

      ndof=SIZE(v%nbas)

!     Get the exponents for modes with different nbas values      
      allocate(pows(ndof))
      div=0.d0
      DO d=1,ndof
         pows(d)=sqrt(REAL(v%nbas(d)))
         div=div+pows(d)
      ENDDO
      pows(:)=pows(:)/div
      
!     Multiply out coefficient
      DO i=1,SIZE(v%coef)
         DO d=1,ndof
            fac=v%coef(i)**pows(d) 
            v%base(v%ibas(d):v%fbas(d),i)=fac*&
            v%base(v%ibas(d):v%fbas(d),i)
         ENDDO
         v%coef(i)=1.d0
      ENDDO

      deallocate(pows)

      end subroutine DistributeCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_1(v,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      TYPE (CP), INTENT(IN) :: v
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real*8, intent(out) :: val

      IF (SIZE(v%coef).ne.1) THEN
         write(*,*) 'Rank of v (',SIZE(v%coef),') must be 1'
         call AbortWithError('Error in GetRank1DominantEntry_1()')
      ENDIF

      call GetRank1DominantEntry_gen(v,1,rowi,coli,val)

      end subroutine GetRank1DominantEntry_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetRank1DominantEntry_gen(v,irk,rowi,coli,val)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Get largest absolute value, indices of term irk in CP object v

      implicit none
      TYPE (CP), INTENT(IN) :: v
      integer, intent(in) :: irk
      integer, allocatable, intent(out) :: rowi(:),coli(:)
      real*8, intent(out) :: val
      integer :: d,ndof
      integer :: imx(1)
      real*8  :: vmx

      IF (irk.lt.1 .or. irk.gt.SIZE(v%coef)) THEN
         write(*,*) 'irk (',irk,') out of range: [1,',SIZE(v%coef),']'
         call AbortWithError('Error in GetRank1DominantEntry_gen()')
      ENDIF

      ndof=SIZE(v%nbas)
      ALLOCATE(rowi(ndof),coli(ndof))

      val=v%coef(irk)
      DO d=1,ndof
         IF (v%sym(d)) THEN
!!!         Careful with sym: 
!!!         1) array indexing is different
!!!         2) off-diag elements have been scaled
            call AbortWithError('GR1DE(): sym not yet implemented')
         ELSE
!           Use imx (range [1:v%nbas(d)]) to extract row,col indices
            imx=MAXLOC(ABS(v%base(v%ibas(d):v%fbas(d),irk)))
            rowi(d)=mod(imx(1)-1,v%rows(d))+1
            coli(d)=(imx(1)-1)/v%rows(d)+1
            vmx=v%base(v%ibas(d)-1+imx(1),irk)
         ENDIF
         val=val*vmx
      ENDDO

      end subroutine GetRank1DominantEntry_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractCPvec(M,indx,getcol) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-vector from CP-matrix M from index list, returns as v.
! Extracts a column by default, set getcol=.FALSE. to extract a row

      implicit none
      TYPE (CP), INTENT(IN)  :: M
      TYPE (CP)  :: v
      integer, intent(in)  :: indx(:)
      logical, intent(in)  :: getcol
      integer, allocatable :: one(:)
      logical, allocatable :: sym(:)
      integer :: i,j,k,l,d,st,fi

!     No symmetry since we want a vector
      allocate(sym(SIZE(M%nbas)),one(SIZE(M%nbas)))
      sym(:)=.FALSE.
      one(:)=1
      IF (getcol) THEN
         v=NewCP(SIZE(M%coef),M%rows,one,sym)
      ELSE
         v=NewCP(SIZE(M%coef),one,M%cols,sym)
      ENDIF
      deallocate(sym,one)

!     Extract base from M      
      DO d=1,SIZE(M%nbas)
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

         IF (M%sym(d)) THEN
            DO i=1,SIZE(M%coef)
               DO j=1,M%cols(d)
                  l=max(j,indx(d))
                  k=l-min(j,indx(d))+1
                  st=M%ibas(d)-1+l+(k-1)*M%rows(d)-k*(k-1)/2
                  fi=v%ibas(d)+j-1
                  v%base(fi,i)=M%base(st,i)
               ENDDO
            ENDDO
         ELSEIF (getcol) THEN
            st=M%ibas(d)+(indx(d)-1)*M%rows(d)
            fi=st+M%rows(d)-1
            DO i=1,SIZE(M%coef)
               v%base(v%ibas(d):v%fbas(d),i)=M%base(st:fi,i)
            ENDDO
         ELSE
            st=M%ibas(d)+indx(d)-1
            DO i=1,SIZE(M%coef)
               fi=0
               DO j=0,M%cols(d)-1
                  v%base(v%ibas(d)+j,i)=M%base(st+fi,i)
                  fi=fi+M%rows(d)
               ENDDO
            ENDDO
         ENDIF
      END DO

!     Coefs copy directly
      v%coef(:)=M%coef(:)

      end function ExtractCPvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractCPsubmatrix(M,irs,irf,ics,icf) result(V)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets CP-submatrix M from list of row, col starting and ending indices

      implicit none
      TYPE (CP), INTENT(IN), TARGET :: M
      TYPE (CP) :: V
      integer, intent(in)  :: irs(:),irf(:),ics(:),icf(:)
      integer, allocatable :: rows(:),cols(:)
      real*8, pointer :: mat(:,:)
      integer :: d,ndof,i,rk

      ndof=SIZE(M%nbas)
      rk=SIZE(M%coef)

!     Loads of error checking
      IF (any(M%sym)) &
         call AbortWithError('ExtractCPsubmatrix(): no sym implemented')
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
      V=NewCP(rk,rows,cols,M%sym)
      deallocate(rows,cols)

!     Copy coefs of M -> V directly, base by reshaping
      V%coef(:)=M%coef(:)
      DO d=1,ndof
         DO i=1,rk
!           Array 'mat' points to base of M
            call C_F_POINTER(C_LOC(M%base(M%ibas(d):M%fbas(d),i)),&
                             mat,[M%rows(d),M%cols(d)])
!           Copy selected portion of 'mat' into base of V
            V%base(V%ibas(d):V%fbas(d),i)=RESHAPE(&
                mat(irs(d):irf(d),ics(d):icf(d)),(/V%nbas(d)/))
         ENDDO
      ENDDO

      end function ExtractCPsubmatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutCPsubmatrix(V,irs,ics,M)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Puts CP-submatrix V into M, overwriting that portion of M, from the 
! row, col starting indices provided

      implicit none
      TYPE (CP), INTENT(IN) :: V
      TYPE (CP), INTENT(INOUT), TARGET :: M
      integer, intent(in)  :: irs(:),ics(:)
      integer, allocatable :: irf(:),icf(:)
      real*8, pointer :: mat(:,:)
      integer :: d,ndof,i,rk

      ndof=SIZE(V%nbas)
      rk=SIZE(V%coef)

!     Error checking extravaganza!
      IF (any(V%sym)) &
         call AbortWithError('PutCPsubmatrix(): no sym implemented (V)')
      IF (any(M%sym)) &
         call AbortWithError('PutCPsubmatrix(): no sym implemented (M)')

      IF (SIZE(M%nbas).ne.ndof) &
         call AbortWithError('PutCPsubmatrix(): M,V must have = ndof')
      IF (SIZE(M%coef).ne.rk) &
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
         DO i=1,rk
!           Array 'mat' points to base of M
            call C_F_POINTER(C_LOC(M%base(M%ibas(d):M%fbas(d),i)),&
                          mat,[M%rows(d),M%cols(d)])
!           Reshape V into selected portion of 'mat'
            mat(irs(d):irf(d),ics(d):icf(d))=&
               RESHAPE(V%base(V%ibas(d):V%fbas(d),i),&
               (/V%rows(d),V%cols(d)/))
         ENDDO
      ENDDO

      deallocate(irf,icf)

      end subroutine PutCPsubmatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function ExtractCPmatrixElement(M,ir,ic)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets element of CP-matrix M from list of row, col indices

      implicit none
      TYPE (CP), INTENT(IN) :: M
      integer, intent(in) :: ir(:),ic(:)
      integer :: d,ndof,i,rk,indx,k,l
      real*8 :: ExtractCPmatrixElement
      real*8, allocatable :: prod(:)
      real*8, parameter   :: s=1/sqrt(2.d0)

!!! TO BE TESTED (for symmetry only, no-sym already works)
      ndof=SIZE(M%nbas)
      rk=SIZE(M%coef)

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
         IF (M%sym(d)) THEN
            l=max(ir(d),ic(d))
            k=l-min(ir(d),ic(d))+1
            indx=M%ibas(d)-1+l+(k-1)*M%rows(d)-k*(k-1)/2
!           De-scale off-diagonal elements since they are stored
!           multiplied by the extra factor of sqrt(2)
            IF (l.ne.k) prod(i)=s*prod(i)
         ELSE        
            indx=M%ibas(d)-1+(ic(d)-1)*M%rows(d)+ir(d)
         ENDIF
         DO i=1,rk
            prod(i)=prod(i)*M%base(indx,i)
         ENDDO
      ENDDO

!     Now accumulate the sum-over-term-products
      ExtractCPmatrixElement=SUM(prod)

      DEALLOCATE(prod)

!!! END TO BE TESTED (symmetry only)
      end function ExtractCPmatrixElement

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPModeJoin(v,w) result(u)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Merges v and w by joining into a tensor of same rank with the sum of
! the dimensions

      implicit none
      TYPE (CP), INTENT(IN) :: v,w
      TYPE (CP) :: u
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: sym(:)
      integer :: rk,ndof,ndofv,ndofw,vs,ws

      IF (SIZE(v%coef).ne.SIZE(w%coef)) THEN
         write(*,*) 'Ranks of v,w (',SIZE(v%coef),',',SIZE(v%coef),&
                    ') must match'
         call AbortWithError('CPModeJoin(): v,w rank mismatch')
      ENDIF

!     Construct u with same rank and sum of mode sizes from v and w
      rk=SIZE(v%coef)
      ndofv=SIZE(v%nbas)
      ndofw=SIZE(w%nbas)
      ndof=ndofv+ndofw
      vs=SIZE(v%base,1)
      ws=SIZE(w%base,1)

      ALLOCATE(rows(ndof),cols(ndof),sym(ndof))

      rows(1:ndofv)=v%rows(:)
      rows(ndofv+1:ndofv+ndofw)=w%rows(:)
      cols(1:ndofv)=v%cols(:)
      cols(ndofv+1:ndofv+ndofw)=w%cols(:)
      sym(1:ndofv)=v%sym(:)
      sym(ndofv+1:ndofv+ndofw)=w%sym(:)

      u=NewCP(rk,rows,cols,sym)
      u%coef(:)=v%coef(:)*w%coef(:)
      u%base(1:vs,:)=v%base(:,:)
      u%base(vs+1:vs+ws,:)=w%base(:,:)

      DEALLOCATE(rows,cols,sym)

      end function CPModeJoin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPEntrywiseCompare(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares F and G, entry by entry.
! WARNING: runtime is exponential in ndof

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
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
         valf=ExtractCPmatrixelement(F,ir,ic)
         valg=ExtractCPmatrixelement(G,ir,ic)
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

      end subroutine CPEntrywiseCompare

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CHECKNBAS(v1,v2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks nbas of 2 CP-format vectors to make sure they are the same

      implicit none
      TYPE (CP), INTENT(IN) :: v1,v2
      LOGICAL :: CHECKNBAS
      INTEGER :: i

      CHECKNBAS=.TRUE.

      IF (SIZE(v1%rows).ne.SIZE(v2%rows)) THEN 
         CHECKNBAS=.FALSE.
      ELSE
         DO i=1,SIZE(v1%nbas)
            IF ((v1%rows(i).ne.v2%rows(i)) .or. &
                (v1%cols(i).ne.v2%cols(i)) .or. &
                (v1%sym(i).neqv.v2%sym(i))) THEN
               CHECKNBAS=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      end function CHECKNBAS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CHECKCOEFS(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks coefficients for good (non-NaN) values

      implicit none
      TYPE (CP), INTENT(IN) :: v
      LOGICAL :: CHECKCOEFS
      INTEGER :: i

      CHECKCOEFS=.TRUE.

      DO i=1,SIZE(v%coef)
         IF (v%coef(i).ne.v%coef(i)) THEN
            CHECKCOEFS=.FALSE.
            EXIT
         ENDIF
      ENDDO

      end function CHECKCOEFS
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE SEPDREPN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
