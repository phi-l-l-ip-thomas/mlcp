!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE GENCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine nextconfig(qns,nbas,maxnmode,nmcfg,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines the next configuration in the sequence:
! 1) Sum-of-indices ranges from low to high, then ...
! 2) number of nonzero indices (nmode) ranges from low to high, then ...
! 3) differences in indices ranges from maximum to minimum, then ...
! 4) the index to migrate ranges from largest to smallest, then ...
! 5) indices are migrated from left to right
! Set nmcfg to .TRUE. to swap the order of conditions 1) and 2)

      implicit none
      integer, intent(in)    :: nbas(:)
      integer, intent(inout) :: qns(:)
      integer, intent(in)    :: maxnmode
      logical, intent(in)    :: nmcfg
      integer, allocatable   :: qt(:),qn(:),lq(:),nb(:),ln(:)
      logical, intent(out)   :: success
      integer :: i,j,k,ndof,nmode,nmmax
      character*60 :: frmt

!     Set parameters
      ndof=SIZE(qns)
      write(frmt,'(A,I0,A)') '(A,',ndof,'(I2,X))'

!     Generated sorted nbas array, 'nb'
      ALLOCATE(qt(ndof),qn(ndof),lq(ndof),nb(ndof),ln(ndof))
      nb=nbas-1
      call sortindlist(nb,ln)

!     Sort qns according to order of descending nb, store in qt
      DO k=1,ndof
         qt(k)=qns(ln(k))-1
      ENDDO
!     The qn array contains qt in descending order
      qn=qt
      call sortindlist(qn,lq)

!     Count the number of coupled DOFs in qn and the number of coupled
!     DOFs allowed by the limits stored in nb
      nmode=0
      nmmax=0
      DO i=1,ndof
         IF (qn(i).gt.0) nmode=nmode+1
         IF (nb(i).gt.0) nmmax=nmmax+1
      ENDDO

!     Error checking
      IF (SIZE(nbas).ne.ndof) &
         call AbortWithError("nextconfig(): qns/nbas size mismatch")
      IF (nmode.gt.maxnmode) &
         call AbortWithError("nextconfig(): nmode > maxnmode")
      IF (nmmax.lt.maxnmode) &
         call AbortWithError("nextconfig(): nmmax < maxnmode")

!     n-mode config version: outermost loop increases nmode
      IF (nmcfg) THEN
         call rearrangeinds(qn,lq,qt,nb,success)
         IF (.not.success) THEN
            qt=qn
            call redistributeinds(qt,nb,success)
            IF (.not.success) THEN
               call increaseorderindsalt(qt,nb,nmode,success)
               IF (.not.success) THEN
                  call increasenmodeinds(qt,nb,maxnmode,nmcfg,success)
               ENDIF
            ENDIF
         ENDIF

!     Regular version: outermost loop increases order
      ELSE
         call rearrangeinds(qn,lq,qt,nb,success)
         IF (.not.success) THEN
            qt=qn
            call redistributeinds(qt,nb,success)
            IF (.not.success) THEN
               call increasenmodeinds(qt,nb,maxnmode,nmcfg,success)
               IF (.not.success) THEN
                  call increaseorderinds(qt,nb,maxnmode,success)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

!     Put qns back into its original order
      DO k=1,ndof
         qns(ln(k))=qt(k)+1
      ENDDO

      DEALLOCATE(qt,qn,lq,nb,ln)

      end subroutine nextconfig

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortindlist(qn,lc)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts qn in descending order of magnitude and then in increasing order
! of position, as stored in the 'lc' array

      implicit none
      integer, intent(inout) :: qn(:)
      integer, allocatable, intent(out) :: lc(:)
      real*8, allocatable :: qnt(:),loc(:)
      integer :: i,ndof,nsame
      real*8  :: val

      ndof=SIZE(qn)

!     Copy the qn values to a temp array and record their locations
      ALLOCATE(loc(ndof),qnt(ndof))
      DO i=1,ndof
         loc(i)=REAL(i)
         qnt(i)=REAL(qn(i))
      ENDDO

!     Sort 'qnt' in descending order of size
      call dsort(qnt,loc,ndof,-2)
!     Since dsort does not sort repeated instances of identical 'qnt'
!     values according to their positions in 'loc', sort them here
      val=nint(qnt(1))
      nsame=1
      DO i=2,ndof
         IF (nint(qnt(i)).eq.val) THEN
            nsame=nsame+1
         ELSE
!           Sort the existing list if it needs it
            IF (nsame.gt.1) &
               call dsort(loc(i-nsame:i-1),qnt(i-nsame:i-1),nsame,2)
!           Reset the nsame counter and value to compare
            val=nint(qnt(i))
            nsame=1
         ENDIF
         IF (i.eq.ndof .and. nsame.gt.1) &
            call dsort(loc(i-nsame+1:i),qnt(i-nsame+1:i),nsame,2)
      ENDDO

!     Copy real arrays from dsort back into integer arrays
      ALLOCATE(lc(ndof))
      DO i=1,ndof
         qn(i)=nint(qnt(i))
         lc(i)=nint(loc(i))
      ENDDO

      DEALLOCATE(loc,qnt)

      end subroutine sortindlist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine rearrangeinds(qn,lc,qns,nb,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tries to find a place in 'qns' after index 'ind' put 'val'

      implicit none
      integer, intent(in)    :: qn(:),nb(:)
      integer, intent(inout) :: lc(:),qns(:)
      logical, intent(out)   :: success
      integer :: i,j,ndof,nmode

!     Set parameters
      ndof=SIZE(qn)
      nmode=0
      DO i=1,ndof
         IF (qn(i).gt.0) nmode=nmode+1
      ENDDO

!     Try to generate the next item in the list by rearranging the order
!     of the indices in the configuration
      success=.FALSE.
      DO i=nmode,1,-1
         call moveindright(qn(i),lc(i),qns,nb,.TRUE.,success)
         IF (success) THEN
            DO j=i+1,nmode
               IF (qn(j).eq.qn(j-1)) THEN
                  lc(j)=lc(j-1)
               ELSE
                  lc(j)=0
               ENDIF
               call moveindright(qn(j),lc(j),qns,nb,.FALSE.,success)
               IF (.not.success) EXIT
            ENDDO
            IF (success) EXIT
         ENDIF
      ENDDO

      end subroutine rearrangeinds

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine moveindright(val,ind,qns,nb,repl,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tries to find a place in 'qns' after index 'ind' put 'val'

      implicit none
      integer, intent(in)    :: nb(:)
      integer, intent(inout) :: qns(:)
      integer, intent(in)    :: val
      integer, intent(inout) :: ind
      logical, intent(in)    :: repl
      logical, intent(out)   :: success
      integer :: j,k,ndof

      ndof=SIZE(qns)

      IF (SIZE(nb).ne.ndof) &
         call AbortWithError("moveright(): qns/nb size mismatch")
      IF (ind.lt.0 .or. ind.gt.ndof) &
         call AbortWithError("moveright(): ind is out of range")

      success=.FALSE.
      DO j=ind+1,ndof
!        val replaces the element in qns
         IF (val.gt.qns(j) .and. val.le.nb(j)) THEN
            qns(j)=val
            IF (ind.gt.0 .and. repl) qns(ind)=0
            DO k=1,ndof
               IF (qns(k).lt.val .or. (qns(k).eq.val .and. k.gt.j)) THEN
                  qns(k)=0
               ENDIF
            ENDDO
            success=.TRUE.
            EXIT
!        Do not let equal values "cross"
         ELSEIF (val.eq.qns(j)) THEN
            EXIT
         ENDIF
      ENDDO

      IF (success) ind=j

      end subroutine moveindright

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine redistributeinds(qn,nb,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Redistributes indices consistent with the number of coupled DOFs
! (nmode) and the index sum of the qn array)

      implicit none
      integer, intent(inout) :: qn(:)
      integer, intent(in)    :: nb(:)
      logical, intent(out)   :: success
      integer :: j,k,l,ndof,kmod,nmode,sumtmp,sumqn

      ndof=SIZE(qn)

!     Make sure the index set is in descending order
      DO j=2,ndof
         IF (qn(j).gt.qn(j-1)) &
            call AbortWithError('redistributeinds(): qn is not ordered')
      ENDDO

!     Determine the number of coupled DOF (nmode) and the index sum
      nmode=0
      sumqn=0
      DO j=1,ndof
         IF (qn(j).gt.0) THEN
            nmode=nmode+1
            sumqn=sumqn+qn(j)
         ENDIF
      ENDDO

      success=.FALSE.
      DO j=nmode,2,-1
!        Look for a qn value that is two more than the j-th
         DO k=j-1,1,-1
            IF (qn(k).gt.qn(j)+1) THEN
               kmod=k
               success=.TRUE.
               EXIT
            ENDIF
         ENDDO
!        Make sure the remaining quanta can fit in the slots after kmod
!        since the kmod-th slot is decremented by 1
         IF (success) THEN
            sumtmp=sumqn+1
            DO k=1,kmod
               sumtmp=sumtmp-qn(k)
            ENDDO
            l=qn(kmod)-1
            DO k=kmod+1,nmode
               l=min(sumtmp-nmode+k,l,nb(k))
               sumtmp=sumtmp-l
            ENDDO
            IF (sumtmp.gt.0) success=.FALSE.
         ENDIF
         IF (success) EXIT
      ENDDO

!     Redistribute the quanta
      IF (success) THEN
         sumtmp=sumqn+1
         DO k=1,kmod
            sumtmp=sumtmp-qn(k)
         ENDDO
         qn(kmod)=qn(kmod)-1
         DO k=kmod+1,nmode
            qn(k)=min(sumtmp-nmode+k,qn(k-1),nb(k))
            sumtmp=sumtmp-qn(k)
         ENDDO
      ENDIF

      end subroutine redistributeinds

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine increasenmodeinds(qn,nb,mxnmode,nmcfg,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Redistributes indices consistent with the number of coupled DOFs
! (nmode) and the index sum of the qn array)

      implicit none
      integer, intent(inout) :: qn(:)
      integer, intent(in)    :: nb(:)
      integer, intent(in)    :: mxnmode
      logical, intent(in)    :: nmcfg
      logical, intent(out)   :: success
      integer :: j,ndof,nmode,sumqn

      ndof=SIZE(qn)
      success=.FALSE.

!     Determine the number of coupled DOF (nmode) and the index sum
      nmode=0
      sumqn=0
      DO j=1,ndof
         IF (qn(j).gt.0) THEN
            nmode=nmode+1
         ENDIF
         sumqn=sumqn+qn(j)
      ENDDO

!     Conditions for failure
      IF (nmode.eq.ndof .or. nmode.ge.mxnmode) THEN
         IF (nmcfg) qn=0
         RETURN
      ENDIF
      IF (.not.nmcfg .and. nmode.ge.sumqn) RETURN

!     Increment nmode and generate the new configuration
      nmode=nmode+1
      qn=0
      qn(1:nmode)=1
      sumqn=sumqn-nmode
!     Alternative version: generate the lowest-order config consistent
!     with the new value of nmode
      IF (nmcfg) sumqn=0
      DO j=1,nmode
         qn(j)=qn(j)+min(sumqn,nb(j)-1)
         sumqn=sumqn-min(sumqn,nb(j)-1)
         IF (sumqn.eq.0) EXIT
      ENDDO
      success=.TRUE.

      end subroutine increasenmodeinds

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine increaseorderinds(qn,nb,mxnmode,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Increases the total order, that is, the sum of the indices in qn, and
! generates the first configuration consistent with the new total order,
! constrained by the qn limits in 'nb' and the maximum nmode order.
! Upon failure the qn list is reset to its minimum values.
! The nb array should be zero-based and in descending order on entry.

      implicit none
      integer, intent(inout) :: qn(:)
      integer, intent(in)    :: nb(:)
      integer, intent(in)    :: mxnmode
      logical, intent(out)   :: success
      integer :: j,ndof,sumqn,sumnb

      ndof=SIZE(qn)
      success=.TRUE.

      IF (SIZE(nb).ne.ndof) call &
         AbortWithError('increaseorderinds(): qn/nb size mismatch')
      IF (mxnmode.gt.ndof) call &
         AbortWithError('increaseorderinds(): mxnmode > ndof')

!     Determine the qn and nb sums
      sumqn=0
      sumnb=0
      DO j=1,ndof
         sumqn=sumqn+qn(j)
         sumnb=sumnb+nb(j)
      ENDDO

!     Error if the qn sum is greater than nb (should never happen)
      IF (sumqn.gt.sumnb) &
         call AbortWithError('increaseorderinds(): sumqn > sumnb')

!     Fill each entry in qn with as many quanta as allowed by nb
      qn=0
      sumqn=sumqn+1
      DO j=1,mxnmode
         qn(j)=min(sumqn,nb(j))
         sumqn=sumqn-min(sumqn,nb(j))
         IF (sumqn.eq.0) EXIT
      ENDDO

!     If not all quanta are placed, reset to the minimum qn values
      IF (sumqn.gt.0) THEN
         success=.FALSE.
         qn(1:mxnmode)=0
      ENDIF

      end subroutine increaseorderinds

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine increaseorderindsalt(qn,nb,mxnmode,success)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Increases the total order, that is, the sum of the indices in qn, and
! generates the first configuration consistent with the new total order,
! constrained by the qn limits in 'nb' and the current nmode order.
! Upon failure the original qn list is retained.
! The nb array should be zero-based and in descending order on entry.

      implicit none
      integer, intent(inout) :: qn(:)
      integer, intent(in)    :: nb(:)
      integer, intent(in)    :: mxnmode
      logical, intent(out)   :: success
      integer, allocatable   :: qnstor(:)
      integer :: j,ndof,sumqn,sumnb,nmode

      ndof=SIZE(qn)
      success=.TRUE.

      IF (SIZE(nb).ne.ndof) call &
         AbortWithError('increaseorderinds(): qn/nb size mismatch')
      IF (mxnmode.gt.ndof) call &
         AbortWithError('increaseorderinds(): mxnmode > ndof')

!     Since the qn array is overwritten, store the indices so that they
!     can be restored if the subroutine does not return success
      ALLOCATE(qnstor(ndof))
      qnstor=qn

!     Determine nmode and the qn and nb sums
      sumqn=0
      sumnb=0
      nmode=0
      DO j=1,ndof
         IF (qn(j).gt.0) THEN
            nmode=nmode+1
         ENDIF
         sumqn=sumqn+qn(j)
         sumnb=sumnb+nb(j)
      ENDDO

!     Error if the qn sum is greater than nb (should never happen)
      IF (sumqn.gt.sumnb) &
         call AbortWithError('increaseorderinds(): sumqn > sumnb')

!     Fill each entry in qn with as many quanta as allowed by nb,
!     retaining the same value of nmode
      qn=0
      qn(1:nmode)=1
      sumqn=sumqn+1-nmode
      DO j=1,mxnmode
         qn(j)=qn(j)+min(sumqn,nb(j)-1)
         sumqn=sumqn-min(sumqn,nb(j)-1)
         IF (sumqn.eq.0) EXIT
      ENDDO

!     If not all quanta are placed, return to the previous configuration
      IF (sumqn.gt.0) THEN
         success=.FALSE.
         qn=qnstor
      ENDIF

      DEALLOCATE(qnstor)

      end subroutine increaseorderindsalt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE GENCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
