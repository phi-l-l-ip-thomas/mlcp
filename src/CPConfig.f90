!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module manages "configuration" states and their conversion to
! and from CP-format

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN

      TYPE Configs
         INTEGER, ALLOCATABLE :: nbas(:),qns(:,:)
         REAL*8,  ALLOCATABLE :: coef(:)
      END TYPE Configs

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE NewConfigs(v,nbas,rnk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Makes a new configuration list. Entries are initialized to zero

      IMPLICIT NONE
      TYPE (Configs), INTENT(OUT) :: v
      INTEGER, INTENT(IN) :: nbas(:)
      INTEGER, INTENT(IN) :: rnk
      INTEGER :: i,ndof

      ndof=SIZE(nbas)

!     Error checking
      IF (rnk.lt.1) THEN
         write(*,*) 'rnk = ',rnk
         call AbortWithError('NewConfigs(): rank < 1')
      ENDIF
      DO i=1,ndof
         IF (nbas(i).lt.1) THEN
            write(*,*) 'i = ',i,'; nbas(i) = ',nbas(i)
            call AbortWithError('NewConfigs(): invalid basis size')
         ENDIF
      ENDDO

      ndof=SIZE(nbas)

      ALLOCATE(v%nbas(ndof),v%qns(rnk,ndof),v%coef(rnk))
      v%nbas=nbas
      v%coef=0.d0
      v%qns=0

      END SUBROUTINE NewConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintConfigs(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Print out list of configurations

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      integer :: i,j,nconf,ndof
      character*64 :: frmt

      IF (.not.ALLOCATED(v%qns)) THEN
          write(*,*) '{ }'
          return
      ENDIF

      nconf=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

      DO i=1,nconf
         write(frmt,'(A,I0,A)') '(X,',ndof,'(I3,X),f26.12)'
         write(*,frmt) (v%qns(i,j),j=1,ndof),v%coef(i)
      ENDDO

      end subroutine PrintConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE FlushConfigs(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes configuration list

      IMPLICIT NONE
      TYPE (Configs) :: v

      IF (ALLOCATED(v%nbas)) DEALLOCATE(v%nbas)
      IF (ALLOCATED(v%qns)) DEALLOCATE(v%qns)
      IF (ALLOCATED(v%coef)) DEALLOCATE(v%coef)

      END SUBROUTINE FlushConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateConfigs(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks configurations for input errors

      IMPLICIT NONE
      TYPE (Configs), INTENT(IN) :: v
      INTEGER :: i,j,ndof,nrk

      ndof=SIZE(v%nbas)
      nrk=SIZE(v%coef)

!     Make sure array sizes match
      IF (SIZE(v%qns,1).ne.nrk) THEN
         write(*,*) 'nrk: from coef = ',nrk,&
                    ', from qns = ',SIZE(v%qns,1)
         call AbortWithError('ValidateConfigs(): inconsistent nrk')
      ENDIF
      IF (SIZE(v%qns,2).ne.ndof) THEN
         write(*,*) 'ndof: from nbas = ',ndof,&
                    ', from qns = ',SIZE(v%qns,2)
         call AbortWithError('ValidateConfigs(): inconsistent ndof')
      ENDIF

!     Make sure basis sizes are positive integers
      DO i=1,ndof
         IF (v%nbas(i).lt.1) THEN
            write(*,*) 'i = ',i,'; nbas(i) = ',v%nbas(i)
            call AbortWithError('ValidateConfigs(): invalid basis size')
         ENDIF
      ENDDO

!     Make sure quantum numbers are in range allowed by nbas
      DO i=1,nrk
         DO j=1,ndof
!           Give a warning if the qns value is not yet set
            IF (v%qns(i,j).eq.0) THEN
               write(*,*) 'i,j = ',i,j,'; qns(i,j) = ',v%qns(i,j)
               call ShowWarning('ValidateConfigs(): unset qn')
            ENDIF
!           Give a warning if the qns value is negative but in range
            IF (v%qns(i,j).lt.0 .and. abs(v%qns(i,j)).le.v%nbas(j)) THEN
               write(*,*) 'i,j = ',i,j,'; qns(i,j) = ',v%qns(i,j)
               call ShowWarning('ValidateConfigs(): negative qn value')
            ENDIF
!           Error out for out-of-range qns value
            IF (abs(v%qns(i,j)).gt.v%nbas(j)) THEN
               write(*,*) 'i,j = ',i,j,'; qns(i,j) = ',v%qns(i,j)
               call AbortWithError('ValidateConfigs(): qn out of range')
            ENDIF
         ENDDO
      ENDDO

      end subroutine ValidateConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ReplaceConfigsVwithW(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, disposing W afterwards

      IMPLICIT NONE
      TYPE (Configs) :: v,w

      call FlushConfigs(v)
      call CopyConfigsWtoV(v,w)
      call FlushConfigs(w)

      END SUBROUTINE ReplaceConfigsVwithW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CopyConfigsWtoV(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies W into V, leaving W intact.

      IMPLICIT NONE
      TYPE (Configs) :: v,w
      INTEGER :: nrk

      nrk=SIZE(w%coef)
      call NewConfigs(v,w%nbas,nrk)
      call GenCopyConfigsWtoV(v,w,1,nrk,1,nrk)

      END SUBROUTINE CopyConfigsWtoV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE GenCopyConfigsWtoV(v,w,vi,ve,wi,we)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies a group of consecutive terms in W to consecutive slots in V, 
! leaving W intact. v must be allocated beforehand.

      IMPLICIT NONE
      TYPE (Configs) :: v,w
      INTEGER, INTENT(in) :: vi,ve,wi,we
      INTEGER :: rkv,rkw

      rkv=SIZE(v%coef)
      rkw=SIZE(w%coef)

      IF (vi.lt.1 .or. ve.gt.rkv .or. vi.gt.ve .or. &
          wi.lt.1 .or. we.gt.rkw .or. wi.gt.we .or. &
          we-wi.ne.ve-vi) THEN
          write(*,'(2A,6(X,I0))') 'Bad rank indices: ',&
          'vi,ve,rkv,wi,we,rkw =',vi,ve,rkv,wi,we,rkw
          CALL AbortWithError('Error in GenCopyConfigsWtoV()')
      ENDIF

      v%qns(vi:ve,:)=w%qns(wi:we,:)
      v%coef(vi:ve)=w%coef(wi:we)

      END SUBROUTINE GenCopyConfigsWtoV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE ResizeConfigList(v,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resizes vector v, either by truncating at a smaller rank or by adding
! space for extra terms. If the size is set to zero, a zero vector of
! rank-1 is generated

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      TYPE (Configs) :: w
      INTEGER, INTENT(IN) :: nrk
      INTEGER :: rkv

      rkv=SIZE(v%coef)

      IF (nrk.lt.0) &
         call AbortWithError('Error in ResizeConfigList(): nrk < 0')

      call ReplaceConfigsVwithW(w,v)
      call NewConfigs(v,w%nbas,MAX(1,nrk))
      IF (nrk.gt.0) &
         call GenCopyConfigsWtoV(v,w,1,MIN(rkv,nrk),1,MIN(rkv,nrk))
      call FlushConfigs(w)

      END SUBROUTINE ResizeConfigList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE TrimZeroConfigs(v,tol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resize vector v by trimming columns with coefficients with absolute
! values smaller than tol

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      TYPE (Configs) :: w
      real*8, intent(in) :: tol
      integer, allocatable :: iok(:)
      integer :: i,nok,rk

      rk=SIZE(v%coef)
      ALLOCATE(iok(rk))

!     Check the cols for nonzero coef
      nok=0
      DO i=1,rk
         IF (abs(v%coef(i)).le.abs(tol)) CYCLE
         nok=nok+1
         iok(nok)=i
      ENDDO

      IF (nok.lt.rk) THEN
         call NewConfigs(w,v%nbas,max(1,nok))
         DO i=1,nok
            call GenCopyConfigsWtoV(w,v,i,i,iok(i),iok(i))
         ENDDO
         call ReplaceConfigsVwithW(v,w)
      ENDIF

      DEALLOCATE(iok)

      END SUBROUTINE TrimZeroConfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE TrimZeroQns(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Resize vector v by trimming columns with containing a zero qn

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      TYPE (Configs) :: w
      integer, allocatable :: iok(:)
      integer :: i,nok,rk

      rk=SIZE(v%coef)
      ALLOCATE(iok(rk))

!     Check the cols for nonzero qn
      nok=0
      DO i=1,rk
         IF (ANY(v%qns(i,:).eq.0)) CYCLE
         nok=nok+1
         iok(nok)=i
      ENDDO

      IF (nok.lt.rk) THEN
         call NewConfigs(w,v%nbas,max(1,nok))
         DO i=1,nok
            call GenCopyConfigsWtoV(w,v,i,i,iok(i),iok(i))
         ENDDO
         call ReplaceConfigsVwithW(v,w)
      ENDIF

      DEALLOCATE(iok)

      END SUBROUTINE TrimZeroQns

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSignOnCoef(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks config list and changes signs to make all qns positive

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      INTEGER :: i,j,ndof,nrk

      ndof=SIZE(v%nbas)
      nrk=SIZE(v%coef)

      DO i=1,nrk
         DO j=1,ndof
            IF (v%qns(i,j).lt.0) THEN
               v%coef(i)=-v%coef(i)
               v%qns(i,j)=-v%qns(i,j)
            ENDIF
         ENDDO
      ENDDO

      end subroutine PutSignOnCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PutSignOnQn(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks config list and changes signs to make all coefficients positive

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      INTEGER :: i,j,ndof,nrk
      real*8, parameter :: tol=1.d-12

      ndof=SIZE(v%nbas)
      nrk=SIZE(v%coef)

!     Change sign on 1st dof if coef is negative
      DO i=1,nrk
         IF (v%coef(i).lt.0) THEN
            v%coef(i)=-v%coef(i)
            v%qns(i,1)=-v%qns(i,1)
         ENDIF
      ENDDO

      end subroutine PutSignOnQn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Config2CP(F,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Constructs CP-vector from list of configurations

      implicit none
      TYPE (CPvec), INTENT(OUT)  :: F
      TYPE (Configs), INTENT(IN) :: v
      integer :: ndof,i,j,rF,gst

!     Set parameters
      ndof=SIZE(v%nbas)
      rF=SIZE(v%coef)

!     If v is a zero config, get a zero CP-vec
      IF (rF.eq.1 .and. ALL(v%qns(1,:).eq.0)) THEN
         call GetZeroCPvec(F,v%nbas)
         RETURN
      ENDIF

!     Reconstruct F
      call NewCPvec(F,v%nbas,rF)
      F%base=1.d-15
      DO i=1,rF
         F%coef(i)=abs(v%coef(i))
         gst=0
         DO j=1,ndof
            IF (v%qns(i,j).ne.0) &
               F%base(gst+abs(v%qns(i,j)),i)=SIGN(1.d0,REAL(v%qns(i,j)))
            gst=gst+v%nbas(j)
         ENDDO
!        If coef is negative, change sign of base for 1st DOF
         IF (v%coef(i).lt.0.d0) &
            F%base(abs(v%qns(i,1)),i)=-F%base(abs(v%qns(i,1)),i)
      ENDDO

      end subroutine Config2CP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CP2Config(F,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Converts rank-1 term in CP-format rep'n of quadratic/cubic/quartic
! tensor to list of indices.

      implicit none
      TYPE (CPvec), INTENT(IN)    :: F
      TYPE (Configs), INTENT(OUT) :: v
      integer :: i,j,k,gst,ndof,rF
      logical :: found
      real*8, parameter :: tol=1.d-12

!     Parameters
      ndof=SIZE(F%nbas)
      rF=SIZE(F%coef)

      call NewConfigs(v,F%nbas,rF)

!     Convert a zero CP-vector to a zero config
      IF (rF.eq.1 .and. F%coef(1).eq.0.d0) RETURN

!     Look for the excited component for each rank and DOF
      DO i=1,rF
         v%coef(i)=F%coef(i)
         gst=0
         DO j=1,ndof
            found=.FALSE.
            DO k=1,F%nbas(j)
!              Base has +1 or -1...
               IF (abs(abs(F%base(gst+k,i))-1.d0).lt.tol) THEN

!                 Error if more than one excited component
                  IF (found) THEN
                     write(*,*) 'F has more than one excited component'
                     call AbortWithError("CP2Config(): bad F")
                  ENDIF
                  v%qns(i,j)=k
                  found=.TRUE.
!                 If term is negative, put sign on coefficient
                  IF (F%base(gst+k,i).lt.0.d0) v%coef(i)=-v%coef(i)

!              Error if base has something other than +1, -1, or zero...
               ELSEIF (abs(F%base(gst+k,i)).gt.tol) THEN
                  write(*,*) 'F is unnormalized or is not a unit vector'
                  call AbortWithError("CP2Config(): bad F")
               ENDIF
            ENDDO

!           Error if no excited components are found
            IF (.not.found) THEN
               write(*,*) 'F has no excited components'
               call AbortWithError("CP2Config(): bad F")
            ENDIF
            gst=gst+F%nbas(j)
         ENDDO
      ENDDO

      end subroutine CP2Config

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ConfigSetOverlap(B,C,F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the overlap of each of a set of configurations in B with
! CP-vector F, returning the overlaps in C

      implicit none
      TYPE (Configs), INTENT(IN)  :: B
      TYPE (Configs), INTENT(OUT) :: C
      TYPE (CPvec), INTENT(IN)    :: F
      integer :: i,rk

      rk=SIZE(B%coef)

!     Copy the configuration list B -> C
      call CopyConfigsWtoV(C,B)

!     Get the new coefficients of C
!$omp parallel
!$omp do private(i)
      DO i=1,rk
         C%coef(i)=PRODVConfig(F,B%qns(i,:))
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine ConfigSetOverlap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODVConfig(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates <F,G>, with G as a single configuration
! Here G is just a list of indices, not an object of TYPE "Configs"

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      integer, intent(in) :: G(:)
      real*8  :: PRODVConfig,prod
      integer :: i,j,ndof,rF,gst

      ndof=SIZE(G)
      rF=SIZE(F%coef)

!     Error checking
      IF (SIZE(F%nbas).ne.ndof) &
         call AbortWithError('PRODVConfig(): ndof_F != ndof_G')
      DO i=1,ndof
         IF (G(i).lt.0 .or. G(i).gt.F%nbas(i)) &
         call AbortWithError('PRODVConfig(): G index is out of range')
      ENDDO

!     Compute the overlap <F,G>
      PRODVConfig=0.d0
      DO i=1,rF
         gst=0
         prod=1.d0
         DO j=1,ndof
            prod=prod*F%base(gst+G(j),i)
            gst=gst+F%nbas(j)
         ENDDO
         PRODVConfig=PRODVConfig+F%coef(i)*prod
      ENDDO

      end function PRODVConfig

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine collectconfigs(v,avg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Combines and sums equivalent terms of a symmetric tensor (in
! configuration format), e.g. 1 1 2 , 1 2 1, and 2 1 1.

      implicit none
      TYPE (Configs), INTENT(INOUT) :: v
      TYPE (Configs) :: w
      logical, intent(in)  :: avg
      integer, allocatable :: tmpq(:),dupct(:)
      integer :: i,j,k,l,jn,ndof,nrk,tmp
      logical :: found,match
      real*8, parameter :: tol=1.d-12

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     Auxiliary arrays
      ALLOCATE(tmpq(ndof),dupct(nrk))
      dupct=1

      call NewConfigs(w,v%nbas,nrk)

!     Error checking
      tmp=v%nbas(1)
      DO i=2,ndof
         IF (v%nbas(i).ne.tmp) THEN
            write(*,*) 'nbas(1) = ',tmp,'; i,nbas(i) = ',i,v%nbas(i)
            call AbortWithError('collectconfigs(): v is not symmetric')
         ENDIF
      ENDDO

!     Copy first term in v to w
      call GenCopyConfigsWtoV(w,v,1,1,1,1)
      jn=1

!     Compare each term in original list to term in temp list
      DO i=2,nrk
         found=.FALSE.

         DO j=1,jn

            tmpq(:)=w%qns(j,:)
            DO k=1,ndof

               match=.FALSE.
               DO l=k,ndof
                  IF (v%qns(i,k).eq.tmpq(l)) THEN
                     tmp=tmpq(l)
                     tmpq(l)=tmpq(k)
                     tmpq(k)=tmp
                     match=.TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (.not.match) EXIT
               IF (k.eq.ndof) found=.TRUE.
            ENDDO

!           Term is found: add coef to existing term. Also keep track of
!           how many times term occurs for averaging
            IF (found) THEN
               w%coef(j)=w%coef(j)+v%coef(i)
               dupct(j)=dupct(j)+1
               EXIT
            ENDIF
         ENDDO
!        If term is not found, add to temp list
         IF (.not.found) THEN
            jn=jn+1
            w%qns(jn,:)=v%qns(i,:)
            w%coef(jn)=v%coef(i)
         ENDIF
      ENDDO

!     Average the coefficient, if requested
      IF (avg) THEN
         DO j=1,jn
            w%coef(j)=w%coef(j)/dupct(j)
         ENDDO
      ENDIF

      call TrimZeroConfigs(w,tol)
      call ReplaceConfigsVwithW(v,w)

      DEALLOCATE(tmpq,dupct)

      end subroutine collectconfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine removeduplicateconfigs(v,typ)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts a configuration list by index and removes duplicate configs:
! typ = 0: duplicates are removed
! typ = 1: duplicates are summed
! typ =-1: duplicates are averaged

      implicit none
      TYPE (Configs), INTENT(INOUT) :: v
      real*8, allocatable :: navg(:)
      integer, intent(in) :: typ
      integer :: i,j,k,nrk,ndof
      logical :: difft

!     Error checking
      IF (abs(typ).gt.1) call &
         AbortWithError('removeduplicateconfigs(): bad "typ" value')

      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     For averaging configs: navg keeps number of occurences
      IF (typ.eq.-1) THEN
         ALLOCATE(navg(nrk))
         navg=1.d0
      ENDIF

!     Sort the list. Once the list of qns is sorted, one only needs to 
!     compare configurations until a different index is encountered
      call SortConfigsByIndex(v)

      DO i=1,nrk-1
!        Skip configurations in i that have been flagged
         IF (v%qns(i,1).eq.0) CYCLE
!        Compare the configuration with those that follow
         DO j=i+1,nrk
!           Skip configurations in j that have been flagged
            IF (v%qns(j,1).eq.0) CYCLE
!           Exit j loop as soon as a different config is found
            difft=.FALSE.
            DO k=ndof,1,-1
               IF (v%qns(i,k).ne.v%qns(j,k)) THEN
                  difft=.TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (difft) EXIT
!           Config j is the same as config i.
            IF (abs(typ).eq.1) THEN
               v%coef(i)=v%coef(i)+v%coef(j)          ! Sum coefs
               IF (typ.eq.-1) navg(i)=navg(i)+navg(j) ! # occurences
            ENDIF
!           Flag config for removal
            v%qns(j,1)=0
         ENDDO
      ENDDO

!     Averaging: divide each coefficient by the number of occurences
      IF (typ.eq.-1) THEN
         DO i=1,nrk
            v%coef(i)=v%coef(i)/navg(i)
         ENDDO
         DEALLOCATE(navg)
      ENDIF

!     Remove all configurations having a zero in v%qns
      call TrimZeroQns(v)

      end subroutine removeduplicateconfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine randomizeconfiglist(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Puts the R terms in a configuration list in random order

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      TYPE (Configs) :: w
      REAL*8, ALLOCATABLE :: tabindx(:),rndnr(:)
      INTEGER :: nrk,i

      nrk=SIZE(v%qns,1)

!     Generate a list of random numbers and an index set
      ALLOCATE(tabindx(nrk),rndnr(nrk))
      DO i=1,nrk
         tabindx(i)=REAL(i)
      ENDDO
      call random_number(rndnr(:))

!     Sort the random array
      call dsort(rndnr,tabindx,nrk,-2)

!     Resort v according to the random index set
      call NewConfigs(w,v%nbas,nrk)
      DO i=1,nrk
         w%qns(i,:)=v%qns(int(tabindx(i)),:)
         w%coef(i)=v%coef(int(tabindx(i)))
      ENDDO

      DEALLOCATE(tabindx,rndnr)
      call ReplaceConfigsVwithW(v,w)

      end subroutine randomizeconfiglist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine swapconfigs(v,i,j)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts configurations in descending order of coefficients

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      INTEGER, INTENT(IN)  :: i,j
      INTEGER, ALLOCATABLE :: qnt(:)
      REAL*8  :: ct
      INTEGER :: ndof,nrk

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     Error checking
      IF (i.lt.0 .or. j.lt.0 .or. i.gt.nrk .or. j.gt.nrk) &
         call AbortWithError('swaponfigs(): i or j out of range')

!     Quick exit if i.eq.j
      IF (i.eq.j) RETURN

!     Swap terms i and j
      ALLOCATE(qnt(ndof))
      qnt=v%qns(i,:)
      ct=v%coef(i)
      v%qns(i,:)=v%qns(j,:)
      v%coef(i)=v%coef(j)
      v%qns(j,:)=qnt
      v%coef(j)=ct
      DEALLOCATE(qnt)

      end subroutine swapconfigs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortConfigsByCoef(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts configurations in descending order of coefficients

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      INTEGER, ALLOCATABLE :: qnt(:,:)
      REAL*8, ALLOCATABLE  :: ord(:)
      INTEGER :: i,ndof,nrk

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     Quick exit for rank-1 list
      IF (nrk.eq.1) THEN
         call PutSignOnCoef(v)
         RETURN
      ENDIF

!     Make sure all coefs are positive
      call PutSignOnQn(v)

!     Copy the quantum number list and generate the order index list
      ALLOCATE(ord(nrk),qnt(nrk,ndof))
      qnt=v%qns
      DO i=1,nrk
         ord(i)=REAL(i)
      ENDDO

!     Sort by coefficients
      call dsort(v%coef,ord,nrk,-2)

      DO i=1,nrk
         v%qns(i,:)=qnt(nint(ord(i)),:)
      ENDDO

!     Make sure qns are positive before exiting
      call PutSignOnCoef(v)

      DEALLOCATE(ord,qnt)

      end subroutine SortConfigsByCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortConfigsByIndex(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts configurations in ascending order of indices, with the first
! index changing the slowest

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      integer, allocatable :: ist(:),iend(:)
      integer :: nrk,ndof,j

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     Make sure qns are positive before sorting
      call PutSignOnCoef(v)

!     Quick exit for rank-1 list
      IF (nrk.eq.1) RETURN

      ALLOCATE(ist(ndof),iend(ndof))
      iend=0
      iend(1)=nrk
      ist(1)=1
      j=1
      DO
!        Sort the configurations by DOF index j
         call SortConfigBlock(v,j,ist(j),iend(j))

!        Update DOF index j
         IF (j.lt.ndof) j=j+1
         DO 
            IF (iend(j).lt.iend(j-1)) EXIT
            j=j-1
            IF (j.eq.1) EXIT
         ENDDO

!        When j returns to 1, the entire vector is sorted
         IF (j.eq.1) EXIT

!        Update the sort ranges
         ist(j)=iend(j)+1
         iend(j)=ibisect(v%qns(ist(j):iend(j-1),j-1),1)+ist(j)-1
      ENDDO

      DEALLOCATE(ist,iend)

      end subroutine SortConfigsByIndex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SortConfigBlock(v,i,ri,rf)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts configurations in v in ascending order of index i
! Negative signs must be on coefficients (not on qns) on entry

      IMPLICIT NONE
      TYPE (Configs), INTENT(INOUT) :: v
      integer, intent(in)  :: i,ri,rf
      integer, allocatable :: qnt(:,:)
      real*8, allocatable  :: ct(:),qit(:),ord(:)
      integer :: j,nrk,ndof,rks

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)
      rks=rf-ri+1

!     Error checking
      IF (i.lt.1 .or. i.gt.ndof) &
         CALL AbortWithError("SortConfigBlock(): i out of range")
      IF (ri.lt.1 .or. rf.gt.nrk .or. ri.gt.rf) &
         CALL AbortWithError("SortConfigBlock(): ri/rf out of range")

!     No sorting needed if ri.eq.rf
      IF (rks.eq.1) RETURN

!     If ri and rf differ by 1, just swap if needed
      IF (rks.eq.2) THEN
         IF (v%qns(ri,i).gt.v%qns(rf,i)) call swapconfigs(v,ri,rf)
         RETURN
      ENDIF

!     Copy the quantum number list and coefs, 
!     and generate the real qn and order index lists
      ALLOCATE(qnt(rks,ndof),ct(rks),qit(rks),ord(rks))
      qnt(1:rks,:)=v%qns(ri:rf,:)
      ct(1:rks)=v%coef(ri:rf)
      qit(1:rks)=REAL(v%qns(ri:rf,i))
      DO j=1,rks
         ord(j)=REAL(j)
      ENDDO

!     Sort by indices
      call dsort(qit,ord,rks,2)

      DO j=1,rks
         v%qns(ri+j-1,:)=qnt(nint(ord(j)),:)
         v%coef(ri+j-1)=ct(nint(ord(j)))
      ENDDO

      DEALLOCATE(qnt,ct,qit,ord)

      end subroutine SortConfigBlock

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function findconfigindex(v,c,rset)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds position of configuration c in v, restricted to index range
! rset(1):rset(2). Returns zero if c is not found.

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      integer, intent(in) :: c(:)
      integer, intent(in) :: rset(2)
      integer :: jset(2),tset(2)
      integer :: i,ndof,nrk,findconfigindex

!     Set parameters
      nrk=SIZE(v%qns,1)
      ndof=SIZE(v%qns,2)

!     Error checking
      IF (rset(1).lt.1 .or. rset(2).gt.nrk .or. rset(1).gt.rset(2)) &
         CALL AbortWithError("findconfigindex(): rset out of range")
      IF (SIZE(c).ne.ndof) &
         CALL AbortWithError("findconfigindex(): c,v ndof mismatch")

      tset=rset
      DO i=1,ndof
         call ifind(v%qns(tset(1):tset(2),i),c(i),jset)
         IF (jset(1).eq.0) THEN
            findconfigindex=0
            RETURN
         ENDIF
         tset(1)=tset(1)+jset(1)-1
         tset(2)=tset(1)+jset(2)-jset(1)
      ENDDO
      findconfigindex=tset(1)

      end function findconfigindex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function zeronorm(v,tol)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||v||_0, the l-0 norm of v (number of nonzero entries)
! tol is the threshold for being "non-zero"

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      real*8, intent(in) :: tol
      integer :: i,rk,zeronorm

      rk=SIZE(v%coef)

      zeronorm=0
      DO i=1,rk
         IF (abs(v%coef(i)).gt.tol) zeronorm=zeronorm+1
      ENDDO

      end function zeronorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function onenorm(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||v||_1, the l-1 norm of v

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      real*8  :: onenorm
      integer :: i,rk

      rk=SIZE(v%coef)

      onenorm=0.d0
      DO i=1,rk
         onenorm=onenorm+abs(v%coef(i))
      ENDDO

      end function onenorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function twonorm(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||v||_2, the l-2 norm of v

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      real*8  :: twonorm
      integer :: i,rk

      rk=SIZE(v%coef)

      twonorm=0.d0
      DO i=1,rk
         twonorm=twonorm+v%coef(i)**2
      ENDDO

      twonorm=sqrt(twonorm)

      end function twonorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function infnorm(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes ||v||_infinity, the l-infinity norm of v

      implicit none
      TYPE (Configs), INTENT(IN) :: v
      real*8  :: infnorm
      integer :: i,rk

      rk=SIZE(v%coef)

      infnorm=0.d0
      DO i=1,rk
         IF (abs(v%coef(i)).gt.infnorm) infnorm=abs(v%coef(i))
      ENDDO

      end function infnorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODCC(v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the dot product of configuration vectors v*w

      implicit none
      TYPE (Configs), INTENT(IN) :: v,w
      real*8  :: PRODCC
      integer :: i,j,k,rkv,rkw

      IF (SIZE(v%nbas).ne.SIZE(w%nbas)) &
         call AbortWithError("PRODCC(): v and w have different # DOFs")

      rkv=SIZE(v%coef)
      rkw=SIZE(w%coef)

      PRODCC=0.d0
      DO i=1,rkv
         DO j=1,rkw
            IF (ALL(v%qns(i,:).eq.w%qns(j,:))) &
               PRODCC=PRODCC+v%coef(i)*w%coef(j)
         ENDDO
      ENDDO

      end function PRODCC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
