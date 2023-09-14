!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Manages the multi-layer wavefunction

      USE ERRORTRAP
      USE UTILS

      implicit none
      TYPE MLtree
          integer, dimension(:,:), allocatable :: modcomb,modstart,&
          whichmod,gdim
          integer, dimension(:), allocatable :: nmode,resort
          integer :: nlayr,ndof
      END TYPE MLtree

      real*8  :: mc_time
      logical :: MC_SETUP = .FALSE.

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_ModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initializes the module

      implicit none

      mc_time=0.d0
      MC_SETUP=.TRUE.

      end subroutine Init_ModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadModeDat(ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads input file ('layers.inp') containing mode combination data and
! fills type 'MLtree'

      implicit none
      TYPE (MLtree) :: ML
      character(len=64), intent(in) :: fnm      
      integer, allocatable :: nmode_tmp(:),blankline(:),res_tmp(:)
      integer, allocatable :: bas_tmp(:,:),layrs_tmp(:,:)
      integer :: il,im,ii,i2,u,InpStat,ReadStat,linelen,maxndof
      integer :: iri,irf,ili,ilf,ibi,ibf,maxlines,ilb,ill
      integer :: iline,nlines,iblk,nblk,rblk,bblk,lblk
      character(len=1024)  :: line
      character(len=32)    :: string
      real*8  :: t1,t2

      IF (.not. MC_SETUP) CALL Init_ModeComb

      CALL CPU_TIME(t1)

!     For line reading (change if working with a larger system)
      maxndof=1024
      maxlines=1024
      iri=0
      irf=0
      ibi=0
      ibf=0
      ili=0
      ilf=0

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         call AbortWithError("Oh, no! Error reading input file")
      ENDIF

      ALLOCATE(blankline(maxlines))
      iline=0
      iblk=0
      DO
!        Read line
         READ(u,"(A1024)",IOSTAT=ReadStat) Line
         IF ( ReadStat /= 0 ) EXIT
         iline=iline+1
!        remove leading spaces...
         line = ADJUSTL(line)
!        get length of line without the spaces at the end...
         linelen = LEN_TRIM(line)
!        keep track of empty or comment lines
         IF (linelen.eq.0 .or. line(1:1).eq.'#') THEN
            iblk=iblk+1
            blankline(iblk)=iline
!        determine where sections begin and end
         ELSEIF (line(1:7) == '$resort') THEN
            iri=iline
         ELSEIF (line(1:11) == '$end-resort') THEN
            irf=iline
         ELSEIF (line(1:6) == '$basis') THEN 
            ibi=iline
         ELSEIF (line(1:10) == '$end-basis') THEN 
            ibf=iline
         ELSEIF (line(1:7) == '$layers') THEN
            ili=iline
         ELSEIF (line(1:11) == '$end-layers') THEN
            ilf=iline
         ELSEIF (line(linelen:linelen).ne.'/') THEN
            CALL AbortWithError('Error: input list must end with "/"')
         ENDIF
      ENDDO
      nblk=iblk
      nlines=iline

!     Detect input errors
      IF (iri.eq.0 .or. irf.eq.0 .or. ibi.eq.0 .or. ibf.eq.0 .or. &
          ili.eq.0 .or. ilf.eq.0 .or. iri.gt.irf .or. &
          ibi.gt.ibf .or. ili.gt.ilf .or. &
          (ili.ge.iri .and. ili.le.irf) .or. &
          (ili.ge.ibi .and. ili.le.ibf) .or. &
          (ilf.ge.iri .and. ilf.le.irf) .or. &
          (ilf.ge.ibi .and. ilf.le.ibf) .or. &
          (ibi.ge.iri .and. ibi.le.irf) .or. &
          (ibi.ge.ili .and. ibi.le.ilf) .or. &
          (ibf.ge.iri .and. ibf.le.irf) .or. &
          (ibf.ge.ili .and. ibf.le.ilf) .or. &
          (iri.ge.ibi .and. iri.le.ibf) .or. &
          (iri.ge.ili .and. iri.le.ilf) .or. &
          (irf.ge.ibi .and. irf.le.ibf) .or. &
          (irf.ge.ili .and. irf.le.ilf)) &
          CALL AbortWithError('unable to read input groups')

!     Make sure layer counts are consistent
      rblk=0
      bblk=0
      lblk=0
      DO iblk=1,nblk
        IF (blankline(iblk).gt.iri .and. blankline(iblk).lt.irf) &
           rblk=rblk+1
        IF (blankline(iblk).gt.ibi .and. blankline(iblk).lt.ibf) &
           bblk=bblk+1
        IF (blankline(iblk).gt.ili .and. blankline(iblk).lt.ilf) &
           lblk=lblk+1
      ENDDO
      IF (irf-iri-rblk.ne.2) &
         CALL AbortWithError('Resort section must have exactly 1 line')
      IF (ibf-ibi-bblk.ne.ilf-ili-lblk+1) &
         CALL AbortWithError('Inconsistent input layer numbers')
      ML%nlayr=ibf-ibi-bblk-1
      IF (ML%nlayr.lt.1) CALL AbortWithError('No layers!')

      REWIND(u)

!     Read sections
      ALLOCATE(res_tmp(maxndof),bas_tmp(ML%nlayr,maxndof))
      IF (ML%nlayr.gt.1) THEN
         ALLOCATE(layrs_tmp(ML%nlayr-1,maxndof))
      ENDIF
      res_tmp=0
      bas_tmp=0
      iblk=1
      ilb=0
      ill=0
      DO iline=1,nlines
         IF (iline.eq.blankline(iblk)) THEN
            READ(u,*,IOSTAT=ReadStat)
            iblk=iblk+1
         ELSEIF (iline.gt.iri .and. iline.lt.irf) THEN
            READ(u,*,err=225) (res_tmp(im),im=1,maxndof)
         ELSEIF (iline.gt.ibi .and. iline.lt.ibf) THEN
            ilb=ilb+1
            READ(u,*,err=225) (bas_tmp(ilb,im),im=1,maxndof)
         ELSEIF (iline.gt.ili .and. iline.lt.ilf) THEN
            ill=ill+1
            READ(u,*,err=225) (layrs_tmp(ill,im),im=1,maxndof)
         ELSE
            READ(u,*,IOSTAT=ReadStat)
         ENDIF
225      continue
      ENDDO

!     Determine nmode from length of temporary basis array
      ALLOCATE(ML%nmode(ML%nlayr))
      DO il=1,ML%nlayr
         ML%nmode(il)=0
         DO im=1,maxndof
            IF (bas_tmp(il,im).gt.0) ML%nmode(il)=ML%nmode(il)+1
         ENDDO
      ENDDO
      ML%ndof=ML%nmode(1)

      ALLOCATE(ML%modcomb(ML%nlayr,ML%ndof),ML%gdim(ML%nlayr,ML%ndof))
      ALLOCATE(ML%modstart(ML%nlayr,ML%ndof),ML%resort(ML%ndof))
      ALLOCATE(ML%whichmod(max(1,ML%nlayr-1),ML%ndof))

!     Fill resort array
      ML%resort(1:ML%ndof)=res_tmp(1:ML%ndof)

!     Fill modcomb array
      DO im=1,ML%ndof
         ML%modcomb(1,im)=1
      ENDDO
      DO il=2,ML%nlayr
         DO im=1,ML%nmode(il)
            ML%modcomb(il,im)=layrs_tmp(il-1,im)
         ENDDO
      ENDDO

!     Fill modstart array (points to the 1st mode on the previous layer
!     in the super-mode on the current layer)
      DO im=1,ML%ndof
         ML%modstart(1,im)=im
      ENDDO
      DO il=2,ML%nlayr
         ML%modstart(il,1)=1
         DO im=2,ML%nmode(il)
            ML%modstart(il,im)=ML%modstart(il,im-1)+ML%modcomb(il,im-1)
         ENDDO
      ENDDO

!     Fill whichmod array (points to the super-mode number in next layer
!     where a given mode is found, i.e. the "reverse" of modstart)
      IF (ML%nlayr.gt.1) THEN
         DO il=1,ML%nlayr-1
            DO im=1,ML%nmode(il)
               DO ii=1,ML%nmode(il+1)
                  IF (im.ge.ML%modstart(il+1,ii))&
                  ML%whichmod(il,im)=ii
               ENDDO
            ENDDO
         ENDDO
      ELSE  ! Single layer case
         DO im=1,ML%nmode(1)
            ML%whichmod(1,im)=im
         ENDDO
      ENDIF

!     Fill gdim array
      DO il=1,ML%nlayr
         DO im=1,ML%nmode(il)
            ML%gdim(il,im)=bas_tmp(il,im)
         ENDDO
      ENDDO

      DEALLOCATE(bas_tmp)
      IF (ML%nlayr.gt.1) THEN
         DEALLOCATE(layrs_tmp)
      ENDIF
      DEALLOCATE(blankline)
      CLOSE(u)

      call CPU_TIME(t2)
      mc_time=mc_time+t2-t1

      end subroutine ReadModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints mode combination data

      implicit none
      TYPE (MLtree)        :: ML
      integer :: il,im

      write(*,'(X,A/)') '** Structure of multilayer CP-format tree **'
      write(*,*) 'Number of DOF    : ',ML%ndof
      write(*,*) 'Number of layers : ',ML%nlayr

      write(*,'(/A)') 'Resort modes in Hamiltonian in this order:'
      write(*,1233) (ML%resort(im),im=1,ML%ndof)

      write(*,'(/A)') 'modcomb: # prev. layer modes in current'
      write(*,'(A)') 'v-layer;mode-> '
      write(*,1233) (im,im=1,ML%ndof)
      write(*,1235) '=====',('----',im=1,ML%ndof)
      DO il=1,ML%nlayr
         write(*,1234) il,(ML%modcomb(il,im),im=1,ML%nmode(il))
      ENDDO

      write(*,'(/A)') 'modstart: prev. layer mode where current begins'
      write(*,'(A)') 'v-layer;mode-> '
      write(*,1233) (im,im=1,ML%ndof)
      write(*,1235) '=====',('----',im=1,ML%ndof)
      DO il=1,ML%nlayr
         write(*,1234) il,(ML%modstart(il,im),im=1,ML%nmode(il))
      ENDDO

      write(*,'(/A)') 'whichmod: next layer mode where current is found'
      write(*,'(A)') 'v-layer;mode-> '
      write(*,1233) (im,im=1,ML%ndof)
      write(*,1235) '=====',('----',im=1,ML%ndof)
      DO il=1,max(1,ML%nlayr-1)
         write(*,1234) il,(ML%whichmod(il,im),im=1,ML%nmode(il))
      ENDDO

      write(*,'(/A)') 'gdim: number of basis functions per mode'
      write(*,'(A)') 'v-layer;mode-> '
      write(*,1233) (im,im=1,ML%ndof)
      write(*,1235) '=====',('----',im=1,ML%ndof)
      DO il=1,ML%nlayr
         write(*,1234) il,(ML%gdim(il,im),im=1,ML%nmode(il))
      ENDDO

      write(*,'(/X,A)') '********************************************'

1233  format(6X,32(I4,X))
1234  format(I4,2X,32(I4,X))
1235  format(A5,1X,32(A4,X))

      end subroutine PrintModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Makes sure mode combination data is not bogus

      implicit none
      TYPE (MLtree) :: ML
      integer :: il,im,k,nbloc,nsubm,mstart,sum,prod

      DO im=1,ML%nmode(1)
         sum=ML%resort(im)
         IF (sum.lt.1 .or. sum.gt.ML%nmode(1)) &
            CALL AbortWithError('Bad resort data: mode improperly set')
         DO il=im+1,ML%nmode(1)
            IF (ML%resort(il).eq.sum) &
               CALL AbortWithError('Bad resort data: duplicated mode')
         ENDDO
      ENDDO

!     Check mode counts
      DO il=2,ML%nlayr
         IF (ML%nmode(il).gt.ML%nmode(il-1)) THEN
            write(*,*) 'Modes in layer ',il,' : ',ML%nmode(il)
            write(*,*) 'Modes in layer ',il-1,' : ',ML%nmode(il-1)
            CALL AbortWithError('The number of modes must not&
                 & increase with increasing layer number')
         ENDIF
      ENDDO

!     Check mode mappings
      DO il=2,ML%nlayr
         sum=0
         DO im=1,ML%nmode(il)
            sum=sum+ML%modcomb(il,im)
         ENDDO
         IF (sum.ne.ML%nmode(il-1)) THEN
            write(*,*) '# modes in layer ',il-1,' : ',ML%nmode(il-1)
            write(*,*) 'Of these, ',sum,' are represented in layer ',il
            CALL AbortWithError('DOF are not mapped 1:1 into modes')
         ENDIF
      ENDDO

!     Make sure basis numbers are within acceptable range
      DO il=1,ML%nlayr
         DO im=1,ML%nmode(il)
            IF (ML%gdim(il,im).lt.1) THEN
               write(*,*) 'Layer: ',il,' Mode: ',im, ', nbasis: ',&
               ML%gdim(il,im)
               CALL AbortWithError('Must have >=1 basis fxn per mode')
            ENDIF
            IF (il.gt.1) THEN
               nbloc=ML%gdim(il,im)
               nsubm=ML%modcomb(il,im)
               mstart=ML%modstart(il,im)
               prod=1
               DO k=1,nsubm
                  prod=prod*ML%gdim(il-1,mstart+k-1)
                  IF (nbloc.lt.prod) EXIT
               ENDDO
               IF (nbloc.gt.prod) THEN
                  write(*,*) 'Layer: ',il,' Mode: ',im,&
                  ' # functions desired: ',nbloc,&
                  ' product basis size: ',prod
                  CALL AbortWithError('Product basis exceeded')
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      end subroutine ValidateModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function uppermodenr(ulayr,llayr,lmodenr,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns uppermodenr in ulayr containing lmodenr in llayr

      implicit none
      TYPE (MLtree)       :: ML
      integer, intent(in) :: ulayr,llayr,lmodenr
      integer :: il,uppermodenr

!     Error checking
      IF (lmodenr.lt.1 .or. lmodenr.gt.ML%nmode(llayr)) &
         call AbortWithError('Error in uppermodenr(): bad mode nr.')
      IF (llayr.lt.1 .or. llayr.gt.ML%nlayr .or. ulayr.lt.1 .or. &
          ulayr.gt.ML%nlayr .or. ulayr.lt.llayr) &
         call AbortWithError('Error in uppermodenr(): bad layer nr.')

!     Trace mode from lower layer to upper layer using whichmod array
      uppermodenr=lmodenr
      DO il=llayr,ulayr-1
         uppermodenr=ML%whichmod(il,uppermodenr)
      ENDDO

      end function uppermodenr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function firstmode(ulayr,llayr,umodenr,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds the index of the first mode in the lower-layer (llayr) that 
! belongs to mode umodenr in upper-layer (ulayr)

      implicit none
      TYPE (MLtree)       :: ML
      integer, intent(in) :: ulayr,llayr,umodenr
      integer :: il,firstmode

!     Error checking
      IF (umodenr.lt.1 .or. umodenr.gt.ML%nmode(ulayr)) &
         call AbortWithError('Error in firstmode(): bad mode nr.')
      IF (llayr.lt.1 .or. llayr.gt.ML%nlayr .or. ulayr.lt.1 .or. &
          ulayr.gt.ML%nlayr .or. ulayr.lt.llayr) &
         call AbortWithError('Error in firstmode(): bad layer nr.')

!     Trace mode from upper layer to lower layer using modstart array
      firstmode=umodenr
      DO il=ulayr,llayr+1,-1
         firstmode=ML%modstart(il,firstmode)
      ENDDO

      end function firstmode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function lastmode(ulayr,llayr,umodenr,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds the index of the last mode in the lower-layer (llayr) that 
! belongs to mode umodenr in upper-layer (ulayr)

      implicit none
      TYPE (MLtree)       :: ML
      integer, intent(in) :: ulayr,llayr,umodenr
      integer :: lastmode

!     Error checking
      IF (umodenr.lt.1 .or. umodenr.gt.ML%nmode(ulayr)) &
         call AbortWithError('Error in lastmode(): bad mode nr.')
      IF (llayr.lt.1 .or. llayr.gt.ML%nlayr .or. ulayr.lt.1 .or. &
          ulayr.gt.ML%nlayr .or. ulayr.lt.llayr) &
         call AbortWithError('Error in lastmode(): bad layer nr.')

!     If the upper mode number is the last in ulayr, then
!     lastmode will be the very last mode in llayr
      IF (umodenr.eq.ML%nmode(ulayr)) THEN
         lastmode=ML%nmode(llayr)
!     Otherwise, just call firstmode() with umodenr+1 and then
!     subtract 1 from the result
      ELSE
         lastmode=firstmode(ulayr,llayr,umodenr+1,ML)-1
      ENDIF

      end function lastmode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_ModeComb(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes the MLtree type

      implicit none
      TYPE (MLtree) :: ML

      IF (.not. MC_SETUP) CALL Init_ModeComb

      MC_SETUP=.FALSE.
      IF (ALLOCATED(ML%modcomb)) DEALLOCATE(ML%modcomb)
      IF (ALLOCATED(ML%modstart)) DEALLOCATE(ML%modstart)
      IF (ALLOCATED(ML%whichmod)) DEALLOCATE(ML%whichmod)
      IF (ALLOCATED(ML%gdim)) DEALLOCATE(ML%gdim)
      IF (ALLOCATED(ML%nmode)) DEALLOCATE(ML%nmode)

      end subroutine Flush_ModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
