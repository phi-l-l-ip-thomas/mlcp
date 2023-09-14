!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODEGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Manages the multi-layer wavefunction

      USE ERRORTRAP
      USE UTILS
      USE MODECOMB

      implicit none
      TYPE GridPar
          real*8, allocatable  :: bounds(:,:)
          integer, allocatable :: npts(:)
          logical, allocatable :: egrid(:)
          character(3), allocatable :: bc(:)
      END TYPE GridPar

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadGridDat(gp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads input file ('grids.inp') and fills type GridPar

      implicit none
      TYPE (GridPar), INTENT(OUT)   :: gp
      character(len=64), intent(in) :: fnm      
      integer, allocatable :: blankline(:),ninpl(:),secbound(:,:)
      integer, allocatable :: sectyp(:),currln(:),inttmp(:,:,:),blks(:)
      real*8, allocatable  :: realtmp(:,:,:)
      integer :: i,j,k,u,InpStat,ReadStat,linelen,maxndof,nsect
      integer :: maxline,iline,nlines,iblk,nblk,ndof
      real*8, parameter   :: initval=12345678.90
      character(len=1024) :: line
      character(len=64), ALLOCATABLE :: secnm(:),sechf(:,:)
      logical :: lineid

!     For line reading (change if working with a larger system)
      maxndof=1024
      maxline=1024

!!!   Input file specific material !!!

      nsect=4

!     The sections to look for are included here
      ALLOCATE(secnm(nsect),sectyp(nsect),ninpl(nsect))
      ninpl=0

!     Names of the sections
      secnm=(/'basis     ',&
              'boundary  ',&
              'gridtype  ',&
              'parameters'/)

!     Section types: 1=integer,2=real
      sectyp=(/1,1,1,2/)

!     Here is the number of input lines expected for each section
      ninpl=(/1,1,1,2/)

!!!   End input file specific material !!!

!     Arrays for the names to search, start/ending lines of sections
      ALLOCATE(secbound(nsect,2),sechf(nsect,2))
      secbound=0
      DO j=1,nsect
         write(sechf(j,1),'(2A)') '$',TRIM(ADJUSTL(secnm(j)))
         write(sechf(j,2),'(2A)') '$end-',TRIM(ADJUSTL(secnm(j)))
      ENDDO

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         call AbortWithError("Oh, no! Error reading input file")
      ENDIF

      ALLOCATE(blankline(maxline))
      iline=0
      iblk=0
      DO
!        Read line
         READ(u,"(A1024)",IOSTAT=ReadStat) Line
         IF ( ReadStat /= 0 ) EXIT
         iline=iline+1
!        remove leading/trailing spaces...
         line = TRIM(ADJUSTL(line))
!        get length of line without leading/trailing spaces...
         linelen = LEN_TRIM(line)
!        keep track of empty or comment lines
         IF (linelen.eq.0 .or. line(1:1).eq.'#') THEN
            iblk=iblk+1
            blankline(iblk)=iline
!        Determine where sections begin and end
         ELSE
            lineid=.FALSE.
            DO j=1,nsect
               DO k=1,2
                  IF (line.eq.TRIM(ADJUSTL(sechf(j,k)))) THEN
                     secbound(j,k)=iline
                     lineid=.TRUE.
                  ENDIF
               ENDDO
            ENDDO
!           Make sure lines within sections are properly terminated
            IF (.not.lineid .and. line(linelen:linelen).ne.'/') &
               call AbortWithError( &
               'ReadGridDat(): input list must end with "/"')
         ENDIF
      ENDDO
      nblk=iblk
      nlines=iline

      DEALLOCATE(sechf)

!     Detect input errors
      IF (ANY(secbound(:,:).eq.0)) &
         call AbortWithError('ReadGridDat(): missing section flag')
      DO i=1,nsect
         IF (secbound(i,1).gt.secbound(i,2)) call &
            AbortWithError('ReadGridDat(): section flag order reversed')
         DO j=1,nsect
            IF (j.ne.i) THEN
               DO k=1,2
                  IF (secbound(i,k).ge.secbound(j,1) .and. &
                     secbound(i,k).le.secbound(j,2)) call &
                   AbortWithError('ReadGridDat(): sections intertwined')
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     Count number of blank lines in section to be subtracted from total
      ALLOCATE(blks(nsect))
      blks=0
      DO iblk=1,nblk
         DO j=1,nsect
            IF (blankline(iblk).gt.secbound(j,1) .and. &
                blankline(iblk).lt.secbound(j,2)) blks(j)=blks(j)+1
         ENDDO
      ENDDO

!     Make sure each input section has correct number of lines
      DO j=1,nsect
         IF (secbound(j,2)-secbound(j,1)-blks(j).ne.ninpl(j)+1) THEN
            write(*,'(A,X,A,X,I0,X,A)') TRIM(ADJUSTL(secnm(j))),&
            'section must have exactly',ninpl(j),'lines'
            call AbortWithError('ReadGridDat(): sections read error')
         ENDIF
      ENDDO

      REWIND(u)

!     Read sections
      ALLOCATE(currln(nsect),inttmp(nsect,maxndof,MAXVAL(ninpl)))
      ALLOCATE(realtmp(nsect,maxndof,MAXVAL(ninpl)))
      currln=1
      inttmp=-1
      realtmp=initval
      iblk=1
      DO iline=1,nlines
         IF (iline.eq.blankline(iblk)) THEN
            READ(u,*,IOSTAT=ReadStat)
            iblk=iblk+1
         ELSE
            lineid=.FALSE.
            DO j=1,nsect
               IF (iline.gt.secbound(j,1) .and. & 
                   iline.lt.secbound(j,2)) THEN
                   IF (sectyp(j).eq.1) THEN
                      READ(u,*,err=225) &
                          (inttmp(j,k,currln(j)),k=1,maxndof)
                   ELSEIF (sectyp(j).eq.2) THEN
                      READ(u,*,err=225) &
                          (realtmp(j,k,currln(j)),k=1,maxndof)
                   ELSE
                      call AbortWithError(&
                           'ReadGridDat(): unrecognized section type')
                   ENDIF
                   currln(j)=currln(j)+1
                   lineid=.TRUE.
               ENDIF
            ENDDO
            IF (.not.lineid) READ(u,*,IOSTAT=ReadStat)
         ENDIF
225      continue
      ENDDO

      DEALLOCATE(secnm,secbound,sectyp,ninpl,blks,blankline,currln)


!!!   Input file specific material !!!

!     Determine number of DOFs from temporary basis array
      ndof=0
      DO i=1,maxndof
         IF (inttmp(1,i,1).eq.-1) EXIT
         ndof=ndof+1
      ENDDO

!     Fill the npts array
      ALLOCATE(gp%npts(ndof))
      gp%npts(1:ndof)=inttmp(1,1:ndof,1)

!     Determine number of DOFs from temporary boundary array
      ndof=0
      DO i=1,maxndof
         IF (inttmp(2,i,1).eq.-1) EXIT
         ndof=ndof+1
      ENDDO

      IF (ndof.ne.SIZE(gp%npts)) call AbortWithError( &
          'ReadGridDat(): ndof mismatch, boundary section')

!     Fill/validate the boundary array
      ALLOCATE(gp%bc(ndof))
      DO i=1,ndof
         IF (inttmp(2,i,1).eq.0) THEN
            gp%bc(i)='fin'
         ELSEIF (inttmp(2,i,1).eq.1) THEN
            gp%bc(i)='sem'
         ELSEIF (inttmp(2,i,1).eq.2) THEN
            gp%bc(i)='inf'
         ELSEIF (inttmp(2,i,1).eq.3) THEN
            gp%bc(i)='per'
         ELSEIF (inttmp(2,i,1).eq.4) THEN
            gp%bc(i)='cos'
         ELSE
            write(*,*) '$boundary, entry:',i
            write(*,*) ' must be one of : 0 (finite),'
            write(*,*) '                  1 (semi-infinite),' 
            write(*,*) '                  2 (infinite)'
            write(*,*) '                  3 (periodic)'
            write(*,*) '                  4 (cosine)'
            call AbortWithError(&
            'ReadGridDat(): unrecognized boundary condition')
         ENDIF
      ENDDO

!     Determine number of DOFs from temporary gridtype array
      ndof=0
      DO i=1,maxndof
         IF (inttmp(3,i,1).eq.-1) EXIT
         ndof=ndof+1
      ENDDO

      IF (ndof.ne.SIZE(gp%npts)) call AbortWithError( &
          'ReadGridDat(): ndof mismatch, gridtype section')

!     Fill/validate the boundary type array
      ALLOCATE(gp%egrid(ndof))
      DO i=1,ndof
         IF (inttmp(3,i,1).eq.0) THEN
            gp%egrid(i)=.FALSE.
         ELSEIF (inttmp(3,i,1).eq.1) THEN
            gp%egrid(i)=.TRUE.
         ELSE
            write(*,*) '$gridtype, entry:',i
            write(*,*) ' must be one of : 0 (nodes grid),'
            write(*,*) '               or 1 (extrema grid)'
            call AbortWithError(&
            'ReadGridDat(): unrecognized grid type')
         ENDIF
      ENDDO

!     Determine number of DOFs from temporary grid type array
      DO k=1,2
         ndof=0
         DO i=1,maxndof
            IF (realtmp(4,i,k).eq.initval) EXIT
            ndof=ndof+1
         ENDDO
         IF (ndof.ne.SIZE(gp%npts)) call AbortWithError( &
            'ReadGridDat(): ndof mismatch, parameters section')
      ENDDO

!     Fill/validate the boundary type array
      ALLOCATE(gp%bounds(ndof,2))
      gp%bounds(1:ndof,1:2)=realtmp(4,1:ndof,1:2)

!!!   End input file specific material !!!

      DEALLOCATE(inttmp,realtmp)

      CLOSE(u)

      end subroutine ReadGridDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FlushGridDat(gp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes the GridPar type

      implicit none
      TYPE (GridPar) :: gp

      IF (ALLOCATED(gp%bounds)) DEALLOCATE(gp%bounds)
      IF (ALLOCATED(gp%npts)) DEALLOCATE(gp%npts)
      IF (ALLOCATED(gp%egrid)) DEALLOCATE(gp%egrid)
      IF (ALLOCATED(gp%bc)) DEALLOCATE(gp%bc)

      end subroutine FlushGridDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODEGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
