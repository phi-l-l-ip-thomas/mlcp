!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Manages the multi-layer wavefunction

      USE ERRORTRAP
      USE UTILS
      USE MYMPI

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE MLtree
          integer, dimension(:,:), allocatable :: modcomb,modstart,&
          whichmod,gdim,truncate
          integer, dimension(:), allocatable :: nmode,resort
          integer :: nlayr,ndof,ntrunc
      END TYPE MLtree

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_ModeComb_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_ModeComb_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_ModeComb_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()
      call Get_MPI_Timings('ModeComb module',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_Modecomb_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Flush_ModeComb(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Disposes the MLtree type

      implicit none
      TYPE (MLtree) :: ML

      IF (ALLOCATED(ML%modcomb)) DEALLOCATE(ML%modcomb)
      IF (ALLOCATED(ML%modstart)) DEALLOCATE(ML%modstart)
      IF (ALLOCATED(ML%whichmod)) DEALLOCATE(ML%whichmod)
      IF (ALLOCATED(ML%gdim)) DEALLOCATE(ML%gdim)
      IF (ALLOCATED(ML%truncate)) DEALLOCATE(ML%truncate)
      IF (ALLOCATED(ML%nmode)) DEALLOCATE(ML%nmode)
      IF (ALLOCATED(ML%resort)) DEALLOCATE(ML%resort)

      end subroutine Flush_ModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine StartModeComb(ML,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main routine for reading and processing mode combination data

      implicit none
      TYPE (MLtree) :: ML
      integer, intent(in) :: verbosity
      character(len=64) :: inpfile

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A/)') 'Setting up mode-combination module...'

      inpfile='layers.inp'
      CALL ReadModeDat(ML,inpfile)
      CALL BcastModeDat(ML)
      CALL ValidateModeDat(ML)
      CALL PrintModeDat(ML,verbosity)

      end subroutine StartModeComb

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
      integer :: iri,irf,ili,ilf,ibi,ibf,iti,itf,maxlines,ilb,ill,ilt
      integer :: iline,nlines,iblk,nblk,rblk,bblk,lblk,tblk
      character(len=1024)  :: line
      character(len=32)    :: string
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

!     Read from mpi_io_rank
      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!     For line reading (change if working with a larger system)
      maxndof=1024
      maxlines=1024
      iri=0
      irf=0
      ibi=0
      ibf=0
      ili=0
      ilf=0
      iti=0
      itf=0

!     Open input file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
      IF (InpStat /= 0) THEN
         write(*,*) TRIM(ADJUSTL(fnm)),' not found'
         call AbortWithError("Error reading input file")
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
         ELSEIF (line(1:9) == '$truncate') THEN
            iti=iline
         ELSEIF (line(1:13) == '$end-truncate') THEN
            itf=iline
         ELSEIF (line(linelen:linelen).ne.'/') THEN
            CALL AbortWithError('Error: input list must end with "/"')
         ENDIF
      ENDDO
      nblk=iblk
      nlines=iline

!     Detect input errors:

!     (required) Basis section
      IF (ibi.eq.0 .or. ibf.eq.0 .or. ibi.gt.ibf .or. &
          (ibi.ge.iri .and. ibi.le.irf) .or. &
          (ibi.ge.ili .and. ibi.le.ilf) .or. &
          (ibi.ge.iti .and. ibi.le.itf) .or. &
          (ibf.ge.iri .and. ibf.le.irf) .or. &
          (ibf.ge.ili .and. ibf.le.ilf) .or. &
          (ibf.ge.iti .and. ibf.le.itf)) &
         CALL AbortWithError('ReadModeDat(): bad $basis section')

!     (required) Layers section
      IF (ili.eq.0 .or. ilf.eq.0 .or. ili.gt.ilf .or. &
          (ili.ge.iri .and. ili.le.irf) .or. &
          (ili.ge.ibi .and. ili.le.ibf) .or. &
          (ili.ge.iti .and. ili.le.itf) .or. &
          (ilf.ge.iri .and. ilf.le.irf) .or. &
          (ilf.ge.ibi .and. ilf.le.ibf) .or. &
          (ilf.ge.iti .and. ilf.le.itf)) &
          CALL AbortWithError('ReadModeDat(): bad $layers section')

!     (optional) Resort section
      IF (.not.(iri.eq.0 .and. irf.eq.0)) THEN
         IF (iri.eq.0 .or. irf.eq.0 .or. iri.gt.irf .or. &
             (iri.ge.ibi .and. iri.le.ibf) .or. &
             (iri.ge.ili .and. iri.le.ilf) .or. &
             (iri.ge.iti .and. iri.le.itf) .or. &
             (irf.ge.ibi .and. irf.le.ibf) .or. &
             (irf.ge.ili .and. irf.le.ilf) .or. &
             (irf.ge.iti .and. irf.le.itf)) &
            CALL AbortWithError('ReadModeDat(): bad $resort section')
      ENDIF

!     (optional) Truncation section
      IF (.not.(iti.eq.0 .and. itf.eq.0)) THEN
         IF (iti.eq.0 .or. itf.eq.0 .or. iti.gt.itf .or. &
             (iti.ge.ili .and. iti.le.ilf) .or. &
             (iti.ge.iri .and. iti.le.irf) .or. &
             (iti.ge.ibi .and. iti.le.ibf) .or. &
             (itf.ge.ili .and. itf.le.ilf) .or. &
             (itf.ge.iri .and. itf.le.irf) .or. &
             (itf.ge.ibi .and. itf.le.ibf)) &
            CALL AbortWithError('ReadModeDat(): bad $truncate section')
      ENDIF

!     Make sure layer counts are consistent
      rblk=0
      bblk=0
      lblk=0
      tblk=0
      DO iblk=1,nblk
        IF (blankline(iblk).gt.iri .and. blankline(iblk).lt.irf) &
           rblk=rblk+1
        IF (blankline(iblk).gt.ibi .and. blankline(iblk).lt.ibf) &
           bblk=bblk+1
        IF (blankline(iblk).gt.ili .and. blankline(iblk).lt.ilf) &
           lblk=lblk+1
        IF (blankline(iblk).gt.iti .and. blankline(iblk).lt.itf) &
           tblk=tblk+1
      ENDDO
      IF (.not.(iri.eq.0 .and. irf.eq.0)) THEN
         IF (irf-iri-rblk.gt.2) CALL AbortWithError(&
            'ReadModeDat(): Resort section must have 1 input line')
      ENDIF
      IF (ibf-ibi-bblk.ne.ilf-ili-lblk+1) CALL AbortWithError(&
            'ReadModeDat(): Inconsistent input layer numbers')
      ML%nlayr=ibf-ibi-bblk-1
      ML%ntrunc=itf-iti-tblk-1
      IF (ML%nlayr.lt.1) CALL AbortWithError('No layers!')

      REWIND(u)

!     Read sections
      ALLOCATE(res_tmp(maxndof),bas_tmp(ML%nlayr,maxndof))
      IF (ML%ntrunc.gt.0) THEN
         ALLOCATE(ML%truncate(ML%ntrunc,5))
         ML%truncate(:,:)=-1
      ENDIF
      IF (ML%nlayr.gt.1) THEN
         ALLOCATE(layrs_tmp(ML%nlayr-1,maxndof))
      ENDIF
      res_tmp=0
      bas_tmp=0
      iblk=1
      ilb=0
      ill=0
      ilt=0
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
         ELSEIF (iline.gt.iti .and. iline.lt.itf) THEN
            ilt=ilt+1
            READ(u,*,err=225) (ML%truncate(ilt,im),im=1,5)
         ELSE
            READ(u,*,IOSTAT=ReadStat)
         ENDIF
225      continue
      ENDDO

      ML%ndof=0
      DO im=1,maxndof
         IF (bas_tmp(1,im).gt.0) ML%ndof=ML%ndof+1
      ENDDO

      ALLOCATE(ML%nmode(ML%nlayr),ML%resort(ML%ndof))
      ALLOCATE(ML%modcomb(ML%nlayr,ML%ndof),ML%gdim(ML%nlayr,ML%ndof))
      ALLOCATE(ML%modstart(ML%nlayr,ML%ndof))
      ALLOCATE(ML%whichmod(max(1,ML%nlayr-1),ML%ndof))

!     Determine nmode from length of temporary basis array
      DO il=1,ML%nlayr
         ML%nmode(il)=0
         DO im=1,maxndof
            IF (bas_tmp(il,im).gt.0) ML%nmode(il)=ML%nmode(il)+1
         ENDDO
      ENDDO

!     Fill resort array
      IF ((iri.eq.0 .and. irf.eq.0) .or. (irf-iri-rblk.eq.1)) THEN
         ! No resort section; modes in default order
         DO im=1,ML%ndof
            ML%resort(im)=im
         ENDDO
      ELSE
         ML%resort(1:ML%ndof)=res_tmp(1:ML%ndof)
      ENDIF

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

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine ReadModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveModeDat(ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('layers.inp') with mode combination data to
! file (i.e. a restart file)

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      character(len=48), intent(in) :: fnm
      character(len=64) :: fname,frmt
      integer :: u,il,im
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

      write(fname,'(2A)') TRIM(ADJUSTL(fnm)),'_layers.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fname)), STATUS="UNKNOWN")

!     Write resort section
      write(u,*)
      write(u,'(A)') '$resort'
      write(frmt,'(A,I0,A)') '(',ML%ndof,'(I0,X),A)'
      write(u,frmt) (ML%resort(im),im=1,ML%ndof),'/'
      write(u,'(A)') '$end-resort'

!     Write basis section
      write(u,*)
      write(u,'(A)') '$basis'
      DO il=1,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X),A)'
         write(u,frmt) (ML%gdim(il,im),im=1,ML%nmode(il)),'/'
      ENDDO
      write(u,'(A)') '$end-basis'

!     Write layers section
      write(u,*)
      write(u,'(A)') '$layers'
      DO il=2,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X),A)'
         write(u,frmt) (ML%modcomb(il,im),im=1,ML%nmode(il)),'/'
      ENDDO
      write(u,'(A)') '$end-layers'
      write(u,*)

!     Write truncate section
      IF (ML%ntrunc.gt.0) THEN
         write(u,*)
         write(u,'(A)') '$truncate'
         DO il=1,ML%ntrunc
            write(frmt,'(A)') '(5(I0,X),A)'
            write(u,frmt) (ML%truncate(il,im),im=1,5),'/'
         ENDDO
         write(u,'(A)') '$end-truncate'
         write(u,*)
      ENDIF

      close(u)

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine SaveModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintModeDat(ML,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints mode combination data

      implicit none
      TYPE (MLtree) :: ML
      integer, intent(in) :: verbosity
      integer :: il,im
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

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

      IF (verbosity.ge.1) THEN
         write(*,'(/A)') &
               'modstart: prev. layer mode where current begins'
         write(*,'(A)') 'v-layer;mode-> '
         write(*,1233) (im,im=1,ML%ndof)
         write(*,1235) '=====',('----',im=1,ML%ndof)
         DO il=1,ML%nlayr
            write(*,1234) il,(ML%modstart(il,im),im=1,ML%nmode(il))
         ENDDO

         write(*,'(/A)') &
               'whichmod: next layer mode where current is found'
         write(*,'(A)') 'v-layer;mode-> '
         write(*,1233) (im,im=1,ML%ndof)
         write(*,1235) '=====',('----',im=1,ML%ndof)
         DO il=1,max(1,ML%nlayr-1)
            write(*,1234) il,(ML%whichmod(il,im),im=1,ML%nmode(il))
         ENDDO
      ENDIF

      write(*,'(/A)') 'gdim: number of basis functions per mode'
      write(*,'(A)') 'v-layer;mode-> '
      write(*,1233) (im,im=1,ML%ndof)
      write(*,1235) '=====',('----',im=1,ML%ndof)
      DO il=1,ML%nlayr
         write(*,1234) il,(ML%gdim(il,im),im=1,ML%nmode(il))
      ENDDO

      IF (ML%ntrunc.gt.0) then
         write(*,'(/A)') 'truncation: trim basis based on criteria'
         write(*,'(A)') 'layer - mode: nmode-max sum-max q.n.-max '
         write(*,'(4(A,X))') '=====   =====','---------','-------',&
                         '--------'
         DO il=1,ML%ntrunc
            write(*,1236) ML%truncate(il,1),'-',ML%truncate(il,2),':',&
                          (ML%truncate(il,im),im=3,5)
         ENDDO
      ENDIF

      write(*,'(/X,A)') '********************************************'

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

1233  format(6X,32(I4,X))
1234  format(I4,2X,32(I4,X))
1235  format(A5,1X,32(A4,X))
1236  format(I5,X,A,X,I4,A,X,I9,X,I7,X,I8)

      end subroutine PrintModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ValidateModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Makes sure mode combination data is not bogus

      implicit none
      TYPE (MLtree) :: ML
      integer :: il,im,k,nbloc,nsubm,mstart,sum,prod
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

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
            IF (mpirank.eq.mpi_prnt_rank) THEN
               write(*,*) 'Modes in layer ',il,' : ',ML%nmode(il)
               write(*,*) 'Modes in layer ',il-1,' : ',ML%nmode(il-1)
            ENDIF
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
            IF (mpirank.eq.mpi_prnt_rank) THEN
               write(*,*) '# modes in layer ',il-1,' : ',ML%nmode(il-1)
               write(*,*) 'Of these, ',sum,' are represented in layer ',il
            ENDIF
            CALL AbortWithError('DOF are not mapped 1:1 into modes')
         ENDIF
      ENDDO

!     Make sure basis numbers are within acceptable range
      DO il=1,ML%nlayr
         DO im=1,ML%nmode(il)
            IF (ML%gdim(il,im).lt.1) THEN
               IF (mpirank.eq.mpi_prnt_rank) THEN
                  write(*,*) 'Layer: ',il,' Mode: ',im, ', nbasis: ',&
                  ML%gdim(il,im)
               ENDIF
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
                  IF (mpirank.eq.mpi_prnt_rank) THEN
                     write(*,*) 'Layer: ',il,' Mode: ',im,&
                     ' # functions desired: ',nbloc,&
                     ' product basis size: ',prod
                  ENDIF
                  CALL AbortWithError('Product basis exceeded')
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!     Truncation node validation
      DO k=1,ML%ntrunc
         il=ML%truncate(k,1)
         im=ML%truncate(k,2)

!        Check for out-of-range entries
         IF (il.lt.2 .or. il.gt.ML%nlayr) THEN
            write(*,'(2(A,I0),A)') 'truncation: layer-mode ',&
                                   il,'-',im,', layer is out of range'
            CALL AbortWithError('ValidateModeDat(): bad input')
         ELSEIF (im.lt.1 .or. im.gt.ML%nmode(il)) THEN
            write(*,'(2(A,I0),A)') 'truncation: layer-mode ',&
                                   il,'-',im,', mode is out of range'
            CALL AbortWithError('ValidateModeDat(): bad input')
         ENDIF

!        Truncation nodes must have only one parent
         IF (firstmode(il,il-1,im,ML).ne.lastmode(il,il-1,im,ML)) THEN
            write(*,'(2(A,I0),A)') 'truncation: layer-mode ',&
                         il,'-',im,' must have exactly 1 parent node'
            CALL AbortWithError('ValidateModeDat(): bad input')
         ENDIF

!        Enforce ordering to prevent duplicate entries
         IF (k.gt.1) THEN
            IF (il.lt.ML%truncate(k-1,1) .or. (il.eq.ML%truncate(k-1,1)&
                .and. im.le.ML%truncate(k-1,2))) THEN
               write(*,'(2(A,I0),A)') &
               'truncation: layer-mode ',il,'-',im,&
               ': entries must be in ascending order and not repeated.'
               CALL AbortWithError('ValidateModeDat(): bad input')
            ENDIF
         ENDIF
      ENDDO

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine ValidateModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts MLtree to all MPI ranks

      implicit none
      TYPE (MLtree) :: ML
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      call bcast(ML%nlayr,mpi_io_rank)
      call bcast(ML%ndof,mpi_io_rank)
      call bcast(ML%ntrunc,mpi_io_rank)

      IF (mpirank.ne.mpi_io_rank) THEN
         ALLOCATE(ML%nmode(ML%nlayr),ML%resort(ML%ndof))
         ALLOCATE(ML%modcomb(ML%nlayr,ML%ndof),ML%gdim(ML%nlayr,ML%ndof))
         ALLOCATE(ML%modstart(ML%nlayr,ML%ndof))
         ALLOCATE(ML%whichmod(max(1,ML%nlayr-1),ML%ndof))
         IF (ML%ntrunc.gt.0) ALLOCATE(ML%truncate(ML%ntrunc,5))
      ENDIF

!     Broadcast arrays
      call bcast(ML%nmode,mpi_io_rank)
      call bcast(ML%resort,mpi_io_rank)
      call bcast(ML%modcomb,mpi_io_rank)
      call bcast(ML%gdim,mpi_io_rank)
      call bcast(ML%modstart,mpi_io_rank)
      call bcast(ML%whichmod,mpi_io_rank)
      IF (ML%ntrunc.gt.0) call bcast(ML%truncate,mpi_io_rank)

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine BcastModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CopyModeDat(MLin,MLout)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts MLtree to all MPI ranks

      implicit none
      TYPE (MLtree), INTENT(IN) :: MLin
      TYPE (MLtree), INTENT(INOUT) :: MLout
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      call Flush_ModeComb(MLout)

      MLout%nlayr=MLin%nlayr
      MLout%ndof=MLin%ndof
      MLout%ntrunc=MLin%ntrunc
      
      ALLOCATE(MLout%nmode(MLin%nlayr),MLout%resort(MLin%ndof))
      ALLOCATE(MLout%modcomb(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%modstart(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%gdim(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%whichmod(max(1,MLin%nlayr-1),MLin%ndof))

      MLout%nmode(:)=MLin%nmode(:)
      MLout%resort(:)=MLin%resort(:)
      MLout%modcomb(:,:)=MLin%modcomb(:,:)
      MLout%modstart(:,:)=MLin%modstart(:,:)
      MLout%gdim(:,:)=MLin%gdim(:,:)
      MLout%whichmod(:,:)=MLin%whichmod(:,:)

      IF (MLin%ntrunc.gt.0) THEN
         ALLOCATE(MLout%truncate(MLin%ntrunc,5))
         MLout%truncate(:,:)=MLin%truncate(:,:)
      ENDIF

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine CopyModeDat

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

      function resortedmode(modenr,ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds the index in the resorted list of mode 'modenr'

      implicit none
      TYPE (MLtree), intent(in)  :: ML
      integer, intent(in)        :: modenr
      integer :: resortedmode,i,n

      n=ML%ndof

!     Error checking
      IF (modenr.lt.1 .or. modenr.gt.n) &
         call AbortWithError('resortedmode(): modenr out of range')

      DO i=1,n
         IF (ML%resort(i).eq.modenr) THEN
            resortedmode=i
            EXIT
         ENDIF
      ENDDO

      end function resortedmode


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getnodenumber(il,im,ML) result (k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets node number in sequence (for testing)

      implicit none
      TYPE (MLtree), intent(in)  :: ML
      integer, intent(in)        :: il,im
      integer :: i,j,k

      k=0
      do i=1,ML%nlayr
         do j=1,ML%nmode(i)
            k=k+1
            if (i.ge.il .and. j.ge.im) exit
         enddo
         if (i.ge.il) exit
      enddo

      end function getnodenumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
