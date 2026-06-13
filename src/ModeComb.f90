!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Manages the multi-layer wavefunction

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE INPUTCP

      implicit none
      real(kind=8), allocatable, private :: module_time(:)
      logical, private :: MODULE_SETUP = .FALSE.

      TYPE MLtree
           TYPE (CPPar) :: dpp
           TYPE (CPPar), allocatable :: cpp(:)
           integer, dimension(:,:), allocatable :: modcomb,modstart,&
           whichmod,gdim
           integer, dimension(:), allocatable :: nmode,resort
           integer :: nlayr,ndof,nnodelist,nkey
           character(len=128), allocatable :: keyfields(:), keyvalues(:)
           ! Control parameters
           character(len=128) :: system,pes_path,resfile,pe_transform
           logical            :: dorestart
           integer            :: rs(33)
           real(kind=8)       :: pe_trans_fac
           character(len=128), dimension(6) :: &
           fieldlist=[character(len=128) :: &
                       'pe_trans_fac',&
                       'pe_transform',&
                       'system',&
                       'pes_path',&
                       'resfile',&
                       'rs'&
                      ]
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
      IF (ALLOCATED(ML%cpp)) DEALLOCATE(ML%cpp)
      IF (ALLOCATED(ML%nmode)) DEALLOCATE(ML%nmode)
      IF (ALLOCATED(ML%resort)) DEALLOCATE(ML%resort)
      IF (ALLOCATED(ML%keyfields)) DEALLOCATE(ML%keyfields)
      IF (ALLOCATED(ML%keyvalues)) DEALLOCATE(ML%keyvalues)

      end subroutine Flush_ModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine StartModeComb(ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main routine for reading and processing mode combination data

      implicit none
      TYPE (MLtree) :: ML
      character(len=128),intent(in) :: fnm

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A/)') 'Setting up mode-combination module...'

      CALL ReadModeDat(ML,fnm)
      CALL BcastModeDat(ML)
      CALL ValidateModeDat(ML)
      CALL PrintModeDat(ML)

      end subroutine StartModeComb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadModeDat(ML,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads input file ('layers.inp') containing mode combination data and
! fills type 'MLtree'

      implicit none
      TYPE (MLtree) :: ML
      TYPE (CPPar)  :: cppref
      character(len=128), intent(in) :: fnm
      character(len=128), dimension(:), allocatable :: tfields, tvalues
      character(len=128), dimension(:,:), allocatable :: mfields, mvalues
      integer, allocatable :: nmode_tmp(:),blankline(:),res_tmp(:)
      integer, allocatable :: bas_tmp(:,:),layrs_tmp(:,:),nfield_tmp(:)
      integer, allocatable :: lm_key(:)
      real(kind=8), allocatable :: lm_tmp(:,:)
      integer :: il,im,ii,i2,u,InpStat,ReadStat,linelen,maxndof
      integer :: ici,icf,iri,irf,ili,ilf,ibi,ibf,iti,itf,maxfields
      integer :: ilb,ill,ilt,maxlines,nread
      integer :: iline,nlines,iblk,nblk,cblk,rblk,bblk,lblk,tblk
      character(len=1024)  :: line
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

!     Read from mpi_io_rank
      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!     For line reading (change if working with a larger system)
      maxndof=1024
      maxlines=1024
      maxfields=128
      ici=0
      icf=0
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
         ELSEIF (line(1:8) == '$control') THEN
            ici=iline
         ELSEIF (line(1:12) == '$end-control') THEN
            icf=iline
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
         ELSEIF (line(1:5) == '$node') THEN
            iti=iline
         ELSEIF (line(1:9) == '$end-node') THEN
            itf=iline
         ENDIF
      ENDDO
      nblk=iblk
      nlines=iline

!     Detect input errors:
!     (required) Control section
      IF (ici.eq.0 .or. icf.eq.0 .or. ici.gt.icf .or. &
          (ici.ge.ibi .and. ici.le.ibf) .or. &
          (ici.ge.iri .and. ici.le.irf) .or. &
          (ici.ge.ili .and. ici.le.ilf) .or. &
          (ici.ge.iti .and. ici.le.itf) .or. &
          (icf.ge.ibi .and. icf.le.ibf) .or. &
          (icf.ge.iri .and. icf.le.irf) .or. &
          (icf.ge.ili .and. icf.le.ilf) .or. &
          (icf.ge.iti .and. icf.le.itf)) &
         CALL AbortWithError('ReadModeDat(): bad $control section')

!     (required) Basis section
      IF (ibi.eq.0 .or. ibf.eq.0 .or. ibi.gt.ibf .or. &
          (ibi.ge.ici .and. ibi.le.icf) .or. &
          (ibi.ge.iri .and. ibi.le.irf) .or. &
          (ibi.ge.ili .and. ibi.le.ilf) .or. &
          (ibi.ge.iti .and. ibi.le.itf) .or. &
          (ibf.ge.ici .and. ibf.le.icf) .or. &
          (ibf.ge.iri .and. ibf.le.irf) .or. &
          (ibf.ge.ili .and. ibf.le.ilf) .or. &
          (ibf.ge.iti .and. ibf.le.itf)) &
         CALL AbortWithError('ReadModeDat(): bad $basis section')

!     (required) Layers section
      IF (ili.eq.0 .or. ilf.eq.0 .or. ili.gt.ilf .or. &
          (ili.ge.ici .and. ili.le.icf) .or. &
          (ili.ge.iri .and. ili.le.irf) .or. &
          (ili.ge.ibi .and. ili.le.ibf) .or. &
          (ili.ge.iti .and. ili.le.itf) .or. &
          (ilf.ge.ici .and. ilf.le.icf) .or. &
          (ilf.ge.iri .and. ilf.le.irf) .or. &
          (ilf.ge.ibi .and. ilf.le.ibf) .or. &
          (ilf.ge.iti .and. ilf.le.itf)) &
          CALL AbortWithError('ReadModeDat(): bad $layers section')

!     (optional) Resort section
      IF (.not.(iri.eq.0 .and. irf.eq.0)) THEN
         IF (iri.eq.0 .or. irf.eq.0 .or. iri.gt.irf .or. &
             (iri.ge.ici .and. iri.le.icf) .or. &
             (iri.ge.ibi .and. iri.le.ibf) .or. &
             (iri.ge.ili .and. iri.le.ilf) .or. &
             (iri.ge.iti .and. iri.le.itf) .or. &
             (irf.ge.ici .and. irf.le.icf) .or. &
             (irf.ge.ibi .and. irf.le.ibf) .or. &
             (irf.ge.ili .and. irf.le.ilf) .or. &
             (irf.ge.iti .and. irf.le.itf)) &
            CALL AbortWithError('ReadModeDat(): bad $resort section')
      ENDIF

!     (optional) Node section
      IF (.not.(iti.eq.0 .and. itf.eq.0)) THEN
         IF (iti.eq.0 .or. itf.eq.0 .or. iti.gt.itf .or. &
             (iti.ge.ici .and. iti.le.icf) .or. &
             (iti.ge.ili .and. iti.le.ilf) .or. &
             (iti.ge.iri .and. iti.le.irf) .or. &
             (iti.ge.ibi .and. iti.le.ibf) .or. &
             (itf.ge.ici .and. itf.le.icf) .or. &
             (itf.ge.ili .and. itf.le.ilf) .or. &
             (itf.ge.iri .and. itf.le.irf) .or. &
             (itf.ge.ibi .and. itf.le.ibf)) &
            CALL AbortWithError('ReadModeDat(): bad $node section')
      ENDIF

!     Make sure layer counts are consistent
      cblk=0
      rblk=0
      bblk=0
      lblk=0
      tblk=0
      DO iblk=1,nblk
        IF (blankline(iblk).gt.ici .and. blankline(iblk).lt.icf) &
           cblk=cblk+1
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
      IF (ML%nlayr.lt.1) CALL AbortWithError('No layers!')
      ML%nnodelist=itf-iti-tblk-1

      REWIND(u)

!     Read sections
      ALLOCATE(res_tmp(maxndof),bas_tmp(ML%nlayr,maxndof))
      IF (ML%nlayr.gt.1) ALLOCATE(layrs_tmp(ML%nlayr-1,maxndof))
      IF (ML%nnodelist.gt.0) THEN
         ALLOCATE(nfield_tmp(ML%nnodelist))
         ALLOCATE(mfields(maxfields,ML%nnodelist))
         ALLOCATE(mvalues(maxfields,ML%nnodelist))
         ALLOCATE(lm_tmp(ML%nnodelist,2))
         ALLOCATE(lm_key(ML%nnodelist))
         ALLOCATE(ML%cpp(ML%nnodelist))
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
         ELSEIF (iline.gt.ici .and. iline.lt.icf) THEN
            READ(u,"(A1024)") Line
            call parse_line(line,ML%keyfields,ML%keyvalues)
         ELSEIF (iline.gt.iri .and. iline.lt.irf) THEN
            READ(u,"(A1024)") Line
            nread=stringcountintegers(TRIM(ADJUSTL(line)))
            if (nread.lt.1 .or. nread.gt.maxndof) &
               call AbortWithError("Error reading $resort data")
            res_tmp(1:nread)=string2integerarray(TRIM(ADJUSTL(line)),nread)
         ELSEIF (iline.gt.ibi .and. iline.lt.ibf) THEN
            ilb=ilb+1
            READ(u,"(A1024)") Line
            nread=stringcountintegers(TRIM(ADJUSTL(line)))
            if (nread.lt.1 .or. nread.gt.maxndof) &
               call AbortWithError("Error reading $basis data")
            bas_tmp(ilb,1:nread)=string2integerarray(TRIM(ADJUSTL(line)),nread)
         ELSEIF (iline.gt.ili .and. iline.lt.ilf) THEN
            ill=ill+1
            READ(u,"(A1024)") Line
            nread=stringcountintegers(TRIM(ADJUSTL(line)))
            if (nread.lt.1 .or. nread.gt.maxndof) &
               call AbortWithError("Error reading $layers data")
            layrs_tmp(ill,1:nread)=string2integerarray(TRIM(ADJUSTL(line)),nread)
         ELSEIF (iline.gt.iti .and. iline.lt.itf) THEN
            ilt=ilt+1
            READ(u,"(A1024)") Line
            call parse_line(line,tfields,tvalues)
            nfield_tmp(ilt)=SIZE(tfields)
            mfields(1:nfield_tmp(ilt),ilt)=tfields(1:nfield_tmp(ilt))
            mvalues(1:nfield_tmp(ilt),ilt)=tvalues(1:nfield_tmp(ilt))
            DEALLOCATE(tfields,tvalues)
         ELSE
            READ(u,*,IOSTAT=ReadStat)
         ENDIF
      ENDDO

      ML%nkey=0
      IF (ALLOCATED(ML%keyfields)) ML%nkey=SIZE(ML%keyfields)

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

!     Set control parameters
      call SetMLDefaults(ML)
      call processmlfields(ML)

!     Process node list and extract layer-mode data
      IF (ML%nnodelist.gt.0) THEN
         DO il=1,ML%nnodelist
            ML%cpp(il)=ML%dpp
!           Remove whitespace due to special characters
            DO im=1,nfield_tmp(il)
               call scrubstring(mfields(im,il))
            ENDDO
            call processcppfields(ML%cpp(il),mfields(1:nfield_tmp(il),il),&
                                             mvalues(1:nfield_tmp(il),il))
            lm_key(il)=il
            lm_tmp(il,1)=REAL(ML%cpp(il)%layer)
            lm_tmp(il,2)=REAL(ML%cpp(il)%mode)
         ENDDO

         call hsort2Drlist(lm_key,lm_tmp)

!        Builded sorted node list
         i2=0
         DO il=1,ML%nnodelist
            im=lm_key(il)
            if (NINT(lm_tmp(il,1)).gt.i2) then ! Reset ref. settings
               cppref=ML%dpp
            endif
            ML%cpp(il)=cppref
            call processcppfields(ML%cpp(il),mfields(1:nfield_tmp(im),im),&
                                             mvalues(1:nfield_tmp(im),im))
            if (ML%cpp(il)%mode.eq.0) then ! Set layer ref. settings
               i2=ML%cpp(il)%layer
               cppref=ML%cpp(il)
            endif
         ENDDO

         DEALLOCATE(nfield_tmp,mfields,mvalues,lm_tmp,lm_key)
      ENDIF

      IF (ML%nlayr.gt.1) DEALLOCATE(layrs_tmp)
      DEALLOCATE(bas_tmp,blankline)

      CLOSE(u)

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine ReadModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('layers.inp') with mode combination data to
! file (i.e. a restart file)

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CPPar) :: cppref
      character(len=128)  :: fnm,frmt
      character(len=1024) :: line
      integer :: u,il,im,i2
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

      write(fnm,'(2A)') TRIM(ADJUSTL(ML%resfile)),'_layers.rst'

!     Open output file
      u = LookForFreeUnit()
      OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

!     Write control section
      write(u,*)
      write(u,'(A)') '$control'
      DO il=1,ML%nkey
         write(line,'(4A,X)') &
         TRIM(ADJUSTL(ML%keyfields(il))),"='",&
         TRIM(ADJUSTL(ML%keyvalues(il))),"'"
         write(u,'(A)') TRIM(ADJUSTL(line))
      ENDDO
      write(u,'(A)') '$end-control'

!     Write resort section
      write(u,*)
      write(u,'(A)') '$resort'
      write(frmt,'(A,I0,A)') '(',ML%ndof,'(I0,X))'
      write(u,frmt) (ML%resort(im),im=1,ML%ndof)
      write(u,'(A)') '$end-resort'

!     Write basis section
      write(u,*)
      write(u,'(A)') '$basis'
      DO il=1,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X))'
         write(u,frmt) (ML%gdim(il,im),im=1,ML%nmode(il))
      ENDDO
      write(u,'(A)') '$end-basis'

!     Write layers section
      write(u,*)
      write(u,'(A)') '$layers'
      DO il=2,ML%nlayr
         write(frmt,'(A,I0,A)') '(',ML%nmode(il),'(I0,X))'
         write(u,frmt) (ML%modcomb(il,im),im=1,ML%nmode(il))
      ENDDO
      write(u,'(A)') '$end-layers'
      write(u,*)

!     Write node section
      IF (ML%nnodelist.gt.0) THEN
         write(u,*)
         write(u,'(A)') '$node'

         i2=0
         DO il=1,ML%nnodelist
            if (ML%cpp(il)%layer.gt.i2) then ! Reset ref. settings
               cppref=ML%dpp
            endif
            call WriteCPPInputs(cppref,ML%cpp(il),line)
            if (ML%cpp(il)%mode.eq.0) then ! Set layer ref. settings
               i2=ML%cpp(il)%layer
               cppref=ML%cpp(il)
            endif
            write(u,'(A)') TRIM(ADJUSTL(line))
         ENDDO


         write(u,'(A)') '$end-node'
         write(u,*)
      ENDIF

      close(u)

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine SaveModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintModeDat(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints mode combination data

      implicit none
      TYPE (MLtree) :: ML
      TYPE (CPPar)  :: cppref
      integer :: il,im,i2
      character(len=1024) :: line
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

      write(*,'(X,A/)') '** Input vars (defaults / set in $control) **'
      write(*,'(X,A,2X,A)') 'Hamiltonian to be set up      (system):',&
                             TRIM(ADJUSTL(ML%system))
      write(*,'(X,A,2X,A)') 'Path to potential constants (pes_path):',&
                             TRIM(ADJUSTL(ML%pes_path))
      write(*,'(X,A,2X,A)') 'Restart file name            (resfile):',&
                             TRIM(ADJUSTL(ML%resfile))
      write(*,'(X,A,2X,A)') 'PES transformation type (pe_transform):',&
                             TRIM(ADJUSTL(ML%pe_transform))
      write(*,'(X,A,X,ES11.4)') &
              'PE transform factor     (pe_trans_fac):',ML%pe_trans_fac
      write(*,'(X,A,2X,33(I0,X))') &
                            'Random seed                       (rs):',&
                             (ML%rs(il),il=1,33)
      call PrintCPPInputs(ML%dpp)

      IF (ML%nnodelist.gt.0) then
         write(*,'(/A/)') &
                        '** Input vars specified in $node namelist  **'
         i2=0
         DO il=1,ML%nnodelist
            if (ML%cpp(il)%layer.gt.i2) then ! Reset ref. settings
               cppref=ML%dpp
            endif
            call WriteCPPInputs(cppref,ML%cpp(il),line)
            if (ML%cpp(il)%mode.eq.0) then ! Set layer ref. settings
               i2=ML%cpp(il)%layer
               cppref=ML%cpp(il)
            endif
            write(*,'(A)') TRIM(ADJUSTL(line))
         ENDDO
      ENDIF

      write(*,'(/X,A/)') '** Structure of multilayer CP-format tree **'
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

      IF (ML%dpp%verbosity.ge.1) THEN
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

      write(*,'(/X,A)') '********************************************'

      ENDIF rank0

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

1233  format(6X,32(I4,X))
1234  format(I4,2X,32(I4,X))
1235  format(A5,1X,32(A4,X))
1236  format(I4,4X,I4,2X,A)

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

!     Check parameters read in $control list
      IF (TRIM(ADJUSTL(ML%system)).seq.'none') &
         call AbortWithError(&
         "ValidateModeDat: 'system' not set in $control")
      IF (ML%dpp%layer.ne.0 .or. ML%dpp%mode.ne.0) &
         call AbortWithError(&
         "ValidateModeDat(): 'layer-mode' must not be set in $control")
!      IF (.not. ML%dpp%donode) call AbortWithError(&
!         "ValidateModeDat: 'donode' must not be set to .F. in $control")
      IF (ML%dpp%max_nmode.ne.-1) call AbortWithError(&
         "ValidateModeDat(): 'max_nmode' must not be set in $control")
      IF (ML%dpp%max_sum.ne.-1) call AbortWithError(&
         "ValidateModeDat(): 'max_sum' must not be set in $control")
      IF (ML%dpp%max_qn.ne.-1) call AbortWithError(&
         "ValidateModeDat(): 'max_qn' must not be set in $control")
      IF (ML%dpp%Etarget.ne.0.d0) call AbortWithError(&
         "ValidateModeDat(): 'Etarget' must not be set in $control")
      IF (ML%dpp%nactivations.ne.1) call AbortWithError(&
         "ValidateModeDat(): 'nactivations' must not be set in $control")

!     Check resort section
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
               IF (mpirank.eq.mpi_prnt_rank) write(*,*) &
               'Layer: ',il,' Mode: ',im, ', nbasis: ',ML%gdim(il,im)
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
                  IF (mpirank.eq.mpi_prnt_rank) &
                     write(*,*) 'Layer: ',il,' Mode: ',im,&
                     ' # functions desired: ',nbloc,&
                     ' product basis size: ',prod
                  CALL AbortWithError('Product basis exceeded')
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!     Node namelist validation
      DO k=1,ML%nnodelist
         il=ML%cpp(k)%layer
         im=ML%cpp(k)%mode

!        Check for out-of-range entries
         IF (il.lt.1 .or. il.gt.ML%nlayr) THEN
            IF (mpirank.eq.mpi_prnt_rank) &
            write(*,'(2(A,I0),A)') 'node list: layer-mode ',&
                                  il,'-',im,', layer is out of range'
            CALL AbortWithError('ValidateModeDat(): bad input')
         ELSEIF (im.lt.0 .or. im.gt.ML%nmode(il)) THEN
            IF (mpirank.eq.mpi_prnt_rank) &
            write(*,'(2(A,I0),A)') 'node list: layer-mode ',&
                                  il,'-',im,', mode is out of range'
            CALL AbortWithError('ValidateModeDat(): bad $node input')
         ENDIF

!        Enforce ordering to prevent duplicate entries
         IF (k.gt.1) THEN
            IF (il.lt.ML%cpp(k-1)%layer .or. (il.eq.ML%cpp(k-1)%layer &
                .and. im.le.ML%cpp(k-1)%mode)) THEN
               IF (mpirank.eq.mpi_prnt_rank) write(*,'(2(A,I0),A)') &
               'node list: layer-mode ',il,'-',im,&
               ': entries must be in ascending order and not repeated.'
               CALL AbortWithError('ValidateModeDat(): bad $node input')
            ENDIF
         ENDIF

!        Since setting ncycle=0 forces a dry run, make sure $node
!        namelist sets this consistently with $control
         IF ((ML%dpp%ncycle.eq.0 .and. ML%cpp(k)%ncycle.ne.0) .or. &
             (ML%dpp%ncycle.ne.0 .and. ML%cpp(k)%ncycle.eq.0)) THEN
             IF (mpirank.eq.mpi_prnt_rank) write(*,'(2(A,I0),A)') &
               'node list: layer-mode ',il,'-',im,&
               ': "ncycle" must be set =0 or !=0 same as $control.' 
             call AbortWithError('ValidateModeDat(): bad $node input')
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
      integer :: i
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      call bcast(ML%nlayr,mpi_io_rank)
      call bcast(ML%ndof,mpi_io_rank)
      call bcast(ML%nnodelist,mpi_io_rank)
      call bcast(ML%nkey,mpi_io_rank)
      call bcast(ML%system,mpi_io_rank)
      call bcast(ML%dorestart,mpi_io_rank)
      call bcast(ML%pes_path,mpi_io_rank)
      call bcast(ML%resfile,mpi_io_rank)
      call bcast(ML%pe_transform,mpi_io_rank)
      call bcast(ML%pe_trans_fac,mpi_io_rank)
      call bcast(ML%rs,mpi_io_rank)

      IF (mpirank.ne.mpi_io_rank) THEN
         ALLOCATE(ML%nmode(ML%nlayr),ML%resort(ML%ndof))
         ALLOCATE(ML%modcomb(ML%nlayr,ML%ndof),ML%gdim(ML%nlayr,ML%ndof))
         ALLOCATE(ML%modstart(ML%nlayr,ML%ndof))
         ALLOCATE(ML%whichmod(max(1,ML%nlayr-1),ML%ndof))
         IF (ML%nkey.gt.0) &
            ALLOCATE(ML%keyfields(ML%nkey),ML%keyvalues(ML%nkey))
         IF (ML%nnodelist.gt.0) ALLOCATE(ML%cpp(ML%nnodelist))
      ENDIF

!     Broadcast arrays
      call bcast(ML%nmode,mpi_io_rank)
      call bcast(ML%resort,mpi_io_rank)
      call bcast(ML%modcomb,mpi_io_rank)
      call bcast(ML%gdim,mpi_io_rank)
      call bcast(ML%modstart,mpi_io_rank)
      call bcast(ML%whichmod,mpi_io_rank)
      call BcastCPP(ML%dpp)
      IF (ML%nkey.gt.0) THEN
         DO i=1,ML%nkey
            call bcast(ML%keyfields(i),mpi_io_rank)
            call bcast(ML%keyvalues(i),mpi_io_rank)
         ENDDO
      ENDIF
      IF (ML%nnodelist.gt.0) THEN
         DO i=1,ML%nnodelist
            call BcastCPP(ML%cpp(i))
         ENDDO
      ENDIF

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine BcastModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CopyModeDat(MLin,MLout)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copies ML data MLin -> MLout

      implicit none
      TYPE (MLtree), INTENT(IN) :: MLin
      TYPE (MLtree), INTENT(INOUT) :: MLout
      integer :: i
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_ModeComb_Module()

      call CPU_TIME(ti1)

      call Flush_ModeComb(MLout)

      MLout%nlayr=MLin%nlayr
      MLout%ndof=MLin%ndof
      MLout%nnodelist=MLin%nnodelist
      MLout%nkey=MLin%nkey
      MLout%system=MLin%system
      MLout%pes_path=MLin%pes_path
      MLout%resfile=MLin%resfile
      MLout%pe_transform=MLin%pe_transform
      MLout%pe_trans_fac=MLin%pe_trans_fac
      MLout%rs(:)=MLin%rs(:)
      MLout%dpp=MLin%dpp
      
      ALLOCATE(MLout%nmode(MLin%nlayr),MLout%resort(MLin%ndof))
      ALLOCATE(MLout%modcomb(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%modstart(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%gdim(MLin%nlayr,MLin%ndof))
      ALLOCATE(MLout%whichmod(max(1,MLin%nlayr-1),MLin%ndof))
      ALLOCATE(MLout%keyfields(MLin%nkey))
      ALLOCATE(MLout%keyvalues(MLin%nkey))

      MLout%nmode(:)=MLin%nmode(:)
      MLout%resort(:)=MLin%resort(:)
      MLout%modcomb(:,:)=MLin%modcomb(:,:)
      MLout%modstart(:,:)=MLin%modstart(:,:)
      MLout%gdim(:,:)=MLin%gdim(:,:)
      MLout%whichmod(:,:)=MLin%whichmod(:,:)
      MLout%keyfields(:)=MLin%keyfields(:)
      MLout%keyvalues(:)=MLin%keyvalues(:)

      IF (MLin%nnodelist.gt.0) THEN
         ALLOCATE(MLout%cpp(MLin%nnodelist))
         DO i=1,MLin%nnodelist
            MLout%cpp(i)=MLin%cpp(i)
         ENDDO
      ENDIF

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine CopyModeDat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetMLDefaults(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set MLCP parameter defaults

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML

      call SetCPPDefaults(ML%dpp,0,0)

!     Real parameters
      ML%pe_trans_fac=1.d0

!     Character parameters
      ML%pe_transform='none'
      ML%system='none'
      ML%pes_path='./'
      ML%resfile='none'
      
!     Array parameters
      ML%rs(:)=0

      end subroutine SetMLDefaults

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine processmlfields(ML)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assigns values to CPpar variables if found in input file

      implicit none
      TYPE (MLtree), INTENT(INOUT) :: ML
      character(len=128), allocatable :: cppfields(:),cppvalues(:)
      character(len=128) :: thevalue
      integer :: i,j,n,nfields,nfieldlist,fcount,cppct
      logical :: valid

      nfields=SIZE(ML%keyfields)
      nfieldlist=SIZE(ML%fieldlist)

      ALLOCATE(cppfields(nfieldlist),cppvalues(nfieldlist))
      cppct=0

!     Make sure each field in the input file is found in 'fieldlist'
      do i=1,nfields
         valid=.FALSE.
         call scrubstring(ML%keyfields(i))
         do j=1,nfieldlist
            if ( TRIM(ADJUSTL(ML%keyfields(i))) .seq. &
                 TRIM(ADJUSTL(ML%fieldlist(j))) ) then
                valid=.true.
                exit
            endif
         enddo

!        If not a valid ML entry, check if this is a cpp entry
         if (.not.valid) then
            cppct=cppct+1
            cppfields(cppct)=TRIM(ADJUSTL(ML%keyfields(i)))
            cppvalues(cppct)=TRIM(ADJUSTL(ML%keyvalues(i)))
         endif
      enddo

!     Process the accumulated list of cpp field/value pairs
      IF (cppct.gt.0) call processcppfields(ML%dpp,cppfields(1:cppct),&
                                                   cppvalues(1:cppct))
      DEALLOCATE(cppfields,cppvalues)

!     Real fields
      call get_field(ML%keyfields,ML%keyvalues,'pe_trans_fac',thevalue,fcount)
      if (fcount.eq.1) ML%pe_trans_fac=string2real8(thevalue)
      
!     String fields
      call get_field(ML%keyfields,ML%keyvalues,'pe_transform',thevalue,fcount)
      if (fcount.eq.1) ML%pe_transform=TRIM(ADJUSTL(thevalue))

      call get_field(ML%keyfields,ML%keyvalues,'system',thevalue,fcount)
      if (fcount.eq.1) ML%system=TRIM(ADJUSTL(thevalue))

      call get_field(ML%keyfields,ML%keyvalues,'pes_path',thevalue,fcount)
      if (fcount.eq.1) ML%pes_path=TRIM(ADJUSTL(thevalue))

      call get_field(ML%keyfields,ML%keyvalues,'resfile',thevalue,fcount)
      if (fcount.eq.1) ML%resfile=TRIM(ADJUSTL(thevalue))

!     Integer array fields
      n=SIZE(ML%rs)
      call get_field(ML%keyfields,ML%keyvalues,'rs',thevalue,fcount)
      if (fcount.eq.1) ML%rs=string2integerarray(thevalue,n)

      end subroutine processmlfields

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getnodecppar(ML,k) result (cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns CPPar object to use for node k

      implicit none
      TYPE (MLtree), INTENT(IN) :: ML
      TYPE (CPpar) :: cpp
      integer, intent(in) :: k
      integer :: i,il,im

      call getlayermode(ML,k,il,im)

      cpp=ML%dpp
      DO i=1,ML%nnodelist
         IF (ML%cpp(i)%layer.gt.il) THEN
            EXIT
         ELSEIF (ML%cpp(i)%layer.eq.il) THEN
            IF (ML%cpp(i)%mode.eq.0) THEN
               cpp=ML%cpp(i) ! Overwrite global settings with layer ref
            ELSEIF (ML%cpp(i)%mode.eq.im) THEN
               cpp=ML%cpp(i) ! Overwrite with layer-mode settings
               EXIT
            ELSEIF (ML%cpp(i)%mode.gt.im) THEN
               EXIT
            ENDIF
         ENDIF
      ENDDO

      end function getnodecppar

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
! Gets node number in sequence

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

      subroutine getlayermode(ML,k,il,im)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gets layer and mode corresponding to node number k

      implicit none
      TYPE (MLtree), intent(in)  :: ML
      integer, intent(in)  :: k
      integer, intent(out) :: il,im
      integer :: i,j

      il=0
      im=0
      j=0
      do i=1,ML%nlayr
         j=j+ML%nmode(i)
         if (j.ge.k) then
            il=i
            im=k-j+ML%nmode(i)
            exit
         endif
      enddo

      if (il.eq.0 .or. im.eq.0) then
          write(*,*) 'Node number',k,' is out of range'
          call AbortWithError("getlayermode(): bad node number")
      endif

      end subroutine getlayermode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODECOMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
