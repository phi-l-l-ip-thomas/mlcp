!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
       
      implicit none

      TYPE CPpar
           integer :: ncycle,npow,lowmem,truncation,algo
           integer :: psirank,hrank,psinals,hnals,verbosity
           integer :: rs(33)
           real(kind=8)  :: etarget,solvtol,alspenalty,pe_trans_fac
           logical :: update,dorestart,dotopnode
           character(len=64) :: resfile
           character(len=64) :: system
           character(len=64) :: solver,calcbounds
           character(len=64) :: red2D,redND
           character(len=64) :: h_sort_alg
           character(len=64) :: als_linsys_alg
           character(len=64) :: pe_transform
           ! 'fieldlist' holds parameters that can be set in input file
           character(len=64), dimension(26) :: &
           fieldlist=(/&
                       'update',&
                       'dotopnode',&
                       'ncycle',&
                       'npow',&
                       'psirank',&
                       'psinals',&
                       'hrank',&
                       'hnals',&
                       'h_sort_alg',&
                       'verbosity',&
                       'alspenalty',&
                       'als_linsys_alg',&
                       'algo',&
                       'lowmem',&
                       'truncation',&
                       'etarget',&
                       'solvtol',&
                       'pe_trans_fac',&
                       'pe_transform',&
                       'system',&
                       'solver',&
                       'calcbounds',&
                       'red2D',&
                       'redND',&
                       'resfile',&
                       'rs'&
                      /)
      END TYPE CPpar

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine StartInputCP(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main routine for reading and processing mode combination data

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64) :: inpfile

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(/X,A/)') 'Reading input file (CP.inp)...'

      inpfile='CP.inp'
      CALL SetMLCPparameterDefaults(cpp)
      CALL ReadMLCPInputs(cpp,inpfile)
      CALL BcastMLCPInputs(cpp)
      CALL PrintMLCPInputs(cpp)

      end subroutine StartInputCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetMLCPparameterDefaults(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set MLCP parameter defaults

      implicit none
      TYPE (CPpar) :: cpp

!     Logical parameters
      cpp%update=.true.
      cpp%dorestart=.false.
      cpp%dotopnode=.true.

!     Integer parameters
      cpp%ncycle=10
      cpp%npow=10
      cpp%psirank=10
      cpp%psinals=10
      cpp%hrank=0
      cpp%hnals=200
      cpp%verbosity=0
      cpp%algo=-1
      cpp%lowmem=2
      cpp%truncation=0

!     Real parameters
      cpp%etarget=0.d0
      cpp%solvtol=1.d-10
      cpp%alspenalty=1.d-10
      cpp%pe_trans_fac=1.d0

!     Character parameters
      cpp%pe_transform='none'
      cpp%system='CpOsc'
      cpp%solver='powr'
      cpp%calcbounds='calc'
      cpp%red2D='SVD'
      cpp%redND='ALS'
      cpp%resfile='none'
      cpp%h_sort_alg='pack'
      cpp%als_linsys_alg='LU'

!     Array parameters
      cpp%rs(:)=0
!      cpp%rs(27:32)=(/14,2,2022,20,21,13/)

      end subroutine SetMLCPparameterDefaults

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadMLCPInputs(cpp,fnm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Reads input file for CP-format code for computing eigenvalues

      implicit none
      TYPE (CPpar) :: cpp
      character(len=64), intent(in) :: fnm
      character(1024) :: line
      character(len=128), dimension(:), allocatable :: fields, values
      integer :: i,u,InpStat

!     Read from mpi_io_rank 
      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

!        Open input file
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="OLD", IOSTAT=InpStat)
         IF (InpStat /= 0) THEN
            write(*,*) TRIM(ADJUSTL(fnm)),' not found'
           call AbortWithError("Error reading input file")
         ENDIF
         rewind(u)

!        Read loop
         DO
            READ(u,"(A1024)",IOSTAT=InpStat) line
            IF (InpStat /= 0) EXIT
            call parse_line(line,fields,values)
         ENDDO

         CLOSE(u)

         IF (.not.allocated(fields)) call &
            AbortWithError('ReadMLCPInputs(): no valid inputs found!')

         call processfieldlist(cpp,fields,values)
         deallocate(fields,values)

      ENDIF rank0

      end subroutine ReadMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine parse_line(line,fields,vals)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parse line for input fields and values

      implicit none
      character(len=*), intent(in) :: line
      character(len=*), dimension(:), allocatable :: fields, vals
      character(len=128), dimension(:), allocatable :: field_tmp, val_tmp
      character(len=128) :: field, val
      character(len=1)   :: q
      integer :: i,fi,ff,qi,qf,m,fct,itn
      integer :: nfields, nvals

      if (allocated(fields)) then
         nfields=SIZE(fields)
      else
         nfields=0
      endif

      if (allocated(vals)) then
         nvals=SIZE(vals)
      else
         nvals=0
      endif

      if (nvals.ne.nfields) then
         write(*,*) 'nfields = ',nfields,'; nvals = ',nvals
         call AbortWithError('parse_line(): nfields, nvals must be =')
      endif

      do itn=1,2
         fct=0
         m=1
         fi=-1
         ff=-1
         q='$'
         qi=-1
         qf=-1
         do while (m.lt.1024)

            if (line(m:m).eq.'=' .and. qi.eq.-1) then
               ff=m-1 ! field ends before equals
               m=m+1
               qi=m   ! val must begin with quote after equals
               q=line(m:m)
               if (q.ne.'"' .and. q.ne."'") then ! select single or double quotes
                  call AbortWithError(&
                  'parse_line(): quoted section must follow =')
               endif

            elseif (line(m:m).eq.' ') then
               if (qi.eq.-1) then
                  fi=m+1 ! advance the field begin tag if not inside a val
               endif

            elseif (line(m:m).eq.q(1:1)) then
               if (q(1:1).ne.'$') then ! finishing a field-val pair
                  qf=m
                  fct=fct+1
                  if (itn.eq.2) then
                     write(fields(nfields+fct),'(A)') line(fi:ff)
                     write(vals(nfields+fct),'(A)') line(qi+1:qf-1)
                  endif
                  fi=-1
                  ff=-1
                  q='$'
                  qi=-1
                  qf=-1
               endif

            endif
            m=m+1
         enddo

!        Extend length of 'fields' by number found in this input line
         if (itn.eq.1) then
            if (fct.eq.0) exit

            if (nfields.gt.0) then
               allocate(field_tmp(nfields),val_tmp(nfields))
               field_tmp(1:nfields)=fields(1:nfields)
               val_tmp(1:nfields)=vals(1:nfields)
               deallocate(fields,vals)
            endif

            allocate(fields(nfields+fct),vals(nfields+fct))

            if (nfields.gt.0) then
               fields(1:nfields)=field_tmp(1:nfields)
               vals(1:nfields)=val_tmp(1:nfields)
               deallocate(field_tmp,val_tmp)
            endif
         endif

      enddo ! itn

      end subroutine parse_line

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine processfieldlist(cpp,fields,values)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assigns values to CPpar variables if found in input file

      implicit none
      TYPE (CPpar), intent(inout) :: cpp
      character(len=*), dimension(:), intent(in) :: fields, values
      character(len=256) :: thevalue
      integer :: i,j,n,nfields,nfieldlist,fcount
      logical :: valid

      nfields=SIZE(fields)
      nfieldlist=SIZE(cpp%fieldlist)

!     Make sure each field in the input file is found in 'fieldlist'
      do i=1,nfields 
         valid=.FALSE.
         do j=1,nfieldlist
            if (TRIM(ADJUSTL(fields(i))).seq.&
                TRIM(ADJUSTL(cpp%fieldlist(j)))) then
                valid=.true.
                exit
            endif
         enddo
         if (.not.valid) then
             write(*,*) 'Input item "',TRIM(ADJUSTL(fields(i))),&
                        '" does not correspond to a valid input flag'
             call AbortWithError('processfieldlist(): bad input item')
         endif 
      enddo

!     Check for items in field list; overwrite default if found

!     Logical fields
      call get_field(fields,values,'update',thevalue,fcount)
      if (fcount.eq.1) cpp%update=string2logical(thevalue)

      call get_field(fields,values,'dotopnode',thevalue,fcount)
      if (fcount.eq.1) cpp%dotopnode=string2logical(thevalue)

!     Integer fields
      call get_field(fields,values,'ncycle',thevalue,fcount)
      if (fcount.eq.1) cpp%ncycle=string2integer(thevalue)

      call get_field(fields,values,'npow',thevalue,fcount)
      if (fcount.eq.1) cpp%npow=string2integer(thevalue)

      call get_field(fields,values,'psirank',thevalue,fcount)
      if (fcount.eq.1) cpp%psirank=string2integer(thevalue)

      call get_field(fields,values,'psinals',thevalue,fcount)
      if (fcount.eq.1) cpp%psinals=string2integer(thevalue)

      call get_field(fields,values,'hrank',thevalue,fcount)
      if (fcount.eq.1) cpp%hrank=string2integer(thevalue)

      call get_field(fields,values,'hnals',thevalue,fcount)
      if (fcount.eq.1) cpp%hnals=string2integer(thevalue)

      call get_field(fields,values,'verbosity',thevalue,fcount)
      if (fcount.eq.1) cpp%verbosity=string2integer(thevalue)

      call get_field(fields,values,'algo',thevalue,fcount)
      if (fcount.eq.1) cpp%algo=string2integer(thevalue)

      call get_field(fields,values,'lowmem',thevalue,fcount)
      if (fcount.eq.1) cpp%lowmem=string2integer(thevalue)

      call get_field(fields,values,'truncation',thevalue,fcount)
      if (fcount.eq.1) cpp%truncation=string2integer(thevalue)

!     Real fields
      call get_field(fields,values,'etarget',thevalue,fcount)
      if (fcount.eq.1) cpp%etarget=string2real8(thevalue)

      call get_field(fields,values,'solvtol',thevalue,fcount)
      if (fcount.eq.1) cpp%solvtol=string2real8(thevalue)

      call get_field(fields,values,'alspenalty',thevalue,fcount)
      if (fcount.eq.1) cpp%alspenalty=string2real8(thevalue)

      call get_field(fields,values,'pe_trans_fac',thevalue,fcount)
      if (fcount.eq.1) cpp%pe_trans_fac=string2real8(thevalue)

!     String fields
      call get_field(fields,values,'pe_transform',thevalue,fcount)
      if (fcount.eq.1) cpp%pe_transform=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'system',thevalue,fcount)
      if (fcount.eq.1) cpp%system=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'solver',thevalue,fcount)
      if (fcount.eq.1) cpp%solver=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'calcbounds',thevalue,fcount)
      if (fcount.eq.1) cpp%calcbounds=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'red2D',thevalue,fcount)
      if (fcount.eq.1) cpp%red2D=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'redND',thevalue,fcount)
      if (fcount.eq.1) cpp%redND=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'h_sort_alg',thevalue,fcount)
      if (fcount.eq.1) cpp%h_sort_alg=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'als_linsys_alg',thevalue,fcount)
      if (fcount.eq.1) cpp%als_linsys_alg=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'resfile',thevalue,fcount)
      if (fcount.eq.1) cpp%resfile=TRIM(ADJUSTL(thevalue))

!     Integer array fields
      n=SIZE(cpp%rs)
      call get_field(fields,values,'rs',thevalue,fcount)
      if (fcount.eq.1) cpp%rs=string2integerarray(thevalue,n)

      end subroutine processfieldlist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine get_field(fields, values, thefield, thevalue, fcount)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Looks for field-value pair in list of fields, values, and extracts the
! value if found

      implicit none
      character(len=*), dimension(:), intent(in) :: fields, values
      character(len=*), intent(in)  :: thefield
      character(len=*), intent(out) :: thevalue
      character(len=1) :: q
      integer :: fcount
      integer :: i,nfields

      nfields=SIZE(fields)
      fcount=0
      thevalue=''
      do i=1,nfields
         if (trim(adjustl(fields(i))).seq.trim(adjustl(thefield))) then
            write(thevalue,'(A)') values(i)
            fcount=fcount+1
         endif
      enddo

      if (fcount.gt.1) then
         write(*,*) 'Input field "',trim(adjustl(thefield)),&
                    '" must appear [0,1]x, appears ',fcount
         call AbortWithError('get_field(): duplicated field')
      endif

      end subroutine get_field

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SaveMLCPInputFile(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Regurgitates input file ('CP.inp') to another file for restart

      implicit none
      TYPE (CPpar), intent(in) :: cpp
      character(len=64) :: fnm,fa,fi,fl,fr,fii
      integer :: u,j,nrs

      rank0 : IF (mpirank.eq.mpi_io_rank) THEN

         nrs=SIZE(cpp%rs)

         fa='(3A)'
         fi='(A,I0,A)'
         fl='(A,L0,A)'
         fr='(A,ES14.6,A)'
         write(fii,'(A,I0,A)') '(A,',nrs,'(I0,X),A)'

         write(fnm,'(2A)') TRIM(ADJUSTL(cpp%resfile)),'_CP.rst'

!        Open output file
         u = LookForFreeUnit()
         OPEN(u, FILE=TRIM(ADJUSTL(fnm)), STATUS="UNKNOWN")

         write(u,*) "* Physical system:"
         write(u,fa) "system='",TRIM(ADJUSTL(cpp%system)),"'"
         write(u,*)
         write(u,*) "* Eigensolver algorithm:"
         write(u,fa) "solver='",TRIM(ADJUSTL(cpp%solver)),"'"
         write(u,*)
         write(u,*) "* Calculate spectral bounds using power method:"
         write(u,fa) "calcbounds='",TRIM(ADJUSTL(cpp%calcbounds)),"'"
         write(u,*)
         write(u,*) "* Rank-reduction for 2D nodes:"
         write(u,fa) "red2D='",TRIM(ADJUSTL(cpp%red2D)),"'"
         write(u,*)
         write(u,*) "* Rank-reduction for >2D nodes:"
         write(u,fa) "redND='",TRIM(ADJUSTL(cpp%redND)),"'"
         write(u,*)
         write(u,*) "* Restart file name (set to 'none' to omit):"
         write(u,fa) "resfile='",TRIM(ADJUSTL(cpp%resfile)),"'"
         write(u,*)
         write(u,*) "* Number of solver cycles:"
         write(u,fi) "ncycle='",cpp%ncycle,"'"
         write(u,*)
         write(u,*) "* Number of power iterations:"
         write(u,fi) "npow='",cpp%npow,"'"
         write(u,*)
         write(u,*) "* Target reduction rank for eigenstates"
         write(u,fi) "psirank='",cpp%psirank,"'"
         write(u,*)
         write(u,*) "* Number of ALS iterations for psi reduction"
         write(u,fi) "psinals='",cpp%psinals,"'"
         write(u,*)
         write(u,*) "* Target reduction rank for Hamiltonian:"
         write(u,*) "  (set to 0 to bypass H reduction)"
         write(u,fi) "hrank='",cpp%hrank,"'"
         write(u,*)
         write(u,*) "* Number of ALS iterations for H reduction:"
         write(u,fi) "hnals='",cpp%hnals,"'"
         write(u,*)
         write(u,*) "* Algorithm used to sort H terms:"
         write(u,fa) "h_sort_alg='",TRIM(ADJUSTL(cpp%h_sort_alg)),"'"
         write(u,*)
         write(u,*) "* Potential energy transformation:"
         write(u,fa) "pe_transform='",&
                      TRIM(ADJUSTL(cpp%pe_transform)),"'"
         write(u,*)
         write(u,*) "* PE transform scaling factor:"
         write(u,fr) "pe_trans_fac='",cpp%pe_trans_fac,"'"
         write(u,*)
         write(u,*) "* Print verbosity (less <-{0,1,2,3}-> more):"
         write(u,fi) "verbosity='",cpp%verbosity,"'"
         write(u,*)
         write(u,*) "* ALS regularization penalty:"
         write(u,fr) "alspenalty='",cpp%alspenalty,"'"
         write(u,*)
         write(u,*) "* ALS linear system algorithm:"
         write(u,fa) "als_linsys_alg='",TRIM(ADJUSTL(cpp%als_linsys_alg)),"'"
         write(u,*)
         write(u,*) "* Numerical algorithm (-1 legacy, 0 CPU, 1 GPU):"
         write(u,fi) "algo='",cpp%algo,"'"
         write(u,*)
         write(u,*) "* Low memory option:"
         write(u,fi) "lowmem='",cpp%lowmem,"'"
         write(u,*)
         write(u,*) "* Basis truncation option:"
         write(u,fi) "truncation='",cpp%truncation,"'"
         write(u,*)
         write(u,*) "* Update psi with solutions of subspace problem:"
         write(u,fl) "update='",cpp%update,"'"
         write(u,*)
         write(u,*) "* Process top node of tree:"
         write(u,fl) "dotopnode='",cpp%dotopnode,"'"
         write(u,*)
         write(u,*) "* Solver target energy:"
         write(u,fr) "Etarget='",cpp%etarget,"'"
         write(u,*)
         write(u,*) "* Solver convergence tolerance:"
         write(u,fr) "solvtol='",cpp%solvtol,"'"
         write(u,*)
         write(u,*) "* Random seed:"
         write(u,fii) "rs='",(cpp%rs(j),j=1,nrs),"'"

         CLOSE(u)

      ENDIF rank0

      end subroutine SaveMLCPInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CP.inp

      implicit none
      TYPE (CPpar),INTENT(IN) :: cpp
      integer :: i

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

      write(*,'(X,A/)') '********** Input parameters read: ***********'
      write(*,'(X,A,2X,A5)') 'The Hamiltonian will be set up for    :',&
                             cpp%system
      write(*,'(X,A,4X,A3)') 'Reduction type for 2-D modes   (red2D):',&
                             cpp%red2D
      write(*,'(X,A,4X,A3)') 'Reduction type for >2-D modes  (redND):',&
                             cpp%redND
      write(*,'(X,A,2X,I5)') 'Wavefunction reduced rank    (psirank):',&
                             cpp%psirank
      write(*,'(X,A,2X,I5)') 'Hamiltonian reduced rank       (hrank):',&
                             cpp%hrank
      write(*,'(X,A,2X,I5)') 'Number of ALS iterations-w.f.(psinals):',&
                             cpp%psinals
      write(*,'(X,A,2X,I5)') 'Number of ALS iterations-H     (hnals):',&
                             cpp%hnals
      write(*,'(X,A,4X,A3)') 'H sorting algorithm       (h_sort_alg):',&
                             cpp%h_sort_alg
      write(*,'(X,A,2X,A5)') 'PES transformation type (pe_transform):',&
                             cpp%pe_transform
      write(*,'(X,A,2X,ES11.4)') &
              'PE transform factor     (pe_trans_fac):',cpp%pe_trans_fac
      write(*,'(X,A,2X,I5)') 'Printout verbosity         (verbosity):',&
                             cpp%verbosity
      write(*,'(X,A,2X,ES11.4)') &
              'ALS regularization        (alspenalty):',cpp%alspenalty
      write(*,'(X,A,4X,A3)') 'ALS solver algorithm  (als_linsys_alg):',&
                             cpp%als_linsys_alg
      write(*,'(X,A,2X,A5)') 'Eigensolver algorithm to use  (solver):',&
                             cpp%solver
      write(*,'(X,A,2X,A5)') 'Calculate spectral bounds (calcbounds):',&
                             cpp%calcbounds
      write(*,'(X,A,2X,I5)') 'Data algorithm for OMP/ACC      (algo):',&
                             cpp%algo
      write(*,'(X,A,2X,I5)') 'Number of solver cycles       (ncycle):',&
                             cpp%ncycle
      write(*,'(X,A,2X,I5)') 'Number of Power iteratons       (npow):',&
                             cpp%npow
      write(*,'(X,A,2X,I5)') 'Low-memory calculation type   (lowmem):',&
                             cpp%lowmem
      write(*,'(X,A,2X,I5)') 'Truncation criterion      (truncation):',&
                             cpp%truncation
      write(*,'(X,A,2X,L5)') 'Use vector updates            (update):',&
                             cpp%update
      write(*,'(X,A,2X,L5)') 'Process top node of tree   (dotopnode):',&
                             cpp%dotopnode
      write(*,'(X,A,2X,ES11.4)') &
                'Solver target energy         (etarget):',cpp%etarget
      write(*,'(X,A,2X,ES11.4)') &
                'Solver convergence criterion (solvtol):',cpp%solvtol
      write(*,'(X,A,2X,A)') 'Restart file name            (resfile):',&
                             cpp%resfile
      write(*,'(X,A,2X,33(I0,X))') &
                            'Random seed                       (rs):',&
                             (cpp%rs(i),i=1,33)
      write(*,'(/X,A)') '*********************************************'

      ENDIF rank0

      end subroutine PrintMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastMLCPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts MLtree to all MPI ranks

      implicit none
      TYPE (CPpar) :: cpp

!     Broadcast variables
      call bcast(cpp%ncycle,mpi_io_rank)
      call bcast(cpp%npow,mpi_io_rank)
      call bcast(cpp%lowmem,mpi_io_rank)
      call bcast(cpp%truncation,mpi_io_rank)
      call bcast(cpp%psirank,mpi_io_rank)
      call bcast(cpp%hrank,mpi_io_rank)
      call bcast(cpp%psinals,mpi_io_rank)
      call bcast(cpp%hnals,mpi_io_rank)
      call bcast(cpp%verbosity,mpi_io_rank)
      call bcast(cpp%alspenalty,mpi_io_rank)
      call bcast(cpp%rs,mpi_io_rank)
      call bcast(cpp%etarget,mpi_io_rank)
      call bcast(cpp%solvtol,mpi_io_rank)
      call bcast(cpp%update,mpi_io_rank)
      call bcast(cpp%dorestart,mpi_io_rank)
      call bcast(cpp%dotopnode,mpi_io_rank)
      call bcast(cpp%algo,mpi_io_rank)
      call bcast(cpp%resfile,mpi_io_rank)
      call bcast(cpp%system,mpi_io_rank)
      call bcast(cpp%pe_transform,mpi_io_rank)
      call bcast(cpp%pe_trans_fac,mpi_io_rank)
      call bcast(cpp%solver,mpi_io_rank)
      call bcast(cpp%calcbounds,mpi_io_rank)
      call bcast(cpp%red2D,mpi_io_rank)
      call bcast(cpp%redND,mpi_io_rank)
      call bcast(cpp%h_sort_alg,mpi_io_rank)
      call bcast(cpp%als_linsys_alg,mpi_io_rank)

      end subroutine BcastMLCPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
