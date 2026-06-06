!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
       
      implicit none

      TYPE CPpar
           integer :: layer,mode,ncycle,npow,lowmem,algo,verbosity
           integer :: max_nmode,max_sum,max_qn
           integer :: psirank,hrank,psinals,hnals,nactivations
           real(kind=8) :: Etarget,solvtol,ovrlpenalty,alspenalty,hcut
           real(kind=8) :: padbounds
           logical :: update,diag,donode,reduceHQ
           character(len=128) :: solver,calcbounds,orthogalg
           character(len=128) :: red2D,redND
           character(len=128) :: h_sort_alg,activation
           character(len=128) :: als_linsys_alg
           ! 'fieldlist' holds parameters that can be set in input file
           character(len=128), dimension(32) :: &
           fieldlist=[character(len=128) :: &
                       'layer-mode',&
                       'donode',&
                       'max_nmode',&  
                       'max_sum',&  
                       'max_qn',&  
                       'Etarget',&  
                       'activation',&
                       'nactivations',&
                       'hcut',&
                       'verbosity',&
                       'algo',&
                       'lowmem',&
                       'h_sort_alg',&
                       'calcbounds',&
                       'padbounds',&
                       'solver',&
                       'solvtol',&
                       'ncycle',&
                       'npow',&
                       'orthogalg',&
                       'diag',&
                       'reduceHQ',&
                       'update',&
                       'ovrlpenalty',&
                       'red2D',&
                       'redND',&
                       'psirank',&
                       'psinals',&
                       'hrank',&
                       'hnals',&
                       'alspenalty',&
                       'als_linsys_alg'&
                      ]
      END TYPE CPpar

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetCPPDefaults(cpp,layer,mode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set MLCP parameter defaults

      implicit none
      TYPE (CPpar) :: cpp
      integer, intent(in) :: layer,mode

!     Integer parameters
      cpp%layer=layer
      cpp%mode=mode
      cpp%verbosity=0
      cpp%algo=-1
      cpp%lowmem=2
      cpp%max_nmode=-1
      cpp%max_sum=-1
      cpp%max_qn=-1
      cpp%nactivations=1
      cpp%ncycle=10
      cpp%npow=10
      cpp%psirank=10
      cpp%psinals=10
      cpp%hrank=0
      cpp%hnals=200

!     Real parameters
      cpp%Etarget=0.d0
      cpp%hcut=0.d0
      cpp%padbounds=0.d0
      cpp%solvtol=1.d-10
      cpp%ovrlpenalty=0.d0
      cpp%alspenalty=1.d-10

!     Character parameters
      cpp%calcbounds='calc'
      cpp%h_sort_alg='pack'
      cpp%activation='scale-linear'
      cpp%solver='powr'
      cpp%orthogalg='gram'
      cpp%red2D='SVD'
      cpp%redND='ALS'
      cpp%als_linsys_alg='LU'

!     Logical parameters
      cpp%diag=.true.
      cpp%reduceHQ=.true.
      cpp%update=.true.
      cpp%donode=.true.

      end subroutine SetCPPDefaults

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
         fi=1 !-1
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
!                     write(*,'(3A,3(X,I0))') '"',line,'"',fi,ff,m !!!
                     write(fields(nfields+fct),'(A)') line(fi:ff)
                     write(vals(nfields+fct),'(A)') line(qi+1:qf-1)
                  endif
                  fi=1 !-1
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

      subroutine processcppfields(cpp,fields,values)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assigns values to CPpar variables if found in input file

      implicit none
      TYPE (CPpar), intent(inout) :: cpp
      character(len=*), dimension(:), intent(in) :: fields, values
      character(len=128) :: thevalue
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
             call AbortWithError('processcppfields(): bad input item')
         endif 
      enddo

!     Check for items in field list; overwrite default if found

!     Logical fields
      call get_field(fields,values,'diag',thevalue,fcount)
      if (fcount.eq.1) cpp%diag=string2logical(thevalue)

      call get_field(fields,values,'reduceHQ',thevalue,fcount)
      if (fcount.eq.1) cpp%reduceHQ=string2logical(thevalue)

      call get_field(fields,values,'update',thevalue,fcount)
      if (fcount.eq.1) cpp%update=string2logical(thevalue)

      call get_field(fields,values,'donode',thevalue,fcount)
      if (fcount.eq.1) cpp%donode=string2logical(thevalue)

!     Integer fields
      call get_field(fields,values,'verbosity',thevalue,fcount)
      if (fcount.eq.1) cpp%verbosity=string2integer(thevalue)

      call get_field(fields,values,'algo',thevalue,fcount)
      if (fcount.eq.1) cpp%algo=string2integer(thevalue)

      call get_field(fields,values,'lowmem',thevalue,fcount)
      if (fcount.eq.1) cpp%lowmem=string2integer(thevalue)

      call get_field(fields,values,'max_nmode',thevalue,fcount)
      if (fcount.eq.1) cpp%max_nmode=string2integer(thevalue)

      call get_field(fields,values,'max_sum',thevalue,fcount)
      if (fcount.eq.1) cpp%max_sum=string2integer(thevalue)

      call get_field(fields,values,'max_qn',thevalue,fcount)
      if (fcount.eq.1) cpp%max_qn=string2integer(thevalue)

      call get_field(fields,values,'nactivations',thevalue,fcount)
      if (fcount.eq.1) cpp%nactivations=string2integer(thevalue)

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

!     Real fields
      call get_field(fields,values,'Etarget',thevalue,fcount)
      if (fcount.eq.1) cpp%Etarget=string2real8(thevalue)

      call get_field(fields,values,'hcut',thevalue,fcount)
      if (fcount.eq.1) cpp%hcut=string2real8(thevalue)

      call get_field(fields,values,'padbounds',thevalue,fcount)
      if (fcount.eq.1) cpp%padbounds=string2real8(thevalue)

      call get_field(fields,values,'solvtol',thevalue,fcount)
      if (fcount.eq.1) cpp%solvtol=string2real8(thevalue)

      call get_field(fields,values,'ovrlpenalty',thevalue,fcount)
      if (fcount.eq.1) cpp%ovrlpenalty=string2real8(thevalue)

      call get_field(fields,values,'alspenalty',thevalue,fcount)
      if (fcount.eq.1) cpp%alspenalty=string2real8(thevalue)

!     String fields
      call get_field(fields,values,'layer-mode',thevalue,fcount)
      if (fcount.eq.1) call parselayermode(thevalue,cpp%layer,cpp%mode)

      call get_field(fields,values,'calcbounds',thevalue,fcount)
      if (fcount.eq.1) cpp%calcbounds=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'h_sort_alg',thevalue,fcount)
      if (fcount.eq.1) cpp%h_sort_alg=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'activation',thevalue,fcount)
      if (fcount.eq.1) cpp%activation=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'solver',thevalue,fcount)
      if (fcount.eq.1) cpp%solver=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'orthogalg',thevalue,fcount)
      if (fcount.eq.1) cpp%orthogalg=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'red2D',thevalue,fcount)
      if (fcount.eq.1) cpp%red2D=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'redND',thevalue,fcount)
      if (fcount.eq.1) cpp%redND=TRIM(ADJUSTL(thevalue))

      call get_field(fields,values,'als_linsys_alg',thevalue,fcount)
      if (fcount.eq.1) cpp%als_linsys_alg=TRIM(ADJUSTL(thevalue))

      end subroutine processcppfields

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
                    '" must appear [0,1]x, appears ',fcount,'x'
         call AbortWithError('get_field(): duplicated field')
      endif

      end subroutine get_field

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine WriteCPPInputs(dpp,cpp,line)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Writes parameters in 'cpp' differing from 'dpp' to line

      implicit none
      TYPE (CPpar), intent(in) :: dpp,cpp
      character(len=1024), intent(out) :: line
      character(len=128) :: tag

      line=""

      IF (cpp%layer.ne.dpp%layer .or. cpp%mode.ne.dpp%mode) THEN
         write(tag,'(A,I0,A,I0,A)') &
         "layer-mode='",cpp%layer,"-",cpp%mode,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%donode.neqv.dpp%donode) THEN
         write(tag,'(A,L0,A)') &
         "donode='",cpp%donode,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%max_nmode.ne.dpp%max_nmode) THEN
         write(tag,'(A,I0,A)') &
         "max_nmode='",cpp%max_nmode,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%max_sum.ne.dpp%max_sum) THEN
         write(tag,'(A,I0,A)') &
         "max_sum='",cpp%max_sum,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%max_qn.ne.dpp%max_qn) THEN
         write(tag,'(A,I0,A)') &
         "max_qn='",cpp%max_qn,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%Etarget.ne.dpp%Etarget) THEN
         write(tag,'(A,ES14.6,A)') &
         "Etarget='",cpp%Etarget,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%activation.seq.dpp%activation)) THEN
         write(tag,'(A,A,A)') &
         "activation='",TRIM(ADJUSTL(cpp%activation)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%nactivations.ne.dpp%nactivations) THEN
         write(tag,'(A,I0,A)') &
         "nactivations='",cpp%nactivations,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%hcut.ne.dpp%hcut) THEN
         write(tag,'(A,ES14.6,A)') &
         "hcut='",cpp%hcut,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%verbosity.ne.dpp%verbosity) THEN
         write(tag,'(A,I0,A)') &
         "verbosity='",cpp%verbosity,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%algo.ne.dpp%algo) THEN
         write(tag,'(A,I0,A)') &
         "algo='",cpp%algo,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%lowmem.ne.dpp%lowmem) THEN
         write(tag,'(A,I0,A)') &
         "lowmem='",cpp%lowmem,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%h_sort_alg.seq.dpp%h_sort_alg)) THEN
         write(tag,'(A,A,A)') &
         "h_sort_alg='",TRIM(ADJUSTL(cpp%h_sort_alg)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%calcbounds.seq.dpp%calcbounds)) THEN
         write(tag,'(A,A,A)') &
         "calcbounds='",TRIM(ADJUSTL(cpp%calcbounds)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%padbounds.ne.dpp%padbounds) THEN
         write(tag,'(A,ES14.6,A)') &
         "padbounds='",cpp%padbounds,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%solver.seq.dpp%solver)) THEN
         write(tag,'(A,A,A)') &
         "solver='",TRIM(ADJUSTL(cpp%solver)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%solvtol.ne.dpp%solvtol) THEN
         write(tag,'(A,ES14.6,A)') &
         "solvtol='",cpp%solvtol,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%ncycle.ne.dpp%ncycle) THEN
         write(tag,'(A,I0,A)') &
         "ncycle='",cpp%ncycle,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%npow.ne.dpp%npow) THEN
         write(tag,'(A,I0,A)') &
         "npow='",cpp%npow,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%orthogalg.seq.dpp%orthogalg)) THEN
         write(tag,'(A,A,A)') &
         "orthogalg='",TRIM(ADJUSTL(cpp%orthogalg)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%diag.neqv.dpp%diag) THEN
         write(tag,'(A,L0,A)') &
         "diag='",cpp%diag,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%reduceHQ.neqv.dpp%reduceHQ) THEN
         write(tag,'(A,L0,A)') &
         "reduceHQ='",cpp%reduceHQ,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%update.neqv.dpp%update) THEN
         write(tag,'(A,L0,A)') &
         "update='",cpp%update,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%ovrlpenalty.ne.dpp%ovrlpenalty) THEN
         write(tag,'(A,ES14.6,A)') &
         "ovrlpenalty='",cpp%ovrlpenalty,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%red2D.seq.dpp%red2D)) THEN
         write(tag,'(A,A,A)') &
         "red2D='",TRIM(ADJUSTL(cpp%red2D)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%redND.seq.dpp%redND)) THEN
         write(tag,'(A,A,A)') &
         "redND='",TRIM(ADJUSTL(cpp%redND)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%psirank.ne.dpp%psirank) THEN
         write(tag,'(A,I0,A)') &
         "psirank='",cpp%psirank,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%psinals.ne.dpp%psinals) THEN
         write(tag,'(A,I0,A)') &
         "psinals='",cpp%psinals,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%hrank.ne.dpp%hrank) THEN
         write(tag,'(A,I0,A)') &
         "hrank='",cpp%hrank,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%hnals.ne.dpp%hnals) THEN
         write(tag,'(A,I0,A)') &
         "hnals='",cpp%hnals,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (cpp%alspenalty.ne.dpp%alspenalty) THEN
         write(tag,'(A,ES14.6,A)') &
         "alspenalty='",cpp%alspenalty,"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      IF (.not.(cpp%als_linsys_alg.seq.dpp%als_linsys_alg)) THEN
         write(tag,'(A,A,A)') &
         "als_linsys_alg='",TRIM(ADJUSTL(cpp%als_linsys_alg)),"'"
         line=TRIM(line) // " " // TRIM(tag)
      ENDIF

      end subroutine WriteCPPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintCPPInputs(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints out parameters read in CP.inp

      implicit none
      TYPE (CPpar),INTENT(IN) :: cpp
      integer :: i

      rank0 : IF (mpirank.eq.mpi_prnt_rank) THEN

      write(*,'(X,A,X,I5)') 'Guess max nr coupled modes (max_nmode):',&
                             cpp%max_nmode
      write(*,'(X,A,X,I5)') 'Guess max sum of quanta      (max_sum):',&
                             cpp%max_sum
      write(*,'(X,A,X,I5)') 'Guess max quantum number      (max_qn):',&
                             cpp%max_qn
      write(*,'(X,A,X,ES11.4)') &
                             'Guess target energy          (Etarget):',&
                             cpp%Etarget
      write(*,'(X,A,2X,A)') 'Solver activation algo    (activation):',&
                             TRIM(ADJUSTL(cpp%orthogalg))
      write(*,'(X,A,X,I5)') 'Nr. solver activations  (nactivations):',&
                             cpp%nactivations
      write(*,'(X,A,X,ES11.4)') &
                             'Hamiltonian term minimum cutoff (hcut):',&
                             cpp%hcut
      write(*,'(X,A,X,I5)') 'Printout verbosity         (verbosity):',&
                             cpp%verbosity
      write(*,'(X,A,X,I5)') 'Data algorithm for OMP/ACC      (algo):',&
                             cpp%algo
      write(*,'(X,A,X,I5)') 'Low-memory calculation type   (lowmem):',&
                             cpp%lowmem
      write(*,'(X,A,2X,A)') 'H sorting algorithm       (h_sort_alg):',&
                             TRIM(ADJUSTL(cpp%h_sort_alg))
      write(*,'(X,A,2X,A)') 'Calculate spectral bounds (calcbounds):',&
                             TRIM(ADJUSTL(cpp%calcbounds))
      write(*,'(X,A,X,ES11.4)') &
                'Spectral bounds pad factor (padbounds):',cpp%padbounds
      write(*,'(X,A,2X,A)') 'Eigensolver algorithm to use  (solver):',&
                             TRIM(ADJUSTL(cpp%solver))
      write(*,'(X,A,X,ES11.4)') &
                'Solver convergence criterion (solvtol):',cpp%solvtol
      write(*,'(X,A,X,I5)') 'Number of solver cycles       (ncycle):',&
                             cpp%ncycle
      write(*,'(X,A,X,I5)') 'Number of Power iteratons       (npow):',&
                             cpp%npow
      write(*,'(X,A,2X,A)') 'Orthogonalization algo.    (orthogalg):',&
                             TRIM(ADJUSTL(cpp%orthogalg))
      write(*,'(X,A,X,L5)') 'Use subspace diagonalization    (diag):',&
                             cpp%diag
      write(*,'(X,A,X,L5)') 'Rank-reduce HQ during diag. (reduceHQ):',&
                             cpp%reduceHQ
      write(*,'(X,A,X,L5)') 'Use vector updates            (update):',&
                             cpp%update
      write(*,'(X,A,X,ES11.4)') &
              'Gen. eigv regularization (ovrlpenalty):',cpp%ovrlpenalty
      write(*,'(X,A,4X,A)') 'Reduction type for 2-D modes   (red2D):',&
                             TRIM(ADJUSTL(cpp%red2D))
      write(*,'(X,A,4X,A)') 'Reduction type for >2-D modes  (redND):',&
                             TRIM(ADJUSTL(cpp%redND))
      write(*,'(X,A,X,I5)') 'Wavefunction reduced rank    (psirank):',&
                             cpp%psirank
      write(*,'(X,A,X,I5)') 'Number of ALS iterations-w.f.(psinals):',&
                             cpp%psinals
      write(*,'(X,A,X,I5)') 'Hamiltonian reduced rank       (hrank):',&
                             cpp%hrank
      write(*,'(X,A,X,I5)') 'Number of ALS iterations-H     (hnals):',&
                             cpp%hnals
      write(*,'(X,A,X,ES11.4)') &
              'ALS regularization        (alspenalty):',cpp%alspenalty
      write(*,'(X,A,2X,A)') 'ALS solver algorithm  (als_linsys_alg):',&
                             TRIM(ADJUSTL(cpp%als_linsys_alg))

      ENDIF rank0

      end subroutine PrintCPPInputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine compareCPP(cpp,cppo,same1,same2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares two cpp objects

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp,cppo
      logical, intent(out) :: same1,same2

!     The following parameters can only be changed on restart if
!     processing for the node has not yet begun
      same1=.TRUE.
      same1=(same1.and.(cpp%donode.eqv.cppo%donode))
      same1=(same1.and.(cpp%max_nmode.eq.cppo%max_nmode))
      same1=(same1.and.(cpp%max_sum.eq.cppo%max_sum))
      same1=(same1.and.(cpp%max_qn.eq.cppo%max_qn))
      same1=(same1.and.(cpp%Etarget.eq.cppo%Etarget))
      same1=(same1.and.(cpp%activation.eq.cppo%activation))
      same1=(same1.and.(cpp%nactivations.eq.cppo%nactivations))

!     The following parameters can be changed on restart even if
!     processing for the node has begun
      same2=.TRUE.
      same2=(same2.and.(cpp%hcut.eq.cppo%hcut))
      same2=(same2.and.(cpp%verbosity.eq.cppo%verbosity))
      same2=(same2.and.(cpp%algo.eq.cppo%algo))
      same2=(same2.and.(cpp%lowmem.eq.cppo%lowmem))
      same2=(same2.and.(cpp%h_sort_alg.eq.cppo%h_sort_alg))
      same2=(same2.and.(cpp%calcbounds.eq.cppo%calcbounds))
      same2=(same2.and.(cpp%padbounds.eq.cppo%padbounds))
      same2=(same2.and.(cpp%solver.eq.cppo%solver))
      same2=(same2.and.(cpp%solvtol.eq.cppo%solvtol))
      same2=(same2.and.(cpp%ncycle.eq.cppo%ncycle))
      same2=(same2.and.(cpp%npow.eq.cppo%npow))
      same2=(same2.and.(cpp%orthogalg.eq.cppo%orthogalg))
      same2=(same2.and.(cpp%diag.eqv.cppo%diag))
      same2=(same2.and.(cpp%reduceHQ.eqv.cppo%reduceHQ))
      same2=(same2.and.(cpp%update.eqv.cppo%update))
      same2=(same2.and.(cpp%ovrlpenalty.eq.cppo%ovrlpenalty))
      same2=(same2.and.(cpp%red2D.eq.cppo%red2D))
      same2=(same2.and.(cpp%redND.eq.cppo%redND))
      same2=(same2.and.(cpp%psirank.eq.cppo%psirank))
      same2=(same2.and.(cpp%psinals.eq.cppo%psinals))
      same2=(same2.and.(cpp%hrank.eq.cppo%hrank))
      same2=(same2.and.(cpp%hnals.eq.cppo%hnals))
      same2=(same2.and.(cpp%alspenalty.eq.cppo%alspenalty))
      same2=(same2.and.(cpp%als_linsys_alg.eq.cppo%als_linsys_alg))

      end subroutine compareCPP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine showCPPcomparisons(cpp,cppo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Summarize changes in current CPPar object vs. old one

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp,cppo
      logical :: same1,same2

      call compareCPP(cpp,cppo,same1,same2)
      if (same1.and.same2) return

!     Print changes to job parameters
      IF (mpirank.eq.mpi_prnt_rank) THEN
         IF (cpp%donode.neqv.cppo%donode) write(*,2232) &
            '  * donode changed from ',cppo%donode,' to ',cpp%donode
         IF (cpp%max_nmode.ne.cppo%max_nmode) write(*,2233) &
            '  * max_nmode changed from ',&
            cppo%max_nmode,' to ',cpp%max_nmode
         IF (cpp%max_sum.ne.cppo%max_sum) write(*,2233) &
            '  * max_sum changed from ',cppo%max_sum,' to ',cpp%max_sum
         IF (cpp%max_qn.ne.cppo%max_qn) write(*,2233) &
            '  * max_qn changed from ',cppo%max_qn,' to ',cpp%max_qn
         IF (cpp%Etarget.ne.cppo%Etarget) write(*,2234) &
            '  * Etarget changed from ',cppo%Etarget,&
            ' to ',cpp%Etarget
         IF (.not.(cpp%activation .seq. cppo%activation)) &
            write(*,2235) '  * activation changed from "',&
            TRIM(ADJUSTL(cppo%activation)),'" to "',&
            TRIM(ADJUSTL(cpp%activation)),'"'
         IF (cpp%nactivations.ne.cppo%nactivations) write(*,2233) &
            '  * nactivations changed from ',&
            cppo%nactivations,' to ',cpp%nactivations
         IF (cpp%hcut.ne.cppo%hcut) write(*,2234) &
            '  * hcut changed from ',cppo%hcut,&
            ' to ',cpp%hcut
         IF (cpp%verbosity.ne.cppo%verbosity) write(*,2233) &
            '  * verbosity changed from ',&
            cppo%verbosity,' to ',cpp%verbosity
         IF (cpp%algo.ne.cppo%algo) write(*,2233) &
            '  * algo changed from ',cppo%algo,' to ',cpp%algo
         IF (cpp%lowmem.ne.cppo%lowmem) write(*,2233) &
            '  * lowmem changed from ',cppo%lowmem,' to ',cpp%lowmem
         IF (.not.(cpp%h_sort_alg .seq. cppo%h_sort_alg)) &
            write(*,2235) '  * h_sort_alg changed from "',&
            TRIM(ADJUSTL(cppo%h_sort_alg)),'" to "',&
            TRIM(ADJUSTL(cpp%h_sort_alg)),'"'
         IF (.not.(cpp%calcbounds .seq. cppo%calcbounds)) &
            write(*,2235) '  * calcbounds changed from "',&
            TRIM(ADJUSTL(cppo%calcbounds)),'" to "',&
            TRIM(ADJUSTL(cpp%calcbounds)),'"'
         IF (cpp%padbounds.ne.cppo%padbounds) write(*,2234) &
            '  * padbounds changed from ',cppo%padbounds,&
            ' to ',cpp%padbounds
         IF (.not.(cpp%solver .seq. cppo%solver)) write(*,2235) &
            '  * solver changed from "',&
            TRIM(ADJUSTL(cppo%solver)),'" to "',&
            TRIM(ADJUSTL(cpp%solver)),'"'
         IF (cpp%solvtol.ne.cppo%solvtol) write(*,2234) &
            '  * solvtol changed from ',cppo%solvtol,' to ',cpp%solvtol
         IF (cpp%ncycle.ne.cppo%ncycle) write(*,2233) &
            '  * ncycle changed from ',cppo%ncycle,' to ',cpp%ncycle
         IF (cpp%npow.ne.cppo%npow) write(*,2233) &
            '  * npow changed from ',cppo%npow,' to ',cpp%npow
         IF (.not.(cpp%orthogalg .seq. cppo%orthogalg)) &
            write(*,2235) '  * orthogalg changed from "',&
            TRIM(ADJUSTL(cppo%orthogalg)),'" to "',&
            TRIM(ADJUSTL(cpp%orthogalg)),'"'
         IF (cpp%diag.neqv.cppo%diag) write(*,2232) &
            '  * diag changed from ',cppo%diag,' to ',cpp%diag
         IF (cpp%reduceHQ.neqv.cppo%reduceHQ) write(*,2232) &
            '  * reduceHQ changed from ',cppo%reduceHQ,' to ',cpp%reduceHQ
         IF (cpp%update.neqv.cppo%update) write(*,2232) &
            '  * update changed from ',cppo%update,' to ',cpp%update
         IF (cpp%ovrlpenalty.ne.cppo%ovrlpenalty) write(*,2234) &
            '  * ovrlpenalty changed from ',cppo%ovrlpenalty,&
            ' to ',cpp%ovrlpenalty
         IF (.not.(cpp%red2D .seq. cppo%red2D)) write(*,2235) &
            '  * red2D changed from "',&
            TRIM(ADJUSTL(cppo%red2D)),'" to "',&
            TRIM(ADJUSTL(cpp%red2D)),'"'
         IF (.not.(cpp%redND .seq. cppo%redND)) write(*,2235) &
            '  * redND changed from "',&
            TRIM(ADJUSTL(cppo%redND)),'" to "',&
            TRIM(ADJUSTL(cpp%redND)),'"'
         IF (cpp%psirank.ne.cppo%psirank) write(*,2233) &
            '  * psirank changed from ',cppo%psirank,' to ',cpp%psirank
         IF (cpp%psinals.ne.cppo%psinals) write(*,2233) &
            '  * psinals changed from ',cppo%psinals,' to ',cpp%psinals
         IF (cpp%hrank.ne.cppo%hrank) write(*,2233) &
            '  * hrank changed from ',cppo%hrank,' to ',cpp%hrank
         IF (cpp%hnals.ne.cppo%hnals) write(*,2233) &
            '  * hnals changed from ',cppo%hnals,' to ',cpp%hnals
         IF (cpp%alspenalty.ne.cppo%alspenalty) write(*,2234) &
            '  * alspenalty changed from ',cppo%alspenalty,&
            ' to ',cpp%alspenalty
         IF (.not.(cpp%als_linsys_alg .seq. cppo%als_linsys_alg)) &
            write(*,2235) '  * als_linsys_alg changed from "',&
            TRIM(ADJUSTL(cppo%als_linsys_alg)),'" to "',&
            TRIM(ADJUSTL(cpp%als_linsys_alg)),'"'
         write(*,*)
      ENDIF

2232  format(X,A,L0,A,L0)
2233  format(X,A,I0,A,I0)
2234  format(X,A,ES11.4,A,ES11.4)
2235  format(X,5A)

      end subroutine showCPPcomparisons

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BcastCPP(cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Broadcasts MLtree to all MPI ranks

      implicit none
      TYPE (CPpar) :: cpp

!     Broadcast variables
      call bcast(cpp%layer,mpi_io_rank)
      call bcast(cpp%mode,mpi_io_rank)
      call bcast(cpp%donode,mpi_io_rank)
      call bcast(cpp%max_nmode,mpi_io_rank)
      call bcast(cpp%max_sum,mpi_io_rank)
      call bcast(cpp%max_qn,mpi_io_rank)
      call bcast(cpp%Etarget,mpi_io_rank)
      call bcast(cpp%activation,mpi_io_rank)
      call bcast(cpp%nactivations,mpi_io_rank)
      call bcast(cpp%hcut,mpi_io_rank)
      call bcast(cpp%verbosity,mpi_io_rank)
      call bcast(cpp%algo,mpi_io_rank)
      call bcast(cpp%lowmem,mpi_io_rank)
      call bcast(cpp%h_sort_alg,mpi_io_rank)
      call bcast(cpp%calcbounds,mpi_io_rank)
      call bcast(cpp%padbounds,mpi_io_rank)
      call bcast(cpp%solver,mpi_io_rank)
      call bcast(cpp%solvtol,mpi_io_rank)
      call bcast(cpp%ncycle,mpi_io_rank)
      call bcast(cpp%npow,mpi_io_rank)
      call bcast(cpp%orthogalg,mpi_io_rank)
      call bcast(cpp%diag,mpi_io_rank)
      call bcast(cpp%reduceHQ,mpi_io_rank)
      call bcast(cpp%update,mpi_io_rank)
      call bcast(cpp%ovrlpenalty,mpi_io_rank)
      call bcast(cpp%red2D,mpi_io_rank)
      call bcast(cpp%redND,mpi_io_rank)
      call bcast(cpp%psirank,mpi_io_rank)
      call bcast(cpp%psinals,mpi_io_rank)
      call bcast(cpp%hrank,mpi_io_rank)
      call bcast(cpp%hnals,mpi_io_rank)
      call bcast(cpp%alspenalty,mpi_io_rank)
      call bcast(cpp%als_linsys_alg,mpi_io_rank)

      end subroutine BcastCPP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine parselayermode(tag,il,im)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Extracts layer and mode numbers from tag

      implicit none
      character(len=128), intent(in) :: tag
      integer, intent(out) :: il,im
      integer :: j,lf,mi

      il=-1
      im=-1

      do j=1,128
         lf=j-1
         mi=j+1
         if (tag(j:j).eq.'-') exit
      enddo
      
      if (lf.gt.0 .and. mi.le.128) then
         il=string2integer(tag(1:lf))
         im=string2integer(tag(mi:128))
      endif

      end subroutine parselayermode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE INPUTCP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
