!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE FFPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads Force Field potential files and parses data into useful formats
! To add a PES to the code, place all potential constant files in the
! 'pes' directory and add an appropriate call to GetPotential() below

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE SEPDREPN
      USE CPCONFIG

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPotential(V,sys,ndof,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Makes a call to the appropriate PES routine or reads potential
! constant files

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in) :: ndof,verbosity
      character(len=*), intent(in) :: sys

!     Get the constants for the PES of choice

!     Dummy PES (always evaluates to 0)
      IF (trim(adjustl(sys)).seq.'dummy') THEN
         call DummyHamiltonian(V,ndof)

!     d-D bi-linearly coupled oscillators
      ELSEIF (trim(adjustl(sys)).seq.'CpOsc') THEN
         call CoupledOscillatorHamiltonian(V,ndof)

!     (modified) Henon-Heiles, with quartic terms to enforce boundedness
      ELSEIF (trim(adjustl(sys)).seq.'Henon') THEN
         call HenonHeilesHamiltonian(V,ndof)

!     Formaldehyde QFF (level of theory and source uncertain)
      ELSEIF (trim(adjustl(sys)).seq.'forma') THEN
         call ReadFFHamiltonian(V,sys,.FALSE.)

!     CH3CN CCSD(T)/cc-pVTZ harmonic + B3LYP/cc-pVTZ cubic/quartic QFF
!     Original: Begue et al, JPCA 109 (2005) 4611.
!     interpreted by: Avila and Carrington, JCP 134 (2011) 054126.
      ELSEIF (trim(adjustl(sys)).seq.'ch3cn') THEN
         call ReadFFHamiltonian(V,sys,.FALSE.)

!     Arbitrary QFF, e.g. Gaussian format
      ELSE
         call ReadFFHamiltonian(V,sys,.TRUE.)
      ENDIF

      call PrintPotentialConstants(V,verbosity)

      IF (ndof.ne.V(1)%nbas(1)) THEN
         write(*,'(X,A,X,I0)') '# DOF from input:',ndof
         write(*,'(X,A,X,I0)') '# DOF from read :',V(1)%nbas(1)
         call AbortWithError('Wrong # DOF for this PES!')
      ENDIF

      end subroutine GetPotential

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PrintPotentialConstants(V,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints potential constants read from file

      implicit none
      TYPE (Configs), INTENT(IN) :: V(:)
      integer, intent(in) :: verbosity
      integer :: i,j,k,ncoef,ncp
      character*64 :: frmt

!     Print out potential constants
      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.1) THEN
         ncp=SIZE(V)
         DO k=1,ncp
            ncoef=SIZE(V(k)%coef)
            IF (ncoef.gt.1 .or. ANY(V(k)%qns(1,:).gt.0)) THEN
               write(*,'(/X,A,I0/)') 'Potential constants, order: ',k
               write(frmt,'(A,I0,A)') '(X,',k,'(I3,X),f26.12)'
               DO i=1,ncoef
                  write(*,frmt) (V(k)%qns(i,j),j=1,k),V(k)%coef(i)
               ENDDO
            ENDIF
         ENDDO
         write(*,*)
      ENDIF

      end subroutine PrintPotentialConstants

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DummyHamiltonian(V,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a dummy Hamiltonian which evaluates to zero. For testing.

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in)  :: ndof
      integer, allocatable :: nbas(:)
      integer :: i,j,k,betalen

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A)') "--> Setting up Dummy Hamiltonian"

!     Allocate configs. This PES has two linear terms which cancel one
!     another. This is to prevent FF2Configs() from detecting a single
!     zero entry which would cause it to exit without generating a PES 
!     config array
      ALLOCATE(V(1),nbas(1))
      nbas(:)=ndof
      call NewConfigs(V(1),nbas,2)
      V(1)%qns(:,1)=1
      V(1)%coef(1)=1.d0
      V(1)%coef(2)=-1.d0
      DEALLOCATE(nbas)

      end subroutine DummyHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CoupledOscillatorHamiltonian(V,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds the (primitive) coupled oscillator Hamiltonian

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in)  :: ndof
      integer, allocatable :: nbas(:)
      real*8, parameter    :: beta=0.1 ! Bilinear coupling constant
      integer :: i,j,k,betalen

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A)') "--> Setting up Coupled Oscillator Hamiltonian"

!     Size of anharmonic term array
      betalen=0
      DO i=ndof-1,1,-1
         betalen=betalen+i
      ENDDO

!     Allocate configs. Generate a zero config for the linear terms
      ALLOCATE(V(2),nbas(1))
      nbas(:)=ndof
      call NewConfigs(V(1),nbas,1)
      DEALLOCATE(nbas)

      ALLOCATE(nbas(2))
      nbas(:)=ndof
      call NewConfigs(V(2),nbas,ndof+betalen)
      DEALLOCATE(nbas)

!     Diagonal Hamiltonian elements: omega values
      k=1
      do i=1,ndof
         V(2)%qns(k,:)=i
         V(2)%coef(k)=0.5*sqrt(i*0.5d0)
         k=k+1
      enddo

!     Off-diagonal Hamiltonian elements: beta values
      do i=1,ndof-1
         do j=i+1,ndof
            V(2)%qns(k,:)=(/i,j/)
            V(2)%coef(k)=beta
            k=k+1
         enddo
      enddo

      end subroutine CoupledOscillatorHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HenonHeilesHamiltonian(V,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds the (primitive) Henon-Heiles Hamiltonian

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in)  :: ndof
      integer, allocatable :: nbas(:)
      real(kind=8) :: omega,beta
      integer :: i,k

      IF (mpirank.eq.mpi_prnt_rank) &
      write(*,'(X,A)') "--> Setting up Henon-Heiles Hamiltonian"

!     Potential constants (same for all DOF)
      omega=1.d0             ! Harmonic frequency
      beta=1.d0/sqrt(80.d0)  ! Anharmonic constant

!     Allocate configs. Generate a zero config for the linear terms
      ALLOCATE(V(3),nbas(1))
      nbas(:)=ndof
      call NewConfigs(V(1),nbas,1)
      DEALLOCATE(nbas)

      ALLOCATE(nbas(2))
      nbas(:)=ndof
      call NewConfigs(V(2),nbas,ndof)
      DEALLOCATE(nbas)

      ALLOCATE(nbas(3))
      nbas(:)=ndof
      call NewConfigs(V(3),nbas,2*(ndof-1))
      DEALLOCATE(nbas)

!     Diagonal Hamiltonian elements: omega values
      do i=1,ndof
         V(2)%qns(i,:)=i
         V(2)%coef(i)=0.5*omega
      enddo

!     Cubic part of Hamiltonian
      k=1
      do i=1,ndof-1
!        q_i^2*q_i+1 term
         V(3)%qns(k,:)=(/i,i,i+1/)
         V(3)%coef(k)=beta
         k=k+1
!        q_i+1^3 term
         V(3)%qns(k,:)=i+1
         V(3)%coef(k)=-beta/3
         k=k+1
      enddo

      end subroutine HenonHeilesHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadFFHamiltonian(W,id,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads Quartic Force Field Hamiltonian from files containing harmonic,
! cubic, quartic, etc. constants

      implicit none
      TYPE (Configs), ALLOCATABLE,INTENT(OUT) :: W(:)
      logical, intent(in) :: divide
      character(len=*), intent(in) :: id
      integer, allocatable :: ncoef(:),qns(:),nbas(:),modpowr(:,:)
      integer :: i,j,k,u,ndof,ndf,InpStat,ReadStat
      integer, parameter :: ncp=99 ! Max number of coupled DOF
      real*8  :: ftmp
      character(len=128) :: fname
      character*64 :: frmt

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,'(/X,A,A)') '--> Reading force field for: ',&
                              trim(adjustl(id))
      ENDIF

      ALLOCATE(ncoef(ncp),W(ncp))
      ncoef(:)=0
      ndof=0

      IF (mpirank.eq.mpi_io_rank) THEN

!        Count the number of potential constants and DOF
         DO k=1,ncp

!           Look for potential file with k coupled DOFs
            write(fname,'(A,I0,A,A)') 'pes/f',k,&
                                        trim(adjustl(id)),'.dat'
            u=LookForFreeUnit()
            open(u,status='old',file=trim(adjustl(fname)),IOSTAT=InpStat)

!           Next k if file cannot be found
            IF (InpStat /= 0) CYCLE

!           Count the successful potential term reads
            ALLOCATE(qns(k))
            DO
               read(u,*,IOSTAT=ReadStat) (qns(j),j=1,k),ftmp
               IF (ReadStat /= 0) EXIT
               ncoef(k)=ncoef(k)+1
!              For ndof to be determined correctly there must be at least 1
!              potential constant for the last DOF (should be always true)
               ndof=MAX(ndof,MAXVAL(qns))
            ENDDO
            DEALLOCATE(qns)
            close(u)
         ENDDO

      ENDIF

      call bcast(ncoef,mpi_io_rank)
      call bcast(ndof,mpi_io_rank)

      IF (mpirank.eq.mpi_prnt_rank) THEN
          DO k=1,ncp
             IF (ncoef(k).eq.0) CYCLE
             write(*,'(5X,2(A,X,I0,X))') &
            'Potential constants of order',k,&
            'read from file :',ncoef(k)
          ENDDO
          write(*,'(5X,A,X,I0/)') &
         'Number of DOF detected in force constant files:',ndof
      ENDIF

!     Read potential constants and store as configurations
      DO k=1,ncp

!        Allocate configs. If ncoefs(k)=0, then generate a zero config
         ALLOCATE(nbas(k))
         nbas(:)=ndof
         call NewConfigs(W(k),nbas,max(1,ncoef(k)))
         DEALLOCATE(nbas)

         IF (ncoef(k).lt.1) CYCLE

         IF (mpirank.eq.mpi_io_rank) THEN

            write(fname,'(A,I0,A,A)') 'pes/f',k,&
                                        trim(adjustl(id)),'.dat'
            u=LookForFreeUnit()
            open(u,status='old',file=trim(adjustl(fname)))

            DO i=1,ncoef(k)
               read(u,*) (W(k)%qns(i,j),j=1,k),W(k)%coef(i)

!              Divide here if needed to account for degeneracy factors
               IF (divide) THEN
                  call DistribModePower(W(k)%qns(i,:),modpowr)
                  ndf=SIZE(modpowr,1)
                  do j=1,ndf
                     W(k)%coef(i)=W(k)%coef(i)/FACRL(modpowr(j,2))
                  enddo
                  deallocate(modpowr)
               ENDIF

            ENDDO
            close(u)
         ENDIF ! rnk0

         call bcast(W(k)%qns,mpi_io_rank)
         call bcast(W(k)%coef,mpi_io_rank)
      ENDDO

      end subroutine ReadFFHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FF2Configs(V,C)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds configuration representation of PES in terms of powers from
! configuration rep'n based on products of x_i

      implicit none
      TYPE (Configs), INTENT(IN)  :: V(:)
      TYPE (Configs), INTENT(OUT) :: C
      integer, allocatable :: nbas(:),modpowr(:,:)
      integer :: i,j,k,l,ncp,ndof,nrk,ndf

!     Set parameters
      ncp=SIZE(V)
      ndof=V(1)%nbas(1)

!     Determine the rank of the PES
      nrk=0
      DO k=1,ncp
!        Skip if there are no terms for this k
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE
         nrk=nrk+SIZE(V(k)%coef)
      ENDDO

!     The values in nbas are temporarily set to 1. These will be
!     corrected after the configs are generated
      ALLOCATE(nbas(ndof))
      nbas(:)=1
      call NewConfigs(C,nbas,nrk)
      DEALLOCATE(nbas)
      C%qns(:,:)=1

      l=1
      DO k=1,ncp

!        If V(k) is a zero vector, skip
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

         DO i=1,SIZE(V(k)%coef)

!           Convert config to product of q(s) to some power(s)
            call DistribModePower(V(k)%qns(i,:),modpowr)
            ndf=SIZE(modpowr,1)

            DO j=1,ndf
               C%qns(l,modpowr(j,1))=C%qns(l,modpowr(j,1))+modpowr(j,2)
            ENDDO
            C%coef(l)=V(k)%coef(i)
            l=l+1

            DEALLOCATE(modpowr)
         ENDDO
      ENDDO

!     Set the values of nbas to the max value for each DOF
      DO j=1,ndof
         C%nbas(j)=MAXVAL(C%qns(:,j))
      ENDDO

      end subroutine FF2Configs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ExtractOmegas(V,omega,verbosity)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Transforms potential into "Morsified" coordinates with asymptotic long
! range behavior

      implicit none
      TYPE (Configs), INTENT(INOUT) :: V(:)
      real(kind=8), allocatable, intent(out) :: omega(:)
      integer, intent(in)  :: verbosity
      integer, allocatable :: modpowr(:,:)
      integer :: i,ndof,ndf,mode

!     Set parameters
      ndof=V(1)%nbas(1)

      ALLOCATE(omega(ndof))
      omega=0.d0

!     Extract the harmonic constant from the quadratic terms
      DO i=1,SIZE(V(2)%coef)
         call DistribModePower(V(2)%qns(i,:),modpowr)
         ndf=SIZE(modpowr,1)
         IF (ndf.eq.1) THEN  ! 1D potential term
            mode=modpowr(1,1)
            omega(mode)=omega(mode)+V(2)%coef(i)
         ENDIF
         deallocate(modpowr)
      ENDDO

      IF (mpirank.eq.mpi_prnt_rank .and. verbosity.ge.1) THEN
         write(*,*) 'Harmonic constants extracted from PES:',&
                    '(used to construct KEO)'
         write(*,*)
         write(*,'(X,A,9X,A))') 'DOF','Omega'
         DO i=1,ndof
            write(*,'(X,I3,X,f22.12)') i,2*omega(i)
         ENDDO
         write(*,*)
      ENDIF

      end subroutine ExtractOmegas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine TransformPES(V,vtype,alpha,trans,afac,opmap,optable)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Transforms potential into "Morsified" coordinates with asymptotic long
! range behavior

      implicit none
      TYPE (Configs), INTENT(INOUT) :: V(:)
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: vtype(:)
      character(len=64), intent(in) :: trans
      real(kind=8), intent(in) :: afac
      real(kind=8), allocatable, intent(out) :: alpha(:)
      integer, intent(in)    :: opmap(:)
      logical, intent(inout) :: optable(:,:,:)
      integer, allocatable :: modpowr(:,:)
      integer :: i,j,k,l,ndof,ndf,ncoup,ot
      real(kind=8), allocatable  :: v1d(:,:)
      logical, allocatable :: sympes(:)
      real(kind=8) :: am1

!     Set parameters
      ncoup=SIZE(V)
      ndof=V(1)%nbas(1)

!     Error checking
      IF ((trans .seq. 'poly-tanh') .or. &
          (trans .seq. 'morse-tanh')) THEN
         IF (ncoup.lt.4) THEN
            write(*,'(A,I0,2A)') "Max order: ",ncoup," for this PES; ",&
            "must be >= 4 to do 'poly-tanh' or 'morse-tanh' transform"
         ENDIF
      ENDIF

      ALLOCATE(v1d(ndof,ncoup),alpha(ndof),sympes(ndof),vtype(ncoup))
      alpha(:)=1.d0
      v1d(:,:)=0.d0
      sympes(:)=.FALSE.

!     Loop over quadratic, cubic, quartic, ... terms in the PES
      DO k=1,ncoup
!        Array for holding transformation types
!        (note: NewConfigs() initializes this to zero)
         call NewConfigs(vtype(k),V(k)%nbas,SIZE(V(k)%coef))

!        If V(k) is a zero vector, skip
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

!        Extract the 1D potentials if transformation is to be made
         IF (.not.(trans .seq. 'none')) THEN
            DO i=1,SIZE(V(k)%coef)
               call DistribModePower(V(k)%qns(i,:),modpowr)
               ndf=SIZE(modpowr,1)
               IF (ndf.eq.1) THEN  ! 1D potential term
                  v1d(modpowr(1,1),modpowr(1,2))=&
                  v1d(modpowr(1,1),modpowr(1,2))+V(k)%coef(i)
               ENDIF
               deallocate(modpowr)
            ENDDO
         ENDIF
      ENDDO

      IF (.not.(trans .seq. 'none')) THEN

!        Compute the value of the alpha parameter and transform the 1D PES
         DO i=1,ndof
!           Symmetric potential: compute alpha for coupling terms only
            IF (v1d(i,3).eq.0.d0) THEN
               sympes(i)=.TRUE.
               alpha(i)=-1.5d0*v1d(i,4)/v1d(i,2)
               alpha(i)=SIGN(sqrt(abs(alpha(i))),alpha(i))

!           Asymmetric potential: compute alpha and morsify 1D terms
            ELSE
               alpha(i)=-v1d(i,3)/v1d(i,2)
               am1=1.d0/alpha(i)
!              Morse series reversal defs
               v1d(i,4) = am1*(am1*(am1*(am1*v1d(i,4) + 1.5d0*v1d(i,3)) + &
                              (11.d0/12.d0)*v1d(i,2)) + 0.25d0*v1d(i,1))
!              Alpha is chosen to give v1d(i,3) = 0. Enforce v1d(i,3) = 0
!              here since the expression below may give a nonzero alpha
!              due to roundoff error
!               v1d(i,3) = am1*(am1*(am1*v1d(i,3) + v1d(i,2)) + &
!                              (1.d0/3.d0)*v1d(i,1))
               v1d(i,3) = 0.d0
               v1d(i,2) = am1*(am1*v1d(i,2) + 0.5d0*v1d(i,1))
               v1d(i,1) = am1*v1d(i,1)
            ENDIF
            alpha(i)=abs(afac*alpha(i)) ! Scaled alpha
         ENDDO

      ENDIF

      IF (mpirank.eq.mpi_prnt_rank .and. &
          (.not.(trans .seq. 'none'))) THEN
         write(*,'(/X,A,A/)') '--> The PES will be transformed ',&
              'into asymptotically-decaying coordinates'
         IF (trans .seq. 'poly-tanh') THEN 
            write(*,'(X,A,A)') 'Asymmetric 1D potentials :',&
                               ' y_i = q_i'
            write(*,'(X,A,A)') ' Symmetric 1D potentials :',&
                               ' y_i = q_i' 
            write(*,'(X,A,A)') ' d-D coupling potentials :',&
                               ' y_i = tanh(alpha_i*q_i)'
         ELSEIF (trans .seq. 'morse-tanh') THEN
            write(*,'(X,A,A)') 'Asymmetric 1D potentials :',&
                               ' y_i = 1-exp(-alpha_i*q_i)'
            write(*,'(X,A,A)') ' Symmetric 1D potentials :',&
                               ' y_i = q_i' 
            write(*,'(X,A,A)') ' d-D coupling potentials :',&
                               ' y_i = tanh(alpha_i*q_i)'
         ENDIF
         write(*,'(/X,A,f10.6,A)') &
                           'DOF Sym Alpha-values (scaled by ',afac,')'
         DO i=1,ndof
            write(*,'(X,I3,2X,L1,2X,ES15.8)') i,sympes(i),alpha(i)
         ENDDO
         write(*,*)
      ENDIF

!     Check which PEO are present and transform the PES, if requested
      DO k=1,ncoup

!        If V(k) is a zero vector, skip
         IF (SIZE(V(k)%coef).eq.1 .and. V(k)%coef(1).eq.0.d0) CYCLE

         DO i=1,SIZE(V(k)%coef)
            call DistribModePower(V(k)%qns(i,:),modpowr)
            ndf=SIZE(modpowr,1)

!           1D potential terms:
            IF (ndf.eq.1) THEN

!              Symmetric 1D potential or no transformation
               IF (sympes(modpowr(1,1)).or.(trans .seq. 'none')) THEN
                  ot=findival(opmap,0) ! Leave as power of q
!              Asymmetric 1D potential, poly-tanh
               ELSEIF (trans .seq. 'poly-tanh') THEN
                  ot=findival(opmap,0) ! Leave as power of q
!              Asymmetric 1D potential, morse-tanh
               ELSEIF (trans .seq. 'morse-tanh') THEN
                  V(k)%coef(i)=v1d(modpowr(1,1),modpowr(1,2))
                  ot=findival(opmap,2) ! morse
               ENDIF

            ELSE
!              Coupling terms
               IF (trans .seq. 'none') THEN
                  ot=findival(opmap,0) ! Leave as power of q
               ELSE
                  ot=findival(opmap,1) ! tanh
                  DO j=1,ndf
                     V(k)%coef(i)=V(k)%coef(i)/&
                                  alpha(modpowr(j,1))**modpowr(j,2)
                  ENDDO
               ENDIF

            ENDIF

!           Record primitive operator as present
            DO j=1,ndf
               IF (.not.optable(modpowr(j,1),modpowr(j,2),ot)) &
                  optable(modpowr(j,1),modpowr(j,2),ot)=.TRUE.
            ENDDO

            vtype(k)%qns(i,1:ndf)=ot
            deallocate(modpowr)
         ENDDO
      ENDDO

      DEALLOCATE(v1d,sympes)

      end subroutine TransformPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DistribModePower(modlist,modpowr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts Quartic Force Field entries into list of modes and powers
! E.g. '1 1 7 8' = q1 * q1 * q7 * q8 --> q1^2 * q7^1 * q8^1

      implicit none
      integer, intent(in)  :: modlist(:)
      integer, allocatable, intent(out):: modpowr(:,:)
      integer, allocatable :: modt(:,:)
      integer :: lm,i,j,ndof

      lm=SIZE(modlist)
      ALLOCATE(modt(lm,2))
      modt=0
      modt(1,1)=modlist(1)
      modt(1,2)=1
      ndof=1

      DO i=2,lm
         DO j=1,ndof
!           If the DOF is repeated, increment its power
            IF (modlist(i).eq.modt(j,1)) THEN
               modt(j,2)=modt(j,2)+1
               EXIT
            ENDIF
!           If the DOF is unique, add to list
            IF (j.eq.ndof) THEN
               ndof=ndof+1
               modt(ndof,1)=modlist(i)
               modt(ndof,2)=1
            ENDIF
         ENDDO
      ENDDO

!     Copy modt to modpowr
      ALLOCATE(modpowr(ndof,2))
      modpowr(1:ndof,:)=modt(1:ndof,:)
      DEALLOCATE(modt)

      end subroutine DistribModePower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE FFPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
