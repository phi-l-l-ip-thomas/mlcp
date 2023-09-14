!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE FFPES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads Force Field potential files and parses data into useful formats
! To add a PES to the code, place all potential constant files in the
! 'pes' directory and add an appropriate call to GetPotential() below

      use ERRORTRAP
      use UTILS
      use SEPDREPN
      use CPCONFIG

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetPotential(V,sys,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Makes a call to the appropriate PES routine or reads potential
! constant files

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in) ::ndof
      character(5), intent(in) :: sys

!     Get the constants for the PES of choice

!     Dummy PES (always evaluates to 0)
      IF (sys(1:5).seq.'dummy') THEN
         call DummyHamiltonian(V,ndof)

!     d-D bi-linearly coupled oscillators
      ELSEIF (sys(1:5).seq.'CpOsc') THEN
         call CoupledOscillatorHamiltonian(V,ndof)

!     (modified) Henon-Heiles, with quartic terms to enforce boundedness
      ELSEIF (sys(1:5).seq.'Henon') THEN
         call HenonHeilesHamiltonian(V,ndof)

!     HO2 (X electronic state) UB3LYP/aug-cc-pVTZ QFF
!     PST (unpublished)
      ELSEIF (sys(1:5).seq.'HO2gs') THEN
         call ReadFFHamiltonian(V,'HO2gs',4,.TRUE.)

!     Formaldehyde QFF (level of theory and source uncertain)
      ELSEIF (sys(1:5).seq.'forma') THEN
         call ReadFFHamiltonian(V,'forma',4,.FALSE.)

!     CH3CN CCSD(T)/cc-pVTZ harmonic + B3LYP/cc-pVTZ cubic/quartic QFF
!     Original: Begue et al, JPCA 109 (2005) 4611.
!     interpreted by: Avila and Carrington, JCP 134 (2011) 054126.
      ELSEIF (sys(1:5).seq.'ch3cn') THEN
         call ReadFFHamiltonian(V,'ch3cn',4,.FALSE.)

!     Vinoxy radical (X electronic state) UB3LYP/aug-cc-pVTZ QFF
!     PST et al, JCP 132 (2010) 114302
      ELSEIF (sys(1:5).seq.'VinOx') THEN
         call ReadFFHamiltonian(V,'VinOx',4,.TRUE.)

!     Vinoxy radical (A electronic state) UB3LYP/aug-cc-pVTZ QFF
!     PST et al, JCP 132 (2010) 114302
      ELSEIF (sys(1:5).seq.'VinOA') THEN
         call ReadFFHamiltonian(V,'VinOA',4,.TRUE.)

!     Ethylene oxide CCSD(T)/cc-pVTZ harmonic + B3LYP/cc-pVTZ
!     cubic/quartic QFF (this version includes the 1 cubic and 7
!     quartic constants omitted from supporting info)
!     Begue et al, JCP 127 (2007) 164115
      ELSEIF (sys(1:5).seq.'EthOx') THEN
         call ReadFFHamiltonian(V,'EthOx',4,.TRUE.)

!     Propargyl peroxy radical (ace-T conformer, X electronic state)
!     B3LYP/cc-pVDZ QFF
!     PST et al, JPCA 114 (2010) 12437
      ELSEIF (sys(1:5).seq.'PglOO') THEN
         call ReadFFHamiltonian(V,'PglOO',4,.TRUE.)

!     Ethyl peroxy radical (G conformer, X electronic state)
!     UB3LYP/aug-cc-pVTZ QFF
!     Melnik et al, JPCA 115 (2011) 13931
      ELSEIF (sys(1:5).seq.'EtPGX') THEN
         call ReadFFHamiltonian(V,'EtPGX',4,.TRUE.)

!     Cyclopentadiene (all proteo) B971/TZ2P QFF
!     Cane and Trombetti, PCCP 11 (2009) 2428
      ELSEIF (sys(1:5).seq.'CPDie') THEN
         call ReadFFHamiltonian(V,'CPDie',4,.TRUE.)

!     Uracil "best estimate" harmonic + MP2/cc-pVTZ cubic/quartic QFF
!     Harmonic terms from Puzzarini et al, JCTC 7 (2011) 3702
!     cubic/quartic terms from Krasnoshchekov et al. JPCA 119 (2015) 6723.
      ELSEIF (sys(1:5).seq.'urSer') THEN
         call ReadFFHamiltonian(V,'urSer',4,.TRUE.)

!     Napthalene-h8 B971/TZ2P QFF
!     Cane et al, JPCA 111 (2007) 8218
      ELSEIF (sys(1:5).seq.'napht') THEN
         call ReadFFHamiltonian(V,'napht',4,.TRUE.)

!     Napthalene B971/TZ2P QFF
!     Mackie et al, JCP 143 (2015) 224314
      ELSEIF (sys(1:5).seq.'nnaph') THEN
         call ReadFFHamiltonian(V,'nnaph',4,.TRUE.)

!     Anthracene B971/TZ2P QFF
!     Mackie et al, JCP 143 (2015) 224314
      ELSEIF (sys(1:5).seq.'anthr') THEN
         call ReadFFHamiltonian(V,'anthr',4,.TRUE.)

!     Tetracene B971/TZ2P QFF
!     Mackie et al, JCP 143 (2015) 224314
      ELSEIF (sys(1:5).seq.'tetra') THEN
         call ReadFFHamiltonian(V,'tetra',4,.TRUE.)

      ELSE
         call AbortWithError('GetPotential(): PES not recognized')
      ENDIF

      IF (ndof.ne.V(1)%nbas(1)) THEN
         write(*,'(X,A,X,I0)') '# DOF from input:',ndof
         write(*,'(X,A,X,I0)') '# DOF from read :',V(1)%nbas(1)
         call AbortWithError('Wrong # DOF for this PES!')
      ENDIF

      end subroutine GetPotential

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DummyHamiltonian(V,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a dummy Hamiltonian which evaluates to zero. For testing.

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in)  :: ndof
      integer, allocatable :: nbas(:)
      integer :: i,j,k,betalen

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

      write(*,'(/X,A,I0/)') 'Potential constants, order: ',1
      call PrintConfigs(V(1))
      write(*,*)

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

      write(*,'(/X,A,I0/)') 'Potential constants, order: ',2
      call PrintConfigs(V(2))
      write(*,*)

      end subroutine CoupledOscillatorHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HenonHeilesHamiltonian(V,ndof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds the (primitive) Henon-Heiles Hamiltonian

      implicit none
      TYPE (Configs), ALLOCATABLE, INTENT(OUT) :: V(:)
      integer, intent(in)  :: ndof
      integer, allocatable :: nbas(:)
      real*8, parameter    :: omega=1.d0 ! Harmonic frequencies
      real*8, parameter    :: beta=0.2   ! Anharmonic constant
      integer :: i,k
      real*8  :: beta2

      write(*,'(X,A)') "--> Setting up Henon-Heiles Hamiltonian"

!     Potential constants (same for all DOF)
      beta2=beta**2/16

!     Allocate configs. Generate a zero config for the linear terms
      ALLOCATE(V(4),nbas(1))
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

      ALLOCATE(nbas(4))
      nbas(:)=ndof
      call NewConfigs(V(4),nbas,3*(ndof-1))
      DEALLOCATE(nbas)

!     Diagonal Hamiltonian elements: omega values
      do i=1,ndof
         V(2)%qns(i,:)=i
         V(2)%coef(i)=0.5*omega
      enddo

!     Cubic part of Hamiltonian (original HH Hamiltonian)
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

!     Quartic part of Hamiltonian (modified HH, to keep PES bound)
      k=1
      do i=1,ndof-1
!        q_i^4 term
         V(4)%qns(k,:)=i
         V(4)%coef(k)=beta2
         k=k+1
!        q_i+1^4 term
         V(4)%qns(k,:)=i+1
         V(4)%coef(k)=beta2
         k=k+1
!        q_i^2*q_i+1^2 term
         V(4)%qns(k,:)=(/i,i,i+1,i+1/)
         V(4)%coef(k)=2*beta2
         k=k+1
      enddo

      DO i=2,4
         write(*,'(/X,A,I0/)') 'Potential constants, order: ',i
         call PrintConfigs(V(i))
      ENDDO
      write(*,*)

      end subroutine HenonHeilesHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ReadFFHamiltonian(W,id,ncp,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads Quartic Force Field Hamiltonian from files containing harmonic,
! cubic, quartic, etc. constants

      implicit none
      TYPE (Configs), ALLOCATABLE,INTENT(OUT) :: W(:)
      integer, intent(in) :: ncp
      logical, intent(in) :: divide
      character(LEN=5), intent(in) :: id
      integer, allocatable :: ncoef(:),qns(:),nbas(:),modpowr(:,:)
      integer :: i,j,k,u,ndof,ndf,InpStat,ReadStat
      real*8  :: ftmp
      character(LEN=20) :: fname
      character*64 :: frmt

      write(*,'(/X,A,A/)')   '--> Reading force field for: ',id
      write(*,'(X,A,X,I0)') 'Max nr of products per term:',ncp

      ALLOCATE(ncoef(ncp),W(ncp))
      ncoef(:)=0
      ndof=0

!     Count the number of potential constants and DOF
      DO k=1,ncp

!        Look for potential file with k coupled DOFs
         write(fname,'(A5,I0,A5,A4)') 'pes/f',k,id,'.dat'
         u=LookForFreeUnit()
         open(u,status='old',file=fname,IOSTAT=InpStat)

!        Next k if file cannot be found
         IF (InpStat /= 0) CYCLE

!        Count the successful potential term reads
         ALLOCATE(qns(k))
         DO
            read(u,*,IOSTAT=ReadStat) (qns(j),j=1,k),ftmp
            IF (ReadStat /= 0) EXIT
            ncoef(k)=ncoef(k)+1
!           For ndof to be determined correctly there must be at least 1
!           potential constant for the last DOF (should be always true)
            ndof=MAX(ndof,MAXVAL(qns))
         ENDDO
         DEALLOCATE(qns)
         close(u)

         IF (ncoef(k).eq.0) CYCLE

         write(*,'(X,2(A,X,I0,X))') &
         'Potential constants of order',k,&
         'read from file :',ncoef(k)
      ENDDO

      write(*,'(X,A,X,I0)') &
            'Number of DOF detected in force constant files:',ndof

!     Read potential constants and store as configurations
      DO k=1,ncp

!        Allocate configs. If ncoefs(k)=0, then generate a zero config
         ALLOCATE(nbas(k))
         nbas(:)=ndof
         call NewConfigs(W(k),nbas,max(1,ncoef(k)))
         DEALLOCATE(nbas)

         IF (ncoef(k).lt.1) CYCLE

         write(*,'(/X,2A,I0/)') 'Potential constants (.dat file',&
                                ' ordering), order: ',k

         write(fname,'(A5,I0,A5,A4)') 'pes/f',k,id,'.dat'
         u=LookForFreeUnit()
         open(u,status='old',file=fname)

         write(frmt,'(A,I0,A)') '(X,',k,'(I3,X),f26.12)'

         DO i=1,ncoef(k)
            read(u,*) (W(k)%qns(i,j),j=1,k),W(k)%coef(i)

!           Divide here if needed to account for degeneracy factors
            IF (divide) THEN
               call DistribModePower(W(k)%qns(i,:),modpowr)
               ndf=SIZE(modpowr,1)
               do j=1,ndf
                  W(k)%coef(i)=W(k)%coef(i)/FACRL(modpowr(j,2))
               enddo
               deallocate(modpowr)
            ENDIF

            write(*,frmt) (W(k)%qns(i,j),j=1,k),W(k)%coef(i)
         ENDDO
         close(u)
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
