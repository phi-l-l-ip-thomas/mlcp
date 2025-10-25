!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module LINSOLVER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Linear solver powered by Alternating Least Squares reduction
! see Beylkin

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE SEPDREPN
      USE MODVECVEC
      USE LINALG
      USE CPMMM
      USE ALSPOW
!!!
      USE REDUCTION
!!!

      implicit none
      real(kind=8), private :: als_penalty=-1.d0
      real(kind=8), allocatable, private :: module_time(:)
      character(len=64), private :: als_solver='uninitialized'
      logical, private :: MODULE_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Init_LinSolver_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      allocate(module_time(mpinodes))
      module_time(:) = 0.d0
      MODULE_SETUP = .TRUE.

      end subroutine Init_LinSolver_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Dispose_LinSolver_Module()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MODULE_SETUP) call Init_LinSolver_Module()
      call Get_MPI_Timings('ALS linear solver',module_time)
      MODULE_SETUP = .FALSE.
      deallocate(module_time)

      end subroutine Dispose_LinSolver_Module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine set_als_settings_LinSolver(penalty,solver)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sets the ALS regularization penalty for the module

      implicit none
      real(kind=8), intent(in) :: penalty
      character(len=64), intent(in) :: solver

      if (penalty.gt.1.d0 .or. penalty.lt.0.d0) then
         write(*,'(A,ES11.4,A)') 'ALS regularization penalty ',&
         penalty,' must be in range: 0 <= penalty <= 1'
         call AbortWithError('set_als_settings_LinSolver(): wrong value')
      endif

      als_penalty=penalty
      als_solver=solver

      end subroutine set_als_settings_LinSolver

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_alg(A,F,G,nitn,ishift,Eshift,which,show)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Select the linear solver algorithm depending on the value of 'lowmem'

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      logical, optional, intent(in) :: show
      integer, intent(in) :: nitn,ishift,which
      real*8, intent(in)  :: Eshift
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_LinSolver_Module()

      call CPU_TIME(ti1)

      IF (which.eq.1) THEN
         call LinSolver_1(A,F,G,nitn,ishift,Eshift,show)
      ELSEIF (which.eq.2) THEN
         call LinSolver_2(A,F,G,nitn,ishift,Eshift,show)
      ELSE
         write(*,*) 'LinSolver_alg(): only which={1,2} implemented'
         call AbortWithError('LinSolver_alg(): invalid algorithm')
      ENDIF

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine LinSolver_alg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_1(A,F,G,nitn,ishift,Eshift,show)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: ATA,ATG,AF
      logical, optional, intent(in) :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,lhs,rhs
      real*8  :: GGprod,AFAFprod,AFGprod,conver,FF,FAF,tmp
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical :: useSVD=.FALSE.

      IF (nitn.eq.0) return

!     For checking convergence later
      GGprod=PRODVV(G)

!     Precompute ATA, ATG, initial AF
      call CPMM(A,ishift,Eshift,.TRUE.,A,ishift,Eshift,.FALSE.,ATA)
      call CPMM(A,ishift,Eshift,.TRUE.,G,0,0.d0,.FALSE.,ATG)
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Set parameters
      rA=SIZE(A%coef)
      IF (ishift.ne.0) rA=rA+1
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=SIZE(G%nbas)

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < G_i^l , (A*F)_i^l' >
      allocate(BB(rAF,rAF),PS(rG,rAF))
      call CONSTPT(AF,0,BB)
      call CONSTPT(AF,G,0,PS)

      IF (present(show)) THEN
         IF (show) THEN
            write(*,*)
            write(*,*) "LinSolver_1() iterations..."
         ENDIF
      ENDIF

      FF=PRODVV(F)
      FAF=PRODVV(F,AF)
      write(*,*) 'ials = ',0,'; RQ = ',(FAF+ishift*Eshift*FF)/FF

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension d
         do d=1,ndof
            m=F%cols(d) ! also equals cols of G
            n=A%cols(d) ! also equals rows of F
            rFn=rF*n

!           Downdate the BB and PS matrices of the linear system
            call UPDATEP(AF,d,BB,.TRUE.)
            call UPDATEP(AF,G,d,PS,.TRUE.)

!           Compute the left-hand sides of the linear system
            allocate(lhs(rFn,rFn))
            lhs=0.d0
            do s=1,rA
               browa=(s-1)*rF ! row group of BB from A^T
               do t=1,rA
                  bcola=(t-1)*rF ! col group of BB from A
                  atark=(s-1)*rA+t ! rank of ATA
                  do k=1,n
                     lrowa=(k-1)*rF ! row group of lhs from basis
                     do l=1,n
                        lcola=(l-1)*rF ! col group of lhs from basis
                        ataind=ATA%ibas(d)-1+(l-1)*n+k ! ATA elem index
                        tmp=ATA%coef(atark)*ATA%base(ataind,atark)
                        do i=1,rF
                           brow=browa+i ! row of BB from AT and F
                           lrow=lrowa+i ! row of lhs from A and F
                           do j=1,rF
                              bcol=bcola+j ! col of BB from A and F
                              lcol=lcola+j ! col of lhs from A and F
                              lhs(lrow,lcol)=&
                              lhs(lrow,lcol)+tmp*BB(brow,bcol)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!           Compute the right-hand sides of the linear system
            allocate(rhs(rFn,m))
            rhs=0.d0
            do s=1,rA
               pcola=(s-1)*rF ! col group of PS from A
               do t=1,rG
                  atgrk=(s-1)*rG+t ! rank of ATG from AT and G
                  do l=1,m
                     rcola=(l-1)*n
                     do k=1,n
                        rrowa=(k-1)*rF ! row group of rhs
                        atgind=ATG%ibas(d)-1+rcola+k ! element of ATG
                        tmp=ATG%coef(atgrk)*ATG%base(atgind,atgrk)
                        do i=1,rF
                           rrow=rrowa+i ! row of rhs
                           pcol=pcola+i ! col of PS
                           rhs(rrow,l)=rhs(rrow,l)+tmp*PS(t,pcol)
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!!! TEST
            IF (itn.eq.nitn) THEN
!               write(*,*) '<AF,AF> matrix, CP-ified:'
!               call truncateCPtest(BB,(/rF,rA/),(/rF,rA/))
!               write(*,*) 'lhs, CP-ified:'
!               call truncateCPtest(lhs,(/rF,n/),(/rF,n/))
            ENDIF
!!! END TEST

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            call SolveLinSys(lhs,rhs,als_penalty,als_solver)

!           Construct improved F
            call UpdateFfromSoln(F,rhs,d)
            deallocate(lhs,rhs)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, we're doomed for now.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'LinSolver_1: NaN on update; itn = ',itn,&
                          '; k = ',d
               call AbortWithError('LinSolver_1 crashed')
            ENDIF

!           Update AF, BB, PS using new Fs
            call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call UPDATEP(AF,d,BB,.FALSE.)
            call UPDATEP(AF,G,d,PS,.FALSE.)

!           Compute ||AF-G|| as a convergence check (always do it since
!           it's cheap compared to building+solving the linear system!)
            AFAFprod=0.d0    
            DO i=1,SIZE(AF%coef)
               DO k=1,SIZE(AF%coef)
                  AFAFprod=AFAFprod+AF%coef(i)*AF%coef(k)*BB(i,k)
               ENDDO
            ENDDO
            AFGprod=0.d0
            DO i=1,SIZE(G%coef)
               DO k=1,SIZE(AF%coef)
                  AFGprod=AFGprod+G%coef(i)*AF%coef(k)*PS(i,k)
               ENDDO
            ENDDO
            conver=sqrt(abs((AFAFprod+GGprod-2*AFGprod)/GGprod))

            IF (present(show)) THEN
               IF (show) &
               write(*,'(2(A,X,I3,X),2(A,X,ES16.8))') &
                          'Itn: ',itn,', d: ',d,', ||AF-G||/||G|| = ',&
                          conver,'; <AF,AF> = ',AFAFprod!PRODVV(F)
            ENDIF
         enddo  ! loop over d
         FF=PRODVV(F)
         FAF=PRODVV(F,AF)
         write(*,*) 'ials = ',itn,'; RQ = ',(FAF+ishift*Eshift*FF)/FF
      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCP(ATA)
      call FlushCP(ATG)
      call FlushCP(AF)

      end subroutine LinSolver_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_2(A,F,G,nitn,ishift,Eshift,show)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver, no normal form

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: AF,Id
      logical, optional, intent(in) :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,lhs,rhs
      real*8  :: GGprod,AFAFprod,AFGprod,conver,tmp
      real*8  :: AFFprod,FGprod,conver2,FF,FAF
      integer :: rA,rAs,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical, allocatable :: sym(:)
      logical :: useSVD=.FALSE.

      IF (nitn.eq.0) return

!     For checking convergence later
      GGprod=PRODVV(G)

!     Precompute initial AF
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Set parameters
      rA=SIZE(A%coef)
      rAs=rA
      IF (ishift.ne.0) rAs=rAs+1
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      rAF=rAs*rF
      ndof=SIZE(G%nbas)

!     Precompute identity matrix for Eshift
      allocate(sym(ndof))
      sym(:)=.FALSE.
      Id=IdentityCPMatrix(A%rows,A%cols,sym)
      call VecScalarMult(Id,-Eshift)
      deallocate(sym)

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < G_i^l , (A*F)_i^l' >
      allocate(BB(rAF,rF),PS(rG,rF))
      call CONSTPT(F,AF,0,BB)
      call CONSTPT(F,G,0,PS)

      IF (present(show)) THEN
         IF (show) THEN
            write(*,*)
            write(*,*) "LinSolver_2() iterations..."
         ENDIF
      ENDIF

      FF=PRODVV(F)
      FAF=PRODVV(F,AF)
      write(*,*) 'ials = ',0,'; RQ = ',(FAF+ishift*Eshift*FF)/FF

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension d
         do d=1,ndof
            m=F%cols(d) ! also equals cols of G
            n=A%cols(d) ! also equals rows of F
            rFn=rF*n

!           Downdate the BB and PS matrices of the linear system
            call UPDATEP(F,AF,d,BB,.TRUE.)
            call UPDATEP(F,G,d,PS,.TRUE.)

!           Compute the left-hand sides of the linear system
            allocate(lhs(rFn,rFn))
            lhs=0.d0
            do t=1,rA ! (does not include Eshift term from Id)
               browa=(t-1)*rF ! row group of BB from A
               do k=1,n
                  lrowa=(k-1)*rF ! row group of lhs from basis
                  do l=1,n
                     lcola=(l-1)*rF ! col group of lhs from basis
                     ataind=A%ibas(d)-1+(l-1)*n+k ! A elem index
                     tmp=A%coef(t)*A%base(ataind,t)
!                     write(*,*) 't = ',t,';A%coef(t) = ',A%coef(t),&
!                     '; A%base(atind,t) = ',A%base(ataind,t)
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
!                        lrow=lrowa+i ! row of lhs from A and F
                        do j=1,rF
                           lrow=lrowa+i
                           lcol=lcola+j
                           lhs(lrow,lcol)=&
                           lhs(lrow,lcol)+tmp*BB(brow,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!           Accumulate energy shift term into LHS (not included above!)
            do t=rA+1,rAs
               browa=(t-1)*rF ! row group of BB from A
               do k=1,n
                  lrowa=(k-1)*rF ! row group of lhs from basis
                  do l=1,n
                     lcola=(l-1)*rF ! col group of lhs from basis
                     ataind=A%ibas(d)-1+(l-1)*n+k ! A elem index
                     tmp=Id%coef(t-rA)*Id%base(ataind,t-rA)
!                     write(*,*) 't = ',t,';Id%coef(1) = ',Id%coef(t-rA),&
!                     '; Id%base(atind,1) = ',Id%base(ataind,t-rA)
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
!                        lrow=lrowa+i ! row of lhs from A and F
                        do j=1,rF
                           lrow=lrowa+i !!! i
                           lcol=lcola+j !!! j
                           lhs(lrow,lcol)=&
                           lhs(lrow,lcol)+tmp*BB(brow,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo

            if (ishift.lt.0) lhs(:,:)=-lhs(:,:)

!           Compute the right-hand sides of the linear system
            allocate(rhs(rFn,m))
            rhs=0.d0
            do t=1,rG
               do l=1,m
                  rcola=(l-1)*n
                  do k=1,n
                     rrowa=(k-1)*rF ! row group of rhs
                     atgind=G%ibas(d)-1+rcola+k
                     tmp=G%coef(t)*G%base(atgind,t)
                     do i=1,rF
                        rrow=rrowa+i ! row of rhs
                        rhs(rrow,l)=rhs(rrow,l)+tmp*PS(t,i)
                     enddo
                  enddo
               enddo
            enddo

!!! TEST
            IF (itn.eq.nitn) THEN
!               write(*,*) '<AF,AF> matrix, CP-ified:'
!               call truncateCPtest(BB,(/rF,rA/),(/rF,rA/))
!               write(*,*) 'lhs, CP-ified:'
!               call truncateCPtest(lhs,(/rF,n/),(/rF,n/))
            ENDIF

!            write(*,*) 'LinSolver_2 B: itn,d = ',itn,d
!            write(*,*) '(',rAF,' x ',rF,')'
!            call PrintMatrix(BB)
!            write(*,*) 'LinSolver_2 P: itn,d = ',itn,d
!            write(*,*) '(',rG,' x ',rF,')'
!            call PrintMatrix(PS)
!            write(*,*) 'LinSolver_2 lhs: itn,d = ',itn,d
!            write(*,*) '(',rFn,' x ',rFn,')'
!            call PrintMatrix(lhs)
!            write(*,*) 'LinSolver_2 rhs: itn,d = ',itn,d
!            write(*,*) '(',rFn,' x ',m,')'
!            call PrintMatrix(rhs)

!!! END TEST

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            call SolveLinSys(lhs,rhs,als_penalty,als_solver)

!           Construct improved F
            call UpdateFfromSoln(F,rhs,d)
            deallocate(lhs,rhs)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, we're doomed for now.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'LinSolver_2: NaN on update; itn = ',itn,&
                          '; k = ',d
               call AbortWithError('LinSolver_2 crashed')
            ENDIF

!           Update AF, BB, PS using new Fs
            call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call UPDATEP(F,AF,d,BB,.FALSE.)
            call UPDATEP(F,G,d,PS,.FALSE.)

!           Compute ||AF-G|| as a convergence check (always do it since
!           it's cheap compared to building+solving the linear system!)
!            AFAFprod=0.d0    
!            DO i=1,SIZE(AF%coef)
!               DO k=1,SIZE(AF%coef)
!                  AFAFprod=AFAFprod+AF%coef(i)*AF%coef(k)*BB(i,k)
!               ENDDO
!            ENDDO
!            AFGprod=0.d0
!            DO i=1,SIZE(G%coef)
!               DO k=1,SIZE(AF%coef)
!                  AFGprod=AFGprod+G%coef(i)*AF%coef(k)*PS(i,k)
!               ENDDO
!            ENDDO

!            AFAFprod=PRODVV(AF,AF) !!! TEMP, for testing
!            AFGprod=PRODVV(AF,G) !!! TEMP, for testing
!            conver=sqrt(abs((AFAFprod+GGprod-2*AFGprod)/GGprod))
            AFFprod=PRODVV(F,AF) !!! TEMP, for testing
            FGprod=PRODVV(F,G) !!! TEMP, for testing
!            conver2=sqrt(abs((AFFprod-2*FGprod)/FGprod))

            IF (present(show)) THEN
               IF (show) &
!               write(*,'(2(A,X,I0),2(A,X,ES16.8))') &
!                          'Itn:',itn,', d:',d,', ||AF-G||/||G|| = ',&
!                          conver,'; <AF,AF> = ',AFAFprod!PRODVV(F)
               write(*,'(2(A,X,I0),2(A,X,ES16.8))') &
                          'Itn:',itn,', d:',d,', <AF,F> = ',&
                         abs(AFFprod),'; <G,F> = ',abs(FGprod)
            ENDIF
         enddo  ! loop over d
         FF=PRODVV(F)
         write(*,*) 'ials = ',itn,'; RQ = ',(AFFprod+ishift*Eshift*FF)/FF
      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCP(Id)
      call FlushCP(AF)

      end subroutine LinSolver_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LintertwinedInvItn_1(A,F,nitn,Eguess)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided intertwined inverse iteration, cause intertwining rocks!

      implicit none
      TYPE (CP), INTENT(IN)    :: A
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: ATA,AF
      real*8, intent(in)   :: Eguess
      integer, intent(in)  :: nitn
!!!   TEST
      integer :: irep,nrep
!!!
      real*8, dimension (:,:), allocatable :: BB,PS,lhs,rhs
      real*8 :: FFprod,AFFprod,rqold,rqnew,tmp
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,pcola,pcol,atgrk,atgind,ataind
      real(kind=8) :: ti1,ti2

      IF (.NOT. MODULE_SETUP) call Init_LinSolver_Module()

      call CPU_TIME(ti1)

      IF (nitn.eq.0) return

!     Rayleigh quotient for initial E-shift
      rqold=Eguess

!     Precompute ATA, AF
      call CPMM(A,1,rqold,.TRUE.,A,1,rqold,.FALSE.,ATA)
      call CPMM(A,1,rqold,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Set parameters
      rA=SIZE(A%coef)+1 ! +1 due to Eshift
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=SIZE(F%nbas)

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < F_i^l , (A*F)_i^l' >
      allocate(BB(rAF,rAF),PS(rF,rAF))
      call CONSTPT(AF,0,BB)
      call CONSTPT(AF,F,0,PS)

!!! TEST: check RQ before calculation to make sure it is what I expect!
            AFFprod=0.d0
            FFprod=0.d0
            DO i=1,rF
!              Calc <F,(A-E1)F>, excluding the Eshift (last rF) terms
               DO k=1,rAF-rF
                  AFFprod=AFFprod+F%coef(i)*AF%coef(k)*PS(i,k)
               ENDDO
!              Calc -<F,-F> from the Eshift terms using F's coefs
               DO k=1,rF
                  FFprod=FFprod-F%coef(i)*F%coef(k)*PS(i,rAF-rF+k)
               ENDDO
            ENDDO
            rqnew=AFFprod/FFprod
            write(*,*) 'On entry, energy from previous iteration is = ',rqold
            write(*,*) 'Itn: ',0,', d: ',0,', <F,AF>/<F,F> = ',rqnew,&
                       '; <F,F> = ',FFprod                       
!!!

!!! Repetitions on a single mode
      nrep=1
!!!

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension d
         do d=1,ndof
!            do irep=1,nrep !!! Reps for a single DOF

            n=A%cols(d) ! also equals rows of F
            rFn=rF*n

!           Downdate the BB and PS matrices of the linear system
            call UPDATEP(AF,d,BB,.TRUE.)
            call UPDATEP(AF,F,d,PS,.TRUE.)

!           Compute the left-hand sides of the linear system
            allocate(lhs(rFn,rFn))
            lhs=0.d0
            do s=1,rA
               browa=(s-1)*rF ! row group of BB from A^T
               do t=1,rA
                  bcola=(t-1)*rF ! col group of BB from A
                  atark=(s-1)*rA+t ! rank of ATA
                  do k=1,n
                     lrowa=(k-1)*rF ! row group of lhs from basis
                     do l=1,n
                        lcola=(l-1)*rF ! col group of lhs from basis
                        ataind=ATA%ibas(d)-1+(l-1)*n+k ! ATA elem index
                        tmp=ATA%coef(atark)*ATA%base(ataind,atark)
                        do i=1,rF
                           brow=browa+i ! row of BB from AT and F
                           lrow=lrowa+i ! row of lhs from A and F
                           do j=1,rF
                              bcol=bcola+j ! col of BB from A and F
                              lcol=lcola+j ! col of lhs from A and F
                              lhs(lrow,lcol)=&
                              lhs(lrow,lcol)+tmp*BB(brow,bcol)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!           Compute the right-hand sides of the linear system
            allocate(rhs(rFn,1))
            rhs=0.d0
            do s=1,rA
               pcola=(s-1)*rF ! col group of PS from A
               do t=1,rF
                  atgrk=(s-1)*rF+t ! rank of ATG from AT and F
                  do k=1,n
                     rrowa=(k-1)*rF ! row group of rhs
                     atgind=AF%ibas(d)-1+k ! element of ATF
                     tmp=AF%coef(atgrk)*AF%base(atgind,atgrk)
                     do i=1,rF
                        rrow=rrowa+i ! row of rhs
                        pcol=pcola+i ! col of PS
                        rhs(rrow,1)=rhs(rrow,1)+tmp*PS(t,pcol)
                     enddo
                  enddo
               enddo
            enddo
            
!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            call SolveLinSys(lhs,rhs,als_penalty,als_solver)

!           Construct improved F
            call UpdateFfromSoln(F,rhs,d)
            deallocate(lhs,rhs)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'LintertwinedInvItn_1: ',&
                             'NaN on update; itn = ',itn,'; k = ',d
                  call AbortWithError('LintertwinedInvItn_1 crashed')
               ENDIF
            ENDDO

!!!         Normalize F to avoid underflow
            FFprod=0.d0
            DO i=1,rF
               DO k=1,rF
                  FFprod=FFprod-F%coef(i)*F%coef(k)*PS(i,rAF-rF+k)*&
                         dot_product(F%base(F%ibas(d):F%fbas(d),i),&
                                     F%base(F%ibas(d):F%fbas(d),k))
               ENDDO
            ENDDO
            FFPROD=1.d0/sqrt(abs(FFPROD))
!            write(*,*) 'Norm factor = ',FFPROD
            F%coef(:)=FFPROD*F%coef(:)
!            call NORMALIZE(F)

!           Update AF, Rayleigh Quotient, BB, PS using new Fs
!!! Watch out: the call below uses the old Rayleigh Quotient
            call CPMM(A,1,rqold,.FALSE.,F,0,0.d0,.FALSE.,AF,d)

!!! Watch out here also: if RQ changes sign, part of PS corresponding to
! Eshift will have wrong sign
            call UPDATEP(AF,d,BB,.FALSE.)
            call UPDATEP(AF,F,d,PS,.FALSE.)

!           Rayleigh quotient update
            AFFprod=0.d0
            FFprod=0.d0
            tmp=0.d0 !!!
            DO i=1,rF
!              Calc <F,(A-E1)F>, excluding the Eshift (last rF) terms
               DO k=1,rAF-rF
                  AFFprod=AFFprod+F%coef(i)*AF%coef(k)*PS(i,k)
               ENDDO
               DO k=1,rF
                  FFprod=FFprod-F%coef(i)*F%coef(k)*PS(i,rAF-rF+k)
               ENDDO
               DO k=1,rF
                  tmp=tmp+F%coef(i)*AF%coef(rAF-rF+k)*PS(i,rAF-rF+k)
               ENDDO
            ENDDO
            rqnew=AFFprod/FFprod

!           Sign change if needed
            IF ((rqold*rqnew).lt.0.d0) THEN
               write(*,*)
               write(*,*) "!!! Changing signs on PS and BB !!!"
               write(*,*)
               PS(:,rAF-rF+1:rAF)=-PS(:,rAF-rF+1:rAF)
               BB(rAF-rF+1:rAF,1:rAF-rF)=-BB(rAF-rF+1:rAF,1:rAF-rF)
               BB(1:rAF-rF,rAF-rF+1:rAF)=-BB(1:rAF-rF,rAF-rF+1:rAF)
            ENDIF
!!!
            write(*,'(X,A,I3,X,A,I2,X,3(A,f16.8))') &
            'Itn: ',itn,', d: ',d,', <F,AF>/<F,F> = ',rqnew,&
            '; <F,(A-EI)F>/<F,F> = ',(AFFprod+tmp)/FFprod,&
            '; <F,F> = ',FFprod
!!!
!           Update the coefs in ATA, AF with the new Rayleigh Quotient
!            call CPMMreshift(A,1,rqold,rqnew,A,1,rqold,rqnew,ATA)
!            call CPMMreshift(A,1,rqold,rqnew,F,0,rqold,rqnew,AF)
!            rqold=rqnew

!            enddo ! loop over irep
         enddo  ! loop over d

!        Alt: Update shift only after loop over d
!         call CPMMreshift(A,1,rqold,rqnew,A,1,rqold,rqnew,ATA)
!         call CPMMreshift(A,1,rqold,rqnew,F,0,rqold,rqnew,AF)
!         rqold=rqnew



      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCP(ATA)
      call FlushCP(AF)

      call CPU_TIME(ti2)
      module_time=module_time+ti2-ti1

      end subroutine LintertwinedInvItn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine truncateCPtest(M,rows,cols)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided intertwined inverse iteration, cause intertwining rocks!

      implicit none
      TYPE (CP) :: MCP
      real*8, intent(in)  :: M(:,:)
      integer, intent(in) :: rows(:),cols(:)
      real*8, allocatable :: Mtr(:,:)
      integer :: rk

      write(*,*) 'SVD untruncated'
      call PrintSVD(M)

      write(*,*) 'Matrix-to-CP...'
      MCP=Matrix2CP(M,rows,cols)
      call PrintCPmat(MCP,.FALSE.)
      rk=MAX(MCP%R()/2,1)
      write(*,*) 'truncating matrix at rank: ',rk,'/',MCP%R()
      call ResizeV(MCP,rk)
      write(*,*) 'truncated MCP:'
!      call MCP%print()
      call PrintCPmat(MCP,.FALSE.)

      write(*,*) 'MCP back to matrix...'
      Mtr=CP2Matrix(MCP)
      write(*,*) 'SVD truncated'
      call PrintSVD(Mtr)

      call FlushCP(MCP)
      DEALLOCATE(Mtr)

      end subroutine truncateCPtest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module LINSOLVER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
