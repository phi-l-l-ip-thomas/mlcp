!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module LINSOLVER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Linear solver powered by Alternating Least Squares reduction
! see Beylkin

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE MODVECVEC
      USE LINALG
      USE CPMMM
      USE ALSPOW
!!!
      USE REDUCTION
!!!

      implicit none
      real*8, private  :: linsolver_time=0.d0
      logical, private :: LINSOLVER_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupLinSolver()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Set up linear solver module

      implicit none

      linsolver_time=0.d0
      LINSOLVER_SETUP=.TRUE.

      end subroutine SetupLinSolver

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeLinSolver()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Dispose linear solver module

      implicit none

      IF (.NOT. LINSOLVER_SETUP) call SetupLinSolver()

!     Set up the module if it was not set up already
      LINSOLVER_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total reduction time (LINSOLVER)  (s)',&
                            linsolver_time

      end subroutine DisposeLinSolver

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_alg(A,F,G,nitn,ishift,Eshift,lowmem,show)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Select the linear solver algorithm depending on the value of 'lowmem'

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      logical, optional, intent(in) :: show
      integer, intent(in) :: nitn,ishift,lowmem
      real*8, intent(in)  :: Eshift


      IF (.NOT. LINSOLVER_SETUP) call SetupLinSolver()

      IF (lowmem.ne.1) THEN
         write(*,*) 'LinSolver_alg(): only lowmem=1 implemented so far'
         call AbortWithError('LinSolver_alg(): invalid lowmem')
      ELSE
         call LinSolver_1(A,F,G,nitn,ishift,Eshift,show)
      ENDIF

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
      real*8  :: valpen,tmp,t1,t2
      real*8  :: GGprod,AFAFprod,AFGprod,conver
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical :: useSVD=.TRUE.

      IF (nitn.eq.0) return

!     For checking convergence later
      GGprod=PRODVV(G)

!     Precompute ATA, ATG, initial AF
      call CPMM(A,ishift,Eshift,.TRUE.,A,ishift,Eshift,.FALSE.,ATA)
      call CPMM(A,ishift,Eshift,.TRUE.,G,0,0.d0,.FALSE.,ATG)
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

      call CPU_TIME(t1)

!     Set parameters
      rA=SIZE(A%coef)
      IF (ishift.ne.0) rA=rA+1
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=SIZE(G%nbas)

!     Penalty to avoid bad conditioning
!      valpen=maxval(F%coef)*1.d-10
      valpen=sqrt(2.d0)*EPSILON(1.d0)

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
            IF (useSVD) THEN
               call SolveLinSysSVD(lhs,rhs,valpen)
            ELSE
               call SolveLinSysLU(lhs,rhs,valpen)
            ENDIF

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
            call CPU_TIME(t2)
            linsolver_time=linsolver_time+t2-t1
            call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call CPU_TIME(t1)
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
               write(*,*) 'Itn: ',itn,', d: ',d,', ||AF-G||/||G|| = ',&
                          conver,'; <F,F> = ',PRODVV(F)
            ENDIF
         enddo  ! loop over d
      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCP(ATA)
      call FlushCP(ATG)
      call FlushCP(AF)

      call CPU_TIME(t2)
      linsolver_time=linsolver_time+t2-t1

      end subroutine LinSolver_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_2(A,F,G,nitn,ishift,Eshift,show)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver, no normal form

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: AF
      logical, optional, intent(in) :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,lhs,rhs
      real*8  :: valpen,tmp,t1,t2
      real*8  :: GGprod,AFAFprod,AFGprod,conver
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical :: useSVD=.TRUE.

      IF (nitn.eq.0) return

!     For checking convergence later
      GGprod=PRODVV(G)

!     Precompute initial AF
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

      call CPU_TIME(t1)

!     Set parameters
      rA=SIZE(A%coef)
      IF (ishift.ne.0) rA=rA+1
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=SIZE(G%nbas)

!     Penalty to avoid bad conditioning
!      valpen=maxval(F%coef)*1.d-10
      valpen=sqrt(2.d0)*EPSILON(1.d0)

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
            do t=1,rA
               browa=(t-1)*rF ! row group of BB from A
               do k=1,n
                  lrowa=(k-1)*rF ! row group of lhs from basis
                  do l=1,n
                     lcola=(l-1)*rF ! col group of lhs from basis
                     ataind=A%ibas(d)-1+(l-1)*n+k ! A elem index
                     tmp=A%coef(t)*A%base(ataind,t)
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
                        lrow=lrowa+i ! row of lhs from A and F
                        do j=1,rF
                           lcol=lcola+j ! col of lhs from A and F
                           lhs(lrow,lcol)=&
                           lhs(lrow,lcol)+tmp*BB(brow,j)
                        enddo
                     enddo
                  enddo
               enddo
            enddo

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
!!! END TEST

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            IF (useSVD) THEN
               call SolveLinSysSVD(lhs,rhs,valpen)
            ELSE
               call SolveLinSysLU(lhs,rhs,valpen)
            ENDIF

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
            call CPU_TIME(t2)
            linsolver_time=linsolver_time+t2-t1
            call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call CPU_TIME(t1)
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
            AFAFprod=PRODVV(AF) !!! TEMP, for testing
            AFGprod=PRODVV(AF,G) !!! TEMP, for testing
            conver=sqrt(abs((AFAFprod+GGprod-2*AFGprod)/GGprod))

            IF (present(show)) THEN
               IF (show) &
               write(*,*) 'Itn: ',itn,', d: ',d,', ||AF-G||/||G|| = ',&
                          conver,'; <F,F> = ',PRODVV(F)
            ENDIF
         enddo  ! loop over d
      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCP(AF)

      call CPU_TIME(t2)
      linsolver_time=linsolver_time+t2-t1

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
      real*8  :: valpen,tmp,t1,t2
      real*8 :: FFprod,AFFprod,rqold,rqnew
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,pcola,pcol,atgrk,atgind,ataind

      IF (nitn.eq.0) return

!     Rayleigh quotient for initial E-shift
      rqold=Eguess

!     Precompute ATA, AF
      call CPMM(A,1,rqold,.TRUE.,A,1,rqold,.FALSE.,ATA)
      call CPMM(A,1,rqold,.FALSE.,F,0,0.d0,.FALSE.,AF)

      call CPU_TIME(t1)

!     Set parameters
      rA=SIZE(A%coef)+1 ! +1 due to Eshift
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=SIZE(F%nbas)

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

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
            call SolveLinSysLU(lhs,rhs,valpen)
!!! TEST: use SVD instead
!            call SolveLinSysSVD(lhs,rhs,valpen)
!!!

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
            call CPU_TIME(t2)
            linsolver_time=linsolver_time+t2-t1

!!! Watch out: the call below uses the old Rayleigh Quotient
            call CPMM(A,1,rqold,.FALSE.,F,0,0.d0,.FALSE.,AF,d)

            call CPU_TIME(t1)

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

      call CPU_TIME(t2)
      linsolver_time=linsolver_time+t2-t1

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
