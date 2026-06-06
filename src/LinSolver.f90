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
!      USE BLOCKUTILS
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

      subroutine LinSolver_1(A,F,G,nitn,ishift,Eshift,showconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: ATA,ATG,AF
      logical, optional, intent(in) :: showconv
      logical :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,B1,B2,PS,lhs,rhs
      real*8  :: GG,AFAF,AFG,FF,FAF,RQ,conver,tmp
      integer :: rA,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind

      IF (nitn.eq.0) return

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

      show=.FALSE.
      IF (present(showconv).and.mpirank.eq.mpi_prnt_rank) show=showconv

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < G_i^l , (A*F)_i^l' >
      allocate(BB(rAF,rAF),PS(rG,rAF),B1(rF,rAF),B2(rF,rF))
      call CONSTPT(AF,0,BB)
      call CONSTPT(AF,G,0,PS)
      call CONSTPT(AF,F,0,B1)
      call CONSTPT(F,0,B2)

      IF (show) THEN
         write(*,*)
         write(*,*) "LinSolver_1() iterations..."
      ENDIF

      AFAF=PVVfromPS(AF,AF,BB)
      AFG=PVVfromPS(AF,G,PS)
      FAF=PVVfromPS(AF,F,B1)
      FF=PVVfromPS(F,F,B2)
      GG=PRODVV(G)

      RQ=(FAF+ishift*Eshift*FF)/FF
      conver=sqrt(abs((AFAF+GG-2*AFG)/GG))

      IF (show) THEN
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,A)') &
        'ials:',0,'d:',0,'Es = ',Eshift,'||AF - G||'
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES15.9)') &
        'ials:',0,'d:',0,'RQ = ',RQ,conver
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
            call UPDATEP(AF,F,d,B1,.TRUE.)
            call UPDATEP(F,d,B2,.TRUE.)

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
            call UPDATEP(AF,F,d,B1,.FALSE.)
            call UPDATEP(F,d,B2,.FALSE.)

!           Compute ||AF-G|| as a convergence check (always do it since
!           it's cheap compared to building+solving the linear system!)
            AFAF=PVVfromPS(AF,AF,BB)
            AFG=PVVfromPS(AF,G,PS)
            FAF=PVVfromPS(AF,F,B1)
            FF=PVVfromPS(F,F,B2)

            RQ=(FAF+ishift*Eshift*FF)/FF
            conver=sqrt(abs((AFAF+GG-2*AFG)/GG))

            IF (show) THEN
               write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES15.9)') &
               'ials:',itn,'d:',d,'RQ = ',RQ,conver
            ENDIF

         enddo  ! loop over d
      ENDDO  ! loop over iterations

      IF (show) write(*,*) 

      deallocate(BB,PS,B1,B2)
      call FlushCP(ATA)
      call FlushCP(ATG)
      call FlushCP(AF)

      end subroutine LinSolver_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LinSolver_2(A,F,G,nitn,ishift,Eshift,showconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver, no normal form

      implicit none
      TYPE (CP), INTENT(IN)    :: A,G
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: AF,Id
      logical, optional, intent(in) :: showconv
      logical :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,B2,lhs,rhs
      real*8  :: conver,tmp,FF,FAF,FG,GG,RQ
      integer :: rA,rAs,rG,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical, allocatable :: sym(:)

      IF (nitn.eq.0) return

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

      show=.FALSE.
      IF (present(showconv).and.mpirank.eq.mpi_prnt_rank) show=showconv

!     Precompute identity matrix for Eshift
      allocate(sym(ndof))
      sym(:)=.FALSE.
      Id=IdentityCPMatrix(A%rows,A%cols,sym)
      call VecScalarMult(Id,-Eshift)
      deallocate(sym)

!     BB(l,l') = Pi_{i=2}^ndof < (F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < (F)_i^l , (G)_i^l' >
      allocate(BB(rAF,rF),PS(rG,rF),B2(rF,rF))
      call CONSTPT(F,AF,0,BB)
      call CONSTPT(F,G,0,PS)
      call CONSTPT(F,F,0,B2)

      IF (show) THEN
         write(*,*)
         write(*,*) "LinSolver_2() iterations..."
      ENDIF

      FAF=PVVfromPS(F,AF,BB)
      FG=PVVfromPS(F,G,PS)
      FF=PVVfromPS(F,F,B2)
      GG=PRODVV(G)

      RQ=(FAF+ishift*Eshift*FF)/FF

      IF (show) THEN
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,A)') &
        'ials:',0,'d:',0,'Es = ',Eshift,'(<F,AF>/<F,F>)**2'
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
        'ials:',0,'d:',0,'RQ = ',RQ,(FAF/FF)**2
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
            call UPDATEP(F,F,d,B2,.TRUE.)

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
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
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

!           Solve linear system B*c_j_k = b_j_k
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
            call UPDATEP(F,F,d,B2,.FALSE.)

!           Compute ||AF-G|| as a convergence check (always do it since
!           it's cheap compared to building+solving the linear system!)
            FAF=PVVfromPS(F,AF,BB)
            FG=PVVfromPS(F,G,PS)
            FF=PVVfromPS(F,F,B2)

            RQ=(FAF+ishift*Eshift*FF)/FF

            IF (show) &
               write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
              'ials:',itn,'d:',d,'RQ = ',RQ,(FAF/FF)**2

         enddo  ! loop over d
      ENDDO  ! loop over iterations

      IF (show) write(*,*) 

      deallocate(BB,PS,B2)
      call FlushCP(Id)
      call FlushCP(AF)

      end subroutine LinSolver_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LintertwinedInvItn_1(A,F,nitn,ishift,Eshift,showconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver

      implicit none
      TYPE (CP), INTENT(IN)    :: A
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: ATA,AF
      logical, optional, intent(in) :: showconv
      logical :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,B2,PS,lhs,rhs
      real*8  :: AFAF,AFF,FF,RQ,rqold,conver,tmp
      integer :: rA,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind

      IF (nitn.eq.0) return

!     Precompute ATA, ATG, initial AF
      call CPMM(A,ishift,Eshift,.TRUE.,A,ishift,Eshift,.FALSE.,ATA)
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Set parameters
      rA=SIZE(A%coef)
      IF (ishift.ne.0) rA=rA+1
      rF=SIZE(F%coef)
      rAF=rA*rF
      ndof=F%D()

      show=.FALSE.
      IF (present(showconv).and.mpirank.eq.mpi_prnt_rank) show=showconv

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (A*F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < F_i^l , (A*F)_i^l' >
      allocate(BB(rAF,rAF),PS(rF,rAF),B2(rF,rF))
      call CONSTPT(AF,0,BB)
      call CONSTPT(AF,F,0,PS)
      call CONSTPT(F,0,B2)

      IF (show) THEN
         write(*,*)
         write(*,*) "LinSolver_1() iterations..."
      ENDIF

      AFAF=PVVfromPS(AF,AF,BB)
      AFF=PVVfromPS(AF,F,PS)
      FF=PVVfromPS(F,F,B2)

      rqold=Eshift
      RQ=(AFF+ishift*Eshift*FF)/FF
      conver=sqrt(abs((AFAF+FF-2*AFF)/FF))

      IF (show) THEN
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,A)') &
        'ials:',0,'d:',0,'Es = ',Eshift,'||AF - F||'
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
        'ials:',0,'d:',0,'RQ = ',RQ,conver
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
            call UPDATEP(AF,F,d,PS,.TRUE.)
            call UPDATEP(F,d,B2,.TRUE.)

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
               do t=1,rF
                  atgrk=(s-1)*rF+t ! rank of ATG from AT and G
                  do l=1,m
                     rcola=(l-1)*n
                     do k=1,n
                        rrowa=(k-1)*rF ! row group of rhs
                        atgind=AF%ibas(d)-1+rcola+k ! element of AF
                        tmp=AF%coef(atgrk)*AF%base(atgind,atgrk)
                        do i=1,rF
                           rrow=rrowa+i ! row of rhs
                           pcol=pcola+i ! col of PS
                           rhs(rrow,l)=rhs(rrow,l)+tmp*PS(t,pcol)
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!           Solve linear system B*c_j_k = b_j_k
!           (B includes all inner products except the kth)
            call SolveLinSys(lhs,rhs,als_penalty,als_solver)

!           Construct improved F
            call UpdateFfromSoln(F,rhs,d)
            deallocate(lhs,rhs)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, we're doomed for now.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'LintertwinedInvItn_1: NaN on update; itn = ',itn,&
                          '; k = ',d
               call AbortWithError('LinSolver_1 crashed')
            ENDIF

!           Update B2 and FF, normalizing to prevent overflow
            call UPDATEP(F,d,B2,.FALSE.)
            FF=PVVfromPS(F,F,B2)
            F%coef=F%coef*1/sqrt(abs(FF))
            FF=1.d0

!           Update ATA, AF, BB, PS using new Fs
            call CPMM(A,ishift,rqold,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call UPDATEP(AF,d,BB,.FALSE.)
            call UPDATEP(AF,F,d,PS,.FALSE.)

            AFF=PVVfromPS(AF,F,PS)
            RQ=(AFF+ishift*rqold*FF)/FF
            call CPMMreshift(A,1,rqold,RQ,F,0,rqold,RQ,AF)
            call CPMMreshift(A,1,rqold,RQ,A,1,rqold,RQ,ATA)
            rqold=RQ

            AFAF=PVVfromPS(AF,AF,BB)
            conver=sqrt(abs((AFAF+FF-2*AFF)/FF))
          
            IF (show) THEN
               write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
               'ials:',itn,'d:',d,'RQ = ',RQ,conver
            ENDIF

         enddo  ! loop over d
      ENDDO  ! loop over iterations

      IF (show) write(*,*) 

      deallocate(BB,PS,B2)
      call FlushCP(ATA)
      call FlushCP(AF)

      end subroutine LintertwinedInvItn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine LintertwinedInvItn_2(A,F,nitn,ishift,Eshift,showconv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided linear solver, no normal form

      implicit none
      TYPE (CP), INTENT(IN)    :: A
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: AF,Id
      logical, optional, intent(in) :: showconv
      logical :: show
      integer, intent(in) :: nitn,ishift
      real*8, intent(in)  :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,lhs,rhs
      real*8  :: conver,tmp,FF,FAF,FG,RQ
      integer :: rA,rAs,rF,rAF,rFn,ndof
      integer :: d,i,j,k,l,m,n,s,t
      integer :: ir,imod,gst,itn
      integer :: browa,brow,bcola,bcol,lrow,lrowa,lcol,lcola,atark
      integer :: rrowa,rrow,rcola,pcola,pcol,atgrk,atgind,ataind
      logical, allocatable :: sym(:)

      IF (nitn.eq.0) return

!     Precompute initial AF
      call CPMM(A,ishift,Eshift,.FALSE.,F,0,0.d0,.FALSE.,AF)

!     Set parameters
      rA=SIZE(A%coef)
      rAs=rA
      IF (ishift.ne.0) rAs=rAs+1
      rF=SIZE(F%coef)
      rAF=rAs*rF
      ndof=SIZE(F%nbas)

      show=.FALSE.
      IF (present(showconv).and.mpirank.eq.mpi_prnt_rank) show=showconv

!     Precompute identity matrix for Eshift
      allocate(sym(ndof))
      sym(:)=.FALSE.
      Id=IdentityCPMatrix(A%rows,A%cols,sym)
      call VecScalarMult(Id,-Eshift)
      deallocate(sym)

!     BB(l,l') = Pi_{i=2}^ndof < (A*F)_i^l , (F)_i^l' >
!     PS(l,l') = Pi_{i=2}^ndof < F_i^l , (F)_i^l' >
      allocate(BB(rAF,rF),PS(rF,rF))
      call CONSTPT(F,AF,0,BB)
      call CONSTPT(F,F,0,PS)

      IF (show) THEN
         write(*,*)
         write(*,*) "LintertwinedInvItn_2() iterations..."
      ENDIF

      FAF=PVVfromPS(F,AF,BB)
      FF=PVVfromPS(F,F,PS)
      RQ=Eshift

      IF (show) THEN
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,A)') &
        'ials:',0,'d:',0,'Es = ',Eshift,'(<F,AF>/<F,F>)**2'
         write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
        'ials:',0,'d:',0,'RQ = ',(FAF+ishift*Eshift*FF)/FF,(FAF/FF)**2
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
            call UPDATEP(F,F,d,PS,.TRUE.)

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
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
                        lrow=lrowa+i ! row of lhs from A and F
                        do j=1,rF
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
                     do i=1,rF
                        brow=browa+i ! row of BB from A and F
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
            do t=1,rF
               do l=1,m
                  rcola=(l-1)*n
                  do k=1,n
                     rrowa=(k-1)*rF ! row group of rhs
                     atgind=F%ibas(d)-1+rcola+k
                     tmp=F%coef(t)*F%base(atgind,t)
                     do i=1,rF
                        rrow=rrowa+i ! row of rhs
                        rhs(rrow,l)=rhs(rrow,l)+tmp*PS(t,i)
                     enddo
                  enddo
               enddo
            enddo

!           Solve linear system B*c_j_k = b_j_k
!           (B includes all inner products except the kth)
            call SolveLinSys(lhs,rhs,als_penalty,als_solver)

!           Construct improved F
            call UpdateFfromSoln(F,rhs,d)
            deallocate(lhs,rhs)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, we're doomed for now.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'LintertwinedInvItn_2: NaN on update; itn = ',itn,&
                          '; k = ',d
               call AbortWithError('LinSolver_2 crashed')
            ENDIF

!           Update PS and FF, normalizing to prevent overflow
            call UPDATEP(F,F,d,PS,.FALSE.)
            FF=PVVfromPS(F,F,PS)
            F%coef=F%coef*1/sqrt(abs(FF))
            FF=1.d0

!           Update AF, BB from new Fs
            call CPMM(A,ishift,RQ,.FALSE.,F,0,0.d0,.FALSE.,AF,d)
            call UPDATEP(F,AF,d,BB,.FALSE.)
            FAF=PVVfromPS(F,AF,BB)
            RQ=(FAF+ishift*RQ*FF)

            IF (show) &
               write(*,'(X,A,X,I5,X,A,I2,X,A,f22.11,3X,ES9.3)') &
              'ials:',itn,'d:',d,'RQ = ',RQ,(FAF/FF)**2

         enddo  ! loop over d
      ENDDO  ! loop over iterations

      IF (show) write(*,*) 

      deallocate(BB,PS)
      call FlushCP(Id)
      call FlushCP(AF)

      end subroutine LintertwinedInvItn_2

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
