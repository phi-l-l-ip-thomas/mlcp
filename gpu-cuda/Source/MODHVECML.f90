!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODHVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use ERRORTRAP
      use UTILS
      use MODVECVEC
      use SEPDREPN
      use HAMILSETUP

      implicit none
      real*8, private  :: mvp_time=0.d0
      logical, private :: MVP_SETUP = .FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializePRODHVModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      mvp_time = 0.d0
      MVP_SETUP = .TRUE.

      end subroutine InitializePRODHVModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposePRODHVModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      MVP_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total matrix-vector product time  (s)',&
                            mvp_time

      end subroutine DisposePRODHVModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHV(F,G,H,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product H*F = G, where H is the Hamiltonian and
! F and G are vectors, all of which are in CP-format

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: F
      TYPE (CPvec), INTENT(OUT) :: G
      integer, intent(in) :: ishift
      real*8, intent(in)  :: Eshift
      integer :: i,j,k,ik,rF,rG,rH,rFH,ndof,gi,gf,npop
      real*8  :: t1,t2

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      call CPU_TIME(t1)

!     Set parameters
      ndof=SIZE(F%nbas)   ! # modes
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      rFH=rF*rH
      rG=rFH
      IF (ishift.ne.0) rG=rG+rF
      npop=SIZE(H%pops)

      call NewCPvec(G,F%nbas,rG)

!     Loop over terms in F
      ik=1
      DO k=1,rH
!        Loop over terms in H
         DO i=1,rF
            G%coef(ik)=F%coef(i)
!           Loop over sub-modes
            gi=1
            DO j=1,ndof
               gf=gi+F%nbas(j)-1
               call HVBaseProd(F%base(gi:gf,i),G%base(gi:gf,ik),H,j,k)
               gi=gi+F%nbas(j)
            ENDDO
            ik=ik+1
         ENDDO
      ENDDO

!     Energy shift if applicable
      IF (ishift.ne.0) THEN
         call GenCopyWtoV(G,F,rFH+1,rG,1,rF)
         call VecScalarMult(G,Eshift,rFH+1,rG)
         IF (ishift.gt.0) THEN
            call VecSignChange(G,rFH+1,rG)
         ELSE
            call VecSignChange(G,1,rFH)
         ENDIF
      ENDIF

      call CPU_TIME(t2)
      mvp_time=mvp_time+t2-t1

      end subroutine PRODHV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHVR1(F,G,H,i,k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product H*F = G, where H is the Hamiltonian and
! F and G are vectors, for only ONE term of H and F. No E-shift is done.
! i and k correspond to terms of F and H to use, respectively.

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: F
      TYPE (CPvec), INTENT(OUT) :: G
      integer, intent(in) :: i,k
      integer :: j,rF,rH,ndof,gi,gf,npop
      real*8  :: t1,t2

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      call CPU_TIME(t1)

!     Set parameters
      ndof=SIZE(F%nbas)   ! # modes
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      npop=SIZE(H%pops)

!     Error checking
      IF (i.lt.1 .or. i.gt.rF) &
         call AbortWithError('PRODHVR1(): i is out of range')
      IF (k.lt.1 .or. k.gt.rH) &
         call AbortWithError('PRODHVR1(): k is out of range')

!     Early exit if H is applied to a zero vector
      IF (rF.eq.1 .and. F%coef(1).eq.0.d0) THEN
         call CopyWtoV(G,F)
         RETURN
      ENDIF

      call NewCPvec(G,F%nbas,1)
      G%coef(1)=F%coef(i)

!     Calculate matrix-vector product for each DOF
      gi=1
      DO j=1,ndof
         gf=gi+F%nbas(j)-1
         call HVBaseProd(F%base(gi:gf,i),G%base(gi:gf,1),H,j,k)
         gi=gi+F%nbas(j)
      ENDDO

      call CPU_TIME(t2)
      mvp_time=mvp_time+t2-t1      

      end subroutine PRODHVR1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine HVBaseProd(F,G,H,j,k)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product H*F = G, where H is the Hamiltonian and
! F and G are vectors. j indexes the DOF label; k indexes the term in H

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      integer, intent(in)   :: j,k
      real*8, intent(in)    :: F(:)
      real*8, intent(inout) :: G(:)
      integer :: n,npop
      real*8  :: fac


      n=SIZE(F)
      npop=SIZE(H%pops)
      fac=1.d0
      IF (j.eq.1) fac=H%opcoef(k)

!     Identity operator (just copy F)
      IF (H%opindx(k,j).eq.0) THEN
         G(:)=fac*F(:)
      ELSEIF (H%opindx(k,j).gt.0) THEN
!        Primitive operator (look in pops array)
         IF (H%opindx(k,j).le.npop) THEN
            call dspmv('U',n,fac,H%pops(H%opindx(k,j))%mat,F,1,&
                       0.d0,G,1)
!        Compound operator (look in cops array)
         ELSE
            call dspmv('U',n,fac,H%cops(H%opindx(k,j)-npop)%mat,F,1,&
                       0.d0,G,1)
         ENDIF
!     Eigenoperator (look in eops array)
      ELSE
         call dspmv('U',n,fac,H%eops(abs(H%opindx(k,j)))%mat,F,1,&
                    0.d0,G,1)
      ENDIF

      end subroutine HVBaseProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MODHVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
