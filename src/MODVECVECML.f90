!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODVECVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Basic operations (add, dot product, normalize) for CP-format vectors

      USE SEPDREPN
      USE LINALG

      implicit none
      real*8, private  :: norm_time=0.d0,pvv_time=0.d0,svv_time=0.d0
      logical, private :: MVV_SETUP = .FALSE.

      INTERFACE PRODVV
         MODULE PROCEDURE PRODVV_AA,PRODVV_AB
      END INTERFACE PRODVV

      INTERFACE ordre
         MODULE PROCEDURE ordre_AA,ordre_AB
      END INTERFACE ordre

      INTERFACE SUMLCVEC
         MODULE PROCEDURE SUMLCVEC_facs,SUMLCVEC_nofacs
      END INTERFACE SUMLCVEC

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializeMVV()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      norm_time = 0.d0
      pvv_time  = 0.d0
      svv_time = 0.d0
      MVV_SETUP = .TRUE.

      end subroutine InitializeMVV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeMVV()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MVV_SETUP) call InitializeMVV()

      MVV_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total vector inner-product time   (s)',&
                            pvv_time
      write(*,'(X,A,X,f20.3)') 'Total vector normalization time   (s)',&
                            norm_time
      write(*,'(X,A,X,f20.3)') 'Total vector-vector addition time (s)',&
                            svv_time

      end subroutine DisposeMVV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecSignChange(F,ri,re)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Changes sign of CP-vec by negating the base for the first DOF. ri and 
! re are the rank indices over which over which to change the sign

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in) :: ri,re
      integer :: n

      n=F%nbas(1)
      F%base(1:n,ri:re)=-F%base(1:n,ri:re)

      end subroutine VecSignChange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecScalarMult(F,fac,ri,re)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies CP-vec by scalar factor, which is distributed over all DOFs
! If the factor is negative then the sign is also changed

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in) :: ri,re
      real*8, intent(in)  :: fac
      integer :: ndof

      ndof=SIZE(F%nbas)
      F%coef(ri:re)=abs(fac)*F%coef(ri:re)
      IF (fac.lt.0.d0) call VecSignChange(F,ri,re)

      end subroutine VecScalarMult

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SUMVECVEC(F,Ffac,G,Gfac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums Ffac*F and Gfac*G in CP-format, with factors distributed over all
! dofs. The summed vector replaces F; G is unchanged.

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec), INTENT(IN)    :: G
      real*8, intent(in) :: Ffac,Gfac
      integer, allocatable :: nbas(:)
      integer :: rF,rG,ndof
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      IF (.NOT.CHECKNBAS(F,G)) THEN
         write(*,*) 'Dimensions or type of F and G do not match'
         CALL AbortWithError('Error in SUMVECVEC()')
      ENDIF

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)
      ndof=SIZE(F%nbas)
      ALLOCATE(nbas(ndof))

!     If one of the two terms is zero, it is replaced by the other
      IF (rF.eq.1 .and. F%coef(1).eq.0.d0) THEN
         call FlushCPvec(F)
         call CopyWtoV(F,G)
         call VecScalarMult(F,Gfac,1,rG)

      ELSE IF (rG.eq.1 .and. G%coef(1).eq.0.d0) THEN
         call VecScalarMult(F,Ffac,1,rF)

!     General case
      ELSE
         call ResizeV(F,rF+rG)
         call GenCopyWtoV(F,G,rF+1,rF+rG,1,rG)
         call VecScalarMult(F,Ffac,1,rF)
         call VecScalarMult(F,Gfac,rF+1,rF+rG)
      ENDIF

      call CPU_TIME(t2)
      svv_time=svv_time+t2-t1

      end subroutine SUMVECVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SUMLCVEC_facs(F,Q,facs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums a linear combination of vectors in Q and stores in F
! Faster than SUMVECVEC for summing many vectors with coefs since only 
! one allocate is performed

      implicit none
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      real*8, intent(in) :: facs(:)
      integer :: i,nbloc,rst,nrki
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      nbloc=SIZE(facs)
      IF (SIZE(Q).ne.nbloc) THEN
         write(*,*) 'Mismatch in dimension of Q, facs'
         call AbortWithError("Error in SUMLCVEC_facs()")
      ENDIF

!     Get the rank of F (=sum of ranks in Q)
      nrki=0
      DO i=1,nbloc
         nrki=nrki+SIZE(Q(i)%coef)
      ENDDO

!     Build F from Q and the factors
      rst=0
      call NewCPvec(F,Q(1)%nbas,nrki)
      DO i=1,nbloc
         nrki=SIZE(Q(i)%coef)
         call GenCopyWtoV(F,Q(i),rst+1,rst+nrki,1,nrki)
         call VecScalarMult(F,facs(i),rst+1,rst+nrki)
         rst=rst+nrki
      ENDDO

      call CPU_TIME(t2)
      svv_time=svv_time+t2-t1

      end subroutine SUMLCVEC_facs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SUMLCVEC_nofacs(F,Q)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums a linear combination of vectors in Q and stores in F
! This version does not require factors, assuming them to be 1

      implicit none
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (CPvec), INTENT(IN)  :: Q(:)
      integer :: i,nbloc,rst,nrki
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      nbloc=SIZE(Q)

!     Get the rank of F (=sum of ranks in Q)
      nrki=0
      DO i=1,nbloc
         nrki=nrki+SIZE(Q(i)%coef)
      ENDDO

!     Build F from Q
      rst=0
      call NewCPvec(F,Q(1)%nbas,nrki)
      DO i=1,nbloc
         nrki=SIZE(Q(i)%coef)
         call GenCopyWtoV(F,Q(i),rst+1,rst+nrki,1,nrki)
         rst=rst+nrki
      ENDDO

      call CPU_TIME(t2)
      svv_time=svv_time+t2-t1

      end subroutine SUMLCVEC_nofacs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ordre_AA(G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts CP-vector H in decreasing order of coefficients, but uses the
! replacement function so that ordre can be called with a single vector

      implicit none
      TYPE (CPvec), INTENT(INOUT)  :: G
      TYPE (CPvec) :: F

      call ordre(G,F)
      call ReplaceVwithW(G,F)

      end subroutine ordre_AA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ordre_AB(G,F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts CP-vector H in decreasing order of coefficients. This only makes
! sense if the base is normalized, so call NORMBASE before this. This 
! version generates a sorted vector in addition to the original

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (CPvec), INTENT(OUT) :: F
      real*8, allocatable :: tabindex(:)
      integer :: rG,i

!     Set parameters       
      rG=SIZE(G%coef)

!     List of ordered indices
      allocate(tabindex(rG))
      do i=1,rG
         tabindex(i)=i
      enddo

!     Sort the coefficients into Fcoef
      call NewCPvec(F,G%nbas,rG)
      F%coef=G%coef
      call dsort(F%coef,tabindex,rG,-2)

!     Reorder the base to match the coefficients
      do i=1,rG
         F%base(:,i)=G%base(:,int(tabindex(i)))
      enddo

      deallocate(tabindex)

      end subroutine ordre_AB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NORMBASE(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Renormalizes Fbase, which modifies Fcoef

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      real*8  :: norm1D,normND,t1,t2
      integer :: ndim,ir,irk,id,imod,rF,gst

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      rF=SIZE(F%coef)
      ndim=SIZE(F%nbas)

      do irk=1,rF
         normND=1.d0
         gst=0
         do id=1,ndim
            norm1D=0.d0
            IF (id.gt.1) gst=gst+F%nbas(id-1)
            do ir=1,F%nbas(id)
               imod=gst+ir
               norm1D=norm1D+F%base(imod,irk)*F%base(imod,irk)
            enddo
            norm1D=sqrt(norm1D)
            do ir=1,F%nbas(id)
               imod=gst+ir
               F%base(imod,irk)=F%base(imod,irk)/norm1D
            enddo
            normND=normND*norm1D
         enddo
         F%coef(irk)=F%coef(irk)*normND
      enddo

      call CPU_TIME(t2)
      norm_time=norm_time+t2-t1

      end subroutine NORMBASE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NORMCOEF(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Renormalizes Fbase with a call to NORMBASE
! then renormalizes Fcoef here

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      real*8  :: normND

      call NORMBASE(F)
      normND=1/sqrt(abs(PRODVV(F)))
      F%coef=F%coef*normND

      end subroutine NORMCOEF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODVV_AA(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the dot product of F with itself, using factor symmetry

      implicit none
      TYPE (CPvec), intent(in) :: F
      real*8  :: PRODVV_AA,t1,t2
      integer :: nrk,i,j

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      nrk=SIZE(F%coef)

      PRODVV_AA=0.d0
!     Upper triangle gets multiplied by 2
      do i=1,nrk-1
         do j=i+1,nrk
            PRODVV_AA=PRODVV_AA+F%coef(i)*F%coef(j)*&
                      PRODND(F%base(:,i),F%base(:,j),F%nbas,0)
         enddo
      enddo
      PRODVV_AA=2*PRODVV_AA

!     Sum of diagonal entries
      do i=1,nrk
         PRODVV_AA=PRODVV_AA+F%coef(i)**2*&
                   PRODND(F%base(:,i),F%base(:,i),F%nbas,0)
      enddo

      call CPU_TIME(t2)
      pvv_time=pvv_time+t2-t1

      end function PRODVV_AA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODVV_AB(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the dot product of F and G
! <F|G> = SUM_if=1,rF SUM_ig=1,rG s_if*s_ig*
!         <F^if_1|G^ig_1>*...*<F^if_ndim|G^ig_ndim >

      implicit none
      TYPE (CPvec), intent(in) :: F,G
      real*8  :: PRODVV_AB,t1,t2
      integer :: rF,rG,i,j

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      IF (.NOT.CHECKNBAS(F,G)) THEN
         write(*,*) 'Error: base dimensions of F and G do not match'
         CALL AbortWithError('Error in PRODVV_AB()')
      ENDIF

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)

      PRODVV_AB=0.d0
      do i=1,rF
         do j=1,rG
            PRODVV_AB=PRODVV_AB+F%coef(i)*G%coef(j)*&
                      PRODND(F%base(:,i),G%base(:,j),F%nbas,0)
         enddo
      enddo

      call CPU_TIME(t2)
      pvv_time=pvv_time+t2-t1

      end function PRODVV_AB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODND(F,G,nbas,kskip)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates dot product between two rank-1 parts of two vectors in
! CP-format, skipping the kskip DOF (set kskip=0 to do full dot product)
! Only the bases are multiplied here

      implicit none
      real*8, intent(in)  :: F(:),G(:)
      integer, intent(in) :: nbas(:)
      integer, intent(in) :: kskip
      real*8  :: PRODND,prod1D
      integer :: i,ndim,gi,j

      ndim=SIZE(nbas)

      PRODND=1.d0
      gi=0
      do i=1,ndim
         IF (i.ne.kskip) THEN
            prod1D=0.d0
            do j=1,nbas(i)
               prod1D=prod1D+F(gi+j)*G(gi+j)
            enddo
            PRODND=PRODND*prod1D
         ENDIF
         gi=gi+nbas(i)
      enddo

      end function PRODND

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTPSTOT(F,G,kskip,PS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSTOT.
! This subroutine calculates PS(l,l')= Pi_{i=1}^ndim <G_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
      
      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in) :: kskip
      integer :: rF,rG,i,j

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)

!     Compute inner products
      PS=0.d0
      DO i=1,rG
         DO j=1,rF
            PS(i,j)=PRODND(F%base(:,j),G%base(:,i),F%nbas,kskip)
         ENDDO
      ENDDO

      end subroutine CONSTPSTOT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBTOT(F,kskip,BB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BB must be allocated before the call to CONSTBBTOT.
! This subroutine calculates BB(l,l')= Pi_{i=1}^ndim <F_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Different from CONSTPSTOT since BB matrix is symmetric
! Skips the kskip-th DOF (set kskip=0 to include all DOF)

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: kskip
      integer :: rF,i,j

      rF=SIZE(F%coef)

!     Compute inner products 
      BB=0.d0
      DO i=1,rF
         DO j=i,rF
            BB(i,j)=PRODND(F%base(:,j),F%base(:,i),F%nbas,kskip)
         ENDDO
      ENDDO

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
         ENDDO
      ENDDO

      end subroutine CONSTBBTOT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTPSk(F,G,k,PSk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PSk must be allocated before the call to CONSTPSk. 
! This subroutine calculates the matrix PSk(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF
! (dividing factor by which to produce the product excluding id = k)
      
      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PSk(:,:)
      integer, intent(in)   :: k
      real*8  :: prod1D
      integer :: rF,rG,i,j,id,ir,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)

!     Calc basis ranges
      gi=1
      DO id=1,k-1
         gi=gi+F%nbas(id)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner products 
      PSk=0.d0
      DO i=1,rG
         DO j=1,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+G%base(ir,i)*F%base(ir,j)
            ENDDO
            PSk(i,j)=prod1D
         ENDDO
      ENDDO

      end subroutine CONSTPSk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBk(F,k,BBk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BBk must be allocated before the call to CONSTBBk. 
! This subroutine calculates the matrix BBk(l,l')=  < F_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < nrk and 1 < l' < nrk
! (dividing factor by which to produce the product excluding id = k)

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BBk(:,:)
      integer, intent(in)   :: k
      real*8  :: prod1D
      integer :: rF,i,j,id,ir,gi,gf,is

      rF=size(F%coef)

!     Calc basis ranges
      gi=1
      DO id=1,k-1
         gi=gi+F%nbas(id)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner products
      BBk=0.d0
      DO i=1,rF
         DO j=i,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+F%base(ir,i)*F%base(ir,j)
            ENDDO
            BBk(i,j)=prod1D
         ENDDO
      ENDDO

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BBk(i,j)=BBk(j,i)
         ENDDO
      ENDDO

      end subroutine CONSTBBk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdatePS(F,G,k,PS,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PSk must be allocated before the call to CONSTPSk. 
! This subroutine calculates the matrix PSk(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF
! (dividing factor by which to produce the product excluding id = k)

      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k
      logical, intent(in)   :: divide
      real*8  :: prod1D
      integer :: rF,rG,i,j,id,ir,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)

!     Calc basis ranges
      gi=1
      DO id=1,k-1
         gi=gi+F%nbas(id)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner products 
      DO i=1,rG
         DO j=1,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+G%base(ir,i)*F%base(ir,j)
            ENDDO
            IF (divide) prod1D=1.d0/prod1D
            PS(i,j)=PS(i,j)*prod1D
         ENDDO
      ENDDO

      end subroutine UpdatePS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateBB(F,k,BB,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine updates matrix BB(l,l')= PI_kp < F_kp^l , F_kp^l' >
! by multiplying by < F_k^l , F_k^l' >, (divide=.FALSE.)
! or dividing by < F_k^l , F_k^l' >, (divide=.TRUE.)

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: k
      logical, intent(in) :: divide
      real*8  :: prod1D
      integer :: rF,i,j,id,ir,gi,gf,is

      rF=size(F%coef)

!     Calc basis ranges
      gi=1
      DO id=1,k-1
         gi=gi+F%nbas(id)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner products
      DO i=1,rF
         DO j=i,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+F%base(ir,i)*F%base(ir,j)
            ENDDO
            IF (divide) prod1D=1.d0/prod1D
            BB(i,j)=BB(i,j)*prod1D
         ENDDO
      ENDDO

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
         ENDDO
      ENDDO

      end subroutine UpdateBB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTPSTOT_OMP(F,G,kskip,PS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSTOT.
! This subroutine calculates PS(l,l')= Pi_{i=1}^ndim <G_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in) :: kskip
      integer :: rF,rG,rT,i,j,l

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)
      rT=rG*rF

!     Compute inner products
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         PS(i,j)=PRODND(F%base(:,j),G%base(:,i),F%nbas,kskip)
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine CONSTPSTOT_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBTOT_OMP(F,kskip,BB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BB must be allocated before the call to CONSTBBTOT.
! This subroutine calculates BB(l,l')= Pi_{i=1}^ndim <F_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Different from CONSTPSTOT since BB matrix is symmetric
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: kskip
      integer :: rF,rT,i,j,l

      rF=SIZE(F%coef)
      rT=rF*rF

!     Compute inner products
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         IF (j.ge.i) &
            BB(i,j)=PRODND(F%base(:,j),F%base(:,i),F%nbas,kskip)
      ENDDO
!$omp enddo
!$omp end parallel

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
         ENDDO
      ENDDO

      end subroutine CONSTBBTOT_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTPSk_OMP(F,G,k,PSk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PSk must be allocated before the call to CONSTPSk.
! This subroutine calculates the matrix PSk(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF
! (dividing factor by which to produce the product excluding id = k)
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PSk(:,:)
      integer, intent(in)   :: k
      integer :: rF,rG,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      rT=rG*rF

!     Calc basis ranges
      gi=1
      DO l=1,k-1
         gi=gi+F%nbas(l)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner product for the k-th DOF
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         PSk(i,j)=dot_product(G%base(gi:gf,i),F%base(gi:gf,j))
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine CONSTPSk_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBk_OMP(F,k,BBk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BBk must be allocated before the call to CONSTBBk.
! This subroutine calculates the matrix BBk(l,l')=  < F_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < nrk and 1 < l' < nrk
! (dividing factor by which to produce the product excluding id = k)
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BBk(:,:)
      integer, intent(in)   :: k
      real*8  :: prod1D
      integer :: rF,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rT=rF*rF

!     Calc basis ranges
      gi=1
      DO l=1,k-1
         gi=gi+F%nbas(l)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner product for the k-th DOF
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         IF (j.ge.i) &
            BBk(i,j)=dot_product(F%base(gi:gf,i),F%base(gi:gf,j))
      ENDDO
!$omp enddo
!$omp end parallel

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BBk(i,j)=BBk(j,i)
         ENDDO
      ENDDO

      end subroutine CONSTBBk_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdatePS_OMP(F,G,k,PS,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSk.
! This subroutine modifies the matrix PS(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF, by multiplying by or dividing out
! PS for the k-th DOF
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k
      logical, intent(in)   :: divide
      real*8  :: prod1D
      integer :: rF,rG,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      rT=rG*rF

!     Calc basis ranges
      gi=1
      DO l=1,k-1
         gi=gi+F%nbas(l)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner product and update or downdate
!$omp parallel
!$omp do private(i,j,l,prod1D)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         prod1D=dot_product(G%base(gi:gf,i),F%base(gi:gf,j))
         IF (divide) prod1D=1.d0/prod1D
         PS(i,j)=PS(i,j)*prod1D
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine UpdatePS_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateBB_OMP(F,k,BB,divide)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine updates matrix BB(l,l')= PI_kp < F_kp^l , F_kp^l' >
! by multiplying by < F_k^l , F_k^l' >, (divide=.FALSE.)
! or dividing by < F_k^l , F_k^l' >, (divide=.TRUE.)
! OpenMP parallel version

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: k
      logical, intent(in) :: divide
      real*8  :: prod1D
      integer :: rF,rT,i,j,l,gi,gf,is

      rF=size(F%coef)
      rT=rF*rF

!     Calc basis ranges
      gi=1
      DO l=1,k-1
         gi=gi+F%nbas(l)
      ENDDO
      gf=gi+F%nbas(k)-1

!     Compute inner product and update or downdate
!$omp parallel
!$omp do private(i,j,l,prod1D)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         IF (j.ge.i) THEN
            prod1D=dot_product(F%base(gi:gf,i),F%base(gi:gf,j))
            IF (divide) prod1D=1.d0/prod1D
            BB(i,j)=BB(i,j)*prod1D
         ENDIF
      ENDDO
!$omp enddo
!$omp end parallel

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
         ENDDO
      ENDDO

      end subroutine UpdateBB_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function COND(F) 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the 'condition number' for multiD vector F
! in separated representation: F = sum _l s_l F^l_1 x ... x F^l_d
! COND = sqrt (sum_l=1^rank s_l) / ||F||

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      real*8  :: COND,num
      integer :: rF,i
      
      rF=SIZE(F%coef)
      num=0.d0
      do i=1,rF
         num=num+F%coef(i)*F%coef(i)
      enddo
      COND=sqrt(num/PRODVV(F))

      end function COND

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function calc_FmG(G,F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes the two norm ||F-G||_2

      implicit none
      TYPE (CPvec), INTENT(IN) :: G,F
      TYPE (CPvec)  :: H
      integer :: rF,rG
      real*8  :: calc_FmG

      rF=SIZE(F%coef)

      call CopyWtoV(H,F)
      call SUMVECVEC(H,1.d0,G,-1.d0)
      calc_FmG=sqrt(abs(PRODVV(H)))
      call FlushCPvec(H)

      end function calc_FmG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function diagnorm(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes sum-of-squares of diagonal elements of a tensor with all 
! dimensions equal

      implicit none
      TYPE (CPvec), INTENT(IN) :: F
      integer :: i,j,k,l,n,rF,ndof
      real*8  :: prod1,sum1,diagnorm

!     Parameters
      ndof=SIZE(F%nbas)
      rF=SIZE(F%coef)
      n=F%nbas(1)

!     Error if dimensions are not all the same
      DO j=2,ndof
         IF (F%nbas(j).ne.n) THEN
            write(*,*) 'n, bad-DOF, n_bad-DOF = ',n,j,F%nbas(j)
            call AbortWithError('diagnorm(): unequal dimensions')
         ENDIF
      ENDDO

!     Compute the multiplicitive trace
      diagnorm=0.d0
      DO k=1,n
         sum1=0.d0
         DO i=1,rF
            prod1=F%coef(i)
            l=k
            DO j=1,ndof
               prod1=prod1*F%base(l,i)
               l=l+n
            ENDDO
            sum1=sum1+prod1
         ENDDO
         diagnorm=diagnorm+sum1**2
      ENDDO

      end function diagnorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CPMatVecProd(M,v,w,trans)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix-vector product: CP-format matrix times CP-format vector:
! M * v = w. Set trans = .TRUE. to do M^T * v = w. M must be square.

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN)  :: M,v
      TYPE (CPvec), INTENT(OUT) :: w
      logical, intent(in)  :: trans
      integer, allocatable :: nbas(:),gst(:,:)
      integer :: gi(3),gf(3)
      integer :: i,j,k,ik,l,rkM,rkv,rkw,ndof,rtot

      rkM=SIZE(M%coef)
      rkv=SIZE(v%coef)
      rkw=rkM*rkv
      ndof=SIZE(M%nbas)
      rtot=rkw*ndof

!     Error checking
      IF (ndof.ne.SIZE(v%nbas)) THEN
         write(*,*) '# DOF in M: ',ndof,'; # DOF in v: ',SIZE(M%nbas)
         call AbortWithError("CPMatVecProd(): M, v have # DOF mismatch")
      ENDIF

!     Determine nbas for w, make sure M, v have consistent sizes
      ALLOCATE(nbas(ndof))
      DO j=1,ndof
         IF (mod(M%nbas(j),v%nbas(j)).ne.0) THEN
            write(*,*) 'DOF ',j,': SIZE(M) = ',M%nbas(j),&
                                '; SIZE(v) = ',v%nbas(j)
            call AbortWithError("CPMatVecProd(): M, v size mismatch")
         ENDIF
         nbas(j)=M%nbas(j)/v%nbas(j)
      ENDDO

!     Compute the matrix-vector product; store in w
      call NewCPvec(w,nbas,rkw)

      ALLOCATE(gst(ndof,3))
      gst(1,:)=1
      DO j=2,ndof
         gst(j,1)=gst(j-1,1)+M%nbas(j-1)
         gst(j,2)=gst(j-1,2)+v%nbas(j-1)
         gst(j,3)=gst(j-1,3)+w%nbas(j-1)         
      ENDDO

!$omp parallel
!$omp do private(l,j,gi,gf,ik,k,i)
      DO l=1,rtot
         j=mod(l-1,ndof)+1
         gi(:)=gst(j,:)
         gf(1)=gi(1)+M%nbas(j)-1
         gf(2)=gi(2)+v%nbas(j)-1
         gf(3)=gi(3)+w%nbas(j)-1
         ik=(l-1)/ndof+1
         k=mod(ik-1,rkM)+1
         i=(ik-1)/rkM+1
         w%coef(ik)=v%coef(i)*M%coef(k)
         call WrappedMatVecProd(M%base(gi(1):gf(1),k),&
                                v%base(gi(2):gf(2),i),&
                                w%base(gi(3):gf(3),ik),trans)
      ENDDO
!$omp end do
!$omp end parallel

      DEALLOCATE(nbas,gst)

      END SUBROUTINE CPMatVecProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CPMatOrdVecProd(M,v,w)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix-vector product: matrix in CP-format times ordinary vector
! M * v = w

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN) :: M
      REAL*8, INTENT(IN)  :: v(:)
      REAL*8, INTENT(OUT) :: w(:)
      INTEGER :: i,j,k,nr,nc,nrk
      REAL*8 :: wtmp

!     Error checking
      IF (SIZE(M%nbas).ne.2) THEN
         write(*,*) 'Must have exactly 2 DOFs in CP-mat-vec-product'
         call AbortWithError('Error in CP2MatOrdVecProd()')
      ENDIF
      IF (M%nbas(2).ne.SIZE(v)) THEN
         write(*,*) 'Mismatch in matrix and vector-in dimensions'
         call AbortWithError('Error in CP2MatOrdVecProd()')
      ENDIF
      IF (M%nbas(1).ne.SIZE(w)) THEN
         write(*,*) 'Mismatch in matrix and vector-out dimensions'
         call AbortWithError('Error in CP2MatOrdVecProd()')
      ENDIF

!     Set parameters
      nrk=SIZE(M%coef)
      nr=M%nbas(1)
      nc=M%nbas(2)

      w=0.d0

!     Matrix-vector product
      DO i=1,nrk
         wtmp=0.d0
         DO j=1,nc
            wtmp=wtmp+M%base(nr+j,i)*v(j)
         ENDDO
         wtmp=M%coef(i)*wtmp
         DO j=1,nr
            w(j)=w(j)+M%base(j,i)*wtmp
         ENDDO
      ENDDO

      END SUBROUTINE CPMatOrdVecProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CPMatMatProd(M1,t1,M2,t2,M3,nr1,nc1,nr2,nc2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix-matrix product: matrices in CP-format M3=M1*M2
! t1,t2 = transpose of 1st,2nd matrices (choose .FALSE. for normal)

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN)  :: M1,M2
      TYPE (CPvec), INTENT(OUT) :: M3
      LOGICAL, INTENT(IN)  :: t1,t2
      INTEGER, INTENT(IN)  :: nr1(:),nc1(:),nr2(:),nc2(:)
      INTEGER, ALLOCATABLE :: nbas(:)
      REAL*8, ALLOCATABLE  :: tmat1(:,:),tmat2(:,:),tmat3(:,:),tvec(:)
      INTEGER :: gi(3),gf(3),ndof,rk1,rk2,rk3,i,j,k,ik

!     Initializations
      ndof=SIZE(M1%nbas)
      rk1=SIZE(M1%coef)
      rk2=SIZE(M2%coef)
      rk3=rk1*rk2

!     Error checking
      IF (SIZE(M2%nbas).ne.ndof .or. &
          SIZE(nr1).ne.ndof .or. SIZE(nr2).ne.ndof .or. &
          SIZE(nc1).ne.ndof .or. SIZE(nc2).ne.ndof) call &
         AbortWithError('CPMatMatProd(): inconsistent # DOF')

!     Determine nbas based on the expected dimensions of M3. If there
!     is a mismatch in dimension, MatrixMult() will trap this error,
!     so it is not tested here
      ALLOCATE(nbas(ndof))
      IF (t1) THEN
         nbas=nc1
      ELSE
         nbas=nr1
      ENDIF
      IF (t2) THEN
         nbas=nbas*nr2
      ELSE
         nbas=nbas*nc2
      ENDIF
!     Use NewCPvec() instead of NewCPmat() since M3 may not be square
      call NewCPvec(M3,nbas,rk3)
      DEALLOCATE(nbas)

!     Multiply the coefficients
      ik=1
      DO k=1,rk2
         DO i=1,rk1
            M3%coef(ik)=M1%coef(i)*M2%coef(k)
            ik=ik+1
         ENDDO
      ENDDO

!     Multiply the matrices stored in M1%base and M2%base
      gi(:)=1
      DO j=1,ndof
         gf(1)=gi(1)+M1%nbas(j)-1
         gf(2)=gi(2)+M2%nbas(j)-1
         gf(3)=gi(3)+M3%nbas(j)-1

!        Loop over terms in H
         ik=1
         DO k=1,rk2
!           Get the matrix in M2 for the j-th DOF
            call Vec2Mat(M2%base(gi(2):gf(2),k),tmat2,nr2(j),nc2(j))
            DO i=1,rk1
!              Get the matrix in M1 for the j-th DOF
               call Vec2Mat(M1%base(gi(1):gf(1),i),tmat1,nr1(j),nc1(j))
!              Matrix multiply tmat1*tmat2=tmat3
               call MatrixMult(tmat1,t1,tmat2,t2,tmat3)
!              Convert product matrix to vector; then put in M3
               call Mat2Vec(tvec,tmat3,.FALSE.)
               M3%base(gi(3):gf(3),ik)=tvec
               DEALLOCATE(tmat1,tmat3,tvec)
               ik=ik+1
            ENDDO
            DEALLOCATE(tmat2)
         ENDDO
         gi(1)=gi(1)+M1%nbas(j)
         gi(2)=gi(2)+M2%nbas(j)
         gi(3)=gi(3)+M3%nbas(j)
      ENDDO

      END SUBROUTINE CPMatMatProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE CPMatMatProd2(M1,tp1,M2,tp2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix-matrix product: matrices in CP-format M3=M1*M2
! tp1,tp2 = transpose of 1st,2nd matrices (choose .FALSE. for normal)

      IMPLICIT NONE
      TYPE (CPvec), INTENT(IN)  :: M1,M2
      TYPE (CPvec), INTENT(OUT) :: M3
      TYPE (CPvec) :: v
      LOGICAL, INTENT(IN) :: tp1,tp2
      INTEGER :: nbas3(2)
      INTEGER :: i,j,k,r1,c1,r2,c2,kr1,kc1,kr2,kc2,nrk1,nrk2
      INTEGER :: nr1,nc1,nr2,nc2
      REAL*8  :: wtmp

!     Error checking
      IF (SIZE(M1%nbas).ne.2 .or. SIZE(M2%nbas).ne.2) THEN
         write(*,*) 'Error: must have 2 DOFs in CP-mat-mat-product'
         call AbortWithError('Error in CP2MatMatProd2()')
      ENDIF

!     Set parameters: some must be swapped if 'tpX' is .TRUE. (matrix
!     transpose case)
      nrk1=SIZE(M1%coef)
      nrk2=SIZE(M2%coef)
      r1=1
      c1=2
      kr1=0
      kc1=M1%nbas(1)
      IF (tp1) THEN
         r1=2
         c1=1
         kr1=M1%nbas(1)
         kc1=0
      ENDIF
      r2=1
      c2=2
      kr2=0
      kc2=M2%nbas(1)
      IF (tp2) THEN
         r2=2
         c2=1
         kr2=M2%nbas(1)
         kc2=0
      ENDIF
      nr1=M1%nbas(r1)
      nc1=M1%nbas(c1)
      nr2=M2%nbas(r2)
      nc2=M2%nbas(c2)

      IF (nc1.ne.nr2) THEN
         write(*,*) 'Error: mismatch in matrix dimensions'
         call AbortWithError('Error in CP2MatMatProd2()')
      ENDIF

      nbas3=(/nr1,nc2/)

      call GetZeroCPvec(M3,nbas3)

!     Matrix-matrix product
      DO i=1,nrk1
         DO j=1,nrk2
            call NewCPvec(v,nbas3,1)
            v%base(1:nr1,1)=M1%base(kr1+1:kr1+nr1,i)
            v%base(nr1+1:nr1+nc2,1)=M2%base(kc2+1:kc2+nc2,j)
            v%coef(1)=M1%coef(i)*M2%coef(j)
            wtmp=0.d0
            DO k=1,nc1
               wtmp=wtmp+M1%base(kc1+k,i)*M2%base(kr2+k,j)
            ENDDO
            call SUMVECVEC(M3,1.d0,v,wtmp)
            call FlushCPvec(v)
         ENDDO
      ENDDO

      END SUBROUTINE CPMatMatProd2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProjectOut(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Removes contributions of G from F, returning F.
! No reduction or normalization is done here

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec), INTENT(IN)    :: G
      real*8  :: vivj

      vivj=PRODVV(F,G)
      call SUMVECVEC(F,1.d0,G,-vivj)

      end subroutine ProjectOut

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AugmentVWithRandom(F,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments F to rank nrk by adding random terms with small coefficients

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec) :: G
      real*8, parameter   :: smallnr=1.d-12
      integer, intent(in) :: nrk
      integer :: rF

      rF=SIZE(F%coef)

!     No need to do anything if the rank of F is already large enough
      IF (rF.ge.nrk) RETURN

      call GetRandomCPvec(G,F%nbas,nrk-rF)
      call NORMBASE(G)
      G%coef=MINVAL(F%coef)*smallnr

      call SUMVECVEC(F,1.d0,G,1.d0)

      end subroutine AugmentVWithRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module MODVECVEC 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
