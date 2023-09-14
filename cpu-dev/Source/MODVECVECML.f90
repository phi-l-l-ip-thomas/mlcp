!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MODVECVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Basic operations (add, dot product, normalize) for CP-format vectors

      USE SEPDREPN
      USE LINALG

      implicit none
      real*8, private  :: norm_time=0.d0,pvv_time=0.d0,svv_time=0.d0
      logical, private :: MVV_SETUP = .FALSE.

      INTERFACE VecScalarMult
         MODULE PROCEDURE VecScalarMult_all,VecScalarMult_gen
      END INTERFACE VecScalarMult

      INTERFACE PRODVV
         MODULE PROCEDURE PRODVV_AA,PRODVV_AB
      END INTERFACE PRODVV

      INTERFACE ordre
         MODULE PROCEDURE ordre_AA,ordre_AB
      END INTERFACE ordre

      INTERFACE CONSTPT
         MODULE PROCEDURE CONSTPSTOT,CONSTBBTOT
         MODULE PROCEDURE CONSTPSTOT_OMP,CONSTBBTOT_OMP
      END INTERFACE

      INTERFACE CONSTPk
         MODULE PROCEDURE CONSTPSk,CONSTBBk
         MODULE PROCEDURE CONSTPSk_OMP,CONSTBBk_OMP
      END INTERFACE

      INTERFACE UPDATEP
         MODULE PROCEDURE UpdatePS,UpdateBB
         MODULE PROCEDURE UpdatePS_OMP,UpdateBB_OMP
      END INTERFACE

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
      TYPE (CP), INTENT(INOUT) :: F
      integer, intent(in) :: ri,re
      integer :: n

      n=F%nbas(1)
      F%base(1:n,ri:re)=-F%base(1:n,ri:re)

      end subroutine VecSignChange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecScalarMult_all(F,fac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies CP-vec by scalar factor, which is distributed over all DOFs
! If the factor is negative then the sign is also changed

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      real*8, intent(in) :: fac

      call VecScalarMult(F,fac,1,F%R())

      end subroutine VecScalarMult_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine VecScalarMult_gen(F,fac,ri,re)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multiplies CP-vec by scalar factor, which is distributed over all DOFs
! If the factor is negative then the sign is also changed

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      integer, intent(in) :: ri,re
      real*8, intent(in)  :: fac
      integer :: ndof

      ndof=SIZE(F%nbas)
      F%coef(ri:re)=abs(fac)*F%coef(ri:re)
      IF (fac.lt.0.d0) call VecSignChange(F,ri,re)

      end subroutine VecScalarMult_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SUMVECVEC(F,Ffac,G,Gfac)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums Ffac*F and Gfac*G in CP-format, with factors distributed over all
! dofs. The summed vector replaces F; G is unchanged.

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP), INTENT(IN)    :: G
      real*8, intent(in) :: Ffac,Gfac
      integer, allocatable :: nbas(:)
      integer :: rF,rG,ndof
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      IF (.NOT.CHECKNBAS(F,G)) THEN
         write(*,*) 'Dimensions or type of F and G do not match'
         write(*,*) 'F:'
         call F%show()
         write(*,*) 'G:'
         call G%show()
         call AbortWithError('Error in SUMVECVEC()')
      ENDIF

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)
      ndof=SIZE(F%nbas)
      ALLOCATE(nbas(ndof))

!     If one of the two terms is zero, it is replaced by the other
      IF (rF.eq.1 .and. F%coef(1).eq.0.d0) THEN
         call FlushCP(F)
         F=CopyCP(G)
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

      subroutine SUMLCVEC(F,Q,facs)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sums a linear combination of vectors in Q and stores in F
! Faster than SUMVECVEC for summing many vectors with coefs since only 
! one allocate is performed

      implicit none
      TYPE (CP), INTENT(OUT) :: F
      TYPE (CP), INTENT(IN)  :: Q(:)
      real*8, intent(in) :: facs(:)
      integer :: i,nbloc,rst,nrki
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      nbloc=SIZE(facs)
      IF (SIZE(Q).ne.nbloc) THEN
         write(*,*) 'Mismatch in dimension of Q, facs'
         call AbortWithError("Error in SUMLCVEC()")
      ENDIF

!     Get the rank of F (=sum of ranks in Q)
      nrki=0
      DO i=1,nbloc
         nrki=nrki+SIZE(Q(i)%coef)
      ENDDO

!     Build F from Q and the factors
      rst=0
      F=NewCP(nrki,Q(1)%nbas)
      DO i=1,nbloc
         nrki=SIZE(Q(i)%coef)
         call GenCopyWtoV(F,Q(i),rst+1,rst+nrki,1,nrki)
         call VecScalarMult(F,facs(i),rst+1,rst+nrki)
         rst=rst+nrki
      ENDDO

      call CPU_TIME(t2)
      svv_time=svv_time+t2-t1

      end subroutine SUMLCVEC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ordre_AA(G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts CP-vector H in decreasing order of coefficients, but uses the
! replacement function so that ordre can be called with a single vector

      implicit none
      TYPE (CP), INTENT(INOUT)  :: G
      TYPE (CP) :: F

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
      TYPE (CP), INTENT(IN)  :: G
      TYPE (CP), INTENT(OUT) :: F
      real*8, allocatable :: tabindex(:)
      integer :: rG,i

!     Set parameters       
      rG=G%R() !SIZE(G%coef)

!     List of ordered indices
      allocate(tabindex(rG))
      do i=1,rG
         tabindex(i)=i
      enddo

!     Sort the coefficients into Fcoef
      F=NewCP(G,rG)
      F%coef=G%coef
      call dsort(F%coef,tabindex,rG,-2)

!     Reorder the base to match the coefficients
      do i=1,rG
         F%base(:,i)=G%base(:,int(tabindex(i)))
      enddo

      deallocate(tabindex)

      end subroutine ordre_AB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NORMBASE1(F,d,ovr)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Renormalizes Fbase, which modifies or overwrites (ovr=F/T) Fcoef

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      logical, intent(in) :: ovr
      integer, intent(in) :: d
      real*8  :: norm1D,t1,t2
      integer :: i,rF

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      rF=F%R()

      do i=1,rF
         norm1D=sqrt(abs(dot_product(F%base(F%ibas(d):F%fbas(d),i),&
                                     F%base(F%ibas(d):F%fbas(d),i))))
         F%base(F%ibas(d):F%fbas(d),i)=&
         F%base(F%ibas(d):F%fbas(d),i)/norm1D
         IF (ovr) THEN
            F%coef(i)=norm1D
         ELSE
            F%coef(i)=F%coef(i)*norm1D
         ENDIF
      enddo

      call CPU_TIME(t2)
      norm_time=norm_time+t2-t1

      end subroutine NORMBASE1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NORMBASE(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Renormalizes Fbase, which modifies Fcoef

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
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
! Normalizes coefs of F. Base should be normalized prior to calling this

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      real*8  :: normND

      normND=1/sqrt(abs(PRODVV(F)))
      F%coef=F%coef*normND

      end subroutine NORMCOEF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine NORMALIZE(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Renormalizes Fbase,Fcoef with calls to NORMBASE,NORMCOEF, respectively

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      real*8  :: normND

      call NORMBASE(F)
      call NORMCOEF(F)

      end subroutine NORMALIZE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateFfromSoln(F,RHS,mode)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fills terms of F for given mode from RHS. New terms are normalized to
! generate new coefficients

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      real*8, intent(in)  :: RHS(:,:)
      integer, intent(in) :: mode
      integer :: i,k,l,row,col
      real*8, parameter :: thresh=1.d-12
      real*8  :: rnorm

      IF (mode.lt.1 .or. mode.gt.SIZE(F%nbas)) THEN
         write(*,*) 'mode (',mode,') must be in [1,',SIZE(F%nbas),']'
         call AbortWithError('UpdateFfromSoln(): mode out of range')
      ENDIF

!     RHS contains multiple solutions from "small" linear system
      IF ((SIZE(RHS,1).eq.SIZE(F%coef)) .and. &
          (SIZE(RHS,2).eq.F%nbas(mode))) THEN
         do i=1,SIZE(F%coef)
            rnorm=sqrt(abs(dot_product(RHS(i,:),RHS(i,:))))
            F%coef(i)=rnorm
            if (rnorm.gt.0.d0) then
               F%base(F%ibas(mode):F%fbas(mode),i)=&
               RHS(i,1:F%nbas(mode))/rnorm
            else
               call RandomNormVec(F%base(F%ibas(mode):F%fbas(mode),i))
            endif
         enddo

!     RHS contains one solution from "large" linear system
!     1st dimension is n vectors of length rF, 2nd dimension is length 1
      ELSEIF ((MOD(SIZE(RHS,1),SIZE(F%coef)).eq.0) .and. & 
              (SIZE(RHS,2).eq.1)) THEN
         do i=1,SIZE(F%coef)
            rnorm=0.d0
            do k=0,F%nbas(mode)-1
               row=k*SIZE(F%coef)+i
               rnorm=rnorm+rhs(row,1)**2
            enddo
            rnorm=sqrt(abs(rnorm))
            F%coef(i)=rnorm
            if (rnorm.gt.0.d0) then
               do k=0,F%nbas(mode)-1
                  row=k*SIZE(F%coef)+i
                  F%base(F%ibas(mode)+k,i)=rhs(row,1)/rnorm
               enddo
            else
               call RandomNormVec(F%base(F%ibas(mode):F%fbas(mode),i))
            endif
         enddo

!     RHS contains multiple solutions from "large" linear system
!     1st dimension is n vectors of length rF, 2nd dimension is length m
      ELSEIF ((MOD(SIZE(RHS,1),SIZE(F%coef)).eq.0) .and. &
              (SIZE(RHS,2).eq.F%cols(mode))) THEN
         do i=1,SIZE(F%coef)
            rnorm=0.d0
            do k=0,F%rows(mode)-1
               row=k*SIZE(F%coef)+i
               do l=1,F%cols(mode)
                  rnorm=rnorm+rhs(row,l)**2
               enddo
            enddo
            rnorm=sqrt(abs(rnorm))
            F%coef(i)=rnorm
            if (rnorm.gt.0.d0) then
               do l=1,F%cols(mode)
                  col=(l-1)*F%rows(mode)
                  do k=0,F%rows(mode)-1
                     row=k*SIZE(F%coef)+i
                     F%base(F%ibas(mode)+col+k,i)=rhs(row,l)/rnorm
                  enddo
               enddo
            else
               call RandomNormVec(F%base(F%ibas(mode):F%fbas(mode),i))
            endif
         enddo

!     Error for incorrect sizes
      ELSE
         write(*,*) 'RHS (',SIZE(RHS,1),' x ',SIZE(RHS,2),') must be ',&
         'either (',SIZE(F%coef),' x ',F%nbas(mode),') or (',&
         SIZE(F%coef)*F%nbas(mode),' x 1) to match F'        
         call AbortWithError('UpdateFfromSoln(): F, RHS size mismatch')

      ENDIF

!     Replace zero coefficients with small values since these may cause
!     zero division in future B,P matrix downdates
      rnorm=thresh*MAXVAL(F%coef)
      do i=1,F%R()
         if (F%coef(i).eq.0.d0) F%coef(i)=rnorm
      enddo

      end subroutine UpdateFfromSoln

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPMatSelfProdTot(M1,M2,sumoverrows) result (F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Product along rows or cols of two CP-matrices
! Returns CP-vector of products as a result
! Set sumoverrows=.TRUE. to sum over rows, .FALSE. to sum over cols

      implicit none
      TYPE (CP), intent(in)  :: M1,M2
      TYPE (CP) :: F
      logical, intent(in)  :: sumoverrows
      integer, allocatable, dimension(:) :: s,nr,nc,mlim,nlim
      integer :: d,i,k,ik,m,ms,n,ns,ndof,rk1,rk2
      real*8  :: prod1d

      IF (.NOT.CHECKNBAS(M1,M2)) THEN
         write(*,*) 'Dimensions or type of M1 and M2 do not match'
         call AbortWithError('Error in CPMatSelfProdTot()')
      ENDIF
      IF (ANY(M1%sym) .or. ANY(M2%sym)) THEN
         write(*,*) 'Symmetry is not yet implemented'
         call AbortWithError('CPMatSelfProdTot(): symmetry not allowed')
      ENDIF

      ndof=SIZE(M1%nbas)
      rk1=SIZE(M1%coef)
      rk2=SIZE(M2%coef)
      ALLOCATE(s(ndof),nr(ndof),nc(ndof),mlim(ndof),nlim(ndof))

      IF (sumoverrows) THEN
         s(:)=M1%rows(:)
         nr(:)=1
         nc(:)=M1%cols(:)
         mlim(:)=M1%cols(:)
         nlim(:)=M1%rows(:)
      ELSE
         s(:)=1
         nr(:)=M1%rows(:)
         nc(:)=1
         mlim(:)=M1%rows(:)
         nlim(:)=M1%cols(:)
      ENDIF

      F=NewCP(rk1*rk2,nr,nc,M1%sym)

      DO d=1,ndof
         ik=1
         DO i=1,rk1
            DO k=1,rk2
               F%coef(ik)=M1%coef(i)*M2%coef(k)
               ms=M1%ibas(d)
               DO m=0,mlim(d)-1
                  prod1D=0.d0
                  ns=ms
                  DO n=0,nlim(d)-1
                     prod1D=prod1D+M1%base(ns,i)*M2%base(ns,k)
                     ns=ns+nr(d)
                  ENDDO
                  F%base(F%ibas(d)+m,ik)=prod1D
                  ms=ms+s(d)
               ENDDO
               ik=ik+1
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(s,nr,nc,mlim,nlim)

      end function CPMatSelfProdTot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function CPHadamardProd(M1,M2) result(M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes element-wise CP-matrix product for M1, M2 with same dims

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP) :: M3
      integer :: i,j,k,ik
      real*8  :: t1,t2

      IF (.NOT. MVV_SETUP) call InitializeMVV()
      call CPU_TIME(t1)

      IF (.NOT.CHECKNBAS(M1,M2)) THEN
         write(*,*) 'Dimensions or type of M1 and M2 do not match'
         write(*,*) 'M1:'
         call M1%show()
         write(*,*) 'M2:'
         call M2%show()
         call AbortWithError('Error in CPHadamardProd()')
      ENDIF

      M3=NewCP(M1,M1%R()*M2%R())

      ik=1
      DO k=1,M1%R()
         DO i=1,M2%R()
            M3%coef(ik)=M1%coef(k)*M2%coef(i)  
            DO j=1,SIZE(M3%base,1)
               M3%base(j,ik)=M1%base(j,k)*M2%base(j,i)
            ENDDO
            ik=ik+1
         ENDDO
      ENDDO

      call CPU_TIME(t2)
      pvv_time=pvv_time+t2-t1

      end function CPHadamardProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function PRODVV_AA(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculates the dot product of F with itself, using factor symmetry

      implicit none
      TYPE (CP), intent(in) :: F
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
      TYPE (CP), intent(in) :: F,G
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

      function CHECKPDIMS(F,G,P)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Checks dimensions of BB or PS inner product matrix to make sure their
! sizes match those of the input CP items

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(in) :: P(:,:)
      logical :: CHECKPDIMS

      CHECKPDIMS=(SIZE(F%coef).eq.SIZE(P,2)) .and. &
                 (SIZE(G%coef).eq.SIZE(P,1))

      end function CHECKPDIMS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTPSTOT(F,G,kskip,PS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSTOT.
! This subroutine calculates PS(l,l')= Pi_{i=1}^ndim <G_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
      
      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in) :: kskip
      integer :: rF,rG,i,j

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')' 
         call AbortWithError('CONSTPSTOT(): dimension mismatch')
      ENDIF

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
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: kskip
      integer :: rF,i,j

      rF=SIZE(F%coef)

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('CONSTBBTOT(): dimension mismatch')
      ENDIF

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

      subroutine CONSTPSk(F,G,k,PS)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PSk must be allocated before the call to CONSTPSk. 
! This subroutine calculates the matrix PSk(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF
! (dividing factor by which to produce the product excluding id = k)
      
      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k
      real*8  :: prod1D
      integer :: rF,rG,i,j,ir,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')'
         call AbortWithError('CONSTPSk(): dimension mismatch')
      ENDIF

!     Compute inner products 
      PS=0.d0
      DO i=1,rG
         DO j=1,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+G%base(ir,i)*F%base(ir,j)
            ENDDO
            PS(i,j)=prod1D
         ENDDO
      ENDDO

      end subroutine CONSTPSk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBk(F,k,BB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BBk must be allocated before the call to CONSTBBk. 
! This subroutine calculates the matrix BBk(l,l')=  < F_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < nrk and 1 < l' < nrk
! (dividing factor by which to produce the product excluding id = k)

      implicit none
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in)   :: k
      real*8  :: prod1D
      integer :: rF,i,j,ir,gi,gf

      rF=size(F%coef)
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('CONSTBBk(): dimension mismatch')
      ENDIF

!     Compute inner products
      BB=0.d0
      DO i=1,rF
         DO j=i,rF
            prod1D=0.d0
            DO ir=gi,gf
               prod1D=prod1D+F%base(ir,i)*F%base(ir,j)
            ENDDO
            BB(i,j)=prod1D
         ENDDO
      ENDDO

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
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
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k
      logical, intent(in)   :: divide
      real*8  :: prod1D
      integer :: rF,rG,i,j,ir,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')'
         call AbortWithError('UpdatePS(): dimension mismatch')
      ENDIF

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
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: k
      logical, intent(in) :: divide
      real*8  :: prod1D
      integer :: rF,i,j,ir,gi,gf

      rF=size(F%coef)
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('UpdateBB(): dimension mismatch')
      ENDIF

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

      subroutine CONSTPSTOT_OMP(F,G,kskip,PS,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSTOT.
! This subroutine calculates PS(l,l')= Pi_{i=1}^ndim <G_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in) :: kskip,pll
      integer :: rF,rG,rT,i,j,l

      rF=SIZE(F%coef)
      rG=SIZE(G%coef)
      rT=rG*rF

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')'
         call AbortWithError('CONSTPSTOT_OMP(): dimension mismatch')
      ENDIF

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

      subroutine CONSTBBTOT_OMP(F,kskip,BB,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BB must be allocated before the call to CONSTBBTOT.
! This subroutine calculates BB(l,l')= Pi_{i=1}^ndim <F_i^l,F_i^l'>
! which appears in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF.
! Different from CONSTPSTOT since BB matrix is symmetric
! Skips the kskip-th DOF (set kskip=0 to include all DOF)
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: kskip,pll
      integer :: rF,rT,i,j,l

      rF=SIZE(F%coef)
      rT=rF*rF

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('CONSTBBTOT_OMP(): dimension mismatch')
      ENDIF

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

      subroutine CONSTPSk_OMP(F,G,k,PS,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PSk must be allocated before the call to CONSTPSk.
! This subroutine calculates the matrix PSk(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF
! (dividing factor by which to produce the product excluding id = k)
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k,pll
      integer :: rF,rG,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      rT=rG*rF
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')'
         call AbortWithError('CONSTPSk_OMP(): dimension mismatch')
      ENDIF

!     Compute inner product for the k-th DOF
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         PS(i,j)=dot_product(G%base(gi:gf,i),F%base(gi:gf,j))
      ENDDO
!$omp enddo
!$omp end parallel

      end subroutine CONSTPSk_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CONSTBBk_OMP(F,k,BB,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BBk must be allocated before the call to CONSTBBk.
! This subroutine calculates the matrix BBk(l,l')=  < F_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < nrk and 1 < l' < nrk
! (dividing factor by which to produce the product excluding id = k)
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in)   :: k,pll
      real*8  :: prod1D
      integer :: rF,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rT=rF*rF
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('CONSTBBk_OMP(): dimension mismatch')
      ENDIF

!     Compute inner product for the k-th DOF
!$omp parallel
!$omp do private(i,j,l)
      DO l=1,rT
         i=(l-1)/rF+1
         j=mod(l-1,rF)+1
         IF (j.ge.i) &
            BB(i,j)=dot_product(F%base(gi:gf,i),F%base(gi:gf,j))
      ENDDO
!$omp enddo
!$omp end parallel

!     Copy elements that are same due to symmetry
      DO i=2,rF
         DO j=1,i-1
            BB(i,j)=BB(j,i)
         ENDDO
      ENDDO

      end subroutine CONSTBBk_OMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdatePS_OMP(F,G,k,PS,divide,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PS must be allocated before the call to CONSTPSk.
! This subroutine modifies the matrix PS(l,l')=  < G_k^l , F_k^l' >
! which occurs in the second part of the normal equations for ALS
! for 1 < l < rG and 1 < l' < rF, by multiplying by or dividing out
! PS for the k-th DOF
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F,G
      real*8, intent(inout) :: PS(:,:)
      integer, intent(in)   :: k,pll
      logical, intent(in)   :: divide
      real*8  :: prod1D
      integer :: rF,rG,rT,i,j,l,gi,gf

      rF=size(F%coef)
      rG=size(G%coef)
      rT=rG*rF
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,G,PS)) THEN
         write(*,*) 'PS is (',SIZE(PS,1),' x ',SIZE(PS,2),&
                    ') but must be (',rG,' x ',rF,')'
         call AbortWithError('UpdatePS_OMP(): dimension mismatch')
      ENDIF

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

      subroutine UpdateBB_OMP(F,k,BB,divide,pll)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine updates matrix BB(l,l')= PI_kp < F_kp^l , F_kp^l' >
! by multiplying by < F_k^l , F_k^l' >, (divide=.FALSE.)
! or dividing by < F_k^l , F_k^l' >, (divide=.TRUE.)
! OpenMP parallel version (pll used to direct interface)

      implicit none
      TYPE (CP), INTENT(IN) :: F
      real*8, intent(inout) :: BB(:,:)
      integer, intent(in) :: k,pll
      logical, intent(in) :: divide
      real*8  :: prod1D
      integer :: rF,rT,i,j,l,gi,gf,is

      rF=size(F%coef)
      rT=rF*rF
      gi=F%ibas(k)
      gf=F%fbas(k)

!     Error checking
      IF (.not.CHECKPDIMS(F,F,BB)) THEN
         write(*,*) 'BB is (',SIZE(BB,1),' x ',SIZE(BB,2),&
                    ') but must be (',rF,' x ',rF,')'
         call AbortWithError('UpdateBB_OMP(): dimension mismatch')
      ENDIF

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
      TYPE (CP), INTENT(IN) :: F
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
      TYPE (CP), INTENT(IN) :: G,F
      TYPE (CP)  :: H
      integer :: rF,rG
      real*8  :: calc_FmG

      rF=SIZE(F%coef)

      H=CopyCP(F)
      call SUMVECVEC(H,1.d0,G,-1.d0)
      calc_FmG=sqrt(abs(PRODVV(H)))
      call FlushCP(H)

      end function calc_FmG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function diagnorm(F)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes sum-of-squares of diagonal elements of a tensor with all 
! dimensions equal

      implicit none
      TYPE (CP), INTENT(IN) :: F
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

      subroutine CPMatVecProd(M,v,w,trans)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matrix-vector product: CP-format matrix times CP-format vector:
! M * v = w. Set trans = .TRUE. to do M^T * v = w. M must be square.

      implicit none
      TYPE (CP), INTENT(IN)  :: M,v
      TYPE (CP), INTENT(OUT) :: w
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
      w=NewCP(rkw,nbas)

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

      end subroutine CPMatVecProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ProjectOut(F,G)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Removes contributions of G from F, returning F.
! No reduction or normalization is done here

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP), INTENT(IN)    :: G
      real*8  :: vivj

      vivj=PRODVV(F,G)
      call SUMVECVEC(F,1.d0,G,-vivj)

      end subroutine ProjectOut

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AugmentVWithRandom(F,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments F to rank nrk by adding random terms with small coefficients

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: G
      real*8, parameter   :: smallnr=1.d-12
      integer, intent(in) :: nrk
      integer :: rF

      rF=SIZE(F%coef)

!     No need to do anything if the rank of F is already large enough
      IF (rF.ge.nrk) RETURN

      G=RandomCP(F,nrk-rF)
      call NORMBASE(G)
      G%coef=MINVAL(F%coef)*smallnr

      call SUMVECVEC(F,1.d0,G,1.d0)

      end subroutine AugmentVWithRandom

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine AugmentVWithRandomN(F,nrk)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Augments F to rank nrk by adding random terms with small coefficients

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP) :: G
      real*8, parameter    :: smallnr=1.d-12
      logical, allocatable :: domode(:)
      integer, intent(in)  :: nrk
      integer :: rF

      rF=SIZE(F%coef)

!     No need to do anything if the rank of F is already large enough
      IF (rF.ge.nrk) RETURN

      G=RandomCP(F,nrk-rF)
      G%base(:,:)=1.d0
      call NORMALIZE(G)
      G%coef=MINVAL(F%coef)*smallnr
      allocate(domode(SIZE(G%nbas)))
      domode(:)=.TRUE.
!      call CPMatrixZeroOffDiag(G,domode)
      deallocate(domode)

      call SUMVECVEC(F,1.d0,G,1.d0)

      end subroutine AugmentVWithRandomN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module MODVECVEC 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
