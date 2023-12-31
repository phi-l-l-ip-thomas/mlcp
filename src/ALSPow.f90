!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module ALSPOW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Power method guided by Alternating Least Squares reduction

      use ERRORTRAP
      use UTILS
      use SEPDREPN
      use MODVECVEC
      use LINALG
      use MODHVEC

      implicit none
      real*8, private  :: alspow_time=0.d0
      logical, private :: ALSPOW_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupALSPow()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      alspow_time=0.d0
      ALSPOW_SETUP=.TRUE.

      end subroutine SetupALSPow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeALSPow()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. ALSPOW_SETUP) call SetupALSPow()

!     Set up the module if it was not set up already
      ALSPOW_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total reduction time (ALSPOW)     (s)',&
                            alspow_time

      end subroutine DisposeALSPow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getALSPOWmem(rG,rF,n,ndof,ncpu,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines memory needed for ALS

      implicit none
      integer, intent(in) :: rG,rF,n,ndof,ncpu,lowmem
      real*8 :: getALSPOWmem

      IF (.NOT. ALSPOW_SETUP) call SetupALSPow()

!     Calculate memory usage depending on algorithm
      IF (lowmem.eq.0) THEN
!        Memory allocated for ndof PS, ndof BB matrices, bjk+IPV
         getALSPOWmem=REAL(ncpu)*rF*(ndof*(rG+rF)+n+1)
      ELSEIF (lowmem.eq.2) THEN
!        Memory allocated for PS, 2 BB matrices, bjk, IPV
         getALSPOWmem=REAL(ncpu)*rF*(rG+2*rF+n+1)
      ELSEIF (lowmem.eq.3) THEN
!        Memory allocated for BB, bjk, IPV
         getALSPOWmem=REAL(ncpu)*rF*(rF+n+1)
      ELSE ! lowmem=1 or invalid choice of lowmem
!        Memory allocated for PS, 2 BB, and either PS or bjk+IPV
         getALSPOWmem=REAL(ncpu)*rF*(rG+2*rF+n+1)
      ENDIF

      end function getALSPOWmem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_POW_alg(H,F,nitn,ishift,Eshift,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Select the intertwined power iteration algorithm depending on the
! value of 'lowmem'

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in) :: nitn,ishift,lowmem
      real*8, intent(in)  :: Eshift

      IF (.NOT. ALSPOW_SETUP) call SetupALSPow()

      IF (lowmem.eq.0) THEN
         call ALS_POW0(H,F,nitn,ishift,Eshift)
      ELSEIF (lowmem.eq.1) THEN
         call ALS_POW1(H,F,nitn,ishift,Eshift)
      ELSEIF (lowmem.eq.2) THEN
         call ALS_POW2(H,F,nitn,ishift,Eshift)
      ELSEIF (lowmem.eq.3.or.lowmem.eq.4) THEN
         call ALS_POW3(H,F,nitn,ishift,Eshift)
      ELSE
         write(*,'(/X,A)') 'ALS_POW_alg(): invalid low-mem choice'
         write(*,'(X,A)')  'Choose: FASTEST <--0,1,2,3--> LEAST MEMORY'
         write(*,'(X,A/)') 'ALS_POW1() selected by default...'
         call ALS_POW1(H,F,nitn,ishift,Eshift)
      ENDIF

      end subroutine ALS_POW_alg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHV_ALS_alg(Fold,F,H,ishift,Eshift,nitn,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Selects the algorithm to compute an ALS guided single power iteration, 
! depending on the value of 'lowmem'

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: Fold
      TYPE (CPvec), INTENT(OUT) :: F
      integer, intent(in) :: nitn,ishift,lowmem
      real*8, intent(in)  :: Eshift

      IF (.NOT. ALSPOW_SETUP) call SetupALSPow()

      IF (lowmem.ge.0 .and. lowmem.le.2) THEN
         call PRODHV_ALS2(Fold,F,H,ishift,Eshift,nitn)
      ELSEIF (lowmem.eq.3) THEN
         call PRODHV_ALS3(Fold,F,H,ishift,Eshift,nitn)
      ELSE
         write(*,'(/X,A)') 'PRODHV_ALS_alg(): invalid low-mem choice'
         write(*,'(X,A)')  'Choose: FASTEST <--0=1=2,3--> LEAST MEMORY'
         write(*,'(X,A/)') 'PRODHV_ALS2() selected by default...'
         call PRODHV_ALS2(Fold,F,H,ishift,Eshift,nitn)
      ENDIF

      end subroutine PRODHV_ALS_alg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_POW0(H,F,nitn,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided power method. Given CP-format Hamiltonian H and vector F, F
! is iteratively refined using the power method (H - Eshift*I)*F. The
! final vector has the same rank as initial vector F. This version
! generates a 'long' G vector and stores up to D rG x rF 'PS' matrices
! (FASTEST, largest memory usage)

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec) :: G
      integer, intent(in)  :: nitn,ishift
      real*8, intent(in)   :: Eshift
      integer, allocatable :: IPV(:)
      real*8, allocatable  :: bjk(:,:),BBmem(:,:)
      real*8, allocatable  :: BB(:,:,:),PS(:,:,:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rG,rF,ndim,i,j,ir,imod,k,l,n,info,gst,itn,kp

      IF (nitn.eq.0) return

!     First, generate G using a matrix-vector product
      call PRODHV(F,G,H,ishift,Eshift)

      call CPU_TIME(t1)

!     Set parameters
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      ndim=SIZE(G%nbas)

      allocate(BB(rF,rF,ndim),PS(rG,rF,ndim))

!     Penalty to avoid bad conditioning
      valpen=maxval(G%coef)*1.d-15

!     BB(l,l',k) = < F_i^l , F_i^l' > for all but 1st DOF
!     PS(l,l',k) = < G_i^l , F_i^l' > for all but 1st DOF
      DO k=2,ndim
         call CONSTBBk(F,k,BB(:,:,k))
         call CONSTPSk(F,G,k,PS(:,:,k))
      ENDDO

!     Main loop over ALS iterations
      kp=0
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           BB = Pi_{i != k} < F_i^l , F_i^l' >
!           PS = Pi_{i != k} < G_i^l , F_i^l' >
            kp=mod(k,ndim)+1
            BB(:,:,k)=BB(:,:,kp)
            PS(:,:,k)=PS(:,:,kp)
            DO l=2,ndim-1
               kp=mod(k+l-1,ndim)+1
               BB(:,:,k)=BB(:,:,k)*BB(:,:,kp)
               PS(:,:,k)=PS(:,:,k)*PS(:,:,kp)
            ENDDO

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i,k)=BB(i,i,k)+valpen
            enddo

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4)
            allocate(bjk(rF,n))
            bjk=0.d0
            do ir=1,n
               imod=gst+ir
               do j=1,rG
                  gtmp=G%coef(j)*G%base(imod,j)
                  do i=1,rF
                     bjk(i,ir)=bjk(i,ir)+gtmp*PS(j,i,k)
                  enddo
               enddo
            enddo

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            allocate(IPV(rF))
            call dgetrf(rF,rF,BB(:,:,k),rF,IPV,info)
            call dgetrs('N',rF,n,BB(:,:,k),rF,IPV,bjk,rF,INFO)
            deallocate(IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

            call CPU_TIME(t2)
            alspow_time=alspow_time+t2-t1

!           Normalize the coefficients of F after each complete mat-vec
!           to avoid overflow, which can occur after several iterations
            IF (k.eq.ndim) THEN
               rnorm=1/sqrt(abs(PRODVV(F)))
               F%coef=F%coef*rnorm
            ENDIF

            call CPU_TIME(t1)

!           Matrix-vector product on the k-th DOF
            IF (itn.lt.nitn .or. k.lt.ndim) &
               call PRODHV1(F,G,H,k,ishift,Eshift)

!           Calculate BB and PS for the k-th DOF
            IF (itn.lt.nitn .or. k.lt.ndim) THEN
               call CONSTBBk(F,k,BB(:,:,k))
               call CONSTPSk(F,G,k,PS(:,:,k))
            ENDIF

            gst=gst+n
            kp=k
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      deallocate(BB,PS)
      call FlushCPvec(G)

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine ALS_POW0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_POW1(H,F,nitn,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided power method. Given CP-format Hamiltonian H and vector F, F
! is iteratively refined using the power method (H - Eshift*I)*F. The
! final vector has the same rank as initial vector F. This version 
! generates a 'long' G vector and stores up to two rG x rF 'PS' matrices

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: F
      TYPE (CPvec) :: G
      integer, intent(in)  :: nitn,ishift
      real*8, intent(in)   :: Eshift
      integer, allocatable :: IPV(:)
      real*8, dimension (:,:), allocatable :: BB,PS,bjk,BBmem
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rG,rF,ndim,i,j,ir,imod,k,l,n,info,gst,itn,kp
      logical :: update

      IF (nitn.eq.0) return

!     First, generate G using a matrix-vector product
      call PRODHV(F,G,H,ishift,Eshift)

      call CPU_TIME(t1)

!     Set parameters
      rG=SIZE(G%coef)
      rF=SIZE(F%coef)
      ndim=SIZE(G%nbas)

!     Update the ALS matrices BB and PS if ndim > 3. For ndim <= 3 it
!     is faster to build BB and PS at each iteration.
!     Building also avoids possible zero division during the update.
      update=.TRUE.
      IF (ndim.le.3)  update=.FALSE.

      allocate(BB(rF,rF),BBmem(rF,rF),PS(rG,rF))

!     Penalty to avoid bad conditioning
      valpen=maxval(G%coef)*1.d-15

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
!     PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!     If the ALS matrices BB and PS are to be updated, initialize
!     them here with the first DOF removed
      IF (update) THEN
         call CONSTBBTOT(F,1,BBmem)
         call CONSTPSTOT(F,G,1,PS)
      ENDIF

!     Main loop over ALS iterations
      kp=0
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Update the BB and PS matrices of the linear system. This
!           also requires copying BBmem <-- BB since LAPACK destroys BB
            IF (update) THEN
               IF (kp.ne.0) then
                  call UpdateBB(F,k,BBmem,.TRUE.)
                  call UpdatePS(F,G,k,PS,.TRUE.)
               ENDIF
               BB=BBmem

!           Alternatively, build BB and PS from scratch
            ELSE
               call CONSTBBTOT(F,k,BB)
               call CONSTPSTOT(F,G,k,PS)
            ENDIF

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4)
            allocate(bjk(rF,n))
            bjk=0.d0
            do ir=1,n
               imod=gst+ir
               do j=1,rG
                  gtmp=G%coef(j)*G%base(imod,j)
                  do i=1,rF
                     bjk(i,ir)=bjk(i,ir)+gtmp*PS(j,i)
                  enddo
               enddo
            enddo

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            ALLOCATE(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            DEALLOCATE(IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

            call CPU_TIME(t2)
            alspow_time=alspow_time+t2-t1

!           Normalize the coefficients of F after each complete mat-vec
!           to avoid overflow, which can occur after several iterations
            IF (k.eq.ndim) THEN
               rnorm=1/sqrt(abs(PRODVV(F)))
               F%coef=F%coef*rnorm
            ENDIF

            call CPU_TIME(t1)

!           Matrix-vector product on the k-th DOF
            call PRODHV1(F,G,H,k,ishift,Eshift)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'ALS_POW1(): NaN on update; itn = ',itn,&
                             '; k = ',k
                  call AbortWithError('ALS_POW1 crashed')
               ENDIF
            ENDDO

!           Update BB and PS with new Fs
            IF (update) THEN
               IF (itn.lt.nitn .or. k.lt.ndim) THEN
                  call UpdateBB(F,k,BBmem,.FALSE.)
                  call UpdatePS(F,G,k,PS,.FALSE.)
               ENDIF
            ENDIF

            gst=gst+n
            kp=k
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      deallocate(BB,BBmem,PS)
      call FlushCPvec(G)

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine ALS_POW1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_POW2(H,F,nitn,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided power method. Given CP-format Hamiltonian H and vector F, F
! is iteratively refined using the power method (H - Eshift*I)*F. The
! final vector has the same rank as initial vector F. This version
! recomputes PSk and bjk from H*F, as needed, to bypass forming the
! "long" G vector, and stores one rG x rF PS matrix

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in)  :: nitn,ishift
      real*8, intent(in)   :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,bjk,BBmem
      integer, allocatable :: nbas(:),IPV(:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rG,rF,ndim,i,j,ir,k,l,n,info,gst,itn,kp

      IF (nitn.eq.0) return

      call CPU_TIME(t1)

!     Set parameters
      ndim=SIZE(F%nbas)
      rF=SIZE(F%coef)
      rG=rF*SIZE(H%opcoef,1)
      IF (ishift.ne.0) rG=rG+rF

      allocate(BB(rF,rF),BBmem(rF,rF),PS(rG,rF))

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
!     PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!     (Initialize with the first DOF removed)
      call CONSTBBTOT(F,1,BBmem)

!     Build PS without storing G
      PS(:,:)=1.d0
      DO k=1,ndim
         call UpdateXPS(F,F,H,PS,k,ishift,Eshift)
      ENDDO

!     Main loop over ALS iterations
      kp=0
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Update the BBmem matrix of the linear system. Copy to BB and
!           pass BB to LAPACK since LAPACK destroys BB
            IF (kp.ne.0) then
               call UpdateBB(F,k,BBmem,.TRUE.)
            ENDIF
            BB=BBmem

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4) and remove
!           inner products for the k-th DOF from PS
            allocate(bjk(rF,n))
            call DowndateXPS(F,F,H,PS,bjk,k,ishift,Eshift)

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            allocate(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            deallocate(IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

            call CPU_TIME(t2)
            alspow_time=alspow_time+t2-t1

!           Normalize the coefficients of F after each complete mat-vec
!           to avoid overflow, which can occur after several iterations
            IF (k.eq.ndim) THEN
               rnorm=1/sqrt(abs(PRODVV(F)))
               F%coef=F%coef*rnorm
            ENDIF

            call CPU_TIME(t1)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'ALS_POW2(): NaN on update; itn = ',itn,&
                             '; k = ',k
                  call AbortWithError('ALS_POW2 crashed')
               ENDIF
            ENDDO

!           Update BB and PS with the new Fs
            IF (itn.lt.nitn .or. k.lt.ndim) THEN
               call UpdateBB(F,k,BBmem,.FALSE.)
               call UpdateXPS(F,F,H,PS,k,ishift,Eshift)
            ENDIF

            gst=gst+n
            kp=k
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      deallocate(BB,BBmem,PS)

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine ALS_POW2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_POW3(H,F,nitn,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided power method. Given CP-format Hamiltonian H and vector F, F
! is iteratively refined using the power method (H - Eshift*I)*F. The
! final vector has the same rank as initial vector F. This version
! recomputes bjk from H*F on each iteration, and avoids both
! constructing the PS matrix and the 'long' G vector. BB is computed
! along with bjk (SLOWEST, least memory usage)

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(INOUT) :: F
      integer, intent(in)  :: nitn,ishift
      real*8, intent(in)   :: Eshift
      real*8, allocatable  :: BB(:,:),bjk(:,:)
      integer, allocatable :: nbas(:),IPV(:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rF,ndim,i,j,k,ir,l,n,info,gst,itn

      IF (nitn.eq.0) return

      call CPU_TIME(t1)

!     Set parameters
      ndim=SIZE(F%nbas)
      rF=SIZE(F%coef)

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4) and BB
            call BuildbjkDirectX(F,F,H,BB,bjk,k,ishift,Eshift,.FALSE.)

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            allocate(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            deallocate(BB,IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

            call CPU_TIME(t2)
            alspow_time=alspow_time+t2-t1

!           Normalize the coefficients of F after each complete mat-vec
!           to avoid overflow, which can occur after several iterations
            IF (k.eq.ndim) THEN
               rnorm=1/sqrt(abs(PRODVV(F)))
               F%coef=F%coef*rnorm
            ENDIF

            call CPU_TIME(t1)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, error out
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'ALS_POW3(): NaN on update; itn = ',itn
                  call AbortWithError('ALS_POW3 crashed')
               ENDIF
            ENDDO

            gst=gst+n
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine ALS_POW3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHV_ALS2(Fold,F,H,ishift,Eshift,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided single power iteration. Given CP-format Hamiltonian H and 
! vector F, (H - Eshift*I)*F is computed AND reduced using ALS. The
! final vector has the same rank as initial vector F. This version
! recomputes PSk and bjk from H*F, as needed, to bypass forming the
! "long" G vector, and stores one rG x rF PS matrix.

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: Fold
      TYPE (CPvec), INTENT(OUT) :: F
      integer, intent(in)  :: ishift,nitn
      real*8, intent(in)   :: Eshift
      real*8, dimension (:,:), allocatable :: BB,PS,bjk,BBmem
      integer, allocatable :: nbas(:),IPV(:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rG,rF,ndim,i,ir,k,l,n,info,gst,itn,kp

      IF (nitn.eq.0) return

      call CPU_TIME(t1)

      call CopyWtoV(F,Fold)

!     Set parameters
      ndim=SIZE(F%nbas)
      rF=SIZE(F%coef)
      rG=rF*SIZE(H%opcoef,1)
      IF (ishift.ne.0) rG=rG+rF

      allocate(BB(rF,rF),BBmem(rF,rF),PS(rG,rF))

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
!     PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!     If the ALS matrices BB and PS are to be updated, initialize
!     them here with the first DOF removed
      call CONSTBBTOT(F,1,BBmem)

!     Build PS without storing G
      PS(:,:)=1.d0
      DO k=1,ndim
         call UpdateXPS(Fold,F,H,PS,k,ishift,Eshift)
      ENDDO

!     Main loop over ALS iterations
      kp=0
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Update the BBmem matrix of the linear system. Copy to BB and
!           pass BB to LAPACK since LAPACK destroys BB
            IF (kp.ne.0) then
               call UpdateBB(F,k,BBmem,.TRUE.)
            ENDIF
            BB=BBmem

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4) and remove
!           inner products for the k-th DOF from PS
            allocate(bjk(rF,n))
            call DowndateXPS(Fold,F,H,PS,bjk,k,ishift,Eshift)

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            allocate(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            deallocate(IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'PRODHV_ALS2(): NaN on update; itn = ',itn
                  call AbortWithError('PRODHV_ALS2 crashed')
               ENDIF
            ENDDO

!           Update BB and PS with the new Fs (except on the last iteration)
            IF (itn.lt.nitn .or. k.lt.ndim) THEN
               call UpdateBB(F,k,BBmem,.FALSE.)
               call UpdateXPS(Fold,F,H,PS,k,ishift,Eshift)
            ENDIF

            gst=gst+n
            kp=k
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      deallocate(BB,BBmem,PS)

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine PRODHV_ALS2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHV_ALS3(Fold,F,H,ishift,Eshift,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided single power iteration. Given CP-format Hamiltonian H and 
! vector F, (H - Eshift*I)*F is computed AND reduced using ALS. The
! final vector has the same rank as initial vector F. This version
! recomputes bjk from H*F, as needed, bypassing the "long" G vector and
! the PS matrix (SLOWEST, least memory usage)

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)  :: Fold
      TYPE (CPvec), INTENT(OUT) :: F
      integer, intent(in)  :: nitn,ishift
      real*8, intent(in)   :: Eshift
      real*8, allocatable  :: BB(:,:),bjk(:,:)
      integer, allocatable :: nbas(:),IPV(:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: rF,ndim,i,j,k,ir,l,n,info,gst,itn

      IF (nitn.eq.0) return

      call CPU_TIME(t1)

      call CopyWtoV(F,Fold)

!     Set parameters
      ndim=SIZE(F%nbas)
      rF=SIZE(F%coef)

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     Main loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         do k=1,ndim
            n=F%nbas(k)

!           Calculate BB and b_j_k ( l', nr) (Beylkin, eq. 3.4)
            call BuildbjkDirectX(Fold,F,H,BB,bjk,k,ishift,Eshift,.TRUE.)

!           Add penalty to avoid ill-conditioning (section 3.2, Beylkin)
            do i=1,rF
               BB(i,i)=BB(i,i)+valpen
            enddo

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
!           Use DGETRF for LU factorization + DGETRS to solve system
            allocate(IPV(rF))
            call dgetrf(rF,rF,BB,rF,IPV,info)
            call dgetrs('N',rF,n,BB,rF,IPV,bjk,rF,INFO)
            deallocate(BB,IPV)

!           Construct improved F
            do i=1,rF
               rnorm=sqrt(abs(dot_product(bjk(i,:),bjk(i,:))))
               F%coef(i)=rnorm
               F%base(gst+1:gst+n,i)=bjk(i,1:n)/rnorm
            enddo
            deallocate(bjk)

!           Check coefs of F for NaN values resulting from zero
!           division. If there are any, restart ALS without updating
            DO i=1,rF
               IF (F%coef(i).ne.F%coef(i)) THEN
                  write(*,*) 'PRODHV_ALS3(): NaN on update; itn = ',itn
                  call AbortWithError('PRODHV_ALS3 crashed')
               ENDIF
            ENDDO

            gst=gst+n
         enddo  ! loop over k
      ENDDO  ! loop over iterations

      call CPU_TIME(t2)
      alspow_time=alspow_time+t2-t1

      end subroutine PRODHV_ALS3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine PRODHV1(F,G,H,idof,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product H*F = G, where H is the Hamiltonian and
! F and G are vectors, all of which are in CP-format. This subroutine
! applies the matrix-vector product only to the 'idof'-th DOF

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN)    :: F
      TYPE (CPvec), INTENT(INOUT) :: G
      integer, intent(in) :: idof,ishift
      real*8, intent(in)  :: Eshift
      integer :: i,j,k,ik,rF,rG,rH,rFH,gi,gf

!     Set parameters
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      rFH=rF*rH
      rG=rFH
      IF (ishift.ne.0) rG=rG+rF

      IF (rG.ne.SIZE(G%coef)) THEN
         write(*,*) 'rG = ',rG,'; SIZE(G%coef) = ',SIZE(G%coef)
         call AbortWithError('PRODHV1(): unexpected rank of G')
      ENDIF

!     Find the index range for DOF 'idof'
      gi=1
      DO j=1,idof-1
         gi=gi+F%nbas(j)
      ENDDO
      gf=gi+F%nbas(idof)-1

!     Loop over terms in H
      ik=1
      DO k=1,rH
!        Loop over terms in F
         DO i=1,rF
!           Replace the coefficients since the calling subroutine
!           changes the coefficients after operating on each DOF
            G%coef(ik)=F%coef(i)
            call HVBaseProd(F%base(gi:gf,i),G%base(gi:gf,ik),H,idof,k)
            ik=ik+1
         ENDDO
      ENDDO

!     Energy shift if applicable
      IF (ishift.ne.0) THEN
         G%base(gi:gf,rFH+1:rG)=F%base(gi:gf,1:rF)

!        As above, replace the coefficients and apply energy shift for
!        all DOFs, not just the first
         G%coef(rFH+1:rG)=F%coef(1:rF)
         call VecScalarMult(G,Eshift,rFH+1,rG)
         IF (idof.eq.1) THEN
            IF (ishift.gt.0) THEN  ! (H-E*1)*v
               call VecSignChange(G,rFH+1,rG)
            ELSE                   ! (E*1-H)*v
               call VecSignChange(G,1,rFH)
            ENDIF
         ENDIF
      ENDIF

      end subroutine PRODHV1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine UpdateXPS(Fold,F,H,PS,idof,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product to compute <F,H*Fold>, where H is the
! Hamiltonian and F,Fold are CP vectors, for coordinate idof. The 
! rF x rF*rH matrix of overlaps ('PSk') for coordinate idof is 
! multiplied by PS. 

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: Fold,F
      integer, intent(in)   :: idof,ishift
      real*8, intent(in)    :: Eshift
      real*8, intent(inout) :: PS(:,:)
      real*8, allocatable   :: T(:)
      integer :: i,j,k,ik,l,rF,rG,rH,rFH,gi,gf
      real*8  :: prod

!     Set parameters
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      rFH=rF*rH
      rG=rFH
      IF (ishift.ne.0) rG=rG+rF

!     Error checking
      IF (rF.ne.SIZE(Fold%coef)) THEN
         write(*,*) 'rF = ',rF,'; rFold = ',SIZE(Fold%coef)
         call AbortWithError('UpdateXPS(): rF must equal rFold')
      ENDIF
      IF (rG.ne.SIZE(PS,1)) THEN
         write(*,*) 'rG = ',rG,'; SIZE(PS,1) = ',SIZE(PS,1)
         call AbortWithError('UpdateXPS(): mismatch between rG,PS')
      ENDIF
      IF (rF.ne.SIZE(PS,2)) THEN
         write(*,*) 'rF = ',rF,'; SIZE(PS,2) = ',SIZE(PS,2)
         call AbortWithError('UpdateXPS(): mismatch between rF,PS')
      ENDIF

!     Find the index range for DOF 'idof'
      gi=1
      DO j=1,idof-1
         gi=gi+F%nbas(j)
      ENDDO
      gf=gi+F%nbas(idof)-1

      ALLOCATE(T(F%nbas(idof)))

!     Loop over terms in H
      ik=1
      DO k=1,rH
!        Loop over terms in Fold
         DO i=1,rF
!           Matrix-vector product: H*Fold=G
            call HVBaseProd(Fold%base(gi:gf,i),T,H,idof,k)
!           Dot product <F,H*Fold>
            DO l=1,rF
               prod=0.d0
               DO j=1,F%nbas(idof)
                  prod=prod+F%base(gi+j-1,l)*T(j)
               ENDDO
               PS(ik,l)=PS(ik,l)*prod
            ENDDO
            ik=ik+1
         ENDDO
      ENDDO

!     Energy shift if applicable
      IF (ishift.ne.0) THEN
!        The usual case: (H-E*1)*Fold
         DO i=1,rF
            T(:)=Fold%base(gi:gf,i)
            IF (idof.eq.1) T(:)=-T(:)
            DO l=1,rF
               prod=0.d0
               DO j=1,F%nbas(idof)
                  prod=prod+F%base(gi+j-1,l)*T(j)
               ENDDO
               PS(ik,l)=PS(ik,l)*prod
            ENDDO
            ik=ik+1
         ENDDO

!        The rare case: (E*1-H)*Fold
         IF (idof.eq.1) THEN
            IF (ishift.lt.0) THEN
               PS=-PS
            ENDIF
         ENDIF
      ENDIF

      DEALLOCATE(T)

      end subroutine UpdateXPS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DowndateXPS(Fold,F,H,PS,bjk,idof,ishift,Eshift)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product to compute <F,H*Fold>, where H is the
! Hamiltonian and F,Fold are CP vectors, for coordinate idof. The 
! rF x rF*rH matrix of overlaps ('PSk') for coordinate idof is 
! divided out of PS. bjk is also computed here.

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: Fold,F
      integer, intent(in)   :: idof,ishift
      real*8, intent(in)    :: Eshift
      real*8, intent(inout) :: PS(:,:),bjk(:,:)
      real*8, allocatable   :: T(:)
      integer :: i,j,k,ik,l,rF,rG,rH,rFH,gi,gf
      real*8  :: prod,tmp,tmp2

!     Set parameters
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      rFH=rF*rH
      rG=rFH
      IF (ishift.ne.0) rG=rG+rF

!     Error checking
      IF (rF.ne.SIZE(Fold%coef)) THEN
         write(*,*) 'rF = ',rF,'; rFold = ',SIZE(Fold%coef)
         call AbortWithError('DowndateXPS(): rF must equal rFold')
      ENDIF
      IF (rG.ne.SIZE(PS,1)) THEN
         write(*,*) 'rG = ',rG,'; SIZE(PS,1) = ',SIZE(PS,1)
         call AbortWithError('DowndateXPS(): mismatch between rG,PS')
      ENDIF
      IF (rF.ne.SIZE(PS,2)) THEN
         write(*,*) 'rF = ',rF,'; SIZE(PS,2) = ',SIZE(PS,2)
         call AbortWithError('DowndateXPS(): mismatch between rF,PS')
      ENDIF
      IF (rF.ne.SIZE(bjk,1)) THEN
         write(*,*) 'rF = ',rF,'; SIZE(bjk,1) = ',SIZE(bjk,1)
         call AbortWithError('DowndateXPS(): mismatch between rF,bjk')
      ENDIF
      IF (F%nbas(idof).ne.SIZE(bjk,2)) THEN
         write(*,*) 'n = ',F%nbas(idof),'; SIZE(bjk,2) = ',SIZE(bjk,2)
         call AbortWithError('DowndateXPS(): mismatch between n,bjk')
      ENDIF

!     Find the index range for DOF 'idof'
      gi=1
      DO j=1,idof-1
         gi=gi+F%nbas(j)
      ENDDO
      gf=gi+F%nbas(idof)-1

      ALLOCATE(T(F%nbas(idof)))

      bjk=0.d0

!     Loop over terms in H
      ik=1
      DO k=1,rH
!        Loop over terms in F
         DO i=1,rF
            call HVBaseProd(Fold%base(gi:gf,i),T,H,idof,k)
!           Dot product <F,H*F>
            DO l=1,rF
!              Downdate PS
               prod=0.d0
               DO j=1,F%nbas(idof)
                  prod=prod+F%base(gi+j-1,l)*T(j)
               ENDDO
               PS(ik,l)=PS(ik,l)/prod
!              Build bjk
               tmp=Fold%coef(i)*PS(ik,l)
               DO j=1,F%nbas(idof)
                  bjk(l,j)=bjk(l,j)+tmp*T(j)
               ENDDO
            ENDDO
            ik=ik+1
         ENDDO
      ENDDO

!     Energy shift if applicable
      IF (ishift.ne.0) THEN
!        The usual case: (H-E*1)*v
         DO i=1,rF
            tmp=Eshift*Fold%coef(i)
            T(:)=Fold%base(gi:gf,i)
            IF (idof.eq.1) T(:)=-T(:)
            DO l=1,rF
!              Downdate PS
               prod=0.d0
               DO j=1,F%nbas(idof)
                  prod=prod+F%base(gi+j-1,l)*T(j)
               ENDDO
               PS(ik,l)=PS(ik,l)/prod
!              Build bjk
               tmp2=tmp*PS(ik,l)
               DO j=1,F%nbas(idof)
                  bjk(l,j)=bjk(l,j)+tmp2*T(j)
               ENDDO
            ENDDO
            ik=ik+1
         ENDDO

!        The rare case: (E*1-H)*v
         IF (idof.eq.1) THEN
            IF (ishift.lt.0) THEN
               PS=-PS
               bjk=-bjk
            ENDIF
         ENDIF
      ENDIF

      DEALLOCATE(T)

      end subroutine DowndateXPS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine BuildbjkDirectX(Fold,F,H,BB,bjk,idof,ishift,Eshift,&
                                 recalcB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies matrix-vector product to compute <F,(H-Es)*F>, where H is the
! Hamiltonian and F is a CP vectors This subroutine builds bjk directly 
! by computing elements of PS on the fly, to bypass storing PS.

      implicit none
      TYPE (Hamiltonian), INTENT(IN) :: H
      TYPE (CPvec), INTENT(IN) :: Fold,F
      logical, intent(in)  :: recalcB
      integer, intent(in)  :: idof,ishift
      integer, allocatable :: ind(:,:)
      real*8, intent(in)   :: Eshift
      real*8, allocatable  :: T(:)
      real*8, allocatable, intent(out) :: BB(:,:),bjk(:,:)
      integer :: i,j,k,l,m,mmod,ndof,rF,rG,rH,rFH,gi,gf
      real*8  :: prod1D,tmp,tmp2

!     Set parameters
      ndof=SIZE(F%nbas)
      rF=SIZE(F%coef)     ! rank of F
      rH=SIZE(H%opcoef,1) ! rank of H
      rFH=rF*rH
      rG=rFH
      IF (ishift.ne.0) rG=rG+rF

!     Error checking
      IF (rF.ne.SIZE(Fold%coef)) THEN
         write(*,*) 'rF = ',rF,'; rFold = ',SIZE(Fold%coef)
         call AbortWithError('BuildbjkDirectX(): rF must equal rFold')
      ENDIF

!     Find the index range for DOF 'idof'
      allocate(ind(ndof,2))
      ind(1,1)=1
      ind(1,2)=F%nbas(1)
      DO m=2,ndof
         ind(m,1)=ind(m-1,2)+1
         ind(m,2)=ind(m-1,2)+F%nbas(m)
      ENDDO

      ALLOCATE(BB(rF,rF),bjk(rF,F%nbas(idof)),T(MAXVAL(F%nbas)))
      bjk=0.d0

!     Loop over terms in H
      DO k=1,rH
!        Loop over terms in F
         DO i=1,rF
!           Dot product <F,H*Fold>: build including all DOF except idof
            BB(i,:)=1.d0
            DO m=1,ndof-1
               mmod=mod(idof+m-1,ndof)+1
               gi=ind(mmod,1)
               gf=ind(mmod,2)
               call HVBaseProd(Fold%base(gi:gf,i),&
                               T(1:F%nbas(mmod)),H,mmod,k)
!              Construct a row of PS and store in BB for now
               DO l=1,rF
                  prod1D=0.d0
                  DO j=1,F%nbas(mmod)
                     prod1D=prod1D+F%base(gi+j-1,l)*T(j)
                  ENDDO
                  BB(i,l)=BB(i,l)*prod1D
               ENDDO
            ENDDO

!           Build bjk (start with T for this DOF)
            gi=ind(idof,1)
            gf=ind(idof,2)
            call HVBaseProd(Fold%base(gi:gf,i),&
                            T(1:F%nbas(idof)),H,idof,k)
            DO l=1,rF
               tmp=Fold%coef(i)*BB(i,l)
               DO j=1,F%nbas(idof)
                  bjk(l,j)=bjk(l,j)+tmp*T(j)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!     Add the energy shift term, if applicable
!     The usual case: (H-E*1)*v
      DO i=1,rF
         BB(i,:)=1.d0
         tmp=-Eshift*Fold%coef(i)
!        Build PS including all DOF except idof
         DO m=1,ndof-1
            mmod=mod(idof+m-1,ndof)+1
            gi=ind(mmod,1)
            gf=ind(mmod,2)
            T(1:F%nbas(mmod))=Fold%base(gi:gf,i)
!           Construct a row of PS
            DO l=1,rF
               prod1D=0.d0
               DO j=1,F%nbas(mmod)
                  prod1D=prod1D+F%base(gi+j-1,l)*T(j)
               ENDDO
               BB(i,l)=BB(i,l)*prod1D
            ENDDO
         ENDDO

!        Build bjk (start with T for this DOF)
         IF (ishift.ne.0) THEN
            gi=ind(idof,1)
            gf=ind(idof,2)
            T(1:F%nbas(idof))=Fold%base(gi:gf,i)
            DO l=1,rF
               tmp2=tmp*BB(i,l)
               DO j=1,F%nbas(idof)
                  bjk(l,j)=bjk(l,j)+tmp2*T(j)
               ENDDO
            ENDDO
         ENDIF
      ENDDO

!     The rare case: (E*1-H)*v
      IF (idof.eq.1) THEN
         IF (ishift.lt.0) THEN
             bjk=-bjk
         ENDIF
      ENDIF

!     The preceding steps produced BB = <F,Fold>. Sometimes BB = <F,F>
!     is desired; if so, then compute it here
      IF (recalcB) THEN
         DO i=1,rF
            BB(i,:)=1.d0
!           Build BB including all DOF except idof
            DO m=1,ndof-1
               mmod=mod(idof+m-1,ndof)+1
               gi=ind(mmod,1)
               gf=ind(mmod,2)
               T(1:F%nbas(mmod))=F%base(gi:gf,i)
!              Construct a row of BB (upper triangle only)
               DO l=i,rF
                  prod1D=0.d0
                  DO j=1,F%nbas(mmod)
                     prod1D=prod1D+F%base(gi+j-1,l)*T(j)
                  ENDDO
                  BB(i,l)=BB(i,l)*prod1D
               ENDDO
            ENDDO
         ENDDO

!        Copy elements of BB that are same due to symmetry
         DO i=2,rF
            DO l=1,i-1
               BB(i,l)=BB(l,i)
            ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(ind,T)

      end subroutine BuildbjkDirectX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module ALSPOW

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
