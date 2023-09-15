!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module ALSUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Utility functions guided by Alternating Least Squares reduction

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE MODVECVEC
      USE LINALG

      implicit none
      real*8, private  :: alsutils_time=0.d0
      logical, private :: ALSUTILS_SETUP=.FALSE.

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SetupALSUtils()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      alsutils_time=0.d0
      ALSUTILS_SETUP=.TRUE.

      end subroutine SetupALSUtils

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposeALSUtils()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Set up the module if it was not set up already
      ALSUTILS_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total reduction time (ALSUTILS)   (s)',&
                            alsutils_time

      end subroutine DisposeALSUtils

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getALSOrthoMem(rG,rF,n,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines memory needed for ALS

      implicit none
      integer, intent(in) :: rG,rF,n,lowmem
      real*8 :: getALSOrthoMem

      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Calculate memory usage
      IF (lowmem.eq.3) THEN
!        Memory allocated for 3 BB, bjk
         getALSOrthoMem=REAL(rF)*(3*rF+n)
      ELSE
!        Memory allocated for PS, max(PS,BB+bjk+IPV)
         getALSOrthoMem=REAL(rF)*(rG+max(rG,rF+n+1))
      ENDIF

      end function getALSOrthoMem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function getALSUpdateMem(rG,rF,n,ncpu,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Determines memory needed for ALS

      implicit none
      integer, intent(in) :: rG,rF,n,ncpu,lowmem
      real*8  :: getALSUpdateMem

      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Calculate memory usage
      IF (lowmem.eq.3) THEN
!        Memory allocated for 2 BB, bjk, IPV
         getALSUpdateMem=REAL(ncpu)*rF*(2*rF+n+1)
      ELSE
!        Memory allocated for PS, 2 BB, bjk, IPV
         getALSUpdateMem=REAL(ncpu)*rF*(rG+2*rF+n+1)
      ENDIF

      end function getALSUpdateMem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_ORTHO_alg(Q,nitn,lowmem,svec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Select the algorithm for intertwined vector updates depending on the 
! value of'lowmem'

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      integer, intent(in)  :: nitn,lowmem
      integer, intent(in), optional :: svec
      integer :: s2

      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Select the first vector in the block to be orthogonalized
      IF (present(svec)) THEN
         s2=svec
      ELSE
         s2=2 ! Orthogonalize beginning with second vector
      ENDIF

      IF (lowmem.ge.0 .and. lowmem.le.2) THEN
         call ALS_ORTHO_2(Q,nitn,s2)
      ELSEIF (lowmem.eq.3) THEN
         call ALS_ORTHO_3(Q,nitn,s2)
      ELSE
         write(*,'(/X,A)') 'ALS_ORTHO_alg(): invalid low-mem choice'
         write(*,'(X,A)')  'Choose: FASTEST <--0=1=2,3--> LEAST MEMORY'
         write(*,'(X,A/)') 'ALS_ORTHO_2() selected by default...'
         call ALS_ORTHO_2(Q,nitn,s2)
      ENDIF

      end subroutine ALS_ORTHO_alg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_SUMLCVEC_alg(F,Q,coefs,nitn,lowmem)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Select the algorithm for intertwined vector updates depending on the 
! value of'lowmem'

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP), INTENT(IN)    :: Q(:)
      real*8, intent(in)   :: coefs(:)
      integer, intent(in)  :: nitn,lowmem

      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

      IF (lowmem.ge.0 .and. lowmem.le.2) THEN
         call ALS_SUMLCVEC_2(F,Q,coefs,nitn)
      ELSEIF (lowmem.eq.3) THEN
         call ALS_SUMLCVEC_3(F,Q,coefs,nitn)
      ELSE
         write(*,'(/X,A)') 'ALS_SUMLCVEC_alg(): invalid low-mem choice'
         write(*,'(X,A)')  'Choose: FASTEST <--0=1=2,3--> LEAST MEMORY'
         write(*,'(X,A/)') 'ALS_SUMLCVEC_2() selected by default...'
         call ALS_SUMLCVEC_2(F,Q,coefs,nitn)
      ENDIF

      end subroutine ALS_SUMLCVEC_alg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_ORTHO_2(Q,nitn,svec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided Gram-Schmidt orthogonalization. Given CP-format vectors
! stored in Q, the vectors are iteratively orthogonalized using ALS.
! This version updates the overlap coefficients after each ALS iteration
! Vectors before svec must be already orthonormalized before calling

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      integer, intent(in)  :: nitn,svec
      integer, allocatable :: nbas(:),ranks(:,:)
      real*8, allocatable  :: coefs(:),BB(:,:),PS(:,:),bjk(:,:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: iv,nbloc,rG,rF,rT,ndim,i,j,ir,imod,k,l,m,n,gst,itn
      integer :: nr
      character(60) :: frmt

      IF (nitn.eq.0) return
      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Set parameters
      nbloc=SIZE(Q)
      ndim=SIZE(Q(1)%nbas)

!     Error checking
      IF (svec.lt.2 .or. svec.gt.nbloc) THEN
         write(*,'(X,2(A,I0),A)') &
         'svec is ',svec,' but must be [2,',nbloc,']'
         call AbortWithError('ALS_ORTHO_2(): svec out of range')
      ENDIF

!     Store a pointer for the ranks of the vectors
      ALLOCATE(coefs(nbloc),ranks(nbloc,3))
      ranks(1,1)=1
      ranks(1,2)=SIZE(Q(1)%coef)
      ranks(1,3)=SIZE(Q(1)%coef)

!     Normalize the first vector
      call NORMALIZE(Q(1))

!     Loop over vectors in the block
      DO iv=2,nbloc

!        Set ranks. rG is the sum of the ranks of the first iv vectors
         coefs(iv)=1.d0
         ranks(iv,1)=ranks(iv-1,2)+1
         ranks(iv,2)=ranks(iv-1,2)+SIZE(Q(iv)%coef)
         ranks(iv,3)=SIZE(Q(iv)%coef)
         rG=ranks(iv,2)
         rF=ranks(iv,3)

         IF (iv.lt.svec) CYCLE

!        Penalty to avoid bad conditioning
         valpen=maxval(Q(iv)%coef)*1.d-10

         ALLOCATE(PS(rG,rF))

!        PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!$omp parallel
!$omp do private(l,t1,t2)
         DO l=1,iv
            call CPU_TIME(t1)            
            call CONSTPT(Q(iv),Q(l),0,PS(ranks(l,1):ranks(l,2),1:rF))
            call CPU_TIME(t2)
            alsutils_time=alsutils_time+t2-t1
         ENDDO
!$omp enddo
!$omp end parallel

!        Main loop over ALS iterations
         DO itn=1,nitn

!           Loop over dimension k
            gst=0
            DO k=1,ndim
               call CPU_TIME(t1)

               n=Q(iv)%nbas(k)

!              Compute the Gram-Schmidt expansion coefficients
!$omp parallel
!$omp do private(l,rT,i,j,t1,t2)
               DO l=1,iv-1
                  call CPU_TIME(t1)
                  coefs(l)=0.d0
                  rT=ranks(l,3)
                  DO i=1,rF
                     DO j=1,rT
                        coefs(l)=coefs(l)-PS(ranks(l,1)+j-1,i)*&
                                 Q(iv)%coef(i)*Q(l)%coef(j)
                     ENDDO
                  ENDDO
                  call CPU_TIME(t2)
                  alsutils_time=alsutils_time+t2-t1
               ENDDO
!$omp enddo
!$omp end parallel

!              Downdate PS (remove the k-th DOF)
               ALLOCATE(BB(rG,rF))
!$omp parallel
!$omp do private(l,rT,t1,t2)
               DO l=1,iv
                  call CPU_TIME(t1)
                  rT=ranks(l,3)
                  call CONSTPk(Q(iv),Q(l),k,&
                                BB(ranks(l,1):ranks(l,2),1:rF))
                  PS(ranks(l,1):ranks(l,2),1:rF)= &
                  PS(ranks(l,1):ranks(l,2),1:rF)/ &
                  BB(ranks(l,1):ranks(l,2),1:rF)
                  call CPU_TIME(t2)
                  alsutils_time=alsutils_time+t2-t1
               ENDDO
!$omp enddo
!$omp end parallel
               call CPU_TIME(t1)
               DEALLOCATE(BB)

!              Copy the last rF x rF block of PS into BB
               ALLOCATE(BB(rF,rF))
               BB(1:rF,1:rF)=PS(ranks(iv,1):ranks(iv,2),1:rF)

!              Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4). This also
!              requires updating the Gram-Schmidt coefficients
               ALLOCATE(bjk(rF,n))
               nr=n*rF
               call CPU_TIME(t2)
               alsutils_time=alsutils_time+t2-t1
!$omp parallel
!$omp do private(m,i,ir,imod,l,rT,j,t1,t2)
               DO m=1,nr
                  call CPU_TIME(t1)
                  i=(m-1)/n+1
                  ir=mod(m-1,n)+1
                  imod=gst+ir
                  bjk(i,ir)=0.d0
                  DO l=1,iv
                     rT=ranks(l,3)
                     DO j=1,rT
                        bjk(i,ir)=bjk(i,ir)+coefs(l)*Q(l)%coef(j)*&
                           Q(l)%base(imod,j)*PS(ranks(l,1)+j-1,i)
                     ENDDO
                  ENDDO
                  call CPU_TIME(t2)
                  alsutils_time=alsutils_time+t2-t1
               ENDDO
!$omp enddo
!$omp end parallel
               call CPU_TIME(t1)

!              Solve linear system B*c_j_k = b_j_k (eq 3.5)
!              (B includes all inner products except the kth)
               call SolveLinSysLU(BB,bjk,valpen)
               DEALLOCATE(BB)

!              Construct improved F
               call UpdateFfromSoln(Q(iv),bjk,k)
               DEALLOCATE(bjk)

!              Check coefs for NaN values resulting from zero division.
               IF (.NOT. CHECKCOEFS(Q(iv))) THEN
                  write(*,*) 'ALS_ORTHO_2(): NaN on update; itn = ',itn,&
                             '; mode = ',k
                  call AbortWithError('ALS_ORTHO_2 crashed')
               ENDIF

!              Update PS (calc inner products with new fs for k-th DOF)
               ALLOCATE(BB(rG,rF))
               call CPU_TIME(t2)
               alsutils_time=alsutils_time+t2-t1

!$omp parallel
!$omp do private(l,rT,t1,t2)
               DO l=1,iv
                  call CPU_TIME(t1)
                  rT=ranks(l,3)
                  call CONSTPk(Q(iv),Q(l),k,&
                                BB(ranks(l,1):ranks(l,2),1:rF))
                  PS(ranks(l,1):ranks(l,2),1:rF)= &
                  PS(ranks(l,1):ranks(l,2),1:rF)* &
                  BB(ranks(l,1):ranks(l,2),1:rF)
                  call CPU_TIME(t2)
                  alsutils_time=alsutils_time+t2-t1
               ENDDO
!$omp enddo
!$omp end parallel
               call CPU_TIME(t1)
               DEALLOCATE(BB)

!              Normalization
               rnorm=0.d0
               DO i=1,rF
                  rnorm=rnorm+PS(ranks(iv,1)+i-1,i)*Q(iv)%coef(i)**2
                  DO j=1,i-1
                     rnorm=rnorm+2*PS(ranks(iv,1)+j-1,i)*&
                              Q(iv)%coef(i)*Q(iv)%coef(j)
                  ENDDO
               ENDDO
               rnorm=1/sqrt(abs(rnorm))
               Q(iv)%coef=rnorm*Q(iv)%coef

               call CPU_TIME(t2)
               alsutils_time=alsutils_time+t2-t1

               gst=gst+n
            ENDDO  ! loop over k
         ENDDO  ! loop over iterations
         DEALLOCATE(PS)
      ENDDO  ! loop over vectors

      DEALLOCATE(ranks,coefs)

      end subroutine ALS_ORTHO_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_ORTHO_3(Q,nitn,svec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided Gram-Schmidt orthogonalization. Given CP-format vectors
! stored in Q, the vectors are iteratively orthogonalized using ALS.
! This version updates the overlap coefficients after each ALS iteration
! and does not store a full-size (nbloc*rF x rF) PS matrix
! Vectors before svec must be already orthonormalized before calling

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      integer, intent(in)  :: nitn,svec
      integer, allocatable :: nbas(:),ind(:,:)
      real*8, allocatable  :: coefs(:)
      real*8, allocatable  :: BB(:,:),BBk(:,:),BBmem(:,:),bjk(:,:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: iv,nbloc,rF,rT,ndim,i,j,k,l,m,n,itn,gik,gfk

      IF (nitn.eq.0) return
      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Set parameters
      nbloc=SIZE(Q)
      ndim=SIZE(Q(1)%nbas)

!     Error checking
      IF (svec.lt.2 .or. svec.gt.nbloc) THEN
         write(*,'(X,2(A,I0),A)') &
         'svec is ',svec,' but must be [2,',nbloc,']'
         call AbortWithError('ALS_ORTHO_3(): svec out of range')
      ENDIF

!     Normalize the first vector
      call NORMALIZE(Q(1))

      call CPU_TIME(t1)

!     Store pointers for the basis indices for each DOF 
      ALLOCATE(ind(ndim,2),coefs(nbloc))
      ind(1,1)=1
      ind(1,2)=Q(1)%nbas(1)
      DO m=2,ndim
         ind(m,1)=ind(m-1,2)+1
         ind(m,2)=ind(m-1,2)+Q(1)%nbas(m)
      ENDDO

!     Loop over vectors in the block
      DO iv=svec,nbloc

!        Penalty to avoid bad conditioning
         valpen=maxval(Q(iv)%coef)*1.d-10

         rF=SIZE(Q(iv)%coef)
         ALLOCATE(BBmem(rF,rF))

!        Main loop over ALS iterations
         DO itn=1,nitn

!           Loop over dimension k
            DO k=1,ndim
               n=Q(iv)%nbas(k)
               gik=ind(k,1)
               gfk=ind(k,2)

!              Calculate b_j_k (l',n) (Beylkin, eq. 3.4). This also
!              requires computing the Gram-Schmidt coefficients
               ALLOCATE(bjk(rF,n))
               bjk=0.d0

!              Loop over vectors in block
               DO l=1,iv

                  rT=SIZE(Q(l)%coef)
                  coefs(l)=0.d0

!                 Calc <Q(iv),Q(l)> for all DOFs except the k-th
                  ALLOCATE(BB(rT,rF))
                  call CONSTPT(Q(iv),Q(l),k,BB,0)

!                 Compute the Gram-Schmidt coefficient. This requires
!                 evaluating the full inner product <Q(iv),Q(l)>, so we
!                 need to do this for the k-th DOF also
                  IF (l.lt.iv) THEN
                     ALLOCATE(BBk(rT,rF))
                     call CONSTPk(Q(iv),Q(l),k,BBk,0)
                     DO j=1,rT
                        DO i=1,rF
                           coefs(l)=coefs(l)-BB(j,i)*BBk(j,i)*&
                                   Q(l)%coef(j)*Q(iv)%coef(i)
                        ENDDO
                     ENDDO
                     DEALLOCATE(BBk)
                  ELSE
!                    Q(iv) is normalized already, so its coef is 1
                     coefs(l)=1.d0
                  ENDIF

!                 Loops over terms in each vector
                  DO j=1,rT
!                    Add contribution from l-th vector to bjk
                     gtmp=coefs(l)*Q(l)%coef(j)
                     DO i=1,rF
                        bjk(i,1:n)=bjk(i,1:n)+&
                                   gtmp*BB(j,i)*Q(l)%base(gik:gfk,j)
                     ENDDO
                  ENDDO  ! terms-in-vector

!                 Keep BB for l=iv since BB is needed to solve the
!                 linear system in the next step
                  IF (l.lt.iv) DEALLOCATE(BB)

               ENDDO  ! vectors

!              Store BBmem <- BB for later normalization
               BBmem(1:rF,1:rF)=BB(1:rF,1:rF)

!              Solve linear system B*c_j_k = b_j_k (eq 3.5)
!              (B includes all inner products except the kth)
               call SolveLinSysLU(BB,bjk,valpen)
               DEALLOCATE(BB)

!              Construct improved F
               call UpdateFfromSoln(Q(iv),bjk,k)
               DEALLOCATE(bjk)

!              Check coefs for NaN values resulting from zero division.
               IF (.NOT. CHECKCOEFS(Q(iv))) THEN
                  write(*,*) 'ALS_ORTHO_3(): NaN on update; itn = ',itn,&
                             '; mode = ',k
                  call AbortWithError('ALS_ORTHO_3 crashed')
               ENDIF

!              Update BBmem with k-th DOF and normalize Q(iv)
               call UpdateP(Q(iv),k,BBmem,.FALSE.,0)

               rnorm=0.d0
               DO i=1,rF
                  rnorm=rnorm+BBmem(i,i)*Q(iv)%coef(i)**2
                  DO j=1,i-1
                     rnorm=rnorm+2*BBmem(j,i)*&
                              Q(iv)%coef(i)*Q(iv)%coef(j)
                  ENDDO
               ENDDO
               rnorm=1/sqrt(abs(rnorm))
               Q(iv)%coef=rnorm*Q(iv)%coef

            ENDDO  ! loop over k
         ENDDO  ! loop over iterations
         DEALLOCATE(BBmem)
      ENDDO  ! loop over vectors

      DEALLOCATE(coefs)

      call CPU_TIME(t2)
      alsutils_time=alsutils_time+t2-t1

      end subroutine ALS_ORTHO_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_ORTHOb(Q,nitn,svec)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ALS guided Gram-Schmidt orthogonalization. Given CP-format vectors
! stored in Q, the vectors are iteratively orthogonalized using ALS.
! This version does not update the overlap coefficients
! Vectors before svec must be already orthonormalized before calling

      implicit none
      TYPE (CP), INTENT(INOUT) :: Q(:)
      integer, intent(in)  :: nitn,svec
      integer, allocatable :: nbas(:),ranks(:,:)
      real*8, allocatable  :: coefs(:),BB(:,:),PS(:,:),bjk(:,:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: iv,nbloc,rG,rF,rT,ndim,i,j,ir,imod,k,l,m,n,gst,itn
      integer :: nr

      IF (nitn.eq.0) return
      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

!     Set parameters
      nbloc=SIZE(Q)
      ndim=SIZE(Q(1)%nbas)

!     Error checking
      IF (svec.lt.2 .or. svec.gt.nbloc) THEN
         write(*,'(X,2(A,I0),A)') &
         'svec is ',svec,' but must be [2,',nbloc,']'
         call AbortWithError('ALS_ORTHOb(): svec out of range')
      ENDIF

!     Store a pointer for the ranks of the vectors
      ALLOCATE(coefs(nbloc),ranks(nbloc,3))
      ranks(1,1)=1
      ranks(1,2)=SIZE(Q(1)%coef)
      ranks(1,3)=SIZE(Q(1)%coef)

!     Normalize the first vector
      call NORMALIZE(Q(1))

      call CPU_TIME(t1)

!     Loop over vectors in the block
      DO iv=2,nbloc

!        Set ranks. rG is the sum of the ranks of the first iv vectors
         coefs(iv)=1.d0
         ranks(iv,1)=ranks(iv-1,2)+1
         ranks(iv,2)=ranks(iv-1,2)+SIZE(Q(iv)%coef)
         ranks(iv,3)=SIZE(Q(iv)%coef)
         rG=ranks(iv,2)
         rF=ranks(iv,3)

         IF (iv.lt.svec) CYCLE

!        Penalty to avoid bad conditioning
         valpen=maxval(Q(iv)%coef)*1.d-10

         ALLOCATE(PS(rG,rF))

!        PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
!$omp parallel
!$omp do private(l)
         DO l=1,iv
            call CONSTPT(Q(iv),Q(l),0,PS(ranks(l,1):ranks(l,2),1:rF))
         ENDDO
!$omp enddo
!$omp end parallel

!        Compute the Gram-Schmidt expansion coefficients
         DO l=1,iv-1
            coefs(l)=0.d0
            rT=ranks(l,3)
            DO i=1,rF
               DO j=1,rT
                  coefs(l)=coefs(l)-PS(ranks(l,1)+j-1,i)*&
                           Q(iv)%coef(i)*Q(l)%coef(j)
               ENDDO
            ENDDO
         ENDDO

!        Main loop over ALS iterations
         DO itn=1,nitn

!           Loop over dimension k
            gst=0
            DO k=1,ndim
               n=Q(iv)%nbas(k)

!              Downdate PS (remove the k-th DOF)
               DO l=1,iv
                  rT=ranks(l,3)
                  ALLOCATE(BB(rT,rF))
                  call CONSTPk(Q(iv),Q(l),k,BB)
                  PS(ranks(l,1):ranks(l,2),1:rF)= &
                  PS(ranks(l,1):ranks(l,2),1:rF)/ &
                  BB(1:rT,1:rF)
                  DEALLOCATE(BB)
               ENDDO

!              Copy the last rF x rF block of PS into BB
               ALLOCATE(BB(rF,rF))
               BB(1:rF,1:rF)=PS(ranks(iv,1):ranks(iv,2),1:rF)

!              Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4). This also
!              requires updating the Gram-Schmidt coefficients
               ALLOCATE(bjk(rF,n))
               bjk=0.d0
               DO ir=1,n
                  imod=gst+ir
                  DO l=1,iv
                     rT=ranks(l,3)
                     DO j=1,rT
                        gtmp=coefs(l)*Q(l)%coef(j)*Q(l)%base(imod,j)
                        DO i=1,rF
                           bjk(i,ir)=bjk(i,ir)+gtmp*PS(ranks(l,1)+j-1,i)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

!              Solve linear system B*c_j_k = b_j_k (eq 3.5)
!              (B includes all inner products except the kth)
               call SolveLinSysLU(BB,bjk,valpen)
               DEALLOCATE(BB)

!              Construct improved F
               call UpdateFfromSoln(Q(iv),bjk,k)
               DEALLOCATE(bjk)

!              Check coefs for NaN values resulting from zero division.
               IF (.NOT. CHECKCOEFS(Q(iv))) THEN
                  write(*,*) 'ALS_ORTHO(): NaN on update; itn = ',itn,&
                             '; mode = ',k
                  call AbortWithError('ALS_ORTHO crashed')
               ENDIF

!              Update PS (calc inner products with new fs for k-th DOF)
               DO l=iv,1,-1
                  rT=ranks(l,3)
                  ALLOCATE(BB(rT,rF))
                  call CONSTPk(Q(iv),Q(l),k,BB)
                  PS(ranks(l,1):ranks(l,2),1:rF)= &
                  PS(ranks(l,1):ranks(l,2),1:rF)* &
                  BB(1:rT,1:rF)
                  DEALLOCATE(BB)
!                 Only update BB (for normalization) on last pass
                  IF (itn.eq.nitn .and. k.eq.ndim) EXIT
               ENDDO

!              Normalization
               rnorm=0.d0
               DO i=1,rF
                  rnorm=rnorm+PS(ranks(iv,1)+i-1,i)*Q(iv)%coef(i)**2
                  DO j=1,i-1
                     rnorm=rnorm+2*PS(ranks(iv,1)+j-1,i)*&
                              Q(iv)%coef(i)*Q(iv)%coef(j)
                  ENDDO
               ENDDO
               rnorm=1/sqrt(abs(rnorm))
               Q(iv)%coef=rnorm*Q(iv)%coef

               gst=gst+n
            ENDDO  ! loop over k
         ENDDO  ! loop over iterations
         DEALLOCATE(PS)
      ENDDO  ! loop over vectors

      DEALLOCATE(ranks,coefs)

      call CPU_TIME(t2)
      alsutils_time=alsutils_time+t2-t1

      end subroutine ALS_ORTHOb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_SUMLCVEC_2(F,Q,coefs,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sum-and-reduce a linear combination of CP-vectors using ALS. This 
! subroutine iteratively refines F, which is a linear combination of
! CP-vectors stored in Q (with coefficients in 'coefs'), without
! generating a long vector

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP), INTENT(IN)    :: Q(:)
      real*8, intent(in)   :: coefs(:)
      integer, intent(in)  :: nitn
      integer, allocatable :: nbas(:),ranks(:,:)
      real*8, allocatable  :: BBmem(:,:),BB(:,:),PS(:,:),bjk(:,:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: nbloc,rG,rF,rT,ndim,i,j,ir,imod,k,l,n,gst,itn

      IF (nitn.eq.0) return
      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

      call CPU_TIME(t1)

!     Set parameters
      nbloc=SIZE(Q)
      ndim=SIZE(Q(1)%nbas)
      rF=SIZE(F%coef)

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     Store a pointer for the ranks of the vectors
      ALLOCATE(ranks(nbloc,3))
      ranks(1,1)=1
      ranks(1,2)=SIZE(Q(1)%coef)
      ranks(1,3)=SIZE(Q(1)%coef)
      DO l=2,nbloc
         ranks(l,1)=ranks(l-1,2)+1
         ranks(l,2)=ranks(l-1,2)+SIZE(Q(l)%coef)
         ranks(l,3)=SIZE(Q(l)%coef)
         rG=ranks(l,2)
      ENDDO
      rG=ranks(nbloc,2)

      ALLOCATE(PS(rG,rF),BBmem(rF,rF))

!     PS(l,l') = Pi_{i=2}^ndim < G_i^l , F_i^l' >
      DO l=1,nbloc
         call CONSTPT(F,Q(l),0,PS(ranks(l,1):ranks(l,2),1:rF))
      ENDDO

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
      call CONSTPT(F,0,BBmem)

!     Loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension k
         gst=0
         DO k=1,ndim
            n=F%nbas(k)

!           Downdate PS (remove the k-th DOF)
            DO l=1,nbloc
               rT=ranks(l,3)
               ALLOCATE(BB(rT,rF))
               call CONSTPk(F,Q(l),k,BB)
               PS(ranks(l,1):ranks(l,2),1:rF)= &
               PS(ranks(l,1):ranks(l,2),1:rF)/ &
               BB(1:rT,1:rF)
               DEALLOCATE(BB)
            ENDDO

!           Downdate BBmem
            call UPDATEP(F,k,BBmem,.TRUE.)
            ALLOCATE(BB(rF,rF))
            BB=BBmem 

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4).
            ALLOCATE(bjk(rF,n))
            bjk=0.d0
            DO ir=1,n
               imod=gst+ir
               DO l=1,nbloc
                  rT=ranks(l,3)
                  DO j=1,rT
                     gtmp=coefs(l)*Q(l)%coef(j)*Q(l)%base(imod,j)
                     DO i=1,rF
                        bjk(i,ir)=bjk(i,ir)+gtmp*PS(ranks(l,1)+j-1,i)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            call SolveLinSysLU(BB,bjk,valpen)
            DEALLOCATE(BB)

!           Construct improved F
            call UpdateFfromSoln(F,bjk,k)
            DEALLOCATE(bjk)

!           Check coefs for NaN values resulting from zero division.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'ALS_SUMLCVEC_2(): NaN on update; itn = ',itn,&
                          '; mode = ',k
               call AbortWithError('ALS_SUMLCVEC_2 crashed')
            ENDIF

!           Update PS with the new Fs (skip on last pass)
            IF (itn.lt.nitn .or. k.lt.ndim) THEN
               DO l=1,nbloc
                  rT=ranks(l,3)
                  ALLOCATE(BB(rT,rF))
                  call CONSTPk(F,Q(l),k,BB)
                  PS(ranks(l,1):ranks(l,2),1:rF)= &
                  PS(ranks(l,1):ranks(l,2),1:rF)* &
                  BB(1:rT,1:rF)
                  DEALLOCATE(BB)
               ENDDO
            ENDIF

!           Update BB (do even on last pass for normalization)
            call UPDATEP(F,k,BBmem,.FALSE.)

            gst=gst+n
         ENDDO  ! loop over k
      ENDDO  ! loop over iterations

!     Normalization
      rnorm=0.d0
      DO i=1,rF
         rnorm=rnorm+BBmem(i,i)*F%coef(i)**2
         DO j=1,i-1
            rnorm=rnorm+2*BBmem(i,j)*F%coef(i)*F%coef(j)
         ENDDO
      ENDDO
      rnorm=1/sqrt(abs(rnorm))
      F%coef=rnorm*F%coef

      DEALLOCATE(PS,BBmem,ranks)

      call CPU_TIME(t2)
      alsutils_time=alsutils_time+t2-t1

      end subroutine ALS_SUMLCVEC_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ALS_SUMLCVEC_3(F,Q,coefs,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sum-and-reduce a linear combination of CP-vectors using ALS. This 
! subroutine iteratively refines F, which is a linear combination of
! CP-vectors stored in Q (with coefficients in 'coefs'), without
! generating a long vector or a PS matrix

      implicit none
      TYPE (CP), INTENT(INOUT) :: F
      TYPE (CP), INTENT(IN)    :: Q(:)
      real*8, intent(in)   :: coefs(:)
      integer, intent(in)  :: nitn
      integer, allocatable :: nbas(:),ranks(:,:),ind(:,:)
      real*8, allocatable  :: BBmem(:,:),BB(:,:),bjk(:,:),PS(:)
      real*8  :: valpen,rnorm,gtmp,t1,t2
      integer :: nbloc,rF,rT,ndim,i,j,ir,k,l,m,n,itn
      integer :: mmod,gi,gf,gik,gfk

      IF (nitn.eq.0) return
      IF (.NOT. ALSUTILS_SETUP) call SetupALSUtils()

      call CPU_TIME(t1)

!     Set parameters
      nbloc=SIZE(Q)
      ndim=SIZE(Q(1)%nbas)
      rF=SIZE(F%coef)

!     Penalty to avoid bad conditioning
      valpen=maxval(F%coef)*1.d-10

!     Store pointers for the ranks of the vectors
      ALLOCATE(ranks(nbloc,3))
      ranks(1,1)=1
      ranks(1,2)=SIZE(Q(1)%coef)
      ranks(1,3)=SIZE(Q(1)%coef)
      DO l=2,nbloc
         ranks(l,1)=ranks(l-1,2)+1
         ranks(l,2)=ranks(l-1,2)+SIZE(Q(l)%coef)
         ranks(l,3)=SIZE(Q(l)%coef)
      ENDDO

!     Store pointers for the basis indices for each DOF 
      ALLOCATE(ind(ndim,2))
      ind(1,1)=1
      ind(1,2)=F%nbas(1)
      DO m=2,ndim
         ind(m,1)=ind(m-1,2)+1
         ind(m,2)=ind(m-1,2)+F%nbas(m)
      ENDDO

      ALLOCATE(BBmem(rF,rF),PS(rF))

!     BB(l,l') = Pi_{i=2}^ndim < F_i^l , F_i^l' >
      call CONSTPT(F,0,BBmem)

!     Loop over ALS iterations
      DO itn=1,nitn

!        Loop over dimension k
         DO k=1,ndim
            n=F%nbas(k)
            gik=ind(k,1)
            gfk=ind(k,2)

!           Calculate b_j_k ( l', nr) (Beylkin, eq. 3.4).
            ALLOCATE(bjk(rF,n))
            bjk=0.d0
!           Loop over vectors in block
            DO l=1,nbloc
               rT=ranks(l,3)
!              Loops over terms in each vector
               DO j=1,rT
                  PS(:)=1.d0
!                 Loop over DOFs except the k-th
                  DO m=1,ndim-1
                     mmod=mod(k+m-1,ndim)+1
                     gi=ind(mmod,1)
                     gf=ind(mmod,2)
!                    Loop over terms in F
                     DO i=1,rF
                        rnorm=0.d0
!                       1D dot product
                        DO ir=gi,gf
                           rnorm=rnorm+F%base(ir,i)*Q(l)%base(ir,j)
                        ENDDO
                        PS(i)=PS(i)*rnorm
                     ENDDO
                  ENDDO

!                 Now construct bjk with current row of PS
                  gtmp=coefs(l)*Q(l)%coef(j)
                  DO i=1,rF
                     bjk(i,1:n)=bjk(i,1:n)+&
                                gtmp*PS(i)*Q(l)%base(gik:gfk,j)
                  ENDDO
               ENDDO  ! terms-in-vector
            ENDDO  ! vectors

!           Downdate BBmem
            call UPDATEP(F,k,BBmem,.TRUE.)
            ALLOCATE(BB(rF,rF))
            BB=BBmem

!           Solve linear system B*c_j_k = b_j_k (eq 3.5)
!           (B includes all inner products except the kth)
            call SolveLinSysLU(BB,bjk,valpen)
            DEALLOCATE(BB)

!           Construct improved F
            call UpdateFfromSoln(F,bjk,k)
            DEALLOCATE(bjk)

!           Check coefs for NaN values resulting from zero division.
            IF (.NOT. CHECKCOEFS(F)) THEN
               write(*,*) 'ALS_SUMLCVEC_3(): NaN on update; itn = ',itn,&
                          '; mode = ',k
               call AbortWithError('ALS_SUMLCVEC_3 crashed')
            ENDIF

!           Update BB (do even on last pass for normalization)
            call UPDATEP(F,k,BBmem,.FALSE.)

         ENDDO  ! loop over k
      ENDDO  ! loop over iterations

!     Normalization
      rnorm=0.d0
      DO i=1,rF
         rnorm=rnorm+BBmem(i,i)*F%coef(i)**2
         DO j=1,i-1
            rnorm=rnorm+2*BBmem(i,j)*F%coef(i)*F%coef(j)
         ENDDO
      ENDDO
      rnorm=1/sqrt(abs(rnorm))
      F%coef=rnorm*F%coef

      DEALLOCATE(BBmem,PS,ranks,ind)

      call CPU_TIME(t2)
      alsutils_time=alsutils_time+t2-t1

      end subroutine ALS_SUMLCVEC_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module ALSUTILS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
