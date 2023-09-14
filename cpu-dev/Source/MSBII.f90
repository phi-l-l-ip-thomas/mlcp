!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MSBII

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains the MSBII solver

      USE ERRORTRAP
      USE UTILS
      USE INPUTFIELDS
      USE SEPDREPN
      USE BLOCKUTILS
      USE ALSPOW
      USE ALSUTILS
      USE ALSDRVR
!!!
      USE LINSOLVER
!!!
      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MSBII_solve(Q,H,eigv,cpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main MSBII solver routine

      implicit none
      TYPE (CPpar), INTENT(IN) :: cpp
      TYPE (CP), INTENT(INOUT) :: Q(:)
      TYPE (CP), INTENT(IN)    :: H
      TYPE (CP), ALLOCATABLE :: Qb(:)
      TYPE (CP) :: G
      real*8, allocatable :: QHQ(:,:),S(:,:),shifts(:)
      real*8, intent(inout) :: eigv(:)
      real*8, parameter     :: tol=1.d-15
      real*8, parameter     :: ranfac=1.d-2
      character(len=18)     :: tag
      integer, parameter :: os=2
      integer :: i,j,k,l,nbloc,conv,ind,nsbloc,nscbloc

      nbloc=SIZE(Q)
      nsbloc=nbloc*os
      nscbloc=nbloc*os*cpp%ncycle
      tag='MSBII'

! Q and eigv have the start vectors and their shift energies
! vectors in Q must be rank-1

      write(*,*) '0th stage: ',0,(eigv(j),j=1,nbloc)

!     Error checking
      do j=1,nbloc
         if (Q(j)%R().gt.1) then
            write(*,*) 'Rank of vector ',j,':',Q(j)%R(),' (must be 1)'
            call AbortWithError('MSBII(): guess rank > 1')
         endif
      enddo

!     Set up arrays
      ALLOCATE(Qb(nscbloc),shifts(nscbloc))
      ALLOCATE(QHQ(nscbloc,nscbloc),S(nscbloc,nscbloc))
      QHQ(:,:)=0.d0
      S(:,:)=0.d0

! First stage: get better shifts

      call prepareshifts(eigv,shifts(1:nsbloc))

      write(*,*) '0th shifts: ',0,(shifts(j),j=1,nsbloc)

      do j=1,nbloc
         do k=1,os
            ind=(j-1)*os+k

!            if (k.eq.1) then
               Qb(ind)=CopyCP(Q(j))
               call AugmentVWithRandom(Qb(ind),cpp%psirank)
               IF (cpp%psirank.gt.1) &
               Qb(ind)%coef(2:cpp%psirank)=ranfac
!            else
!               Qb(ind)=RandomCP(Q(j),cpp%psirank)
!            endif
!!! TEST
!            write(*,*) 'Vector before filter: ',ind
!            call Qb(ind)%printvec()
!!!
!           Apply the filter
            do l=1,cpp%npow
               G=CopyCP(Qb(ind))
               conv=ALS_solve(H,Qb(ind),G,cpp%psinals,1,&
                    shifts(ind),tag)
!               call LinSolver_1(H,Qb(ind),G,cpp%psinals,1,&
!                    shifts(ind),.TRUE.)

!                write(*,*) 'just testing linsolver_2 (Espig)'
!                call LinSolver_2(H,Qb(ind),G,cpp%psinals,0,&
!                    shifts(ind),.TRUE.)
!                call AbortWithError("done")

!               write(*,*)
               call NORMALIZE(Qb(ind))
               call FlushCP(G)
            enddo
!!! TEST
!            write(*,*) 'Vector  after filter:',ind
!            call Qb(ind)%printvec()
!!!
         enddo
      enddo

!     Solve the generalized eigenvalue problem & update shifts
      call GetQHQ2a(Qb(1:nsbloc),H,QHQ(1:nsbloc,1:nsbloc))
      call GetOverlaps(Qb(1:nsbloc),S(1:nsbloc,1:nsbloc))
!!! TEST
!     write(*,*) '(1) QHQ:'
!     call ReflectTriangle(QHQ(1:nsbloc,1:nsbloc),.TRUE.)
!     call PrintMatrix(QHQ(1:nsbloc,1:nsbloc))
!     write(*,*) '(1) S:'
!     call ReflectTriangle(S(1:nsbloc,1:nsbloc),.TRUE.)
!     call PrintMatrix(S(1:nsbloc,1:nsbloc))
!!!
      call SolveGenEigval(shifts(1:nsbloc),S(1:nsbloc,1:nsbloc),&
                          QHQ(1:nsbloc,1:nsbloc),'V')
      eigv(1:nbloc)=shifts(1:nbloc)
!!! TEST
!     write(*,*) '(1) eigenvectors:'
!     call ReflectTriangle(QHQ(1:nsbloc,1:nsbloc),.TRUE.)
!     call PrintMatrix(QHQ(1:nsbloc,1:nsbloc))
!     write(*,*) '(1) eigenvaues:'
!     call PrintVector(shifts(1:nsbloc))
!!!

      write(*,*) '1st stage: ',1,(eigv(j),j=1,nbloc)

!      call AbortWithError('Done stage 1')
! Second stage: calc eigenvalues
     
      call prepareshifts(eigv,shifts(1:nsbloc))
!      write(*,*) '1st shifts: ',0,(shifts(j),j=1,nsbloc)
!      write(*,*)

      do i=2,cpp%ncycle
         shifts((i-1)*nsbloc+1:i*nsbloc)=shifts(1:nsbloc)
         do j=1,nbloc
            do k=1,os
               ind=(i-1)*nsbloc+(j-1)*os+k

               Qb(ind)=CopyCP(Q(j))
               call AugmentVWithRandom(Qb(ind),cpp%psirank)
               IF (cpp%psirank.gt.1) &
                  Qb(ind)%coef(2:cpp%psirank)=ranfac
!               Qb(ind)=RandomCP(Q(j),cpp%psirank)

!              Apply the filter
               do l=1,cpp%npow
                  G=CopyCP(Qb(ind))
                  conv=ALS_solve(H,Qb(ind),G,cpp%psinals,1,&
                       shifts(ind),tag)
                  call NORMALIZE(Qb(ind))
!                  write(*,*)
                 call FlushCP(G)
               enddo
            enddo
         enddo
      enddo

!     Solve the generalized eigenvalue problem for the full block
      call GetQHQ2a(Qb,H,QHQ)
      call GetOverlaps(Qb,S)

!!! TEST
!     write(*,*) '(2) QHQ:'
!     call PrintMatrix(QHQ(1:nscbloc,1:nscbloc))
!     write(*,*) '(2) S:'
!     call PrintMatrix(S(1:nscbloc,1:nscbloc))
!!!

      call SolveGenEigval(shifts,S,QHQ,'V')
      eigv(1:nbloc)=shifts(1:nbloc)

!!! TEST
!     write(*,*) '(2) eigenvectors:'
!     call PrintMatrix(QHQ(1:nscbloc,1:nscbloc))
!     write(*,*) '(2) eigenvaues:'
!     call PrintVector(shifts(1:nscbloc))
!!!

      write(*,*) '2nd stage: ',2,(shifts(j),j=1,nbloc)
      call AbortWithError('Done stage 2')

!     Vectors in Q are still rank-1, so augment to rank R here
      call AugmentQWithRandom(Q,cpp%psirank)

!     Update vectors in Q for the lowest nbloc eigenvalues
      do j=1,nbloc
         call ALS_SUMLCVEC_alg(Q(j),Qb,QHQ(:,j),cpp%psinals,2)
      enddo

      end subroutine MSBII_solve

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine prepareshifts(eigv,shifts)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fill shifts array with shifts at evenly spaced intervals between each
! pair of values in eigv

      real*8, intent(inout) :: shifts(:)
      real*8, intent(in) :: eigv(:)
      integer :: mbloc,nbloc,os,j,k
      real*8  :: shift1
      logical :: shiftfirst=.TRUE.

      mbloc=SIZE(shifts)
      nbloc=SIZE(eigv)
      os=mbloc/nbloc

      IF (os.lt.1 .or. mod(mbloc,nbloc).ne.0) THEN
         write(*,'(X,2(A,I0),A)') &
         'mbloc (',mbloc,'must be a multiple of nbloc(',nbloc,')'
         call AbortWithError('prepareshifts(): bad shifts size')
      ENDIF

      IF (shiftfirst) THEN ! Put shifts below lowest eigenvalue
!        Shifts for first eigenvalue
         shift1=abs(eigv(nbloc)-eigv(1))/mbloc
         do k=1,os
            shifts(k)=eigv(1)-shift1*(os-k)
         enddo

!        Shifts for later eigenvalues
         do j=2,nbloc
            shift1=abs(eigv(j)-eigv(j-1))/os
            do k=1,os
               shifts((j-1)*os+k)=eigv(j)-shift1*(os-k)
            enddo
         enddo

      ELSE ! Put shifts above last eigenvalue
!        Shifts for earlier eigenvalues
         do j=1,nbloc-1
            shift1=abs(eigv(j+1)-eigv(j))/os
            do k=1,os
               shifts((j-1)*os+k)=eigv(j)+shift1*(k-1)
            enddo
         enddo

!        Shift for last eigenvalue
         shift1=abs(eigv(nbloc)-eigv(1))/mbloc
         do k=1,os
            shifts(k)=eigv(nbloc)+shift1*(k-1)
         enddo
      ENDIF

      end subroutine prepareshifts

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE MSBII

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
