!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE FEAST8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines implementing the FEAST solver

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE CHEBLIB
      USE INPUTCP
      USE CPr8
      USE ALS8DRVR
      USE BLOCKUTILS

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetWindows_bycluster8(eigv,rclus,windows)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Partitions energy region defined by eigenvalues in 'eigv' into windows
! determined by the cluster radius

      implicit none
      real(kind=8), intent(in) :: eigv(:)
      real(kind=8), intent(in) :: rclus
      real(kind=8), allocatable, intent(out) :: windows(:,:)
      real(kind=8), allocatable :: wtmp(:,:)
      integer :: i,nev,nwind

      if (rclus.lt.0.d0) then
          write(*,*) 'cluster radius = ',rclus,' must not be negative'
          call AbortWithError('GetWindows_bycluster(): bad value rclus')
      endif

      nev=SIZE(eigv)

      allocate(wtmp(nev,2))
      wtmp(1,1)=eigv(1)-rclus
      wtmp(1,2)=eigv(1)+rclus
      nwind=1

!     Check if each eigenvalue is within the cluster radius of the
!     previous; if so then expand the window, otherwise create a new
!     window
      do i=2,nev
         if (eigv(i).gt.wtmp(nwind,2)) then
            wtmp(nwind,2)=0.5*(eigv(i)+eigv(i-1))
            nwind=nwind+1
            wtmp(nwind,1)=0.5*(eigv(i)+eigv(i-1))
         endif
         wtmp(nwind,2)=eigv(i)+rclus
      enddo

!     Resize the window array
      allocate(windows(nwind,2))
      windows(1:nwind,:)=wtmp(1:nwind,:)
      deallocate(wtmp)

      end subroutine GetWindows_bycluster8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FEAST_iterate_CP8(Q,H,W,eigvo,eigv,cpp,nconv,icyc,tiledH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FEAST iteration master routine
!!! UNFINISHED, UNDER CONSTRUCTION AS OF 6/13/2025

      implicit none
      TYPE (CPpar), INTENT(IN)  :: cpp
      TYPE (CP8), INTENT(INOUT) :: Q(:)
      TYPE (CP8), INTENT(IN)    :: H,W
      TYPE (CP8), ALLOCATABLE   :: bigQ(:),littleQ(:)
      TYPE (CP8) :: Qdummy
      integer, intent(inout) :: nconv
      integer, intent(in)    :: icyc
      logical, intent(in)    :: tiledH
      real(kind=8), intent(inout) :: eigv(:)
      real(kind=8), intent(in)    :: eigvo(:)
      real(kind=8), parameter   :: tol=1.d-15
      real(kind=8), allocatable :: windows(:,:)
      real(kind=8), allocatable :: pts(:),wts(:),shifts(:),bigeigv(:)
      integer, allocatable :: statesinwindow(:),exstatesinwindow(:)
      integer, allocatable :: vecindx(:)
      character(len=18) :: tag
      real(kind=8)  :: Eshift,rq,rclus,exfac,avdf(2),val
      integer :: i,j,k,l,m,nbloc,nblocb,sz,os,szmx,nwind
      integer :: thevec,theshift,thewindow,conv

!     Easy exit for zero iterations
      IF (cpp%npow.lt.1) RETURN

      tag='FEAST      cycle: '
      nbloc=SIZE(Q)

!!!   These parameters will eventually be defined in CP.inp
      rclus=250.0
      exfac=1.2d0
      IF (mpirank.eq.mpi_prnt_rank) write(*,*) '--- FEAST ---'
!!!

!     Quadrature points and weights
      call ChebNodes(cpp%npow,pts) !!! Need complex vsn
      allocate(wts(cpp%npow))
      wts(:)=PI/cpp%npow

!     Get the window boundaries
      call GetWindows_bycluster8(eigv(nconv+1:),rclus,windows)
      nwind=SIZE(windows,1)

!     Count the states inside each window
      allocate(statesinwindow(nwind),exstatesinwindow(nwind))
      statesinwindow(:)=0
      do j=nconv+1,nbloc
         thewindow=rbisectL(windows(:,2),eigv(j))
         statesinwindow(thewindow)=statesinwindow(thewindow)+1
      enddo

!     Expand the number of vecs included in the window to allow
!     for missed states
      do k=1,nwind
         exstatesinwindow(k)=ceiling(exfac*statesinwindow(k))
      enddo

!!!   TEST
      IF (mpirank.eq.mpi_prnt_rank) THEN
        write(*,*) 'window boundaries:'
         do j=1,nwind
            write(*,*) '[',windows(j,1),',',windows(j,2),']'
         enddo
      ENDIF
!!!

!     Determine the expanded block size
      nblocb=0
      do k=1,nwind
         nblocb=nblocb+exstatesinwindow(k)
      enddo
      nblocb=nblocb*cpp%npow

!     Create big block of vectors with shifts determined by quadrature
      ALLOCATE(bigQ(nblocb),shifts(nblocb),vecindx(nblocb))

!!!   TEST
      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,*)
         write(*,*) 'Big block: vecs and shifts for ',nblocb,' vecs'
      ENDIF
!!!
      vecindx(:)=0
      l=nconv
      m=0
      do k=1,nwind
         avdf=averdif(windows(k,:))
         do j=1,exstatesinwindow(k)
            if (j.le.statesinwindow(k)) then
               l=l+1
               thevec=l
               val=eigv(thevec)
            else
               thevec=0
               val=0.d0
            endif
            do i=1,cpp%npow
               m=m+1
               vecindx(m)=thevec
               shifts(m)=mapx(-pts(i),avdf,'fin',.FALSE.)
               IF (mpirank.eq.mpi_prnt_rank) &
                  write(*,'(4(I4,X),2(f14.6,X))') &
                  m,i,thevec,k,val,shifts(m)
            enddo
         enddo
      enddo
      write(*,*)

!     Assign vectors to this MPI rank
      call calc_mpi_partition(nblocb,sz,os,szmx)

!     Solve (Eshift*I - H)*F = G for each vector in the block
!!!   Complex solver needed for this section
      call nvtx_start('Solver Iterations')
      DO j=1,sz
         k=os+j
         IF (vecindx(k).gt.0) THEN
            call bigQ(k)%copyfrom(Q(vecindx(k)))
         ELSE
            call bigQ(k)%clonerand(Q(1))
            if (cpp%algo.eq.1) call bigQ(k)%copyintodevice()
            call Normalize_CP8(bigQ(k),cpp%algo)
         ENDIF
!!!      Replace call to solver with call to complex version
!         call ALS_invpow_normal_CP8(H,-1,shifts(k),bigQ(k),cpp%psinals,1,cpp%algo)
!         call ALS_invpow_fast_CP8(H,-1,shifts(k),bigQ(k),cpp%psinals,1,cpp%algo)
         if (cpp%algo.eq.1) call bigQ(k)%updatefromdevice()
      ENDDO

!     If using the H-tiled algorithm, every rank must participate in
!     cycling portions of H. However, if the number of vectors to 
!     iterate doesn't evenly divide the number of MPI ranks, then some
!     ranks will make fewer calls to the cycling routine. In this case,
!     include a additional vector to preserve the number of MPI calls.
      IF (sz.lt.szmx .and. tiledH) THEN
         call Qdummy%copyfrom(Q(1))
!         call ALS_invpown_normal_CP8(H,-1,shifts(k),Qdummy,cpp%psinals,1,cpp%algo)
!         call ALS_invpow_fast_CP8(H,-1,shifts(k),Qdummy,cpp%psinals,1,cpp%algo)
         call Qdummy%flush()
      ENDIF
      call nvtx_stop()

!!! Looks like full diag of bigQ hits 'error in the diagonalization'
!      call Diagonalize_CP8(bigQ,H,bigeigv,cpp%psinals,cpp%algo,tiledH)
!      write(*,*) 'big eigv:'
!      do j=1,nblocb
!         write(*,*) j,bigeigv(j)
!      enddo
!      write(*,*)
!!! Diag of all calcs in window also has 'error in the diagonalization
!      m=0
!      do k=1,nwind
!         ALLOCATE(littleQ(exstatesinwindow(k)*cpp%npow))
!         ALLOCATE(bigeigv(exstatesinwindow(k)*cpp%npow))
!         l=0
!         do j=1,exstatesinwindow(k)
!            do i=1,cpp%npow
!               l=l+1
!               m=m+1
!               call littleQ(l)%copyfrom(bigQ(m))
!            enddo
!         enddo
!         call Diagonalize_CP8(littleQ,H,bigeigv,cpp%psinals,cpp%algo,tiledH)
!         write(*,*) 'little eigv, window ',k,':'
!         do j=1,exstatesinwindow(k)*cpp%npow
!            write(*,*) j,bigeigv(j)
!         enddo
!         write(*,*)
!         DEALLOCATE(littleQ)
!         DEALLOCATE(bigeigv)
!      enddo
!!! Once complex solver is implemented, select states, sum vecs for each
!!! window, and solve generalized eigenvalue problem 
!      m=0
!      do k=1,nwind
!         ALLOCATE(littleQ(exstatesinwindow(k)),bigeigv(exstatesinwindow(k)))
!        Sum filtered vecs and reduce rank
!         do j=1,exstatesinwindow(k)
!            call littleQ(j)%copyfrom() ! copy original Q vec or random vec
!            call Qdummy%sumlccp(bigQ(m+1:m+cpp%npow),wts)
!            write(*,*) 'Qdummy(',j,')'
!            call Qdummy%printvec
!            call littleQ(j)%copyfrom(Q(1))
!            conv=ALS_reduce_CP8(littleQ(j),Qdummy,cpp%psinals,cpp%algo)
!            write(*,*) 'littleQ(',j,')'
!            call littleQ(j)%printvec
!            call Qdummy%flush()
!            m=m+cpp%npow
!         enddo
!         write(*,*) 'diag:'
!        Solve Generalized eigenvalue problem for this window
!         call pGRAMORTHO_CP8(littleQ,0,cpp%psinals,cpp%algo)
!         call Diagonalize_CP8(littleQ,H,bigeigv,cpp%psinals,cpp%algo,tiledH)
!         write(*,*) 'little eigv, window ',k,':'
!         do j=1,exstatesinwindow(k)
!            write(*,*) j,bigeigv(j)
!         enddo
!         write(*,*)
!         DEALLOCATE(littleQ,bigeigv)
!      enddo

      deallocate(bigQ,shifts,vecindx,pts,wts,windows)
      deallocate(statesinwindow,exstatesinwindow)

      end subroutine FEAST_iterate_CP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module FEAST8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
