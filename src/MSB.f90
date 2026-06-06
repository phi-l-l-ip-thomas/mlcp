!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE MSB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Contains subroutines implementing the MSB solver

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
!      USE CHEBLIB
!      USE INPUTCP
!      USE CPr8
!      USE ALS8DRVR
      USE BLOCKUTILS

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MSB_expand_block(Qr,Q,eigvr,eigv,nblocr,nbloc,&
                                    nshift,esig)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Qr(:),Q(:)
      real(kind=8), allocatable, intent(inout) :: eigvr(:),eigv(:)
      integer, intent(inout) :: nblocr,nbloc
      integer, intent(in) :: nshift
      real(kind=8), intent(in) :: esig
      integer :: j,k,k2
      real(kind=8), parameter :: tol=1.d-12

!      write(*,*) 'MSB expand block'

!     Copy Q -> Qr, eigv -> eigvr, nbloc -> nblocr
!      write(*,*) '...first copy'
      IF (ALLOCATED(Qr)) DEALLOCATE(Qr)
      IF (ALLOCATED(eigvr)) DEALLOCATE(eigvr)
      nblocr=nbloc
      ALLOCATE(Qr(nblocr),eigvr(nblocr))
      DO j=1,nblocr
         Qr(j)=CopyCP(Q(j))
      ENDDO
      eigvr(:)=eigv(:)

!     Create expanded Q,eigv with augmented block size
      DEALLOCATE(Q,eigv)

      nbloc=nblocr+nshift
!      write(*,*) '...get shifts'
!      write(*,*) 'Original shifts'
!      do j=1,nblocr
!         write(*,*) j,eigvr(j)
!      enddo

      call MSB_get_shifts(eigvr,eigv,nshift,esig)

!      write(*,*) '...second copy'
      ALLOCATE(Q(nbloc))
      DO j=1,nbloc
!        Copy vec associated with the energy nearest the shift
         k=max(rbisectH(eigvr,eigv(j)),1)
         k2=min(k+1,nblocr)
         if (abs(eigv(j)-eigvr(k2)).lt.abs(eigv(j)-eigvr(k))) k=k2
!         write(*,*) '   copy j -> k',j,k
         if (abs(eigv(j)-eigvr(k)).lt.tol) then
            Q(j)=CopyCP(Qr(k))
         else
            Q(j)=RandomCP(Qr(k))
         endif
      ENDDO

!!!
!      write(*,*) 'Expanded shifts'
!      do j=1,nbloc
!         write(*,*) j,eigv(j)
!      enddo
!      write(*,*) 'MSB expand block done'
!      call AbortWithError('done')
!!!

      end subroutine MSB_expand_block

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MSB_get_shifts(eigvr,eigv,nshift,esig)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the shifts used in MSB by adding shifts at midpoints of
! largest gaps

      implicit none
      real(kind=8), intent(in) :: eigvr(:)
      real(kind=8), allocatable, intent(out) :: eigv(:)
      integer, intent(in) :: nshift
      real(kind=8), intent(in) :: esig
      real(kind=8), allocatable :: gaps(:),tabindex(:),eigt(:)
      real(kind=8) :: gt(2),et(2),gh
      integer :: i,j,p,nblocr,nbloc
      integer :: l

      nblocr=SIZE(eigvr)
      nbloc=nblocr+nshift

      ALLOCATE(eigv(nbloc))

!     Add shifts
      select case(nshift)
      case(0) ! no shifts: just copy array

        eigv(1:nblocr)=eigv(1:nblocr)

      case(1) ! add upper limit

        eigv(1:nblocr)=eigvr(1:nblocr)
        eigv(nblocr+1)=eigvr(nblocr)+abs(esig)

      case(2) ! add upper and lower limit

        eigv(1)=eigvr(1)-abs(esig)
        eigv(2:nblocr+1)=eigvr(1:nblocr)
        eigv(nblocr+2)=eigvr(nblocr)+abs(esig)

      case default ! add midpoint shifts

        eigv(1)=eigvr(1)-abs(esig)
        eigv(2:nblocr+1)=eigvr(1:nblocr)
        eigv(nblocr+2)=eigvr(nblocr)+abs(esig)
        p=nblocr+2

!       Create array of gaps
        ALLOCATE(gaps(nbloc),tabindex(nbloc),eigt(nbloc))
        gaps(:)=0
        do i=2,p
           gaps(i-1)=eigv(i)-eigv(i-1)
        enddo

!       Sort gaps by increasing size
        do i=1,p
           tabindex(i)=i
        enddo
        call dsort(gaps(1:p),tabindex(1:p),p,2)

!       Reorder energies to match sorted gaps
        do i=1,p
           eigt(i)=eigv(int(tabindex(i)))
        enddo

!       Add shifts
        do while (p.lt.nbloc)
!!!
!        write(*,*) 'sorted gaps and start energies'
!        do l=1,p
!           write(*,*) l,gaps(l),eigt(l)
!        enddo
!!!

!          Add the new shift and gap, bisect the current largest gap
           gh=0.5d0*gaps(p)
           eigt(p+1)=eigt(p)+gh

!           write(*,*) '  Add shift at: ',eigt(p+1)

           gaps(p+1)=gh
           gaps(p)=gh

!          Reorder arrays of gaps,energies in ascending order of gap size
           gt(1:2)=gaps(p:p+1)
           et(1:2)=eigt(p:p+1)
           j=rbisectH(gaps(:p-1),gh)+1
           do i=p-1,j,-1
              gaps(i+2)=gaps(i)
              eigt(i+2)=eigt(i)
           enddo
           gaps(j:j+1)=gt(1:2)
           eigt(j:j+1)=et(1:2)
           p=p+1
        enddo

!       Sort final list of shifts in ascending order
        do i=1,nbloc
           tabindex(i)=i
        enddo
        call dsort(eigt,tabindex,nbloc,2)
        eigv(:)=eigt(:)
        DEALLOCATE(gaps,tabindex,eigt)
      end select

      end subroutine MSB_get_shifts

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine MSB_extract_block(Qr,Q,eigvr,eigv,nblocr,nbloc,esig)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assign vectors, eigenvals in Q,eigv as "matches" to those in Qr,eigvr
! and extracts these, returning them in Q,eigv.
! Qr,eigvr are destroyed on exit.

      implicit none
      TYPE (CP), ALLOCATABLE, INTENT(INOUT) :: Qr(:),Q(:)
      real(kind=8), allocatable, intent(inout) :: eigvr(:),eigv(:)
      integer, intent(inout) :: nblocr,nbloc
      integer, allocatable :: matches(:)
      real(kind=8), allocatable :: ix(:)
      real(kind=8), intent(in) :: esig
      integer :: j

!     Assign vecs in Qr to 'most similar' ones in Q
      call MatchVecs(Qr,Q,eigvr,eigv,esig,matches)

      nblocr=SIZE(Qr)
      DO j=1,nblocr
         call ReplaceVwithW(Qr(j),Q(matches(j)))
         eigvr(j)=eigv(matches(j))
      ENDDO
      DEALLOCATE(Q,eigv)

!     Sort energies in ascending order
      ALLOCATE(ix(nblocr))
      DO j=1,nblocr
         ix(j)=j
      ENDDO
      call dsort(eigvr,ix,nblocr,2)

!     Copy retained vectors from Qr,eigvr -> Q,eigv
      nbloc=nblocr
      ALLOCATE(Q(nbloc),eigv(nbloc))
      eigv(:)=eigvr(:)
      DO j=1,nbloc
         Q(j)=CopyCP(Qr(int(ix(j))))
      ENDDO
      DEALLOCATE(Qr,eigvr)
      nblocr=0

      end subroutine MSB_extract_block

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module MSB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
