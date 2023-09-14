!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE TARGETEDSTATES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS

      implicit none

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetStatesinWindow(nbloc,evalsND,qns,evals1Dr,nbas,Etarget)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts eigenvalues of a direct product separable wavefunction

      implicit none
      integer, allocatable, intent(out) :: qns(:,:)
      real*8, allocatable, intent(inout):: evalsND(:)
      integer, allocatable, intent(in)  :: nbas(:)
      real*8, allocatable, intent(in)   :: evals1Dr(:,:)
      integer, intent(in)  :: nbloc
      real*8, intent(in) :: Etarget
      integer, allocatable :: iranges(:,:),nranges(:),indx(:)
      real*8, allocatable  :: Eranges(:,:),evals1D(:,:),Edifs(:)
      real*8 :: Esofar,eps
      integer :: i,j,k,l,ndof,prodND,nguess,ii,jj,ind,mstate,istate
      integer :: nmax,imin,imax,nmsofar,nexsofar
      character(len=64) :: frmt
!!! TEMP: nmode and nexci arrays
      integer, allocatable :: nmranges(:,:),neranges(:,:)
      integer, allocatable :: nmode(:,:),nexci(:,:)
      integer :: nmtarget(2),netarget(2)
!!!

      ndof=SIZE(nbas)
      nmax=SIZE(evals1Dr,1)

      IF (nbloc.ne.PRODUCT(nbas)) THEN
         write(*,'(2(A,I0),A)') 'nbloc (',nbloc,&
         ') must equal size of DP basis (',PRODUCT(nbas),')'
      ENDIF

      ALLOCATE(qns(nbloc,ndof),Edifs(nbloc))

      write(frmt,'(A,I0,A)') '(',ndof,'(X,I2),2(X,f16.8))'

!!!   TEMP: nmode and nexci arrays
      ALLOCATE(nmode(nmax,ndof),nexci(nmax,ndof))
      DO j=1,ndof
         DO i=1,nbas(j)
            nmode(i,j)=1
            nexci(i,j)=i-1
         ENDDO
         nmode(1,j)=0
      ENDDO

      nmtarget(1)=0
      nmtarget(2)=4
      netarget(1)=0
      netarget(2)=3
!!!

!     ZPE-corrected evals1D
      ALLOCATE(evals1D(nmax,SIZE(evals1Dr,2)))
      DO i=1,ndof
         DO j=1,nbas(i)
            evals1D(i,j)=evals1Dr(i,j)-evals1Dr(i,1)
         ENDDO
      ENDDO

      write(*,*) 'evals1D:'
      DO i=1,ndof
         DO j=1,nbas(i)
            write(*,*) 'i = ',i,'; j = ',j,&
            '; evals1D(j,i) = ',evals1D(i,j),nmode(j,i),nexci(j,i)
         ENDDO
         write(*,*)
      ENDDO

      write(*,*) 'Initializin...'

      ALLOCATE(indx(ndof),nranges(ndof),iranges(nmax,ndof),Eranges(ndof,2))
      ALLOCATE(nmranges(ndof,2),neranges(ndof,2))
      indx(:)=1
      indx(1)=0

!     Esofar holds sum of energies of modes whose indices vary more
!     slowly than or equal to that of mode i,
      Esofar=0

!     Eranges(i) holds the maximum and minimum possible summed energy
!     values for modes whose indices vary more rapidly than mode i
      Eranges(1,:)=0.d0
      DO i=2,ndof
         Eranges(i,1)=Eranges(i-1,1)+evals1D(i-1,1)
         Eranges(i,2)=Eranges(i-1,2)+evals1D(i-1,nbas(i-1))
      ENDDO

!     Window initially covers full energy range
      eps=Eranges(ndof,2)+evals1D(ndof,2)-Eranges(ndof,1)-evals1D(ndof,1)

!     nmsofar, nmranges, analogous to energy case above, but for nmode
      write(*,*) 'nmranges'
      nmsofar=0
      nmranges(1,:)=0
      DO i=2,ndof
         nmranges(i,1)=nmranges(i-1,1)+minval(nmode(1:nbas(i),i))
         nmranges(i,2)=nmranges(i-1,2)+maxval(nmode(1:nbas(i),i))
         write(*,*) i, nmranges(i,1), nmranges(i,2)
      ENDDO

!     nexsofar, neranges, analogous to energy case above, but for nexci
      write(*,*) 'neranges'
      nexsofar=0
      neranges(1,:)=0
      DO i=2,ndof
         neranges(i,1)=neranges(i-1,1)+minval(nexci(1:nbas(i),i))
         neranges(i,2)=neranges(i-1,2)+maxval(nexci(1:nbas(i),i))
         write(*,*) i, neranges(i,1), neranges(i,2)
      ENDDO

!     Number of states so far
      mstate=0

      ! Initialize loop limits for outermost mode
      ind=ndof
      call GetIndexRanges(evals1D(ndof,:nbas(ndof)),Esofar,&
                          Eranges(ndof,1),Eranges(ndof,2),Etarget,eps,&
                          imin,imax)
      nranges(ndof)=0
      DO i=imin,imax
         if (withinranges(nmode(i,ndof),nmsofar,nmtarget,nmranges(ndof,:))) then
            if (withinranges(nexci(i,ndof),nexsofar,netarget,neranges(ndof,:))) then
               nranges(ndof)=nranges(ndof)+1
               iranges(nranges(ndof),ndof)=i
            endif
         endif
      ENDDO
      indx(ndof)=0

      write(*,*) 'Everything initialized!'

      DO

!        Increment this mode and update energy, nmode, nex
         indx(ind)=indx(ind)+1
         Esofar=Esofar+evals1D(ind,iranges(indx(ind),ind))
         nmsofar=nmsofar+nmode(iranges(indx(ind),ind),ind)
         nexsofar=nexsofar+nexci(iranges(indx(ind),ind),ind)

!        Reset earlier modes
         DO j=ind-1,1,-1
            call GetIndexRanges(evals1D(j,:nbas(j)),Esofar,&
                                Eranges(j,1),Eranges(j,2),Etarget,eps,&
                                imin,imax)
            nranges(j)=0
            DO i=imin,imax
            ! Add nmode, nexec conditions here
               if (withinranges(nmode(i,j),nmsofar,nmtarget,nmranges(j,:))) then
                  if (withinranges(nexci(i,j),nexsofar,netarget,neranges(j,:))) then
                    nranges(j)=nranges(j)+1
                    iranges(nranges(j),j)=i
!                    write(*,*) 'adding',i
                  endif
               endif
            ENDDO
            indx(j)=imin
            Esofar=Esofar+evals1D(j,iranges(indx(j),j))
            nmsofar=nmsofar+nmode(iranges(indx(j),j),j)
            nexsofar=nexsofar+nexci(iranges(indx(j),j),j)

!               if (iranges(j,1).gt.iranges(j,2)) &
!                  write(*,*) 'IRANGES: ', iranges(j,1), ' > ',iranges(j,2),'; mode',j,' !!!'
!               write(*,*) Esofar,' = ',Esofar-evals1D(j,indx(j)),' + ',evals1D(j,indx(j))
         ENDDO

         IF (nranges(1).gt.0) THEN
!           Add the element to the list of qns
            call InsertArrayItem(evalsND,Esofar,Edifs,&
                                 abs(Esofar-Etarget),&
                                 qns,indx,istate,mstate)
!           Once nbloc states have been captured, contract the window so
!           that only states closer to the target than the worst state
!           are added
            IF (mstate.eq.nbloc) THEN
               eps=abs(evalsND(mstate)-Etarget)
            ENDIF
            write(*,frmt) (iranges(indx(i),i)-1,i=1,ndof),Esofar,Esofar-Etarget
               write(*,*) 'mstate = ',mstate,'; eps = ',eps
         ENDIF

!        Find the mode with the index to increment
         DO i=1,ndof
            ind=i
!           This mode: downdate energy, nmode, nexci
            Esofar=Esofar-evals1D(i,iranges(indx(i),i))
            nmsofar=nmsofar-nmode(iranges(indx(i),i),i)
            nexsofar=nexsofar-nexci(iranges(indx(i),i),i)

            IF (indx(i).lt.nranges(i)) EXIT
         ENDDO

!        Exit when outermost mode exceeds final index value
         IF (ind.eq.ndof .and. indx(ind).eq.nranges(ind)) EXIT
      ENDDO

      DEALLOCATE(indx,iranges,nranges,Eranges,nmranges,neranges,evals1D)

      write(*,*) 'Final list:'
      do j=1,nbloc
         write(*,frmt) (qns(j,i)-1,i=1,ndof),evalsND(j),Edifs(j)
      enddo

      call AbortWithError("Done")

      end subroutine GetStatesinWindow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InsertArrayItem(evalsND,Esofar,Edifs,dif,qns,indx,i,m)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Inserts item p into array v, in position to keep v in ascending order. 
! On input, elements in v must be in ascending order.
! v is of length m on entry and min(n,m+1) on exit.
! i is the position of element v

      implicit none
      integer, intent(in)    :: indx(:)
      integer, intent(inout) :: qns(:,:)
      integer, intent(inout) :: m
      integer, intent(out)   :: i
      real*8, intent(in)     :: Esofar,dif
      real*8, intent(inout)  :: evalsND(:),Edifs(:)
      integer :: j,n

      IF (m.lt.1) THEN
         i=1
         Edifs(i)=dif
         evalsND(i)=Esofar
         qns(i,:)=indx(:)
         m=1
         RETURN
      ELSE
         n=SIZE(Edifs)
         i=rbisectL(Edifs(1:m),dif)
         m=min(m+1,n)

         DO j=min(m,n),i+1,-1
            Edifs(j)=Edifs(j-1)
            evalsND(j)=evalsND(j-1)
            qns(j,:)=qns(j-1,:)
         ENDDO
         IF (i.le.n) THEN
            Edifs(i)=dif
            evalsND(i)=Esofar
            qns(i,:)=indx(:)
         ENDIF
      ENDIF

      end subroutine InsertArrayItem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetIndexRanges(evals1D,Esofar,Emin,Emax,Etarget,eps,&
                                ilo,ihi)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      integer, intent(inout) :: ilo, ihi
      real*8, intent(in) :: Esofar, Emin, Emax, Etarget, eps
      real*8, intent(in) :: evals1D(:)
      real*8 :: limitlo, limithi

!     Esofar+Emin+largestallowed <= Etarget+eps 
!     Esofar+Emax+smallestallowed >= Etarget-eps
      limitlo = Etarget-eps-Esofar-Emax
      limithi = Etarget+eps-Esofar-Emin
      ilo=rbisectL(evals1D,limitlo)
      ihi=rbisectH(evals1D,limithi)

      end subroutine GetIndexRanges

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function withinranges(val,sofar,targ,ranges) RESULT(within)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
      logical :: within
      integer, intent(in) :: val,sofar,targ(2),ranges(2)
      integer :: limitlo,limithi

      limitlo = targ(1)-sofar-ranges(2)
      limithi = targ(2)-sofar-ranges(1)

      within=.FALSE.
      IF (val.ge.limitlo .and. val.le.limithi) within=.TRUE.

      end function withinranges

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function rbisectL(v,i) RESULT(ju)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds index of first item of v that is >= i, by bisection
! Data values must be in monotonically-ascending order

      implicit none
      real*8, intent(in) :: v(:)
      real*8, intent(in) :: i
      integer :: n,ju,jl,t

      n=SIZE(v)

      jl=0
      ju=n+1
      IF (i.le.v(1)) THEN
         ju=1
      ELSEIF (i.gt.v(n)) THEN
         jl=n
      ENDIF

      DO
        IF (ju-jl.le.1) EXIT
        t=(ju+jl)/2
        IF (v(t).lt.i) THEN
           jl=t
        ELSE
           ju=t
        ENDIF
      ENDDO

!      write(*,*) 'L: i,jl,ju=',i,jl,ju

      end function rbisectL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function rbisectH(v,i) RESULT(jl)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finds index of last item of v that is >= i, by bisection
! Data values must be in monotonically-ascending order

      implicit none
      real*8, intent(in) :: v(:)
      real*8, intent(in) :: i
      integer :: n,ju,jl,t

      n=SIZE(v)

      jl=0
      ju=n+1
      IF (i.lt.v(1)) THEN
         ju=1
      ELSEIF (i.ge.v(n)) THEN
         jl=n
      ENDIF

      DO
        IF (ju-jl.le.1) EXIT
        t=(ju+jl)/2
        IF (v(t).le.i) THEN
           jl=t
        ELSE
           ju=t
        ENDIF
      ENDDO

!      write(*,*) 'H: i,jl,ju=',i,jl,ju

      end function rbisectH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TARGETEDSTATES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
