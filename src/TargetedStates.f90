!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE TARGETEDSTATES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS

      implicit none

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetStatesinWindow(nbloc,evalsND,qns,evals1Dr,nbas,&
                  nmode,nexci,nmtarget,netarget,Etarget,mstate)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts eigenvalues of a direct product separable wavefunction

      implicit none
      integer, intent(inout) :: qns(:,:)
      real*8, intent(inout)  :: evalsND(:)
      integer, intent(in)    :: nbloc
      integer, intent(in)    :: nbas(:)
      integer, intent(in)    :: nmode(:,:),nexci(:,:)
      integer, intent(in)    :: nmtarget(2),netarget(2)
      integer, intent(out)   :: mstate
      real*8, intent(in)     :: Etarget
      real*8, intent(in)     :: evals1Dr(:,:)
      integer, allocatable   :: iranges(:,:),nranges(:),indx(:),jndx(:)
      integer, allocatable   :: nmranges(:,:),neranges(:,:)
      real*8, allocatable    :: Eranges(:,:),evals1D(:,:),Edifs(:)
      real*8, parameter      :: smallnr=1.d-8
      real*8  :: ZPE,Esofar,eps
      integer :: i,j,k,l,ndof,prodND,nguess,ii,jj,ind,istate
      integer :: nmax,imin,imax,nmsofar,nexsofar
      character(len=64) :: frmt

      ndof=SIZE(nbas)
      nmax=SIZE(evals1Dr,2)

      ALLOCATE(Edifs(nbloc))

!     ZPE-corrected evals1D
      ALLOCATE(evals1D(ndof,nmax))
      ZPE=0.d0
      DO i=1,ndof
         DO j=1,nbas(i)
            evals1D(i,j)=evals1Dr(i,j)-evals1Dr(i,1)
         ENDDO
         ZPE=ZPE+evals1Dr(i,1)
      ENDDO

      ALLOCATE(indx(ndof),nranges(ndof),jndx(ndof),iranges(nmax,ndof),Eranges(ndof,2))
      ALLOCATE(nmranges(ndof,2),neranges(ndof,2))

!     Esofar holds sum of energies of modes whose indices vary more
!     slowly than or equal to that of mode i,
      Esofar=0

!     Eranges(i) holds the maximum and minimum possible summed energy
!     values for modes whose indices vary more rapidly than mode i
      Eranges(1,:)=0.d0
      DO i=2,ndof
         !!! Should look at min and max evals1D, which might not
         !!! be same as first and last val? Or should we enforce
         !!! sorting evals1D in ascending order of energy?
         Eranges(i,1)=Eranges(i-1,1)+evals1D(i-1,1)
         Eranges(i,2)=Eranges(i-1,2)+evals1D(i-1,nbas(i-1))
      ENDDO

!     Window initially covers full energy range (plus room for roundoff)
      eps=Eranges(ndof,2)+evals1D(ndof,nbas(ndof))-Eranges(ndof,1)-evals1D(ndof,1)
      eps=eps*(1.d0+smallnr)

!     nmsofar, nmranges, analogous to energy case above, but for nmode
      nmsofar=0
      nmranges(1,:)=0
      DO i=2,ndof
         nmranges(i,1)=nmranges(i-1,1)+minval(nmode(1:nbas(i),i))
         nmranges(i,2)=nmranges(i-1,2)+maxval(nmode(1:nbas(i),i))
      ENDDO

!     nexsofar, neranges, analogous to energy case above, but for nexci
      nexsofar=0
      neranges(1,:)=0
      DO i=2,ndof
         neranges(i,1)=neranges(i-1,1)+minval(nexci(1:nbas(i),i))
         neranges(i,2)=neranges(i-1,2)+maxval(nexci(1:nbas(i),i))
      ENDDO

      ! Initialize loop limits for outermost mode
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
      ind=ndof
      indx(:)=-1
      indx(ndof)=0
      jndx(:)=-1

!     Number of states so far
      mstate=0

      DO
!        Exit when outermost mode exceeds final index value
         IF (ind.eq.ndof .and. indx(ind).eq.nranges(ind)) EXIT

!        Increment this mode and update energy, nmode, nex
         indx(ind)=indx(ind)+1
         jndx(ind)=iranges(indx(ind),ind)
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

            indx(j)=1
            jndx(j)=iranges(indx(j),j)
            Esofar=Esofar+evals1D(j,iranges(indx(j),j))
            nmsofar=nmsofar+nmode(iranges(indx(j),j),j)
            nexsofar=nexsofar+nexci(iranges(indx(j),j),j)

         ENDDO

         IF (nranges(1).gt.0) THEN

!           Add the element to the list of qns
            call InsertArrayItem(evalsND,Esofar,Edifs,&
                                 abs(Esofar-Etarget),&
                                 qns,jndx,istate,mstate)
!           Once nbloc states have been captured, contract the window so
!           that only states closer to the target than the worst state
!           are added
            IF (mstate.eq.nbloc) THEN
               eps=abs(evalsND(mstate)-Etarget)
            ENDIF

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
      ENDDO

!     Add the ZPE back in
      evalsND(:)=evalsND(:)+ZPE

      DEALLOCATE(indx,jndx,iranges,nranges,Eranges,nmranges,neranges,evals1D)

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

      end function rbisectH

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TARGETEDSTATES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
