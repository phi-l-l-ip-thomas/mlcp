!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module HGORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CP-matrix orthogonalization via Alternating Least Squares

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE MODVECVEC
      USE LINALG
      USE CPMMM
      USE REDUCTION
      USE ALSDRVR
      USE CPMATH

      INTERFACE OrthogCPmat
         MODULE PROCEDURE OrthogCPmat_noweights,OrthogCPmat_weights
      END INTERFACE OrthogCPmat

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine OrthogCPmat_noweights(U,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Orthogonalizes CP-format matrix U using an ALS-type algorithm. This
! version builds the columns and reduces their ranks

      implicit none
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP) :: V,W,Y,Z
      logical, allocatable :: domode(:)
      integer, intent(in)  :: nitn
      integer, parameter   :: nals=10
      integer, allocatable :: irs(:),irf(:),ics(:),icf(:)
      integer :: itn,d,ndof,ic,nc,rk,j,redstat
      integer :: stu,fiu,rki,rke
      real*8  :: norm
      character(len=64) :: tag

      IF (nals.lt.1) return

      ndof=U%D()
      rk=U%R()

      allocate(irs(ndof),irf(ndof),ics(ndof),icf(ndof),domode(ndof))

!     Main loop over ALS iterations      
      DO itn=1,nitn
!        Loop over modes
         DO d=1,ndof

            nc=U%cols(d)
            stu=U%ibas(d)
            fiu=U%ibas(d)+U%rows(d)-1

            irs(:)=1
            irf(:)=U%rows(:)
            ics(:)=1
            icf(:)=U%cols(:)

!           Initialize W as a zero tensor
            W=NewCP(U,1)
            W%coef(:)=0.d0

!           Loop over cols for mode d
            DO ic=1,nc

!              Extract the next col of mode d of U
               ics(d)=1
               icf(d)=1
               V=ExtractCPsubmatrix(U,irs,irf,ics,icf)

!              Normalization: calc norm tensor, mult by col, ALS
               call CPMatNormA(V,nals)

!              Put normalized col ic of U into W as next rk terms
!              Y contains V, resized for W               
               Y=NewCP(W,rk) 
               Y%base(:,:)=0.d0
               ics(d)=ic
               call PutCPsubmatrix(V,irs,ics,Y)
               call FlushCP(V)
               call SUMVECVEC(W,1.d0,Y,1.d0)
               call FlushCP(Y)

!              Exit upon reaching last col since nothing to project out
               IF (ic.eq.nc) EXIT

!              Resize U so projector only operates on cols > ic
               ics(d)=2
               icf(d)=U%cols(d)
               V=ExtractCPsubmatrix(U,irs,irf,ics,icf)
               call ReplaceVwithW(U,V)

!              Project col ic out of other cols using ic-th block of W
!              U = [I - W*W^T]*U (1st step): Y = W^T*U
               rki=(ic-1)*rk+1
               rke=ic*rk
!               call CPMM(W,rki,rke,0,0.d0,.TRUE.,&
               call CPMM(W,1,rke,0,0.d0,.TRUE.,&
                         U,1,rk,0,0.d0,.FALSE.,Y)

!              Projection operator: diag(Y) for modes > d
               domode(:)=.TRUE.
               domode(1:d)=.FALSE.
               call CPMatrixZeroOffDiag(Y,domode)

!              Projection operator: reduce rank Z <-ALS-- Y
               Z=RandomCP(Y,rk)
               write(tag,'(A,3(I0,X))') &
               'OrthogCPmat(): Z<-Y, itn,d,ic = ',itn,d,ic
               redstat=ALS_reduce(Z,Y,nals,tag)
               call ReplaceVwithW(Y,Z)

!              Projection operator (2nd step): V = -W*Y + U
               call CPMM(W,rki,rke,0,0.d0,.FALSE.,&
                         Y,1,Y%R(),0,0.d0,.FALSE.,V)
               call FlushCP(Y)
               call SUMVECVEC(V,-1.d0,U,1.d0)

!              Reduce rank of result, U <-ALS--V
               write(tag,'(A,3(I0,X))') &
               'OrthogCPmat(): U<-V, itn,d,ic = ',itn,d,ic
               redstat=ALS_reduce(U,V,nals,tag)
               call FlushCP(V)

               stu=stu+U%rows(d)
               fiu=fiu+U%rows(d)
            ENDDO ! cols loop (ic)

!           Reduce rank and replace : U <-ALS-- W.
            call FlushCP(U)
            U=RandomCP(W,rk)
            redstat=ALS_reduce(U,W,10*nals,'OrthogCPmat(): U<-W')
            call FlushCP(W)

         ENDDO ! modes loop (d)
      ENDDO ! main ALS loop (itn)

      deallocate(irs,irf,ics,icf,domode)

      end subroutine OrthogCPmat_noweights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine OrthogCPmat_weights(U,Wt,nitn)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Orthogonalizes CP-format matrix U using an ALS-type algorithm. This
! version builds the columns and reduces their ranks

      implicit none
      TYPE (CP), INTENT(IN) :: Wt
      TYPE (CP), INTENT(INOUT) :: U
      TYPE (CP) :: V,W,Y,Z,Wtcol
      logical, allocatable :: domode(:)
      integer, intent(in)  :: nitn
      integer, parameter   :: nals=10
      integer, allocatable :: irs(:),irf(:),ics(:),icf(:)
      integer, allocatable :: wrs(:),wrf(:),wcs(:),wcf(:)
      integer :: itn,d,ndof,ic,nc,rk,j,redstat
      integer :: stu,fiu,rki,rke
      real*8  :: norm
      character(len=64) :: tag

      IF (nals.lt.1) return

      ndof=U%D()
      rk=U%R()

      allocate(irs(ndof),irf(ndof),ics(ndof),icf(ndof),domode(ndof))
      allocate(wrs(ndof),wrf(ndof),wcs(ndof),wcf(ndof))

!     Main loop over ALS iterations      
      DO itn=1,nitn
!        Loop over modes
         DO d=1,ndof

            nc=U%cols(d)
            stu=U%ibas(d)
            fiu=U%ibas(d)+U%rows(d)-1

            irs(:)=1
            irf(:)=U%rows(:)
            ics(:)=1
            icf(:)=U%cols(:)

            wrs(:)=1
            wrf(:)=Wt%rows(:)
            wcs(:)=1
            wcf(:)=Wt%cols(:)

!           Initialize W as a zero tensor
            W=NewCP(U,1)
            W%coef(:)=0.d0

!            write(*,*) 'Wtsub, initial'
!            call Wt%print()

!           Loop over cols for mode d
            DO ic=1,nc

!              Extract the next col of mode d of U
               ics(d)=1
               icf(d)=1
               V=ExtractCPsubmatrix(U,irs,irf,ics,icf)

!              Extract weights for col extracted above
!              (Wt is square, so use col index also for rows)
               wrs(d)=ic
               wrf(d)=ic
               wcs(d)=ic
               wcf(d)=ic
               Wtcol=ExtractCPsubmatrix(Wt,wrs,wrf,wcs,wcf)

!               write(*,*) 'Wtsub, ic-th col, itn,d,ic = ',itn,ic,d
!               call Wtcol%print()

!              Normalization: calc norm tensor, mult by col, ALS
               call CPMatNormA(V,Wtcol,nals)
               call FlushCP(Wtcol)

!              Put normalized col ic of U into W as next rk terms
!              Y contains V, resized for W               
               Y=NewCP(W,rk) 
               Y%base(:,:)=0.d0
               ics(d)=ic
               call PutCPsubmatrix(V,irs,ics,Y)
               call FlushCP(V)
               call SUMVECVEC(W,1.d0,Y,1.d0)
               call FlushCP(Y)

!              Exit upon reaching last col since nothing to project out
               IF (ic.eq.nc) EXIT

!              Resize U --> V so projector only operates on cols > ic
               ics(d)=2
               icf(d)=U%cols(d)
               V=ExtractCPsubmatrix(U,irs,irf,ics,icf)
               call ReplaceVwithW(U,V)

!              Extract weights for cols remaining in V
               wrs(d)=ic+1
               wrf(d)=nc
               wcs(d)=ic+1
               wcf(d)=nc
               Wtcol=ExtractCPsubmatrix(Wt,wrs,wrf,wcs,wcf)

!               write(*,*) 'Wtsub, other cols, itn,d,ic = ',itn,ic,d
!               call Wtcol%print()

!              Project col ic out of other cols using ic-th block of W
!              U = [I - W*W^T]*U (1st step): Y = W^T*U
               rki=(ic-1)*rk+1
               rke=ic*rk
!               call CPMM(W,rki,rke,0,0.d0,.TRUE.,&
               call CPMM(W,1,rke,0,0.d0,.TRUE.,&
                         U,1,rk,0,0.d0,.FALSE.,Y)

!              Projection operator: diag(Y) for modes > d
               domode(:)=.TRUE.
               domode(1:d)=.FALSE.
               call CPMatrixZeroOffDiag(Y,domode)

!              Projection operator: reduce rank Z <-ALS-- Y
               Z=RandomCP(Y,rk)
               write(tag,'(A,3(I0,X))') &
               'OrthogCPmat(): Z<-Y, itn,d,ic = ',itn,d,ic
               call CPMatrixTranspose(Y)
               call CPMatrixTranspose(Z)
               redstat=ALS_reduce(Z,Y,Wtcol,nals,tag)
               call CPMatrixTranspose(Z)
!!! TEST
!              call CPMatrixTranspose(Y)
!              call CPEntrywiseCompare(Z,Y)
!              call ShowRank1Diff(Z,Y)
!!!

               call ReplaceVwithW(Y,Z)
!              Projection operator (2nd step): V = -W*Y + U
               call CPMM(W,rki,rke,0,0.d0,.FALSE.,&
                         Y,1,Y%R(),0,0.d0,.FALSE.,V)
               call FlushCP(Y)
               call SUMVECVEC(V,-1.d0,U,1.d0)

!              Reduce rank (column-weighted) of result, U <-ALS--V
               write(tag,'(A,3(I0,X))') &
               'OrthogCPmat(): U<-V, itn,d,ic = ',itn,d,ic
               call CPMatrixTranspose(V)
               call CPMatrixTranspose(U)
               redstat=ALS_reduce(U,V,Wtcol,nals,tag)
               call CPMatrixTranspose(U)
!!! TEST     
!              call CPMatrixTranspose(V)
!              call CPEntrywiseCompare(U,V)
!              call ShowRank1Diff(U,V)
!!!

               call FlushCP(Wtcol)
               call FlushCP(V)

               stu=stu+U%rows(d)
               fiu=fiu+U%rows(d)

!               call AbortWithError('Done with ic=1')
            ENDDO ! cols loop (ic)

!           Reduce rank (column-weighted) and replace : U <-ALS-- W.
            call FlushCP(U)

            U=RandomCP(rk,W%cols,W%rows,W%sym) ! Initalize as U^T
            call CPMatrixTranspose(W)
!            U=ALS_trials(rk,W,Wt,20,10*nals,'OrthogCPmat(): U<-W')
            redstat=ALS_reduce(U,W,Wt,10*nals,'OrthogCPmat(): U<-W')
            call CPMatrixTranspose(U) ! U <- U^T

!!! TEST     
!            call CPMatrixTranspose(W)
!            call CPEntrywiseCompare(U,W)
!            call ShowRank1Diff(U,W)
!!!
            call FlushCP(W)
         ENDDO ! modes loop (d)
      ENDDO ! main ALS loop (itn)

      deallocate(irs,irf,ics,icf,domode)
      deallocate(wrs,wrf,wcs,wcf)

      end subroutine OrthogCPmat_weights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module HGORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
