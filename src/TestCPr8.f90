!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE TESTCPR8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test the CPr8 module

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE CPr8

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine maintestcpr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Run some simple tests on CPr8 types

      implicit none
      TYPE (CP8) :: F,G,H
      integer, allocatable :: rows(:),cols(:)
      integer, allocatable :: rowi(:),coli(:)
      integer, allocatable :: rowf(:),colf(:)
      integer :: rk,ndof
      integer :: d,i,j,k,l,rtest
      real(kind=8)  :: val

      IF (mpirank.eq.0) THEN

      ndof=3
      rk=2
      allocate(rows(ndof),cols(ndof))
      allocate(rowi(ndof),coli(ndof))
      allocate(rowf(ndof),colf(ndof))
      rows(:)=3
      cols(:)=3
      rows(2)=4
      cols(3)=4

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'RandomRef_CP8 test:'
      write(*,*) '****************************************************'
      F=Random_CP(rk,rows,cols) 
      write(*,'(/X,A/)') 'ShowStats_CP8 test:'
      call F%show
      write(*,'(/X,A/)') 'PrintMat_CP8 test:'
      call F%printmat
      write(*,'(/X,A/)') 'PrintVec_CP8 test:'
      call F%printvec
      call Flush_CP8(F) !call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'IdentityMatrix_CP8 test'
      write(*,*) '****************************************************'
      F=IdentityMatrix(rows)
      call F%printmat
      call Flush_CP8(F) !call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'VectoDiagMatrix_CP8 test:'
      write(*,*) '****************************************************'
      rows(:)=2
      cols(:)=2
      rows(1)=3
      cols(2)=3
      F=Random_CP(rk,rows,cols)
      write(*,*) 'F:'
      call F%printmat
      G=F%diag()
      !G=VectoDiagMatrix(F)
      write(*,*) 'G:'
      call G%printmat
      call Flush_CP8(G)       !call G%flush
      call Flush_CP8(F)       !call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'GenCopyWtoV_CP8 test:'
      write(*,*) '****************************************************'
      call F%new(rk,rows,cols)
      call G%new(2*rk,rows,cols)
      F%coef(1)=1.0
      F%coef(2)=2.0
      F%base(:)=0.2d0
      call F%printmat
      call G%zero
      call GenCopyWtoV_CP8(G,F,2,3,1,2)
      call G%printmat
!      call F%flush
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ReplaceVwithW_CP8 test:'
      write(*,*) '****************************************************'
      call F%replace(G)
      call F%printmat
!      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'Resize_CP8 test:'
      write(*,*) '****************************************************'
      call F%resize(3)
      call F%printmat

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'TrimZeros_CP8 test:'
      write(*,*) '****************************************************'
      val=1.d-8
      call F%trim(val) !TrimZeros(F,val)
      write(*,*) 'F:'
      call F%printmat
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'MatrixTranspose_CP8 test:'
      write(*,*) '****************************************************'
      rows=(/4,3,2/)
      cols=(/1,1,1/)
      F=Random_CP(rk,rows,cols)
      call F%printmat
      call F%transpose() !MatrixTranspose_CP(F)
      call F%printmat
      call F%transpose() !MatrixTranspose_CP(F)
      call F%printmat
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'MatrixTranspose_CP8 test 2:'
      write(*,*) '****************************************************'
      rows=(/4,3,2/)
      cols=(/2,3,4/)
      F=Random_CP(rk,rows,cols)
      call F%printmat
      call F%transpose() !MatrixTranspose_CP(F)
      call F%printmat

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractVec_CP8:'
      write(*,*) '****************************************************'
      rowi(:)=1
      coli(:)=2
      write(*,*) 'G: second cols of F'
      G=ExtractVec(F,coli,.true.)
      call G%printmat
      call G%Flush
      write(*,*) 'G: first rows of F'
      G=ExtractVec(F,rowi,.false.)
      call G%printmat
      call G%Flush
      DEALLOCATE(rowi,coli)
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'MatrixZeroOffDiag_CP8 test:'
      write(*,*) '****************************************************'
      call MatrixZeroOffDiag(F)
      call F%printmat
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'MultOutCoefSmallest_CP8 test:'
      write(*,*) '****************************************************'
      rows=(/5,4,3/)
      cols=(/1,1,1/)
      call F%new(rk,rows,cols)
      F%coef(1)=2.0
      F%coef(2)=4.0
      F%base(:)=0.1d0
      write(*,*) 'F before:'
      call F%printvec
      call MultOutCoefSmallest_CP8(F)
      write(*,*) 'F after:'
      call F%printvec
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'DistributeCoef_CP8 test:'
      write(*,*) '****************************************************'
      rows=(/5,4,3/)
      cols=(/1,1,1/)
      call F%new(rk,rows,cols)
      F%coef(1)=2.0
      F%coef(2)=4.0
      F%base(:)=0.1d0
      write(*,*) 'F before:'
      call F%printvec
      call DistributeCoef(F)
      write(*,*) 'F after:'
      call F%printvec
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'GetRank1DominantEntry_gen_CP8:'
      write(*,*) '****************************************************'
      rows(:)=3
      cols(:)=3
      call F%new(rk,rows,cols)
      F%coef(1)=2.0
      F%coef(2)=4.0
      F%base(:)=0.1d0
      call F%put(2,2,2,1,0.5)  ! (2,2) in r=2, d=1
      call F%put(3,1,2,2,0.25) ! (3,1) in r=2, d=2
      call F%put(1,2,2,3,4.0)  ! (1,2) in r=2, d=3
      write(*,*) 'F:'
      call F%printmat
      rtest=2
      call GetRank1DominantEntry_gen_CP8(F,rtest,rowi,coli,val)
      do d=1,ndof
         write(*,'(A,I0,A,I0,A)') '(',rowi(d),',',coli(d),')'
      enddo
      write(*,*) 'largest entry for r = ',rtest,' is ',val
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractDiagfromMatrix_CP8:'
      write(*,*) '****************************************************'
      rows(:)=4
      F=Random_CP(rk,rows,rows)
      write(*,*) 'F:'
      call F%printmat
      write(*,*) 'G: F diagonal as column vec:'
      G=ExtractDiagfromMatrix(F,.false.)
      call G%printmat
      call G%flush
      write(*,*) 'G: F diagonal as row vec:'
      G=ExtractDiagfromMatrix(F,.true.)
      call G%printmat
      call G%flush
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractSubmatrix_CP8:'
      write(*,*) '****************************************************'
      rows=(/5,4,3/)
      cols=(/3,4,5/)
      F=Random_CP(rk,rows,cols)
      write(*,*) 'F:'
      call F%printmat
      F%coef=(/1.d0,2.d0/)
      rowi=(/4,2,1/)
      rowf=(/5,3,3/)
      coli=(/2,2,2/)
      colf=(/3,3,2/)
      G=F%submatrix(rowi,rowf,coli,colf)
      write(*,*) 'G:'
      call G%printmat
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'PutSubmatrix_CP8:'
      write(*,*) '****************************************************'
      G%base(:)=0.1d0
      call PutSubmatrix(G,rowi,coli,F)
      write(*,*) 'F: part extracted as G overwritten'
      call F%printmat
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractMatrixElement_CP8:'
      write(*,*) '****************************************************'
      rowi=(/2,2,3/)
      coli=(/1,1,4/)
      do d=1,ndof
         write(*,'(A,I0,A,I0,A)') '(',rowi(d),',',coli(d),')'
      enddo      
      val=ExtractMatrixElement(F,rowi,coli)
      write(*,*) 'Full entry is ',val
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ModeJoin_CP8:'
      write(*,*) '****************************************************'
      H=ModeJoin(F,G) 
      write(*,*) 'H: F joined to G'
      call H%printmat
      call H%flush
      call G%flush
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'Sum_CP8:'
      write(*,*) '****************************************************'
      rows(:)=4
      cols(:)=1
      write(*,*) 'F:'
      F=Random_CP(rk,rows,cols)
      call F%printvec
      write(*,*) 'G:'
      call G%new0(rows,cols) !G=Zero_CP(rows,cols)
      call G%printvec
      call F%sumcp(2.0,G,10.0)
      write(*,*) 'F = 2.0*F + 10.0*G:'
      call F%printvec
      write(*,*) 'G = 3.0*G - 0.25*F:'
      call G%sumcp(3.0,F,-0.25)
      call G%printvec
      write(*,*) 'F=-F+G (with new G)'
      call F%sumcp(-1.0,G,1.0)
      call F%printvec


      call G%flush
      call F%flush

      deallocate(rows,cols,rowi,rowf,coli,colf)


      ENDIF ! MPI
      call AbortWithError("Done testing CPr8")

      end subroutine maintestcpr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TESTCPR8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
