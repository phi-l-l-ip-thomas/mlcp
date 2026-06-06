!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE TESTCPR8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Test the CPr8 module

      USE ERRORTRAP
      USE UTILS
      USE MYMPI
      USE SEPDREPN
      USE CPr8
      USE MODVECVEC
      USE ALSOO
      USE ALSOO8
      USE ALSDRVR
      USE ALS8DRVR
      USE CPMMM
      USE LINSOLVER
      USE CPMM8
      USE CPVV8
      USE CPLS8
      USE LinAlg8
      USE ALSDRVR
#if ACC_ENABLED
      USE MyACC
#endif

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine maintestcpr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Run some simple tests on CPr8 types

      implicit none
      TYPE (CP8) :: F,G,H
      TYPE (CP)  :: Z
      integer, allocatable :: rows(:),cols(:)
      integer, allocatable :: rowi(:),coli(:)
      integer, allocatable :: rowf(:),colf(:)
      logical, allocatable :: modes(:)
      integer :: rk,ndof,ish1,ish2
      integer :: d,i,j,k,l,rtest,s1,s2,lim1,lim2,sz1,sz2
      integer :: it1,it2,szf1,szf2,szg1,szg2
      logical :: trans1,trans2,passtest
      real(kind=8)  :: val,Esh1,Esh2,stime,ftime

      IF (mpirank.eq.0) THEN

      ndof=3
      rk=2
      allocate(modes(ndof))
      allocate(rows(ndof),cols(ndof))
      allocate(rowi(ndof),coli(ndof))
      allocate(rowf(ndof),colf(ndof))
      modes(:)=.true.
      rows(:)=3
      cols(:)=3
      rows(2)=4
      cols(3)=4

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'RandomRef_CP8 test:'
      write(*,*) '****************************************************'
      call F%newrand(rk,rows,cols) 
      write(*,'(/X,A/)') 'ShowStats_CP8 test:'
      call F%show
      write(*,'(/X,A/)') 'PrintMat_CP8 test:'
      call F%printmat
      write(*,'(/X,A/)') 'PrintVec_CP8 test:'
      call F%printvec
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'IdentityMatrix_CP8 test'
      write(*,*) '****************************************************'
      call F%identity(rows,rows)
      call F%printmat
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'VectoDiagMatrix_CP8 test:'
      write(*,*) '****************************************************'
      rows(:)=2
      cols(:)=2
      rows(1)=3
      cols(2)=3
      call F%newrand(rk,rows,cols)
      write(*,*) 'F:'
      call F%printmat
      G=F%diag()
      !G=VectoDiagMatrix(F)
      write(*,*) 'G:'
      call G%printmat
      call G%flush
      call F%flush

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
      call G%copy_terms(F,2,3,1,2)
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
      call F%newrand(rk,rows,cols)
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
      call F%newrand(rk,rows,cols)
      call F%printmat
      call F%transpose() !MatrixTranspose_CP(F)
      call F%printmat

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractVec_CP8:'
      write(*,*) '****************************************************'
      rowi(:)=1
      coli(:)=2
      write(*,*) 'G: second cols of F'
      call G%ExtractVec(F,coli,.true.)
      call G%printmat
      call G%Flush
      write(*,*) 'G: first rows of F'
      call G%ExtractVec(F,rowi,.false.)
      call G%printmat
      call G%Flush
      DEALLOCATE(rowi,coli)
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'MatrixZeroOffDiag_CP8 test:'
      write(*,*) '****************************************************'
!      call MatrixZeroOffDiag(F)
      call F%zerooffdiagonalall
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
      call F%multoutcoef
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
      call F%distributecoef
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
      call F%largestentryinterm(rtest,rowi,coli,val)
      do d=1,ndof
         write(*,'(A,I0,A,I0,A)') '(',rowi(d),',',coli(d),')'
      enddo
      write(*,*) 'largest entry for r = ',rtest,' is ',val
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractDiagfromMatrix_CP8:'
      write(*,*) '****************************************************'
      rows(:)=4
      call F%newrand(rk,rows,rows)
      write(*,*) 'F:'
      call F%printmat
      write(*,*) 'G: F diagonal as column vec:'
      call G%extractdiagonal(F,.false.)
      call G%printmat
      call G%flush
      write(*,*) 'G: F diagonal as row vec:'
      call G%extractdiagonal(F,.true.)
!      G=ExtractDiagfromMatrix(F,.true.)
      call G%printmat
      call G%flush
      call F%flush

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractSubmatrix_CP8:'
      write(*,*) '****************************************************'
      rows=(/5,4,3/)
      cols=(/3,4,5/)
      call F%newrand(rk,rows,cols)
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
      write(*,*) 'G: will overwrite part of F'
      call G%printmat
!      call F%copyto(H)
      call F%putsubmatrix(G,rowi,coli)
      write(*,*) 'F: part extracted as G overwritten'
      call F%printmat
!      write(*,*) '****************************************************'
!      write(*,'(/X,A/)') 'EntrywiseCompare_CP8:'
!      write(*,*) '****************************************************'
!      call F%compareentries(H)
!      call H%flush
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ExtractMatrixElement_CP8:'
      write(*,*) '****************************************************'
      rowi=(/2,2,3/)
      coli=(/1,1,4/)
      do d=1,ndof
         write(*,'(A,I0,A,I0,A)') '(',rowi(d),',',coli(d),')'
      enddo      
      val=F%tensorelement(rowi,coli)
      write(*,*) 'Full entry is ',val
      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ModeJoin_CP8:'
      write(*,*) '****************************************************'
      call H%modejoin(F,G) 
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
      call F%newrand(rk,rows,cols)
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

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'CP8toCP() and CPtoCP8():'
      write(*,*) '****************************************************'

      write(*,*) 'G:'
      call G%printvec

      call G%toCP(Z)
      write(*,*) 'G (CP8) -> Z (CP); Z:'
      call PrintCPvec(Z)
      passtest=compareCPandCP8(Z,G)

      call H%fromCP(Z)
      write(*,*) 'Z (CP) -> H (CP8); H:'
      call H%printvec

      call G%flush
      call F%flush
      call H%Flush

!      call testcpmm()
!      call testpvv()
!       call testpvvperf()
!      call testls
!      call testaux
!      call testals
      call testals8all()

      call AbortWithError("Done testing CPr8")

      ENDIF ! MPI rank 0

      end subroutine maintestcpr8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testcpmm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Run some simple tests on CPr8 types

      implicit none
      TYPE (CP8) :: F,G,H
      TYPE (CP)  :: ZF,ZG,ZH
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: modes(:)
      integer :: rk,ndof,ish1,ish2,algo
      integer :: s1,s2,lim1,lim2,sz1,sz2
      integer :: it1,it2,szf1,szf2,szg1,szg2
      logical :: trans1,trans2,passtest,t3
      real(kind=8)  :: val,Esh1,Esh2,stime,ftime

      algo=1
      t3=.TRUE.
      ndof=3
      rk=2
      allocate(modes(ndof))
      allocate(rows(ndof),cols(ndof))
      modes(:)=.true.

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'CPMM8_CP8 test:'
      write(*,*) '****************************************************'
!      rows(:)=4
!      cols(:)=3
!      call F%newrand(rk,rows,cols)
!      F%coef(1:2)=(/3.0,5.0/)
!      write(*,*) 'F:'
!      call F%printmat

!      rows(:)=3
!      cols(:)=4
!      call G%newrand(rk,rows,cols)
!      G%coef(1:2)=(/2.25,1.25/)
!      write(*,*) 'G:'
!      call G%printmat

      do it1=1,2
         trans1=(it1.eq.1)
         do it2=1,2
            trans2=(it2.eq.1)
            szf1=4
            szf2=3
            szg1=3
            szg2=4
            if (trans1) then
               if (.not.trans2) call swap(szg1,szg2)
            endif
            if (trans2) then
               if (.not.trans1) call swap(szf1,szf2)
            endif
            rows(:)=szf1
            cols(:)=szf2
            call F%newrand(rk+1,rows,cols)
            F%coef(1:3)=(/3.0,5.0,7.7/)
            rows(:)=szg1
            cols(:)=szg2
            call G%newrand(rk,rows,cols)
            G%coef(1:2)=(/2.25,1.25/)

            write(*,*)
            write(*,*) 'F transposed: ',trans1,'; G transposed: ',trans2
            write(*,*)

!           Test all possible combinations of sign and ishift values
            write(*,*) 'Checking all ish and sign-of-Esh combinations'
            do ish2=-1,1,1
               lim2=-1
               if (ish2.eq.0) lim2=1
               do ish1=-1,1,1
                  lim1=-1
                  if (ish1.eq.0) lim1=1
                  do s1=lim1,1,2
                     do s2=lim2,1,2
                        Esh1=s1*4.4
                        Esh2=s2*2.9

                        write(*,*) ish1,ish2,Esh1,Esh2 
                        if (algo.eq.1) then
                           call F%copyintodevice()
                           call G%copyintodevice()
                        endif
                        call CPMM_CP8(F,ish1,Esh1,trans1,&
                                      G,ish2,Esh2,trans2,H,t3,algo)
                        if (algo.eq.1) then
                           call F%deletefromdevice()
                           call G%deletefromdevice()
                           call H%updatefromdevice()
                        endif
!                  write(*,*) 'H (result of matrix mult):'
!                   call H%printmat

                        write(*,*) 'Compare with existing vsn...'
                        call F%toCP(ZF)
                        call G%toCP(ZG)
                        call CPMM(ZF,ish1,Esh1,trans1,ZG,ish2,Esh2,trans2,ZH)
                        passtest=compareCPandCP8(ZH,H)

                        call H%flush()
                        call FlushCP(ZF)
                        call FlushCP(ZG)
                        call FlushCP(ZH)
                     enddo
                  enddo
               enddo
            enddo
            call F%flush()
            call G%flush()
         enddo
      enddo

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'CPMM8_CP8 perf test:'
      write(*,*) '****************************************************'
      rows(:)=256
      cols(:)=384
      call F%newrand(rk,rows,cols)
      F%coef(1:3)=(/3.0,5.0,7.7/)
!      write(*,*) 'F:'
!      call F%printmat

      rows(:)=384
      cols(:)=256
      call G%newrand(rk,rows,cols)
      G%coef(1:2)=(/2.25,1.25/)
!      write(*,*) 'G:'
!      call G%printmat

      trans1=.FALSE.
      trans2=.FALSE.
      ish1=1
      ish2=1
      Esh1=4.4
      Esh2=2.9

      write(*,*) 'Reference calc on host...'
      call F%toCP(ZF)
      call G%toCP(ZG)
      call CPMM(ZF,ish1,Esh1,trans1,ZG,ish2,Esh2,trans2,ZH)

      write(*,*) 'CPMM8, new host version...'
      call CPU_TIME(stime)
      call CPMM_CP8(F,ish1,Esh1,trans1,&
                    G,ish2,Esh2,trans2,H,modes,t3,0)
      call CPU_TIME(ftime)
      write(*,*) ' -> Time for host CPMM: ',ftime-stime
      passtest=compareCPandCP8(ZH,H)
      call H%flush()

!      write(*,*) 'CPMM8, CUBLAS version...' ! (Removed)
!      call CPU_TIME(stime)
!      call CPMM_CP8(F,ish1,Esh1,trans1,&
!                    G,ish2,Esh2,trans2,H,modes,t3,1)
!      call CPU_TIME(ftime)
!      write(*,*) ' -> Time for host CPMM: ',ftime-stime
!      passtest=compareCPandCP8(ZH,H)
!      call H%flush()

      write(*,*) 'CPMM8, CUTENSOR version...'
      call CPU_TIME(stime)
      if (algo.eq.1) then
         call F%copyintodevice()
         call G%copyintodevice()
      endif
      call CPMM_CP8(F,ish1,Esh1,trans1,&
                    G,ish2,Esh2,trans2,H,modes,t3,1)
      if (algo.eq.1) then
         call F%deletefromdevice()
         call G%deletefromdevice()
         call H%updatefromdevice()
      endif

      call CPU_TIME(ftime)
      write(*,*) ' -> Time for host CPMM: ',ftime-stime
      passtest=compareCPandCP8(ZH,H)
      call H%flush()

      call FlushCP(ZF)
      call FlushCP(ZG)
      call FlushCP(ZH)
      call F%flush()
      call G%flush()

      deallocate(rows,cols)

      end subroutine testcpmm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      logical function compareCPandCP8(F,F8) result(pass)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares all entries in F and F8

      implicit none
      TYPE (CP8), intent(in) :: F8
      TYPE (CP), intent(in) :: F
      integer :: d,ndof,r,rk,i,m,k,n,l
      real(kind=8), parameter :: tol=1.d-13
      real(kind=8) :: fval,f8val
      logical :: same
      
      pass=.true.

      if (F8%R().ne.F%R()) then
         write(*,*) 'Ranks not equal: F8:',F8%R(),'; F:',F%R()
         pass=.false.
      elseif (F8%D().ne.F%D()) then
         write(*,*) 'ndof not equal: F8:',F8%D(),'; F:',F%D()
         pass=.false.
      else
         ndof=F%D()
         same=.true.
         do d=1,ndof
            if (F8%M(d).ne.F%M(d) .or. F8%N(d).ne.F%N(d)) then
               same=.false.
               pass=.false.
               write(*,*) 'mode',d,': F8(',F8%M(d),' x ',F8%N(d),&
                          ') != F(',F%M(d),' x ',F%N(d),')'
            endif
         enddo
         if (same) then
            rk=F%R()
            write(*,*) 'Comparing coefs...'
            do r=1,rk       
               if (abs(F8%coef(r)-F%coef(r)).gt.tol) then
                  write(*,*) r,F8%coef(r),F%coef(r),F8%coef(r)-F%coef(r)
                  pass=.false.
               endif
            enddo
            write(*,*)
            write(*,*) 'Comparing base...'
            do d=1,ndof
               m=F%M(d)
               n=F%N(d)
               do r=1,rk
                  l=1
                  do k=1,n
                     do i=1,m
                        fval=F%base(F%ibas(d)-1+l,r)
                        f8val=F8%get(i,k,r,d)
                        if (abs(f8val-fval).gt.tol) then
                           write(*,'(4(I5,X),3(ES11.4,x))') & 
                           r,d,i,k,f8val,fval,f8val-fval
                           pass=.false.
                        endif
                        l=l+1
                     enddo
                  enddo
               enddo
            enddo
            write(*,*) 'Compare done!'
         endif
      endif

      end function compareCPandCP8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testpvv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Run some simple tests on CPr8 types

      implicit none
      TYPE (CP8) :: F,G
      TYPE (CP)  :: ZF,ZG
      integer, allocatable :: rows(:),cols(:)
      real(kind=8), allocatable :: Ptot(:),Pk(:),Pnk(:),Pfull(:)
      real(kind=8), allocatable :: Qtot(:),Qk(:),Qnk(:),Qfull(:)
      real(kind=8), allocatable :: Rtot(:,:),Rk(:,:),Rnk(:,:),Rfull(:,:)
      logical, allocatable :: modes(:)
      integer :: rkf,rkg,ndof,i,j
      character*64 :: frmt
      logical :: passtest
      real(kind=8) :: stime,ftime,pvvP,pvvQ,pvvR

      ndof=3
      allocate(rows(ndof),cols(ndof),modes(ndof))
      rows(:)=3
      cols(:)=3
      rows(2)=5
      cols(2)=1
      rows(3)=1
      cols(3)=10
      rkf=5
      rkg=4

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'CPVV8_CP8 tests:'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
!      F%coef(1:5)=(/3.0,5.0,7.7,0.3,2.9/)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)
      write(*,*) 'F:'
!      call F%printvec

      call G%newrand(rkg,rows,cols)
!      G%coef(1:4)=(/2.25,1.25,3.55,0.11/)
      call random_number(G%coef)
      G%coef(:)=12*G%coef(:)
      write(*,*) 'G:'
!      call G%printvec

      allocate(Rtot(F%R(),G%R()),Rk(F%R(),G%R()))
      allocate(Rnk(F%R(),G%R()),Rfull(F%R(),G%R()))
      allocate(Qtot(F%R()*G%R()),Qk(F%R()*G%R()))
      allocate(Qnk(F%R()*G%R()),Qfull(F%R()*G%R()))
      allocate(Ptot(F%R()*G%R()),Pk(F%R()*G%R()))
      allocate(Pnk(F%R()*G%R()),Pfull(F%R()*G%R()))

      Rtot=0.d0
      Rk=0.d0
      Rnk=0.d0
      Rfull=0.d0
      Qtot=0.d0
      Qk=0.d0
      Qnk=0.d0
      Qfull=0.d0
      Ptot=0.d0
      Pk=0.d0
      Pnk=0.d0
      Pfull=0.d0


      write(*,*) 'Reference calcs on host...'
      call F%toCP(ZF)
      call G%toCP(ZG)

      call CONSTPT(ZG,ZF,0,Rtot)
      call CONSTPT_CP8(F,G,0,Qtot,0)
      call checkRPtest(Rtot,Qtot,'Qtot (new host)')
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Ptot)
      call CONSTPT_CP8(F,G,0,Ptot,1)
      !$acc exit data copyout(Ptot)
      call F%deletefromdevice()
      call G%deletefromdevice()
      call checkRPtest(Rtot,Ptot,'Ptot (new cutensor)')
#endif

      call UPDATEP(ZG,ZF,2,Rtot,.TRUE.)
      call UPDATEP_CP8(F,G,2,Qtot,.TRUE.,0)
      call checkRPtest(Rtot,Qtot,'Qtot downdate (new host)')
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Ptot)
      call UPDATEP_CP8(F,G,2,Ptot,.TRUE.,1)
      !$acc exit data copyout(Ptot)
      call F%deletefromdevice()
      call G%deletefromdevice()
      call checkRPtest(Rtot,Ptot,'Ptot downdate (new cutensor)')
#endif

      call UPDATEP(ZG,ZF,2,Rtot,.FALSE.)
      call UPDATEP_CP8(F,G,2,Qtot,.FALSE.,0)
      call checkRPtest(Rtot,Qtot,'Qtot update (new host)')
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Ptot)
      call UPDATEP_CP8(F,G,2,Ptot,.FALSE.,1)
      !$acc exit data copyout(Ptot)
      call F%deletefromdevice()
      call G%deletefromdevice()
      call checkRPtest(Rtot,Ptot,'Ptot update (new cutensor)')
#endif

      call CONSTPT(ZG,ZF,1,Rnk)
      call CONSTPT_CP8(F,G,1,Qnk,0)
      call checkRPtest(Rnk,Qnk,'Qnk (new host)')
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Pnk)
      call CONSTPT_CP8(F,G,1,Pnk,1)
      !$acc exit data copyout(Pnk)
      call F%deletefromdevice()
      call G%deletefromdevice()
      call checkRPtest(Rnk,Pnk,'Pnk (new cutensor)')
#endif

      call CONSTPk(ZG,ZF,3,Rk)
      call CONSTPk_CP8(F,G,3,Qk,0)
      call checkRPtest(Rk,Qk,'Qk (new host)')
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Pk)
      call CONSTPk_CP8(F,G,3,Pk,1) 
      !$acc exit data copyout(Pk)
      call F%deletefromdevice()
      call G%deletefromdevice()
      call checkRPtest(Rk,Pk,'Pk (new cutensor)')
#endif

      pvvR=PRODVV(ZF,ZG)
      call UPDATEPCoef_CP8(F,G,Qtot,.FALSE.,0)
      pvvQ=SUM(Qtot)
#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      !$acc enter data copyin(Ptot)
      call UPDATEPCoef_CP8(F,G,Ptot,.FALSE.,1)
      !$acc exit data copyout(Ptot)
      call F%deletefromdevice()
      call G%deletefromdevice()
      pvvP=SUM(Ptot)
#endif
      write(*,*) 'ProdVV(F,G) host (original algo)   = ',pvvR
      write(*,*) 'ProdVV(F,G) from new host PTOT     = ',pvvQ

#if ACC_ENABLED
      write(*,*) 'ProdVV(F,G) from new cutensor PTOT = ',pvvP
#endif

      call FlushCP(ZF)
      call FlushCP(ZG)
      call F%flush()
      call G%flush()

      deallocate(rows,cols)
      deallocate(Ptot,Pk,Pnk,Pfull)
      deallocate(Rtot,Rk,Rnk,Rfull)

      end subroutine testpvv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testpvvperf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Run some simple tests on CPr8 types

      implicit none
      TYPE (CP8) :: F,G
      TYPE (CP)  :: ZF,ZG
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: modes(:)
      integer :: rkf,rkg,ndof
      real(kind=8) :: ti1,ti2,pvv

      ndof=6
      allocate(rows(ndof),cols(ndof),modes(ndof))
      rows(:)=256
      cols(:)=1
      rkf=100
      rkg=30000

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'PVV performance tests:'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)

      call G%newrand(rkg,rows,cols)
      call random_number(G%coef)
      G%coef(:)=12*G%coef(:)

      call F%toCP(ZF)
      call G%toCP(ZG)

      call CPU_TIME(ti1)
      pvv=PRODVV(ZF,ZG)
      call CPU_TIME(ti2)
      write(*,*) 'Original alg time: ',ti2-ti1

      call CPU_TIME(ti1)
      pvv=PRODVV_CP8(F,G,0)
      call CPU_TIME(ti2)
      write(*,*) 'New cpu alg time: ',ti2-ti1

#if ACC_ENABLED
      call F%copyintodevice()
      call G%copyintodevice()
      call CPU_TIME(ti1)
      pvv=PRODVV_CP8(F,G,1)
      call CPU_TIME(ti2)
      write(*,*) 'New gpu alg time: ',ti2-ti1
#else
         call AbortWithError('testals(): OpenACC alg not implemented')
#endif

      call F%flush()
      call G%flush()
      call FlushCP(ZF)
      call FlushCP(ZG)

      end subroutine testpvvperf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testls

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests rhs module

      implicit none
      TYPE (LS8) :: hU,tU
      TYPE (ALS) :: X
      TYPE (CP8) :: H,F,F2,G,G2,HF,HG,HWG,WHF,V,W,WH,WF,WG,HTWH,Wd,Wd2,WG2
      TYPE (CP8) :: dummyH,dummyW
      TYPE (CP)  :: ZH,ZF,ZG,ZHF,ZHG,ZV,ZW,ZWFG,ZWG,ZWF,ZWHF,ZWd
      integer, allocatable :: rows(:),cols(:)
      real(kind=8), allocatable :: tRHS(:),tLHS(:)
      real(kind=8), allocatable :: hPnk(:),hBnk(:),hRHS(:),hLHS(:)
      real(kind=8), allocatable :: aPnk(:),aBnk(:),aRHS(:),aLHS(:)
      real(kind=8), allocatable :: oRHS(:,:),oLHS(:,:)
      logical, allocatable :: modes(:)
      logical :: passtest,t3
      integer :: rkf,rkg,rkh,rkhs,rkw,ndof,i,j,ish,tstmode
      real(kind=8) :: stime,ftime,Esh,reg

      t3=.TRUE.
      ndof=3
      tstmode=3
      allocate(rows(ndof),cols(ndof),modes(ndof))
      modes(:)=.TRUE.
      rows(:)=5
      cols(:)=1
      rows(2)=7
      cols(2)=1
      rows(3)=4
      cols(3)=2
      rkf=4
      rkg=8
      rkh=3
      rkw=2
      ish=1
      Esh=11.3
      rkhs=rkh
      if (ish.ne.0) rkhs=rkhs+1
      reg=0.0000001

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'CPLS8_CP8 LinSolv and ALS tests:'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)
      write(*,*) 'F:'
      call F%printvec
      call F%toCP(ZF)

      call G%newrand(rkg,rows,cols)
      call random_number(G%coef)
      G%coef(:)=12*G%coef(:)
      write(*,*) 'G:'
      call G%printvec
      call G%toCP(ZG)
      call G2%copyfrom(G)
      call PrepG_CP8(G2,0)
      write(*,*) 'G2:'
      call G2%printvec

      call H%newrand(rkh,rows,rows) ! Square operator matrix
      call random_number(H%coef)
      H%coef(:)=5*H%coef(:)
      write(*,*) 'H:'
      call H%printvec
      call H%toCP(ZH)

      call W%newrand(rkw,rows,rows) ! Square weights matrix
      W%base=abs(W%base)
      call random_number(W%coef)
      W%coef(:)=5*W%coef(:)
      write(*,*) 'W:'
      call W%printvec
      call W%toCP(ZW)

      call Wd%copyfrom(W) 
      call Wd%zerooffdiagonalall()
      write(*,*) 'diagonal W (for wALS test):'
      call Wd%printmat
      call Wd%toCP(ZWd)
      call Wd2%copyfrom(Wd)
      call PrepG_CP8(Wd2,0)
      write(*,*) 'Wd2:'
      call Wd2%printvec

      write(*,*)
      write(*,*) ' --- ALS test (no weights): A->I, F, G ---'
      write(*,*)

      allocate(hBnk(F%R()*F%R()),hPnk(F%R()*G%R()))

      write(*,*) ' 1) Build B and P matrices (only on host)...'
      call CONSTPT_CP8(F,F,tstmode,hBnk,0)
      call CONSTPT_CP8(G,F,tstmode,hPnk,0)

      write(*,*) ' 2) Create U,X types'
      call hU%new(dummyH,F,G,dummyW,0,0.0,0)
#if ACC_ENABLED
      call tU%new(dummyH,F,G,dummyW,0,0.0,1)
#endif
      call NewALS(X,ZF,ZG)
      call ProdMatsDD(X,ZF,ZG,tstmode)

      write(*,*) ' 3) Build LHS,mode: ',tstmode
      oLHS=GetLHSforALS(X,tstmode)
      call GetLHS_CP8(hU,hBnk,hLHS) 
      call checkRPtest(oLHS,hLHS,'hLHS (new host)')
#if ACC_ENABLED
      !$acc enter data copyin(hBnk)
      call GetLHS_CP8(tU,hBnk,tLHS)
      !$acc exit data delete(hBnk) copyout(tLHS)
      call checkRPtest(oLHS,tLHS,'tLHS (new cutensor)')
#endif

      write(*,*) ' 4) Build RHS, mode: ',tstmode
      oRHS=GetRHSforALS(X,ZG,tstmode)
      call GetRHS_CP8(hU,hPnk,G2,hRHS,tstmode)     
      call checkRPtest(oRHS,hRHS,'hRHS (new host)')
#if ACC_ENABLED
      !$acc enter data copyin(hPnk)
      call G2%copyintodevice()
      call GetRHS_CP8(tU,hPnk,G2,tRHS,tstmode)
      call G2%deletefromdevice()
      !$acc exit data delete(hPnk) copyout(tRHS)
      call checkRPtest(oRHS,tRHS,'tRHS (new cutensor)')
#endif

      write(*,*) ' 5) Clean up'
      call FlushALS(X)
      call hU%flush()
#if ACC_ENABLED
      call tU%flush()
      deallocate(tLHS,tRHS)
#endif
      deallocate(oLHS,oRHS,hLHS,hRHS,hBnk,hPnk)

      write(*,*)
      write(*,*) ' --- ALS test (weights): A->I, F, G ---'
      write(*,*)

      allocate(hBnk(rkW*rkF*rkF))
      allocate(hPnk(rkW*rkG*rkF))

      write(*,*) '1) CPMM ->W*F, ->W*G'
      ! Z versions use the old host algorithm
      call NewALS(X,ZF,ZG,ZWd)
      call CPMM(ZWd,.FALSE.,ZF,.FALSE.,ZWF)
      call CPMM(ZWd,.FALSE.,ZG,.FALSE.,ZWG)

      call hU%new(dummyH,F,G,Wd,ish,Esh,0)
#if ACC_ENABLED
      call tU%new(dummyH,F,G,Wd,ish,Esh,1)
#endif

!     Create WF, WG, and rearranged-WG
      call CPMM_CP8(Wd,.FALSE.,F,.FALSE.,WF,.FALSE.,0)
      call CPMM_CP8(Wd,.FALSE.,G,.FALSE.,WG,t3,0)
      call WG2%copyfrom(WG)
      call PrepG_CP8(WG2,0)

      write(*,*) '2) Calc B matrix'
      call ProdMatsDD(X,ZF,ZWF,ZG,ZWG,tstmode)
      call CONSTPT_CP8(WF,F,tstmode,hBnk,0)
      call CONSTPT_CP8(WG,F,tstmode,hPnk,0)

      write(*,*) '3) Build LHS,mode: ',tstmode
      oLHS=GetLHSforALS(X,ZWd,tstmode)
      call GetLHS_CP8(hU,hBnk,Wd2,hLHS,tstmode)
      call checkRPtest(oLHS,hLHS,'hLHS (new host)')
      write(*,*) ' ...hLHS size is ',SIZE(hLHS)
#if ACC_ENABLED
      !$acc enter data copyin(hBnk)
      call Wd2%copyintodevice()
      call GetLHS_CP8(tU,hBnk,Wd2,tLHS,tstmode)
      call Wd2%deletefromdevice()
      !$acc exit data delete(hBnk) copyout(tLHS)
      call checkRPtest(oLHS,tLHS,'tLHS (new cutensor)')
#endif

      write(*,*) '4) Build RHS, mode: ',tstmode
      oRHS=GetRHSforALS(X,ZWG,tstmode)
      call GetRHS_CP8(hU,hPnk,WG2,hRHS,tstmode)
      call checkRPtest(oRHS,hRHS,'hRHS (new host)')
      write(*,*) ' ...hRHS size is ',SIZE(hRHS)
#if ACC_ENABLED
      !$acc enter data copyin(hPnk)
      call WG2%copyintodevice()
      call GetRHS_CP8(tU,hPnk,WG2,tRHS,tstmode)
      call WG2%deletefromdevice()
      !$acc exit data delete(hPnk) copyout(tRHS)
      call checkRPtest(oRHS,tRHS,'tRHS (new cutensor)')
#endif

      call X%flush()
      call WF%flush()
      call WG%flush()
      call WG2%flush()
      call hU%flush()
#if ACC_ENABLED
      call tU%flush()
      deallocate(tLHS,tRHS)
#endif
      call FlushCP(ZWF)
      call FlushCP(ZWG)

      deallocate(hBnk,hPnk,oLHS,hlhs,orhs,hrhs)

      write(*,*)
      write(*,*) ' --- LinSolv test: A->(H-E*I), F, (H-E*I)^TG ---'
      write(*,*)

      allocate(hBnk(rkW*(rkF*rkHs)**2))
      allocate(hPnk(rkHs*rkF*rkW*rkG))

      write(*,*) '1) CPMM ->H*F, ->H^T*G'
      ! Z versions use the old host algorithm
      call NewLinSolver(X,ZH,ZF,ZG,ZW,.FALSE.,.FALSE.,.FALSE.,ish,Esh)
      ! The call above makes ZHTWH and ZHWG
      call CPMM(ZH,ish,Esh,.FALSE.,ZF,0,0.0,.FALSE.,ZHF)
      call CPMM(ZW,.FALSE.,ZHF,.FALSE.,ZWHF)
      call CPMM(ZW,.FALSE.,ZG,.FALSE.,ZWG)

      call hU%new(H,F,G,W,ish,Esh,0)
#if ACC_ENABLED
      call tU%new(H,F,G,W,ish,Esh,1)
#endif

!     Create HTWH
      call CPMM_CP8(W,0,0.0,.FALSE.,H,ish,Esh,.FALSE.,WH,modes,.FALSE.,0)
      call CPMM_CP8(H,ish,Esh,.TRUE.,WH,0,0.0,.FALSE.,HTWH,modes,.TRUE.,0)
      call PrepG_CP8(HTWH,0)

!     Create WHF
      call CPMM_CP8(H,ish,Esh,.FALSE.,F,0,0.0,.FALSE.,HF,modes,.FALSE.,0)
      call CPMM_CP8(W,.FALSE.,HF,.FALSE.,WHF,.FALSE.,0)

!     Create HTWG
      call CPMM_CP8(W,.FALSE.,G,.FALSE.,WG,t3,0)
      call CPMM_CP8(H,ish,Esh,.TRUE.,WG,0,0.0,.FALSE.,HWG,modes,t3,0)
      call PrepG_CP8(HWG,0)

      write(*,*) '2) Calc B matrix'
      call RecomputeNormalEquations(X,ZH,ish,Esh,ZWG,ZW,tstmode)
      call ProdMatsDD(X,ZHF,ZWHF,ZG,ZWG,tstmode)
      call CONSTPT_CP8(WHF,HF,tstmode,hBnk,0)
      call CONSTPT_CP8(WG,HF,tstmode,hPnk,0)

      write(*,*) '3) Build LHS, mode: ',tstmode
      oLHS=GetLHSforLinSys(X,ZF,tstmode)
      call GetLHS_CP8(hU,hBnk,HTWH,hLHS,tstmode)
      call checkRPtest(oLHS,hLHS,'hLHS (new host)')
#if ACC_ENABLED
      !$acc enter data copyin(hBnk)
      call HTWH%copyintodevice()
      call GetLHS_CP8(tU,hBnk,HTWH,tLHS,tstmode)
      call HTWH%deletefromdevice()
      !$acc exit data delete(hBnk) copyout(tLHS)
      call checkRPtest(oLHS,tLHS,'tLHS (new cutensor)')
#endif

      write(*,*) '4) Build RHS, mode: ',tstmode
      oRHS=GetRHSforLinSys(X,tstmode)
      call GetRHS_CP8(hU,hPnk,HWG,hRHS,tstmode)
      call checkRPtest(oRHS,hRHS,'hRHS (new host)')
#if ACC_ENABLED
      !$acc enter data copyin(hPnk)
      call HWG%copyintodevice()
      call GetRHS_CP8(tU,hPnk,HWG,tRHS,tstmode)
      call HWG%deletefromdevice()
      !$acc exit data delete(hPnk) copyout(tRHS)
      call checkRPtest(oRHS,tRHS,'tRHS (new cutensor)')
#endif

      write(*,*)
      write(*,*) 'Now to call the solver with LHS,RHS'
      call SolveLinSysLU(oLHS,oRHS,reg)
      call SolveLinSysLU8(hLHS,hRHS,SIZE(oRHS,1),SIZE(oRHS,2),reg,0)
      call checkRPtest(oRHS,hRHS,'hRHS (new host)')
#if ACC_ENABLED
      !$acc enter data copyin(tLHS,tRHS)
      call SolveLinSysLU8(tLHS,tRHS,SIZE(oRHS,1),SIZE(oRHS,2),reg,1)
      !$acc exit data copyout(tRHS) delete(tLHS)
      call checkRPtest(oRHS,tRHS,'tRHS (new cutensor)')
#endif

      call F2%copyfrom(F)
      write(*,*)
      write(*,*) 'Updating F with the result: original code'
      call UpdateFfromSoln(ZF,oRHS,tstmode)

      write(*,*) 'Updating F with the result: CPU'
      call PutSolninF_CP8(F,hRHS,tstmode,0)
      call NormBaseD_CP8(F,tstmode,.TRUE.,0)
      passtest=compareCPandCP8(ZF,F)
#if ACC_ENABLED
      write(*,*)
      write(*,*) 'Updating F with the result: GPU'
      call F2%copyintodevice()
      !$acc data copyin(hRHS)
      call PutSolninF_CP8(F2,hRHS,tstmode,1)
      !$acc end data
      call NormBaseD_CP8(F2,tstmode,.TRUE.,1)
      call F2%copyoutfromdevice()
      passtest=compareCPandCP8(ZF,F2)
#endif

      write(*,*)
      write(*,*) 'Normalizing G test: original code'
      call G2%flush()
      call G2%copyfrom(G)
      call FlushCP(ZG)
      call G%toCP(ZG)
      call Normalize(ZG)
      write(*,*)
      write(*,*) 'Normalizing G test: CPU'
      call Normalize_CP8(G,0)
      passtest=compareCPandCP8(ZG,G)
#if ACC_ENABLED
      write(*,*)
      write(*,*) 'Normalizing G test: GPU'
      call G2%copyintodevice()
      call Normalize_CP8(G2,1)
      call G2%copyoutfromdevice()
      passtest=compareCPandCP8(ZG,G2)
#endif

      call FlushALS(X)
      call FlushCP(ZHF)
      call FlushCP(ZWHF)
      call FlushCP(ZWG)
      call HTWH%flush()
      call HWG%flush()
      call HF%flush()
      call WHF%flush()
      call WG%flush()

      call FlushCP(ZF)
      call FlushCP(ZG)
      call FlushCP(ZH)
      call FlushCP(ZW)
      call FlushCP(ZWd)
      call F%flush()
      call F2%flush()
      call G%flush()
      call H%flush()
      call W%flush()
      call Wd%flush()

      deallocate(rows,cols,modes)

      end subroutine testls

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testaux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests auxiliary kernels

      implicit none
      TYPE (CP8) :: F,G
      TYPE (CP)  :: ZF,ZG
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: modes(:)
      logical :: passtest,t3
      integer :: rkf,rkg,rkh,rkw,ndof,i,j,ish,tstmode
      real(kind=8) :: stime,ftime,Esh,reg

      t3=.TRUE.
      ndof=3
      tstmode=3
      allocate(rows(ndof),cols(ndof),modes(ndof))
      modes(:)=.TRUE.
      rows(:)=5
      cols(:)=1
      rows(2)=7
      cols(2)=1
      rows(3)=4
      cols(3)=2
      rkf=4
      rkg=8

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'NORMBASE tests:'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)
      write(*,*) 'F:'
      call F%printvec
      call F%toCP(ZF)
      call G%copyfrom(F)

      write(*,*) 'Normalized F (existing alg)'
      call NORMBASE(ZF)
      call PrintCPvec(ZF) 

      write(*,*) 'CPU alg:'
      call NormBase_CP8(F,0)
      passtest=compareCPandCP8(ZF,F)

#if ACC_ENABLED
      write(*,*) 'GPU alg:'
      call G%copyintodevice()
      call NormBase_CP8(G,1)
      call G%copyoutfromdevice()
      passtest=compareCPandCP8(ZF,G)
#endif

      call F%flush()
      call G%flush()
      call FlushCP(ZF)

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'PrepG tests (new F generated):'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)
      call G%copyfrom(F)

      write(*,*) 'CPU alg: here is F'
      call PrepG_CP8(F,0)
      call F%printvec()
#if ACC_ENABLED
      write(*,*) 'GPU alg: here is G (should be same as F)'
      call G%copyintodevice()
      call PrepG_CP8(G,1)
      call G%copyoutfromdevice()
      call G%printvec
#endif

      end subroutine testaux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests rhs module

      implicit none
      TYPE (ALS) :: X
      TYPE (CP8) :: F,F1,F2,F3,G
      TYPE (CP)  :: ZF,ZG
      integer, allocatable :: rows(:),cols(:)
      logical, allocatable :: modes(:)
      logical :: passtest,t3
      integer :: rkf,rkg,ndof,i,j,ish,nals,conv
      real(kind=8) :: stime,ftime,Esh,reg

      nals=5

      ndof=4
      allocate(rows(ndof),cols(ndof),modes(ndof))
      modes(:)=.TRUE.
      rows(:)=25
      cols(:)=1
      rows(2)=17
      cols(2)=1
      rows(3)=14
      cols(3)=12
      rows(4)=4
      cols(4)=4

      rkf=10
      rkg=20

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ALS full sweep test:'
      write(*,*) '****************************************************'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      F%coef(:)=10*F%coef(:)
!      write(*,*) 'F (before ALS):'
!      call F%printvec
      call F1%copyfrom(F)
      call F2%copyfrom(F)
      call F3%copyfrom(F)
      call F%toCP(ZF)

      call G%newrand(rkg,rows,cols)
      call random_number(G%coef)
      G%coef(:)=12*G%coef(:)
!      write(*,*) 'G:'
!      call G%printvec
      call G%toCP(ZG)

      conv=ALS_reduce(ZF,ZG,nals)
      write(*,*) 'ZF, post-ALS (old alg):'
!      call PrintCPvec(ZF) 

      write(*,*) 'New als, host alg'
      call myals(F,G,nals,0)
      passtest=compareCPandCP8(ZF,F)
!      write(*,*) 'F on host, post-ALS:'
!      call F%printvec

      write(*,*) 'New als, device alg'
      call myals(F1,G,nals,1)
      passtest=compareCPandCP8(ZF,F1)
!      write(*,*) 'F on device, post-ALS:'
!      call F1%printvec

      write(*,*) 'New alsoo, host alg'
      call myalsoo(F2,G,nals,0)
      passtest=compareCPandCP8(ZF,F2)
!      write(*,*) 'F on host, post-ALS:'
!      call F2%printvec

      write(*,*) 'New alsoo, device alg'
      call myalsoo(F3,G,nals,1)
      passtest=compareCPandCP8(ZF,F3)
!      write(*,*) 'F on device, post-ALS:'
!      call F3%printvec

      call F%flush()
      call F1%flush()
      call F2%flush()
      call F3%flush()
      call G%flush()
      call FlushCP(ZF)
      call FlushCP(ZG)

      end subroutine testals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine myalsoo(F,G,nals,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests auxiliary kernels

      implicit none
      TYPE (CP8), intent(inout) :: F,G
      TYPE (CP8)  :: dummyA,dummyW,dummyF
      TYPE (ALS8) :: X
      integer, intent(in) :: nals,algo
      integer :: i,d,ndof,conv

      if (algo.eq.1) then
#if ACC_ENABLED
         call F%copyintodevice()
         call G%copyintodevice()
#else
         call AbortWithError('testals(): OpenACC alg not implemented')
#endif
      endif

      call X%new(dummyA,0,0.d0,F,G,dummyW,.FALSE.,.FALSE.,.FALSE.,algo)
      call X%setname('TestALS')
      call X%setlogical('chkconv',.TRUE.)
      call X%setlogical('prtconv',.TRUE.)
      ndof=F%D()
      conv=0
      do i=1,nals
         do d=1,ndof
            call X%resetconv()
            call X%accumprods(F,F,dummyF,G)
            conv=X%checkconv(F,G,G)
            call X%downdateBP(F,F,dummyF,G,d)
            call X%constls(d)
            call X%solvels(F,d)
            call X%updateBP(F,F,dummyF,G,d)
         enddo
      enddo

      call X%resetconv()
      call X%accumprods(F,F,dummyF,G)
      conv=X%checkconv(F,G,G)

      if (algo.eq.1) then
#if ACC_ENABLED
         call F%copyoutfromdevice()
         call G%deletefromdevice()
#endif
      endif

      call X%flush()

      end subroutine myalsoo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine myals(F,G,nals,algo)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests auxiliary kernels

      implicit none
      TYPE (CP8), intent(inout) :: F,G
      TYPE (CP8) :: G2,dummyH,dummyW
      TYPE (VV8) :: VB,VP
      TYPE (LS8) :: LS
      integer, intent(in) :: nals,algo
      real(kind=8), allocatable :: B(:),P(:),lhs(:),rhs(:)
      logical :: passtest,t3
      integer :: i,d,ndof,rF,rG,m,n
      real(kind=8) :: stime,ftime,reg

      ndof=F%D()
      rF=F%R()
      rG=G%R()
      reg=1.d-15

      call G2%copyfrom(G)

      call LS%new(dummyH,F,G,dummyW,0,0.0,algo)
      call VB%new(F,F,algo)
      call VP%new(G,F,algo)

      ALLOCATE(B(rF*rF),P(rG*rF))

      if (algo.eq.0) then
         
      elseif (algo.eq.1) then
#if ACC_ENABLED
         call F%copyintodevice()
         call G%copyintodevice()
         call G2%copyintodevice()
         !$acc enter data create(B,P)
#else
         call AbortWithError('testals(): OpenACC alg not implemented')

#endif
      endif

      call PrepG_CP8(G2,algo)
      call CONSTPT_CP8(VB,F,F,0,B)
      call CONSTPT_CP8(VP,G,F,0,P)

      do i=1,nals

         do d=1,ndof
            m=F%M(d)
            n=F%N(d)
            ALLOCATE(rhs(rF*m*n))
            ! acc enter data create(rhs)

!           Downdate B,P
            call UPDATE_P_CP8(VB,F,F,d,B,.TRUE.)
            call UPDATE_P_CP8(VP,G,F,d,P,.TRUE.)

!           Build LHS,RHS
            call GetLHS_CP8(LS,B,lhs)
            call GetRHS_CP8(LS,P,G2,rhs,d)

!           Solve linear system and update F
            call SolveLinSysLU8(lhs,rhs,rF,m*n,reg,algo)
            call PutSolninF_CP8(F,rhs,d,algo)
            call NormBaseD_CP8(F,d,.TRUE.,algo)

!           Update B,P
            call UPDATE_P_CP8(VB,F,F,d,B,.FALSE.)
            call UPDATE_P_CP8(VP,G,F,d,P,.FALSE.)

            !$acc exit data delete(lhs,rhs)
            DEALLOCATE(lhs,rhs)
         enddo

      enddo

      if (algo.eq.1) then
#if ACC_ENABLED
         !$acc exit data delete(B,P)
         call F%copyoutfromdevice()
         call G%deletefromdevice()
         call G2%deletefromdevice()
#endif
      endif

      DEALLOCATE(B,P)
      call VB%flush()
      call VP%flush()
      call LS%flush()
      call G2%flush()

      end subroutine myals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testals8all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests rhs module

      implicit none
      TYPE (ALS) :: X
      TYPE (CP8) :: A,G,W,F,F1,pA,AF
      TYPE (CP)  :: ZA,ZG,ZW,ZF,ZF1
      TYPE (MM8) :: tAF
      integer, allocatable :: rows(:),cols(:)
      logical :: passtest
      integer :: rkpA,rkw,rkf,rkg,ndof,i,j,ish,nals,conv
      real(kind=8) :: stime,ftime,Esh,reg,rq
      character(len=64) :: solver

      nals=3

      ndof=4
      allocate(rows(ndof),cols(ndof))
      cols(:)=1
      rows(:)=3
      rows(2)=7
      rows(3)=8
      rows(4)=4

      rkpA=2 ! rkA=rkpA**2 (+ shifted term, if applicable)
      rkw=2
      rkf=3
      rkg=6

      write(*,*) '****************************************************'
      write(*,'(/X,A/)') 'ALS8 module test:'
      write(*,*) '****************************************************'

      write(solver,*) 'lu'
      call set_als_settings_ALS(1.d-12,solver)
      call set_als_settings_ALS8(1.d-12,solver)
      call set_als_settings_LinSolver(1.d-12,solver)

      write(*,*) 'Creating Operator A...'
      call pA%newrand(rkpA,rows,rows)
      A%coef(:)=7*A%coef(:)
      call CPMM_CP8(pA,.FALSE.,pA,.TRUE.,A,.FALSE.,0)
      call pA%flush()
      call A%toCP(ZA)
!      write(*,*) 'Here is A (unshifted):'
!      call A%printmat()
      ish=1
      Esh=5.5

      write(*,*) 'Creating weights W...'
      call W%newrand(rkW,rows,rows)
      W%coef(:)=3*W%coef(:)
      call W%zerooffdiagonalall()
      W%base=abs(W%base)
      call W%toCP(ZW)
!      write(*,*) 'Here is W:'
!      call W%printmat()

      write(*,*) 'Creating vector F...'
      call F%newrand(rkf,rows,cols)
      call random_number(F%coef)
      call Normalize_CP8(F,0)
!      write(*,*) 'F (before ALS):'
!      call F%printvec

      write(*,*) 'Creating vector G...'
      call G%newrand(rkg,rows,cols)
      call random_number(G%coef)
      call Normalize_CP8(G,0)
      G%coef(:)=12*G%coef(:)
      call G%toCP(ZG)
!      write(*,*) 'G:'
!      call G%printvec

      write(*,*) '**************************'
      write(*,'(X,A)') 'unweighted ALS:'
      write(*,*) '**************************'

      write(*,*) 'ALS test (no wts, orig):'
      call F%toCP(ZF1)
      conv=ALS_reduce(ZF1,ZG,nals,'ALS test (no wts orig)')

      write(*,*) 'ALS test (unweighted, host):'
      call F1%copyfrom(F)
      conv=ALS_reduce_CP8(F1,G,nals,0,'ALS test (unweighted, host)')
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()

      write(*,*) 'ALS test (unweighted, device):'
      call F1%copyfrom(F)
      call F1%copyintodevice()
      call G%copyintodevice()
      conv=ALS_reduce_CP8(F1,G,nals,1,'ALS test (unweighted, device)')
      call G%deletefromdevice()
      call F1%copyoutfromdevice()
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()   
      call FlushCP(ZF1)

      write(*,*) '**************************'
      write(*,'(X,A)') 'weighted ALS:'
      write(*,*) '**************************'

      write(*,*) 'ALS test (weighted, orig):'
      call F%toCP(ZF1)
      conv=ALS_reduce(ZF1,ZG,ZW,nals,'ALS test (weighted orig)')

      write(*,*) 'ALS test (weighted, host):'
      call F1%copyfrom(F)
      conv=ALS_reduce_CP8(F1,G,W,nals,0,'ALS test (weighted, host)')
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()

      write(*,*) 'ALS test (weighted, device):'
      call F1%copyfrom(F)
      call F1%copyintodevice()
      call G%copyintodevice()
      call W%copyintodevice()
      conv=ALS_reduce_CP8(F1,G,W,nals,1,'ALS test (weighted, device)')
      call W%deletefromdevice()
      call G%deletefromdevice()
      call F1%copyoutfromdevice()
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()
      call FlushCP(ZF1)

      write(*,*) '**************************'
      write(*,'(X,A)') 'unweighted LS:'
      write(*,*) '**************************'

      write(*,*) 'LS test (unweighted, orig):'
      call F%toCP(ZF1)
      conv=ALS_solve(ZA,ZF1,ZG,nals,ish,Esh,'LS test (unweighted orig)')

      write(*,*) 'LS test (unweighted, host):'
      call F1%copyfrom(F)
      conv=ALS_solve_CP8(A,ish,Esh,F1,G,nals,0,&
                         'LS test (unweighted, host)')
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()

      write(*,*) 'LS test (unweighted, device):'
      call F1%copyfrom(F)
      call F1%copyintodevice()
      call G%copyintodevice()
      call A%copyintodevice()
      conv=ALS_solve_CP8(A,ish,Esh,F1,G,nals,1,&
                          'LS test (unweighted, device)')
      call A%deletefromdevice()
      call G%deletefromdevice()
      call F1%copyoutfromdevice()
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()
      call FlushCP(ZF1)

      write(*,*) '**************************'
      write(*,'(X,A)') 'weighted LS:'
      write(*,*) '**************************'

      write(*,*) 'LS test (weighted, orig):'
      call F%toCP(ZF1)
      conv=ALS_solve(ZA,ZF1,ZG,ZW,nals,ish,Esh,'LS test (weighted orig)')

      write(*,*) 'LS test (weighted, host):'
      call F1%copyfrom(F)
      conv=ALS_solve_CP8(A,ish,Esh,F1,G,W,nals,0,&
                         'LS test (weighted, host)')
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()

      write(*,*) 'LS test (weighted, device):'
      call F1%copyfrom(F)
      call F1%copyintodevice()
      call G%copyintodevice()
      call A%copyintodevice()
      call W%copyintodevice()
      conv=ALS_solve_CP8(A,ish,Esh,F1,G,W,nals,1,&
                          'LS test (weighted, device)')
      call W%deletefromdevice()
      call A%deletefromdevice()
      call G%deletefromdevice()
      call F1%copyoutfromdevice()
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()
      call FlushCP(ZF1)

      write(*,*) '**************************'
      write(*,'(X,A)') 'fast unweighted LS:'
      write(*,*) '**************************'

      write(*,*) 'LS test (unweighted, orig):'
      call F%toCP(ZF1)
!      conv=ALS_solve(ZA,ZF1,ZG,nals,ish,Esh,'LS test (unweighted orig)')
!      call LinSolver_1(ZA,ZF1,ZG,nals,ish,Esh,.TRUE.)
      call LinSolver_2(ZA,ZF1,ZG,nals,ish,Esh,.TRUE.)

      write(*,*) 'LS test (fast unweighted, host):'
      call F1%copyfrom(F)
      conv=ALS_fastsolve_CP8(A,ish,Esh,F1,G,nals,0,&
                         'LS test (fast unweighted, host)')
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()

      write(*,*) 'LS test (fast unweighted, device):'
      call F1%copyfrom(F)
      call F1%copyintodevice()
      call G%copyintodevice()
      call A%copyintodevice()
      conv=ALS_fastsolve_CP8(A,ish,Esh,F1,G,nals,1,&
                          'LS test (fast unweighted, device)')
      call A%deletefromdevice()
      call G%deletefromdevice()
      call F1%copyoutfromdevice()
      passtest=compareCPandCP8(ZF1,F1)
      call F1%flush()
      call FlushCP(ZF1)

!      write(*,*) '**************************'
!      write(*,'(X,A)') 'fast weighted LS:'
!      write(*,*) '**************************'
!
!      write(*,*) 'LS test (weighted, orig):'
!      call F%toCP(ZF1)
!      conv=ALS_solve(ZA,ZF1,ZG,ZW,nals,ish,Esh,'LS test (weighted orig)')
!
!      write(*,*) 'LS test (fast weighted, host):'
!      call F1%copyfrom(F)
!      conv=ALS_fastsolve_CP8(A,ish,Esh,F1,G,W,nals,0,&
!                         'LS test (fast weighted, host)')
!      passtest=compareCPandCP8(ZF1,F1)
!      call F1%flush()
!
!      write(*,*) 'LS test (fast weighted, device):'
!      call F1%copyfrom(F)
!      call F1%copyintodevice()
!      call G%copyintodevice()
!      call A%copyintodevice()
!      call W%copyintodevice()
!      conv=ALS_fastsolve_CP8(A,ish,Esh,F1,G,W,nals,1,&
!                          'LS test (fast weighted, device)')
!      call W%deletefromdevice()
!      call A%deletefromdevice()
!      call G%deletefromdevice()
!      call F1%copyoutfromdevice()
!      passtest=compareCPandCP8(ZF1,F1)
!      call F1%flush()
!      call FlushCP(ZF1)

!      write(*,*) '**************************'
!      write(*,'(X,A)') 'ALS_pow:'
!      write(*,*) '**************************'

!      write(*,*) 'ALS_pow test (orig):'
!      call F%toCP(ZF1)
!      call ALS_pow(ZA,ZF1,nals,ish,Esh)

!      write(*,*) 'ALS_pow test (host):'
!      call F1%copyfrom(F)
!      rq=ALS_pow_CP8(A,ish,Esh,F1,nals,0)
!      passtest=compareCPandCP8(ZF1,F1)
!      call F1%flush()

!      write(*,*) 'ALS_pow test (device):'
!      call F1%copyfrom(F)
!      call F1%copyintodevice()
!      call A%copyintodevice()
!      rq=ALS_pow_CP8(A,ish,Esh,F1,nals,1)
!      call A%deletefromdevice()
!      call F1%copyoutfromdevice()
!      passtest=compareCPandCP8(ZF1,F1)
!      call F1%flush()
!      call FlushCP(ZF1)

      end subroutine testals8all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine checkRPtest(R,P,tag)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Runs test to compare R and P:

      implicit none
      real(kind=8), intent(in) :: R(:,:)
      real(kind=8), intent(in) :: P(:)
      character(len=*) :: tag
      integer :: rk1,rk2
      logical :: passtest

      rk1=SIZE(R,1)
      rk2=SIZE(R,2)

      write(*,*) '-------------'
      write(*,*) ' Test: ',TRIM(ADJUSTL(tag))
      write(*,*) '-------------'

      write(*,*) 'Host R array:'
      call printRarray(R)

      write(*,*) 'test P array for: ',TRIM(ADJUSTL(tag))
      call printParray(P,rk1,rk2)
      passtest=comparepvvarrays(R,P)
      if (passtest) then
         write(*,*) TRIM(ADJUSTL(tag)),': test passed!'
      else
         write(*,*) TRIM(ADJUSTL(tag)),': entries differ!'
      endif

      end subroutine checkRPtest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function comparepvvarrays(R,P) result(passed)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares all entries in arrays R (reference) and P

      implicit none
      real(kind=8), intent(in) :: R(:,:)
      real(kind=8), intent(in) :: P(:)
      integer i,j,rk1,rk2
      real(kind=8), parameter :: tol=1.d-12
      real(kind=8) :: rval,pval
      logical :: passed

      rk1=SIZE(R,1)
      rk2=SIZE(R,2)
      passed=.TRUE.

      do j=1,rk2
         do i=1,rk1
            rval=R(i,j)
            pval=P((j-1)*rk1+i)
            if (abs(rval-pval).gt.tol) then
               write(*,'(2(I5,X),3(ES11.4,x))') &
               i,j,rval,pval,rval-pval
               passed=.false.
            endif
         enddo
      enddo

      end function comparepvvarrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine printRarray(R)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares all entries in arrays R (reference) and P

      implicit none
      real(kind=8), intent(in) :: R(:,:)
      integer i,j,rk1,rk2
      character*64 :: frmt

      rk1=SIZE(R,1)
      rk2=SIZE(R,2)
      write(frmt,'(A,I0,A)') '(',rk2,'(X,f14.6))'
      do i=1,rk1
         write(*,frmt) (R(i,j),j=1,rk2)
      enddo

      end subroutine printRarray

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine printParray(P,rk1,rk2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares all entries in arrays R (reference) and P

      implicit none
      real(kind=8), intent(in) :: P(:)
      integer, intent(in) :: rk1,rk2
      integer i,j
      character*64 :: frmt

      write(frmt,'(A,I0,A)') '(',rk2,'(X,f14.6))'
      do i=1,rk1
         write(*,frmt) (P((j-1)*rk1+i),j=1,rk2)
      enddo

      end subroutine printParray

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine testmpicycle()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests CP8 cycling over MPI ranks

      implicit none
      TYPE (CP8) :: v,w
      integer :: reqs(4),stats(4),rows(2),cols(2)
      integer :: i,j,npass,rk,torank,fromrank,ierr,nreqs
      integer :: tst
      logical :: swapem

      nreqs=4
      rk=3
      rows(:)=2
      cols(:)=2
      call v%new(rk,rows,cols)
      v%coef(:)=10+REAL(mpirank+1)
      v%base(:)=REAL(mpirank+1)
      call w%copyfrom(v)

      npass=5
      swapem=.FALSE.
      do i=1,npass
         do j=1,mpinodes
!           Source and destination ranks
            torank=mod(mpirank+mpinodes-1,mpinodes)
            fromrank=mod(mpirank+1,mpinodes)

            if (swapem) then
               call MPI_ISendRecv_CP8(w,v,torank,fromrank,reqs,0)
            else
               call MPI_ISendRecv_CP8(v,w,torank,fromrank,reqs,0)
            endif

            IF (mpirank.eq.mpi_prnt_rank) THEN
               if (swapem) then
                  write(*,*) 'pass ',i,', cycle ',j,': using w on rank ',mpirank
                  call w%printvec()                  
               else
                  write(*,*) 'pass ',i,', cycle ',j,': using v on rank ',mpirank
                  call v%printvec()
               endif
            ENDIF

#if MPI_ENABLED
            call MPI_WAITALL(nreqs,reqs,MPI_STATUSES_IGNORE,ierr)
#endif
            if (ierr.ne.0) write(*,*) 'rank ',mpirank,': ierr = ',ierr

            swapem=(.not.swapem)
         enddo
      enddo

!     Self w->v in case odd ranks,passes leaves initial vec in w
      if (swapem) then
         call MPI_ISendRecv_CP8(w,v,mpirank,mpirank,reqs,0)
#if MPI_ENABLED
         call MPI_WAITALL(nreqs,reqs,MPI_STATUSES_IGNORE,ierr)
#endif
      endif

      IF (mpirank.eq.mpi_prnt_rank) THEN
         write(*,*) 'final v on rank ',mpirank
         call v%printvec()
      ENDIF
      call w%flush()

      call AbortWithError('done testing')

      end subroutine testmpicycle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TESTCPR8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
