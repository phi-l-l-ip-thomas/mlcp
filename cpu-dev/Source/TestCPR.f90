!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module TESTCPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests CP vectors with CP-in-rank representation

      USE ERRORTRAP
      USE UTILS
      USE LINALG
      USE SEPDREPN
      USE MODVECVEC
      USE CPMMM
      USE REDUCTION
      USE LINSOLVER
      USE ALSOO
      USE ALSDRVR
      USE CPMATH
      USE CP2CP

      contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Test_calcPk()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Tests CPR_calcPk()

      implicit none
      TYPE (CPR) :: v,v2
      TYPE (CP)  :: w,w2,P
      real*8, allocatable  :: wtmp(:,:),Pv(:,:),Pw(:,:)
      integer, allocatable :: rows(:),cols(:),srows(:),scols(:),srk(:)
      integer, allocatable :: srows2(:),scols2(:)
      integer, parameter   :: ndof=3,pdof=4,pdof2=3,crk=2
      integer :: d,k,i,j,rk,rk2,rprod,cprod,rprod2,cprod2
      logical :: dow2coef,dowcoef

      write(*,*) 'Test_calcPk()...'

      dow2coef=.TRUE.
      dowcoef=.TRUE.

      ALLOCATE(rows(ndof),cols(ndof),srk(ndof),srows(pdof),scols(pdof))
      ALLOCATE(srows2(pdof2),scols2(pdof2))
      rows=(/3,5,4/)
      cols(:)=1
      srk=(/3,3,3,4,4,4,4,4,5,5,5,5/)
      srows=(/3,2,2,3/)
      scols(:)=1
      srows2=(/2,4,3/)
      scols2(:)=1

!     Generate random CPR vec v, v2
!!! Terms are all very nearly orthogonal to one another...
!!! TO DO: synthesize a data set of terms with more realistic overlaps
      v=NewCPRgen(rows,cols,srows,scols,srk,crk,.TRUE.)
      v2=NewCPRgen(rows,cols,srows2,scols2,srk,crk,.TRUE.)

!     Calculate (ordinary) CP-repn of v,v2
      rprod=PRODUCT(srows)
      cprod=PRODUCT(scols)
      rk=rprod*cprod
      rprod2=PRODUCT(srows2)
      cprod2=PRODUCT(scols2)
      rk2=rprod2*cprod2

      w=NewCP(rk,rows,cols)
      w2=NewCP(rk2,rows,cols)

      wtmp=CP2Matrix(v%coef)
      IF (SIZE(wtmp,1).ne.rprod .or. &
          SIZE(wtmp,2).ne.cprod) THEN
          write(*,*) 'wtmp(',SIZE(wtmp,1),',',SIZE(wtmp,2),&
          ') dims must be (',rprod,',',cprod,&
          ') for coefs'
         call AbortWithError('Test_calcPk(): bad wtmp dims')
      ENDIF
      w%coef(:)=RESHAPE(wtmp,(/rk/))

      DEALLOCATE(wtmp)
      DO d=1,ndof
         DO k=v%ibas(d),v%fbas(d)
            wtmp=CP2Matrix(v%base(k))
            IF (SIZE(wtmp,1).ne.rprod .or. &
                SIZE(wtmp,2).ne.cprod) THEN
               write(*,*) 'wtmp(',SIZE(wtmp,1),',',SIZE(wtmp,2),&
               ') dims must be (',rprod,',',cprod,&
               ') for base d = ',d,', k = ',k
               call AbortWithError('Test_calcPk(): bad wtmp dims')
            ENDIF
            w%base(k,:)=RESHAPE(wtmp,(/rk/))
            DEALLOCATE(wtmp)
         ENDDO
      ENDDO

      wtmp=CP2Matrix(v2%coef)
      IF (SIZE(wtmp,1).ne.rprod2 .or. &
          SIZE(wtmp,2).ne.cprod2) THEN
          write(*,*) 'wtmp(',SIZE(wtmp,1),',',SIZE(wtmp,2),&
          ') dims must be (',rprod2,',',cprod2,&
          ') for coefs'
         call AbortWithError('Test_calcPk(): bad wtmp2 dims')
      ENDIF
      w2%coef(:)=RESHAPE(wtmp,(/rk2/))

      DEALLOCATE(wtmp)
      DO d=1,ndof
         DO k=v2%ibas(d),v2%fbas(d)
            wtmp=CP2Matrix(v2%base(k))
            IF (SIZE(wtmp,1).ne.rprod2 .or. &
                SIZE(wtmp,2).ne.cprod2) THEN
               write(*,*) 'wtmp(',SIZE(wtmp,1),',',SIZE(wtmp,2),&
               ') dims must be (',rprod2,',',cprod2,&
               ') for base d = ',d,', k = ',k
               call AbortWithError('Test_calcPk(): bad wtmp dims')
            ENDIF
            w2%base(k,:)=RESHAPE(wtmp,(/rk2/))
            DEALLOCATE(wtmp)
         ENDDO
      ENDDO

!     Calc Bk of each for a mode and compare
      DO k=1,ndof
         write(*,*) 'CPR_calcPk(), k = ',k
!         P=CPR_calcPk(v,v2,k) ! Bk = Pk = <v,v2>_k
         P=CPR_calcPtot(v,v2,k,.FALSE.,.FALSE.,100) ! ALT: P != k

         write(*,*) 'CP2Matrix(), k = ',k
         Pv=CP2Matrix(P) ! Bk as an ordinary matrix for v

         write(*,*) 'CONSTPk, k = ',k
         ALLOCATE(Pw(rk,rk2))
!         call CONSTPk(w2,w,k,Pw) ! as an ordinary matrix for v,w
         call CONSTPT(w2,w,k,Pw) ! ALT: calc P != k

         write(*,'(3(A,10X),A,I0)') 'Pv','Pw','Rel. diff',' for d = ',k
         call comparematrices(Pw,Pv)

         call FlushCP(P)
         DEALLOCATE(Pv,Pw)
      ENDDO

!     Compute Ptot, with or without coefs
      write(*,*) 'CPR_calcPtot(), all modes'
      P=CPR_calcPtot(v,v2,0,dowcoef,dow2coef,100)

      write(*,*) 'CP2Matrix(), all modes'
      Pv=CP2Matrix(P) ! Bk as an ordinary matrix for v

      write(*,*) 'CONSTPTOT()'
      ALLOCATE(Pw(rk,rk2))
      call CONSTPT(w2,w,0,Pw) ! ALT: calc P != k

!     For comparing coefs:
      IF (dowcoef) THEN
         write(*,*) 'multiply w coefs'
         DO i=1,rk
            Pw(i,:)=w%coef(i)*Pw(i,:)
         ENDDO
      ENDIF

      IF (dow2coef) THEN
         write(*,*) 'multiply w2 coefs'
         DO i=1,rk2
            Pw(:,i)=w2%coef(i)*Pw(:,i)
         ENDDO
      ENDIF

      write(*,'(3(A,10X),A,I0)') 'Pv','Pw','Rel. diff',' (all modes)'
      call comparematrices(Pw,Pv)

      write(*,*) 'Pw(1,1) Pw(1,1)[8]  Pw(1,1)[4]'
      write(*,*) Pw(1,1),REAL(Pw(1,1),8),REAL(Pw(1,1),4)

      call FlushCP(P)
      DEALLOCATE(Pv,Pw)

      DEALLOCATE(rows,cols,srk,srows,scols)

      call AbortWithError('Done Test_calcPk()')

      end subroutine Test_calcPk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function GetSyntheticCPR(rows,cols,srows,scols,srk,crk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates "synthetic" data set

      implicit none
      TYPE (CPR) :: v
      integer, intent(in) :: rows(:),cols(:),srows(:),scols(:),srk(:)
      integer, intent(in) :: crk
      integer :: d,ndof,i
      real*8  :: norm,norms

!     Initialize the CPR type but do not fill with data yet
      v=NewCPRgen(rows,cols,srows,scols,srk,crk)

!     Coefs are still random
      call random_number(v%coef%coef(:))

!     Generate initially parallel terms; make small, random variations
      ndof=V%D()
      DO d=1,ndof
         norm=1.d0/sqrt(REAL(v%rows(i)*v%cols(i)))
         DO i=v%ibas(d),v%fbas(d)
!!! UC
!            norms=
            v%base(i)%coef(:)=1.d0
!            v%base(i)%base(, )=norms
!!! EC
         ENDDO
      ENDDO

      end function GetSyntheticCPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Test_randomnorm(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generate a random vec w (with dims of v) and check the norm

      implicit none
      TYPE (CP), INTENT(IN) :: v
      TYPE (CP) :: w
      real*8, allocatable :: BB(:,:)
      integer :: i,d,ntrials,rk
      real*8  :: accum,maxip,minip,tmp

      ntrials=1000
      accum=0.d0
      maxip=0.d0
      minip=1.d99
      rk=14

      ALLOCATE(BB(rk,rk))

      DO i=1,ntrials
         w=RandomCPref_A(v,rk)
         tmp=sqrt(abs(PRODVV(w)))
         accum=accum+tmp
         maxip=MAX(maxip,tmp)
         minip=MIN(minip,tmp)
         IF (i.eq.ntrials) THEN
            call CONSTPT(w,0,BB)
            write(*,*) 'BBtot (last vector)'
            call PrintMatrix(BB)
         ENDIF
         call FlushCP(w)
      ENDDO

      write(*,*) 'Test_randomnorm(): mean inner prod is ',&
      accum/REAL(ntrials),' after ',ntrials,' trials.'
      write(*,*) 'Test_randomnorm(): max  inner prod is ',&
      maxip,' after ',ntrials,' trials.'
      write(*,*) 'Test_randomnorm(): min  inner prod is ',&
      minip,' after ',ntrials,' trials.'


      DEALLOCATE(BB)

      call AbortWithError('Done with Test_randomnorm()')

      end subroutine Test_randomnorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function RandomCPref_A(w,rk) result(v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Builds a CP item with random entries, with dimensions of reference w

      implicit none
      TYPE (CP) :: v
      TYPE (CP), INTENT(IN) :: w
      INTEGER, INTENT(IN), OPTIONAL :: rk
      INTEGER :: i,d,rv
      REAL*8  :: fac

      IF (present(rk)) THEN
         rv=rk
      ELSE
         rv=w%R()
      ENDIF

!     Generate v with random entries and equal coefs for all terms
      v=NewCP(rv,w%rows,w%cols,w%sym)
      v%coef(:)=1.d0/sqrt(REAL(rv))
      call random_number(v%base(:,:))

!     Shift, scale entries for each mode to make rms norm ~ unity
      DO d=1,v%D()
         fac=sqrt(12.d0/REAL(v%nbas(d)))
         v%base(v%ibas(d):v%fbas(d),:)=fac*&
         (v%base(v%ibas(d):v%fbas(d),:)-0.5d0)
      ENDDO

      end function RandomCPref_A


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine comparematrices(M1,M2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares matrices entrywise. M1 is treated as the reference for
! relative difference

      implicit none
      real*8, intent(in) :: M1(:,:),M2(:,:)
      integer :: i,j,r,c

      r=SIZE(M1,1)
      c=SIZE(M1,2)

      IF (SIZE(M2,1).ne.r .or. SIZE(M2,2).ne.c) THEN
         write(*,*) 'M2 is [',SIZE(M2,1),' x ',SIZE(M2,2),&
         '] but must be [',r,' x ',c,'] (same as M1)'
         call AbortWithError('comparematrices(): M1, M2 size mismatch')
      ENDIF

      write(*,'(2(9X,A),2(24X,A),2(18X,A))') &
           'col','row','M1','M2',&
           'Abs. dif','Rel. dif'
      DO j=1,c
         DO i=1,r
            IF (M1(i,j).ne.0.d0) THEN
               write(*,*) j,i,M1(i,j),M2(i,j),&
                          M1(i,j)-M2(i,j),&
                          abs((M1(i,j)-M2(i,j))/M1(i,j))
            ELSE
               write(*,*) j,i,M1(i,j),M2(i,j),&
                          M1(i,j)-M2(i,j)
            ENDIF
         ENDDO
      ENDDO

      end subroutine comparematrices

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module TESTCPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
