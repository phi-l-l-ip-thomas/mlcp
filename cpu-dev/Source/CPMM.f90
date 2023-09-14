!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE CPMMM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      USE ERRORTRAP
      USE UTILS
      USE MODVECVEC
      USE SEPDREPN

      implicit none
      real*8, private  :: mvp_time=0.d0
      logical, private :: MVP_SETUP = .FALSE.

      INTERFACE CPMM
         MODULE PROCEDURE CPMM_noshiftall
         MODULE PROCEDURE CPMM_shiftedall
         MODULE PROCEDURE CPMM_noshiftsub
         MODULE PROCEDURE CPMM_allmodes
         MODULE PROCEDURE CPMM_1mode_allterms
         MODULE PROCEDURE CPMM_1mode
         MODULE PROCEDURE CPMM_fillM3
         MODULE PROCEDURE CPMM_general
      END INTERFACE CPMM

      INTERFACE CPMM_vec
         MODULE PROCEDURE CPMM_vec_noshiftall
         MODULE PROCEDURE CPMM_vec_shiftedall
         MODULE PROCEDURE CPMM_vec_noshiftsub
         MODULE PROCEDURE CPMM_vec_allmodes
         MODULE PROCEDURE CPMM_vec_1mode_allterms
         MODULE PROCEDURE CPMM_vec_1mode
         MODULE PROCEDURE CPMM_vec_fillM3
         MODULE PROCEDURE CPMM_vec_general
      END INTERFACE CPMM_vec

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine InitializePRODHVModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      mvp_time = 0.d0
      MVP_SETUP = .TRUE.

      end subroutine InitializePRODHVModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DisposePRODHVModule()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      MVP_SETUP = .FALSE.
      write(*,'(X,A,X,f20.3)') 'Total matrix-vector product time  (s)',&
                            mvp_time

      end subroutine DisposePRODHVModule

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_noshiftall(M1,t1,M2,t2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      logical, intent(in) :: t1,t2
      integer, parameter :: i=0
      real*8, parameter  :: E=0.d0

      call CPMM(M1,i,E,t1,M2,i,E,t2,M3)

      end subroutine CPMM_noshiftall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_shiftedall(M1,ishift1,Eshift1,t1,&
                                 M2,ishift2,Eshift2,t2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Each of M1, M2 can
! be shifted by E*I.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ishift1,ishift2
      logical, intent(in) :: t1,t2

      call CPMM(M1,1,SIZE(M1%coef),ishift1,Eshift1,t1,&
                M2,1,SIZE(M2%coef),ishift2,Eshift2,t2,M3)

      end subroutine CPMM_shiftedall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_noshiftsub(M1,ri1,re1,t1,M2,ri2,re2,t2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Subsets of ranks
! of M1 and M2 can be selected

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      integer, intent(in) :: ri1,re1,ri2,re2
      logical, intent(in) :: t1,t2
      integer, parameter :: i=0
      real*8, parameter  :: E=0.d0

      call CPMM(M1,ri1,re1,i,E,t1,M2,ri2,re2,i,E,t2,M3)

      end subroutine CPMM_noshiftsub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_allmodes(M1,ri1,re1,ishift1,Eshift1,t1,&
                               M2,ri2,re2,ishift2,Eshift2,t2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. All modes are
! multiplied. Subsets of ranks of M1 and M2 can be selected.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ri1,re1,ishift1,ri2,re2,ishift2
      logical, intent(in)  :: t1,t2
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=MAX(M1%D(),M2%D())

      allocate(modes(ndof))
      modes(:)=.TRUE.
      call CPMM(M1,ri1,re1,ishift1,Eshift1,t1,&
                M2,ri2,re2,ishift2,Eshift2,t2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_allmodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_1mode_allterms(M1,ishift1,Eshift1,t1,&
                                     M2,ishift2,Eshift2,t2,&
                                     M3,idof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Only mode idof is
! multiplied. All R terms are multiplied.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ishift1,ishift2,idof
      logical, intent(in)  :: t1,t2

      call CPMM(M1,1,SIZE(M1%coef),ishift1,Eshift1,t1,&
                M2,1,SIZE(M2%coef),ishift2,Eshift2,t2,&
                M3,idof)

      end subroutine CPMM_1mode_allterms

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_1mode(M1,ri1,re1,ishift1,Eshift1,t1,&
                            M2,ri2,re2,ishift2,Eshift2,t2,&
                            M3,idof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Only mode idof is
! multiplied. Subsets of ranks of M1 and M2 can be selected.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ri1,re1,ishift1,ri2,re2,ishift2,idof
      logical, intent(in)  :: t1,t2
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=MAX(M1%D(),M2%D())

      IF (idof.lt.1 .or. idof.gt.ndof) THEN
         write(*,*) 'idof (',idof,') must be in range [1,',ndof,']'
         call AbortWithError('Error in CPMM()')
      ENDIF

      allocate(modes(ndof))
      modes(:)=.FALSE.
      modes(idof)=.TRUE.
      call CPMM(M1,ri1,re1,ishift1,Eshift1,t1,&
                M2,ri2,re2,ishift2,Eshift2,t2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_1mode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_fillM3(M1,ri1,re1,ishift1,Eshift1,t1,&
                             M2,ri2,re2,ishift2,Eshift2,t2,&
                             M3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created (or
! overwritten, if already present) as a result. Set tx (x={1,2}) to 
! .TRUE. to do transpose; add 'modes' argument to select a mode subset.
! This version assumes that M3 has the same size as M1*M2 (where M1,M2
! may be shifted)

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ri1,re1,ishift1,ri2,re2,ishift2
      logical, intent(in) :: t1,t2
      logical, intent(in) :: modes(:)

      call CPMM(M1,ri1,re1,ishift1,Eshift1,t1,&
                M2,ri2,re2,ishift2,Eshift2,t2,&
                M3,1,modes) 

      end subroutine CPMM_fillM3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_general(M1,ri1,re1,ishift1,Eshift1,t1,&
                              M2,ri2,re2,ishift2,Eshift2,t2,&
                              M3,ri3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created (or
! overwritten, if already present) as a result. Set tx (x={1,2}) to 
! .TRUE. to do transpose; add 'modes' argument to select a mode subset.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      TYPE (CP) :: I1,I2
      integer, allocatable, dimension(:) :: ld1,ld2,r1,c1,r2,c2,&
                                            subm,resm,r1sub,c2sub
      integer :: d,h,i,j,k,ik,rk1,rk2,rk3,re3,ndof,nsubm
      real*8  :: alpha1,alpha2,afac,ti1,ti2
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ri1,re1,ishift1,ri2,re2,ishift2,ri3
      logical, intent(in) :: t1,t2
      logical, intent(in) :: modes(:)
      logical, allocatable :: symsub(:)
      character(1) :: tN1,tN2            

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      call CPU_TIME(ti1)

!     Set parameters: # modes, ranks
      ndof=MAX(M1%D(),M2%D())
      nsubm=COUNT(modes)
      rk1=re1-ri1+1
      rk2=re2-ri2+1
      rk3=rk1*rk2

!     Alpha multiplies the matrix product and is always opposite in sign
!     to the energy shift
      alpha1=1.d0
      alpha2=1.d0

!     Modify ranks for shifts and generate shift terms
      IF (ishift1.ne.0) THEN
         rk3=rk3+rk2
         IF (ishift1.lt.0) alpha1=-alpha1
         I1=IdentityCPMatrix(M1%rows,M1%cols,M1%sym)
         call VecScalarMult(I1,-alpha1*Eshift1,1,1)
      ENDIF
      IF (ishift2.ne.0) THEN
         rk3=rk3+rk1
         IF (ishift2.lt.0) alpha2=-alpha2
         I2=IdentityCPMatrix(M2%rows,M2%cols,M2%sym)
         call VecScalarMult(I2,-alpha2*Eshift2,1,1)
      ENDIF
      IF (ishift1.ne.0 .and. ishift2.ne.0) rk3=rk3+1
      re3=ri3+rk3-1

!     Error checking
      IF (M1%D().lt.1 .or. M2%D().lt.1) THEN
         write(*,*) 'M1, M2 have ',M1%D(),', ',M2%D(),&
                    ' modes, respectively; must both be > 0'
         call AbortWithError('CPMM(): # modes in M1, M2 must be > 0')
      ENDIF

      IF (M3%R().gt.0) THEN
         IF (M3%D().ne.ndof) THEN
            write(*,*) '(existing) M3 has ',M3%D(),&
            ' modes but must have ',ndof,' (max of M1,M2) modes'
            call AbortWithError('CPMM(): wrong mode size for M3')
         ENDIF
         IF (M3%R().lt.re3) THEN
            write(*,*) 'M3 has rank ',M3%R(),&
            ', but must be at least ',re3
            call AbortWithError('CPMM(): rk_M3 < re3')
         ENDIF
      ENDIF

      IF (SIZE(modes).ne.ndof) THEN
            write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',ndof
            call AbortWithError('Error in CPMM()')
      ENDIF

      IF (nsubm.eq.0) THEN
         write(*,*) 'At least one mode must be selected'
         call AbortWithError('Error in CPMM()')
      ENDIF

      IF (ANY(M1%sym) .or. ANY(M2%sym)) &
         call AbortWithError('CPMM(): symmetry not yet implemented')

      IF (ri1.lt.1 .or. re1.gt.M1%R() .or. ri1.gt.re1 .or. &
          ri2.lt.1 .or. re2.gt.M2%R() .or. ri2.gt.re2 .or. &
          ri3.lt.1) THEN
          write(*,*) 'Rank index out-of-range: ',&
          'ri1,re1 = ',ri1,',',re1,' for M1(1:',M1%R(),'); ',&
          'ri2,re2 = ',ri2,',',re2,' for M2(1:',M2%R(),');' ,&
          'ri3,re3 = ',ri3,',',re3,' for M3(1:',rk3,')'
          call AbortWithError('CPMM(): bad rank indices')
      ENDIF

!     Prepare DGEMM arguments
!     Set entries initially to 1 in case number of modes mismatches
      ALLOCATE(ld1(ndof),ld2(ndof),r1(ndof),c1(ndof),r2(ndof),c2(ndof))
      ld1(:)=1
      r1(:)=1
      c1(:)=1
      ld2(:)=1
      r2(:)=1
      c2(:)=1
      ld1(1:M1%D())=M1%rows(1:M1%D())
      ld2(1:M2%D())=M2%rows(1:M2%D())
      IF (t1) THEN
         tN1='T'
         r1(1:M1%D())=M1%cols(1:M1%D())
         c1(1:M1%D())=M1%rows(1:M1%D())
      ELSE
         tN1='N'
         r1(1:M1%D())=M1%rows(1:M1%D())
         c1(1:M1%D())=M1%cols(1:M1%D())
      ENDIF
      IF (t2) THEN
         tN2='T'
         r2(1:M2%D())=M2%cols(1:M2%D())
         c2(1:M2%D())=M2%rows(1:M2%D())
      ELSE
         tN2='N'
         r2(1:M2%D())=M2%rows(1:M2%D())
         c2(1:M2%D())=M2%cols(1:M2%D())
      ENDIF

!     Map from modes to subset-of-modes, and subset dimensions
      ALLOCATE(subm(nsubm),resm(nsubm))
      ALLOCATE(r1sub(nsubm),c2sub(nsubm),symsub(nsubm))
      symsub(:)=.FALSE.
      d=1
      DO j=1,ndof
         IF (modes(j)) THEN
            subm(d)=j
            IF (M3%R().gt.0) THEN
               resm(d)=j
            ELSE
               resm(d)=d
            ENDIF
            r1sub(d)=r1(j)
            c2sub(d)=c2(j)
            d=d+1
         ENDIF
      ENDDO

      IF (M3%R().eq.0) M3=NewCP(re3,r1sub,c2sub,symsub)
      deallocate(r1sub,c2sub,symsub)

!     Loop over modes
      DO d=1,nsubm
         j=subm(d)
         h=resm(d)

         IF (c1(j).ne.r2(j)) THEN
            write(*,*) 'Mode ',j,': M1 cols (',c1(j),&
                       '), M2 rows (',r2(j),') must match'
            call AbortWithError('Error in CPMM()')
         ENDIF

         IF (r1(j).ne.M3%rows(h) .or. c2(j).ne.M3%cols(h)) THEN
            write(*,*) 'Mode ',j,': M3 must be (',r1(j),' x ',c2(j),&
            '), but is (',M3%rows(h),' x ',M3%cols(h),')'
            call AbortWithError('Error in CPMM()')
         ENDIF

!        Loop over rank terms M1
         ik=ri3
         DO k=ri1,re1 ! k=1,rk1
!           Loop over rank terms in M2
            DO i=ri2,re2 ! i=1,rk2
               IF (d.eq.1) M3%coef(ik)=M1%coef(k)*M2%coef(i)
               IF (j.eq.1) THEN
                  afac=alpha1*alpha2
               ELSE
                  afac=1.d0
               ENDIF

!              j references a nonexistant mode in M1, so copy M2->M3
!              afac is not needed since j must be greater than 1
               IF (j.gt.M1%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=& 
                  M2%base(M2%ibas(j):M2%fbas(j),i)
!              j references a nonexistant mode in M2, so copy M1->M3
!              afac is not needed since j must be greater than 1
               ELSEIF (j.gt.M2%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=&
                  M1%base(M1%ibas(j):M1%fbas(j),k)
               ELSE
                  call DGEMM(tN1,tN2,r1(j),c2(j),c1(j),afac,&
                          M1%base(M1%ibas(j):M1%fbas(j),k),ld1(j),&
                          M2%base(M2%ibas(j):M2%fbas(j),i),ld2(j),0.d0,&
                          M3%base(M3%ibas(h):M3%fbas(h),ik),r1(j))
               ENDIF
               ik=ik+1
            ENDDO

!           Shift on M2
            IF (ishift2.ne.0) THEN
               IF (d.eq.1) M3%coef(ik)=M1%coef(k)*I2%coef(1)
               IF (j.eq.1) THEN
                  afac=alpha1
               ELSE
                  afac=1.d0
               ENDIF

               IF (j.gt.M1%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=&
                  I2%base(I2%ibas(j):I2%fbas(j),1)
               ELSEIF (j.gt.M2%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=&
                  M1%base(M1%ibas(j):M1%fbas(j),k)
               ELSE
                  call DGEMM(tN1,tN2,r1(j),c2(j),c1(j),afac,&
                          M1%base(M1%ibas(j):M1%fbas(j),k),ld1(j),&
                          I2%base(I2%ibas(j):I2%fbas(j),1),ld2(j),0.d0,&
                          M3%base(M3%ibas(h):M3%fbas(h),ik),r1(j))
               ENDIF
               ik=ik+1
            ENDIF        
         ENDDO

!        Shift on M1
         IF (ishift1.ne.0) THEN
            DO i=ri2,re2 ! i=1,rk2
               IF (d.eq.1) M3%coef(ik)=I1%coef(1)*M2%coef(i)
               IF (j.eq.1) THEN
                  afac=alpha2
               ELSE
                  afac=1.d0
               ENDIF

               IF (j.gt.M1%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=&
                  M2%base(M2%ibas(j):M2%fbas(j),i)
               ELSEIF (j.gt.M2%D()) THEN
                  M3%base(M3%ibas(h):M3%fbas(h),ik)=&
                  I1%base(I1%ibas(j):I1%fbas(j),1)
               ELSE
                  call DGEMM(tN1,tN2,r1(j),c2(j),c1(j),afac,&
                          I1%base(I1%ibas(j):I1%fbas(j),1),ld1(j),&
                          M2%base(M2%ibas(j):M2%fbas(j),i),ld2(j),0.d0,&
                          M3%base(M3%ibas(h):M3%fbas(h),ik),r1(j))
               ENDIF
               ik=ik+1
           ENDDO
         ENDIF

!        Shift on both M1 and M2
         IF (ishift1.ne.0 .and. ishift2.ne.0) THEN
            IF (d.eq.1) M3%coef(ik)=I1%coef(1)*I2%coef(1)

            IF (j.gt.M1%D()) THEN
               M3%base(M3%ibas(h):M3%fbas(h),ik)=&
               I2%base(I2%ibas(j):I2%fbas(j),1)
            ELSEIF (j.gt.M2%D()) THEN
               M3%base(M3%ibas(h):M3%fbas(h),ik)=&
               I1%base(I1%ibas(j):I1%fbas(j),1)
            ELSE
               call DGEMM(tN1,tN2,r1(j),c2(j),c1(j),1.d0,&
                    I1%base(I1%ibas(j):I1%fbas(j),1),ld1(j),&
                    I2%base(I2%ibas(j):I2%fbas(j),1),ld2(j),0.d0,&
                    M3%base(M3%ibas(h):M3%fbas(h),ik),r1(j))
            ENDIF
         ENDIF
      ENDDO ! modes loop (j)

      DEALLOCATE(ld1,ld2,r1,c1,r2,c2,subm,resm)
      call FlushCP(I1)
      call FlushCP(I2)

      call CPU_TIME(ti2)
      mvp_time=mvp_time+ti2-ti1

      end subroutine CPMM_general

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_noshiftall(M1,M2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. M1,M2,M3 are passed as vectors assumed to be the diagonal 
! entries of these matrices, respectively.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      integer, parameter :: i=0
      real*8, parameter  :: E=0.d0

      call CPMM_vec(M1,i,E,M2,i,E,M3)

      end subroutine CPMM_vec_noshiftall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_shiftedall(M1,ishift1,Eshift1,&
                                     M2,ishift2,Eshift2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Each of M1, M2 can be shifted by E*I. M1,M2,M3 are passed as 
! vectors assumed to be the diagonal entries of these matrices, resp.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ishift1,ishift2

      call CPMM_vec(M1,1,SIZE(M1%coef),ishift1,Eshift1,&
                    M2,1,SIZE(M2%coef),ishift2,Eshift2,M3)

      end subroutine CPMM_vec_shiftedall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_noshiftsub(M1,ri1,re1,M2,ri2,re2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Subsets of ranks of M1 and M2 can be selected. M1,M2,M3 are
! passed as vectors assumed to be the diagonal entries of these matrices
! respectively.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      integer, intent(in) :: ri1,re1,ri2,re2
      integer, parameter :: i=0
      real*8, parameter  :: E=0.d0

      call CPMM_vec(M1,ri1,re1,i,E,M2,ri2,re2,i,E,M3)

      end subroutine CPMM_vec_noshiftsub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_allmodes(M1,ri1,re1,ishift1,Eshift1,&
                                   M2,ri2,re2,ishift2,Eshift2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. All modes are
! multiplied. Subsets of ranks of M1 and M2 can be selected.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ri1,re1,ishift1,ri2,re2,ishift2
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=M1%D()

      allocate(modes(ndof))
      modes(:)=.TRUE.
      call CPMM_vec(M1,ri1,re1,ishift1,Eshift1,&
                    M2,ri2,re2,ishift2,Eshift2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_vec_allmodes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_1mode_allterms(M1,ishift1,Eshift1,&
                                         M2,ishift2,Eshift2,&
                                         M3,idof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Only mode idof is
! multiplied. All R terms are multiplied.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ishift1,ishift2,idof

      call CPMM_vec(M1,1,SIZE(M1%coef),ishift1,Eshift1,&
                    M2,1,SIZE(M2%coef),ishift2,Eshift2,&
                    M3,idof)

      end subroutine CPMM_vec_1mode_allterms

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_1mode(M1,ri1,re1,ishift1,Eshift1,&
                                M2,ri2,re2,ishift2,Eshift2,&
                                M3,idof)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. Set tx (x={1,2}) to .TRUE. to do transpose. Only mode idof is
! multiplied. Subsets of ranks of M1 and M2 can be selected.

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)   :: Eshift1,Eshift2
      integer, intent(in)  :: ri1,re1,ishift1,ri2,re2,ishift2,idof
      logical, allocatable :: modes(:)
      integer :: ndof

      ndof=M1%D()

      IF (idof.lt.1 .or. idof.gt.ndof) THEN
         write(*,*) 'idof (',idof,') must be in range [1,',ndof,']'
         call AbortWithError('Error in CPMM()')
      ENDIF

      allocate(modes(ndof))
      modes(:)=.FALSE.
      modes(idof)=.TRUE.
      call CPMM_vec(M1,ri1,re1,ishift1,Eshift1,&
                    M2,ri2,re2,ishift2,Eshift2,M3,modes)
      deallocate(modes)

      end subroutine CPMM_vec_1mode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_fillM3(M1,ri1,re1,ishift1,Eshift1,&
                                 M2,ri2,re2,ishift2,Eshift2,&
                                 M3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. This is the Hadamard product version which computes products  
! of diagonal-matrix-diagonal-matrix (so no transposing done here).
! This version assumes that M3 has the same size as M1*M2 (where M1,M2
! may be shifted)

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ri1,re1,ishift1,ri2,re2,ishift2
      logical, intent(in) :: modes(:)

      call CPMM_vec(M1,ri1,re1,ishift1,Eshift1,&
                    M2,ri2,re2,ishift2,Eshift2,&
                    M3,1,modes)

      end subroutine CPMM_vec_fillM3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMM_vec_general(M1,ri1,re1,ishift1,Eshift1,&
                                  M2,ri2,re2,ishift2,Eshift2,&
                                  M3,ri3,modes)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Applies CP-format matrix-matrix product M1*M2 = M3; M3 is created as a
! result. This is the Hadamard product version which computes products  
! of diagonal-matrix-diagonal-matrix (so no transposing done here)

      implicit none
      TYPE (CP), INTENT(IN)  :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      TYPE (CP) :: I1,I2
      real*8, intent(in)  :: Eshift1,Eshift2
      integer, intent(in) :: ri1,re1,ishift1,ri2,re2,ishift2,ri3
      logical, intent(in) :: modes(:)
      logical, allocatable :: symsub(:)
      integer, allocatable :: subm(:),resm(:),rsub(:),csub(:)
      integer :: d,h,i,j,k,ik,rk1,rk2,rk3,re3,ndof,nsubm
      real*8  :: alpha1,alpha2,afac,ti1,ti2

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      call CPU_TIME(ti1)

!     Set parameters: # modes, ranks
      ndof=M1%D()
      nsubm=COUNT(modes)
      rk1=re1-ri1+1
      rk2=re2-ri2+1
      rk3=rk1*rk2

!     Alpha multiplies the matrix product and is always opposite in sign
!     to the energy shift
      alpha1=1.d0
      alpha2=1.d0

!     Modify ranks for shifts and generate shift terms
      IF (ishift1.ne.0) THEN
         rk3=rk3+rk2
         IF (ishift1.lt.0) alpha1=-alpha1
         I1=NewCP(1,M1%nbas)
         I1%coef(:)=1.d0
         I1%base(:,:)=1.d0
         call VecScalarMult(I1,-alpha1*Eshift1,1,1)
      ENDIF
      IF (ishift2.ne.0) THEN
         rk3=rk3+rk1
         IF (ishift2.lt.0) alpha2=-alpha2
         I2=NewCP(1,M2%nbas)
         I2%coef(:)=1.d0
         I2%base(:,:)=1.d0
         call VecScalarMult(I2,-alpha2*Eshift2,1,1)
      ENDIF
      IF (ishift1.ne.0 .and. ishift2.ne.0) rk3=rk3+1
      re3=ri3+rk3-1

!     Error checking
      IF (.NOT.CHECKNBAS(M1,M2)) &
         call AbortWithError('CPMM_vec(): M1,M2 dimension mismatch')
      IF (allocated(M3%coef)) THEN
         IF (.NOT.CHECKNBAS(M1,M3)) &
            call AbortWithError('CPMM_vec(): M1,M3 dimension mismatch')
         IF (SIZE(M3%coef).lt.re3) THEN
            write(*,*) 'M3 has rank ',SIZE(M3%coef),&
            ', but must be at least ',re3
            call AbortWithError('CPMM_vec(): rk_M3 < re3')
         ENDIF
      ENDIF
      IF (SIZE(modes).ne.ndof) THEN
            write(*,*) 'SIZE(modes) = ',SIZE(modes),' must be ',ndof
            call AbortWithError('Error in CPMM_vec()')
      ENDIF
      IF (nsubm.eq.0) THEN
         write(*,*) 'At least one mode must be selected'
         call AbortWithError('Error in CPMM_vec()')
      ENDIF
      IF (ANY(M1%sym) .or. ANY(M2%sym)) &
         call AbortWithError('CPMM_vec(): symmetry not yet implemented')
      IF (ri1.lt.1 .or. re1.gt.SIZE(M1%coef) .or. ri1.gt.re1 .or. &
          ri2.lt.1 .or. re2.gt.SIZE(M2%coef) .or. ri2.gt.re2 .or. &
          ri3.lt.1) THEN
          write(*,*) 'Rank index out-of-range: ',&
          'ri1,re1 = ',ri1,',',re1,' for M1(1:',SIZE(M1%coef),'); ',&
          'ri2,re2 = ',ri2,',',re2,' for M2(1:',SIZE(M2%coef),'); ',&
          'ri3,re3 = ',ri3,',',re3,' for M3(1:',rk3,')'
          call AbortWithError('CPMM_vec(): bad rank indices')
      ENDIF

!     Map from modes to subset-of-modes, and subset dimensions
      ALLOCATE(subm(nsubm),resm(nsubm))
      ALLOCATE(rsub(nsubm),csub(nsubm),symsub(nsubm))
      d=1
      DO j=1,ndof
         IF (modes(j)) THEN
            subm(d)=j
            IF (allocated(M3%coef)) THEN
               resm(d)=j
            ELSE
               resm(d)=d
            ENDIF
            rsub(d)=M1%rows(j)
            csub(d)=M1%cols(j)
            symsub=M1%sym(j)
            d=d+1
         ENDIF
      ENDDO

      IF (.not.allocated(M3%coef)) M3=NewCP(rk3,rsub,csub,symsub)
      deallocate(rsub,csub,symsub)

!     Loop over modes
      DO d=1,nsubm
         j=subm(d)
         h=resm(d)

         IF (M1%nbas(j).ne.M2%nbas(j)) THEN
            write(*,*) 'Mode ',j,': M1 length (',M1%nbas(j),&
                       '), M2 length (',M2%nbas(j),') must match'
            call AbortWithError('Error in CPMM_vec()')
         ENDIF
         IF (M1%nbas(j).ne.M3%nbas(h)) THEN
            write(*,*) 'Mode ',j,': M1 length (',M1%nbas(j),&
                       '), M3 length (',M3%nbas(h),') must match'
            call AbortWithError('Error in CPMM_vec()')
         ENDIF

!        Loop over rank terms M1
         ik=ri3
         DO k=ri1,re1 ! k=1,rk1
!           Loop over rank terms in M2
            DO i=ri2,re2 ! i=1,rk2
               IF (d.eq.1) M3%coef(ik)=M1%coef(k)*M2%coef(i)
               IF (j.eq.1) THEN
                  afac=alpha1*alpha2
               ELSE
                  afac=1.d0
               ENDIF               
               M3%base(M3%ibas(h):M3%fbas(h),ik)=afac*&
               M1%base(M1%ibas(j):M1%fbas(j),k)*&
               M2%base(M2%ibas(j):M2%fbas(j),i)
               ik=ik+1
            ENDDO

!           Shift on M2
            IF (ishift2.ne.0) THEN
               IF (d.eq.1) M3%coef(ik)=M1%coef(k)*I2%coef(1)
               IF (j.eq.1) THEN
                  afac=alpha1
               ELSE
                  afac=1.d0
               ENDIF
               M3%base(M3%ibas(h):M3%fbas(h),ik)=afac*&
               M1%base(M1%ibas(j):M1%fbas(j),k)*&
               I2%base(I2%ibas(j):I2%fbas(j),1)
               ik=ik+1
            ENDIF        
         ENDDO

!        Shift on M1
         IF (ishift1.ne.0) THEN
            DO i=ri2,re2 ! i=1,rk2
               IF (d.eq.1) M3%coef(ik)=I1%coef(1)*M2%coef(i)
               IF (j.eq.1) THEN
                  afac=alpha2
               ELSE
                  afac=1.d0
               ENDIF
               M3%base(M3%ibas(h):M3%fbas(h),ik)=afac*&
               I1%base(I1%ibas(j):I1%fbas(j),1)*&
               M2%base(M2%ibas(j):M2%fbas(j),i)
               ik=ik+1
           ENDDO
         ENDIF

!        Shift on both M1 and M2
         IF (ishift1.ne.0 .and. ishift2.ne.0) THEN
            IF (d.eq.1) M3%coef(ik)=I1%coef(1)*I2%coef(1)
            M3%base(M3%ibas(h):M3%fbas(h),ik)=&
            I1%base(I1%ibas(j):I1%fbas(j),1)*&
            I2%base(I2%ibas(j):I2%fbas(j),1)
         ENDIF
      ENDDO ! loop over j

      deallocate(subm,resm)

      call CPU_TIME(ti2)
      mvp_time=mvp_time+ti2-ti1

      end subroutine CPMM_vec_general

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine CPMMreshift(M1,ishift1,Eold1,Enew1,&
                             M2,ishift2,Eold2,Enew2,M3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Computes new coefs of M3 resulting from changing only an Eshift
! Old and new Eshifts must be provided on input (to make sure signs did 
! not change)

      implicit none
      TYPE (CP), INTENT(IN)    :: M1,M2
      TYPE (CP), INTENT(INOUT) :: M3
      integer, intent(in) :: ishift1,ishift2
      real*8, intent(in)  :: Eold1,Enew1,Eold2,Enew2
      logical :: sc1,sc2
      integer :: i,k,ik,rk1,rk2,rk3
      real*8  :: ti1,ti2

      IF (.NOT. MVP_SETUP) call InitializePRODHVModule()

      call CPU_TIME(ti1)

!     Set parameters: # modes, ranks
      rk1=SIZE(M1%coef)
      rk2=SIZE(M2%coef)
      rk3=rk1*rk2

!     Check if sign changes
      sc1=(Eold1*Enew1).lt.0
      sc2=(Eold2*Enew2).lt.0
      
!     Modify ranks for shifts
      IF (ishift1.ne.0) THEN
         rk3=rk3+rk2
      ENDIF
      IF (ishift2.ne.0) THEN
         rk3=rk3+rk1
      ENDIF
      IF (ishift1.ne.0 .and. ishift2.ne.0) rk3=rk3+1

!     Error checking
      IF (SIZE(M3%coef).lt.rk3) THEN
         write(*,*) 'M3 has rank ',SIZE(M3%coef),&
                    ', but must be at least ',rk3
         call AbortWithError('CPMMreshift(): insufficient rk_M3')
      ENDIF

!     Shift on M2
      ik=1
      IF (ishift2.ne.0) THEN
         DO k=1,rk1
            ik=ik+rk2
            M3%coef(ik)=M1%coef(k)*abs(Enew2)
            IF (sc2) call VecSignChange(M3,ik,ik)
            ik=ik+1
         ENDDO
      ELSE
         ik=ik+rk1*rk2
      ENDIF

!     Shift on M1
      IF (ishift1.ne.0) THEN
         M3%coef(ik:ik+rk2-1)=abs(Enew1)*M2%coef(1:rk2)
         IF (sc1) call VecSignChange(M3,ik,ik+rk2-1)
         ik=ik+rk2

!        Shift on both M1 and M2
         IF (ishift2.ne.0) THEN
            M3%coef(ik)=abs(Enew1)*abs(Enew2)
            IF (sc1.neqv.sc2) call VecSignChange(M3,ik,ik)
         ENDIF
      ENDIF

      call CPU_TIME(ti2)
      mvp_time=mvp_time+ti2-ti1

      end subroutine CPMMreshift

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE CPMMM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
