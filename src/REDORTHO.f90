!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE REDORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates initial guess wavefunctions

      USE ERRORTRAP
      USE UTILS
      USE SEPDREPN
      USE CPCONFIG
      USE MODVECVEC

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_orthog(G,F,rF)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Projects G onto an orthogonal basis of unit vectors. The largest rF 
! terms are retained

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (CPvec), INTENT(OUT) :: F
      TYPE (Configs) :: v
      integer, intent(in)  :: rF

      call GetConfigList(G,rF,v)
      call Config2CP(F,v)
      call FlushConfigs(v)

      end subroutine reduc_orthog

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetConfigList(G,rF,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This wrapper which GetConfigList algorithm to use. 
! GetConfigList1() has quadratic complexity in rF, but is faster than
! GetConfigList1(), an O(rF*log_2(rF)) algorithm, if rF is <~100 or so.
! Due to the way the config list is updated, the two algorithms can
! give slightly different results.

      implicit none
      TYPE (CPvec), INTENT(IN)    :: G
      TYPE (Configs), INTENT(OUT) :: v
      integer, intent(in) :: rF
      integer, parameter  :: switch=100

      IF (rF.lt.switch) THEN
         call GetConfigList1(G,rF,v)
      ELSE
         call GetConfigList2(G,rF,v)
      ENDIF

      end subroutine GetConfigList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetConfigList1(G,rF,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Beginning with CP-vector G, this extracts the dominant configurations
! in terms of the product configurations of unit vectors. rF is the max
! number of configurations to keep

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (Configs), intent(out) :: v
      TYPE (Configs) :: w
      integer, intent(in)  :: rF
      real*8, allocatable  :: bas1D(:,:),tabindx(:,:)
      real*8  :: cuttol,biggest,mincoef,t1,t2
      integer :: rG,rrx,ndof,mxbas,i,j,k,l,n,imod,gst,minl,rterm,nfound
      logical :: found
      real*8, parameter  :: tol=1.d-15

      call CPU_TIME(t1)

!     Set parameters
      rG=SIZE(G%coef)
      ndof=SIZE(G%nbas)
      mxbas=MAXVAL(G%nbas)

!     Truncate rank if > direct-product size (unlikely if > 2 DOF)
      rrx=1
      DO j=1,ndof
         rrx=rrx*G%nbas(j)
         IF (rrx.ge.rF) EXIT
      ENDDO
      rrx=min(rrx,rF)

      call NewConfigs(v,G%nbas,rrx)

      ALLOCATE(bas1D(ndof,mxbas),tabindx(ndof,mxbas))
      tabindx=0.d0

!     Loop over terms in G
      DO i=1,rG

!        Get the list of 1D basis functions and sort by decreasing
!        coefficient magnitude
         bas1D=0.d0
         gst=0
         DO j=1,ndof
            n=G%nbas(j)
!           Copy coefficient into bas1D and set the table index
            DO k=1,n
               imod=gst+k
               bas1D(j,k)=abs(G%base(imod,i))
               tabindx(j,k)=SIGN(REAL(k),G%base(imod,i))
            ENDDO
            call dsort(bas1D(j,1:n),tabindx(j,1:n),G%nbas(j),-2)
            gst=gst+G%nbas(j)
         ENDDO

!        First iteration: set the biggest value
         IF (i.eq.1) THEN
            biggest=abs(G%coef(1))
            DO j=1,ndof
               biggest=biggest*bas1D(j,1)
            ENDDO
         ENDIF
         cuttol=tol*biggest/abs(G%coef(i))

!        Set the new cutoff tolerence and generate the configuration list
         call NewConfigs(w,G%nbas,rrx)

         IF (abs(G%coef(i)).le.tol*biggest/1.d300) cuttol=1.d300
         call sortDPbasis(w,bas1D,cuttol,rterm)

         biggest=max(biggest,abs(G%coef(i))*w%coef(1))

!        Replace the quantum numbers with those in tabindx
!        Scale the coefficients in w by the coefficient in G
         w%coef=G%coef(i)*w%coef
         DO k=1,rterm
            DO j=1,ndof
               IF (tabindx(j,w%qns(k,j)).lt.0.d0) w%coef(k)=-w%coef(k)
               w%qns(k,j)=abs(tabindx(j,w%qns(k,j)))
            ENDDO
         ENDDO

!        Update the total list of configurations
         IF (i.gt.1) THEN

!           Check each term in the new list against the total list
            DO k=1,rterm

               minl=0
               mincoef=9.d99
!               nfound=1
!               found=.FALSE.
               DO l=1,rrx
!                 Term matches: modify an existing coefficient
                  IF (ALL(w%qns(k,:).eq.v%qns(l,:))) THEN
                     v%coef(l)=v%coef(l)+w%coef(k)
!                     found=.TRUE.
                     EXIT

!                 First "empty" slot in v: add new term here
                  ELSEIF (v%qns(l,1).eq.0) THEN
                     call GenCopyConfigsWtoV(v,w,l,l,k,k)
                     EXIT

!                 Term does not match: keep track of which existing
!                 coef has the smallest magnitude...
                  ELSE
                     IF (abs(v%coef(l)).lt.mincoef) THEN
                        minl=l
                        mincoef=abs(v%coef(l))
                     ENDIF
                  ENDIF

!                 Term not found in v: if term is larger than smallest 
!                 term in v, replace the smallest term in v
                  IF (l.eq.rrx .and. abs(w%coef(k)).gt.mincoef) THEN
                     call GenCopyConfigsWtoV(v,w,minl,minl,k,k)
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            IF (rterm.gt.0) call GenCopyConfigsWtoV(v,w,1,rterm,1,rterm)
         ENDIF
         call FlushConfigs(w)
      ENDDO

!     Sort the configuration list and resize
      call SortConfigsByCoef(v)
      k=0
      DO i=1,rrx
         IF (v%qns(i,1).ne.0 .and. &
             abs(v%coef(i)).gt.abs(v%coef(1))*tol) k=k+1
      ENDDO
      call ResizeConfigList(v,k)

      DEALLOCATE(bas1D,tabindx)

      end subroutine GetConfigList1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GetConfigList2(G,rF,v)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Beginning with CP-vector G, this extracts the dominant configurations
! in terms of the product configurations of unit vectors. rF is the max
! number of configurations to keep per term

      implicit none
      TYPE (CPvec), INTENT(IN)  :: G
      TYPE (Configs), intent(out) :: v
      TYPE (Configs) :: w
      integer, intent(in)  :: rF
      real*8, allocatable  :: bas1D(:,:),tabindx(:,:)
      real*8  :: cuttol,biggest,t1,t2
      integer :: rG,rrx,ndof,mxbas,i,j,k,l,n,imod,gst,rterm
      real*8, parameter  :: tol=1.d-15

!     Set parameters
      rG=SIZE(G%coef)
      ndof=SIZE(G%nbas)
      mxbas=MAXVAL(G%nbas)

!     Truncate rank if > direct-product size (unlikely if > 2 DOF)
      rrx=1
      DO j=1,ndof
         rrx=rrx*G%nbas(j)
         IF (rrx.ge.rF) EXIT
      ENDDO
      rrx=min(rrx,rF)

      call NewConfigs(v,G%nbas,rrx)

      ALLOCATE(bas1D(ndof,mxbas),tabindx(ndof,mxbas))
      tabindx=0.d0

!     Loop over terms in G
      DO i=1,rG

!        Get the list of 1D basis functions and sort by decreasing
!        coefficient magnitude
         bas1D=0.d0
         gst=0
         DO j=1,ndof
            n=G%nbas(j)
!           Copy coefficient into bas1D and set the table index
            DO k=1,n
               imod=gst+k
               bas1D(j,k)=abs(G%base(imod,i))
               tabindx(j,k)=SIGN(REAL(k),G%base(imod,i))
            ENDDO
            call dsort(bas1D(j,1:n),tabindx(j,1:n),G%nbas(j),-2)
            gst=gst+G%nbas(j)
         ENDDO

!        First iteration: set the biggest value
         IF (i.eq.1) THEN
            biggest=abs(G%coef(1))
            DO j=1,ndof
               biggest=biggest*bas1D(j,1)
            ENDDO
         ENDIF
         cuttol=tol*biggest/abs(G%coef(i))

!        Set the new cutoff tolerence and generate the config list
         call NewConfigs(w,G%nbas,rrx)

         IF (abs(G%coef(i)).le.tol*biggest/1.d300) cuttol=1.d300
         call sortDPbasis(w,bas1D,cuttol,rterm)

         biggest=max(biggest,abs(G%coef(i))*w%coef(1))

!        Replace the quantum numbers with those in tabindx
!        Scale the coefficients in w by the coefficient in G
         w%coef=G%coef(i)*w%coef
         DO k=1,rterm
            DO j=1,ndof
               IF (tabindx(j,w%qns(k,j)).lt.0.d0) w%coef(k)=-w%coef(k)
               w%qns(k,j)=abs(tabindx(j,w%qns(k,j)))
            ENDDO
         ENDDO

!        Update the total list of configurations
         IF (i.gt.1) THEN
            DO k=1,rterm
!              See if config in w is already present in v
               l=findconfigindex(v,w%qns(k,:),(/1,rrx/))
!              Config in w is in v: add coef to v, erase coef in w
               IF (l.ne.0) THEN
                  v%coef(l)=v%coef(l)+w%coef(k)
                  w%coef(k)=0.d0
               ENDIF
            ENDDO

!           Sort the lists in v and w by magnitude. Coefs in w that are
!           larger than ones in v will replace ones in v
            call SortConfigsByCoef(v)
            call SortConfigsByCoef(w)
         ENDIF

         DO k=1,rterm
            l=rrx-k+1
!           Coef in w replaces smaller one in v
            IF (abs(w%coef(k)).gt.abs(v%coef(l))) THEN
               call GenCopyConfigsWtoV(v,w,l,l,k,k)
!           All remaining coefs in w are smaller than all coefs in v
            ELSE
               EXIT
            ENDIF
         ENDDO

         call FlushConfigs(w)

!        Sort the master list by index for the next term
         IF (i.lt.rG) call SortConfigsByIndex(v)
      ENDDO

!     Sort the configuration list and resize
      call SortConfigsByCoef(v)
      call TrimZeroConfigs(v,v%coef(1)*tol)

      DEALLOCATE(bas1D,tabindx)

      end subroutine GetConfigList2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine sortDPbasis(v,bas1D,cuttol,rterm)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sorts product basis in decreasing order of dominant coefficients

      implicit none
      TYPE (Configs), intent(inout) :: v
      real*8,  intent(in)  :: bas1D(:,:)
      integer, allocatable :: qnlist(:,:),tmpqn(:)
      real*8, allocatable  :: prods(:)
      integer, intent(out) :: rterm
      real*8, intent(in)   :: cuttol
      integer :: i,j,k,l,ndof,prodND,nguess,ii,jj
      integer :: minind(1),nallowed,nrk
      real*8  :: prodt,btmp
      logical :: allowed,add2end

      ndof=SIZE(v%nbas)
      nrk=SIZE(v%coef)

!     Error if more states are requested than allowed by the basis
      prodND=1
      DO j=1,ndof
         prodND=prodND*v%nbas(j)
         IF (prodND.ge.nrk) EXIT
      ENDDO
      IF (nrk.gt.prodND) THEN
         write(*,*) 'Error: nrk is too large for product basis size'
         call AbortWithError('Error in sortDPbasis()')
      ENDIF

      ALLOCATE(qnlist(ndof*nrk,ndof),prods(ndof*nrk),tmpqn(ndof))

!     Initialize the list of q.n. list and the prods with the ground
!     state
      nguess=1
      qnlist(1,:)=1
      prods(1)=1.d0
      DO j=1,ndof
         prods(1)=prods(1)*bas1D(j,1)
      ENDDO

!     Build the list of states
      rterm=0
      DO i=1,nrk

!        Get the state in the list with the largest coefficient
         btmp=MAXVAL(prods(1:nguess))
         IF (btmp.lt.cuttol) EXIT
         rterm=rterm+1
         v%coef(i)=btmp
         minind=MAXLOC(prods(1:nguess))
         v%qns(i,:)=qnlist(minind(1),:)

!        Find the allowed excitations out of the previously-found
!        largest coef states and add to qnlist
         add2end=.FALSE.
         nallowed=0
         DO j=1,ndof
            prodt=v%coef(i)
            tmpqn=v%qns(i,:)
            call ExcitepQN(j,bas1D(j,:),v%nbas(j),tmpqn,prodt)
            allowed=chkallowd(qnlist(1:nguess,:),tmpqn,minind(1),&
                    v%nbas)
            IF (allowed) THEN
               nallowed=nallowed+1
               IF (add2end) THEN ! Allowed new q.n.: add to end
                  nguess=nguess+1
                  qnlist(nguess,:)=tmpqn
                  prods(nguess)=prodt
               ELSE  ! First new allowed q.n.: replace old one
                  qnlist(minind(1),:)=tmpqn
                  prods(minind(1))=prodt
                  add2end=.TRUE.
               ENDIF
            ENDIF
         ENDDO
!        If none of the excitations are allowed, remove from the list
         IF (nallowed.eq.0) THEN
            qnlist(minind(1),:)=v%nbas(:)+1
            prods(minind(1))=0.d0
         ENDIF
      ENDDO

      DEALLOCATE(qnlist,prods,tmpqn)

      end subroutine sortDPbasis

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ExcitepQN(dofnr,bas,nbas,qnlist,ddprod)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Excites a quantum number in a multi-D list, and replaces the prod 
! resulting from the old value with the prod resulting from the new value

      implicit none
      integer, intent(in) :: dofnr,nbas
      integer, intent(inout) :: qnlist(:)
      real*8, intent(in) :: bas(:)
      real*8, intent(inout) :: ddprod
      real*8 :: mult,div

!     If qnlist already exceeds the basis, do not excite
      IF (qnlist(dofnr).lt.nbas) THEN ! excite and update prod
         div=bas(qnlist(dofnr))
         mult=bas(qnlist(dofnr)+1)
         ddprod=ddprod*mult/div
         qnlist(dofnr)=qnlist(dofnr)+1
      ELSE
         ddprod=0.d0
         qnlist(dofnr)=nbas+1
      ENDIF

      end subroutine ExcitepQN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function chkallowd(qnlist,tmpqn,iskip,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Return true if excitation of quantum number is allowed, i.e. if there
! are no states in the q.n. list with ALL q.n.s <= the ones in tmpqn

      implicit none
      integer, intent(in) :: qnlist(:,:),tmpqn(:),nbas(:)
      integer, intent(in) :: iskip
      integer :: k,l,ndof,nterms
      logical :: chkallowd

      ndof=SIZE(tmpqn)
      nterms=SIZE(qnlist,1)

      chkallowd=.TRUE.

      DO k=1,nterms
         IF (k.eq.iskip) CYCLE
         DO l=1,ndof
            IF (tmpqn(l).gt.nbas(l)) THEN
               chkallowd=.TRUE.
               EXIT
            ENDIF
            IF (tmpqn(l).lt.qnlist(k,l)) EXIT
            IF (l.eq.ndof) chkallowd=.FALSE.
         ENDDO
         IF (.not. chkallowd) EXIT
      ENDDO

      end function chkallowd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine reduc_bysorting(G,rrG)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reduces rank of G analytically by sorting and combining common factors
! The base of G should be normalized before passing to this subroutine
! After sorting, the rank is truncated at rrG if it exceeds this value

      implicit none
      TYPE (CPvec), INTENT(INOUT) :: G
      TYPE (CPvec) :: F
      integer, intent(in) :: rrG
      integer :: rF,rG,i,j,k,l,ndof,ndiff,idiff,gst,pos
      real*8  :: norm1D,cuttol
      logical :: found
      REAL*8, PARAMETER :: tol=1.d-12

      ndof=SIZE(G%nbas)
      rG=SIZE(G%coef)
      cuttol=MAXVAL(G%coef)*tol
      call NewCPvec(F,G%nbas,rG)

      DO
!        Copy the first term from G to F
         rF=1
         call GenCopyWtoV(F,G,rF,rF,rF,rF)

!        Now sweep through the terms in G and F and compare the bases
         DO i=2,rG
!           Filter out terms with small coefs
            IF (abs(G%coef(i)).lt.cuttol) CYCLE

            found=.FALSE.
            DO j=1,rF
!              Compare the base for each DOF for equality
               call comparebases(G,i,F,j,ndiff,idiff,pos)
               IF (ndiff.lt.2) THEN

!                 No DOFs different: just add coef
                  IF (ndiff.eq.0) THEN
                     F%coef(j)=F%coef(j)+pos*G%coef(i)
                     IF (F%coef(j).lt.0.d0) THEN
                        F%coef(j)=abs(F%coef(j))
                        call VecSignChange(F,j,j)
                     ENDIF

!                 One DOF different: sum different bases and renormalize
                  ELSEIF (ndiff.eq.1) THEN
                     gst=0
                     DO k=2,idiff
                        gst=gst+G%nbas(k-1)
                     ENDDO
                     norm1D=0.d0
                     DO k=1,G%nbas(idiff)
                        F%base(gst+k,j)=F%base(gst+k,j)+&
                        G%base(gst+k,i)*pos*G%coef(i)/F%coef(j)
                        norm1D=norm1D+F%base(gst+k,j)*F%base(gst+k,j)
                     ENDDO
                     norm1D=sqrt(abs(norm1D))

                     F%coef(j)=norm1D*F%coef(j)
                     DO k=1,G%nbas(idiff)
                        F%base(gst+k,j)=F%base(gst+k,j)/norm1D
                     ENDDO
                  ENDIF
                  found=.TRUE.
                  EXIT
               ENDIF
            ENDDO
!           If G isn't already in F, add it
            IF (.not.found) THEN
               rF=rF+1
               call GenCopyWtoV(F,G,rF,rF,i,i)
            ENDIF
         ENDDO

!        Cannot reduce G further: get out of here
         IF (rF.eq.rG) EXIT

         call GenCopyWtoV(G,F,1,rF,1,rF)
         rG=rF
      ENDDO

!     At this point the sorted terms are the first rF terms in F. Post-
!     process (normalize/sort/truncate) and transfer terms to G
      call FlushCPvec(G)
      call NewCPvec(G,F%nbas,rF)
      call GenCopyWtoV(G,F,1,rF,1,rF)
      call FlushCPvec(F)
      call NORMBASE(G)
      call ordre(G,F)
      call FlushCPvec(G)
      rF=min(rF,rrG)
      call NewCPvec(G,F%nbas,rF)
      call GenCopyWtoV(G,F,1,rF,1,rF)
      call FlushCPvec(F)

      end subroutine reduc_bysorting

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine comparebases(G,gr,F,fr,ndiff,idiff,pos)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compares the gr-th term of G and the fr-th term of F for equality

      implicit none
      TYPE (CPvec), INTENT(IN) :: G,F
      integer, intent(in)  :: gr,fr
      integer, intent(out) :: ndiff,idiff,pos
      integer :: ndof,gst,gi,k,l,bsign
      logical :: same
      REAL*8, PARAMETER :: tol=1.d-12
      ndof=SIZE(G%nbas)
      ndiff=0
      idiff=0
      pos=1
      gst=0

      DO k=1,ndof

         same=.TRUE.
         bsign=0
         DO l=1,G%nbas(k)
            gi=gst+l
!           Compare terms that are farther than 0.5*tol from zero
            IF (abs(G%base(gi,gr)).gt.0.5*tol .or. &
                abs(F%base(gi,fr)).gt.0.5*tol) THEN

!              Terms are same
               IF (abs(G%base(gi,gr)-F%base(gi,fr)).lt.tol) THEN
                  IF (bsign.eq.0) THEN
                     bsign=1
                  ELSEIF (bsign.eq.-1) THEN
                     same=.FALSE.
                     EXIT
                  ENDIF

!              Same with opposite sign
               ELSEIF (abs(G%base(gi,gr)+F%base(gi,fr)).lt.tol) THEN
                  IF (bsign.eq.0) THEN
                     bsign=-1
                  ELSEIF (bsign.eq.1) THEN
                     same=.FALSE.
                     EXIT
                  ENDIF

!              Terms are different
               ELSE
                  bsign=1
                  same=.FALSE.
                  EXIT
               ENDIF
            ENDIF
         ENDDO

         IF (.not.same) THEN
            ndiff=ndiff+1
            IF (ndiff.gt.1) EXIT
            idiff=k
         ELSE
            pos=pos*bsign
         ENDIF

         gst=gst+G%nbas(k)
      ENDDO

      end subroutine comparebases

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE REDORTHO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
