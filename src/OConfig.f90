!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      MODULE OCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module generates 2-orthogonal configurations

      USE ERRORTRAP
      USE UTILS
      USE CPCONFIG

      CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Gen2OrthogConfigList(v,nrk,evals1D,nbas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Generates 2-orthogonal configurations

      implicit none
      TYPE (Configs), INTENT(OUT) :: v
      TYPE (Configs) :: wacc,wadd
      real(kind=8), intent(in) :: evals1D(:,:)
      real(kind=8), parameter  :: smallnr=1.d-14
      integer, intent(in) :: nbas(:)
      integer, intent(in) :: nrk
      integer :: ndof,maxbas,irk,ark,i,j,k,nadd,ndiff

      ndof=SIZE(evals1D,1)
      maxbas=SIZE(evals1D,2)
      nadd=ndof*(ndof-1)/2

!     Initial config
      call NewConfigs(v,nbas,nrk)
      call NewConfigs(wacc,nbas,1)
      wacc%qns(1,:)=1
      wacc%coef(1)=product(evals1D(:,1))

      write(*,*) 'The weights'
      call PrintMatrix(evals1D)

      DO irk=1,nrk
!       Add the config at the top of the list to v
        call GenCopyConfigsWtoV(v,wacc,irk,irk,1,1)

        IF (irk.eq.nrk) EXIT

        write(*,*) 'Iteration ',irk,': candidate config list'
        call PrintConfigs(wacc)

!       Generate new configs that are orthogonal to the top one
        call NewConfigs(wadd,nbas,nadd)
        k=0
        DO i=1,ndof-1
           IF (wadd%qns(i,1).eq.nbas(i)) CYCLE
           DO j=i+1,ndof
              IF (wadd%qns(j,1).eq.nbas(j)) CYCLE
              k=k+1

              call GenCopyConfigsWtoV(wadd,wacc,k,k,1,1)
              wadd%coef(k)=wadd%coef(k)/evals1D(i,wadd%qns(k,i))
              wadd%coef(k)=wadd%coef(k)/evals1D(j,wadd%qns(k,j))
              wadd%qns(k,i)=wadd%qns(k,i)+1
              wadd%qns(k,j)=wadd%qns(k,j)+1
              wadd%coef(k)=wadd%coef(k)*evals1D(i,wadd%qns(k,i))
              wadd%coef(k)=wadd%coef(k)*evals1D(j,wadd%qns(k,j))
           ENDDO
        ENDDO

        write(*,*) 'New configs added to candidate list:'
        call PrintConfigs(wadd)

!       Zero out configs that are less than 2-orthogonal to selected one
        ark=SIZE(wacc%coef)
        DO i=1,ark
           ndiff=0
           DO j=1,ndof
              IF (wacc%qns(i,j).ne.wacc%qns(1,j)) ndiff=ndiff+1
              IF (ndiff.eq.2) EXIT
           ENDDO
           IF (ndiff.lt.2) wacc%coef(i)=0.d0
        ENDDO

!       Combine the new and old lists, resort, and remove any duplicates
        call ResizeConfigList(wacc,ark+k)
        call GenCopyConfigsWtoV(wacc,wadd,ark+1,ark+k,1,k)
        call removeduplicateconfigs(wacc,0)
        call SortConfigsByCoef(wacc)
        call TrimZeroConfigs(wacc,smallnr)
        call FlushConfigs(wadd)
      ENDDO

      call FlushConfigs(wacc)
      write(*,*) 'Final selected config list:'
      call PrintConfigs(v)

      end subroutine Gen2OrthogConfigList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE OCONFIG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

