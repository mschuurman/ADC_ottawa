!######################################################################
! density_matrix: routines for the calculation of 1-electron reduced
!                 density and transition density matrices
!######################################################################

module density_matrix

  use constants
  use parameters
  use iomod
  use vpqrsmod
  
  implicit none

  contains

!######################################################################
! adc2_trden_gs: calculation of the ADC(2) ground-to-excited state
!                transition density matrix
!######################################################################
    
    subroutine adc2_trden_gs(trdens,ndimf,kpqf,rvec,nstates)

      use mp2
      
      implicit none

      integer                                   :: i
      integer                                   :: ndimf,nstates
      integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
      real(d), dimension(nbas,nbas,nstates)     :: trdens
      real(d), dimension(ndimf,nstates)         :: rvec
      real(d), allocatable                      :: rhomp2(:,:)
      
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
      trdens=0.0d0

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
      allocate(rhomp2(nbas,nbas))
      rhomp2=0.0d0
      
!----------------------------------------------------------------------
! Precalculation of the MP2 correction to the ground state density
! matrix
!----------------------------------------------------------------------
      ! MP2 ground state density matrix
      call rho_mp2(rhomp2)

      ! Subtraction of the zeroth-order contribution to obtain the
      ! MP2 correction
      do i=1,nocc
         rhomp2(i,i)=rhomp2(i,i)-2.0d0
      enddo

!----------------------------------------------------------------------
! Calculation of the occupied-occupied block
!----------------------------------------------------------------------
      call adc2_trdens_gs_occ_occ(trdens,ndimf,kpqf,rvec,nstates,&
           rhomp2)
      
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
      deallocate(rhomp2)
      
      return
      
    end subroutine adc2_trden_gs

!######################################################################
! adc2_trden_gs_occ_occ: calculation of the occupied-occupied block of
!                        the ADC(2) ground-to-excited state transition
!                        density matrix
!######################################################################
    
    subroutine adc2_trdens_gs_occ_occ(trdens,ndimf,kpqf,rvec,nstates,&
         rhomp2)

      implicit none

      integer                                   :: ndimf,nstates
      integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
      integer                                   :: cnt,nlim1,nlim2
      integer                                   :: i,j,a
      real(d), dimension(nbas,nbas,nstates)     :: trdens
      real(d), dimension(ndimf,nstates)         :: rvec
      real(d), dimension(nbas,nbas)             :: rhomp2

!----------------------------------------------------------------------
! 1h1p contributions
!----------------------------------------------------------------------
      nlim1=1
      nlim2=kpqf(1,0)

      do cnt=nlim1,nlim2
         j=kpqf(3,cnt)
         a=kpqf(5,cnt)
         do i=1,nocc
            trdens(i,j,:)=trdens(i,j,:)+rhomp2(i,a)*rvec(cnt,:)
         enddo
      enddo

      trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:)*2.0d0*sqrt(2.0d0)
      
!----------------------------------------------------------------------
! 2h2p a=b, i=j contributions
!----------------------------------------------------------------------
      nlim1=nlim2+1
      nlim2=nlim2+kpqf(2,0)

      do cnt=nlim1,nlim2
         
      enddo
         
      STOP
      
      return
      
    end subroutine adc2_trdens_gs_occ_occ
      
!######################################################################
    
end module density_matrix
