!#######################################################################
! propagate_adc1: Routines to perform ADC(1) wavepacket propagations
!                 including the interaction of the molecule with a
!                 laser pulse
!#######################################################################

module propagate_adc1

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                               :: matdim
  complex(d), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser_adc1: ADC(1) wavepacket propagation including the
!                       molecule-laser interaction.
!                       The wavepacket is represented in a basis
!                       consisting of the intermediate state basis plus
!                       the HF ground state.
!########################################################################
! IMPORTANT:            For ease of implementation of the Hamiltonian
!                       and dipole matrix-vector products, the HF ground
!                       state is taken to be the last basis function in
!                       the set.
!#######################################################################
  
  subroutine propagate_laser_adc1(ndimf,kpqf)

    use tdsemod
    
    implicit none

    integer, intent(in)                       :: ndimf
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(d)                                   :: tw1,tw2,tc1,tc2

    print*,"HERE"
    
    return
    
  end subroutine propagate_laser_adc1

!#######################################################################
  
end module propagate_adc1
