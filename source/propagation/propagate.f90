!#######################################################################
! molecule_laser: Routines to perform wavepacket propagations
!                 including the interaction of the molecule with a
!                 laser pulse
!#######################################################################

module propagate

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                               :: wpdim
  complex(d), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser: Wavepacket propagation including the molecule-laser
!                  interaction.
!                  The wavepacket is represented in a basis consisting
!                  of the intermediate state basis plus the MP2 ground
!                  state.
!########################################################################
! IMPORTANT:       For ease of implementation of the Hamiltonian and
!                  dipole matrix-vector products, the MP2 ground state
!                  is taken to be the last basis function in the set.
!#######################################################################
  
  subroutine propagate_laser(ndimf,noffdf)

    use tdsemod
    
    implicit none

    integer, intent(in)   :: ndimf
    integer*8, intent(in) :: noffdf
    integer               :: k
    real(d)               :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Output some information about the calculation to be performed
!----------------------------------------------------------------------
    call wrinfo

!----------------------------------------------------------------------
! Set up the initial wavefunction vector
!----------------------------------------------------------------------
    call init_wavefunction(ndimf)

!----------------------------------------------------------------------
! Peform the wavepacket propagation
!----------------------------------------------------------------------
    errmsg='Finish writing the wavepacket propagation code!'
    call error_control
    
!----------------------------------------------------------------------
! Finalise and deallocate arrays
!----------------------------------------------------------------------
    call finalise
    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine propagate_laser

!#######################################################################

  subroutine wrinfo

    implicit none

    integer :: k

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Wavepacket propagation'
    write(ilog,'(72a)') ('-',k=1,72)

!----------------------------------------------------------------------
! SIL parameters
!----------------------------------------------------------------------         
    if (lcap) then
       write(ilog,'(2x,a,/)') 'Wavepacket propagation performed &
            using the CSIL method'
    else
       write(ilog,'(2x,a,/)') 'Wavepacket propagation performed &
            using the SIL method'
    endif
       
    write(ilog,'(2x,a,x,i2,/)') 'Maximum Krylov subspace dimension:',&
         kdim

!    write(ilog,'(2x,a,x,ES15.8)') 'Error tolerance:',autotol

!!----------------------------------------------------------------------
!! Matrix-vector multiplication algorithm
!!----------------------------------------------------------------------
!    if (hincore) then
!       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
!            will proceed in-core'
!    else
!       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
!            will proceed out-of-core'
!    endif
    
    return
    
  end subroutine wrinfo

!#######################################################################

  subroutine init_wavefunction(ndimf)

    implicit none

    integer :: ndimf
    
!----------------------------------------------------------------------
! Wavepacket dimension: IS basis plus the ground state
!----------------------------------------------------------------------
    wpdim=ndimf+1

!----------------------------------------------------------------------
! Allocate the wavepacket array
!----------------------------------------------------------------------
    allocate(psi(wpdim))
    psi=0.0d0
    
!----------------------------------------------------------------------
! Set the initial wavepacket
!
! For now, we will only support excitation/ionisation from the ground
! state, corresponding to the vector (0,...,0,1)^T
!----------------------------------------------------------------------
    psi=0.0d0
    psi(wpdim)=cone
    
    return
    
  end subroutine init_wavefunction

!#######################################################################

  subroutine finalise

    implicit none

    deallocate(psi)
    
    return
    
  end subroutine finalise
    
!#######################################################################
  
end module propagate
