!######################################################################
! fdstates: Routines for the calculation of filter diagonalisation
!           states via the real-time propagation of the f-vector
!           f_J = < Psi_J | D | Psi_0 >
!######################################################################
module fdstates

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  implicit none

  integer                               :: matdim
  complex(d), dimension(:), allocatable :: psi0
  
contains

!######################################################################

  subroutine calc_fdstates(fvec,ndimf)

    implicit none

    integer, intent(in)                   :: ndimf
    integer                               :: k
    real(d), dimension(ndimf), intent(in) :: fvec
    real(d)                               :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Output where we are at and what we are doing
!----------------------------------------------------------------------
    call wrinfo

!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine calc_fdstates

!######################################################################

  subroutine wrinfo

    implicit none

    integer :: k

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of filter diagonalisation &
         eigenstates'
    write(ilog,'(70a,/)') ('-',k=1,70)

!----------------------------------------------------------------------
! SIL parameters
!----------------------------------------------------------------------         
    write(ilog,'(2x,a,/)') 'Wavepacket propagation performed using &
         the SIL method'
    write(ilog,'(2x,a,x,i2,/)') 'Maximum Krylov subspace dimension:',&
         kdim
    write(ilog,'(2x,a,x,ES15.8,/)') 'Error tolerance:',autotol

    return

  end subroutine wrinfo
  
!######################################################################
  
end module fdstates
