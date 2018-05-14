  module rixsmod
    
    use constants

    implicit none

    save

    ! Spectrum parameters
    integer                :: nener1,nener2
    real(dp), dimension(2) :: ener1,ener2
    real(dp)               :: de1,de2
    real(dp)               :: gammaint,gammaf
    real(dp)               :: theta
    character(len=70)      :: datfile

    ! Transition matrix elements
    integer                                 :: istate,nval
    integer                                 :: nlanc
    real(dp), dimension(:), allocatable     :: enerval,enerlanc
    real(dp), dimension(:,:), allocatable   :: transsq
    real(dp), dimension(:,:,:), allocatable :: tdm,zeta

  end module rixsmod

