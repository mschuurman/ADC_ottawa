  module rixsmod
    
    use constants

    implicit none

    save

    ! Spectrum parameters
    integer                              :: nener1,nener2
    real(d), dimension(2)                :: ener1,ener2
    real(d)                              :: de1,de2
    real(d)                              :: gammaint,gammaf
    character(len=70)                    :: datfile

    ! Transition matrix elements
    integer                              :: istate,nlanc,nval
    real(d), dimension(:), allocatable   :: enerval,enerlanc
    real(d), dimension(:,:), allocatable :: tdm,transsq

  end module rixsmod

