  module rixsmod
    
    use constants

    implicit none

    save

    ! Spectrum parameters
    integer                              :: nener1,nener2
    real(d), dimension(2)                :: ener1,ener2
    real(d)                              :: de1,de2
    real(d)                              :: gammaint,gammaf
    character(len=70), dimension(3)      :: datfile

    ! Transition matrix elements
    integer                                :: istate,nval,maxlanc
    integer, dimension(3)                  :: nlanc
    real(d), dimension(:), allocatable     :: enerval
    real(d), dimension(:,:), allocatable   :: enerlanc,transsq
    real(d), dimension(:,:,:), allocatable :: tdm

  end module rixsmod

