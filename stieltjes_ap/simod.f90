  module simod

    save

    integer                           :: npoints,nord,ntrial,method
    real*8, dimension(2)              :: erange
    real*8, dimension(:), allocatable :: ener,osc
    real*8, dimension(:), allocatable :: si_e,si_f,si_osc,si_cosc1,&
                                         si_cosc2
    real*8                            :: ade
    character(len=120)                :: asiinp,aosc

  end module simod
