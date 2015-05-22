  module simod

    save

    integer                              :: npoints,nord
    integer, parameter                   :: ntrial=100
    real*8, dimension(2)                 :: erange
    real*8, dimension(:), allocatable    :: ener,osc
    real*16, dimension(:), allocatable   :: si_e,si_f,si_osc,si_cosc1,&
                                            si_cosc2
    character(len=120)                   :: asiinp,aosc

  end module simod
