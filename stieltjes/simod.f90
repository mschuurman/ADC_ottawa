  module simod

    save

    integer                              :: npoints,nord
    integer, parameter                   :: ntrial=100
    real*8, dimension(2)                 :: erange
    real*8, dimension(:), allocatable    :: ener,osc
    real*16, dimension(:), allocatable   :: e,f,a,b,si_e,si_F,si_osc
    real*16, dimension(:,:), allocatable :: Q
    character(len=120)                   :: asiinp,aosc

  end module simod
