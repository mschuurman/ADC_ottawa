  module simod

    save

    integer                           :: npoints,nord,ntrial
    real*8, dimension(2)              :: erange
    real*8, dimension(:), allocatable :: ener,osc
    character(len=120)                :: asiinp,aosc

  end module simod
