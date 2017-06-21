  module mcspmod

    save

    integer                                      :: nfiles,maxdat,&
                                                    maxintvl,npoints
    integer, dimension(:), allocatable           :: ndat,nintvl
    real*8, dimension(:,:), allocatable          :: dat,x,deriv,s,dx
    real*8, dimension(:,:,:), allocatable        :: coeff
    real*8, dimension(2)                         :: erange
    character(len=80), dimension(:), allocatable :: datfile

  end module mcspmod
