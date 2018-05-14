  module hessmod
    
    use constants
    
    implicit none

    save

    integer                                        :: nlines,geomline,&
                                                      natm,ncoo,ijob,&
                                                      ngeom,nsta,nsta_i,&
                                                      nsta_f
    real(dp), dimension(:), allocatable            :: xcoo0,xcoo,ref
    real(dp), dimension(:,:), allocatable          :: pos,neg
    real(dp), dimension(:,:,:), allocatable        :: pospos,negneg,&
                                                      posneg,negpos
    real(dp), parameter                            :: dx=0.01d0
    character(len=2), dimension(:), allocatable    :: aatm
    character(len=120)                             :: infile,listfile
    character(len=120), dimension(:), allocatable  :: aline,adcout

  end module hessmod
