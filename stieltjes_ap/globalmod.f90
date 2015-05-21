  module globalmod

    implicit none

    integer, parameter                 :: maxbas=2048
    integer                            :: nbas,nocc,adc2
    integer, dimension(8,8)            :: mt
    integer, dimension(maxbas)         :: sym,epsilon,occ
    real*8, dimension(maxbas,maxbas)   :: fmatrix
    real*8, dimension(maxbas,maxbas,2) :: d_mo
    real*8                             :: eps_min,eps_max,coup_threshold
    character(len=80)                  :: aout

  end module globalmod
