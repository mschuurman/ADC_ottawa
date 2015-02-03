  module stieltjesmod

    implicit none

    save

    integer, parameter                 :: maxbas=2048
    integer                            :: mt(8,8)
    integer                            :: nocc,adc2
    real*8, dimension(maxbas)          :: epsilon,sym,occ
    real*8                             :: nbas
    real*8                             :: eps_min,eps_max,coup_threshold
    real*8, dimension(maxbas,maxbas,2) :: d_mo
    real*8, dimension(maxbas,maxbas)   :: fmatrix

  end module stieltjesmod
