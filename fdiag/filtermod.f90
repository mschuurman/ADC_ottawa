module filtermod

  use constants
  
  implicit none

  save

  ! Energy window
  integer                               :: nener
  real(d), dimension(2)                 :: ebound
  real(d)                               :: de
  
  ! Window function
  integer                               :: iwfunc

  ! Autocorrelation functions
  integer                               :: ntauto,ntprop
  real(d)                               :: t0,dt,proptime,&
                                           autotime
  complex(d), dimension(:), allocatable :: auto,auto1,auto2

  ! Hamiltonian and overlap matrices
  integer                               :: nrbas
  real(d), dimension(:,:), allocatable  :: hfbas,hrbas,sfbas,&
                                           h2fbas,h2rbas,ubar,&
                                           transmat,normfac,&
                                           rvec
  real(d), dimension(:), allocatable    :: rener

  ! Intensities
  real(d), dimension(:), allocatable    :: dvec,avec,intens
  
  ! Unit conversion factors
  real(d), parameter                    :: fs2au=41.3413745758d0
  real(d), parameter                    :: au2fs=0.02418884254d0
  real(d), parameter                    :: ev2eh=0.0367493d0
  real(d), parameter                    :: eh2ev=27.2113845d0
  
end module filtermod
