module filtermod

  use constants
  
  implicit none

  save

  ! Energy window
  real(d), dimension(2)                 :: ebound
  integer                               :: nener

  ! Window function
  integer                               :: iwfunc

  ! Autocorrelation functions
  integer                               :: nt
  real(d)                               :: t0,dt
  complex(d), dimension(:), allocatable :: auto,auto1,auto2

  ! Unit conversion factors
  real(d), parameter                    :: fs2au=41.3413745758d0
  real(d), parameter                    :: au2fs=0.02418884254d0
  real(d), parameter                    :: ev2eh=0.0367493d0
  real(d), parameter                    :: eh2ev=27.2113845d0
  
end module filtermod
