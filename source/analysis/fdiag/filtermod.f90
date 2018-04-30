module filtermod

  use constants
  
  implicit none

  save

  ! Energy window
  integer                                :: nener
  real(dp), dimension(2)                 :: ebound
  real(dp)                               :: de
  
  ! Window function
  integer                                :: iwfunc

  ! Autocorrelation functions
  integer                                :: ntauto,ntprop
  real(dp)                               :: t0,dt,proptime,&
                                            autotime,tcutoff
  complex(dp), dimension(:), allocatable :: auto,auto1,auto2

  ! Hamiltonian and overlap matrices
  integer                                :: nrbas
  real(dp), dimension(:,:), allocatable  :: hfbas,hrbas,sfbas,&
                                            h2fbas,h2rbas,&
                                            transmat,normfac,&
                                            rvec
  real(dp), dimension(:), allocatable    :: rener
  real(dp)                               :: ovrthrsh
  
  ! Intensities
  real(dp), dimension(:), allocatable    :: dvec,avec,intens

  ! Error estimates
  real(dp), dimension(:), allocatable    :: error
  real(dp)                               :: errthrsh
  
  ! Unit conversion factors
  real(dp), parameter                    :: fs2au=41.3413745758d0
  real(dp), parameter                    :: au2fs=0.02418884254d0
  real(dp), parameter                    :: ev2eh=0.0367493d0
  real(dp), parameter                    :: eh2ev=27.2113845d0

  ! Output
  integer                                :: idat
  character(len=120)                     :: adat
  
end module filtermod
