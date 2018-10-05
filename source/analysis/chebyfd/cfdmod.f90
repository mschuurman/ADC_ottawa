module cfdmod

  use constants
  
  implicit none
  
  save

  ! Order domain autocorrelation function
  integer                :: order,kfinal,Kdim
  real(dp), allocatable  :: auto(:)

  ! Spectral bounds
  real(dp), dimension(2)  :: bounds

  ! Filter functions
  integer :: ifilter
  integer :: nfsbas

  !----------------------------------------------------------
  ! DPSS information
  !
  ! No. quadrature points
  integer, parameter  :: npts=5001
  !
  ! NxW
  real(dp)            :: fw
  !
  ! DPSSs
  real(dp), dimension(:,:), allocatable :: v
  !
  ! Eigenvalues
  real(dp), dimension(:), allocatable :: lambda
  !
  ! Tabulated optimal time-bandwidth products
  integer, parameter        :: nopt=50
  real(dp), dimension(nopt) :: optfw
  !
  ! Variable time-bandwidth products
  integer, parameter :: nprecalc=50
  logical            :: varfw
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Gaussian information
  !
  ! Width
  real(dp) :: sigma
  !
  ! Scaled width
  real(dp) :: sigmabar
  !----------------------------------------------------------
  
  ! Energy interval
  real(dp) :: Ea,Eb
  
  ! Scaled energy interval
  real(dp) :: Eabar,Ebbar

  ! Expansion coefficients
  real(dp), dimension(:,:), allocatable :: fkn

  ! Filtered state overlap and Hamiltonian matrices
  real(dp), dimension(:,:), allocatable :: smat,hmat

  ! Eigenvectors and eigenvalues of the filtered-state
  ! Hamiltonian matrix
  real(dp), dimension(:,:), allocatable :: eigvec
  real(dp), dimension(:), allocatable   :: eigval
  real(dp), dimension(:), allocatable   :: ener
  
  ! Projector  onto the orthogonal complement of the
  ! null space
  real(dp), dimension(:,:), allocatable :: transmat

  ! Dimension of the orthogonal complement of the null space
  integer :: nrbas

  ! Transition dipoles and oscillator strengths
  real(dp), dimension(:), allocatable :: tdm,osc
  
  ! Output
  integer            :: idat
  character(len=120) :: adat

  ! Unit conversion factors
  real(dp), parameter :: eh2ev=27.2113845d0
  real(dp)            :: convfac
  logical             :: lau
  
end module cfdmod
