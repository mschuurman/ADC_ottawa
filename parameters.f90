module parameters

  use constants
  use integrals_mo2e

!!$*************************************
!!$*******MOLCAS provided variables*****
!!$************************************* 
!!$int*4 nIrr - number of irreps
!!$int*4 nBas - number of basis functions
!!$int*4 nCen - number of atomic centres
!!$int*4 nOcc - number of occupied orbitals
!!$int*4 nVirt - number of virtual orbitals
!!$int*4 Ehf - Hartree-Fock energy of the ground state
!!$orbSym - int*4 array(nBas) containing irrep labels of MO's
!!$e - real*8 array(nBas) containing MO's energies
!!$occNum - real*8 array(nBas) containing MO's occupation numbers 
!!$x,y,z-dipole real*8 array(nBas,nBas) containing dipole moment matrix elements 
  
  integer*4                            :: nelec,nIrr,nBas,nCen
  integer                              :: nOcc,nVirt
  real(d)                              :: Ehf
  real(d)                              :: E_MP2
  integer*4, dimension(:), allocatable :: orbSym
  character*2,dimension(3)             :: labSym
  real(d), dimension(:), allocatable   :: e,occNum
  real(d), dimension(:,:), allocatable :: x_dipole,y_dipole,z_dipole,dpl
  real(d), dimension(:,:), allocatable :: density_matrix
  type(moint2e_cache)                  :: moIntegrals  ! Currently active MO integrals context

!!$**************************************************
!!$*********User provided variables******************
!!$**************************************************

!!$logical debug - debug levels of output
!!$int hinit - initial hole orbital number
!!$int nirrep - number of irrep of the initial excitation
!!$int idiag - diagonalisation proc. for the init. state:1-Lapack,2-Davidson; activated if chrun is not 'direct'
!!$int fdiag - diagonalisation proc. for the final state:1-Lapack,2-Lanczos; activated if chrun is not 'direct' 
!!$int method - comp. method: (0) - fano, (1) - tda, (2) - adc2 (3) - adc2e
!!$int fmethod - fano method: (1) - tda, (2) - adc2 (3) - adc2e
!!$int davstates - number of requested Davidson initial states 
!!$int lancstates - number of lanczos computed final states; Number of printed states in simple ADC2 calculations
!!$int numinista - number of initial states in Fano calculations;
!!$int array(1:numinista) ninista - the initial states for which appropriate Lanczos generated final state manifold 
!!$is to be computed  
!!$int stiprilev - print level in the Stieltjes subroutine
!!$int norder - Stieltjes order
!!$real*8 minc - minimum value of accepted adc matr. elem.
!!$real*8 mspacewi - minimum weight of single exc. in an allowed initial state
!!$real*8 mspacewf - minimum weight of single exc. in an allowed final state
!!$logical readband - a flag requesting the final eigenvectors in an energy band (eupper,elower)
!!$real*8 eupper,elower - energy band borders 
!!$roccnum - int array(nBas) rearranged MO's in the order occ -> virt.
!!$int array(8,8) MT -  multiplication table for the abelian groups d2h, etc.
!!$int array(nhcentre) hcentre - contains labels of the occupied orbitals on an atom carrying the initial hole
!!$char(4) chrun - 'direct'-Lapack diagonalization of both matrices, 'save'-save matrices 
!!$needed for either Lanczos or Davidson routine, 'read' read vectors produced by the external 
!!$diag. routine and get lifetimes 
!!$char(36) davname,lancname  - the names of the davidson and lanczos vector files
!!$char(1) tranmom  - dipole moment index 'x','y', or 'z' corresp. to the nirrep.
!!$integer array(1:lmain) stvc_lbl  - Damit wirder der Lanc-Startblock festgelegt
!!$integer info if 1 stops execution after printing the configuration tables
!!$integer ninista gives the number of the fanostate among davidson eigenvectors

!!$logical ladc1guess - greater than zero if the initial vectors for the Davidson diagonalisation 
!!                      are to be generated from an ADC(1) calculation
!!$real*8 davtol - error tolerance for the block-Davidson eigensolver, default of 10^-8
!!$logical lcvs - flag to switch on the CVS approximation
!!$integer array(1:nhcentre) icore - array of indices indexing the core orbitals
!!$integer ncore - no. core orbitals
!!$logical lfakeip - true if we are going to include a basis function
!!                      with an extermely small exponent to mimic ionization
!!$integer ifakeorb - index of the 'continuum' orbital for a fake ip calculation
!!$ifakeex array(1:dmain) - array of of indices indexing the 1h1p
!!                          configs to be taken as guesses for the
!!                          Davidson diagonalisation in the case of a fake
!!                          IP calculation
!!$expfakeip - value of the exponent of the diffuse function used in a
!!             fake IP calculation

  logical                            :: debug
  character(1)                       :: tranmom,tranflag,tranmom2
  character(4)                       :: chrun
  character(4)                       :: chrun2
  character(4)                       :: WHAT
  integer                            :: matvec
  character(36)                      :: lancname,davname
  integer                            :: hinit,nirrep,nirrep2,idiag,fdiag,method,&
                                        fmethod,davstates,lancstates,stiprilev,&
                                        numinista,norder,info,statenumber,denord
  integer, parameter                 :: nhcentre=40
  integer, dimension(0:nhcentre)     ::  hcentre
  real(d)                            :: minc,mspacewi,mspacewf,eupper,elower
  real(d)                            :: dlim
  integer, dimension(:), allocatable :: roccnum
  integer                            :: mgvdim
  real(d), dimension(:), allocatable :: mgvec
  integer, dimension(8,8)            :: MT
  logical                            :: readband
  integer, dimension(400)            :: stvc_lbl
  integer                            :: ninista
  logical                            :: ladc1guess
  real(d)                            :: davtol
  logical                            :: lcvs
  integer, dimension(nhcentre)       :: icore
  integer                            :: ncore
  logical                            :: lfakeip
  integer                            :: ifakeorb
  integer, dimension(:), allocatable :: ifakeex
  real(d)                            :: expfakeip

!!$************************************************
!!$**********Physical Cobnstants*******************
!!$************************************************

!!$ abohr - Bohr radius [cm]
!!$ fsconstinv - inverted fine structure constant
!!$ os2cs - oscillator strength [a.u.] to cross-section [Mb] conversion factor
!!$ omega - photon energy, required by Stieltjes_phi, photoionisation routine.

!  real(d), parameter :: abohr=5.2918e-9
  real(d), parameter :: fsconstinv=137._d
  real(d), parameter :: os2cs=4.0347443
  real(d), parameter :: omega=3.0_d
  
!!$************************************************
!!$**********Lanczos Parameters********************
!!$************************************************  

  character*4           :: mtxidl
  integer               :: ncycles,maxmem,memx,lmain
  integer               :: mode,nprint
  integer               :: maxiter
  integer, dimension(2) :: iparm
  real(d)               :: wthr
  real(d), dimension(2) :: erange
  real(d)               :: unit
  real(d),dimension(5)  :: fparm

!!$************************************************
!!$**********Davidson Parameters********************
!!$************************************************  

  character*4 :: mtxidd
  integer     :: dmain
  logical     :: myb0,transp
  character*2 :: GO
  character*3 :: POLARIZATION

INTEGER                            :: NSYMA
INTEGER                            :: NSYMA_PROP
INTEGER                            :: HAM_PIECES
INTEGER                            :: KLPDTOT
integer, dimension(:), allocatable :: SYM_MAP
integer, dimension(:), allocatable :: DIM_PROP
integer, dimension(:), allocatable :: NDIV
integer, dimension(:), allocatable :: NDIV_DIP
integer, dimension(:), allocatable :: NBUF_SYM
integer, dimension(:), allocatable :: SIMMETRIA
integer, dimension(:), allocatable :: NLPD
integer, dimension(:), allocatable :: LPSYM
integer, dimension(:), allocatable :: MU1D
integer, dimension(:), allocatable :: MU2D
integer, dimension(:), allocatable :: KAPPAD
integer, dimension(:), allocatable :: DIPOLE_BLOCK
integer, dimension(:), allocatable :: NUM
integer, dimension(:), allocatable :: NUM_DIAG
integer, dimension(:), allocatable :: NREC_VECTOR
integer, dimension(:), allocatable :: NRECTOT_VECT
integer, dimension(:), allocatable :: NREC_VECTOR_BIS
integer, dimension(:), allocatable :: RECINI_VECT
integer, dimension(3)              :: DIPOLESYM
REAL*8 , dimension(3)              :: ELECTRIC_FIELD
integer                            :: CHECK_dip

end module parameters
