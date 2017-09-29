module parameters

  use constants
  use integrals_mo2e

!!$*************************************
!!$*******Hartree-Fock MO arrays********
!!$************************************* 
!!$int*4 nIrr - number of irreps
!!$int*4 nBas - number of basis functions
!!$int*4 nCen - number of atomic centres
!!$int*4 nOcc - number of occupied orbitals
!!$int*4 nVirt - number of virtual orbitals
!!$real*8 Ehf - Hartree-Fock energy of the ground state
!!$orbSym - int*4 array(nBas) containing irrep labels of MO's
!!$e - real*8 array(nBas) containing MO's energies
!!$occNum - real*8 array(nBas) containing MO's occupation numbers 
!!$x,y,z-dipole real*8 array(nBas,nBas) containing dipole moment matrix elements
!!$hcore - real*8 array(nOcc,nOcc) containing the elements h_ij of the
!!         core Hamiltonian in the MO basis
  integer*4                            :: nelec,nIrr,nBas,nCen
  integer                              :: nOcc,nVirt
  real(d)                              :: Ehf
  real(d)                              :: E_MP2
  integer*4, dimension(:), allocatable :: orbSym
  character*3,dimension(8)             :: labSym
  real(d), dimension(:), allocatable   :: e,occNum
  real(d), dimension(:,:), allocatable :: x_dipole,y_dipole,z_dipole,dpl
  real(d), dimension(:,:), allocatable :: density_matrix
  type(moint2e_cache)                  :: moIntegrals  ! Currently active MO integrals context
  character*100                        :: moType ! Either 'incore' or 'disk'
  integer                              :: imotype
  real(d), dimension(:,:), allocatable :: hcore
  
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
!!
!!$logical ld2 - flag to turn on the calculation of the MP2 D2 diagnostic
!!
!!$logical lscf - flag used to determine whether to perform an SCF calculation
!!                to ensure that the MOs as read from the GAMESS output are converged
!!
!!$logical ltdm_gs2i - flag used to determine whether to calculate the
!!                     TDM between the ground state and the initial
!!                     space states
!!
!!$logical lifrzcore - flag used to determine whether to freeze the
!!                     core orbitals in the calculation of the initial
!!                     space states
!!
!!$logical lffrzcore - flag used to determine whether to freeze the
!!                     core orbitals in the calculation of the final
!!                     space states
!!
!!$logical ldiagfinal - flag used to determine whether to perform
!!                      diagonalisation in the final space -
!!                      used when considering excitation from a
!!                      valence-excited neutral state to a
!!                      core-excited neutral state
!!
!!$integer lancguess - integer value used to determine how the initial
!!                     Lanczos vectors are set:
!!
!!                     lancguess=1 <-> 1h1p (1h1p, 2h2p) unit vectors
!!                                     for ADC(2)-s (ADC(2)-x)
!!                               2 <-> ADC(1) eigenvectors for
!!                                     ADC(2)-s (plus 2h2p unit
!!                                     vectors for ADC(s)-x)
!!                               3 <-> linear combinations of 1h1p and
!!                                     2h2p unit vectors for ADC(2)-s
!!                               4 <-> linear combinations of ADC(1)
!!                                     eigenvectors and 2h2p unit
!!                                     vectors for ADC(2)-s
!!
!!$real maxmem - memory threshold (in Mb) before disk-based algorithms
!!               are used
!!
!!$logical ldipole - flag used to determine whether to calculate the
!!                   excited-state dipole moments
!!
!!$real dipmom - dipole moments for the initial space states
!!$real dipmom_f - dipole moments for the final space states

  logical                              :: debug
  character(1)                         :: tranmom,tranmom2
  character(36)                        :: lancname,davname,davname_f
  integer                              :: hinit,nirrep,nirrep2,method,&
                                          method_f,davstates,lancstates,&
                                          numinista,norder,info,&
                                          statenumber,denord,davstates_f
  integer, parameter                   :: nhcentre=40
  integer, dimension(0:nhcentre)       :: hcentre
  real(d)                              :: mspacewi,mspacewf
  real(d), parameter                   :: minc=1e-12_d
  real(d)                              :: dlim
  integer, dimension(:), allocatable   :: roccnum
  integer, dimension(8,8)              :: MT
  integer, dimension(400)              :: stvc_lbl
  integer, dimension(:), allocatable   :: stvc_mxc
  integer                              :: ninista
  logical                              :: ladc1guess,ladc1guess_f
  logical                              :: lnoise,lnoise_f
  logical                              :: lsubdiag,lsubdiag_f
  real(d), dimension(:), allocatable   :: adc1en,adc1en_f
  real(d)                              :: davtol,davtol_f
  logical                              :: lcvs,lcvsfinal
  integer, dimension(nhcentre)         :: icore,iexpfrz
  integer                              :: ncore
  logical                              :: lfakeip
  integer                              :: ifakeorb
  integer, dimension(:), allocatable   :: ifakeex
  real(d)                              :: expfakeip
  logical                              :: ld2
  logical                              :: lscf
  integer                              :: ivpqrs
  real(d), dimension(:,:), allocatable :: density
  logical                              :: ltdm_gs2i
  logical                              :: lifrzcore,lffrzcore
  logical                              :: ldiagfinal
  real(d)                              :: maxmem
  logical                              :: ldipole
  real(d), dimension(:), allocatable   :: dipmom,dipmom_f

!!$************************************************
!!$**********Physical Constants********************
!!$************************************************

!!$ abohr - Bohr radius [cm]
!!$ fsconstinv - inverted fine structure constant
!!$ os2cs - oscillator strength [a.u.] to cross-section [Mb] conversion factor
!!$ omega - photon energy, required by Stieltjes_phi, photoionisation routine.
!!$ eh2ev - Hartree to eV conversion factor

!  real(d), parameter :: abohr=5.2918e-9
  real(d), parameter :: fsconstinv=137.0d0
  real(d), parameter :: os2cs=4.0347443d0
  real(d), parameter :: omega=3.0d0
  real(d), parameter :: eh2ev=27.2113845d0

!!$************************************************
!!$**********Lanczos Parameters********************
!!$************************************************  
  integer :: ncycles,lmain,maxblock
  integer :: lancguess
  integer :: orthotype
  real(d) :: tdtol
  logical :: ldynblock

!!$************************************************
!!$**********Davidson Parameters*******************
!!$************************************************
  integer     :: maxiter,dmain,maxiter_f,dmain_f,&
                 ndavcalls,eigentype,solver,solver_f,&
                 precon,precon_f,maxsubdim,maxsubdim_f
  real(d)     :: davtarg
  
  logical     :: lextdiag,ldfl,ldfl_f

!!$************************************************
!!$**********Relaxation Parameters*****************
!!$************************************************
  integer :: kdim,kdim_f,guessdim,guessdim_f
  real(d) :: stepsize,stepsize_f,siltol,siltol_f


  
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

!!$************************************************
!!$**********SCF Parameters************************
!!$************************************************
integer :: scfiter

!!$************************************************
!!$**********GAMESS Parameters*********************
!!$************************************************
integer, parameter                          :: maxao=100
integer, parameter                          :: maxprim=100
integer                                     :: ncoo,difftype,contcent
integer, dimension(5)                       :: ndiff
integer, dimension(:), allocatable          :: naogms
integer, dimension(:,:), allocatable        :: nprim,aotype,ilquant
real(d), dimension(:,:,:), allocatable      :: aocoeff,aoexp
real(d), dimension(:,:), allocatable        :: diffexp
character(len=60)                           :: basname
character(len=2), dimension(:), allocatable :: atlbl
logical                                     :: lrungamess,ldiffcom

!!$************************************************
!!$**********Dyson orbital parameters**************
!!$************************************************
integer                              :: nbas_ao,ccent,cshell
integer                              :: dysirrep,dysdiag
integer, dimension(5)                :: dysout
real(d)                              :: dyslim
real(d), dimension(:,:), allocatable :: ao2mo
real(d), dimension(10,10)            :: ftransmat
character(len=1), dimension(0:3)     :: shlbl
logical                              :: ldyson,lrmatom

! ezdyson input file parameters
integer :: lmax,nelen,ngrdpnts
real(d) :: zcore,eleni,elenf,grdi,grdf

!!$************************************************
!!$**********Target State Matching*****************
!!$************************************************
real(d)                          :: detthrsh,ovrthrsh
character(len=120)               :: detfile,mofile
logical                          :: ltarg

!!$************************************************
!!$**********RIXS Calculations*********************
!!$************************************************
real(d), dimension(:,:,:), allocatable :: dpl_all
logical                                :: lrixs

!!$************************************************
!!$**********TPA Calculations**********************
!!$************************************************
integer, dimension(2)                  :: tpblock
real(d), dimension(3)                  :: tdmgsi
real(d), dimension(:,:), allocatable   :: travec_ic,&
                                          travec_iv,&
                                          tdmgsf
real(d), dimension(:,:,:), allocatable :: travec_fc,&
                                          travec_fv
real(d), dimension(:), allocatable     :: edavf
logical                                :: ltpa

!!$************************************************
!!$**********Autospec Calculations*****************
!!$************************************************
real(d)                               :: tfinal,tout,&
                                         autotol
integer                               :: autoord
logical                               :: lautospec

!!$************************************************
!!$**********Filter Diagonalisation State**********
!!$******************Calculations******************
!!$************************************************
integer                               :: nfbas,neig,&
                                         nsel,iwfunc
integer, dimension(:), allocatable    :: isel
real(d), dimension(:,:), allocatable  :: fbas2eig
real(d), dimension(2)                 :: ebound
character(len=60)                     :: fdiagdat,&
                                         fdiagsel
logical                               :: lfdstates

!!$************************************************
!!$************External Electric Field*************
!!$************************************************
integer*8, dimension(3)      :: nel_dip
integer, dimension(3)        :: nbuf_dip
integer, parameter           :: mxenvpar=20
integer                      :: nenvpar
integer                      :: ipulse
integer                      :: ienvelope
real(d), dimension(3)        :: d00
real(d), allocatable         :: d0j(:,:)
real(d), allocatable         :: dij(:,:,:)
real(d), dimension(3)        :: pulse_vec
real(d)                      :: proptol
real(d)                      :: freq
real(d)                      :: strength
real(d), dimension(mxenvpar) :: envpar
logical                      :: lpropagation

!!$************************************************
!!$***************CAP Parameters*******************
!!$************************************************
integer                       :: icap
integer                       :: capord
integer*8                     :: nel_cap
integer                       :: nbuf_cap
integer, dimension(2)         :: nrad,nang
real(d)                       :: capstr
real(d), dimension(3)         :: boxpar
real(d)                       :: densthrsh
real(d)                       :: w00
real(d), allocatable          :: w0j(:)
real(d), allocatable          :: wij(:,:)
logical                       :: lcap
logical                       :: lprojcap
logical                       :: lautobox

!!$************************************************
!!$***************ADC(1) Hamiltonian***************
!!$************************************************
real(d), allocatable :: h1(:,:)

!!$************************************************
!!$**********Other Parameters**********************
!!$************************************************
integer                                     :: natm
integer(d), dimension(:), allocatable       :: nrec_omp
real(d), dimension(:), allocatable          :: xcoo
real(d), dimension(:,:,:,:), allocatable    :: fvpqrs
character(len=2), dimension(:), allocatable :: aatm
character(len=3)                            :: pntgroup

!!$************************************************
!!$**********Hamiltonian matrices******************
!!$************************************************
character*1        :: hamflag

!!$************************************************
!!$*******Record size for disk-based storage*******
!!$************************************************
! 10 Mb per record
integer, parameter :: buf_size=655360

! 1 Mb per record
!integer, parameter :: buf_size=65536

! 128 Kb per record
!integer, parameter :: buf_size=8192

end module parameters
