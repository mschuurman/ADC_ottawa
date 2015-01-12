!
!  Direct evaluation of quantum 1-particle partition functions on a grid using
!  matrix exponentiation. Matrix exponentiation is performed in real space, and
!  relies critically on the sparseness of the kinetic energy operator.
!
!  The partition function is given by the trace of matrix exponential of 
!  (-beta*H), where H is the Hamiltonian matrix. We calculate the exponent
!  using explicit expansion of the exponential into series of matrix powers.
!  As a result, a complete density matrix is also available at the end of the
!  process. If any observables are of interest, these can be calculated from
!  traced of products of the density matrix and the property matrix. For the
!  time being, we are simply dumping the diagonal of the density matrix, for
!  graphical illustration purposes.
!
!  This example is hydrogen absorption in a variety of materials, using both
!  quantized ideal gas and quantized-liquid density functional theory, in 
!  orthorhombicc cells.
!
!  A number of interaction potentials is built-in, including (see comments
!  for centre_types() array below):
!
!  1. C-H2 dispersion-only exp-6 interaction for graphitic carbons ('CARB')
!  2. A 6-12 dispersion-only interaction, with user-controllable
!     potential, suitable for non-polar hosts. ('LJ612')
!  3. For polar hosts, charge-polarizability and charge-quadrupole
!     interactions can be handled approximately, assuming rapid
!     rotational equilibration at each adsorption site (see h2_rotation.f90
!     for comments).A ('Q612')
!
!  The 6-12 parameters are for the potential in the form:
!
!    v_{612}(r) = 4\epsilon (  -(\sigma/r)^6 + (\sigma/r)^12 )
!
!  with \epsilon in kcal/mol and \sigma in Angstroms
!
!  The electrostatic contribution is described by atom-centered effective
!  charges, given in the units of positron charge (absolute electron 
!  charge |e|).
!
  module partition
    use accuracy
    use math
    use vector
    use sparse
    use lapack
    use timer
    use convolution
    use h2_thermo
    use h2_potential
    use h2_rotation
    implicit none
    !
    !  ==== Fixed parameters, changing which does not make any sense ====
    !
    integer(ik), parameter :: unit_struct = 35
    integer(ik), parameter :: unit_dens   = 36
    integer(ik), parameter :: unit_mol    = 37
    integer(ik), parameter :: max_exclude = 60                  ! Max number of excluded volume 
                                                                ! elements.
    real(rk), parameter    :: bond_cc     = 1.421_rk / abohr    ! Experimental C=C bond in graphite
    real(rk), parameter    :: dens_fact   = 1e30_rk/(abohr**3 * N_Avogadro) 
                                                                ! conversion factor: 1/Bohr^3 -> mole/m^3
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik) :: verbose       = 2                           ! Verbosity level
    real(rk)    :: mass          = (1.0_rk+1838.7_rk) * 2.0_rk ! Particle mass, in e.m.u
    real(rk)    :: mass_guess    = 1e5_rk                      ! Particle mass used for "classical" guess
    integer(ik) :: mass_steps    = 3                           ! Number of mass steps for the 'graded' guess
    real(rk)    :: temperature   = 300._rk                     ! Temperature at which partition function
                                                               ! is required, in Kelvin
    real(rk)    :: v_clip        = 5000._rk                    ! Clip energy for the potential evaluation,
                                                               ! in Kelvin
    integer(ik) :: max_terms     = 180                         ! Maximum number of terms to include in
                                                               ! matrix exponentiation series
    real(rk)    :: power_screen  = 1e-3_rk                     ! Screen out smaller terms in matrix powers
    real(rk)    :: dens_cutoff   = 0.01_rk/dens_fact           ! Zero densities make our F[] go crazy; do not
                                                               ! let molar volume drop below 100 m^3/mole
    integer(ik) :: n_points(3)   = (/ 25, 26, 26 /)            ! Number of sampling points along each
    integer(ik) :: sum_cells(3)  = (/ 5, 5, 5 /)               ! Number of cells to include in potential
    integer(ik) :: plot_cells(3) = (/ 1, 1, 1 /)               ! Cells to include in the background image
    integer(ik) :: oversample_vads = 2                         ! Values greater than 1 request oversampling
                                                               ! of adsorption potential over volume elements,
                                                               ! to minimize spurious potential variations.
                                                               ! Cost rises cubically with oversampling
    !
    !  ==== Parameters relates to thermodynamic ensemble ====
    !
    character(len=10) :: mode     = 'ideal'                    ! Calculation mode. Allowed values:
                                                               ! 'ideal' - use the ideal gas approximation
                                                               ! 'dft'   - use liquid density functional theory
    character(len=10) :: guess    = 'classical'                ! Initial guess for dft calculation, either of:
                                                               ! 'classical' - use classical-liquid DFT
                                                               ! 'graded'    - use gradual mass switch-on from
                                                               !               the large (classical) limit to
                                                               !               the final quantum mass. 'classical'  
                                                               !               is the special case of 'graded' with
                                                               !               one step.
                                                               ! 'zero'      - use zero Vxc
                                                               ! 'read'      - read initial guess for the
                                                               !               density from file.
    character(len=10) :: scf_mode = 'both'                     ! Primary variable for SCF convergence driver, either:
                                                               ! 'density'   - Guest density
                                                               ! 'potential' - Effective potential
                                                               ! 'both'      - Use mixing for both density and potential
    real(rk)    :: molar_volume   = 0.0224_rk                  ! Desired average molar volume, in m^3 per mole
                                                               ! The value refers to the unscaled unit cell.
    real(rk)    :: rho_mix        = 0.10_rk                    ! Density mixing parameter
    real(rk)    :: rho_mix_guess  = 0.03_rk                    ! Density mixing parameter for the initial guess
    real(rk)    :: veff_mix       = 0.10_rk                    ! Effective potential mixing parameter
    real(rk)    :: veff_mix_guess = 0.03_rk                    ! Effective potential mixing parameter for the initial guess
    real(rk)    :: eps_f          = 1e-1_rk                    ! Desired convergence in free energy.
    real(rk)    :: eps_f_guess    = 1.0_rk                     ! Desired convergence in free energy for the initial guess
    real(rk)    :: eps_rho        = 1e-4_rk                    ! Desired convergence in density.
    real(rk)    :: eps_rho_guess  = 1e-3_rk                    ! Desired convergence in density for the initial guess
    integer(ik) :: rho_iterations = 80                         ! Number of density iterations
    integer(ik) :: rho_iterations_guess = 160                  ! Number of density iterations for the initial guess
    !
    !  ==== Excluded spherical volume elements ====
    !
    integer(ik) :: n_exclude     = 0               ! Number of excluded VEs (spherical)
    real(rk)    :: xyzr_exclude(4,max_exclude)     ! X/Y/Z/R (Bohr) of excluded VEs
    !
    !  ==== Hack: place a hard wall at the origin ====
    !
    real(rk)    :: hard_wall_width = -1._rk        ! Width < 0 means no hard wall.
                                                   ! The wall is centered at the origin, and is
                                                   ! perpendicular to Z direction. I told you this is
                                                   ! a hack!
    !
    !  ==== Hack: place a hydrogen-like potential at the origin ====
    !
    logical     :: central_probe   = .false.       ! Activate for calculating radial distribution profiles
    !
    !  ==== Structure data ====
    !
    character(len=80)     :: coord_file = ' '      ! Name of the file containing structure
                                                   ! (blank if inline)
    real(rk)              :: lat_vec(3)            ! Orthorhombic periodic lattice parameters (Bohr)
                                                   ! Note that the unit cell MUST be centered at
                                                   ! the origin - a lot of things won't work otherwise.
    real(rk)              :: lat_clip(3) = (/ 0._rk, 0._rk, 0._rk /)
                                                   ! Symmetrically reduce the box from a given direction 
                                                   ! by lat_clip(:) amount. The reduction affects only 
                                                   ! the wavefunction, and should be used when we -know-
                                                   ! that the potential in the excluded area is so large
                                                   ! that the wavefunction -will- be negligible.
    real(rk)              :: cell_scale = 1.0_rk   ! Uniform cell scaling factor (for pressure calculation)
    integer(ik)           :: n_centres             ! Number of centres within the unit cell
    real(rk), allocatable :: centres(:,:)          ! Coordinates of the atoms (Bohr)
    character(len=10), allocatable &
                          :: centre_labels(:)      ! Atom labels, purely for cosmetics
    character(len=10), allocatable &
                          :: centre_types(:)       ! Potential type on the atom. Can be one of:
                                                   !   'CARB'  - Hard-coded graphitic carbon (this is the defalt)
                                                   !   'LJ612' - 6-12 Lennard-Jones interaction, requires two
                                                   !             numbers to follow (sigma and epsilon)
                                                   !   'Q612'  - Electrostatic + 6-12 Lennard-Jones interaction,
                                                   !             requires three numbers to follow (sigma, epsilon, and
                                                   !             the effective charge). Note that the rotational
                                                   !             electrostatic terms are treated approximately!
    real(rk), allocatable :: centre_params(:,:)    ! Additional parameters defining the potential
    integer(ik), parameter:: max_centre_params = 3 ! Maximum number of numerical parameters for a single centre
    character(len=80)     :: vads_read_file = ' '  ! Name of the file, from which adsorption potential (in Kelvin)
                                                   ! should be loaded; blank means compute from the structure.
    character(len=80)     :: vads_file = 'vads.dx' ! Name of the file for the adsorption potential (in Kelvin)
                                                   ! (blank if plot is not needed)
    character(len=80)     :: veff_file = 'veff.dx' ! Name of the file for the effective potential (blank
                                                   ! if the effective potential is not needed)
    character(len=80)     :: dens_file = 'rho.dx'  ! Name of the file for the density (blank if density
                                                   ! need not be saved)
    character(len=80)     :: image_file = 'graphene.dx' ! Name of the structure image
    character(len=80)     :: ascii_file = ' '      ! Name of the ascii dump of the density & energy density
    character(len=80)     :: guess_file = 'guess.dx'
                                                   ! Name of the file containing density guess. The density
                                                   ! need not not be normalized, but the grids must match.
                                                   ! Data is expected in the format produced by dump_dx_field,
                                                   ! but most of the header field will be simply skipped.
    character(len=80)     :: v12_form = 'diep-sph' ! Choice of the inter-particle interaction potential
                                                   ! Possibilities are:
                                                   !  'none'     = LIE-QLDFT approximation, no v12 at all
                                                   !  'diep-'    = Fit to the 0-0-0 part of Diep & Johnson's
                                                   !               H2-H2 potential, flatlined for r<v12_cutpoint
                                                   !  'diep-sph' = As 'diep-', but matched to an arc with the
                                                   !               maximum height of v12_sphmax at the origin.
                                                   !               Matching point is at the zero crossing.
    real(rk)              :: v12_scale = 1.0_rk    ! Overall scale for the v12 interaction. Useful in
                                                   ! cases of severe convergence problems.
    real(rk)              :: v12_cutpoint = 2.850_rk  ! Best for v12 functional term dependent on n
                                                   ! Flat-line the repulsive part of the v12 potential
                                                   ! at r<v12_cutpoint.
    real(rk)              :: v12_sphmax  = 180._rk ! The maximum height of the capping sphere for 
                                                   ! v12_form=='diep-sph', in Kelvin
    real(rk)              :: v12_epsilon = 1e-9_rk ! Cut-off level for v12 potential (atomic units)
    character(len=80)     :: wda_shape = 'Fermi-Dirac'
                                                   ! Possibilities are:
                                                   !   'Fermi-Dirac' = 1/(1+Exp((r-r0)/d))
                                                   !   'Gaussian'    = Exp(-r**2/r0**2) / ( pi**1.5 * r0**3 )
                                                   !   'Delta'       = Dirac's delta function
    real(rk)              :: wda_r0    = 2.500_rk  ! Characteristic withth of the filter
    real(rk)              :: wda_d     = 0.001_rk  ! For Fermi-Dirac, rate of switch-off of the filter
    real(rk)              :: wda_epsilon = 1e-7_rk ! Cut-off for the density averaging filter
    logical               :: debug_convergence = .false.
                                                   ! If true, dump density and effective potential on each iteration
    character(len=240)    :: comment = ' '         ! Descriptive comment, will be chopped up at 80-byte boundaries
    character(len=80)     :: functional = 'user'   ! Set up an internally-parameterized functional. Possible values:
                                                   !   'user'   = The user knows that they are doing; do not
                                                   !              adjust the parameters
                                                   !   'lie-0'  = Purely local "LIE-0" functional
                                                   !   'lie-1'  = WDA-style "LIE-1" functional
                                                   ! A number of variant LIE-1 functionals is recognized in 
                                                   ! process_input_shortcuts - do not use them unless you know
                                                   ! what you are doing
    logical               :: let_it_be = .false.   ! Allow calculation to continie even if some of the input 
                                                   ! appears to be non-sensical.
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    real(rk)                 :: beta            ! 1/(k T) in atomic units
    real(rk)                 :: g_box (2,3)     ! Periodic box extent: min and max in each direction
    real(rk)                 :: g_step(3)       ! This is only the part of the box which contains grid
                                                ! points - excluded volume is not gridded!
    integer(ik)              :: n_states        ! Dimension of the Hamiltonian matrix
    integer(ik), allocatable :: neighbour(:,:)  ! Neighbour list - indices of non-zero elements in
                                                ! each row of the kinetic matrix
    real(rk)                 :: kmat(0:6)       ! Values of non-zero elements in the kinetic matrix
                                                ! Displacements 1:6 in kmat must match directions 1:6
                                                ! in the neighbour table
    type(SparseMatrixT)      :: hpow            ! Powers of the Hamiltonian matrix
    real(rk), allocatable    :: vads(:)         ! External potential matrix
    real(rk), allocatable    :: hdia(:)         ! Diagonal of the Hamiltonian matrix' powers
    real(rk), allocatable    :: v0dia(:)        ! Diagonal of the external potential matrix
    real(rk), allocatable    :: vdia(:)         ! Diagonal of the effective potential matrix (incl. XC potential)
    real(rk), allocatable    :: old_vdia(:)     ! Diagonal of the effective potential from previous iteration
    real(rk), allocatable    :: v12(:)          ! Mean-field potential
    real(rk), allocatable    :: eps_xc(:)       ! Exchange-correlation energy density
    real(rk), allocatable    :: vxc(:)          ! Exchange-correlation potential
    real(rk), allocatable    :: vxc_n(:)        ! Scratch pad needed in evaluation of vxc() - density part
    real(rk), allocatable    :: vxc_nave(:)     ! Scratch pad needed in evaluation of vxc() - smoothed density part
    real(rk), allocatable    :: dens(:)         ! Diagonal of the density matrix
    real(rk), allocatable    :: old_dens(:)     ! Diagonal of the density matrix from previous iteration
    real(rk), allocatable    :: smooth_dens(:)  ! Smoothed density for evaluation of the exchange-correlation potential
    real(rk), allocatable    :: xyz (:,:)       ! Coordinates of the "grid" points - only needed during
                                                ! the initial set-up phase.
    real(rk)                 :: particle_count  ! Number of particles in the simulation box
                                                ! needed to achieve the target density.
    real(rk)                 :: q_potential     ! Partition function with the potential included
    real(rk)                 :: q_accuracy      ! Estimate of the numerical uncertaintly in Q
    real(rk)                 :: q_free          ! Partition function of the empty box - control
    real(rk)                 :: q_free_integral ! Analytical result for q_free - control, continuous limit
    real(rk)                 :: q_free_limit    ! Analytical result for q_free - control, exact limit
    real(rk)                 :: cell_volume     ! Simulation cell volume, in Bohr^3
    real(rk)                 :: p_min, p_max    ! Highest and lowest "local pressure"
    real(rk)                 :: max_dens_delta  ! Maximum change in density on the last iteration
    real(rk)                 :: max_veff_delta  ! Maximum change in effective potential on the last iteration
    integer(ik)              :: debug_sequence  ! Incremented by 1 for each SCF loop call
    type(ConvolutionT)       :: cnv_v12         ! Convolution for v12 potential
    type(ConvolutionT)       :: cnv_wda         ! Density-smoothing convolution
    real(rk)                 :: v12_uniform     ! Integral of v12 over the entire space
                                                ! Components of the total free energy, in Hartree/particle
    real(rk)                 :: energy_v12      ! Mean-field energy of the density distribution
    real(rk)                 :: energy_vext     ! Energy in the external potential
    real(rk)                 :: energy_fxc      ! Exchange-correlation energy
    real(rk)                 :: energy_vxc      ! "Energy" of the density in the exchange-correlation potential
    real(rk)                 :: energy_kin      ! Kinetic energy
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /partfun/ verbose, &
                       mode, functional, v12_form, wda_shape,      &
                       guess, molar_volume,                        &
                       rho_iterations, rho_iterations_guess,       &
                       eps_f, eps_f_guess, eps_rho, eps_rho_guess, &
                       scf_mode,                                   &
                       rho_mix, rho_mix_guess,                     &
                       veff_mix, veff_mix_guess,                   &
                       mass, mass_guess, mass_steps,               &
                       temperature, v_clip,                        &
                       max_terms, power_screen, dens_cutoff,       &
                       n_points, sum_cells, oversample_vads,       &
                       lat_vec, lat_clip, cell_scale,              &
                       coord_file,                                 &
                       n_exclude, xyzr_exclude,                    &
                       plot_cells,                                 &
                       vads_file, veff_file, dens_file,            &
                       image_file, ascii_file, guess_file,         &
                       vads_read_file,                             &
                       hard_wall_width, central_probe,             &
                       v12_scale, v12_epsilon, v12_cutpoint,       &
                       v12_sphmax,                                 &
                       wda_r0, wda_d, wda_epsilon,                 &
                       debug_convergence,                          &
                       comment, let_it_be
    !
    !  ==== End of global data ====
    !
    contains
    !   
    !  Pairwise potential for C-H2 interaction, where H2 is taken as a
    !  feature-less particle at the centre of mass of H2.
    !
    function vdw_potential(r) result(v)
      real(rk), intent(in) :: r(3)
      real(rk)             :: v
      !
      real(rk)             :: rx
      !
      rx = sqrt(sum(r**2)) * abohr  ! Bohr -> Angstroms
      rx = max(1.8_rk,rx)           ! Cut off hard repulsive tail
      !
      !  Potential C - fit to MP2/cc-pvtz in coronene
      ! 
      v = -400.37_rk/rx**6 + 25352.0_rk * exp(-3.5763_rk * rx)
      !
      v = v / ( 627.51_rk ) ! kcal/mol -> Hartree
    end function vdw_potential
    !
    !  Lennard-Jones 6-12 potential, parameterized by the sigma and epsilon
    !
    function lj612_potential(r,par) result(v)
      real(rk), intent(in) :: r  (3)  ! Coordinates, in Bohr
      real(rk), intent(in) :: par(2)  ! LJ parameters: 1=sigma, in Angstrom; 2=epsilon, in kcal/mol
      real(rk)             :: v       ! Interaction potential, in Hartree
      !
      real(rk)             :: rx
      !
      rx = sqrt(sum(r**2)) * abohr             ! Bohr -> Angstroms
      rx = max(0.5_rk,rx)                      ! Cut off the singularity
      !
      rx = par(1)/rx                           ! Convert into units of sigma
      v  = 4._rk * par(2) * ( rx**12 - rx**6 ) ! kcal/mol
      v  = v / ( 627.51_rk )                   ! kcal/mol -> Hartree
    end function lj612_potential
    !
    !  First and second derivatives of the Coulomb potential
    !
    !  We will stop evaluation of the denominator in rotational corrections 
    !  at the point where the repulsive vdW energy will reach about 1750 times 
    !  the maximum vdW attraction.
    !
    subroutine accumulate_ef_efg(r,rc,q,ef,efg)
      real(rk), intent(in) :: r(3)        ! Coordinates of the probe (charge at origin), in Bohrs
      real(rk), intent(in) :: rc          ! 6-12 distance parameter, in Angstroms
      real(rk), intent(in) :: q           ! Charge, in |e| units
      real(rk), intent(inout) :: ef(3)    ! Gradient of the potential, NOT THE FIELD
      real(rk), intent(inout) :: efg(3,3) ! Second derivative of the potential
      !
      real(rk)    :: r2   ! Square of the distance from the centre, for use in the NUMERATOR
      real(rk)    :: rd   ! Distance from the centre, for use in the DENOMINATOR
      integer(ik) :: i, j
      !
      r2 = sum(r**2)
    ! rd = max(0.75_rk*rc/abohr,sqrt(r2))
      rd = max(0.60_rk*rc/abohr,sqrt(r2))
      ef = ef - q * r/rd**3
      do i=1,3
        do j=1,3
          efg(i,j) = efg(i,j) + q*3*r(i)*r(j)/rd**5
        end do
        efg(i,i) = efg(i,i) - q*r2/rd**5
      end do
    end subroutine accumulate_ef_efg
    !
    !  Sample total interaction potential. Input coordinate is assumed to be
    !  inside the unit cell (no checking!)
    !
    !  NOTE: This routine should not be called directly - pot_create()
    !        is responsible for oversampling the potential
    !
    function pot_create_sample(x) result(vr)
      real(rk), intent(in) :: x(3)        ! Coordinates of the point
      real(rk)             :: vr          ! Interaction potential
      integer(ik)          :: ixc         ! Exclusion zone counter
      real(rk)             :: push_x(3)   ! Adjusted coordinates, pushed "outside"
                                          ! of exclusion zones if necessary
      real(rk)             :: dx(3)       ! Relative position wrt exclusion zone centre
      real(rk)             :: lat_off(3)  ! Replicated cell offset
      real(rk)             :: rx          ! Displacement from the exclusion centre
      real(rk)             :: rs          ! Exclusion radius
      integer(ik)          :: ia, ib, ic  ! Replication counters
      integer(ik)          :: iat         ! Centres within the unit cell
      logical              :: have_efield ! Activated if electric field parameters are present
      real(rk)             :: ef(3)       ! Electric field (actually, gradient of the potential)
      real(rk)             :: efg(3,3)    ! Electric field gradient (actually, second derivative)
      real(rk)             :: vrot        ! Rotational correction
      !
      ! vr = -1/beta ; return !
      !
      !  Check for the exclusion zones. The check here works only for non-overlapping
      !  exclusion zones - otherwise, there is no guarantee we'll clear all exclusion
      !  zones in the system.
      !
      push_x = fold_vector(x/cell_scale)
      !
      exclusion: do ixc=1,n_exclude
        dx = fold_vector(push_x - xyzr_exclude(1:3,ixc))
        !
        rs = xyzr_exclude(4,ixc)
        rx = sqrt(sum(dx**2))
        if ( rx>=rs ) cycle exclusion
        !
        !  Point is within the exclusion zone, push it out. If point happens to
        !  be at the origin, arbitraliry push it along the third coordinate ("up")
        !
        if ( rx>=0.01_rk ) then
          push_x = xyzr_exclude(1:3,ixc) + rs*dx/rx
        else
          push_x = xyzr_exclude(1:3,ixc) + (/ 0._rk, 0._rk, rs /)
        end if
        !
        !  Pushing could lead to a vector outside the unit cell, renormalize
        !
        push_x = fold_vector(push_x)
      end do exclusion
      !
      !  push_x(:) coordinate should now be either outside or on the boundary of
      !            the exclusion zones. 
      !            
      !  Sum the potential over enough surrounding zones
      !
      vr  = 0
      ef  = 0
      efg = 0
      have_efield = .false.
      sum_a: do ia=-sum_cells(1),sum_cells(1)      ! Replicate the cell in x direction
        sum_b: do ib=-sum_cells(2), sum_cells(2)   ! Replicate the cell in y direction
          sum_c: do ic=-sum_cells(3),sum_cells(3)  ! Replicate the cell in z direction
            lat_off = lat_vec * (/ ia, ib, ic /)
            do iat=1,n_centres
              dx = push_x - (centres(:,iat) + lat_off)
              select case (centre_types(iat))
                case default
                  write (out,"('Atom ',i6,' has unrecognized potential type ',a)") &
                         iat, trim(centre_types(iat))
                  stop 'partition%pot_create_sample - bad potential type'
                case ('CARB')
                  vr = vr + vdw_potential(dx)
                case ('LJ612')
                  vr = vr + lj612_potential(dx,centre_params(:2,iat))
                case ('Q612')
                  vr = vr + lj612_potential(dx,centre_params(:2,iat))
                  call accumulate_ef_efg(dx,centre_params(1,iat),centre_params(3,iat),ef,efg)
                  have_efield = .true.
              end select
            end do
          end do sum_c
        end do sum_b
      end do sum_a
      !
      !  If electric field is present, calculate rotational correction
      !
      if (have_efield) then
        vrot = h2_vrot(temperature,ef,efg)
        ! write (out,"('Rotational term = ',g12.6,' rest = ',g12.6,' total = ',g12.6)") &
        !        vrot, vr, vrot + vr
        vr = vr + vrot
      end if
      !
      !  Add central probe potential if it was requested
      !
      if (central_probe) then
        vr = vr + h2_v12(sqrt(sum(push_x**2)))
      end if
      !
    end function pot_create_sample
    !
    !  Return abscissas and weights for Gauss-Legendre integration; naive
    !  version
    !
    subroutine gauss_legendre (order,x,w)
      integer(ik), intent(in) :: order ! Desired integration order
      real(rk), intent(out)   :: x(:)  ! Abscissas in the [-1:1] range
      real(rk), intent(out)   :: w(:)  ! Weights
      !
      select case (order)
        case default
          write (out,"('Gauss-Legendre integration order ',i5,' is not implemented')") order
          stop 'partition-in-small-box%gauss_legendre - invalid order'
        case (2) ! 2 points
          x(1:2) = (/ -0.5773502691896257645091487805019575_rk, &
                       0.5773502691896257645091487805019575_rk /)
          w(1:2) = (/ 1.0_rk, 1.0_rk /)
        case (3) ! 3 points
          x(1:3) = (/ -0.7745966692414833770358530799564799_rk, 0.0_rk, &
                       0.7745966692414833770358530799564799_rk /)
          w(1:3) = (/ 0.5555555555555555555555555555555556_rk, &
                      0.8888888888888888888888888888888889_rk, &
                      0.5555555555555555555555555555555556_rk /)
        case (4) ! 4 points
          x(1:4) = (/ -0.8611363115940525752239464888928095_rk, &
                      -0.3399810435848562648026657591032447_rk, &
                       0.3399810435848562648026657591032447_rk, &
                       0.8611363115940525752239464888928095_rk /)
          w(1:4) = (/ 0.3478548451374538573730639492219994_rk, &
                      0.6521451548625461426269360507780006_rk, &
                      0.6521451548625461426269360507780006_rk, &
                      0.3478548451374538573730639492219994_rk /)
      end select
    end subroutine gauss_legendre
    !
    !  Calculate adsorption potential for a volume element, possibly
    !  oversampling the analytical potential.
    !
    function pot_create(x) result(vr)
      real(rk), intent(in) :: x(3)       ! Coordinates of the point
      real(rk)             :: vr         ! Interaction potential
      !
      integer(ik) :: sx, sy, sz ! Sampling indices
      real(rk)    :: xt(3)      ! Current coordinate
      real(rk)    :: dx(4)      ! 1D displacements; We support up to 4-point formulae
      real(rk)    :: wg(4)      ! 1D weights
      !
      !  Hard-wall term, if present. It will never get triggered if 
      !  hard_wall_width is < 0.
      !
      if (abs(x(3))<=0.5_rk*hard_wall_width) then
        vr = 0.1_rk ! In Hartrees - this is impossibly huge for adsorption
        return
      end if
      !
      !  Quick return if no oversampling requested
      !
      if (oversample_vads<=1) then
        vr = pot_create_sample(x)
        return
      end if
      !
      !  Sample potential through the volume element. We use product of
      !  Gauss formulae to do the integration. The offsets and weights
      !  are tabulated for the [-1:1] interval; we'll rescale to the
      !  [-1/2:1/2] interval we need
      !
      dx = 0 ; wg = 0
      call gauss_legendre(oversample_vads,dx,wg)
      dx = 0.5_rk * dx
      wg = 0.5_rk * wg
      !
      !  Calculate the average
      !
      vr = 0
      sample_z: do sz=1,oversample_vads
        xt(3) = x(3) + dx(sz) * g_step(3)
        sample_y: do sy=1,oversample_vads
          xt(2) = x(2) + dx(sy) * g_step(2)
          sample_x: do sx=1,oversample_vads
            xt(1) = x(1) + dx(sx) * g_step(1)
            vr = vr + wg(sx)*wg(sy)*wg(sz)*pot_create_sample(xt)
          end do sample_x
        end do sample_y
      end do sample_z
    end function pot_create
    !
    !  Fold a vector within a zero-centered unit cell 
    !
    function fold_vector(xin) result(x)
      real(rk), intent(in) :: xin(:)  ! Arbitrary vector
      real(rk)             :: x  (3)  ! Vector folded within the origin cell
      !
      x = xin
      normalize_x: do while (any(abs(x)>lat_vec/2))
        where (x>lat_vec/2) 
          x = x - lat_vec
        end where
        where (x<-lat_vec/2)
          x = x + lat_vec
        end where
      end do normalize_x
    end function fold_vector
    !
    !  Read unit cell structure, formatted as an extended XYZ file.
    !  The input accepts a parameterized adsorption potential as a hack:
    !  the standard element X Y Z entries can be optionally followed by
    !  a character string defining the potential type, together with zero
    !  or more numerical parameters defining the potential.
    !
    subroutine read_xyz_structure
      integer(ik)           :: inp     ! I/O channel used for reading the structure data
      integer(ik)           :: iat     ! Atom number
      character(len=300)    :: buf     ! Temporary buffer - needed for making sure our entries
                                       ! obey the line boundaries, while still allowing the use
                                       ! of optional free-format parts. 
      integer(ik)           :: ios
      !
      !  If the input file name was provided in the namelist, use it.
      !  Otherwise, continue reading the in-line input.
      !
      inp = input
      if (coord_file/=' ') then
        inp = unit_struct
        open(inp,file=trim(coord_file),status='old',action='read')
      end if
      !
      !  Read the XYZ file. The first line contains the number of atoms.
      !
      n_centres = 0
      read(inp,*) n_centres
      if (n_centres==0) return
      allocate(centres(3,n_centres),centre_labels(n_centres), &
               centre_types(n_centres),centre_params(max_centre_params,n_centres))
      !
      !  Defaults for the potential, if none is seen in the input
      !
      centre_types   (:) = 'CARB'
      centre_params(:,:) = 0._rk
      !
      !  Skip the XYZ comment line
      !
      read(inp,"(a)") buf
      write (out,"('Reading ',i5,' centres. XYZ comment is: ',a)") n_centres, trim(buf)
      !
      ! Read the atoms. Coordinates are expected to be in Angstroms,
      ! do that we'll do the conversion later.
      !
      read_atoms: do iat=1,n_centres
        !
        !  Input is split by lines; there is very litle error checking - so beware
        !
        read(inp,"(a)") buf
        read(buf,*,iostat=ios) centre_labels(iat), centres(:,iat), centre_types(iat), centre_params(:,iat)
        if (any(abs(centres(:,iat))>lat_vec/2)) then
          write (out,"('Atom at coordinates ',3f14.7,' is outside unit cell')") &
            centres(:,iat)
          stop 'bad geometry'
        end if
      end do read_atoms
      !
      !  If the structure was coming from an external file, close it
      !
      if (coord_file/=' ') then
        close(inp,status='keep')
      end if
    end subroutine read_xyz_structure
    !
    !  Report the structure, in case anybody cares ...
    !
    subroutine echo_xyz_structure
      integer(ik) :: iat     ! Atom number
      integer(ik) :: npar    ! Number of parameters
      real(rk)    :: qtot    ! Total charge of the structure
      !
      if (n_centres<=0) return 
      qtot = 0._rk
      !
      !  Report the coordinates and potential parameters. This routine is called before the
      !  conversion to atomic units, so the coordinates are still in Angstrom.
      !
      write (out,"(/t10,'Potential-carrying centres (could be atoms)'/)")
      write (out,"((1x,a10,2x,3(a15,1x),2x,a10))") &
          ' Label ', '    X    ', '    Y    ', '    Z    ', '   Type   ', &
          ' ----- ', '  -----  ', '  -----  ', '  -----  ', '  ------  '
      echo_atoms: do iat=1,n_centres
        select case (centre_types(iat))
          case default
            write (out,"('Atom ',i6,' has an unrecognized potential type ',a)") iat, trim(centre_types(iat))
            stop 'partition%echo_xyz_structure - Bad potential type'
          case ('CARB') ;  npar = 0 
          case ('LJ612') ; npar = 2
          case ('Q612') ; npar = 3
            qtot = qtot + centre_params(3,iat)
        end select
        write (out,"(1x,a10,2x,3(f15.7,1x),2x,a10,10(1x,g14.7))") &
               centre_labels(iat), centres(:,iat), centre_types(iat), centre_params(1:npar,iat)
      end do echo_atoms
      write (out,"()")
      write (out,"(' Total charge = ',f12.6)") qtot
      write (out,"()")
    end subroutine echo_xyz_structure
    !
    !  Calculate number of particles expected in this box
    !
    subroutine prepare_density_parameters
      real(rk) :: number_density  ! In particles/Bohr^3
      real(rk) :: p, f
      !
      number_density = 1._rk/(dens_fact * molar_volume)
      cell_volume    = product(lat_vec)
      particle_count = cell_volume * number_density
      call h2_tv_p(temperature,molar_volume,p)
      call h2_tp_f(temperature,p,f)
      !
      if (verbose>=0) then
        write (out,"()")
        write (out,"(' Average molar volume [m^3/mole] = ',g12.6)") molar_volume
        write (out,"('         Free gas pressure [bar] = ',g12.6)") p / bar
        write (out,"('    Free gas Helmholtz F [J/mol] = ',g12.6)") f
        write (out,"('Desired number density [Bohr^-3] = ',g12.6)") number_density
        write (out,"(' Simulation cell volume [Bohr^3] = ',g12.6)") cell_volume
        write (out,"('          Particle normalization = ',g12.6)") particle_count
        write (out,"()")
      end if
      !
      !  Report scaling parameters, useful for taking derivatives
      !
      write (out,"()")
      write (out,"('          Cell scaling parameter = ',f20.12)") cell_scale
      write (out,"('Change in molar volume [m^3/mol] = ',g20.12)") &
             (cell_scale**3-1._rk)*molar_volume
      write (out,"()")
    end subroutine prepare_density_parameters
    !
    !  Calculate gridded box extent and grid spacing
    !
    subroutine prepare_grid_parameters
      g_box(2,:) = 0.5_rk*lat_vec - lat_clip
      if (any(g_box(2,:)<=0)) then
        write (out,"(' Encountered negative upper dimension for the box: '/t8,3f14.6)") &
               g_box(2,:)
        stop 'bad box - too much clipping?'
      end if
      !
      g_box(1,:) = -g_box(2,:)
      g_box      =  g_box * cell_scale
      g_step(:)  = (g_box(2,:)-g_box(1,:))/n_points
      !
      if (verbose>=1) then
        write (out,"(//t12,'Grid-carrying part of the box (Bohr)'/)")
        write (out,"(t8,'   ',3(a14  ,2x))") '  Min  ', '  Max  ', '  Step  '
        write (out,"(t8,'X: ',3(f14.6,2x))") g_box(:,1), g_step(1)
        write (out,"(t8,'Y: ',3(f14.6,2x))") g_box(:,2), g_step(2)
        write (out,"(t8,'Z: ',3(f14.6,2x))") g_box(:,3), g_step(3)
        write (out,"()")
        write (out,"(t8,'(Scaling factor of ',f16.10,' was used.)')") cell_scale
        write (out,"()")
      end if
      !
    end subroutine prepare_grid_parameters
    !
    !  Calculate neighbour list and "grid" coordinates. We use periodic boundary
    !  conditions.
    !
    subroutine build_neighbour_list
      integer(ik) :: ix, iy, iz   ! X/Y/Z grid point
      real(rk)    :: x, y, z      ! Coordinates of a current grid point
      integer(ik) :: ipt
      !
      !$omp parallel do private(iz,z,iy,y,ix,x,ipt)
      loop_z: do iz=1,n_points(3)
        z = g_step(3)*(iz-0.5_rk*(1+n_points(3)))
        loop_y: do iy=1,n_points(2)
          y = g_step(2)*(iy-0.5_rk*(1+n_points(2)))
          loop_x: do ix=1,n_points(1)
            x = g_step(1)*(ix-0.5_rk*(1+n_points(1)))
            ipt = point_index(ix,iy,iz)
            xyz      (:,ipt) = (/ x, y, z /)
            neighbour(1,ipt) = point_index(ix-1,iy,  iz  )
            neighbour(2,ipt) = point_index(ix+1,iy,  iz  )
            neighbour(3,ipt) = point_index(ix,  iy-1,iz  )
            neighbour(4,ipt) = point_index(ix,  iy+1,iz  )
            neighbour(5,ipt) = point_index(ix,  iy,  iz-1)
            neighbour(6,ipt) = point_index(ix,  iy,  iz+1)
            if (verbose>=5) then
              write (out,"('IX= ',i4,' IY= ',i4,' IZ= ',i4,' ISEQ= ',i10,' XYZ = ',3f12.6)") &
                     ix, iy, iz, ipt, xyz(:,ipt)
              write (out,"((t15,6i12))") neighbour(:,ipt)
            end if
          end do loop_x
        end do loop_y
      end do loop_z
      !$omp end parallel do
      !
      contains 
        integer(ik) function point_index(ix,iy,iz)
          integer(ik), intent(in) :: ix, iy, iz
          integer(ik)             :: ixt, iyt, izt
          !
          ixt = ix ; iyt = iy ; izt = iz
          if (ixt<1) ixt = n_points(1)
          if (iyt<1) iyt = n_points(2)
          if (izt<1) izt = n_points(3)
          if (ixt>n_points(1)) ixt = 1
          if (iyt>n_points(2)) iyt = 1
          if (izt>n_points(3)) izt = 1
          point_index = ((izt-1)*n_points(2) + (iyt-1))*n_points(1) + ixt
        end function point_index
    end subroutine build_neighbour_list
    !
    !  Calculate non-zero elements of the kinetic energy matrix
    !
    subroutine build_kinetic_matrix
      !
      !  Negative off-diagonal displacements
      !
      kmat(1:5:2) = -0.5_rk / (mass*g_step(:)**2)
      !
      !  Positive off-diagonal displacements
      !
      kmat(2:6:2) = kmat(1:5:2)
      !
      !  Central position
      !
      kmat(0) = -sum(kmat(1:6))
      !
      if (verbose>=1) then
        write (out,"(//t12,'Kinetic matrix:'/)") 
        write (out,"('Center: ', g14.6)") kmat(0)
        write (out,"('     X: ',2g14.6)") kmat(1:2)
        write (out,"('     Y: ',2g14.6)") kmat(3:4)
        write (out,"('     Z: ',2g14.6)") kmat(5:6)
        write (out,"()")
      end if
    end subroutine build_kinetic_matrix
    !
    !  Flat-lined v12 interaction potential
    !
    function flatline_h2_v12(r) result(v)
      real(rk), intent(in) :: r ! Distance
      real(rk)             :: v ! Interaction potential
      !
      v = h2_v12(max(r,v12_cutpoint))
    end function flatline_h2_v12
    !
    !  Spherically capped v12 interaction potential. The constants are valid ONLY
    !  for our fit to the Diep-Johnson interaction potential.
    !
    function spherical_cap_h2_v12(r) result(v)
      real(rk), intent(in) :: r ! Distance
      real(rk)             :: v ! Interaction potential
      !
      real(rk), parameter  :: zp  = 3.024129479132871_rk ! r in Angstrom where v12 crosses zero
      real(rk), parameter  :: zpg = -312.042887128799_rk ! gradient (Kelvin/Angstrom) at the crossing point
      real(rk)             :: r_ang                      ! r in Angstrom
      real(rk)             :: av0
      !
      if (r>=zp/abohr) then
        v = h2_v12(r)
      else
        !
        !  Cap. Very few calls will end up here, so being efficient is not a priority
        !
        r_ang = r * abohr
        av0   = sqrt( (zp*(v12_sphmax+zp*zpg)**2)/(zpg*(2*v12_sphmax+zp*zpg)) )
        v     = v12_sphmax * ( sqrt(av0**2-r_ang**2) - sqrt(av0**2 - zp**2) ) &
                           / ( av0 - sqrt(av0**2-zp**2) )
        v     = v * k_Boltzmann
      end if
      ! write (out,"('+++ ',f14.7,1x,g15.9)") r, v
    end function spherical_cap_h2_v12
    !
    !  Filter for density smoothing
    !
    function wda_filter(r) result(w)
      real(rk), intent(in) :: r ! Distance from the origin
      real(rk)             :: w ! Smoothing weight
      real(rk)             :: efact
      !
      select case (wda_shape)
        case default
          write (out,"('partition%wda_filter is called for wda_shape = ',a)") trim(wda_shape)
          stop 'partition%wda_filter - invalid call'
        case ('Gaussian')
          efact = (r/wda_r0)**2
          if (efact<max_exp) then
            w = exp(-efact)/(pi**1.5_rk * wda_r0**3)
          else
            w = 0._rk
          end if
        case ('Fermi-Dirac')
          efact = (r-wda_r0)/wda_d
          if (efact<max_exp) then
            !
            !  Normalization factor is approximate - convolution generator 
            !  will fix it up later. The exact expression is:
            !    1/(-8*\pi*wda_d**3*Li_3(-Exp(wda_r0/wda_d))
            !  where Li_3 is the 3-rd order polylog function. Unfortunatly,
            !  numerical evaluation of polylog is way too ugly.
            !
            w = 1._rk/(8._rk*pi*wda_r0**2*wda_d*(1._rk + exp(efact)))
          else
            w = 0._rk
          end if
      end select
    end function wda_filter
    !
    !  Prepare data tables for v12 evaluation. This is the only routine
    !  which cares about the specific form of the potential or the WDA
    !  filters.
    !
    subroutine initialize_v12
      real(rk) :: wda_norm
      !
      !  v12 potential
      !
      select case (v12_form)
        case default
          write (out,"('Interaction potential ',a,' is not recognized')") trim(v12_form)
          stop 'partition%initialize_v12 - bad v12_form'
        case ('none')
          v12_uniform = 0._rk
        case ('diep-')
          if (verbose>=0) then
            write (out,"(/'Constructing V12 convolution:')")
            write (out,"( 'v12 flat-lined below ',f12.5,' Bohr = ',f12.5,' Angstrom')") &
                   v12_cutpoint, v12_cutpoint * abohr
          end if
          call convolution_initialize(cnv_v12,v12=flatline_h2_v12,eps=v12_epsilon, &
                  step=g_step,npoints=n_points,scale=v12_scale,verbose=verbose)
          v12_uniform = convolution_uniform_shift(cnv_v12)
        case ('diep-sph')
          if (verbose>=0) then
            write (out,"(/'Constructing V12 convolution:')")
            write (out,"( 'v12 capped by an arc segment. v12(0) is ',f12.8,' Hartree = ',f12.5,' Kelvin')") &
                   v12_sphmax * k_Boltzmann, v12_sphmax
          end if
          call convolution_initialize(cnv_v12,v12=spherical_cap_h2_v12,eps=v12_epsilon, &
                  step=g_step,npoints=n_points,scale=v12_scale,verbose=verbose)
          v12_uniform = convolution_uniform_shift(cnv_v12)
      end select
      !
      if (verbose>=0) then
        write (out,"(/' v12 factor for the uniform density = ',g14.7,' Hartree-Bohr^3')") &
               v12_uniform
      end if
      !
      if (v12_uniform>0._rk) then
        write (out,"(//1x,30('='),' E R R O R ',30('='))")
        write (out,"(' The functional is not variationally stable. The outcome of this calculation')")
        write (out,"(' is not meaningful, even if no variational collapse occurs for the current')")
        write (out,"(' starting guess.')")
        write (out,"( 1x,30('='),' E R R O R ',30('=')//)")
        call conditional_abort
      end if
      !
      !  Density smoothing filters
      !
      select case (wda_shape)
        case default
          write (out,"('WDA convolution shape ',a,' is not recognized')") trim(wda_shape)
          stop 'partition%initialize_v12 - bad wda_shape'
        case ('Delta')
          ! Do nothing
        case ('Gaussian','Fermi-Dirac')
          if (verbose>=0) then
            write (out,"(/'Constructing density weighting convolution:')")
            write (out,"( 'Weighting filter shape: ',a)") trim(wda_shape)
            write (out,"( 'Filter width: ',f12.5,' Bohr = ',f12.5,' Angstrom')") &
                   wda_r0, wda_r0 * abohr
            if (wda_shape/='Gaussian') &
            write (out,"( '  Edge width: ',f12.5,' Bohr = ',f12.5,' Angstrom')") &
                   wda_d, wda_d * abohr
          end if
          call convolution_initialize(cnv_wda,v12=wda_filter,eps=wda_epsilon, &
                  step=g_step,npoints=n_points,scale=1.0_rk,verbose=verbose)
          wda_norm = convolution_uniform_shift(cnv_wda)
          call convolution_rescale(cnv_wda,1.0_rk/wda_norm)
          if (verbose>=0) then
            write (out,"(/'Convolution filter norm before renormalization = ',g14.8)") wda_norm
          end if
      end select
    end subroutine initialize_v12
    !
    subroutine destroy_v12
      select case (v12_form)
        case default
          write (out,"('Interaction potential ',a,' is not recognized')") trim(v12_form)
          stop 'partition%initialize_v12 - bad v12_form'
        case ('none')
        case ('diep-','diep-sph')
          call convolution_destroy(cnv_v12)
      end select
      !
      if (wda_shape/='Delta') then
        call convolution_destroy(cnv_wda)
      end if
    end subroutine destroy_v12
    !
    !  Build the (diagonal) potential matrix
    !
    subroutine build_potential_matrix
      integer(ik) :: ipt
      !
      if (verbose>=1) then
        write (out,"(/' Adsorption potential is sampled over ',i0,' points for each volume element'/)") &
               oversample_vads**3
      end if
      !
      if (verbose>=0 .and. oversample_vads<=1) then
        write (out,"(/' WARNING: Potential oversampling has been disabled. If simulation volume includes'"// &
                    "/' WARNING: rapidly varying adsorption potentials, accuracy of the results may decrease'/)")
      end if
      !
      !$omp parallel do private(ipt)
      point_loop: do ipt=1,n_states
        v0dia(ipt) = pot_create(xyz(:,ipt))
      end do point_loop
      !$omp end parallel do
    end subroutine build_potential_matrix
    !
    !  Calculate and report accessible volume. We'll have three
    !  different versions of it:
    !   a) v <= 0
    !   b) v <= k*T
    !   c) v <= k*v_clip
    !
    subroutine report_accessible_volume
      integer(ik) :: n_acc        ! Number of accessible sites
      real(rk)    :: f_acc        ! Fraction of accessible sites
      real(rk)    :: e_min, e_max ! Minimum and maximum of the potential
      !
      write (out,"(/'Grid sites where adsorption potential is ...')")
      call count_accessible(0._rk,n_acc,f_acc)
      write (out,"(t5,'  ... less than zero: ',i10,' (',f10.5,' %)')") n_acc, 100._rk*f_acc
      call count_accessible(k_Boltzmann*temperature,n_acc,f_acc)
      write (out,"(t5,'    ... less than kT: ',i10,' (',f10.5,' %)')") n_acc, 100._rk*f_acc
      call count_accessible(k_Boltzmann*v_clip,n_acc,f_acc)
      write (out,"(t5,'... less than v_clip: ',i10,' (',f10.5,' %)')") n_acc, 100._rk*f_acc
      call vector_minmax(v0dia,e_min,e_max)
      write (out,"(t5,'The minimum of the potential is ',g12.5,' Hartree = ',g12.5,' Kelvin')") &
             e_min, e_min/k_Boltzmann
      !
      contains
        subroutine count_accessible(v,n,f)
          real(rk), intent(in)     :: v ! Max value of the potential included in the count
          integer(ik), intent(out) :: n ! Number of sites with potential <= v
          real(rk), intent(out)    :: f ! Fraction of sites with potential <= v
          !
          integer(ik) :: ipt
          !
          n = 0
          !$omp parallel do private(ipt) reduction(+:n)
          count_sites: do ipt=1,n_states
            if (v0dia(ipt)<=v) n = n+1
          end do count_sites
          !$omp end parallel do
          f = real(n,kind=rk)/n_states
        end subroutine count_accessible
    end subroutine report_accessible_volume
    !
    !  Evaluate LDA exchange-correlation energy density and its derivative wrt density
    !
    subroutine evaluate_exc(n,nave,p,exc,ddn,ddnave)
      real(rk), intent(in)  :: n      ! Density at a grid point, in bohr^-3. ZERO IS NOT ALLOWED!
      real(rk), intent(in)  :: nave   ! Smoothed density at a grid point, in bohr^-3. ZERO IS NOT ALLOWED!
      real(rk), intent(out) :: p      ! Effective pressure at a grid point, in Pascal
      real(rk), intent(out) :: exc    ! Exchange-correlation energy density, in Hartree/particle
      real(rk), intent(out) :: ddn    ! Derivative of exc with respect to density, in Hartree/bohr^3
      real(rk), intent(out) :: ddnave ! Derivative of exc with respect to smoothed density, in Hartree/bohr^3
      !
      real(rk) :: vmol, dvmoldn ! Molar volume, in m^3/mole, and its derivative wrt particle density
      real(rk) :: dpdv          ! d P/d V, at constant T, in mole-Pa/m^3
      real(rk) :: f, dfdp       ! Helmholtz free energy, in J/mole and it's pressure derivative at constant T
      real(rk) :: f_au, dfdp_au ! ditto, in atomic units [the derivative is still wrt Pascals]
      real(rk) :: dfdn_au       ! Derivative of the Helmholtz energy wrt particle density, at constant T
      real(rk) :: q_particle    ! Kinetic partition function per particle - needed for size-extensivity
      !
      vmol    =  1._rk/(dens_fact*nave)
      dvmoldn = -1._rk/(dens_fact*nave**2)
      call h2_tv_p(temperature,vmol,p,dpdv)    ! Pressure and derivative, from experimental data
      call h2_tp_f(temperature,p,f,dfdp)       ! Free energy and derivative, ditto
      f_au    = f / (Hartree * N_Avogadro)     ! Free energy in atomic units
      dfdp_au = dfdp / (Hartree * N_Avogadro)  ! derivative in mixed units: Hartree/(particle-Pascal)
      dfdn_au = dfdp_au * dpdv * dvmoldn       ! Use chain rule to differentiate free energy wrt density
      !
!     q_particle = q_free / (n*cell_volume)    ! Indistinguishability correction - local density version
      q_particle = q_free / (nave*cell_volume) ! Indistinguishability correction - weighted density version
      !
      exc = f_au + (1_rk/beta)*log(q_particle) ! Subtract contribution due to partition function of 
                                               ! free ideal gas of INdistinguishable particles
      exc = exc - n*0.5_rk*v12_uniform         ! Add correction to the LDA form due to v12 potential - physical density form
      ddn    = -0.5_rk*v12_uniform             ! Derivative of vxc with respect to density - physical density form
!     ddn    = ddn - 1._rk/(beta*n)            ! ... and the indistinguishability correction - version with dependence on n
      ddnave = dfdn_au                         ! Derivative of vxc with respect to smoothed density
      ddnave = ddnave - 1._rk/(beta*nave)      ! ... and the indistinguishability correction - version with dependence on nave
      !
      if (verbose>=5) then
        write (out,"('EXC: n= ',g12.6,' nave= ',g12.6,' exc= ',g12.6,' ddn= ',g12.6,' ddnave= ',g12.6)") &
               n, nave, exc, ddn, ddnave
      end if
    end subroutine evaluate_exc
    !
    !  We calculate two quantities here:
    !   eps_xc () is the exchange-correlation energy density, at the _smoothed_ density
    !   vxc_n  () is the derivative of eps_xc wrt density, times density itself
    !   vxc_nave() is the derivative of eps_xc wrt the smoothed density, times the real
    !             density. It still needs to be convoluted with the wda smoothing function
    !             to form the vxc. (Of course, for pure LDA convolution is a delta function).
    !
    subroutine evaluate_vxc_components(n_wda,n,vxc_n,vxc_nave)
      real(rk), intent(in)  :: n_wda   (:)  ! Smoothed density, at which LDA functional is evaluated
      real(rk), intent(in)  :: n       (:)  ! Actual density, needed to construct the vxc 
      real(rk), intent(out) :: vxc_n   (:)  ! Derivative wrt density times the actual density
      real(rk), intent(out) :: vxc_nave(:)  ! Derivative wrt smoothed density times the actual density
      !
      integer(ik) :: ipt
      real(rk)    :: ncut_wda    ! Weighted density with cut-off applied
      real(rk)    :: ncut        ! Density with cut-off applied
      real(rk)    :: p           ! Effective LDA "pressure", in Pascal
      real(rk)    :: exc         ! LDA exchange-correlation energy
      real(rk)    :: ddn, ddnave ! ... and its derivatives wrt density and smoothed density
      real(rk)    :: pl, ph      ! Lowest and highest pressures
      !
      pl = 1e20 ; ph = -1e20 ;
      !$omp parallel do private(ipt,ncut,ncut_wda,p,exc,ddn,ddnave)  &
      !$omp& reduction(min:pl) reduction(max:ph)
      point_loop: do ipt=1,n_states
        ncut     = max(dens_cutoff,n    (ipt))
        ncut_wda = max(dens_cutoff,n_wda(ipt))
        call evaluate_exc(ncut,ncut_wda,p,exc,ddn,ddnave)
        pl  = min(pl,p) ; ph  = max(ph,p)
        eps_xc  (ipt) = exc
        vxc_n   (ipt) = ddn*ncut
        vxc_nave(ipt) = ddnave*ncut
      end do point_loop
      !$omp end parallel do
      p_min = pl ; p_max = ph
    end subroutine evaluate_vxc_components
    !
    !  For density given in dens(), construct the exchange-correlation
    !  potential in vxc(). We have two principal forks, depending on
    !  whether the weighted-density approximation is invoked or not.
    !  As a side-effect, we'll also leave the exchange-correlation
    !  energy density in eps_xc(), so that the post-SCF correction
    !  can be constructed later.
    !
    subroutine construct_exchange_correlation_potential
      if (wda_shape=='Delta') then
        !
        !  In the absence of the WDA averaging, actual and smoothed densities
        !  coincide, and the final convolution step is not needed
        !
        call evaluate_vxc_components(dens,dens,vxc_n,vxc_nave)
        call vector_copy(dst=vxc,src=vxc_n)
        call vector_add_to(dst=vxc,src=vxc_nave)
        call vector_add_to(dst=vxc,src=eps_xc)
      else
        !
        !  Weighted-density expression for vxc has to be evaluated.
        !
        !  Step 1: Calculate averaged density.
        !
        call convolution_evaluate(cnv_wda,dens,smooth_dens)
        !
        !  Step 2: Evaluate eps_xc and its derivative times instantaneous
        !          density at the smoothed density grid. The results are
        !          left in eps_xc() and vxc_pad()
        !
        call evaluate_vxc_components(smooth_dens,dens,vxc_n,vxc_nave)
        !
        !  Step 3: Convolution of the derivative term with the weighting
        !          function.
        !
        call convolution_evaluate(cnv_wda,vxc_nave,vxc)
        !
        !  Step 4: Add exchange-correlation energy density to vxc(), and
        !          we are done here.
        !
        call vector_add_to(dst=vxc,src=vxc_n)
        call vector_add_to(dst=vxc,src=eps_xc)
      end if
    end subroutine construct_exchange_correlation_potential
    !
    !  Add mean-field (v12) term to the effective potential.
    !  We will also calculate and store the mean-field energy.
    !
    subroutine add_v12_potential
      integer(ik) :: ipt
      real(rk)    :: ncut  ! Density with cut-off applied
      real(rk)    :: wgt   ! Volume element of the grid
      !
      energy_v12 = 0
      if (v12_form=='none') return ! No mean-field potential - nothing to do
      !
      !  Evaluate mean-field potential using previously computed
      !  Fourier image of the 2-particle interaction potential.
      !
      call convolution_evaluate(cnv_v12,dens,v12)
      !
      wgt = product(g_step)
      !$omp parallel do private(ipt,ncut) reduction(+:energy_v12)
      point_loop: do ipt=1,n_states
        ncut = max(dens_cutoff,dens(ipt))
        vdia(ipt)  = vdia(ipt) - beta*v12(ipt)
        energy_v12 = energy_v12 + ncut*v12(ipt)
      end do point_loop
      !$omp end parallel do
      !
      !  Total mean-field energy includes the second volume element
      !  Additionally, pro-rate it to the number of particles
      !
      energy_v12 = 0.5_rk * energy_v12 * wgt / particle_count
      ! write (out,"('energy_v12 = ',g14.7)") energy_v12
    end subroutine add_v12_potential
    !
    !  Evaluate the post-SCF contribution to the DFT energy, which has to be
    !  added to the sum of the eigenvalues to get the total energy. In both the
    !  LDA and WDA approximations, this contribution is given simply by the 
    !  density-weighted integral of the difference between exchange-correlation 
    !  energy density and exchange-correlation potential.
    !
    subroutine dft_energy_correction(edft)
      real(rk), intent(out) :: edft ! Post-SCF contribution to the energy
      !
      integer(ik) :: ipt
      real(rk)    :: ncut         ! Density with cut-off applied
      real(rk)    :: n_sum, d_sum ! Sum of the sanitized and raw densities
      !
      edft  = 0      ! In principle, this should be =energy_fxc-energy_vxc,
                     ! but error accumulation properties are better this way
      n_sum = 0      ! Sum of the sanitized densities
      d_sum = 0      ! Sum of the true densities
      energy_fxc = 0 ! Exchange-correlation energy
      energy_vxc = 0 ! Exchange-correlation potential "energy"
      !$omp parallel do private(ipt,ncut) reduction(+:edft,n_sum,d_sum,energy_fxc,energy_vxc)
      point_loop: do ipt=1,n_states
        ncut  = max(dens_cutoff,dens(ipt))
        n_sum = n_sum + ncut
        d_sum = d_sum + dens(ipt)
        edft  = edft + ncut*(eps_xc(ipt)-vxc(ipt))  ! Grid weight factor is included below
        energy_fxc = energy_fxc + ncut*eps_xc(ipt)
        energy_vxc = energy_vxc + ncut*vxc(ipt)
      end do point_loop
      !$omp end parallel do
      !
      !  The correction is now for a system of particle_count particles
      !  rescale it to per-particle and include grid weight factor
      !
      !write (out,"(' d_sum = ',g15.8,' n_sum = ',g15.8,' correction = 1 + ',g15.8)") &
      !       d_sum, n_sum, d_sum/n_sum - 1.0_rk
      edft = edft*(d_sum/n_sum)*product(g_step)/particle_count
      energy_fxc = energy_fxc*(d_sum/n_sum)*product(g_step)/particle_count
      energy_vxc = energy_vxc*(d_sum/n_sum)*product(g_step)/particle_count
    end subroutine dft_energy_correction
    !
    !  Renormalize total density to match the particle count in the simulation box
    !
    subroutine rescale_density(shout)
      logical, intent(in), optional :: shout ! Report total particle counts
      !
      integer(ik) :: ipt
      real(rk)    :: rho_total ! Sum of the densities - must match q_potential
      real(rk)    :: scl       ! Scale factor needed to give correct particle count
      !
      !  Start by sanitizing the density, removing negative densities.
      !
      rho_total = 0
      !$omp parallel do private(ipt) reduction(+:rho_total)
      sanitize_density: do ipt=1,n_states
        if (dens(ipt)<0) dens(ipt) = 0
        rho_total = rho_total + dens(ipt)
      end do sanitize_density
      !$omp end parallel do
      !
      if (verbose>=2) then
        write (out,"('Sanitized density sums up to ',g16.8,' (expected ',g16.8,')')") &
               rho_total, q_potential
      end if
      !
      scl = particle_count / (rho_total * product(g_step))
      !
      !$omp parallel do private(ipt)
      renormalize_density: do ipt=1,n_states
        dens(ipt) = scl * dens(ipt)
      end do renormalize_density
      !
      !  Output relevant for guess densities, but not for anything else
      !
      if (present(shout)) then
        if (shout) then
          write (out,"(/' Guess density was normalized to ',g14.5,' particles.')") rho_total * product(g_step)
          write (out,"( '                 renormalized to ',g14.5,' particles.'/)") particle_count
        end if
      end if
    end subroutine rescale_density
    !
    !  Calculate hydrogen gas pressure matching given free energy.
    !  Because this operation is of negligible cost overall, we'll 
    !  use simple bisection.
    !
    subroutine external_pressure(f,p)
      real(rk), intent(in)  :: f  ! Desired free energy
      real(rk), intent(out) :: p  ! Pressure giving this free energy
      !
      real(rk) :: p_low, p_high, p_mid, f_mid
      !
      p_low  =  0.001_rk * bar ! 0.001 bar - No point in looking below this!
      p_high = 5e4 * bar       ! 50 kbar - No point in looking aboce this!
      !
      !  There is no point in looking for a solution more accurate than
      !  a millibar!
      !
      bisect: do while( p_high-p_low > max(1e-3_rk*bar,100._rk*spacing(p_low+p_high)) )
        p_mid = 0.5_rk*(p_high+p_low)
        call h2_tp_f(temperature,p_mid,f_mid)
        if (f_mid<f) then
          p_low  = p_mid
        else
          p_high = p_mid
        end if
      end do bisect
      !
      p = 0.5_rk * (p_low + p_high)
    end subroutine external_pressure
    !
    !  Analytical result for eigenvalues of a cyclic matrix. The eigenvalues
    !  are not sorted.
    !
    subroutine kinetic_eigenvalues_1D(a,b,npts,eps)
      real(rk), intent(in)    :: a      ! Diagonal matrix elements
      real(rk), intent(in)    :: b      ! Coupling matrix elements
      integer(ik), intent(in) :: npts   ! Number of grid points
      real(rk), intent(out)   :: eps(:) ! Eigenvalues - npts of them
      !
      integer(ik) :: k
      !
      eigencircle: do k=0,npts-1
        eps(1+k) = a + 2*b*cos((twopi*k)/npts)
      end do eigencircle
    end subroutine kinetic_eigenvalues_1D
    !
    !  Evaluate kinetic partition function for an inifinitely dense 
    !  1D cyclic cluster grid. The limit happens to be given by:
    !
    !    t3(0,Exp(-0.5*(beta/mass)*(2*pi/L)**2))
    !
    !  where t3 is the elliptic theta function of the third kind,
    !  and L is the box extent. We evaluate the t3 using a rapidly
    !  converging special-case series expansion, valid for the zero 
    !  first argument.
    !
    function cyclic_sum_limit(arg) result(s)
      real(rk), intent(in) :: arg ! the value of 0.5*(beta/mass)*(2*pi/L)**2
      real(rk)             :: s
      !
      real(rk)    :: x    ! exp(-arg)
      real(rk)    :: term ! exp(-arg)
      integer(ik) :: n    ! sum iterator
      !
      if (arg>=max_exp) then
        s = 1._rk ! The low-temperature limit
        return
      end if
      !
      if (arg<=0._rk) then
        write (out,"('cyclic_sum_limit: asked for a limit of a divergent series. arg = ',g20.10)") arg
        stop 'partition-in-small-box%cyclic_sum_limit - argument domain error'
      end if
      !
      x = exp(-arg)
      s = 0
      t3_series_sum: do n=1_ik,huge(1_ik)
        term = x**(n**2)
        s    = s + term
        if (term<spacing(s)) exit t3_series_sum
      end do t3_series_sum
      !
      s = 1._rk + 2._rk*s
      !
    end function cyclic_sum_limit
    !
    !  Evaluate partition function for the empty box, using analytical
    !  results for the eigenspectrum
    !
    subroutine calculate_q_free
      real(rk)    :: eps_x(n_points(1))  ! Kinetic energy eigenvalues along the X axis
      real(rk)    :: eps_y(n_points(2))  ! ... Y axis
      real(rk)    :: eps_z(n_points(3))  ! ... Z axis 
      real(rk)    :: eps
      integer(ik) :: kx, ky, kz
      !
      call TimerStart('Empty box Q')
      !
      !  The continuous distribution limit is -very- easy
      !
      q_free_integral = product(g_box(2,:)-g_box(1,:))*(mass/(twopi*beta))**1.5_rk
      !
      !  The cyclic-cluster limit for an infinitely dense grid is also
      !  fairly easy:
      !
      q_free_limit = cyclic_sum_limit((beta/(2*mass))*(twopi/(g_box(2,1)-g_box(1,1)))**2) &
                   * cyclic_sum_limit((beta/(2*mass))*(twopi/(g_box(2,2)-g_box(1,2)))**2) &
                   * cyclic_sum_limit((beta/(2*mass))*(twopi/(g_box(2,3)-g_box(1,3)))**2)
      !
      !  The grid-based result is a little harder, as the eigenspectrum is
      !  discrete, and must be calculated. First, calculate eigenvalues for
      !  each orthogonal spatial direction
      !
      call kinetic_eigenvalues_1D(-2*kmat(1),kmat(1),n_points(1),eps_x)
      call kinetic_eigenvalues_1D(-2*kmat(3),kmat(3),n_points(2),eps_y)
      call kinetic_eigenvalues_1D(-2*kmat(5),kmat(5),n_points(3),eps_z)
      !
      !  Combine all grid-based eigenvalues in the partition function
      !
      q_free = 0 
      !$omp parallel do reduction(+:q_free) private (kx,ky,kz,eps)
      eigenloop_z: do kz=1,n_points(3)
        eigenloop_y: do ky=1,n_points(2)
          eigenloop_x: do kx=1,n_points(1)
            eps    = eps_x(kx) + eps_y(ky) + eps_z(kz)
            q_free = q_free + exp(-beta*eps)
            ! write (out,"(1x,'# ',3i5,f20.10)") kx, ky, kz, eps
          end do eigenloop_x
        end do eigenloop_y
      end do eigenloop_z
      !$omp end parallel do
      call TimerStop('Empty box Q')
    end subroutine calculate_q_free
    !
    !  Save final density for visualization in OpenDX
    !
    subroutine dump_dx_field(file,field,comment)
      character(len=*), intent(in) :: file     ! Name of OpenDX file
      real(rk), intent(in)         :: field(:) ! Data field to dump
      character(len=*), intent(in) :: comment  ! Descriptive comment
      !
      !  Create OpenDX file - the header was shamelessly stolen from 
      !  OpenDX dump, so it should (mostly) "just work"
      !
      open (unit_dens,form='formatted',status='replace',file=file)
      write (unit_dens,"(" // &
         " 'object 1 class array type float rank 0 items ',i12,' data 0'" // &
         "/'attribute ""dep"" string ""positions""'"                      // &
         "/'object 2 class gridpositions counts ',3(1x,i5)"               // &
         "/' origin  ',3(1x,f14.7)"                                       // &
         "/' delta   ',3(1x,f14.7)"                                       // &
         "/' delta   ',3(1x,f14.7)"                                       // &
         "/' delta   ',3(1x,f14.7)"                                       // &
         "/'attribute ""dep"" string ""positions""'"                      // &
         "/'object 3 class gridconnections counts ',3(1x,i5)"             // &
         "/'attribute ""element type"" string ""cubes""'"                 // &
         "/'attribute ""dep"" string ""connections""'"                    // &
         "/'attribute ""ref"" string ""positions""'"                      // &
         "/'object ""default"" class field'"                              // &
         "/'component ""data"" value 1'"                                  // &
         "/'component ""positions"" value 2'"                             // &
         "/'component ""connections"" value 3'"                           // &
         "/'attribute ""name"" string ""default""'"                       // &
         "/'end')") n_states, n_points(3:1:-1), &
                    g_box(1,3:1:-1) + 0.5_rk * g_step(3:1:-1), &
                    g_step(3), 0._rk,     0._rk,     &
                    0._rk,     g_step(2), 0._rk,     &
                    0._rk,     0._rk,     g_step(1), &
                    n_points(3:1:-1)
       write (unit_dens,"(6(g12.6,1x))") field
       close (unit_dens)
    end subroutine dump_dx_field
    !
    !  Load starting density produced by dump_dx_field. Only minimal 
    !  checking is done, so beware.
    !
    subroutine load_dx_field(file,field)
      character(len=*), intent(in) :: file     ! Name of OpenDX file
      real(rk), intent(out)        :: field(:) ! Data field to fill
      !
      character(len=256) :: linebuf
      character(len=256) :: errmsg
      integer(ik)        :: ios, counter, pos, ngrid(3)
      !
      read_block: do
        errmsg = 'Error opening '//trim(file)//' for reading'
        linebuf = ' '
        counter = 0
        open (unit_dens,form='formatted',status='old',action='read',file=file,iostat=ios)
        errmsg = 'Error reading '//trim(file)//' while skipping header records'
        skip_header: do
          counter = counter + 1
          read (unit_dens,"(a)",iostat=ios) linebuf
          if (ios/=0) exit read_block
          if (linebuf=='end') exit skip_header  ! End of header
          !
          !  Minimal sanity test - grid counts must match
          !
          pos = index(linebuf,' gridpositions counts ')
          if (pos>0) then
            pos = pos + len(' gridpositions counts ')
            read (linebuf(pos:),*,iostat=ios) ngrid(3:1:-1)
            if (ios/=0) then
              errmsg = 'Error parsing grid dimensions in file '//trim(file)
              exit read_block
            end if
            if (any(ngrid/=n_points)) then
              write (errmsg,"('Grid in file ',a,' (',3i5,') does not match active grid (',3i5,')')") &
                     trim(file), ngrid, n_points
              exit read_block
            end if
          end if
        end do skip_header
        counter = counter + 1
        errmsg = 'Error reading data field from file '//trim(file)
        read  (unit_dens,"(6(g12.6,1x))",iostat=ios) field
        if (ios/=0) exit read_block
        errmsg = 'Error closing file '//trim(file)
        close (unit_dens,iostat=ios)
        if (ios/=0) exit read_block
        return
      end do read_block
      !
      !  We only come here if there was an error
      ! 
      write (out,"(/1x,a)") trim(errmsg)
      write (out,"(' Last input line was:'/1x,a)") trim(linebuf)
      write (out,"(' iostat       = ',i8)") ios
      write (out,"(' line counter = ',i8)") counter
      stop 'load_dx_field'
    end subroutine load_dx_field
    !
    !  Generate image of the structure, suitable for plotting with OpenDX
    !
    subroutine plotStructure
      integer(ik)              :: iatom, jatom, natoms, ibond, nbonds,k,j,i,n
      real(rk)                 :: rij 
      real(rk),allocatable     :: atoms (:,:)
      integer(ik), allocatable :: bonds(:,:)
      !
      ! Count number of atoms in the replicated structure
      !
      natoms = n_centres*product(2*plot_cells+1)
      !
      allocate (atoms(3,natoms),bonds(2,(natoms*(natoms-1))/2))
      !
      ! Explicitly replicate atoms
      !
      iatom = 0
      !
      sum_a: do k = -plot_cells(1),plot_cells(1)       ! Replicate the cell in x direction
         sum_b: do j = -plot_cells(2),plot_cells(2)    ! Replicate the cell in y direction
            sum_c: do n = -plot_cells(3),plot_cells(3) ! replicate the cell in z direction
               do i = 1,n_centres                      ! sum over all atoms
                  iatom = iatom+1
                  atoms(:,iatom) = centres(:,i) + lat_vec * (/k,j,n/)
               end do
            end do sum_c
         end do sum_b
      end do sum_a

      if (iatom/=natoms) then
         write (out,"(' iatom = ',i5,' natoms = ',i5,'. Oops.')") iatom, natoms
         stop 'atom count error'
      end if
      !
      !  Count and record bonds
      !
      ibond = 0
      bond_scan_i: do iatom=1,natoms
        bond_scan_j: do jatom=iatom+1,natoms
          rij = sqrt(sum( (atoms(:,iatom)-atoms(:,jatom))**2 ))
          if (rij>bond_cc*1.25_rk) cycle bond_scan_j
          ibond = ibond + 1
          bonds(1,ibond) = iatom
          bonds(2,ibond) = jatom
        end do bond_scan_j
      end do bond_scan_i
      nbonds = ibond
      !
      !  Write OpenDX output file
      !
      open (unit_mol,form='formatted',status='replace',file=image_file)
      write (unit_mol,"('object 1 class array type float rank 0 "// &
                 "items ',i6,' data follows')") natoms
      write (unit_mol,"(8(1x,i5))") (6,iatom=1,natoms)
      write (unit_mol,"('attribute ""dep"" string ""positions""')")
      !
      !  Our data fields go out to DX in the reverse order - so we'll
      !  need to swap the XYZ coordinates as well.
      !
      write (unit_mol,"('object 2 class array type float rank 1 shape 3"// &
                 " items ',i6,' data follows')") natoms
      write (unit_mol,"(3(1x,f15.8))") atoms(3:1:-1,1:natoms)
      write (unit_mol,"('attribute ""dep"" string ""positions""')")
      !
      write (unit_mol,"('object 3 class array type int rank 1 shape 2"// &
                 " items ',i6,' data follows')") nbonds
      write (unit_mol,"(2(1x,i7))") bonds(:,1:nbonds)-1
      write (unit_mol,"('attribute ""element type"" string ""lines""')")
      write (unit_mol,"('attribute ""ref"" string ""positions""')")
      !
      write (unit_mol,"('object ""molecule"" class field')")
      write (unit_mol,"('component ""data"" value 1')")
      write (unit_mol,"('component ""positions"" value 2')")
      write (unit_mol,"('component ""connections"" value 3')")
      write (unit_mol,"('attribute ""name"" string ""graphene""')")
      !
      close (unit_mol)
      !
      deallocate (atoms,bonds)
      return
    end subroutine plotStructure
    !
    !  Dump density and potential, useful (?) for debugging convergence
    !
    subroutine dump_density_and_potential(iter)
      integer(ik), intent(in) :: iter ! Iteration number, to make locating the right dump easier
      character (len=200)     :: file ! File name to use for dumps
      !
      write (file,"('iter-',i2.2,'-',i4.4,'-dens.dx')") debug_sequence, iter
      call dump_dx_field(file,dens,'Final density, particle/Bohr^3')
      write (file,"('iter-',i2.2,'-',i4.4,'-veff.dx')") debug_sequence, iter
      call dump_dx_field(file,vdia,'Effective potential, Kelvin')
    end subroutine dump_density_and_potential
    !
    !  Comprehensive ASCII dump of the results. For all grid points, dump
    !  their coordinates, final SCF density, final local pressure, external 
    !  potential, exchange-correlation energy density, and exchange-
    !  correlation potential.
    !
    !  This is a very simple-minded routine, which assumes that the external 
    !  potential is in the form pot_create() returns it
    !
    subroutine dump_ascii_density
      integer(ik) :: ipt          ! Grid point
      real(rk)    :: vext         ! External potential at a point
      !
      call TimerStart('Final ASCII dump')
      open (unit_dens,form='formatted',status='replace',file=ascii_file)
      write (unit_dens,"('#')")
      write (unit_dens,"('# Values and units:')")
      write (unit_dens,"('#   X, Y, Z - Grid coordinate [Bohr]')")
      write (unit_dens,"('#   N       - Final density [particles/Bohr^3]')")
      write (unit_dens,"('#   N_WDA   - Smoothed final density [particles/Bohr^3]')")
      write (unit_dens,"('#   VEXT    - External potential [Hartree/particle]')")
      write (unit_dens,"('#   EPS     - Exchange-correlation energy density [Hartree/particle]')")
      write (unit_dens,"('#   VXC     - Exchange-correlation potential [Hartree/particle]')")
      write (unit_dens,"('#  ')")
      write (unit_dens,"('# Grid extent:')")
      write (unit_dens,"('#   X: ',f14.7,' - ',f14.7,' ',i6,' pts')") g_box(:,1), n_points(1)
      write (unit_dens,"('#   Y: ',f14.7,' - ',f14.7,' ',i6,' pts')") g_box(:,2), n_points(2)
      write (unit_dens,"('#   Z: ',f14.7,' - ',f14.7,' ',i6,' pts')") g_box(:,3), n_points(3)
      write (unit_dens,"('# Volume of each grid point: ',e14.7,' bohr^3')") product(g_step)
      write (unit_dens,"('#  ')")
      write (unit_dens,"('#',a12,2(1x,a12),6(1x,a14))") &
             ' X ', ' Y ', ' Z ', ' N ', ' N_WDA ', ' VEXT ', ' EPS ', ' VXC '
      point_loop: do ipt=1,n_states
        vext = pot_create(xyz(:,ipt))
        write (unit_dens,"(3(1x,f12.6),6(1x,e14.7))") &
               xyz(:,ipt), dens(ipt), smooth_dens(ipt), vext, eps_xc(ipt), vxc(ipt)
      end do point_loop
      close (unit_dens)
      call TimerStop('Final ASCII dump')
    end subroutine dump_ascii_density
    !
    !  Fold and print the comment line.
    !
    subroutine print_comment
      integer(ik)            :: pos, end, max_split, split
      integer(ik), parameter :: max_str = 78, max_slack = 10
      !
      if (comment==' ') return
      write (out,"(/1x,78('='))")
      pos = 1
      end = len_trim(comment)
      cut_comment: do while(pos<=end)
        max_split = min(end,pos+max_str)
        if (max_split-pos+1<=max_str) then
          split = max_split
        else
          split = pos + index(comment(pos:max_split),' ',back=.true.) - 1
          if (split<=max_split-max_slack) split = max_split
        end if
        write (out,"(1x,a)") trim(comment(pos:split))
        pos = split + 1
      end do cut_comment
      write (out,"( 1x,78('=')/)")
    end subroutine print_comment
    !
    !  Something happened to make calculation meaningless. Abort unless user
    !  insists that we should continue
    !
    subroutine conditional_abort
      if (let_it_be) then
        write (out,"(/' You have specified ''LET_IT_BE = .TRUE.'' in the input. The calculation will')")
        write (out,"( ' continue, even though the results are likely to be meaningless.'/)")
      else
        write (out,"(/' Aborting the calculation. Although the results will likely be meaningless,')")
        write (out,"( ' you can choose to continue the calculation by adding ''LET_IT_BE = .TRUE.''')")
        write (out,"( ' to the input block.'/)")
        stop 'conditional abort'
      end if
    end subroutine conditional_abort
    !
    !  Initialize parameters for the SCF loop
    !
    subroutine scf_initialize
      real(rk) :: grid_err
      !
      !  Increment the number of calls to scf_loop, to make convergence dumps
      !  easier to locate.
      !
      debug_sequence = debug_sequence + 1
      !
      !  Neighbour list and kinetic matrix
      !
      call TimerStart('Kinetic matrix')
      call build_neighbour_list  ! Also fills xyz array, to make sure it is 
                                 ! consistent between kinetic and potential
                                 ! matrices
      call build_kinetic_matrix
      call TimerStop('Kinetic matrix')
      write (out,"('Finished with the kinetic matrix')")
      call flush(out)
      !
      !  External potential matrix. We use Hartrees internally - this is what
      !  build_potential_matrix generates. However, the external potential is
      !  is stored in Kelvin (this makes visualization easier), so there is a
      !  bit of constant jiggery below.
      !
      call TimerStart('Potential matrix')
      if (vads_read_file/=' ') then
        !
        !  External file is expected to be in the format produced for
        !  vads_file= below
        !
        if (verbose>=1) then
          write (out,"(/' Loading adsorption potential from ',a)") trim(vads_read_file)
        end if
        call load_dx_field(vads_read_file,v0dia)
        call vector_scale(dst=v0dia,a=0._rk,b=k_Boltzmann)
      else
        call build_potential_matrix
      end if
      call vector_copy(dst=vads,src=v0dia)
      write (out,"('Finished with the potential matrix')")
      call flush(out)
      call TimerStop('Potential matrix')
      !
      if (verbose>=5) then
        write (out,"(/' v0dia = '/(5(1x,g18.12))/)") v0dia
      end if
      !
      !  Save adsorption potential for plotting, if needed. Plot should
      !  be in Kelvin, so some rescaling will be done before and after.
      !
      if (vads_file/=' ') then
        call vector_scale(dst=v0dia,a=0._rk,b=1._rk/k_Boltzmann)
        call dump_dx_field(vads_file,v0dia,'Adsorption potential, Kelvin')
        call vector_scale(dst=v0dia,a=0._rk,b=k_Boltzmann)
      end if
      !
      !  Calculate various versions of accessible volume
      !
      call report_accessible_volume
      !
      !  Temperature stuff
      !
      call TimerStart('Temperature')
      write (out,"()")
      write (out,"('Temperature [Kelvin] = ',g14.7)") temperature
      write (out,"('   beta [Hartree^-1] = ',g14.7)") beta
      write (out,"('  Vpot clip [Kelvin] = ',g14.7)") v_clip
      write (out,"(' Vpot clip [Hartree] = ',g14.7)") v_clip*k_Boltzmann
      !
      !  Calculate analytical result for Q(Free), to test for grid
      !  convergence. We depend on having kinetic matrix avalable.
      !
      call calculate_q_free
      write (out,"(/'Analytical results for 1-particle partition function in empty box')")
      write (out,"(t5,'     Continuos limit = ',g16.8)") q_free_integral
      write (out,"(t5,'Converged grid limit = ',g16.8)") q_free_limit
      write (out,"(t5,'        Current grid = ',g16.8)") q_free
      write (out,"()")
      grid_err = q_free/q_free_limit
      if ( verbose>=1 .and. (grid_err<0.5_rk .or. grid_err>2.0_rk) ) then
        write (out,"(t5,'WARNING: A significant difference between the grid limit and the current'"// &
                   "/t5,'WARNING: grid was detected, indicating that the current grid is approaching'"// &
                   "/t5,'WARNING: saturation. Note that for mildly deficient grids, it is not unusual'"// &
                   "/t5,'WARNING: to obtain one-particle partition functions in excess of the converged'"// &
                   "/t5,'WARNING: result. Proceed with extreme caution!'//)")
      end if
      call flush(out)
      !
      !  Apply temperature to kinetic and potential matrices, and add the
      !  diagonal part of the kinetic matrix to the potential.
      !
      kmat  = -beta*kmat
      !
      !  We'll clip later, once the XC potential has been added. 
      !
      call vector_scale(dst=v0dia,a=kmat(0),b=-beta) ! dst = a + b*dst
      !
      if (verbose>=5) then
        write (out,"(/'-beta*kmat ='/7(1x,g16.8))") kmat
        write (out,"(/'-beta*(v0dia+kmat(0) = '/(5(1x,g18.12))/)") v0dia
      end if
      !
      if (verbose>=2) then
        write (out,"(/'Diagonal of (-beta*H):')")
        write (out,"(t5,'Smallest value = ',g14.7)") minval(v0dia)
        write (out,"(t5,' Largest value = ',g14.7)") maxval(v0dia)
        write (out,"()")
      end if
      call TimerStop('Temperature')
      !
      call print_comment
      !
      call TimerReport
      write (out,"(/'Ready to begin iterations'/)")
      call flush(out)
    end subroutine scf_initialize
    !
    !  SCF mixing and exchange-correlation potential
    !
    subroutine xc_potential(rho_iter)
      integer(ik) :: rho_iter ! For rho_iter>0, effective potential mixing makes sense
      !
      !  Change in density and density mixing
      !
      call TimerStart('SCF mixing')
      max_dens_delta = vector_max_change(old_dens,dens)
      if (scf_mode=='density' .or. scf_mode=='both') then
        call vector_merge_to(wgt_dst=rho_mix,dst=dens,wgt_src=1.0_rk-rho_mix,src=old_dens)
      end if
      call vector_copy(dst=old_dens,src=dens)
      call TimerStop('SCF mixing')
      !
      !  Calculate v12 contribution
      !
      call TimerStart('Mean-field potential')
      call add_v12_potential
      call TimerStop('Mean-field potential')
      !
      !  Calculate new exchange-correlation potential
      !
      call TimerStart('Exchange-correlation potential')
      call construct_exchange_correlation_potential
      call vector_scale_and_add_to(dst=vdia,scale=-beta,src=vxc)
      call TimerStop('Exchange-correlation potential')
      if (verbose>=5) then
        write (out,"(/' -beta*veff = '/(5(1x,g18.12))/)") vdia
      end if
      !
      !  Change in the effective potential and effective potential mixing
      !
      call TimerStart('SCF mixing')
      max_veff_delta = vector_max_change(old_vdia,vdia)
      if (scf_mode=='potential' .or. scf_mode=='both') then
        if (rho_iter>0) then
          call vector_merge_to(wgt_dst=veff_mix,dst=vdia,wgt_src=1.0_rk-veff_mix,src=old_vdia)
        end if
      end if
      call vector_copy(dst=old_vdia,src=vdia)
      call TimerStop('SCF mixing')
    end subroutine xc_potential
    !
    !  Massage the potential to improve convergence of the exponential
    !  expansion, but leave physics (almost) unchanged
    !
    subroutine clip_and_shift_potential(v_shift)
      real(rk), intent(out) :: v_shift ! Overall potential shift we applied
      !
      real(rk)    :: v_min, v_max ! Range of potential values
      real(rk)    :: v_limit      ! Cut-off value for clipping
      integer(ik) :: n_clip       ! Number of points clipped
      !
      !  Clip the repulsive part of the potential. On physical
      !  grounds, we know that it will have minimal impact on
      !  the total partition function, but is a major cause of
      !  slow convergence of the exponential series. Note that
      !  our potential is multiplied by (-beta), so that the
      !  strongest repulsion corresponds to the minimum of vdia,
      !  while the strongest attraction is given by its maximum.
      !
      call vector_minmax(src=vdia,min_src=v_min,max_src=v_max)
      v_limit = v_max - beta*v_clip*k_Boltzmann
      call vector_clip_below(dst=vdia,v_clip=v_limit,n_clip=n_clip)
      if (verbose>=2 .or. (n_clip>n_states/10 .and. verbose>=1) ) then
        write (out,"('-V(eff) : max = ',g14.7,' (',g12.6,' K)')") -v_min, -v_min/(beta*k_Boltzmann)
        write (out,"('          min = ',g14.7,' (',g12.6,' K)')") -v_max, -v_max/(beta*k_Boltzmann)
        write (out,"('Clip above -V = ',g14.7,' (',g12.6,' K)')") -v_limit, -v_limit/(beta*k_Boltzmann)
        write (out,"('Clipped ',f8.2,'% of sites (',i8,' out of ',i8,')')") &
               (100._rk*n_clip)/n_states, n_clip, n_states
      end if
      !
      !  Add potential shift, to minimize max(abs(vdia))
      !
      call vector_minmax(src=vdia,min_src=v_min,max_src=v_max)
      v_shift = 0.5_rk*(v_min+v_max)
      call vector_scale(dst=vdia,a=-v_shift,b=1.0_rk)
      if (verbose>=1) then
        write (out,"('Applying uniform potential shift = ',g12.6)") v_shift
      end if
    end subroutine clip_and_shift_potential
    !
    !  Single exponentiation loop: calculate density matrix given guess
    !  potential
    !
    subroutine exponentiate_hamiltonian(rho_iter,v_shift,used,allocated)
      integer(ik), intent(in) :: rho_iter        ! DFT: current SCF iteration number (affects printing)
      real(rk), intent(in)    :: v_shift         ! Overall potential shift (affects printing)
      real(ark), intent(out)  :: allocated, used ! Memory usage statistics
      !
      integer(ik)             :: iter     ! Power expansion term
      integer(ik)             :: i        ! Debugging counter
      real(rk)                :: q_term
      real(rk)                :: q_max_term ! Largest contribution to Q ever seen 
      !
      q_potential = 0 
      q_max_term  = 0
      q_term      = 1e5  ! Needed to shut off gfortran warning
      hpow        = SparseCreateUnitMatrix(n_states)
      call vector_zero(dens)
      partition_iterations: do iter=0,max_terms
        call TimerStart('Iterations')
        call SparseExtractDiagonal(hpow,hdia)
        call vector_add_to(dst=dens,src=hdia)
        !
        if (verbose>=5) then
          write (out,"(1x,a10,1x,a14  ,1x,a14  ,1x,a14  )") &
                'POINT', 'DENSITY', 'H^ITER', 'V0'
          write (out,"(1x,i10,1x,g14.7,1x,g14.7,1x,g14.7)") &
                (i, dens(i), hdia(i), vdia(i),i=1,n_states)
          write (out,"()")
        end if
        !
        q_term = vector_sum(hdia)
        q_potential = q_potential + q_term
        q_max_term = max(q_max_term,abs(q_term),abs(q_potential),vector_absmax(hdia))
        call TimerStop('Iterations')
        !
        if ((mode=='ideal' .and. verbose>=1) .or. &
            (mode=='dft' .and. rho_iter==0 .and. verbose>=1) .or. verbose>=2) then
          write (out,"(/'Iteration ',i5,'/',i5,' Q(Potential) = ',g14.7,' K = ',g14.7)") &
                 rho_iter, iter, q_potential*exp(v_shift), q_potential*exp(v_shift)/q_free
        end if
        !
        call flush(out)
        !
        if (iter<max_terms) then
          call SparseHamiltonianMultiply(hpow,neighbour,kmat(1:6)/(iter+1),vdia/(iter+1),power_screen)
          call SparseStatistics(hpow,max_used=used,max_allocated=allocated)
        end if
      end do partition_iterations
      call SparseDestroyMatrix(hpow)
      !
      !  Finish construction of the diagonal density - it's needed to complete
      !  evaluation of the free energy and to calculate vxc for the next iteration
      !
      call rescale_density
      !
      !  See how much accuracy we have left in Q
      !
      q_accuracy = spacing(q_max_term)/abs(q_potential)
      if ((mode=='ideal'.and.verbose>=1) .or. verbose>=2 .or. q_accuracy>=1e-3_rk) then
        write (out,"(t5,'estimated number of significant digits in Q = ',f5.1)") -log10(q_accuracy)
        if (q_accuracy>=3e-1_rk) then
          write (out,"(//8(' WARNING ')/t8,'TOTAL LOSS OF SIGNIFICANCE IN Q'/" // & 
                     "t8,'THE RESULTS ARE MEANINGLESS'/8(' WARNING ')//)")
          call conditional_abort
        end if
      end if
      if (abs(q_term)>=1e-5_rk*abs(q_potential)) then
        write (out,"('Exponentiation has not converged. Increase MAX_TERMS')")
        call conditional_abort
      end if
    end subroutine exponentiate_hamiltonian
    !
    !  Sanity check on the final density - are local pressures within the
    !  valid range for the LDA expression we use?
    !
    subroutine local_pressure_check
      if (p_max>=4e4_rk*bar) then
        write (out,"("// &
          "/t5,'********************************************************************',"// &
          "/t5,'**** WARNING: MAXIMUM SELF-CONSISTENT LOCAL PRESSURE IS OUTSIDE ****',"// &
          "/t5,'****   THE SAFE EXTRAPOLATION REGION. THE RESULTS ARE INVALID.  ****',"// &
          "/t5,'********************************************************************'/)")
      end if
      if (p_min<=-bar) then
        write (out,"("// &
          "/t5,'********************************************************************',"// &
          "/t5,'**** WARNING: MIMIMUM SELF-CONSISTENT LOCAL PRESSURE IS OUTSIDE ****',"// &
          "/t5,'****   THE SAFE EXTRAPOLATION REGION. THE RESULTS ARE INVALID.  ****',"// &
          "/t5,'********************************************************************'/)")
      end if
    end subroutine local_pressure_check
    !
    !  SCF procedure. The initial guess for the density is expected in
    !  the dens array.
    !
    subroutine scf_loop(have_guess)
      logical, intent(in) :: have_guess ! Initial guess is present
      integer(ik)         :: rho_iter
      real(rk)            :: f_old, f_new, p_ext, f_ext, f_dftextra
      real(rk)            :: v_shift
      real(ark)           :: allocated, used
      !
      !  All dirty work to initialize the SCF: 
      !    build the grid and connectivity lists
      !    calculate kinetic energy matrix
      !    calculate adsorption potential
      !    do some sanity checking on grid parameters
      !
      call scf_initialize
      !
      !  At this point: 
      !    v0dia contains adsorption potential + diagonal of the kinetic matrix
      !          it has been pre-multiplied by (-beta)
      !    kmat  contains kinetic matrix. It has been pre-multiplied by (-beta)
      !
      !  Let iterations begin. 
      !
      f_old = 0
      call vector_copy(dst=old_dens,src=dens)
      call vector_copy(dst=old_vdia,src=v0dia)
      density_iterations: do rho_iter=0,rho_iterations
        !
        !  Exchange-correlation potential evaluation
        !
        call vector_copy(dst=vdia,src=v0dia)
        if (mode=='dft' .and. (have_guess .or. rho_iter>0)) then
          call xc_potential(rho_iter) 
        end if
        !
        !  Clip the repulsive part of the potential, and shift whatever
        !  remains. This improves convergence dramatically!
        !
        call clip_and_shift_potential(v_shift) 
        !
        !  Calculate density matrix by exponentiating the shifted Hamiltonian
        !
        call exponentiate_hamiltonian(rho_iter,v_shift,used,allocated)
        !
        !  Evaluate the free energy.
        !
        if (mode=='dft') then
          !
          !  DFT requires taking an additional integral over the density. Here, we include
          !  only the leading statistical factor (-RT N Log(N)), which is consistent with
          !  the deifinition of the exchange-correlation functional.
          !
          !  We could use either the density from the previous iteration, or the
          !  current density to evaluate the correction. (At convergence, the choice
          !  should not matter.) Prior to May 25, 2008, we were using the final
          !  density here. With WDA-type functionals it is more sensible to use
          !  the density from the previous iteration - which is done now. If the
          !  old behaviour is essential, calling construct_exchange_correlation_potential()
          !  before dft_energy_correction() should recover it.
          !
          call dft_energy_correction(f_dftextra)
          f_new = -R_gas*temperature*log(q_potential*exp(v_shift)/particle_count) &
                + (f_dftextra-energy_v12)*(Hartree*N_Avogadro)
          energy_vext = vector_dot_product(dens,vads)*product(g_step)/particle_count
          energy_kin  = -(R_gas/(Hartree*N_Avogadro))*temperature*log(q_potential*exp(v_shift)/particle_count) &
                        - 2*energy_v12 - energy_vxc - energy_vext
        else 
          !
          !  Ideal gas case: include all statistical factors here.
          !
          f_new = -R_gas*temperature*(log(q_potential*exp(v_shift)) &
                - (log(particle_count)-1))
        end if
        !
        call external_pressure(f_new,p_ext)
        call h2_tp_f(temperature,p_ext,f_ext)
        !
        !  Remove potential shift from vdia, and convert it to Kelvin, in case
        !  it is plotted later on
        !
        call vector_scale(dst=vdia,a=-(v_shift-kmat(0))/(beta*k_Boltzmann),b=-1._rk/(beta*k_Boltzmann))
        if (debug_convergence) then
          call dump_density_and_potential(rho_iter)
        end if
        !
        if (verbose>=1) then
          call report_progress
        else 
          write (out,"(1x,'Iteration ',i4,', F= ',g14.7,' J/K-mol, Pext= ',g12.5,' bar')") &
                 rho_iter, f_new, p_ext/bar
        end if
        !
        if (mode/='dft') then
          if (verbose<1) call report_progress
          exit density_iterations
        end if
        !
        if (rho_iter>0 .and. abs(f_new-f_old)<eps_f .and. &
                       (max_dens_delta/abohr**3)<eps_rho ) then
          if (verbose<1) call report_progress
          exit density_iterations
        end if
        f_old = f_new
      end do density_iterations
      !
      if (mode=='dft') then
        call local_pressure_check
      end if
      !
      if (verbose>=0) then
        call print_comment
        call TimerReport
      end if
      !
      contains
        subroutine report_progress
          write (out,"(/t10,'At the end of density iteration ',i4)") rho_iter
          write (out,"(t4,'            Q(Potential)  = ',g14.7)") q_potential*exp(v_shift)
          write (out,"(t4,'            Q(Free,Grid)  = ',g14.7)") q_free
          write (out,"(t4,'            Q(Free,Limit) = ',g14.7)") q_free_limit
          write (out,"(t4,'            Equilibrium K = ',g14.7)") q_potential*exp(v_shift)/q_free
          write (out,"()")
          write (out,"(t4,'              F [J/mole]  = ',g14.7)") f_new
          if (rho_iter>0) then
            write (out,"(t4,'                 Delta F  = ',g14.7)") f_new - f_old
          end if
          write (out,"(t4,'Matching external P [bar] = ',g14.7)") p_ext / bar
          write (out,"(t4,'      External F [J/mole] = ',g14.7)") f_ext
          if ( (rho_iter>0 .or. have_guess) .and. (mode=='dft') ) then
            write (out,"(t4,'Min. local pressure [bar] = ',g14.7)") p_min / bar
            write (out,"(t4,'Max. local pressure [bar] = ',g14.7)") p_max / bar
            write (out,"(t4,'Max. dens. chg. [Angs^-3] = ',g14.7)") max_dens_delta / abohr**3
            write (out,"(t4,'     Max. veff change [K] = ',g14.7)") max_veff_delta / (beta*k_Boltzmann)
          end if
          if (mode=='dft') then
            write (out,"()")
            write (out,"(t4,' Free energy contributions: ',a18,1x,a9)") &
                   '[Hartree/particle]','[J/mole]'
            write (out,"(t4,'                  kinetic = ',g16.8,' = ',g16.8)") &
                   energy_kin,  energy_kin *(Hartree*N_Avogadro)
            write (out,"(t4,'       external potential = ',g16.8,' = ',g16.8)") &
                   energy_vext, energy_vext*(Hartree*N_Avogadro)
            write (out,"(t4,'               mean-field = ',g16.8,' = ',g16.8)") &
                   energy_v12,  energy_v12 *(Hartree*N_Avogadro)
            write (out,"(t4,'     exchange-correlation = ',g16.8,' = ',g16.8)") &
                   energy_fxc,  energy_fxc *(Hartree*N_Avogadro)
            write (out,"(t4,'               vxc energy = ',g16.8,' = ',g16.8)") &
                   energy_vxc,  energy_vxc *(Hartree*N_Avogadro)
            write (out,"()")
          end if
          write (out,"()")
          write (out,"(t4,' Max. memory used [Mbytes] = ',f10.3)") used/(1024._rk**2)
          write (out,"(t4,'    Max. memory allocated  = ',f10.3)") allocated/(1024._rk**2)
          write (out,"()")
       end subroutine report_progress
    end subroutine scf_loop
    !
    !  Initial guess driver, from increasingly quantum DFT
    !
    subroutine graded_density_guess
      real(rk)     :: mass_save, rho_mix_save, veff_mix_save
      real(rk)     :: eps_f_save, eps_rho_save
      integer(ik)  :: rho_iterations_save
      integer(ik)  :: mass_step
      real(rk)     :: f_guess, f_quant
      !
      !  Save parameters which will be adjusted during preliminary iterations
      !
      mass_save           = mass
      rho_mix_save        = rho_mix
      veff_mix_save       = veff_mix
      eps_f_save          = eps_f
      eps_rho_save        = eps_rho
      rho_iterations_save = rho_iterations
      !
      guess_iterations: do mass_step=1,mass_steps
        !
        !  Mass interpolation is a little tricky: I'd like to have linear
        !  steps in the logarithm of the de Broglie wavelength of the guest.
        !  For everything else, interpolate linearly between the limiting
        !  values.
        !
        f_quant = real(mass_step-1,kind=rk)/mass_steps
        f_guess = 1._rk - f_quant
        !
        mass           =     (mass_guess         **f_guess) * (mass_save         **f_quant)
        rho_mix        =      rho_mix_guess       *f_guess  +  rho_mix_save       *f_quant
        veff_mix       =      veff_mix_guess      *f_guess  +  veff_mix_save      *f_quant
        eps_f          =      eps_f_guess         *f_guess  +  eps_f_save         *f_quant
        eps_rho        =      eps_rho_guess       *f_guess  +  eps_rho_save       *f_quant
        rho_iterations = nint(rho_iterations_guess*f_guess  +  rho_iterations_save*f_quant)
        !
        write (out,"(/1x,80('*'))")
        write (out,"(15x,'Initial guess: guest mass step ',i4)") mass_step
        write (out,"( 1x,80('*')/)")
        !
        write (out,"('For the current guess iteration:'/)")
        write (out,"(t5,'                 Guest mass [EMU] = ',e12.5)") mass
        write (out,"(t5,'                   Density mixing = ',e12.5)") rho_mix
        write (out,"(t5,'                 Potential mixing = ',e12.5)") veff_mix
        write (out,"(t5,'Free energy convergence [J/mol-K] = ',e12.5)") eps_f
        write (out,"(t5,'Density convergence [Angstrom^-3] = ',e12.5)") eps_rho
        write (out,"(t5,'                  Iteration limit = ',i4)") rho_iterations
        write (out,"()")
        !
        !  Do an SCF loop with these parameters
        !
        call scf_loop(have_guess=mass_step>1)
      end do guess_iterations
      !
      !  Restore parameters to their final values
      !
      mass           = mass_save
      rho_mix        = rho_mix_save
      veff_mix       = veff_mix_save
      eps_f          = eps_f_save
      eps_rho        = eps_rho_save
      rho_iterations = rho_iterations_save
    end subroutine graded_density_guess
    !
    !  "Summary" input keywords
    !
    subroutine process_input_shortcuts
      select case (functional)
        case default
          write (out,"(/'Functional '',a,'' is not recognized.')") trim(functional)
          stop 'bad functional specification'
        case ('user','USER','User',' ')
        case ('lie-0','LIE-0')
          v12_form     = 'none'
          wda_shape    = 'Delta'
          write (out,"(/' Setting up LIE-0 functional:')")
          write (out,"('    V12_FORM  = ',a)") trim(v12_form)
          write (out,"('    WDA_SHAPE = ',a)") trim(wda_shape)
          write (out,"()")
        case ('lie-1-old','LIE-1-OLD')
          v12_form     = 'diep-'
          wda_shape    = 'Fermi-Dirac'
          v12_cutpoint = 2.827_rk
          wda_r0       = 2.197_rk
          wda_d        = 0.229_rk
          write (out,"(/' Setting up LIE-1 (old) functional:')")
          write (out,"('    V12_FORM     = ',a)") trim(v12_form)
          write (out,"('    WDA_SHAPE    = ',a)") trim(wda_shape)
          write (out,"('    V12_CUTPOINT = ',f14.8)") v12_cutpoint
          write (out,"('    WDA_R0       = ',f14.8)") wda_r0
          write (out,"('    WDA_D        = ',f14.8)") wda_d
          write (out,"()")
        case ('lie-1-flat','LIE-1-FLAT')
          v12_form     = 'diep-'
          wda_shape    = 'Fermi-Dirac'
          v12_cutpoint = 2.850_rk
          wda_r0       = 2.310_rk
          wda_d        = 0.0141_rk
          write (out,"(/' Setting up LIE-1 (flatline) functional:')")
          write (out,"('    V12_FORM     = ',a)") trim(v12_form)
          write (out,"('    WDA_SHAPE    = ',a)") trim(wda_shape)
          write (out,"('    V12_CUTPOINT = ',f14.8)") v12_cutpoint
          write (out,"('    WDA_R0       = ',f14.8)") wda_r0
          write (out,"('    WDA_D        = ',f14.8)") wda_d
          write (out,"()")
        case ('lie-1','LIE-1')
          v12_form   = 'diep-sph'
          wda_shape  = 'Fermi-Dirac'
          v12_sphmax = 180._rk
          wda_r0     = 2.500_rk
          wda_d      = 0.001_rk
          write (out,"(/' Setting up LIE-1 functional:')")
          write (out,"('    V12_FORM   = ',a)") trim(v12_form)
          write (out,"('    WDA_SHAPE  = ',a)") trim(wda_shape)
          write (out,"('    V12_SPHMAX = ',f14.8)") v12_sphmax
          write (out,"('    WDA_R0     = ',f14.8)") wda_r0
          write (out,"('    WDA_D      = ',f14.8)") wda_d
          write (out,"()")
      end select
    end subroutine process_input_shortcuts
    !
    !  Problem driver
    !
    subroutine pfun
      integer(ik)        :: info
      logical            :: have_guess
      real(rk)           :: rdum
      !
      call TimerStart('pfun')
      call accuracyInitialize
      rdum = MathLogFactorial(200)
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk 
      write (out,"('DFT statistical factors are included')")
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=partfun,iostat=info)
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=partfun)
      write (out,"()")
      call print_comment
      !
      !  Additional structure input - comes as an XYZ file. Note that the
      !  structure must be consistent with the positions and radii of the
      !  exclusion centres
      !
      call read_xyz_structure
      if (verbose>=0) call echo_xyz_structure
      !
      !  Input short-cuts: if summary parameters such as "functional"
      !  are present, set up the dependent parameters.
      !
      call process_input_shortcuts
      !
      !  Units conversion: Angstrom -> Bohr
      !
      beta            = 1/(k_Boltzmann*temperature)
      lat_vec         = lat_vec      / abohr
      lat_clip        = lat_clip     / abohr
      hard_wall_width = hard_wall_width / abohr
      v12_cutpoint    = v12_cutpoint / abohr
      wda_r0          = wda_r0       / abohr
      wda_d           = wda_d        / abohr
      if (n_exclude>0) xyzr_exclude = xyzr_exclude / abohr
      if (n_centres>0) centres      = centres      / abohr
      !
      !  Density parameters
      !
      call prepare_density_parameters
      !
      !  Initialize grid parameters
      !
      call prepare_grid_parameters
      !
      !  Build background plot
      !
      if (image_file/=' ') call plotStructure
      !
      !  How many states do we have?
      !
      n_states     = product(n_points)
      write (out,"('Dimension of the Hamiltonian is ',i10)") n_states
      write (out,"('Static memory overhead is ',f12.3,' MBytes')") &
            (n_states*(rk_bytes*15._rk+ik_bytes*6._rk))/1024._rk**2
      call flush(out)
      !
      allocate (xyz(3,n_states),neighbour(6,n_states),hdia(n_states), vads(n_states), &
                v0dia(n_states),vdia(n_states),old_vdia(n_states),dens(n_states), &
                old_dens(n_states),v12(n_states),vxc(n_states),vxc_n(n_states), &
                vxc_nave(n_states), eps_xc(n_states),smooth_dens(n_states),stat=info)
      if (info/=0) then
        write (out,"('Error ',i10,' allocating diagonal arrays')") info
        stop 'partition - no memory for diagonal arrays'
      end if
      !
      call initialize_v12
      !
      dens       = 0
      have_guess = .false.
      if (mode=='dft') then
        call TimerStart('Initial guess')
        select case (guess)
          case default
            write (out,"('Guess = ',a,' is not valid in a DFT calculation')") trim(guess)
            stop 'partition - invalid guess choice'
          case ('zero') ! Do nothing
          case ('read')  
            write (out,"(/'Reading initial density from ',a/)") trim(guess_file)
            call load_dx_field(guess_file,dens)
            call rescale_density(shout=.true.)
            have_guess = .true.
          case ('classical','graded')
            !
            !  During the first pass, we will solve classical DFT problem,
            !  to obtain a (rather good) starting guess for the density.
            !  The cost of the classical solution is essentially negligible,
            !  but it reduces the number of iterations by a large factor.
            !
            if (guess=='classical') mass_steps = 1
            verbose    = verbose - 1
            call graded_density_guess
            verbose    = verbose + 1
            have_guess = .true.
        end select
        call TimerStop('Initial guess')
        !
        write (out,"(/1x,80('*'))") 
        write (out,"(20x,'Quantum-liquid DFT calculation starts now')")
        write (out,"( 1x,80('*')/)") 
      end if ! mode=='dft'
      !
      call scf_loop(have_guess=have_guess)
      !
      if (dens_file/=' ') then
        call dump_dx_field(dens_file,dens,'Final density, particle/Bohr^3')
      end if
      !
      if (ascii_file/=' ') then
        call construct_exchange_correlation_potential
        call dump_ascii_density
      end if
      !
      if (veff_file/=' ') then
        call dump_dx_field(veff_file,vdia,'Effective potential, Kelvin')
      end if
      !
      call destroy_v12
      !
      call TimerStop('pfun')
      call print_comment
      call TimerReport
    end subroutine pfun
  end module partition
  !
  program driver
    use partition
    
    call pfun
  end program driver
