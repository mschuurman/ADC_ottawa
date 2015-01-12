!
!  Calculation of tunneling lifetimes in intense static fields, using
!  FON-SCF and MR-CI wavefunctions with absorbing boundary conditions.
!  Note that both FON-SCF and MR-CI are implemented internally in the 
!  basis of spin-orbitals, so that UHF is also available if needed.
!  Although the code is currently non-relativistic, it would be trivial
!  to convert it to two-component.
!
!  WARNING: Only FON-SCF is implemented right now!
!  WARNING: In fact, not even that.
!
!  This is not a fully stand-alone code: we lack sophisticated SCF and
!  CI solvers, and instead rely on simple iterations (SCF) and iterative
!  refinement (MR-CI) to calculate the solutions. As the result, we 
!  expect to be given a converged wavefunction in the absence of the
!  external field as the starting guess.
!
!  WARNING: This code does not support truncation of variational space
!           to spherical functions (GAMESS option ISPHER=+1). Using MOs
!           coming from such truncated calculation will likely cause
!           trouble!
!
!  WARNING: OpenMP code in this module is miscompiled by Intel Fortran 9.1
!  WARNING: It seems to work OK for Intel Fortran 12.1
!
!  WARNING: ECP implementation is based on atom-based resolution of identity.
!           Unless basis sets used in RI are saturated, the results may
!           differ from the exact integrals used by quantum-chemistry codes.
!
!  Some notes on the basis function and orbital layout
!  ---------------------------------------------------
!
!  Atomic basis sets are taken to be spin-orbitals. All alpha-spin AOs
!  are coming first, at positions (1:nao). The beta-spin AOs follow, at
!  positions (nao+1:nao_spin). An exception is made for the 2e AO integrals:
!  these are kept internally in the spin-less form to reduce storage 
!  requirements.
!
!  MOs are expanded over the spin-AOs and are intrinsically of mixed spin.
!  The guess MOs will be typically taken from GAMESS; these orbitals are
!  arranged with spin-alpha at odd orbital indices, and spin-beta at even
!  indices.
!
!  During the FON-SCF, we are rotating the degenerate MO pairs to maintain
!  maximum alignment to the initial guess orbitals. If the guess orbitals
!  are of pure spin, and we do not have any spin-orbital interactions present,
!  the final FON-SCF orbitals should also remain pure-spin, and should be
!  arranged in the odd-even sequence. However, this property should not
!  be relied upon.
!
  module static_tunnel
    use accuracy
    use timer
    use math
    use import_gamess
    use gamess_internal
    use integral_tools
    use integrals_mo2e
    use ci_tools
    use fock_tools
    use mp2_tools
    use ecp_convert_gamess
    use sort_tools
    use basis_cap
    use lapack
    use block_diag
    use matrix_tools
    use diis
    use biorthogonal_tools
    use scf_tools
    implicit none
    private
    public start, xk
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: iu_temp    = 36          ! I/O unit for subroutine-local I/O
    integer(ik), parameter :: iu_2e_ao   = 37          ! I/O unit used for storing 2e integrals over the atomic bfs
    integer(ik), parameter :: iu_mos     = 38          ! I/O unit used for punching MOs
    integer(ik), parameter :: iu_2e_mo   = 39          ! I/O unit used for storing 2e integrals over the MOs
    integer(ik), parameter :: max_fields = 100         ! Max number of field scaling factors
    integer(ik), parameter :: max_ormas  = 10          ! Max numver of ORMAS partitions
    integer(ik), parameter :: max_isa    = 20          ! Max number of intruder-state avoidance shifts
    !
    !  ==== Intermediate accuracy parameter ====
    !
    integer, parameter     :: xk              = xrk    ! Real kind used for handling potentially linearly-
                                                       ! dependent quantities. Should be used sparingly, and
                                                       ! hopefully never in a critical path!
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose         = 1              ! Level of output
    character(len=clen) :: molecule_file   = 'mol.dat'      ! Name of GAMESS checkpoint file, containing
                                                            ! the structure and (optionally) converged orbitals. Both
                                                            ! spin-restricted (RHF, MCSCF) and spin-unrestricted
                                                            ! (UHF) inputs are OK; the choice between RHF and UHF
                                                            ! guess is based on the number of MOs found in the file.
                                                            ! Note that internally, our FON-SCF is -always- UHF
    character(len=clen) :: orbital_guess   = ' '            ! Name of the file containing starting complex spin-orbitals 
                                                            ! (e.g. as punched in mos_file). Blank will use real orbitals
                                                            ! from molecule_file instead. Finally, 'project filename'
                                                            ! will tread 'filename' as a GAMESS checkpoint file, possibly
                                                            ! using a different basis set. Orbitals in this file will be
                                                            ! projected, one-by-one, onto the currently active basis.
    real(rk)            :: charge          = 0              ! Overall molecular charge; non-integer values 
                                                            ! are OK and will trigger FON-SCF. Must be consistent 
                                                            ! with mo_occ(:) below
    real(rk)            :: efield(3)       = 0              ! Static electric field in a.u.
    real(rk)            :: efield_current(3)                ! Static electric field on the current iteration
    integer(ik)         :: cnt_efield_scale= 2              ! Number of scaling factors to apply to efield in succession
                                                            ! The default is to do a calculation with a rather weak field
                                                            ! to get symmetry-adapted orbitals, then continue with the 
                                                            ! desired calculation.
    integer(ik)         :: constructor_i                    ! Stupid syntax rules ...
    real(rk)            :: efield_scale(max_fields) = (/ 0.001_rk, 1.0_rk, (0._rk, constructor_i=1,max_fields-2) /)
    integer(ik)         :: grid_nrad       = 120_ik         ! Basic number of radial points in atomic spheres; actual number of points 
                                                            ! may depend on this value and atom types
    integer(ik)         :: grid_nang       = 770_ik         ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik)         :: grid_outer_nrad =1200_ik         ! Basic number of radial points in the outer sphere
    integer(ik)         :: grid_outer_nang = 770_ik         ! Angular points in the outer sphere
    integer(ik)         :: max_scf_iter    = 100            ! Max number of SCF iterations
    character(len=clen) :: scf_converger   = 'diis'         ! Can be either 'diis' or 'damp'
    integer(ik)         :: max_diis_nvec   = 50             ! Maximum number of DIIS vectors allowed
    real(rk)            :: max_diis_coeff  = 20._rk         ! Restart DIIS if any of the coefficients exceed this threshold
    real(rk)            :: scf_mixing      = 0.1_rk         ! Mixing coefficient to use for the Fock matrix
    real(rk)            :: scf_eps_energy  = 1e-8_rk        ! Desired SCF convergence for the total energy
    real(rk)            :: scf_eps_rho     = 1e-6_rk        ! Desired SCF convergence for the density matrix
                                                            ! Default (1e-6) may be difficult to achieve for large, diffuse
                                                            ! basis sets, but is essential for getting sensible MP2 results!
    real(xk)            :: eps_smat        = 1e-7_rk        ! Threshold for linear dependence in AO overlap.
                                                            ! Zero will use near-machine accuracy, which almost certainly will cause
                                                            ! SCF to fail. The default (1e-7) is barely adquate for use with Kaufmann
                                                            ! basis sets with n_max=10; 1e-8 or 1e-9 would have been better, but
                                                            ! requires quadruple precision
    real(xk)            :: eps_geev        = 1e-5_rk        ! Threshold for declaring ZGEEV eigenvalues degenerate;
                                                            ! Zero means to go for machine accuracy; this should
                                                            ! be safe when it works, but may cause diagonalize_fmat to fail
    logical             :: follow_null     = .false.        ! Include null-space orbitals when trying to match orbitals in SCF.
    logical             :: absorb_fieldfree= .true.         ! Include absorbing boundary in the field-free case
    logical             :: skip_2e         = .false.        ! Omit all 2-electron terms; this flag is strictly for debugging!
    logical             :: skip_fieldfree  = .false.        ! Do not perform field-free calculation first; not recommended!
    integer(hik)        :: iosize_2e       =220000000_hik   ! Integral I/O buffer size in words for conventional SCF
    character(len=clen) :: scf_type        = 'conventional' ! Can be one of:
                                                            ! 'direct'        - compute all 2e AO integrals on the fly
                                                            ! 'incore'        - keep 2e integrals on disk
                                                            ! 'conventional'  - store integrals on disk
    character(len=clen) :: all_math        = 'as-is'        ! Override for all numerical accuracy settings. Can be one of:
                                                            ! 'as-is' - Use the default settings
                                                            ! 'real'  - Set all accuracy settings to 'real'
                                                            ! 'quad'  - Set all accuracy settings to 'quad'
    character(len=clen) :: diag_math       = 'real'         ! Numerical accuracy to use in diagonalizing the Fock matrix
                                                            ! Changing accuract of diag_math will also affect DIIS and
                                                            ! orbital following
    character(len=clen) :: ints_2e_math    = 'real'         ! Numerical accuracy to use for calculating 2-electron integrals.
                                                            ! Can be one of:
                                                            ! 'real' - Use standard accuracy for the integrals
                                                            ! 'quad' - Use extended accuracy for the integrals
    character(len=clen) :: fock_2e_math    = 'real'         ! Numerical accuracy to use for calculating the 2-electron
                                                            ! part of the Fock matrix. Can be one of:
                                                            ! 'real' - Use standard accuracy for fock matrix and density
                                                            ! 'quad' - Use extended accuracy for fock matrix and density
    character(len=clen) :: mp2_math        = 'real'         ! Orbital coefficients accuracy to use for MP2 step. Integral accuracy
                                                            ! is controlled separately.
    character(len=clen) :: mp2_storage     = 'real'         ! Storage accuracy for the transformed MO integrals in MP2.
    character(len=clen) :: scf_sort_mos    = 'never'        ! Can be one of: 
                                                            ! 'never'      - don't sort the MOs; try to maintain the guess 
                                                            !                orbitals order regardless of the energy
                                                            ! 'field-free' - sort MOs in the order of ascending real part of 
                                                            !                the energy after the field-free SCF converges
                                                            ! It's a bad idea to sort MOs once the field has been applied!
    character(len=clen) :: moint_mode      = 'incore'       ! Can be one of:
                                                            ! 'incore'        - keep transformed 2e MO integrals in memory
                                                            ! 'disk'          - keep them on disk
    logical             :: use_block_diag  = .true.         ! .True. if block-diagonalization routines should be used;
                                                            ! Otherwise, we will use straight LAPACK routines
    character(len=clen) :: mos_file        = 'field_mos.dat'! Quasi-GAMESS punch file for molecular orbitals
                                                            ! Appending this file to 'mol.dat' will make it
                                                            ! readable by our OpenDX routines (once I get around to fixing them, that is).
                                                            ! This output can also be used for getting the starting guess; it is better
                                                            ! than reading GAMESS MOS since we punch the complex spin-MO coefficients
                                                            ! Blank disables the output.
    !
    character(len=clen) :: ecp_file        = ' '            ! File containing the ECP and the RI basis; blank
                                                            ! means use ECP and basis in mol_file.
                                                            ! WARNING: standard molecular basis sets are likely 
                                                            !          not flexible enough to match GAMESS ECP energies!
    character(len=clen) :: occ_file        = 'guess'        ! File containing list of -all- spin-orbital occupation numbers; 
                                                            ! blank means that occupations follow the main namelist in the
                                                            ! input file. 'guess' is special, and will use aufbau principle to
                                                            ! fill in occupation numbers.
    real(rk)            :: ecp_eps_min     = 1e-6_rk        ! Small shift value cut-off (absolute) in ECP
    real(rk)            :: ecp_eps_max     = 1e+6_rk        ! Large shift value cut-off (positive) in ECP
    real(rk)            :: ecp_eps_grid    = 1e-6_rk        ! Characteristic grid spacing cut-off in ECP
    character(len=clen) :: ecp_report      = ' '            ! File to report ECP level-shift projectors to; empty name
                                                            ! suppresses the output.
    character(len=clen) :: correlation_type= 'none'         ! Treatment of correlation; Can be one of:
                                                            ! 'none'      - Just do SCF
                                                            ! 'ormas'     - Perform ORMAS-CI, no Sz constraint
                                                            ! 'sz ormas'  - ORMAS-CI, restricted to a specified Sz
                                                            ! 'mp2'       - MP2 energy
                                                            ! 'mp2 final' - MP2 energy, calculated at the full field strength only
    character(len=clen) :: ci_solver       = 'direct'       ! CI solver to use. Can be either of:
                                                            ! 'direct'   - Only useful for very small active spaces; done not
                                                            !              require a guess state vector
                                                            ! 'arpack'   - Iterative refinement of the guess vector; guess
                                                            !              must be supplied!
    integer(ik)         :: ci_max_iter     = 1000_ik        ! Maximum number of Arnoldi iterations in iterative solver
    integer(ik)         :: ci_n_root       = 3              ! Number of CI roots to find; we'll choose the right one based
                                                            ! on Cartesian distance of the eigenvector from the guess
    real(rk)            :: ci_eps_root     = 1e-6_rk        ! Desired convergence of the CI roots
    complex(rk)         :: ci_energy       = 0._rk          ! Guess for the energy of the CI state; MUST be given if ci_solver=='arpack'
    character(len=clen) :: ci_vector_file  = ' '            ! File containing guess for the CI state vector. Black will use a random guess
    integer(ik)         :: ci_nelectrons   = 0              ! Number of active electrons in CI
    real(rk)            :: ci_sz           = 0._rk          ! Restrict CI calculation to a given Sz
    integer(ik)         :: ormas_npart     = 0              ! Number of ORMAS partitions; can not exceed max_ormas????
    integer(ik)         :: ormas_orbitals(max_ormas+1)      ! Starting orbital index of each orbital in an ORMAS partitions
                                                            ! The last entry [ormas_orbitals(ormas_npart)] marks the orbital
                                                            ! one past the last active partition
    integer(ik)         :: ormas_mine    (max_ormas)        ! Minimum number of electrons allowed in a partition
    integer(ik)         :: ormas_maxe    (max_ormas)        ! Maximum number ... etc
    character(len=clen) :: mp2_ref_file    = 'guess'        ! Name of the file to read integer occupation numbers of the MP2
                                                            ! reference determinant from. 'guess' is special, and will occupy 
                                                            ! the lowest ci_nelectrons (or nel_scf of ci_nelectrons is non-positive)
                                                            ! orbitals; ' ' is also special, and will cause reading from the 
                                                            ! stanard input (compare to occ_file)
    character(len=clen) :: mp2_mode        = 'incore'       ! MP2 mode must be 'incore'. The 'stupid' version has been deleted.
                                                            ! Except for very small inputs, 'incore' is the only reasonable choice;
                                                            ! it will however need 2*Ne**2*(Norb-Nel)**2 complex words for the integrals
    real(rk)            :: mp2_isa_de(max_isa) = reshape((/ 0.02_rk, 0.04_rk, 0.08_rk, 0.16_rk, 0.32_rk /), (/max_isa/), (/-1._rk/))
                                                            ! Intruder-state avoidance parameter in MP2. Energy denominator
                                                            ! in MP2 will never drop below 2*mp2_isa_de (Hartree)
                                                            ! Setting mp2_isa_de to zero will disable ISA
                                                            ! The list terminates when a negative value is reached 
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type(gam_structure)         :: gam               ! Structure descriptor; for use with both the field-free and static-field calculations
    type(ecp_molecule), target  :: ecp               ! ECPs (if any)
    type(diis_state)            :: diis_st           ! DIIS state, see diis.f90
    integer(ik)                 :: natoms            ! Number of atoms in the structure
    integer(ik)                 :: nao               ! Number of spin-less atomic basis functions.
    integer(ik)                 :: nao_spin          ! Number of atomic basis function with spin factors included
    integer(ik)                 :: nmo               ! Number of molecular spin-orbitals. 
    integer(ik)                 :: nmo_null          ! Number of molecular spin-orbitals within the null-space of
                                                     ! of the overlap matrix. These MOs will be kept at the end
                                                     ! of the vector table.
    real(rk)                    :: nel_scf           ! Number of electrons in FON-SCF; does not have to be an integer 
    integer(ik)                 :: nocc_scf          ! Number of occupied MOs in FON-SCF.
    real(xrk), allocatable      :: mo_occ(:)         ! Occupation numbers of molecular orbitals in FON-SCF. The numbers 
                                                     ! are expected to be in the [0:1] range, with odd entries being
                                                     ! the alpha-spin orbitals. Even entries are the betas.
    real(rk), allocatable       :: mo_occ_rk(:)      ! ditto, standard real kind
    complex(xk), allocatable    :: mo_energy(:)      ! Energies of molecular orbitals
    complex(rk), allocatable    :: mo_energy_rk(:)   ! Energies of molecular orbitals, standard precision (for diag_math='real')
    complex(xk), allocatable    :: mo_energy0(:)     !   ... field-free
    complex(xk), allocatable    :: mo_energyf(:)     !   ... in the presence of the field
    complex(xk), allocatable    :: mos (:,:,:)       ! Molecular orbitals for the currently active SCF. 
                                                     ! For all MO arrays, the indices are:
                                                     !   1 = spin-AO
                                                     !   2 = spin-MO
                                                     !   3 = 1 for the left eigenvectors; 2 for the right eigenvectors
    complex(rk), allocatable    :: mos_rk(:,:,:)     ! Molecular orbitals, standard precision (for diag_math='real')
    complex(xk), allocatable    :: mosg(:,:,:)       ! Guess molecular orbitals
    complex(rk), allocatable    :: mosg_rk(:,:,:)    ! ditto, standard precision
    complex(xk), allocatable    :: mos0(:,:,:)       ! Converged (or converging ;) molecular orbitals in the absence of the field
    complex(xk), allocatable    :: mosf(:,:,:)       ! ditto, in the field
    complex(xk), allocatable    :: rho (:,:)         ! Electronic density matrix (AO basis)
    complex(rk), allocatable    :: rho_rk (:,:)      ! Same as rho(), but using rk-kind complex numbers
    complex(xk), allocatable    :: rho_old(:,:)      ! Electronic density matrix from previous SCF iteration
    real(xk), allocatable       :: smat(:,:)         ! Overlap matrix (AO basis), null-space projected out
    real(rk), allocatable       :: smat_rk(:,:)      ! ditto, standard precision
    real(xk), allocatable       :: sphalf(:,:)       ! S^{+1/2}, null-space is projected out
    real(rk), allocatable       :: sphalf_rk(:,:)    ! ditto, standard precision
    real(xk), allocatable       :: smhalf(:,:)       ! S^{-1/2}, null-space is projected out
    real(rk), allocatable       :: smhalf_rk(:,:)    ! ditto, standard precision
    real(xk), allocatable       :: h0  (:,:)         ! Field-free part of the current Hamiltonial, including the ECP terms
    complex(xk), allocatable    :: h0f (:,:)         ! External field contributions
    complex(rk), allocatable    :: h0c (:,:)         ! CAP contributions
    complex(xk), allocatable    :: hmat(:,:)         ! Current 1-electron Hamiltonian matrix
    complex(xk), allocatable    :: gmat(:,:)         ! 2-electron contribution to the Fock matrix
    complex(rk), allocatable    :: gmat_rk(:,:)      ! 2-electron contribution to the Fock matrix, in standard real precision
    complex(xk), allocatable    :: fmat(:,:)         ! Fock matrix
    complex(rk), allocatable    :: fmat_rk(:,:)      ! Fock matrix, standard precision
    complex(xk), allocatable    :: fmat_old(:,:)     ! Fock matrix from the previous iteration
    real(xk)                    :: enuc              ! Nuclear repulsion energy
    complex(xk)                 :: escf              ! SCF electronic energy; since Hamiltonian is non-self-adjoint, energy is complex
    integer(ik), allocatable    :: active_mos(:)     ! List of the MO indices needed for electron correlation step
                                                     ! Only meaningful if correlation_type/='none'
    real(rk), allocatable       :: orbital_sz(:)     ! Per-orbital Sz value, defined as <L|Sz|R>
    integer(sik), allocatable   :: mp2_ref_occ(:)    ! Occupation numbers of the MP2 reference determinat; only
                                                     ! relevant of correlation_type == 'mp2'
    complex(rk), allocatable    :: ci_vector_guess(:)! Guess for the right CI eigenvector
    complex(rk), allocatable    :: ci_vector(:)      ! Converged right CI eigenvector
    !
    !  Parameters related to the absorbing boundary
    !
    real(rk), parameter         :: ma_c = 2.622057554292119810464840_rk ! The position of the CAP singularity
    real(rk), parameter         :: ma_a = 1._rk - 16._rk*ma_c**(-3)
    real(rk), parameter         :: ma_b = (1._rk - 17._rk*ma_c**(-3)) * ma_c**(-2)
    real(rk)                    :: cap_width         ! Overall width of the CAP
    real(rk)                    :: cap_scale         ! Scaling parameter
    !
    !  Variables related to 2-electron integral handling
    !
    type(int2e_cache)   :: int2e              ! Currently active integrals context
    type(moint2e_cache) :: moint2e            ! Currently active MO integrals context
    type(ci_data)       :: ci                 ! Currently active CI
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /tunnel_data/ verbose, molecule_file, orbital_guess, mos_file, charge, occ_file, &
                           max_scf_iter, scf_mixing, scf_eps_energy, scf_eps_rho, &
                           scf_converger, max_diis_nvec, max_diis_coeff, &
                           eps_smat, eps_geev, follow_null, absorb_fieldfree, &
                           ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
                           scf_type, scf_sort_mos, iosize_2e, &
                           efield, cnt_efield_scale, efield_scale, &
                           cap_type, cap_centre, cap_r0, &
                           cap_strength, cap_order, &
                           cap_lambda, cap_theta, cap_mpole, cap_efield, cap_diff_step, &
                           cap_kmin, cap_delta, cap_limit, &
                           grid_nrad, grid_nang, grid_outer_nrad, grid_outer_nang, &
                           skip_2e, skip_fieldfree, use_block_diag, &
                           correlation_type, ci_solver, ci_nelectrons, ci_sz, &
                           ci_max_iter, ci_n_root, ci_eps_root, ci_energy, ci_vector_file, &
                           ormas_npart, ormas_orbitals, ormas_mine, ormas_maxe, moint_mode, &
                           mp2_ref_file, mp2_mode, mp2_isa_de, &
                           all_math, diag_math, fock_2e_math, mp2_math, mp2_storage, ints_2e_math
    !
    contains
    !
    subroutine print_geometry(gam)
      type(gam_structure), intent(in) :: gam
      integer(ik)  :: iat
      real(rk)     :: xyz(3), q
      !
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,gam%natoms
        xyz = real(gam%atoms(iat)%xyz,kind=rk)
        q   = real(gam%atoms(iat)%znuc,kind=rk)
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") q, xyz/abohr, xyz
      end do print_atoms
      write (out,"()")
    end subroutine print_geometry
    !
    subroutine load_gamess_data
      integer(ik)  :: iat, jat
      integer(ik)  :: nv
      real(xk)     :: xyz(3), q, r
      !
      !  Load the structure and copy some basic parameters to the global arrays
      !
      call gamess_load_orbitals(file=trim(molecule_file),structure=gam)
      write (out,"(/'Loaded GAMESS checkpoint file ',a/)") trim(molecule_file)
      natoms   = gam%natoms
      nao      = gam%nbasis
      nao_spin = 2*nao
      nmo      = 2*nao
      nel_scf  = sum(real(gam%atoms(:)%znuc,kind=kind(nel_scf))) - charge
      enuc     = 0.0_rk
      nuclear_repulsion: do iat=2,natoms
        xyz = real(gam%atoms(iat)%xyz,kind=kind(xyz))
        q   = real(gam%atoms(iat)%znuc,kind=kind(q))
        do jat=1,iat-1
          r = sqrt(sum((real(gam%atoms(jat)%xyz,kind=kind(xyz))-xyz)**2)) / abohr
          enuc = enuc + q*real(gam%atoms(jat)%znuc,kind=kind(q)) / r
        end do
      end do nuclear_repulsion
      !
      !  A bit of sanity checking
      !
      if (orbital_guess==' ') then
        !
        !  We'll be using real orbitals from GAMESS-US file as the starting guess,
        !  make sure the input is at least comprehensible
        !
        if (gam%nvectors/=nao .and. gam%nvectors/=nmo) then
          write (out,"('Number of MOs (',i0,') does not match the number of basis functions (',i0,') or MOs (',i0,')')") &
                 gam%nvectors, nao, nmo
          if (gam%nvectors>nao .and. gam%nvectors<nmo .and. mod(gam%nvectors,2)==0) then
            write (out,"('Assuming calculation is from linearly-dependent spin-unrestricted run. Padding with zero MOs.'/)")
            nv = gam%nvectors / 2
            gam%vectors(:,nao:nao+nv-1) = gam%vectors(:,nv+1:2*nv)
            gam%vectors(:,nv+1:nao) = 0
            gam%vectors(:,nao+nv:nmo) = 0
            gam%nvectors = nmo
          else if (gam%nvectors<nao) then
            write (out,"('Assuming calculation is from linearly-dependent spin-restricted run. Padding with zero MOs.'/)")
            gam%vectors(:,gam%nvectors+1:nao) = 0
            gam%nvectors = nao
          else
            write (out,"('Can''t figure out what to do with this set of MOs.')")
            call stop('static_tunnel%load_gamess_data - incomplete')
          end if
        end if
      end if
      if (nel_scf<0.0_rk .or. nel_scf>2.0_rk*nao) then
        write (out,"('Number of electrons (',f0.5,') is strange.')") nel_scf
        call stop('static_tunnel%load_gamess_data - bad electron count')
      end if
      if (verbose<0) return
      !
      !  Tell a bit more!
      !
      write (out,"('                    Number of atoms = ',i0)") natoms
      write (out,"(' Number of spinless basis functions = ',i0)") nao
      write (out,"('    Number of MOs in GAMESS-US file = ',i0)") gam%nvectors
      write (out,"('                Number of electrons = ',f0.5)") nel_scf
      write (out,"('                       Total charge = ',f0.5)") charge
      call print_geometry(gam)
      write (out,"('  Nuclear repulsion energy = ',f20.12)") enuc
      write (out,"()")
    end subroutine load_gamess_data
    !
    !  Convert GAMESS ECPs to the projector form. 
    !
    subroutine prepare_ecps
      type(gam_structure) :: gamx
      integer(ik)         :: iat
      !
      call TimerStart('Prepare ECPs')
      if (ecp_file==' ') then
        if (any(gam%atoms(:)%ecp_nterms>0)) then
          write (out,"(/'WARNING: Using molecular basis for ECP resolution of identity'/)")
        end if
        call ecp_convert(gam,ecp,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report)
      else
        if (any(gam%atoms(:)%ecp_nterms>0)) then
          write (out,"('Loading ECP RI basis from ',a)") trim(ecp_file)
        end if
        call gamess_load_orbitals(file=ecp_file,structure=gamx)
        if (gamx%natoms/=gam%natoms) call stop('ECP atom count mismatch')
        !
        forall (iat=1:gam%natoms) gamx%atoms(iat)%xyz = gam%atoms(iat)%xyz
        call ecp_convert(gamx,ecp,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report)
        !
        call gamess_destroy(gamx)
      end if
      call TimerStop('Prepare ECPs')
    end subroutine prepare_ecps
    !
    !  Memory allocated by allocate_dynamic_data() is never released; this is by design.
    !
    subroutine allocate_dynamic_data
      integer(ik) :: alloc
      real(rk)    :: rnao
      integer(ik) :: xk_bytes
      !
      xk_bytes = rk_bytes
      if (xk/=rk) xk_bytes = xrk_bytes
      if (verbose>=0) then
        rnao = nao_spin ! To force calculation to floating point!
        write (out,"(/'Allocating about ',f0.6,' GBytes in quadratic arrays.'/)") &
               (1.0_rk/1024.0_rk**3) * ( &
               rk_bytes * (11*rnao**2) &
             + xk_bytes * (4*2*2*rnao*nmo + 2*2*rnao**2 + 4*rnao**2 + 6*2*rnao) )
      end if
      !
      allocate (mo_occ(nmo), mo_occ_rk(nmo), &
                mo_energy(nmo), mo_energy_rk(nmo), mo_energy0(nmo), mo_energyf(nmo), &
                mos(nao_spin,nmo,2), mos_rk(nao_spin,nmo,2), &
                mosg(nao_spin,nmo,2), mosg_rk(nao_spin,nmo,2), mos0(nao_spin,nmo,2), mosf(nao_spin,nmo,2), &
                rho(nao_spin,nao_spin), rho_rk(nao_spin,nao_spin), rho_old(nao_spin,nao_spin), &
                smat(nao_spin,nao_spin), smat_rk(nao_spin,nao_spin), &
                sphalf(nao_spin,nao_spin), sphalf_rk(nao_spin,nao_spin), &
                smhalf(nao_spin,nao_spin), smhalf_rk(nao_spin,nao_spin), &
                h0(nao_spin,nao_spin), h0f(nao_spin,nao_spin), h0c(nao_spin,nao_spin), &
                hmat(nao_spin,nao_spin), gmat(nao_spin,nao_spin), gmat_rk(nao_spin,nao_spin), &
                fmat(nao_spin,nao_spin), fmat_rk(nao_spin,nao_spin), &
                fmat_old(nao_spin,nao_spin), stat=alloc)
      if (alloc/=0) then
        write (out,"('static_tunnel%allocate_dynamic_data: Error ',i0,' allocating quadratic arrays. nao = ',i0)") &
               alloc, nao
        call stop('static_tunnel%allocate_dynamic_data - out of memory (1)')
      end if
      select case (scf_converger)
        case ('diis')
          call diis_initialize(diis_st,nao_spin,max_diis_nvec,max_diis_coeff,diag_math)
      end select
    end subroutine allocate_dynamic_data
    !
    subroutine fill_occupations
      integer(ik)         :: iu      ! Can be either input or iu_temp, depending on whether occ_file is blank or not
      integer(ik)         :: ios
      integer(ik)         :: imo, kmo
      real(rk)            :: e_left  ! Number of electrons left to assign
      real(rk)            :: e_this  ! Number of electrons to go to the current MO
      !
      mo_occ = 0
      if (occ_file/='guess') then
        !
        !  We are given explicit occupation numbers; read them in.
        !
        iu = input
        if (occ_file/=' ') then
          iu = iu_temp
          open(iu,file=trim(occ_file),status='old',action='read',iostat=ios)
          if (ios/=0) then
            write (out,"('Error ',i0,' opening MO occupation file ',a)") ios, trim(occ_file)
            call stop('static_tunnel%fill_occupations - bad open')
          end if
        end if
        read (iu,*) mo_occ
        if (occ_file/=' ') close(iu)
        !
        !  Bit of sanity checking ...
        !
        if (any(mo_occ>1.0) .or. any(mo_occ<0.0)) then
          write (out,"('MO occupation numbers must be in the range [0:1]')")
          write (out,"((10(1x,f10.5)))") mo_occ
          call stop('static_tunnel%fill_occupations - bad occupations (1)')
        end if
        nocc_scf = 0
        find_nocc: do imo=nmo,1,-1
          if (mo_occ(imo)/=0) then
            nocc_scf = imo
            exit find_nocc
          end if
        end do find_nocc
        if (abs(sum(mo_occ(:nocc_scf))-nel_scf)>10*spacing(nel_scf)) then
          write (out,"('Sum of orbital occupations (',g24.14,') does not match number of electrons (',g24.14,')')") &
                 sum(mo_occ(:nocc_scf)), nel_scf
          call stop('static_tunnel%fill_occupations - bad occupations (2)')
        end if
      else
        !
        !  No explicit occupations were given; use aufbau principle to assign occupations
        !
        e_left   = nel_scf
        nocc_scf = 0
        fill_occ: do imo=1,nmo
          e_this = min(e_left,1.0_rk)
          mo_occ(imo) = e_this
          nocc_scf = nocc_scf + 1
          e_left = e_left - e_this
          if (e_left<=0.0_rk) exit fill_occ
        end do fill_occ
        if (e_left>0.0_rk) call stop('static_tunnel%fill_occupations - unallocated electrons')
      end if
      if (verbose<0) return
      !
      write (out,"(t5,a)") 'MO occupation numbers', &
                           '---------------------'
      print_occ: do imo=1,nmo,10
        kmo = min(imo+9,nmo)
        write (out,"(1x,i5,2x,10(1x,f7.5))") imo, mo_occ(imo:kmo)
      end do print_occ
      write (out,"()")
    end subroutine fill_occupations
    !
    subroutine fetch_initial_mos_real
      mosg(:,:,1) = 0
      if (gam%nvectors==nao) then
        !
        !  spin-restricted orbitals on input; double them up
        !
        mosg(    1:nao     ,1:nmo-1:2,1) = gam%vectors(1:nao,1:nao)
        mosg(nao+1:nao_spin,2:nmo  :2,1) = gam%vectors(1:nao,1:nao)
      else if (gam%nvectors==nmo) then
        !
        !  spin-unrestricted orbitals on input
        !
        mosg(    1:nao     ,1:nmo-1:2,1) = gam%vectors(1:nao,    1:nao)
        mosg(nao+1:nao_spin,2:nmo  :2,1) = gam%vectors(1:nao,nao+1:nao_spin)
      else
        call stop('static_tunnel%fetch_initial_mos_real - unexpected orbital count')
      end if
      !
      !  Set right eigenvectors to match the left ones (ie assume that guess
      !  orbitals come from a self-adjoint Hamiltonian).
      !
      mosg(:,:,2) = mosg(:,:,1)
    end subroutine fetch_initial_mos_real
    !
    subroutine fetch_initial_mos_complex
      integer(ik) :: ios
      integer(ik) :: file_line, imo
      integer(ik) :: iao, ifield
      integer(ik) :: chk_mo, chk_line, line
      complex(rk) :: loc_val(4)
      !
      write (out,"(/'Loading initial complex orbitals from: ',a/)") trim(orbital_guess)
      open (gam_file,file=trim(orbital_guess),action='read',position='rewind',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"('static_tunnel%fetch_initial_mos_complex - error ',i8,' opening ',a)") ios, trim(orbital_guess)
        call stop('static_tunnel%fetch_initial_mos_complex - nofile')
      end if
      file_line = 0
      scan_lines: do 
        call gam_readline ; file_line = file_line + 1
        if (gam_line_buf(1:6)==' $CVEC') exit scan_lines
      end do scan_lines
      !
      !  This section is ripped of import_gamess%parse_orbital_data
      !
      iao = 1 ; imo = 0 ; line = 0
      read_vectors: do
        call gam_readline ; file_line = file_line + 1
        line = line + 1
        if (gam_line_buf==' $END') exit read_vectors
        if (imo+1>2*nmo) then
          write (out,"('static_tunnel%fetch_initial_mos_complex: Too many lines in $CVEC at line ',i0,':'/1x,a)") &
                 file_line, trim(gam_line_buf)
          call stop('static_tunnel%fetch_initial_mos_complex - too many lines')
        end if
        !
        read (gam_line_buf,"(i2,i3,8g20.13)",iostat=ios) chk_mo, chk_line, loc_val
        if (ios/=0) then
          write (out,"('static_tunnel%fetch_initial_mos_complex: format error ',i0,' at line ',i0,':'/1x,a)") &
                 ios, file_line, trim(gam_line_buf)
          call stop('static_tunnel%fetch_initial_mos_complex - format error')
        end if
        ! 
        !  Line parsed, check it for validity
        !
        if (modulo(imo+1,100)/=chk_mo) then
          write (out,"('static_tunnel%fetch_initial_mos_complex: Expecting MO ',i0,', got ',i0,' at line ',i0,':'/1x,a)") &
                   imo+1, chk_mo, file_line, trim(gam_line_buf)
          call stop('static_tunnel%fetch_initial_mos_complex - MO mismatch')
        end if
        if (line/=chk_line) then
          write (out,"('static_tunnel%fetch_initial_mos_complex: Expecting line ',i0,', got ',i0,' at line ',i0,':'/1x,a)") &
                   line, chk_line, file_line, trim(gam_line_buf)
          call stop('import_gamess%parse_orbital_data - line mismatch')
        end if
        !
        !  Stuff vector in
        !
        stuff_c: do ifield=1,4
          mosg(iao,imo/2+1,modulo(imo,2)+1) = loc_val(ifield)
          iao = iao + 1
          if (iao>nao_spin) then
            imo  = imo + 1 ; iao  = 1 ; line = 0
            exit stuff_c
          end if
        end do stuff_c
      end do read_vectors
      if (imo/=2*nmo) then
        write (out,"('static_tunnel%fetch_initial_mos_complex: Wanted ',i0,' vectors, but got ',i0)") 2*nmo, imo
        call stop('static_tunnel%fetch_initial_mos_complex - bad vector count')
      end if
      close (gam_file)
    end subroutine fetch_initial_mos_complex
    !
    subroutine project_initial_orbitals(guess_file)
      character(len=*), intent(in) :: guess_file ! GAMESS checkpoint file containing the orbitals
      !
      integer(ik)           :: alloc
      type (gam_structure)  :: guess
      integer(ik)           :: guess_nmo, guess_nao
      real(xk), allocatable :: s_tt(:,:), s_tg(:,:), rhs(:,:)
      !
      call TimerStart('Project initial guess')
      write (out,"(/'Projecting initial guess orbitals from ',a)") trim(guess_file)
      call gamess_load_orbitals(file=guess_file,structure=guess)
      guess_nao = guess%nbasis
      guess_nmo = guess%nvectors
      if (guess%nvectors<=0) call stop('static_tunnel%project_initial_orbitals - no MOs in the guess file')
      !
      allocate (s_tt(nao,nao),s_tg(nao,guess_nao),rhs(nao,guess_nmo),stat=alloc)
      if (alloc/=0) call stop('static_tunnel%project_initial_orbitals - allocation failed (1)')
      call gamess_1e_integrals('AO OVERLAP',s_tt,bra=gam,ket=gam  )
      call gamess_1e_integrals('AO OVERLAP',s_tg,bra=gam,ket=guess)
      !
      rhs = mt_matmul(s_tg,guess%vectors(:,:guess_nmo))
      call lapack_gelss(s_tt,rhs)
      !
      if (allocated(gam%vectors)) deallocate (gam%vectors)
      allocate (gam%vectors(nao,nao),stat=alloc)
      if (alloc/=0) call stop('static_tunnel%project_initial_orbitals - allocation failed (2)')
      gam%nvectors = nao
      gam%vectors  = 0
      gam%vectors(:,1:min(nao,guess_nmo)) = rhs(:,1:min(nao,guess_nmo))
      !
      !  We've replaced the vectors in the gam% structure with the projected guess MOs;
      !  let fetch_initial_mos_real() to finish the job.
      !
      call fetch_initial_mos_real
      !
      deallocate (s_tt,s_tg,rhs)
      call gamess_destroy(guess)
      call TimerStop('Project initial guess')
    end subroutine project_initial_orbitals
    !
    !
    subroutine fetch_initial_mos
      !
      if (orbital_guess==' ') then
        call fetch_initial_mos_real
      else if (orbital_guess(1:8)=='project ') then
        call project_initial_orbitals(orbital_guess(9:))
      else 
        call fetch_initial_mos_complex
      end if
    end subroutine fetch_initial_mos
    !
    !  Since our integrals package operates with spin-free AOs,
    !  we need to do a bit of juggling around to accommodate the
    !  spin-AOs we use here.
    !
    subroutine evaluate_1e_hamiltonian
      integer(ik)                 :: iat, ic, alloc
      real(xk)                    :: xyz(3), q
      real(rk)                    :: maxs, maxserror
      character(len=1), parameter :: icc(3) = (/ 'X', 'Y', 'Z' /)
      real(xk), allocatable       :: ssf(:,:)                                 ! Quantities needed in high accuracy
      real(xk), allocatable       :: hsf(:,:), tmp_xk(:,:)
      real(rk), allocatable       :: tmp_rk(:,:), ssfn(:,:), lsfn(:,:)        ! SF stands for "spin-free", not "science fiction"
      complex(xk), allocatable    :: fsf(:,:)
      complex(rk), allocatable    :: csf(:,:)  
      !
      call TimerStart('1e Hamiltonian')
      !
      allocate (hsf(nao,nao),ssf(nao,nao),tmp_xk(nao,nao),tmp_rk(nao,nao), &
                ssfn(nao,nao),lsfn(nao,nao),fsf(nao,nao),csf(nao,nao),stat=alloc)
      if (alloc/=0) then
        write (out,"('static_tunnel%field_free_1e_hamiltonian: Error ',i0,' allocating spin-free temporaries')") alloc
        call stop('static_tunnel%field_free_1e_hamiltonian - no memory')
      end if
      !
      if (verbose>=0) then
        write (out,"(/'ECP and CAP terms are evaluated using kind-',i0,' arithmetics.')") rk
        write (out,"( 'All other 1-electron terms are calcuated using kind-',i0,' floating point.'/)") xk
      end if
      !
      call gamess_1e_integrals('AO KINETIC',hsf,bra=gam,ket=gam)
      if (verbose>=3) then
        write (out,"(/t5,'KINETIC ENERGY INTEGRALS'/)")
        call gamess_print_1e_integrals(hsf,bra=gam,ket=gam)
      end if
      !
      !  Fill nuclear attraction integrals.
      !
      nuclear_attraction: do iat=1,natoms
        xyz = real(gam%atoms(iat)%xyz,kind=kind(xyz)) / abohr
        q   = real(gam%atoms(iat)%znuc,kind=kind(q))
        call gamess_1e_integrals('AO 3C 1/R',tmp_xk,bra=gam,ket=gam,op_xyz=xyz)
        tmp_xk = -q * tmp_xk
        if (verbose>=4) then
          write (out,"(/t5,'NUCLEAR ATTRACTION TO ATOM ',i3,' Z= ',f12.5,' XYZ= ',3f12.5/)") iat, q, xyz
          call gamess_print_1e_integrals(tmp_xk,bra=gam,ket=gam)
        end if
        hsf = hsf + tmp_xk
      end do nuclear_attraction
      !
      tmp_rk = 0
      call ecp_evaluate_matrix_elements(gam,ecp,tmp_rk)  ! Put ECP terms in tmp
      hsf = hsf + tmp_rk
      !
      fsf = 0
      efield_loop: do ic=1,3
        call gamess_1e_integrals('AO DIPOLE '//icc(ic),tmp_xk,bra=gam,ket=gam)
        fsf = fsf + tmp_xk * efield(ic) ! Electrons are negative; since E = -Grad Scalar_potential, we get overall plus
      end do efield_loop
      call cap_evaluate(gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang,csf,ssfn,lsfn)  ! CAP terms are done numerically
      ! 
      !  For a non-local CAP, check Laplacian matrix elements
      !
      if (cap_type=='moiseyev' .or. cap_type=='atom moiseyev') then
        call gamess_1e_integrals('AO KINETIC',tmp_rk,bra=gam,ket=gam)
        tmp_rk = -2.0_rk * tmp_rk ! Convert them back to a Laplacian
        maxs = maxval(abs(tmp_rk))
        maxserror = maxval(abs(lsfn-tmp_rk))
        write (out,"(/'Largest difference between analytical and numerical AO Laplacian was ',g20.12)") maxserror
        write (out,"( '                                 Largest analytical AO Laplacian was ',g20.12/)") maxs
        if (maxserror>1e-4_rk*maxs) then
          write (out,"('Numerical integration accuracy is insufficient for sensible results. Threshold = ',g20.12)") &
                 1e-4_rk * maxs
          call stop('static_tunnel%evaluate_1e_hamiltonian - bad numerical integration (1)')
        end if
      end if
      !
      !  Sanity check on quad-precision math; costs us next to nothing.
      !
      call gamess_1e_integrals('AO OVERLAP',ssf,   bra=gam,ket=gam)
      call gamess_1e_integrals('AO OVERLAP',tmp_rk,bra=gam,ket=gam)
      !
      !  Now do the check of the numerical integration grid for the overlap matrix elements
      !  Also check the difference in regular and quad-precision overlap intergrals; these
      !  better agree as well!
      !
      maxs = maxval(abs(tmp_rk))
      maxserror = maxval(abs(ssfn-tmp_rk))
      write (out,"(/'Largest difference between analytical and numerical AO overlaps was ',g20.12)") maxserror
      maxserror = real(maxval(abs(ssf-tmp_rk)),kind=kind(maxserror))
      write (out,"( ' Largest difference between regular and quad-precision overlaps was ',g20.12)") maxserror
      write (out,"( '                                  Largest analytical AO overlap was ',g20.12/)") maxs
      if (maxserror>1e-4_rk*maxs) then
        write (out,"('Numerical integration accuracy is insufficient for sensible results. Threshold = ',g20.12)") &
               1e-4_rk * maxs
        call stop('static_tunnel%evaluate_1e_hamiltonian - bad numerical integration (2)')
      end if
      !
      maxserror = real(maxval(abs(ssf-transpose(ssf))),kind=kind(maxserror))
      write (out,"( ' Largest deviation of analytical AO overlap from index symmetry was ',g20.12/)") maxserror
      !
      ssf = 0.5_rk * (ssf + transpose(ssf))
      !
      if (verbose>=2) then
        write (out,"(/t5,'SPIN-FREE, FIELD-FREE 1-ELECTRON HAMILTONIAN'/)")
        call gamess_print_1e_integrals(hsf,bra=gam,ket=gam)
        write (out,"(/t5,'SPIN-FREE OVERLAP MATRIX'/)")
        call gamess_print_1e_integrals(ssf,bra=gam,ket=gam)
        write (out,"(/t5,'SPIN-FREE FIELD-INTERACTION (REAL)'/)")
        call gamess_print_1e_integrals(real(fsf,kind=rk),bra=gam,ket=gam)
        write (out,"(/t5,'SPIN-FREE CAP (COMPLEX SYMMETRIC)'/)")
        call gamess_print_1e_integrals(csf,bra=gam,ket=gam,symmetry='SYMMETRIC')
      end if
      !
      !  Expand 1-electron matrix elements to spin-AO basis
      !
      h0   = 0 ; h0f = 0 ; h0c = 0 ; smat = 0
      h0  (    1:nao     ,    1:nao     ) = hsf
      h0  (nao+1:nao_spin,nao+1:nao_spin) = hsf
      h0f (    1:nao     ,    1:nao     ) = fsf
      h0f (nao+1:nao_spin,nao+1:nao_spin) = fsf
      h0c (    1:nao     ,    1:nao     ) = csf
      h0c (nao+1:nao_spin,nao+1:nao_spin) = csf
      smat(    1:nao     ,    1:nao     ) = ssf
      smat(nao+1:nao_spin,nao+1:nao_spin) = ssf
      deallocate (hsf,ssf,lsfn,ssfn,fsf,tmp_xk,tmp_rk)
      !
      call TimerStop('1e Hamiltonian')
    end subroutine evaluate_1e_hamiltonian
    !
    subroutine density_matrix
      select case (fock_2e_math)
        case default
          call stop('static_tunnel%density_matrix - numerical mode '//trim(fock_2e_math)//' not implemented')
        case ('real')
          mo_occ_rk = real(mo_occ,kind=kind(mo_occ_rk))
          mos_rk    = cmplx(mos,kind=kind(mos_rk))
          call st_density_matrix(mo_occ_rk,mos_rk,rho_rk)
          rho       = cmplx(rho_rk,kind=kind(rho))
        case ('quad')
          call st_density_matrix(mo_occ,mos,rho)
      end select
    end subroutine density_matrix
    !
    subroutine g_matrix
      select case (fock_2e_math)
        case default
          call stop('static_tunnel%g_matrix - numerical mode '//trim(fock_2e_math)//' not implemented')
        case ('real')
          if (verbose>=0) then
            write (out,"(/'2e contributions are evaluates with kind-',i0,' density'/)") rk
          end if
          rho_rk = cmplx(rho,kind=rk)
          call fock_g_matrix(int2e,rho_rk,gmat_rk)
          gmat = gmat_rk
        case ('quad')
          if (verbose>=0) then
            write (out,"(/'2e contributions are evaluates with kind-',i0,' density'/)") xk
          end if
          call fock_g_matrix(int2e,rho,gmat)
      end select
    end subroutine g_matrix
    !
    !  Solve generalized eigenvalue problem
    !
    subroutine diagonalize_fmat
      select case (diag_math)
        case default; call stop('static_tunnel%diagonalize_fmat - bad diag_math')
        case ('quad')
          call st_diagonalize_fmat(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy,use_block_diag=use_block_diag)
        case ('real')
          smhalf_rk = real(smhalf,kind=kind(smhalf_rk))
          sphalf_rk = real(sphalf,kind=kind(sphalf_rk))
          fmat_rk   = cmplx(fmat,kind=kind(fmat_rk))
          call st_diagonalize_fmat(smhalf_rk,sphalf_rk,fmat_rk,nmo_null,real(eps_geev,kind=rk), &
                                   mos_rk,mo_energy_rk,use_block_diag=use_block_diag)
          mo_energy = cmplx(mo_energy_rk,kind=kind(mo_energy))
          mos       = cmplx(mos_rk,kind=kind(mos))
      end select
    end subroutine diagonalize_fmat
    !
    subroutine follow_mos
      integer(ik) :: nmo_act
      !
      if (follow_null) then
        nmo_act = nmo
        write (out,"('Following all MOs, including null-space')")
      else
        nmo_act = nmo - nmo_null
        write (out,"('Following only the MOs outside of the null-space')")
      end if
      !
      select case (diag_math)
        case default; call stop('follow_mos: bad diag_math')
        case ('real')
          sphalf_rk    = real(sphalf,kind=kind(sphalf_rk))
          mosg_rk      = cmplx(mosg,kind=kind(mosg_rk))
          mo_energy_rk = cmplx(mo_energy,kind=kind(mo_energy_rk))
          mos_rk       = cmplx(mos,kind=kind(mos_rk))
          call bt_follow_mos(nmo_act,real(mo_occ,kind=rk),real(eps_geev,kind=rk),sphalf_rk,mosg_rk,mo_energy_rk,mos_rk)
          mo_energy    = cmplx(mo_energy_rk,kind=kind(mo_energy))
          mos          = cmplx(mos_rk,kind=kind(mos))
        case ('quad')
          call bt_follow_mos(nmo_act,mo_occ,eps_geev,sphalf,mosg,mo_energy,mos)
      end select
    end subroutine follow_mos
    !
    !  Hartree-Fock electronic energy
    !
    subroutine total_energy
      complex(xk) :: efock, eg, elen

      efock = sum(rho * fmat)
      eg    = sum(rho * gmat)
      elen  = efock - 0.5_xk * eg
      escf  = elen + enuc
      if (verbose>=0) then
        write (out,"(/'Total SCF energy =',2(1x,g20.12)/)") escf
      end if
    end subroutine total_energy
    !
    !  Transform initial Fock matrix into the MO basis, giving an estimate
    !  of the eigenvalues
    !
    subroutine initial_eigenvalues
      integer(ik) :: imo
      !
      !$omp parallel do default(none) private(imo) shared(nmo,mo_energy,mos,fmat)
      fock_diagonal: do imo=1,nmo
        mo_energy(imo) = dot_product(mos(:,imo,1),matmul(fmat,mos(:,imo,2)))
      end do fock_diagonal
      !$omp end parallel do
      !
      where (mo_occ>0) mo_energy = mo_energy * mo_occ
      write (out,"(/t5,'Initial Fock matrix diagonal:')")
      write (out,"(5(1x,f12.6,2x,f12.6))") mo_energy
      write (out,"()")
    end subroutine initial_eigenvalues
    !
    subroutine scf_loop
      integer(ik) :: iter
      complex(xk) :: escf_old, de
      real(xk)    :: drho
      logical     :: converged
      !
      call TimerStart('SCF loop')
      converged = .false.
      repeat_scf: do iter=1,max_scf_iter
        rho_old  = rho
        fmat_old = fmat
        escf_old = escf
        !
        write (out,"('Beginning SCF cycle ',i4)") iter
        call flush(out)
        call density_matrix
        !
        gmat = 0
        if (.not.skip_2e) then
          call g_matrix
        end if
        fmat = hmat + gmat
        !
        if (verbose>=0 .and. iter==1) then
          call initial_eigenvalues
        end if
        !
        call total_energy
        call flush(out)
        !
        select case (scf_converger)
          case default
            write (out,"('SCF converger ',a,' is not recognized')") trim(scf_converger)
            call stop('static_tunnel%scf_loop - bad SCF converger choice')
          case ('mixing')
            if (iter>1) then
              fmat = scf_mixing * fmat + (1._rk-scf_mixing) * fmat_old
            end if
          case ('diis')
            select case (diag_math)
              case default; call stop('static_tunnel%scf_loop - bad diag_math')
              case ('real')
                smat_rk = real (smat,kind=kind(smat_rk))
                rho_rk  = cmplx(rho, kind=kind(rho_rk))
                fmat_rk = cmplx(fmat,kind=kind(fmat_rk))
                call diis_extrapolate(diis_st,iter,smat_rk,rho_rk,fmat_rk)
                fmat    = cmplx(fmat_rk,kind=kind(fmat))
              case ('quad')
                call diis_extrapolate(diis_st,iter,smat,rho,fmat)
            end select
        end select
        !
        call diagonalize_fmat
        call flush(out)
        !
        call follow_mos
        !
        !  Check for convergence; should only check after one full cycle!
        !
        if (iter<=1) cycle repeat_scf
        de   = escf - escf_old
        drho = maxval(abs(rho-rho_old))
        if (verbose>=0) then
          write (out,"('Iteration ',i4,' escf=',2(1x,g20.12),' de=',2(1x,g12.5),' drho= ',g12.5)") iter, escf, de, drho
          call flush(out)
        end if
        if (abs(de)<=scf_eps_energy .and. drho<=scf_eps_rho) then
          converged = .true.
          write (out,"('SCF converged for efield =',3(1x,f12.8))") efield_current
          write (out,"('Final SCF energy = ',g24.14,1x,g24.14)") escf
          call flush(out)
          exit repeat_scf
        end if
      end do repeat_scf
      !
      if (.not.converged) then
        call stop('static_tunnel%scf_loop - SCF convergence failure')
      end if
      !
      call TimerStop('SCF loop')
      call TimerReport
    end subroutine scf_loop
    !
    subroutine sort_scf_eigenvalues
      integer(ik), allocatable :: order(:)
      integer(ik)              :: alloc
      integer(ik)              :: nmo_act
      !
      call TimerStart('Sort SCF eigenvalues')
      !
      write (out,"(/'WARNING: Sorting SCF eigenvalues in the ascending order of Re(eps)')")
      write (out,"( 'WARNING: Correlation to the guess MOs will no longer be maintained'/)")
      !
      if (follow_null) then
        nmo_act = nmo
      else
        nmo_act = nmo - nmo_null
      end if
      !
      allocate (order(nmo_act),stat=alloc)
      if (alloc/=0) then
        call stop('static_tunnel%sort_scf_eigenvalues - allocation failed for order array')
      end if
      !
      call order_keys(real(mo_energy(:nmo_act),kind=rk),order)
      mo_energy(:nmo_act) = mo_energy(order)
      mos(:,:nmo_act,:)   = mos(:,order,:)
      !
      deallocate (order)
      call TimerStop('Sort SCF eigenvalues')
    end subroutine sort_scf_eigenvalues
    !
    subroutine fill_active_mos
      integer(ik) :: ipart, alloc
      integer(ik) :: nactive, imo
      integer(ik) :: non_null_nmo
      !
      non_null_nmo = nmo - nmo_null
      select case (correlation_type)
        case default
          write (out,"('static_tunnel%fill_active_mos: Correlation type ',a,' is not recognized')") trim(correlation_type)
          call stop('static_tunnel%fill_active_mos - bad correlation_type')
        case ('mp2','mp2 final')
          allocate (active_mos(non_null_nmo),stat=alloc)
          if (alloc/=0) then
            call stop('static_tunnel%fill_active_mos - allocation failed (A)')
          end if
          fill_active_mp2: do imo=1,non_null_nmo
            active_mos(imo) = imo
          end do fill_active_mp2
        case ('ormas','sz ormas')
          !
          ! Do a little bit of sanity checking; this is the first time we are using this data
          !
          if (ormas_npart<=0 .or. ormas_npart>max_ormas) then
            call stop('static_tunnel%fill_active_mos - bad max_ormas')
          end if
          if (any(ormas_orbitals(:ormas_npart)<=0)) then
            call stop('static_tunnel%fill_active_mos - non-positive ormas_orbitals')
          end if
          if (any(ormas_orbitals(:ormas_npart+1)>non_null_nmo+1)) then
            call stop('static_tunnel%fill_active_mos - ormas_orbitals too large')
          end if
          if (any(ormas_orbitals(2:ormas_npart+1)-ormas_orbitals(:ormas_npart)<=0)) then
            call stop('static_tunnel%fill_active_mos - empty partitions')
          end if
          !
          if (verbose>=0) then
            write (out,"(/'ORMAS partitions: ',i0,' CI electrons: ',i0/)") ormas_npart, ci_nelectrons
            write (out,"(1x,a3,2x,a5,1x,a5,2x,a4,1x,a4)") &
                   ' I ', 'MO_1', 'MO_n', 'MINE', 'MAXE', &
                   '---', '----', '----', '----', '----'
            print_ormas_partitions: do ipart=1,ormas_npart
              write (out,"(1x,i3,2x,i5,1x,i5,2x,i4,1x,i4)") &
                     ipart, ormas_orbitals(ipart), ormas_orbitals(ipart+1)-1, ormas_mine(ipart), ormas_maxe(ipart)
            end do print_ormas_partitions
            write (out,"()")
          end if
          !
          nactive = ormas_orbitals(ormas_npart+1)-ormas_orbitals(1)
          allocate (active_mos(nactive),stat=alloc)
          if (alloc/=0) then
            call stop('static_tunnel%fill_active_mos - allocation failed (B)')
          end if
          fill_active_ormas: do imo=ormas_orbitals(1),ormas_orbitals(ormas_npart+1)-1
            active_mos(imo-ormas_orbitals(1)+1) = imo
          end do fill_active_ormas
      end select
    end subroutine fill_active_mos
    !
    !  Calculate Sz for all MOs; thanks to orbital tracking, we do not expect
    !  the Sz value to change after this point (except if we are dealing with
    !  spin-orbit operators - which is an excercize for later).
    !
    subroutine fill_orbital_sz
      integer(ik) :: non_null_nmo
      integer(ik) :: alloc, imo
      complex(xk) :: tmp_l(nao_spin), tmp_r(nao_spin)
      complex(xk) :: sz
      !
      non_null_nmo = nmo - nmo_null
      if (allocated(orbital_sz)) deallocate (orbital_sz)
      allocate (orbital_sz(non_null_nmo),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i0,' allocating orbital_sz()')") alloc
        call stop('static_tunnel%fill_orbital_sz - allocation failed')
      end if
      !
      !  Evaluation of Sz below depends on sphalf separating into the
      !  alpha-alpha and beta-beta blocks, with no cross-spin terms.
      !  This happens to be true for our choice of basis set.
      !
      mo_loop: do imo=1,non_null_nmo
        tmp_l = matmul(sphalf,mos(:,imo,1))
        tmp_r(     :nao) =  0.5_xk * mos(     :nao,imo,2)
        tmp_r(nao+1:   ) = -0.5_xk * mos(nao+1:   ,imo,2)
        tmp_r = matmul(sphalf,tmp_r)
        sz    = real(sum(tmp_l*tmp_r),kind=rk)
        if (abs(aimag(sz))>1e3_rk*spacing(abs(sz))) then
          write (out,"('WARNING: Sz for orbital ',i0,' is not real. Sz = ',2g20.12)") imo, sz
        end if
        orbital_sz(imo) = real(sz,kind=rk) ! We are dropping the precision - linear dependencies should no longer happen
        write (out,"(' MO = ',i4,' Sz = ',f10.8,1x,f10.8,' eps = ',f14.10,1x,g14.7)") imo, sz, mo_energy(imo)
      end do mo_loop
    end subroutine fill_orbital_sz
    !
    subroutine read_mp2_reference
      integer(ik)  :: iu      ! Can be either input or iu_temp, depending on whether occ_file is blank or not
      integer(ik)  :: ios, alloc
      integer(ik)  :: imo, kmo
      integer(ik)  :: e_left  ! Number of electrons left to assign
      !
      allocate (mp2_ref_occ(nmo),stat=alloc)
      if (alloc/=0) call stop('static_tunnel%read_mp2_reference - allocation failed')
      !
      mp2_ref_occ = 0
      if (mp2_ref_file/='guess') then
        !
        !  We are given explicit occupation numbers; read them in.
        !
        iu = input
        if (mp2_ref_file/=' ') then
          iu = iu_temp
          open(iu,file=trim(mp2_ref_file),status='old',action='read',iostat=ios)
          if (ios/=0) then
            write (out,"('Error ',i0,' opening MP2 reference occupation file ',a)") ios, trim(mp2_ref_file)
            call stop('static_tunnel%read_mp2_reference - bad open')
          end if
        end if
        read (iu,*) mp2_ref_occ
        if (mp2_ref_file/=' ') close(iu)
        !
        !  Bit of sanity checking ...
        !
        if (any(mp2_ref_occ>1) .or. any(mp2_ref_occ<0)) then
          write (out,"('MP2 reference occupation numbers must be either 0 or 1')")
          write (out,"((20(1x,i3)))") mo_occ
          call stop('static_tunnel%read_mp2_reference - bad occupations (1)')
        end if
      else
        !
        !  No explicit occupations were given; use aufbau principle to assign occupations
        !
        e_left = ci_nelectrons
        if (e_left<=0) e_left = nint(nel_scf)
        if (e_left>size(mp2_ref_occ)) call stop('static_tunnel%read_mp2_reference - too many electrons')
        mp2_ref_occ(:e_left) = 1
      end if
      if (verbose<0) return
      !
      write (out,"(t5,a)") 'MP2 reference occupation numbers', &
                           '--------------------------------'
      print_occ: do imo=1,nmo,40
        kmo = min(imo+39,nmo)
        write (out,"(1x,i5,2x,40i2)") imo, mp2_ref_occ(imo:kmo)
      end do print_occ
      write (out,"()")
      write (out,"('Total number of electrons in MP2 reference = ',i0)") sum(mp2_ref_occ)
      write (out,"('               Charge of the MP2 reference = ',f0.5)") charge + nel_scf - sum(mp2_ref_occ)
      write (out,"()")

    end subroutine read_mp2_reference
    !
    subroutine read_ci_state_guess
      integer(sik) :: det(ci%nelectrons)
      complex(rk)  :: amp
      integer(ik)  :: ios
      integer(hik) :: irec, bra
      !
      call TimerStart('CI Read guess')
      ci_vector_guess = 0
      if (ci_vector_file/=' ') then
        open(iu_temp,file=trim(ci_vector_file),status='old',action='read',iostat=ios)
        if (ios/=0) then
          write (out,"('Error ',i0,' opening CI vector guess file ',a)") ios, trim(ci_vector_file)
          call stop('static_tunnel%read_ci_state_guess - bad open')
        end if
        irec = 0
        read_determinants: do
          irec = irec + 1
          read (iu_temp,*,iostat=ios) det, amp
          if (ios<0) exit read_determinants ! End-of-file
          if (ios>0) then
            write (out,"('Error ',i0,' reading record ',i0,' of ',a)") ios, irec, trim(ci_vector_file)
            call stop('static_tunnel%read_ci_state_guess - bad data')
          end if
          find_determinant: do bra=1,ci%ndets
            if (any(ci%dets(:,bra)/=det)) cycle find_determinant
            ci_vector_guess(bra) = amp
            cycle read_determinants
          end do find_determinant
          write (out,"('Determinant in record ',i0,' of ',a,' is not part of the CI')") irec, trim(ci_vector_file)
          write (out,"((20(1x,i5)))") det
          write (out,"('Amplitude = ',2g24.14)") amp
          call stop('static_tunnel%read_ci_state_guess - bad determinant')
        end do read_determinants
        close (iu_temp)
        write (out,"('Loaded in ',i0,' CI amplitudes from ',a)") irec-1, trim(ci_vector_file)
      end if
      call TimerStop('CI Read guess')
    end subroutine read_ci_state_guess
    !
    subroutine initialize_correlation
      integer(ik) :: alloc
      !
      call TimerStart('Initialize correlation')
      !
      select case (correlation_type)
        case ('mp2','mp2 final')
          call read_mp2_reference
          call fill_active_mos
        case default
          call fill_orbital_sz
          call fill_active_mos
          call ci_initialize_determinants(ci,mode=correlation_type,norbs=size(active_mos),nelectrons=ci_nelectrons, &
                 occ_ref=real(mo_occ(active_mos),kind=rk), &
                 fix_sz=ci_sz,orbital_sz=orbital_sz(active_mos), &
                 ormas_orbitals=ormas_orbitals(:ormas_npart+1)-ormas_orbitals(1)+1, &
                 ormas_mine=ormas_mine(:ormas_npart),ormas_maxe=ormas_maxe(:ormas_npart))
          if (ci_solver=='arpack') then
            allocate (ci_vector(ci%ndets),ci_vector_guess(ci%ndets),stat=alloc)
            if (alloc/=0) call stop('static_tunnel%initialize_correlation - allocation failure')
            call read_ci_state_guess
          end if
      end select
      call TimerStop('Initialize correlation')
      call TimerReport
    end subroutine initialize_correlation
    !
    subroutine finalize_correlation
      call TimerStart('Finalize correlation')
      select case (correlation_type)
        case ('mp2','mp2 final')
        case default
          call ci_destroy(ci)
      end select
      call TimerStop('Finalize correlation')
    end subroutine finalize_correlation
    !
    subroutine electron_correlation(final)
      logical, intent(in) :: final  ! Full electric field strength has been reached; this is the final call
      !
      logical     :: success
      complex(rk) :: final_e
      !
      call TimerStart('Electron correlation')
      !
      ! We only need transformed integrals over active MOs
      !
      correlate: select case (correlation_type)
        case ('mp2','mp2 final')
          if (final .or. correlation_type=='mp2') then
            select case (mp2_mode)
              case default
                write (out,"('mp2_mode ',a,' is not recognized')") trim(mp2_mode)
                call stop('static_tunnel%electron_correlation - bad mp2_mode')
              case ('incore')
                write (out,"(/'Performing MP2 using ',a,' math with ',a,' AO integrals and ',a,' MO integrals.'/)") &
                       trim(mp2_math), trim(ints_2e_math), trim(mp2_storage)
                select case (mp2_math)
                  case default
                    call stop('static_tunnel%electron_correlation - bad mp2_math '//trim(mp2_math))
                  case ('real')
                    call mp2_energy_incore(      mp2_ref_occ(active_mos),                    &
                                           real (mo_occ(active_mos),kind=rk),                &
                                           cmplx(escf,kind=rk),                              &
                                           cmplx(mo_energy(active_mos),kind=rk),             &
                                                 int2e,                                      &
                                           cmplx(mos(:,active_mos,:),kind=rk),               &
                                           real(pack(mp2_isa_de,mp2_isa_de>=0._rk),kind=rk), &
                                           storage_mode=mp2_storage)
                  case ('quad')
                    call mp2_energy_incore(      mp2_ref_occ(active_mos),                    &
                                           real (mo_occ(active_mos),kind=xk),                &
                                           cmplx(escf,kind=xk),                              &
                                           cmplx(mo_energy(active_mos),kind=xk),             &
                                                 int2e,                                      &
                                           cmplx(mos(:,active_mos,:),kind=xk),               &
                                           real(pack(mp2_isa_de,mp2_isa_de>=0._rk),kind=xk), &
                                           storage_mode=mp2_storage)
                end select
            end select
          end if
        case default
          select case (ci_solver)
            case default
              write (out,"('CI solver ',a,' is not recognized.')") trim(ci_solver)
              call stop('static_tunnel%electron_correlation - bad ci_solver')
          ! case ('direct')
          !   call transform_moint2e(int2e,moint_mode, &
          !           mos(:,active_mos,1),mos(:,active_mos,2),mos(:,active_mos,1),mos(:,active_mos,2),moint2e,iu_2e_mo)
          !   call ci_direct_diagonalization(ci,escf,mo_energy(active_mos),moint2e)
          !   call destroy_moint2e(moint2e)
          ! case ('arpack')
          !   call transform_moint2e(int2e,moint_mode, &
          !           mos(:,active_mos,1),mos(:,active_mos,2),mos(:,active_mos,1),mos(:,active_mos,2),moint2e,iu_2e_mo)
          !   call ci_build_hamiltonian(ci,mo_energy(active_mos),moint2e)
          !   call ci_diagonalization_arpack(ci,eref=escf,guess_e=ci_energy,guess_c=ci_vector_guess, &
          !           eps_root=ci_eps_root,max_iter=ci_max_iter,n_root=ci_n_root, &
          !           success=success,final_e=final_e,final_c=ci_vector)
          !   if (.not.success) call stop('static_tunnel%electron_correlation - CI convergence failure')
          !   ci_energy       = final_e
          !   ci_vector_guess = ci_vector
          !   write (out,"('Final CI energy is ',f24.14,1x,e24.14)") ci_energy
          !   write (out,"('Using final CI state vector as stating guess for the next CI')")
          !   call destroy_moint2e(moint2e)
          end select
      end select correlate
      call TimerStop('Electron correlation')
      call TimerReport
    end subroutine electron_correlation
    !
    !  We'll always punch the latest set of converged MOs; this way, we can restart
    !  calculation from the latest point it got to, if necessary.
    !
    subroutine punch_final_mos(ifield)
      integer(ik), intent(in) :: ifield
      integer(ik)             :: ios, iat, imo
      real(rk)                :: xyz(3), q
      !
      if (mos_file==' ') return
      !
      call TimerStart('Punch final MOs')
      open(iu_mos,file=trim(mos_file),form='formatted',action='write',status='replace',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' opening file ',a,' for the final MOs.')") ios, trim(mos_file)
        call stop('static_tunnel%punch_final_mos - Can''t create file')
      end if
      if (ifield==0) then
        write (iu_mos,"('Field-free MOs for ',a)") trim(molecule_file)
      else
        write (iu_mos,"('Field-perturbed MOs for ',a)") trim(molecule_file)
        write (iu_mos,"('External electric field increment ',i0,' scale ',f20.10)") ifield, efield_scale(ifield)
        write (iu_mos,"('External electric field F= ',3(1x,f12.8))") efield * efield_scale(ifield)
      end if
      write (iu_mos,"('Nuclear coordinates (atomic units): ')")
      write (iu_mos,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,natoms
        xyz = real(gam%atoms(iat)%xyz,kind=rk)
        q   = real(gam%atoms(iat)%znuc,kind=rk)
        write (iu_mos,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") q, xyz/abohr
      end do print_atoms
      write (iu_mos,"('Final energy (complex!) =',2(1x,f25.15))") escf
      write (iu_mos,"('Final orbital eigenvalues (complex!):')")
      write (iu_mos,"(10(1x,f12.6))") mo_energy
      !
      write (iu_mos,"(' Complex spin-orbitals follow for each MO index, left eigenvector')")
      write (iu_mos,"(' is punched before the right eigenvector.')")
      write (iu_mos,"(' Calculation was done using real kind ',i0)") xk
      write (iu_mos,"(' Eigenvectors are converted to real kind ',i0,' for export')") rk
      write (iu_mos,"(' $CVEC  ')")
      punch_mos_1: do imo=1,nao_spin
        call punch_vector_complex(iu_mos,2*imo-1,cmplx(mos(:,imo,1),kind=rk))  ! Left vector
        call punch_vector_complex(iu_mos,2*imo-0,cmplx(mos(:,imo,2),kind=rk))  ! Right vector
      end do punch_mos_1
      write (iu_mos,"(' $END   ')")
      !
      close (iu_mos)
      call TimerStop('Punch final MOs')
    end subroutine punch_final_mos
    !
    subroutine choose_accuracy
      select case (all_math)
        case default
          write (out,"('Global accuracy override all_math=',a,' is not recognized')") trim(all_math)
        case ('as-is')
        case ('real')
          ints_2e_math = 'real'
          fock_2e_math = 'real'
          diag_math    = 'real'
          mp2_math     = 'real'
          mp2_storage  = 'real'
        case ('quad')
          ints_2e_math = 'quad'
          fock_2e_math = 'quad'
          diag_math    = 'quad'
          mp2_math     = 'quad'
          mp2_storage  = 'quad'
      end select
      !
      if (verbose>-1) then
        write (out,"(/'    2-electron AO integrals precision: ',a)")  trim(ints_2e_math)
        write (out,"( '   Fock matrix construction precision: ',a)")  trim(fock_2e_math)
        write (out,"( '            Diagonalization precision: ',a)")  trim(diag_math)
        write (out,"( 'MP2 integral transformation precision: ',a)")  trim(mp2_math)
        write (out,"( '   MP2 MO integrals storage precision: ',a/)") trim(mp2_storage)
      end if
    end subroutine choose_accuracy
    !
    subroutine start
      integer(ik) :: info, ifield
      !
      call TimerStart('start')
      !
      read (input,nml=tunnel_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=tunnel_data)
      write (out,"()")
      !
      call choose_accuracy
      !
      call load_gamess_data
      !
      call allocate_dynamic_data
      !
      call prepare_ecps
      !
      call fill_occupations
      !
      call fetch_initial_mos
      !
      call TimerReport
      call flush(out)
      !
      write (out,"(/'Preparing the integrals'/)")
      call flush(out)
      !
      call evaluate_1e_hamiltonian
      !
      !  st_invert_smat() will also condition smat for linear dependencies
      !
      call st_invert_smat(nmo_null,smat,smhalf,sphalf,use_block_diag=use_block_diag,eps_smat=eps_smat)
      call flush(out)
      !
      if (.not.skip_2e) call prepare_2e(int2e,gam,scf_type,iu_2e_ao,iosize_2e,ints_math=ints_2e_math)
      !
      mos  = mosg
      if (skip_fieldfree) then
        write (out,"(/'WARNING: Using guess orbitals to stabilize calculation in the field'/)")
      else
        write (out,"(/'Starting field-free SCF iterations'/)")
        hmat = h0
        if (absorb_fieldfree) hmat = hmat + h0c
        efield_current = 0
        call scf_loop
        if (scf_sort_mos=='field-free') then
          call sort_scf_eigenvalues
        end if
        mos0 = mos ; mo_energy0 = mo_energy
        call punch_final_mos(0_ik)
        write (out,"(/'Using converged field-free orbitals as the new reference'/)")
      end if
      if (correlation_type/='none') then 
        !
        !  Correlation part can benefit from having converged field-free MOs;
        !  delay CI initialization until these are available
        !
        call initialize_correlation
        if (.not.skip_fieldfree) call electron_correlation(.false.)
      end if
      !
      !
      scan_through_fields: do ifield=1,cnt_efield_scale
        write (out,"(/'Starting SCF iterations in the external field. Scale factor = ',f20.10/)") efield_scale(ifield)
        mosg = mos
        hmat = h0 + h0c + h0f * efield_scale(ifield)
        efield_current = efield * efield_scale(ifield)
        call scf_loop
        mosf = mos ; mo_energyf = mo_energy
        call punch_final_mos(ifield)
        if (correlation_type/='none') call electron_correlation(ifield==cnt_efield_scale)
      end do scan_through_fields
      !
      if (.not.skip_2e) call clear_2e(int2e)
      if (correlation_type/='none') call finalize_correlation
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
    !
    subroutine stop(message)
      character(len=*), intent(in) :: message
      !
      write (out,"('STOP: ',a)") trim(message)
      write (0,"('STOP: ',a)") trim(message)
      call flush(out)
      call flush(0)
      stop 'error in static_tunnel.f90'
    end subroutine stop
  end module static_tunnel
!
  program driver
    use static_tunnel
    use accuracy
    use math

    real(rk) :: hole
    real(xk) :: hole_xk

    write (out,"('Version: ',a/)") __BUILD_ID__
    !
    call accuracyInitialize
    !
    !  We have to issue calls below in case we run under OpenMP.
    !  See comments in math.f90
    !
    hole    = MathFactorial(80)
    hole_xk = MathFactorial(80,1._xk)
    hole    = MathDoubleFactorial(80)
    hole_xk = MathDoubleFactorial(80,1._xk)
    hole    = MathLogFactorial(80)
    hole_xk = MathLogFactorial(80,1._xk)
    !
    !$ write (out,"(/'WARNING: Some compilers are known to miscompile OpenMP constructs in this program')")
    !$ write (out,"( 'WARNING: Please make sure the compiler you are using yields correct results!'/)")
    !
    call start
  end program driver
