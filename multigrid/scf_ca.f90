!
!  This is an example of a very simplistic closed-shell single-determinant 
!  SCF code, using Nafie's Complete Adiabatic Hamiltonian. The only observable
!  we calculate here is the electronic current density.
!
!  It is not intended as a stand-alone electronic structure treatment: 
!  in all sane application scenarios, one would converge Born-Oppenheimer
!  SCF in GAMESS, then use this solution to work with the (very slightly
!  different) Complete Adiabatic SCF.
!
!  The only property we are currently interested in are the electronic
!  currents.
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
!  2013 Feb 04 - Added ECP support using resolution into projectors 
!                This is intended to test ECP implementation in multigrid 
!                more than anything else.
!
  module scf_ca
    use accuracy
    use timer
    use math
    use import_gamess
    use gamess_internal
    use ecp_convert_gamess
    use opendx
    use sort_tools
    use lapack
    implicit none
    private
    public start
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: iu_velo = 35          ! I/O unit used for fetching atomic velocities data
    integer(ik), parameter :: iu_mos  = 36          ! I/O unit used for punching final MOs
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose        = 1            ! Level of output
    character(len=clen) :: molecule_file  = 'mol.dat'    ! Name of GAMESS checkpoint file, containing
                                                         ! the structure and converged orbitals
    character(len=clen) :: velocity_file  = 'vel.dat'    ! Name of a free-format file containing Cartesian
                                                         ! velocities of all atoms in the system
    real(rk)            :: charge         = 0            ! Overall molecular charge; non-integer values 
                                                         ! will trigger half-electron treatment
    integer(ik)         :: max_scf_iter   = 20           ! Max number of SCF iterations
    integer(ik)         :: max_ca_iter    = 15           ! Max number of double-SCF iterations
    real(rk)            :: scf_eps_energy = 1e-7_rk      ! Desired SCF convergence for the total energy
    real(rk)            :: scf_eps_rho    = 1e-7_rk      ! Desired SCF convergence for the density matrix
    real(rk)            :: ca_eps_energy  = 1e-6_rk      ! Desired double-SCF convergence for the total energy
    real(rk)            :: ca_eps_rho     = 1e-3_rk      ! Desired double-SCF convergence for the density matrix
    real(rk)            :: accuracy_2e    = 1e-10_rk     ! Accuracy required of 2-e contributions. Only applies
                                                         ! for scf_type = 'direct'
    character(len=clen) :: scf_type       = 'auto'       ! Can be one of:
                                                         ! 'direct'        - compute all integrals on the fly
                                                         ! 'in-memory'     - keep integrals in memory for each SCF cycle
                                                         ! 'all in-memory' - maintain integrals in memory between cycles
                                                         ! 'auto'          - choose 'all in-memory' for small molecules, 
                                                         !                   and 'direct' otherwise
    real(rk)            :: scf_threshold  = 8._rk        ! Memory which can be used for integral storage, in Gbytes
                                                         !
    real(rk)            :: epsd_include   = 1e-4_rk      ! Orbital eigenvalue degeneracy tolerance for inclusion
    real(rk)            :: epsd_tolerate  = 0.1_rk       ! Orbital eigenvalue degeneract tolerance for mixing-coeff.
                                                         ! inclusion
    real(rk)            :: mix_degenerate = 0.1_rk       ! Threshould for mixing coefficient to be considered a
                                                         ! degenerate rotation
    logical             :: use_ca         = .true.       ! Include CA terms in the Hamiltonian
    real(rk)            :: ca_time        = 0.1_rk       ! Time differerence used to calculate CA term in the
                                                         ! Fock operator.
    real(rk)            :: hca_mixing     = 0.3_rk       ! Mixing factor used for the old and new HCA, to improve
                                                         ! stability of the double-SCF process
    character(len=clen) :: mos_file       = 'ca_mos.dat' ! Quasi-GAMESS punch file for molecular orbitals
                                                         ! Appending this file to 'mol.dat' will make it
                                                         ! readable by our OpenDX routines
                                                         ! Blank disables the output.
    character(len=clen) :: odx_file       = 'ca_grid.dx' ! OpenDX formatted file to write out the grid cube
    real(rk)            :: grid_origin(3)   = (/-5._rk, -5._rk, -5._rk /)
                                                         ! Origin of the 3D uniform grid used to evaluate the currents
    real(rk)            :: grid_step  (3,3) = reshape( &
                                              (/ 1._rk, 0._rk, 0._rk, &
                                                 0._rk, 1._rk, 0._rk, &
                                                 0._rk, 0._rk, 1._rk /), & 
                                              (/ 3, 3 /))
                                                         ! Grid increments
    integer(ik)         :: grid_counts(3)   = (/ 10_ik, 10_ik, 10_ik /)
                                                         ! Number of grid points along each increment
    real(rk)            :: grid_dx          = 1e-5_rk    ! Numerical differentiation step; used for div j evaluation
    logical             :: print_grid       = .false.    ! Force verbose printout of current distribution on grid
    !
    character(len=clen) :: ecp_file        = ' '         ! File containing the ECP and the RI basis; blank
                                                         ! means use ECP and basis in mol_file.
                                                         ! WARNING: standard molecular basis sets are likely 
                                                         !          not flexible enough to match GAMESS ECP energies!
    real(rk)            :: ecp_eps_min     = 1e-6_rk     ! Small shift value cut-off (absolute) in ECP
    real(rk)            :: ecp_eps_max     = 1e+6_rk     ! Large shift value cut-off (positive) in ECP
    real(rk)            :: ecp_eps_grid    = 1e-6_rk     ! Characteristic grid spacing cut-off in ECP
    character(len=clen) :: ecp_report      = ' '         ! File to report ECP level-shift projectors to
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type(gam_structure), target :: gamd              ! Structure descriptor used at a dispaced geometry
    type(gam_structure), target :: gam0              ! Structure descriptor used at a central geometry
    type(ecp_molecule), target  :: ecpd              ! ECP used at displaced geometry
    type(ecp_molecule), target  :: ecp0              ! ECP used at the central geometry
    integer(ik)                 :: natoms            ! Number of atoms in the structure
    integer(ik)                 :: nbasis            ! Number of basis functions.
    real(rk)                    :: nelectrons        ! Number of electrons in the molecule; does not have
                                                     ! to be an integer - although only integer numbers are
                                                     ! physical. Note that the CA term us only defined for
                                                     ! even values of nelectrons at the moment.
    integer(ik)                 :: nocc              ! Number of occupied MOs
    real(rk), allocatable       :: mo_energy(:)      ! Energies of molecular orbitals
    real(rk), allocatable       :: mo_energy0(:)     !   ... at the central geometry (degeneracy checking)
    real(rk), allocatable       :: mo_energyd(:)     !   ... at the displaced geometry
    real(rk), allocatable       :: mo_occ   (:)      ! Occupation numbers of molecular orbitals
    complex(rk), allocatable    :: mos (:,:)         ! Molecular orbitals at the currentry active geometry
    complex(rk), allocatable    :: mos0(:,:)         ! Molecular orbitals at the central geometry
    complex(rk), allocatable    :: mosd(:,:)         ! Molecular orbitals at the displaced geometry
    complex(rk), allocatable    :: rho (:,:)         ! Electronic density matrix (AO basis)
    complex(rk), allocatable    :: rho_old(:,:)      ! Electronic density matrix from previous iteration
    complex(rk), allocatable    :: rho0(:,:)         ! Electronic density matrix from previous CA iteration
    real(rk), allocatable       :: smat(:,:)         ! Overlap matrix (AO basis)
    real(rk), allocatable       :: ddrmat(:,:,:)     ! Matrix elements of d/dx, d/dy, d/dz (the last index is
                                                     ! the Cartesian component of the gradient)
    real(rk), allocatable       :: shalf(:,:)        ! S^{-1/2}
    real(rk), allocatable       :: hbo (:,:)         ! Born-Oppenheimer part of 1-electron Hamiltonian
    real(rk), allocatable       :: hecp(:,:)         ! ECP contribution to 1-electron Hamiltonian
    complex(rk), allocatable    :: hca (:,:)         ! "Complete adiabatic" contribution to 1-electron Hamiltonian
    complex(rk), allocatable    :: hca_old(:,:)      ! "old" HCA, used for stability mixing
    complex(rk), allocatable    :: gmat(:,:)         ! 2-electron contribution to the Fock matrix
    complex(rk), allocatable    :: fmat(:,:)         ! Fock matrix
    real(rk), allocatable       :: ca_ao_smat(:,:)   ! AO overlap matrix, central and displaced geometry
    complex(rk), allocatable    :: ca_mo_smat(:,:)   ! MO overlap matrix, central and displaced geometry
    real(rk)                    :: enuc              ! Nuclear repulsion energy
    real(rk)                    :: etot              ! Total energy
    real(rk)                    :: energy_ca         ! "Complete adiabatic" contribution to the energy
    real(rk), allocatable       :: vatom(:,:)        ! Atomic velocities, in Bohr/au[t]. The first index is the
                                                     ! Cartesian velocity component; the second index is the atomic
                                                     ! index, which must match the ordering in the GAMESS input file.
    integer(ik)                 :: grid_points       ! Number of property grid points
    real(rk), allocatable       :: grid_xyz(:,:)     ! Coordinates of grid points
    real(rk), allocatable       :: grid_property(:,:)! Calculated properties at grid points. The first index is:
                                                     ! 0     = electron density rho
                                                     ! 1,2,3 = components of the electron current j
                                                     ! 4     = \Nabla . j  a.k.a. time derivative of rho
    !
    !  Variables related to 2-electron integral handling
    !
    type int2e_block
      real(rk), allocatable          :: ints(:,:,:,:)      ! Integrals cache
    end type int2e_block
    type int2e_cache
      integer(ik)                    :: nbatch_2e          ! Number of 2e batches
      integer(ik)                    :: maxbatch_2e        ! Largest 2e batch
      integer(ik), allocatable       :: batch_indices(:,:) ! Index 1: first orbital in batch; # of orbitals in batch; atom batch belongs to
                                                           ! Index 2: the batch 
      type(int2e_block), allocatable :: blocks(:,:,:,:)    ! Integral block cache
    end type int2e_cache
    !
    type(int2e_cache), pointer :: int2e => NULL()    ! Currently active integrals context
    type(int2e_cache), target  :: int2e0             ! Integrals cache at the target geometry
    type(int2e_cache), target  :: int2ed             ! Integrals cache at the displaced geometry
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /scf_ca_data/ verbose, use_ca, ca_time, &
                           molecule_file, velocity_file, mos_file, odx_file, &
                           charge, &
                           max_scf_iter, scf_eps_energy, scf_eps_rho, &
                           max_ca_iter, epsd_include, epsd_tolerate, mix_degenerate, &
                           ca_eps_energy, ca_eps_rho, hca_mixing, &
                           grid_origin, grid_step, grid_counts, grid_dx, print_grid, &
                           ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
                           accuracy_2e, scf_type, scf_threshold
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
        xyz = real(gam%atoms(iat)%xyz,kind=kind(xyz))
        q   = real(gam%atoms(iat)%znuc,kind=kind(xyz))
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") q, xyz/abohr, xyz
      end do print_atoms
      write (out,"()")
    end subroutine print_geometry
    !
    subroutine load_gamess_data
      integer(ik)  :: iat, jat
      real(rk)     :: xyz(3), q, r
      !
      !  Load the structure and copy some basic parameters to the global arrays
      !
      call gamess_load_orbitals(file=trim(molecule_file),structure=gam0)
      write (out,"(/'Loaded GAMESS checkpoint file ',a/)") trim(molecule_file)
      natoms     = gam0%natoms
      nbasis     = gam0%nbasis
      nelectrons = sum(real(gam0%atoms(:)%znuc,kind=kind(xyz))) - charge
      enuc       = 0.0_rk
      nuclear_repulsion: do iat=2,natoms
        xyz = real(gam0%atoms(iat)%xyz,kind=kind(xyz))
        q   = real(gam0%atoms(iat)%znuc,kind=kind(xyz))
        do jat=1,iat-1
          r = sqrt(sum((real(gam0%atoms(jat)%xyz,kind=kind(xyz))-xyz)**2)) / abohr
          enuc = enuc + q*real(gam0%atoms(jat)%znuc,kind=kind(xyz)) / r
        end do
      end do nuclear_repulsion
      !
      !  A bit of sanity checking
      !
      if (gam0%nvectors/=nbasis) then
        write (out,"('Number of MOs (',i0,') does not match the number of basis functions (',i0,')')") gam0%nvectors, nbasis
        stop 'scf_ca%load_gamess_data - incomplete'
      end if
      if (nelectrons<0.0_rk .or. nelectrons>2.0_rk*nbasis) then
        write (out,"('Number of electrons (',f0.5,') is strange.')") nelectrons
        stop 'scf_ca%load_gamess_data - bad electron count'
      end if
      if (verbose<0) return
      !
      !  Tell a bit more!
      !
      write (out,"('           Number of atoms = ',i0)") natoms
      write (out,"(' Number of basis functions = ',i0)") nbasis
      write (out,"('       Number of electrons = ',f0.5)") nelectrons
      write (out,"('              Total charge = ',f0.5)") charge
      call print_geometry(gam0)
      write (out,"('  Nuclear repulsion energy = ',f20.12)") enuc
      write (out,"()")
    end subroutine load_gamess_data
    !
    !  Load atomic velocities; needed for the construction of CA Hamiltonian
    !
    subroutine load_velocities
      integer(ik) :: iat, ios
      real(rk)    :: xyz(3), q
      !
      if (use_ca) then
        open(iu_velo,file=trim(velocity_file),form='formatted',action='read',position='rewind',iostat=ios)
        if (ios/=0) then
          write (out,"('Error ',i0,' opening velocities file ',a)") ios, trim(velocity_file)
          stop 'scf_ca%load_velocies - no velocities?!'
        end if
        velo_read: do iat=1,natoms
          read(iu_velo,*,iostat=ios) vatom(:,iat)
          if (ios/=0) then
            write (out,"('Error ',i0,' reading velocity of atom ',i0,' from ',a)") ios, iat, trim(velocity_file)
            stop 'scf_ca%load_velocities - error reading velocities'
          end if
        end do velo_read
        close (iu_velo)
      else
        vatom = 0
        write (out,"(/'WARNING: Complete Adiabatic terms are neglected; setting atomic velocities to zero'/)")
      end if
      !
      if (verbose<0) return
      !
      !  Tell a bit more!
      !
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Velocity (Bohr/au[time])  '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '------------------------  '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,natoms
        xyz = real(gam0%atoms(iat)%xyz,kind=kind(xyz))
        q   = real(gam0%atoms(iat)%znuc,kind=kind(xyz))
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") q, xyz/abohr, vatom(:,iat)
      end do print_atoms
      write (out,"()")
    end subroutine load_velocities
    !
    !  Calculation of the Complete Adiabatic contribution to the Fock operator requires 
    !  an additional SCF calculation at a displaced geometry. 
    !
    subroutine prepare_displaced_geometry
      integer(ik) :: iat
      real(rk)    :: rmax
      !
      !  First, load the struture we used for the central geometry
      !
      call gamess_load_orbitals(file=trim(molecule_file),structure=gamd)
      !
      !  Now, move all atoms along the velocity vector
      !
      write (out,"('Displacing the atoms to time ',f12.5,' au[t]')") ca_time
      rmax = 0
      move_atoms: do iat=1,natoms
        gamd%atoms(iat)%xyz = gam0%atoms(iat)%xyz + vatom(:,iat)*ca_time*abohr
        rmax = max(rmax,sum(vatom(:,iat)**2))
      end do move_atoms
      rmax = ca_time * sqrt(rmax)
      write (out,"('Max. coordinate displacement is ',f12.5,' Bohr')") rmax
      if (verbose<0) return
      !
      write (out,"(/t5,'Displaced geometry for the CA term')")
      call print_geometry(gamd)
    end subroutine prepare_displaced_geometry
    !
    !  Convert GAMESS ECPs to the projector form. This has to be done at both the
    !  original and the displaced geometries.
    !
    subroutine prepare_ecps
      type(gam_structure) :: gam
      integer(ik)         :: iat
      !
      call TimerStart('Prepare ECPs')
      if (ecp_file==' ') then
        if (any(gam0%atoms(:)%ecp_nterms>0)) then
          write (out,"(/'WARNING: Using molecular basis for ECP resolution of identity'/)")
        end if
        call ecp_convert(gam0,ecp0,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report)
        call ecp_convert(gamd,ecpd,ecp_eps_min,ecp_eps_max,ecp_eps_grid,report_file=' ')
      else
        if (any(gam0%atoms(:)%ecp_nterms>0)) then
          write (out,"('Loading ECP RI basis from ',a)") trim(ecp_file)
        end if
        call gamess_load_orbitals(file=ecp_file,structure=gam)
        if (gam%natoms/=gam0%natoms) stop 'ECP atom count mismatch'
        !
        forall (iat=1:gam%natoms) gam%atoms(iat)%xyz = gam0%atoms(iat)%xyz
        call ecp_convert(gam,ecp0,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report)
        !
        forall (iat=1:gam%natoms) gam%atoms(iat)%xyz = gamd%atoms(iat)%xyz
        call ecp_convert(gam,ecpd,ecp_eps_min,ecp_eps_max,ecp_eps_grid,report_file=' ')
        !
        call gamess_destroy(gam)
      end if
      call TimerStop('Prepare ECPs')
    end subroutine prepare_ecps
    !
    !  Memory allocated by allocate_dynamic_data() is never released; this is by design.
    !
    subroutine allocate_dynamic_data
      integer(ik) :: alloc
      !
      allocate (mo_energy(nbasis), mo_energy0(nbasis), mo_energyd(nbasis), mo_occ(nbasis), &
                mos(nbasis,nbasis), mos0(nbasis,nbasis), mosd(nbasis,nbasis), &
                rho(nbasis,nbasis), smat(nbasis,nbasis), hbo(nbasis,nbasis), &
                hecp(nbasis,nbasis), &
                gmat(nbasis,nbasis), hca(nbasis,nbasis), hca_old(nbasis,nbasis), &
                fmat(nbasis,nbasis), ddrmat(nbasis,nbasis,3), &
                shalf(nbasis,nbasis), rho_old(nbasis,nbasis), rho0(nbasis,nbasis), &
                ca_ao_smat(nbasis,nbasis), ca_mo_smat(nbasis,nbasis), &
                vatom(3,natoms), &
                stat=alloc)
      if (alloc/=0) then
        write (out,"('scf_ca%allocate_dynamic_data: Error ',i0,' allocating quadratic arrays. nbasis = ',i0)") &
               alloc, nbasis
        stop 'scf_ca%allocate_dynamic_data - out of memory'
      end if
    end subroutine allocate_dynamic_data
    !
    subroutine fill_occupations
      integer(ik) :: imo, kmo
      real(rk)    :: e_left  ! Number of electrons left to assign
      real(rk)    :: e_this  ! Number of electrons to go to the current MO
      !
      mo_occ = 0.0_rk 
      e_left = nelectrons
      nocc   = 0
      fill_occ: do imo=1,nbasis
        e_this = min(e_left,2.0_rk)
        mo_occ(imo) = e_this
        nocc = nocc + 1
        e_left = e_left - e_this
        if (e_left<=0.0_rk) exit fill_occ
      end do fill_occ
      if (e_left>0.0_rk) stop 'scf_ca%fill_occupations - unallocated electrons'
      if (verbose<0) return
      !
      write (out,"(t5,a)") 'MO occupation numbers', &
                           '---------------------'
      print_occ: do imo=1,nbasis,10
        kmo = min(imo+9,nbasis)
        write (out,"(1x,i5,2x,10(1x,f7.5))") imo, mo_occ(imo:kmo)
      end do print_occ
      write (out,"()")
    end subroutine fill_occupations
    !
    subroutine evaluate_ecps(gam,ecp)
      type(gam_structure), intent(inout), target :: gam
      type(ecp_molecule), intent(inout), target  :: ecp
      !
      integer(ik)                  :: iecp     ! Current ECP index
      integer(ik)                  :: nprj     ! Projector count in a ECP
      integer(ik)                  :: prj_nbas ! Number of basis functions in a projector
      integer(ik)                  :: alloc
      type(ecp_atom), pointer      :: aecp     ! Current ECP
      type(gam_structure), pointer :: prj      ! Projector orbitals for the current ECP
      integer(ik)                  :: iprj, ia, ib
      !
      real(rk), allocatable :: ecp_prj_ao(:,:), ecp_prj(:,:)
      !
      if (ecp%necps<=0) return
      !
      call TimerStart('ECP Hamiltonian')
      hecp = 0
      apply_ecps: do iecp=1,ecp%necps
        aecp     => ecp%ecps(iecp)
        nprj     =  size(aecp%vshift)
        prj      => aecp%projectors
        prj_nbas =  prj%nbasis
        allocate (ecp_prj_ao(nbasis,prj_nbas),ecp_prj(nbasis,nprj),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i0,' allocating buffers for ECP projectors')") alloc
          stop 'scf_ca%evaluate_ecps - memory allocation failed'
        end if
        !
        !  Evaluate projectors
        !
        call gamess_1e_integrals('AO OVERLAP',ecp_prj_ao,bra=gam,ket=prj)
        ecp_prj = matmul(ecp_prj_ao,prj%vectors(:,1:nprj))
        !
        !  Accumulate the ECP contribution; this can be done more efficiently
        !  using DSYR2K, but we won't bother - this routine is not on a critical path (?)
        !
        apply_projectors: do iprj=1,nprj
          aos_right: do ib=1,nbasis
            aos_left: do ia=1,nbasis
              hecp(ia,ib) = hecp(ia,ib) + aecp%vshift(iprj)*ecp_prj(ia,iprj)*ecp_prj(ib,iprj)
            end do aos_left
          end do aos_right
        end do apply_projectors
        !
        deallocate (ecp_prj_ao,ecp_prj)
      end do apply_ecps
      if (verbose>=3) then
        write (out,"(/t5,'ECP INTEGRALS'/)")
        call gamess_print_1e_integrals(hecp,bra=gam,ket=gam)
      end if
      hbo = hbo + hecp
      call TimerStop('ECP Hamiltonian')
    end subroutine evaluate_ecps
    !
    subroutine one_electron_hamiltonian(gam,ecp)
      type(gam_structure), intent(inout), target :: gam
      type(ecp_molecule), intent(inout), target  :: ecp
      integer(ik) :: iat
      real(rk)    :: xyz(3), q
      !
      call TimerStart('1e Hamiltonian')
      call gamess_1e_integrals('AO KINETIC',hbo,bra=gam,ket=gam)
      if (verbose>=3) then
        write (out,"(/t5,'KINETIC ENERGY INTEGRALS'/)")
        call gamess_print_1e_integrals(hbo,bra=gam,ket=gam)
      end if
      !
      !  Fill nuclear attraction integrals; the overlap matrix has not
      !  been initialized yet, so that we can use that as a scratch area 
      !
      nuclear_attraction: do iat=1,natoms
        xyz = real(gam%atoms(iat)%xyz,kind=kind(xyz)) / abohr
        q   = real(gam%atoms(iat)%znuc,kind=kind(xyz))
        call gamess_1e_integrals('AO 3C 1/R',smat,bra=gam,ket=gam,op_xyz=xyz)
        smat = -q * smat
        if (verbose>=4) then
          write (out,"(/t5,'NUCLEAR ATTRACTION TO ATOM ',i3,' Z= ',f12.5,' XYZ= ',3f12.5/)") iat, q, xyz
          call gamess_print_1e_integrals(smat,bra=gam,ket=gam)
        end if
        hbo = hbo + smat
      end do nuclear_attraction
      !
      call evaluate_ecps(gam,ecp)
      !
      call gamess_1e_integrals('AO OVERLAP',smat,bra=gam,ket=gam)
      shalf = smat
      call lapack_ginverse(shalf,0.25_rk) ! power of 1/4 gives S^{-1/2}
      !
      if (verbose>=2) then
        write (out,"(/t5,'BORN-OPPENHEIMER 1-ELECTRON HAMILTONIAN'/)")
        call gamess_print_1e_integrals(hbo,bra=gam,ket=gam)
        write (out,"(/t5,'OVERLAP MATRIX'/)")
        call gamess_print_1e_integrals(smat,bra=gam,ket=gam)
      end if
      if (verbose>=3) then
        write (out,"(/t5,'S**-0.5 MATRIX'/)")
        call gamess_print_1e_integrals(shalf,bra=gam,ket=gam)
      end if
      !
      call TimerStop('1e Hamiltonian')
    end subroutine one_electron_hamiltonian
    !
    !  The subroutine below is correct ONLY for the case of uniform motion:
    !  It neglects the orbital relaxation contribution in the nuclear derivative.
    !
    subroutine complete_adiabatic_hamiltonian_guess
      integer(ik) :: iat, jat            ! Atom indices
      integer(ik) :: ib,  jb             ! Batch indices
      integer(ik) :: ip1, ipn, jp1, jpn  ! Orbital ranges
      integer(ik) :: ic                  ! Cartesian component
      real(rk)    :: hmax, hdev          ! Largest matrix element and max deviation from Hermiticity
      !
      call TimerStart('CA Hamiltonian guess')
      !
      !  We'll need matrix elements of nuclear gradient operator. These happen
      !  to be the same as electronic gradient, taken with negative sign
      !
      call gamess_1e_integrals('AO D/DX',ddrmat(:,:,1),bra=gam0,ket=gam0)
      call gamess_1e_integrals('AO D/DY',ddrmat(:,:,2),bra=gam0,ket=gam0)
      call gamess_1e_integrals('AO D/DZ',ddrmat(:,:,3),bra=gam0,ket=gam0)
      !
      !  We'll also need the (closely related) oveerlap integrals between two geometries later on
      !
      call gamess_1e_integrals('AO OVERLAP',ca_ao_smat,bra=gam0,ket=gamd)
      !
      !
      if (verbose>=3) then
        write (out,"(/t5,'(D/DX)'/)")
        call gamess_print_1e_integrals(ddrmat(:,:,1),bra=gam0,ket=gam0,symmetry='ANTI-HERMITIAN')
        write (out,"(/t5,'(D/DY)'/)")
        call gamess_print_1e_integrals(ddrmat(:,:,2),bra=gam0,ket=gam0,symmetry='ANTI-HERMITIAN')
        write (out,"(/t5,'(D/DZ)'/)")
        call gamess_print_1e_integrals(ddrmat(:,:,3),bra=gam0,ket=gam0,symmetry='ANTI-HERMITIAN')
        write (out,"(/t5,'DISPLACED AO OVERLAP MATRIX'/)")
        call gamess_print_1e_integrals(ca_ao_smat,bra=gam0,ket=gamd,symmetry='NONE')
      end if
      !
      !  Build the CA 1-electron term. For each nucleus, we must multiply by
      !  the nuclear velocity. This takes into account the change of sign 
      !  due to substitution of the electronic-coordinate derivative.
      !
      !  We'll use data structures prepared for 2-e integrals here, so we
      !  need to initialize the tables
      !
      call prepare_2e(gam0,int2e0)
      hca = 0
      right_shells: do jb=1,int2e%nbatch_2e
        jp1 = int2e%batch_indices(1,jb)
        jpn = jp1 + int2e%batch_indices(2,jb) - 1
        jat = int2e%batch_indices(3,jb)
        left_shells: do ib=1,int2e%nbatch_2e
          ip1 = int2e%batch_indices(1,ib)
          ipn = ip1 + int2e%batch_indices(2,ib) - 1
          iat = int2e%batch_indices(3,ib)
          do ic=1,3
            hca(ip1:ipn,jp1:jpn) = hca(ip1:ipn,jp1:jpn) + (0._rk,1.0_rk) * ddrmat(ip1:ipn,jp1:jpn,ic) * vatom(ic,jat)
          end do
        end do left_shells
      end do right_shells
      call clear_2e
      !
      if (verbose>=2) then
        write (out,"(/t5,'(GUESS) COMPLETE ADIABATIC CONTRIBUTION BEFORE HERMITIAN PROJECTION'/)")
        call gamess_print_1e_integrals(hca,bra=gam0,ket=gam0)
      end if
      !
      !  We require the Hermitian part of the matrix element
      !
      hdev = 0.5_rk * maxval(abs(hca-transpose(conjg(hca))))
      hca  = 0.5_rk * (hca + transpose(conjg(hca)))
      hca_old = hca
      hmax = maxval(abs(hca))
      if (verbose>=0) then
        write (out,"('     Largest matrix element magnitude of the Hermitian part of guess H(CA): ',g14.6)") hmax
        write (out,"('Largest matrix element magnitude of the anti-Hermitian part of guess H(CA): ',g14.6/)") hdev
      endif
      !
      if (verbose>=2) then
        write (out,"(/t5,'(GUESS) COMPLETE ADIABATIC CONTRIBUTION TO 1-ELECTRON HAMILTONIAN'/)")
        call gamess_print_1e_integrals(hca,bra=gam0,ket=gam0)
      end if
      !
      call TimerStop('CA Hamiltonian guess')
    end subroutine complete_adiabatic_hamiltonian_guess
    !
    !  Initialize basic data on 2e integrals; we'll pull (and possibly cache) the
    !  actual integrals later on
    !
    subroutine prepare_2e(gam,ints)
      type(gam_structure), intent(inout), target :: gam
      type(int2e_cache), intent(inout), target   :: ints
      integer(ik) :: info_2e(4), alloc, ib
      real(rk)    :: mem_size
      !
      call TimerStart('Prepare 2E Buffer')
      int2e => ints
      !
      if (.not.allocated(int2e%batch_indices)) then
        call gamess_2e_info('GET INFO',gam=gam,op_iparam=info_2e)
        int2e%nbatch_2e   = info_2e(1)
        int2e%maxbatch_2e = info_2e(2)
        !
        if (verbose>=2) then
          write (out,"('         Number of 2e blocks: ',i4)") int2e%nbatch_2e
          write (out,"('Orbitals in largest 2e block: ',i4)") int2e%maxbatch_2e
        end if
        !
        allocate (int2e%batch_indices(3,int2e%nbatch_2e),stat=alloc)
        if (alloc/=0) stop 'scf_ca%prepare_2e - no memory for batch table'
        !
        remember_batches: do ib=1,int2e%nbatch_2e
          info_2e(1) = ib
          call gamess_2e_info('GET BATCH INFO',gam=gam,op_iparam=info_2e)
          int2e%batch_indices(1,ib) = info_2e(3)  ! Initial position
          int2e%batch_indices(2,ib) = info_2e(2)  ! Size
          int2e%batch_indices(3,ib) = info_2e(4)  ! Atom batch belongs to
        end do remember_batches
      end if
      !
      if (scf_type=='auto') then
        mem_size = (rk_bytes*(1._rk/8._rk)*real(nbasis,kind=rk)**4)/1024._rk**3
        write (out,"('2e integralls at a single geometry require ',f10.3,' Gbytes')") mem_size
        write (out,"('                   2e storage available is ',f10.3,' Gbytes')") scf_threshold
        if (     2*mem_size<=scf_threshold) then
           scf_type = 'all in-memory'
           write (out,"('2e integrals will be kept in memory between SCF cycles')") 
        else if (  mem_size<=scf_threshold) then
           scf_type = 'in-memory'
           write (out,"('2e integrals will be kept in memory during SCF cycle')")
        else
           scf_type = 'direct'
           write (out,"('2e integrals will be recomputed for each cycle')")
        end if
      end if
      select case (scf_type)
        case default
          write (out,"('SCF type ',a,' is not recognized')") trim(scf_type)
          stop 'scf_ca%prepare_2e - bad SCF type'
        case ('direct')
          write (out,"(/'Will recompute 2-e integrals. Target accuracy of 2-e terms: ',g14.6,' Hartree'/)") accuracy_2e
        case ('in-memory','all in-memory')
          if (.not.allocated(int2e%blocks)) then
            mem_size = (rk_bytes*(1._rk/8._rk)*real(nbasis,kind=rk)**4)/1024._rk**3
            write (out,"(/'Will use about ',f10.3,' Gbytes for SCF integral storage'/)") mem_size
            allocate (int2e%blocks(int2e%nbatch_2e,int2e%nbatch_2e,int2e%nbatch_2e,int2e%nbatch_2e),stat=alloc)
            if (alloc/=0) stop 'scf_ca%prepare_2e - no memory for integral cache'
          end if
      end select
      call TimerStop('Prepare 2E Buffer')
    end subroutine prepare_2e
    !
    subroutine clear_2e
      integer(ik) :: i, j, k, l
      !
      if (scf_type=='all in-memory') return
      if (.not.associated(int2e)) return
      !
      call TimerStart('Clear 2E buffer')
      if (scf_type=='in-memory' .and. allocated(int2e%blocks)) then
        do l=1,int2e%maxbatch_2e
          do k=1,int2e%maxbatch_2e
            do j=1,int2e%maxbatch_2e
              do i=1,int2e%maxbatch_2e
                if (allocated(int2e%blocks(i,j,k,l)%ints)) deallocate (int2e%blocks(i,j,k,l)%ints)
              end do
            end do
          end do
        end do
        deallocate (int2e%blocks)
        write (out,"('Released integrals at this geometry.')")
      end if
      !
      if (allocated(int2e%batch_indices)) deallocate (int2e%batch_indices)
      call TimerStop('Clear 2E buffer')
    end subroutine clear_2e
    !
    !  This is not the most efficient implementation. The "correct" way
    !  is to pre-multiply by square root of occupation numbers, and call 
    !  ZHER2K BLAS3 routine.
    !
    subroutine density_matrix
      integer(ik) :: mu, nu, i
      !
      call TimerStart('AO Density matrix')
      rho = 0
      rho_mos: do i=1,nbasis
        if (mo_occ(i)<=0.0_rk) cycle rho_mos
        !$omp parallel do default(none) private(mu,nu) shared(mos,mo_occ,rho,i,nbasis)
        rho_nu: do nu=1,nbasis
          rho_mu: do mu=1,nbasis
            rho(mu,nu) = rho(mu,nu) + mo_occ(i)*conjg(mos(mu,i))*mos(nu,i)
          end do rho_mu
        end do rho_nu
        !$omp end parallel do
      end do rho_mos
      if (verbose>=4) then
        write (out,"(/t5,'AO density matrix')")
        call gamess_print_1e_integrals(rho,bra=gam0,ket=gam0)
      end if
      call TimerStop('AO Density matrix')
    end subroutine density_matrix
    !
    subroutine get_integrals(gam,bi,p0,sz,a2e)
      type(gam_structure), intent(inout), target :: gam
      integer(ik), intent(in)  :: bi(:)         ! Block indices
      integer(ik), intent(out) :: p0(:)         ! First orbital in the block
      integer(ik), intent(out) :: sz(:)         ! Number of orbitals
      real(rk), intent(out)    :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      !
      integer(ik) :: pe(4)
      real(rk)    :: rhomax, accuracy
      integer(ik) :: alloc
      !
      if (.not.associated(int2e)) stop 'scf_ca%get_integrals - called before initialization'
      if (any(bi<=0).or.any(bi>int2e%nbatch_2e)) stop 'scf_ca%get_integrals - bad batch indices'
      !
      p0 = int2e%batch_indices(1,bi)
      sz = int2e%batch_indices(2,bi)
      !
      !  Try to find this integral block in the cache
      !
      if (allocated(int2e%blocks)) then
        if (allocated(int2e%blocks(bi(1),bi(2),bi(3),bi(4))%ints)) then
          a2e(:sz(1),:sz(2),:sz(3),:sz(4)) = int2e%blocks(bi(1),bi(2),bi(3),bi(4))%ints
          return
        end if
      end if
      !
      if (scf_type=='direct') then
        !
        !  Figure out largest density matrix element we'll contract with
        !  This involves quite a few permutations of indices ... but at
        !  least we know rho is Hermitian ...
        !
        pe = p0 + sz - 1
        rhomax = spacing(1._rk)
        rhomax = max(rhomax,maxval(abs(rho(p0(1):pe(1),p0(2):pe(2)))))
        rhomax = max(rhomax,maxval(abs(rho(p0(1):pe(1),p0(3):pe(3)))))
        rhomax = max(rhomax,maxval(abs(rho(p0(1):pe(1),p0(4):pe(4)))))
        rhomax = max(rhomax,maxval(abs(rho(p0(2):pe(2),p0(3):pe(3)))))
        rhomax = max(rhomax,maxval(abs(rho(p0(2):pe(2),p0(4):pe(4)))))
        rhomax = max(rhomax,maxval(abs(rho(p0(3):pe(3),p0(4):pe(4)))))
        accuracy = accuracy_2e/rhomax
      else
        accuracy = -1 ! Non-direct methods need full integral accuracy
      end if
      !
      if (verbose>=5) then
        write (out,"('Retrieving integral block: ',4i8)") bi
        write (out,"('       r Initial AO index: ',4i8)") p0
        write (out,"('          r Number of AOs: ',4i8)") sz
        write (out,"('     Max absolute density: ',g15.6)") rhomax
      end if
      !
      !  Direct SCF for the moment 
      !
      call gamess_2e_integrals('AO 4C 1/R',a2e,bi,a=gam,b=gam,c=gam,d=gam,accuracy=accuracy)
      !
      !  Fill the cache
      !
      if (scf_type/='direct') then
        allocate(int2e%blocks(bi(1),bi(2),bi(3),bi(4))%ints(sz(1),sz(2),sz(3),sz(4)),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i8,' allocating storage for integral block ',4i5,' of size ',4i4)") alloc, bi, sz
          stop 'scf_ca%get_integrals - out of memory'
        end if
        int2e%blocks(bi(1),bi(2),bi(3),bi(4))%ints = a2e(:sz(1),:sz(2),:sz(3),:sz(4))
      end if
    end subroutine get_integrals
    !
    subroutine accumulate_g(bi,p0,sz,a2e,glocal)
      integer(ik), intent(in)    :: bi(:)         ! Block indices; debug only
      integer(ik), intent(in)    :: p0(:)         ! First orbital in the block
      integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
      real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      integer(ik) :: pe(4)          ! Last global index
      integer(ik) :: i, j, k, l     ! Index in the local integral buffer
      integer(ik) ::     gj, gk, gl ! Index in the global array
      complex(rk) :: s
      complex(rk) :: bcl(int2e%maxbatch_2e,int2e%maxbatch_2e)  ! Coulomb accumulation bufffer
      complex(rk) :: bex(int2e%maxbatch_2e,int2e%maxbatch_2e)  ! Exchange accumulation buffer
      !
      if (verbose>=5) then
        write (out,"('    Adding integral block: ',4i8)") bi
        write (out,"('       a Initial AO index: ',4i8)") p0
        write (out,"('          a Number of AOs: ',4i8)") sz
        write (out,"(10(1x,f16.10))") a2e(:sz(1),:sz(2),:sz(3),:sz(4))
      end if
      !
      !  Coulomb contribution first
      !
      coul_j: do j=1,sz(2)
        coul_i: do i=1,sz(1)
          s = 0
          coul_l: do l=1,sz(4)
            gl = p0(4)+l-1
            coul_k: do k=1,sz(3)
              gk = p0(3)+k-1
              s = s + rho(gk,gl) * a2e(i,j,k,l)
            end do coul_k
          end do coul_l
          bcl(i,j) = s
        end do coul_i
      end do coul_j
      !
      !  Now the exchange term
      !
      xchg_l: do l=1,sz(4)
        xchg_i: do i=1,sz(1)
          s = 0
          xchg_j: do j=1,sz(2)
            gj = p0(2)+j-1
            xchg_k: do k=1,sz(3)
              gk = p0(3)+k-1
              s = s + rho(gk,gj) * a2e(i,j,k,l)
            end do xchg_k
          end do xchg_j
          bex(i,l) = -0.5_rk*s
        end do xchg_i
      end do xchg_l
      !
      !  Update local copy of G matrix
      !
      pe = p0 + sz - 1
      glocal(p0(1):pe(1),p0(2):pe(2)) = glocal(p0(1):pe(1),p0(2):pe(2)) + bcl(:sz(1),:sz(2))
      glocal(p0(1):pe(1),p0(4):pe(4)) = glocal(p0(1):pe(1),p0(4):pe(4)) + bex(:sz(1),:sz(4))
    end subroutine accumulate_g
    !
    !  A very naive routine for computing 2-electron part of the closed-shell Fock matrix
    !
    subroutine g_matrix_block(gam,ci,cj,ck,cl,glocal)
      type(gam_structure), intent(inout), target :: gam
      !
      integer(ik), intent(in)    :: ci, cj, ck, cl ! Integral block to work on
      complex(rk), intent(inout) :: glocal(:,:)    ! Pre-thread copy of gmat; safe to update 
                                                   ! without locking.
      integer(ik) :: bi(4), b2(4)                  ! Batch indices; copy of ci,cj,ck,cl
      integer(ik) :: p0(4), p2(4)                  ! Initial position of the integral block
      integer(ik) :: sz(4), s2(4)                  ! Size of the integral block
      real(rk)    :: a2e(int2e%maxbatch_2e,int2e%maxbatch_2e,int2e%maxbatch_2e,int2e%maxbatch_2e)
      real(rk)    :: a22(int2e%maxbatch_2e,int2e%maxbatch_2e,int2e%maxbatch_2e,int2e%maxbatch_2e)
      logical     :: go                            ! Set to true if swap is allowed
      !
      bi = (/ ci, cj, ck, cl /)
      call get_integrals(gam,bi,p0,sz,a2e)
      !
      !  This integral block may appear several times in G matrix construction
      !
                               call accumulate_g(bi,p0,sz,a2e,glocal)
      call swap_jikl ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_ijlk ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_jilk ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_lkij ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_klij ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_klji ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      call swap_lkji ; if (go) call accumulate_g(b2,p2,s2,a22,glocal) 
      return
      !
      contains
      !
      !  The routines below apply index permutations to integral blocks, so
      !  that we do not need to compute any of the blocks more than once.
      !  The only tricky bit are the conditional expressions used to cull
      !  redundant permutations.
      !
      subroutine swap_jikl
        integer(ik), parameter :: code(4) = (/2,1,3,4/)
        integer(ik)            :: i, j, k, l
        !
        go = (ci/=cj)
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(j,i,k,l) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_jikl
      !
      subroutine swap_ijlk
        integer(ik), parameter :: code(4) = (/1,2,4,3/)
        integer(ik)            :: i, j, k, l
        !
        go = (ck/=cl)
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(i,j,l,k) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_ijlk
      !
      subroutine swap_jilk
        integer(ik), parameter :: code(4) = (/2,1,4,3/)
        integer(ik)            :: i, j, k, l
        !
        go = (ci/=cj) .and. (ck/=cl)
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(j,i,l,k) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_jilk
      !
      subroutine swap_klij
        integer(ik), parameter :: code(4) = (/3,4,1,2/)
        integer(ik)            :: i, j, k, l
        !
        go = ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(k,l,i,j) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_klij
      !
      subroutine swap_klji
        integer(ik), parameter :: code(4) = (/3,4,2,1/)
        integer(ik)            :: i, j, k, l
        !
        go = (ci/=cj) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(k,l,j,i) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_klji
      !
      subroutine swap_lkij
        integer(ik), parameter :: code(4) = (/4,3,1,2/)
        integer(ik)            :: i, j, k, l
        !
        go = (ck/=cl) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(l,k,i,j) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_lkij
      !
      subroutine swap_lkji
        integer(ik), parameter :: code(4) = (/4,3,2,1/)
        integer(ik)            :: i, j, k, l
        !
        go = (ci/=cj) .and. (ck/=cl) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
        if (.not.go) return
        !
        b2 = bi(code) ; p2 = p0(code) ; s2 = sz(code)
        do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
          a22(l,k,j,i) = a2e(i,j,k,l)
        end do ; end do ; end do ; end do
      end subroutine swap_lkji
    end subroutine g_matrix_block
    !
    subroutine g_matrix(gam)
      type(gam_structure), intent(inout), target :: gam
      integer(ik) :: clk, nclk        ! Combined index; needed to get a decent number of blocks
                                      ! for OpenMP to work on
      integer(ik) :: ci, cj, ck, cl   ! Unfortunately, Fortran does not allow array elements
                                      ! as iterators. A shame, really
      integer(ik) :: cimax            ! We have to be a little careful with this index!
      integer(ik) :: alloc
      integer(ik), allocatable :: clktab(:,:)
      complex(rk), allocatable :: glocal(:,:)
      !
      call TimerStart('2e Fock contribution')
      !
      !  Prepare list of (l,k) index pairs; the explicit list is not actually
      !  necessary (it can be replaced by floating-point math on the combined index),
      !  but the overhead of the list is negligible, and the code is easier to
      !  understand.
      !
      nclk = (int2e%nbatch_2e*(int2e%nbatch_2e+1))/2
      allocate (clktab(2,nclk),stat=alloc)
      if (alloc/=0) stop 'scf_ca%g_matrix - Can''t allocate pair table'
      clk = 0
      batch_l: do cl=1,int2e%nbatch_2e
        batch_k: do ck=1,cl
          clk = clk + 1
          clktab(:,clk) = (/cl,ck/)
        end do batch_k
      end do batch_l
      if (clk/=nclk) stop 'scf_ca%g_matrix - Count error'
      !
      !  Construct the G matrix
      !
      gmat = 0
      !$omp parallel default(none) private(clk,ci,cj,ck,cl,cimax,glocal,alloc) &
      !$omp& shared(clktab,nclk,nbasis,gmat,gam)
      allocate (glocal(nbasis,nbasis),stat=alloc)
      if (alloc/=0) stop 'scf_ca%g_matrix - no memory for OpenMP buffers'
      glocal = 0
      !$omp do schedule(dynamic,1)
      batch_lk: do clk=1,nclk
        ! write (out,"('Block ',i0,' of ',i0)") clk, nclk
        cl = clktab(1,clk)
        ck = clktab(2,clk)
        batch_j: do cj=1,cl
          cimax = cj
          if (cj==cl) cimax = ck
          batch_i: do ci=1,cimax
            call g_matrix_block(gam,ci,cj,ck,cl,glocal)
          end do batch_i
        end do batch_j
      end do batch_lk
      !$omp end do
      !
      !  Accumulate the total 2-electron term
      !
      !$omp critical
      gmat = gmat + glocal
      !$omp end critical
      deallocate (glocal)
      !$omp end parallel
      deallocate (clktab)
      !
      if (verbose>=4) then
        write (out,"(/t5,'2e contribution to the Fock matrix')")
        call gamess_print_1e_integrals(gmat,bra=gam,ket=gam)
      end if
      call TimerStop('2e Fock contribution')
    end subroutine g_matrix
    !
    !  Solve generalized eigenvalue problem
    !
    subroutine diagonalize_fmat
      call TimerStart('Eigenvalue problem')
      !
      !  Construct S^{-1/2} F S^{-1/2}
      !
      mos = matmul(shalf,matmul(fmat,shalf))
      if (verbose>=3) then
        write (out,"(/t5,'S^{-1/2} F S^{-1/2}')")
        call gamess_print_1e_integrals(mos,bra=gam0,ket=gam0)
      end if
      call lapack_heev(mos,mo_energy)
      if (verbose>=0) then
        write (out,"(/t5,'MO energies:')")
        write (out,"(10(1x,f12.6))") mo_energy
        write (out,"()")
      end if
      !
      !  Construct eigenvectors
      !
      mos = matmul(shalf,mos)
      call TimerStop('Eigenvalue problem')
    end subroutine diagonalize_fmat
    !
    !  Hartree-Fock total electronic energy. Since we are dealing with complex
    !  operators and density matrices, the expression is Tr(P*H^\dagger)
    !
    subroutine total_energy
      complex(rk) :: efock, eg, elen

      efock     =      sum(rho * fmat)
      eg        =      sum(rho * gmat)
      energy_ca = real(sum(rho * hca ),kind=rk)
      elen      = efock - 0.5_rk * eg
      ! write (out,*) ' efock = ',efock
      ! write (out,*) ' eg    = ',eg
      if (aimag(elen)>100._rk*spacing(real(elen,kind=rk))) then
        write (out,"('WARNING: Imaginary part of the electronic energy is too large: ',2g24.14)") elen
      end if
      etot = real(elen,kind=rk) + enuc
      if (verbose>=0) then
        write (out,"(/'Total energy = ',f20.12/)") etot
      end if
    end subroutine total_energy
    !
    !  Transform initial Fock matrix into the MO basis, giving an estimate
    !  of the eigenvalues
    !
    subroutine initial_eigenvalues
      integer(ik) :: imo

      fock_diagonal: do imo=1,nbasis
        mo_energy(imo) = real(dot_product(mos(:,imo),matmul(fmat,mos(:,imo))),kind=rk)
      end do fock_diagonal
      !
      write (out,"(/t5,'Initial Fock matrix diagonal:')")
      write (out,"(10(1x,f12.6))") mo_energy
      write (out,"()")
    end subroutine initial_eigenvalues
    !
    subroutine scf_loop(gam)
      type(gam_structure), intent(inout), target :: gam
      !
      integer(ik) :: iter
      real(rk)    :: etot_old
      real(rk)    :: de, drho
      logical     :: converged
      !
      call TimerStart('SCF loop')
      converged = .false.
      repeat_scf: do iter=1,max_scf_iter
        rho_old  = rho
        etot_old = etot
        !
        write (out,"('Beginning SCF cycle ',i4)") iter
        call density_matrix
        call g_matrix(gam)
        fmat = hbo + gmat + hca
        if (verbose>=0 .and. iter==1) then
          call initial_eigenvalues
        end if
        if (verbose>=3) then
          write (out,"(/t5,'Fock matrix')")
          call gamess_print_1e_integrals(fmat,bra=gam,ket=gam)
        end if
        call total_energy
        call diagonalize_fmat
        !
        !  Check for convergence; should only check after one full cycle!
        !
        if (iter<=1) cycle repeat_scf
        de   = etot - etot_old
        drho = maxval(abs(rho-rho_old))
        if (verbose>=0) then
          write (out,"('Iteration ',i4,' etot= ',f20.12,' de= ',g12.5,' drho= ',g12.5)") iter, etot, de, drho
        end if
        if (abs(de)<=scf_eps_energy .and. drho<=scf_eps_rho) then
          converged = .true.
          write (out,"('SCF converged')")
          exit repeat_scf
        end if
      end do repeat_scf
      !
      if (.not.converged) stop 'scf_ca%scf_loop - SCF convergence failure'
      !
      if (verbose>=2) then
        write (out,"(/t5,'Final MOs')")
        call gamess_print_1e_integrals(mos,bra=gam,ket=gam,heading='LEFT',symmetry='NONE')
      end if
      call TimerStop('SCF loop')
    end subroutine scf_loop
    !
    !  Perform SCF at the desired geometry
    !
    subroutine perform_scf(gam,ecp,ints)
      type(gam_structure), intent(inout), target :: gam
      type(ecp_molecule), intent(inout), target  :: ecp
      type(int2e_cache), intent(inout), target   :: ints
      !
      call one_electron_hamiltonian(gam,ecp)
      !
      call prepare_2e(gam,ints)
      !
      call scf_loop(gam)
      !
      call clear_2e
      !
      call TimerReport
    end subroutine perform_scf
!   !
!   !  Classify orbitals into degenerate subsets - too simplistic
!   !
!   subroutine classify_degeneracy(energy,cnt,ind)
!     real(rk), intent(in)     :: energy(:)  ! MO energies
!     integer(ik), intent(out) :: cnt        ! Number of degenerate blocks
!     integer(ik), intent(out) :: ind(:)     ! First orbital of each degenerate block;
!                                            ! one extra entry at ind(cnt+1)
!     !
!     integer(ik) :: nmo, imo, ibl, sbl
!     !
!     nmo = size(energy)
!     cnt = 1 ; ind(1) = 1
!     scan_mos: do imo=2,nmo
!       if (abs(energy(imo)-energy(imo-1))<=eps_degenerate) cycle scan_mos
!       cnt = cnt + 1
!       if (cnt+1>size(ind)) stop 'scf_ca%classify_degeneracy - blown ind buffer'
!       ind(cnt) = imo
!     end do scan_mos
!     ind(cnt+1) = nmo + 1
!     if (verbose<1) return
!     write (out,"('Found ',i0,' degenerate blocks within ',i0,' MOs')") cnt, nmo
!     write (out,"(1x,a5,1x,a5,1x,a4,2x,a12,1x,a12)") &
!            'BLOCK', ' 1ST ', ' CNT ', ' Emin ', ' Emax ', &
!            '-----', '-----', '-----', '------', '------'
!     print_blocks: do ibl=1,cnt
!       sbl = ind(ibl+1)-ind(ibl)
!       if (sbl>1) then
!         write (out,"(1x,i5,1x,i5,1x,i4,2x,f12.6,1x,f12.6)") &
!                ibl, ind(ibl), sbl, energy(ind(ibl)), energy(ind(ibl+1)-1)
!       else
!         write (out,"(1x,i5,1x,i5,1x,i4,2x,f12.6,1x,f12.6)") &
!                ibl, ind(ibl), sbl, energy(ind(ibl))
!       end if
!     end do print_blocks
!     write (out,"()")
!   end subroutine classify_degeneracy
    !
    !  Align a diagonal MO block to maximize the diagonal part of the overlap
    !  We will also update the MOs in mos0() and mosd() and overlap matrix in ca_mo_smat
    !
    subroutine align_mo_block(rb0,rbd)
      integer(ik), intent(in) :: rb0(:)    ! List of the left-hand MOs in the block (mos0)
      integer(ik), intent(in) :: rbd(:)    ! List of the right-hand MOs in the block (mosd)
      !
      complex(rk) :: sblk(size(rb0),size(rbd))       ! Diagonal sub-block of the overlap matrix
      complex(rk) :: su  (size(rb0),size(rb0))       ! Left singular vectors
      complex(rk) :: svth(size(rbd),size(rbd))       ! Right singular vectors
      real(rk)    :: sing(min(size(rb0),size(rbd))) ! Singular values
      !
      sblk = ca_mo_smat(rb0,rbd)
      ! write (out,"('  sblk = ',16f14.6)") sblk
      call lapack_svd(sblk,sing,su,svth)
      ! write (out,"('  sing = ', 4f14.6)") sing
      ! write (out,"('  su   = ',16f14.6)") su
      ! write (out,"('  svth = ',16f14.6)") svth
      !
      !  Transform the MOs and the overlap matrix
      !
      mosd(:,rbd) = matmul(mosd(:,rbd),conjg(transpose(svth)))
      mos0(:,rb0) = matmul(mos0(:,rb0),su)
      !
      ca_mo_smat(:,rbd) = matmul(ca_mo_smat(:,rbd),conjg(transpose(svth)))
      ca_mo_smat(rb0,:) = matmul(conjg(transpose(su)),ca_mo_smat(rb0,:))
    end subroutine align_mo_block
    !
    !  Bring target and displaced MOs into the maximal alignment, by rotating within the
    !  degenerate subsets, and swapping orbitals as needed.
    !
    subroutine align_mos
      integer(ik) :: mo0, mod, mx0
      integer(ik) :: id0, idd
      integer(ik) :: rb0(nbasis), nrb0 ! List of orbitals in the degenerate block on the left
      integer(ik) :: rbd(nbasis), nrbd ! List of orbitals in the degenerate block on the right
      logical     :: grown
      complex(rk) :: tmp(nbasis)
      !
      !  While rotating and swapping, we must make sure that the input MOs are 
      !  sensible. We'll require that matrix elements large enough to be 
      !  considered for a stabilizing rotation only occur between orbitals
      !  which could be considered degenerate. We will also absolutely forbid
      !  occupied-virtual rotations.
      !
      !  Since each rotation block can grow from the "seed" orbital, we'll probably
      !  sweep some blocks more than once. That's OK, since our alignment procedure
      !  is idempotent.
      !
      scan_seed_left: do mo0=1,nbasis
        rb0(1) = mo0 ; nrb0 = 1 ; nrbd = 0 ;
        grown= .true.
        !
        !  We may need several passes before we get a fully consistent block
        !
        grow_block: do while(grown)
          grown = .false.
          call add_all_degenerate(nrb0,rb0,mo_energy0,grown)
          scan_overlaps_left: do id0=1,nrb0
            scan_right: do mod=1,nbasis
              if (abs(ca_mo_smat(rb0(id0),mod))<mix_degenerate) cycle scan_right
              if (any(mod==rbd(1:nrbd))) cycle scan_right
              nrbd = nrbd + 1
              rbd(nrbd) = mod
              grown = .true.
            end do scan_right
          end do scan_overlaps_left
          call add_all_degenerate(nrbd,rbd,mo_energyd,grown)
          scan_overlaps_right: do idd=1,nrbd
            scan_left: do mx0=1,nbasis
              if (abs(ca_mo_smat(mx0,rbd(idd)))<mix_degenerate) cycle scan_left
              if (any(mx0==rb0(1:nrb0))) cycle scan_left
              nrb0 = nrb0 + 1
              rb0(nrb0) = mx0
              grown = .true.
            end do scan_left
          end do scan_overlaps_right
        end do grow_block
        !
        call sort(rb0(1:nrb0))
        call sort(rbd(1:nrbd))
        call check_for_missing_couplings
        !
        !  We now have a consistent block of degenerate orbitals
        !
        ! write (out,"(' Left = ',10i12 )") rb0(1:nrb0)
        ! write (out,"('        ',10f12.6)") mo_energy0(rb0(1:nrb0))
        ! write (out,"('Right = ',10i12 )") rbd(1:nrbd)
        ! write (out,"('        ',10f12.6)") mo_energyd(rbd(1:nrbd))
        !
        call align_mo_block(rb0(1:nrb0),rbd(1:nrbd))
      end do scan_seed_left
      !
      !  We now have a maximally-diagonal overlap matrix; however, the largest
      !  overlaps are not guaranteed to be on the diagonal of the overlap matrix,
      !  so that we may need to swap a few orbitals. 
      !
      swap_sweep: do mo0=1,nbasis
        mod = maxloc(abs(ca_mo_smat(mo0,:)),dim=1)
        if (mod==mo0) cycle swap_sweep
          write (out,"('Need to swap displaced MOs ',i6,' and ',i6)") mo0, mod
        tmp = mosd(:,mo0) ; mosd(:,mo0) = mosd(:,mod) ; mosd(:,mod) = tmp
        tmp = ca_mo_smat(:,mo0) ; ca_mo_smat(:,mo0) = ca_mo_smat(:,mod) ; ca_mo_smat(:,mod) = tmp
      end do swap_sweep
      return
      !
      contains
      !
      subroutine add_all_degenerate(cnt,lst,energy,grown)
        integer(ik), intent(inout) :: cnt
        integer(ik), intent(inout) :: lst(:)
        real(rk), intent(in)       :: energy(:)
        logical, intent(inout)     :: grown
        !
        real(rk)    :: en_min, en_max ! Min and max eigenvalue already included
        integer(ik) :: imo
        !
        en_min = minval(energy(lst(1:cnt)))
        en_max = maxval(energy(lst(1:cnt)))
        if (en_max>en_min+epsd_tolerate) then
          write (out,"('Trying to grow degenerate list containing non-degenerate MOs')")
          write (out,"(5(1x,i16))") lst(1:cnt) 
          write (out,"(5(1x,f16.10))") energy(lst(1:cnt))
          write (out,"('Try decreasing the value of CA_TIME')")
          stop 'scf_ca%align_mos%add_all_degenerate - non-degenerate MOs in degenerate list'
        end if
        !
        grow_list: do imo=1,size(energy)
          if (all(abs(energy(imo)-energy(lst(1:cnt)))>epsd_include)) cycle grow_list
          if (any(imo==lst(1:cnt))) cycle grow_list
          cnt = cnt + 1
          if (cnt>size(lst)) stop 'scf_ca%align_mos%add_all_degenerate - blown buffer'
          lst(cnt) = imo
          grown  = .true.
        end do grow_list
      end subroutine add_all_degenerate
      !
      subroutine check_for_missing_couplings
        integer(ik) :: idd
        integer(ik) :: mo0, mod
        !
        scan_right: do idd=1,nrbd
          mod = rbd(idd)
          scan_left: do mo0=1,nbasis
            if (abs(ca_mo_smat(mo0,mod))<mix_degenerate) cycle scan_left
            if (any(mo0==rb0(1:nrb0))) cycle scan_left
            write (out,"('Found an inconsistent degenerate block mapping')")
            write (out,"('Target orbitals:')")
            write (out,"(5(1x,i16))") rb0(1:nrb0)
            write (out,"(5(1x,f16.10))") mo_energy0(rb0(1:nrb0))
            write (out,"('Displaced orbitals:')")
            write (out,"(5(1x,i16))") rbd(1:nrbd)
            write (out,"(5(1x,f16.10))") mo_energyd(rbd(1:nrbd))
            write (out,"('Block overlap:')")
            write (out,"(5(1x,g13.6,1x,g13.6,1x))") ca_mo_smat(rb0(1:nrb0),rbd(1:nrbd))
            write (out,"('External coupling:')")
            write (out,"(1x,i6,1x,f16.10,' to ',1x,i6,1x,f16.10,' = ',g13.6,1x,g13.6)") &
                   mo0, mo_energy0(mo0), mod, mo_energy(mod), ca_mo_smat(mo0,mod)
            stop 'scf_ca%align_mos%add_all_degenerate - inconsistent degenerate block'
          end do scan_left
        end do scan_right
      end subroutine check_for_missing_couplings
    end subroutine align_mos
    !
    !  The CA term in the Hamiltonian is calculated using numerical differentiation.
    !  The procedure is somewhat elaborate, but in principle simple:
    !    1. We maximally align the MOs at the two geometries, by rotating 
    !       them within the degenerate subsets. This step gets a little
    !       involved because the degeneracies do not need to have the same
    !       dimentionality at the two points, and some MO reordering is
    !       likely to occur
    !    2. We then recalculate the overlaps using the rotated MOs
    !    3. We calculate the CA contribution, by zeroing out the diagonal 
    !       and multiplying with (-I/ca_time)
    !    3. We back-transform the CA contribution to the AO basis.
    !
    subroutine calculate_complete_adiabatic_term
      integer(ik) :: mo0
      complex(rk) :: s0
      real(rk)    :: hdev, hmax
      !
      call TimerStart('CA Hamiltonian')
      ca_mo_smat = matmul(conjg(transpose(mos0)),matmul(ca_ao_smat,mosd))
      if (verbose>=2) then
          write (out,"(/t5,'DISPLACED MO OVERLAP MATRIX BEFORE ADJUSTMENT'/)")
        call gamess_print_1e_integrals(ca_mo_smat,bra=gam0,ket=gamd,symmetry='NONE',heading='NONE')
      end if
      !
      call align_mos
      !
      ca_mo_smat = matmul(conjg(transpose(mos0)),matmul(ca_ao_smat,mosd))
      if (verbose>=2) then
          write (out,"(/t5,'DISPLACED MO OVERLAP MATRIX AFTER ADJUSTMENT'/)")
        call gamess_print_1e_integrals(ca_mo_smat,bra=gam0,ket=gamd,symmetry='NONE',heading='NONE')
      end if
      !
      !  The MOs are in the maximal coincidence; go through them to make sure
      !  the correlations are of acceptable quality. Once we are satisfied,
      !  the diagonal will have to go.
      !
      displaced_mos: do mo0=1,nbasis
        s0 = ca_mo_smat(mo0,mo0)
        if (abs(s0)<=0.90_rk) then
          write (out,"('WARNING: Displaced MO ',i0,' overlap is small (',2g14.6,'). Try decreasing ca_time')") &
                 mo0, s0
          if (abs(s0)<=0.75_rk) stop 'scf_ca%calculate_complete_adiabatic_term - lost MO correlation'
        end if
        ca_mo_smat(mo0,mo0) = 0 
      end do displaced_mos
      !
      !
      !  Finite-difference Complete Adiabatic term, MO basis
      !
      ca_mo_smat = (-(0,1)/ca_time) * ca_mo_smat
      !
      !  Transform to the AO basis:
      !
      !    HAO = S C HMO C^\dagger S
      ! 
      hca_old = hca
      hca = matmul(matmul(smat,mos0),matmul(ca_mo_smat,matmul(transpose(conjg(mos0)),smat)))
      hca = hca_mixing * hca + (1.0_rk-hca_mixing) * hca_old
      !
      if (verbose>=2) then
        write (out,"(/t5,'COMPLETE ADIABATIC CONTRIBUTION BEFORE HERMITIAN PROJECTION'/)")
        call gamess_print_1e_integrals(hca,bra=gam0,ket=gam0)
      end if
      !
      !  We require the Hermitian part of the matrix element
      !
      hdev = 0.5_rk * maxval(abs(hca-transpose(conjg(hca))))
      hca  = 0.5_rk * (hca + transpose(conjg(hca)))
      hmax = maxval(abs(hca))
      if (verbose>=0) then
        write (out,"('     Largest matrix element magnitude of the Hermitian part of H(CA): ',g14.6)") hmax
        write (out,"('Largest matrix element magnitude of the anti-Hermitian part of H(CA): ',g14.6/)") hdev
      endif
      !
      if (verbose>=2) then
        write (out,"(/t5,'COMPLETE ADIABATIC CONTRIBUTION TO 1-ELECTRON HAMILTONIAN'/)")
        call gamess_print_1e_integrals(hca,bra=gam0,ket=gam0)
      end if
      !
      call TimerStop('CA Hamiltonian')
    end subroutine calculate_complete_adiabatic_term
    !
    subroutine punch_final_mos
      integer(ik) :: ios, iat, imo
      real(rk)    :: xyz(3), q
      !
      if (mos_file==' ') return
      !
      call TimerStart('Punch final MOs')
      open(iu_mos,file=trim(mos_file),form='formatted',action='write',status='replace',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' opening file ',a,' for the final MOs.')") ios, trim(mos_file)
        stop 'scf_ca%punch_final_mos - Can''t create file'
      end if
      write (iu_mos,"('""Complete Adiabatic"" MOs for ',a)") trim(molecule_file)
      write (iu_mos,"('Nuclear coordinates and velocities (atomic units): ')")
      write (iu_mos,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  VX  ', '  VY  ', '  VZ  '
      print_atoms: do iat=1,natoms
        xyz = real(gam0%atoms(iat)%xyz,kind=kind(xyz))
        q   = real(gam0%atoms(iat)%znuc,kind=kind(xyz))
        write (iu_mos,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") q, xyz/abohr, vatom(:,iat)
      end do print_atoms
      write (iu_mos,"('Final energy = ',f25.15)") etot
      write (iu_mos,"('Final orbital eigenvalues:')")
      write (iu_mos,"(5(1x,f12.6))") mo_energy
      write (iu_mos,"('For each MO, real part precedes the imaginary part')")
      write (iu_mos,"(' $VECCA ')")
      punch_mos: do imo=1,nbasis
        call punch_vector(iu_mos,2*imo-1, real(mos(:,imo),kind=rk))
        call punch_vector(iu_mos,2*imo-0,aimag(mos(:,imo)))
      end do punch_mos
      write (iu_mos,"(' $END   ')")
      close (iu_mos)
      call TimerStop('Punch final MOs')
    end subroutine punch_final_mos
    !
    !  Allocate data structures for the property evaluation grid
    !
    subroutine prepare_currents_grid
      integer(ik) :: alloc
      integer(ik) :: ia, ib, ic, ipt
      !
      call TimerStart('Prepare currents grid')
      if (any(grid_counts<=0)) then
        stop 'scf_ca%prepare_currents_grid - prooperty grid vanishes.'
      end if
      grid_points = product(grid_counts)
      allocate (grid_xyz(3,grid_points),grid_property(0:4,grid_points),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i0,' allocating property grid with ',i0,' points.')") alloc, grid_points
        stop 'scf_ca%prepare_currents_grid - no memory for properties grid.'
      end if
      !
      if (verbose>=0) then
        write (out,"(/'Evaluating electron density and currents on a grid:')")
        write (out,"( '    Origin: ',3(1x,f15.8))") grid_origin
        write (out,"( '         a: ',3(1x,f15.8),' x ',i6,' pts')") grid_step(:,1), grid_counts(1)
        write (out,"( '         b: ',3(1x,f15.8),' x ',i6,' pts')") grid_step(:,2), grid_counts(2)
        write (out,"( '         c: ',3(1x,f15.8),' x ',i6,' pts')") grid_step(:,3), grid_counts(3)
        write (out,"( '  Total number of points: ',i0/)") grid_points
      end if
      !
      ipt = 0
      fill_c: do ic=1,grid_counts(3)
        fill_b: do ib=1,grid_counts(2)
          fill_a: do ia=1,grid_counts(1)
            ipt = ipt + 1
            grid_xyz(:,ipt) = grid_origin + matmul(grid_step,real((/ia,ib,ic/)-1,kind=rk))
          end do fill_a
        end do fill_b
      end do fill_c
      if (ipt/=grid_points) stop 'scf_ca%prepare_currents_grid - count error'
      call TimerStop('Prepare currents grid')
    end subroutine prepare_currents_grid
    !
    !  Evaluate electron density and currents on a single grid points
    !
    subroutine evaluate_currents(xyz,vr)
      real(rk), intent(in)  :: xyz(3)         ! Coordinates of the point
      real(rk), intent(out) :: vr(0:3)        ! Density (0) and current (1:3) at this point
      !
      real(ark)    :: basval(0:3,gam0%nbasis)  ! Buffer for basis function values
      complex(ark) :: v(0:3)                  ! Temporary accumulation buffer, to make sure imaginary part vanishes
      integer(ik)  :: ja, jb
      real(rk)     :: eps
      !
      call gamess_evaluate_functions(xyz,basval,structure=gam0)
      !
      !  If this routine enters the critical path, we can cut the cost by 1/2 
      !  by using hermitian property of the density matrix.
      !
      v = 0
      basis_right: do jb=1,nbasis
        basis_left: do ja=1,nbasis
          v(0) = v(0) + basval(0,ja)*basval(0,jb) * rho(ja,jb)
          v(1:3) = v(1:3) + rho(ja,jb) * ( basval(1:3,ja)*basval(0,jb) - basval(0,ja)*basval(1:3,jb) )
        end do basis_left
      end do basis_right
      v(1:3) = 0.5_rk * (0,1) * v(1:3)
      !
      !  Make sure the imaginary part of the sum disappeared
      !
      eps = 1e3_rk * spacing(maxval(abs(real(v,kind=rk))))
      if (any(abs(aimag(v))>eps)) then
        write (out,"('Problem evaluating currents at point:'/3(1x,g25.16))") xyz
        write (out,"('Imaginary parts of the observables failed to vanish:'/((1x,g25.16,1x,g25.16)))") v
        stop 'scf_ca%evaluate_currents - unreal'
      end if
      !
      vr = real(v,kind=rk)
    end subroutine evaluate_currents
    !
    !  Evaluate electronic currents on a grid
    !
    subroutine evaluate_currents_on_grid
      integer(ik) :: ipt, ic
      real(rk)    :: xyz(3), up(0:3), down(0:3)
      !
      call TimerStart('Evaluate currents')
      !$omp parallel do default(none) private(ipt,ic,xyz,up,down) &
      !$omp&            shared(grid_points,grid_xyz,grid_property,grid_dx)
      point_loop: do ipt=1,grid_points
        !
        !  First, density and current at the central position
        !
        xyz = grid_xyz(:,ipt)
        call evaluate_currents(xyz,grid_property(0:3,ipt))
        !
        !  Now calculate -div j using numerical differentiation.
        !
        grid_property(4,ipt) = 0
        div_steps: do ic=1,3
          xyz(ic) = xyz(ic) + grid_dx
          call evaluate_currents(xyz,up)
          xyz(ic) = xyz(ic) - 2*grid_dx
          call evaluate_currents(xyz,down)
          xyz(ic) = grid_xyz(ic,ipt)
          !
          grid_property(4,ipt) = grid_property(4,ipt) - (up(ic)-down(ic))/(2*grid_dx)
        end do div_steps
      end do point_loop
      !$omp end parallel do
      !
      if (odx_file/=' ') then
        call odx_write_cube(trim(odx_file),grid_property,grid_origin,grid_step,grid_counts)
        write (out,"('Wrote density and currents data to ',a)") trim(odx_file)
      end if
      !
      !  All done; generate text output if needed.
      !
      if (verbose>=2 .or. print_grid) then
        write (out,"(/t5,'Electronic densities and currents [atomic units]'/)")
        write (out,"((3(1x,a12),1x,a13,1x,3(1x,a13),2x,a13))") &
               '  X    ', '  Y    ', '  Z    ', '  RHO    ', '  JX    ', '  JY    ', '  JZ    ', '  -DIV J    ', &
               '-----  ', '-----  ', '-----  ', '-------  ', '------  ', '------  ', '------  ', '----------  '
        print_loop: do ipt=1,grid_points
          write (out,"('@',3(1x,f12.6),1x,g13.6,1x,3(1x,g13.6),2x,g13.6)") grid_xyz(:,ipt), grid_property(:,ipt)
        end do print_loop
        write (out,"()")
      end if
      call TimerStop('Evaluate currents')
    end subroutine evaluate_currents_on_grid
    !
    subroutine start
      integer(ik) :: info, ca_iter
      real(rk)    :: etot_old, de, drho
      !
      call TimerStart('start')
      !
      read (input,nml=scf_ca_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=scf_ca_data)
      write (out,"()")
      !
      call load_gamess_data
      !
      call allocate_dynamic_data
      !
      call load_velocities
      !
      call prepare_displaced_geometry
      !
      call prepare_ecps
      !
      call fill_occupations
      !
      call prepare_currents_grid
      !
      !  Calculation of the Complete Adiabatic term requires double SCF, involving
      !  two geometries along the displacement vector. Start with the CA part of
      !  the Hamiltonian set to zero, and work from there.
      !
      hca = 0
      if (use_ca) call complete_adiabatic_hamiltonian_guess
      mos0 = gam0%vectors(:,:nbasis)
      call TimerReport
      !
      etot_old = 0
      ca_loop: do ca_iter=1,max_ca_iter
        write (out,"(/'Starting CA iteration ',i4,' (target geometry)')") ca_iter
        !
        !  SCF at the target geometry
        !
        mos  = mos0
        call perform_scf(gam0,ecp0,int2e0)
        mos0 = mos ; mo_energy0 = mo_energy
        if (.not.use_ca) exit ca_loop
        if (ca_iter>=2) then
          de   = etot - etot_old
          drho = maxval(abs(rho-rho0))
          write (out,"('CA iteration ',i4,' etot= ',f20.12,' energy(ca)= ',f20.12,' de= ',g12.5,' drho= ',g12.5)") &
                 ca_iter, etot, energy_ca, de, drho
          if (abs(de)<=ca_eps_energy .and. drho<=ca_eps_rho) then
            write (out,"('Double-SCF converged')")
            exit ca_loop
          end if
          if (ca_iter>=max_ca_iter) stop 'scf_ca%start - Double-SCF failed'
        end if
        rho0     = rho
        etot_old = etot
        !
        !  SCF at the displaced geometry
        !
        write (out,"(/'         CA iteration ',i4,' (displaced geometry)')") ca_iter
        if (ca_iter>=2) mos = mosd
        call perform_scf(gamd,ecpd,int2ed)
        mosd = mos ; mo_energyd = mo_energy
        call calculate_complete_adiabatic_term
      end do ca_loop
      !
      call punch_final_mos
      !
      call evaluate_currents_on_grid
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
  end module scf_ca
!
  program driver
    use scf_ca
    use accuracy
    use math

    real(rk) :: hole

    call accuracyInitialize
    !
    !  We have to issue calls below in case we run under OpenMP.
    !  See comments in math.f90
    !
    hole = MathFactorial(80)
    hole = MathDoubleFactorial(80)
    hole = MathLogFactorial(80)
    !
    !$ write (out,"(/'WARNING: Some compilers are known to miscompile OpenMP constructs in this program')")
    !$ write (out,"( 'WARNING: Please make sure the compiler you are using yields correct results!'/)")
    !
    call start
  end program driver
