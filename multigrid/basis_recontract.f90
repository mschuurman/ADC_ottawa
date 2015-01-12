!
!  Implementation of atomic SCF intended for basis-set recontraction.
!  This is essentially a cut-down, specialized version of static_tunnel.f90
!  In particular, most of the error- and sanity-checking has been removed
!  - so the code will happily accept totally senseless input parameters.
!
!  We are trying to reduce the severity of linear dependence problems in 
!  molecular SCF, by recontracting the complete basis (ie standard atomic
!  + scattering Kaufmann, for example). This is basically the same idea
!  as ANOs, but with a somewhat different aim in mind.
!
!  We proceed as follows:
!  1. Solve atomic SCF problem, with a suitable choice of (fractional) occupations
!  2. For each angular momentum, identify the individual orbital blocks
!  3. Produce the required recontractions
!
!  In preparing the recontractions, we tread all previosly-contracted functions
!  (e.g. from a contracted atomic basis set) as a block. 
!
!  Everything happens in quad precision - we do expect numerical trouble, and 
!  must be prepared to deal with it.
!
!  Anyways, the resulting basis sets are more trouble than it's worth:
!  they require quad accuracy -inside- integral evaluation part to converge
!  at all. Oops.
!
  module basis_recontract
    use accuracy
    use timer
    use math
    use import_gamess
    use gamess_internal
    use integral_tools
    use fock_tools
    use ecp_convert_gamess
    use sort_tools
    use basis_cap
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
    !
    !  ==== Intermediate accuracy parameter ====
    !
    integer, parameter     :: xk           = xrk            ! Real kind used for SCF
    character(len=20), parameter :: global_math = 'quad'    ! Better match the choice maid by xk, or things will go seriously wrong
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
    real(rk)            :: charge          = 0              ! Overall molecular charge; non-integer values 
                                                            ! are OK and will trigger FON-SCF. Must be consistent 
                                                            ! with mo_occ(:) below
    integer(ik)         :: grid_nrad       = 120_ik         ! Basic number of radial points in atomic spheres; actual number of points 
                                                            ! may depend on this value and atom types
    integer(ik)         :: grid_nang       = 770_ik         ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik)         :: grid_outer_nrad =1200_ik         ! Basic number of radial points in the outer sphere
    integer(ik)         :: grid_outer_nang = 770_ik         ! Angular points in the outer sphere
    integer(ik)         :: max_scf_iter    = 100            ! Max number of SCF iterations
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
    logical             :: absorb_fieldfree= .true.         ! Include absorbing boundary in the field-free case
    logical             :: skip_2e         = .false.        ! Omit all 2-electron terms; this flag is strictly for debugging!
    integer(hik)        :: iosize_2e       =220000000_hik   ! Integral I/O buffer size in words for conventional SCF
    character(len=clen) :: scf_type        = 'conventional' ! Can be one of:
                                                            ! 'direct'        - compute all 2e AO integrals on the fly
                                                            ! 'incore'        - keep 2e integrals on disk
                                                            ! 'conventional'  - store integrals on disk
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
    character(len=clen) :: occ_file        = ' '            ! File containing list of -all- spin-orbital occupation numbers; 
                                                            ! blank means that occupations follow the main namelist in the
                                                            ! input file. 
    real(rk)            :: ecp_eps_min     = 1e-6_rk        ! Small shift value cut-off (absolute) in ECP
    real(rk)            :: ecp_eps_max     = 1e+6_rk        ! Large shift value cut-off (positive) in ECP
    real(rk)            :: ecp_eps_grid    = 1e-6_rk        ! Characteristic grid spacing cut-off in ECP
    character(len=clen) :: ecp_report      = ' '            ! File to report ECP level-shift projectors to; empty name
                                                            ! suppresses the output.
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type(gam_structure)         :: gam               ! Structure descriptor; for use with both the field-free and static-field calculations
    type(gam_structure)         :: recon             ! Structure descriptor for the recontracted basis
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
    real(xk), allocatable       :: mo_occ(:)         ! Occupation numbers of molecular orbitals in FON-SCF. The numbers 
                                                     ! are expected to be in the [0:1] range, with odd entries being
                                                     ! the alpha-spin orbitals. Even entries are the betas.
    complex(xk), allocatable    :: mo_energy(:)      ! Energies of molecular orbitals
    complex(xk), allocatable    :: mos (:,:,:)       ! Molecular orbitals for the currently active SCF. 
                                                     ! For all MO arrays, the indices are:
                                                     !   1 = spin-AO
                                                     !   2 = spin-MO
                                                     !   3 = 1 for the left eigenvectors; 2 for the right eigenvectors
    complex(xk), allocatable    :: mosg(:,:,:)       ! Guess molecular orbitals
    complex(xk), allocatable    :: mos0(:,:,:)       ! Converged (or converging ;) molecular orbitals in the absence of the field
    complex(xk), allocatable    :: mosf(:,:,:)       ! ditto, in the field
    complex(xk), allocatable    :: rho (:,:)         ! Electronic density matrix (AO basis)
    complex(xk), allocatable    :: rho_old(:,:)      ! Electronic density matrix from previous SCF iteration
    real(xk), allocatable       :: smat_raw(:,:)     ! Overlap matrix (AO basis), no null-space projection
    real(xk), allocatable       :: smat(:,:)         ! Overlap matrix (AO basis), null-space projected out
    real(xk), allocatable       :: sphalf(:,:)       ! S^{+1/2}, null-space is projected out
    real(xk), allocatable       :: smhalf(:,:)       ! S^{-1/2}, null-space is projected out
    real(xk), allocatable       :: h0  (:,:)         ! Field-free part of the current Hamiltonial, including the ECP terms
    complex(rk), allocatable    :: h0c (:,:)         ! CAP contributions
    complex(xk), allocatable    :: hmat(:,:)         ! Current 1-electron Hamiltonian matrix
    complex(xk), allocatable    :: gmat(:,:)         ! 2-electron contribution to the Fock matrix
    complex(xk), allocatable    :: fmat(:,:)         ! Fock matrix
    complex(xk)                 :: escf              ! SCF electronic energy; since Hamiltonian is non-self-adjoint, energy is complex
    !
    !  Variables related to 2-electron integral handling
    !
    type(int2e_cache)   :: int2e              ! Currently active integrals context
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /recon_data/ verbose, molecule_file, mos_file, charge, occ_file, &
                          max_scf_iter, scf_mixing, scf_eps_energy, scf_eps_rho, &
                          max_diis_nvec, max_diis_coeff, &
                          eps_smat, eps_geev, absorb_fieldfree, &
                          ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
                          scf_type, iosize_2e, skip_2e, &
                          cap_type, cap_centre, cap_r0, &
                          cap_strength, cap_order, &
                          cap_lambda, cap_theta, cap_mpole, cap_efield, cap_diff_step, &
                          grid_nrad, grid_nang, grid_outer_nrad, grid_outer_nang
    !
    contains
    !
    subroutine load_gamess_data
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
      !
      !  A bit of sanity checking
      !
      if (nel_scf<0.0_rk .or. nel_scf>2.0_rk*nao) then
        write (out,"('Number of electrons (',f0.5,') is strange.')") nel_scf
        call stop('basis_recontract%load_gamess_data - bad electron count')
      end if
      !
      !  Tell a bit more!
      !
      write (out,"('                    Number of atoms = ',i0)") natoms
      write (out,"(' Number of spinless basis functions = ',i0)") nao
      write (out,"('                Number of electrons = ',f0.5)") nel_scf
      write (out,"('                       Total charge = ',f0.5)") charge
      write (out,"()")
      if (natoms/=1) call stop('basis_recontract: must have an atom')
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
      !
      allocate (mo_occ(nmo), mo_energy(nmo), mos(nao_spin,nmo,2), &
                rho(nao_spin,nao_spin), rho_old(nao_spin,nao_spin), &
                smat(nao_spin,nao_spin), smat_raw(nao_spin,nao_spin), &
                sphalf(nao_spin,nao_spin), smhalf(nao_spin,nao_spin), &
                h0(nao_spin,nao_spin), h0c(nao_spin,nao_spin), &
                hmat(nao_spin,nao_spin), gmat(nao_spin,nao_spin), &
                fmat(nao_spin,nao_spin), stat=alloc)
      if (alloc/=0) then
        write (out,"('basis_recontract%allocate_dynamic_data: Error ',i0,' allocating quadratic arrays. nao = ',i0)") &
               alloc, nao
        call stop('basis_recontract%allocate_dynamic_data - out of memory (1)')
      end if
      call diis_initialize(diis_st,nao_spin,max_diis_nvec,max_diis_coeff,diis_math=global_math)
    end subroutine allocate_dynamic_data
    !
    subroutine fill_occupations
      integer(ik)         :: iu      ! Can be either input or iu_temp, depending on whether occ_file is blank or not
      integer(ik)         :: ios
      integer(ik)         :: imo, kmo
      !
      mo_occ = 0
      !
      !  We are given explicit occupation numbers; read them in.
      !
      iu = input
      if (occ_file/=' ') then
        iu = iu_temp
        open(iu,file=trim(occ_file),status='old',action='read',iostat=ios)
        if (ios/=0) then
          write (out,"('Error ',i0,' opening MO occupation file ',a)") ios, trim(occ_file)
          call stop('basis_recontract%fill_occupations - bad open')
        end if
      end if
      read (iu,*,iostat=ios) mo_occ
      if (occ_file/=' ') close(iu)
      !
      !  Bit of sanity checking ...
      !
      if (any(mo_occ>1.0) .or. any(mo_occ<0.0)) then
        write (out,"('MO occupation numbers must be in the range [0:1]')")
        write (out,"((10(1x,f10.5)))") mo_occ
        call stop('basis_recontract%fill_occupations - bad occupations (1)')
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
        call stop('basis_recontract%fill_occupations - bad occupations (2)')
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
    !  Since our integrals package operates with spin-free AOs,
    !  we need to do a bit of juggling around to accommodate the
    !  spin-AOs we use here.
    !
    subroutine evaluate_1e_hamiltonian
      integer(ik)                 :: iat, alloc
      real(xk)                    :: xyz(3), q
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
        write (out,"('basis_recontract%field_free_1e_hamiltonian: Error ',i0,' allocating spin-free temporaries')") alloc
        call stop('basis_recontract%field_free_1e_hamiltonian - no memory')
      end if
      !
      if (verbose>=0) then
        write (out,"(/'ECP and CAP terms are evaluated using kind-',i0,' arithmetics.')") rk
        write (out,"( 'All other 1-electron terms are calcuated using kind-',i0,' floating point.'/)") xk
      end if
      !
      call gamess_1e_integrals('AO KINETIC',hsf,bra=gam,ket=gam)
      !
      !  Fill nuclear attraction integrals.
      !
      nuclear_attraction: do iat=1,natoms
        xyz = real(gam%atoms(iat)%xyz,kind=kind(xyz)) / abohr
        q   = real(gam%atoms(iat)%znuc,kind=kind(q))
        call gamess_1e_integrals('AO 3C 1/R',tmp_xk,bra=gam,ket=gam,op_xyz=xyz)
        tmp_xk = -q * tmp_xk
        hsf = hsf + tmp_xk
      end do nuclear_attraction
      !
      tmp_rk = 0
      call ecp_evaluate_matrix_elements(gam,ecp,tmp_rk)  ! Put ECP terms in tmp
      hsf = hsf + tmp_rk
      !
      call cap_evaluate(gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang,csf,ssfn,lsfn)  ! CAP terms are done numerically
      ! 
      !  For a non-local CAP, check Laplacian matrix elements
      !
      call gamess_1e_integrals('AO OVERLAP',ssf,   bra=gam,ket=gam)
      !
      ssf = 0.5_rk * (ssf + transpose(ssf))
      !
      !  Expand 1-electron matrix elements to spin-AO basis
      !
      h0   = 0 ; h0c = 0 ; smat = 0
      h0  (    1:nao     ,    1:nao     ) = hsf
      h0  (nao+1:nao_spin,nao+1:nao_spin) = hsf
      h0c (    1:nao     ,    1:nao     ) = csf
      h0c (nao+1:nao_spin,nao+1:nao_spin) = csf
      smat(    1:nao     ,    1:nao     ) = ssf
      smat(nao+1:nao_spin,nao+1:nao_spin) = ssf
      smat_raw = smat
      deallocate (hsf,ssf,lsfn,ssfn,fsf,tmp_xk,tmp_rk)
      !
      call TimerStop('1e Hamiltonian')
    end subroutine evaluate_1e_hamiltonian
    !
    !  Hartree-Fock electronic energy
    !
    subroutine total_energy
      complex(xk) :: efock, eg

      efock = sum(rho * fmat)
      eg    = sum(rho * gmat)
      escf  = efock - 0.5_xk * eg
      if (verbose>=0) then
        write (out,"(/'Total SCF energy =',2(1x,g20.12)/)") escf
      end if
    end subroutine total_energy
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
        escf_old = escf
        !
        write (out,"('Beginning SCF cycle ',i4)") iter
        call flush(out)
        call st_density_matrix(mo_occ,mos,rho)
        !
        gmat = 0
        if (iter==1) then
          write (out,"(/'First iteration, starting with a bare nuclear Hamiltonian guess'/)") 
        else 
          if (.not.skip_2e) call fock_g_matrix(int2e,rho,gmat)
        end if
        fmat = hmat + gmat
        !
        call total_energy
        call flush(out)
        !
        if (iter>1) then
          !
          !  Since we are starting with bare nucleus Hamiltonian, the first iteration
          !  simply prepares the orbital guess. If we feed it into diis_extrapolate,
          !  we'll get false convergence. Ooops.
          !
          call diis_extrapolate(diis_st,iter-1,smat,rho,fmat)
        end if
        !
        call st_diagonalize_fmat(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy)
        call flush(out)
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
          write (out,"('Final SCF energy = ',g24.14,1x,g24.14)") escf
          call flush(out)
          exit repeat_scf
        end if
      end do repeat_scf
      !
      if (.not.converged) then
        call stop('basis_recontract%scf_loop - SCF convergence failure')
      end if
      !
      call TimerStop('SCF loop')
      call TimerReport
    end subroutine scf_loop
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
        call stop('basis_recontract%punch_final_mos - Can''t create file')
      end if
      write (iu_mos,"('Field-free MOs for ',a)") trim(molecule_file)
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
    subroutine merge_primitives(cnt,z,c)
      integer(ik), intent(out) :: cnt  ! Number of surviving primitives
      real(xrk), intent(inout) :: z(:) ! Exponents
      real(xrk), intent(inout) :: c(:) ! Coefficients
      !
      integer(ik) :: order(size(z))
      integer(ik) :: inp
      !
      !  First, sort the basis functions
      !
      call order_keys(-z,order)
      z = z(order)
      c = c(order)
      !
      cnt = min(1,size(z))
      merge_inputs: do inp=2,size(z)
        if(abs(z(inp)-z(cnt))<=1e3*spacing(z(inp))) then
          !
          !  Same exponent; merge the coefficient
          !
          c(cnt) = c(cnt) + c(inp)
        else
          !
          !  New exponent; copy to the output
          !
          cnt = cnt + 1
          z(cnt) = z(inp)
          c(cnt) = c(inp)
        end if
      end do merge_inputs
    end subroutine merge_primitives
    !
    subroutine remember_contraction(src,dst,lv,c)
      type(gam_atom), intent(in)    :: src
      type(gam_atom), intent(inout) :: dst
      integer(ik), intent(in)       :: lv      ! Current value of the angular momentum
      real(xk), intent(in)          :: c(:)    ! Contraction coefficient of the new basis function
      !
      integer(ik) :: ic      ! Current contraction
      integer(ik) :: src_sh  ! Shell index, input atom
      integer(ik) :: dst_sh  ! ditto, output 
      integer(ik) :: src_p   ! Primitive index, input atom
      integer(ik) :: dst_p   ! ditto, output
      integer(ik) :: dst_p0  ! Starting primitive on the output atom
      integer(ik) :: dst_cnt ! Number of primitives after optimization
      !
      dst_sh = dst%nshell + 1
      dst_p0 = dst%sh_p(dst_sh)
      dst_p  = dst_p0
      if (dst_sh>size(dst%sh_l)) call stop('basis_recontract%remember_contraction - too many shells')
      ic = 0
      scan_source_shells: do src_sh=1,src%nshell
        if (src%sh_l(src_sh)/=lv) cycle scan_source_shells
        !
        ic = ic + 1
        scan_source_primitives: do src_p=src%sh_p(src_sh),src%sh_p(src_sh+1)-1
          if (dst_p>size(dst%p_zet)) call stop('basis_recontract%remember_contraction - too many primitives')
          dst%p_zet(dst_p) = src%p_zet(src_p)
          dst%p_c  (dst_p) = src%p_c  (src_p) * c(ic)
          dst_p = dst_p + 1
        end do scan_source_primitives
      end do scan_source_shells
      if (ic/=size(c)) call stop('basis_recontract%remember_contraction - count error')
      call merge_primitives(dst_cnt,dst%p_zet(dst_p0:dst_p-1),dst%p_c(dst_p0:dst_p-1))
      dst_p = dst_p0 + dst_cnt
      dst%sh_l(dst_sh)   = lv
      dst%sh_p(dst_sh+1) = dst_p
      dst%nshell         = dst_sh
    end subroutine remember_contraction
    !
    !  Assign orbitals to degenerate sub-blocks
    !
    subroutine classify_degeneracy(eval,nset,set1,setn)
      real(xk), intent(in)     :: eval(:) ! Eigenvalues
      integer(ik), intent(out) :: nset    ! Number of degenerate sets of orbitals found
      integer(ik), intent(out) :: set1(:) ! First orbital of each set
      integer(ik), intent(out) :: setn(:) ! Last orbital of each set
      !
      integer(ik) :: iev1, iev2, iev3, nev
      !
      nset = 0
      iev1 = 1
      nev  = size(eval)
      scan_for_degeneracies: do while (iev1<=nev)
        iev3 = iev1
        iev2 = iev3 ! This assignment is not needed, since the loop below
                    ! will always execute at least once; however, it helps
                    ! to shut off a spurious warning from gfortran
        find_identical_eigenvalues: do while(iev3<=nev)
          if (abs(eval(iev3)-eval(iev1))>eps_geev) exit find_identical_eigenvalues
          iev2 = iev3
          iev3 = iev3 + 1
        end do find_identical_eigenvalues
        !
        nset = nset + 1
        set1(nset) = iev1
        setn(nset) = iev2
        !
        iev1 = iev2 + 1
      end do scan_for_degeneracies
    end subroutine classify_degeneracy
    !
    subroutine build_one_contraction(lv,ncntr,s,mos)
      integer(ik), intent(in) :: lv             ! Current value of the angular momentum
      integer(ik), intent(in) :: ncntr          ! Expected contraction length
      real(xk), intent(in)    :: s(:,:)         ! Overlap matrix for sub-block
      real(xk), intent(in)    :: mos(:,:)       ! Orbital coefficients for the degenerate set
      !
      integer(ik) :: ic, ib, ic2, ib2
      integer(ik) :: nc          ! Number of Cartesian functions for this lv
      real(xk)    :: cs(2*ncntr) ! Contraction coefficients for spin-full functions
      real(xk)    :: c (ncntr)   ! Contraction coefficients for spin-less functions
      real(xk)    :: w
      !
      nc = gam_orbcnt(lv)
      unnormalized_contractions: do ic=1,2*ncntr
        ib = 1 + (ic-1)*nc
        cs(ic) = sum(mos(ib:ib+nc-1,:))
      end do unnormalized_contractions
      c = cs(1:ncntr) + cs(ncntr+1:) ! ?????
      ! write (out,"('Unnormalized contractions:',10(1x,g16.8))") c
      !
      !  It does not matter which Cartesian component is used for normalization; choose the first.
      !
      w = 0
      contraction_norm: do ic2=1,ncntr
        ib2 = 1 + (ic2-1)*nc
        do ic=1,ncntr
          ib = 1 + (ic-1)*nc
          w  = w + s(ib,ib2) * c(ic) * c(ic2)
        end do
      end do contraction_norm
      ! write (out,"('Raw norm: ',g16.8)") w
      c = c / sqrt(w)
      !
      !  Make the largest coefficient positive, and we are done
      !
      ib = maxloc(abs(c),dim=1)
      if (c(ib)<0) c = -c
      ! write (out,"('  Normalized contractions:',10(1x,g16.8))") c
      call remember_contraction(gam%atoms(1),recon%atoms(1),lv,c)
    end subroutine build_one_contraction
    !
    subroutine find_degenerate_sets(lv,ncntr,s,en,mos)
      integer(ik), intent(in) :: lv             ! Current value of the angular momentum
      integer(ik), intent(in) :: ncntr          ! Expected contraction length
      real(xk), intent(in)    :: s(:,:)         ! Overlap matrix for sub-block
      real(xk), intent(in)    :: en(:)          ! Eigenvalues of the orbitals
      real(xk), intent(in)    :: mos(:,:)       ! Orbital coefficients
      !
      integer(ik) :: iset, mult, lx
      integer(ik) :: nev                 ! Number of non-null eigenvectors
      integer(ik) :: nset                ! Number of degenerate sets
      integer(ik) :: set1(nao_spin)      ! First orbital of each set
      integer(ik) :: setn(nao_spin)      ! Last orbital of each set
      !
      nev = size(en)
      call classify_degeneracy(en,nset,set1,setn)
      !
      write (out,"()")
      write (out,"(1x,a5,1x,a22,2x,a7,1x,a7,1x,a7)") &
          ' Set ', '  Re(eps)  ', ' First ', ' Last ', ' Mult ', &
          '-----', '-----------', '-------', '------', '------'
      !
      print_sets: do iset=1,nset
        mult = setn(iset)-set1(iset)+1
        write (out,"(1x,i5,1x,f22.12,2x,i7,1x,i7,1x,i7)") &
               iset, en(set1(iset)), set1(iset), setn(iset), mult
        lx = (mult-2)/4
        if (2*(2*lx+1)/=mult) then
          call stop('basis_recontract%build_recontractions - Degeneracy is not in the 2*(2*L+1) form')
        end if
        if (lx>lv .or. lx<0 .or. mod(lv-lx,2)/=0) then
          call stop('basis_recontract%build_recontractions - Angular momentum out of range')
        end if
        if (lx==lv) then
          !
          !  This degenerate set has the correct multiplicity; build contraction coefficients
          !
          call build_one_contraction(lv,ncntr,s,mos(:,set1(iset):setn(iset)))
        end if
      end do print_sets
      write (out,"()")
      call flush(out)
    end subroutine find_degenerate_sets
    !
    subroutine classify_basis(atm,basis_l)
      type(gam_atom), intent(in) :: atm        ! Atom descriptor
      integer(ik), intent(out)   :: basis_l(:) ! Nominal L value of each basis function
      !
      integer(ik) :: ish, ib1, ibn, lv
      !
      ibn = 0
      scan_shells: do ish=1,atm%nshell
        lv  = atm%sh_l(ish)
        ib1 = ibn + 1
        ibn = ib1 + gam_orbcnt(lv) - 1
        basis_l(ib1:ibn) = lv
      end do scan_shells
      if (ibn/=size(basis_l)) call stop('basis_recontract%classify_basis - count error')
    end subroutine classify_basis
    !
    !  Our strategy for recontracting the basis is to:
    !
    !  1. Separate the final Fock matrix into blocks corresponding to nominal
    !     L values of the basis functions. Since we use Cartesian gaussians,
    !     the separation is not perfect; the remaining off-diagonal terms are
    !     quietly ignored.
    !  2. For each sub-block, drop the non-real, non-symmetric part and diagonalize
    !  3. Choose the eigenvalues of the nominal degeneracy; these define our
    !     recontraction coefficients
    !
    subroutine recontract_basis
      integer(ik)              :: lv, io, iev
      integer(ik)              :: basis_l(nao_spin)    ! Nominal L value of each basis function
      integer(ik)              :: lsize, lnull, lev
      integer(ik)              :: ncntr                ! Contraction length
      integer(ik), allocatable :: ao_index(:)          ! List of AOs which belong to a given nominal L
      complex(xk), allocatable :: fblock_cmplx(:,:)    ! Sub-block of the Fock matrix
      real(xk), allocatable    :: fblock_real (:,:)    ! Symmetrized sub-block of the Fock matrix
      real(xk), allocatable    :: sblock(:,:)          ! Sub-block of the overlap
      real(xk), allocatable    :: sblock_ph(:,:)       ! sblock**(+0.5)
      real(xk), allocatable    :: sblock_mh(:,:)       ! sblock**(-0.5)
      complex(xk), allocatable :: moblock_cmplx(:,:,:) ! Eigenvectors of the sub-block
      complex(xk), allocatable :: moen_cmplx(:)        ! Eigenvalues of the sub-block
      real(xk), allocatable    :: moblock_real(:,:)    ! Eigenvectors of the symmetrized sub-block
      real(xk), allocatable    :: moen_real(:)         ! Eigenvalues of the symmetrized sub-block
      !
      !  Prepare output structure by copying (gam) and deleting the basis
      !
      allocate (recon%atoms(1))
      recon%natoms          = 1
      recon%atoms(1)        = gam%atoms(1)
      recon%atoms(1)%nshell = 0
      !
      !  Assign nominal L value to each basis function
      !
      call classify_basis(gam%atoms(1),basis_l(:nao))
      basis_l(nao+1:) = basis_l(:nao)
      !
      scan_l_values: do lv=0,ubound(gam_orbcnt,dim=1)
        lsize = count(basis_l==lv)
        ncntr = lsize / (2*gam_orbcnt(lv)) ! We are working with spin-orbitals, hence the factor of 2
        write (out,"(/'Number of basis functions with L=',i0,' is ',i0)") lv, lsize
        write (out,"( '  Expected contraction length is ',i0)") ncntr
        if (ncntr*2*gam_orbcnt(lv)/=lsize) call stop('basis_recontract%recontract_basis - function count error')
        allocate (ao_index(lsize),fblock_cmplx(lsize,lsize),fblock_real(lsize,lsize), &
                  sblock(lsize,lsize),sblock_ph(lsize,lsize),sblock_mh(lsize,lsize), &
                  moblock_cmplx(lsize,lsize,2),moen_cmplx(lsize), &
                  moblock_real(lsize,lsize),moen_real(lsize))
        ao_index = pack((/(io,io=1,size(basis_l))/),basis_l==lv)
        !
        !  Extract sub-block of the Fock operator
        !
        fblock_cmplx = fmat(ao_index,ao_index)
        fblock_real  = real(fblock_cmplx,kind=xk)
        fblock_real  = 0.5_xk * (fblock_real + transpose(fblock_real))
        !
        write (out,"('Maximum change upon extraction of real, symmetric Fock matrix subblock is: ',g20.10)") &
               maxval(abs(fblock_real-fblock_cmplx))
        !
        !  Prepare sub-block overlaps and needed powers
        !
        sblock = smat_raw(ao_index,ao_index)
        call st_invert_smat(lnull,sblock,sblock_mh,sblock_ph,eps_smat=eps_smat)
        lev = lsize - lnull
        !
        !  Diagonalize both the original and symmetrized sub-blocks
        !
        call st_diagonalize_fmat(sblock_mh,sblock_ph,fblock_cmplx,lnull,eps_geev,moblock_cmplx,moen_cmplx)
        call st_diagonalize_fmat(sblock_mh,sblock_ph,fblock_real,lnull,eps_geev,moblock_real,moen_real)
        !
        write (out,"('Block and symmetrized-block eigenvalues:')")
        write (out,"(1x,a6,2x,a22,1x,a22,2x,a22,2x,a22)") &
                   ' iev ', '  Re(eps)  ', '  Im(eps)  ',  '  eps-symm  ', '  abs(diff)  ', &
                   '-----', '-----------', '-----------',  '------------', '-------------'
        print_block_eigenvalues: do iev=1,lev
          write (out,"(1x,i6,2x,f22.14,1x,f22.14,2x,f22.14,2x,e22.14)") &
                 iev, moen_cmplx(iev), moen_real(iev), abs(moen_cmplx(iev)-moen_real(iev))
        end do print_block_eigenvalues 
        write (out,"()")
        !
        call find_degenerate_sets(lv,ncntr,sblock,moen_real(:lev),moblock_real(:,:lev))
        !
        deallocate (ao_index,fblock_cmplx,fblock_real,sblock,sblock_ph,sblock_mh)
        deallocate (moblock_cmplx,moen_cmplx,moblock_real,moen_real)
      end do scan_l_values
      !
      !  Finish construction of the recontracted basis and report
      !  The gamess structure we prepare is not fully internally consistent
      !  (do not use it for integrals!), but it should be good enough for printing.
      !
      call report_gamess_data(out,recon,'Recontracted basis set (remember the 80-column limit in GAMESS!)')
    end subroutine recontract_basis
    !
    subroutine start
      integer(ik) :: info
      !
      call TimerStart('start')
      !
      read (input,nml=recon_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=recon_data)
      write (out,"()")
      !
      call load_gamess_data
      !
      call allocate_dynamic_data
      !
      call prepare_ecps
      !
      call fill_occupations
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
      call st_invert_smat(nmo_null,smat,smhalf,sphalf,eps_smat=eps_smat)
      call flush(out)
      !
      if (.not.skip_2e) call prepare_2e(int2e,gam,scf_type,iu_2e_ao,iosize_2e,ints_math=global_math)
      !
      hmat = h0
      if (absorb_fieldfree) hmat = hmat + h0c
      call scf_loop
      call punch_final_mos(0_ik)
      !
      call recontract_basis
      !
      if (.not.skip_2e) call clear_2e(int2e)
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
  end module basis_recontract
!
  program driver
    use basis_recontract
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
