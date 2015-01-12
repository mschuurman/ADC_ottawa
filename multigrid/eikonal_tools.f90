!
!  Common subroutines used generation of the eikonal wavefunctions in
!  adiabatic approximation.
!
!  Added Feb 2010 (Michael Spanner): Support for building transition densities and associated transition potentials.
!
!
  module eikonal_tools
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use dft
    use math
    implicit none
    private
    public eikonal_build_potential, eikonal_build_transition_potential, eikonal_build_function
    !
    integer(ik), parameter      :: verbose = 2
    !
    !  Local static data
    !
    integer(ik), save           :: core_nnuc      = 0          ! Number of nuclei in the molecular core
    real(rk), save              :: core_electrons = 0.0_rk     ! Total electronic charge from natural orbitals
    real(rk), allocatable, save :: core_xyzq(:,:)              ! Coordinates and charges of the nuclei in the core
    real(rk), save              :: total_charge                ! Total charge of the ion
    !
    real(rk), save     :: multipoles_nuclear(9)                ! Multipole moments - needed for asymptotic form of
    real(rk), save     :: multipoles_electronic(9)             !                the potential and scattering phase
    real(rk), save     :: multipoles_total(0:9)                ! Total multipole moment
    !
    contains
    !
    !  Build electron density on grid from the natural orbitals on file.
    !
    subroutine build_electron_density(f_rho,f_table,natural_file,natural_occ,rot)
      integer(ik), intent(in)      :: f_rho          ! Total electron density
      integer(ik), intent(in)      :: f_table(:)     ! Scratch fields
      character(len=*), intent(in) :: natural_file   ! Natural orbitals
      real(rk), intent(in)         :: natural_occ(:) ! Natural occupations
      real(ark), intent(in)        :: rot(:,:)       ! Rotation matrix
      !
      real(rk)                 :: norm
      integer(ik)              :: natural_count
      integer(ik)              :: n_scr, ino_1, ino_n, ino, no_count, islot, alloc
      integer(ik), allocatable :: no_ind(:)  ! Indices of the natural orbitals in this batch
      !
      call TimerStart('Build electron density')
      !
      !  There may be a large number of natural orbitals, and we do not want
      !  to load them all at the same time. For this reason, we'll process
      !  NOs in batches, using all available memory as scratch
      !
      natural_count = size(natural_occ)
      n_scr         = size(f_table)
      !
      allocate (no_ind(n_scr),stat=alloc)
      if (alloc/=0) then
        write (out,"('eikonal_tools: Error ',i8,' allocating ',i8,'-element integer index array')") alloc, n_scr
        stop 'eikonal_tools%build_electron_density'
      end if
      call FieldZero(f_rho)
      !
      process_no_batch: do ino_1=1,natural_count,n_scr
        ino_n    = min(ino_1+n_scr-1,natural_count)
        no_count = ino_n - ino_1 + 1
        fill_no_ind: do ino=ino_1,ino_n
          no_ind(ino-ino_1+1) = ino
        end do fill_no_ind
        !
        !  Load the next batch of the NOs 
        !
        if (verbose>=0) then
          write (out,"('Loading natural orbitals ',i0,' through ',i0)") ino_1, ino_n
        end if
        call FieldImport('GAMESS',natural_file,f_table(:no_count),no_ind(:no_count),rot=rot)
        !
        !  Normalize the NOs to 1 to reduce numerical noise, and accumulate total
        !  density.
        !
        accumulate_density: do ino=ino_1,ino_n
          islot = f_table(ino-ino_1+1)
          call QMNormalize(wf=islot,norm_desired=1._rk,norm_original=norm)
          if (verbose>=1) then
            write (out,"(t5,'Natural orbital ',i4,' had norm of ',f20.10)") ino, norm**2
          end if
          call FieldRhoAccumulate(alpha=natural_occ(ino),src=islot,dst=f_rho)
        end do accumulate_density
      end do process_no_batch
      deallocate (no_ind)
      call TimerStop('Build electron density')
    end subroutine build_electron_density
    !
    !  Build electron transition density on grid from the RDM orbitals on file.
    !
    !  Added Feb 2, 2010 (Michael Spanner)
    !
    subroutine build_electron_transition_density(f_rho,f_table,rdm_file,rdm_sv,rot,transition)
      integer(ik), intent(in)      :: f_rho          ! Total electron (pseudo-)density
      integer(ik), intent(in)      :: f_table(:)     ! Scratch fields
      character(len=*), intent(in) :: rdm_file       ! File containing vectors for reduced density matrix deomposition 
      real(rk), intent(in)         :: rdm_sv(:)      ! singular values of the 1-RDM decomposition
      real(ark), intent(in)        :: rot(:,:)       ! Rotation matrix
      logical, intent(in)          :: transition     ! True if we are calculating transition density
      !
      real(rk)    :: norm
      integer(ik) :: rdm_count              ! Number of left/right singular vector pairs in the RDM
      integer(ik) :: n_scr                  ! Number of scratch fields available for loading RDM
      integer(ik) :: rdm_batch              ! Max. size of the RDM vectors batch
      integer(ik) :: ir1, irn               ! First and last RDM vector pairs in a batch
      integer(ik) :: nrn                    ! Number of RDM vector pairs in a batch
      integer(ik) :: irdm                   ! RDM vector pair being processed now
      integer(ik) :: nn, f_left, f_right
      integer(ik) :: rdm_ind(2*size(rdm_sv))  ! Indices of the RDM orbitals for current batch. The size is
                                              ! likely an overkill, but we don't really care at this point.
      !
      call TimerStart('Build electron transition density')
      !
      !  There may be a large number of RDM orbitals, and we do not want
      !  to load them all at the same time. For this reason, we'll process
      !  RDMs in batches, using all available memory as scratch
      !
      rdm_count = size(rdm_sv)
      n_scr     = size(f_table)
      !
      if (n_scr<2) then
        write (out,"('eikonal_tools%build_electron_transition_density: Error, not enough scratch fields, need at least 2')")
        stop 'eikonal_tools%build_electron_transition_density'
      end if
      !
      rdm_batch = n_scr / 2
      !
      call FieldZero(f_rho)
      if (verbose>=0) write (out,"()")
      !
      load_rdm_batch: do ir1=1,rdm_count,rdm_batch
        irn = min(rdm_count,ir1+rdm_batch-1)
        nrn = irn - ir1 + 1
        !
        !  Load batch of the RDM vectors in pairs (ie as many left/right vectors as we can fit in memory)
        !  Orbital loading has a large up-front cost, and batching helps to amortize it.
        !
        fill_rdm_indices: do irdm=ir1,irn
          rdm_ind(2*(irdm-ir1)+1) = 2*irdm-1
          rdm_ind(2*(irdm-ir1)+2) = 2*irdm
        end do fill_rdm_indices
        if (verbose>=0) write(out,"('Loading 1-RDM singular vectors',i5,' through ',i5)") ir1, irn
        call FieldImport('GAMESS',rdm_file,f_table(:2*nrn),rdm_ind(:2*nrn),rot=rot)
        !
        !  Data is in memory, accumulate the 1-RDM from the singular vectors
        !
        process_rdm_batch: do irdm=ir1,irn
          !
          !  Normalize singular vectors to 1 to reduce numerical noise, and 
          !  accumulate total density.
          !
          normalize_pair: do nn=1,2
            call QMNormalize(wf=f_table(2*(irdm-ir1)+nn),norm_desired=1._rk,norm_original=norm)
            if (verbose>=1) then
              write (out,"(t5,'1-RDM singular vector',i4,'.',i1,' had norm of ',f20.10)") irdm, nn, norm**2
            end if
          end do normalize_pair
          !
          f_left  = f_table(2*(irdm-ir1)+1)
          f_right = f_table(2*(irdm-ir1)+2)
          call FieldScale(dst=f_left,con=cmplx(rdm_sv(irdm),0._rk,kind=rk))
          call FieldMulAdd(src_a=f_left,src_b=f_right,dst=f_rho)
          !
        end do process_rdm_batch
      end do load_rdm_batch
      !
      if (transition) then
        !
        !  Sanity check on transition density
        !
        norm = abs(FieldNorm1(src=f_rho)) ! FieldNorm1 returns a complex result
        if (verbose>=1) then
          write (out,"()")
          write (out,"(t5,a,': Integrated transition density = ',f20.10)") trim(rdm_file), norm
          write (out,"()")
        end if
        if (norm>1e-4) then
          write (out,"()")
          write (out,"('******** WARNING ******** Integrated transition density is kinda big...  Maybe a problem? ')")
          write (out,"()")
        end if
        if (norm>1._rk) then
          stop 'build_electron_transition_density - integrated density too inaccurate'
        end if
      end if
      call TimerStop('Build electron transition density')
    end subroutine build_electron_transition_density
    !
    !  Store and report nuclear coordinates and charges used for the potential
    !
    subroutine get_nuclear_coordinates
      integer(ik) :: iat
      !
      call gamess_report_nuclei(core_nnuc)
      write (out,"('Density data file contained ',i5,' nuclei')") core_nnuc
      if (.not.allocated(core_xyzq)) allocate (core_xyzq(4,core_nnuc))
      call gamess_report_nuclei(core_nnuc,core_xyzq)
      total_charge = sum(core_xyzq(4,:))-core_electrons
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,core_nnuc
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") core_xyzq(4,iat), core_xyzq(1:3,iat), core_xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
      write (out,"('       Total core charge = ',f16.8)") sum(core_xyzq(4,:))
      write (out,"(' Total charge of the ion = ',f16.8)") total_charge
      write (out,"()")
      !* fixed: S.P. Jan 16, 2008
      ! if (abs(total_charge-1.0_rk)>1e-4_rk) then
      !   write (out,"('FIXME: Eikonal boundary conditions are not implemented for ions of charge/=1')")
      !   stop 'boundary conditions not implemented'
      ! end if
    end subroutine get_nuclear_coordinates
    !
    !  Build total multiplicative potential of the ion. At the moment, we include
    !  potential of the point charges, plus the Hartree potential of the ion core
    !
    subroutine build_hartree_potential(f_rho,f_bare,f_core_pot)
      integer(ik), intent(in) :: f_rho      ! Input: Field containing the density
      integer(ik), intent(in) :: f_bare     ! Output: Bare nuclei potential
      integer(ik), intent(in) :: f_core_pot ! Output: Nuclear + Hartree potential
      !
      integer(ik) :: in1, in2
      real(rk)    :: e_core, e_hartree, e_chh, e_nuc
      complex(rk) :: mult_temp(9)
      !
      call TimerStart('Build Hartree potential')
      !
      !  Potential due to the nuclei will be left in f_scr (probe charge of -1 already in)
      !
      call fock_external_potential(v=f_bare,f_scr=(/f_core_pot/),xyzq=core_xyzq)
      call FLnuclearMultipoles(multipoles_nuclear)
      e_core = real(FieldProductIntegrate(src_a=f_rho,src_b=f_bare),kind=rk)
      !
      !  Potential due to the electrons will be left in f_core_pot
      !
      call fock_hartree_potential(v=f_core_pot,rho=f_rho)
      e_hartree = 0.5_rk * real(FieldProductIntegrate(src_a=f_rho,src_b=f_core_pot),kind=rk)
      call FieldNorm1Multipoles(f_rho,mult_temp)
      multipoles_electronic = real(mult_temp,kind=rk)
      multipoles_total(0)   = total_charge
      multipoles_total(1:9) = multipoles_nuclear - multipoles_electronic ! Electrons are positive - oops.
      !
      !  Combine the two potentials
      !
      call FieldAXPY(alpha=(1._rk,0._rk),src=f_bare,dst=f_core_pot)
      e_chh = real(FieldProductIntegrate(src_a=f_rho,src_b=f_core_pot),kind=rk)
      !
      !  Also calculate nuclear repulsion energy for completeness
      !
      e_nuc = 0
      nuclear_repulsion: do in1=1,core_nnuc
        do in2=1,in1-1
          e_nuc = e_nuc + core_xyzq(4,in1)*core_xyzq(4,in2)/ &
                          sqrt(sum((core_xyzq(1:3,in1)-core_xyzq(1:3,in2))**2))
        end do
      end do nuclear_repulsion
      !
      !  Report integrals for sanity checking
      !
      if (verbose>=0) then
        write (out,"(/' Nuclear repulsion energy = ',f14.7)") e_nuc
        write (out,"( '              Core energy = ',f14.7)") e_core
        write (out,"( '           Hartree energy = ',f14.7)") e_hartree
      ! write (out,"( '    Core+2*Hartree energy = ',f14.7)") e_chh
        write (out,"(/' Dipole moment of the charge distribution [electron-Bohr]')")
        write (out,"( '             ',3(1x,a15))") '   X   ', '   Y   ', '   Z   '
        write (out,"( '    Nuclear: ',3(1x,f15.8))") multipoles_nuclear   (1:3)
        write (out,"( ' Electronic: ',3(1x,f15.8))") multipoles_electronic(1:3)
        write (out,"( '      Total: ',3(1x,f15.8))") multipoles_total     (1:3)
        write (out,"(/' Quadrupole moment of the charge distribution [electron-Bohr^2]')")
        write (out,"( '             ',6(1x,a15))") '   XX   ', '   YY   ', '   ZZ   ', &
                                                   '   XY   ', '   XZ   ', '   YZ   '
        write (out,"( '    Nuclear: ',6(1x,f15.8))") multipoles_nuclear   (4:9)
        write (out,"( ' Electronic: ',6(1x,f15.8))") multipoles_electronic(4:9)
        write (out,"( '      Total: ',6(1x,f15.8))") multipoles_total     (4:9)
      end if
      !
      call TimerStop('Build Hartree potential')
    end subroutine build_hartree_potential
    !
    subroutine build_exchange_correlation(v_xc,f_rho,f_vxc)
      character(len=*), intent(in) :: v_xc  ! Name of the exchange-correlation functional
      integer(ik), intent(in)      :: f_rho ! Input: Field containing total electron density
      integer(ik), intent(in)      :: f_vxc ! Output: Field containing exchange-correlation potential
      !
      real(rk) :: exc
      !
      call dft_set_functional(v_xc)
      !
      !  First, a sanity check: calculate exchange-correlation
      !  energy for the input density.
      !
      if (verbose>=0) then
        call FieldCopy(src=f_rho,dst=f_vxc)
        call FieldProcess(dst=f_vxc,func=dft_epsilon_xc)
        exc = real(FieldNorm1(src=f_vxc),kind=rk)
        write (out,"(/' DFT (',a,') exchange-correlation energy = ',f15.8,' Hartree'/)") &
               trim(v_xc), exc
      end if
      call FieldCopy(src=f_rho,dst=f_vxc)
      call FieldProcess(dst=f_vxc,func=dft_v_xc)
    end subroutine build_exchange_correlation
    !
    !  Constructs multiplicative potential used for building eikonal
    !  phase corrections. The potential is the sum of the core and
    !  Hartree potential, obtained from the natural orbitals produced
    !  by a separate ab initio run.
    !
    subroutine eikonal_build_potential(natural_file,natural_occ,f_table,f_core_pot,f_bare,ion_charge,v_xc, &
                                       multipoles,electronic_multipoles,rot,use_rdms)

      character(len=*), intent(in)           :: natural_file    ! GAMESS-format MO coefficients for the natural orbitals
      real(rk), intent(in)                   :: natural_occ(:)  ! Occupation numbers for the orbitals
      integer(ik), intent(in)                :: f_table(:)      ! List of fields which can be used for scratch
      integer(ik), intent(in)                :: f_core_pot      ! Field where eikonal potential will be stored
      integer(ik), intent(in), optional      :: f_bare          ! Field for the bare-nuclei potential
      real(rk), intent(out), optional        :: ion_charge      ! Total ion charge
      complex(rk), intent(out), optional     :: multipoles(0:9) ! Multipole moments matching total potential,
                                                                ! can be used for the long-range expansion.
      complex(rk), intent(out), optional     :: electronic_multipoles(0:9)
                                                                ! Electonic part of the multipole moments - no
                                                                ! nuclear contributions are included.
      character(len=*), intent(in), optional :: v_xc            ! Exchange-correlation potential to use. Can be:
                                                                ! ' ' or 'none' - no exchange potential (default)
                                                                ! 'Slater' - Dirac-Slater exchange
                                                                ! 'SVWN'   - Slater exchange + VWN local correlation
      real(ark), intent(in), optional         :: rot(:,:)       ! Rotation matrix to use
      logical, intent(in), optional           :: use_rdms       ! Input file is in superdyson 1-RDM format, not GAMESS checkpoint
                                                                ! format
      !
      integer(ik) :: f_rho, f_bare_temp, f_vxc_temp
      real(ark)   :: rotmat(3,3)
      logical     :: rdm
      real(rk)    :: norm
      !
      call TimerStart('EikonalBuildPotential')
      !
      !  Temporarily, use f_core_pot to store the total electron density
      !
      call MathSetUnitMatrix(rotmat)
      if (present(rot)) rotmat = rot
      rdm = .false.
      if (present(use_rdms)) rdm = use_rdms
      if (rdm) then
        call build_electron_transition_density(f_core_pot,f_table,natural_file,natural_occ,rotmat,transition=.false.)
      else
        call build_electron_density(f_core_pot,f_table,natural_file,natural_occ,rotmat)
      end if
      !
      !  Sanity check on the total density
      !
      core_electrons = sum(natural_occ)
      norm           = real(FieldNorm1(src=f_core_pot),kind=rk) ! FieldNorm1 returns a complex result
      write (out,"(/'    Core electrons from NO occupations: ',f20.10)") core_electrons
      write (out,"( 'Core electrons from integrated density: ',f20.10)") norm
      if (abs(norm-core_electrons)>1e-4) then
        stop 'build_electron_density - integrated density too inaccurate'
      end if
      !
      call get_nuclear_coordinates
      if (present(ion_charge)) ion_charge = total_charge
      !
      !  Copy density elsewhere, and allocate scratch field for the bare nuclei
      !  potential if it was not supplied on input.
      !
      if (size(f_table)<2) then
        write (out,"('Not enough scratch fields to construct eikonal potntial: ',i3)") size(f_table)
        stop 'eikonal_tools%eikonal_build_potential'
      end if
      f_rho       = f_table(1)
      f_bare_temp = f_table(2)
      if (present(f_bare)) f_bare_temp = f_bare
      call FieldCopy(src=f_core_pot,dst=f_rho)
      !
      call build_hartree_potential(f_rho=f_rho,f_bare=f_bare_temp,f_core_pot=f_core_pot)
      !
      if (present(v_xc)) then
        if (v_xc/=' ' .and. v_xc/='none') then
          f_vxc_temp = f_table(2)
          call build_exchange_correlation(v_xc,f_rho,f_vxc_temp)
          call FieldAXPY(alpha=(1._rk,0._rk),src=f_vxc_temp,dst=f_core_pot)
        end if
      end if
      !
      !  If caller asked for the multipole moments, return them. After this point,
      !  we switch to trace-less moments - and the long-range potential needs 
      !  straight Cartesian moments.
      !
      if (present(multipoles)) multipoles = multipoles_total
      if (present(electronic_multipoles)) then
        electronic_multipoles(0)   = core_electrons
        electronic_multipoles(1:9) = multipoles_electronic
      end if
      !
      !  Convert total multipole moment: Cartesian -> Spherical
      !
      multipoles_total(4:6) = 3*multipoles_total(4:6) - sum(multipoles_total(4:6))
      write (out,"(' Traceless quadrupole diagonal [electron-Bohr^2]: ',3(1x,f15.7))") &
             multipoles_total(4:6)
      call FLsetMultipolesPhase(mult=multipoles_total)
      call TimerStop('EikonalBuildPotential')
    end subroutine eikonal_build_potential
    !
    !  Constructs transition potentials for coulomb transition between
    !  ionic states caused by an external electron
    !        v(r) = \int \rho(r',r') \frac{1}{|r-r'|} d r'
    !  where
    !        \rho(r',r) = \sum \phi^L_i(r') n_i \phi^R_i(r)
    !  
    !  This gives the matrix elements <I_L|(1/|r-r'|)|I_R> where |I_L> and |I_R> are two ionic states.
    !
    !  Description from Serguei Patchkovskii:
    !  ---------------------------------------------------------------------------  
    !  The RDMs are factored as:
    !  
    !   \rho(r',r) = \sum \phi^L_i(r') n_i \phi^R_i(r)
    ! 
    !  where n_i are singular values of the 1-RDM, and \phi^L and \phi^R are
    !  respectively left and right singular vectors. The singular values are
    !  in the data section $RDMPCE of the each .dat file. Pairs of singular
    !  vectors (\phi^L_1, \phi^R_1, \phi^L_2, \phi^R_2, ...) are in the $VECRDM
    !  group.
    !  
    !  To get the "interaction" potential:
    !    v(r) = \int \rho(r',r') \frac{1}{|r-r'|} d r'
    ! 
    !  you'll need to slightly modify function "eikonal_build_potential" in
    !  eikonal_tools.f90. The only change which should be needed is to replace
    !  the call to build_electron_density with an alternative version, which
    !  instead of FieldRhoAccumulate (scaled square modulus of an orbital) will 
    !  calculate a scaled product of the left and right eigenvectors.
    !  
    !  Oh, you'll also have to change the sanity check on the total density
    !  (it should integrate to zero), and replace call to build_hartree_potential
    !  (which also does the nuclear part) with (a possibly scaled) call to
    !  fock_hartree_potential - see fock_exchange_potential in fock.f90 for
    !  something very similar.
    !  ---------------------------------------------------------------------------  
    !
    !  Added Feb 2, 2010 (Michael Spanner)
    !
    subroutine eikonal_build_transition_potential(rdm_file,rdm_sv,f_table,f_trans_pot,rot)
      character(len=*), intent(in)     :: rdm_file      ! GAMESS-format MO coefficients for the reduced density matrix
      real(rk), intent(in)             :: rdm_sv(:)     ! singular values of the 1-RDM decomposition
      integer(ik), intent(in)          :: f_table(:)    ! List of fields which can be used for scratch
      integer(ik), intent(in)          :: f_trans_pot   ! Field where transition potential will be stored
      real(ark), intent(in), optional  :: rot(:,:)      ! Rotation matrix
      !
      integer(ik) :: f_rho
      real(ark)   :: rotmat(3,3)
      !
      call TimerStart('EikonalBuildTransitionPotential')
      !
      !  Temporarily, use f_trans_pot to store the total electron density
      !
      call MathSetUnitMatrix(rotmat)
      if (present(rot)) rotmat = rot
      call build_electron_transition_density(f_trans_pot,f_table,rdm_file,rdm_sv,rotmat,transition=.true.)
      !
      !  Copy density elsewhere, and allocate scratch field for the bare nuclei
      !  potential if it was not supplied on input.
      !
      if (size(f_table)<1) then
        write (out,"('Not enough scratch fields to construct eikonal transition potential: ',i3)") size(f_table)
        stop 'eikonal_tools%eikonal_build_transition_potential'
      end if
      f_rho = f_table(1)
      call FieldCopy(src=f_trans_pot,dst=f_rho)
      !
      call fock_hartree_potential(v=f_trans_pot,rho=f_rho)
      !
      call TimerStop('EikonalBuildTransitionPotential')
    end subroutine eikonal_build_transition_potential
    !
    !  Construct an eikonal function for the specified momentum and scattering potential
    !
    subroutine eikonal_build_function(dst,kvec,potential,prefactor,scratch,norm,phase_only)
      integer(ik), intent(in)                :: dst         ! Output: Field filled with eikonal function
      real(rk), intent(in)                   :: kvec(3)     ! Input: Cartesian momentum vector 
      integer(ik), intent(in)                :: potential   ! Input: Field containing scattering potential
      logical, intent(in), optional          :: prefactor   ! Input: Include the prefactor?
      integer(ik), intent(in), optional      :: scratch     ! Input: Scratch field. Must be provided if
                                                            !        prefactor = .true.
      character(len=*), intent(in), optional :: norm        ! Input: Normalization of the eikonal function
                                                            !        'natural'  - normalize flux to abs(kvec)
                                                            !        'momentum' - normalize to delta-function in momentum
                                                            !        'energy'   - normalize to delta-function in energy
                                                            ! 'natural' is the default.
      logical, intent(in), optional          :: phase_only  ! Input: If .true., only evaluate the eikonal-Volkov
                                                            !        phase correction. phase_only==.true. implies
                                                            !        prefactor=.false. and norm='natural'
      !
      logical     :: do_prefactor, do_exponentiate
      complex(rk) :: factor ! Normalization factor
      !
      call TimerStart('EikonalBuildFunction')
      do_exponentiate = .true.
      if (present(phase_only)) do_exponentiate = .not.phase_only
      do_prefactor = .false.
      if (present(prefactor)) do_prefactor = prefactor
      !
      !  Minimal sanity checking
      !
      if (dst<=0 .or. potential<=0 .or. dst==potential) then
        stop 'eikonal_tools%eikonal_build_function - dst and/or potentials are bad'
      end if
      if (do_prefactor) then
        if (present(scratch)) then
          if (scratch<=0 .or. scratch==dst .or. scratch==potential) then
            stop 'eikonal_tools%eikonal_build_function - scratch is bad'
          end if
        else
          stop 'eikonal_tools%eikonal_build_function - scratch is not supplied'
        end if
      end if
      !
      !  Fill output array with analytical solution for quadrupole potential,
      !  then integrate from the appropriate boundary.
      !
      call FLsetPlanewave(k=kvec)
      call FieldInit(dst=dst,func=FLexactphase)
      call FieldBuildPhase(pot=potential,phase=dst,kvec=kvec)
      !
      if (do_exponentiate) then
        !
        !  Eikonal prefactor, will internally construct gradient contributions
        !
        if (do_prefactor) then
          call FieldBuildPref(phase=dst,pref=scratch,kvec=kvec)
        end if
        !
        !  Construct eikonal wavefunction from the eikonal correction factor
        !
        call FieldProcess(dst=dst,func=FLplanewavewithphase)
        !
        if (do_prefactor) then
          call FieldMul(src=scratch,dst=dst)
        end if
        !
        if (present(norm)) then
          select case (norm)
            case default
              write (out,"('eikonal_build_function: normalization ',a,' is not recognized')") &
                     trim(norm)
              stop 'eikonal_tools%eikonal_build_function - bad normalization'
            case ('natural')
              ! Do nothing
            case ('momentum')
              factor = twopi**(-1.5_rk)
              call FieldScale(dst=dst,con=factor)
            case ('energy')
              factor = twopi**(-1.5_rk) * sum(kvec**2)**(-0.25_rk)
              call FieldScale(dst=dst,con=factor)
          end select
        end if ! present(norm)
      end if ! do_exponentiate
      !
      call TimerStop('EikonalBuildFunction')
    end subroutine eikonal_build_function
    !
  end module eikonal_tools
