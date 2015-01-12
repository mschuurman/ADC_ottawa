!
!  Multigrid test - calculation of matrix elements for the eikonal-Volkov
!                   continuum solutions. No adiabatic approximation is
!                   assumed. Adiabatic solutions are easier, and are handled
!                   in amplitudes.f90
!
  module amplitudes_eva
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use eikonal_tools
    use eikonal_tools_eva
    implicit none
    private
    public run_amplitudes
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_mos      = 4                   ! Max number of MOs we support: Dyson + 3x exchange
    integer(ik), parameter :: max_fields   = 3+max_mos           ! 4 = 3x RDF components + planewave + phase correction
    integer(ik), parameter :: max_naturals = 200                 ! Max number of natural orbitals allowed
    integer(ik), parameter :: unit_dump    = 34                  ! A more or less random unit
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_points(3)     = (/ 100, 100, 100 /)   ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    character(len=100) :: dyson_file      = 'dyson.dat'           ! Name of the file containing MO coefficients
                                                                  ! for the Dyson and exchange orbitals
    character(len=100) :: natural_file    = 'natural.dat'         ! Name of the file containing MO coefficients
                                                                  ! for the natural orbitals of the cation
                                                                  ! (used for calculating Hartree potential)
    character(len=100) :: potential_file  = ' '                   ! File for the potential dump - can be HUGE
    integer(ik)        :: spin_parts(4)   = (/ 2, 4, 6, 8 /)      ! Orbital components to load from "dyson_file"
    integer(ik)        :: natural_count   = 0                     ! Number of natural orbitals in the density matrix
    real(rk)           :: natural_occ(max_naturals)               ! Natural orbital occupation numbers
    real(rk)           :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                  ! and exchange potential(s)
    real(rk)           :: sor_rate        = 1.95_rk               ! Over-relaxation rate in solving Poisson equations
    logical            :: plot_phases     = .false.               ! Plot numerical phases
    integer(ik)        :: plot_each       = 10                    ! Reduce plotting frequency to each (plot_each)-th step
    logical            :: eikonal_pref    = .true.                ! Include eikonal pre-factor. NOT IMPLEMENTED
    logical            :: node_factors    = .false.               ! Include integrals necessary for calculating emission
                                                                  ! from orbitals with nodes
    character(len=100) :: v_xc            = ' '                   ! Exchange-correlation potential to use in construction of the
                                                                  ! eikonal functions. Possible choices are:
                                                                  ! ' ' or 'none' - no exchange potential
                                                                  ! 'Slater' - Dirac-Slater exchange
                                                                  ! 'SVWN'   - Slater exchange + VWN local correlation
    real(rk)           :: field_e         = 0.000_rk              ! Electric field magnitude
    real(rk)           :: electron_k(3)   = (/ 0., 1., 0. /)      ! Electron momentum
    real(rk)           :: max_timestep    = 0.5_rk                ! Maximum allowed time step in EVA trajectory integral
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_dyson                                 ! Dyson orbital
    integer(ik)        :: f_exchange(3)                           ! Exchange orbitals
    integer(ik)        :: f_rdf(3)                                ! Recombination dipole field
    integer(ik)        :: f_core_pot                              ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_eva_phase                             ! Current eikonal-Volkov phase
    integer(ik)        :: f_eva_phase2                            ! New eikonal-Volkov phase
    !
    real(rk)           :: dipole_dyson(3)                         ! Dipole moment (e-Bohr) of the Dyson orbital
    complex(rk)        :: multipoles(0:9)                         ! Multipole moments of the ion core
    type(EVAdataT)     :: trajectory                              ! Current state of trajectory 
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /ampdata_eva/ verbose, &
                       n_points, box_extent, &
                       dyson_file, spin_parts, &
                       eps_hartree, sor_rate, &
                       natural_count, natural_file, natural_occ, &
                       potential_file, plot_phases, plot_each, &
                       eikonal_pref, node_factors, v_xc, &
                       electron_k, max_timestep, &
                       field_e
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Prepare simulation box and allocate a sufficient number of data fields
    !
    subroutine initialize_grid
      real(rk)    :: box(2,3)
      integer(ik) :: field
      !
      call TimerStart('Grid initialization')
      !
      box(1,:) = -0.5_rk*box_extent
      box(2,:) =  0.5_rk*box_extent
      !
      write (out,"(//t5,'Simulation box size (Bohrs)'/)")
      write (out,"(t8,i4,' pts for X: ',2f14.6)") n_points(1), box(:,1)
      write (out,"(t8,i4,' pts for Y: ',2f14.6)") n_points(2), box(:,2)
      write (out,"(t8,i4,' pts for Z: ',2f14.6)") n_points(3), box(:,3)
      write (out,"()")
      !
      call MultiGridInit(max_grids=1,max_fields=max_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate_fields: do field=1,max_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      f_free = max_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Allocate fields for and load Dyson and exchange correction orbitals
    !
    subroutine load_gamess_mos
      real(rk)              :: norm
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      f_free = f_free - 4
      if (f_free<0) stop 'build_electron_density - out of fields'
      call FieldImport('GAMESS',dyson_file,f_table(f_free+1:f_free+4),spin_parts)
      f_dyson         = f_table(f_free+1)
      f_exchange(1:3) = f_table(f_free+2:f_free+4)
      write (out,"('Loaded orbital components ',4i3)") spin_parts
      norm = FieldNorm(f_dyson)**2
      write (out,"('<psid|psid> = ',f20.10)") norm
      !
      call gamess_report_nuclei(nnuc)
      write (out,"('Data file contained ',i5,' nuclei')") nnuc
      allocate (xyzq(4,nnuc))
      call gamess_report_nuclei(nnuc,xyzq)
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,nnuc
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
      deallocate (xyzq)
    end subroutine load_gamess_mos
    !
    !  Build RDF (Recombination dipole field) from the Dyson and exchange orbitals
    !
    subroutine build_rdf
      integer(ik) :: ic, scr, scr_rdf
      real(rk)    :: ef(3)
      !
      scr_rdf = f_table(f_free) ; f_free = f_free - 1
      scr     = f_table(f_free) ; f_free = f_free - 1
      if (f_free<0) stop 'build_rdf - not enough fields!'
      !
      build_component: do ic=1,3
        ef = 0 ; ef(ic) = 1._rk
        call FLsetField(ef)
        call FieldCopy(src=f_exchange(ic),dst=scr_rdf)
        call FieldInit(dst=scr,func=FLelectricField)
        call FieldMulAdd(src_a=scr,src_b=f_dyson,dst=scr_rdf)
        call FieldCopy(src=scr_rdf,dst=f_exchange(ic)) ! Exchange component is no longer needed
      end do build_component
      !
      !  Copy handles from f_exchange to f_rdf, and kill values in f_exchange
      !  to prevent reuse
      !
      f_rdf      = f_exchange
      f_exchange = -1
      f_free     = f_free + 2
    end subroutine build_rdf
    !
    !  Fill input parameters for the EVA persistent state structure
    !
    subroutine fill_eva_parameters
      trajectory%s_fK         = electron_k                   ! Electron momentum
      trajectory%s_fF         = field_e                      ! Electric field magnitude
      trajectory%s_fw         = 0.057_rk                     ! Laser frequency
      trajectory%s_fAn        = (/ 1.0_rk, 0.0_rk, 0.0_rk /) ! Vector-potential direction
      trajectory%s_fMinT      = 0.0_rk*pi/trajectory%s_fw    ! Initial time
      trajectory%s_fMaxT      = 10.0_rk*pi/trajectory%s_fw   ! Final time
      trajectory%max_timestep = max_timestep
      trajectory%verbose      = verbose
    end subroutine fill_eva_parameters
    !
    !  Very simple visualization routine
    !
    subroutine visualize_wavefunction(text,psi)
      character(len=*), intent(in) :: text
      integer(ik), intent(in)      :: psi
      !
      call FieldVisualize(slot=0,src=psi,text=trim(text))
    end subroutine visualize_wavefunction
    !
    subroutine dump_potential
      call TimerStart('Dump core potential')
      open (unit=unit_dump,form='formatted',status='replace',file=potential_file)
      call FieldDump(unit_dump,src=f_core_pot)
      close (unit=unit_dump)
      call TimerStop('Dump core potential')
    end subroutine dump_potential
    !
    !  Calculate dipole amplitudes for a continuum wavefunction in
    !
    subroutine generate_amplitudes(tag,label,wf)
      character(len=1), intent(in) :: tag    ! Distinctive tag to use in the printout
      character(len=*), intent(in) :: label  ! Descriptive label for the continnum wf parameters
      integer(ik), intent(in)      :: wf     ! Field containing the continuum wavefunction
      !
      integer(ik) :: ic
      complex(rk) :: amp(3), r_amp(3,3)             ! Amplitudes and nodal corrections
      complex(rk) :: amp_overlap, r_amp_overlap(3)  ! Overlap with the Dyson orbital
      complex(rk) :: amp_ort(3), r_amp_ort(3,3)     ! Orthogonalized amplitudes
      !
      character(len=1), parameter :: xyz(3) = (/ 'X', 'Y', 'Z' /)
      !
      !  Exact scattering solutions would have zero overlap, but ...
      !
      amp_overlap = FieldConjgIntegrate(left=f_dyson,right=wf)
      amplitude_components: do ic=1,3
        amp(ic)     = FieldConjgIntegrate(left=f_rdf(ic),right=wf)
        amp_ort(ic) = amp(ic) - dipole_dyson(ic)*amp_overlap
      end do amplitude_components
      !
      write (out,"(1x,a1,'  ',1x,a44,1x,3(1x,f10.5,1x,f10.5,1x),2(f10.5,1x))") &
             tag, label, amp, amp_overlap
      write (out,"(1x,a1,'O ',1x,a44,1x,3(1x,f10.5,1x,f10.5,1x))") &
             tag, label, amp_ort
      !
      !  Now the nodal corrections
      !
      if (node_factors) then
        call FieldNorm2Multipoles(left=f_dyson,right=wf,mult=r_amp_overlap)
        amplitude_components_nodes: do ic=1,3
          call FieldNorm2Multipoles(left=f_rdf(ic),right=wf,mult=r_amp(:,ic))
          r_amp_ort(:,ic) = r_amp(:,ic) - r_amp_overlap(:)*dipole_dyson(ic)
        end do amplitude_components_nodes
        !
        !  Printing is a separate loop: we need to transpose the result.
        !
        amplitude_nodes_print: do ic=1,3
          write (out,"(1x,a1,'N',a1,1x,a44,1x,3(1x,f10.5,1x,f10.5,1x),2(f10.5,1x))") &
                 tag, xyz(ic), label, r_amp(ic,:), r_amp_overlap(ic)
          write (out,"(1x,a1,a1,'O',1x,a44,1x,3(1x,f10.5,1x,f10.5,1x))") &
                 tag, xyz(ic), label, r_amp_ort(ic,:)
        end do amplitude_nodes_print
      end if
    end subroutine generate_amplitudes
    !
    !  Swap field numbers
    !
    subroutine swap(i1,i2)
      integer(ik), intent(inout) :: i1, i2
      integer(ik)                :: t
      !
      t = i1 ; i1 = i2 ; i2 = t
    end subroutine swap
    !
    !  Problem driver
    !
    subroutine run_amplitudes
      integer(ik)        :: info
      character(len=200) :: comment
      !
      call TimerStart('Amplitudes EVA')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=ampdata_eva,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=ampdata_eva)
      write (out,"()")
      !
      call initialize_grid
      !
      call fill_eva_parameters
      call EVA_Initialize(trajectory=trajectory)
      !
      !  Do we have natural orbitals to construct the density and the Hartree potential?
      !
      if (natural_count>max_naturals) then
        write (out,"('Too many natural orbitals. Increase ''max_naturals''" &
                   // " to at least ',i0,' and recompile')") natural_count
        stop 'amplitudes - too many natural orbitals'
      end if
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      if (f_free<1) stop 'amplitudes - not enough fields for core potential'
      f_core_pot = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_potential(natural_file,natural_occ(:natural_count), &
                                   f_table(:f_free),f_core_pot,v_xc=v_xc,multipoles=multipoles)
      if (potential_file/=' ') then
        call dump_potential ! This call is likely to produce a HUGE file
      end if
      !
      !  Orbitals
      !
      call load_gamess_mos
      !
      !  Note that build_rdf will reuse fields originally allocated for f_exchange
      !
      call build_rdf
      !
      !  Multipole moments for the Dyson orbital - we need the dipole
      !
      call FieldNormMultipoles(f_dyson,dipole_dyson)
      write (out,"('Dipole moment of the Dyson orbital [electron-Bohr]: ',3f12.6)") dipole_dyson
      !
      !  Boundary conditions: we start from the adiabatic eikonal-Volkov solution
      !
      if (f_free<2) stop 'amplitudes - not enough fields for the initial w.f.'
      f_eva_phase  = f_table(f_free) ; f_free = f_free - 1
      f_eva_phase2 = f_table(f_free) ; f_free = f_free - 1
      !
      !  eikonal_build_function() will internally call FLsetPlanewave() and 
      !  FLsetMultipolesPhase(), so that FLexactphase() will yield the long-
      !  range adiabatic eikonal solution constent with the T=0 boundary
      !  conditions.
      !
      !  Adiabatic eikonal solutions are for a scattering solution with the
      !  flat _leaving_ front. The boundary condition we need presently is
      !  for the flat _incoming_ front. The easiest solution seems to be
      !  to change the direction of the electron movement in the adiabatic
      !  solution, then change the sign of the adiabatic solution.
      !
      call eikonal_build_function(dst=f_eva_phase,kvec=-electron_k,potential=f_core_pot,phase_only=.true.)
      call FieldScale(dst=f_eva_phase,con=(-1._rk,0._rk))
      call FieldZero(dst=f_eva_phase2)
      !
      if (plot_phases) then
        call visualize_wavefunction('Adiabatic eikonal boundary conditions',f_eva_phase)
      end if
      !
      if (verbose>=0) then
        write (out,"('Constructed field-free solution at the initial time')")
      end if
      !
      !  Set up multipole moments for long-range expansion of the ion core potential
      !  The width is entirely arbitrary here. The sign change is needed to make the
      !  potential attractive.
      !
      call FLsetLREmultipoles(-multipoles,1.0_rk)
      !
      write (out,"(/'Done with preliminaries'/)")
      call TimerReport
      !
      !  Ready to calculate eikonal-Volkov wavefunctions
      !
      eva_trajectory_loop: do while(trajectory%s_nt<trajectory%s_nNumT)
        write (out,"('# ',i6)") trajectory%s_nt
        !
        !  Warning: EVA_Advance depends on correct initialization of FLsetLREmultipoles,
        !           FLsetMultipolesPhase, and FLsetPlanewave. This state is (currently)
        !           not maintained by EVA_Advance.
        !
        call EVA_Advance(trajectory,potential=f_core_pot,src=f_eva_phase,dst=f_eva_phase2)
        call swap(f_eva_phase,f_eva_phase2)
        if (plot_phases .and. mod(trajectory%s_nt,plot_each)==0) then
          write (comment,"('EVA phase for time step ',i5,' time ',f12.3)") &
                 trajectory%s_nt, trajectory%s_pft(trajectory%s_nt)
          call visualize_wavefunction(trim(comment),f_eva_phase)
        end if
        !
!         call generate_amplitudes(tag='R',label=label,wf=f_plane)
      end do eva_trajectory_loop
      !
      call EVA_Destroy(trajectory)
      !
      call TimerStop('Amplitudes EVA')
      call TimerReport
    end subroutine run_amplitudes

  end module amplitudes_eva
!
  subroutine driver
    use amplitudes_eva

    call run_amplitudes
  end subroutine driver
