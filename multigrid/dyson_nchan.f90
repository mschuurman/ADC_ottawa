!
!  Propagation of many-electron wavefunction in the basis of ionic states.
!  The theory is described in: M. Spanner and S. Patchkovskii, PRA 80, 063411 (2009).
!
!  This module implements the general multi-channel case of the paper - see section 
!  IID and the appendix.
!
!  Input preparation for this program is rather non-trivial.
!
  module dyson_nchan
    use dyson_tools
    use lapack

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !       Keep in mind that some of the relevant constants are in dyson_tools.f90
    !
    integer(ik), parameter :: max_channels = 20                  ! Max number of channels allowed (as of March 2010, 20 channels is _very_ large)
    integer(ik), parameter :: max_rdm      = 200                 ! Max number of RDM orbitals allowed in a single rdm file
    integer(ik), parameter :: max_rdm_fields = max_channels*(max_channels-1)/2   
                                                                 ! Max number of RDM fields needed for chosen value of max_channels
                                                                 ! (number of unique off-digonal matrix elements of ionic basis)
    integer(ik), parameter :: max_fields   = 13 * max_channels + max_rdm_fields + 1  
                                                                 ! max_fields = max_channels * [ dyson + 3*(exchange or radf) + 3*calT + core + 
                                                                 !      total_pot + psi + psi2 + Hpsi + phiIm ] + max_rdm_fields + 1*abs_pot
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: num_channels    = 2                     ! Number of ionic channels to use

    character(len=clen) :: dyson_files(max_channels) = ' '         ! Name of the files containing MO coefficients
                                                                   ! for the Dyson and exchange orbitals
    integer(ik)         :: spin_parts(4,max_channels)              ! Orbital components to load from "dyson_file"
    logical             :: naturals_are_rdms = .false.             ! Set to true to read natural orbitals in superdyson 1-RDM format
    character(len=clen) :: natural_files(max_channels) = ' '       ! Name of the files containing MO coefficients
                                                                   ! for the natural orbitals of the cation
                                                                   ! (used for calculating Hartree potential)
    character(len=clen) :: rdm_files( max_rdm_fields ) = ' '       ! Name of the files containing rdm orbital coefficients
    logical             :: use_interchannel_coupling = .true.      ! Flag to toggle inter-channel coupling. Can be further controlled by
                                                                   ! the two options below.
    logical             :: use_rdm_coupling          = .true.      ! Include inter-channel coupling terms due to transition potentials
    logical             :: use_dipole_coupling       = .true.      ! Include inter-channel coupling terms due to transition dipoles
    real(rk)            :: neutral_ener           =  0._rk         ! Neutral groundstate energy, note: zero is at ionic ground energy
    real(rk)            :: ion_ener(max_channels) =  0._rk         ! energy of current ionic state, note: zero is at ionic ground energy
                       
    real(rk)            :: neutral_dipole(3)          = 0._rk      ! IMPORTANT: The dipole moments include electron charge (-1)
    real(rk)            :: ion_dipole(3,max_channels) = 0._rk      ! They are otherwise in atomic units - NOT in Debye
    real(rk)            :: transition_dipole(3,max_rdm_fields) = 0._rk
                                                                   ! all the unique offdiagonal dipole moments of ion: -|e| <I1|sum_r|I2>, 
                                                                   !                                                   -|e| <I1|sum_r|I3>,
                                                                   !                                                   ...
                                                                   !                                                   -|e| <I_(n-1)|sum_r|In>
                                                                   ! That is, the transition dipoles include the negative electron charge
                                                                   ! in them (eq. A4). This is consistent with the definition of the permanent
                                                                   ! dipoles.
    character(len=clen) :: ions_preview = ' '                      ! Non-empty file name here will trigger a preview of ion dynamics
    real(rk)            :: Tkill_bound_state = 1e50_rk             ! HACK: Delete the ground-state population at this time, to produce cleaner
                                                                   ! HACK: photoelectron spectra
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: total_rdm_fields                        ! Total number of actual RDM fields requires for the input value of num_channels
    integer(ik)        :: f_dyson(max_channels)                   ! Dyson orbital. Eventually renormalized to the "source"
                                                                  ! form, so beware.
    integer(ik)        :: f_exchange(3,max_channels)              ! Exchange orbitals (shares storage with f_radf)
    integer(ik)        :: f_radf(3,max_channels)                  ! Recombination dipole field or acceleration dipole field,
                                                                  ! depending on the value of harmonic_operator. Shares storage
                                                                  ! with f_exchange, so only one would be in use at any given time
    integer(ik)        :: f_calT(3,max_channels)                  ! "transfer" state (1='main', 2='F part', 3= full)
                                                                  ! See eq. 36. The f_calT(1) part does not depend on the
                                                                  ! electric field. The f_calT(2) part is linear in the magnitude
                                                                  ! of the field. This partitioning is for the special case of
                                                                  ! the linear polarization of the electric field.
                                                                  ! IMPORTANT: Note that the f_calT(2) part is actually linear in (-CurrentElaser),
                                                                  ! IMPORTANT: which may cause some confusion if you are not careful!
    integer(ik)        :: f_core_pot(max_channels)                ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_trans_pot(max_rdm_fields)             ! Coulomb transition density potentials for ionic coupling
    integer(ik)        :: f_total_pot(max_channels)               ! total potential (laser and core and ion energy)
    integer(ik)        :: f_psi(max_channels)                     ! present wavefunction 
    integer(ik)        :: f_psi2(max_channels)                    ! and 'next' wavefunction
    integer(ik)        :: f_Hpsi(max_channels)                    !  = H |psi>
    integer(ik)        :: f_phiIm(max_channels)                   !  = |phi^I_m> = ampD |dyson> + |chi^I_m>

    integer(ik), allocatable  &
                       :: f_flux_momentum(:,:)                    ! Fields used to accumulate the photoelectron spectrum in momentum space

    integer(ik)        :: trans_id(max_channels,max_channels)     ! Table of indices relating the transition matrix element <n|0|m> to the actual storage index
                                                                  ! Used for transition_dipole and trans_pot
    integer(ik)        :: natural_count   = 0                     ! Number of natural orbitals in the density matrix
    real(rk)           :: natural_occ(max_naturals)               ! Natural orbital occupation numbers
                                                                  ! natural_count and natural_occ are used only for loading density matrices, and are
                                                                  ! discarded/overwritten afterwards
    integer(ik)        :: rdm_count   = 0                         ! Number of 1-rdm singular orbital pairs in the transition density
    real(rk)           :: rdm_sv(max_rdm)                         ! Transition 1-RDM decomposition singular values

    real(rk)           :: dipole_dyson(3,max_channels)            ! Dipole moments (e-Bohr) of the "source" (renormalized Dyson) orbitals 
                                                                  ! <dyson|r|dyson>. Note that the sign of electron charge is NOT included.
!   real(rk)           :: acceleration_dyson(3,max_channels)      ! Acceleration of the Dyson orbitals - only defined 
!                                                                 ! if harmonic_operator='acceleration'. NOT IMPLEMENTED
    complex(rk)        :: multipoles_ion(0:9,max_channels)        ! Multipole moments of the ion core (electrons only)
                                                                  ! Charge of an electron is taken as +1 for this field
                                                                  ! multipoles_ion(0,nc) is the total number of electrons in the cation nc
    real(rk)           :: N_electrons                             ! Number of electrons in the neutral species
    real(rk)           :: eta(max_channels)                       ! amplitudes of normalized dyson state contained in the neutral
                                                                  !  -->  eta = (<normalized dyson|<I_m|) |neutral>
    real(rk)           :: calN                                    ! Normalization factor for |\tilde N> state

    complex(rk)        :: d_Hm_d(max_channels)                    ! = <dyson|H_m|dyson>
    complex(rk)        :: HtildeN_main                            ! = {\cal H}^{\tilde N}  --- 'main' part, see eq. 38
    complex(rk)        :: HtildeN_F                               ! = {\cal H}^{\tilde N}  --- 'F(t)' part, see eq. 38 (special case of
                                                                  !   the linear polarization)

    complex(rk)        :: AmpD(max_channels)  = 0._rk             ! Amplitudes of the dyson states
    complex(rk)        :: AmpD2(max_channels) = 0._rk             ! Amplitudes of the dyson states at 'next' time step
                                                                  ! In the paper, AmpD/AmpD2 is given by a_m(t)
                                                                  
    complex(rk)        :: AmpB           = 0._rk                  ! Amplitude of the |\tilde N> state
    complex(rk)        :: AmpB2          = 0._rk                  ! Amplitude of the |\tilde N> state at 'next' time step
                                                                  
    real(rk)           :: AbsorbedNorm (max_channels)   = 0._rk   ! Amount of absorbed population (norm)
    real(rk)           :: recombination(3)              = 0._rk   ! Total recombination dipole (harmonic_operator='dipole')
                                                                  ! or dipole acceleration (harmonic_operator='acceleration') at the
                                                                  ! last time step. Accelertion form is currently not implemented
    real(rk)           :: rec_channel(3,max_channels,max_channels) = 0._rk   
                                                                  ! Per-channel contributions to recombination dipole
                                                                  ! The second index is the "to" channel; 
                                                                  ! the third index is the "from" channel
    real(rk)           :: ChiNorm(max_channels) = 0._rk           ! Norm of the excited/continuum parts of the wavefunction

    real(rk)           :: VecsD(max_channels,max_channels)        ! Amplitudes of dressed-state vectors
    real(rk)           :: EnerD(max_channels)                     ! Energies of dressed-state vectors
    logical            :: bound_state_killed = .false.            ! HACK: Set to .true. once the ground-state has been deleted

    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /dyson_nchan_data/ verbose, &
                       n_points, box_extent, &
                       num_channels, &
                       dyson_files, spin_parts, &
                       nspin, euler_angles, &
                       eps_hartree, sor_rate, &
                       naturals_are_rdms, natural_files, rdm_files, &
                       use_interchannel_coupling, use_rdm_coupling, use_dipole_coupling, &
                       v_xc, &
                       ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
                       neutral_ener, ion_ener, &
                       tMin, tMax, dt, &
                       ElectricFieldShape, &
                       Elaser, Edir, Wlaser, Plaser, TlaserOn, TlaserOff, RampTime, FlatTime, & 
                       FWHM, Tlaser, &
                       neutral_dipole, ion_dipole, transition_dipole, abs_bound_kmin, abs_bound_type, &
                       harmonic_operator, &
                       output_prefix, detailed_output, report_each, &
                       potential_file, odx_plot_each, odx_plot_momentum, efield_preview, &
                       ions_preview, &
                       checkpoint_each, checkpoint_out, checkpoint_in, &
                       oversample_orbitals, &
                       potentials_in, potentials_out, &
                       continue_flux, flux_guard, flux_sensors, &
                       flux_samples, flux_max_p, flux_restart,  &
                       fft_thread_threshold, Tkill_bound_state
    !
    !  ==== All parameters needed for checkpoint/restart are in the namelist below ====
    !
    namelist /dyson_nchan_checkpoint/ ampD, ampD2, ampB, ampB2, AbsorbedNorm, &
                                      recombination, ChiNorm, time, tstep, rec_channel, &
                                      CurrentAField, bound_state_killed
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Rotate static input data - at the moment, these are only the dipoles
    !
    subroutine rotate_dipoles
      logical             :: go
      character(len=clen) :: buf
      integer(ik)         :: i, j, ij
      ! 
      go = verbose>=0 .and. .not.MathIsUnitMatrix(rotmat)
      !
      if (go) then
        write (out,"()")
        write (out,"(1x,12x,     a39 ,6x,     a39   )") '   Rotated dipoles [e-Bohr]    ', '  Original dipoles [e-Bohr]    '
        write (out,"(1x,12x,     a39 ,6x,     a39   )") '-----------------------------  ', '-----------------------------  '
        write (out,"(1x,a12,3(1x,a12),6x,3(1x,a12  ))") '       ', '   X   ', '   Y   ', '   Z   ', '   X   ', '   Y   ', '   Z   '
        write (out,"(1x,a12,3(1x,a12),6x,3(1x,a12  ))") '       ', '-------', '-------', '-------', '-------', '-------', '-------'
      end if
      !
      call do_dipole(' neutral ',neutral_dipole)
      permanent_dipoles: do i=1,num_channels
        write (buf,"('ion ',i3)") i
        call do_dipole(buf,ion_dipole(:,i))
      end do permanent_dipoles
      ij = 0
      transition_dipoles_l: do i=1,num_channels-1
        transition_dipoles_r: do j=i+1,num_channels
          ij = ij + 1
          write (buf,"('tr ',i3,'<>',i3)") i, j
          call do_dipole(buf,transition_dipole(:,ij))
        end do transition_dipoles_r
      end do transition_dipoles_l
      !
      if (go) write (out,"()")
      !
      contains
      subroutine do_dipole(name,dip)
        character(len=*), intent(in) :: name
        real(rk), intent(inout)      :: dip(3)
        !
        real(rk) :: dip_save(3)
        !
        dip_save = dip
        dip      = matmul(rotmat,dip)
        if (go) write (out,"(1x,a12,3(1x,f12.6),6x,3(1x,f12.6))") name, dip, dip_save
      end subroutine do_dipole
    end subroutine rotate_dipoles
    !
    !  Prepare ion-specific potentials: core+Hartree and transition potentials
    !
    subroutine load_ion_potentials
      integer(ik) :: nc
      !
      call TimerStart('Load ion potentials')
      !     
      !  Do we have natural orbitals to construct the density and the Hartree potential?
      !
      load_cores: do nc=1,num_channels
        write (out,"(/'Loading ',a,' naturals_are_rdms = ',l3)") trim(natural_files(nc)), naturals_are_rdms
        if (naturals_are_rdms) then
          call gamess_load_rdmsv(trim(natural_files(nc)),natural_occ,natural_count)
        else
          call gamess_load_natocc(trim(natural_files(nc)),natural_occ,natural_count)
        end if
        if (f_free<1) stop 'amplitudes - not enough fields for core potentials'
        call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
        f_core_pot(nc) = f_table(f_free) ; f_free = f_free - 1
        call eikonal_build_potential(trim(natural_files(nc)),natural_occ(:natural_count), &
                                     f_table(:f_free),f_core_pot(nc),v_xc=v_xc, &
                                     electronic_multipoles=multipoles_ion(:,nc),rot=rotmat,use_rdms=naturals_are_rdms)
      end do load_cores
      !
      !  Make sure all ion cores have the same electron count - if the user made a simple
      !  mistake here, we'd rather know now than in three weeks when the run is done.
      !
      scan_electron_counts: do nc=2,num_channels
        if (abs(multipoles_ion(0,nc-1)-multipoles_ion(0,nc))<1e-4_rk) cycle
        write (out,"('Ion core ',i3,' has ',g21.14,' electrons, while core ',i3,' has ',g21.14,'. Oops.')") &
               nc-1, real(multipoles_ion(0,nc-1),kind=rk), nc, real(multipoles_ion(0,nc),kind=rk)
        stop 'dyson_nchan%load_ion_potentials - inconsistent electron counts'
      end do scan_electron_counts
      !
      !  Load transition potential - this is something we didn't have in single-channel
      !
      if (use_interchannel_coupling .and. use_rdm_coupling) then
        call TimerStart('Load transition potentials')
        load_transition_potentials: do nc=1,total_rdm_fields
          write (out,"(/'Loading ',a)") trim(rdm_files(nc))
          call gamess_load_rdmsv(trim(rdm_files(nc)),rdm_sv,rdm_count)
          call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
          if (f_free<1) stop 'setup_fields - not enough fields for transition potentials'
          f_trans_pot(nc) = f_table(f_free) ; f_free = f_free - 1
          call eikonal_build_transition_potential(trim(rdm_files(nc)),rdm_sv(:rdm_count), &
                                                  f_table(:f_free),f_trans_pot(nc),rot=rotmat)
        end do load_transition_potentials
        call TimerStop('Load transition potentials')
      end if
      call TimerStop('Load ion potentials')
    end subroutine load_ion_potentials
    !
    !  Save ion potentials and associated data to disk for later re-use
    !
    subroutine checkpoint_ion_potentials
      integer(ik)         :: nc
      integer(ik)         :: seq
      character(len=clen) :: tmp_name
      !
      call TimerStart('Checkpoint ion potentials')
      !     
      open(unit_ckpt,form='formatted',status='replace',file=trim(potentials_out)//'_multipoles.dat')
      write(unit_ckpt,*) multipoles_ion(:,:num_channels)
      close(unit_ckpt)
      !
      seq = 1
      dump_cores: do nc=1,num_channels
        write (tmp_name,"(a,i5.5,'.bin')") trim(potentials_out), seq
        seq = seq + 1
        call FieldExport('binary',f_core_pot(nc),trim(tmp_name))
      end do dump_cores
      !
      if (use_interchannel_coupling .and. use_rdm_coupling) then
        dump_transitions: do nc=1,total_rdm_fields
          write (tmp_name,"(a,i5.5,'.bin')") trim(potentials_out), seq
          seq = seq + 1
          call FieldExport('binary',f_trans_pot(nc),trim(tmp_name))
        end do dump_transitions
      end if
      call TimerStop('Checkpoint ion potentials')
    end subroutine checkpoint_ion_potentials
    !
    !  Restart ion potentials from disk. No checking of any kind!
    !
    subroutine restart_ion_potentials
      integer(ik)         :: nc
      integer(ik)         :: seq
      character(len=clen) :: tmp_name
      !
      call TimerStart('Restart ion potentials')
      !
      open(unit_ckpt,form='formatted',action='read',position='rewind',file=trim(potentials_in)//'_multipoles.dat')
      read(unit_ckpt,*) multipoles_ion(:,:num_channels)
      close(unit_ckpt)
      !
      seq = 1
      fetch_cores: do nc=1,num_channels
        if (f_free<1) stop 'amplitudes - not enough fields for core potentials'
        f_core_pot(nc) = f_table(f_free) ; f_free = f_free - 1
        write (tmp_name,"(a,i5.5,'.bin')") trim(potentials_in), seq
        seq = seq + 1
        call FieldImport('binary',trim(tmp_name),(/f_core_pot(nc)/),(/1/))
      end do fetch_cores
      !
      if (use_interchannel_coupling .and. use_rdm_coupling) then
        fetch_transitions: do nc=1,total_rdm_fields
          if (f_free<1) stop 'setup_fields - not enough fields for transition potentials'
          f_trans_pot(nc) = f_table(f_free) ; f_free = f_free - 1
        write (tmp_name,"(a,i5.5,'.bin')") trim(potentials_in), seq
        seq = seq + 1
        call FieldImport('binary',trim(tmp_name),(/f_trans_pot(nc)/),(/1/))
        end do fetch_transitions
      end if
      call TimerStop('Restart ion potentials')
    end subroutine restart_ion_potentials
    !
    !  Allocate fields for and load Dyson and exchange correction orbitals
    !
    subroutine load_gamess_mos
      real(rk)              :: norm
      integer(ik)           :: nc
      !
      call TimerStart('Load Gamess MOs')
      !
      !  Load all dyson and "exchange correction" orbitals
      !
      load_orbitals: do nc=1,num_channels
        f_free = f_free - 4
        if (f_free<0) stop 'dyson_nchan%load_gamess_mos - out of fields'
        write (out,"('Loading ',a)") trim(dyson_files(nc))
        call FieldImport('GAMESS',trim(dyson_files(nc)),f_table(f_free+1:f_free+4),spin_parts(1:4,nc),rot=rotmat)
        f_dyson(nc)        = f_table(f_free+1)
        f_exchange(1:3,nc) = f_table(f_free+2:f_free+4)
        write (out,"('Loaded orbital components ',4i3)") spin_parts(:,nc)
        norm = FieldNorm(f_dyson(nc))**2
        write (out,"('<psid|psid> = ',f20.10)") norm
        !
        !  Rotate exchange correction vector fields
        !
        call FieldRotateVectorComponents(rotmat,f_exchange(1,nc),f_exchange(2,nc),f_exchange(3,nc))
        !
        !  As a very simple sanity check, report number of nuclei in each file; for the
        !  final structure, also print the coordinates.
        !
        if (nc< num_channels) call report_structure(2_ik,trim(dyson_files(nc)))
        if (nc==num_channels) call report_structure(1_ik,trim(dyson_files(nc)))
      end do load_orbitals
      call TimerStop('Load Gamess MOs')
    end subroutine load_gamess_mos
    !
    !  Prepare "source" orbitals - these are eqs. 14, 15, and 19.
    !
    subroutine prepare_dyson_fields_source
      integer(ik) :: nc
      real(rk)    :: calNtmp
      !
      call TimerStart('Dyson fields: Source')
      !
      !  Remove the sqrt(N) normalization factor from Dyson and "exchange" orbitals
      !
      N_electrons = real(multipoles_ion(0,1),kind=rk) + 1._rk
      !
      calNtmp = 0._rk
      !
      do nc=1,num_channels
        call FieldScale(f_dyson(nc),     cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
        call FieldScale(f_exchange(1,nc),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
        call FieldScale(f_exchange(2,nc),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
        call FieldScale(f_exchange(3,nc),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))

        eta(nc) = FieldNorm(f_dyson(nc))
        if (abs(eta(nc))<=1e-10_rk) then
          write (out,"('WARNING: Norm of Dyson orbital for channel ',i3,' is very small or zero: ',g16.8)") nc, eta(nc)
          write (out,"('WARNING: Skipping renormalization of the channel ',i3,' and resetting the norm to hard zero')") nc
          eta(nc) = 0._rk
        else
          call FieldScale(f_dyson(nc),cmplx(1._rk/eta(nc),0._rk,kind=rk))
        end if
 
        calNtmp = calNtmp + abs(eta(nc))**2
      end do

      calN = 1._rk / sqrt(1._rk-nspin*calNtmp)
      !
      write (out,"()")
      write (out,"(' Number of electrons (N_electrons)    = ',g14.7)") N_electrons
      write (out,"(' Spin factor (nspin)                  = ',i7   )") nspin
      write (out,"(' Neutral residue normalization (calN) = ',g14.7)") calN
      if (mod(nint(N_electrons)+nspin,2)/=0) then
        write (out,"(//'  WARNING: Number of electrons and number of spin components are not consistent!'//)")
      end if
      write (out,"()")
      !
      !  Multipole moments for the Dyson orbital - we need the dipole
      !  These are used in eq. A19 later on.
      !
      do nc=1,num_channels
        call FieldNormMultipoles(f_dyson(nc),dipole_dyson(:,nc))
        write (out,"('           ""Source"" orbital # ',i3,' normalization (eta)   = ',g14.7)") nc, eta(nc)
        write (out,"('Average coordinate of the ""source"" orbital # ',i3,' [Bohr] = ',3f12.6)") nc, dipole_dyson(:,nc)
        write (out,"(' Overlaps of renormalized dyson and ""craddle"" orbitals # ',i3,':')") nc
        write (out,"('    <D|C_x> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson(nc),f_exchange(1,nc))
        write (out,"('    <D|C_y> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson(nc),f_exchange(2,nc))
        write (out,"('    <D|C_z> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson(nc),f_exchange(3,nc))
        write (out,"()")
      end do
      call TimerStop('Dyson fields: Source')
    end subroutine prepare_dyson_fields_source
    !
    !  Build calT 'main' - see Eq. 36. (Note that the paper does not define calT for the
    !  multichannel case; the same-channel part is identical to eq. 36; the channel-coupling
    !  terms are obtained by inspection of A7/A19 and A14/A23, together with eq. 28
    !
    !  The best way of thinking of calT is probably as of the matrix element:
    !
    !    <I_m|H^F|{\tilde N}>, where I_m is the (n-1)-electron wavefunction of the
    !    ion core in the channel m; H^F is the total Hamiltonian of the system,
    !    including the laser field contribution; and {\tilde N} is the renormalized
    !    residue of the neutral wavefunction, left after projecting out all Dyson
    !    channels.
    !
    subroutine prepare_dyson_fields_calTmain
      integer(ik) :: nc, ic, ncc
      integer(ik) :: f_tmp
      complex(rk) :: w_tmp
      complex(rk) :: orb_e, orb_ekh
      !
      call TimerStart('Dyson fields: calT(main)')
      build_calT_loop_channels: do nc=1,num_channels
        f_free = f_free - 3
        if (f_free<0) stop 'dyson_nchan%prepare_dyson_fields - out of fields for calT storage'
        f_calT(1:3,nc) = f_table(f_free+1:f_free+3)
        !
        ! 'main' part (i.e. part with no laser field prefactor)
        !
        call QMHpsi(1._rk,f_core_pot(nc),f_dyson(nc),f_calT(1,nc))                     ! 
        orb_ekh = FieldConjgIntegrate(f_dyson(nc),f_calT(1,nc))
        call ecp_apply(src=f_dyson(nc),dst=f_calT(1,nc),ecp=ecp)                       ! ECP is a non-local operator; it can't be folded into f_core_pot
        orb_e   = FieldConjgIntegrate(f_dyson(nc),f_calT(1,nc))
        if (verbose>=0) then
          write (out,"('      Orbital eigenvalue (T+V_nuc+V_ee) for channel ',i3,' source orbital: ',f14.7,1x,f14.7)") nc, orb_ekh
          write (out,"('Orbital eigenvalue (T+V_nuc+V_ee+V_ecp) for channel ',i3,' source orbital: ',f14.7,1x,f14.7)") nc, orb_e
        end if
        call FieldAXPY(cmplx(ion_ener(nc),0._rk,kind=rk),f_dyson(nc),f_calT(1,nc))     !  Operator of Eq. 33 acting on the source orbital
        d_Hm_d(nc) = FieldConjgIntegrate(f_dyson(nc),f_calT(1,nc))                     ! <dyson|H_m|dyson> is used below for HtildeN_main (eq. 38)
        call FieldScale(f_calT(1,nc),cmplx(-1._rk,0._rk,kind=rk))
        call FieldAXPY(cmplx(neutral_ener,0._rk,kind=rk),f_dyson(nc),f_calT(1,nc))
        call FieldScale(f_calT(1,nc),cmplx(eta(nc)*calN,0._rk,kind=rk))                ! Finished construction of the field-independent part of eq. 36
        !
        !  'F' part (i.e. part _with_ laser field prefactor) - see Eq. 36. Note that the construction below depends on electric field
        !  being linearly polarized, even though the theory itself allows general polarization.
        !
        call FieldCopy(src=f_dyson(nc),dst=f_calT(2,nc))
        call FieldScale(f_calT(2,nc),cmplx(dot_product(Edir,ion_dipole(:,nc))*eta(nc)*calN,0._rk,kind=rk)) ! PS dipole change
        ! The equation above _appears_ to differ from the paper, but is not:
        ! the overall sign becomes - once the electric 
        ! magnitude is included in performDysonTimeStep.
        !
        do ic=1,3
          call FieldAXPY(cmplx(calN*Edir(ic),0._rk,kind=rk),f_exchange(ic,nc),f_calT(2,nc))
          !       Published eq. 36 has a minus sign above - however, this is OK, since the
          !       sign is recovered upon multiplication with electric field magnitude in 
          !       performDysonTimeStep.
        end do 
        !
        !  Inter-channel coupling.
        !
        if (use_interchannel_coupling) then
          if (f_free<1) stop 'dyson_nchan%prepare_dyson_fields_calTmain - not enough fields for interchannel coupling temp field'
          f_tmp = f_table(f_free) ; f_free = f_free - 1
          couple_ncc: do ncc=1,num_channels
            if (ncc==nc) cycle couple_ncc
            if (use_dipole_coupling) then
              !
              !  F part - see the last term in eq. A19. This is non-adiabatic transtion within the ion core.
              !           Note that eq. A19 has a minus sign for this term; there is no error, since the
              !           sign in compensated in performDysonTimeStep, analogously to the permanent dipoles.
              !
              w_tmp = cmplx(calN*eta(ncc)*dot_product(Edir,transition_dipole(:,trans_id(nc,ncc))),0._rk,kind=rk)
              call FieldAXPY(w_tmp,f_dyson(ncc),f_calT(2,nc))
            end if
            if (use_rdm_coupling) then
              !
              !  Main part - the next to last term in eq. A19. This is the inter-channel electron correlation.
              !
              call FieldCopy(src=f_dyson(ncc),dst=f_tmp)
              call FieldMul(src=f_trans_pot(trans_id(nc,ncc)),dst=f_tmp)
              w_tmp = cmplx(-calN*eta(ncc),0._rk,kind=rk)
              call FieldAXPY(w_tmp,f_tmp,f_calT(1,nc))
            end if
          end do couple_ncc
          f_tmp = -1; f_free = f_free + 1;
        end if
      end do build_calT_loop_channels
      call TimerStop('Dyson fields: calT(main)')
    end subroutine prepare_dyson_fields_calTmain
    !
    !  Build HtildeN, including spin factor - see Eq. 38. Again, the intra-channel part
    !                 is very similar to the uncoupled case; the remaining terms are
    !                 from eq. A21
    !
    subroutine prepare_dyson_fields_HtildeN
      integer(ik) :: ic, nc, ncc, f_tmp
      complex(rk) :: cx_tmp  ! just a temp complex value...
      !
      call TimerStart('Dyson fields: H(tildeN)')
      !
      !  The first term in A21
      !
      HtildeN_main = abs(calN)**2 * neutral_ener       ! A21, first part of A12
      HtildeN_F = -dot_product(Edir,neutral_dipole)    ! A21, second part of A12 - sign reversal compensated later
      !
      !  The last sum in A21 - terms from A7 and A9
      !
      build_HtildeN: do nc=1,num_channels
        !
        !  Build HtildeN_main
        !  d_Hm_d is from A7 (first term); the 2*neutral_ener is the first terms in A9 and its conjugate
        HtildeN_main = HtildeN_main + abs(calN)**2 * nspin*abs(eta(nc))**2*(d_Hm_d(nc) - 2._rk*neutral_ener) 
        !
        !  Build HtildeN_F
        !  The ion_dipole term is the last part of A7.
        !  The dipole_dyson term is from A7 (+) and twice from A9 (-), giving overall (-)
        HtildeN_F = HtildeN_F - nspin*abs(eta(nc))**2 * dot_product(Edir,ion_dipole(:,nc)) &
                              - nspin*abs(eta(nc))**2 * dot_product(Edir,dipole_dyson(:,nc))  ! PS dipole change
                              ! Together with the minus sign for CurrentElaser in performDysonTimeStep,
                              ! this gives correct sign overall
        !
        !  The loop below is the last term in A9
        ! 
        do ic=1,3
          cx_tmp = FieldConjgIntegrate(f_exchange(ic,nc),f_dyson(nc)) * eta(nc)
          ! Published eq. 38 has a plus sign here?!
          ! See comment immediately above
          HtildeN_F = HtildeN_F - nspin*( cx_tmp + conjg(cx_tmp) ) * Edir(ic)
        end do 
      end do build_HtildeN
      !
      !  inter-channel coupling for HtildeN
      !
      if (use_interchannel_coupling) then
        if (f_free<1) stop 'dyson_nchan%prepare_dyson_fields_HtildeN - not enough fields for interchannel coupling temp field'
        f_tmp = f_table(f_free) ; f_free = f_free - 1
        !
        !  This is the second term in eq. A21 - terms from eq. A8
        !  Note that HtildeN_main is accumulated with calN**2 included, while HtildeN_F is not.
        !  I know this is confusing.
        !
        do nc=1,num_channels
        do ncc=nc+1,num_channels
          if (use_dipole_coupling) then
            ! F part - last term in A8. Extra factor of 2 is because we are have both nc-ncc and ncc-nc channels
            HtildeN_F = HtildeN_F - nspin*2._rk*real( eta(nc)*eta(ncc) * dot_product(Edir,transition_dipole(:,trans_id(nc,ncc))) * &
                                            FieldConjgIntegrate(f_dyson(nc),f_dyson(ncc)) )
                                                ! Together with the minus sign for CurrentElaser in performDysonTimeStep,
                                                ! this gives correct sign overall
          end if
          if (use_rdm_coupling) then
            ! main part - first term in A8. Again, extra factor of two.
            call FieldCopy(src=f_dyson(ncc),dst=f_tmp)
            call FieldMul(src=f_trans_pot(trans_id(nc,ncc)),dst=f_tmp)
            HtildeN_main = HtildeN_main + nspin*2._rk*real( eta(nc)*eta(ncc)*FieldConjgIntegrate(f_dyson(nc),f_tmp) )*abs(calN)**2
          end if
        end do
        end do
        f_tmp = -1; f_free = f_free + 1;
      end if
      !
      HtildeN_F = HtildeN_F * abs(calN)**2  ! Overall normalization factor in eq. 38
      !
      call TimerStop('Dyson fields: H(tildeN)')
    end subroutine prepare_dyson_fields_HtildeN
    !
    !  Prepare "transfer" orbitals needed for Hamiltonian evaluation.
    !  This routine was getting unwieldingly large, so I had to split it.
    !
    subroutine prepare_dyson_fields
      call TimerStart('Dyson Fields')
      call prepare_dyson_fields_source
      !
      call prepare_dyson_fields_calTmain
      !
      call prepare_dyson_fields_HtildeN
      call TimerStop('Dyson Fields')
    end subroutine prepare_dyson_fields
    !
    !  Fills the trans_id table
    !
    subroutine fill_trans_id
      integer(ik) :: m,n,counter
      !
      counter = 0_rk
      !
      do n=1,num_channels
      do m=(n+1),num_channels
        !
        counter = counter + 1_ik
        !
        trans_id(n,m) = counter;    ! This assumes that the ionic states are real
        trans_id(m,n) = counter;    ! which implies symmetric matrix elements
        !
      end do
      end do
      !
      if (counter .ne. total_rdm_fields) stop 'transition_index - problem with number of off-diagonal elements...'
      !
    end subroutine fill_trans_id
    !
    !  Build RDF (Recombination dipole field) from the Dyson and exchange orbitals
    !  All we do here is simply make a copy of the "craddle" orbitals.
    !
    subroutine build_rdf
      f_radf     = f_exchange
      f_exchange = -1
    end subroutine build_rdf
!   !
!   !  Build RAF (Recombination acceleration field) from the Dyson and exchange orbitals
!   !  The routine below misses channel-coupling contributions. Until somebody derives
!   !  and implements them, it should not be used.
!   !
!   subroutine build_raf
!     integer(ik) :: ic, nc
!     !
!     if (f_core_pot(1)<=0) then
!       write (out,"('Effective potential is not available; can''t take gradients!')") 
!       stop 'build_raf - Missing the potential'
!     end if
!     !
!     !  "Exchange" terms are not defined for this matrix elements - wipe out and reuse
!     !  the "exchange" fields. Also calculate acceleration of the dyson orbital. It should
!     !  be zero for bound states, but one never knows what we are going to get...
!     !
!     !  Also keep in mind that the field in f_dyson needs to be rescaled by eta*sqrt(N_electrons)
!     !  to get the true Dyson orbital
!     !
!     loop_channels: do nc=1,num_channels
!       build_component: do ic=1,3
!         call FieldGradientComponent(dir=ic,src=f_core_pot(nc),dst=f_exchange(ic,nc))
!         call FieldMul(dst=f_exchange(ic,nc),src=f_dyson(nc))
!         call FieldScale(dst=f_exchange(ic,nc),con=cmplx(eta(nc)*sqrt(N_electrons),0._rk,kind=rk))
!         acceleration_dyson(ic,nc) = FieldConjgIntegrate(left=f_exchange(ic,nc),right=f_dyson(nc))
!       end do build_component
!       write (out,"('Acceleration of the Dyson orbital ',i3,' [Bohr/au[t]^2]:',3(1x,g12.6))") nc, acceleration_dyson(:,nc)
!     end do loop_channels
!     !
!     !  Copy handles from f_exchange to f_radf, and kill values in f_exchange
!     !  to prevent reuse
!     !
!     f_radf     = f_exchange
!     f_exchange = -1
!   end subroutine build_raf
    !
    subroutine build_rdaf
      call TimerStart('Build RDAF')
      select case (harmonic_operator)
        case default
          write (out,"('Harmonic probe operator ""',a,'"" is not recognized. Oops.')") trim(harmonic_operator)
          stop 'build_rdaf - bad operator (1)'
        case ('none')
        case ('dipole')
          write (out,"(/'Using dipole form of the recombination matrix elements')")
          call build_rdf
        case ('acceleration')
          write (out,"(/'Using acceleration form of the recombination matrix elements')")
        ! call build_raf
          stop 'dyson_nchan%build_rdaf - acceleration of dipole is not implemented'
      end select
      write (out,"()")
      call TimerStop('Build RDAF')
    end subroutine build_rdaf
    !
    !  An input file for a multi-channel simulation can be *very* complicated; making a
    !  mistake is way too easy. Let's echo the input in a slightly different format,
    !  which could make it easier to spot mistakes.
    !
    subroutine summarize_inputs
       integer(ik)        :: nc, ncc, id
       integer(ik)        :: max_dyson_len  ! Longest name of the dyson_file among all channels
       character(len=100) :: fmt_buf
       character(len=100), parameter :: lb = ' ' ! A long blank string
       !
       max_dyson_len = maxval(len_trim(dyson_files(:num_channels)))
       !
       write (out,"()")
       write (fmt_buf,"('(1x,a2,1x,a12,1x,3(1x,a8),1x,a12,2x,a',i0,',2x,a)')") max_dyson_len
       write (out,fmt_buf) 'nc', 'Energy', 'Dip-X', 'Dip-Y', 'Dip-Z', 'Dys. ind.', 'Dyson orbital'//lb, 'Natural orbitals'
       write (out,fmt_buf) '--', '------', '-----', '-----', '-----', '---------', '-------------'//lb, '----------------'
       !
       write (fmt_buf,"('(1x,i2,1x,f12.6,1x,3(1x,f8.4),1x,4(1x,i2),2x,a',i0,',2x,a)')") max_dyson_len
       write (out,fmt_buf) 0, neutral_ener, neutral_dipole(:), 0, 0, 0, 0, 'n/a: neutral species'
       print_channels: do nc=1,num_channels
         write (out,fmt_buf) nc, ion_ener(nc), ion_dipole(:,nc), spin_parts(:,nc), dyson_files(nc), trim(natural_files(nc))
       end do print_channels
       write (out,"()")
       !
       if (.not.use_interchannel_coupling) return
       !
       if (use_rdm_coupling) then
         write (out,"()")
         write (fmt_buf,"('(1x,a2,1x,a2,1x,a4,2x,a',i0,',1x,a',i0,',1x,a)')") max_dyson_len, max_dyson_len
         write (out,fmt_buf) 'nc', 'cc', 'id', 'Target channel'//lb, 'Source channel'//lb, 'Transition 1-density'
         write (out,fmt_buf) '--', '--', '--', '--------------'//lb, '--------------'//lb, '--------------------'
         write (fmt_buf,"('(1x,i2,1x,i2,1x,i4,2x,a',i0,',1x,a',i0,',1x,a)')") max_dyson_len, max_dyson_len
         print_rdm_couplings: do nc=1,num_channels
           print_rdm_cpl_inner: do ncc=1,num_channels
             if (ncc==nc) cycle print_rdm_cpl_inner
             !
             id = trans_id(nc,ncc)
             write (out,fmt_buf) nc, ncc, id, dyson_files(nc), dyson_files(ncc), trim(rdm_files(id))
           end do print_rdm_cpl_inner
           write (out,"()")
         end do print_rdm_couplings
       end if
       !
       if (use_dipole_coupling) then
         write (out,"()")
         write (fmt_buf,"('(1x,a2,1x,a2,1x,a4,1x,3(1x,a8),2x,a',i0,',1x,a',i0,')')") max_dyson_len, max_dyson_len
         write (out,fmt_buf) 'nc', 'cc', 'id', 'Dipole-X', 'Dipole-Y', 'Dipole-Z', 'Target channel'//lb, 'Source channel'//lb
         write (out,fmt_buf) '--', '--', '--', '--------', '--------', '--------', '--------------'//lb, '--------------'//lb
         write (fmt_buf,"('(1x,i2,1x,i2,1x,i4,1x,3(1x,f8.4),2x,a',i0,',1x,a',i0,')')") max_dyson_len, max_dyson_len
         print_dipole_couplings: do nc=1,num_channels
           print_dipole_cpl_inner: do ncc=1,num_channels
             if (ncc==nc) cycle print_dipole_cpl_inner
             !
             id = trans_id(nc,ncc)
             write (out,fmt_buf) nc, ncc, id, transition_dipole(:,id), dyson_files(nc), dyson_files(ncc)
           end do print_dipole_cpl_inner
           write (out,"()")
         end do print_dipole_couplings
       end if
       !
    end subroutine summarize_inputs
    !
    !  Initializes the fields
    !
    subroutine setup_fields
      integer(ik) :: nc, ips
      !
      call TimerStart('setup_fields')
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call fftw_activate_threads(fft_thread_threshold)
      !
      !  Set number of required fields
      !
      !  total_rdm_fields: number of unique off-digonal matrix elements of ionic basis)
      !  total_fields: num_channels * [ dyson + 3*(exchange or radf) + 3*calT + core + total_pot + psi + psi2 + Hpsi + phiIm ] + total_rdm_fields + 1*abs_pot
      !
      write (out,"(/'Running an ',i4,'-channel simulation')") num_channels
      total_rdm_fields = num_channels*(num_channels-1)/2   
      total_fields     = 13 * num_channels + 1
      if (odx_plot_momentum) total_fields = total_fields + 1 ! Extra scratch field for FFT
      if (use_interchannel_coupling .and. use_rdm_coupling) then
        total_fields = total_fields + total_rdm_fields
      end if
      if (continue_flux) then
        total_fields = total_fields + product(flux_samples) * num_channels + 2
      end if
      !
      call initialize_grid
      !
      call initialize_rotation
      !
      !  All input data related to the molecular frame will need to be rotated;
      !  orbitals and potentials are handled elsewhere, which leaves us with the
      !  dipoles here.
      !
      call rotate_dipoles
      !     
      call fill_trans_id
      !
      if (verbose>=2) then
        call summarize_inputs
      end if
      !
      !  Load the ECP from the first channel; it's up to the user to make sure
      !  all channels use the same ECP.
      !
      call load_ion_ecp
      !
      !  Ion potentials: this step can be very expensive for a many-channel simulation
      !
      if (potentials_in/=' ') then
        call restart_ion_potentials
      else
        call load_ion_potentials
        if (potentials_out/=' ') then
          call checkpoint_ion_potentials
        end if
      end if
      !
      !  Orbitals - load Dyson and "craddle"/"exchange" orbital
      !
      call load_gamess_mos
      !
      !  Prepare for Hamiltonian evaluation
      !
      call prepare_dyson_fields
      !
      !  Build the RDF/RAF for harmonics calculations
      ! 
      call build_rdaf
      !
      !  Set potential fields
      !
      do nc=1,num_channels
        if (f_free<1) stop 'dyson_nchan%setup_fields - not enough fields for the laser potential'
        f_total_pot(nc) = f_table(f_free) ; f_free = f_free - 1
      end do
      !
      !  Set absorbing boundary
      !
      if (f_free<1) stop 'dyson_nchan%setup_fields - not enough fields for the absorbing boundary'
      f_abs_pot   = f_table(f_free) ; f_free = f_free - 1
      call buildCAP(f_abs_pot)
      !
      !  Set propagation wavefunction fields
      !
      do nc=1,num_channels
        if (f_free<4) stop 'dyson_nchan%setup_fields - not enough fields for the propagation w.f. fields'
        f_psi(nc)   = f_table(f_free) ; f_free = f_free - 1
        f_psi2(nc)  = f_table(f_free) ; f_free = f_free - 1
        f_Hpsi(nc)  = f_table(f_free) ; f_free = f_free - 1
        f_phiIm(nc) = f_table(f_free) ; f_free = f_free - 1
      end do
      !
      !  Prepare for calculation of photoelectron spectra. Photoelectron spectra are
      !  per-channel.
      !
      if (continue_flux) then
        write (out,"('Using ',i6,' fields for flux continuation')") product(flux_samples) * num_channels
        allocate (f_flux_momentum(product(flux_samples),num_channels))
        flux_channels: do nc=1,num_channels
          flux_momentum: do ips=1,size(f_flux_momentum,dim=1)
            if (f_free<1) stop 'dyson_nchan%setup_fields - no fields for momentum continuation'
            f_flux_momentum(ips,nc) = f_table(f_free) ; f_free = f_free - 1
            call FieldZero(f_flux_momentum(ips,nc))
          end do flux_momentum
        end do flux_channels
        call FluxInitialize(mode='3D Cartesian',tMin=tMin,tMax=tMax,dt=dt, &
                            a0=(/0._ark,0._ark,0._ark/),evaluateF=ElectricField3D, &
                            flux_guard=flux_guard,sensors=flux_sensors,samples=flux_samples,max_p=flux_max_p, &
                            grid_limits=flux_grid_limits)
      end if
      !
      write (out,"()")
      write (out,"('[Grid fields free at end of setup_fields: ',i3,']')") f_free
      write (out,"()")
      !
      call TimerStop('setup_fields')
    end subroutine setup_fields
    !
    !  Finds the adiabatic field-dressed ionic states
    !
    subroutine get_dressed_states
      integer(ik) :: ne, ne2
      !
      call TimerStart('get_dressed_states')
      !
      ! Setup hamiltonian matrix using VecsD as storage
      !
      VecsD = 0._rk
  
      do ne=1,num_channels
        VecsD(ne,ne) = ion_ener(ne) + CurrentElaser*dot_product(ion_dipole(:,ne),Edir)
      end do
  
      do ne=1,num_channels
        do ne2=(ne+1),num_channels
          VecsD(ne,ne2) = + CurrentElaser*dot_product(Edir,transition_dipole(:,trans_id(ne,ne2)))
          VecsD(ne2,ne) = + CurrentElaser*dot_product(Edir,transition_dipole(:,trans_id(ne2,ne)))
        end do
      end do
      !
      ! Get eigenstates
      !
      call lapack_syev(VecsD(1:num_channels,1:num_channels),EnerD(1:num_channels))
  
      call TimerStop('get_dressed_states')
    end subroutine get_dressed_states
    ! 
    !  Calculate recombination matrix element. The correct choice for the observables
    !  is not obvious in this treatment. The expression we use is not final. Treat 
    !  all recombination results as preliminary and qualitative only!
    !
    !  We assume that that the wavefunction is determined by f_psi2, AmpD2, and AmpB2
    ! 
    subroutine EvaluateRecombination
      complex(rk) :: cd(3)        ! Complex recombination matrix element
      integer(ik) :: ic, nc, ncc
      !
      call TimerStart('EvaluateRecombination')
      cd = 0
      select case (harmonic_operator) 
        case default
          stop 'Recombination - unknown harmonic_operator value'
        case ('none')
!       case ('acceleration')
!         !
!         !  The form of this operator is most certainly NOT CORRECT. At the moment,
!         !  we simply take the integral <\phi_{Dyson}|\partial{V}{dr}|\chi>
!         !
!         do nc=1,num_channels
!           acceleration: do ic=1,3
!             cd(ic) = cd(ic) + FieldConjgIntegrate(left=f_radf(ic,nc),right=f_psi2(nc))
!           end do acceleration
!         end do 
        case ('dipole')
          !
          !  The dipole form is the result we get by taking expectation value of the
          !  proxy wavefunction, then multiplying it by square root of the number of
          !  electrons. It is not quite correct, as we should have applied the 
          !  antisymmetrizer instead of simply scaling - but there is a bit of
          !  difficulty with working out the correct antisymmetrizer.
          !
          recombination_dipole: do nc=1,num_channels
            !
            !  Recombination to the Dyson part within the same channel
            !
            call FieldNorm2Multipoles(left=f_dyson(nc),right=f_psi2(nc),mult=cd)
            cd = cd * conjg(ampD2(nc))
            !
            !  Recombination to the neutral residue, still within the same channel
            !
            dipole_neutral: do ic=1,3
              cd(ic) = cd(ic) + conjg(calN*AmpB2)*FieldConjgIntegrate(left=f_radf(ic,nc),right=f_psi2(nc))
            end do dipole_neutral
            rec_channel(:,nc,nc) = real(cd,kind=rk) * sqrt(N_electrons*nspin) * 2._rk
            !
            !  Now the cross-channel terms
            !
            cross_dipole: do ncc=1,num_channels
              if (nc==ncc) cycle cross_dipole
              !
              !  Note the change is sign: the rest of our recombination dipoles does not
              !  include (-1) for the electron charge, so flit the sign here as well.
              !
              cd = -transition_dipole(:,trans_id(nc,ncc))
              cd = cd * conjg(ampD2(nc)-calN*AmpB2*eta(nc)) * FieldConjgIntegrate(left=f_dyson(nc),right=f_psi2(ncc))
              rec_channel(:,nc,ncc) = real(cd,kind=rk) * sqrt(N_electrons*nspin) * 2._rk
            end do cross_dipole
          end do recombination_dipole
      end select
      !
      recombination = sum(sum(rec_channel(:,1:num_channels,1:num_channels),dim=3),dim=2)
      call TimerStop('EvaluateRecombination')
    end subroutine EvaluateRecombination
    !
    !  Prepare the local part of the operator in the square brackets in eq. (32)
    !  LaserPlusIon includes E^I_m (eq. 33) and F(t) . (r_n - d^I_{mm}) (eq. 32)
    !  The rest of eq. 33 (V_{nuc} and V^H_{mm}) are in f_core_pot, and are handled by AXPY.
    !
    subroutine prepare_channel_potentials
      integer(ik) :: nc
      !
      call TimerStart('Per-channel potentials')
      do nc=1,num_channels
        IonEnergyShift = CurrentElaser * dot_product(Edir, ion_dipole(:,nc)) + ion_ener(nc)  ! (IonEnergyShift used in LaserPlusIon)
        call FieldInit(f_total_pot(nc),LaserPlusIon,mask=f_psi(nc))
        call FieldAXPY(cmplx(1._rk,0._rk,kind=rk),src=f_core_pot(nc),dst=f_total_pot(nc))
      end do
      call TimerStop('Per-channel potentials')
    end subroutine prepare_channel_potentials
    !
    !  Propagator from the 
    !     Psi(t) = b(t)|\tilde N> + \Sum_m [ a_m(t) |D> + |\chi_m(t)>|I_m> ]
    !  formulation
    !
    !  This routine implements eq. 28 of the paper.
    !
    !  Per-channel potentials at this time step are expected in f_total_pot
    !
    !  At the entry to performDysonTimeStep, f_psi2 contains wavefunction at the current time
    !  (which corresponds to the electric field in CurrentElaser). f_psi contains wavefunction
    !  at the previous time step (time-dt).
    !
    !  Upon exit, f_psi contains wavefunction at the current time step, while f_psi2 is now
    !  the wavefunction at the next time step (time+dt).
    ! 
    subroutine performDysonTimeStep(dt,norm,bFirstStep)
      real(rk), intent(in)    :: dt                       ! Desired time step
      real(rk), intent(out)   :: norm(max_channels)       ! norm of the propagated W.F.
      logical, intent(in)     :: bFirstStep               ! Flag to use 'first step' mode (i.e. euler, not leap frog)
      !
      real(rk)                :: norm0(max_channels)      ! Temp for normalization calculations
      complex(rk)             :: dyson_amp(max_channels)  ! Dyson component after applying the Hamiltonian 
      complex(rk)             :: calT_projections         ! Projection of |Phi^I_m> onto calT
      integer(ik)             :: temp   
      complex(rk)             :: tempcx   
      integer(ik)             :: nc, ncc
      real(rk)                :: dnorm    ! Temp for absorbed norm calculations
      !
      call TimerStart('performDysonTimeStep')
      !
      ! Make |Phi^I_m>  - Eq. (31). AmpD2 is a_m(t).
      !
      do nc=1,num_channels
        call FieldCopy(f_psi2(nc),f_phiIm(nc))
        call FieldAXPY(AmpD2(nc),f_dyson(nc),f_phiIm(nc))
      end do
      !
      calT_projections = 0._rk
      do nc=1,num_channels
        !
        ! Make current calT and get <calT|Phi^I_m> projection - Eq. (36)
        ! The ingredients in f_calT(1:2) have been constructed in prepare_dyson_fields_calTmain
        ! The f_cal(2) part is proportional to (-CurrentElaser), so all is fine.
        !
        call FieldCopy(f_calT(1,nc),f_calT(3,nc))
        call FieldAXPY(cmplx(-CurrentElaser,0._rk,kind=rk),f_calT(2,nc),f_calT(3,nc))
        !
        ! The second and third parts of Eq. 28a, rewritten a-la eq. 30a.
        ! Since f_calT_k amounts to <{\tilde N}|H^F*t)|I_k>, acting on |phiIm> gives us the
        !       two sums we are interested in.
        !
        calT_projections = calT_projections + nspin*FieldConjgIntegrate(f_calT(3,nc),f_phiIm(nc))
      end do
      !
      ! Apply Hamiltonian - Eq. 32. This is the first contribution in Eqs. 30b/30c
      ! Add multi-channel coupling
      ! And add calT to Hpsi - second contribution in Eqs. 30b/30c
      !
      do nc=1,num_channels
        !
        !  Intra-channel contribution: H^F (|S_m>a_m + |X_m>) from eqs. 28b and c
        !
        call QMHpsi(mass=1._rk,pot=f_total_pot(nc),psi=f_phiIm(nc),Hpsi=f_Hpsi(nc))
        call ecp_apply(src=f_phiIm(nc),dst=f_Hpsi(nc),ecp=ecp)
        !
        ! inter-channel coupling - off-diagonal terms in 28b/c
        !
        if (use_interchannel_coupling) then  
          do ncc=1,num_channels
            if (ncc == nc) cycle
            if (use_dipole_coupling) then
              !
              !  Non-adiabatic coupling between ion cores - the last term in A8 and A17
              !
              tempcx = cmplx(CurrentElaser*dot_product(Edir,transition_dipole(:,trans_id(nc,ncc))),0_rk,kind=rk)
              call FieldAXPY(tempcx,f_phiIm(ncc),f_Hpsi(nc))
            end if
            if (use_rdm_coupling) then
              !
              !  Ion-core transitions due to interactions with the continuum - the first term in A8 and A17
              !
              call FieldMulAdd(f_trans_pot(trans_id(nc,ncc)),f_phiIm(ncc),f_Hpsi(nc))
            end if
          end do
        end if
        !
        ! Add calT to Hpsi - this is the first term in eqs. 28b/c
        !
        call FieldAXPY(AmpB2,f_calT(3,nc),f_Hpsi(nc))
      end do
      !
      ! At this point, we've constructed the right-hand side of both eqs. 28b and 28c.
      ! The result has to be separated into the "dyson" and "continuum" channels.
      !
      do nc=1,num_channels
        ! Get dyson projection and remove from HPsi
        ! Finish construction of the RHS in Eq. 28b
        dyson_amp(nc) = FieldConjgIntegrate(f_dyson(nc),f_Hpsi(nc)) 
        ! Finish construction of the RHS in Eq. 28c
        call FieldAXPY(-dyson_amp(nc),f_dyson(nc),f_Hpsi(nc))
      end do
      ! Propagate the wavefunction at T=t-dt
      ! The equations are 28c, 28b, and 28a (in this order)
      if (bFirstStep) then
        !
        !  Simple first-order step
        !
        do nc=1,num_channels
          call FieldAXPY(cmplx(0._rk,-1._rk,kind=rk)*dt,f_Hpsi(nc),f_psi(nc))
          AmpD(nc) = AmpD2(nc) + cmplx(0._rk,-1._rk,kind=rk)*dt * dyson_amp(nc)
        end do
        ! Note that eq. 38 in the paper has a + sign for Elaser here.
        AmpB = AmpB2 + cmplx(0._rk,-1._rk,kind=rk)*dt * ((HtildeN_main-CurrentElaser*HtildeN_F)*AmpB2 + calT_projections) 
      else
        !
        !  Leap-frog
        !
        do nc=1,num_channels
          call FieldAXPY(cmplx(0._rk,-2._rk,kind=rk)*dt,f_Hpsi(nc),f_psi(nc))
          AmpD(nc) = AmpD(nc) + cmplx(0._rk,-2._rk,kind=rk)*dt * dyson_amp(nc)
        end do
        AmpB = AmpB + cmplx(0._rk,-2._rk,kind=rk)*dt * ((HtildeN_main-CurrentElaser*HtildeN_F)*AmpB2 + calT_projections) 
      end if
                                                         
      ! Apply the absorbing boundary
      do nc=1,num_channels
        norm0(nc) = real(FieldConjgIntegrate(f_psi(nc),f_psi(nc)),kind=rk)
        call FieldMul(f_abs_pot, f_psi(nc) )
        call FieldMul(f_abs_pot, f_psi2(nc))
        norm(nc) = real(FieldConjgIntegrate(f_psi(nc),f_psi(nc)),kind=rk)
        dnorm = (norm0(nc)-norm(nc))
        AbsorbedNorm(nc) = AbsorbedNorm(nc) + dnorm
      end do
      
      !  Swap field for the wavefunctions
      do nc=1,num_channels
        temp = f_psi(nc); f_psi(nc) = f_psi2(nc); f_psi2(nc) = temp
        tempcx = AmpD(nc); AmpD(nc) = AmpD2(nc); AmpD2(nc) = tempcx
      end do
      tempcx = AmpB; AmpB = AmpB2; AmpB2 = tempcx
      call TimerStop('performDysonTimeStep')
    end subroutine performDysonTimeStep
    !
    !  performFluxTimeStep expects to find current wavefunction in f_psi,
    !  with CurrentEField and CurrentAField corresponding to the current
    !  time as well.
    !
    subroutine performFluxTimeStep
      integer(ik)        :: f_flux_scr1  ! Fields for calculating the current step's
      integer(ik)        :: f_flux_scr2  ! contributions the photoelectron spectrum
      integer(ik)        :: nc           ! current channel
      !
      if (.not.continue_flux) return
      call TimerStart('performFluxTimeStep')
      if (f_free<2) stop 'dyson_nchan%performFluxTimeStep - out of scratch for the flux'
      f_flux_scr1 = f_table(f_free) ; f_free = f_free - 1
      f_flux_scr2 = f_table(f_free) ; f_free = f_free - 1
      !
      channel_flux: do nc=1,num_channels
        call FluxTimeStep(mode='3D Cartesian',f_psi=f_psi(nc),tstep=tstep,ctime=time,aField=CurrentAField, &
                          f_momentum=f_flux_momentum(:,nc),f_scr1=f_flux_scr1,f_scr2=f_flux_scr2)
      end do channel_flux
      !
      f_free = f_free + 2
      call TimerStop('performFluxTimeStep')
    end subroutine performFluxTimeStep
    !
    !  Job checkpoint - should match the restart routine below. The checkpoint file
    !  format is a bit complicated: In addition to the chi wavefunction fields, we also
    !  require amplitudes of other components. As a result, the complete checkpoint 
    !  consists of an ASCII file with brief description of the simulation, and two
    !  binary files per channel, containing the wavefunctions. It would be impractical
    !  to use FieldCheckpoint - there are too many scratch fields present.
    !
    subroutine checkpoint_simulation
      character(len=2*clen) :: filename
      character(len=2*clen) :: filename2
      character(len=clen)   :: fluxext
      integer(ik)           :: nc, ips
      !
      call TimerStart('Create checkpoint')
      write (filename,fmt=checkpoint_out) tstep
      if (filename==checkpoint_in) then
        write (out,"('Cowardly refusing to overwrite ',a,', which we just restarted from.')") &
               trim(filename)
      else
        open (unit_ckpt,form='formatted',status='replace',file=trim(filename))
        write (unit_ckpt,nml=dyson_nchan_checkpoint)
        close (unit_ckpt)
        write (out,"('Created checkpoint ',a)") trim(filename)
        write (out,"('Associated wavefunctions: ')")
        save_channels: do nc=1,num_channels
          write (filename2,"(a,'_c',i0)") trim(filename), nc
          write (out,"(t5,'channel ',i3,': ',a,', ',a)") nc, trim(filename2)//'_psi', trim(filename2)//'_psi2'
          call FieldExport('binary',f_psi (nc),trim(filename2)//'_psi')
          call FieldExport('binary',f_psi2(nc),trim(filename2)//'_psi2')
          if (continue_flux) then
            save_fluxes: do ips=1,size(f_flux_momentum,dim=1)
              write (fluxext,"('_flux',i0)") ips
              call FieldExport('binary',f_flux_momentum(ips,nc),trim(filename2)//trim(fluxext))
            end do save_fluxes
          end if
        end do save_channels
        write (out,"()")
      end if
      call TimerStop('Create checkpoint')
    end subroutine checkpoint_simulation
    !
    !  Absolutely no checking is done here. If anything goes wrong, it's your pigeon
    !
    subroutine restart_simulation
      character(len=2*clen) :: filename2
      character(len=clen)   :: fluxext
      integer(ik)           :: nc, ips
      !
      call TimerStart('Load checkpoint')
      write (out,"('Reading checkpoint ',a)") trim(checkpoint_in)
      open (unit_ckpt,form='formatted',action='read',position='rewind',file=trim(checkpoint_in))
      read (unit_ckpt,nml=dyson_nchan_checkpoint)
      close (unit_ckpt)
      write (out,"('Associated wavefunctions: ')")
      load_channels: do nc=1,num_channels
        write (filename2,"(a,'_c',i0)") trim(checkpoint_in), nc
        write (out,"(t5,'channel ',i3,': ',a,', ',a)") nc, trim(filename2)//'_psi', trim(filename2)//'_psi2'
        call FieldImport('binary',trim(filename2)//'_psi', (/f_psi (nc)/),(/1/))
        call FieldImport('binary',trim(filename2)//'_psi2',(/f_psi2(nc)/),(/1/))
        if (continue_flux .and. flux_restart) then
          load_fluxes: do ips=1,size(f_flux_momentum,dim=1)
            write (fluxext,"('_flux',i0)") ips
            call FieldImport('binary',trim(filename2)//trim(fluxext),(/f_flux_momentum(ips,nc)/),(/1/))
          end do load_fluxes
        end if
      end do load_channels
      write (out,"()")
      call TimerStop('Load checkpoint')
    end subroutine restart_simulation
    !
    !  Preview of adiabatic ion core dynamics
    !
    subroutine preview_ion_dynamics
      integer(ik) :: ir
      !
      if (ions_preview==' ') return
      if (num_channels>20) then
        write (out,"(/'WARNING: Too many states; preview_ion_dynamics will not report eigenvectors!'/)")
      end if
      !
      call TimerStart('Preview ion dynamics')
      !
      !  Let's just simulate the laser field, and dump it to a file
      !
      open (unit=unit_efield,form='formatted',status='replace',file=trim(ions_preview))
      time      = tMin
      tstep     = 0
      CurrentElaser = ElectricField(time)
      time_loop: do while(time<tMax)
        CurrentElaser = ElectricField(time)
        call get_dressed_states
        write (unit_efield,"(1x,i10,(t12,20(1x,g16.9)))") tstep+1, time, CurrentElaser, EnerD(1:num_channels)
        write_eigenvectors: do ir=1,num_channels
          write (unit_efield,"(' +',t12,20(1x,g16.9))") vecsD(ir,1:num_channels)
        end do write_eigenvectors
        time  = time + dt
        tstep = tstep + 1
      end do time_loop
      close (unit_efield)
      write (out,"(/'Adiabatic dynamics of bare ion core stored in file ',a)") trim(ions_preview)
      call TimerStop('Preview ion dynamics')
    end subroutine preview_ion_dynamics
    !
    subroutine report_header
      character(len=2*clen) :: filename
      !
      write (out,"('# [DATA] = Sum across all channels')")
      write (out,"('# [Cnnn] = Intra-channel')")
      write (out,"('#      ',a10,1x,a12,7(1x,a12))") &
             'step', 'time, au', '|AmpB|^2', '|AmpDyson|^2', 'AbsorbedNorm', &
             'ChiNorm', 'TotalNorm', 'Recombination', 'CurrentElaser'
      !
      if (detailed_output/=' ') then
        open (unit_detail,form='formatted',status='replace',recl=255,file=detailed_output)
        write (unit_detail,"('#',a4,1x,a8,11(1x,a16))") &
              ' 0 ', 'tstep', 'time', 'Re[AmpB2]', 'Im[AmpB2]', 'AbsorbedNorm', 'ChiNorm', &
              'Rec_d_x', 'Rec_d_y', 'Rec_d_z', 'CurrentElaser'
        write (unit_detail,"('#',a4,1x,a8,11(1x,a16))") &
              'chan', 'tstep', 'time', 'Re[AmpD2]', 'Im[AmpD2]', 'AbsorbedNorm', 'ChiNorm', &
              'Rec_self_x', 'Rec_self_y', 'Rec_self_z', 'Rec_intra_x', 'Rec_intra_y', 'Rec_intra_z'
      end if
      !
      if (output_prefix/=' ') then
        write(filename,"(a,'.ampd')") trim(output_prefix); open (unit_ampd,status='replace',file=filename)
        write(filename,"(a,'.chi ')") trim(output_prefix); open (unit_chi ,status='replace',file=filename)
        write(filename,"(a,'.abs ')") trim(output_prefix); open (unit_abs ,status='replace',file=filename)
        write(filename,"(a,'.rec ')") trim(output_prefix); open (unit_rec ,status='replace',file=filename)
      end if
      !
    end subroutine report_header
    !
    subroutine report_end
      if (detailed_output/=' ') then
        close (unit_detail)
      end if
      if (output_prefix/=' ') then
        close(unit_ampd)
        close(unit_chi )
        close(unit_abs )
        close(unit_rec )
      end if
    end subroutine report_end
    !
    subroutine report_results(final)
      logical, intent(in) :: final ! True if this is the final call in this simulation
      integer(ik)         :: nc
      real(rk)            :: sAmpD2, sAbsorbedNorm, sChiNorm, totNorm
      real(rk)            :: rec_projection
      character(len=256)  :: tag, tag2
      character(len=256)  :: formatstring
      !
      !  Prepare descriptive tag, even if no output will be produced.
      !  For the summary, it's enough to report projection of the recombination dipole on the laser field direction
      !
      rec_projection = dot_product(Edir,Recombination)
      sAmpD2         = sum(abs(AmpD2(1:num_channels))**2)
      sAbsorbedNorm  = sum(AbsorbedNorm(1:num_channels))
      sChiNorm       = sum(ChiNorm(1:num_channels))
      totNorm        = abs(AmpB2)**2 + nspin*(sChiNorm + sAmpD2 + sAbsorbedNorm)
      write (tag,"('[DATA] ',i10,1x,f12.4,7(1x,g12.6))") &
                 tstep+1, time, abs(AmpB2)**2, sAmpD2, sAbsorbedNorm, &
                 sChiNorm, totNorm, rec_projection, CurrentElaser
      !
      !  Report summary of the results
      !
      if (mod(tstep+1,report_each)==0) then
        write (out,"(a)") trim(tag)
        report_channels: do nc=1,num_channels
          rec_projection = dot_product(Edir,sum(rec_channel(:,nc,1:num_channels),dim=2))
          write (out,"('[C',i3.3,'] ',i10,1x,f12.4,1x,a12,6(1x,g12.6))") &
                 nc, tstep+1, time, 'n/a', abs(AmpD2(nc))**2, AbsorbedNorm(nc), ChiNorm(nc), &
                 nspin*(ChiNorm(nc)+abs(AmpD2(nc))**2+AbsorbedNorm(nc)), rec_projection, CurrentElaser
        end do report_channels
        !
        !  And a bit more detail if requested
        !
        if (output_prefix/=' ') then
          !
          write(formatstring,"('(f12.5,',i2,'(1x,g15.6e3))')") num_channels+3
          write(unit_ampd,formatstring) time, abs(AmpB2)**2, abs(AmpD2(1:num_channels))**2, totNorm, CurrentElaser
          !
          write(formatstring,"('(f12.5,',i2,'(1x,g15.6e3))')") num_channels+2
          write(unit_chi ,formatstring) time, ChiNorm(1:num_channels), totNorm, CurrentElaser 
          !
          write(formatstring,"('(f12.5,',i2,'(1x,g15.6e3))')") num_channels+2
          write(unit_abs ,formatstring) time, AbsorbedNorm(1:num_channels), totNorm, CurrentElaser
          !
          write(formatstring,"('(f12.5,1x,i3,',i2,'(1x,g15.6e3))')") 3*num_channels+2
          write_rec: do nc=1,num_channels
            write(unit_rec,formatstring) time, nc, rec_channel(1:3,nc,1:num_channels), totNorm, CurrentElaser
          end do write_rec
        end if
      end if 
      !
      !  Report detailed results at each time step. This file is for machine reading only.
      !
      if (detailed_output/=' ') then
        write (unit_detail,"(i3,1x,i10,9(1x,g16.10))") &
               0, tstep+1, time, AmpB2, sAbsorbedNorm, sChiNorm, Recombination, CurrentElaser
        detail_channels: do nc=1,num_channels
          write (unit_detail,"(i3,1x,i10,11(1x,g16.10))") &
                 nc, tstep+1, time, AmpD2(nc), AbsorbedNorm(nc), ChiNorm(nc), rec_channel(:,nc,nc), &
                 sum(rec_channel(:,nc,1:num_channels),dim=2)
        end do detail_channels
      end if
      !
      !  Plot if necessary
      !
      if (odx_plot_each>0 .or. (odx_plot_each==-1 .and. final)) then
        if (mod(tstep+1,odx_plot_each)==0 .or. final) then
          visualize_loop: do nc=1,num_channels
            write (tag2,"('C',i3.3,': ',a)") nc, trim(tag)
            call visualize_wavefunction(trim(tag2),f_psi2(nc))
            if (continue_flux) then
              call visualize_photoelectron(trim(tag2),f_flux_momentum(:,nc))
            end if
          end do visualize_loop
        end if
      end if
    end subroutine report_results
    !
    !  Time propagation. Expensive!
    !
    subroutine propagate
      integer(ik) :: nc
      !
      !  Start propagation
      !
      time      = tMin 
      tstep     = 0
      ChiNorm   = 0
      call FieldSetOuterWall('REFLECTING')
      !
      if (checkpoint_in/=' ') then
        !
        ! initial conditions - start from a checkpoint
        !
        call restart_simulation
      else
        !
        ! initial conditions - start from scratch
        !
        do nc=1,num_channels
          AmpD2(nc) = eta(nc)
          AmpD(nc)  = 0._rk
          call FieldZero(f_psi(nc) )
          call FieldZero(f_psi2(nc))
        end do
        AmpB2 = 1._rk / calN
        AmpB  = 0._rk
        !
        ! first step
        !
        CurrentElaser = ElectricField(time)
        CurrentEField = ElectricField3D(time)
        CurrentAField = 0._rk
        call prepare_channel_potentials
        call performDysonTimeStep(dt,ChiNorm,bFirstStep=.true.)
        call performFluxTimeStep
        call EvaluateRecombination
        AbsorbedNorm = 0._rk
      end if
      !
      !  Main loop - this is where we spend most of the time
      !
      write (out,"(/'Beginning time evolution'/)") 
      call TimerStart('Time evolution')
      call report_header
      time_loop: do while(time<tMax)
        !
        !  Periodic checkpoints
        !
        if (checkpoint_each>0) then
          if (mod(tstep+1,checkpoint_each)==0) then
            call checkpoint_simulation
          end if
        end if
        call report_results(time+dt>=tMax) ! Must call this BEFORE the new electric field gets set.
                                           ! The flag is intended to signal the last call to report_results
        !
        !  Time step
        !
        CurrentElaser = ElectricField(time)
        CurrentEField = ElectricField3D(time)
        call prepare_channel_potentials
        call performDysonTimeStep(dt,ChiNorm,bFirstStep=.false.)
        call performFluxTimeStep
        call EvaluateRecombination
        !
        !  Advance time, and repeat
        !
        time  = time + dt 
        tstep = tstep + 1
        CurrentAField = CurrentAField - CurrentEField*dt
        !
        !  HACK: Eliminate the ground state
        !
        if (time>=Tkill_bound_state .and. .not.bound_state_killed) then
          write (out,"(//10('WARNING: ELIMINATING THE GROUND STATE NOW'/)/)")
          ampD = 0 ; ampD2 = 0 ; ampB = 0 ; ampB2 = 0
          bound_state_killed = .true.
        end if
      end do time_loop
      call report_end
      !
      !  Write the final checkpoint
      !
      if (checkpoint_each>=0) then
        call checkpoint_simulation
      end if
      call TimerStop('Time evolution')
    end subroutine propagate
    !
    !
    !
    subroutine start
      integer(ik)         :: info
      
      call TimerStart('start')

      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=dyson_nchan_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=dyson_nchan_data)
      write (out,"()")
      !
      ! normalize polarization vector (just to be sure...)
      !
      Edir = Edir / sqrt(sum(Edir**2))
      !
      ! sanity control for inter-channel coupling
      !
      if (use_interchannel_coupling .and. .not.(use_rdm_coupling.or.use_dipole_coupling)) then
        use_interchannel_coupling = .false.
        write (out,"(/'WARNING: Setting use_interchannel_coupling = .FALSE.'/)")
      end if
      !
      call setup_fields

      call preview_efield

      call preview_ion_dynamics

      call TimerReport

      call propagate

      call TimerStop('start')
      
      call TimerReport
    end subroutine start
    
  end module dyson_nchan
  !
  !
  !
  subroutine driver
    use dyson_nchan
    use accuracy

    call accuracyInitialize

    call start

  end subroutine driver
