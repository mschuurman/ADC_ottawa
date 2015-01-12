!
!  Propagation of many-electron wavefunction in the basis of ionic states.
!  The theory is described in: M. Spanner and S. Patchkovskii, PRA 80, 063411 (2009).
!
!  This module implements a special case of a singlet/doublet molecule with a single
!  uncoupled ionic state (See section IIIA of the paper).
!
!  Input preparation for this program is rather non-trivial.
!
  module dyson_1chan
    use dyson_tools

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !       Keep in mind that some of the relevant constants are in dyson_tools.f90
    !
    integer(ik), parameter :: max_mos      = 7                   ! Max number of MOs we support: Dyson + 3x exchange + 3x calT ('transfer state')
    integer(ik), parameter :: max_fields   = max_mos + 8         ! max_mos + core + absorbing + laser + psi + psi2 + Hpsi + phiIm + visualization scr
    !
    !  ==== User-adjustable parameters =====
    !
    character(len=clen) :: dyson_file      = 'dyson.dat'           ! Name of the file containing MO coefficients
                                                                   ! for the Dyson and exchange orbitals
    character(len=clen) :: natural_file    = 'natural.dat'         ! Name of the file containing MO coefficients
                                                                   ! for the natural orbitals of the cation
                                                                   ! (used for calculating Hartree potential)
    integer(ik)         :: spin_parts(4)   = (/ 2, 4, 6, 8 /)      ! Orbital components to load from "dyson_file"
    integer(ik)         :: natural_count   =-1                     ! Number of natural orbitals in the density matrix
                                                                   ! Value of -1 (the default) requests that natural
                                                                   ! occupations are loaded from the same file as the
                                                                   ! natural orbitals themselves. Use of this option
                                                                   ! is strongly recommended!
    real(rk)            :: natural_occ(max_naturals)               ! Natural orbital occupation numbers. If natural_count
                                                                   ! is positive, natural_count values should be supplied
                                                                   ! on input.
    real(rk)            :: neutral_ener  =  0.0                    ! Neutral groundstate energy, note: zero is at ionic ground energy
    real(rk)            :: ion_ener      =  0.0                    ! energy of current ionic state, note: zero is at ionic ground energy
                        
                                                                   ! IMPORTANT: The dipole moments include electron charge (-1)
                                                                   ! They are otherwise in atomic units - NOT in Debye
    real(rk)            :: neutral_dipole(3) = (/ 0._rk, 0._rk, 0._rk /) 
    real(rk)            :: ion_dipole(3)     = (/ 0._rk, 0._rk, 0._rk /)
                        
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)         :: f_dyson                                 ! Dyson orbital. Eventually renormalized to the "source"
                                                                   ! form, so beware.
    integer(ik)         :: f_exchange(3)                           ! Exchange orbitals (shares storage with f_radf)
    integer(ik)         :: f_radf(3)                               ! Recombination dipole field or acceleration dipole field,
                                                                   ! depending on the value of harmonic_operator. Shares storage
                                                                   ! with f_exchange, so only one would be in use at any given time
    integer(ik)         :: f_calT(3)                               ! "transfer" state (1='main', 2='F part', 3= full)
                                                                   ! See eq. 36. The f_calT(1) part does not depend on the
                                                                   ! electric field. The f_calT(2) part is linear in the magnitude
                                                                   ! of the field. This partitioning is for the special case of
                                                                   ! the linear polarization of the electric field.
                                                                   ! IMPORTANT: Note that the f_calT(2) part is actually linear in (-CurrentElaser),
                                                                   ! IMPORTANT: which may cause some confusion if you are not careful!
    integer(ik)         :: f_core_pot                              ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)         :: f_total_pot                             ! total potential (laser and core and ion energy)
    integer(ik)         :: f_psi                                   ! present wavefunction 
    integer(ik)         :: f_psi2                                  ! and 'next' wavefunction
    integer(ik)         :: f_Hpsi                                  !  = H |psi>
    integer(ik)         :: f_phiIm                                 !  = |phi^I_m> = ampD |dyson> + |chi^I_m>
                        
    real(rk)            :: dipole_dyson(3)                         ! Dipole moment (e-Bohr) of the "source" (renormalized Dyson) orbital
                                                                   ! <dyson|r|dyson>. Note that the sign of electron charge is NOT included.
    real(rk)            :: acceleration_dyson(3)                   ! Acceleration of the Dyson orbital - only defined 
                                                                   ! if harmonic_operator='acceleration'
    complex(rk)         :: multipoles_ion(0:9)                     ! Multipole moments of the ion core (electrons only)
                                                                   ! Charge of an electron is taken as +1 for this field
                                                                   ! multipoles_ion(0) is the total number of electrons in the cation
    real(rk)            :: N_electrons                             ! Number of electrons in the neutral species
    real(rk)            :: eta                                     ! amplitude of normalized dyson state contained in the neutral
                                                                   !  -->  eta = (<normalized dyson|<I_m|) |neutral>
    real(rk)            :: calN                                    ! Normalization factor for |\tilde N> state
                        
    complex(rk)         :: d_Hm_d                                  ! = <dyson|H_m|dyson>
    complex(rk)         :: HtildeN_main                            ! = {\cal H}^{\tilde N}  --- 'main' part, see eq. 38
    complex(rk)         :: HtildeN_F                               ! = {\cal H}^{\tilde N}  --- 'F(t)' part, see eq. 38 (special case of
                                                                   !   the linear polarization)
                        
    complex(rk)         :: AmpD           = 1._rk                  ! Amplitude of the dyson state
    complex(rk)         :: AmpD2          = 1._rk                  ! Amplitude of the dyson state at 'next' time step
                                                                   ! In the paper, AmpD/AmpD2 is given by a_m(t)
                                                                   
    complex(rk)         :: AmpB           = 1._rk                  ! Amplitude of the |\tilde N> state
    complex(rk)         :: AmpB2          = 1._rk                  ! Amplitude of the |\tilde N> state at 'next' time step
                                                                   
    real(rk)            :: AbsorbedNorm   = 0._rk                  ! Amount of absorbed population (norm)
    real(rk)            :: recombination(3)                        ! Recombination dipole (harmonic_operator='dipole')
                                                                   ! or dipole acceleration (harmonic_operator='acceleration') at the
                                                                   ! last time step.
    real(rk)            :: ChiNorm                                 ! Norm of the excited/continuum part of the wavefunction
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /dyson_1chan_data/ verbose, &
                       n_points, box_extent, &
                       dyson_file, spin_parts, &
                       nspin, euler_angles, &
                       eps_hartree, sor_rate, &
                       natural_count, natural_file, natural_occ, &
                       v_xc, &
                       ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
                       tMin, tMax, dt, &
                       neutral_ener, ion_ener, &
                       ElectricFieldShape, &
                       Elaser, Edir, Wlaser, Plaser, TlaserOn, TlaserOff, RampTime, FlatTime, & 
                       FWHM, Tlaser, &
                       neutral_dipole, ion_dipole, abs_bound_kmin, abs_bound_type, &
                       harmonic_operator, &
                       detailed_output, report_each, &
                       potential_file, odx_plot_each, odx_plot_momentum, efield_preview, &
                       checkpoint_each, checkpoint_out, checkpoint_in, &
                       oversample_orbitals
    !
    !  ==== All parameters needed for checkpoint/restart are in the namelist below ====
    !
    namelist /dyson_1chan_checkpoint/ ampD, ampD2, ampB, ampB2, AbsorbedNorm, &
                                      recombination, ChiNorm, time, tstep
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Rotate static input data - at the moment, these are only the dipoles
    !
    subroutine rotate_dipoles
      logical :: go
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
      call do_dipole('   ion   ',    ion_dipole)
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
    !  Prepare ion-specific potential: core+Hartree
    !
    subroutine load_ion_potentials
      !
      call TimerStart('Load ion potentials')
      !
      !  Do we have natural orbitals to construct the density and the Hartree potential?
      !
      if (natural_count>max_naturals) then
        write (out,"('Too many natural orbitals. Increase ''max_naturals''" &
                   // " to at least ',i0,' and recompile')") natural_count
        stop 'amplitudes - too many natural orbitals'
      end if
      if (natural_count<0) then
        call gamess_load_natocc(trim(natural_file),natural_occ,natural_count)
      end if
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      if (f_free<1) stop 'amplitudes - not enough fields for core potential'
      f_core_pot = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_potential(natural_file,natural_occ(:natural_count), &
                                   f_table(:f_free),f_core_pot,v_xc=v_xc, &
                                   electronic_multipoles=multipoles_ion,rot=rotmat)
      if (potential_file/=' ') then
        call dump_potential(f_core_pot,potential_file) ! This call is likely to produce a HUGE file
      end if
      call TimerStop('Load ion potentials')
    end subroutine load_ion_potentials
    !
    !  Allocate fields for and load Dyson and exchange correction orbitals
    !
    subroutine load_gamess_mos
      real(rk)              :: norm
      !
      call TimerStart('Load Gamess MOs')
      !
      !  Load dyson and "exchange correction" orbitals
      !
      f_free = f_free - 4
      if (f_free<0) stop 'dyson_1chan%load_gamess_mos - out of fields'
      call FieldImport('GAMESS',dyson_file,f_table(f_free+1:f_free+4),spin_parts,rot=rotmat)
      f_dyson         = f_table(f_free+1)
      f_exchange(1:3) = f_table(f_free+2:f_free+4)
      write (out,"('Loaded orbital components ',4i3)") spin_parts
      norm = FieldNorm(f_dyson)**2
      write (out,"('<psid|psid> = ',f20.10)") norm
      !
      !  Rotate exchange correction vector fields
      !
      call FieldRotateVectorComponents(rotmat,f_exchange(1),f_exchange(2),f_exchange(3))
      !
      !  Report nuclei info
      !
      call report_structure(1_ik,trim(dyson_file))
      call TimerStop('Load Gamess MOs')
    end subroutine load_gamess_mos
    !
    !  Prepare "source" orbitals
    !
    subroutine prepare_dyson_fields_source
      !
      call TimerStart('Dyson fields: Source')
      !
      !  Remove the sqrt(N) normalization factor from Dyson and "exchange" orbitals
      !
      N_electrons = real(multipoles_ion(0),kind=rk) + 1._rk
      !
      call FieldScale(f_dyson,      cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
      call FieldScale(f_exchange(1),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
      call FieldScale(f_exchange(2),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))
      call FieldScale(f_exchange(3),cmplx(1._rk/sqrt(N_Electrons),0._rk,kind=rk))

      eta = FieldNorm(f_dyson)
      call FieldScale(f_dyson,cmplx(1._rk/eta,0._rk,kind=rk))

      calN = 1._rk / sqrt(1._rk-nspin*abs(eta)**2)
      !
      write (out,"()")
      write (out,"(' Number of electrons (N_electrons)    = ',g14.7)") N_electrons
      write (out,"(' Spin factor (nspin)                  = ',i7   )") nspin
      write (out,"(' ""Source"" orbital normalization (eta)   = ',g14.7)") eta
      write (out,"(' Neutral residue normalization (calN) = ',g14.7)") calN
      if (mod(nint(N_electrons)+nspin,2)/=0) then
        write (out,"(//'  WARNING: Number of electrons and number of spin components are not consistent!'//)")
      end if
      write (out,"(' Overlaps of renormalized dyson and ""craddle"" orbitals:')")
      write (out,"('    <D|C_x> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson,f_exchange(1))
      write (out,"('    <D|C_y> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson,f_exchange(2))
      write (out,"('    <D|C_z> = ',2(1x,g14.7))") FieldConjgIntegrate(f_dyson,f_exchange(3))
      write (out,"()")
      !
      !  Multipole moments for the Dyson orbital - we need the dipole
      !
      call FieldNormMultipoles(f_dyson,dipole_dyson)
      write (out,"('Average coordinate of the ""source"" orbital [Bohr]: ',3f12.6)") dipole_dyson
      call TimerStop('Dyson fields: Source')
    end subroutine prepare_dyson_fields_source
    !
    !  Build calT 'main' - see Eq. 36
    !
    subroutine prepare_dyson_fields_calTmain
      integer(ik) :: ic
      !
      call TimerStart('Dyson fields: calT(main)')
      f_free = f_free - 3
      if (f_free<0) stop 'dyson_1chan%prepare_dyson_fields - out of fields for calT storage'
      f_calT(1:3) = f_table(f_free+1:f_free+3)
      !
      !  Build calT 'main' - see Eq. 36
      !
      call QMHpsi(1._rk,f_core_pot,f_dyson,f_calT(1))                     ! 
      call ecp_apply(src=f_dyson,dst=f_calT(1),ecp=ecp)                   ! ECP is a non-local operator; it can't be folded into f_core_pot
      call FieldAXPY(cmplx(ion_ener,0._rk,kind=rk),f_dyson,f_calT(1))     !  Operator of Eq. 33 acting on the source orbital
      d_Hm_d = FieldConjgIntegrate(f_dyson,f_calT(1))                     ! <dyson|H_m|dyson> is used below for HtildeN_main (eq. 38)
      call FieldScale(f_calT(1),cmplx(-1._rk,0._rk,kind=rk))
      call FieldAXPY(cmplx(neutral_ener,0._rk,kind=rk),f_dyson,f_calT(1))
      call FieldScale(f_calT(1),cmplx(eta*calN,0._rk,kind=rk))            ! Finished construction of the field-independent part of eq. 36
      !
      !  Build calT 'F part' - see Eq. 36. Note that the construction below depends on electric field
      !  being linearly polarized, even though the theory itself allows general polarization.
      !
      call FieldCopy(src=f_dyson,dst=f_calT(2))
      call FieldScale(f_calT(2),cmplx(dot_product(Edir,ion_dipole)*eta*calN,0._rk,kind=rk)) ! PS dipole change
      ! The equation above _appears_ to differ from the paper, but is not:
      ! the overall sign becomes - once the electric 
      ! magnitude is included in performDysonTimeStep.
      !
      do ic=1,3
        call FieldAXPY(cmplx(calN*Edir(ic),0._rk,kind=rk),f_exchange(ic),f_calT(2))
        ! ????? Published eq. 36 has a minus sign above - however, this is OK, since the
        !       sign is recovered upon multiplication with electric field magnitude in 
        !       performDysonTimeStep.
      end do 
      call TimerStop('Dyson fields: calT(main)')
    end subroutine prepare_dyson_fields_calTmain
    !
    !  Build HtildeN, including spin factor - see Eq. 38
    !
    subroutine prepare_dyson_fields_HtildeN
      integer(ik) :: ic
      complex(rk) :: cx_tmp  ! just a temp complex value...
      !
      call TimerStart('Dyson fields: H(tildeN)')
      ! Build HtildeN_main
      HtildeN_main = abs(calN)**2 * ( neutral_ener + nspin*abs(eta)**2*(d_Hm_d - 2._rk*neutral_ener) )
      !
      ! Build HtildeN_F
      HtildeN_F = -dot_product(Edir,neutral_dipole) - nspin*abs(eta)**2 * dot_product(Edir,ion_dipole) &
                           - nspin*abs(eta)**2 * dot_product(Edir,dipole_dyson)  ! PS dipole change
                           ! Together with the minus sign for CurrentElaser in performDysonTimeStep,
                           ! this gives correct sign overall
      do ic=1,3
        cx_tmp = FieldConjgIntegrate(f_exchange(ic),f_dyson) * eta
        ! ?????? Published eq. 38 has a plus sign here?!
        ! See comment immediately above
        HtildeN_F = HtildeN_F - nspin*( cx_tmp + conjg(cx_tmp) ) * Edir(ic)
      end do 

      HtildeN_F = HtildeN_F * abs(calN)**2  ! Overall normalization factor in eq. 38
      call TimerStop('Dyson fields: H(tildeN)')
    end subroutine prepare_dyson_fields_HtildeN
    !
    !  Prepare "transfer" orbitals needed for Hamiltonian evaluation.
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
    !  Build RDF (Recombination dipole field) from the Dyson and exchange orbitals
    !  All we do here is simply make a copy of the "craddle" orbitals.
    !
    subroutine build_rdf
      f_radf     = f_exchange
      f_exchange = -1
    end subroutine build_rdf
    !
    !  Build RAF (Recombination acceleration field) from the Dyson and exchange orbitals
    !
    subroutine build_raf
      integer(ik) :: ic
      !
      if (f_core_pot<=0) then
        write (out,"('Effective potential is not available; can''t take gradients!')") 
        stop 'build_raf - Missing the potential'
      end if
      if (ecp%active) then
        write (out,"('Acceleration form of recombination dipole is not implemented for ECPs')")
        stop 'build_raf - Can''t handle ECPs!'
      end if
      !
      !  "Exchange" terms are not defined for this matrix elements - wipe out and reuse
      !  the "exchange" fields. Also calculate acceleration of the dyson orbital. It should
      !  be zero for bound states, but one never knows what we are going to get...
      !
      !  Also keep in mind that the field in f_dyson needs to be rescaled by eta*sqrt(N_electrons)
      !  to get the true Dyson orbital
      !
      build_component: do ic=1,3
        call FieldGradientComponent(dir=ic,src=f_core_pot,dst=f_exchange(ic))
        call FieldMul(dst=f_exchange(ic),src=f_dyson)
        call FieldScale(dst=f_exchange(ic),con=cmplx(eta*sqrt(N_electrons),0._rk,kind=rk))
        acceleration_dyson(ic) = real(FieldConjgIntegrate(left=f_exchange(ic),right=f_dyson),kind=rk)
      end do build_component
      !
      !  Copy handles from f_exchange to f_radf, and kill values in f_exchange
      !  to prevent reuse
      !
      f_radf     = f_exchange
      f_exchange = -1
      write (out,"('Acceleration of the Dyson orbital [Bohr/au[t]^2]:',3(1x,g12.6))") acceleration_dyson
    end subroutine build_raf
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
          call build_raf
      end select
      write (out,"()")
      call TimerStop('Build RDAF')
    end subroutine build_rdaf
    !
    !  Initializes the fields
    !
    subroutine setup_fields
      call TimerStart('setup_fields')
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      write (out,"(/'Running a single-channel simulation')")
      total_fields = max_fields
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
      !  Load the ECP
      !
      call load_ion_ecp
      !
      call load_ion_potentials
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
      if (f_free<2) stop 'dyson_1chan%setup_fields - not enough fields for the laser potential'
      f_total_pot = f_table(f_free) ; f_free = f_free - 1
      f_abs_pot   = f_table(f_free) ; f_free = f_free - 1
      !
      call buildCAP(f_abs_pot)
      !
      !  Set propagation wavefunction fields
      !
      if (f_free<4) stop 'dyson_1chan%setup_fields - not enough fields for the propagation w.f. fields'
      f_psi   = f_table(f_free) ; f_free = f_free - 1
      f_psi2  = f_table(f_free) ; f_free = f_free - 1
      f_Hpsi  = f_table(f_free) ; f_free = f_free - 1
      f_phiIm = f_table(f_free) ; f_free = f_free - 1
      !
      call TimerStop('setup_fields')
    end subroutine setup_fields
    ! 
    !  Calculate recombination matrix element. The correct choice for the observables
    !  is not obvious in this treatment. The expression we use is not final. Treat 
    !  all recombination results as preliminary and qualitative only!
    !
    !  We assume that that the wavefunction is determined by f_psi2, AmpD2, and AmpB2
    ! 
    subroutine EvaluateRecombination
      complex(rk) :: cd(3) ! Complex recombination matrix element
      integer(ik) :: ic
      !
      call TimerStart('EvaluateRecombination')
      cd = 0
      select case (harmonic_operator) 
        case default
          stop 'Recombination - unknown harmonic_operator value'
        case ('none')
        case ('acceleration')
          !
          !  The form of this operator is most certainly NOT CORRECT. At the moment,
          !  we simply take the integral <\phi_{Dyson}|\partial{V}{dr}|\chi>
          !
          acceleration: do ic=1,3
            cd(ic) = FieldConjgIntegrate(left=f_radf(ic),right=f_psi2)
          end do acceleration
        case ('dipole')
          !
          !  The dipole form is the result we get by taking expectation value of the
          !  proxy wavefunction, then multiplying it by square root of the number of
          !  electrons. It is not quite correct, as we should have applied the 
          !  antisymmetrizer instead of simply scaling - but there is a bit of
          !  difficulty with working out the correct antisymmetrizer.
          !
          call FieldNorm2Multipoles(left=f_dyson,right=f_psi2,mult=cd)
          cd = cd * conjg(ampD2)
          dipole: do ic=1,3
            cd(ic) = cd(ic) + conjg(calN*AmpB2)*FieldConjgIntegrate(left=f_radf(ic),right=f_psi2)
          end do dipole
          !
          !  Since this is a single-channel implementation, there are no cross-channel terms
          !
          cd = sqrt(N_electrons)*cd
      end select
      !
      !  The observable is (cd + Conjg(cd)), once for each spin component
      !
      recombination = 2._rk * sqrt(real(nspin,kind=rk)) * real(cd,kind=rk)
      call TimerStop('EvaluateRecombination')
    end subroutine EvaluateRecombination
    !
    !  Prepare the local part of the operator in the square brackets in eq. (32)
    !  LaserPlusIon includes E^I_m (eq. 33) and F(t) . (r_n - d^I_{mm}) (eq. 32)
    !  The rest of eq. 33 (V_{nuc} and V^H_{mm}) are in f_core_pot, and are handled by AXPY.
    !
    subroutine prepare_channel_potentials
      call TimerStart('Per-channel potentials')
      IonEnergyShift = CurrentElaser * dot_product(Edir, ion_dipole) + ion_ener  ! (IonEnergyShift used in LaserPlusIon)
      call FieldInit(f_total_pot,LaserPlusIon)
      call FieldAXPY(cmplx(1._rk,0._rk,kind=rk),src=f_core_pot,dst=f_total_pot)
      call TimerStop('Per-channel potentials')
    end subroutine prepare_channel_potentials
    !
    !  Propagator from the 
    !     Psi(t) = b(t)|\tilde N> + a_m(t) |D> + |\chi_m(t)>|I_m> 
    !  formulation
    !
    !  This routine implements eq. 29 of the paper.
    !
    !  The potential is in f_total_pot, to keep things superficially consistent with the
    !  multi-channel version.
    ! 
    subroutine performDysonTimeStep(dt,norm,bFirstStep)
      real(rk), intent(in)    :: dt                ! Desired time step
      real(rk), intent(out)   :: norm              ! norm of the propagated W.F.
      logical, intent(in)     :: bFirstStep        ! Flag to use 'first step' mode (i.e. euler, not leap frog)
      !
      real(rk)                :: norm0             ! Temp for normalization calculations
      complex(rk)             :: dyson_amp         ! Dyson component after applying the Hamiltonian
      complex(rk)             :: calT_projection   ! Projection of |Phi^I_m> onto calT
      integer(ik)             :: temp   
      complex(rk)             :: tempcx   
      
      call TimerStart('performDysonTimeStep')
      !
      ! Make |Phi^I_m>  - Eq. (31). AmpD2 is a_m(t).
      !
      call FieldCopy(f_psi2,f_phiIm)
      call FieldAXPY(AmpD2,f_dyson,f_phiIm)
      !
      ! Make current calT and get <calT|Phi^I_m> projection - Eq. (36)
      ! The ingredients in f_calT(1:2) have been constructed in prepare_dyson_fields_calTmain
      ! The f_cal(2) part is proportional to (-CurrentElaser), so all is fine.
      !
      call FieldCopy(f_calT(1),f_calT(3))
      call FieldAXPY(cmplx(-CurrentElaser,0._rk,kind=rk),f_calT(2),f_calT(3))
      !
      ! The second part of Eq. 30a
      calT_projection = nspin*FieldConjgIntegrate(f_calT(3),f_phiIm)

      ! Apply Hamiltonian - Eq. 32. This is the first contribution in Eqs. 30b/30c
      call QMHpsi(mass=1._rk,pot=f_total_pot,psi=f_phiIm,Hpsi=f_Hpsi)
      call ecp_apply(src=f_phiIm,dst=f_Hpsi,ecp=ecp)

      ! Add calT to Hpsi - second contribution in Eqs. 30b/30c
      call FieldAXPY(AmpB2,f_calT(3),f_Hpsi)

      ! Get dyson projection and remove from HPsi
      ! Finish construction of the RHS in Eq. 30b
      dyson_amp = FieldConjgIntegrate(f_dyson,f_Hpsi) 
      ! Finish construction of the RHS in Eq. 30c
      call FieldAXPY(-dyson_amp,f_dyson,f_Hpsi)

      ! Propagate the wavefunction at T=t-dt
      ! The equations are 30c, 30b, and 30a (in this order)
      if (bFirstStep) then
        call FieldAXPY(cmplx(0._rk,-1._rk,kind=rk)*dt,f_Hpsi,f_psi)
        AmpD = AmpD2 + cmplx(0._rk,-1._rk,kind=rk)*dt * dyson_amp
        ! ?????? Note that eq. 38 in the paper has a + sign for Elaser here
        AmpB = AmpB2 + cmplx(0._rk,-1._rk,kind=rk)*dt * ((HtildeN_main-CurrentElaser*HtildeN_F)*AmpB2 + calT_projection) 
      else
        call FieldAXPY(cmplx(0._rk,-2._rk,kind=rk)*dt,f_Hpsi,f_psi)
        AmpD = AmpD + cmplx(0._rk,-2._rk,kind=rk)*dt * dyson_amp
        AmpB = AmpB + cmplx(0._rk,-2._rk,kind=rk)*dt * ((HtildeN_main-CurrentElaser*HtildeN_F)*AmpB2 + calT_projection) 
      end if
                                                         
      ! Apply the absorbing boundary
      norm0 = FieldNorm(f_psi)**2
      call FieldMul(f_abs_pot, f_psi )
      call FieldMul(f_abs_pot, f_psi2)
      norm = FieldNorm(f_psi)**2
      AbsorbedNorm = AbsorbedNorm + (norm0-norm)
      
      !  Swap field for the wavefunctions
      temp = f_psi; f_psi = f_psi2; f_psi2 = temp
      tempcx = AmpD; AmpD = AmpD2; AmpD2 = tempcx
      tempcx = AmpB; AmpB = AmpB2; AmpB2 = tempcx
      !
      call TimerStop('performDysonTimeStep')
    end subroutine performDysonTimeStep
    !
    !  Job checkpoint - should match the restart routine below. The checkpoint file
    !  format is a bit complicated: In addition to the chi wavefunction fields, we also
    !  require amplitudes of other components. As a result, the complete checkpoint 
    !  consists of an ASCII file with brief description of the simulation, and two
    !  binary files, containing the wavefunctions. It would be impractical to use 
    !  FieldCheckpoint - there are too many scratch fields present.
    !
    subroutine checkpoint_simulation
      character(len=2*clen) :: filename
      !
      call TimerStart('Create checkpoint')
      write (filename,fmt=checkpoint_out) tstep
      if (filename==checkpoint_in) then
        write (out,"('Cowardly refusing to overwrite ',a,', which we just restarted from.')") &
               trim(filename)
      else
        open (unit_ckpt,form='formatted',status='replace',file=trim(filename))
        write (unit_ckpt,nml=dyson_1chan_checkpoint)
        close (unit_ckpt)
        call FieldExport('binary',f_psi, trim(filename)//'_psi')
        call FieldExport('binary',f_psi2,trim(filename)//'_psi2')
        write (out,"('Created checkpoint ',a)") trim(filename)
        write (out,"('Associated wavefunctions: ',a,', ',a)") &
               trim(filename)//'_psi', trim(filename)//'_psi2'
      end if
      call TimerStop('Create checkpoint')
    end subroutine checkpoint_simulation
    !
    !  Absolutely no checking is done here. If anything goes wrong, it's your pigeon
    !
    subroutine restart_simulation
      call TimerStart('Load checkpoint')
      write (out,"('Reading checkpoint ',a)") trim(checkpoint_in)
      write (out,"('Associated wavefunctions: ',a,', ',a)") &
             trim(checkpoint_in)//'_psi', trim(checkpoint_in)//'_psi2'
      open (unit_ckpt,form='formatted',action='read',position='rewind',file=trim(checkpoint_in))
      read (unit_ckpt,nml=dyson_1chan_checkpoint)
      close (unit_ckpt)
      call FieldImport('binary',trim(checkpoint_in)//'_psi', (/f_psi /),(/1/))
      call FieldImport('binary',trim(checkpoint_in)//'_psi2',(/f_psi2/),(/1/))
      call TimerStop('Load checkpoint')
    end subroutine restart_simulation
    !
    subroutine report_header
      write (out,"('#      ',a10,1x,a12,7(1x,a12))") &
             'step', 'time, au', '|AmpB|^2', '|AmpDyson|^2', 'AbsorbedNorm', &
             'ChiNorm', 'TotalNorm', 'Recombination', 'CurrentElaser'
      if (detailed_output/=' ') then
        open (unit_detail,form='formatted',status='replace',recl=255,file=detailed_output)
        write (unit_detail,"('#',1x,a8,11(1x,a16))") &
              'tstep', 'time', 'Re[AmpB2]', 'Im[AmpB2]', 'Re[AmpD2]', 'Im[AmpD2]', &
              'AbsorbedNorm', 'ChiNorm', 'Rec_d_x', 'Rec_d_y', 'Rec_d_z', 'CurrentElaser'
      end if
    end subroutine report_header
    !
    subroutine report_end
      if (detailed_output/=' ') then
        close (unit_detail)
      end if
    end subroutine report_end
    !
    subroutine report_results
      real(rk)           :: rec_projection
      character(len=256) :: tag
      !
      !  Prepare descriptive tag, even if no output will be produced.
      !  For the summary, it's enough to report projection of the recombination dipole on the laser field direction
      !
      rec_projection = dot_product(Edir,Recombination)
      write (tag,"('[DATA] ',i10,1x,f12.4,7(1x,g12.6))") tstep+1, time, abs(AmpB2)**2, abs(AmpD2)**2, AbsorbedNorm, &
                                 ChiNorm, abs(AmpB2)**2 + nspin*(ChiNorm+abs(AmpD2)**2+AbsorbedNorm), rec_projection, CurrentElaser
      !
      !  Report summary of the results
      !
      if (mod(tstep+1,report_each)==0) then
        write (out,"(a)") trim(tag)
      end if 
      !
      !  Report detailed results
      !
      if (detailed_output/=' ') then
        write (unit_detail,"(i10,11(1x,g16.10))") tstep+1, time, AmpB2, AmpD2, AbsorbedNorm, ChiNorm, Recombination, CurrentElaser
      end if
      !
      !  Plot if necessary
      !
      if (odx_plot_each>0) then
        if (mod(tstep+1,odx_plot_each)==0) then
          call visualize_wavefunction(trim(tag),f_psi2)
        end if
      end if
    end subroutine report_results
    !
    !  Time propagation. Expensive!
    !
    subroutine propagate
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
        AmpD2 = eta
        AmpD  = 0._rk
        AmpB2 = sqrt(1-nspin*abs(eta)**2)
        AmpB  = 0._rk
        !
        call FieldZero(f_psi )
        call FieldZero(f_psi2)
        !
        ! first step
        !
        CurrentElaser = ElectricField(time)
        call prepare_channel_potentials
        call performDysonTimeStep(dt,ChiNorm,bFirstStep=.true.)
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
        call report_results ! Must call this BEFORE the new electric field gets set.
        !
        !  Time step
        !
        CurrentElaser = ElectricField(time)
        call prepare_channel_potentials
        call performDysonTimeStep(dt,ChiNorm,bFirstStep=.false.)
        call EvaluateRecombination
        !
        !  Advance time, and repeat
        !
        time  = time + dt 
        tstep = tstep + 1
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
      read (input,nml=dyson_1chan_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=dyson_1chan_data)
      write (out,"()")
      !
      ! normalize polarization vector (just to be sure...)
      !
      Edir = Edir / sqrt(sum(Edir**2))
      
      call setup_fields

      call preview_efield

!     call WriteSlice(2,0,f_dyson,"slices/dyson.dat")
!     call WriteSlice(2,0,f_core_pot,"slices/core.dat")
!     call WriteSlice(2,0,f_abs_pot,"slices/abs.dat")

      call TimerReport

      call propagate

      call TimerStop('start')
      
      call TimerReport
    end subroutine start
    
  end module dyson_1chan
  !
  !
  !
  subroutine driver
    use dyson_1chan
    use accuracy

    call accuracyInitialize

    call start

  end subroutine driver
