!
!  Multigrid test - calculation of relative HHG intensitions of chiral isomers
!
!  WARNING: The version prior to Apr. 10, 2013 was missing the Jacobian in the
!  WARNING: angular integration part. Oops. We were also integrating over [0:2Pi]
!  WARNING: for beta. Oops again.
!
!  WARNING: The version prior to June 20, 2013 was missing the ionization phase
!  WARNING: factor, which must change sign across nodal planes in the Dyson
!  WARNING: orbital. Oops.
!
  module chiral
    use accuracy
    use timer
    use math
    use sfa
    use import_gamess
    implicit none
    private
    public run_chiral
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: iu_input        = 34_ik                 ! Unit number to use for reading ionization/recombination fits
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose            = 0                     ! Verbosity level
    real(rk)            :: omega              = 0.057_rk              ! Laser field frequency
    real(rk)            :: e0x                = 0.05_rk               ! Laser field amplitude, X direction
    real(rk)            :: e0y                = 0.00_rk               ! Laser field amplitude, Y direction
    real(rk)            :: ip                 = 0.50_rk               ! Ionization potential
    real(rk)            :: ion_dipole(3)      = 0.00_rk               ! Dipole moment of the residual ion in atomic units,
                                                                      ! using the same origin used for calculating
                                                                      ! ionization yield and recombination cross-sections
                                                                      ! Nuclear conributions must be included.
    real(rk)            :: ion_alpha(3,3)     = 0.00_rk               ! Ion polarizability tensor, in atomic units
    real(rk)            :: phr_min            = 1.5707963268_rk       ! Recollision phase, min
    real(rk)            :: phr_max            = 7.85398163398_rk      ! Recollision phase, max
    real(rk)            :: phr_step           = 0.01_rk               ! Recollision phase, step
    real(rk)            :: phi_min            =-1.5707963268_rk       ! Minimum real part of the ionization phase
    real(rk)            :: phi_max            = 1.5707963268_rk       ! Maxsimum real part of the ionization phase
    logical             :: sfa_magnetic       = .true.                ! Include 1/c terms in the Hamiltonian
    logical             :: trajectories_only  = .false.               ! True if only the SFA trajectories are to be determined
    integer(ik)         :: n_return           = 1_ik                  ! Max. number of trajectories returning at given phase. 
                                                                      ! -1 means all within the phase window
    character(len=clen) :: ionization_file    = 'ionization.dat'      ! Data file containing ionization amplitudes (see read_ion() below)
    character(len=clen) :: recombination_file = 'recombination.dat'   ! Data file containing ionization amplitudes (see read_rec() below)
    !
    integer(ik)         :: ionization_l_max   = 0                     ! Order of harmonic expansion of ionization yields
    real(rk)            :: ionization_e_peak  = 0                     ! Peak laser electric field used for ionization calculation
    real(rk)            :: ionization_duration= 0                     ! Duration of ionization pulse; used to calculate rates from yields
    character(len=clen) :: ionization_comment = ' '
    !
    integer(ik)         :: recombination_l_max   = 0                  ! Order of harmonic expansion of recombination dipoles
    integer(ik)         :: recombination_n_ekin  = 0                  ! Number of different kinetic energies in the file
    character(len=clen) :: recombination_comment = ' '
    !
    integer(ik)         :: euler_npoints      =  0                    ! Number of discretization per Euler angle in molecular 
                                                                      ! integration. The total number of integration points is
                                                                      ! 0.5_rk*euler_npoints**3 (we only use half of the number
                                                                      ! points for the beta angle, which goes from 0 to Pi, instead
                                                                      ! of 0 to 2Pi for alpha and gamma
                                                                      ! Note that it is pointles to exceed
                                                                      !   euler_npoints=3+ionization_l_max+recombination_l_max
                                                                      ! since integration becomes exact at this point.
                                                                      ! Setting euler_npoints to a 0 will use the exact integration
                                                                      ! For convenience, euler_npoints has to be even.
    logical             :: use_dyson_phase    = .false.               ! Apply phase factor taken from the Dyson orbital to ionization
                                                                      ! amplitudes
    real(rk)            :: phase_origin(3)    = 0._rk                 ! Origin used to determine the phase factor
    real(rk)            :: phase_radius       = 10._rk                ! Radius of the sphere used to determine the phase factor
    character(len=clen) :: dyson_file         = 'dyson.dat'           ! Data file (gamess-like .chk format) containing the Dyson orbital
                                                                      ! used in ionization calculations. Only used to determine the
                                                                      ! ionization phase.
    integer(ik)         :: dyson_index        = 2                     ! Index of the Dyson orbital in the input file; most likely 2
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type(sfa_data)           :: sfa_par
    complex(rk), allocatable :: ion_clm(:,:)                          ! (m,l) Spherical harmonics expansion of ionization yields
    real(rk), allocatable    :: rec_p(:)                              ! Asymptotic recombination momentum, atomic units.
                                                                      ! The input file has kinetic energy in eV; the conversion is done
                                                                      ! while reading the file. The energies/momenta must be in the 
                                                                      ! ascending order; there is no checking.
    complex(rk), allocatable :: rec_clm(:,:,:,:)                      ! (m,l,ic,ie) Spherical harmonics expansion of ionization yields;
                                                                      ! The third component is 1/2/3, standing for X/Y/Z part of the dipole
                                                                      ! The fourth component is energy grid
    type(gam_structure)      :: gam_dyson                             ! Structure/MO descriptor for the Dyson orbital
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /chiral_data/ verbose, &
                           omega, e0x, e0y, ip, &
                           phr_min, phr_max, phr_step, phi_min, phi_max, &
                           n_return, sfa_magnetic, trajectories_only, &
                           ionization_file, recombination_file, &
                           euler_npoints, ion_dipole, ion_alpha, &
                           use_dyson_phase, dyson_file, dyson_index, phase_origin, phase_radius
    !
    namelist /ionization_fit/ ionization_l_max, &
                              ionization_e_peak, &
                              ionization_duration, &
                              ionization_comment
    !
    namelist /recombination_fit/ recombination_l_max, &
                                 recombination_n_ekin, &
                                 recombination_comment
    !
    !  ==== End of global data ====
    !
    contains
    !
    subroutine read_harmonic_fit(iu,clm,l_max)
      integer(ik), intent(in)  :: iu   ! I/O unit containing the data
      integer(ik), intent(in)  :: l_max
      complex(rk), intent(out) :: clm(-l_max:l_max,0:l_max)
      !
      integer(ik) :: mv, lv, m, l
      real(rk)    :: re, im
      !
      clm = 0
      read_l: do l=0,l_max
        read_m: do m=-l,l
          read (iu,*) lv, mv, re, im
          if (lv/=l .or. mv/=m) stop 'read_harmonic_fit - format mismatch'
          clm(m,l) = cmplx(re,im,kind=rk)
        end do read_m
      end do read_l
    end subroutine read_harmonic_fit
    !
    subroutine read_ion
      open (iu_input,form='formatted',action='read',position='rewind',file=trim(ionization_file))
      read (iu_input,nml=ionization_fit)
      write (out,"(/'Reading ionization fit data: ',a)") trim(ionization_file)
      write (out,nml=ionization_fit)
      write (out,"(80('='))")
      allocate (ion_clm(-ionization_l_max:ionization_l_max,0:ionization_l_max))
      call read_harmonic_fit(iu_input,ion_clm,ionization_l_max)
      close (iu_input)
    end subroutine read_ion
    !
    subroutine read_rec
      integer(ik) :: ien, ic
      real(rk)    :: en ! Asymptotic kinetic energy, in eV
      !
      open (iu_input,form='formatted',action='read',position='rewind',file=trim(recombination_file))
      read (iu_input,nml=recombination_fit)
      write (out,"(/'Reading recombination fit data: ',a)") trim(recombination_file)
      write (out,nml=recombination_fit)
      write (out,"(80('='))")
      allocate (rec_p(recombination_n_ekin))
      allocate (rec_clm(-recombination_l_max:recombination_l_max,0:recombination_l_max,3,recombination_n_ekin))
      read_energy_points: do ien=1,recombination_n_ekin
        read (iu_input,*) en
        rec_p(ien) = sqrt(2._rk*en/h2ev)
        read_dipole_components: do ic=1,3
          call read_harmonic_fit(iu_input,rec_clm(:,:,ic,ien),recombination_l_max)
        end do read_dipole_components
      end do read_energy_points
      close (iu_input)
    end subroutine read_rec
    !
    function interpolate_ionization(ei) result(ai)
      real(rk), intent(in) :: ei(3) ! Electric field direction (not normalized)
      complex(rk)          :: ai    ! Ionization amplitude
      !
      integer(ik) :: l, m
      complex(rk) :: ylm(-ionization_l_max:ionization_l_max,0:ionization_l_max)
      real(rk)    :: coord(3,1,1,1)
      complex(rk) :: moval(1,1,1,1)
      ! 
      !  Calculate the yield in the specified direction
      !
      call MathAllYLM(ionization_l_max,ei,ylm)
      ai = 0._rk
      contract_l: do l=0,ionization_l_max
        contract_m: do m=-l,l
          ai = ai + ion_clm(m,l) * ylm(m,l)
          ! ai = ai + ion_clm(m,l) * MathYLM(l,m,ei)
        end do contract_m
      end do contract_l
      !
      !  Convert to the rate
      !
      ai = (0.5_rk/ionization_duration) * ai
      !
      if (.not.use_dyson_phase) return
      !
      !  We have to add a +/-1 phase factor, depending on the phase of the Dyson orbital
      !  We need to sample Dyson orbital in the direction electron will be leaving in,
      !  far enough from the origin to make it more or less into the asymptotic region
      !
      coord(:,1,1,1) = phase_origin + electron_q*ei*phase_radius/(abs(electron_q)*sqrt(sum(ei**2)))
      call gamess_load_orbitals(mos=(/dyson_index/),dst=(/1_ik/),coord=coord,grid=moval,structure=gam_dyson)
      ! write (out,"('At point ',3f12.5,' MO = ',2g14.6)") coord, moval
      if (real(moval(1,1,1,1))<0._rk) ai = -ai
    end function interpolate_ionization
    !
    function interpolate_recombination(pr) result(ar)
      real(rk), intent(in) :: pr(3) ! Kinetic momentum at recombination
      complex(rk)          :: ar(3) ! Recombination dipole
      !
      integer(ik)   :: l, m
      real(rk)      :: p       ! Target momentum
      integer(ik)   :: r1, r2  ! Indices of the energy grid we interpolate between
      real(rk)      :: c1, c2  ! Interpolation coefficients
      logical, save :: warn_min = .true.
      logical, save :: warn_max = .true.
      complex(rk)   :: ylm(-recombination_l_max:recombination_l_max,0:recombination_l_max)
      !
      !  Choose two points on the recombination energy grid we'll be interpolating between
      !
      p = sqrt(sum(pr**2))
      if (p<rec_p(1)) then
        if (warn_min) then
          warn_min = .false.
          write (out,"(/'WARNING: Requested momentum (',g16.7,') is below the grid minimum (',g16.7,')')") p, rec_p(1)
          write (out,"( 'WARNING: Substituting grid minimum. Further warnings are suppressed.'/)")
        end if
        r1 = 1 ; r2 = 1 ;
        c1 = 1.0_rk ; c2 = 0.0_rk
      else if (p>rec_p(size(rec_p))) then
        if (warn_max) then
          warn_max = .false.
          write (out,"(/'WARNING: Requested momentum (',g16.7,') is above the grid maximum (',g16.7,')')") p, rec_p(size(rec_p))
          write (out,"( 'WARNING: Substituting grid maximum. Further warnings are suppressed.'/)")
        end if
        r1 = size(rec_p) ; r2 = size(rec_p)
        c1 = 1.0_rk ; c2 = 0.0_rk
      else
        scan_momenta: do r1=1,size(rec_p)-1
          if (rec_p(r1+1)>=p) exit scan_momenta
        end do scan_momenta
        r2 = r1 + 1
        if (r2>size(rec_p)) stop 'interpolate_recombination - logic error'
        c2 = (p-rec_p(r1))/(rec_p(r2)-rec_p(r1))
        c1 = (rec_p(r2)-p)/(rec_p(r2)-rec_p(r1))
      end if
      ! 
      !  Calculate the yield in the specified direction
      !
      ar = 0._rk
      call MathAllYLM(recombination_l_max,pr,ylm)
      contract_l: do l=0,recombination_l_max
        contract_m: do m=-l,l
          ar = ar + (c1*rec_clm(m,l,:,r1)+c2*rec_clm(m,l,:,r2)) * ylm(m,l)
          ! ar = ar + (c1*rec_clm(m,l,:,r1)+c2*rec_clm(m,l,:,r2)) * MathYLM(l,m,pr)
        end do contract_m
      end do contract_l
    end function interpolate_recombination
    !
    !  Calculate exponential ionization factor for the current field intensity and initial momentum
    !  We'll go with the tunneling exponent; since this factor is common to the entire molecular
    !  integral, it should be good enough.
    !
    function ion_factor(ei,pi) result(c_ion)
      real(rk), intent(in) :: ei(3)  ! Electric field vector, atomic units
      real(rk), intent(in) :: pi(3)  ! Initial momentum, atomic units
      real(rk)             :: c_ion  ! Exponential scale factor
      !
      real(rk) :: ef
      real(rk) :: exp_ref, exp_curr
      !
      ef     = sqrt(sum(ei**2))
      if (ef<=0._rk) then
        c_ion = 0._rk
        return
      end if
      !
      exp_curr = -(1._rk/3._rk) * (2*ip+sum(pi**2))**1.5_rk / ef
      exp_ref  = -(1._rk/3._rk) * (2*ip)**1.5_rk / ionization_e_peak
      !
      c_ion = exp(exp_curr-exp_ref)
    end function ion_factor
    !
    !  Additional ionization momentum contribution induced by the molecular potential
    !  We are using perturbation-theory contribution, calculated along the SFA ougoing
    !  trajectory. 
    !
    function ionization_momentum_correction(ei) result (dp)
      real(rk), intent(in)     :: ei(3)  ! Electric field at the time of ionization, molecular frame
      real(rk)                 :: dp(3)  ! Potential-induced correction to ionization momentum
      !
      real(rk) :: r0(3)             ! Direction of the ionizing field
      real(rk) :: f2                ! Instantanous intensity of the field
      real(rk) :: instant_dipole(3) ! Adiabatic dipole moment - permanent + induced
      !
      instant_dipole = ion_dipole + matmul(ion_alpha,ei)
      !
      f2 = sum(ei**2)
      r0 = ei / sqrt(f2)
      dp = 3._rk*sum(instant_dipole*r0)*r0 - instant_dipole
      dp = dp * ((3._rk*pi)/(8._rk*sqrt(2._rk))) * electron_q**3 * f2 / ip**2.5_rk
      !
      !  The expression above gives the moment induced by the ion's dipole while the electron is leaving;
      !  we now need to flip the sign, to give us the extra moment we must have at ionization time to 
      !  compensate for the shift
      !
      dp = -dp
    end function ionization_momentum_correction
    !
    !  Calculation of molecular part of the SFA amplitudes for two stereoisomers
    !
    subroutine molecular_amplitude(ei,pi,pr,amol1,amol2)
      real(rk), intent(in)     :: ei(3)    ! Electric field at the time of ionization, lab frame
      real(rk), intent(in)     :: pi(3)    ! Kinetic momentum at the time of ionization, lab frame
      real(rk), intent(in)     :: pr(3)    ! Kinetic momentum at the time of recombination, lab frame
      complex(rk), intent(out) :: amol1(3) ! Molecular amplitude averaged over the ensemble; original isomer
      complex(rk), intent(out) :: amol2(3) ! Molecular amplitude averaged over the ensemble; mirror isomer
      !
      integer(ik) :: ialp, ibet, igam
      real(rk)    :: alp, bet, gam     ! Euler angles
      real(rk)    :: jacobian          ! Jacobian of the angular averaging integral
      real(rk)    :: step
      real(rk)    :: c_ion             ! Ionization scaling factor (Keldysh exponent)
      real(rk)    :: c_ion_mir         ! Ionization scaling factor, mirror isomer
      real(rk)    :: u(3,3)
      real(rk)    :: pii(3)            ! Lab frame ionization direction, reflected in the lab (xy) axis
      real(rk)    :: pri(3)            ! Lab frame recombination direction, reflected in the lab (xy) axis
      real(rk)    :: ei_mol(3)         ! Ionizing field direction, molecular frame
      real(rk)    :: pr_mol(3)         ! Recombination direction, molecular frame
      real(rk)    :: pi_mol(3)         ! Ionization direction, molecular frame
      real(rk)    :: tot_norm          ! Sanity test of the Jacobian
      complex(rk) :: tot_ion           ! Integrated ionization amplitude
      complex(rk) :: tot_ion_mir       ! Integrated ionization amplitude, mirror isomer
      complex(rk) :: a_ion             ! Ionization amplitude, before scaling
      complex(rk) :: a_rec(3)          ! Recombination amplitude
      complex(rk) :: a_tot(3)          ! Total amplitude, before scaling and normalization
      complex(rk) :: aloc1(3), aloc2(3)! Intel Fortran 9.1 and 10.1 have a bug with handling OpenMP reductions
                                       ! passed as subroutine arguments. Hence a local copy here
      real(rk)    :: ion_p0(3)         ! Dipole-potential correction to ionizaton momentum
      !
      ! c_ion = ion_factor(ei,pi) ! Now depends on the molecular orientation
      step  = twopi/euler_npoints
      pii   = (/ pi(1), pi(2), -pi(3) /)
      pri   = (/ pr(1), pr(2), -pr(3) /)
      !
      if (verbose>=2) then
        write (out,"('        Lab frame ionizing field: ',3f15.6)") ei
        write (out,"('      Lab frame initial momentum: ',3f15.6)") pi
        write (out,"('Lab frame recombination momentum: ',3f15.6)") pr
        write (out,"('    Euler angle integration step: ',g15.6)") step
      ! write (out,"('  Ionization rate scaling factor: ',g15.6)") c_ion
      end if
      !
      tot_ion     = 0
      tot_ion_mir = 0
      tot_norm    = 0
      aloc1       = 0
      aloc2       = 0
      if (mod(euler_npoints,2)/=0) stop 'chiral%molecular_amplitude - euler_npoints must be even!'
      !
      !$omp parallel do default(none) &
      !$omp&  reduction(+:aloc1,aloc2,tot_ion,tot_ion_mir,tot_norm) &
      !$omp&  private(ialp,ibet,igam,alp,bet,gam,u,ei_mol,a_ion,jacobian) &
      !$omp&  private(c_ion,c_ion_mir,pi_mol) &
      !$omp&  private(pr_mol,a_rec,a_tot,ion_p0) &
      !$omp&  shared(step,euler_npoints,ei,pi,pii,pr,pri)
      euler_gamma: do igam=1,euler_npoints
        gam = step * (igam-0.5_rk)
        euler_beta: do ibet=1,euler_npoints/2
          bet = step * (ibet-0.5_rk)
          jacobian = sin(bet)
          if (jacobian<0._rk) stop 'chiral%molecular_amplitude - jacobian became negative'
          euler_alpha: do ialp=1,euler_npoints
            alp = step * (ialp-0.5_rk)
            call MathRotationMatrix(euler_angles=(/alp,bet,gam/),mat=u)
            !
            !  Sanity test
            !
            tot_norm = tot_norm + jacobian
            !
            !  Ionization scaling factor is now a function of the isomer, but only
            !  if ion_dipole is not zero.
            !
            ei_mol      = matmul(u,ei)
            ion_p0      = ionization_momentum_correction(ei_mol)
            pi_mol      = matmul(u,pi)
            c_ion       = ion_factor(ei_mol,pi_mol+ion_p0)
            ! ???
            ! write (out,"('pi normal = ',3f14.7,' corr = ',3f14.7,' c_ion = ',g14.6)") pi_mol, ion_p0, c_ion
            pi_mol      = matmul(u,pii)
            c_ion_mir   = ion_factor(ei_mol,pi_mol+ion_p0)
            ! write (out,"('   mirror = ',3f14.7,'    . = ',3f14.7,'     . = ',g14.6)") pi_mol, ion_p0, c_ion_mir
            a_ion       = interpolate_ionization(ei_mol)
            tot_ion     = tot_ion     + a_ion * jacobian * c_ion
            tot_ion_mir = tot_ion_mir + a_ion * jacobian * c_ion_mir
            !
            !  Recombination: original isomer first
            !
            pr_mol = matmul(u,pr)
            a_rec  = interpolate_recombination(pr_mol)
            a_tot  = matmul(transpose(u),a_rec) * a_ion * c_ion
            aloc1  = aloc1 + a_tot * jacobian
            !
            !  Now the mirror isomer. 
            !
            pr_mol = matmul(u,pri)
            a_rec  = interpolate_recombination(pr_mol)
            a_tot  = matmul(transpose(u),a_rec) * a_ion * c_ion_mir
            a_tot(3) = -a_tot(3)   ! Reflection in the XY plane
            aloc2   = aloc2 + a_tot * jacobian
          end do euler_alpha
        end do euler_beta
      end do euler_gamma
      !$omp end parallel do
      !
      !  Scale by the overall ionization factor and integral prefactors
      !
      amol1       = aloc1       * step**3 / (2._rk*twopi**2)
      amol2       = aloc2       * step**3 / (2._rk*twopi**2)
      tot_ion     = tot_ion     * step**3 / (2._rk*twopi**2)
      tot_ion_mir = tot_ion_mir * step**3 / (2._rk*twopi**2)
      tot_norm    = tot_norm    * step**3 / (2._rk*twopi**2)
      if (verbose>=1) then
        write (out,"('             Molecular amplitude: ',3(g15.6,1x,g15.6,2x))") amol1
        write (out,"('    Molecular amplitude (mirror): ',3(g15.6,1x,g15.6,2x))") amol2
        write (out,"(' Ionization amplitude (included): ',  g15.6,1x,g15.6,2x )") tot_ion
        write (out,"('   Ionization amplitude (mirror): ',  g15.6,1x,g15.6,2x )") tot_ion_mir
        write (out,"('           Norm test (must be 1): ',  g15.6             )") tot_norm
      end if
      if (abs(tot_norm-1.0_rk)>=1e-2) stop 'chiral%molecular_amplitude - bad norm in angular integral'
    end subroutine molecular_amplitude
    !
    !  Problem driver
    !
    subroutine run_chiral
      integer(ik)        :: info
      complex(rk)        :: phr     ! Current recollision phase
      complex(rk)        :: phi     ! Current ionization phase
      complex(rk)        :: phi_old ! Previosly found ionization phase
      complex(rk)        :: phs     ! Real part of the ionization phase
      complex(rk)        :: p0(3), pi(3), pr(3)
      complex(rk)        :: et
      complex(rk)        :: phase
      complex(rk)        :: aprop, amol1(3), amol2(3)
      real(rk)           :: ei(3), er(3)
      real(rk)           :: d1(3), d2(3), rat
      real(rk)           :: ad1(3), ad2(3), arat
      logical            :: success
      integer(ik)        :: i_return
      real(rk)           :: dummy
      !
      call TimerStart('Chiral')
      call accuracyInitialize
      dummy = MathFactorial(80_ik)
      !
      write (out,"(/'Version with the (hopefully) correct Jacobian')")
      write (out,"( 'Version with the (hopefully) correct ionization phase'/)")
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=chiral_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=chiral_data)
      write (out,"()")
      !
      call read_ion
      !
      call read_rec
      !
      write (out,"(/'Ion dipole moment = ',3f16.12)") ion_dipole
      write (out,"( 'Ion polarizability tensor:')")
      write (out,"(5x,3f16.10)") transpose(ion_alpha)
      write (out,"()")
      !
      if (use_dyson_phase) then
        write (out,"('Loading Dyson orbital from: ',a)") trim(dyson_file)
        call gamess_load_orbitals(file=dyson_file,structure=gam_dyson)
      end if
      !
      if (euler_npoints<=0) then
        euler_npoints = 3 + ionization_l_max + recombination_l_max
        euler_npoints = euler_npoints + mod(euler_npoints,2)
      end if
      !*HACK!
      ! euler_npoints = 2*euler_npoints
      ! write (out,"('HACK ALERT: DOUBLING euler_npoints')")
      !*HACK!
      write (out,"('Using ',i4,' points per Euler angle.')") euler_npoints
      write (out,"('Exact integration requires ',i4,' points per angle.')") 3+ionization_l_max+recombination_l_max
      write (out,"('Total number of grid points = ',i12)") (euler_npoints**3)/2
      !
      call sfa_init(sfa_par,omega=omega,ex=e0x,ey=e0y,ip=ip,magnetic=sfa_magnetic)
      !
      write (out,"(/'  The ''#@'' line gives the result of the stationary-phase analysis, including:')")
      write (out,"( '   2 Re(Ion.p)  = Real part of the stationary ionization phase (aka phase of birth), Radian')")
      write (out,"( '   3 Im(Ion.p)  = Imaginary part of the stationary ionization phase, Radian')")
      write (out,"( '   4 Ion. Px    = X component of the momentum at the phase of birth [Re(Ion.p)], atomic units')")
      write (out,"( '   5 Ion. Py    = Y component of the same')")
      write (out,"( '   6 Ion. Pz    = Z component of the same')")
      write (out,"( '   7 Rec. phas  = Stationary recombination phase, Radian')")
      write (out,"( '   8 Rec. Px    = X component of the electron momentum at the time of recombination, atomic units')")
      write (out,"( '   9 Rec. Py    = Y component of the same')")
      write (out,"( '  10 Rec. Pz    = Z component of the same')")
      write (out,"( '  11 Re(phase)  = Propagation phase between times of birth [Re(ion.p)/omega] and recollision')")
      write (out,"( '                  This phase includes the IP*(tr-ts) term, but not the ionization phase.')")
      write (out,"( '  12 Im(phase)  = Must be zero')")
      write (out,"( '          The ionization and recombination phases are for the laser field at the origin in the form:')")
      write (out,"( '             Ex = E0x * Cos(p)')")
      write (out,"( '             Ey = E0y * Sin(p)')")
      write (out,"( '             Ez = 0')")
      !
      write (out,"(/'  The ''#+'' line gives recombination dipole for normal and mirror')")
      write (out,"( '  species at the time of recombination, including the phase factor:')")
      write (out,"( '   2 Ion.ph.    = Real part of the stationary ionization phase (aka phase of birth), Radian')")
      write (out,"( '   3 Rec.ph.    = Stationary recombination phase, Radian')")
      write (out,"( '   4 Energy,H   = Emission energy (IP+0.5*rec.p**2), Hartree')")
      write (out,"( '   5 difference = (d(mir),x**2+d(mir),y**2)/(d(reg),x**2+d(reg),y**2)-1.0')")
      write (out,"( '                  Note that only the X and Y components of radiating dipole survive propagation')")
      write (out,"( '   6 d(reg),x   = X component of radiating dipole [2*Re(d)], normal isomer')")
      write (out,"( '   7 d(reg),y   = Y component ')")
      write (out,"( '   8 d(reg),z   = Z component')")
      write (out,"( '   9 d(mir),x   = X component of radiating dipole, mirror isomer')")
      write (out,"( '  10 d(mir),y   = Y component')")
      write (out,"( '  11 d(mir),z   = Z component')")
      !
      write (out,"(/'  The ''#-'' line omits all the phase factors, for clearer comparisons:')")
      write (out,"( '   2 Ion.ph.    = Real part of the stationary ionization phase (aka phase of birth), Radian')")
      write (out,"( '   3 Rec.ph.    = Stationary recombination phase, Radian')")
      write (out,"( '   4 Energy,H   = Emission energy (IP+0.5*rec.p**2), Hartree')")
      write (out,"( '   5 difference = (d(mir),x**2+d(mir),y**2)/(d(reg),x**2+d(reg),y**2)-1.0')")
      write (out,"( '   6 d(reg),x   = X component of radiating dipole magnitude [Abs(d)], normal isomer')")
      write (out,"( '   7 d(reg),y   = Y component')")
      write (out,"( '   8 d(reg),z   = Z component')")
      write (out,"( '   9 d(mir),x   = X component of radiating dipole magnitude, mirror isomer')")
      write (out,"( '  10 d(mir),y   = Y component')")
      write (out,"( '  11 d(mir),z   = Z component')")
      write (out,"()")
      !
      write (out,"('#@',a10,1x,a10,1x,3(1x,a12),2x,a10,1x,3(1x,a12),1x,2(1x,a12))") &
          ' Re(Ion.p.) ', ' Im(Ion.p.)',  ' Ion. Px ', ' Ion. Py ', ' Ion. Pz ', &
                          ' Rec. phase ', ' Rec. Px ', ' Rec. Py ', ' Rec. Pz ', &
          ' Re(phase) ', ' Im(phase) '
      write (out,"('#+',a10,1x,a10,1x,a10,1x,a12,2(1x,3(1x,a14)))") &
          ' Ion.ph. ', ' Rec.ph. ', ' Energy,H ', ' difference ', &
          ' d(reg),x ', ' d(reg),y ', ' d(reg),z ', ' d(mir),x ', ' d(mir),y ', ' d(mir),z '
      write (out,"('#-',a10,1x,a10,1x,a10,1x,a12,2(1x,3(1x,a14)))") &
          ' Ion.ph. ', ' Rec.ph. ', ' Energy,H ', ' difference ', &
          '|d(reg),x|', '|d(reg),y|', '|d(reg),z|', '|d(mir),x|', '|d(mir),y|', '|d(mir),z|'
      !
      !  Recombination-point scan moves forwards in time
      !
      phr = phr_min
      scan_recollision_phases: do while(real(phr,kind=rk)<=phr_max)
        !
        !  Ionization-point search proceeds backwards in time,
        !  starting with the largest point within the search windown
        !
        i_return = 0
        phi_old  = phi_max
        scan_ionization_phases: do 
          success = sfa_stationary_t1(sfa_par,ph2=phr,ph1_max=phi_old,ph1=phi)
          !
          !  Bail out if the final search phase is outside the ionization 
          !  window, even if the search was successful otherwise. We must
          !  bail out if optimization gives a too-early point, too: otherise,
          !  we'll hit an infinite loop here.
          !
          if (real(phi,kind=rk)<phi_min) exit scan_ionization_phases
          if (real(phi,kind=rk)>phi_max) exit scan_ionization_phases
          !
          if (success) then
            ! 
            !  Make sure we got the stationary ionization time correctly
            !
            p0 = sfa_stationary_momentum(sfa_par,ph2=phr,ph1=phi)
            pi = sfa_kvec(sfa_par,p=p0,ph=phi) * sfa_par%s_p0
            et = ip + 0.5_rk * sum(pi**2)
            if (abs(et)>1e-5_rk) then
              stop 'chral%run_chiral - Bad solution for ionization phase!'
            end if
            !
            !  Extract time of exit from under the barrier, and construct the
            !  corresponding trajectory.
            !
            phs   = real(phi,kind=rk)
            p0    = sfa_stationary_momentum(sfa_par,ph2=phr,ph1=phs)
            pi    = sfa_kvec(sfa_par,p=p0,ph=phs) * sfa_par%s_p0
            pr    = sfa_kvec(sfa_par,p=p0,ph=phr) * sfa_par%s_p0
            ei    = real(sfa_efield(sfa_par,ph=phs),kind=rk) * sfa_par%s_e0
            er    = real(sfa_efield(sfa_par,ph=phr),kind=rk) * sfa_par%s_e0
            phase = sfa_phase(sfa_par,p=p0,ph1=phs,ph2=phr)
            write (out,"('@ ',f10.6,1x,f10.6,1x,3(1x,f12.6),2x,f10.6,1x,3(1x,f12.6),1x,2(1x,f12.6))") &
                   real(phs,kind=rk), aimag(phi), real(pi,kind=rk), real(phr,kind=rk), real(pr,kind=rk), phase
            !
            if (.not.trajectories_only) then
              aprop = (twopi/(sfa_par%s_t0*(phr-phs)))**1.5_rk * exp((0,1)*phase)
              call molecular_amplitude(ei,real(pi,kind=rk),real(pr,kind=rk),amol1,amol2)
              if (verbose>=1) then
                write (out,"('           Propagation amplitude: ',2g15.6)") aprop
              end if
              d1    = 2.0_rk * real(aprop*amol1/sqrt((0._rk,1._rk)))
              d2    = 2.0_rk * real(aprop*amol2/sqrt((0._rk,1._rk)))
              rat   = sum(d2(1:2)**2)/sum(d1(1:2)**2) - 1.0_rk
              write (out,"('+ ',f10.6,1x,f10.6,1x,f10.6,1x,e12.5,2(1x,3(1x,g14.6)))") &
                     real(phs,kind=rk), real(phr,kind=rk), ip+0.5_rk*real(sum(pr**2),kind=rk), rat, d1, d2
              ad1   = abs(aprop*amol1/sqrt((0._rk,1._rk)))
              ad2   = abs(aprop*amol2/sqrt((0._rk,1._rk)))
              arat  = sum(ad2(1:2)**2)/sum(ad1(1:2)**2) - 1.0_rk
              write (out,"('- ',f10.6,1x,f10.6,1x,f10.6,1x,e12.5,2(1x,3(1x,g14.6)))") &
                     real(phs,kind=rk), real(phr,kind=rk), ip+0.5_rk*real(sum(pr**2),kind=rk), arat, ad1, ad2
            end if
            !
            !  Count a successful ionization point; if we reached the limit on
            !  the number of ionization points desired, bail out.
            !
            i_return = i_return + 1
            if (i_return>=n_return .and. n_return>0) exit scan_ionization_phases
          end if
          phi_old = phi
        end do scan_ionization_phases
        phr = phr + phr_step
      end do scan_recollision_phases
      !
      call TimerStop('Chiral')
      call TimerReport
    end subroutine run_chiral

  end module chiral
!
  program driver
    use chiral
    
    call run_chiral
  end program driver
