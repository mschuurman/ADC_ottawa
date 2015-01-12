!
!  Common routines and definitions, shared by the various versions of the 
!  "dyson" package (dyson_1chan.f90 and dyson_nchan.f90). Since these routines
!  and variables are a tightly coupled part of the parent module, everything
!  here is public.
!
  module dyson_tools
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use eikonal_tools
    use ecp_gamess
    use caps
    use math
    use flux
    use fftw

    implicit none

    public
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_naturals = 200                 ! Max number of natural orbitals allowed
    integer(ik), parameter :: unit_dump    = 34                  ! A more or less random unit
    integer(ik), parameter :: unit_ampd    = 35                  ! Output unit for amplitudes of dyson (source) states and ampitude b(t)
    integer(ik), parameter :: unit_chi     = 36                  ! Output unit for absorbed populations in each channel
    integer(ik), parameter :: unit_abs     = 37                  ! Output unit for absorbed populations in each channel
    integer(ik), parameter :: unit_rec     = 38                  ! Output unit for recombination matrix elements
    integer(ik), parameter :: unit_detail  = 39                  ! A more or less random unit
    integer(ik), parameter :: unit_efield  = 40                  ! A more or less random unit
    integer(ik), parameter :: unit_ckpt    = 41                  ! A more or less random unit
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose         = 2                     ! Verbosity level
    integer(ik)         :: n_points(3)     = (/ 100, 100, 100 /)   ! Number of sampling points along each direction
    real(rk)            :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    !
    real(rk)            :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                   ! and exchange potential(s)
    real(rk)            :: sor_rate        = 1.95_rk               ! Over-relaxation rate in solving Poisson equations
    character(len=clen) :: v_xc            = ' '                   ! Exchange-correlation potential to use in construction of the
                                                                   ! eikonal functions. Possible choices are:
                                                                   ! ' ' or 'none' - no exchange potential
                                                                   ! 'Slater' - Dirac-Slater exchange
                                                                   ! 'SVWN'   - Slater exchange + VWN local correlation
    character(len=clen) :: ecp_file        = ' '                   ! Name of a file containing ECP
    real(rk)            :: ecp_eps_min     = 1e-2_rk               ! Small shift value cut-off (absolute) in ECP
    real(rk)            :: ecp_eps_max     = 1e3_rk                ! Large shift value cut-off (positive) in ECP
    real(rk)            :: ecp_eps_grid    = 0.2_rk                ! Characteristic grid spacing cut-off in ECP
    character(len=clen) :: ecp_report      = ' '                   ! File to report ECP level-shift projectors to
    !
    integer(ik)         :: nspin           = 2                     ! Number of (degenerate) spin components in the continuum.
                                                                   ! Note that the code in its present form can't explicitly handle 
                                                                   ! continuum channels with different spins; as implemented, our
                                                                   ! equations are correct for a) all channels being the same spin;
                                                                   ! or b) each channel being doubly degenerate in spin.
                                                                   ! Use 2 for the overall singlet state, and 1 for the doublet
    real(ark)           :: euler_angles(3) = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
                                                                   ! Euler rotation angles - used to reorient the molecule.
                                                                   ! See MathRotationMatrix in math.f90 for definitions,
                                                                   ! and for some useful special cases.
    !
    real(rk)            :: tMin    = 0._rk                         ! Initial time, atomic units
    real(rk)            :: tMax    = 200._rk                       ! Final time, atomic units
    real(rk)            :: dt                                      ! Time step, atomic units

    character(len=clen) :: ElectricFieldShape = ' '                ! Electric field shape - see ElectricField() below for the possibilities
                                                                   ! Note that the interpretation of the entries below will be affected by
                                                                   ! the choice of ElectricFieldShape.
    real(rk)            :: Elaser    = 0.1_rk                      ! Maximum electric field strength
    real(rk)            :: Edir(3)   = (/ 0._rk, 0._rk, 1._rk /)   ! Direction of laser polarization
    real(rk)            :: Wlaser    = 2._rk                       ! Frequency of the laser, atomic units
    real(rk)            :: TlaserOn  = 0._rk                       ! Time for the pulse turn-on, in cycles 
    real(rk)            :: TlaserOff = 1._rk                       ! Time for the pulse turn-off, in cycles
    real(rk)            :: RampTime  = 1._rk                       ! Time for ramp on and off of the flat-top pulse, in cycles
    real(rk)            :: FlatTime  = 1._rk                       ! Time for the flat part of the flat-top pulse, in cycles
    real(rk)            :: Plaser    = 0._rk                       ! Laser phase at t=0, fractions of the cycle
    real(rk)            :: FWHM    = 1000._rk                      ! Full width at half maximum
    real(rk)            :: Tlaser  = 2.000                         ! Center of pulse in time
    !
    real(rk)            :: abs_bound_kmin   = 1.5_rk               ! Kmin parameter for absorbing boundary
    character(len=clen) :: abs_bound_type  = ' '                   ! Absorbing boundary edges:
                                                                   ! ' '   = no abs. boundaries
                                                                   ! 'all' = abs. boundaries on all edges
                                                                   ! 'xy' = abs. boundaries on x and y edges
                                                                   ! 'xz' = abs. boundaries on x and z edges
                                                                   ! 'yz' = abs. boundaries on y and z edges

    character(len=clen) :: harmonic_operator = 'dipole'            ! One-electron operator for the emission matrix elements. Can be:
                                                                   !  'none'         - Do not calculate the emission (minor speed gain)
                                                                   !  'dipole'       - Cartesian dipole operator. "Exchange"
                                                                   !                   corrections will be included, as long
                                                                   !                   as the non-orthogonal Dyson term is
                                                                   !                   present.
                                                                   !  'acceleration' - Gradient of the potential. External
                                                                   !                   potential (natural_file) must be
                                                                   !                   supplied for this form to work.
                                                                   !   Note that the acceleration form is not internally consistent in
                                                                   !   the present formulation, and is not recommended.
    character(len=clen) :: potential_file  = ' '                   ! File for the potential dump - can be HUGE
    character(len=clen) :: output_prefix   = ' '                   ! Prefix used for output files
    character(len=clen) :: detailed_output = ' '                   ! If non-emply, dump the detailed output on each time step
    integer(ik)         :: report_each     = 100                   ! Reduce report frequency in the main output
    integer(ik)         :: odx_plot_each   = -1                    ! Reduce plotting frequency to each (plot_each)-th step
                                                                   ! Set to -1 to disable OpenDX plotting during the simulation,
                                                                   !     but still produce a prot at the final time
                                                                   ! Set to -2 to disable OpenDX plotting altogether
    logical             :: odx_plot_momentum = .false.             ! Enabling this one could cause segfaults for large arrays,
                                                                   ! due to a bug in FFTW.
    character(len=clen) :: efield_preview  = ' '                   ! File for the e-field "preview" for this simulation
                                                                   ! Empty string suppresses the output
    character(len=clen) :: checkpoint_out  = "('ckpt',i10.10,'.dat')"
                                                                   ! Format string used to produce the checkpoint file name. The
                                                                   ! current time step will be passed as a parameter.
    character(len=clen) :: checkpoint_in   = ' '                   ! Name of the checkpoint file to load
    character(len=clen) :: potentials_out  = ' '                   ! Base name for checkpointing Hartree and transition potentials
                                                                   ! Actual names will be formed by appending unique strings to it.
                                                                   ! See checkpoint_ion_potentials() in dyson_nchan.f90 for the
                                                                   ! exact conventions. 
                                                                   ! Empty string disables checkpoints altogether.
    character(len=clen) :: potentials_in   = ' '                   ! Reload checkpointed potentials. Absolutely no sanitu checking!
    integer(ik)         :: checkpoint_each = 1000                  ! Generate a binary checkpoint file after these many time steps,
                                                                   ! and at the end of the run. Set to -1 to disable checkpoints
                                                                   ! altogether. Setting checkpoint_each to zero will produce a
                                                                   ! single checkpoint at the normal completion.
    logical             :: oversample_orbitals = .true.            ! Oversample orbitals on grid; should give smoother grid variatons
    
    !
    !  Analytical continuation of photoelectron momentum, a-la "t-SURFF" method 
    !  of Tao and Scrinzi (arXiv:1109.4053v1). This choice is only meaniningful
    !  if complex absorbing potential is active as well. This option is noticeably
    !  expensive; only use it of you do need photoelectron spectrum.
    !
    logical                :: continue_flux    = .false.
    logical                :: flux_restart     = .false.           ! Try to restart flux simulation. NOT TESTED IN ANY WAY
    real(rk)               :: flux_guard(3)    = 10._rk
                                                                   ! The distance between the flux sensor and start of 
                                                                   ! the absorber, each of the axes. The buffer 
                                                                   ! should be wide enough to allow field-driven
                                                                   ! re-crossing before absorption takes place
    integer(ik)            :: flux_sensors(3)  =  2                ! Number of sensor plates to use along each of the axes
    integer(ik)            :: flux_samples(3)  =  1                ! Number of k samples to use for each nominal k vector
                                                                   ! Long simulations may require oversampling to resolve
                                                                   ! spectral features.
    real(rk)               :: flux_max_p  (3)  = -1._rk            ! Max absolute value of momentum to include in flux
                                                                   ! continuation, along each of the axes. In a strong-field 
                                                                   ! calculation, large values of final momentum will
                                                                   ! correlate to very large kinetic momentum during the
                                                                   ! simulation, and will be severely distorted during
                                                                   ! propagation.
                                                                   ! Zero or negative means include all momenta
    integer(ik)            :: flux_grid_limits(2,3)                ! Extent of the meaningful part of the continued momentum,
                                                                   ! calculated from flux_max_p(:)
    integer(hik)           :: fft_thread_threshold = 10000         ! Minimum number of elements required to activate FFTW threads
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)              :: total_fields                       ! Total number of actual fields required for the input value of num_channels
    integer(ik)              :: f_free                             ! Index of the last unused position in
                                                                   ! f_table. All fields in f_table(1:f_free) are
                                                                   ! available for scratch use.
    integer(ik), allocatable :: f_table (:)                        ! List of all fields allocated
    integer(ik)              :: f_abs_pot                          ! Absorbing boundary potential
    !
    real(rk)                 :: CurrentElaser    = 0._rk           ! Current value of the laser field, used in LaserPlusIon energy function
    real(rk)                 :: CurrentEField(3) = 0._rk           ! Current value of the electric field vector
    real(rk)                 :: CurrentAField(3) = 0._rk           ! Current value of the vector potential
    real(rk)                 :: IonEnergyShift   = 0._rk           ! Value of ion-channel-dependent energy shift used in LaserPlusIon() energy function
    !
    real(rk)                 :: time                               ! Current time
    integer(ik)              :: tstep                              ! Current time step
    real(ark)                :: rotmat(3,3)                        ! Rotation matrix for the molecule, calculated from euler_angles()
    type(ecp_molecule)       :: ecp                                ! The effective core potential of this molecule; must be shared
                                                                   ! between all channels if more than one channel is present
    !
    contains
    !
    !  Prepare rotation matrix. Note that we'll be rotating the object, leaving the
    !  axes alone - so take the transpose. 
    !
    subroutine initialize_rotation
      !
      !  Rotation matrix
      !
      call MathRotationMatrix(euler_angles,rotmat)
      rotmat = transpose(rotmat)
      !
      if (verbose>=0 .and. .not.MathIsUnitMatrix(rotmat)) then
        write (out,"(/' Euler angles [Radian]:')")
        write (out,"( '    alpha = ',f12.6)") euler_angles(1)
        write (out,"( '    beta  = ',f12.6)") euler_angles(2)
        write (out,"( '    gamma = ',f12.6)") euler_angles(3)
        write (out,"( ' The Euler angles conventions here are:')")
        write (out,"( '   1. Rotate object by alpha around Z')")
        write (out,"( '   2. Rotate object by beta  around Y')")
        write (out,"( '   3. Rotate object by gamma around Z')")
        write (out,"( ' Rotation matrix:')")
        write (out,"(3x,3(1x,f12.8))") rotmat(1,:)
        write (out,"(3x,3(1x,f12.8))") rotmat(2,:)
        write (out,"(3x,3(1x,f12.8))") rotmat(3,:)
        write (out,"()")
      end if
    end subroutine initialize_rotation
    !
    !  Prepare simulation box and allocate a sufficient number of data fields
    !
    subroutine initialize_grid
      real(rk)    :: box(2,3)
      integer(ik) :: field
      !
      call TimerStart('Grid initialization')
      !
      write (out,"('Total number of grid fields required = ',i10)") total_fields
      write (out,"('      Estimated number of time steps = ',i10)") 1+nint((tMax-tMin)/dt)
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
      call MultiGridInit(max_grids=1,max_fields=total_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate (f_table(total_fields))
      allocate_fields: do field=1,total_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      f_free = total_fields
      !
      !  Set oversampling flag
      !
      write (out,"('Orbital oversampling = ',l8)") oversample_orbitals
      call gamess_oversample(oversample_orbitals)
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Load and convert the ECP
    !
    subroutine load_ion_ecp
      type(gam_structure) :: gam
      !
      if (ecp_file==' ') return
      !
      call TimerStart('Prepare ECP')
      call gamess_load_orbitals(file=ecp_file,structure=gam)
      call ecp_convert(gam,ecp,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report,rot=rotmat)
      call gamess_destroy(gam)
      call TimerStop('Prepare ECP')
    end subroutine load_ion_ecp
    !
    !  Electric fields represented by differentiating vector potentials
    !
    function VectorPotential( ft ) result(vp)
      real(rk), intent(in) :: ft  ! Time at which vector potential is needed
      real(rk)             :: vp
      !
      real(rk) :: StartTime, EndTime
      real(rk) :: t0
      real(rk) :: sigma
      real(rk) :: A0
  
      real(rk) :: sfTau, sfTurnOn, sfPlateau, fPre

      select case (ElectricFieldShape)
        case default
          write (out,"('Electric field shape ',a,' is not implemented')") trim(ElectricFieldShape)    
          stop 'dyson_tools%VectorPotential - in bad shape'
        case ('gaussian vp')
          !
          !  Gaussian vector potential
          ! 
          StartTime = TlaserOn  * 2._rk * pi / Wlaser
          EndTime   = TlaserOff * 2._rk * pi / Wlaser
          t0        = ft - 0.33_rk * (StartTime + EndTime)
          sigma     = pi / Wlaser / 2.0_rk
          A0        = Elaser * sigma * 1.6487_rk
          vp        = A0 * exp( -0.5_rk * (t0/sigma)**2 )
        case ('sine carrier gaussian envelope vp')
          !
          !  Gaussian envelope, sine carrier
          !
          t0 = ft - Tlaser
          A0 = Elaser/Wlaser
          vp = A0 * exp( -4._rk*log(2._rk) * (t0/FWHM)**2 ) * sin(Wlaser*t0)
        case ('flat-top vp')
          !
          !  Flat-top pulse
          !
          sfTau = twopi/Wlaser
          sfTurnOn = RampTime*sfTau
          sfPlateau = FlatTime*sfTau
          if (ft<0._rk) then
            fPre = 0._rk
          else if (ft<sfTurnOn) then
            fPre = (sin(pi*(ft/sfTurnOn)/2._rk))**2
          else if (ft<(sfTurnOn+sfPlateau)) then
            fPre = 1._rk
          else if (ft<(2.0*sfTurnOn+sfPlateau)) then
            fPre = (sin(pi*((ft-sfPlateau)/sfTurnOn)/2._rk))**2
          else
            fPre = 0._rk
          end if
          vp = fPre*sin(Wlaser*ft+Plaser*twopi)*Elaser/Wlaser
      end select
    end function VectorPotential
    !
    !
    !
    function ElectricField( ft )
      real(rk) :: ElectricField
      real(rk), intent(in) :: ft  ! Current time
      real(rk)             :: phi ! Current field phase
      real(rk)             :: ton
      !
      select case (ElectricFieldShape)
        case default
          !
          !  Assume that anything we do not recognize here is given by its
          !  vector-potential. 
          !
          ElectricField =  - ( (VectorPotential(ft+0.1_rk*dt)-VectorPotential(ft-0.1_rk*dt))/(0.2_rk*dt) )
        case ('static')
          ElectricField = Elaser
        case ('sine')
          !
          !  Sine wave, instant on and off
          !
          phi = Wlaser * ft 
          if (phi < TlaserOn*twopi .or. phi>= TlaserOff*twopi) then
            ElectricField = 0._rk
          else
            ElectricField = Elaser * sin(phi+Plaser*twopi)
          end if
        case ('half pulse sin^2')
          !
          !  half-pulse sin^2. The shape of the pulse is chosen to mimic the high-intensity
          !  part of the profile of a half-cycle sin pulse of frequency Wlaser
          !
          ton = pi / Wlaser / sqrt(2._rk)
          ElectricField = Elaser * sin(Wlaser*ft/sqrt(2._rk))**2
          if (ft .gt. (2._rk*ton)) then
            ElectricField = 0._rk
          end if
        case ('sin^2 envelope')
          !
          !  Sin^2 envelope, sine carrier
          !
          if ( ft<0 .or. ft>2*FWHM ) then
            ElectricField = 0._rk
          else
            ElectricField = Elaser * sin( 0.5_rk * pi * ft / FWHM )**2 * sin(Wlaser*(ft-FWHM)+Plaser*twopi)
          end if
        end select
    end function ElectricField
    !
    function ElectricField3D(time) result(e)
      real(ark), intent(in) :: time
      real(ark)             :: e(3)
      !
      e = edir * ElectricField(time)
    end function ElectricField3D
    !
    !  E-field preview - helps to avoid stupid mistakes with the input setup
    !
    subroutine preview_efield
      real(rk) :: afield
      !
      if (efield_preview==' ') return
      !
      call TimerStart('Preview E field')
      !
      !  Let's just simulate the laser field, and dump it to a file
      !
      open (unit=unit_efield,form='formatted',status='replace',file=trim(efield_preview))
      time      = tMin 
      tstep     = 0
      CurrentElaser = ElectricField(time)
      afield        = 0._rk
      time_loop: do while(time<tMax)
        CurrentElaser = ElectricField(time)
        write (unit_efield,"(1x,i10,3(1x,g17.10))") tstep+1, time, CurrentElaser, afield
        time   = time + dt 
        tstep  = tstep + 1
        afield = afield - CurrentElaser * dt
      end do time_loop
      close (unit_efield)
      write (out,"(/'Laser electric field evolution stored to file ',a)") trim(efield_preview)
      call TimerStop('Preview E field')
    end subroutine preview_efield
    !
    !  Very simple visualization routine
    !
    subroutine visualize_wavefunction(text,psi)
      character(len=*), intent(in) :: text
      integer(ik), intent(in)      :: psi
      !
      integer(ik) :: f_scr
      !
      call TimerStart('OpenDX visualization')
      call FieldVisualize(slot=0,src=psi,text=trim(text))
      !
      if (odx_plot_momentum) then
        f_scr = f_table(f_free) ; f_free = f_free - 1
        if (f_free<0) stop 'dyson_tools%visualize_wavefunction - not enough scratch fields'
        call FieldFFT(src=psi,dst=f_scr)
        call FieldVisualize(slot=1,src=f_scr,text=trim(text))
        f_free = f_free + 1
      end if
      !
      call FieldShow
      call TimerStop('OpenDX visualization')
    end subroutine visualize_wavefunction
    !
    !  Photoelectron spectrum visualization routine
    !
    subroutine visualize_photoelectron(text,f_flux)
      character(len=*), intent(in) :: text
      integer(ik), intent(in)      :: f_flux(:)  ! Array of momentum continuation fields
      !
      integer(ik) :: f_scr
      !
      call TimerStart('Visualization (PE)')
      !
      f_scr = f_table(f_free) ; f_free = f_free - 1
      if (f_free<0) stop 'dyson_tools%visualize_photoelectron - not enough scratch fields'
      !
      call FluxMergeSamples(f_flux,f_scr)
      call FieldVisualize(slot=3,src=f_scr,text=trim(text),limits=flux_grid_limits)
      f_free = f_free + 1
      !
      call TimerStop('Visualization (PE)')
    end subroutine visualize_photoelectron
    !
    !  An even simpler visualization routine
    !  Writes out the XY plane of 'src' field at Z=0 as an ascii array.
    ! 
    subroutine WriteSlice(dir,shift,src,filename)
      integer(ik), intent(in)      :: dir        ! Cartesian axis perpendicular to the slice
      integer(ik), intent(in)      :: shift      ! Shift of the slab position relative to mid-axis
      integer(ik), intent(in)      :: src        ! Data field for slice
      character(len=*), intent(in) :: filename
      !
      integer(ik)                  :: npts(3), nslb(3)
      integer(ik)                  :: pos
      complex(rk), allocatable     :: slice(:,:,:)
      !
      call TimerStart('WriteSlice')
      !
      npts      = FieldGridNPoints()
      nslb      = npts
      nslb(dir) = 1    ! Only one point along defining direction
      allocate (slice(nslb(1),nslb(2),nslb(3)))
      pos = npts(dir)/2_ik + shift
      call FieldFetchSlab(src,dir=dir,n1=pos,n2=pos,data=slice)
      !
      open(unit=unit_dump,form='formatted',status='replace',file=filename)
      write(unit_dump,"(2(1x,f20.10))") slice
      close(unit_dump)
      !
      deallocate(slice)
      call TimerStop('WriteSlice')
    end subroutine WriteSlice
    !
    !  Construct absorbing potential for the boundary
    !
    subroutine buildCAP(f_abs_pot)
      integer(ik), intent(in) :: f_abs_pot  ! Field for the absorbing potential
      !
      integer(ik) :: axis
      !
      call TimerStart('Build CAP')
      !
      !  More boundary types are available - see CAPsetUp in caps.f90 if needed.
      !
      select case (abs_bound_type)
        case default
          write (out,"('Absorbing boundary geometry ',a,' is not implemented')") trim(abs_bound_type)
          stop 'dyson_tools%buildCAP - bad CAP type'
        case ('all','ALL')
          axis = 0
        case ('xy','XY')
          axis = 13
        case ('xz','XZ')
          axis = 12
        case ('yz','YZ')
          axis = 11
      end select
      !
      call CAPsetUp(type='Manolopoulos',dt=dt,axis=axis,kmin=abs_bound_kmin)
      !
      !  It's not a good idea to use mask= here: the absorbing bondary is very likely to be outside
      !  the spatial region where wavefunction is significant, so it may be silently discarded. Oops.
      !
      call FieldInit(f_abs_pot,CAPsplitPotential)
      call TimerStop('Build CAP')
    end subroutine buildCAP
    !
    !  Report molecular geometry
    !
    subroutine report_structure(level,name)
      integer(ik), intent(in)      :: level ! Minimal verbosity level to print the structure;
                                            ! Report just the number of atoms otherwise
      character(len=*), intent(in) :: name  ! Name of the structure/file
      !
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call gamess_report_nuclei(nnuc)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      if (level>verbose) return
      !
      !  Tell a bit more!
      !
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
    end subroutine report_structure
    !
    !  Constructs part of the eq. 32, which is constant or linear in the laser electric field
    !
    function LaserPlusIon(r)
      complex(rk) :: LaserPlusIon
      real(rk), intent(in) :: r(3)
      LaserPlusIon = -CurrentElaser * dot_product(Edir, r) + IonEnergyShift
    end function LaserPlusIon
    !
    subroutine dump_potential(f_core_pot,potential_file)
      integer(ik), intent(in)      :: f_core_pot
      character(len=*), intent(in) :: potential_file
      !
      call TimerStart('Dump core potential')
      open (unit=unit_dump,form='formatted',status='replace',file=potential_file)
      call FieldDump(unit_dump,src=f_core_pot)
      close (unit=unit_dump)
      call TimerStop('Dump core potential')
    end subroutine dump_potential
  end module dyson_tools
