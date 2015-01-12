!
!  Auxiliary routines for analytical continuation of long-range photoelectron
!  distributions, a-la the "t-SURFF" method of Tao and Scrinzi (arXiv:1109.4053v1)
!
!  This routine will likely be incomprehensible without the notes in:
!
!      cylindrical-momentum-continuation_v1.pdf
!      cylindrical-momentum-continuation_v2.pdf
!      cylindrical-momentum-continuation_technical.pdf
!
!  IMPORTANT: The flux method is not appropriate for very sharp spectral features
!  IMPORTANT: However, it is perfect for strong-field processes, where widths of
!  IMPORTANT: the photoelectron lines are broder than the momentum grid resolution.
!
!  The code is necessarily specific to the grid topology. Currently, two grid
!  topologies are implemented: 'cylindrical' and '3D Cartesian'
!
!  WARNING: Prior to 21 Feb 2013, implementation of oversampling was (very slightly)
!  WARNING: wrong, with a correction term containing momentum offset missing from 
!  WARNING: the surface-sensing integrals. Oops. Unfortunately, the correct implementation
!  WARNING: is quite a bit slower than the wrong one, since sensing now has to be
!  WARNING: repeated for each sampling point.
!
 module flux
   use accuracy
   use bessel
   use fftw
   use math
   use multigrid
   use caps
   use timer
   implicit none

   private
   public FluxInitialize, FluxTimeStep, FluxMergeSamples

   integer(ik), parameter       :: verbose          =  2    ! Level of verboseness, what else?
   integer(ik), parameter       :: electron_q       =  1    ! Electron charge; choose either +1 or -1
                                                            ! The choice will affect interpretation of the
                                                            ! vector-potential.
   !
   !  Flux sensor location
   !
   integer(ik), save            :: flux_sensors(3)          ! Number of sensor plates in each direction; in 2D only the
                                                            ! first two entries are valid
   !
   !  The constants below are specific to the cylindrical grid topology
   !
   integer(ik), save            :: flux_rhopos              ! Positions of the first sensor plane in the rho direction 
   integer(ik), save            :: flux_zminpos             ! Two sensor planes are needed; the second one is one grid
   integer(ik), save            :: flux_zmaxpos             ! point further away from the origin. 
   integer(ik), save            :: cyl_nrho                 ! Grid spacing along the rho direction
   integer(ik), save            :: cyl_nz                   ! Grid spacing along the z direction
   real(rk), save               :: cyl_drho                 ! Grid spacing along the rho direction
   real(rk), save               :: cyl_dz                   ! Grid spacing along the z direction
   real(rk), save               :: cyl_dkrho                ! Momentum grid spacing along rho
   real(rk), save               :: cyl_dkz                  ! Momentum grid spacing along z
   real(rk), allocatable, save  :: cyl_krho_tab(:)          ! Values of momentum resolved on our grid, rho direction
   real(rk), allocatable, save  :: cyl_rho_tab (:)          ! Values of coordinate resolved on our grid
   real(rk), allocatable, save  :: cyl_kz_tab(:)            ! Values of momentum resolved on our grid, z direction
   real(rk), allocatable, save  :: cyl_z_tab (:)            ! Values of coordinate resolved on our grid
   real(rk), allocatable, save  :: bessel_tab(:,:,:)        ! Nominally, this is a table BesselJ(0,krho*rho)*rho at
                                                            ! all grid points. In practice, see comments below, in
                                                            ! fill_bessel_tab. The last index is the "subsample" index,
                                                            ! with the krho adjusted by the corresponding momentum shift
                                                            ! from p_offset(1,:)*cyl_dkrho
                                                            ! First index is krho index; second index is rho
   real(rk), allocatable, save  :: struve_tab(:,:,:)        ! Same as bessel_tab, except for StruveH
   !
   !  Constants below are specific to the Cartesian-product 3D grids
   !
   integer(ik), save            :: flux_pos(2,3)            ! Cartesian 3D: Same meaning as flux_zminpos/flux_zmaxpos above;
                                                            ! however, there are paired sensor surfaces on both sides now.
                                                            ! First index: 1 is minpos; 2 is maxpos
   integer(ik), save            :: crt_np(3)                ! Number of grid points along each direction
   real(rk), save               :: crt_dx(3)                ! Spatial grid spacing along each direction
   real(rk), save               :: crt_dk(3)                ! Momentum grid spacing along each direction
   type real_array
     real(rk), allocatable :: x(:)
   end type real_array
   type(real_array), save       :: crt_r_tab(3)             ! Coordinates of spatial grid points; since dimensions of arrays
                                                            ! can be different for X/Y/Z, the declaration is a bit elaborate
   type(real_array), save       :: crt_k_tab(3)             ! Coordinates of momentum grid points
   !
   !  Accumulation limits in momentum space. Momentum positions outside these limits are treated as zero
   !
   integer(ik), save            :: p_limits(2,3)            ! The first index is the lower(=1)/upper(=2) limit
                                                            ! The second index is the coordinate 
                                                            !   Cylindrical:  Z=1, Rho=2
                                                            !   Cartesian 3D: X=1, Y=2, Z=3
   !
   !  Sub-sampling in momentum space
   !
   integer(ik), save            :: p_samples                ! Number of p samples
   real(rk), allocatable, save  :: p_offset(:,:)            ! p_offset(:,i) is the offset of i-th sample
   real(rk), allocatable, save  :: p_weight  (:)            ! integration weights
   !
   !  Vector-potential and it's integrals; needed for tracking contributions to the 
   !  appropriate final-momentum state.
   !
   integer(ik), save            :: lastTimeStep             ! Total number of time steps is one higher (we count from 0)
   real(ark), allocatable, save :: time      (:)            ! Time at a given time step
   real(ark), allocatable, save :: timeA   (:,:)            ! Vector-potential at each time step, 0-based
   real(ark), allocatable, save :: timeAint(:,:)            ! Integral of vector-potential from the starting time to the
                                                            ! time step
   real(ark), allocatable, save :: timeA2    (:)            ! Integral of the square of the vector-potential
   !
   !  Parameters needed to evaluate the phases of the current step's contribution
   !  to the total integral. The phases are a quadratic function of the meomentum,
   !  so this is easy enough ...
   !
 contains
   !
   !  Externally visible interfaces
   !
   subroutine FluxInitialize(mode,tMin,tMax,dt,a0,evaluateF,flux_guard,sensors,samples,max_p,grid_limits)
     character(len=*), intent(in) :: mode           ! Must be 'cylindrical'
     real(ark), intent(in)        :: tMin           ! Initial time, time step 0
     real(ark), intent(in)        :: tMax           ! Final time
     real(ark), intent(in)        :: dt             ! Time step
     real(ark), intent(in)        :: a0(3)          ! Initial vector-potential
!    external                     :: evaluateF      ! Function for evaluating electric field as a function
!                                                   ! of time; takes real(ark) time as an argument
     real(rk), intent(in)         :: flux_guard(:)  ! Extra space between flux sensor and absorber, to
                                                    ! accommodate field-driven oscillations.
     integer(ik), intent(in)      :: sensors(:)     ! Number of sensor plates along each boundary
     integer(ik), intent(in)      :: samples(:)     ! Number of sampling points to use for each p grid
                                                    ! direction. Accumulation of the photoelectrion spectrum
                                                    ! will require product(samples) fields.
     real(rk), intent(in)         :: max_p(:)       ! Maximum absolute momentum to include; zero or negative
                                                    ! means all momenta
     integer(ik), intent(out), optional :: grid_limits(:,:) ! Copy of p_limits(:,:), in case caller wants to visualize 
                                                            ! just the meaniningful part of the momentum grids
     !
     interface 
       function evaluateF(time) result(e)
         use accuracy
         real(ark), intent(in) :: time
         real(ark)             :: e(3)
       end function evaluateF
     end interface
     !
     integer(ik) :: ic   ! Cartesian component
     integer(ik) :: ndim ! Effective dimensionality of the problem
     !
     call TimerStart('FluxInitialize')
     flux_sensors(1:size(sensors)) = sensors
     select case (mode)
       case default
         stop 'flux%FluxInitialize - unknown grid choice'
       case ('cylindrical')
         ndim = 2
         !
         !  Sanity check on array parameters
         !
         if (size(flux_guard)/=2 .or. size(sensors)/=2 .or. size(samples)/=2 .or. size(max_p)/=2) then
           stop 'flux%FluxInitialize - cylindrical case called with array size(s) /=2'
         end if
         !
         !  First, I need to determine the position of the sensor planes
         !     
         flux_zminpos = locate_sensor(2_ik,-1_ik,flux_guard(1))
         flux_zmaxpos = locate_sensor(2_ik, 1_ik,flux_guard(1))
         flux_rhopos  = locate_sensor(2_ik, 2_ik,flux_guard(2))
         !
         !  Determine accumulation range for the momentum grid. This is a
         !  performance-optimization parameter.
         !
         p_limits(1,1) = locate_p_limit(2_ik,-1_ik,max_p(1))
         p_limits(2,1) = locate_p_limit(2_ik, 1_ik,max_p(1))
         p_limits(1,2) = 1
         p_limits(2,2) = locate_p_limit(2_ik, 2_ik,max_p(2))
         !
         !  Prepare p sub-sampling points
         !
         call subsample_p(2_ik,samples)
         !
         !  Fill table of Bessel functions, we'll need those for Fourier-
         !  Bessel transforms at each time step.
         !
         call fill_bessel_tab(2_ik)
       case ('3D Cartesian')
         ndim = 3
         !
         !  Sanity check on array parameters
         !
         if (size(flux_guard)/=3 .or. size(sensors)/=3 .or. size(samples)/=3 .or. size(max_p)/=3) then
           stop 'flux%FluxInitialize - 3D Cartesian case called with array size(s) /=3'
         end if
         locate_cartesian_sensors: do ic=1,3
           !
           !  First, I need to determine the position of the sensor planes.
           !     
           flux_pos(1,ic) = locate_sensor(3_ik,-ic,flux_guard(ic))
           flux_pos(2,ic) = locate_sensor(3_ik, ic,flux_guard(ic))
           !
           !  Determine accumulation range for the momentum grid. This is a
           !  performance-optimization parameter; passing a negative number for the
           !  last argument should return the same results _within_ the -/+max_p range
           !
           p_limits(1,ic) = locate_p_limit(3_ik,-ic,max_p(ic))
           p_limits(2,ic) = locate_p_limit(3_ik, ic,max_p(ic))
         end do locate_cartesian_sensors
         !
         !  Prepare p sub-sampling points
         !
         call subsample_p(3_ik,samples)
         !
         call fill_cartesian_tables(3_ik)
     end select
     !
     !  I will need an integral of the vector potential and its square
     !  at every time point.
     !
     call fill_vp_integrals(tMin,tMax,dt,a0,evaluateF)
     !
     if (present(grid_limits)) then
       if (size(grid_limits,dim=1)/=2 .or. size(grid_limits,dim=2)/=ndim) then
         stop 'flux%FluxInitialize - grid_limits argument has incorrect size'
       end if
       grid_limits = p_limits(:,:ndim)
     end if
     !
     call TimerStop('FluxInitialize')
   end subroutine FluxInitialize
   !
   subroutine FluxTimeStep(mode,f_psi,tstep,ctime,aField,f_momentum,f_scr1,f_scr2)
     character(len=*)        :: mode          ! Must be 'cylindrical'
     integer(ik), intent(in) :: f_psi         ! Current wavefunction handle
     integer(ik), intent(in) :: tstep         ! Current time step
     real(ark), intent(in)   :: ctime         ! Caller's idea of the current time; sanity check
     real(ark), intent(in)   :: aField(:)     ! Caller's idea of the vector-poential at the current time; sanity check
     integer(ik), intent(in) :: f_momentum(:) ! Final momentum in the continuum, will be updated 
                                              ! for the current time step
     integer(ik), intent(in) :: f_scr1        ! Scratch fields, used to compose current time step's
     integer(ik), intent(in) :: f_scr2        ! contribution to the momentum
     !
     complex(rk), allocatable :: wgt_rho(:), wgt_z(:) ! Weights along each coordinate
     complex(rk), allocatable :: wgt_x(:), wgt_y(:)
     integer(ik)              :: alloc
     integer(ik)              :: ndim         ! Dimensionality of the problem
     integer(ik)              :: ic           ! Cartesian dimension
     integer(ik)              :: is           ! Sensor ordinal, 1 .. flux_sensors
     integer(ik)              :: sample       ! Momentum sample ordinal, 1 .. p_samples
     real(rk)                 :: sample_dk(3) ! Momentum shift for all the grid points during oversampling,
                                              ! in fractions of the corresponding momentum grid resolution.
     real(ark)                :: AFtemp(3)    ! Temporary for "shifted" vector-potential during oversampling
     real(rk)                 :: wgt_s        ! Sensor weight
     real(ark)                :: wgt_a        ! Coefficient in front of k**2
     real(ark)                :: wgt_b(3)     ! Coefficient in front of -2*k
     real(ark)                :: wgt_c        ! Free coefficient
     complex(ark)             :: wgt_scale    ! Overall scale coefficient
     !
     call TimerStart('FluxTimeStep')
     !
     !  Do a bit of sanity checking before we do anything else. The caller's idea
     !  of the current time and vector-potential must match ours - otherwise, the
     !  things will go wrong in mysterious ways.
     !
     if (tstep<0 .or. tstep>lastTimeStep-1) stop 'flux%FluxTimeStep - bad time step index'
     if (abs(ctime-time(tstep))>10._ark*spacing(ctime)) then
       write (out,"('At time step ',i15,' caller''s time is ',g25.15,' while ours is ',g25.15)") &
              tstep, ctime, time(tstep)
       stop 'flux%FluxTimeStep - bad time accounting'
     end if
     if (any(abs(aField-timeA(:,tstep))>10._ark*spacing(maxval(abs(aField))))) then
       write (out,"('At time step ',i15,' caller''s VP is '/3g25.15/' while ours is '/3g25.15)") &
              tstep, aField, timeA(:,tstep)
       stop 'flux%FluxTimeStep - bad VP accounting'
     end if
     if (size(f_momentum)/=p_samples) stop 'flux%FluxTimeStep - bad number of accumulators'
     !
     !  Topology-specific setup
     !
     select case (mode)
       case default; stop 'flux%FluxTimeStep - unknown grid choice (0)'
       case ('cylindrical')
         ndim = 2
         !
         !  Sanity check - make sure fields are compatible with cylindrical symmetry
         !
         if (any(aField(2:)/=0)) stop 'flux%FluxTimeStep - aField is not cylindrical'
         !
         allocate (wgt_rho(cyl_nrho),wgt_z(cyl_nz),stat=alloc)
       case ('3D Cartesian')
         ndim = 3
         allocate (wgt_x(crt_np(1)),wgt_y(crt_np(2)),wgt_z(crt_np(3)),stat=alloc)
     end select
     if (alloc/=0) stop 'flux%FluxTimeStep - bad weight allocation'
     !
     !  Phase of the length-gauge Volkov solutions is a quadratic function of the
     !  momentum. Let's calculate the coefficients, using tabulated integrals of
     !  the vector-potential.
     !
     wgt_a     = time(lastTimeStep) - time(tstep)
     wgt_b     = electron_q * (timeAint(:,lastTimeStep) - timeAint(:,tstep))
     wgt_c     = electron_q**2 * (timeA2(lastTimeStep) - timeA2(tstep))
     wgt_scale = (0,-1)*(time(tstep+1)-time(tstep))/(sqrt2pi**3)
     !
     if (verbose>=3) then
       write (out,"(/' Volkov phase parameters [-0.5*I*(a*K**2 - 2*b.K +c)]:')")
       write (out,"( '         a = ', f25.15)") wgt_a
       write (out,"( '         b = ',3f25.15)") wgt_b
       write (out,"( '         c = ', f25.15)") wgt_c
       write (out,"( ' prefactor = ',2f25.15/)") wgt_scale
     end if
     !
     !  This quadratic form can be evaluated as a sum of component phases,
     !  which makes things much faster:
     !
     !   -0.5*I*[ (c-|b|^2/a) + a*(k_1-b_1/a)^2 + ... ]
     !
     if (abs(wgt_a)<=spacing(10._rk)) then
        ! This could only happen if t=t_{final}, and the overall phase correction is zero
        wgt_b = 0 
        wgt_c = 0
     else
        wgt_b = wgt_b / wgt_a
     end if
     wgt_c = wgt_c - sum(wgt_b**2)*wgt_a
     !
     !  Ready to sense and accumulate at each momentum sub-sample
     !
     momentum_samples: do sample=1,p_samples
       sample_dk(1:ndim) = p_offset(:,sample)
       !
       !  Sensing part: we evaluate the g-integrals from our notes.
       !  These integrals are specific to the boundary topology.
       !
       call TimerStart('FluxTimeStep: Sense')
       call FieldZero(f_scr1,limits=p_limits(:,1:ndim))
       select case (mode)
         case default ; stop 'flux%FluxTimeStep - unknown grid choice (1)'
         case ('cylindrical')
           !
           !  Our surface is a cylinder; we have two "lids" (the top and the bottom), and the side surface.
           !  For the cylindrical grid, oversampling has to be handled differently for the kz and krho
           !  components. For kz, we can simply shift the the vector-potential by (-dk/electron_q).
           !  For krho, we have to pre-calculate an extra Bessel and Struve table for each shift value,
           !  and use the correct table in evaluating the integrals.
           !
           AFtemp(1) = aField(1) - cyl_dkz*sample_dk(1)/electron_q
           call TimerStart('FluxTimeStep: Sense: Lid')
           wgt_s = 1.0_rk/flux_sensors(1)
           measure_flux_lids: do is=0,flux_sensors(1)-1
             call surface_cylindrical_lid (sample,f_psi,f_scr1,flux_zmaxpos-is,flux_zmaxpos+1-is,aFtemp(1),wgt_s)
             call surface_cylindrical_lid (sample,f_psi,f_scr1,flux_zminpos+is,flux_zminpos-1+is,aFtemp(1),wgt_s)
           end do measure_flux_lids
           call TimerStop('FluxTimeStep: Sense: Lid')
           !
           call TimerStart('FluxTimeStep: Sense: Side')
           wgt_s = 1.0_rk/flux_sensors(2)
           measure_flux_sides: do is=0,flux_sensors(2)-1
             call surface_cylindrical_side(sample,f_psi,f_scr1,flux_rhopos -is,flux_rhopos +1-is,aFtemp(1),wgt_s)
           end do measure_flux_sides
           call TimerStop('FluxTimeStep: Sense: Side')
           !
         case ('3D Cartesian')
           !
           !  Our surface is a cube; we have two sides per Cartesian dimension.
           !  For the Cartesian product grid, shifting all momentum points by (dk) in the g integral
           !  is equivalent to shifting the vector-potential by (-dk/electron_q).
           !
           AFtemp = aField - crt_dk*sample_dk(1:ndim)/electron_q
           measure_flux_cartesian: do ic=1,3
             wgt_s = 1.0_rk/flux_sensors(ic)
             measure_flux_cartesian_side: do is=0,flux_sensors(ic)-1
               call surface_cartesian(ic,f_psi,f_scr1,flux_pos(1,ic)+is,flux_pos(1,ic)-1+is,aFtemp,wgt_s)
               call surface_cartesian(ic,f_psi,f_scr1,flux_pos(2,ic)-is,flux_pos(2,ic)+1-is,aFtemp,wgt_s)
             end do measure_flux_cartesian_side
           end do measure_flux_cartesian
           !
       end select
       call TimerStop('FluxTimeStep: Sense')
       !
       !  We can now calculate the Volkov-phase correction to the surface integrals, 
       !  and to accumulate the current contribution to the time integral
       !
       call TimerStart('FluxTimeStep: Accumulate')
       call FieldZero(f_scr2,limits=p_limits(:,:ndim))
       select case (mode)
         case default ; stop 'flux%FluxTimeStep - unknown grid choice (2)'
         case ('cylindrical')
           call weight_1d(wgt_a,wgt_b(1),cyl_kz_tab  +p_offset(1,sample)*cyl_dkz,  wgt_z)
           call weight_1d(wgt_a,wgt_b(2),cyl_krho_tab+p_offset(2,sample)*cyl_dkrho,wgt_rho)
           wgt_rho = wgt_rho * wgt_scale * exp(-(0._rk,0.5_rk)*wgt_c)
           call FieldGER(f_scr2,wgt_z,wgt_rho,limits=p_limits(:,:ndim))
         case ('3D Cartesian')
           call weight_1d(wgt_a,wgt_b(1),crt_k_tab(1)%x+p_offset(1,sample)*crt_dk(1),wgt_x)
           call weight_1d(wgt_a,wgt_b(2),crt_k_tab(2)%x+p_offset(2,sample)*crt_dk(2),wgt_y)
           call weight_1d(wgt_a,wgt_b(3),crt_k_tab(3)%x+p_offset(3,sample)*crt_dk(3),wgt_z)
           wgt_x = wgt_x * wgt_scale * exp(-(0._rk,0.5_rk)*wgt_c)
           call FieldGER(f_scr2,wgt_x,wgt_y,wgt_z,limits=p_limits(:,:ndim))
       end select
       call FieldMulAdd(dst=f_momentum(sample),src_a=f_scr1,src_b=f_scr2,limits=p_limits(:,:ndim))
       call TimerStop('FluxTimeStep: Accumulate')
     end do momentum_samples
     !
     !  Topology-specific termination
     !
     select case (mode)
       case default; stop 'flux%FluxTimeStep - unknown grid choice (3)'
       case ('cylindrical')
         deallocate(wgt_rho,wgt_z)
       case ('3D Cartesian')
         deallocate(wgt_x,wgt_y,wgt_z)
     end select
     call TimerStop('FluxTimeStep')
   end subroutine FluxTimeStep
   !
   !  Perform averaging over momentum sub-sample grid
   !
   subroutine FluxMergeSamples(f_momentum,f_ave)
     integer(ik), intent(in) :: f_momentum(:)  ! Time integrals for momentum sub-samples
     integer(ik), intent(in) :: f_ave          ! Average of the time integrals
     !
     integer(ik) :: is
     !
     call TimerStart('FluxMergeSamples')
     if (size(f_momentum)/=p_samples) stop 'flux%FluxTimeStep - bad number of accumulators'
     !
     if (p_samples==1) then
       call FieldCopy(src=f_momentum(1),dst=f_ave)
     else
       !
       if (verbose>=0) then
         write (out,"('Photoelectron distribution averaged over ',i0,' k sub-samples')") p_samples
         write (out,"('Phases have been arbitralily set to zero')")
       end if
       call FieldZero(f_ave)
       momentum_samples: do is=1,p_samples
         call FieldRhoAccumulate(alpha=p_weight(is),src=f_momentum(is),dst=f_ave)
       end do momentum_samples
       call FieldProcess(f_ave,field_sqrt)
     end if
     !
     call TimerStop('FluxMergeSamples')
   end subroutine FluxMergeSamples
   !
   function field_sqrt(coord,val) result(f)
     real(rk), intent(in)    :: coord(*)  ! Coordinates of the point; ignored
     complex(rk), intent(in) :: val       ! Field value; must be real
     complex(rk)             :: f
     !
     f = sqrt(abs(val))
   end function field_sqrt
   !
   !  Internal routines below this point
   !
   subroutine weight_1d(a,b,k,wgt)
     real(rk), intent(in) :: a, b       ! weight factor: exp(-0.5*I*A*(K-B)**2)
     real(rk), intent(in) :: k(:)       ! table of K values
     complex(rk), intent(out) :: wgt(:) ! Table of phase factors
     !
     if (size(k)/=size(wgt)) stop 'flux%weight_id - oops'
     !
     wgt(:) = exp(-(0._rk,0.5_rk)*a*(k(:)-b)**2)
   end subroutine weight_1d
   !
   function locate_sensor(dim,axis,guard) result(ind)
     integer(ik), intent(in) :: dim      ! Dimensionality of the grid
     integer(ik), intent(in) :: axis     ! Axis perpendicular to the sensor surface; negative means looking down
     real(rk), intent(in)    :: guard    ! Flux guard for this direction
     !
     integer(ik)             :: nd(dim)
     integer(ik)             :: ind      ! Sensor plane closest to the origin; the second plane is at ind + sign(1,axis)
     real(rk)                :: rcap     ! Position where absorbing potential starts
     integer(ik)             :: ngrid    ! Number of grid points
     real(rk), allocatable   :: coord(:) ! Grid coordinates
     integer(ik)             :: alloc
     integer(ik)             :: ic, dir, pos
     !
     if (.not.CAPgetZone(axis,rcap)) stop 'flux%locate_sensor - missing a CAP'
     nd    = FieldGridNPoints()
     ngrid = nd(abs(axis))
     allocate (coord(ngrid),stat=alloc)
     if (alloc/=0) stop 'flux%locat_sensor - allocation failure'
     call FieldGridCoordinates(abs(axis),coord)
     if (axis>0) then
       scan_up: do ind=1,ngrid-2,1
         if (coord(ind+2)<rcap-guard) cycle scan_up
         exit scan_up
       end do scan_up
     else
       scan_down: do ind=ngrid,3,-1
         if (coord(ind-2)>rcap+guard) cycle scan_down
         exit scan_down
       end do scan_down
     end if
     !
     if (verbose>=0) then
       dir = sign(1_ik,axis)
       print_sensors: do ic=1,flux_sensors(abs(axis))
         pos = ind - dir*(ic-1)
         write (out,"('For axis ',i2,' flux sensor ',i3,' is at ',i8,' r1,r2 = ',2f12.5,' rcap = ',f12.5)") &
                axis, ic, pos, coord(pos), coord(pos+dir), rcap
       end do print_sensors
     end if
     !
     deallocate (coord)
   end function locate_sensor
   !
   function locate_p_limit(dim,axis,p_max) result(ind)
     integer(ik), intent(in) :: dim      ! Dimensionality of the grid
     integer(ik), intent(in) :: axis     ! Momentum axis; negative means looking at negative momenta
     real(rk), intent(in)    :: p_max    ! Maximum absolute value of momentum to handle
                                         ! p_max<=0 means treat all values
     integer(ik)             :: ind      ! First (axis<0) or last (axis>0) point to include
     !
     integer(ik)             :: nd(dim)
     integer(ik)             :: ngrid    ! Number of grid points
     real(rk), allocatable   :: coord(:) ! Grid coordinates
     integer(ik)             :: alloc
     !
     nd    = FieldGridNPoints()
     ngrid = nd(abs(axis))
     allocate (coord(ngrid),stat=alloc)
     if (alloc/=0) stop 'flux%locate_p_limit - allocation failure'
     !
     call FieldGridCoordinates(-abs(axis),coord)  ! Ie ask for the momenta
     if (p_max<=0._rk) then
       ind = 1
       if (axis>0) ind = ngrid
     else if (axis<0) then
       scan_up: do ind=1,ngrid
         if (abs(coord(ind))<=p_max) exit scan_up
       end do scan_up
     else
       scan_down: do ind=ngrid,1,-1
         if (abs(coord(ind))<=p_max) exit scan_down
       end do scan_down
     end if
     !
     if (verbose>=0) then
       if (axis>0) then 
         write (out,"('For axis ',i2,' upper limit for flux momentum is ',i6,' = ',f14.7)") abs(axis), ind, coord(ind)
       else
         write (out,"('For axis ',i2,' lower limit for flux momentum is ',i6,' = ',f14.7)") abs(axis), ind, coord(ind)
       end if
     end if
     !
     deallocate (coord)
   end function locate_p_limit
   !
   !  We will need to tabulate values of BesselJ(0,krho*rho)*rho and 
   !  StruveH(0,krho*rho)*rho.
   !
   !  These functions are significantly variable over our grid volumes.
   !  To get reasonable results, we have to average them over the volume
   !  elements; see bes_smooth() and str_smooth() below.
   !
   !  IMPORTANT: These function should NOT be averaged over the momentum!
   !
   subroutine fill_bessel_tab(dim)
     integer(ik), intent(in) :: dim  ! Dimensionality of the current space. This _must_ be 2!
     integer(ik)             :: nd(dim), ikrho, irho, alloc, sample
     real(rk)                :: dx(dim), krho, rho, rhou, rhod, val, dkrho
     !
     nd   = FieldGridNPoints()
     dx   = FieldGridSpacing()
     cyl_dz   = dx(1)
     cyl_drho = dx(2)
     cyl_nz   = nd(1)
     cyl_nrho = nd(2)
     allocate (cyl_krho_tab(cyl_nrho),cyl_rho_tab(cyl_nrho), &
               cyl_kz_tab(cyl_nz),cyl_z_tab(cyl_nz), &
               bessel_tab(cyl_nrho,cyl_nrho,p_samples), &
               struve_tab(cyl_nrho,cyl_nrho,p_samples), stat=alloc)
     if (alloc/=0) stop 'flux%fill_bessel_tab - allocation failure'
     call FieldGridCoordinates( 1,   cyl_z_tab)
     call FieldGridCoordinates(-1,  cyl_kz_tab)
     call FieldGridCoordinates( 2, cyl_rho_tab)
     call FieldGridCoordinates(-2,cyl_krho_tab)
     !
     cyl_dkrho = cyl_krho_tab(2)-cyl_krho_tab(1)
     cyl_dkz   = cyl_kz_tab(2)-cyl_kz_tab(1)
     !
     !  First index of bessel_tab is ikrho
     !
     fill_sample_offsets: do sample=1,p_samples
       dkrho = p_offset(2,sample)*cyl_dkrho
       write (out,"('sample = ',i5,' dkrho = ',g25.15)") sample, dkrho
       !$omp parallel do default(none) shared(cyl_nrho,bessel_tab,cyl_drho,cyl_rho_tab,cyl_krho_tab) &
       !$omp&            shared(cyl_dkrho,struve_tab,sample,dkrho) &
       !$omp&            private(irho,ikrho,krho,rho,val,rhou,rhod)
       fill_bessel_rho: do irho=1,cyl_nrho
         rho  = cyl_rho_tab(irho) + dkrho
         rhod = rho - 0.5_rk * cyl_drho
         rhou = rho + 0.5_rk * cyl_drho
         !
         fill_bessel_k: do ikrho=1,cyl_nrho
           krho  = cyl_krho_tab(ikrho)
           val   = bes_smooth(rhod,rhou,krho)
           bessel_tab(ikrho,irho,sample) = val
           val   = str_smooth(rhod,rhou,krho)
           struve_tab(ikrho,irho,sample) = val
         end do fill_bessel_k
       end do fill_bessel_rho
       !$omp end parallel do
     end do fill_sample_offsets
   end subroutine fill_bessel_tab
   !
   !  Average the function: 
   !
   !     f(r,k) = r BesselJ[0,r k]
   !
   !  over a surface element in both r and k. The integral over r can be done
   !  analytically:
   !
   !     /
   !     | f(r,k) d k = (r/k) BesselJ[1,r k]
   !     /
   !
   !  The integral over k is known, but ugly - so that we choose to integrate
   !  numerically.
   !
   !  The same procedure applies to the StruveH (except the second integral
   !  is even uglier)
   !
   !  WARNING: Averaging over the moments is NOT a correct procedure. The
   !  WARNING: spikes in the momentum spectra come from the time integral,
   !  WARNING: and should be dealt with using a finer momentum grid at that
   !  WARNING: step. The integration over k has been moved to flux_rejects.f90
   !
   !
   !  Bessel function, smoothed over a spatial grid patch
   !
   function bes_smooth(r0,r1,k) result(bes)
     real(rk), intent(in) :: r0, r1 ! Volume element limits in rho
     real(rk), intent(in) :: k      ! Desired k
     real(rk)             :: bes
     !
     bes = (fl(r1,k) - fl(r0,k))/(r1-r0)
     !
     contains 
       !
       !  (r/k)*BesselJ(1,k*r), numerically accurate at small k
       !
       real(rk) function fl(r,k)
         real(rk), intent(in) :: r, k
         !
         if (abs(k*r)<=(spacing(10._rk)**0.25_rk)) then
           fl = 1.0_rk - (r*k)**2/8._rk + (r*k)**4/192._rk
           fl = 0.5_rk * r**2 * fl
         else
           fl = (r/k)*BesselJ(1,r*k)
         end if
       end function fl
   end function bes_smooth
   !
   !  Struve function, smoothed over a spatial grid patch
   !
   function str_smooth(r0,r1,k) result(str)
     real(rk), intent(in) :: r0, r1 ! Volume element limits in rho
     real(rk), intent(in) :: k      ! Desired k
     real(rk)             :: str
     !
     str = (fl(r1,k) - fl(r0,k))/(r1-r0)
     !
     contains 
       !
       !  (r/k)*StruveH(1,k*r), numerically accurate at small k
       !
       real(rk) function fl(r,k)
         real(rk), intent(in) :: r, k
         !
         if (abs(k*r)<=(spacing(10._rk)**0.25_rk)) then
           fl = 1.0_rk - (r*k)**2/15._rk + (r*k)**4/525._rk
           fl = (2._rk*k*r**3/(3._rk*pi)) * fl
         else
           fl = (r/k)*StruveH(1,r*k)
         end if
       end function fl
   end function str_smooth
   !
   !  Prepare data tables for Cartesian flux calculation, fori anything needed at each 
   !  time step.
   !
   subroutine fill_cartesian_tables(dim)
     integer(ik), intent(in) :: dim  ! Dimensionality of the current space. This _must_ be 3!
     !
     integer(ik) :: nd(dim), alloc
     real(rk)    :: dx(dim)
     integer(ik) :: ic
     !
     nd   = FieldGridNPoints()
     dx   = FieldGridSpacing()
     fill_coordinates: do ic=1,3
       allocate (crt_r_tab(ic)%x(nd(ic)), crt_k_tab(ic)%x(nd(ic)), stat=alloc)
       if (alloc/=0) stop 'flux%fill_cartesian_tables - allocation failure'
       call FieldGridCoordinates( ic,crt_r_tab(ic)%x)
       call FieldGridCoordinates(-ic,crt_k_tab(ic)%x)
       crt_dk(ic) = crt_k_tab(ic)%x(2) - crt_k_tab(ic)%x(1)
     end do fill_coordinates
     !
     crt_np = nd
     crt_dx = dx
     !
   end subroutine fill_cartesian_tables
   !
   subroutine fill_vp_integrals(tMin,tMax,dt,a0,evaluateF)
     real(ark), intent(in)             :: tMin       ! Initial time, time step 0
     real(ark), intent(in)             :: tMax       ! Final time
     real(ark), intent(in)             :: dt         ! Time step
     real(ark), intent(in)             :: a0(3)      ! Initial vector-potential
!    external                          :: evaluateF  ! Function for evaluating electric field as a function
!                                                    ! of time; takes real(ark) time as an argument
     !
     interface 
       function evaluateF(time) result(f)
         use accuracy
         real(ark), intent(in) :: time
         real(ark)             :: f(3)
       end function evaluateF
     end interface
     !
     integer(ik) :: it, alloc
     real(ark)   :: f(3)
     !
     lastTimeStep = 2+int((tMax-tMin+dt)/dt,kind=ik) ! "2" is the safety margin for round-off error;
                                                     ! if we ever hit the "lastTimeStep", we are in trouble.
     !
     allocate (time(0:lastTimeStep),timeA(3,0:lastTimeStep),timeAint(3,0:lastTimeStep), &
               timeA2(0:lastTimeStep),stat=alloc)
     if (alloc/=0) stop 'flux%fill_vp_integrals - allocation failure'
     !
     time      (0) = tMin
     timeA   (:,0) = a0
     timeAint(:,0) = 0
     timeA2    (0) = 0
     time_integral: do it=1,lastTimeStep
       f       (:)    = evaluateF(time(it-1))
       time      (it) = time      (it-1) + dt
       timeA   (:,it) = timeA   (:,it-1) - dt*f(:)
       timeAint(:,it) = timeAint(:,it-1) + dt*timeA(:,it-1)
       timeA2    (it) = timeA2    (it-1) + dt*sum(timeA(:,it-1)**2)
     end do time_integral
     !
     if (verbose>=0) then
       write (out,"('flux:        Time step first ',i15,  ' last ',i15  )")            0,             lastTimeStep
       write (out,"('flux:             Time first ',f15.6,' last ',f15.6)") time      (0), time      (lastTimeStep)
       write (out,"('flux: Vector-potential first ',g15.7,' last ',g15.7)") timeA   (1,0), timeA   (1,lastTimeStep)
       write (out,"('flux:                        ',g15.7,'      ',g15.7)") timeA   (2,0), timeA   (2,lastTimeStep)
       write (out,"('flux:                        ',g15.7,'      ',g15.7)") timeA   (3,0), timeA   (3,lastTimeStep)
       write (out,"('flux: VP time integral first ',g15.7,' last ',g15.7)") timeAint(1,0), timeAint(1,lastTimeStep)
       write (out,"('flux:                        ',g15.7,'      ',g15.7)") timeAint(2,0), timeAint(2,lastTimeStep)
       write (out,"('flux:                        ',g15.7,'      ',g15.7)") timeAint(3,0), timeAint(3,lastTimeStep)
       write (out,"('flux:    VP^2 integral first ',g15.7,' last ',g15.7)") timeA2    (0), timeA2    (lastTimeStep)
     end if
     write (out,"('Photoelectron spectra will be calculated at time step ',i15,' time ',f16.6)") &
            lastTimeStep, time(lastTimeStep)
   end subroutine fill_vp_integrals
   !
   !  Prepare table of p sample offsets from the central momentum
   !
   subroutine subsample_p(dims,samples)
     integer(ik) :: dims       ! Dimensionality of space
     integer(ik) :: samples(:) ! Number of samples along each dimension
     !
     integer(ik) :: id, ns
     integer(ik) :: is, ist(dims)
     real(rk)    :: poly_x(maxval(samples),dims)  ! 1D polynomial roots
     real(rk)    :: poly_w(maxval(samples),dims)  ! 1D polynomial weights
     integer(ik) :: alloc
     !
     if (dims/=size(samples)) stop 'flux%subsample_p - bad number of entries in samples(:)'
     if (any(samples<=0)) stop 'flux%subsample_p - bad sample counts in samples(:)'
     !
     p_samples = product(samples)
     if (verbose>=0) then
       write (out,"('Will use ',i6,' samples for each final momentum')") p_samples
       if (p_samples>1) then
         write (out,"(/'WARNING: momentum averaging removes phase information in the final spectrum'/)") 
       end if
     end if
     ! 
     allocate (p_offset(dims,p_samples),p_weight(p_samples),stat=alloc)
     if (alloc/=0) stop 'flux%subsample_p - allocation failed'
     !
     !  Fill the necessary polynomial roots and integration weights
     !  Since we know nothing about the integrand, go for Gauss-Legendre.
     !  MathGetQuadrature returns the roots in the [-1:+1] range; we need them
     !  in the [-0.5:0.5] range, so rescale the roots at the end.
     !
     fill_chebychev: do id=1,dims
       ns = samples(id)
       call MathGetQuadrature('Legendre',ns,poly_x(1:ns,id),poly_w(1:ns,id))
       !
       !  Fold in the weight function in the integration weights
       !
       poly_x(1:ns,id) = 0.5_rk * poly_x(1:ns,id)
       poly_w(1:ns,id) = 0.5_rk * poly_w(1:ns,id)
     end do fill_chebychev
     !
     ist = 1
     fill_samples: do is=1,p_samples
       p_weight(is) = 1._rk
       fill_dims: do id=1,dims
         p_offset(id,is) = poly_x(ist(id),id)
         p_weight   (is) = poly_w(ist(id),id) * p_weight(is)
       end do fill_dims
       if (verbose>=1) then
         write (out,"('Sample ',i6,' with weight ',f14.10,' at ',3f14.10)") is, p_weight(is), p_offset(:,is)
       end if
       advance_state: do id=1,dims
         ist(id) = ist(id) + 1
         if (ist(id)<=samples(id)) exit advance_state
         ist(id) = 1  ! And we'll move the next coordinate ...
       end do advance_state
     end do fill_samples
     if (verbose>=1) then
       write (out,*) ' sum of the weights = ', sum(p_weight)
     end if
   end subroutine subsample_p
   !
   !  Contribution from the upper/lower "lid" surfaces 
   !
   subroutine surface_cylindrical_lid (sample,f_psi,f_scr,sensor,next,aField,wgt)
     integer(ik), intent(in) :: sample ! Sub-sample index to use for the Bessel and Struve tables
     integer(ik), intent(in) :: f_psi  ! Wavefunction
     integer(ik), intent(in) :: f_scr  ! Field to update with the surface integral contribution
     integer(ik), intent(in) :: sensor ! Position of the sensor plane
     integer(ik), intent(in) :: next   ! Position of the next plane along the outward-looking normal
     real(rk), intent(in)    :: aField ! Vector-potential along the cylindrical "Z" axis.
     real(rk), intent(in)    :: wgt    ! Weight of this sensor
     !
     complex(rk), allocatable :: psi_sn(:,:)  ! Wavefunction within the sensor plane and the "next" plane
     complex(rk), allocatable :: int_drho(:)  ! Buffer for the inner part of the integral, which depends
                                              ! only on the value of krho
     complex(rk), allocatable :: int_phase(:) ! Buffer for the outer part of the integral, which depends
                                              ! only on the momentum
     integer(ik)              :: ikrho        ! Index of krho
     integer(ik)              :: ikz          ! Index of kz
     real(rk)                 :: kin_k        ! Instantaneous kinetic momentum 
     real(rk)                 :: zs           ! Sensor Z position
     !
     allocate (psi_sn(cyl_nrho,2),int_drho(cyl_nrho),int_phase(cyl_nz))
     !
     call FieldFetchSlice(f_psi,1_ik,sensor,psi_sn(:,1))
     call FieldFetchSlice(f_psi,1_ik,next,  psi_sn(:,2))
     !
     !  Evaluate the derivative of psi_sn along the normal vector.
     !  Remember that the wavefunction table actually contains 
     !  chi(z,rho), where Psi(rho,z) = exp(i m phi) chi(z,rho) / Sqrt(2*pi*rho)
     !  Additionally, pre-multiply by the volume element for the int_drho
     !  integral below.
     !
     psi_sn(:,1) = (cyl_drho/cyl_dz)*(psi_sn(:,2)-psi_sn(:,1))/sqrt(cyl_rho_tab)
     ! write (out,*) ' Lid gradient '
     ! write (out,*) psi_sn(:,1)
     !
     !  For the last point, we should only keep half: the surface passes through
     !  the middle of the volume element. 
     !
     psi_sn(flux_rhopos,1) = 0.5_rk * psi_sn(flux_rhopos,1)
     psi_sn(flux_rhopos+1:,1) = 0
     !
     !  The inner integral:
     !
     !    rho_max
     !      /                               d
     !      |  d rho J_0(krho*rho) * rho * ---  Psi(rho,z)
     !      /                              d z
     !      0
     !
     !  where rho_max is the radius of the sensor surface lid.
     !  The integral of J_0 * rho over each surface patch is
     !  already tabulated in "bessel_tab", so this whole integral
     !  becomes just a dot product.
     !
     !$omp parallel do default(none) private(ikrho) shared(int_drho,bessel_tab,flux_rhopos,psi_sn,cyl_nrho,sample)
     fill_inner: do ikrho=1,cyl_nrho
       int_drho(ikrho) = dot_product(bessel_tab(ikrho,1:flux_rhopos,sample),psi_sn(1:flux_rhopos,1))
     end do fill_inner
     !$omp end parallel do
     ! write (out,*) ' Bessel transform '
     ! write (out,*) int_drho
     !
     !  The outer part of the integral. All kz pointing away from the sensor surface should
     !  be excluded: electrons leaving to the detector in the sensor direction must have
     !  asymptotic momentum pointing towards the sensor. The kinetic momentum could still
     !  point the other way; this is simply the free electron oscillation in the field.
     !
     !  We have to be careful to only include kz=0 once (otherwise, we get very strange
     !  interference patterns from data on the two sensors); in theory, these electrons 
     !  could never cross the sensor, and should be discarded anyways ...
     !
     zs = cyl_z_tab(sensor)
     int_phase = 0
     fill_outer: do ikz=1,cyl_nz
       if (next>sensor .and. cyl_kz_tab(ikz)< 0) cycle fill_outer
       if (next<sensor .and. cyl_kz_tab(ikz)>=0) cycle fill_outer
       !
       kin_k = cyl_kz_tab(ikz)-electron_q*aField
       int_phase(ikz) = exp((0,-1)*zs*kin_k)
     end do fill_outer
     ! write (out,*) ' Phase correction '
     ! write (out,*) int_phase
     !
     !  The rest of the operation is a rank-1 update of the target matrix
     !
     int_phase = wgt * sqrt2pi * int_phase
     call FieldGER(f_scr,int_phase,int_drho,limits=p_limits(:,1:2))
     !
     deallocate (psi_sn,int_drho,int_phase)
   end subroutine surface_cylindrical_lid
   !
   !  Contribution from the "side" (constant rho) surface.
   !
   subroutine surface_cylindrical_side(sample,f_psi,f_scr,sensor,next,aField,wgt)
     integer(ik), intent(in) :: sample ! Sub-sample index to use for the Bessel and Struve tables
     integer(ik), intent(in) :: f_psi  ! Wavefunction
     integer(ik), intent(in) :: f_scr  ! Field to update with the surface integral contribution
     integer(ik), intent(in) :: sensor ! Position of the sensor "plane" (ring)
     integer(ik), intent(in) :: next   ! Position of the next "plane" along the outward-looking normal
     real(rk), intent(in)    :: aField ! Vector-potential along the cylindrical "Z" axis.
     real(rk), intent(in)    :: wgt    ! Weight of this sensor
     !
     complex(rk), allocatable :: psi_sn(:,:)  ! Wavefunction within the sensor plane and the "next" plane
     complex(rk), allocatable :: int_phase(:) ! Buffer for the outer part of the integral, which depends
                                              ! only on the value of krho
     real(rk)                 :: rhos         ! Position of the sensor plane
     real(rk)                 :: rhon         ! Position of the sensor plane
     !
     allocate (psi_sn(cyl_nz,2),int_phase(cyl_nrho))
     !
     call FieldFetchSlice(f_psi,2_ik,sensor,psi_sn(:,1))
     call FieldFetchSlice(f_psi,2_ik,next,  psi_sn(:,2))
     !
     rhos = cyl_rho_tab(sensor)
     rhon = cyl_rho_tab(next)
     !
     !  Evaluate the derivative of psi_sn along the normal vector.
     !  Remember that the wavefunction table actually contains 
     !  chi(z,rho), where Psi(rho,z) = exp(i m phi) chi(z,rho) / Sqrt(2*pi*rho)
     !
     psi_sn(:,1) = (psi_sn(:,2)/sqrt(rhon)-psi_sn(:,1)/sqrt(rhos)) / cyl_drho
     ! write (out,*) ' side gradient '
     ! write (out,*) psi_sn(:,1)
     !
     !  Multiply by the extra phase factor: exp(I*q*A(t)*z)
     !
     psi_sn(:,1) = psi_sn(:,1) * exp((0,1)*electron_q*aField*cyl_z_tab(:))
     ! write (out,*) ' extra phase '
     ! write (out,*) psi_sn(:,1)
     !
     !  For the first and the last points, we should only keep half: 
     !  the surface passes through the middle of the volume element.
     !  All elements beyond that should be zero, or FFT will get unhappy.
     !
     psi_sn(flux_zminpos,1) = 0.5_rk * psi_sn(flux_zminpos,1)
     psi_sn(flux_zmaxpos,1) = 0.5_rk * psi_sn(flux_zmaxpos,1)
     psi_sn(:flux_zminpos-1,1) = 0
     psi_sn(flux_zmaxpos+1:,1) = 0
     !
     !  The inner integral:
     !
     !    zmax
     !      /                                    d
     !      |  dz Exp(-I kz z) Exp(I q A(t) z) ----- Psi
     !      /                                  d rho
     !    zmin
     !
     !  Once we have multiplied the derivative with the A-dependent phase
     !  factor, the rest is simply an FFT. We need to be careful with
     !  phase factors, to match the FFTW's definition of the transform
     !  to ours, since the result will be accumulated coherently over
     !  the entire simulation.
     !  We never had to worry about these phase factors in FieldFFT (which 
     !  was was for the final observable!), so that phase factors we get
     !  here do not have to match FieldFFT. However, the amplitudes must
     !  of couse match.
     !
     call fftw_1d(cyl_nz,1_ik,psi_sn(:,1:1),.false.)
     psi_sn(:,1) = psi_sn(:,1) * exp((0,-1)*cyl_z_tab(1)*cyl_kz_tab(:)) * cyl_dz
     ! write (out,*) ' FFTW '
     ! write (out,*) psi_sn(:,1)
     !
     !  The krho-dependent part of the integral:
     !
     !    sqrt(2*pi)*rhos*0.5*(J_0(krho*rho)-I*H_0(krho*rho))
     !
     !  We will substitute the average of the rho*J_0(krho*rho) over the 
     !  spatial grid volume belonging to the boundary, but otherwise this 
     !  is a trivial contribution.
     !
     int_phase(:) = bessel_tab(:,sensor,sample) - (0,1) * struve_tab(:,sensor,sample)
     int_phase(:) = wgt * sqrt2pi * 0.5_rk * int_phase(:)
     !
     ! write (out,*) ' int_phase '
     ! write (out,*) int_phase
     !
     !  The rest of the operation is a rank-1 update of the target matrix
     !
     call FieldGER(f_scr,psi_sn(:,1),int_phase,limits=p_limits(:,1:2))
     !
     deallocate (psi_sn,int_phase)
   end subroutine surface_cylindrical_side
   !
   !  Contribution from 3D Cartesian-product surface
   !
   subroutine surface_cartesian(ic,f_psi,f_scr,sensor,next,aField,wgt)
     integer(ik), intent(in) :: ic        ! Axis perperndicular to the sensor surface (1=X, etc)
     integer(ik), intent(in) :: f_psi     ! Wavefunction
     integer(ik), intent(in) :: f_scr     ! Field to update with the surface integral contribution
     integer(ik), intent(in) :: sensor    ! Position of the sensor "plane" (ring)
     integer(ik), intent(in) :: next      ! Position of the next "plane" along the outward-looking normal
     real(rk), intent(in)    :: aField(3) ! Vector-potential
     real(rk), intent(in)    :: wgt       ! Weight of this sensor
     !
     integer(ik)              :: nic(2)         ! Cartesian dimensions within the sensor plane
     integer(ik)              :: np(2)          ! Number of grid points along each direction
     integer(ik)              :: alloc          ! Allocation status
     integer(ik)              :: j
     complex(rk), allocatable :: psi_sn(:,:,:)  ! Wavefunction within the sensor plane and the "next" plane
     complex(rk), allocatable :: ac1(:), ac2(:) ! Phase factors due to vector-potential components along
                                                ! the "first" and "second" plane directions; also used to
                                                ! shift FFT origin to the grid centre from the edge.
     complex(rk), allocatable :: ac3(:)         ! Buffer for the outer part of the integral, which depends
                                                ! only on the value of k along the normal of the sensor plane
     real(rk)                 :: rs             ! Position of the sensor plane
     !
     !  We'll need to know the axes corresponding to the sensor plane, and
     !  a few associated parameters
     !
     select case (ic)
       case default; stop 'flux%surface_cartesian - bad ic parameter'
       case (1); nic = (/ 2, 3 /)
       case (2); nic = (/ 1, 3 /)
       case (3); nic = (/ 1, 2 /)
     end select
     np = crt_np(nic)
     !
     allocate (psi_sn(np(1),np(2),2),ac1(np(1)),ac2(np(2)), &
               ac3(crt_np(ic)),stat=alloc)
     if (alloc/=0) stop 'flux%surface_cartesian - allocation failure'
     !
     call FieldFetchSlice(f_psi,ic,sensor,psi_sn(:,:,1))
     call FieldFetchSlice(f_psi,ic,next,  psi_sn(:,:,2))
     !
     rs = crt_r_tab(ic)%x(sensor)
     !
     !  Evaluate the derivative of psi_sn along the normal vector.
     !
     psi_sn(:,:,1) = (psi_sn(:,:,2)-psi_sn(:,:,1)) / crt_dx(ic)
     ! write (out,*) ' side gradient '
     ! write (out,*) psi_sn(:,1)
     !
     !  Multiply by the extra phase factor: exp(I*q*(A.r))
     !  Fortunately, it factorizes for the two directions within the plane,
     !  so that the number of complex exponentials remains small.
     !
     ac1 = exp((0,1)*electron_q*aField(nic(1))*crt_r_tab(nic(1))%x)
     ac2 = exp((0,1)*electron_q*aField(nic(2))*crt_r_tab(nic(2))%x)
     !
     add_spatial_vp_planewave: do j=1,np(2)
       psi_sn(:,j,1) = psi_sn(:,j,1) * ac1(:) * ac2(j)
     end do add_spatial_vp_planewave
     ! write (out,*) ' extra phase '
     ! write (out,*) psi_sn(:,1)
     !
     !  For edge points, we should only keep half:
     !  the surface passes through the middle of the volume element.
     !  All elements beyond that should be zero, or FFT will get unhappy.
     !
     psi_sn(  flux_pos(1,nic(1)),:,1) = 0.5_rk * psi_sn(  flux_pos(1,nic(1)),:,1)
     psi_sn(  flux_pos(2,nic(1)),:,1) = 0.5_rk * psi_sn(  flux_pos(2,nic(1)),:,1)
     psi_sn(:,flux_pos(1,nic(2)),  1) = 0.5_rk * psi_sn(:,flux_pos(1,nic(2)),  1)
     psi_sn(:,flux_pos(2,nic(2)),  1) = 0.5_rk * psi_sn(:,flux_pos(2,nic(2)),  1)
     !
     psi_sn(  :flux_pos(1,nic(1))-1 ,:,1) = 0
     psi_sn(   flux_pos(2,nic(1))+1:,:,1) = 0
     psi_sn(:,:flux_pos(1,nic(2))-1 ,  1) = 0
     psi_sn(:, flux_pos(2,nic(2))+1:,  1) = 0
     !
     !  The inner integral:
     !
     !      /                                          d
     !      | dr1 dr2 Exp(-I k . r) Exp(I q A(t) . r) --- Psi
     !      /                                         d z
     !
     !  Once we have multiplied the derivative with the A-dependent phase
     !  factor, the rest is simply an FFT. We need to be careful with
     !  phase factors, to match the FFTW's definition of the transform
     !  to ours, since the result will be accumulated coherently over
     !  the entire simulation.
     !  We never had to worry about these phase factors in FieldFFT (which 
     !  was for the final observable!), so that phase factors we get
     !  here do not have to match FieldFFT. However, the amplitudes must
     !  of couse match.
     !
     call TimerStart('FluxTimeStep: Sense: FFT')
     call fftw_2d(np(1),np(2),psi_sn(:,:,1),.false.,good_plan=.true.)
     call TimerStop('FluxTimeStep: Sense: FFT')
     !
     ac1 = exp((0,-1)*crt_r_tab(nic(1))%x(1)*crt_k_tab(nic(1))%x(:)) * crt_dx(nic(1))
     ac2 = exp((0,-1)*crt_r_tab(nic(2))%x(1)*crt_k_tab(nic(2))%x(:)) * crt_dx(nic(2))
     !
     origin_shift: do j=1,np(2)
       psi_sn(:,j,1) = psi_sn(:,j,1) * ac1(:) * ac2(j)
     end do origin_shift
     ! write (out,*) ' FFTW '
     ! write (out,*) psi_sn(:,1)
     !
     !  The part of the integral along the sensor normal:
     !
     !    exp(-I (k-q A)_z z0)
     !
     ac3(:) = wgt * exp((0,-1)*rs*(crt_k_tab(ic)%x(:)-electron_q*aField(ic)))
     !
     ! write (out,*) ' int_phase '
     ! write (out,*) int_phase
     !
     !  The rest of the operation is a rank-2+1 update of the target tensor
     !
     call FieldGER(f_scr,ic,psi_sn(:,:,1),ac3,limits=p_limits(:,1:3))
     !
     deallocate (psi_sn,ac1,ac2,ac3)
     !
   end subroutine surface_cartesian
   !
 end module flux
