!
!  Selection of complex absorbing potentials.
!
!  This module implements:
!  1. general complex absorbing potentials with no linear term, described in:
!     B. Poirier and T. Carrington, JCP 119, 77 (2003)
!  2. "transmission-free" absorbind potentials of:
!     D.E. Manolopoulos, JCP 117, 9552 (2002)
!     We'll be using an approximate version of the potential
!
!  Moiseyev's RF-CAP is non-local, and can't be accessed through this interface.
!
!  The two CAPs have rather different properties - see the original publications.
!
!  The Manolopoulos CAP and the numerically stable way of including them in
!  leap-frog propagation has been suggested by Michael Spanner
!
!  At the moment, only the Manolopoulos version is accessible through the "new"
!  interface. The "correct" way of using CAPsplitPotential is (quoting Michael's
!  e-mail of 28th September 2009):
!
!    At some point in the program, I set a field (f_abs_pot) for the
!    absorbing boundary using:
!         call FieldInit(f_abs_pot,AbsorbPot,mask=f_psi)
!    before I can use the absorbing boundary field
!    
!    performTimeStep: Demontrates how I use the absorbing boundary in
!                     conjunction with simple leap-frog algorithm.
!    
!                     I found that a stable way to add the absorbing boundary
!                     is to use a split-operator-type strategy.  The unitary
!                     Hamiltonian evolution coming from the system Hamiltonian
!                     and the laser field are applied using leap-frog.  Then
!                     the absorbing boundary is applied
!                     using the simplest method:
!                             Psi(t+dt) = exp(-dt*AbsrobingPotential) * Psi(t)
!                     However, for consistency, this needs to be done to both
!                     time-components of the leap-frog algoriothm, thus in the
!                     code I do:
!                             call FieldMul(f_abs_pot, f_psi )
!                             call FieldMul(f_abs_pot, f_psi2)
!                     Although multiplying the same operator on both timesteps
!                     looks weird, I've checked against 1D
!                     Fourier-split-opterator methods, where
!                     adding the Abs.Pot. directly into the main propagation
!                     algorithm is stable, and the results agree.
!    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        subroutine performTimeStep(f_pot,dt,eval,norm)
!          integer(ik), intent(in) :: f_pot  ! Field index containing the
!                                            ! potential
!          real(rk), intent(in)    :: dt     ! Desired time step
!          real(rk), intent(out)   :: eval   ! Total energy of f_psi
!                                            ! wavefunction
!          real(rk), intent(out)   :: norm   ! Initial norm of the propagated
!                                            ! W.F.
!          real(rk)                :: norm0  ! Temp for nomalization
!                                            ! calculations
!          integer(ik)   :: temp   ! temp
!    
!          ! Apply Hamiltonian
!          call QMHpsi(1._rk,f_pot,f_psi2,f_Hpsi)
!    
!          !  Propagate the wavefunction at T=t-dt
!          call FieldAXPY(cmplx(0._rk,-2._rk,kind=rk)*dt,f_Hpsi,f_psi)
!    
!          !  Apply the absorbing boundary
!          norm0 = FieldConjgIntegrate(f_psi,f_psi)
!          call FieldMul(f_abs_pot, f_psi )
!          call FieldMul(f_abs_pot, f_psi2)
!          norm  = FieldConjgIntegrate(f_psi,f_psi)
!          AbsorbedNorm = AbsorbedNorm + (norm0-norm)
!    
!          !  Swap field for the wavefunctions
!          temp = f_psi; f_psi = f_psi2 ; f_psi2 = temp
!    
!        end subroutine performTimeStep
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  module caps
    use accuracy
    use multigrid
    implicit none
    private
    public CAPsetUpPolynomial1D, CAPpolynomial1D  ! Legacy entry points, don't use.
    public CAPsetUp, CAPsplitPotential, CAPgetZone
    !
    integer(ik), parameter :: verbose = 0
    !
    !  General parameters
    !
    integer(ik), save :: zone_count   ! Number of absorbing zones in the box - could be up to 6 for 3D
    integer(ik), save :: zone_axis(6) ! Cartesian axis perpendicular to the zone
    real(rk), save    :: zone_r(2,6)  ! Beginning and end of each zone. The end would normally 
                                      ! correspond to the outer boundary of the grid.
    !
    !  Parameters for the Poirier's CAP
    !
    integer(ik), save  :: po_order    ! Polynomial order. Only 6, 10, and 14 are implemented
    integer(ik), save  :: po_axis     ! Axis perpendicular to which CAP is to be set up
    real(rk), save     :: po_en       ! Characteristic absorption energy
    real(rk), save     :: po_xmin     ! Beginning of the CAP zone
    real(rk), save     :: po_xmax     ! End of the CAP zone (normally coincides with cell boundary)
    real(rk), save     :: po_ac(2:14) ! Coefficients of f(z), eq. 11 / Table II
    !
    !  Parameters for the Manolopoulos CAP
    !
    real(rk), save      :: ma_scale
    real(rk), parameter :: ma_c = 2.622057554292119810464840_rk ! The position of the CAP singularity
    real(rk), parameter :: ma_a = 1._rk - 16._rk*ma_c**(-3)
    real(rk), parameter :: ma_b = (1._rk - 17._rk*ma_c**(-3)) * ma_c**(-2)
    !
    !  Routines
    !
    contains
    !
    subroutine CAPsetUpPolynomial1D(order_,axis_,en_,xmin_,xmax_)
      integer(ik), intent(in) :: order_ ! Polynomial order. Only 6, 10, and 14 are implemented
      integer(ik), intent(in) :: axis_  ! Axis perpendicular to which CAP is to be set up
      real(rk), intent(in)    :: en_    ! Characteristic absorption energy
      real(rk), intent(in)    :: xmin_  ! Beginning of the CAP zone
      real(rk), intent(in)    :: xmax_  ! End of the CAP zone (normally coincides with cell boundary)
      !
      po_axis  = axis_
      po_en    = en_
      po_xmin  = xmin_
      po_xmax  = xmax_
      po_order = order_
      select case (po_order)
        case default
          write (out,"('caps%CAPsetUpPolynomial1D: potential of order ',i6,' is not implemented')") po_order
          stop 'caps%CAPsetUpPolynomial1D - bad order'
        case  (6_ik)
          po_ac(2: 6) = (/ 3.94_rk, -7.38_rk,11.12_rk, -9.25_rk, 5.05_rk /)
        case (10_ik) 
          po_ac(2:10) = (/ 4.48_rk,-11.03_rk,23.22_rk,-35.94_rk,46.44_rk, -46.04_rk, 36.09_rk, -19.16_rk,  6.51_rk /)
        case (14_ik) 
          po_ac(2:14) = (/ 4.72_rk,-12.95_rk,31.02_rk,-58.48_rk,96.36_rk,-134.79_rk,164.07_rk,-169.26_rk,148.35_rk, &
                        -105.66_rk,59.97_rk,-23.92_rk,6.02_rk /)
      end select
    end subroutine CAPsetUpPolynomial1D
    !
    function CAPpolynomial1D(r) result (v)
      real(rk), intent(in) :: r(:) ! Coordinates - try to make it work for both 2D and 3D
      complex(rk)          :: v    ! Potential
      !
      real(rk)    :: x, z, fz
      complex(rk) :: wz
      integer(ik) :: ip
      !
      x = r(po_axis)
      if ( (po_xmin<po_xmax .and. (x<po_xmin .or. x>po_xmax)) .or. &
           (po_xmin>po_xmax .and. (x>po_xmin .or. x<po_xmax)) ) then
        !
        !  Not inside the absorbing zone
        !
        v = 0
        return
      end if
      !
      !  We are inside the absorbing zone.
      !
      z  = (x-po_xmin)/(po_xmax-po_xmin)
      !
      !  Evaluate auxiliary real polynomial f(z) - eq. 11
      !
      fz = po_ac(po_order)
      real_polynomial: do ip=po_order-1,2,-1
        fz = fz*z + po_ac(ip)
      end do real_polynomial
      fz = fz*z**2
      !
      !  Evaluate renormalized absorbing potential w(z) - eq. 13
      !
      wz = cmplx(2*fz,fz**2,kind=rk)
      !
      !  Evaluate problem-scale absorbing potential - eq. 1
      !
      v  = (0._rk,-1._rk)*po_en*wz
      ! write (out,"(' r = ',3f10.5,' v = ',2g16.7)") r, v
    end function CAPpolynomial1D
    !
    !  The JWKB reflection-free CAP shape, in renormalized coordinates
    !
    function maPotential(x) result(v)
      real(rk), intent(in) :: x ! Dimensionless penetration coordinate
      real(rk)             :: v ! The potential
      real(rk)             :: xt
      !
      xt = ma_c * min(max(0._rk,x),1._rk-spacing(100._rk))
      v  = ma_a * xt - ma_b * xt**3 + 4._rk * (ma_c-xt)**(-2) - 4._rk * (ma_c+xt)**(-2)
    end function maPotential
    !
    !  This "potential" should be used in a split-operator like propagation scheme 
    !  (suggested by Michael Spanner). Using it in any other way will likely cause
    !  no end of trouble.
    !
    function CAPsplitPotential(r) result(v)
      real(rk), intent(in) :: r(*) ! Coordinates of the grid point
      complex(rk)          :: v    ! exp(-dt*vcap(r)), where vcap is the imaginary part
                                   ! of the absorbing potential.
      integer(ik) :: ic, zone
      real(rk)    :: vs, x
      !
      vs = 0._rk
      check_zones: do zone=1,zone_count
        ic = zone_axis(zone)
        x  = (r(ic) - zone_r(1,zone))/(zone_r(2,zone)-zone_r(1,zone))
        if (x>0) vs = vs + maPotential(x)
      end do check_zones
      !
      v = 1._rk
      if (vs/=0._rk) v = exp(-ma_scale*vs)
    end function CAPsplitPotential
    !
    subroutine CAPsetUp(type,dt,axis,kmin,delta,mass)
      character(len=*), intent(in)      :: type  ! Type of the CAP to set up. Possible values are:
                                                 !  'Manolopoulos' 
      real(rk), intent(in)              :: dt    ! Time step used in propagation
      integer(ik), intent(in), optional :: axis  ! Cartesian absorption direction. Possible values are:
                                                 !   1 (-1) - Upper (lower) boundary of X
                                                 !   2 (-2) - Upper (lower) boundary of Y
                                                 !   3 (-3) - Upper (lower) boundary of Z
                                                 !  11      - All boundaries in Y and Z directions
                                                 !  12      - All boundaries in X and Z directions
                                                 !  13      - All boundaries in X and Y directions
                                                 !   0      - All boundaries. This is the default
      !
      !  The following parameters apply to Manolopoulos CAP 
      !
      real(rk), intent(in), optional    :: kmin  ! Minimal momentum to be absorbed
      real(rk), intent(in), optional    :: delta ! JWKB scaling constant. Will affect the 
                                                 ! efficiency with with kmin momentum is absorbed
                                                 ! Default is 0.2, corresponding to 1% efficiency
      real(rk), intent(in), optional    :: mass  ! Particle mass. Default is 1 (electron)
      !
      integer(ik) :: axis_
      real(rk)    :: mass_, delta_
      !
      axis_ = 0 ; if (present(axis)) axis_ = axis
      select case (type)
        case default
          write (out,"('CAP type ',a,' is not recognized')") trim(type)
          stop 'caps%CAPsetUp - bad CAP type'
        case ('Manolopoulos')
          if (.not.present(kmin)) then
            write (out,"('kmin must be specified for the Manolopoulos CAP')") 
            stop 'caps%CAPsetUp - missing required argument'
          end if
          mass_  = 1.0_rk ; if (present(mass )) mass_  = mass
          delta_ = 0.2_rk ; if (present(delta)) delta_ = delta
          call setupManolopoulos(axis_,dt,kmin,delta_,mass_)
      end select
    end subroutine CAPsetUp
    !
    subroutine setupManolopoulos(axis,dt,kmin,delta,mass)
      integer(ik), intent(in) :: axis  
      real(rk), intent(in)    :: dt    
      real(rk), intent(in)    :: kmin  
      real(rk), intent(in)    :: delta 
      real(rk), intent(in)    :: mass  
      !
      integer(ik) :: ndim     ! Dimensionality of the grid
      integer(ik) :: ic
      real(rk)    :: width    ! Width of the absorbing potential
      real(rk)    :: ext(2,3) ! Grid extent in each direction
      real(rk)    :: dx(3)    ! Grid spacing in each direction
      !
      ndim  = FieldComponentCount()
      if (ndim<2 .or. ndim>3) stop 'caps%setupManolopoulos - unsupported grid dimensionality'
      width    = ma_c / ( 2._rk * delta * kmin )
      ma_scale = dt * 0.5_rk * kmin**2 / mass
      if (verbose>=0) then
        write (out,"('cap%setupManolopoulos: Width of the CAP     = ',g14.6,' Bohr')") width
        write (out,"('cap%setupManolopoulos: Magnitude of the CAP = ',g14.6,' Hartree')") ma_scale/dt
        write (out,"('cap%setupManolopoulos: Magnitude of the CAP = ',g14.6,' Hartree-Jiffies')") ma_scale
      end if
      dx(1:ndim) = FieldGridSpacing()
      get_extent: do ic=1,ndim
        call FieldGridExtent(ic,ext(:,ic))
        if ( (ext(2,ic)-ext(1,ic)) <= 2*width ) then
          stop 'cap%setupManolopoulos: Absorbing potential reaches the origin'
        end if
        if ( width/dx(ic) < 5 ) then
          stop 'cap%setupManolopoulos: Absorbing potential has insufficient grid coverage'
        end if
      end do get_extent
      zone_count = 0
      if (              any(axis==(/0, 1,12,13/)) ) then  ! upper X (3D) or Z (cyl 2D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: positive X')")
          if (ndim==2) write (out,"('cap%setupManolopoulos: Adding CAP zone: positive Z')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 1
        zone_r(1,zone_count) = ext(2,1) - width
        zone_r(2,zone_count) = ext(2,1)
      end if
      if (              any(axis==(/0,-1,12,13/)) ) then  ! lower X (3D) or Z (cyl 2D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: negative X')")
          if (ndim==2) write (out,"('cap%setupManolopoulos: Adding CAP zone: negative Z')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 1
        zone_r(1,zone_count) = ext(1,1) + width
        zone_r(2,zone_count) = ext(1,1)
      end if
      if (              any(axis==(/0, 2,11,13/)) ) then  ! upper Y (3D) or Rho (cyl 2D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: positive Y')")
          if (ndim==2) write (out,"('cap%setupManolopoulos: Adding CAP zone: positive Rho')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 2
        zone_r(1,zone_count) = ext(2,2) - width
        zone_r(2,zone_count) = ext(2,2)
      end if
      if (ndim>=3 .and. any(axis==(/0,-2,11,13/)) ) then  ! lower Y (3D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: negative Y')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 2
        zone_r(1,zone_count) = ext(1,2) + width
        zone_r(2,zone_count) = ext(1,2)
      end if
      if (ndim>=3 .and. any(axis==(/0, 3,11,12/)) ) then  ! upper Z (3D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: positive Z')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 3
        zone_r(1,zone_count) = ext(2,3) - width
        zone_r(2,zone_count) = ext(2,3)
      end if
      if (ndim>=3 .and. any(axis==(/0,-3,11,12/)) ) then  ! lower Z (3D)
        if (verbose>=0) then
          if (ndim==3) write (out,"('cap%setupManolopoulos: Adding CAP zone: negative Z')")
        end if
        zone_count = zone_count + 1
        zone_axis(zone_count) = 3
        zone_r(1,zone_count) = ext(1,3) + width
        zone_r(2,zone_count) = ext(1,3)
      end if
      if (verbose>=1) then
        write (out,*) ' zone_count = ', zone_count
        write (out,*) ' zone_axis  = ', zone_axis
        write (out,*) ' zone_r     = ', zone_r
      end if
      if (zone_count==0) then
        write (out,"('Manolopoulos CAP axis specification ',i9,' did not generate any zones.')") axis
        stop 'caps%setupManolopoulos - invalid axis specification'
      end if
    end subroutine setupManolopoulos
    !
    function CAPgetZone(axis,rmin) result(OK)
      integer(ik), intent(in) :: axis  ! Can be +/- 1, 2, and 3.
      real(rk), intent(out)   :: rmin  ! Beginning of the zone; coordinate closest to the origin
      logical                 :: OK    ! .true. if zone found
      !
      integer(ik)             :: iz
      !
      OK = .false.
      scan_zones: do iz=1,zone_count
        if (zone_axis(iz)/=abs(axis)) cycle scan_zones
        if (axis<0 .and. zone_r(2,iz)>zone_r(1,iz)) cycle scan_zones
        rmin = zone_r(1,iz)
        OK   = .true.
        return
      end do scan_zones
    end function CAPgetZone
  end module caps
