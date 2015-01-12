module sfa
!
!  Subroutines useful for doing Strong-Field Approximation calculations
!  For simplicity, we'll use a fixed field geometry:
!
!   Laser field propagation is towards positive Z direction. The field
!   is monochromatic, elliptically polarized:
!
!     E(z,t) = [ Ex*Cos(eta), Ey*Sin(eta), 0 ]
!     eta    = omega*(t-z/c)
!
!   The electric field is related to the vector-potential as:
!
!     E      = -dA/dt 
!     A      = [ (-Ex/omega) Sin[eta], (Ey/omega) Cos[eta], 0 ]
!
!   So that the magnetic field is:
!
!     B      = \Nabla\cross A
!            = [ (-Ey/c) Sin[eta], (Ex/c) Cos[eta], 0 ]
!
!   We are working in non-relativistic length gauge, including 
!   non-dipole terms through order (1/c). The Hamiltonian is:
!
!         1   2           q z
!    H = --- p  - q E.r - --- E.p
!        2 m              m c
!
!   The Volkov solution for this Hamiltonian is [see Chirila et at,
!   PRA 66, 063411 (2002)]
!
!                   I
!    psi_V = Exp[- ---- phi_V]
!                  hbar
!
!                         1  / t   2
!    phi_V = - r.pi(t) + --- |   pi (t') d t'
!                        2 m /
!                                   q                   2
!    pi(t) = p - q A(t) - [0,0,1] ----- ( 2 p.A(t) - q A (t) )
!                                 2 m c
!
!   Unternally, all calculations use SFA units, defined as:
!
!                         Atomic                     SFA
!    time:                t                        omega*t (aka phase)
!    electric field:      E0*e_s                     e_s
!    vector-potential:    (E0/omega)*a_s             a_s
!    momentum:            (E0/omega)*p_s             p_s
!    energy:              (E0**2/m*omega**2)*e_s     e_s
!
  use accuracy
  use math
  implicit none
  private
  public sfa_data
  public sfa_init, sfa_destroy
  public sfa_stationary_momentum, sfa_kvec
  public sfa_stationary_t1, sfa_phase
  public sfa_efield, electron_q, electron_m
  !
  integer(ik), parameter :: verbose    = 1_ik   ! Verbosity level
  real(rk), parameter    :: electron_q = -1._rk ! Electron charge
  real(rk), parameter    :: electron_m =  1._rk ! Electron mass
  !
  type sfa_data
    !
    real(rk) :: s_ex   ! Relative amplitude of the peak electric field in the X direction
    real(rk) :: s_ey   ! Relative amplitude of the peak electric field in the Y direction
    real(rk) :: s_ip   ! Ionization potential; SFA units
    real(rk) :: s_e0   ! Conversion factor: SFA electric field to a.u. electric field
    real(rk) :: s_p0   ! Conversion factor: SFA momentum to a.u. momentum
    real(rk) :: s_u0   ! Conversuon factor: SFA energy to Hartree
    real(rk) :: s_t0   ! Conversion factor: SFA time (phase) to a.u. time
    real(rk) :: s_mcm1 ! (electron mass times speed of light in SFA units)**-1
    !
    real(rk), allocatable :: s_grid_x(:)  ! Numerical integration grid, used for calculating phases
    real(rk), allocatable :: s_grid_w(:)
  end type sfa_data
  !
  contains
  !
  !  External interfaces
  !
  subroutine sfa_init(s,omega,ex,ey,ip,magnetic)
    type(sfa_data), intent(inout) :: s        ! Data structure describing current laser field
    real(rk), intent(in)          :: omega    ! Laser field frequency, atomic units
    real(rk), intent(in)          :: ex, ey   ! Laser electric field amplitudes along x and y directions
    real(rk), intent(in)          :: ip       ! Target ionization potential, Hartree
    logical, intent(in), optional :: magnetic ! Include magnetic (1/c) terms. The default is .true.
    !
    integer(ik), parameter :: grid_order = 1000_ik
    !
    s%s_e0    = sqrt(ex**2+ey**2)
    s%s_p0    = s%s_e0/omega
    s%s_u0    = s%s_p0**2/(electron_m)
    s%s_t0    = 1.0_rk/omega
    !
    s%s_ex    = ex/s%s_e0
    s%s_ey    = ey/s%s_e0
    s%s_ip    = ip/s%s_u0
    s%s_mcm1  = s%s_p0/(electron_m * vlight)
    !
    if (present(magnetic)) then
      if (.not.magnetic) s%s_mcm1 = 0._rk
    end if
    !
    allocate (s%s_grid_x(grid_order),s%s_grid_w(grid_order))
    call MathGetQuadrature('Legendre',grid_order,s%s_grid_x,s%s_grid_w)
    !
    if (verbose>=1) then
      write (out,"()")
      write (out,"(' SFA electric field unit = ',g25.14,' a.u.')") s%s_e0
      write (out,"('  SFA momentum/V.P. unit = ',g25.14,' a.u.')") s%s_p0
      write (out,"('         SFA energy unit = ',g25.14,' Hartree')") s%s_u0
      write (out,"('           SFA time unit = ',g25.14,' a.u.')") s%s_t0
      write (out,"()")
      write (out,"(' peak electric field (X) = ',g25.14,' SFA units')") s%s_ex
      write (out,"(' peak electric field (Y) = ',g25.14,' SFA units')") s%s_ey
      write (out,"('    ionization potential = ',g25.14,' SFA units')") s%s_ip
      write (out,"('            1/(m*vlight) = ',g25.14,' SFA units')") s%s_mcm1
      write (out,"()")
    end if
  end subroutine sfa_init
  !
  subroutine sfa_destroy(s)
    type(sfa_data), intent(inout) :: s      ! Data structure describing current laser field
    !
    if (allocated(s%s_grid_x)) deallocate (s%s_grid_x)
    if (allocated(s%s_grid_w)) deallocate (s%s_grid_w)
  end subroutine sfa_destroy
  !
  !  Calculate stationary momentum, given times of ionization and recollision
  !
  function sfa_stationary_momentum(s,ph2,ph1) result(p)
    type(sfa_data), intent(inout) :: s      ! Data structure describing current laser field
    complex(rk), intent(in)       :: ph2    ! Field phase at the time of recollision, radian
    complex(rk), intent(in)       :: ph1    ! Field phase at the time of ionization, radian
    complex(rk)                   :: p(3)   ! Stationary momentum, in SFA units
    !
    !  Stationary solutions for the momemtum are OK in the order O(1/c); this is
    !  the same accuracy as the Hamiltonian, so there is no need to try any harder than this.
    !
    if (abs(ph2-ph1)<=sqrt(spacing(max(abs(ph2),abs(ph1))))) then
      !
      !  If ionization and recollision times become similar, we have to use the L'Hopital rule
      !
      p(1) = s%s_ex * (-sin(ph1) - 0.5_rk * cos(ph1)*(ph2-ph1))
      p(2) = s%s_ey * ( cos(ph1) - 0.5_rk * sin(ph1)*(ph2-ph1))
      p(3) = s%s_ex**2 * (-4.0_rk * sin(ph1)**2 - 2.0_rk * sin(2*ph1)*(ph2-ph1)) &
           + s%s_ey**2 * (-4.0_rk * cos(ph1)**2 + 2.0_rk * sin(2*ph1)*(ph2-ph1))
    else
      p(1) = s%s_ex * (cos(ph2)-cos(ph1))/(ph2-ph1)
      p(2) = s%s_ey * (sin(ph2)-sin(ph1))/(ph2-ph1)
      p(3) = s%s_ex**2 * (-2._rk + (sin(2*ph2)-sin(2*ph1))/(ph2-ph1)) &
           + s%s_ey**2 * (-2._rk - (sin(2*ph2)-sin(2*ph1))/(ph2-ph1))
    end if
    !
    !  Include the electron charge
    !
    p(1:2) = electron_q    * p(1:2)
    p(3)   = electron_q**2 * p(3)
    !
    !  Add the p^2 term to the Z momentum, and apply the overall factor
    !
    p(3) = ( p(1)**2 + p(2)**2 + 0.125_rk * p(3) ) * s%s_mcm1
  end function sfa_stationary_momentum
  !
  !  Calculate instanteneous kinetic momentum, given stationary momentum and field phase
  !
  function sfa_kvec(s,ph,p) result(pi)
    type(sfa_data), intent(inout) :: s      ! Data structure describing current laser field
    complex(rk), intent(in)       :: ph     ! Field phase 
    complex(rk), intent(in)       :: p(3)   ! Stationary momentum, in SFA units
    complex(rk)                   :: pi(3)  ! Kinetic momentum
    !
    complex(rk) :: a(3)
    !
    a(1) = s%s_ex * (-sin(ph))
    a(2) = s%s_ey *   cos(ph)
    a(3) = 0
    ! write (out,*) ' a = ',a
    !
    !  Dipole contribution
    !
    pi = p - electron_q * a
    !
    !  Extra (1/c) part
    !
    pi(3) = pi(3) - 0.5_rk * electron_q * s%s_mcm1 * sum((2._rk*p-electron_q*a)*a)
  end function sfa_kvec
  !
  !  Find stationary ionization time closest to the desired recollision time
  !
  function sfa_stationary_t1(s,ph2,ph1_max,ph1) result(success)
    type(sfa_data), intent(inout) :: s       ! Data structure describing current laser field
    complex(rk), intent(in)       :: ph2     ! Field phase at recollistion
    complex(rk), intent(in)       :: ph1_max ! Start the search at Re(ph1_max)
    complex(rk), intent(out)      :: ph1     ! success==.True.:  Field phase at ionization; Re(ph1)<ph1_max
                                             ! success==.False.: Field phase at which the search was terminated
    logical                       :: success ! Matching ionization time was found 
    !
    integer(ik), parameter :: steps_real      = 200_ik  ! Number of steps in [0:2pi] for the initial search
    integer(ik), parameter :: steps_imaginary = 200_ik
    real(rk), parameter    :: dp              = 1e-9_rk ! Numerical differentiation step & convergence criterion
    integer(ik), parameter :: max_iterations  =1000_ik  ! Max number of refinement iterations
    real(rk), parameter    :: max_step        = 1.0_rk/steps_real
    real(rk)               :: test_disp       
    !
    integer(ik)            :: sr, si
    complex(rk)            :: cn, cg
    complex(rk)            :: phr
    !
    !  At the stationary ionization time, we must have:
    !
    !    IP + pi**2/(2*m) == 0
    !
    !  We will begin by performing a grid search for a starting guess.
    !  The criterion for terminating the search is the sign change in
    !  both the real and imaginary components of the stationary condition
    !
    test_disp = 5.0_rk * twopi / min(steps_real,steps_imaginary)
    success = .false.
    scan_real: do sr=1,steps_real
      ph1 = real(ph1_max,kind=rk) - (twopi*sr)/steps_real
      scan_imag: do si=0,steps_imaginary
        ph1 = cmplx(real(ph1,kind=rk),(twopi*si)/steps_imaginary,kind=rk)
        call stationary_condition_and_gradient(ph1,cn,cg)
        if ( test_disp*abs(cg)>=abs(cn) ) then
          !
          !  We are close to a stationary point, try refining it
          !
          if (refine_guess(ph1,phr)) then
            if (abs(phr-ph1_max)>10._rk*dp) then
              !
              !  This is a new solution; return
              !
              ph1 = phr
              success = .true.
              exit scan_real
            end if
          end if
          !
          !  Refinement produced the same solution we started the search from; keep trying
          !
        end if
      end do scan_imag
    end do scan_real
    !
    if (.not.success) then
      write (out,"('sfa_stationary_t1: Can''t find a solution for ph2 = ',2g20.11)") ph2
    end if
    !
    contains

    function stationary_condition(ph1) result(cn)
      complex(rk), intent(in) :: ph1 ! Ionization phase
      complex(rk)             :: cn  ! Current value of the stationary ionization condition
      !
      complex(rk) :: pst(3), pi(3)
      !
      pst = sfa_stationary_momentum(s,ph2,ph1)
      pi  = sfa_kvec(s,ph1,pst)
      cn  = s%s_ip + 0.5_rk*sum(pi**2)
      ! write (out,*) ' ph1 = ', ph1, ' cn = ',cn
    end function stationary_condition
    !
    !  Calculate gradient of the stationary condition by numerical differentiation
    !
    subroutine stationary_condition_and_gradient(ph1,cn,gr)
      complex(rk), intent(in)  :: ph1 ! Ionization phase
      complex(rk), intent(out) :: cn  ! Current value of the stationary ionization condition
      complex(rk), intent(out) :: gr  ! ,,, and its gradient
      !
      complex(rk) :: fp, fm
      !
      cn = stationary_condition(ph1)
      fm = stationary_condition(ph1-dp)
      fp = stationary_condition(ph1+dp)
      gr = (fp-fm)/(2*dp)
    end subroutine stationary_condition_and_gradient
    !
    !  Refine the initial guess for the stationary point, using simple gradient steps.
    !
    function refine_guess(ph1_guess,ph1) result(success)
      complex(rk), intent(in)  :: ph1_guess ! Initial guess
      complex(rk), intent(out) :: ph1       ! Refined guess; only valid if success == .true.
      logical                  :: success   ! Refinement converged
      !
      integer(ik) :: iter
      complex(rk) :: f0, fg
      complex(rk) :: delta
      !
      ph1 = ph1_guess
      refine: do iter=1,max_iterations
        call stationary_condition_and_gradient(ph1,f0,fg)
        if (abs(fg)==0.0_rk) then
          stop 'sfa%sfa_stationary_t1 - gradient vanished identically'
        end if
        delta = -f0/fg
        if (abs(delta)>=max_step) then
          delta = delta*max_step/abs(delta)
        end if
        if (verbose>=2) then
          write (out,"('Phase = ',2f14.7,' step = ',2g16.5,' residue = ',2g16.5,' grad = ',2g16.5)") ph1, delta, f0, fg
        end if
        ph1     = ph1 + delta
        if (abs(delta)<=dp) then
          !
          !  Refinement worked; we can leave now.
          !
          success = .true.
          return
        end if
      end do refine
      success = .false.
    end function refine_guess
  end function sfa_stationary_t1
  !
  !  Return SFA electric field at a given time
  !
  function sfa_efield(s,ph) result(e)
    type(sfa_data), intent(inout) :: s       ! Data structure describing current laser field
    complex(rk), intent(in)       :: ph      ! Phase when electric field is needed
    complex(rk)                   :: e(3)
    !
    e(1) = s%s_ex * cos(ph)
    e(2) = s%s_ey * sin(ph)
    e(3) = 0._rk
  end function sfa_efield
  !
  !  Although I can write an analytical expression for this integral, it is much
  !  simpler to take it numerically.
  !
  function sfa_phase(s,p,ph1,ph2) result(fi)
    type(sfa_data), intent(inout) :: s      ! Data structure describing current laser field
    complex(rk), intent(in)       :: p(3)   ! Momentum, in SFA units
    complex(rk), intent(in)       :: ph1    ! Initial field phase
    complex(rk), intent(in)       :: ph2    ! Final field phase
    complex(rk)                   :: fi     ! Phase accumulated in the continuum is exp(i*fi)
    !
    integer(ik) :: ipt
    complex(rk) :: a, b  ! Linear transformation from [-1:1] to [ph1:ph2]
    complex(rk) :: ph    ! Phase at an integration point
    complex(rk) :: val   ! Integrand at an integration point
    complex(rk) :: pi(3) ! Momentum at integration point
    !
    !  Sanity check: we'd like at least 100 points per laser cycle
    !
    if (100*abs(ph2-ph1)/twopi>size(s%s_grid_x)) then
      stop 'sfa%sfa_stationary_phase - integration order is insufficient to handle phase integration'
    end if
    !
    a = 0.5_rk*(ph2-ph1)
    b = 0.5_rk*(ph2+ph1)
    !
    fi = 0.0_rk
    phase_integrate: do ipt=1,size(s%s_grid_x)
      ph  = a*s%s_grid_x(ipt) + b
      pi  = sfa_kvec(s,p=p,ph=ph)
      val = s%s_ip + 0.5_rk * sum(pi**2)
      fi = fi + s%s_grid_w(ipt) * val
    end do phase_integrate
    fi = -a * fi
    !
    !  Convert to atomic units
    !
    fi = fi * (s%s_u0*s%s_t0)
  end function sfa_phase
end module sfa
