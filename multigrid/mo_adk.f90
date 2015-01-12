!
!  Multigrid test - MO-ADK ionization rates, using data from partial_waves_sp.f90
!
!  This is a literal implementation of Zhao et al, J Phys B 44, 035601 (2011).
!  The only significant change is in the definition of the factor Q(l,m) in eq. 4
!  Since we use Landau&Lifshitz conventions for the spherical harmonics, we have
!  to replace:
!
!    (-1)**((|m|+m)/2)
!
!  by:
!
!    (-1)**((|m|-m)/2) (0,1)**L
!
!
  module mo_adk
    use accuracy
    use timer
    use math
    use opendx
    implicit none
    private
    public run_mo_adk
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: dim_l     = 10_ik          ! Highest order or spherical harmonics in wf
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)            :: verbose = 0_ik
    real(rk)               :: match_r = 1.0_rk           ! Radius of the matching shere; used for debugging
    real(rk)               :: zc      = 1.0_rk           ! Charge of the cation, |e| units
    real(rk)               :: ip      = 0.5_rk           ! Ionization potential, atomic units
    real(rk)               :: f0      = 0.1_rk           ! Electric field intensity, atomic units
    integer(ik)            :: l_max   = 0                ! Max. order of YLM coefficients on input
    integer(ik)            :: m_limit = 100              ! Do not include tunneling contributions with m' higher than this
    complex(rk)            :: clm(0:dim_l,-dim_l:dim_l)  ! YLM coefficients
    integer(ik)            :: grid_theta = 200           ! Number of grid points in theta angle [0:Pi]
    integer(ik)            :: grid_phi   = 400           ! Number of grid points in phi angle [0:2Pi]
    character(len=100)     :: phase      = 'L&L'         ! Choice of the spherical harmonics phase.
                                                         ! 'L&L' is Landau and Lifshitz
                                                         ! 'CDL' is Chi-Dong Lin.
    character(len=100)     :: directions = 'grid'        ! Laser field polarization directions:
                                                         ! 'grid' = use a product grid
                                                         ! 'read' = read theta, phi pairs from input
    character(len=100)     :: out_file   = ' '           ! File for the angular distribution; blank
                                                         ! suppresses the output
    character(len=100)     :: sph_file   = ' '           ! Debug: output MO on grid
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    real(rk)               :: kappa                      ! Calculated from ip
    real(rk), allocatable  :: ion_grid(:,:)              ! Ionization data on the spherical product grid
    real(rk), allocatable  :: mo_grid (:,:)              ! Reconstructed orbital on the grid
    real(rk)               :: step_theta                 ! Calculated from grid_theta
    real(rk)               :: step_phi                   ! Calculated from grid_phi
    real(rk)               :: total_rate                 ! Total ionization rate
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /moadk/ verbose, &
                     match_r, zc, ip, f0, &
                     l_max, m_limit, &
                     directions, &
                     out_file, sph_file, &
                     grid_theta, grid_phi, &
                     phase, clm
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Eq. 4, corrected for the difference in phase factors
    !
    complex(rk) function qlm(l,m)
      integer(ik), intent(in) :: l, m
      integer(ik) :: am
      !
      am  = abs(m)
      qlm = sqrt((2*l+1)*MathFactorial(l+am)/(2*MathFactorial(l-am)))
      qlm = qlm * (-1)**((am-m)/2) * (0,1)**l
    end function qlm
    !
    !  Eq. 10 - Static ionization rate in a given direction
    !
    function static_rate(theta,phi,f,mov) result(w)
      real(rk), intent(in)            :: theta, phi ! Direction of the laser field polarization 
                                                    ! in the molecular frame. 
      real(rk), intent(in)            :: f          ! Electric field strength
      real(rk), intent(out), optional :: mov        ! Debug: Orbital value on grid 
      real(rk)                        :: w          ! MO-ADK ionization rate
      !
      real(rk)    :: euler(3)                       ! Euler angles for field-adapted coordinate system
      complex(rk) :: xlm(0:l_max,-l_max:l_max)      ! YLM coefficients in field-adapted C.S.
      complex(rk) :: dj (-l_max:l_max,-l_max:l_max) ! Rotation matrix for spherical harmonics
      complex(rk) :: bm(-l_max:l_max)               ! Axially-symmetric components of the wf in field direction
      real(rk)    :: pre, prem                      ! Pre-exponential factor in the ionization rate
      integer(ik) :: l, m, am
      !
      if (verbose>=1) then
        write (out,"('Calculating static rate for theta = ',f12.8,' phi = ',f12.8,' f = ',f12.8)") theta, phi, f
      end if
      !
      !  Rotate harmonics to align quantization axis with the field
      !
      euler(1) =  phi
      euler(2) = -theta
      euler(3) = 0._rk ! This angle is arbitrary
      rotate_harmonics: do l=0,l_max
        call MathYJMRotationMatrix(euler,2*l+1,dj(-l:l,-l:l))
        xlm(l,-l:l) = matmul(dj(-l:l,-l:l),clm(l,-l:l))
      end do rotate_harmonics
      !
      if (verbose>=2) then
        write (out,"(/'Rotated spherical harmonics: ')")
        call echo_angular_coefficients(xlm)
        write (out,"()")
      end if
      !
      !  Calculate axial components along field direction
      !
      axial_components: do m=-l_max,l_max
        bm(m) = 0
        axial_l: do l=abs(m),l_max
          bm(m) = bm(m) + qlm(l,m) * xlm(l,m)
        end do axial_l
      end do axial_components
      !
      !  Sanity check: bm(0) must be real
      !
      if (abs(imag(bm(0)))>spacing(1e6_rk*max(1.0_rk,abs(bm(0))))) then
        write (out,"('Fatal: bm(0) at phi= ',g14.6,' theta= ',g14.6,' is not real:',2(1x,g25.14))") &
               phi, theta, bm(0)
        stop 'mo_adk%static_rate - internal check failed'
      end if
      !
      if (verbose>=1) then
        write (out,"('On-axis amplitudes:')")
        write (out,"(1x,a4,2x,a20,1x,a20)") ' M ', ' Re(b(m)) ', ' Im(b(m)) ', &
                                            '---', '----------', '----------'
        print_bm: do m=-l_max,l_max
          write (out,"(1x,i4,2x,f20.12,1x,f20.12)") m, bm(m)
        end do print_bm
        write (out,"()")
      end if
      !
      !  Calculate pre-exponential part of the rate
      !
      if (verbose>=1) then
        write(out,"(1x,a5,1x,a20)") ' M ', ' Prefactor ', &
                                    '---', '-----------'
      end if
      pre = 0
      rate_preexponent: do m=-min(l_max,m_limit),min(l_max,m_limit)
        am   = abs(m)
        prem = + (abs(bm(m))**2/(2**am * MathFactorial(am))) &
               * (kappa**(2*zc/kappa-1))**(-1) &
               * (2*kappa**3/f)**(3*zc/kappa-am-1)
        pre  = pre + prem
        if (verbose>=1) then
          write(out,"(1x,i5,1x,g20.12)") m, prem
        end if
      end do rate_preexponent
      if (verbose>=1) then
        write(out,"(1x,a5,1x,g20.12/)") 'Total', pre
      end if
      !
      !  Add the Keldysh factor, and off we go
      !
      w = pre * exp(-2*kappa**3/(3*f))
      !
      !  Debugging: bm(0) is a rescaled MO value on grid; save it
      ! 
      if (present(mov)) then
        mov = real(bm(0),kind=rk) * match_r**(zc/kappa-1) * exp(-kappa*match_r) / sqrt2pi
      end if
    end function static_rate
    !
    !  Calculate ionization rates on spherical product grid. This
    !  gives us the total ionization yield, and can be used for plotting
    !
    subroutine calculate_total_rate
      integer(ik) :: itheta, iphi
      real(rk)    :: theta, phi, wgt, w, wtot, mov
      !
      step_phi   = twopi/grid_phi
      step_theta =    pi/grid_theta
      wtot       = 0._rk
      total_rate = 0._rk
      allocate (ion_grid(grid_theta,grid_phi),mo_grid(grid_theta,grid_phi))
      !$omp parallel default(none) private(iphi,phi,itheta,theta,wgt,w,mov) &
      !$omp&         shared(ion_grid,mo_grid,grid_theta,grid_phi,step_theta,step_phi,f0) &
      !$omp&         reduction(+:total_rate,wtot)
      wtot       = 0._rk
      total_rate = 0._rk
      !$omp do
      phi_points: do iphi=1,grid_phi  ! [0:2Pi]
        phi = step_phi*(iphi-0.5_rk)
        theta_points: do itheta=1,grid_theta
          theta = step_theta*(itheta-0.5_rk)
          wgt   = sin(theta) * step_phi * step_theta
          w     = static_rate(theta,phi,f0,mov)
          ion_grid(itheta,iphi) = w
          mo_grid (itheta,iphi) = mov
          total_rate = total_rate + wgt * w
          wtot       = wtot + wgt
        end do theta_points
      end do phi_points
      !$omp end do
      !$omp end parallel
      total_rate = total_rate / wtot
    end subroutine calculate_total_rate
    !
    !  Transform phases of spherical harmonics from CDL to L&L convention
    !
    subroutine adjust_phases
      integer(ik) :: l, m, am
      !
      adjust_l: do l=0,l_max
        adjust_m: do m=-l,l
          am = abs(m)
          clm(l,m) = clm(l,m) &
                   * ( (-1)**((am-m)/2) * (0,1)**l ) &
                   / ( (-1)**((am+m)/2) )
        end do adjust_m
      end do adjust_l
    end subroutine adjust_phases
    !
    subroutine echo_angular_coefficients(clm)
      complex(rk), intent(in) :: clm(0:l_max,-l_max:l_max)
      integer(ik) :: l, m
      !
      write (out,"((1x,a3,1x,a4,2x,a14,1x,a14))") &
             ' L ', ' M ', ' Prj. YLM ', ' ', &
             '---', '---', '----------', ' '
      loop_l2: do l=0,l_max
        loop_m2: do m=-l,l
          if (abs(clm(l,m))<1e-5_rk) cycle loop_m2
          write (out,"(1x,i3,1x,i4,2x,g14.7,1x,g14.7)") l, m, clm(l,m)
        end do loop_m2
      end do loop_l2
      write (out,"()")
    end subroutine echo_angular_coefficients
    !
    !  Problem driver
    !
    subroutine run_mo_adk
      integer(ik) :: info
      real(rk)    :: theta, phi, w
      !
      call TimerStart('MO ADK')
      call accuracyInitialize
      w = MathFactorial(80)
      w = MathLogFactorial(80)
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=moadk,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=moadk)
      write (out,"()")
      !
      kappa      = sqrt(2._rk*ip)
      !
      write (out,"('Long-range exponent (kappa) = ',g16.7)") kappa
      !
      if (m_limit<l_max) then
        write (out,"(/'WARNING: Tunneling contributions truncated at m''= ',i4/)") m_limit
      end if
      !
      if (phase=='CDL') call adjust_phases
      !
      write (out,"(/t5,'Long-range wavefunction expansion; L&L phase conventions.'/)")
      call echo_angular_coefficients(clm(0:l_max,-l_max:l_max))
      !
      write (out,"(/'All rates reported below are for static field of the peak intensity. To convert'" &
                //"/'to the cycle-average rates, multiply by the factor below.'/)")
      write (out,"(' fcycle = ',g18.9/)") sqrt(3*f0/(pi*kappa**3))
      !
      select case (directions)
        case default ; stop 'mo_adk%run_mo_adk - invalid directions parameter'
        case ('grid','GRID')
          call calculate_total_rate
          write (out,"('Total MO-ADK ionization rate = ',g16.7)") total_rate
          if (out_file/=' ') call odx_write_spherical(out_file,ion_grid)
          if (sph_file/=' ') call odx_write_spherical(sph_file, mo_grid)
        case ('read','READ')
          write (out,"(3x,a12,1x,a12,2x,a18)") ' theta ', ' phi ', ' rate ', &
                                               '-------', '-----', '------'
          read_loop: do
            read (input,*,iostat=info) theta, phi
            if (info/=0) exit read_loop
            w = static_rate(theta,phi,f0)
            write (out,"(1x,'#',1x,f12.8,1x,f12.8,2x,g18.9)") theta, phi, w
          end do read_loop
      end select
      !
      call TimerStop('MO ADK')
      call TimerReport
    end subroutine run_mo_adk

  end module mo_adk
!
  subroutine driver
    use mo_adk

    call run_mo_adk
  end subroutine driver
