!
!  Data interpolation from a spherical Lebedev grid to a uniform theta-phi product grid
!  We first project the data onto spherical harmonics, using the fact that Lebedev grids
!  integrate harmonics up to a given order exactly. We then evaluate the resulting fit
!  on the desired product grid.
!
!  The input is the &si namelist, immediately followed by the grid data, in the format:
!    X, Y, Z, W, F
!    ...
!
!  The output file is in (a human-readable flavour of) the OpenDX format.
!
  module spherical_interpolate
    use accuracy
    use fields
    use opendx
    use timer
    implicit none
    private
    public run_si
    !
    integer(ik), parameter :: iu_output   = 34
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 0                     ! Verbosity level
    integer(ik)        :: grid_theta      = 200                   ! Number of grid points in theta angle [0:Pi]
    integer(ik)        :: grid_phi        = 400                   ! Number of grid points in phi angle [0:2Pi]
    integer(ik)        :: l_max           = 4                     ! Max. L value to scan in the expansion
    integer(ik)        :: n_lebedev       = 38                    ! 9-th order Lebedev grid has 38 points
    character(len=256) :: out_file        = ' '                   ! Output file; standard output if blank
    character(len=256) :: harmonics_file  = 'fit_harmonics.dat'   ! Output file containing harmonics coefficients
                                                                  ! output disabled if blank
    character(len=256) :: comment         = ' '                   ! Arbitrary text comment to be added to the output
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    real(rk), allocatable    :: leb_grid (:,:)    ! Lebedev grid
                                                  ! Indices are: (X:Y:Z:W,IPT)
    real(rk), allocatable    :: leb_data (:)      ! Data points on the Lebedev grid
    complex(rk), allocatable :: prj_coeff(:,:)    ! Harmonic amplitudes from projection
    real(rk)                 :: step_phi          ! Product grid spacing
    real(rk)                 :: step_theta        ! Product grid spacing
    real(rk), allocatable    :: prod_grid(:,:,:)  ! Spherical product grid of a unit radius at zero
                                                  ! Indices are: (X:Y:Z:W,ITH,IPH)
    complex(rk), allocatable :: prod_ylm   (:,:)  ! Spherical harmonics on the product grid
    complex(rk), allocatable :: prod_data  (:,:)  ! Interpolated data on the product grid
    !
    real(rk), parameter      :: r_zero(3) = 0._rk ! This must be zero
    complex(rk), parameter   :: v_zero    = 0._rk ! This must be zero
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /si/ verbose, &
                  grid_theta, grid_phi, &
                  n_lebedev, l_max,     &
                  out_file, harmonics_file, &
                  comment
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Allocate all data fields used by the code
    !
    subroutine allocate_arrays
      integer(ik) :: alloc, nfit
      !
      nfit = (l_max+1)**2
      allocate (leb_grid(4,n_lebedev), leb_data(n_lebedev),      &
                prj_coeff(-l_max:l_max,0:l_max),                 &
                prod_grid(4,grid_theta,grid_phi),                &
                prod_ylm   (grid_theta,grid_phi),                &
                prod_data  (grid_theta,grid_phi),                &
         stat=alloc)
      if (alloc/=0) then
        write (out,"('Memory allocation failed, error code = ',i8)") alloc
        stop 'spherical_interpolate%allocate_arrays - allocation failed'
      end if
    end subroutine allocate_arrays
    !
    !  Read in data to interpolate
    !
    subroutine read_points
      integer(ik) :: info, ipt
      !
      slurp_points: do ipt=1,n_lebedev
        read(input,*,iostat=info) leb_grid(:,ipt), leb_data(ipt)
        if (info/=0) then
          write (out,"('Error ',i8,' reading data point ',i8)") info, ipt
          stop 'spherical_interpolate%read_point - I/O problem'
        end if
      end do slurp_points
    end subroutine read_points
    !
    !  Assume that the input grid is a valid Lebedev grid of a sufficient order
    !  to handle all needed harmonics; perform the projection.
    !
    subroutine harmonics_project
      integer(ik) :: lv, mv, ipt
      complex(rk) :: hv, acc
      !
      !  Perform projection, and fill the linear system matrix for the fit
      !
      prj_coeff = 0
      loop_l1: do lv=0,l_max
        loop_m1: do mv=-lv,lv
          call FLharmonicsSetParameters(lv,mv,r_zero,0._rk,safe_max)
          acc = 0
          loop_pt: do ipt=1,n_lebedev
            hv = FLharmonics(leb_grid(1:3,ipt),v_zero)
            acc = acc + leb_grid(4,ipt) * conjg(hv) * leb_data(ipt)
          end do loop_pt
          prj_coeff(mv,lv) = acc
        end do loop_m1
      end do loop_l1
      !
      if (verbose>=0) then
        write (out,"((1x,a3,1x,a4,2x,a14,1x,a14))") &
               ' L ', ' M ', ' Prj. YLM ', ' ', &
               '---', '---', '----------', ' '
        loop_l2: do lv=0,l_max
          loop_m2: do mv=-lv,lv
            write (out,"(1x,i3,1x,i4,2x,g14.7,1x,g14.7)") lv, mv, prj_coeff(mv,lv)
          end do loop_m2
        end do loop_l2
        write (out,"()")
        write (out,"('Total yield from grid = ',g14.7)") sum(leb_data*leb_grid(4,:))
        write (out,"('Total yield from fit  = ',g14.7)") real(prj_coeff(0,0),kind=rk) * 2._rk * sqrtpi
        write (out,"()")
      end if
    end subroutine harmonics_project
    !
    !  Save harmonics coefficients, for use by other codes 
    !
    subroutine harmonics_dump
      integer(ik) :: lv, mv
      integer(ik) :: ionization_l_max
      real(rk)    :: ionization_e_peak, ionization_duration
      !
      namelist /ionization_fit/ ionization_l_max, ionization_e_peak, ionization_duration
      !
      if (harmonics_file==' ') return
      !
      open(iu_output,file=trim(harmonics_file),form='formatted',status='replace')
      !
      ionization_l_max    = l_max
      ionization_e_peak   = -1._rk ! To be filled by hand
      ionization_duration = -1._rk ! To be filled by hand
      write (iu_output,nml=ionization_fit)
      !
      loop_l2: do lv=0,l_max
        loop_m2: do mv=-lv,lv
          write (iu_output,"(1x,i3,1x,i4,2x,g24.15,1x,g24.15)") lv, mv, prj_coeff(mv,lv)
        end do loop_m2
      end do loop_l2
      close (iu_output)
    end subroutine harmonics_dump
    !
    !  Prepare basic angular grid, for a sphere of a unit radius.
    !
    subroutine initialize_angular_grid
      integer(ik) :: itheta, iphi
      real(rk)    :: theta, phi
      real(rk)    :: x, y, z, w
      !
      call TimerStart('Initialize grid')
      step_phi   = twopi/grid_phi
      step_theta =    pi/grid_theta
      !
      write (out,"(/' Angular product grid:')")
      write (out,"(' phi [0:2Pi]: pts = ',i5,' step = ',f12.9,' Rad')") grid_phi, step_phi
      write (out,"(' theta[0:Pi]: pts = ',i5,' step = ',f12.9,' Rad')") grid_theta, step_theta
      !
      !$omp parallel do private(iphi,phi,itheta,theta,x,y,z,w)
      phi_points: do iphi=1,grid_phi  ! [0:2Pi]
        phi = step_phi*(iphi-0.5_rk)
        theta_points: do itheta=1,grid_theta
          theta = step_theta*(itheta-0.5_rk)
          x = sin(theta)*cos(phi)
          y = sin(theta)*sin(phi)
          z = cos(theta)
          w = sin(theta) * step_phi * step_theta
          prod_grid(:,itheta,iphi) = (/ x, y, z, w /)
        end do theta_points
      end do phi_points
      !$omp end parallel do
      call TimerStop('Initialize grid')
    end subroutine initialize_angular_grid
    !
    !  Fill data field with a spherical harmonic
    !
    subroutine fill_harmonic(l,m,ylm)
      integer(ik), intent(in)  :: l, m     ! Angular parameters
      complex(rk), intent(out) :: ylm(:,:) ! Desired spherical harmonic
      !
      integer(ik) itheta, iphi
      !
      call TimerStart('Spherical harmonic')
      call FLharmonicsSetParameters(l,m,r_zero,0._rk,safe_max)
      !$omp parallel do private(iphi,itheta)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          ylm(itheta,iphi) = FLharmonics(prod_grid(1:3,itheta,iphi),v_zero)
        end do theta_loop
      end do phi_loop
      !$omp end parallel do
      call TimerStop('Spherical harmonic')
    end subroutine fill_harmonic
    !
    !  Add two kets: |dst> = |dst> + wgt*|src>
    !
    subroutine add_kets(wgt,src,dst)
      complex(rk), intent(in)    :: wgt      ! Weight of the source
      complex(rk), intent(in)    :: src(:,:) ! Vector to add
      complex(rk), intent(inout) :: dst(:,:) ! Destination vector
      !
      integer(ik) itheta, iphi
      !
      call TimerStart('Add two kets')
      !$omp parallel do private (iphi,itheta)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          dst(itheta,iphi) = dst(itheta,iphi) + wgt*src(itheta,iphi)
        end do theta_loop
      end do phi_loop
      !$omp end parallel do
      call TimerStop('Add two kets')
    end subroutine add_kets
    !
    !  Interpolate data
    !
    subroutine interpolate
      integer(ik) :: lv, mv
      !
      prod_data = 0
      loop_l: do lv=0,l_max
        loop_m: do mv=-lv,lv
          call fill_harmonic(lv,mv,prod_ylm)
          call add_kets(prj_coeff(mv,lv),src=prod_ylm,dst=prod_data)
        end do loop_m
      end do loop_l
      !
      if (verbose>=0) then
        write (out,"()")
        write (out,"('Maximum interpolated yield = ',g14.7)") maxval(real(prod_data,kind=rk))
        write (out,"('Minimum interpolated yield = ',g14.7)") minval(real(prod_data,kind=rk))
        write (out,"()")
      end if
      !
      !  Sanity check: interpolated data must be real
      !
      if (maxval(abs(aimag(prod_data)))>=100._rk*spacing(maxval(abs(real(prod_data,kind=rk))))) then
        stop 'spherical_interpolate%interpolate - the result is not real!'
      end if
      !
    end subroutine interpolate
    !
    !  Problem driver
    !
    subroutine run_si
      integer(ik)        :: info
      !
      call TimerStart('Spherical Interpolation')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=si,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=si)
      write (out,"()")
      !
      call allocate_arrays
      !
      call read_points
      !
      call harmonics_project
      !
      call harmonics_dump
      !
      call initialize_angular_grid
      !
      call interpolate
      !
      call odx_write_spherical(out_file,real(prod_data,kind=rk),trim(comment))
      !
      call TimerStop('Spherical Interpolation')
      call TimerReport
    end subroutine run_si

  end module spherical_interpolate
!
  subroutine driver
    use spherical_interpolate

    call run_si
  end subroutine driver
