!
!  Multigrid test - expansion of a 3D orbital in spherical partial waves
!                   (this is an essential step for an MO-ADK calculation)
!                   Spherical-grid version.
!
  module partial_waves_sp
    use accuracy
    use fields
    use math
    use opendx
    use timer
    use import_gamess
    implicit none
    private
    public run_partial
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_fields   = 2                    ! Max. number of fields used in the code
    integer(ik), parameter :: iu_out       = 30                   ! I/O unit used by report_mo_adk
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 0                     ! Verbosity level
    integer(ik)        :: grid_theta      = 200                   ! Number of grid points in theta angle [0:Pi]
    integer(ik)        :: grid_phi        = 400                   ! Number of grid points in phi angle [0:2Pi]
    character(len=100) :: orbital_file    = 'orbital.dat'         ! Name of the file containing MO coefficients
    integer(ik)        :: orbital_index   = 1                     ! Orbital components to load from "orbital_file"
    real(rk)           :: euler_angles(3) = 0.0_rk                ! Euler angles used to rotate the molecule
    real(rk)           :: orbital_ip      = 0.5_rk                ! Ionization potential - for the long-range tail
    real(rk)           :: ion_charge      = 1.0_rk                ! Charge of the ion core - for the long-range tail
    integer(ik)        :: wave_l_max      = 4                     ! Max. L value to scan in the expansion
    real(rk)           :: wave_r_min      = 5._rk                 ! Min. distance cut-off for the long-range tail
    real(rk)           :: wave_r_max      = 10._rk                ! Max. distance cut-off for the long-range tail
    real(rk)           :: wave_r_step     = 0.2_rk                ! Steps to take in the radial grid
    real(rk)           :: wave_r0(3)      = (/ 0._rk, 0._rk, 0._rk /)
    real(rk)           :: report_fraction = 1e-3_rk               ! Smallest relative expansion coefficients to report
    logical            :: mo_adk          = .false.               ! Generate input data for mo_adk.f90
    character(len=100) :: mo_adk_out      = "('mo_adk_',f0.3,'.inp')" ! Name template for the mo_adk.f90 input file.
                                                                  ! Actual file name will be generated using 
                                                                  !    write(fmt=mo_adk_out) r
                                                                  ! where r is the matching-sphere radius
    character(len=100) :: mo_adk_dx       = "('dist_',f0.3,'.dx')" ! Name template for the mo_adk.f90 visialization file
    real(rk)           :: f0              = 0.057_rk              ! Only usef to generate mo_adk.f90 input
    logical            :: plot_sphere     = .false.
    character(len=100) :: sphere_dx       = "('sph_',f0.3,'.dx')" ! Debug: Dump file for orbital sampled on matched sphere
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    real(rk)                 :: wave_r            ! Current radius of the matching sphere
    real(rk), allocatable    :: base_grid(:,:,:)  ! Reference spherical grid of a unit radius at zero
                                                  ! Indices are: (X:Y:Z:W,ITH,IPH)
    real(rk), allocatable    :: grid (:,:,:)      ! Current working grid
    complex(rk), allocatable :: bra    (:,:)      ! Wavefunctions
    complex(rk), allocatable :: ket    (:,:)      ! 
    complex(rk), allocatable :: ylm_wgt(:,:)      ! Table of spherical wave weights
    real(rk)                 :: rotmat (3,3)      ! Rotation matrix for the molecule
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /wavedata/ verbose, &
                       grid_theta, grid_phi, &
                       euler_angles, &
                       orbital_file, orbital_index, &
                       orbital_ip, ion_charge, &
                       wave_r0, wave_l_max, wave_r_min, wave_r_max, wave_r_step, &
                       report_fraction, &
                       mo_adk, mo_adk_out, mo_adk_dx, f0, &
                       plot_sphere, sphere_dx
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Allocate all data fields used by the code
    !
    subroutine allocate_arrays
      integer(ik) :: alloc
      !
      write (out,"(/' Need about ',f12.3,' Mbytes of memory'/)") &
             rk_bytes*(4+4+2+2)*real(grid_theta,kind=rk)*grid_phi/(1024._rk**2)
      !
      allocate (base_grid(4,grid_theta,grid_phi), &
                     grid(4,grid_theta,grid_phi), &
                      bra(  grid_theta,grid_phi), &
                      ket(  grid_theta,grid_phi), &
            ylm_wgt(-wave_l_max:wave_l_max,0:wave_l_max), stat=alloc)
      if (alloc/=0) then
        write (out,"('Memory allocation failed, error code = ',i8)") alloc
        stop 'partial_waves_sp%allocate_arrays - allocation failed'
      end if
    end subroutine allocate_arrays
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
    !  Prepare basic angular grid, for a sphere of a unit radius.
    !
    subroutine initialize_angular_grid
      integer(ik) :: itheta, iphi
      real(rk)    :: step_theta, step_phi
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
          base_grid(:,itheta,iphi) = (/ x, y, z, w /)
        end do theta_points
      end do phi_points
      !$omp end parallel do
      call TimerStop('Initialize grid')
    end subroutine initialize_angular_grid
    !
    !  Scale reference angular grid to reflect the desired expansion centre and radius
    !
    subroutine scale_angular_grid
      integer(ik) itheta, iphi
      !
      call TimerStart('Rescale grid')
      !$omp parallel do private(iphi,itheta)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          grid(1:3,itheta,iphi) = wave_r*base_grid(1:3,itheta,iphi) + wave_r0
          grid(  4,itheta,iphi) =        base_grid(  4,itheta,iphi)
        end do theta_loop
      end do phi_loop
      !$omp end parallel do
      call TimerStop('Rescale grid')
    end subroutine scale_angular_grid
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
      call FLharmonicsSetParameters(l=l,m=m,r0=wave_r0,rmin=0._rk,rmax=safe_max)
      !$omp parallel do private(iphi,itheta)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          ylm(itheta,iphi) = FLharmonics(grid(1:3,itheta,iphi),(0._rk,0._rk))
        end do theta_loop
      end do phi_loop
      !$omp end parallel do
      call TimerStop('Spherical harmonic')
    end subroutine fill_harmonic
    !
    !  Calculate a braket
    !
    function calculate_braket(bra,ket) result(v)
      complex(rk), intent(in) :: bra(:,:) ! Left-hand vector
      complex(rk), intent(in) :: ket(:,:) ! Right-hand vector
      complex(rk)             :: v        ! Overlap integral aka braket
      !
      integer(ik) itheta, iphi
      !
      call TimerStart('Evaluate braket')
      v = 0
      !$omp parallel do private(iphi,itheta) reduction(+:v)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          v = v + grid(4,itheta,iphi)*conjg(bra(itheta,iphi))*ket(itheta,iphi)
        end do theta_loop
      end do phi_loop
      !$omp end parallel do
      call TimerStop('Evaluate braket')
    end function calculate_braket
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
    !  Check quality of an angular grid by integrating spherical harmonics
    !
    subroutine check_angular_grid
      integer(ik) :: bra_vl, bra_vm, ket_vl, ket_vm
      complex(rk) :: ovr
      real(rk)    :: exact, error, max_error
      integer(ik) :: e_bra_vl, e_bra_vm, e_ket_vl, e_ket_vm
      !
      call TimerStart('Check angular grid')
      wave_r = 1.0_rk
      call scale_angular_grid
      !
      max_error = -1
      write (out,"(/' Testing angular grid quality')")
      bra_l: do bra_vl=0,wave_l_max
        bra_m: do bra_vm=-bra_vl,bra_vl
          call fill_harmonic(bra_vl,bra_vm,bra)
          ket_l: do ket_vl=0,wave_l_max
            ket_m: do ket_vm=-ket_vl,ket_vl
              call fill_harmonic(ket_vl,ket_vm,ket)
              ovr = calculate_braket(bra,ket)
              !
              exact = 0.0_rk
              if (bra_vl==ket_vl .and. bra_vm==ket_vm) exact = 1.0_rk
              error = abs(ovr-exact)
              if (error>max_error) then
                max_error = error
                e_bra_vl = bra_vl ; e_bra_vm = bra_vm
                e_ket_vl = ket_vl ; e_ket_vm = ket_vm
              end if
              !
              if (verbose>=1) then
                write (out,"(' <YLM(',i2,',',i3,')|YLM(',i2,',',i3,')> = ',2g20.12)") &
                       bra_vl, bra_vm, ket_vl, ket_vm, ovr
              end if
            end do ket_m
          end do ket_l
        end do bra_m
      end do bra_l
      !
      if (verbose>=1) write (out,"()")
      write (out,"(' Largest error was ',g13.6,' for  <YLM(',i2,',',i3,')|YLM(',i2,',',i3,')>')") &
             max_error, e_bra_vl, e_bra_vm, e_ket_vm, e_ket_vm
      write (out,"()")
      call TimerStop('Check angular grid')
    end subroutine check_angular_grid
    !
    !  Load target Gamess orbital
    !
    subroutine load_gamess_mos(wf)
      complex(rk), intent(out) :: wf(:,:) ! Data field for the orbital
      !
      integer(ik)              :: nnuc, iat, alloc
      real(rk), allocatable    :: xyzq(:,:)
      logical, save            :: first_run = .true.
      real(rk), allocatable    :: tmp_coord(:,:,:,:)
      complex(rk), allocatable :: tmp_data(:,:,:,:)
      !
      !  The interface to GAMESS import routine is for a 3D grid, so this
      !  will get a little messy. 
      !
      allocate (tmp_coord(3,grid_theta,grid_phi,1),tmp_data(grid_theta,grid_phi,1,1),stat=alloc)
      if (alloc/=0) then
        write (out,"('load_gamess_mos: Error ',i8,' allocating temporary I/O buffers')") alloc
        stop 'partial_waves_sp%load_gamess_mos - allocation problem'
      end if
      !
      tmp_coord(:,:,:,1) = grid(1:3,:,:)
      call gamess_load_orbitals(orbital_file,(/orbital_index/),(/1_ik/),tmp_coord,tmp_data,rot=rotmat)
      wf = tmp_data(:,:,1,1)
      !
      deallocate (tmp_coord,tmp_data,stat=alloc)
      if (alloc/=0) then
        write (out,"('load_gamess_mos: Error ',i8,' releasing temporary I/O buffers')") alloc
        stop 'partial_waves_sp%load_gamess_mos - deallocation problem'
      end if
      !
      if (first_run .or. verbose>=1) then
        first_run = .false.
        write (out,"('Loaded orbital ',i4,' from ',a)") orbital_index, trim(orbital_file)
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
      end if
    end subroutine load_gamess_mos
    !
    !  Divide the orbital by the expected long-range radial function
    !
    subroutine set_long_range_tail
      real(rk) :: kappa
      !
      kappa = sqrt(2._rk*orbital_ip)
      write (out,"(' Long-range exponent (kappa) = ',f18.12)") kappa
      call FLasymptoticSetParameters(z=ion_charge,kappa=kappa,r0=wave_r0)
    end subroutine set_long_range_tail
    !
    function get_asymptotic_factor() result(v)
      real(rk) :: v
      !
      v = real(FLasymptotic(grid(1:3,1,1),(0._rk,0._rk)),kind=rk)
    end function get_asymptotic_factor
    !
    subroutine plot_mo_sphere(r,data)
      real(rk), intent(in)    :: r         ! Shere radius
      complex(rk), intent(in) :: data(:,:) ! Orbital data
      character(len=100)      :: file
      !
      write (file,fmt=sphere_dx) r
      call odx_write_spherical(trim(file),real(data,kind=rk))
    end subroutine plot_mo_sphere
    !
    !  Descriptive message
    !
    subroutine simulation_banner
      write (out,"()")
      write (out,"('  ===================================================================')")
      write (out,"('  One-centre expansion of the long-range tail of a molecular orbital ')")
      write (out,"('  in terms of spherical harmonics:')")
      write (out,"('  ')")
      write (out,"('     Psi(r,t,p) = Rad(r) Sum Sum C(L,M) Y(L,M,t,p)')")
      write (out,"('                          L   M')")
      write (out,"('  ')")
      write (out,"('  where the asymptotic radial dependence is in the form:')")
      write (out,"('  ')")
      write (out,"('     Rad(r) = r**(Z/kappa-1) exp(-kappa r)')")
      write (out,"('  ')")
      write (out,"('  where Z is the charge of the remaining ion, kappa in (2*Ip)**0.5,')")
      write (out,"('  and the phase of spherical harmonics is follows Landau&Lifshitz')")
      write (out,"('  conventions. The coefficients C(L,M) are calculated as:')")
      write (out,"('  ')")
      write (out,"('              /      *                      /  *')")
      write (out,"('     C(L,M) = | Psi Y (L,M) / Rad(r) d r  / | Y (L,M) Y(L,M) d r')")
      write (out,"('              /                             /')")
      write (out,"('  ')")
      write (out,"('  and both 2D integrals are taken over a spherical shell of a radius')")
      write (out,"('  between Rmin and Rmax.')")
      write (out,"('  ===================================================================')")
      write (out,"()")
    end subroutine simulation_banner
    !
    !  Report angular weights in a compact form
    !
    subroutine report_weights(radial)
      real(rk), intent(in) :: radial  ! Asymptotic radial wavefunction
      !
      integer(ik)  :: vl, vm
      real(rk)     :: threshold
      !
      threshold = report_fraction*maxval(abs(ylm_wgt))
      !
      write (out,"(/t5,'Significant contributions at R=',f12.5,' Rad(r)= ',g16.8/)") &
             wave_r, radial
      !
      write (out,"(t5,a4,1x,a5,3x,a16,1x,a16)") &
             ' L ', ' M ', '    Re(CLM)    ', '   Im(CLM)   ', &
             '---', '---', '   ---------   ', '  ---------  '
      !
      scan_l: do vl=0,wave_l_max
        scan_m: do vm=-vl,vl
          if (abs(ylm_wgt(vm,vl))<threshold) cycle scan_m
          write (out,"(t5,i4,1x,i5,3x,f16.6,1x,f16.6)") vl, vm, ylm_wgt(vm,vl)
        end do scan_m
      end do scan_l
      write (out,"()")
    end subroutine report_weights
    !
    !  Generate input for mo_adk.f90
    !  It would have been easier to use namelist here, but this produces a rather
    !  ugly output with most compilers. I'd rather make it easier for a human to 
    !  spot possible mistakes.
    !
    subroutine report_mo_adk(r)
      real(rk), intent(in) :: r
      character(len=100)   :: out_file, grid_file
      integer(ik)          :: vl, vm
      !
      write (out_file,fmt=mo_adk_out) r
      write (grid_file,fmt=mo_adk_dx) r
      open (iu_out,form='formatted',status='replace',file=trim(out_file))
      write (iu_out,"(' &MOADK'/"// &
                     "'  MATCH_R    = ',f12.8,', ! Bohr'/"// &
                     "'  VERBOSE    = 0,'/"// &
                     "'  ZC         = ',f12.8,','/"// &
                     "'  IP         = ',f12.8,', ! = ',f12.5,' eV'/"// &
                     "'  F0         = ',f12.8,', ! = ',g12.4,' W/cm^2'/"// &
                     "'  DIRECTIONS = ''GRID'','/"// &
                     "'  OUT_FILE   = ''',a,''','/"// &
                     "'  SPH_FILE   = '' '','/"// &
                     "'  GRID_THETA = 200,'/"// &
                     "'  GRID_PHI   = 400,'/"// &
                     "'  PHASE      = ''L&L'','/"// &
                     "'  L_MAX      = ',i4,',')") &
             r, ion_charge, orbital_ip, orbital_ip * h2ev, f0, f0**2 * au2wcm2, trim(grid_file), wave_l_max
      scan_l: do vl=0,wave_l_max
        scan_m: do vm=-vl,vl
          if (abs(ylm_wgt(vm,vl))<=spacing(100._rk)) cycle scan_m
          write (iu_out,"('  CLM(',i4,',',i5,') = (',f22.14,',',f22.14,'),')") vl, vm, ylm_wgt(vm,vl)
        end do scan_m
      end do scan_l
      write (iu_out,"(' /')")
      close (iu_out)
    end subroutine report_mo_adk
    !
    !  Problem driver
    !
    subroutine run_partial
      integer(ik)        :: info, vl, vm
      real(rk)           :: total_rho, residue_rho, asymptotic_factor
      complex(rk)        :: ovr, nrm
      !
      call TimerStart('Partial Waves')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=wavedata,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=wavedata)
      write (out,"()")
      !
      call TimerStart('Preliminaries')
      !
      call allocate_arrays
      !
      call initialize_rotation
      !
      call initialize_angular_grid
      !
      if (verbose>=0) call check_angular_grid
      !
      call set_long_range_tail
      !
      call TimerStop('Preliminaries')
      !
      !  Expand it!
      !
      call simulation_banner
      !
      call TimerStart('Expansion')
      wave_r = wave_r_min
      scan_spheres: do while (wave_r<=wave_r_max)
        !
        !  Perform angular integration for the current sphere
        !
        call scale_angular_grid
        !
        asymptotic_factor = get_asymptotic_factor()
        call load_gamess_mos(ket)
        if (plot_sphere) then
          call plot_mo_sphere(wave_r,ket)
        end if
        total_rho = real(calculate_braket(ket,ket),kind=rk)
        write (out,"(/' R= ',f12.8,' Bohr, radial density= ',g13.6,' Bohr^-1 = ',g13.6,' x asymptote')") &
               wave_r, total_rho * wave_r**2, total_rho * asymptotic_factor**2
        !
        ylm_wgt = 0
        scan_spherical_l: do vl=0,wave_l_max
          scan_spherical_m: do vm=-vl,vl
            call fill_harmonic(vl,vm,bra)
            ovr = calculate_braket(bra,ket)
            nrm = calculate_braket(bra,bra)
            ylm_wgt(vm,vl) = asymptotic_factor*ovr/nrm 
            !
            call add_kets(-ovr/nrm,bra,ket)
            !
            if (verbose>=1) then
              write (out,"(' L= ',i2,' M= ',i3,' C= ',2f18.8,' (ovr= ',2g12.5,' nrm= ',2g12.5,')')") &
                     vl, vm, ylm_wgt(vm,vl), ovr, nrm
            end if
          end do scan_spherical_m
        end do scan_spherical_l
        !
        residue_rho = real(calculate_braket(ket,ket),kind=rk)
        write (out,"(' R= ',f12.8,' Bohr: ',f10.4,' % of radial density is at L>',i2)") &
               wave_r, 100._rk*residue_rho/total_rho, wave_l_max
        call report_weights(1._rk/asymptotic_factor)
        if (mo_adk) then
          call report_mo_adk(wave_r)
        end if
        !
        wave_r = wave_r + wave_r_step
      end do scan_spheres
      call TimerStop('Expansion')
      !
      call TimerStop('Partial Waves')
      call TimerReport
    end subroutine run_partial

  end module partial_waves_sp
!
  subroutine driver
    use partial_waves_sp

    call run_partial
  end subroutine driver
