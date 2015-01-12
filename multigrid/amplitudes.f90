!
!  Multigrid test - calculation of plane-wave and eikonal recombination
!                   amplitudes from Dyson orbitals and exchage corrections
!                   to them
!
  module amplitudes
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use eikonal_tools
    use lebedev
    use opendx
    implicit none
    private
    public run_amplitudes
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_mos        = 4                   ! Max number of MOs we support: Dyson + 3x exchange
    integer(ik), parameter :: max_fields     = 3+max_mos           ! 4 = 3x RDF components + planewave + phase correction
    integer(ik), parameter :: max_naturals   = 200                 ! Max number of natural orbitals allowed
    integer(ik), parameter :: unit_dump      = 34                  ! A more or less random unit
    integer(ik), parameter :: unit_fit       = 35                  ! Unit for writing angular fit data
    integer(ik), parameter :: maxAmplitudes  = 4                   ! Max number of possible choices of continuum wavefunctions,
                                                                   ! relevant for averaging
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_points(3)     = (/ 200, 200, 200 /)   ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    character(len=100) :: dyson_file      = 'dyson.dat'           ! Name of the file containing MO coefficients
                                                                  ! for the Dyson and exchange orbitals
    character(len=100) :: natural_file    = 'natural.dat'         ! Name of the file containing MO coefficients
                                                                  ! for the natural orbitals of the cation
                                                                  ! (used for calculating Hartree potential)
    character(len=100) :: potential_file  = ' '                   ! File for the potential dump - can be HUGE
    integer(ik)        :: spin_parts(4)   = (/ 2, 4, 6, 8 /)      ! Orbital components to load from "dyson_file"
    integer(ik)        :: natural_count   = -1                    ! Number of natural orbitals in the density matrix
                                                                  ! Value of -1 (which is the default) requests that natural
                                                                  ! occupations are loaded from the same file as the
                                                                  ! natural orbitals themselves.
    real(rk)           :: natural_occ(max_naturals)               ! Natural orbital occupation numbers
    real(rk)           :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                  ! and exchange potential(s)
    real(rk)           :: sor_rate        = 1.95_rk               ! Over-relaxation rate in solving Poisson equations
    logical            :: plot_phases     = .false.               ! Plot numerical phases
    logical            :: eikonal_pref    = .true.                ! Include eikonal pre-factor
    logical            :: node_factors    = .false.               ! Include integrals necessary for calculating emission
                                                                  ! from orbitals with nodes
    character(len=100) :: v_xc            = ' '                   ! Exchange-correlation potential to use in construction of the
                                                                  ! eikonal functions. Possible choices are:
                                                                  ! ' ' or 'none' - no exchange potential
                                                                  ! 'Slater' - Dirac-Slater exchange
                                                                  ! 'SVWN'   - Slater exchange + VWN local correlation
    character(len=100) :: momenta         = 'input'               ! Choice of asymptotic momenta. Can be:
                                                                  ! 'input'  - Energy and direction are specified explicitly
                                                                  ! 'grid'   - Use a regular grid. In the latter case, we'll
                                                                  !            also perform incoherent angular averaging
    character(len=100) :: fit_file        = 'recombination.dat'   ! Data file to write spherical harmonics fit of recombination
                                                                  ! amplitudes. Has no effect unless momenta == 'grid'.
    character(len=20)  :: fit_component   = 'RO'                  ! Amplitude component to fit; can be 'P', 'PO', 'R', or 'RO'
    integer(ik)        :: fit_index       = 0                     ! Derived from fit_component as appropriate
    integer(ik)        :: recombination_l_max = 8                 ! Max angular momentum in matrix element fit
    character(len=100) :: operator        = 'dipole'              ! One-electron operator for the matrix elements. Can be:
                                                                  !  'dipole'       - Cartesian dipole operator. "Exchange" 
                                                                  !                   corrections will be included, as long 
                                                                  !                   as the non-orthogonal Dyson term is 
                                                                  !                   present.
                                                                  !  'acceleration' - Gradient of the potential. External 
                                                                  !                   potential (natural_file) must be 
                                                                  !                   supplied for this form to work.
    character(len=100) :: angular_grid    = 'lg17'                ! Choice of angular grid. Can be one of:
                                                                  ! 'lg9'    - Order-9 Lebedev grid
                                                                  ! 'lg17'   - Order-17 Lebedev grid
                                                                  ! 'lg29'   - Order-29 Lebedev grid
                                                                  ! 'lg47'   - Order-47 Lebedev grid
    real(rk)           :: ekin_min        =  0.5_rk               ! Kinetic energy of the electron, in eV
    real(rk)           :: ekin_max        = 20.0_rk               ! 
    real(rk)           :: ekin_step       =  0.1_rk               ! 
    integer(ik)        :: grid_theta      = 200_ik                ! Spherical product grid, used to resample spherical fit for opendx vis.
    integer(ik)        :: grid_phi        = 400_ik
    character(len=100) :: odx_template    = "('rec_',f0.3,'.dx')" ! Name template for recombination matrix element magnitude
    character(len=100) :: table_template  = "('rec_',f0.3,'.table')" ! Name template for interpolated recombination matrix elements
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_dyson                                 ! Dyson orbital
    integer(ik)        :: f_exchange(3)                           ! Exchange orbitals
    integer(ik)        :: f_rdf(3)                                ! Recombination dipole/acceleration field
    integer(ik)        :: f_core_pot = -1                         ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_plane                                 ! Probing planewave
    integer(ik)        :: f_pref                                  ! Prefactor of the eikonal function
    !
    real(rk)           :: dipole_dyson(3)                         ! Dipole moment (e-Bohr) of the Dyson orbital (operator="dipole")
                                                                  ! Acceleration of the Dyson orbital (operator="acceleration")
    real(rk), pointer  :: angularGrid(:,:) => NULL()              ! Angular integration grid for observation direction
    !
    !  Data for calculating spherical harmonics expansions
    !
    complex(rk), allocatable :: angularData (:,:)                 ! Complex vector amplitudes at Lebedev grid positions
    complex(rk), allocatable :: angularCoeff(:,:,:)               ! Spherical harmonics fit coefficients
    !
    !  Visualization grid
    !
    real(rk), allocatable    :: prod_grid(:,:,:)                  ! Spherical product grid of a unit radius at zero
                                                                  ! Indices are: (X:Y:Z:W,ITH,IPH)
    complex(rk), allocatable :: prod_ylm   (:,:)                  ! Spherical harmonics on the product grid
    complex(rk), allocatable :: prod_data  (:,:,:)                ! Interpolated data on the product grid
    real(rk)                 :: step_phi
    real(rk)                 :: step_theta
    !
    real(rk), parameter      :: r_zero(3) = 0._rk                 ! This must be zero
    complex(rk), parameter   :: v_zero    = 0._rk                 ! This must be zero
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /ampdata/ verbose, &
                       n_points, box_extent, &
                       dyson_file, spin_parts, &
                       eps_hartree, sor_rate, &
                       natural_count, natural_file, natural_occ, &
                       potential_file, plot_phases, eikonal_pref, &
                       node_factors, v_xc, &
                       operator, &
                       momenta, angular_grid, ekin_min, ekin_max, ekin_step, &
                       recombination_l_max, fit_file, fit_component, &
                       grid_theta, grid_phi, odx_template, table_template
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Prepare simulation box and allocate a sufficient number of data fields
    !
    subroutine initialize_grid
      real(rk)    :: box(2,3)
      integer(ik) :: field
      !
      call TimerStart('Grid initialization')
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
      call MultiGridInit(max_grids=1,max_fields=max_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate_fields: do field=1,max_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      f_free = max_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Allocate fields for and load Dyson and exchange correction orbitals
    !
    subroutine load_gamess_mos
      real(rk)              :: norm
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      f_free = f_free - 4
      if (f_free<0) stop 'load_gamess_mos - out of fields'
      call FieldImport('GAMESS',dyson_file,f_table(f_free+1:f_free+4),spin_parts)
      f_dyson         = f_table(f_free+1)
      f_exchange(1:3) = f_table(f_free+2:f_free+4)
      write (out,"('Loaded orbital components ',4i3)") spin_parts
      norm = FieldNorm(f_dyson)**2
      write (out,"('<psid|psid> = ',f20.10)") norm
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
    end subroutine load_gamess_mos
    !
    !  Build RDF (Recombination dipole field) from the Dyson and exchange orbitals
    !
    subroutine build_rdf
      integer(ik) :: ic, scr, scr_rdf
      real(rk)    :: ef(3)
      !
      scr_rdf = f_table(f_free) ; f_free = f_free - 1
      scr     = f_table(f_free) ; f_free = f_free - 1
      if (f_free<0) stop 'build_rdf - not enough fields!'
      !
      build_component: do ic=1,3
        ef = 0 ; ef(ic) = 1._rk
        call FLsetField(ef)
        call FieldCopy(src=f_exchange(ic),dst=scr_rdf)
        call FieldInit(dst=scr,func=FLelectricField)
        call FieldMulAdd(src_a=scr,src_b=f_dyson,dst=scr_rdf)
        call FieldCopy(src=scr_rdf,dst=f_exchange(ic)) ! Exchange component is no longer needed
      end do build_component
      !
      !  Copy handles from f_exchange to f_rdf, and kill values in f_exchange
      !  to prevent reuse
      !
      f_rdf      = f_exchange
      f_exchange = -1
      f_free     = f_free + 2
      !
      !  Multipole moments for the Dyson orbital - we need the dipole
      !
      call FieldNormMultipoles(f_dyson,dipole_dyson)
    end subroutine build_rdf
    !
    !  Build RAF (Recombination acceleration field) from the Dyson and exchange orbitals
    !
    subroutine build_raf
      integer(ik) :: ic
      !
      if (f_core_pot<=0) then
        write (out,"('Effective potential is not available; can''t take gradients!')") 
        stop 'amplitudes%build_raf - Missing the potential'
      end if
      !
      !  "Exchange" terms are not defined for this matrix elements - wipe out and reuse
      !  the "exchange" fields. Also calculate acceleration of the dyson orbital. It should
      !  be zero for bound states, but one never knows what we are going to get...
      !
      build_component: do ic=1,3
        call FieldGradientComponent(dir=ic,src=f_core_pot,dst=f_exchange(ic))
        call FieldMul(dst=f_exchange(ic),src=f_dyson)
        dipole_dyson(ic) = real(FieldConjgIntegrate(left=f_exchange(ic),right=f_dyson),kind=rk)
      end do build_component
      !
      !  Copy handles from f_exchange to f_rdf, and kill values in f_exchange
      !  to prevent reuse
      !
      f_rdf      = f_exchange
      f_exchange = -1
    end subroutine build_raf
    !
    !  Very simple visualization routine
    !
    subroutine visualize_wavefunction(text,psi)
      character(len=*), intent(in) :: text
      integer(ik), intent(in)      :: psi
      !
      call FieldVisualize(slot=0,src=psi,text=trim(text))
    end subroutine visualize_wavefunction
    !
    !  Prepare grid for electron momenta (if requected)
    !
    subroutine initialize_momenta
      select case (momenta)
        case default
          write (out,"('Asymtotic momentum choice ""',a,'"" is not recognized.')") trim(momenta)
          stop 'amplitudes - Bad choice for asymptotic momenta'
        case ('input')
          return
        case ('grid')
      end select
      !
      !  We only get here if momenta == 'grid'
      !
      select case (angular_grid)
        case default
          write (out,"('Angular grid ""',a,'"" is not recognized')") trim(angular_grid)
          stop 'amplitudes - bad angular grid'
        case ('lg9','LG9')
          angularGrid => lebedev_gr9
        case ('lg17','LG17')
          angularGrid => lebedev_gr17
        case ('lg29','LG29')
          angularGrid => lebedev_gr29
        case ('lg47','LG47')
          angularGrid => lebedev_gr47
      end select
      !
      if (verbose>=-1) then
        write (out,"('Using ',i4,'-point observation direction grid')") size(angularGrid,dim=2)
      end if
    end subroutine initialize_momenta
    !
    subroutine dump_potential
      call TimerStart('Dump core potential')
      open (unit=unit_dump,form='formatted',status='replace',file=potential_file)
      call FieldDump(unit_dump,src=f_core_pot)
      close (unit=unit_dump)
      call TimerStop('Dump core potential')
    end subroutine dump_potential
    !
    subroutine print_legend_rdf
      write (out,"(/'  P   = plane wave amplitude, with exchange corrections included')")
      write (out,"( '  PO  = plane wave, with exchange and orthogonalization corrections')")
      write (out,"( '  R   = eikonal amplitude (length gauge), incl. exchange corrections')")
      write (out,"( '  RO  = eikonal (length gauge), incl. exchange and orho. corrections'/)")
      write (out,"( '  PNX = nodal corrections to the plane-wave amplitudes (X), ditto PNY/PNZ')")
      write (out,"( '  PXO = plane wave nodal correction, orthogonalized, ditto PYO/PZO')")
      write (out,"( '  RNX = nodal corrections eikonal amplitudes, ditto RNY/RNZ')")
      write (out,"( '  RXO = eikonal nodal corrections, orthogonalized, ditto RYO/RZO')")
      write (out,"(/'  dx         = <Psi^d|x|k> + <Psi^x|k>, ditto for dy, dz')")
      write (out,"( '  dx-ortho   = dx - <Psi^d|x|Psi^d><Psi^d|k>, ditto for dy-o, dz-o')")
      write (out,"( '  dx-nodal-Z = <Psi^d|x|z k> + <Psi^x|z k>, ditto for dx-X, dx-Y etc.')")
      write (out,"(/'  where Psi^d is the Dyson orbital')")
      write (out,"( '        Psi^x, Psi^y, and Psi^z are exchange correction orbitals')")
      write (out,"( '        |k> are either plane wave or eikonal solutions for a given momentum')")
      write (out,"()")
    end subroutine print_legend_rdf
    !
    subroutine print_legend_raf
      write (out,"(/'  P   = plane wave amplitude')")
      write (out,"( '  PO  = plane wave, with orthogonalization correction')")
      write (out,"( '  R   = eikonal amplitude (length gauge)')")
      write (out,"( '  RO  = eikonal (length gauge), with orhogonalization corrections'/)")
      write (out,"(/'  dx         = <Psi^d|dV/dx|k>, ditto for dy, dz')")
      write (out,"( '  dx-ortho   = dx - <Psi^d|dV/dx|Psi^d><Psi^d|k>, ditto for dy-o, dz-o')")
      write (out,"(/'  where Psi^d is the Dyson orbital')")
      write (out,"( '        |k> are either plane wave or eikonal solutions for a given momentum')")
      write (out,"()")
    end subroutine print_legend_raf
    !
    !  Calculate dipole amplitudes for a continuum wavefunction in
    !
    subroutine generate_amplitudes(tag,label,wf,amptable)
      character(len=1), intent(in) :: tag    ! Distinctive tag to use in the printout
      character(len=*), intent(in) :: label  ! Descriptive label for the continnum wf parameters
      integer(ik), intent(in)      :: wf     ! Field containing the continuum wavefunction
      complex(rk), intent(out)     :: amptable(3,2) ! Results for the caller
      !
      integer(ik) :: ic
      complex(rk) :: amp(3), r_amp(3,3)             ! Amplitudes and nodal corrections
      complex(rk) :: amp_overlap, r_amp_overlap(3)  ! Overlap with the Dyson orbital
      complex(rk) :: amp_ort(3), r_amp_ort(3,3)     ! Orthogonalized amplitudes
      !
      character(len=1), parameter :: xyz(3) = (/ 'X', 'Y', 'Z' /)
      !
      !  Exact scattering solutions would have zero overlap, but ...
      !
      amp_overlap = FieldConjgIntegrate(left=f_dyson,right=wf)
      amplitude_components: do ic=1,3
        amp(ic)     = FieldConjgIntegrate(left=f_rdf(ic),right=wf)
        amp_ort(ic) = amp(ic) - dipole_dyson(ic)*amp_overlap
      end do amplitude_components
      amptable(:,1) = amp
      amptable(:,2) = amp_ort
      !
      if (verbose>=0) then
        write (out,"(1x,a1,'  ',1x,a44,1x,3(1x,g12.6,1x,g12.6,1x),2(f10.5,1x))") &
               tag, label, amp, amp_overlap
        write (out,"(1x,a1,'O ',1x,a44,1x,3(1x,g12.6,1x,g12.6,1x))") &
               tag, label, amp_ort
      end if
      !
      !  Now the nodal corrections
      !
      if (node_factors) then
        call FieldNorm2Multipoles(left=f_dyson,right=wf,mult=r_amp_overlap)
        amplitude_components_nodes: do ic=1,3
          call FieldNorm2Multipoles(left=f_rdf(ic),right=wf,mult=r_amp(:,ic))
          r_amp_ort(:,ic) = r_amp(:,ic) - r_amp_overlap(:)*dipole_dyson(ic)
        end do amplitude_components_nodes
        !
        !  Printing is a separate loop: we need to transpose the result.
        !
        if (verbose>=0) then
          amplitude_nodes_print: do ic=1,3
            write (out,"(1x,a1,'N',a1,1x,a44,1x,3(1x,g12.6,1x,g12.6,1x),2(f10.5,1x))") &
                   tag, xyz(ic), label, r_amp(ic,:), r_amp_overlap(ic)
            write (out,"(1x,a1,a1,'O',1x,a44,1x,3(1x,g12.6,1x,g12.6,1x))") &
                   tag, xyz(ic), label, r_amp_ort(ic,:)
          end do amplitude_nodes_print
        end if
      end if
    end subroutine generate_amplitudes
    !
    !  Evaluate all requested combinations of the matrix elements for
    !  a single choice of electron energy ands observation direction
    !
    subroutine generate_amplitudes_for_a_single_momentum(eval,edir,amptable)
      real(rk), intent(in)     :: eval, edir(3)
      complex(rk), intent(out) :: amptable(3,maxAmplitudes) 
                                  ! The last index of amptable corresponds to:
                                  !   1 = plane wave
                                  !   2 = plane wave, orthogonalized
                                  !   3 = eikonal
                                  !   4 = eikonal, orthogonalized
      real(rk)                 :: kvec(3), lambda
      character(len=200)       :: comment, label
      !
      !  Normalize the direction, and convert harmonic number to the k vector
      !
      kvec   = edir
      kvec   = kvec / sqrt(sum(kvec**2))
      kvec   = kvec * sqrt(2*eval/h2ev)
      lambda = twopi / sqrt(sum(kvec**2))
      !
      !  Prepate descriptive label for later use
      !
      write (label,"(f9.4,1x,3(1x,f7.4),2x,f8.3)") eval, edir, lambda
      amptable = 0._rk
      !
      call FLsetPlanewave(k=kvec)
      !
      !  Pure plane wave for reference
      !
      call FieldInit(dst=f_plane,func=FLplanewave)
      call generate_amplitudes(tag='P',label=label,wf=f_plane,amptable=amptable(:,1:2))
      !
      !  Eikonal correction factor
      !
      if (natural_count>0) then
        !
        !  Phase-corrected plane wave. The first call will construct correction
        !  phase factor, which is then incorporated into the plane wave
        !
        call eikonal_build_function(dst=f_plane,kvec=kvec,potential=f_core_pot, &
                                    prefactor=eikonal_pref,scratch=f_pref,norm='natural')
        !
        if (plot_phases) then
          write (comment,"('Eikonal wavefunction for k= ',3(1x,f15.7))") kvec
          call visualize_wavefunction(trim(comment),f_plane)
        end if
        !
        call generate_amplitudes(tag='R',label=label,wf=f_plane,amptable=amptable(:,3:4))
      end if
    end subroutine generate_amplitudes_for_a_single_momentum
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
      complex(rk), intent(in)    :: wgt(:)     ! Weight of the source
      complex(rk), intent(in)    :: src(:,:)   ! Vector to add
      complex(rk), intent(inout) :: dst(:,:,:) ! Destination vector
      !
      integer(ik) itheta, iphi
      !
      call TimerStart('Add two kets')
      !$omp parallel do private (iphi,itheta)
      phi_loop: do iphi=1,grid_phi
        theta_loop: do itheta=1,grid_theta
          dst(itheta,iphi,:) = dst(itheta,iphi,:) + wgt*src(itheta,iphi)
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
      loop_l: do lv=0,recombination_l_max
        loop_m: do mv=-lv,lv
          call fill_harmonic(lv,mv,prod_ylm)
          call add_kets(AngularCoeff(mv,lv,:),src=prod_ylm,dst=prod_data)
        end do loop_m
      end do loop_l
      !
    end subroutine interpolate
    !
    !  Report an ASCII table of the interpolated matrix elements (complex 3-vector)
    !
    subroutine print_recombination_table(filename)
      character(len=*), intent(in) :: filename
      !
      integer(ik) :: itheta, iphi
      real(rk)    :: phi, theta
      !
      call TimerStart('Print recombination table')
      open (unit=unit_dump,form='formatted',status='replace',file=trim(filename))
      !
      write (unit_dump,"('#',a12,1x,a12,3(2x,a16,1x,a16))") &
             ' phi,Rad ', ' theta,Rad ', ' Re(<dys|x|k>) ', ' Im(<dys|x|k>) ', ' ... '
      phi_points: do iphi=1,grid_phi  ! [0:2Pi]
        phi = step_phi*(iphi-0.5_rk)
        theta_points: do itheta=1,grid_theta
          theta = step_theta*(itheta-0.5_rk)
          write (unit_dump,"(1x,f12.8,1x,f12.8,3(2x,e16.8,1x,e16.8))") phi, theta, prod_data(itheta,iphi,:)
        end do theta_points
      end do phi_points
      !
      close (unit=unit_dump)
      call TimerStop('Print recombination table')
    end subroutine print_recombination_table
    !
    !  Spherical harmonics fitting
    !
    subroutine fit_initialize
      integer(ik)        :: recombination_n_ekin
      character(len=256) :: recombination_comment
      !
      namelist /recombination_fit/ recombination_l_max, recombination_n_ekin, recombination_comment
      !
      if (fit_file==' ') return
      open (unit_fit,form='formatted',status='replace',file=trim(fit_file))
      !
      recombination_n_ekin  = int((ekin_max-ekin_min+(ekin_step-spacing(ekin_step)))/ekin_step,kind=ik)
      recombination_comment = trim(fit_component)//' '//trim(dyson_file)//' '//trim(dyson_file)//' '//trim(potential_file)
      write (unit_fit,nml=recombination_fit)
      !
      select case (fit_component)
        case default
          write (out,"('fit_initialize - fit component designation ',a,' is not recognized')") trim(fit_component)
          stop 'fit_initialize - bad component name'
        case ('P', 'p' ) ; fit_index = 1
        case ('PO','po') ; fit_index = 2
        case ('R' ,'r' ) ; fit_index = 3
        case ('RO','ro') ; fit_index = 4
      end select
      !
      allocate (angularData(3,size(angularGrid,dim=2)))
      allocate (angularCoeff(-recombination_l_max:recombination_l_max,0:recombination_l_max,3))
      allocate (prod_grid(4,grid_theta,grid_phi))
      allocate (prod_ylm   (grid_theta,grid_phi))
      allocate (prod_data  (grid_theta,grid_phi,3))
      !
      call initialize_angular_grid
    end subroutine fit_initialize
    !
    subroutine fit_start(eval)
      real(rk), intent(in) :: eval
      !
      if (fit_file==' ') return
      angularData(:,:) = 0
      write (unit_fit,"(g25.14)") eval
    end subroutine fit_start
    !
    subroutine fit_accumulate(idir,amptable)
      integer(ik), intent(in) :: idir          ! Position within angularGrid
      complex(rk), intent(in) :: amptable(:,:) ! Amplitudes
      !
      if (fit_file==' ') return
      angularData(:,idir) = amptable(:,fit_index)
    end subroutine fit_accumulate
    !
    !  The fit_report routine is modified version of harmonics_project from spherical_interpolate.f90
    !
    subroutine fit_report(eval)
      real(rk), intent(in)   :: eval
      !
      integer(ik)        :: lv, mv, ipt, ic
      complex(rk)        :: hv, acc(3)
      character(len=256) :: visual_file
      !
      if (fit_file==' ') return
      !
      !  Perform projection, and fill the linear system matrix for the fit
      !
      angularCoeff = 0
      loop_l1: do lv=0,recombination_l_max
        loop_m1: do mv=-lv,lv
          call FLharmonicsSetParameters(lv,mv,r_zero,0._rk,safe_max)
          acc = 0
          loop_pt: do ipt=1,size(angularGrid,dim=2)
            hv  = FLharmonics(angularGrid(1:3,ipt),v_zero)
            acc = acc + angularGrid(4,ipt) * conjg(hv) * angularData(:,ipt)
          end do loop_pt
          angularCoeff(mv,lv,:) = acc
        end do loop_m1
      end do loop_l1
      !
      !  Write out the fit
      !
      fit_component: do ic=1,3
        fit_l2: do lv=0,recombination_l_max
          fit_m2: do mv=-lv,lv
            write (unit_fit,"(1x,i3,1x,i4,2x,g24.15,1x,g24.15)") lv, mv, angularCoeff(mv,lv,ic)
          end do fit_m2
        end do fit_l2
      end do fit_component
      !
      !  OpenDX visualization part
      !
      call interpolate
      if (table_template/=" ") then
        write (visual_file,table_template) eval
        call print_recombination_table(visual_file)
      end if
      write (visual_file,odx_template) eval
      call odx_write_spherical(visual_file,sqrt(real(sum(prod_data*conjg(prod_data),dim=3),kind=rk)))
    end subroutine fit_report
    !
    subroutine fit_finalize
      if (fit_file==' ') return
      close (unit_fit)
      deallocate (angularData,angularCoeff)
    end subroutine fit_finalize
    !
    !  Problem driver
    !
    subroutine run_amplitudes
      integer(ik)        :: info, idir
      real(rk)           :: eval, edir(3)
      complex(rk)        :: amptable   (3,maxAmplitudes)  ! Amplitudes for each direction
      real(rk)           :: probability(3,maxAmplitudes)  ! Angular average of the probability
      !
      call TimerStart('Amplitudes')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=ampdata,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=ampdata)
      write (out,"()")
      !
      call initialize_momenta
      !
      call initialize_grid
      !
      !  Do we have natural orbitals to construct the density and the Hartree potential?
      !
      if (natural_count<0) then
        call gamess_load_natocc(natural_file,natural_occ,natural_count)
      end if
      if (natural_count>0) then
        if (natural_count>max_naturals) then
          write (out,"('Too many natural orbitals. Increase ''max_naturals''" &
                     // " to at least ',i0,' and recompile')") natural_count
          stop 'amplitudes - too many natural orbitals'
        end if
        call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
        if (f_free<1) stop 'amplitudes - not enough fields for core potential'
        f_core_pot = f_table(f_free) ; f_free = f_free - 1
        call eikonal_build_potential(natural_file,natural_occ(:natural_count),f_table(:f_free),f_core_pot,v_xc=v_xc)
        if (potential_file/=' ') then
          call dump_potential ! This call is likely to produce a HUGE file
        end if
      end if
      !
      !  Orbitals
      !
      call load_gamess_mos
      !
      !  Note that build_rdf/build_raf will reuse fields originally allocated for f_exchange
      !
      select case (operator)
        case default
          write (out,"('Operator ""',a,'"" is not recognized. Oops.')") trim(operator)
          stop 'amplitudes - bad operator (1)'
        case ('dipole')
          write (out,"(/'Using dipole form of the matrix elements'/)")
          call build_rdf
          write (out,"('Dipole moment of the Dyson orbital [electron-Bohr]:',3(1x,g12.6))") dipole_dyson
        case ('acceleration')
          write (out,"(/'Using acceleration form of the matrix elements'/)")
          if (node_factors) then
            write (out,"('Nodal corrections are not available for the acceleration matrix element'/)")
            node_factors = .false.
          end if
          call build_raf
          write (out,"('Acceleration of the Dyson orbital [Bohr/au[t]^2]:',3(1x,g12.6))") dipole_dyson
      end select
      !
      write (out,"(/'Done with preliminaries'/)")
      !
      !  The rest of the input is a free-format sequence of:
      !
      !  e(eV), ex, ey, ez
      !
      !  where e is the harmonic energy in eV; ex, ey, and ez are the
      !  planewave propagation direction.
      !
      if (f_free<2) stop 'amplitudes - not enough scratch fields'
      f_plane = f_table(f_free) ; f_free = f_free - 1
      f_pref  = f_table(f_free) ; f_free = f_free - 1
      !
      select case (operator)
        case default
          write (out,"('Operator ""',a,'"" is not recognized. Oops.')") trim(operator)
          stop 'amplitudes - bad operator (1)'
        case ('dipole')
          call print_legend_rdf
        case ('acceleration')
          call print_legend_raf
      end select
      !
      if (verbose>=0) then
        write (out,"(4x,a9,1x,3(1x,a7),2x,a8,1x,3(1x,a12,1x,a12,1x),2x,a10,1x,a10)") &
           'E,eV', ' X ', ' Y ', ' Z ', ' L,Bohr ', ' Re(dx) ', ' Im(dx) ', &
                   ' Re(dy) ', ' Im(dy) ', ' Re(dz) ', ' Im(dz) ', &
                   ' Re(<d|k>) ', ' Im(<d|k>) '
      end if
      if (momenta=='grid') then
        write (out,"('@   ',a9,4(1x,a12))") 'E,eV', 'Probability', ' Prob.X ', ' Prob.Y ', ' Prob.Z '
      end if
      !
      select case (momenta)
        case ('input')
          !
          !  Explicit specification of the energies and observation directions
          !
          amplitudes_input_loop: do
            read(input,*,iostat=info) eval, edir
            if (info/=0) exit amplitudes_input_loop
            call generate_amplitudes_for_a_single_momentum(eval,edir,amptable)
          end do amplitudes_input_loop
        case ('grid')
          !
          !  Regular choice of the directions - also do incoherent angular averaging
          !  and prepare spherical harmonics fit.
          !
          call fit_initialize
          !
          eval = ekin_min
          energy_grid: do while(eval<=ekin_max)
            probability = 0._rk
            call fit_start(eval)
            direction_grid: do idir=1,size(angularGrid,dim=2)
              edir = angularGrid(1:3,idir)
              call generate_amplitudes_for_a_single_momentum(eval,edir,amptable)
              !
              call fit_accumulate(idir,amptable)
              !
              !  We are interested in total ionization probabilities for an aligined
              !  molecule. If there is no prefered orientation of the electric field,
              !  then just take 1/3 of the sum across the first dimension.
              !
              probability = probability + angularGrid(4,idir)*abs(amptable)**2
            end do direction_grid
            !
            !  Report averaged ionization probabilities
            !
            write (out,"('@P  ',f9.4,4(1x,e12.6))") eval, sum(probability(:,1))/3._rk, probability(:,1)
            write (out,"('@PO ',f9.4,4(1x,e12.6))") eval, sum(probability(:,2))/3._rk, probability(:,2)
            if (natural_count>0) then
              write (out,"('@R  ',f9.4,4(1x,e12.6))") eval, sum(probability(:,3))/3._rk, probability(:,3)
              write (out,"('@RO ',f9.4,4(1x,e12.6))") eval, sum(probability(:,4))/3._rk, probability(:,4)
            end if
            !
            call fit_report(eval)
            eval = eval + ekin_step
          end do energy_grid
          !
          call fit_finalize
      end select
      !
      call TimerStop('Amplitudes')
      call TimerReport
    end subroutine run_amplitudes

  end module amplitudes
!
  subroutine driver
    use amplitudes

    call run_amplitudes
  end subroutine driver
