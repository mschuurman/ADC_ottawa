!
!  * April 26, 2013: Added support for integrating a (single) square-root 
!                    ring singularity, which arises for certain imaginary
!                    displacements.
!
!  Evaluation of analytical continuation of the "correlation potential" into
!  the complex plane - see e-mail from Lisa Torina on July 27th, 2012.
!
!  WARNING: This is NOT a superset of correlation_potential_v2, even though
!  WARNING: there is some overlap in functionality.
!
!  This is the the third version of the code, using brute-force numerical 
!  integration for the integral. The complex correlation potential is:
!
!         /
!    v(R)=| dr rho(r) [(r-R).(r-R)]**-0.5
!         /
!
!  where R is allowed to be complex, while r is real. The square root phase
!  is chosen such that real part is always non-negative. (?)
!
!  Unfortunately, this operator is not suitable for the standard recursion 
!  relations, and has to be evaluated numerically.
!
!  We'll use Becke's grid construction procedure. As an accuracy check,
!  we'll use the same grid to evaluate Coulomb potential at Re[R] and
!  the norm of the density - both of these quantities are available
!  analytically.
!
  module correlation_potential_v3
    use accuracy
    use timer
    use import_gamess
    use lebedev
    use math
    use molecular_grid
    !$ use OMP_LIB

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)           :: verbose      = 0          ! Verbosity level
    integer(ik)           :: max_rdm_sv   = 1000       ! Should be good enough ...
    character(len=100)    :: rdm_file     = ' '        ! Name of the file containing rdm orbital coefficients
    integer(ik)           :: angular_npts = 302        ! Angular grid used in numerical integration. The number of
                                                       ! points must match one of the Lebedev grids in lebedev.f90,
                                                       ! currently 110, 302, or 770.
    integer(ik)           :: radial_npts  = 100        ! Number of radial grid points at each centre
    logical               :: ring_singularity = .true. ! Try to integrate the ring singularity
    integer(ik)           :: ring_nphi    =  50        ! Number of angular grid points in spheroidal grid
    integer(ik)           :: ring_nrad    =  50        ! Number of "radial" grid points in spheroidal grid
    real(rk)              :: ring_threshold = 0.05_rk  ! Very small rings are intergated by the atomic spheres perfectly well,
                                                       ! so don't bother with them
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)           :: rdm_count                 ! Number of singular values
    real(rk), allocatable :: rdm_sv(:)                 ! Singular values of the transition 1-RDM
    real(rk), allocatable :: rdm_ao(:,:)               ! Reduced density matrix in the AO basis
    type(gam_structure)   :: mol                       ! Gamess basis set and orbital data
    integer(ik)           :: mol_natoms                ! Number of atoms in the molecule
    real(rk), allocatable :: mol_xyzq(:,:)             ! Coordinates and charges of the nuclei; extra entry at the end
                                                       ! for the operator position
    character(len=20), allocatable :: mol_labels(:)    ! Atom labels
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /correlation/ verbose, rdm_file, max_rdm_sv, radial_npts, angular_npts, &
                           ring_singularity, ring_nphi, ring_nrad, ring_threshold
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Report molecular geometry - lifted from dyson_tools.f90
    !
    subroutine slurp_structure(name)
      character(len=*), intent(in) :: name  ! Name of the structure/file
      !
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call gamess_report_nuclei(nnuc,structure=mol)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      !
      !  Tell a bit more!
      !
      allocate (xyzq(4,nnuc))
      call gamess_report_nuclei(nnuc,xyzq,structure=mol)
      !
      !  Make copy of the data for grid construction
      !
      mol_natoms = nnuc
      allocate (mol_xyzq(4,nnuc+1),mol_labels(nnuc+1))
      mol_xyzq(:,:nnuc)  = xyzq
      mol_labels(:nnuc)  = mol%atoms(:nnuc)%name
      mol_xyzq(:,nnuc+1) = 0.0_rk
      mol_labels(nnuc+1) = "Xx"
      !
      if (verbose>=0) then
        write (out,"()")
        write (out,"(      t12,a36,t52,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
        write (out,"(      t12,a36,t52,a36)") '------------------    ', '----------------------    '
        write (out,"(t6,a5,t12,3a12,t52,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
        print_atoms: do iat=1,nnuc
          write (out,"(t2,a3,t6,f5.2,t12,3f12.5,t52,3f12.5)") mol_labels(iat)(1:3), xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
        end do print_atoms
        write (out,"()")
      end if
      deallocate (xyzq)
    end subroutine slurp_structure
    !
    !  Prepare ion-specific potentials: core+Hartree and transition potentials
    !
    subroutine load_molecular_data
      real(rk)    :: tmp_rdm_sv(max_rdm_sv)
      !
      call TimerStart('Load molecular data')
      write (out,"(/'Loading ',a)") trim(rdm_file)
      call gamess_load_rdmsv(trim(rdm_file),tmp_rdm_sv,rdm_count)
      write (out,"( 'Found ',i4,' singular values')") rdm_count
      write (out,"( 'Values are: '/)")
      write (out,"(10(1x,f12.8))") tmp_rdm_sv(:rdm_count)
      !
      allocate (rdm_sv(rdm_count))
      rdm_sv = tmp_rdm_sv(:rdm_count)
      !
      call gamess_load_orbitals(file=trim(rdm_file),structure=mol)
      !
      call slurp_structure(trim(rdm_file))
      call TimerStop('Load molecular data')
    end subroutine load_molecular_data
    !
    !  Transform 1-RDM to the AO basis
    !
    subroutine transform_1rdm
      integer(ik) :: nao, nmo, ird, imo, mu, nu
      real(rk)    :: tmp
      !
      call TimerStart('1-RDM: MO->AO')
      !
      nao = mol%nbasis
      nmo = mol%nvectors
      allocate (rdm_ao(nao,nao))
      !
      !$omp parallel do default(none) private(nu,mu,imo,tmp) &
      !$omp& shared(nao,nmo,mol,rdm_ao,rdm_count,rdm_sv)
      rdm_ao_nu: do nu=1,nao
        rdm_ao_mu: do mu=1,nao
          tmp = 0._rk
          rdm_ao_sv: do ird=1,rdm_count
            imo = 2*ird-1
            tmp = tmp + rdm_sv(ird) * mol%vectors(mu,imo) * mol%vectors(nu,imo+1)
          end do rdm_ao_sv
          rdm_ao(mu,nu) = tmp
        end do rdm_ao_mu
      end do rdm_ao_nu
      !$omp end parallel do
      !
      call TimerStop('1-RDM: MO->AO')
    end subroutine transform_1rdm
    !
    function evaluate_property(what,coord) result(v)
      character(len=*), intent(in) :: what     ! Operator to evaluate
      real(rk), intent(in)         :: coord(:) ! X,Y,Z of the operator
      real(rk)                     :: v        ! Desired integral
      !
      integer(ik)           :: nao
      real(rk), allocatable :: ao_ints(:,:)
      !
      call TimerStart('Property '//trim(what))
      nao = mol%nbasis
      allocate (ao_ints(nao,nao))
      call gamess_1e_integrals('AO 3C '//what,ao_ints,bra=mol,ket=mol,op_xyz=coord)
      v = sum(ao_ints*rdm_ao)
      deallocate (ao_ints)
      call TimerStop('Property '//trim(what))
    end function evaluate_property
    !
    !  Evaluate transition density on a grid of points. For efficiency reasons, we'll
    !  block the evaluation, and use BLAS3-like transformation.
    !
    subroutine evaluate_transition_density(xyz,rho)
      real(rk), intent(in)  :: xyz(:,:) ! XYZ coordinates of the grid points
      real(rk), intent(out) :: rho(:)   ! Transition density at these grid points
      !
      integer(ik)            :: ipt, ird, imo, npt, nmo, nbas
      real(ark), allocatable :: basval(:,:,:)
      real(rk), allocatable  :: moval (:,:)
      !
      if (size(xyz,dim=1)/=3 .or. size(xyz,dim=2)/=size(rho)) stop 'evaluate_transition_density - bad dimensions'
      !
      npt  = size(rho)
      nmo  = mol%nvectors
      nbas = mol%nbasis
      !
      allocate (basval(1,nbas,npt),moval(nmo,npt))
      !
      !  First, evaluate basis functions
      !
      evaluate_basis_functions: do ipt=1,npt
        call gamess_evaluate_functions(xyz(:,ipt),basval(:,:,ipt),mol)
      end do evaluate_basis_functions
      !
      !  Transform AOs to the MOs, for all grid points simultaneously
      !
      moval = matmul(transpose(mol%vectors(:,:nmo)),basval(1,:,:))
      !
      !  Finally, evaluate the transition density at grid points
      !
      rho(:) = 0
      evaluate_rdm: do ird=1,rdm_count
        imo = 2*ird - 1
        rho = rho + rdm_sv(ird) * moval(imo,:) * moval(imo+1,:)
      end do evaluate_rdm
      !
      deallocate (basval,moval)
    end subroutine evaluate_transition_density
    !
    subroutine grid_integrate(coords,vcomplex,vreal,overlap)
      complex(rk), intent(in)  :: coords(:) ! Coordinates of the complex point
      complex(rk), intent(out) :: vcomplex  ! Coulomb potential at the point
      real(rk), intent(out)    :: vreal     ! Coulomb potential at Re[coords]
      real(rk), intent(out)    :: overlap   ! Overlap of the density
      !
      type(mol_grid)        :: grid
      integer(ik)           :: grid_natoms   ! Number of "atomic" centres in the grid
      integer(ik)           :: nbatch, ib    ! Number of grid point batches
      integer(ik)           :: npts          ! Number of grid points in a batch
      integer(ik)           :: ipt           ! Grid point within the batch
      real(rk)              :: r_real        ! Distance in real space
      complex(rk)           :: rc2           ! Distance in complex space, squared
      real(rk)              :: a, b          ! Real and imaginary parts of rc2
      complex(rk)           :: r_complex     ! Distance in complex space
      real(rk), pointer     :: xyzw(:,:)     ! Coordinates and weights of the grid points
      real(rk), allocatable :: rho(:)        ! Transition density on grid
      logical               :: have_ring     ! True if any of the probe coordinates have imaginary component
      !
      call TimerStart('Grid Integrate')
      !
      !  Prepare integration grid. Don't forget to insert a point at the
      !  operator position - but only if the operator is away from all atoms
      !
      grid_natoms = mol_natoms
      !
      have_ring = any(abs(aimag(coords))>ring_threshold) .and. ring_singularity
      !
      if (.not.have_ring) then
        if (all(sqrt(sum((mol_xyzq(1:3,:mol_natoms)-spread(real(coords,kind=rk),dim=2,ncopies=mol_natoms))**2,dim=1))>1e-6_rk)) then
          mol_xyzq(1:3,mol_natoms+1) = real(coords,kind=rk)
          grid_natoms = grid_natoms + 1
        end if
      end if
      !
      !
      call GridInitialize(grid,radial_npts,angular_npts,mol_xyzq(1:3,:grid_natoms),mol_labels(:grid_natoms), &
                          ring=have_ring, &
                          ring_rc=real(coords,kind=rk),ring_rn=aimag(coords), &
                          ring_nphi=ring_nphi,ring_nrad=ring_nrad)
      call GridPointsBatch(grid,'Batches count',count=nbatch)
      !
      !  Numerical integration loop
      !
      vreal    = 0 
      vcomplex = 0
      overlap  = 0
      !$omp parallel default(none) &
      !$omp& shared(grid,nbatch,coords) &
      !$omp& private(ib,r_real,rc2,a,b,r_complex,xyzw,rho,ipt,npts) &
      !$omp& reduction(+:vreal,overlap,vcomplex)
      nullify(xyzw)
      !$omp do
      grid_batches: do ib=1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
        !
        !  If the batch size changed, reallocate rho()
        !
        npts = size(xyzw,dim=2)
        !
        if (allocated(rho)) then
          if (size(rho)/=npts) deallocate (rho)
        end if
        if (.not.allocated(rho)) allocate (rho(npts))
        !
        !  Evaluate density at grid points
        !
        call evaluate_transition_density(xyzw(1:3,:),rho)
        !
        !  Ready to integrate!
        !
        integrate: do ipt=1,npts
          ! write (out,"('#pt ',4g20.12)") xyzw(:,ipt)
          overlap = overlap + xyzw(4,ipt) * rho(ipt)
          r_real  = sqrt(sum((xyzw(1:3,ipt)-real(coords,kind=rk))**2))
          if (abs(r_real)>=spacing(100._rk)) then
            vreal = vreal + xyzw(4,ipt) * rho(ipt) / r_real
          end if
          !
          !  We need to be a little careful with the phase of complex distance,
          !  so it's not a good idea to use the intrinsic complex sqrt()
          !
          rc2 = sum((xyzw(1:3,ipt)-coords)**2)
          a   = real(rc2,kind=rk)
          b   = aimag(rc2)
          r_complex = cmplx(     sqrt(sqrt(a**2+b**2)+a), &
                            sign(sqrt(sqrt(a**2+b**2)-a),b),kind=rk) / sqrt(2._rk)
          ! ???? The line below should give the equivalent result!
          ! r_complex = sqrt(rc2)
          !
          if (abs(r_complex)<=1e-5) then
            write (out,*) ' xyzw = ',xyzw(:,ipt), ' r_complex = ', r_complex
          end if
          if (abs(r_complex)>=spacing(100._rk)) then
            vcomplex = vcomplex + xyzw(4,ipt) * rho(ipt) / r_complex
          end if
        end do integrate
      end do grid_batches
      !$omp end do
      if (associated(xyzw)) deallocate (xyzw)
      if (allocated (rho) ) deallocate (rho)
      !$omp end parallel
      !
      call GridDestroy(grid)
      call TimerStop('Grid Integrate')
    end subroutine grid_integrate
    !
    subroutine start
      integer(ik) :: info
      complex(rk) :: coords(3)  ! Complex coordinates of the trajectory point
      real(rk)    :: overlap_a  ! Grid accuracy check: Norm of the input density, analytical result
      real(rk)    :: vreal_a    ! Grid accuracy check: Coulomb potential at Re[R], analytical result
      complex(rk) :: vcomplex_n ! Coulomb potential at complex R, numerical
      real(rk)    :: overlap_n  ! Norm of the input density, numerical
      real(rk)    :: vreal_n    ! Coulomb potential at Re[R], numerical
      !
      call TimerStart('start')
      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=correlation,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=correlation)
      write (out,"()")
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call load_molecular_data
      !
      call transform_1rdm
      !
      write (out,"(('#',3(1x,a12,1x,a12,1x),2x,2(a18,1x),4(2x,a18)))") &
          '  Re[X]  ', '  Im[X]  ', '  Re[Y]  ', '  Im[Y]  ', '  Re[Z]  ', '  Im[Z]  ', '  Re[V]  ', '  Im[V]  ', &
                       ' V[Re(R)], num ', ' V[Re(R)], anal ', ' Norm, num ', ' Norm, anal ', &
          '---------', '---------', '---------', '---------', '---------', '---------', '---------', '---------', &
                       '---------------', '----------------', '-----------', '------------'
      run_trajectory: do
        read (input,*,iostat=info) coords
        if (info/=0) exit run_trajectory
        !
        vreal_a   = evaluate_property('1/R',real(coords,kind=rk))
        overlap_a = evaluate_property('ONE',real(coords,kind=rk))
        call grid_integrate(coords,vcomplex_n,vreal_n,overlap_n)
        write (out,"('@',3(1x,f12.6,1x,f12.6,1x),2x,2(g18.11,1x),4(2x,g18.11))") &
               coords, vcomplex_n, vreal_n, vreal_a, overlap_n, overlap_a
      end do run_trajectory
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
    
  end module correlation_potential_v3
  !
  !
  !
  subroutine driver
    use accuracy
    use math
    use correlation_potential_v3
    !
    real(rk) :: dummy
    !
    call accuracyInitialize
    dummy = MathDoubleFactorial(100_ik)
    dummy = MathFactorial(100_ik)
    dummy = MathLogDoubleFactorial(1000_ik)
    dummy = MathLogFactorial(1000_ik)

    call start

  end subroutine driver
