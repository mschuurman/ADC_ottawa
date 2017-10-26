!
!  * 2013 Oct 25 - Added support for integrating outer-sphere operators
!  * 2013 Apr 26 - Added support for ring singularity of a square-root type.
!
!  Simple non-uniform spatial grid, used for evaluating integrals of bound
!  orbitals and densities numerically. 
!
!  The implementation essentially follows:
!
!   A.D. Becke, J Chem Phys 88, 2547 (1988)
!
!  More sophisticated schemes do exist, especially for larger molecules.
!  However, the original Becke's scheme is simple to implement, efficient
!  enough for small systems, and is accurate enough.
!
!  The routines in this module are intended for "disposable" single-use
!  grids. If the same grid is used multiple times, consider adding 
!  caching of the integration points and weights - perhaps using an
!  external file.
!
!  It is safe and useful to invoke GridPointsBatch('Next') from a parallel
!  region. Each call will receive a unique set of points.
!
  module molecular_grid
    use accuracy
    use atoms
    use timer
    use lebedev
    use math
    !$ use OMP_LIB
    implicit none
    private
    public mol_grid
    public GridInitialize, GridDestroy, GridPointsBatch
    !
    integer(ik), parameter     :: step_order = 3      ! There does not seem to be a good reason to change this one
    logical, parameter         :: randomize  = .true. ! Randomize orientation of the spherical shell; should give more
                                                      ! accurate, but non-deterministic results
    !
    !  Radial grid type is used to cache integration points and weights on a
    !  [-1:1] grid. They'll be scaled as appropriate for each atom.
    !
    type radial_grid
      real(rk), allocatable :: rw(:,:)  ! Grid positions (1,:) and weights (2,:)
    end type radial_grid
    !
    type atom_grid
      real(rk)    :: xyz(3)   ! Coordinates of the centres
      real(rk)    :: rat      ! Parameter used for scaling the radial grid - eq. 25
      integer(ik) :: nrad     ! Number of radial grid points per centre
      integer(ik) :: nang     ! Number of points per angular shell
    end type atom_grid
    !
    type ring_grid
      real(rk)    :: rc(3)    ! ring centre
      real(rk)    :: rn(3)    ! Normal; |ring_rn| = 1
      real(rk)    :: rr       ! Radius of the ring
      integer(ik) :: nphi     ! Number of elliptical slices (phi)
      integer(ik) :: nrad     ! Number of points along "radial" coordinates (mu, nu)
    end type ring_grid
    !
    type outer_grid
      real(rk)    :: rc(3)    ! Coordinates of the centre of the outer grid
      real(rk)    :: rout     ! Characteristic size of the outer grid
      integer(ik) :: nrad     ! Number of radial grid points per centre
      integer(ik) :: nang     ! Number of points per angular shell
    end type outer_grid
    !
    !  We may maintain multiple grids; as the result, all grid parameters must be collected in mol_grid type.
    !
    type mol_grid
      private
      !
      !  Partition management; this part uses generic "partitions", which have to
      !  be combined according to the Becke's blending formula to give the complete grid.
      !
      integer(ik)                    :: npartitions   ! Number of partitions, what else?
      character(len=10), allocatable :: p_type(:)     ! Partition types; the values currently can be:
                                                      ! 'Atom'  - spherical partition, intended for atomic integration
                                                      ! 'Ring'  - a toroidal partition, for integrating complex-coulomb operator
                                                      ! 'Outer' - an outer-shell partition, for integrating CAPs
      integer(ik), allocatable       :: p_shells(:)   ! Number of shells per partition; the actual definition of a shell
                                                      ! is partition-specific
      integer(ik), allocatable       :: p_shellsz(:)  ! Shell size within each partition (so we assume that partitions are made 
                                                      ! of equally-sized shells)
      integer(ik), allocatable       :: p_index(:)    ! Sequential number of each partition within it's partition-specific table
      !
      !  Position of the next point batch within the grid. We always return an entire shell
      !
      integer(ik)                    :: next_part     ! Next partition to be processed; 1 to npartitions
      integer(ik)                    :: next_shell    ! Next shell witin the partion; 1 to p_shells(next_part)
      !$ integer(OMP_LOCK_KIND)      :: update_lock   ! Lock must be acquired before modifying nexr_part, next_shell
      !
      !  Spherical partitions, typically centered on atoms
      !
      integer(ik)                    :: natoms        ! Number of sphere-like centres
      type (atom_grid), allocatable  :: atoms(:)      ! "atom" partitions - these are as described by Becke
      !
      !  Oblate spheroidal grids, needed for the ring singularity
      !
      integer(ik)                    :: nrings        ! Number of toroidal grids; currently must be 0 or 1
      type (ring_grid), allocatable  :: rings(:)      ! "ring" partitions - needed for complex Coulomb continuation
      !
      !  Outer grids, needed for consistent integration in the outer region
      !
      integer(ik)                    :: nouter        ! Number of outer grids; must be 0 or 1. Can't be combined with a toroidal grid
      type (outer_grid), allocatable :: outer(:)      ! "outer" partitions - needed for CAPs integration
      !
      !  Items below are used for caching quantities, which may be a little expensive to evaluate otherwise
      !
      type(radial_grid), allocatable :: radial(:)     ! Radial grid cache. Grid of order N is kept
                                                      ! at index N of radial(). Some (most) elements
                                                      ! of the cache will be empty.
    end type mol_grid
    !
    contains
    !
    !  Externally visible interfaces
    !
    subroutine GridInitialize(grid,nrad,nang,xyz,types, &
                              ring,ring_rc,ring_rn,ring_nrad,ring_nphi, &
                              outer,outer_rc,outer_rout,outer_nrad,outer_nang)
      type(mol_grid), intent(out)       :: grid        ! Grid to initialize
      integer(ik), intent(in)           :: nrad        ! Basic number of radial points; actual number of points 
                                                       ! may depend on this value and atom types
      integer(ik), intent(in)           :: nang        ! Number of angular points; can be 110, 302, or 770,
                                                       ! giving the corresponding Lebedev's grid
      real(rk), intent(in)              :: xyz(:,:)    ! Coordinates of centres, in Bohr
      character(len=*), intent(in)      :: types(:)    ! Types of centres
      ! Special-purpose grid parameters - ring singularity
      logical, intent(in), optional     :: ring        ! Allow for the ring singularity; absense of ring is treated at ring=.false.
      real(rk), intent(in), optional    :: ring_rc(3)  ! Centre of the ring
      real(rk), intent(in), optional    :: ring_rn(3)  ! Normal direction to the ring plane; |rn| is the ring radius
      integer(ik), intent(in), optional :: ring_nrad   ! Number of radial grid points in ring grid
      integer(ik), intent(in), optional :: ring_nphi   ! Number of angular grid points in ring grid
      ! Special-purpose grid parameters - outer sphere
      logical, intent(in), optional     :: outer       ! Add outer-sphere grid; absense is treated as outer=.false.
      real(rk), intent(in), optional    :: outer_rc(3) ! Coordinates of the grid centre
      real(rk), intent(in), optional    :: outer_rout  ! Characteristic size
      integer(ik), intent(in), optional :: outer_nrad  ! Number of radial points
      integer(ik), intent(in), optional :: outer_nang  ! Number of angular points
      !
      integer(ik) :: alloc
      integer(ik) :: natom, iat, jat
      integer(ik) :: irad, rad_max
      real(rk)    :: rat, rij
      integer(ik) :: nring, npart, ipart
      integer(ik) :: nouter
      !
      call TimerStart('MolGrid Initialize')
      !
      if (size(xyz,dim=1)/=3 .or. size(xyz,dim=2)/=size(types)) then
        stop 'molecular_grid%GridInitialize - bad input array sizes'
      end if
      !
      !  Spherical partitions first
      !
      natom = size(xyz,dim=2)
      grid%natoms   = natom
      allocate (grid%atoms(natom),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating array for an ',i8,'-atom grid.')") alloc, natom
        stop 'molecular_grid%GridInitialize - out of memory (1)'
      end if
      !
      !  Check for collisions: this will cause problems with integrals!
      !
      collision_test: do iat=2,natom
        do jat=1,iat-1
          rij = sqrt(sum( (xyz(:,iat)-xyz(:,jat))**2 ))
          if (rij>=spacing(1e5_rk)) cycle
          write (out,"(/'molecular_grid%GridInitialize - collision of atomic centres is not allowed')")
          write (out,"(' Atom ',i5,' at ',3g25.16)") iat, xyz(:,iat)
          write (out,"(' Atom ',i5,' at ',3g25.16)") jat, xyz(:,jat)
          stop 'molecular_grid%GridInitialize - atom collision'
        end do
      end do collision_test
      !
      !  Fill atom tables
      !
      scan_atom_types: do iat=1,natom
        grid%atoms(iat)%xyz  = xyz(:,iat)
        grid%atoms(iat)%nang = nang        ! Checking will occur later
        grid%atoms(iat)%nrad = nrad        ! We'll modify later as appropriate
        grid%atoms(iat)%rat  = 1._rk       ! We'll modify below as appropriate
        rat = AtomCovalentR(types(iat))    ! Returns radius in Angstrom; we use Bohr here
        if (rat<=0) cycle scan_atom_types
        grid%atoms(iat)%rat  = 0.5_rk * rat / abohr
      end do scan_atom_types
      !
      !  Now the toroidal grids, provided that they are present
      !
      nring = 0
      if (present(ring)) then
        if (ring) nring = 1
      end if
      grid%nrings = nring
      allocate (grid%rings(nring),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating array for an ',i8,'-ring grid.')") alloc, nring
        stop 'molecular_grid%GridInitialize - out of memory (2)'
      end if
      if ( nring>0 ) then
        if (.not.present(ring_rc) .or. .not.present(ring_rn) .or. &
            .not.present(ring_nphi) .or. .not.present(ring_nrad)) then
          stop 'molecular_grid%GridInitialize - missing ring parameters'
        end if
        grid%rings(1)%rc   = ring_rc
        grid%rings(1)%rr   = sqrt(sum(ring_rn**2))
        if (grid%rings(1)%rr<=0._rk) stop 'molecular_grid%GridInitialize - requested zero-radius ring?!'
        grid%rings(1)%rn   = ring_rn / grid%rings(1)%rr
        grid%rings(1)%nphi =  ring_nphi
        grid%rings(1)%nrad =  ring_nrad
        ! write (out,"('Added ring at rc = ',3f14.7)") grid%rings(1)%rc 
        ! write (out,"('          Normal = ',3f14.7)") grid%rings(1)%rn 
        ! write (out,"('          Radius = ', f14.7)") grid%rings(1)%rr
      end if
      !
      !  Finally, the optional outer grid
      !
      nouter = 0
      if (present(outer)) then
        if (outer) nouter = 1
      end if
      grid%nouter = nouter
      allocate (grid%outer(nouter),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating array for an ',i8,'-sphere outer grid.')") alloc, nouter
        stop 'molecular_grid%GridInitialize - out of memory (2a)'
      end if
      if ( nouter>0 ) then
        if (.not.present(outer_rc)   .or. .not.present(outer_rout) .or. &
            .not.present(outer_nrad) .or. .not.present(outer_nang)) then
          stop 'molecular_grid%GridInitialize - missing outer grid parameters'
        end if
        grid%outer(1)%rc   = outer_rc
        grid%outer(1)%rout = outer_rout
        if (grid%outer(1)%rout<=0._rk) stop 'molecular_grid%GridInitialize - requested zero-size outer sphere?!'
        grid%outer(1)%nrad = outer_nrad
        grid%outer(1)%nang = outer_nang
      end if
      !
      !  Either toroidal or outer grids can be present, but not both
      !
      if (grid%nrings+grid%nouter>1) stop 'molecular_grid%GridInitialize - can''t combine rings and outer-sphere partitions'
      !
      !  Fill radial grid cache - recalculating radial grids on each call may
      !  be quite expensive!
      !
      rad_max = max(maxval(grid%atoms(:)%nrad),maxval(grid%rings(:)%nrad),maxval(grid%outer(:)%nrad))
      allocate (grid%radial(rad_max),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating radial grid cache. rad_max = ',i8)") alloc, rad_max
        stop 'molecular_grid%GridInitialize - out of memory (3)'
      end if
      !
      fill_radial_cache_atoms: do iat=1,natom
        irad = grid%atoms(iat)%nrad
        call add_to_cache(irad)
      end do fill_radial_cache_atoms
      !
      fill_radial_cache_rings: do iat=1,nring
        irad = grid%rings(iat)%nrad
        call add_to_cache(irad)
      end do fill_radial_cache_rings
      !
      fill_radial_cache_outer: do iat=1,nouter
        irad = grid%outer(iat)%nrad
        call add_to_cache(irad)
      end do fill_radial_cache_outer
      !
      !  Overall partition control - runs over all types of grid present
      !
      npart = grid%natoms + grid%nrings + grid%nouter
      grid%npartitions = npart
      allocate (grid%p_type(npart),grid%p_shells(npart),grid%p_shellsz(npart),grid%p_index(npart),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating array for an ',i8,'-element partition list.')") alloc, npart
        stop 'molecular_grid%GridInitialize - out of memory (6)'
      end if
      !
      ipart = 0
      fill_atoms: do iat=1,grid%natoms
        ipart = ipart + 1
        grid%p_type   (ipart) = 'Atom'
        grid%p_shells (ipart) = grid%atoms(iat)%nrad
        grid%p_shellsz(ipart) = grid%atoms(iat)%nang
        grid%p_index  (ipart) = iat
      end do fill_atoms
      !
      fill_rings: do iat=1,grid%nrings
        ipart = ipart + 1
        grid%p_type   (ipart) = 'Ring'
        grid%p_shells (ipart) = grid%rings(iat)%nphi
        grid%p_shellsz(ipart) = grid%rings(iat)%nrad**2
        grid%p_index  (ipart) = iat
      end do fill_rings
      !
      fill_outer: do iat=1,grid%nouter
        ipart = ipart + 1
        grid%p_type   (ipart) = 'Outer'
        grid%p_shells (ipart) = grid%outer(iat)%nrad
        grid%p_shellsz(ipart) = grid%outer(iat)%nang
        grid%p_index  (ipart) = iat
      end do fill_outer
      if (ipart/=grid%npartitions) stop 'molecular_grid%GridInitialize - count error'
      !
      grid%next_part  = 1
      grid%next_shell = 1
      !$ call omp_init_lock(grid%update_lock)
      call TimerStop('MolGrid Initialize')
      !
      contains
      subroutine add_to_cache(irad)
        integer(ik), intent(in) :: irad
        !
        if (allocated(grid%radial(irad)%rw)) return
        allocate (grid%radial(irad)%rw(2,irad),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i8,' allocating radial grid ',i8)") alloc, irad
          stop 'molecular_grid%GridInitialize - out of memory (4)'
        end if
        call MathGetQuadrature('Legendre',irad,grid%radial(irad)%rw(1,:),grid%radial(irad)%rw(2,:))
      end subroutine add_to_cache
    end subroutine GridInitialize
    !
    subroutine GridDestroy(grid)
      type(mol_grid), intent(inout) :: grid 
      !
      integer(ik) :: iord
      !
      call TimerStart('MolGrid Destroy')
      free_radial_cache: do iord=1,size(grid%radial)
        if (.not.allocated(grid%radial(iord)%rw)) cycle free_radial_cache
        deallocate (grid%radial(iord)%rw)
      end do free_radial_cache
      !
      deallocate (grid%p_type,grid%p_shells,grid%p_shellsz,grid%p_index)
      if (allocated(grid%atoms)) deallocate (grid%atoms)
      if (allocated(grid%rings)) deallocate (grid%rings)
      if (allocated(grid%outer)) deallocate (grid%outer)
      !$ call omp_destroy_lock(grid%update_lock)
      call TimerStop('MolGrid Destroy')
    end subroutine GridDestroy
    !
    subroutine GridPointsBatch(grid,action,xyzw,count,done)
      type(mol_grid), intent(inout)                  :: grid       ! Grid
      character(len=*), intent(in)                   :: action     ! What to do
!     real(rk), pointer, intent(inout), optional     :: xyzw(:,:)  ! Next batch of points
!    Declaration above causes Intel compiler 9.1 to barf at the invocation point of fill_spherical_grid_shell. Compiler bug?
      real(rk), pointer,                optional     :: xyzw(:,:)  ! Next batch of points
      integer(ik), intent(out), optional             :: count      ! Total number of batches
      logical, intent(out), optional                 :: done       ! Fill be set to .true. if out of points
      !
      integer(ik) :: ipart, ishell ! Local copies of next_part and next_shell; needed to avoid holding
                                   ! a lock for too long.
      integer(ik) :: iatom        
      integer(ik) :: shellsz       ! Size of the points batch
      integer(ik) :: alloc       
      !
      call TimerStart('MolGrid PointsBatch')
      batch_case: select case (action)
        case default
          write (out,"('molecular_grid%GridPointsBatch: unknown action ',a)") trim(action)
          stop 'molecular_grid%GridPointsBatch - Bad action'
        case ('Batches count')
          if (.not.present(count)) stop 'molecular_grid%GridPointsBatch - required argument is missing (1)'
          count = sum(grid%p_shells)
        case ('Reset')
          !$ call omp_set_lock(grid%update_lock)
          grid%next_part  = 1
          grid%next_shell = 1
          !$omp flush
          !$ call omp_unset_lock(grid%update_lock)
        case ('Next batch')
          if (.not.present(xyzw)) stop 'molecular_grid%GridPointsBatch - required argument is missing (2)'
          !$ call omp_set_lock(grid%update_lock)
          !$omp flush
          ipart  = grid%next_part
          ishell = grid%next_shell
          if (ipart>grid%npartitions) then
            if (present(done)) then
              done = .true.
            else
              stop 'molecular_grid%GridPointsBatch - run out of points, but "done" is not present!'
            end if
            !$ call omp_unset_lock(grid%update_lock)
            !
            call TimerStop('MolGrid PointsBatch')
            return
          end if
          grid%next_shell = ishell + 1
          if (grid%next_shell > grid%p_shells(ipart)) then
            grid%next_shell = 1
            grid%next_part  = ipart + 1
          end if
          !$omp flush
          !$ call omp_unset_lock(grid%update_lock)
          !
          !  We are commited to a particular batch of points; the rest of this routine will not
          !  modify anything within the grid structure.
          !
          shellsz = grid%p_shellsz(ipart)
          if (associated(xyzw)) then
            if (size(xyzw,dim=2)/=shellsz) deallocate (xyzw)
          end if
          if (.not.associated(xyzw)) then
            allocate (xyzw(4,shellsz),stat=alloc)
            if (alloc/=0) then
              write (out,"('molecular_grid%GridPointsBatch: Error ',i8,' allocating batch buffer for ',i8,' points')") &
                     alloc, shellsz
              stop 'molecular_grid%GridPointsBatch - out of memory'
            end if
          end if
          !
          !  Done with all preliminaries, fill the grid points
          !
          iatom = grid%p_index(ipart)
          partition_type: select case (grid%p_type(ipart))
            case default
              stop 'molecular_grid%GridPointsBatch - unknown partition type '
            case ('Atom')
              call fill_spherical_grid_shell(grid,ipart,iatom,ishell,xyzw(:,:))
            case ('Ring')
              call fill_toroidal_grid_shell (grid,ipart,iatom,ishell,xyzw(:,:))
            case ('Outer')
              call fill_outer_grid_shell    (grid,ipart,iatom,ishell,xyzw(:,:))
          end select partition_type
      end select batch_case
      call TimerStop('MolGrid PointsBatch')
    end subroutine GridPointsBatch
    !
    !  Internal routines beyond this point
    !
    subroutine fill_spherical_grid_shell(grid,ipart,iatom,irad,xyzw)
      type(mol_grid), intent(in) :: grid      ! Everything we want to know about the grid
      integer(ik), intent(in)    :: ipart     ! Partition index
      integer(ik), intent(in)    :: iatom     ! Atom index; this is the same as grid%p_index(ipart)
      integer(ik), intent(in)    :: irad      ! Radial index; ditto
      real(rk), intent(out)      :: xyzw(:,:) ! Grid positions and weights
      !
      !  Step one: fill grid positions and unscaled weights.
      !
      integer(ik)       :: nrad          ! Number of radial points on this atom
      integer(ik)       :: nang          ! Number of angular points on this atom
      integer(ik)       :: iang
      real(rk)          :: rc(3)         ! Center of the atomic sphere
      real(rk)          :: rad           ! Radius of the current sphere
      real(rk)          :: w_rad         ! Integration weight associated with the radial point
      real(rk)          :: r_mid         ! Characteristic radius of the atom
      real(rk), pointer :: ang_xyzw(:,:) ! Our angular grid
      real(rk)          :: p_wgt         ! Becke partitioning weight
      real(rk)          :: euler(3)      ! Random Euler angles
      real(rk)          :: rm(3,3)       ! Randomized rotation matrix
      !
      nrad  = grid%atoms(iatom)%nrad
      nang  = grid%atoms(iatom)%nang
      rc    = grid%atoms(iatom)%xyz(:)
      r_mid = grid%atoms(iatom)%rat
      rad   = grid%radial(nrad)%rw(1,irad) ! Unscaled grid position, [-1:+1] range
      w_rad = grid%radial(nrad)%rw(2,irad)
      !
      !  Scale grid position and weight. We'll use Becke's scaling function:
      !    r_mid * (1+x)/(1-x)
      !  The weight gets scaled by the derivative of the transformation function,
      !  which is:
      !    r_mid * 2/(1-x)**2
      !
      w_rad = w_rad * r_mid * 2._rk / (1-rad)**2
      rad   = r_mid * (1+rad)/(1-rad)
      !
      !  Scale weight by the spherical volume element, since our angular grid 
      !  is normalized to unity
      !
      w_rad = w_rad * 4._rk * pi * rad**2
      !
      !  Choose angular grid
      !
      select case (nang)
        case default
          write (out,"('molecular_grid%fill_spherical_grid_shell: Angular grid ',i6,' is not recognized.')") nang
          stop 'molecular_grid%fill_spherical_grid_shell - bad angular grid'
        case (110) ; ang_xyzw => lebedev_gr17
        case (302) ; ang_xyzw => lebedev_gr29
        case (770) ; ang_xyzw => lebedev_gr47
      end select
      !
      !  Randomize shell orientation, to prevent accumulation or errors
      !
      if (randomize) then
        call random_number(euler)
        euler = twopi * euler
        call MathRotationMatrix(euler,rm)
      else
        call MathSetUnitMatrix(rm)
      end if
      !
      !  Figure out coordinates and weights of atomic grid points in this shell
      !
      xyzw(1:3,:) = spread(rc,dim=2,ncopies=nang) + rad * matmul(rm,ang_xyzw(1:3,:))
      xyzw(  4,:) = w_rad * ang_xyzw(4,:)
      !
      !  Magic sauce: Becke's partitioning weights
      !
      partitioning_weights: do iang=1,nang
        p_wgt        = partitioning_weight(grid,ipart,xyzw(1:3,iang))
        xyzw(4,iang) = p_wgt * xyzw(4,iang)
      end do partitioning_weights
      !*????
      ! write (out,"(('0',4(1x,g20.12)))") xyzw
    end subroutine fill_spherical_grid_shell
    !
    !  Oblate spheroidal grid
    !
    subroutine fill_toroidal_grid_shell(grid,ipart,iring,iphi,xyzw)
      type(mol_grid), intent(in) :: grid      ! Everything we want to know about the grid
      integer(ik), intent(in)    :: ipart     ! Partition index
      integer(ik), intent(in)    :: iring     ! Ring index; this is the same as grid%p_index(ipart)
      integer(ik), intent(in)    :: iphi      ! Angular index [along phi]
      real(rk), intent(out)      :: xyzw(:,:) ! Grid positions and weights
      !
      !  Step one: fill grid positions and unscaled weights.
      !
      real(rk)              :: rc(3)       ! Center of the ring sphere
      real(rk)              :: rn(3)       ! Normal of the ring
      real(rk)              :: rr          ! Radius of the ring
      integer(ik)           :: nphi        ! Number of elliptical slices (phi) - uniform grid
      integer(ik)           :: nrad        ! Number of points along "radial" coordinates (mu, nu)
      real(rk)              :: p_wgt       ! Becke partitioning weight
      real(rk), allocatable :: g_mu(:,:)   ! Grid positions and weights; mu grid [0:inf)
      integer(ik)           :: alloc
      integer(ik)           :: imu, inu    ! Grid sub-indices
      integer(ik)           :: igrid       ! Composite index
      real(rk)              :: mu, nu, phi ! Oblate spheroidal coordinates
      real(rk)              :: rm(3,3)     ! Rotation matrix for turning the toroid to the desired orientation
      !
      rc   = grid%rings(iring)%rc
      rn   = grid%rings(iring)%rn
      rr   = grid%rings(iring)%rr
      nphi = grid%rings(iring)%nphi
      nrad = grid%rings(iring)%nrad
      !
      !  The nu and phi grids are equallt spaced; the mu grid, however, requires a bit of
      !  more careful handling.
      !
      allocate (g_mu(2,nrad),stat=alloc)
      if (alloc/=0) stop 'molecular_grid%fill_toroidal_grid_shell - out of memory'
      g_mu = grid%radial(nrad)%rw
      !
      !  Scale nu grid position and weight. Since the distance from the origin
      !  grows exponentially with mu, we'll use:
      !
      !    Log[1+Pi*(1+x)/(1-x)]
      !
      !  The weight gets scaled by the derivative of the transformation function,
      !  which is:
      !
      !    2*Pi/((x-1)(Pi+1+(Pi-1)*x))
      !
      g_mu(2,:) = g_mu(2,:)*twopi/((1-g_mu(1,:))*(pi+1+(pi-1)*(g_mu(1,:))))
      g_mu(1,:) = log(1+pi*(1+g_mu(1,:))/(1-g_mu(1,:)))
      !
      !  Construct the product grid of (mu,nu), in the standard orientation
      !  and at the origin.
      !
      phi   = pi * real(2*iphi-nphi-1,kind=rk) / nphi
      igrid = 0
      fill_grid_nu: do inu=1,nrad
        nu = 0.5_rk * pi * real(2*inu-nrad-1,kind=rk) / nrad
        fill_grid_mu: do imu=1,nrad
          mu = g_mu(1,imu)
          igrid = igrid + 1
          xyzw(1,igrid) = rr * cosh(mu) * cos(nu) * cos(phi)
          xyzw(2,igrid) = rr * cosh(mu) * cos(nu) * sin(phi)
          xyzw(3,igrid) = rr * sinh(mu) * sin(nu)
          xyzw(4,igrid) = rr**3 * cosh(mu) * cos(nu) * (sinh(mu)**2 + sin(nu)**2) &
                        * g_mu(2,imu) * (pi/nrad) * (twopi/nphi)
        end do fill_grid_mu
      end do fill_grid_nu
      if (igrid/=size(xyzw,dim=2)) stop 'molecular_grid%fill_toroidal_grid_shell - count error'
      !
      !  Translate and rotate the shell to the desired position
      !
      call MathRotationMatrix((/0._rk,-acos(rn(3)),-atan2(rn(2),rn(1))/),rm)
      ! write (out,*) ' rn = ', rn
      ! write (out,*) ' xx = ', matmul(rm,(/0._rk,0._rk,1._rk/))
      move_grid: do igrid=1,size(xyzw,dim=2)
        xyzw(1:3,igrid) = rc + matmul(rm,xyzw(1:3,igrid))
      end do move_grid
      !
      !  Magic sauce: Becke's partitioning weights
      !
      partitioning_weights: do igrid=1,size(xyzw,dim=2)
        p_wgt         = partitioning_weight(grid,ipart,xyzw(1:3,igrid))
        xyzw(4,igrid) = p_wgt * xyzw(4,igrid)
      end do partitioning_weights
      !* ????
      ! write (out,"(('1',4(1x,g20.12)))") xyzw
    end subroutine fill_toroidal_grid_shell
    !
    subroutine fill_outer_grid_shell(grid,ipart,isphere,irad,xyzw)
      type(mol_grid), intent(in) :: grid      ! Everything we want to know about the grid
      integer(ik), intent(in)    :: ipart     ! Partition index
      integer(ik), intent(in)    :: isphere   ! Outer sphere index; this is the same as grid%p_index(ipart)
      integer(ik), intent(in)    :: irad      ! Radial index; ditto
      real(rk), intent(out)      :: xyzw(:,:) ! Grid positions and weights
      !
      !  Step one: fill grid positions and unscaled weights.
      !
      integer(ik)       :: nrad          ! Number of radial points 
      integer(ik)       :: nang          ! Number of angular points 
      integer(ik)       :: iang
      real(rk)          :: rc(3)         ! Center of the sphere
      real(rk)          :: rad           ! Radius of the current sphere
      real(rk)          :: w_rad         ! Integration weight associated with the radial point
      real(rk)          :: r_mid         ! Characteristic radius of the grid
      real(rk), pointer :: ang_xyzw(:,:) ! Our angular grid
      real(rk)          :: p_wgt         ! Becke partitioning weight
      real(rk)          :: euler(3)      ! Random Euler angles
      real(rk)          :: rm(3,3)       ! Randomized rotation matrix
      !
      nrad  = grid%outer(isphere)%nrad
      nang  = grid%outer(isphere)%nang
      rc    = grid%outer(isphere)%rc(:)
      r_mid = grid%outer(isphere)%rout
      rad   = grid%radial(nrad)%rw(1,irad) ! Unscaled grid position, [-1:+1] range
      w_rad = grid%radial(nrad)%rw(2,irad)
      !
      !  Scale grid position and weight. We'll use Becke's scaling function:
      !    r_mid * (1+x)/(1-x)
      !  The weight gets scaled by the derivative of the transformation function,
      !  which is:
      !    r_mid * 2/(1-x)**2
      !
      w_rad = w_rad * r_mid * 2._rk / (1-rad)**2
      rad   = r_mid * (1+rad)/(1-rad)
      !
      !  Scale weight by the spherical volume element, since our angular grid 
      !  is normalized to unity
      !
      w_rad = w_rad * 4._rk * pi * rad**2
      !
      !  Choose angular grid
      !
      select case (nang)
        case default
          write (out,"('molecular_grid%fill_outer_grid_shell: Angular grid ',i6,' is not recognized.')") nang
          stop 'molecular_grid%fill_outer_grid_shell - bad angular grid'
        case (110) ; ang_xyzw => lebedev_gr17
        case (302) ; ang_xyzw => lebedev_gr29
        case (770) ; ang_xyzw => lebedev_gr47
      end select
      !
      !  Randomize shell orientation, to prevent accumulation or errors
      !
      if (randomize) then
        call random_number(euler)
        euler = twopi * euler
        call MathRotationMatrix(euler,rm)
      else
        call MathSetUnitMatrix(rm)
      end if
      !
      !  Figure out coordinates and weights of atomic grid points in this shell
      !
      xyzw(1:3,:) = spread(rc,dim=2,ncopies=nang) + rad * matmul(rm,ang_xyzw(1:3,:))
      xyzw(  4,:) = w_rad * ang_xyzw(4,:)
      !
      !  Magic sauce: Becke's partitioning weights
      !
      partitioning_weights: do iang=1,nang
        p_wgt        = partitioning_weight(grid,ipart,xyzw(1:3,iang))
        xyzw(4,iang) = p_wgt * xyzw(4,iang)
      end do partitioning_weights
      !*????
      ! write (out,"(('0',4(1x,g20.12)))") xyzw
    end subroutine fill_outer_grid_shell
    !
    !  Calculate Becke's partitioning weight for a grid point
    !
    function partitioning_weight(grid,ipart,xyz) result(w)
      type(mol_grid), intent(in) :: grid
      integer(ik), intent(in)    :: ipart  ! Parent partition
      real(rk), intent(in)       :: xyz(:) ! Grid point
      real(rk)                   :: w      ! Partitioning weight
      !
      real(rk)    :: w_tot
      integer(ik) :: jpart
      !
      w     = raw_weight(grid,ipart,xyz)
      w_tot = w
      all_claims: do jpart=1,grid%npartitions
        if (jpart==ipart) cycle all_claims
        w_tot = w_tot + raw_weight(grid,jpart,xyz)
      end do all_claims
      w = w / w_tot
    end function partitioning_weight
    !
    !  Calculate Becke's primitive cell weight
    !
    function raw_weight(grid,ipart,xyz) result(w)
      type(mol_grid), intent(in) :: grid
      integer(ik), intent(in)    :: ipart  ! Parent partition
      real(rk), intent(in)       :: xyz(:) ! Grid point
      real(rk)                   :: w      ! Partitioning weight
      !
      integer(ik) :: jpart
      integer(ik) :: iatom, jatom
      real(rk)    :: ri, rj   ! Distance from grid point to partitions I and J
      real(rk)    :: rij      ! Distance between partitions I and J
      real(rk)    :: muij     ! Elliptical distance difference coordinate
      real(rk)    :: nuij     ! "Adjusted" cell boundary
      real(rk)    :: aij      ! Adjustment parameter
      real(rk)    :: chi      ! Ratio of atom sizes
      real(rk)    :: uij      ! 
      real(rk)    :: sf
      !
      iatom = grid%p_index(ipart)
      ri = distance_point_to_partition(grid,ipart,xyz)
      w = 1.0_rk
      cell_accumulate: do jpart=1,grid%npartitions
        if (jpart==ipart) cycle cell_accumulate
        jatom = grid%p_index(jpart)
        rj    = distance_point_to_partition(grid,jpart,xyz)
        rij   = distance_partition_to_partition(grid,ipart,jpart,xyz)
        muij  = (ri-rj)/rij
        !
        !  "Adjust" the cell boundary to account for heteronuclear bonding - Becke's appendix A
        !  The adjustment is only applicable if both cells are atoms!
        !
        if (grid%p_type(ipart)=='Atom' .and. grid%p_type(jpart)=='Atom') then
          chi  = grid%atoms(iatom)%rat/grid%atoms(jatom)%rat
          uij  = (chi-1._rk)/(chi+1._rk)
          aij  = uij/(uij**2-1._rk)
          aij  = min(0.5_rk,max(-0.5_rk,aij))
        else
          aij  = 0._rk
        end if
        !
        nuij = muij + aij * (1._rk-muij**2)
      ! sf   = step_function(muij)           !???? This was my original code; it looks like a bug!
        sf   = step_function(nuij)
        w    = w * sf
      end do cell_accumulate
    end function raw_weight
    !
    function distance_point_to_partition(grid,ipart,xyz) result(ri)
      type(mol_grid), intent(in) :: grid 
      integer(ik), intent(in)    :: ipart
      real(rk), intent(in)       :: xyz(:)  ! Point for which we need the distance
      real(rk)                   :: ri
      !
      integer(ik) :: iatom
      !
      iatom = grid%p_index(ipart)
      select case (grid%p_type(ipart))
        case default
          stop 'molecular_grid%distance_point_to_partition - unknown partion type'
        case ('Atom')  
          ri = sqrt(sum( (xyz-grid%atoms(iatom)%xyz)**2 ))
        case ('Ring')
          ri = pt_ring(xyz,grid%rings(iatom))
        case ('Outer')
          ri = abs(sqrt(sum((xyz-grid%outer(iatom)%rc)**2))-grid%outer(iatom)%rout)
       end select
    end function distance_point_to_partition
    !
    function distance_partition_to_partition(grid,ipart,jpart,xyz) result(rij)
      type(mol_grid), intent(in) :: grid 
      integer(ik), intent(in)    :: ipart  ! Partitions we need the distance between
      integer(ik), intent(in)    :: jpart
      real(rk), intent(in)       :: xyz(:) ! Point we are looking at; it may affect the partition-partition distance
      real(rk)                   :: rij
      !
      integer(ik) :: iatom, jatom
      !
      iatom = grid%p_index(ipart)
      jatom = grid%p_index(jpart)
      select case (trim(grid%p_type(ipart))//' '//trim(grid%p_type(jpart)))
        case default
          stop 'molecular_grid%distance_partition_to_partition - unknown partion type combination'
        case ('Atom Atom')  
          rij = sqrt(sum( (grid%atoms(jatom)%xyz-grid%atoms(iatom)%xyz)**2 ))
        case ('Atom Ring')
          rij = sphere_ring(xyz,grid%rings(jatom),grid%atoms(iatom)%xyz)
        case ('Ring Atom')
          rij = sphere_ring(xyz,grid%rings(iatom),grid%atoms(jatom)%xyz)
        case ('Atom Outer')
          rij = sphere_outer(xyz,grid%outer(jatom),grid%atoms(iatom)%xyz)
        case ('Outer Atom')
          rij = sphere_outer(xyz,grid%outer(iatom),grid%atoms(jatom)%xyz)
      end select
    end function distance_partition_to_partition
    !
    !  Distance from a point to the nearest point on a ring
    !
    function pt_ring(xyz,ring) result (ri)
      real(rk), intent(in)        :: xyz(:)
      type(ring_grid), intent(in) :: ring
      real(rk)                    :: ri
      !
      real(rk) :: r12(3)   ! Vector from the centre of the ring to the point
      real(rk) :: rn (3)   ! Projection of r12 to ring normal
      real(rk) :: rp (3)   ! r12 - rn
      real(rk) :: lrn, lrp ! Lengths of rn and rp
      !
      r12 = xyz - ring%rc
      rn  = ring%rn * sum(ring%rn*r12)
      rp  = r12 - rn
      lrn = sqrt(sum(rn**2))
      lrp = sqrt(sum(rp**2))
      ri  = sqrt(lrn**2 + (lrp-ring%rr)**2)
    end function pt_ring
    !
    !  Distance from a centre to the point on a ring nearest to a given point
    !
    function sphere_ring(xyz,ring,sph_xyz) result (ri)
      real(rk), intent(in)        :: xyz(:)      ! Grid point
      type(ring_grid), intent(in) :: ring        ! Ring
      real(rk), intent(in)        :: sph_xyz(:)  ! Sphere centre
      real(rk)                    :: ri
      !
      real(rk) :: tmp_xyz(3) ! Fallback for direction choice
      real(rk) :: rim(3)     ! Point on the ring; set by get_rim_direction()
      !
      if (.not.get_rim_direction(xyz)) then
        !
        !  Grid point is on-axis; try to choose direction towards the sphere centre
        !
        if (.not.get_rim_direction(sph_xyz)) then
          !
          !  Sphere centre is on-axis as well; choose direction towards smallest 
          !  component of ring%rn
          !
          tmp_xyz = 0
          tmp_xyz(minloc(ring%rn,dim=1)) = 1._rk
          if (.not.get_rim_direction(tmp_xyz)) then
            stop 'molecular_grid%sphere_ring - an impossible stop reached'
          end if
        end if
      end if
      !
      ri  = sqrt(sum((sph_xyz-rim)**2))
      !
      contains 
      function get_rim_direction(pt_xyz) result(happy)
        real(rk), intent(in) :: pt_xyz(:) ! Point which defines the required rim direction
        real(rk)             :: r12(3)    ! Vector from the centre of the ring to the point
        real(rk)             :: rn (3)    ! Projection of r12 to ring normal
        real(rk)             :: rp (3)    ! r12 - rn
        real(rk)             :: lrp       ! Lengths of rn and rp
        logical              :: happy     ! Set to .true. if rim direction found
        !
        r12 = pt_xyz - ring%rc
        rn  = ring%rn * sum(ring%rn*r12)
        rp  = r12 - rn
        lrp = sqrt(sum(rp**2))
        if (lrp>=spacing(1e3_rk*ring%rr)) then
          rim   = ring%rc + rp * (ring%rr/lrp)
          happy = .true.
        else
          happy = .false.
        end if
      end function get_rim_direction
    end function sphere_ring
    !
    !  Distance from a centre to the point on a sphere nearest to a given point
    !
    function sphere_outer(xyz,outer,sph_xyz) result (ri)
      real(rk), intent(in)         :: xyz(:)      ! Grid point
      type(outer_grid), intent(in) :: outer       ! Outer sphere
      real(rk), intent(in)         :: sph_xyz(:)  ! Atomic sphere centre
      real(rk)                     :: ri
      !
      real(rk) :: rel_xyz(3) ! Position of the grid point relative to the outer sphere
      real(rk) :: rrel       ! Length of rel_xyz(:)
      real(rk) :: rim(3)     ! Point on the sphere
      !
      !  try to find direction from centre of the outer sphere towards the grid point
      !
      rel_xyz = xyz - outer%rc
      rrel    = sqrt(sum(rel_xyz**2))
      if (rrel>spacing(1e4_rk*outer%rout)) then
        rel_xyz = rel_xyz / rrel
      else
        !
        !  grid point is very close to the center of the outer sphere; chose direction towards the atomic sphere
        !
        rel_xyz = sph_xyz - outer%rc
        rrel    = sqrt(sum(rel_xyz**2))
        if (rrel>spacing(1e4_rk*outer%rout)) then
          rel_xyz = rel_xyz / rrel
        else
          !
          !  atomic sphere is also very close to the outer sphere; then it does not matter what we choose
          !
          rel_xyz = (/ 0, 0, 1 /)
        end if
      end if
      rim = outer%rc + outer%rout * rel_xyz
      ri  = sqrt(sum((sph_xyz-rim)**2))
    end function sphere_outer
    !
    function step_function(mu) result(sf)
      real(rk), intent(in)    :: mu
      real(rk)                :: sf
      !
      integer(ik) :: iord
      !
      sf = mu
      recurse_function: do iord=1,step_order
        sf = 1.5_rk * sf - 0.5_rk * sf**3
      end do recurse_function
      sf = 0.5_rk * (1.0_rk - sf)
    end function step_function
  end module molecular_grid
