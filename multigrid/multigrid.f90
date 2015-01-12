module multigrid
  use accuracy
  use complex
  use interpolate
  use parallel
  use fftw
  use kinetic
  use timer
  use import_gamess
  use phase_interpolation
  implicit none
  private
  private verbose, chkptIO
  public MultiGridInit, SimpleGridNew, FieldNew
  public FieldInit, FieldCopy, FieldConjugate, FieldScale, FieldSum
  public FieldLaplacian, FieldDump, FieldNorm, FieldNorm1
  public FieldProductIntegrate, FieldConjgIntegrate
  public FieldMul, FieldMulAdd, FieldAXPY, FieldZero, FieldIO, FieldVisualize, FieldShow
  public FieldSanitize, FieldFFT, FieldCheckOut, FieldProcess, FieldRotateVectorComponents
  public FieldNormMultipoles, FieldCheckpoint, FieldSetOuterWall, FieldGetLeakage
  public FieldShrinkBox, FieldScalarIntegrate, FieldBraVKet, FieldSetWavefunction
  public FieldInvertHamiltonian, FieldDiffuseLaplacian
  public FieldGradientComponent, FieldGradientComponentRight, FieldComponentCount
  public FieldGridSpacing, FieldGridNPoints, FieldGridCoordinates, FieldGridExtent
  public FieldExport, FieldImport, FieldIterationSOR, FieldRhoAccumulate, FieldNorm1Multipoles
  public FieldNorm2Multipoles, FieldBuildPhase,FieldBuildPref
  public FieldEVAInteriorStep, FieldEVAExteriorStep, FieldFetchSlab, FieldSetSlab
  public FieldFetchSlice, FieldGER
  public FieldECPProject, FieldECPApply, FieldECPGetExtent
!
  interface FieldFetchSlice
    module procedure FieldFetchSlice2D
    module procedure FieldFetchSlice3D
  end interface FieldFetchSlice
  interface FieldGER
    module procedure FieldGER2D
    module procedure FieldGER3D21
    module procedure FieldGER3D111
  end interface FieldGER
!
  integer(ik), parameter :: cl          = 80                     ! Max character string length
  integer(ik), parameter :: verbose     = 1                      ! Verbosity level
  integer(ik), parameter :: gridIO      = 55                     ! External I/O unit for swapping
  integer(ik), parameter :: chkptIO     = 56                     ! External I/O unit for checkpointing
!
!  zeroField parameter is etremely important - anything smaller than this in
!  absolute magnitude will be considered zero, and not included in the calculations.
!  Setting zeroField to a negative value will force the multigrid to use complete
!  meshes, even when the contain zero sections.
!
  real(rk), parameter    :: zeroField          = -1              ! Anything smaller is zero
  logical, parameter     :: allowPolar         = .false.         ! Enable norm-conserving grid reconciliation
  integer(ik), parameter :: interpolationOrder = 3               ! Interpolation order, used for stitching
  integer(ik), parameter :: extrapolationOrder = 1               ! Extrapolation order, used for flux evaluation
!
!  DX is limited in the size of 3D grids it could handle. In any event, it
!  will never be able to handle (100x100x100) grids, or some such - so we
!  must provide a mechanism for downsampling the grid.
!
! integer(ik), parameter :: maxDXgridSize = 45*45*90             ! Anything bigger must be downsampled
  integer(ik), parameter :: maxDXgridSize = 200*200*200          ! Anything bigger must be downsampled
!
!  Simple, rectangular grid type.
!
  type SimpleGridT
    integer(ik)            :: npoints (  3)       ! Number of grid patches in each dimension.
    real(rk)               :: range   (2,3)       ! Coordinate ranges for a given dimension - min and max
                                                  ! The range corresponds to the extreme outside points
                                                  ! of the integration grid
    real(rk)               :: step    (  3)       ! Grid spacing for a given dimension. This is guaranteed to
                                                  ! be: step(:) = ( range(2,:) - range(1,:) ) / npoints(:)
    real(rk)               :: weight              ! Uniform integration weight for the patches
                                                  ! This is really = product(step(:)), or the volume of each
                                                  ! cell, if you wish.
    integer(ik)            :: up_grid             ! Index of the containing grid (0 = none)
    integer(ik)            :: up_index(2,3)       ! Indices within the parent grid, corresponding to
                                                  ! the range(:,:)
    integer(ik), pointer   :: down_index(:,:,:)   ! Indices within this grid, corresponding to each
                                                  ! parent cell index. The first index is low, high,
                                                  ! the second index is parent cell number, and the
                                                  ! third index is the spatial coordinate
    integer(sik), pointer  :: down_grid  (:,:,:)  ! Indices of the subdividing grids (0 = none)
    real(rk), pointer      :: coords   (:,:,:,:)  ! Coordinates and weights of the grid points
                                                  ! Coordinates (the first three values for the first index)
                                                  ! refer to the -centre- of the corridponding grid volume
    complex(rk), pointer   :: fields   (:,:,:,:)  ! Fields on grid - the last index is the field index
    integer(ik), pointer   :: active   (:,:,:)    ! Range of active indices for each field. Data outside
                                                  ! the active range is (implicitly) zero, and should be
                                                  ! ignored. The last index is field. Field of zero is special,
                                                  ! and covers the whole mesh; Second last index
                                                  ! is the coordinate (1=x, etc); the first index gives
                                                  ! the low (1) and the upper (2) bounds
    !
    !  Alternative distribution parameters, used to describe fields in the
    !  momentum space. Currently, the only function to produce such distributions
    !  is FieldFFT. The only way to use momentum distributions is currently to
    !  plot it (with FieldVisualize) or run analysis on it with FieldCheckOut
    !
    !  Unlike real-space grids, momentum grids are not hierarchical - they are
    !  simply the velocity/momentum distributions, corresponding to a spacial
    !  wavefunction.
    !
    real(rk)               :: p_range (2,3)
    real(rk)               :: p_step    (3)
    real(rk)               :: p_weight
    real(rk), pointer      :: p_coords(:,:,:,:)
  end type SimpleGridT
!
!  Multigrid type definition
!
  type MultiGridT
    integer(ik)                :: ngrids_max      ! Number of simple grids requested at initialization
    integer(ik)                :: nfields_max     ! Number of scalar fields requested at initialization
    integer(ik)                :: nborder         ! Number of border points for each subgrid
    integer(ik)                :: ngrids          ! Number of simple grids active
    integer(ik)                :: nfields         ! Number of scalar fields present
    character(len=cl), pointer :: field_names (:) ! Identifying names for the fields
    character(len=cl), pointer :: grid_names  (:) ! Identifying names for the simple grids
    logical, pointer           :: scratch     (:) ! True if the field does not need to be checkpointed
    logical, pointer           :: wavefunction(:) ! Field contains a wavefunction, and has to use
                                                  ! norm-concerving reconciliation.
    type(SimpleGridT), pointer :: grids       (:) ! Simple rectangular grids
  end type MultiGridT
!
!  The grid itself
!
  type(MultiGridT), save  :: grid
  character(len=20), save :: wallType = 'REFLECTING'
  real(rk), save          :: chargeFlux
!
  contains
!
  subroutine MultiGridInit(max_grids,max_fields,nborder)
    integer(ik), intent(in) :: max_grids   ! Max number of simple grid in the multigrid
    integer(ik), intent(in) :: max_fields  ! Max number of fields on grid - must include enough
                                           ! for scratch fields!
    integer(ik), intent(in) :: nborder     ! Border size

    integer(ik) :: i, alloc

    call TimerStart('MultiGridInit')
    grid%ngrids_max  = max_grids
    grid%nfields_max = max_fields
    grid%nborder     = max(1,nborder)      ! At least one border point is necessary, to keep
                                           ! gradient and laplacian continuous.
    grid%ngrids      = 0
    grid%nfields     = 0

    if (grid%nborder>1) then
      write (out,"(' Pad regions with > 1 (',i2,') points are not supported')") grid%nborder
      stop ' MultiGridInit - overpad'
    end if

    allocate (grid%field_names(max_fields),grid%scratch(max_fields), &
              grid%wavefunction(max_fields), &
              grid%grid_names(max_grids),grid%grids(max_grids),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i8,' initializing the multigrid, max_grids = ',i8,' max_fields = ',i8)") &
              alloc, max_grids, max_fields
       stop 'MultiGridInit - alloc'
    end if

    do i=1,max_grids
      write(grid%grid_names(i),"(' Empty grid ',i8)") i
    end do

    do i=1,max_fields
      write(grid%field_names(i),"(' Scratch field ',i8)") i
    end do

    grid%scratch     (:) = .true.
    grid%wavefunction(:) = .false.

    call TimerStop('MultiGridInit')
  end subroutine MultiGridInit

  subroutine SimpleGridNew(name,npoints,range)
    character(len=*), intent(in) :: name         ! Identifying name for this grid
    integer(ik), intent(in)      :: npoints(  3) ! Number of points in each dimenstion
    real(rk), intent(in)         :: range  (2,3) ! Point ranges for this grid. Will be taken
                                                 ! as-is for the outer grid, but rounded to
                                                 ! the nearest containing point for the inner
                                                 ! grids

    integer(ik)                :: igrid          ! Index for the new grid
    integer(ik)                :: container      ! Index for the containing grid
    type(SimpleGridT), pointer :: gr             ! New grid
    type(SimpleGridT), pointer :: gc             ! Containing grid

    call TimerStart('SimpleGridNew')
    call SimpleGridNewFindContainer(range,container,name) ! Also makes sure the topology is right
    if (container/=0) gc => grid%grids(container)
    !
    !  Find the spot for the new simple grid
    !
    igrid = grid%ngrids + 1
    if (igrid>grid%ngrids_max) then
      write (out,"(' New subgrid ',a,' is one too many - the limit is ',i8)") trim(name), grid%ngrids_max
      stop 'SimpleGridNew - too many'
    end if
    grid%ngrids = igrid
    grid%grid_names(igrid) = name
    gr => grid%grids(igrid)
    gr%up_grid = container
    !
    !  Initialize grid range, and other scalars
    !
    call SimpleGridNewInitRange   (igrid,gr,gc,range,npoints,name)
    call SimpleGridNewInitMomenta (gr,name)
    !
    !  Fill the required arrays
    !
    call SimpleGridNewInitArrays  (gr,name)
    call SimpleGridNewInitDownlink(gr,gc,name)
    !
    !  Be verbose!
    !
    call PrintSimpleGridInfo (igrid)
    call TimerStop('SimpleGridNew')
  end subroutine SimpleGridNew

  subroutine SimpleGridNewFindContainer(range,container,name)
    real(rk), intent(in)         :: range(2,3)  ! Cartesian size of the new grid
    integer(ik), intent(out)     :: container   ! Index of the smallest enclosing grid, or 0
    character(len=*), intent(in) :: name        ! Grid name - for error messages

    integer(ik) :: igrid
    !
    !  First of all, we have to make sure this simple grid does not intersect
    !  the boundary of any of the existing grids - we can only do strictly
    !  hierarchical multigrids
    !
    container = 0
    do igrid=1,grid%ngrids
      if ( disjoint (range,grid%grids(igrid)%range) ) cycle
      if ( contained(range,grid%grids(igrid)%range) ) then
        container = igrid ! This will always pick up the smallest containing grid
        cycle
      end if
      write (out,"(' New subgrid ',a,'(',5(g12.5,','),g12.5," // &
                 "') intersects or contains one of the existing subgrids')") &
            trim(name), range
      stop 'SimpleGridNew - topology'
    end do
  end subroutine SimpleGridNewFindContainer

  subroutine SimpleGridNewInitRange(igrid,gr,gc,range,npoints,name)
    integer(ik), intent(in)      :: igrid
    type(SimpleGridT), pointer   :: gr, gc
    real(rk)                     :: range(2,3)
    integer(ik), intent(in)      :: npoints(3)   ! Number of points in each dimenstion
    character(len=*), intent(in) :: name         ! Identifying name for this grid

    integer(ik)   :: up_npoints(3)               ! Number of patches in the containing section
                                                 ! of the higher-level grid
    integer(ik)   :: match_ratio(3)              ! Number of points per coarser patch
    integer(ik)   :: xl, xu, yl, yu, zl, zu
    !
    !  Initialize new grid, taking care to align it to the nearest point of the
    !  containing grid, and that the spacing of the containing grid is an integer
    !  multiple of the spacing of the enclosed grid. None of this applies to the
    !  top-level grid, obviously.
    !
    if (gr%up_grid==0) then
      gr%npoints   = npoints
      gr%range     = range
      gr%up_index  = 0
    else
      !
      !  For the lower side of the range, we should more the boundary down, to
      !  match the containing grid. For the higher side, we must move up - this
      !  way, the fine grid is guaranteed to enclose all the requested volume
      !
      gr%up_index(1,:) = floor  ((range(1,:)-gc%range(1,:)+gc%step(:)/4)/gc%step(:)) + 1
      gr%up_index(2,:) = ceiling((range(2,:)-gc%range(1,:)-gc%step(:)/4)/gc%step(:))
      if (verbose>=2) then
        write (out,"(' The new grid maps to cells: X = ',2i4,'  Y = ',2i4,'  Z = ',2i4)") gr%up_index
      end if
      if (any(gr%up_index(2,:) - gr%up_index(1,:)<=0)) then
        write (out,"(' Subgrid extent is zero or less: ',6i4)") gr%up_index
        stop 'SimpleGridNew - zero extent'
      end if
      !
      !  Magic "2" below is needed to guarantee a continuous Laplacian
      !
      if ( any(gr%up_index(1,:) < max(2,grid%nborder)) .or. &
           any(gr%up_index(2,:) > gc%npoints-max(2,grid%nborder)) ) then
        write (out,"(' Simple grid ',a,' is placed too close to the boundary of the containing grid ',a)") &
               trim(name), trim(grid%grid_names(gr%up_grid))
        stop 'SimpleGridNew - boundary too close'
      end if
      !
      !  Initialize downward pointers on the containing grid, and kill the
      !  integration weights
      !
      xl = gr%up_index(1,1) ; xu = gr%up_index(2,1)
      yl = gr%up_index(1,2) ; yu = gr%up_index(2,2)
      zl = gr%up_index(1,3) ; zu = gr%up_index(2,3)
      if (any(gc%down_grid(xl:xu,yl:yu,zl:zu)/=0)) then
        write (out,"(' Adjusted subgrids intersect - either join or move apart')")
        stop 'SimpleGridNew - subgrids too close'
      end if
      gc%down_grid(xl:xu,yl:yu,zl:zu) = int(igrid,kind=sik)  ! downlink
      gc%coords (4,xl:xu,yl:yu,zl:zu) = 0                    ! Integration weigth
      !
      !  Get the numerical grid range - don't forget that gc%coords stores coordinates for
      !                                 the patch centres, not for the corners!
      !
      gr%range(1,:) = gc%coords(1:3,xl,yl,zl)-gc%step/2.0_rk
      gr%range(2,:) = gc%coords(1:3,xu,yu,zu)+gc%step/2.0_rk
      if (verbose>=2) then
        write (out,"(' Subgrid extent : ',6(1x,a12))") '  XL  ', '  XU  ', '  YL  ', '  YU  ', '  ZL  ', '  ZU  '
        write (out,"('        Original: ',6(1x,f12.6))") range
        write (out,"('        Adjusted: ',6(1x,f12.6))") gr%range
      end if
      !
      !  Adjust the number of points, making sure the two grids are matched,
      !  calculate the actual grid spacing
      !
      up_npoints   = gr%up_index(2,:) - gr%up_index(1,:) + 1
      match_ratio  = (npoints + up_npoints - 1)/up_npoints
      gr%npoints   = up_npoints * match_ratio
      if (verbose>=2) then
        write (out,"(' Subgrid points : ',3(1x,a4))") ' X ', ' Y ', ' Z '
        write (out,"('        Original: ',3(1x,i4))") npoints
        write (out,"('        Adjusted: ',3(1x,i4))") gr%npoints
      end if
    end if
    !
    !  This part is common for both the top-level and contained grids
    !
    gr%step    = (gr%range(2,:) - gr%range(1,:))/gr%npoints
    gr%weight  = product(gr%step)

  end subroutine SimpleGridNewInitRange

  subroutine SimpleGridNewInitMomenta(gr,name)
    type(SimpleGridT), pointer   :: gr
    character(len=*), intent(in) :: name         ! Identifying name for this grid

    gr%p_step   = twopi / (gr%range(2,:) - gr%range(1,:))
    gr%p_weight = product(gr%p_step)
    !
    !  In our conventions, zero momentum corresponds to the point 1+(N-1)/2
    !  Because we also know the step, this fixes the box dimensions.
    !
    gr%p_range(1,:) = - gr%p_step * ( 1 + (gr%npoints-1)/2 - 0.5_rk )
    gr%p_range(2,:) = gr%p_range(1,:) + gr%p_step * gr%npoints
  end subroutine SimpleGridNewInitMomenta

  subroutine SimpleGridNewInitArrays(gr,name)
    type(SimpleGridT), pointer   :: gr
    character(len=*), intent(in) :: name         ! Identifying name for this grid
    !
    integer(ik)              :: xl, xu, yl, yu, zl, zu
    integer(ik)              :: ix, iy, iz, slot
    integer(ik)              :: alloc
    !
    xl = 1 - grid%nborder ; xu = gr%npoints(1) + grid%nborder
    yl = 1 - grid%nborder ; yu = gr%npoints(2) + grid%nborder
    zl = 1 - grid%nborder ; zu = gr%npoints(3) + grid%nborder
    if (verbose>=1) then
      write (out,"(' Subgrid ',a,' needs ',f9.3,' Mbytes of memory (plus a bit)')") &
             trim(name), rk_bytes*(4+2*grid%nfields_max)*real(xu-xl+1,kind=rk)* &
                         real(yu-yl+1,kind=rk)*real(zu-zl+1,kind=rk)/(1024.0_rk**2)
    end if
    allocate (gr%down_grid(  xl:xu,yl:yu,zl:zu), &
              gr%coords   (4,xl:xu,yl:yu,zl:zu), &
              gr%p_coords (4,xl:xu,yl:yu,zl:zu), &
              gr%active   (2,3,0:grid%nfields_max), &
              gr%fields   (  xl:xu,yl:yu,zl:zu,grid%nfields_max),stat=alloc)
    if (alloc/=0) then
      write (out,"(' Error ',i8,' allocating simple grid ',a)") alloc, trim(name)
      stop 'SimpleGridNew - alloc 1'
    end if
    !
    !  Initialize the arrays
    !
    gr%down_grid = 0
    !
    !  NEC SX compiler generates bad code for fields initialization otherwise!
    !
    zero_fields: do slot=1,grid%nfields_max
      !$omp parallel do private(iz)
      zero_fl_iz: do iz=zl,zu
        gr%fields(:,:,iz,slot) = 0
      end do zero_fl_iz
      !$omp end parallel do
    end do zero_fields
    !
    gr%coords    = 0
    gr%p_coords  = 0
    gr%coords  (4,1:gr%npoints(1),1:gr%npoints(2),1:gr%npoints(3)) = gr%weight
    gr%p_coords(4,1:gr%npoints(1),1:gr%npoints(2),1:gr%npoints(3)) = gr%p_weight
    gr%active  (1,:,:) = spread(spread(1-grid%nborder,1,3),2,grid%nfields_max+1)
    gr%active  (2,:,:) = spread(gr%npoints+grid%nborder,   2,grid%nfields_max+1)
    do iz=zl,zu
      do iy=yl,yu
        do ix=xl,xu
          gr%coords(1:3,ix,iy,iz)   = gr%range  (1,:) + &
                  ( ( (gr%range  (2,:)-gr%range  (1,:)) * ( (/ ix, iy, iz /) - 0.5_rk ) ) / gr%npoints )
          gr%p_coords(1:3,ix,iy,iz) = gr%p_range(1,:) + &
                  ( ( (gr%p_range(2,:)-gr%p_range(1,:)) * ( (/ ix, iy, iz /) - 0.5_rk ) ) / gr%npoints )
        end do
      end do
    end do

  end subroutine SimpleGridNewInitArrays

  subroutine SimpleGridNewInitDownlink(gr,gc,name)
    type(SimpleGridT), pointer   :: gr, gc
    character(len=*), intent(in) :: name         ! Identifying name for this grid
    !
    integer(ik)                  :: idim         ! Dimension - X, Y, Z
    integer(ik)                  :: parent_l     ! Index range within the parent grid
    integer(ik)                  :: parent_u
    integer(ik)                  :: pi           ! Index within the parent grid
    integer(ik)                  :: child_l      ! Fine grid index range, corresponding to parent patch
    integer(ik)                  :: child_u
    integer(ik)                  :: alloc
    integer(ik)                  :: p_step       ! Step of the fine grid for each point of the coarse grid
    !
    !  Build the correspondence table for each affected grid point higher up
    !
    if (gr%up_grid==0) return

    parent_l = minval(gr%up_index(1,:) - grid%nborder) ! It may give a few more points than strictly
    parent_u = maxval(gr%up_index(2,:) + grid%nborder) ! necessary (e.g. if grid%border is more than one, but
                                                       ! grids have identical spacing), but there is no harm
                                                       ! in this.
    allocate (gr%down_index(2,parent_l:parent_u,3),stat=alloc)
    if (alloc/=0) then
      write (out,"(' Error ',i8,' allocating simple grid downlink ',a)") alloc, trim(name)
      stop 'SimpleGridNew - alloc downlink'
    end if
    gr%down_index = 0
    do idim=1,3
      parent_l = gr%up_index(1,idim)
      parent_u = gr%up_index(2,idim)
      p_step   = gr%npoints(idim) / ( parent_u - parent_l + 1 )
      do pi=parent_l-grid%nborder, parent_u+grid%nborder
        child_l = (pi-parent_l)*p_step + 1
        child_u = child_l + p_step - 1
        gr%down_index(1,pi,idim) = max(child_l,               1-grid%nborder)
        gr%down_index(2,pi,idim) = min(child_u,gr%npoints(idim)+grid%nborder)
      end do
      !
      !  Sanify check - the whole lower-level grid must be covered, and all spacings
      !                 must be identical
      !
      if (gr%down_index(1,gr%up_index(1,idim),idim)/=1 .or.  &
          gr%down_index(2,gr%up_index(2,idim),idim)/=gr%npoints(idim) ) then
        write (out,"(' Internal logic error constructing grid down_index: '"//   &
                   ",i5,' must be ',i5,'; ',i5,' must be ',i5)")                 &
               gr%down_index(1,gr%up_index(1,idim),idim), 1,                     &
               gr%down_index(2,gr%up_index(2,idim),idim), gr%npoints(idim)
        stop 'SimpleGridNew - down_index insane 1'
      end if
      if (any(gr%down_index(2,parent_l:parent_u,idim)-gr%down_index(1,parent_l:parent_u,idim)+1-p_step/=0)) then
        write (out,"(' Internal login error constructing grid down_index - unexpected spacings!')")
        write (out,*) gr%down_index(2,parent_l:parent_u,idim)-gr%down_index(1,parent_l:parent_u,idim)+1
        stop 'impleGridNew - down_index insane 2'
      end if
    end do
  end subroutine SimpleGridNewInitDownlink

  subroutine PrintSimpleGridInfo (igrid)
    integer(ik), intent(in)     :: igrid
    !
    type (SimpleGridT), pointer :: gr
    integer(ik)                 :: idim
    integer(ik)                 :: ix, xl, xu
    integer(ik)                 :: iy, yl, yu
    integer(ik)                 :: iz, zl, zu
    !
    if (verbose<=0) return
    !
    !  Report grid parameters
    !
    gr => grid%grids(igrid)
    write (out,"(' Added simple grid #',i8,' (',a,')')") igrid, trim(grid%grid_names(igrid))
    write (out,"(' npoints   = ',3i6)")       gr%npoints
    write (out,"(' range     = ',6f10.5)")    gr%range
    write (out,"(' step      = ',3f10.5)")    gr%step
    write (out,"(' weight    = ', f10.5)")    gr%weight
    write (out,"(' parent    = ', i4)")       gr%up_grid
    write (out,"(' up_range  = ', 6i5)")      gr%up_index
    write (out,"(' p_range   = ',6f10.5)")    gr%p_range
    write (out,"(' p_step    = ',3f10.5)")    gr%p_step
    write (out,"(' p_weight  = ',3f10.5)")    gr%p_weight
    write (out,"(' active(0) = ',3(2i5,2x))") gr%active(:,:,0)
    if (verbose>=3 .and. gr%up_grid/=0) then
      do idim=1,3
        xl = gr%up_index(1,idim) - grid%nborder
        xu = gr%up_index(2,idim) + grid%nborder
        write (out,"(' Correspondence table for dimension ',i1)") idim
        write (out,"(('   Upper = ',i4,' Lower = ',i4,' to ',i4))") ( ix, gr%down_index(:,ix,idim), ix=xl,xu)
      end do
    end if
    if (verbose>=4) then
      xl = 1 - grid%nborder ; xu = gr%npoints(1) + grid%nborder
      yl = 1 - grid%nborder ; yu = gr%npoints(2) + grid%nborder
      zl = 1 - grid%nborder ; zu = gr%npoints(3) + grid%nborder
      write (out, "(3(1x,a3),2x,a2,1x,3(1x,a10  ),1x,a15   )") &
            'ix', 'iy', 'iz', 'dn', ' X ', ' Y ', ' Z ', ' wgt '
      write (out,"((3(1x,i3),2x,i2,2(1x,3(1x,f10.5),1x,f15.10)))") &
            (((ix,iy,iz,gr%down_grid(ix,iy,iz),gr%coords(:,ix,iy,iz), &
                        gr%p_coords(:,ix,iy,iz),ix=xl,xu),iy=yl,yu),iz=zl,zu)
    end if
  end subroutine PrintSimpleGridInfo

  !
  !  Allocate a new field
  !
  subroutine FieldNew(name,index,scratch,wavefunction)
    character(len=*), intent(in) :: name
    integer(ik), intent(out)     :: index
    logical, intent(in)          :: scratch      ! True if the grid is used for scratch, and
                                                 ! should not be included in checkpoints
    logical, intent(in)          :: wavefunction ! True if grid contains a wavefunction

    if (grid%nfields>=grid%nfields_max) then
      write (out,"(' Can''t allocate new field ',a,' - the maximum (',i3,') is already reached')") &
             trim(name), grid%nfields_max
      stop 'FieldNew - too many'
    end if
    grid%nfields = grid%nfields + 1
    index        = grid%nfields
    grid%field_names (index) = name
    grid%scratch     (index) = scratch
    grid%wavefunction(index) = wavefunction .and. allowPolar
  end subroutine FieldNew
  !
  !  Set type of the field, used to control the behaviour of interpolation
  !  routines.
  !
  subroutine FieldSetWavefunction(index,wavefunction)
    integer(ik), intent(in) :: index
    logical, intent(in)     :: wavefunction
    !
    if (verbose>=2) write (out,"(' = FieldSetWavefunction ',i3,' = ',l6)") &
                           index, wavefunction
    !
    grid%wavefunction(index) = wavefunction .and. allowPolar
  end subroutine FieldSetWavefunction
  !
  !  The job of FieldDump is to write out the complete definition of a field,
  !  together with information on the grid. Dumps can be -very- large!
  !
  subroutine FieldDump(cout,src)
    integer(ik), intent(in)    :: cout ! Output I/O unit to use
    integer(ik), intent(in)    :: src  ! Field to dump; grid is dumped as well

    integer(ik)                :: xu,  yu,  zu       ! Point ranges, with no borders
    integer(ik)                :: xub, yub, zub      ! Point ranges, with borders included
    integer(ik)                :: nb                 ! Shorthand for grid%nborder
    integer(ik)                :: igrid, ix, iy, iz
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: itype              ! 0 for major point; n for average; -n for border
    !
    call TimerStart('FieldDump')
    if (verbose>=2) write (out,"(' = FieldDump ',i3)") src
    !
    nb = grid%nborder
    write (cout,"('# Field ',i3,': ',a)") src, trim(grid%field_names(src))
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      xu = gr%npoints(1) ; xub = xu + nb
      yu = gr%npoints(2) ; yub = yu + nb
      zu = gr%npoints(3) ; zub = zu + nb
      write (cout,"('# Grid ',i3,': ',a)") igrid, trim(grid%grid_names(igrid))
      write (cout,"('#               ',3a15)") '    X    ','    Y    ','    Z    '
      write (cout,"('#               ',3a15)") '---------','---------','---------'
      write (cout,"('#   Points:     ',3i15  )") gr%npoints
      write (cout,"('#   Min. coord: ',3f15.8)") gr%range(1,:)
      write (cout,"('#   Max. coord: ',3f15.8)") gr%range(2,:)
      write (cout,"('#   Step:       ',3f15.8)") gr%step
      write (cout,"('#',3(1x,a10),2x,a12,2x,2(1x,a15),2(1x,a3))") &
             ' X ', ' Y ', ' Z ', ' WGT ', ' RE ', ' IMAG ', 'TYP', 'GRD'
      do iz=1-nb,zub
        do iy=1-nb,yub
          do ix=1-nb,xub
            if (ix<1 .or. ix>xu .or. iy<1 .or. iy>yu .or. iz<1 .or. iz>zu ) then
              itype = -(1+gr%up_grid)
            else
              itype = gr%down_grid(ix,iy,iz)
            end if
            write (cout,"(3(1x,f10.5),2x,f12.6,2x,2(1x,f15.8),2(1x,i3))") &
                   gr%coords(:,ix,iy,iz), gr%fields(ix,iy,iz,src), itype, igrid
          end do
        end do
      end do
    end do
    call TimerStop('FieldDump')
  end subroutine FieldDump

  !
  !  This chooses either reflecting, or adsorbing outer wall
  !
  subroutine FieldSetOuterWall(type)
    character(len=*), intent(in) :: type  ! This can be 'REFLECTING' or 'ADSORBING'

    select case (type)
      case default
        write (out,"(' FieldSetOuterWall - bad wall type ',a)") trim(type)
        stop 'FieldSetOuterWall - bad wall type'
      case ('REFLECTING')
      case ('ADSORBING')
    end select
    !
    wallType = type
    chargeFlux = 0
  end subroutine FieldSetOuterWall

  !
  !  When wallType is 'ADSORBING', FieldGetLeakage will return
  !  total charge flux through the outer wall of the box.
  !  The flux is determined at time of grid reconciliation
  !  (which means FieldLaplacian, for the time being)
  !
  function FieldGetLeakage() result(v)
    real(rk) :: v

    v = chargeFlux
  end function FieldGetLeakage

  !
  !  The job of FieldVisualize is to communicate specified field to an
  !  external visualization program (currently, OpenDX). If we are not
  !  running with visualization turned on, this is a no-op. On it's own,
  !  FieldVisualize does -not- display anything; it simply prepares
  !  data to be sent over to OpenDX. The actual visualization is done
  !  at the call to FieldShow.
  !
  !  To prevent DX running out of memory, we have to be careful with
  !  the amount of data we send over. Therefore, we'll send the probality
  !  density instead of the wavefunction. We'll also downsample the densities
  !  (taking care to preserve normalizations) if that appears necessary.
  !
  subroutine FieldVisualize(slot,src,text,level,field,replicate,limits)
    integer(ik), intent(in)                :: slot        ! Output channel to show this field on
                                                          ! For even channels (0, 2, etc), the field is assumed to be spatial
                                                          ! For odd channels (1, 3, etc), it is in the momentum representation
    integer(ik), intent(in)                :: src         ! Field to visualize
    character(len=*), intent(in)           :: text        ! Descriptive name of the field
    integer(ik), intent(in), optional      :: level       ! Grid level to visualize. If omitted,
                                                          ! the outermost (coarsest) grid is used.
    real(rk), intent(in), optional         :: field(3)    ! Vector "external driver" value to plot
    character(len=*), intent(in), optional :: replicate   ! Ignored; needed for compatibility with 2D case.
    integer(ik), intent(in), optional      :: limits(:,:) ! Optional grid limits

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    real(rk), pointer          :: cr(:,:,:,:)
    integer(ik)                :: active, alloc
    integer(ik)                :: ix, iy, iz, ox, oy, oz
    integer(ik)                :: oxu, oyu, ozu
    integer(ik)                :: xl, yl, zl     ! Grid extent, including possibility of dumping 
    integer(ik)                :: xu, yu, zu     ! just a part of the grid
    integer(ik)                :: div            ! Downsampling amount
    real(rk)                   :: scale          ! Downscaling coefficient
    integer(ik)                :: o_pnt (  3)    ! Original size of the data mesh
    integer(ik)                :: c_pnt (  3)    ! Original size of the data mesh, commesurate with
                                                 ! the downsampled mesh
    integer(ik)                :: m_pnt (  3)    ! Downsampled mesh dimensions
    real(srk)                  :: m_corners(2,3) ! Extreme positions of the simulation mesh
    real(srk)                  :: m_field(3)     ! Field vlaue to plot
                                                 ! m_corners and m_field are communicated to C, so they
                                                 ! must be kind- less real
!   real(srk), allocatable     :: m_data(:,:,:)  ! (Possibly downsampled) probability density to plot
    real(srk), allocatable     :: m_data(:,:,:,:)! (Possibly downsampled) probability density to plot
    !
    call TimerStart('FieldVisualize')
    !
    !  Get optional grid level parameter, or take the default (1) if it
    !  is absent.
    !
    igrid = 1
    if (present(level)) igrid = level
    m_field = 0
    if (present(field)) m_field = real(field) ! Convert to default real kind
    !
    if (verbose>=2) then
      write (out,"(' = FieldVisualize ',i3,1x,i3,'(',i3,')',a)") slot, src, igrid, trim(text)
    end if
    if (igrid>grid%ngrids) then
      write (out,"(' FieldVisualize - asked to plot non-existent grid ',i8)") igrid
      stop 'FieldVisualize - non-existent grid'
    end if

    active = 0
    call drvGridHead(active)
    if (active==0) then
      call TimerStop('FieldVisualize')
      return
    end if
    !
    !  We'll end up here only if there is a connection to the visualization
    !  program, so that we can meaningfully send some data over. We'll send
    !  the complete simulation mesh, including the borders. Because our data
    !  structures record mesh parameters without the borders, we'll need to
    !  fill out the full mesh parameters first.
    !
    !  Note that our data are in Fortran order - first index changes fastest.
    !  DX works in C, so it expects the last index to change fastest. Therefore,
    !  we have to reverse the order of the X, Y, Z mesh and corners here. The
    !  same applies for vector driving field (m_field) below.
    !
    gr => grid%grids(igrid)
    !
    xl = 1 ; xu = gr%npoints(1)
    yl = 1 ; yu = gr%npoints(2)
    zl = 1 ; zu = gr%npoints(3)
    if (present(limits)) then
      xl = max(xl,limits(1,1)) ; xu = min(xu,limits(2,1))
      yl = max(yl,limits(1,2)) ; yu = min(yu,limits(2,2))
      zl = max(zl,limits(1,3)) ; zu = min(zu,limits(2,3))
    end if
    !
    !  Decide on the number of points and downsampling
    !
    div = 1
    o_pnt = (/ xu-xl, yu-yl, zu-zl /) + 1      ! Total number of points
    do while(product(o_pnt/div)>maxDXgridSize)
      div = div + 1
    end do
    !
    m_pnt (3:1:-1)   = o_pnt / div             ! Downsampled number of points
                                               ! (in the reverse order - damn C!)
    c_pnt (3:1:-1)   = m_pnt * div             ! Commesurate number of points in
                                               ! the main grid
    if (verbose>=1 .and. div/=1) then
      write (out,"(' FieldVisualize: Grid downsampled by: ',i3)") div
      write (out,"(' FieldVisualize: Downsampled grid size: ',3i4)") m_pnt(3:1:-1)
      write (out,"(' FieldVisualize: Commesurate grid size: ',3i4)") c_pnt
    end if
    !
    cr => gr%coords(1:3,xl:xl+c_pnt(1)-1,yl:yl+c_pnt(2)-1,zl:zl+c_pnt(3)-1)   ! Real-space grids
    if (mod(slot,2)==1) &
    cr => gr%p_coords(1:3,xl:xl+c_pnt(1)-1,yl:yl+c_pnt(2)-1,zl:zl+c_pnt(3)-1) ! Momentum-space grids
    !
    !  Coordinates of the corners are centres of the downsampled patches,
    !  which is identical to the averages of the coordinates of the centres
    !  of the constituent patches. Note the reverse order of the coordinates
    !  (damn C again!)
    !
    m_corners(1,3) = real(sum(cr(1,                 1:1+div-1, 1,1))/div)   ! Convert to default real kind
    m_corners(1,2) = real(sum(cr(2,1,               1:1+div-1, 1  ))/div)
    m_corners(1,1) = real(sum(cr(3,1,1,             1:1+div-1     ))/div)
    m_corners(2,3) = real(sum(cr(1,    c_pnt(1)-div+1:c_pnt(1),1,1))/div)
    m_corners(2,2) = real(sum(cr(2,1,  c_pnt(2)-div+1:c_pnt(2),1  ))/div)
    m_corners(2,1) = real(sum(cr(3,1,1,c_pnt(3)-div+1:c_pnt(3)    ))/div)
    if (verbose>1) then
      write (out,"(' Corners of the reduced grid: '/'  X: ',2f10.5/'  Y: ',2f10.5/'  Z: ',2f10.5)") &
             m_corners(:,3:1:-1)
    end if
    !
    !  For real-space grids, data must be reconciled with finer grids first
    !
    if (mod(slot,2)==0) call reconcileGrids(src)
    !
    !  Convert wavefunction to densities, and downsample it
    !
!   allocate (m_data(1:m_pnt(3),1:m_pnt(2),1:m_pnt(1)),stat=alloc)
    allocate (m_data(2,1:m_pnt(3),1:m_pnt(2),1:m_pnt(1)),stat=alloc)
    if (alloc/=0) then
      write (out,"(' FieldVisualize - can''t allocate ',i4,'-word buffer')") &
             product(m_pnt)
      stop 'FieldVisualize - no buffer memory'
    end if
    !
    scale = 1.0_rk / div**3
    !$omp parallel do private(ix,iy,iz,oz,ozu,oy,oyu,ox,oxu)
    do iz=1,m_pnt(1)
      !
      !  First of all, zero the output slice - it may never be
      !  initialized otherwise.
      !
!     m_data(:,:,iz) = 0
      m_data(:,:,:,iz) = 0
      oz  = (iz-1)*div + zl ; ozu = oz+div-1
      oz  = max(gr%active(1,3,src),oz )  ! Anything smaller is zero
      ozu = min(gr%active(2,3,src),ozu)  ! Anything larger is zero
      if (ozu.lt.oz) cycle
      do iy=1,m_pnt(2)
        oy  = (iy-1)*div + yl ; oyu = oy+div-1
        oy  = max(gr%active(1,2,src),oy )  ! Anything smaller is zero
        oyu = min(gr%active(2,2,src),oyu)  ! Anything larger is zero
        if (oyu.lt.oy) cycle
        do ix=1,m_pnt(3)
          ox  = (ix-1)*div + xl ; oxu = ox+div-1
          ox  = max(gr%active(1,1,src),ox )  ! Anything smaller is zero
          oxu = min(gr%active(2,1,src),oxu)  ! Anything larger is zero
          if (oxu.lt.ox) cycle
!         m_data(ix,iy,iz) = scale * sum(abs(gr%fields(ox:oxu,oy:oyu,oz:ozu,src))**2)
          m_data(1,ix,iy,iz) = real(scale * sum( real(gr%fields(ox:oxu,oy:oyu,oz:ozu,src))))
          m_data(2,ix,iy,iz) = real(scale * sum(aimag(gr%fields(ox:oxu,oy:oyu,oz:ozu,src))))
        end do
      end do
    end do
    !$omp end parallel do
    !
    !  Send data over to DX, and release the buffer
    !
    call drvGridData(slot,m_pnt,m_corners,m_field(3:1:-1),m_data,len_trim(text),text)
    deallocate (m_data)
    call TimerStop('FieldVisualize')

  end subroutine FieldVisualize

  !
  !  Tell OpenDX that the data is ready for pick-up
  !
  subroutine FieldShow
    integer(ik)                :: active

    call TimerStart('FieldShow')
    active = 0
    call drvGridHead(active)
    if (active==0) then
      call TimerStop('FieldShow')
      return
    end if
    call drvGridShow
    call TimerStop('FieldShow')

  end subroutine FieldShow

  !
  !  Returns square root of the integral of a field with its conjugate.
  !  Very handy for computing 2-norms :-)
  !
  function FieldNorm(src) result(v)
    integer(ik), intent(in) :: src
    real(rk)                :: v

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu
    real(rk)                   :: vp

    call TimerStart('FieldNorm')
    if (verbose>=2) write (out,"(' = FieldNorm ',i3)") src
    !
    v = 0
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = gr%active(1,3,src) ; izu = gr%active(2,3,src)
      xl  = gr%active(1,1,src) ; xu  = gr%active(2,1,src)
      yl  = gr%active(1,2,src) ; yu  = gr%active(2,2,src)
      vp  = 0
      !$omp parallel do reduction(+:vp) private(iz)
      do iz=izl,izu
        vp = vp + sum(gr%coords(4,xl:xu,yl:yu,iz)     * &
                 real(gr%fields(  xl:xu,yl:yu,iz,src) * &
                conjg(gr%fields(  xl:xu,yl:yu,iz,src)),kind=rk))
      end do
      !$omp end parallel do
      v = v + vp
    end do
    !
    if (v<0) v = 0
    v = sqrt(v)
    call TimerStop('FieldNorm')
    !
  end function FieldNorm
  !
  !  Returns integral of the field elements
  !
  function FieldNorm1(src) result(v)
    integer(ik), intent(in) :: src
    complex(rk)             :: v

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu
    complex(rk)                :: vp

    call TimerStart('FieldNorm1')
    if (verbose>=2) write (out,"(' = FieldNorm1 ',i3)") src
    !
    v = 0
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = gr%active(1,3,src) ; izu = gr%active(2,3,src)
      xl  = gr%active(1,1,src) ; xu  = gr%active(2,1,src)
      yl  = gr%active(1,2,src) ; yu  = gr%active(2,2,src)
      vp  = 0
      !$omp parallel do reduction(+:vp) private(iz)
      do iz=izl,izu
        vp = vp + sum(gr%coords(4,xl:xu,yl:yu,iz) * &
                      gr%fields(  xl:xu,yl:yu,iz,src))
      end do
      !$omp end parallel do
      v = v + vp
    end do
    !
    call TimerStop('FieldNorm1')
    !
  end function FieldNorm1
  !
  !  Returns multipole moments of the squared absolute values of the
  !  field. If the field happens to be a waefunction, this will give
  !  moments of the electron density distribution.
  !
  !  The order of the multipoles is:
  !    1 = X ; 2 = Y ; 3 = Z ;
  !    4 = XX ; 5 = YY ; 6 = ZZ ; 7 = XY ; 8 = XZ ; 9 = YZ
  !
  !  For nuclear multipoles, see FLnuclearMultipoles
  !
  subroutine FieldNormMultipoles(src,mult)
    integer(ik), intent(in) :: src
    real(rk),intent(out)    :: mult(:)

    integer(ik)                :: nmult
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, zl, zu
    integer(ik)                :: ix, xl, xu
    integer(ik)                :: iy, yl, yu
    real(ark)                  :: vd(3), ve(3), vf(3)
    real(rk)                   :: wgl
    !
    call TimerStart('FieldNormMultipoles')
    if (verbose>=2) write (out,"(' = FieldNormMultipoles ',i3)") src
    !
    mult  = 0
    nmult = size(mult)
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      zl = gr%active(1,3,src) ; zu = gr%active(2,3,src)
      xl = gr%active(1,1,src) ; xu = gr%active(2,1,src)
      yl = gr%active(1,2,src) ; yu = gr%active(2,2,src)
      !
      !  The goal of the loop structure is to minimize the number of passes
      !  over memory, while simultaneously eliminating inner-loop branches.
      !  A (good) optimizing compiler should be able to do the transformation
      !  itself, but ...
      !
      if (nmult<6) then
        !
        !  Case 1: Dipoles only
        !
        vd = 0
        !$omp parallel do reduction(+:vd) private(ix,iy,iz,wgl)
        dipoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * real(gr%fields(ix,iy,iz,src) * conjg(gr%fields(ix,iy,iz,src)),kind=rk)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
            end do
          end do
        end do dipoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
      else if (nmult<9) then
        !
        !  Case 2: Dipoles and diagonal of quadrupoles
        !
        vd = 0 ; ve = 0
        !$omp parallel do reduction(+:vd,ve) private(ix,iy,iz,wgl)
        diagonal: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * real(gr%fields(ix,iy,iz,src) * conjg(gr%fields(ix,iy,iz,src)),kind=rk)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
            end do
          end do
        end do diagonal
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
      else
        !
        !  Case 3: Dipoles and full quadripoles
        !
        vd = 0 ; ve = 0 ; vf = 0
        !$omp parallel do reduction(+:vd,ve,vf) private(ix,iy,iz,wgl)
        quadrupoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * real(gr%fields(ix,iy,iz,src) * conjg(gr%fields(ix,iy,iz,src)),kind=rk)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
              vf  = vf + wgl * (/ gr%coords(1,ix,iy,iz)*gr%coords(2,ix,iy,iz),   &
                                  gr%coords(1,ix,iy,iz)*gr%coords(3,ix,iy,iz),   &
                                  gr%coords(2,ix,iy,iz)*gr%coords(3,ix,iy,iz) /)
            end do
          end do
        end do quadrupoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
        mult(7:9) = mult(7:9) + vf(1:3)
      end if
    end do
    call TimerStop('FieldNormMultipoles')
    !
  end subroutine FieldNormMultipoles
  !
  !  Same as FieldNormMultipoles, except that we already have the density
  !  The density is allowed to be complex (think exchange-type integrals!),
  !  and so are the multipoles.
  !
  subroutine FieldNorm1Multipoles(src,mult)
    integer(ik), intent(in) :: src
    complex(rk),intent(out) :: mult(:)
    !
    integer(ik)                :: nmult
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, zl, zu
    integer(ik)                :: ix, xl, xu
    integer(ik)                :: iy, yl, yu
    complex(ark)               :: vd(3), ve(3), vf(3)
    complex(rk)                :: wgl
    !
    call TimerStart('FieldNorm1Multipoles')
    if (verbose>=2) write (out,"(' = FieldNorm1Multipoles ',i3)") src
    !
    mult  = 0
    nmult = size(mult)
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      zl = gr%active(1,3,src) ; zu = gr%active(2,3,src)
      xl = gr%active(1,1,src) ; xu = gr%active(2,1,src)
      yl = gr%active(1,2,src) ; yu = gr%active(2,2,src)
      !
      !  The goal of the loop structure is to minimize the number of passes
      !  over memory, while simultaneously eliminating inner-loop branches.
      !  A (good) optimizing compiler should be able to do the transformation
      !  itself, but ...
      !
      if (nmult<6) then
        !
        !  Case 1: Dipoles only
        !
        vd = 0
        !$omp parallel do reduction(+:vd) private(ix,iy,iz,wgl)
        dipoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * gr%fields(ix,iy,iz,src)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
            end do
          end do
        end do dipoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
      else if (nmult<9) then
        !
        !  Case 2: Dipoles and diagonal of quadrupoles
        !
        vd = 0 ; ve = 0
        !$omp parallel do reduction(+:vd,ve) private(ix,iy,iz,wgl)
        diagonal: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * gr%fields(ix,iy,iz,src)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
            end do
          end do
        end do diagonal
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
      else
        !
        !  Case 3: Dipoles and full quadripoles
        !
        vd = 0 ; ve = 0 ; vf = 0
        !$omp parallel do reduction(+:vd,ve,vf) private(ix,iy,iz,wgl)
        quadrupoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl = gr%coords(4,ix,iy,iz) * gr%fields(ix,iy,iz,src)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
              vf  = vf + wgl * (/ gr%coords(1,ix,iy,iz)*gr%coords(2,ix,iy,iz),   &
                                  gr%coords(1,ix,iy,iz)*gr%coords(3,ix,iy,iz),   &
                                  gr%coords(2,ix,iy,iz)*gr%coords(3,ix,iy,iz) /)
            end do
          end do
        end do quadrupoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
        mult(7:9) = mult(7:9) + vf(1:3)
      end if
    end do
    call TimerStop('FieldNorm1Multipoles')
    !
  end subroutine FieldNorm1Multipoles
  !
  !  Calculate integrals in the form:
  !
  !    <bra|M|ket>
  !
  !  where M is a dipole or quadrupole operator.
  !
  subroutine FieldNorm2Multipoles(left,right,mult)
    integer(ik), intent(in) :: left, right
    complex(rk),intent(out) :: mult(:)
    !
    integer(ik)                :: nmult
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, zl, zu
    integer(ik)                :: ix, xl, xu
    integer(ik)                :: iy, yl, yu
    complex(ark)               :: vd(3), ve(3), vf(3)
    complex(rk)                :: wgl
    !
    call TimerStart('FieldNorm2Multipoles')
    if (verbose>=2) write (out,"(' = FieldNorm2Multipoles ',i3,1x,i3)") left, right
    !
    mult  = 0
    nmult = size(mult)
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      xl = minval(gr%active(1,1,(/left,right/))) 
      yl = minval(gr%active(1,2,(/left,right/))) 
      zl = minval(gr%active(1,3,(/left,right/))) 
      xu = maxval(gr%active(2,1,(/left,right/)))
      yu = maxval(gr%active(2,2,(/left,right/)))
      zu = maxval(gr%active(2,3,(/left,right/)))
      !
      call FieldAugmentBox(left,igrid,xl,xu,yl,yu,zl,zu)
      call FieldAugmentBox(right,igrid,xl,xu,yl,yu,zl,zu)
      !
      !  The goal of the loop structure is to minimize the number of passes
      !  over memory, while simultaneously eliminating inner-loop leftnches.
      !  A (good) optimizing compiler should be able to do the transformation
      !  itself, but ...
      !
      if (nmult<6) then
        !
        !  Case 1: Dipoles only
        !
        vd = 0
        !$omp parallel do reduction(+:vd) private(ix,iy,iz,wgl)
        dipoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl =       gr%coords(4,ix,iy,iz)      &
                  * conjg(gr%fields(  ix,iy,iz,left)) &
                  *       gr%fields(  ix,iy,iz,right)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
            end do
          end do
        end do dipoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
      else if (nmult<9) then
        !
        !  Case 2: Dipoles and diagonal of quadrupoles
        !
        vd = 0 ; ve = 0
        !$omp parallel do reduction(+:vd,ve) private(ix,iy,iz,wgl)
        diagonal: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl =       gr%coords(4,ix,iy,iz)      &
                  * conjg(gr%fields(  ix,iy,iz,left)) &
                  *       gr%fields(  ix,iy,iz,right)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
            end do
          end do
        end do diagonal
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
      else
        !
        !  Case 3: Dipoles and full quadripoles
        !
        vd = 0 ; ve = 0 ; vf = 0
        !$omp parallel do reduction(+:vd,ve,vf) private(ix,iy,iz,wgl)
        quadrupoles: do iz=zl,zu
          do iy=yl,yu
            do ix=xl,xu
              wgl =       gr%coords(4,ix,iy,iz)      &
                  * conjg(gr%fields(  ix,iy,iz,left)) &
                  *       gr%fields(  ix,iy,iz,right)
              vd  = vd + wgl * gr%coords(1:3,ix,iy,iz)
              ve  = ve + wgl * gr%coords(1:3,ix,iy,iz)**2
              vf  = vf + wgl * (/ gr%coords(1,ix,iy,iz)*gr%coords(2,ix,iy,iz),   &
                                  gr%coords(1,ix,iy,iz)*gr%coords(3,ix,iy,iz),   &
                                  gr%coords(2,ix,iy,iz)*gr%coords(3,ix,iy,iz) /)
            end do
          end do
        end do quadrupoles
        !$omp end parallel do
        mult(1:3) = mult(1:3) + vd(1:3)
        mult(4:6) = mult(4:6) + ve(1:3)
        mult(7:9) = mult(7:9) + vf(1:3)
      end if
    end do
    call FieldShrinkBox(left)
    call FieldShrinkBox(right)
    call TimerStop('FieldNorm2Multipoles')
    !
  end subroutine FieldNorm2Multipoles
  !
  !  Calculates an integral of the product of two functions
  !
  function FieldProductIntegrate(src_a,src_b) result(v)
    integer(ik), intent(in) :: src_a  ! First field in the product
    integer(ik), intent(in) :: src_b  ! Second field in the product
    complex(rk)             :: v

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu
    complex(rk)                :: vp

    !
    call TimerStart('FieldProductIntegrate')
    if (verbose>=2) write (out,"(' = FieldProductIntegrate ',i3,' with ',i3)") src_a, src_b
    !
    v = 0
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = maxval(gr%active(1,3,(/src_a,src_b/)))
      izu = minval(gr%active(2,3,(/src_a,src_b/)))
      xl  = maxval(gr%active(1,1,(/src_a,src_b/)))
      xu  = minval(gr%active(2,1,(/src_a,src_b/)))
      yl  = maxval(gr%active(1,2,(/src_a,src_b/)))
      yu  = minval(gr%active(2,2,(/src_a,src_b/)))
      vp  = 0
      !$omp parallel do reduction(+:vp) private(iz)
      do iz=izl,izu
        vp = vp + sum(gr%coords(4,xl:xu,yl:yu,iz)* &
                      gr%fields(xl:xu,yl:yu,iz,src_b)* &
                      gr%fields(xl:xu,yl:yu,iz,src_a))
      end do
      !$omp end parallel do
      v = v + vp
    end do
    call TimerStop('FieldProductIntegrate')
    !
  end function FieldProductIntegrate
  !
  !  Calculates <left|right> braket.
  !
  function FieldConjgIntegrate(left,right) result(v)
    integer(ik), intent(in) :: left   ! Field to be conjugated
    integer(ik), intent(in) :: right  ! Field to appear as-is
    complex(rk)             :: v

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu
    complex(rk)                :: vp

    !
    call TimerStart('FieldConjgIntegrate')
    if (verbose>=2) write (out,"(' = FieldConjgIntegrate ',i3,' with ',i3)") left, right
    !
    v = 0
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = maxval(gr%active(1,3,(/left,right/)))
      izu = minval(gr%active(2,3,(/left,right/)))
      xl  = maxval(gr%active(1,1,(/left,right/)))
      xu  = minval(gr%active(2,1,(/left,right/)))
      yl  = maxval(gr%active(1,2,(/left,right/)))
      yu  = minval(gr%active(2,2,(/left,right/)))
      vp  = 0
      !$omp parallel do reduction(+:vp) private(iz)
      do iz=izl,izu
        vp = vp + sum(gr%coords(4,xl:xu,yl:yu,iz)* &
                      gr%fields(xl:xu,yl:yu,iz,right)* &
                conjg(gr%fields(xl:xu,yl:yu,iz,left)))
      end do
      !$omp end parallel do
      v = v + vp
    end do
    call TimerStop('FieldConjgIntegrate')
    !
  end function FieldConjgIntegrate
  !
  !  Integrate equation in the form:
  !
  !     k . Grad P + U = 0
  !
  !   where k is a 3-vector, P and U are scalar fields.
  !   Boundary conditions for P should already be set at the two points
  !   on each side of the extreme index values along each Cartesian direction.
  !
  !   WARNING: this a very, very, very long routine
  !
  subroutine FieldBuildPhase(pot,phase,kvec)
    integer(ik), intent(in)    :: pot     ! Field containing U; not modified
    integer(ik), intent(in)    :: phase   ! Field containing P; boundary counditions must be set;
                                          ! the rest of the field will be filled on output.
    real(rk), intent(in)       :: kvec(3) ! K-vector
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz, ind, indm(1)
    integer(ik)                :: xl, xu, yl, yu, zl, zu, incr
    integer(ik)                :: r_order(3), Wu_order(3) ! Ordering sequence for the intersection
                                                          ! point and integration weights matrix
    real(rk)                   :: r (3)  ! Intersection point in the "lower" plane.
    real(rk)                   :: r1(3)  ! Intersection point, reordered to shift integration direction
                                         ! to the first two positions.
    real(rk)                   :: Wg(3,3)   ! Interpolation weights
    real(rk)                   :: Wu(3,3,3) ! Integration weights
    !
    call TimerStart('FieldBuildPhase')
    if (verbose>=2) write (out,"(' = FieldBuildPhase pot=',i3,' phase = ',i3,' kvec = ',3g12.5)") pot, phase, kvec
    !
    if (grid%ngrids/=1) stop 'FieldBuildPhase - multigrids not supported'
    if (pot==phase) stop 'FieldBuildPhase - oops: input = output'
    !
    gr => grid%grids(1)
    xl = min(gr%active(1,1,phase),gr%active(1,1,pot))
    xu = max(gr%active(2,1,phase),gr%active(2,1,pot))
    yl = min(gr%active(1,2,phase),gr%active(1,2,pot))
    yu = max(gr%active(2,2,phase),gr%active(2,2,pot))
    zl = min(gr%active(1,3,phase),gr%active(1,3,pot))
    zu = max(gr%active(2,3,phase),gr%active(2,3,pot))
    !
    call FieldAugmentBox(pot  ,1_ik,xl,xu,yl,yu,zl,zu)
    call FieldAugmentBox(phase,1_ik,xl,xu,yl,yu,zl,zu)
    !
    !  The algorithm changes depending on the k-vector direction.
    !  The integration is always along the Cartesian direction most aligned
    !  with the k-vector. Additionally, integration starts from the flat-front
    !  side, and proceeds towards the scattering potential - either upwards
    !  or downwards, depending on the sign of the largest component of k.
    !
    !  Find the largest component of the k-vector
    !
    indm = maxloc(abs(kvec)) ; ind=indm(1)
    !
    !  Find the intersection point in the "lower" plane - we need to know it
    !  to set up the interpolation and integration weights.
    !
    r = -kvec*gr%step(ind)/kvec(ind)
    !
    !  Choose the appropriate reordering sequence for the integration weights
    !
    select case (ind)
      case default
        stop 'multigrid%FieldBuildPhase - impossible ind'
      case (1) !  Integration along x, plane Y/Z
        r_order  = (/2,3,1/)
        Wu_order = (/3,2,1/)
      case (2) ! Integration along y, plane X/Z
        r_order  = (/1,3,2/)
        Wu_order = (/1,3,2/)
      case (3) ! Integration along Z, plane X/Y
        r_order  = (/1,2,3/)
        Wu_order = (/1,2,3/)
    end select
    incr = -nint(sign(1._rk,kvec(ind)))
    !
    !  Prepare interpolation and integration weights. Depending on the integration
    !  direction, the Wu (integration) matrix may need to be reordered
    !
    r1 = r(r_order)
    call GetMatrixWg(incr*r1/gr%step,Wg)
    call GetMatrixWu(incr*r1/gr%step,Wu)
    Wu = incr*Wu*gr%step(ind)/kvec(ind)
    if (incr<0) Wu = Wu(:,:,3:1:-1)  ! Swap orientation of the interpolation matrix
    Wu = reshape (Wu,(/3,3,3/),order=Wu_order)
    !
    ! write (out,"('Wg = '/3(1x,3f20.10/))") Wg
    ! write (out,"('Wu = '/3(3(1x,3f20.10/)/))") Wu
    !
    !  Branch into the appropriate integration loop - there are six cases in all
    !
    select case (3*ind+incr)
      case default
        stop 'multigrid%FieldBuildPhase - impossible integration direction'
      case (3*1-1) ! Integrate down along X
        kx_down: do ix=xu-2,xl+1,-1
          !$omp parallel do private(iy,iz)
          kx_down_z: do iz=zl+2,zu-2
            kx_down_y: do iy=yl+2,yu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(     ix+1,iy-1:iy+1,iz-1:iz+1,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do kx_down_y
          end do kx_down_z
          !$omp end parallel do
        end do kx_down
      case (3*1+1) ! Integrate up along X
        kx_up: do ix=xl+2,xu-1,1
          !$omp parallel do private(iy,iz)
          kx_up_z: do iz=zl+2,zu-2
            kx_up_y: do iy=yl+2,yu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(ix-1     ,iy-1:iy+1,iz-1:iz+1,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do kx_up_y
          end do kx_up_z
          !$omp end parallel do
        end do kx_up
      case (3*2-1) ! Integrate down along Y
        ky_down: do iy=yu-2,yl+1,-1
          !$omp parallel do private(ix,iz)
          ky_down_z: do iz=zl+2,zu-2
            ky_down_x: do ix=xl+2,xu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(ix-1:ix+1,     iy+1,iz-1:iz+1,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do ky_down_x
          end do ky_down_z
          !$omp end parallel do
        end do ky_down
      case (3*2+1) ! Integrate up along Y
        ky_up: do iy=yl+2,yu-1,1
          !$omp parallel do private(ix,iz)
          ky_up_z: do iz=zl+2,zu-2
            ky_up_x: do ix=xl+2,xu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(ix-1:ix+1,iy-1     ,iz-1:iz+1,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do ky_up_x
          end do ky_up_z
          !$omp end parallel do
        end do ky_up
      case (3*3-1)  ! Integrate down along Z
        kz_down: do iz=zu-2,zl+1,-1
          !$omp parallel do private(ix,iy)
          kz_down_y: do iy=yl+2,yu-2
            kz_down_x: do ix=xl+2,xu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(ix-1:ix+1,iy-1:iy+1,     iz+1,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do kz_down_x
          end do kz_down_y
          !$omp end parallel do
        end do kz_down
      case (3*3+1)  ! Integrate up along Z
        kz_up: do iz=zl+2,zu-1,1
          !$omp parallel do private(ix,iy)
          kz_up_y: do iy=yl+2,yu-2
            kz_up_x: do ix=xl+2,xu-2
              gr%fields(ix,iy,iz,phase) = sum(Wg * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1     ,phase),kind=rk)) &
                                        + sum(Wu * real(gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,pot  ),kind=rk))
            end do kz_up_x
          end do kz_up_y
          !$omp end parallel do
        end do kz_up
    end select
    !
    call FieldShrinkBox(pot)
    call FieldShrinkBox(phase)
    call TimerStop('FieldBuildPhase')
    !
  end subroutine FieldBuildPhase
  !
  !  Calculate the prefactor (2nd-order contribution) to the eikonal wavefunction
  !
  subroutine FieldBuildPref(phase,pref,kvec)
    integer(ik), intent(in)    :: phase     ! Field containing the phase; not modified
    integer(ik), intent(in)    :: pref      ! Field containing the prefactor; filled on output
    real(rk), intent(in)       :: kvec(3)   ! K-vector
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    real(rk)                   :: absk, gradw(3), gradp(3), pg, eps_pg
    !
    call TimerStart('FieldBuildPref')
    if (verbose>=2) then
      write (out,"(' = FieldBuildPref phase = ',i3,' pref = ',i3,' kvec = ',3g12.5)") phase, pref, kvec
    end if
    !
    if (grid%ngrids/=1) stop 'FieldBuildPref - multigrids not supported'
    if (phase==pref) stop 'FieldBuildPref - oops: input = output'
    !
    gr => grid%grids(1)
    xl  = max(gr%active(1,1,phase)-1,1) ; xu  = min(gr%active(2,1,phase)+1,gr%npoints(1))
    yl  = max(gr%active(1,2,phase)-1,1) ; yu  = min(gr%active(2,2,phase)+1,gr%npoints(2))
    zl  = max(gr%active(1,3,phase)-1,1) ; zu  = min(gr%active(2,3,phase)+1,gr%npoints(3))
    call FieldAugmentBox(pref, 1_ik,xl,  xu,  yl,  yu,  zl,  zu)
    call FieldAugmentBox(phase,1_ik,xl-1,xu+1,yl-1,yu+1,zl-1,zu+1)
    !
    absk   = sqrt(sum(kvec**2))
    gradw  = 0.5_rk / gr%step
    eps_pg = 4._rk * spacing(absk**2)
    !
    !$omp parallel do private(ix,iy,iz,gradp,pg)
    p_z: do iz=zl,zu
      p_y: do iy=yl,yu
        p_x: do ix=xl,xu
          gradp(1) = gradw(1)*real(gr%fields(ix+1,iy,iz,phase)-gr%fields(ix-1,iy,iz,phase),kind=rk)
          gradp(2) = gradw(2)*real(gr%fields(ix,iy+1,iz,phase)-gr%fields(ix,iy-1,iz,phase),kind=rk)
          gradp(3) = gradw(3)*real(gr%fields(ix,iy,iz+1,phase)-gr%fields(ix,iy,iz-1,phase),kind=rk)
          pg = sum( (kvec+gradp)**2 ) + eps_pg
          gr%fields(ix,iy,iz,pref) = sqrt( absk/sqrt(pg) )
        end do p_x
      end do p_y
    end do p_z
    call FieldShrinkBox(phase)
    call FieldShrinkBox(pref)
    call TimerStop('FieldBuildPref')
    !
  end subroutine FieldBuildPref
  !
  !  Return a 2D slab of the field
  !
  subroutine FieldFetchSlab(src,dir,n1,n2,data)
    integer(ik), intent(in)    :: src  ! Input: the data field
    integer(ik), intent(in)    :: dir  ! Input: Direction of the "thin" axis of the slab
    integer(ik), intent(in)    :: n1   ! Input: First layer to extract. Negative values count from the top
    integer(ik), intent(in)    :: n2   ! Input: Last layer to extract. Negative values count from the top
    complex(rk), intent(out)   :: data(:,:,:) ! Output: where to dump the data
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xu, yu, zu
    integer(ik)                :: cl, cu
    !
    call TimerStart('FieldFetchSlab')
    if (verbose>=2) then
      write (out,"(' = FieldFetchSlab src = ',i3,' dir = ',i2,' range = ',2i6)") src, dir, n1, n2
    end if
    !
    if (grid%ngrids/=1) stop 'FieldFetchSlab - multigrids not supported'
    if (dir<=0 .or. dir>3) stop 'FieldFetchSlab - bad dir argument'
    !
    gr => grid%grids(1)
    xu  = gr%npoints(1)
    yu  = gr%npoints(2)
    zu  = gr%npoints(3)
    call FieldAugmentBox(src,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    !
    !  Get slab dimensions
    !
    if (n1<0) then ; cl = gr%npoints(dir) + 1 + n1 ; else ; cl = n1 ; end if
    if (n2<0) then ; cu = gr%npoints(dir) + 1 + n2 ; else ; cu = n2 ; end if
    if (cl<1 .or. cl>gr%npoints(dir) .or. cu<1 .or. cu>gr%npoints(dir) ) then
      stop 'FieldFetchSlab - bad slab dimensions'
    end if
    !
    select case (dir)
      case default ; stop 'FieldFetchSlab - invalid slab direction'
      case (1)
        if (any(ubound(data)/=(/cu-cl+1,yu,zu/))) stop 'FieldFetchSlab - slab size mismatch'
        data = gr%fields(cl:cu,1:yu,1:zu,src)
      case (2)
        if (any(ubound(data)/=(/xu,cu-cl+1,zu/))) stop 'FieldFetchSlab - slab size mismatch'
        data = gr%fields(1:xu,cl:cu,1:zu,src)
      case (3)
        if (any(ubound(data)/=(/xu,yu,cu-cl+1/))) stop 'FieldFetchSlab - slab size mismatch'
        data = gr%fields(1:xu,1:yu,cl:cu,src)
    end select
    !
    call FieldShrinkBox(src)
    call TimerStop('FieldFetchSlab')
    !
  end subroutine FieldFetchSlab
  !
  !  Replace a 2D slab of the field
  !
  subroutine FieldSetSlab(dst,dir,n1,n2,data)
    integer(ik), intent(in)    :: dst  ! Input/Output: the data field
    integer(ik), intent(in)    :: dir  ! Input: Direction of the "thin" axis of the slab
    integer(ik), intent(in)    :: n1   ! Input: First layer to extract. Negative values count from the top
    integer(ik), intent(in)    :: n2   ! Input: Last layer to extract. Negative values count from the top
    complex(rk), intent(in)    :: data(:,:,:) ! Output: where to dump the data
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xu, yu, zu
    integer(ik)                :: cl, cu
    !
    call TimerStart('FieldSetSlab')
    if (verbose>=2) then
      write (out,"(' = FieldSetSlab dst = ',i3,' dir = ',i2,' range = ',2i6)") dst, dir, n1, n2
    end if
    !
    if (grid%ngrids/=1) stop 'FieldSetSlab - multigrids not supported'
    if (dir<=0 .or. dir>3) stop 'FieldSetSlab - bad dir argument'
    !
    gr => grid%grids(1)
    xu  = gr%npoints(1)
    yu  = gr%npoints(2)
    zu  = gr%npoints(3)
    call FieldAugmentBox(dst,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    !
    !  Get slab dimensions
    !
    if (n1<0) then ; cl = gr%npoints(dir) + 1 + n1 ; else ; cl = n1 ; end if
    if (n2<0) then ; cu = gr%npoints(dir) + 1 + n2 ; else ; cu = n2 ; end if
    if (cl<1 .or. cl>gr%npoints(dir) .or. cu<1 .or. cu>gr%npoints(dir) ) then
      stop 'FieldSetSlab - bad slab dimensions'
    end if
    !
    select case (dir)
      case default ; stop 'FieldSetSlab - invalid slab direction'
      case (1)
        if (any(ubound(data)/=(/cu-cl+1,yu,zu/))) stop 'FieldSetSlab - slab size mismatch'
        gr%fields(cl:cu,1:yu,1:zu,dst) = data
      case (2)
        if (any(ubound(data)/=(/xu,cu-cl+1,zu/))) stop 'FieldSetSlab - slab size mismatch'
        gr%fields(1:xu,cl:cu,1:zu,dst) = data
      case (3)
        if (any(ubound(data)/=(/xu,yu,cu-cl+1/))) stop 'FieldSetSlab - slab size mismatch'
        gr%fields(1:xu,1:yu,cl:cu,dst) = data
    end select
    !
    call FieldShrinkBox(dst)
    call TimerStop('FieldSetSlab')
    !
  end subroutine FieldSetSlab
  !
  !  Dummy entry, to make sure cylindrical modules can be linked against this module
  !
  subroutine FieldFetchSlice2D(src,dir,n,data)
    integer(ik), intent(in)    :: src     ! Input: the data field
    integer(ik), intent(in)    :: dir     ! Input: Direction of the "thin" axis of the slab
    integer(ik), intent(in)    :: n       ! Input: Layer to extract. Negative values count from the top
    complex(rk), intent(out)   :: data(:) ! Output: where to dump the data
    !
    stop 'multigrid%FieldFetchSlice2D - nonsensical entry point called'
  end subroutine FieldFetchSlice2D
  !
  !  Return a 2D slice of the field
  !
  subroutine FieldFetchSlice3D(src,dir,n,data)
    integer(ik), intent(in)    :: src       ! Input: the data field
    integer(ik), intent(in)    :: dir       ! Input: Direction of the "thin" axis of the slab
    integer(ik), intent(in)    :: n         ! Input: First layer to extract. Negative values count from the top
    complex(rk), intent(out)   :: data(:,:) ! Output: where to dump the data
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xu, yu, zu
    integer(ik)                :: i
    integer(ik)                :: cl
    !
    call TimerStart('FieldFetchSlice')
    if (verbose>=2) then
      write (out,"(' = FieldFetchSlice src = ',i3,' dir = ',i2,' slice = ',i6)") src, dir, n
    end if
    !
    if (grid%ngrids/=1) stop 'FieldFetchSlice - multigrids not supported'
    if (dir<=0 .or. dir>3) stop 'FieldFetchSlice - bad dir argument'
    !
    gr => grid%grids(1)
    xu  = gr%npoints(1)
    yu  = gr%npoints(2)
    zu  = gr%npoints(3)
    !
    !  Get slab dimensions
    !
    if (n<0) then ; cl = gr%npoints(dir) + 1 + n ; else ; cl = n ; end if
    if (cl<1 .or. cl>gr%npoints(dir)) then
      stop 'FieldFetchSlice - bad slice dimensions'
    end if
    !
    select case (dir)
      case default ; stop 'FieldFetchSlice - invalid slice direction'
      case (1)
        if (any(ubound(data)/=(/yu,zu/))) stop 'FieldFetchSlice - buffer size mismatch'
        call FieldAugmentBox(src,1_ik,    cl,cl,  1_ik,yu,  1_ik,zu)
        !$omp parallel do private(i)
        do i=1,zu
          data(1:yu,i) = gr%fields(cl,1:yu,i,src)
        end do
      case (2)
        if (any(ubound(data)/=(/xu,zu/))) stop 'FieldFetchSlice - buffer size mismatch'
        call FieldAugmentBox(src,1_ik,  1_ik,xu,    cl,cl,  1_ik,zu)
        !$omp parallel do private(i)
        do i=1,zu
          data(1:xu,i) = gr%fields(1:xu,cl,i,src)
        end do
      case (3)
        if (any(ubound(data)/=(/xu,yu/))) stop 'FieldFetchSlice - buffer size mismatch'
        call FieldAugmentBox(src,1_ik,  1_ik,xu,  1_ik,yu,    cl,cl)
        !$omp parallel do private(i)
        do i=1,yu
          data(1:xu,i) = gr%fields(1:xu,i,cl,src)
        end do
    end select
    !
    call FieldShrinkBox(src)
    call TimerStop('FieldFetchSlice')
    !
  end subroutine FieldFetchSlice3D
  !
  !  Update a field by an inner product of two vectors.
  !  This is a dummy entry for compatibility with cylindrical code
  !
  subroutine FieldGER2D(dst,a,b,limits)
    integer(ik), intent(in)           :: dst         ! Field to be updated
    complex(rk), intent(in)           :: a(:)        ! First-index factor
    complex(rk), intent(in)           :: b(:)        ! Second-index factor
    integer(ik), intent(in), optional :: limits(:,:) ! Optional grid limits
    !
    stop 'multigrid%FieldGER2D - this call is nonsensical in 3D'
  end subroutine FieldGER2D
  !
  !  Two versions of GER operation are possible in 3D: update from a 2D + 1D fields,
  !  and an update from three 1D fields. The 2+1 case is a little more complicated
  !  (it has three distinct versions) and goes first.
  !
  subroutine FieldGER3D21(dst,dir,a,b,limits)
    integer(ik), intent(in)           :: dst         ! Field to be updated
    integer(ik), intent(in)           :: dir         ! Axis for 1-index component
    complex(rk), intent(in)           :: a(:,:)      ! Two-index factor
    complex(rk), intent(in)           :: b(:)        ! One-index factor
    integer(ik), intent(in), optional :: limits(:,:) ! Optional grid limits
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: xl, yl, zl
    integer(ik)                :: xu, yu, zu
    !
    call TimerStart('FieldGER3D21')
    if (verbose>=2) then
      write (out,"(' = FieldGER3D21 dst = ',i3,' dir = ',i2)") dst, dir
    end if
    !
    if (grid%ngrids/=1) stop 'FieldGER3D21 - multigrids not supported'
    if (dir<=0 .or. dir>3) stop 'FieldGER3D21 - bad dir argument'
    !
    gr => grid%grids(1)
    xl = 1 ; xu = gr%npoints(1)
    yl = 1 ; yu = gr%npoints(2)
    zl = 1 ; zu = gr%npoints(3)
    !
    !  Input arrays a and b must match full dimensions of the data fields; check now,
    !  before we clip the ranges
    !
    select case (dir)
      case (1) ; if (size(b)/=xu .or. any(ubound(a)/=(/yu,zu/))) stop 'multigrid%FieldGER3D21 - dimension mismatch (1)'
      case (2) ; if (size(b)/=yu .or. any(ubound(a)/=(/xu,zu/))) stop 'multigrid%FieldGER3D21 - dimension mismatch (2)'
      case (3) ; if (size(b)/=zu .or. any(ubound(a)/=(/xu,yu/))) stop 'multigrid%FieldGER3D21 - dimension mismatch (3)'
    end select
    !
    if (present(limits)) then
      xl = max(xl,limits(1,1))
      xu = min(xu,limits(2,1))
      yl = max(yl,limits(1,2))
      yu = min(yu,limits(2,2))
      zl = max(zl,limits(1,3))
      zu = min(zu,limits(2,3))
    end if
    call FieldAugmentBox(dst,1_ik,  xl,xu,  yl,yu,  zl,zu)
    !
    !  The rest of the code is a bunch of special cases
    !
    select case (dir)
      case default ; stop 'FieldGER3D21 - invalid slice direction'
      case (1)
        !$omp parallel do private(ix,iy,iz)
        fill_1z: do iz=zl,zu
          fill_1y: do iy=yl,yu
            fill_1x: do ix=xl,xu
              gr%fields(ix,iy,iz,dst) = gr%fields(ix,iy,iz,dst) + a(iy,iz) * b(ix)
            end do fill_1x
          end do fill_1y
        end do fill_1z 
      case (2)
        !$omp parallel do private(ix,iy,iz)
        fill_2z: do iz=zl,zu
          fill_2y: do iy=yl,yu
            fill_2x: do ix=xl,xu
              gr%fields(ix,iy,iz,dst) = gr%fields(ix,iy,iz,dst) + a(ix,iz) * b(iy)
            end do fill_2x
          end do fill_2y
        end do fill_2z 
      case (3)
        !$omp parallel do private(ix,iy,iz)
        fill_3z: do iz=zl,zu
          fill_3y: do iy=yl,yu
            fill_3x: do ix=xl,xu
              gr%fields(ix,iy,iz,dst) = gr%fields(ix,iy,iz,dst) + a(ix,iy) * b(iz)
            end do fill_3x
          end do fill_3y
        end do fill_3z 
    end select
    !
    call FieldShrinkBox(dst)
    call TimerStop('FieldGER3D21')
    !
  end subroutine FieldGER3D21
  !
  subroutine FieldGER3D111(dst,a,b,c,limits)
    integer(ik), intent(in)           :: dst         ! Field to be updated
    complex(rk), intent(in)           :: a(:)        ! Factor for the X axis
    complex(rk), intent(in)           :: b(:)        ! Factor for the Y axis
    complex(rk), intent(in)           :: c(:)        ! Factor for the Z axis
    integer(ik), intent(in), optional :: limits(:,:) ! Optional grid limits
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: xl, yl, zl
    integer(ik)                :: xu, yu, zu
    !
    call TimerStart('FieldGER3D111')
    if (verbose>=2) then
      write (out,"(' = FieldGER3D111 dst = ',i3)") dst
    end if
    !
    if (grid%ngrids/=1) stop 'FieldGER3D111 - multigrids not supported'
    !
    gr => grid%grids(1)
    xl = 1 ; xu = gr%npoints(1)
    yl = 1 ; yu = gr%npoints(2)
    zl = 1 ; zu = gr%npoints(3)
    !
    !  Input arrays a, b, and c must match full dimensions of the data fields; check now,
    !  before we clip the ranges
    !
    if (size(a)/=xu .or. size(b)/=yu .or. size(c)/=zu) stop 'multigrid%FieldGER3D111 - dimension mismatch'
    !
    if (present(limits)) then
      xl = max(xl,limits(1,1))
      xu = min(xu,limits(2,1))
      yl = max(yl,limits(1,2))
      yu = min(yu,limits(2,2))
      zl = max(zl,limits(1,3))
      zu = min(zu,limits(2,3))
    end if
    call FieldAugmentBox(dst,1_ik,  xl,xu,  yl,yu,  zl,zu)
    !
    !$omp parallel do private(ix,iy,iz)
    fill_1z: do iz=zl,zu
      fill_1y: do iy=yl,yu
        fill_1x: do ix=xl,xu
          gr%fields(ix,iy,iz,dst) = gr%fields(ix,iy,iz,dst) + a(ix) * b(iy) * c(iz)
        end do fill_1x
      end do fill_1y
    end do fill_1z 
    !
    call FieldShrinkBox(dst)
    call TimerStop('FieldGER3D111')
    !
  end subroutine FieldGER3D111
  !
  !  Perform one time integration step for the eikonal-Volkov continuum solution.
  !  The line integrals are taken over the potential supplied in src_u.
  !  This routine is meant for the interior points only.
  !  The convolution matrices have been calculated elsewhere.
  !  The two layers of the boundary points are not updated.
  !  Update of the inner layer of the boundary is implemented in FieldEVAExteriorStep.
  !  The outer boundary layer is set to zero.
  !
  subroutine FieldEVAInteriorStep(src_u,src_g,dst_g,wgt_u,wgt_g)
    integer(ik), intent(in)    :: src_u     ! Potential; not modified
    integer(ik), intent(in)    :: src_g     ! Initial value of the phase; not modifield
    integer(ik), intent(in)    :: dst_g     ! Final value of the phase; interior points overwritten
    real(rk), intent(in)       :: wgt_u(-1:1,-1:1,-1:1) ! Stencil for the line integral of u
    real(rk), intent(in)       :: wgt_g(-1:1,-1:1,-1:1) ! Stencil for the initial-point terpolation in g
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: xu, yu, zu
    !
    call TimerStart('FieldEVAInteriorStep')
    if (verbose>=2) then
      write (out,"(' = FieldEVAInteriorStep src_u = ',i3,' src_g = ',i3,' dst_g = ',i3)") src_u, src_g, dst_g
    end if
    !
    if (grid%ngrids/=1) stop 'FieldEVAInteriorStep - multigrids not supported'
    if (src_u==src_g .or. src_u==dst_g .or. src_g==dst_g) stop 'FieldEVAInteriorStep - overlapping fields'
    !
    gr => grid%grids(1)
    xu  = gr%npoints(1)
    yu  = gr%npoints(2)
    zu  = gr%npoints(3)
    call FieldAugmentBox(src_u,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    call FieldAugmentBox(src_g,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    call FieldAugmentBox(dst_g,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    !
    !$omp parallel do private(ix,iy,iz)
    p_z: do iz=3,zu-2
      p_y: do iy=3,yu-2
        p_x: do ix=3,xu-2
          gr%fields(ix,iy,iz,dst_g) = sum(wgt_u*gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,src_u)) + &
                                      sum(wgt_g*gr%fields(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,src_g))
        end do p_x
      end do p_y
    end do p_z
    call FieldShrinkBox(src_u)
    call FieldShrinkBox(src_g)
    call FieldShrinkBox(dst_g)
    call TimerStop('FieldEVAInteriorStep')
    !
  end subroutine FieldEVAInteriorStep
  !
  !  Update the inner boundary layer of the eikonal-Volkov continuum solution.
  !  For each edge point, two possibilities exist. It is either 
  !   (a) an "inside" edge, with the trajectory arriving from the interior of the grid, or
  !   (b) an "outside" edge, with the trajectory coming into the grid interior.
  !  The "inside" edges can be handled through one of a finite number of modified
  !  interior cases. The "outside" edges are done elsewhere
  !
  !  Overall, eight possibilities exist for slassifying a point as "inside" or
  !  "outside". These are encoded by the bits of the topology argument:
  !
  !                 0              1
  !   LSB+0:     ix==2    ix==gr%npoints(1)-1   
  !   LSB+1:     iy==2    iy==gr%npoints(2)-1   
  !   LSB+2:     iz==2    iz==gr%npoints(3)-1  
  !
  subroutine FieldEVAExteriorStep(src_u,src_g,dst_g,topology,wgt_u,wgt_g)
    integer(ik), intent(in)    :: src_u     ! Potential; not modified
    integer(ik), intent(in)    :: src_g     ! Initial value of the phase; not modifield
    integer(ik), intent(in)    :: dst_g     ! Final value of the phase; inner boundary overwritten
    integer(ik), intent(in)    :: topology  ! Topology of the trajectory - determines which points
                                            ! are "inside", and which are "outside" (see above).
    real(rk), intent(in)       :: wgt_u(-1:1,-1:1,-1:1)                ! Stencil for the line integral of u
                                     !   dx   dy   dz   bx   by   bz   ! The last three indices enumerate possible cases
    real(rk), intent(in)       :: wgt_g(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1) ! Stencil for the initial-point terpolation in g
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xu, yu, zu  ! Upper bound for grid points
    integer(ik)                :: gx, gy, gz  ! Position of the slice of the phase array
    integer(ik)                :: ux, uy, uz  ! Position of the slice of the potential array
    integer(ik)                :: ix, iy, iz  ! Counters
    !
    call TimerStart('FieldEVAExteriorStep')
    if (verbose>=2) then
      write (out,"(' = FieldEVAExteriorStep src_u = ',i3,' src_g = ',i3,' dst_g = ',i3)") src_u, src_g, dst_g
    end if
    !
    if (grid%ngrids/=1) stop 'FieldEVAExteriorStep - multigrids not supported'
    if (src_u==src_g .or. src_u==dst_g .or. src_g==dst_g) stop 'FieldEVAExteriorStep - overlapping fields'
    if (topology<0 .or. topology>7) stop 'FieldEVAExteriorStep - impossible topology'
    !
    gr => grid%grids(1)
    xu  = gr%npoints(1)
    yu  = gr%npoints(2)
    zu  = gr%npoints(3)
    call FieldAugmentBox(src_u,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    call FieldAugmentBox(src_g,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    call FieldAugmentBox(dst_g,1_ik,  1_ik,xu,  1_ik,yu,  1_ik,zu)
    !
    !  Handle the "inside" points first. We try to eliminate all 
    !  conditional expressions from the inner loops. Hopefully, this
    !  should result in better vectorization.
    !
    select case (topology)
      case (0+0+0) ; gx = 3    ; gy = 3    ; gz = 3    ; ux = 2    ; uy = 2    ; uz = 2   
      case (0+0+1) ; gx = xu-2 ; gy = 3    ; gz = 3    ; ux = xu-1 ; uy = 2    ; uz = 2   
      case (0+2+0) ; gx = 3    ; gy = yu-2 ; gz = 3    ; ux = 2    ; uy = yu-1 ; uz = 2   
      case (0+2+1) ; gx = xu-2 ; gy = yu-2 ; gz = 3    ; ux = xu-1 ; uy = yu-1 ; uz = 2   
      case (4+0+0) ; gx = 3    ; gy = 3    ; gz = zu-2 ; ux = 2    ; uy = 2    ; uz = zu-1
      case (4+0+1) ; gx = xu-2 ; gy = 3    ; gz = zu-2 ; ux = xu-1 ; uy = 2    ; uz = zu-1
      case (4+2+0) ; gx = 3    ; gy = yu-2 ; gz = zu-2 ; ux = 2    ; uy = yu-1 ; uz = zu-1
      case (4+2+1) ; gx = xu-2 ; gy = yu-2 ; gz = zu-2 ; ux = xu-1 ; uy = yu-1 ; uz = zu-1
    end select
    call inside_plane_yz(wgt_g(:,:,:,ux-gx,    0,    0))
    call inside_plane_xz(wgt_g(:,:,:,    0,uy-gy,    0))
    call inside_plane_xy(wgt_g(:,:,:,    0,    0,uz-gz))
    call inside_edge_x  (wgt_g(:,:,:,    0,uy-gy,uz-gz))
    call inside_edge_y  (wgt_g(:,:,:,ux-gx,    0,uz-gz))
    call inside_edge_z  (wgt_g(:,:,:,ux-gx,uy-gy,    0))
    call inside_corner  (wgt_g(:,:,:,ux-gx,uy-gy,uz-gz))
    !
    call FieldShrinkBox(src_u)
    call FieldShrinkBox(src_g)
    call FieldShrinkBox(dst_g)
    call TimerStop('FieldEVAExteriorStep')
    !
    contains
    !
    !  Handle side planes. There are three versions, one for each face of the cube.
    !
    subroutine inside_plane_yz(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: YZ plane, UX = ',i6,' GX = ',i6)") ux, gx
      !$omp parallel do private(iy,iz)
      p_z: do iz=3,zu-2
        p_y: do iy=3,yu-2
          gr%fields(ux,iy,iz,dst_g) = sum(wgt_u*gr%fields(ux-1:ux+1,iy-1:iy+1,iz-1:iz+1,src_u)) + &
                                      sum(wgt_g*gr%fields(gx-1:gx+1,iy-1:iy+1,iz-1:iz+1,src_g))
        end do p_y
      end do p_z
      !$omp end parallel do
    end subroutine inside_plane_yz
    !
    subroutine inside_plane_xz(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: XZ plane, UY = ',i6,' GY = ',i6)") uy, gy
      !$omp parallel do private(ix,iz)
      p_z: do iz=3,zu-2
        p_x: do ix=3,xu-2
          gr%fields(ix,uy,iz,dst_g) = sum(wgt_u*gr%fields(ix-1:ix+1,uy-1:uy+1,iz-1:iz+1,src_u)) + &
                                      sum(wgt_g*gr%fields(ix-1:ix+1,gy-1:gy+1,iz-1:iz+1,src_g))
        end do p_x
      end do p_z
      !$omp end parallel do
    end subroutine inside_plane_xz
    !
    subroutine inside_plane_xy(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: XY plane, UZ = ',i6,' GZ = ',i6)") uz, gz
      !$omp parallel do private(ix,iy)
      p_y: do iy=3,yu-2
        p_x: do ix=3,xu-2
          gr%fields(ix,iy,uz,dst_g) = sum(wgt_u*gr%fields(ix-1:ix+1,iy-1:iy+1,uz-1:uz+1,src_u)) + &
                                      sum(wgt_g*gr%fields(ix-1:ix+1,iy-1:iy+1,gz-1:gz+1,src_g))
        end do p_x
      end do p_y
      !$omp end parallel do
    end subroutine inside_plane_xy
    !
    !  Halndle grid edges. Again, three versions - one for each direction
    !
    subroutine inside_edge_x(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: X edge, UY/UZ = ',2i6,' GY/GZ = ',2i6)") uy, uz, gy, gz
      !$omp parallel do private(ix)
      p_x: do ix=3,xu-2
        gr%fields(ix,uy,uz,dst_g) = sum(wgt_u*gr%fields(ix-1:ix+1,uy-1:uy+1,uz-1:uz+1,src_u)) + &
                                    sum(wgt_g*gr%fields(ix-1:ix+1,gy-1:gy+1,gz-1:gz+1,src_g))
      end do p_x
      !$omp end parallel do
    end subroutine inside_edge_x
    !
    subroutine inside_edge_y(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: Y edge, UX/UZ = ',2i6,' GX/GZ = ',2i6)") ux, uz, gx, gz
      !$omp parallel do private(iy)
      p_y: do iy=3,yu-2
        gr%fields(ux,iy,uz,dst_g) = sum(wgt_u*gr%fields(ux-1:ux+1,iy-1:iy+1,uz-1:uz+1,src_u)) + &
                                    sum(wgt_g*gr%fields(gx-1:gx+1,iy-1:iy+1,gz-1:gz+1,src_g))
      end do p_y
      !$omp end parallel do
    end subroutine inside_edge_y
    !
    subroutine inside_edge_z(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: Z edge, UX/UY = ',2i6,' GX/GY = ',2i6)") ux, uy, gx, gy
      !$omp parallel do private(iz)
      p_z: do iz=3,zu-2
        gr%fields(ux,uy,iz,dst_g) = sum(wgt_u*gr%fields(ux-1:ux+1,uy-1:uy+1,iz-1:iz+1,src_u)) + &
                                    sum(wgt_g*gr%fields(gx-1:gx+1,gy-1:gy+1,iz-1:iz+1,src_g))
      end do p_z
      !$omp end parallel do
    end subroutine inside_edge_z
    !
    !  Finally, handle the grid corner
    !
    subroutine inside_corner(wgt_g)
      real(rk), intent(in)    :: wgt_g(:,:,:) ! Interpolation stencil for the initial phase
      ! write (out,"('FieldEVAExteriorStep: Corner, UX/UY/UZ = ',3i6,' GX/GY/GZ = ',3i6)") ux, uy, uz, gx, gy, gz
      gr%fields(ux,uy,uz,dst_g) = sum(wgt_u*gr%fields(ux-1:ux+1,uy-1:uy+1,uz-1:uz+1,src_u)) + &
                                  sum(wgt_g*gr%fields(gx-1:gx+1,gy-1:gy+1,gz-1:gz+1,src_g))
    end subroutine inside_corner
  end subroutine FieldEVAExteriorStep
  !
  !  Calculates the integral:
  !    < psi | v | psi* >
  !  where v is a real scalar field, defined by an external
  !  function.
  !
  !  The function must take one (vector) argument - coordinates
  !  of the point where it is evaluated.
  !
  function FieldScalarIntegrate(src,func) result(val)
    integer(ik), intent(in)     :: src   ! Field containing the wavefunction
    real(rk), external          :: func  ! Function
    real(rk)                    :: val   ! expectational value of the scalar
                                         ! field func
    integer(ik)                 :: igrid
    type(SimpleGridT), pointer  :: gr
    integer(ik)                 :: ix, iy, iz
    integer(ik)                 :: zl, zu
    integer(ik)                 :: xl, xu, yl, yu

    !
    call TimerStart('FieldScalarIntegrate')
    if (verbose>=2) write (out,"(' FieldScalarIntegrate ',i3)") src
    !
    val = 0.0_rk
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      zl = gr%active(1,3,src) ; zu = gr%active(2,3,src)
      xl = gr%active(1,1,src) ; xu = gr%active(2,1,src)
      yl = gr%active(1,2,src) ; yu = gr%active(2,2,src)
      !
      !  SGI and Pathscale compilers requires "func" to be declared as "shared"
      !
!     !$omp parallel do reduction(+:val) private(ix,iy,iz) shared(func)
      !$omp parallel do reduction(+:val) private(ix,iy,iz)
      do iz=zl,zu
        do iy=yl,yu
          do ix=xl,xu
             val = val +  func(gr%coords(1:3,ix,iy,iz))  &
                        * gr%coords(4,ix,iy,iz)          &
                        * real(gr%fields(ix,iy,iz,src) * conjg(gr%fields(ix,iy,iz,src)),kind=rk)
          end do
        end do
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldScalarIntegrate')
  end function FieldScalarIntegrate
  !
  !  Calculates the integral:
  !    < psi | v | psi* >
  !  where v is a real scalar field, defined by an external
  !  function.
  !
  !  The function must take one (vector) argument - coordinates
  !  of the point where it is evaluated.
  !
  function FieldBraVKet(bra,v,ket) result(val)
    integer(ik), intent(in)     :: bra   ! "Left" wavefunction field 
    integer(ik), intent(in)     :: v     ! Multiplicative operator field
    integer(ik), intent(in)     :: ket   ! "Right" wavefunction field 
    complex(rk)                 :: val   ! Integrated matrix element
    integer(ik)                 :: igrid
    type(SimpleGridT), pointer  :: gr
    integer(ik)                 :: iz, zl, zu
    integer(ik)                 :: xl, xu, yl, yu

    !
    call TimerStart('FieldBraVKet')
    if (verbose>=2) write (out,"(' FieldBraVKet',3i4)") bra,v,ket
    !
    val = 0
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      xl = maxval(gr%active(1,1,(/bra,v,ket/)))
      xu = minval(gr%active(2,1,(/bra,v,ket/)))
      yl = maxval(gr%active(1,2,(/bra,v,ket/)))
      yu = minval(gr%active(2,2,(/bra,v,ket/)))
      zl = maxval(gr%active(1,3,(/bra,v,ket/)))
      zu = minval(gr%active(2,3,(/bra,v,ket/)))
      !
      !$omp parallel do reduction(+:val) private(iz)
      do iz=zl,zu
        val = val + sum(gr%coords(4,xl:xu,yl:yu,iz)       &
                * conjg(gr%fields(  xl:xu,yl:yu,iz,bra))  &
                *       gr%fields(  xl:xu,yl:yu,iz,v  )   &
                *       gr%fields(  xl:xu,yl:yu,iz,ket))
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldBraVKet')
  end function FieldBraVKet
  !
  !  Initialize the field
  !
  !DEC$ATTRIBUTES FORCEINLINE::FieldInit
  subroutine FieldInit(dst,func,mask,grace)
    integer(ik), intent(in)           :: dst
    complex(rk), external             :: func
    integer(ik), intent(in), optional :: mask  ! Field boundaries to use as a template
    integer(ik), intent(in), optional :: grace ! Extra border points to add to the template
    !DEC$ATTRIBUTES FORCEINLINE::func

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: izl, izu, xl, xu, yl, yu
    integer(ik)                :: plus

    !
    call TimerStart('FieldInit')
    if (verbose>=2) write (out,"(' = FieldInit ',i3)") dst
    !
    plus = 0
    if (present(grace)) plus = grace
    !
    !  This would have been a perfect match to an ELEMENTAL function.
    !  Unfortunately, there does not seem to be a way of passing an
    !  ELEMENTAL as an argument?
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  If mask was not specified, initialize the whole field
      !  Otherwise, initialize just the area of interest
      !
      gr%active(1,:,dst) = spread(1,1,3)
      gr%active(2,:,dst) = gr%npoints
      if (present(mask)) then
        gr%active(1,:,dst) = max(gr%active(1,:,mask)-plus,gr%active(1,:,dst))
        gr%active(2,:,dst) = min(gr%active(2,:,mask)+plus,gr%active(2,:,dst))
      end if
      !
      izl = gr%active(1,3,dst) ; izu = gr%active(2,3,dst)
      xl  = gr%active(1,1,dst) ; xu  = gr%active(2,1,dst)
      yl  = gr%active(1,2,dst) ; yu  = gr%active(2,2,dst)
      !
      !  SGI and Pathscale compilers require "func" to be declared as "shared"
      !
!     !$omp parallel do private(ix,iy,iz) shared(func)
      !$omp parallel do private(ix,iy,iz)
      do iz=izl,izu
        do iy=yl,yu
          do ix=xl,xu
            gr%fields(ix,iy,iz,dst) = func(gr%coords(1:3,ix,iy,iz))
          end do
        end do
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(dst)
    call TimerStop('FieldInit')
    !
  end subroutine FieldInit
  !
  !  Pass through the field through a function filter. The function must take
  !  two arguments - coordinate and the field value.
  !
  subroutine FieldProcess(dst,func)
    integer(ik), intent(in)  :: dst
    complex(rk), external    :: func
    !
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: zl, zu, xl, xu, yl, yu
    !
    call TimerStart('FieldProcess')
    if (verbose>=2) write (out,"(' = FieldProcess ',i3)") dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      zl = gr%active(1,3,dst) ; zu = gr%active(2,3,dst)
      xl = gr%active(1,1,dst) ; xu = gr%active(2,1,dst)
      yl = gr%active(1,2,dst) ; yu = gr%active(2,2,dst)
      !
      !  SGI and Pathscale compilers require "func" to be declared as "shared"
      !
!     !$omp parallel do private(ix,iy,iz) shared(func)
      !$omp parallel do private(ix,iy,iz)
      do iz=zl,zu
        do iy=yl,yu
          do ix=xl,xu
            gr%fields(ix,iy,iz,dst) = func(gr%coords(1:3,ix,iy,iz), &
                                           gr%fields(ix,iy,iz,dst))
          end do
        end do
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(dst)
    call TimerStop('FieldProcess')
    !
  end subroutine FieldProcess
  !
  !  Rotates the components of the vetor field (dstx,dsty,dstz)
  !
  !  [dstx']   [                      ] [dstx]
  !  [dsty'] = [  3x3 rotation matrix ] [dsty]
  !  [dstz']   [                      ] [dstz]
  !
  !  (Added March 2010, Michael Spanner)
  !
  subroutine FieldRotateVectorComponents(rotmat,dstx,dsty,dstz)
    real(ark), intent(in)             :: rotmat(3,3)
    integer(ik), intent(in)           :: dstx, dsty, dstz
    !
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: zl, zu, xl, xu, yl, yu
    !
    call TimerStart('FieldRotateVectorComponents')
    if (verbose>=2) write (out,"(' FieldRotateVectorComponents ',3(1x,i3))") dstx, dsty, dstz
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      xl = minval(gr%active(1,1,(/dstx,dsty,dstz/)))
      xu = maxval(gr%active(2,1,(/dstx,dsty,dstz/)))
      yl = minval(gr%active(1,2,(/dstx,dsty,dstz/)))
      yu = maxval(gr%active(2,2,(/dstx,dsty,dstz/)))
      zl = minval(gr%active(1,3,(/dstx,dsty,dstz/)))
      zu = maxval(gr%active(2,3,(/dstx,dsty,dstz/)))
      !
      call FieldAugmentBox(dstx,igrid,xl,xu,yl,yu,zl,zu)
      call FieldAugmentBox(dsty,igrid,xl,xu,yl,yu,zl,zu)
      call FieldAugmentBox(dstz,igrid,xl,xu,yl,yu,zl,zu)
      !
      !$omp parallel do private(ix,iy,iz)
      do iz=zl,zu
        do iy=yl,yu
          do ix=xl,xu
            gr%fields(ix,iy,iz,(/dstx,dsty,dstz/)) = matmul(rotmat,gr%fields(ix,iy,iz,(/dstx,dsty,dstz/)))
          end do
        end do
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(dstx)
    call FieldShrinkBox(dsty)
    call FieldShrinkBox(dstz)
    call TimerStop('FieldRotateVectorComponents')
    !
  end subroutine FieldRotateVectorComponents
  !
  !  Initialize field to zero
  !
  subroutine FieldZero(dst,limits)
    integer(ik), intent(in)           :: dst
    integer(ik), intent(in), optional :: limits(:,:)  ! Optional range of grid points to update
                                                      ! Only applies to the outermost grid

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: ixl, ixu, iyl, iyu

    !
    call TimerStart('FieldZero')
    if (verbose>=2) write (out,"(' = FieldZero ',i3)") dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      ixl = lbound(gr%fields,dim=1)
      ixu = ubound(gr%fields,dim=1)
      iyl = lbound(gr%fields,dim=2)
      iyu = ubound(gr%fields,dim=2)
      izl = lbound(gr%fields,dim=3)
      izu = ubound(gr%fields,dim=3)
      if (present(limits) .and. igrid==1) then
        ixl = max(ixl,limits(1,1))
        ixu = min(ixu,limits(2,1))
        iyl = max(iyl,limits(1,2))
        iyu = min(iyu,limits(2,2))
        izl = max(izl,limits(1,3))
        izu = min(izu,limits(2,3))
      end if
      gr%active(2:1:-1,:,dst) = gr%active(1:2,:,0)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(ixl:ixu,iyl:iyu,iz,dst) = 0
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldZero')
    !
  end subroutine FieldZero

  subroutine FieldCopy(src,dst)
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    !
    call TimerStart('FieldCopy')
    if (verbose>=2) write (out,"(' = FieldCopy ',i3,' to ',i3)") src, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      gr%active(:,:,dst) = gr%active(:,:,src)
      izl = gr%active(1,3,src) ; izu = gr%active(2,3,src)
      xl  = gr%active(1,1,src) ; xu  = gr%active(2,1,src)
      yl  = gr%active(1,2,src) ; yu  = gr%active(2,2,src)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,src)
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldCopy')
    !
  end subroutine FieldCopy

  subroutine FieldConjugate(src,dst)
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    !
    call TimerStart('FieldConjugate')
    if (verbose>=2) write (out,"(' = FieldConjugate ',i3,' to ',i3)") src, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      gr%active(:,:,dst) = gr%active(:,:,src)
      izl = gr%active(1,3,src) ; izu = gr%active(2,3,src)
      xl  = gr%active(1,1,src) ; xu  = gr%active(2,1,src)
      yl  = gr%active(1,2,src) ; yu  = gr%active(2,2,src)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = conjg(gr%fields(xl:xu,yl:yu,iz,src))
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldConjugate')
    !
  end subroutine FieldConjugate

  subroutine FieldScale(dst,con)
    integer(ik), intent(in) :: dst
    complex(rk), intent(in) :: con

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    !
    call TimerStart('FieldScale')
    if (verbose>=2) write (out,"(' = FieldScale ',i3)") dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = gr%active(1,3,dst) ; izu = gr%active(2,3,dst)
      xl  = gr%active(1,1,dst) ; xu  = gr%active(2,1,dst)
      yl  = gr%active(1,2,dst) ; yu  = gr%active(2,2,dst)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = con * gr%fields(xl:xu,yl:yu,iz,dst)
      end do
      !$omp end parallel do
    end do
    call TimerStop('FieldScale')
    !
  end subroutine FieldScale

  !
  !  Removes unsafely large and special IEEE (nan/inv) values from the grid
  !
  subroutine FieldSanitize(dst)
    integer(ik), intent(in) :: dst  ! Field to sanitize

!   integer(ik)                :: igrid
!   type(SimpleGridT), pointer :: gr
!   integer(ik)                :: iz, izl, izu
!   integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldSanitize')
    if (verbose>=2) write (out,"(' = FieldSanitize ',i3)") dst
    !
!   do igrid=1,grid%ngrids
!     gr => grid%grids(igrid)
!     izl = gr%active(1,3,dst) ; izu = gr%active(2,3,dst)
!     xl  = gr%active(1,1,dst) ; xu  = gr%active(2,1,dst)
!     yl  = gr%active(1,2,dst) ; yu  = gr%active(2,2,dst)
!     !$omp parallel do private(iz)
!     do iz=izl,izu
!       where (isnan(abs(gr%fields(xl:xu,yl:yu,iz,dst))))
!         gr%fields(xl:xu,yl:yu,iz,dst) = 0
!       end where
!       where (abs(gr%fields(xl:xu,yl:yu,iz,dst))>safe_max)
!         gr%fields(xl:xu,yl:yu,iz,dst) = safe_max
!       end where
!     end do
!     !$omp end parallel do
!   end do
    call TimerStop('FieldSanitize')
    !
  end subroutine FieldSanitize

  subroutine FieldSum(src,dst)
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldSum')
    if (verbose>=2) write (out,"(' = FieldSum ',i3,' + ',i3,' -> ',i3)") src,dst,dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = minval(gr%active(1,3,(/src,dst/)))
      izu = maxval(gr%active(2,3,(/src,dst/)))
      xl  = minval(gr%active(1,1,(/src,dst/)))
      xu  = maxval(gr%active(2,1,(/src,dst/)))
      yl  = minval(gr%active(1,2,(/src,dst/)))
      yu  = maxval(gr%active(2,2,(/src,dst/)))
      call FieldAugmentBox(src,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(dst,igrid,xl,xu,yl,yu,izl,izu)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,dst) &
                                      + gr%fields(xl:xu,yl:yu,iz,src)
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldSum')
    !
  end subroutine FieldSum
  !
  !  |dst> = |dst> * |src_a>
  !
  subroutine FieldMul(src,dst)
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldMul')
    if (verbose>=2) write (out,"(' = FieldMul ',i3,' * ',i3,' -> ',i3)") dst, src, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  The product field is non-zero within the intersection of the sources
      !
      izl = maxval(gr%active(1,3,(/src,dst/)))
      izu = minval(gr%active(2,3,(/src,dst/)))
      xl  = maxval(gr%active(1,1,(/src,dst/)))
      xu  = minval(gr%active(2,1,(/src,dst/)))
      yl  = maxval(gr%active(1,2,(/src,dst/)))
      yu  = minval(gr%active(2,2,(/src,dst/)))
      !
      call FieldAugmentBox(src,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(dst,igrid,xl,xu,yl,yu,izl,izu)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,dst)   &
                                      * gr%fields(xl:xu,yl:yu,iz,src)
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldMul')
    !
  end subroutine FieldMul
  !
  !  |dst> = |dst> + |src_a> * |src_b>
  !
  subroutine FieldMulAdd(src_a,src_b,dst,limits)
    integer(ik), intent(in)           :: src_a, src_b, dst
    integer(ik), intent(in), optional :: limits(:,:)        ! Grid summation limits for the outermost grid

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldMulAdd')
    if (verbose>=2) write (out,"(' = FieldMulAdd ',i3,' + ',i3,' * ',i3,' -> ',i3)") dst, src_a, src_b, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  The product field is non-zero within the intersection of the sources
      !
      izl = maxval(gr%active(1,3,(/src_a,src_b/)))
      izu = minval(gr%active(2,3,(/src_a,src_b/)))
      xl  = maxval(gr%active(1,1,(/src_a,src_b/)))
      xu  = minval(gr%active(2,1,(/src_a,src_b/)))
      yl  = maxval(gr%active(1,2,(/src_a,src_b/)))
      yu  = minval(gr%active(2,2,(/src_a,src_b/)))
      !
      !  The sum extent is the union of the extents
      !
      izl = min(gr%active(1,3,dst),izl)
      izu = max(gr%active(2,3,dst),izu)
      xl  = min(gr%active(1,1,dst),xl)
      xu  = max(gr%active(2,1,dst),xu)
      yl  = min(gr%active(1,2,dst),yl)
      yu  = max(gr%active(2,2,dst),yu)
      !
      if (present(limits).and.igrid==1) then
        xl = max( xl,limits(1,1))
        xu = min( xu,limits(2,1))
        yl = max( yl,limits(1,2))
        yu = min( yu,limits(2,2))
       izl = max(izl,limits(1,3))
       izu = min(izu,limits(2,3))
      end if
      !
      call FieldAugmentBox(src_a,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(src_b,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(dst,  igrid,xl,xu,yl,yu,izl,izu)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,dst)   &
                                      + gr%fields(xl:xu,yl:yu,iz,src_a) &
                                      * gr%fields(xl:xu,yl:yu,iz,src_b)
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src_a)
    call FieldShrinkBox(src_b)
    call FieldShrinkBox(dst)
    call TimerStop('FieldMulAdd')
    !
  end subroutine FieldMulAdd
  !
  !  |dst> = |dst> + alpha * |src>
  !
  subroutine FieldAXPY(alpha,src,dst)
    complex(rk), intent(in)    :: alpha
    integer(ik), intent(in)    :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldAXPY')
    if (verbose>=2) write (out,"(' = FieldAXPY ',i3,' + ',2g12.5,' * ',i3,' -> ',i3)") dst, alpha, src, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = minval(gr%active(1,3,(/src,dst/)))
      izu = maxval(gr%active(2,3,(/src,dst/)))
      xl  = minval(gr%active(1,1,(/src,dst/)))
      xu  = maxval(gr%active(2,1,(/src,dst/)))
      yl  = minval(gr%active(1,2,(/src,dst/)))
      yu  = maxval(gr%active(2,2,(/src,dst/)))
      call FieldAugmentBox(src,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(dst,igrid,xl,xu,yl,yu,izl,izu)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,dst) &
                              + alpha * gr%fields(xl:xu,yl:yu,iz,src)
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldAXPY')
    !
  end subroutine FieldAXPY
  !
  !  |dst> = |dst> + alpha * Abs(|src>)**2
  !
  subroutine FieldRhoAccumulate(alpha,src,dst)
    real(rk), intent(in)       :: alpha
    integer(ik), intent(in)    :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: iz, izl, izu
    integer(ik)                :: xl, xu, yl, yu

    call TimerStart('FieldRhoAccumulate')
    if (verbose>=2) write (out,"(' = FieldRhoAccumulate ',i3,' + ',2g12.5,' * ',i3,' -> ',i3)") dst, alpha, src, dst
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      izl = minval(gr%active(1,3,(/src,dst/)))
      izu = maxval(gr%active(2,3,(/src,dst/)))
      xl  = minval(gr%active(1,1,(/src,dst/)))
      xu  = maxval(gr%active(2,1,(/src,dst/)))
      yl  = minval(gr%active(1,2,(/src,dst/)))
      yu  = maxval(gr%active(2,2,(/src,dst/)))
      call FieldAugmentBox(src,igrid,xl,xu,yl,yu,izl,izu)
      call FieldAugmentBox(dst,igrid,xl,xu,yl,yu,izl,izu)
      !$omp parallel do private(iz)
      do iz=izl,izu
        gr%fields(xl:xu,yl:yu,iz,dst) = gr%fields(xl:xu,yl:yu,iz,dst) &
                              + alpha * abs(gr%fields(xl:xu,yl:yu,iz,src))**2
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldRhoAccumulate')
    !
  end subroutine FieldRhoAccumulate
  !
  !  Calculate approximate local inverse of the Hamiltonian:
  !
  !       1
  !   (- --- \Delta + v(r) - eps) \psi(r) = -R(r)
  !      2 m
  !
  !  Given R(r), we want to estimate \psi(r). To do this, we
  !  replace the laplacian with an error-diffusion multiplicative
  !  stencil C, making the problem invertible:
  !
  !                           R(r)
  !     \psi(r) = - --------------------------
  !                 - C(0,0)/(2m) - eps + v(r)
  !
  !  This function must be used together with FieldDiffuseLaplacian
  !  (which will appy the stencil to the local inverse) to get
  !  meaningful results. See preconditioning code in "liu.f90"
  !  for an example.
  !
  subroutine FieldInvertHamiltonian(src,dst,eps,mass)
    integer(ik), intent(in)    :: src   ! Field, containing multiplicative
                                        ! potential
    integer(ik), intent(in)    :: dst   ! Input: Vector to invert
                                        ! Output: Approximate local inverse
    real(rk), intent(in)       :: eps   ! Input: Approximate state energy
    real(rk), intent(in)       :: mass  ! Input: Particle mass
    !
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    real(rk)                   :: ps(-1:1,-1:1,-1:1)     ! Error diffusion stencil
    real(rk)                   :: unit_lap               ! Laplacian for normalized w.f.
    real(rk)                   :: base_shift             ! -C(0,0)/(2m) - eps
    real(rk)                   :: shift
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    integer(ik)                :: ix, iy, iz
    !
    call TimerStart('FieldInvertHamiltonian')
    if (verbose>=2) then
      write (out,"(' = FieldInvertHamiltonian ',i3,' -> ',i3,"// &
                 "' eps = ',e10.3,' mass = ',e10.3)") src, dst, eps, mass
    end if
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      xl  = gr%active(1,1,dst) ; xu  = gr%active(2,1,dst)
      yl  = gr%active(1,2,dst) ; yu  = gr%active(2,2,dst)
      zl  = gr%active(1,3,dst) ; zu  = gr%active(2,3,dst)
      call FieldAugmentBox(src,igrid,xl,xu,yl,yu,zl,zu)
      !
      !  Build the stencil for this grid level. Our stencil is
      !  normalized to give unit Laplacian. In computing the
      !  inverse, we actually need the value of Laplacian
      !  for unit wavefunction norm.
      !
      call kineticBuildStencil(gr%step,ps)
      unit_lap   = -1.0_rk/sqrt(gr%weight*sum(ps**2))
      base_shift = -unit_lap/(2*mass) - eps
      ! write (out,"('FieldInvertHamiltonian: base_shift = ',f15.8)") base_shift
      !
      !$omp parallel do private(ix,iy,iz,shift)
      z: do iz=zl,zu
        y: do iy=yl,yu
          x: do ix=xl,xu
            shift = base_shift + real(gr%fields(ix,iy,iz,src),kind=rk)
            if (abs(shift)<1e-2_rk) then
              shift = sign(1e-2_rk,shift)
            end if
            gr%fields(ix,iy,iz,dst) = - gr%fields(ix,iy,iz,dst) / shift
          end do x
        end do y
      end do z
      !$omp end parallel do
    end do
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldInvertHamiltonian')
  end subroutine FieldInvertHamiltonian
  !
  !  FieldDiffuseLaplacian applies laplacian error diffusion map
  !  to src. The only meaningful use is together with FieldInvertHamiltonian
  !
  subroutine FieldDiffuseLaplacian(src,dst)
    integer(ik), intent(in)    :: src   ! Input: Approximate local inverse from
                                        !        FieldInvertHamiltonian
    integer(ik), intent(in)    :: dst   ! Output: Error-diffused inverse, with
                                        !         laplacian stencil applied
    !
    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    real(rk)                   :: ps(-1:1,-1:1,-1:1) ! Error diffusion stencil
    integer(ik)                :: xl, xu, yl         ! Source dimensions
    integer(ik)                :: yu, zl, zu
    integer(ik)                :: xlt, xut, ylt      ! Target dimensions
    integer(ik)                :: yut, zlt, zut      ! Never learn French!
    integer(ik)                :: ix, iy, iz         ! Source indices
    integer(ik)                :: dx, dy, dz         ! Relative target indices
    integer(ik)                :: tx, ty, tz         ! Absolute target indices
    complex(rk)                :: c
    !
    call TimerStart('FieldDiffuseLaplacian')
    if (verbose>=2) then
      write (out,"(' = FieldDiffuseLaplacian ',i3,' -> ',i3)") src, dst
    end if
    !
    call FieldZero(dst)
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  Stencilled box will be one point wider than the local box.
      !  We'll need to make sure the dimenstions don't overflow.
      !
      xl  = gr%active(1,1,src) ; xu  = gr%active(2,1,src)
      yl  = gr%active(1,2,src) ; yu  = gr%active(2,2,src)
      zl  = gr%active(1,3,src) ; zu  = gr%active(2,3,src)
      !
      xlt = max(xl-1,1) ; xut = min(xu+1,gr%npoints(1))
      ylt = max(yl-1,1) ; yut = min(yu+1,gr%npoints(2))
      zlt = max(zl-1,1) ; zut = min(zu+1,gr%npoints(3))
      call FieldAugmentBox(dst,igrid,xlt,xut,ylt,yut,zlt,zut)
      !
      !  Build the stencil for this grid level. Because our
      !  stencils are normalized to unit Laplacian, we'll
      !  need to renormalize here.
      !
      call kineticBuildStencil(gr%step,ps)
      ps = ps / sqrt(gr%weight*sum(ps**2))
      !
      !$omp parallel do private(ix,iy,iz,c,dx,dy,dz,tx,ty,tz)
      src_z: do iz=zl,zu
        src_y: do iy=yl,yu
          src_x: do ix=xl,xu
            c = gr%fields(ix,iy,iz,src)
            !
            !  Applying the stencil. If we blow the grid boundary,
            !  we'll assign the density to the nearest point inside.
            !
            dst_z: do dz=-1,1
              tz = min(max(zlt,iz+dz),zut)
              dst_y: do dy=-1,1
                ty = min(max(ylt,iy+dy),yut)
                dst_x: do dx=-1,1
                  tx = min(max(xlt,ix+dx),xut)
                  !
                  gr%fields(tx,ty,tz,dst) = gr%fields(tx,ty,tz,dst) &
                                          + c * ps(dx,dy,dz)
                end do dst_x
              end do dst_y
            end do dst_z
            !
          end do src_x
        end do src_y
      end do src_z
      !$omp end parallel do
    end do
    call FieldShrinkBox(dst)
    call TimerStop('FieldDiffuseLaplacian')
  end subroutine FieldDiffuseLaplacian
  !
  !  Laplacian can also by synthesized out of a series of gradient
  !  calls. However, this is more efficient (one pass) and more
  !  accurate (it only has to reconcile grids once)
  !
  subroutine FieldLaplacian(src,dst)
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    integer(ik)                :: iz

    !
    call TimerStart('FieldLaplacian')
    if (verbose>=2) write (out,"(' = FieldLaplacian ',i3,' -> ',i3)") src,dst
    !
    call reconcileGrids(src)
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  In laplacian evaluation, the active part of the grid will expand
      !  by 1 point in each direction. This means we may have to augment
      !  the source grid by up to -two- points, to allow general-case
      !  derivative code to work. Unlike the target grid, the source grid
      !  is allowed to expand into the grid border.
      !
      xl  = max(gr%active(1,1,src)-1,1) ; xu  = min(gr%active(2,1,src)+1,gr%npoints(1))
      yl  = max(gr%active(1,2,src)-1,1) ; yu  = min(gr%active(2,2,src)+1,gr%npoints(2))
      zl  = max(gr%active(1,3,src)-1,1) ; zu  = min(gr%active(2,3,src)+1,gr%npoints(3))
      call FieldAugmentBox(dst,igrid,xl,  xu,  yl,  yu,  zl,  zu)
      call FieldAugmentBox(src,igrid,xl-1,xu+1,yl-1,yu+1,zl-1,zu+1)
      !
      !$omp parallel do private(iz)
      do iz=zl,zu
        gr%fields(xl:xu,yl:yu,iz,dst) =                                                     &
              (          gr%fields(xl+1:xu+1,  yl:yu,    iz,   src) +                       &
                         gr%fields(xl-1:xu-1,  yl:yu,    iz,   src) -                       &
                2.0_rk * gr%fields(  xl:xu,    yl:yu,    iz,   src) ) / (gr%step(1)**2) +   &
!
              (          gr%fields(  xl:xu,   yl+1:yu+1, iz,   src) +                       &
                         gr%fields(  xl:xu,   yl-1:yu-1, iz,   src) -                       &
                2.0_rk * gr%fields(  xl:xu,     yl:yu,   iz,   src) ) / (gr%step(2)**2) +   &
!
              (          gr%fields(  xl:xu,     yl:yu,   iz+1, src) +                       &
                         gr%fields(  xl:xu,     yl:yu,   iz-1, src) -                       &
                2.0_rk * gr%fields(  xl:xu,     yl:yu,   iz,   src) ) / (gr%step(3)**2)
      end do
      !$omp end parallel do
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldLaplacian')
    !
  end subroutine FieldLaplacian
  !
  !  An iteration of succesive over-relaxation Poisson solver.
  !  This code WILL NOT WORK on multi-grids. Sorry.
  !
  subroutine FieldIterationSOR(rho,pot,delta,omega,boundary)
    integer(ik), intent(in)           :: rho      ! Field containing the density
    integer(ik), intent(in)           :: pot      ! Input: guess for the potential
                                                  ! Output: (hopefully) improved potential
    real(rk), intent(out)             :: delta    ! Maximum change in the potential
    real(rk), intent(in), optional    :: omega    ! Over-relaxation parameter. If omega is
                                                  ! not given, we'll use Gauss-Seidel method.
    integer(ik), intent(in), optional :: boundary ! Number of pixels in pot containing
                                                  ! boundary conditions.
    !
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    integer(ik)                :: ix, iy, iz
    integer(ik)                :: edge
    real(rk)                   :: wx, wy, wz, wc, wr
!   real(rk)                   :: step
    complex(rk)                :: step
    !
    if (grid%ngrids/=1) then
      stop 'multigrid%FieldIterationSOR - multigrid case not implemented. Sorry.'
    end if
    !
    call TimerStart('FieldIterationSOR')
    if (verbose>=2) write (out,"(' = FieldIterationSOR ',i3,' -> ',i3)") rho, pot
    !
    gr => grid%grids(1)
    !
    !  The potential grid could (in principle) extend one point
    !  beyond the density grid OR the input potential grid, as
    !  long as we do not overwrite the boundary elements.
    !
    edge = 0
    if (present(boundary)) edge = boundary
    xl = max(min(gr%active(1,1,pot),gr%active(1,1,rho))-1,1+edge)
    yl = max(min(gr%active(1,2,pot),gr%active(1,2,rho))-1,1+edge)
    zl = max(min(gr%active(1,3,pot),gr%active(1,3,rho))-1,1+edge)
    !
    xu = min(max(gr%active(2,1,pot),gr%active(2,1,rho))+1,gr%npoints(1)-edge)
    yu = min(max(gr%active(2,2,pot),gr%active(2,2,rho))+1,gr%npoints(2)-edge)
    zu = min(max(gr%active(2,3,pot),gr%active(2,3,rho))+1,gr%npoints(3)-edge)
    !
    call FieldAugmentBox(rho,1_ik,xl,xu,yl,yu,zl,zu)
    call FieldAugmentBox(pot,1_ik,xl,xu,yl,yu,zl,zu)
    !
    !  SOR weights for each directional displacement, the centre, and the density
    !  We'll begin by constructing the Laplacian coefficients first.
    !
    wx = gr%step(1)**(-2) ; wy = gr%step(2)**(-2) ; wz = gr%step(3)**(-2)
    wc = 2._rk*(wx+wy+wz)
    !
    !  Scale the coefficients to form Gauss-Seidel weights
    !
    wx = wx / wc ; wy = wy / wc ; wz = wz / wc
    wr = fourpi / wc
    wc = -1._rk
    !
    !  Apply the over-relaxation scale factor, if present
    !
    if (present(omega)) then
      wx = omega * wx
      wy = omega * wy
      wz = omega * wz
      wc = omega * wc
      wr = omega * wr
    end if
    ! write (out,"(' SOR: ',5f20.12)") wx, wy, wz, wc, wr
    !
    !  The code below is not -quite- the SOR, as the order of iterations
    !  over Z slices is not deterministic. However, the convergence limit
    !  should be the same.
    !
    delta = 0
    !$omp parallel do private(ix,iy,iz,step) reduction(max:delta)
    slice_z: do iz=zl,zu
      slice_y: do iy=yl,yu
        slice_x: do ix=xl,xu
!         step = wx*real(gr%fields(ix+1,iy  ,iz  ,pot),kind=rk) + &
!                wx*real(gr%fields(ix-1,iy  ,iz  ,pot),kind=rk) + &
!                wy*real(gr%fields(ix  ,iy+1,iz  ,pot),kind=rk) + &
!                wy*real(gr%fields(ix  ,iy-1,iz  ,pot),kind=rk) + &
!                wz*real(gr%fields(ix  ,iy  ,iz+1,pot),kind=rk) + &
!                wz*real(gr%fields(ix  ,iy  ,iz-1,pot),kind=rk) + &
!                wc*real(gr%fields(ix  ,iy  ,iz  ,pot),kind=rk) + &
!                wr*real(gr%fields(ix  ,iy  ,iz  ,rho),kind=rk)
!         delta = max(delta,abs(step))
!         gr%fields(ix,iy,iz,pot) = real(gr%fields(ix,iy,iz,pot),kind=rk) + step
          step = wx*gr%fields(ix+1,iy  ,iz  ,pot) + &
                 wx*gr%fields(ix-1,iy  ,iz  ,pot) + &
                 wy*gr%fields(ix  ,iy+1,iz  ,pot) + &
                 wy*gr%fields(ix  ,iy-1,iz  ,pot) + &
                 wz*gr%fields(ix  ,iy  ,iz+1,pot) + &
                 wz*gr%fields(ix  ,iy  ,iz-1,pot) + &
                 wc*gr%fields(ix  ,iy  ,iz  ,pot) + &
                 wr*gr%fields(ix  ,iy  ,iz  ,rho)
          delta = max(delta,abs(step))
          gr%fields(ix,iy,iz,pot) = gr%fields(ix,iy,iz,pot) + step
        end do slice_x
      end do slice_y
    end do slice_z
    !$omp end parallel do
    !
    call FieldShrinkBox(rho)
    call FieldShrinkBox(pot)
    call TimerStop('FieldIterationSOR')
    !
  end subroutine FieldIterationSOR
  !
  !  Preparation for evaluating grid representation of an ECP: returns the extend
  !  and grid parameters, which should be used for calculating the projectors.
  !
  subroutine FieldECPGetExtent(prj,rmax,np)
    type(gam_structure), intent(in)  :: prj      ! GAMESS wavefunction representing the projectors
    real(rk), intent(in)             :: rmax     ! Max distance from the first atom of prj, where 
                                                 ! the projectors are still non-zero
    integer(ik),intent(out)          :: np(2,3)  ! Grid extent in each direction, in points
    !
    type(SimpleGridT), pointer  :: gr
    !
    if (grid%ngrids/=1) stop 'FieldECPGetExtent - multigrids not implemented'
    gr => grid%grids(1)
    call get_projector_grid_range(np,gr,prj,rmax)
  end subroutine FieldECPGetExtent
  !
  !  Calculates overlap integrals between a set of projectors (bras) and the
  !  wavefunction (the ket). This is 1/2 of the implementation of ECPs, with
  !  FieldECPApply being the other half. For simplicity, this function does
  !  not work on multi-grids. The reason for having a specialized function is
  !  that the projectors are known to be very short-range; using general routines
  !  will lead to enourmous computational cost.
  !
  subroutine FieldECPProject(ket,prj,cp,rmax,ovr)
    integer(ik), intent(in)         :: ket         ! "Right" wavefunction field
    type(gam_structure), intent(in) :: prj         ! GAMESS wavefunction representing the projectors
    real(rk), intent(in)            :: cp(:,:,:,:) ! Cartesian representation of the projectors
    real(rk), intent(in)            :: rmax        ! Max distance from the first atom of prj, where 
                                                   ! the projectors are still non-zero
    complex(rk), intent(out)        :: ovr(:)      ! Overlap integrals. The size of the buffer must match
                                                   ! number of the MOs in the projector.
    !
    type(SimpleGridT), pointer  :: gr
    integer(ik)                 :: nprj, ip
    integer(ik)                 :: np(2,3)
    !
    call TimerStart('FieldECPProject')
    if (verbose>=2) write (out,"(' FieldECPProject',3i4)") ket
    if (grid%ngrids/=1) stop 'multigrid%FieldECPProject - multigrids not implemented'
    !
    gr => grid%grids(1)
    !
    nprj = size(ovr)
    call get_projector_grid_range(np,gr,prj,rmax)
    !
    !  Sanity check: the extent and number of overlaps must match the size of the projector
    !
    if (any(shape(cp)/=(/np(2,:)-np(1,:)+1,nprj/))) then
      stop 'multigrid%FieldECPProject - inconsistent shapes'
    end if
    !
    call FieldAugmentBox(ket,1_ik,np(1,1),np(2,1),np(1,2),np(2,2),np(1,3),np(2,3))
    !
    !$omp parallel do default(none) private(ip) shared(ovr,gr,cp,nprj,np,ket)
    do_projectors: do ip=1,nprj
      ovr(ip) = sum(gr%coords(4,np(1,1):np(2,1),np(1,2):np(2,2),np(1,3):np(2,3))  &
                 *  gr%fields(  np(1,1):np(2,1),np(1,2):np(2,2),np(1,3):np(2,3),ket) &
                 *         cp(:,:,:,ip))
    end do do_projectors
    !$omp end parallel do
    call FieldShrinkBox(ket)
    call TimerStop('FieldECPProject')
  end subroutine FieldECPProject
  !
  !  |ket> = |ket> + Sum wgt(i) |prj(i)>
  !  This is a companion to FieldECPProject, primarily useful for implementing ECPs
  !
  subroutine FieldECPApply(ket,prj,cp,rmax,wgt)
    integer(ik), intent(in)         :: ket         ! Wavefunction field to be modified
    type(gam_structure), intent(in) :: prj         ! GAMESS wavefunction representing the projectors
    real(rk), intent(in)            :: cp(:,:,:,:) ! Cartesian representation of the projectors
    real(rk), intent(in)            :: rmax        ! Max distance from the first atom of prj, where 
                                                   ! the projectors are still non-zero
    complex(rk), intent(in)         :: wgt(:)      ! Projector weights. The size of the buffer must match
                                                   ! number of the MOs in the projector.
    !
    type(SimpleGridT), pointer  :: gr
    integer(ik)                 :: nprj, ix, iy, iz
    integer(ik)                 :: np(2,3)
    !
    call TimerStart('FieldECPApply')
    if (verbose>=2) write (out,"(' FieldECPApply',3i4)") ket
    if (grid%ngrids/=1) stop 'multigrid%FieldECPApply - multigrids not implemented'
    !
    gr => grid%grids(1)
    !
    nprj = size(wgt)
    call get_projector_grid_range(np,gr,prj,rmax)
    !
    !  Sanity check: the extent and number of overlaps must match the size of the projector
    !
    if (any(shape(cp)/=(/np(2,:)-np(1,:)+1,nprj/))) then
      stop 'multigrid%FieldECPApply - inconsistent shapes'
    end if
    !
    call FieldAugmentBox(ket,1_ik,np(1,1),np(2,1),np(1,2),np(2,2),np(1,3),np(2,3))
    !
    !$omp parallel do default(none) private(ix,iy,iz) shared(wgt,gr,cp,np,ket)
    add_z: do iz=np(1,3),np(2,3)
      add_y: do iy=np(1,2),np(2,2)
        add_x: do ix=np(1,1),np(2,1)
          gr%fields(ix,iy,iz,ket) = gr%fields(ix,iy,iz,ket) + sum(wgt*cp(ix-np(1,1)+1,iy-np(1,2)+1,iz-np(1,3)+1,:))
        end do add_x
      end do add_y
    end do add_z
    !$omp end parallel do
    !
    call FieldShrinkBox(ket)
    call TimerStop('FieldECPApply')
  end subroutine FieldECPApply
  ! Some shared code for FieldECPProject and FieldECPApply - figure out the integration
  ! extent in points. Very specialized!
  subroutine get_projector_grid_range(np,gr,prj,rmax)
    integer(ik), intent(out)        :: np(2,3)  ! Extent of the integration region in points
    type(SimpleGridT), pointer      :: gr       ! Current grid
    type(gam_structure), intent(in) :: prj      ! Data structure defining the projectors
    real(rk), intent(in)            :: rmax     ! Cut-off radius
    !
    integer(ik) :: ic, iat, ipt
    real(rk)    :: rpt
    real(rk)    :: rotatm(3)
    !
    check_directions: do ic=1,3
      !
      !  Establish the lower bound. If any of the atoms cause the point to enter
      !  the integration zone, then the boundary was at the previous point.
      !
      scan_up: do ipt=2,gr%npoints(ic)
        if (ic==1) rpt = gr%coords(1,ipt,1,1)
        if (ic==2) rpt = gr%coords(2,1,ipt,1)
        if (ic==3) rpt = gr%coords(3,1,1,ipt)
        scan_atoms_up: do iat=1,prj%natoms
          !
          !  We need to rotate the molecule into the desired orientation before looking at the
          !  grid limits. When we rotate the grid (see code involving rotmat in import_gamess),
          !  we use transpose(rotmat). To get the equivalent transformation for the molecule,
          !  use straight rotmat here.
          !
          rotatm = matmul(real(prj%rotmat,kind=kind(rotatm)),real(prj%atoms(iat)%xyz,kind=kind(rotatm)))
          if (rpt>rotatm(ic)/abohr-rmax) exit scan_up
        end do scan_atoms_up
      end do scan_up
      np(1,ic) = ipt-1
      !
      !  Now the upper bound
      !
      scan_down: do ipt=gr%npoints(ic)-1,1,-1
        if (ic==1) rpt = gr%coords(1,ipt,1,1)
        if (ic==2) rpt = gr%coords(2,1,ipt,1)
        if (ic==3) rpt = gr%coords(3,1,1,ipt)
        scan_atoms_down: do iat=1,prj%natoms
          rotatm = matmul(real(prj%rotmat,kind=kind(rotatm)),real(prj%atoms(iat)%xyz,kind=kind(rotatm)))
          if (rpt<rotatm(ic)/abohr+rmax) exit scan_down
        end do scan_atoms_down
      end do scan_down
      np(2,ic) = ipt+1
    end do check_directions
    !
    if (verbose>=2) then
      write (out,"('Projector point range: ',3(i5,'-',i5,'; '))") np
      write (out,"('Coordinates (min): ',3(1x,f16.6))") gr%coords(1:3,np(1,1),np(1,2),np(1,3))
      write (out,"('Coordinates (max): ',3(1x,f16.6))") gr%coords(1:3,np(2,1),np(2,2),np(2,3))
    end if
  end subroutine get_projector_grid_range
  !
  !  Tell caller how many components are available
  !
  function FieldComponentCount() result(c)
    integer(ik) :: c
    !
    c = 3
  end function FieldComponentCount
  !
  !  Tell caller the grid spacing at a given level (finest by default)
  !
  function FieldGridSpacing(level) result(d)
    integer(ik), intent(in), optional :: level
    real(rk)                          :: d(3)
    !
    if (present(level)) then
      d = grid%grids(level)%step
    else
      d = grid%grids(grid%ngrids)%step
    end if
  end function FieldGridSpacing
  !
  !  Return the number of points for a given grid level (finest by default)
  !
  function FieldGridNPoints(level) result(n)
    integer(ik), intent(in), optional :: level
    integer(ik)                       :: n(3)
    !
    if (present(level)) then
      n = grid%grids(level)%npoints
    else
      n = grid%grids(grid%ngrids)%npoints
    end if
  end function FieldGridNPoints
  !
  !  Return coordinates of grid positions for a given grid level (finest by default)
  !
  subroutine FieldGridCoordinates(component,coord,level)
    integer(ik)                       :: component !  1 = X,   2 = Y,   3 = Z
                                                   ! -1 = PX, -2 = PY, -3 = PZ
    integer(ik), intent(in), optional :: level
    real(rk)                          :: coord(:)
    integer(ik)                       :: lv
    !
    lv = grid%ngrids
    if (present(level)) lv = level
    !
    if (abs(component)<1 .or. abs(component)>3) then
      stop 'FieldGridCoordinates - bad component index'
    end if
    if (size(coord)/=grid%grids(lv)%npoints(abs(component))) then
      stop 'FieldGridCoordinates - bad coordinate array size'
    end if
    select case (component)
      case ( 1) ; coord = grid%grids(lv)%  coords(1,    1:grid%grids(lv)%npoints(1),1,1)
      case ( 2) ; coord = grid%grids(lv)%  coords(2,1,  1:grid%grids(lv)%npoints(2)  ,1)
      case ( 3) ; coord = grid%grids(lv)%  coords(3,1,1,1:grid%grids(lv)%npoints(3)    )
      case (-1) ; coord = grid%grids(lv)%p_coords(1,    1:grid%grids(lv)%npoints(1),1,1)
      case (-2) ; coord = grid%grids(lv)%p_coords(2,1,  1:grid%grids(lv)%npoints(2)  ,1)
      case (-3) ; coord = grid%grids(lv)%p_coords(3,1,1,1:grid%grids(lv)%npoints(3)    )
    end select
  end subroutine FieldGridCoordinates
  !
  !  Return grid extent for a given grid level (finest by default)
  !
  subroutine FieldGridExtent(component,coord,level)
    integer(ik)                       :: component ! 1 = X, 2 = Y, 3 = Z
    integer(ik), intent(in), optional :: level
    real(rk), intent(out)             :: coord(2)
    integer(ik)                       :: lv
    !
    lv = grid%ngrids
    if (present(level)) lv = level
    !
    if (component<1 .or. component>3) then
      stop 'FieldGridExtent - bad component index'
    end if
    select case (component)
      case (1) ; coord = grid%grids(lv)%coords(1,    (/1,grid%grids(lv)%npoints(1)/),1,1)
      case (2) ; coord = grid%grids(lv)%coords(2,1,  (/1,grid%grids(lv)%npoints(2)/)  ,1)
      case (3) ; coord = grid%grids(lv)%coords(3,1,1,(/1,grid%grids(lv)%npoints(3)/)    )
    end select
  end subroutine FieldGridExtent
  !
  !  Evaluate a spatial component of the Nabla (gradient) operator
  !
  subroutine FieldGradientComponent(dir,src,dst)
    integer(ik), intent(in) :: dir
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    integer(ik)                :: iz

    !
    call TimerStart('FieldGradientComponent')
    if (verbose>=2) write (out,"(' = FieldGradientComponent',i1,':',i3,' -> ',i3)") dir,src,dst
    !
    call reconcileGrids(src)
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  In gradient evaluation, the active part of the target grid will
      !  expand by 1 point in the direction of the operator component. To
      !  keep things simple, we'll expand the grid in all directions
      !  (just as we did for the Laplacian).
      !
      xl  = max(gr%active(1,1,src)-1,1) ; xu  = min(gr%active(2,1,src)+1,gr%npoints(1))
      yl  = max(gr%active(1,2,src)-1,1) ; yu  = min(gr%active(2,2,src)+1,gr%npoints(2))
      zl  = max(gr%active(1,3,src)-1,1) ; zu  = min(gr%active(2,3,src)+1,gr%npoints(3))
      call FieldAugmentBox(dst,igrid,xl,  xu,  yl,  yu,  zl,  zu)
      call FieldAugmentBox(src,igrid,xl-1,xu+1,yl-1,yu+1,zl-1,zu+1)
      !
      select case (dir)
        case default
          write (out,"(' FieldGradientComponent: bad component ',i8)") dir
          stop 'FieldGradientComponent - bad component'
        case (1)
          !$omp parallel do private(iz)
          grad_x: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl+1:xu+1,yl:yu,iz,src) &
                                            - gr%fields(xl-1:xu-1,yl:yu,iz,src) ) / (2*gr%step(1))
          end do grad_x
          !$omp end parallel do
        case (2)
          !$omp parallel do private(iz)
          grad_y: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl:xu,yl+1:yu+1,iz,src) &
                                            - gr%fields(xl:xu,yl-1:yu-1,iz,src) ) / (2*gr%step(2))
          end do grad_y
          !$omp end parallel do
        case (3)
          !$omp parallel do private(iz)
          grad_z: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl:xu,yl:yu,iz+1,src) &
                                            - gr%fields(xl:xu,yl:yu,iz-1,src) ) / (2*gr%step(3))
          end do grad_z
          !$omp end parallel do
      end select
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldGradientComponent')
    !
  end subroutine FieldGradientComponent
  !
  !  Evaluate a spatial component of the Nabla (gradient) operator,
  !  at right half-point positions. This operator is intended strictly
  !  for use within the Davidson diagonalizer, and may produce unexpected
  !  results otherwise. The main advantage over the central gradient
  !  generated by FieldGradientComponent is that convolution of two
  !  right gradients will give results consistent with taking expectation
  !  of the nabla operator.
  !
  subroutine FieldGradientComponentRight(dir,src,dst)
    integer(ik), intent(in) :: dir
    integer(ik), intent(in) :: src, dst

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    integer(ik)                :: iz
    real(rk)                   :: wgt
    !
    call TimerStart('FieldGradientComponentRight')
    if (verbose>=2) write (out,"(' = FieldGradientComponent',i1,':',i3,' -> ',i3)") dir,src,dst
    !
    call reconcileGrids(src)
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  In gradient evaluation, the active part of the target grid will
      !  expand by 1 point in the direction of the operator component. To
      !  keep things simple, we'll expand the grid in all directions
      !  (just as we did for the Laplacian).
      !
      xl  = max(gr%active(1,1,src)-1,1) ; xu  = min(gr%active(2,1,src)+1,gr%npoints(1))
      yl  = max(gr%active(1,2,src)-1,1) ; yu  = min(gr%active(2,2,src)+1,gr%npoints(2))
      zl  = max(gr%active(1,3,src)-1,1) ; zu  = min(gr%active(2,3,src)+1,gr%npoints(3))
      call FieldAugmentBox(dst,igrid,xl,  xu,  yl,  yu,  zl,  zu)
      call FieldAugmentBox(src,igrid,xl-1,xu+1,yl-1,yu+1,zl-1,zu+1)
      !
      wgt = 1._rk / gr%step(dir)
      !
      select case (dir)
        case default
          write (out,"(' FieldGradientComponent: bad component ',i8)") dir
          stop 'FieldGradientComponent - bad component'
        case (1)
          !$omp parallel do private(iz)
          grad_x: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl+1:xu+1,yl:yu,iz,src) &
                                            - gr%fields(xl  :xu  ,yl:yu,iz,src) ) * wgt
          end do grad_x
          !$omp end parallel do
        case (2)
          !$omp parallel do private(iz)
          grad_y: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl:xu,yl+1:yu+1,iz,src) &
                                            - gr%fields(xl:xu,yl  :yu  ,iz,src) ) * wgt
          end do grad_y
          !$omp end parallel do
        case (3)
          !$omp parallel do private(iz)
          grad_z: do iz=zl,zu
            gr%fields(xl:xu,yl:yu,iz,dst) = ( gr%fields(xl:xu,yl:yu,iz+1,src) &
                                            - gr%fields(xl:xu,yl:yu,iz  ,src) ) * wgt
          end do grad_z
          !$omp end parallel do
      end select
    end do
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    call TimerStop('FieldGradientComponentRight')
    !
  end subroutine FieldGradientComponentRight
  !
  !  The routine FieldFFT computes the FFTW of the field src and writes to the field dst
  !  The temporary array fftw_arr is needed for storing the input src and output dst fields
  !
  !  FieldFFT is a bit special, in the sense that it does -not- operate on the multigrid.
  !  The transform is done of the specified grid only. All other grids simply get zeros.
  !
  subroutine FieldFFT(src,dst,igrd_,inverse,normalize)
    integer(ik), intent(in)                :: src       ! Source field
    integer(ik), intent(in)                :: dst       ! Destination field
    integer(ik), intent(in), optional      :: igrd_     ! Grid level (1 by default)
    logical, intent(in), optional          :: inverse   ! Requests inverse transform
    character(len=*), intent(in), optional :: normalize ! Choice of normalization, one of:
                                                        ! 'WAVEFUNCTION' (default) or
                                                        ! 'POTENTIAL'
    integer(ik)                :: igrd
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: n1, n2, n3
    real(rk)                   :: norm
    logical                    :: invert
    character(len=20)          :: norm_choice
    integer(ik)                :: status
    complex(rk), allocatable   :: data(:,:,:) ! Unfortunately, FFTW -hates- array sections.
                                              ! To prevent creation of a stack temporary,
                                              ! we'll have to do this. Sod it.
    !
    call TimerStart('FieldFFT')
    if (verbose>=2) write (out,"(' = FieldFFT ',i3,' -> ',i3)") src, dst
    !
    igrd = 1
    if (present(igrd_)) igrd = igrd_
    !
    gr => grid%grids(igrd)
    !
    invert = .false.
    if (present(inverse)) invert = inverse
    !
    norm_choice = 'WAVEFUNCTION'
    if (present(normalize)) norm_choice = normalize
    !
    !  Because we can only look at a single grid level for FFT, the grids
    !  must be reconciled before we begin, unless we do an inverse transform
    !
    if (.not.invert) then
      call reconcileGrids(src)
    end if
    !
    !  Expand the source field, and zero the destination
    !
    n1 = gr%npoints(1) ; n2 = gr%npoints(2) ; n3 = gr%npoints(3)
    call FieldAugmentBox(src,igrd,1,n1,1,n2,1,n3)
    call FieldZero(dst)
    gr%fields(:,:,:,dst) = 0
    gr%active(:,:,dst)   = gr%active(:,:,0)
    !
    !  Copy data to a temporary buffer
    !
    allocate (data(n1,n2,n3),stat=status)
    if (status/=0) then
      write (out,"('FieldFFT: allocation of ',i6,' x ',i6,' x ',i6,' complex failed with code ',i9)") &
             n1, n2, n3, status
      stop 'multigrid%FieldFFT - out of memory'
    end if
    data = gr%fields(1:n1,1:n2,1:n3,src)
    !
    !  Calculate the transform
    !
    call fftw_3d(n1,n2,n3,data,invert)
    !
    !  Copy-out
    !
    gr%fields(1:n1,1:n2,1:n3,dst) = data
    deallocate (data)
    !
    !  The output from FFTW is not normalized. We need our transform
    !  normalized to number of particle in the simulation box.
    !
    select case (norm_choice)
      case default
        write (out,"('FieldFFT: normalization ',a,' is not known')") &
               trim(norm_choice)
      case ('WAVEFUNCTION')
        if (invert) then
          norm = gr%p_weight / sqrt(twopi**3)
        else
          norm = gr%weight / sqrt(twopi**3)
        end if
      case ('POTENTIAL')
        if (invert) then
          norm = 1.0_rk
        else
          norm = 1.0_rk / product(gr%npoints)
        end if
    end select
    !
    gr%fields(1:n1,1:n2,1:n3,dst) = norm * gr%fields(1:n1,1:n2,1:n3,dst)
    !
    !  Shrink back zero parts
    !
    call FieldShrinkBox(src)
    call FieldShrinkBox(dst)
    !
    !  Zero frequncies are found at element 1+(N-1)/2
    !  The momentum grid we set up earlier should give this
    !
    if (verbose>=1) then
      if (.not.invert) then
        write(out,"(' Stationary component: V = ',3f12.6,' Psi = ',2f12.6)") &
                  gr%p_coords(1:3,1+(n1-1)/2,1+(n2-1)/2,1+(n3-1)/2), &
                  gr%fields(      1+(n1-1)/2,1+(n2-1)/2,1+(n3-1)/2,dst)
      end if
    !
    !  Check normalization
    !
      if (norm_choice=='WAVEFUNCTION') then
        if (invert) then
          write(out,"(' Total probability in the velocity field   = ',f12.7)") &
                gr%p_weight*sum(abs(gr%fields(1:n1,1:n2,1:n3,src))**2)
          write(out,"(' Total probability in the coordinate field = ',f12.7)") &
                gr%weight  *sum(abs(gr%fields(1:n1,1:n2,1:n3,dst))**2)
        else
          write(out,"(' Total probability in the coordinate field = ',f12.7)") &
                gr%weight  *sum(abs(gr%fields(1:n1,1:n2,1:n3,src))**2)
          write(out,"(' Total probability in the velocity field   = ',f12.7)") &
                gr%p_weight*sum(abs(gr%fields(1:n1,1:n2,1:n3,dst))**2)
        end if
      end if
    end if
    call TimerStop('FieldFFT')
  end subroutine FieldFFT

  !
  !  Returns a pointer to the actual field data, at the requested level of
  !  grid. This is useful for specialized analysis routines - e.g. partial
  !  integrations or convolutions.
  !
  subroutine FieldCheckOut(typ,src,grd,npts,step,coord,field,expand)
    character(len=*), intent(in)  :: typ            ! Field type, either 'COORDINATE', or 'MOMENTUM'
    integer(ik), intent(in)       :: src            ! Field to return
    integer(ik), intent(in)       :: grd            ! Grid level
    integer(ik), intent(out)      :: npts(3)        ! Number of points at this grid level
    real(rk), intent(out)         :: step(3)        ! Grid spacing
    real(rk), pointer             :: coord(:,:,:,:) ! Coordinates and integration weights
    complex(rk), pointer          :: field(:,:,:)   ! Field values
    logical, intent(in), optional :: expand

    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu

    !
    call TimerStart('FieldCheckOut')
    if (verbose>=2) write (out,"(' = FieldCheckOut ',i3)") src
    !
    gr => grid%grids(grd)
    !
    if (present(expand)) then
      call FieldAugmentBox(src,grd,1,gr%npoints(1),1,gr%npoints(2),1,gr%npoints(3))
    end if
    !
    xl = gr%active(1,1,src) ; xu = gr%active(2,1,src)
    yl = gr%active(1,2,src) ; yu = gr%active(2,2,src)
    zl = gr%active(1,3,src) ; zu = gr%active(2,3,src)
    npts = (/ xu - xl + 1, yu - yl + 1, zu - zl + 1 /)
    !
    field => gr%fields(xl:xu,yl:yu,zl:zu,src)
    !
    select case (typ)
      case default
        write (out,"(' FieldCheckOut: grid type ',a,' is unknown')") typ
        stop 'FieldCheckOut - bad grid type'
      case ('COORDINATE')
        step = gr%step
        coord => gr%coords(:,xl:xu,yl:yu,zl:zu)
      case ('MOMENTUM')
        step = gr%p_step
        coord => gr%p_coords(:,xl:xu,yl:yu,zl:zu)
    end select
    call TimerStop('FieldCheckOut')
  end subroutine FieldCheckOut
  !
  !  FieldShrinkBox is responsible for detecting zero parts of the wave
  !  function at the box boundary, and removing them from the active
  !  part of the box, updating the corresponding active lists.
  !
  subroutine FieldShrinkBox(ind)
    integer(ik), intent(in) :: ind ! Field to shrink boxes for

    integer(ik)                :: igrid
    type(SimpleGridT), pointer :: gr
    integer(ik)                :: xl, xu, yl, yu, zl, zu

    if (zeroField<=0._rk) return ! This is a no-op
    !
    call TimerStart('FieldShrinkBox')
    if (verbose>=2) write (out,"(' = FieldShrinkBox ',i8)") ind
    !
    do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      xl = gr%active(1,1,ind) ; xu = gr%active(2,1,ind)
      yl = gr%active(1,2,ind) ; yu = gr%active(2,2,ind)
      zl = gr%active(1,3,ind) ; zu = gr%active(2,3,ind)
      !
      !$omp parallel
      !
      !  Reduce the box along Z direction, coming from both sides
      !
      !$omp sections
      !$omp section
      !  Going up along Z
        up_z: do while( zl <= zu )
          if (any( abs(gr%fields(xl:xu,yl:yu,zl,ind)) >= zeroField )) exit up_z
          !$omp atomic
          zl = zl + 1
        end do up_z
      !$omp section
      !  Going down along Z
        down_z: do while( zl <= zu )
          if (any( abs(gr%fields(xl:xu,yl:yu,zu,ind)) >= zeroField )) exit down_z
          !$omp atomic
          zu = zu - 1
        end do down_z
      !$omp end sections
      !
      !  Reduce the box along Y direction. We can already use the shrunk
      !  extent along Z, which cuts down on the cost
      !
      !$omp sections
      !$omp section
      !  Going up along Y
        up_y: do while( yl <= yu )
          if (any( abs(gr%fields(xl:xu,yl,zl:zu,ind)) >= zeroField )) exit up_y
          !$omp atomic
          yl = yl + 1
        end do up_y
      !$omp section
      !  Going down along Y
        down_y: do while( yl <= yu )
          if (any( abs(gr%fields(xl:xu,yu,zl:zu,ind)) >= zeroField )) exit down_y
          !$omp atomic
          yu = yu - 1
        end do down_y
      !$omp end sections
      !
      !  Finally, reduce the box along X
      !
      !$omp sections
      !$omp section
      !  Going up along X
        up_x: do while( xl <= xu )
          if (any( abs(gr%fields(xl,yl:yu,zl:zu,ind)) >= zeroField )) exit up_x
          !$omp atomic
          xl = xl + 1
        end do up_x
      !$omp section
      !  Going down along X
        down_x: do while( xl <= xu )
          if (any( abs(gr%fields(xu,yl:yu,zl:zu,ind)) >= zeroField )) exit down_x
          !$omp atomic
          xu = xu - 1
        end do down_x
      !$omp end sections
      !$omp end parallel
      !
      if (verbose>=2) then
        write (out,"(' Shrunk grid ',i2,' of field ',i2)") igrid, ind
        write (out,"(' Original extent : ',3(4x,i4,1x,i4))") gr%active(:,:,ind)
        write (out,"(' Reduced extent  : ',3(4x,i4,1x,i4))") xl, xu, yl, yu, zl, zu
      end if
      !
      !  Store the new limits for this box
      !
      gr%active(1,1,ind) = xl ; gr%active(2,1,ind) = xu
      gr%active(1,2,ind) = yl ; gr%active(2,2,ind) = yu
      gr%active(1,3,ind) = zl ; gr%active(2,3,ind) = zu
    end do
    call TimerStop('FieldShrinkBox')
    !
  end subroutine FieldShrinkBox
  !
  !  FieldAugmentBox guarantees that the active part of the box is
  !  at least as large as the specified dimensions. If necessary, it
  !  will fill the extra section(s) of the box with zeros.
  !
  subroutine FieldAugmentBox(ind,igrid,dxl,dxu,dyl,dyu,dzl,dzu)
    integer(ik), intent(in) :: ind        ! Field to be augmented
    integer(ik), intent(in) :: igrid      ! Grid level to work on
    integer(ik), intent(in) :: dxl, dxu   ! Desired box size
    integer(ik), intent(in) :: dyl, dyu
    integer(ik), intent(in) :: dzl, dzu

    type(SimpleGridT), pointer :: gr
    integer(ik)                :: ix, iy, iz
    integer(ik)                ::  xl,  xu,  yl,  yu,  zl,  zu  ! "New" box extent values
    integer(ik)                :: oxl, oxu, oyl, oyu, ozl, ozu  ! "Old" box extent values
    !
    if (verbose>=2) write (out,"(' = FieldAugmentBox ',2i8)") ind, igrid
    !
    xl = dxl ; yl = dyl ; zl = dzl
    xu = dxu ; yu = dyu ; zu = dzu
    !
    gr => grid%grids(igrid)
    !
    !  If active box already encompasses the requested box, we don't
    !  have to do anything.
    !
    if ( all( gr%active(1,:,ind) <= (/ xl, yl, zl /) ) .and. &
         all( gr%active(2,:,ind) >= (/ xu, yu, zu /) ) ) then
      return
    end if
    !
    call TimerStart('FieldAugmentBox')
    !
    !  Make sure requested boundaries are cosher
    !
    if ( any( gr%active(1,:,0) > (/ xl, yl, zl /) ) .or. &
         any( gr%active(2,:,0) < (/ xu, yu, zu /) ) ) then
      write (out,"(' FieldAugmentBox - requested to grow to ',3(2x,i4,1x))") &
             xl, xu, yl, yu, zl, zu
      write (out,"(' The limit for grid ',i2,' is ',3(2x,i4,1x))") &
             igrid, gr%active(:,:,0)
      stop 'FieldAugmentBox - can''t grow active region outside the box!'
    end if
    !
    !  Make sure the new active region is at least not smaller than the old
    !  one along all directions.
    !
    oxl = gr%active(1,1,ind) ; oxu = gr%active(2,1,ind)
    oyl = gr%active(1,2,ind) ; oyu = gr%active(2,2,ind)
    ozl = gr%active(1,3,ind) ; ozu = gr%active(2,3,ind)
    xl = min(xl,oxl) ; xu = max(xu,oxu)
    yl = min(yl,oyl) ; yu = max(yu,oyu)
    zl = min(zl,ozl) ; zu = max(zu,ozu)
    !
    if (verbose>=2) then
      write (out,"(' Growing grid ',i2,' of field ',i2)") igrid, ind
      write (out,"(' Old size : ',3(2x,i4,1x,i4))") oxl, oxu, oyl, oyu, ozl, ozu
      write (out,"(' New size : ',3(2x,i4,1x,i4))") xl,  xu,  yl,  yu,  zl,  zu
    end if
    !
    !  Fill the additional region with zeros. Filling large contiguos
    !  sections is more efficient, so that we'll take the largest
    !  possible slab in the XY plane, then fill the remaining chunks.
    !
    !  Most of the time, we'll be filling rather small chunks along
    !  each direction, so it's not worth it parallelizing the ix/iy/iz
    !  loops.
    !
    !$omp parallel private(ix,iy,iz)
    !$omp sections
    !$omp section
    !  Higher-Z region
      do iz=ozu+1,zu
        gr%fields(xl:xu,yl:yu,iz,ind) = 0
      end do
    !$omp section
    !  Lower-Z region
      do iz=zl,ozl-1
        gr%fields(xl:xu,yl:yu,iz,ind) = 0
      end do
    !$omp section
    !  Higher-Y region
      do iy=oyu+1,yu
        gr%fields(xl:xu,iy,ozl:ozu,ind) = 0
      end do
    !$omp section
    !  Lower-Y region
      do iy=yl,oyl-1
        gr%fields(xl:xu,iy,ozl:ozu,ind) = 0
      end do
    !$omp section
    !  Higher-X region
      do ix=oxu+1,xu
        gr%fields(ix,oyl:oyu,ozl:ozu,ind) = 0
      end do
    !$omp section
    !  Lower-X region
      do ix=xl,oxl-1
        gr%fields(ix,oyl:oyu,ozl:ozu,ind) = 0
      end do
    !$omp end sections
    !$omp end parallel
    !
    !  Store the new active extents, and we are done
    !
    gr%active(1,1,ind) = xl ; gr%active(2,1,ind) = xu
    gr%active(1,2,ind) = yl ; gr%active(2,2,ind) = yu
    gr%active(1,3,ind) = zl ; gr%active(2,3,ind) = zu
    !
    call TimerStop('FieldAugmentBox')
  end subroutine FieldAugmentBox
  !
  !  Checkpoint routines - dump or recover the whole state of the multigrip
  !
  subroutine FieldCheckpoint(action,name)
    character(len=*), intent(in) :: action ! 'SAVE' or 'RESTORE'
    character(len=*), intent(in) :: name   ! File name

    integer(ik)                :: igrid, ifield
    type(SimpleGridT), pointer :: gr

    call TimerStart('FieldCheckpoint')
    select case (action)
      case default
        write (out,"(' FieldCheckpoint - action ',a,' is not valid')") trim(action)
        stop 'FieldCheckpoint - bogus command'
      case ('SAVE')
        call checkpointSave
      case ('RESTORE')
        call checkpointRestore
    end select
    call TimerStop('FieldCheckpoint')

    contains

      subroutine checkpointSave
        integer(ik) :: xl, xu, yl, yu, zl, zu
        !
        open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=name)
          write(chkptIO) 'Start Multigrid'
          write(chkptIO) grid%ngrids_max, grid%nfields_max, grid%nborder, grid%ngrids, grid%nfields
          write(chkptIO) grid%field_names
          write(chkptIO) grid%grid_names
          do igrid=1,grid%ngrids
            gr => grid%grids(igrid)
            write(chkptIO) gr%npoints, gr%range, gr%step, gr%weight, gr%up_grid, gr%up_index
            write(chkptIO) gr%p_range, gr%p_step, gr%p_weight
            if (igrid/=1) &
            write(chkptIO) gr%down_index
            write(chkptIO) gr%down_grid
            write(chkptIO) gr%coords
            write(chkptIO) gr%p_coords
            write(chkptIO) gr%active
!
!          New-style checkpoint
!
            do ifield=1,grid%nfields
              if (grid%scratch(ifield)) cycle
              xl = gr%active(1,1,ifield) ; xu = gr%active(2,1,ifield) ;
              yl = gr%active(1,2,ifield) ; yu = gr%active(2,2,ifield) ;
              zl = gr%active(1,3,ifield) ; zu = gr%active(2,3,ifield) ;
              write(chkptIO) gr%fields(xl:xu,yl:yu,zl:zu,ifield)
            end do
          end do
          write(chkptIO) 'End Multigrid'
        close(chkptIO,status='keep')
        !
      end subroutine checkpointSave

      subroutine checkpointRestore
        character(len=15) :: buf
        integer(ik)       :: t_ngrids_max, t_nfields_max, t_nborder, t_ngrids, t_nfields
        integer(ik)       :: t_npoints(3)
        integer(ik)       :: xl, xu, yl, yu, zl, zu
        !
        open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=name)
          read(chkptIO) buf
          if (buf/='Start Multigrid') then
            write (out,"(' Checkpoint file ',a,' has bogus header: ',a)") name, buf
            stop 'FieldCheckpoint - bogus file format (1)'
          end if
          read(chkptIO) t_ngrids_max, t_nfields_max, t_nborder, t_ngrids, t_nfields
          if (t_ngrids_max/=grid%ngrids_max .or. t_nfields_max/=grid%nfields_max .or. &
              t_nborder   /=grid%nborder    .or. t_ngrids     /=grid%ngrids      .or. &
              t_nfields   /=grid%nfields) then
            write(out,"(' FieldCheckpoint parameter mismatch: ')")
            write(out,"(' ngrids_max = ',i8,'/',i8,' nfields_max = ',i8,'/',i8)") &
                  t_ngrids_max, grid%ngrids_max, t_nfields_max, grid%nfields_max
            write(out,"(' nborder    = ',i8,'/',i8,' ngrids      = ',i8,'/',i8)") &
                  t_nborder,    grid%nborder,    t_ngrids,      grid%ngrids
            write(out,"(' nfields    = ',i8,'/',i8)") t_nfields, grid%nfields
            stop 'FieldCheckpoint - parameter mismatch'
          end if
          read(chkptIO) grid%field_names
          read(chkptIO) grid%grid_names
          do igrid=1,grid%ngrids
            gr => grid%grids(igrid)
            t_npoints = gr%npoints
            read(chkptIO) gr%npoints, gr%range, gr%step, gr%weight, gr%up_grid, gr%up_index
            if (any(t_npoints/=gr%npoints)) then
              write(out,"(' FieldCheckpoint - grid ',i8,' in memory is ',3i8,'; on file ',3i8)") &
                    igrid, t_npoints, gr%npoints
              stop 'FieldCheckpoint - grid size mismatch'
            end if
            read(chkptIO) gr%p_range, gr%p_step, gr%p_weight
            if (igrid/=1) &
            read(chkptIO) gr%down_index
            read(chkptIO) gr%down_grid
            read(chkptIO) gr%coords
            read(chkptIO) gr%p_coords
            read(chkptIO) gr%active
!
!          New-style checkpoint
!
            do ifield=1,grid%nfields
              if (grid%scratch(ifield)) cycle
              xl = gr%active(1,1,ifield) ; xu = gr%active(2,1,ifield) ;
              yl = gr%active(1,2,ifield) ; yu = gr%active(2,2,ifield) ;
              zl = gr%active(1,3,ifield) ; zu = gr%active(2,3,ifield) ;
              read(chkptIO) gr%fields(xl:xu,yl:yu,zl:zu,ifield)
            end do
          end do
          buf = ' '
          read(chkptIO) buf(1:13)
          if (buf/='End Multigrid') then
            write (out,"(' Checkpoint file ',a,' has bogus final header: ',a)") name, buf
            stop 'FieldCheckpoint - bogus file format (2)'
          end if
        close(chkptIO,status='keep')
        !
      end subroutine checkpointRestore

  end subroutine FieldCheckpoint
  !
  !  External I/O
  !
  subroutine FieldIO(action,field,slot)
    character(len=*), intent(in)      :: action ! 'OPEN', 'CLOSE', 'READ', or 'WRITE'
    integer(ik), intent(in), optional :: field  ! Field to read/write
    integer(ik), intent(in), optional :: slot   ! Slot in the direct access file

    integer(ik) :: record                       ! Current direct access record number
    integer(ik) :: rec_len, max_rec_len
    integer(ik) :: igrid

    call TimerStart('FieldIO')
    if (verbose>=2) write (out,"(' FieldIO ',a)") action
    select case (action)
      case default
        write (out,"(' FieldIO: unknown action ',a)") action
        stop 'FieldIO - bad command'
      case ('OPEN')
        max_rec_len = 0
        do igrid=1,grid%ngrids
          inquire(iolength=rec_len) grid%grids(igrid)%active(:,:,1), &
                                    grid%grids(igrid)%fields(:,:,:,1)
          max_rec_len = max(max_rec_len,rec_len)
        end do
        open (GridIO,access='DIRECT',status='SCRATCH',recl=max_rec_len)
      case ('CLOSE')
        close (GridIO)
      case ('READ')
        call parseArgs
        do igrid=1,grid%ngrids
          read (GridIO,rec=record) grid%grids(igrid)%active(:,:,field), &
                                   grid%grids(igrid)%fields(:,:,:,field)
          record = record + 1
        end do
      case ('WRITE')
        call parseArgs
        do igrid=1,grid%ngrids
          write (GridIO,rec=record) grid%grids(igrid)%active(:,:,field), &
                                    grid%grids(igrid)%fields(:,:,:,field)
          record = record + 1
        end do
    end select
    call TimerStop('FieldIO')

    contains
      subroutine parseArgs
        if (.not.present(field) .or. .not.present(slot)) then
          write (out,"(' FieldIO action ',a,' requires field and slot arguments.')") action
          stop 'FieldIO - missing argument'
        end if
        record = 1 + grid%ngrids*(slot-1)
      end subroutine parseArgs
      !
  end subroutine FieldIO
  !
  !  Simple export routines. The format of the 'binary' dump is a little
  !  involved: it is designed to allow importing data from a smaller grid,
  !  as long matching coordinates are present in the bigger grid.
  !
  subroutine FieldExport(type,src,file)
    character(len=*), intent(in) :: type  ! Type of the data file to export
    integer(ik), intent(in)      :: src   ! Multigrid field index to export
    character(len=*), intent(in) :: file  ! Name of the data file to export
    !
    integer(ik)                :: igrid, alloc
    integer(ik)                :: xl, xu, yl, yu, zl, zu
    type(SimpleGridT), pointer :: gr
    complex(rk), allocatable   :: tmp(:,:,:)
    !
    call TimerStart('FieldExport')
    if (verbose>=2) then
      write (out,"(' FieldExport ',a,' src= ',i4,' to ',a)") trim(type), src, trim(file)
    end if
    !
    !  Sanity checking
    !
    if ( src<1 .or. src>grid%nfields_max ) then
      stop 'multigrid%FieldExport - src indices are out of range'
    end if
    if (type/='binary') stop 'multigrid%FieldExport - type is not binary'
    !
    open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=file)
    write (chkptIO) grid%ngrids
    grid_loop: do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      xl = gr%active(1,1,src) ; xu = gr%active(2,1,src)
      yl = gr%active(1,2,src) ; yu = gr%active(2,2,src)
      zl = gr%active(1,3,src) ; zu = gr%active(2,3,src)
      !
      write (chkptIO) xu-xl+1, yu-yl+1, zu-zl+1
      write (chkptIO) gr%step
      write (chkptIO) gr%coords(1:3,xl,yl,zl)
      allocate (tmp(xl:xu,yl:yu,zl:zu),stat=alloc)
      if (alloc/=0) then
        write (out,"('Allocating array: ',3(1x,i5,':',i5),' - code ',i9)") xl,xu,yl,yu,zl,zu,alloc
        stop 'multigrid%FieldExport - allocation failed'
      end if
      tmp = gr%fields(xl:xu,yl:yu,zl:zu,src)
      write (chkptIO) tmp
      deallocate (tmp)
    end do grid_loop
    !
    close(chkptIO)
    !
    call TimerStop('FieldExport')
  end subroutine FieldExport
  !
  !  Import routines
  !
  subroutine FieldImport(type,file,dst,indices,rot)
    character(len=*), intent(in)    :: type        ! Type of the data file to import
    character(len=*), intent(in)    :: file        ! Name of the data file to import
    integer(ik), intent(in)         :: dst(:)      ! Multigrid field indices to import
    integer(ik), intent(in)         :: indices(:)  ! Data object indices to import, one per dst() field
    real(ark), intent(in), optional :: rot(:,:)    ! Rotation matrix; only relevant for type=GAMESS
    !
    integer(ik)                :: idst, igrid, alloc
    integer(ik)                :: lg(3), ug(3)
    integer(ik)                :: ckpt_ngrid, ckpt_n(3)
    real(rk)                   :: ckpt_dx(3), ckpt_r0(3)
    type(SimpleGridT), pointer :: gr
    complex(rk), allocatable   :: tmp(:,:,:)
    !
    call TimerStart('FieldImport')
    if (verbose>=2) then
      write (out,"(' FieldImport ',a,' from ',a)") trim(type), trim(file)
      write (out,"('           dst = ',30i3)") dst
      write (out,"('       indices = ',30i3)") indices
    end if
    !
    !  Sanity checking
    !
    if (size(dst)/=size(indices)) then
      stop 'multigrid%FieldImport - dst/indices size mismatch'
    end if
    if ( any(dst<1) .or. any(dst>grid%nfields_max) ) then
      stop 'multigrid%FieldImport - dst indices are out of range'
    end if
    !
    if (type=='binary') then
      open (chkptIO,form='unformatted',action='read',position='rewind',status='old',file=file)
      read (chkptIO) ckpt_ngrid
      if (ckpt_ngrid/=grid%ngrids) then
        write (out,"(' FieldImport(',a,'): checkpoint has ',i3,' grids; we use ',i3)") &
               trim(file), ckpt_ngrid, grid%ngrids
        stop 'multigrid%FieldImport - incompatible grid layout'
      end if
    end if
    !
    !  Prepare the grid boundary
    !
    grid_loop: do igrid=1,grid%ngrids
      gr => grid%grids(igrid)
      !
      !  Choose the appropriate import subroutine. We'll fill only
      !  the fines-level grid. Because import routines only need to
      !  to "know" about single-level grids, we may call them more
      !  than once.
      !
      select case (type)
        case default
          write (out,"('multigrid%FieldImport: data type ',a,' not recognized')") trim(type)
          stop 'multigrid%FieldImport - bad data type'
        case ('binary')
          read (chkptIO) ckpt_n
          read (chkptIO) ckpt_dx
          read (chkptIO) ckpt_r0
          !
          !  Check grids for conformity, and determine the range of data points
          !
          if ( any(abs(ckpt_dx-gr%step)>10*spacing(gr%step)) ) then
            write (out,"(' FieldImport(',a,'): grid spacing is ',3g16.9,'; we use ',3g16.9)") &
                   trim(file), ckpt_dx, gr%step
            stop 'multigrid%FieldImport - grid steps not compatible'
          end if
          lg = 1 + nint((ckpt_r0 - gr%coords(1:3,1,1,1))/gr%step(:))
          if ( any(lg<1-grid%nborder) .or. any(lg>gr%npoints+grid%nborder) ) then
            write (out,"(' FieldImport(',a,'): matching lower bound is ',3i9)") trim(file), lg
            stop 'multigrid%FieldImport - grid is not large enough (1)'
          end if
          if ( any(abs(ckpt_r0-gr%coords(1:3,lg(1),lg(2),lg(3)))>10*spacing(ckpt_r0)) ) then
            write (out,"(' FieldImport(',a,'): lower corner at ',3g16.9,'; we expect ',3g16.9)") &
                   trim(file), ckpt_r0, gr%coords(1:3,lg(1),lg(2),lg(3))
            stop 'multigrid%FieldImport - grid origins not compatible'
          end if
          ug = lg + ckpt_n(:) - 1
          if ( any(ug<1-grid%nborder) .or. any(ug>gr%npoints+grid%nborder) ) then
            write (out,"(' FieldImport(',a,'): matching upper bound is ',3i9)") trim(file), ug
            stop 'multigrid%FieldImport - grid is not large enough (2)'
          end if
          gr%active(1,:,dst(1)) = lg
          gr%active(2,:,dst(1)) = ug
          !
          write (out,"(' FieldImport(',a,'): loading to range ',3i9,' through ',3i9)") &
                 trim(file), lg, ug
          !
          allocate (tmp(lg(1):ug(1),lg(2):ug(2),lg(3):ug(3)),stat=alloc)
          if (alloc/=0) then
            write (out,"('Allocating array: ',3(1x,i5,':',i5),' - code ',i9)") &
                   lg(1),ug(1),lg(2),ug(2),lg(3),ug(3),alloc
            stop 'multigrid%FieldImport - allocation failed'
          end if
          read (chkptIO) tmp
          gr%fields(lg(1):ug(1),lg(2):ug(2),lg(3):ug(3),dst(1)) = tmp
          deallocate (tmp)
        case ('GAMESS')
          lg = spread(1,1,3)
          ug = gr%npoints
          extent_loop: do idst=1,size(dst)
            gr%active(1,:,dst(idst)) = lg
            gr%active(2,:,dst(idst)) = ug
          end do extent_loop
          if (present(rot)) then
            call gamess_load_orbitals(file,indices,dst, &
                   gr%coords(1:3,lg(1):ug(1),lg(2):ug(2),lg(3):ug(3)), &
                   gr%fields    (lg(1):ug(1),lg(2):ug(2),lg(3):ug(3),:), &
                   rot=rot,dx=gr%step)
          else
            call gamess_load_orbitals(file,indices,dst, &
                   gr%coords(1:3,lg(1):ug(1),lg(2):ug(2),lg(3):ug(3)), &
                   gr%fields    (lg(1):ug(1),lg(2):ug(2),lg(3):ug(3),:), &
                   dx=gr%step)
          end if
      end select
    end do grid_loop
    !
    if (type=='binary') then
      close(chkptIO)
    end if
    !
    !  Reconcile coarse grids. This will have a side effect of
    !  shrinking bounding boxes to their proper values.
    !
    reconcile_loop: do idst=1,size(dst)
      call reconcileGrids(dst(idst))
    end do reconcile_loop
    call TimerStop('FieldImport')
  end subroutine FieldImport
  !
  !  Useful auxiliary routines
  !

  !
  !  disjoint() determines whether two 3D boxes are disjoint
  !
  logical function disjoint(r1,r2)
    real(rk), intent(in)  :: r1(2,3)  ! 3D spatial extents to be compared
    real(rk), intent(in)  :: r2(2,3)

    integer(ik) :: i

    disjoint = .true.
    do i=1,3
      if (r1(2,i) < r2(1,i)) return
      if (r1(1,i) > r2(2,i)) return
    end do
    disjoint = .false.

  end function disjoint

  !
  !  contained() determines whether one of two 3D boxes is inside the other one
  !
  logical function contained(inner,outer)
    real(rk), intent(in)  :: inner(2,3)  ! 3D spatial extents to be compared
    real(rk), intent(in)  :: outer(2,3)

    integer(ik) :: i

    contained = .false.
    do i=1,3
      if (inner(1,i) <= outer(1,i)) return
      if (inner(2,i) >= outer(2,i)) return
    end do
    contained = .true.

  end function contained

  !
  !  reconcileGrids() makes sure different grid levels agree for the specified
  !                   field. If wallType is set to 'ADSORBING', this will also
  !                   calculate the amount of charge, which leaked through the
  !                   wall.
  !
  subroutine reconcileGrids(index)
    integer(ik), intent(in) :: index ! Field to reconcile between the grids

    call TimerStart('reconcileGrids')
    if (verbose>=2) write (out,"(' ReconcileGrids ',i3)") index
    chargeFlux = 0
    call reconcileGridsUpwards  (index)
    call reconcileGridsDownwards(index)
    if (wallType=='ADSORBING') then
      call reconcileOuterWall(index)
    end if
    call calculateChargeFlux(index)
    if (verbose>=2) write (out,"(' Charge flux = ',f12.8)") chargeFlux
    if (wallType/='ADSORBING') chargeFlux = 0
    call FieldShrinkBox(index)
    call TimerStop('reconcileGrids')
  end subroutine reconcileGrids

  subroutine reconcileGridsUpwards(index)
    integer(ik), intent(in) :: index ! Field to reconcile between the grids
  !
  !  Pass 1 - bottom to top. For each subdivided point in the parent grid, copy
  !           the total value from the corresponding section of the contained
  !           grid.
  !
    integer(ik)                :: isrc           ! Grid index for the source grid
    integer(ik)                :: idst           ! Grid index for the destination grid
    type(SimpleGridT), pointer :: src            ! Source grid
    type(SimpleGridT), pointer :: dst            ! Destination grid
    integer(ik)                :: ix1, iy1, iz1  ! First and last indices within the contained grid
    integer(ik)                :: ix2, iy2, iz2  !
    integer(ik)                :: ix, iy, iz     ! Indices within the contained grid
    integer(ik)                :: ixdl, ixdu     ! Range of affected indices within the containing grid (X)
    integer(ik)                :: iydl, iydu     ! Range of affected indices within the containing grid (Y)
    integer(ik)                :: izdl, izdu     ! Range of affected indices within the containing grid (Z)
    integer(ik)                :: ixd, iyd, izd  ! Indices within the containing grid
    integer(ik)                :: npx, npy, npz  ! Number of points in each direction along the fine grid
    integer(ik)                :: npts           ! Total number of points in a section
    complex(rk)                :: sum_field      ! Temporary for Cartesian stitching
    logical                    :: stitchPolar    ! True if norm-conserving stitching is needed
    real(rk)                   :: field_abs      ! Temporaries for polar stitching
    real(rk)                   :: field_arg
    real(rk)                   :: rho
    complex(rk)                :: field_pol
    !
    call TimerStart('reconcileGridsUpwards')
    stitchPolar = grid%wavefunction(index) .and. allowPolar
    !
    loop_level: do isrc=grid%ngrids,1,-1
      src  => grid%grids(isrc)
      idst = src%up_grid
      if (idst==0) cycle ! No parent grid
      dst  => grid%grids(idst)
      ixdl = src%up_index(1,1) ; ixdu = src%up_index(2,1)
      iydl = src%up_index(1,2) ; iydu = src%up_index(2,2)
      izdl = src%up_index(1,3) ; izdu = src%up_index(2,3)
      call FieldAugmentBox(index,idst,ixdl,ixdu,iydl,iydu,izdl,izdu)
      !$omp parallel do private(ixd,iyd,izd,iz1,iz2,npz,iy1,iy2,npy,ix1,ix2,npx,sum_field) &
      !$omp&            private(npts,field_pol,field_abs,field_arg,ix,iy,iz,rho)
      loop_z: do izd=izdl,izdu
        !
        !  Set average for this plane to zero, in case the underlying box
        !  is implicitly zero, and this segment gets skipped.
        !
        dst%fields(ixdl:ixdu,iydl:iydu,izd,index) = 0
        !
        iz1 = src%down_index(1,izd,3)
        iz2 = src%down_index(2,izd,3)
        npz = iz2 - iz1 + 1
        iz1 = max(iz1,src%active(1,3,index))
        iz2 = min(iz2,src%active(2,3,index))
        if (iz2<iz1) cycle loop_z ! There are no non-zero points in this slice
        !
        loop_y: do iyd=iydl,iydu
          iy1 = src%down_index(1,iyd,2)
          iy2 = src%down_index(2,iyd,2)
          npy = iy2 - iy1 + 1
          iy1 = max(iy1,src%active(1,2,index))
          iy2 = min(iy2,src%active(2,2,index))
          if (iy2<iy1) cycle loop_y ! There are no non-zero points in this column
          !
          loop_x: do ixd=ixdl,ixdu
            ix1 = src%down_index(1,ixd,1)
            ix2 = src%down_index(2,ixd,1)
            npx = ix2 - ix1 + 1
            ix1 = max(ix1,src%active(1,1,index))
            ix2 = min(ix2,src%active(2,1,index))
            if (ix2<ix1) cycle loop_x ! There are no non-zero points in this cell
            npts = npx*npy*npz
            !
            !  We have two possible ways of stitching our data: using Cartesian
            !  or polar complex numbers. Cartesian is faster, but polar preserves
            !  the total particle count, and should be more meaningful physically
            !  FOR THE WAVEFUNCTION FIELDS.
            !
            if (stitchPolar) then
              field_abs = 0 ; field_arg = 0
              inner_z: do iz=iz1,iz2
                inner_y: do iy=iy1,iy2
                  inner_x: do ix=ix1,ix2
                    field_pol = cmplx_c2p(src%fields(ix,iy,iz,index))
                    rho       = real(field_pol,kind=rk)**2
                    field_abs = field_abs + rho
                    field_arg = field_arg + rho*aimag(field_pol)
                  end do inner_x
                end do inner_y
              end do inner_z
              !
              !  At this point, field_abs contains total density in this element,
              !  while field_arg is density-weighted phase average
              !
              if (field_abs/=0) field_arg = field_arg/field_abs
              field_abs = sqrt(field_abs/npts)
              sum_field = cmplx_p2c(cmplx(field_abs,field_arg,kind=rk))
            else
              sum_field  = sum( src%fields(ix1:ix2,iy1:iy2,iz1:iz2,index) )/npts
            end if  ! stitchPolar
            ! if (dst%down_grid(ixd,iyd,izd)/=isrc) then
            !   stop 'reconcileGridsUpwards - wrong grid'
            ! end if
            dst%fields(ixd,iyd,izd,index) = sum_field
            !
          end do loop_x
        end do loop_y
      end do loop_z
      !$end omp parallel do
    end do loop_level
    !
    call TimerStop('reconcileGridsUpwards')

  end subroutine reconcileGridsUpwards

  subroutine reconcileGridsDownwards(index)
    integer(ik), intent(in) :: index ! Field to reconcile between the grids
  !
  !  Pass 2 - top to bottom. For each contained grid, calculate the buffer
  !           region by interpolating the values from the coarser grid
  !
    integer(ik)                :: isrc            ! Grid index for the source grid
    integer(ik)                :: idst            ! Grid index for the destination grid
    type(SimpleGridT), pointer :: src             ! Source grid
    type(SimpleGridT), pointer :: dst             ! Destination grid
    integer(ik)                :: sxl, sxu        ! Cube withing the outer grid, containing
    integer(ik)                :: syl, syu        ! the inner grid. This is coming from up_index
    integer(ik)                :: szl, szu
    integer(ik)                :: sx, sy, sz      ! Outer cube we are looking at now
    type(InterpolationT)       :: it              ! Interpolation function
    logical                    :: stitchPolar     ! Use norm-conserving routines
    !
    call TimerStart('reconcileGridsDownwards')
    stitchPolar = grid%wavefunction(index) .and. allowPolar
    !
    !  Shrinking the box is necessary at this point, in order to know whether
    !  we can avoid reconciling some of the box boundaries altogether.
    !
    call FieldShrinkBox(index)
    !
    grid_level: do idst=1,grid%ngrids
      dst  => grid%grids(idst)
      isrc = dst%up_grid
      if (isrc==0) cycle grid_level ! No parent grid
      src  => grid%grids(isrc)
      !
      !  Calculate values for the padding elements. This is a somewhat tricky
      !  business, as we must ensure smooth continuation, by taking some reference
      !  points from the fine grid, and some from the coarse grid. Because we are
      !  lazy, we'll cheat, and do the following:
      !  1. For each outer cube, containing a border point, we'll calculate 4th
      !     order interpolating polynomial, by fitting to all the neighbour points
      !     (this relies on doing bottom-up propagation first)
      !  2. For each border point, we then find the nearest inner grid point, and
      !     adjust the interpolating polynomial's constant term to match field
      !     value in this point
      !  3. We calculate function value at border point from this adjusted fit
      !
      sxl = dst%up_index(1,1) ; sxu = dst%up_index(2,1)
      syl = dst%up_index(1,2) ; syu = dst%up_index(2,2)
      szl = dst%up_index(1,3) ; szu = dst%up_index(2,3)
      !
      !  Possible shortcut: if the non-zero part of the containing box is smaller
      !  than the footprint we just computed, minus the border area (2 pels), the
      !  boundary is known to be zero, and does not need to be reconciled. Otherwise,
      !  we must bite the bullet, augment the containing grid up to the boundary,
      !  and calculate the interpolation region explicitly.
      !
      if ( all( src%active(1,:,index) > (/ sxl, syl, szl /) + 2 ) .and. &
           all( src%active(2,:,index) < (/ sxu, syu, szu /) - 2 ) ) cycle
      call FieldAugmentBox(index,isrc,sxl-2,sxu+2,syl-2,syu+2,szl-2,szl+2)
      call FieldAugmentBox(index,idst,dst%active(1,1,0),dst%active(2,1,0), &
                                      dst%active(1,2,0),dst%active(2,2,0), &
                                      dst%active(1,3,0),dst%active(2,3,0))
      !
      !  Prepare for interpolation
      !
      call prepareInterpolation(src,sxl,syl,szl)
      !
      !$omp parallel do private(it,sx,sy,sz)
      do sz=szl-grid%nborder,szu+grid%nborder
        do sy=syl-grid%nborder,syu+grid%nborder
          do sx=sxl-grid%nborder,sxu+grid%nborder
            if (sx>=sxl .and. sx<=sxu .and. &
                sy>=syl .and. sy<=syu .and. &
                sz>=szl .and. sz<=szu) cycle
            if (verbose>=4) write (out,"(' Outer cube : ',3i4)") sx, sy, sz
            call interpolateFast (src,sx,sy,sz,index,it,stitchPolar)
            call stitchInner     (dst,sx,sy,sz,index,it,stitchPolar)
          end do
        end do
      end do
      !$omp end parallel do
      !
      !  Cleanup after the interpolation routines
      !
      call endInterpolation

    end do grid_level
    !
    call TimerStop('reconcileGridsDownwards')
  end subroutine reconcileGridsDownwards

  !
  !  Extrapolate the outer grid one level out, so that there is a smooth connection
  !  with the inside. In this case, the wavefunction will continue through the wall,
  !  as if it does not exist. We also can calculate the total leakage on the way
  !
  !  Strictly speaking, the wall should be implemented using Fourier transform
  !  close to the box boundary. However, through the first order, it is OK to
  !  simply use linear interpolation.
  !
  subroutine reconcileOuterWall(index)
    integer(ik), intent(in)    :: index           ! Field to work on

    type(SimpleGridT), pointer :: gr              ! Outermost grid
    integer(ik)                :: xl, yl, zl      ! Lower limit of the box, including border
    integer(ik)                :: xu, yu, zu      ! Upper limit of the box, INcluding border
    integer(ik)                :: xb, yb, zb      ! Upper limit of the box, EXcluding border
    integer(ik)                :: sx, sy, sz      ! Point we are extrapolating to
    integer(ik)                :: cx, cy, cz      ! Centre of the cube we are extrapolating from
    integer(ik)                :: ixl, iyl, izl   ! Coordinates of the allowed center points,
    integer(ik)                :: ixu, iyu, izu   ! for the extrapolation cube
    type(InterpolationT)       :: it              ! Interpolation function
    complex(rk)                :: exf
    !
    call TimerStart('reconcileOuterWall')
    if (wallType/='ADSORBING') stop 'reconcileOuterWall - bad arguments'
    !
    gr  => grid%grids(1)
    !
    !  Limits of the box
    !
    xb = gr%npoints(1) ; xl = 1 - grid%nborder ; xu = xb + grid%nborder
    yb = gr%npoints(2) ; yl = 1 - grid%nborder ; yu = yb + grid%nborder
    zb = gr%npoints(3) ; zl = 1 - grid%nborder ; zu = zb + grid%nborder
    !
    !  Limits of the allowed centre points for the extrapolation
    !  "1" is the maximum displacement from the central point, used
    !  in extrapolation.
    !
    ixl = 1 + 1 ; ixu = xb - 1
    iyl = 1 + 1 ; iyu = yb - 1
    izl = 1 + 1 ; izu = zb - 1
    if (ixl>ixu .or. iyl>iyu .or. izl>izu) then
      write (out,"(' reconcileOuterWall - grid too small: ',6i5)") ixl, ixu, iyl, iyu, izl, izu
      stop 'reconcileOuterWall - outer grid too small'
    end if
    !
    !  Possible shortcut - if the wavefunction had not reached the outer wall yet,
    !  we can bail out without doing any work here. Otherwise, we'll have to construct
    !  the full box, and keep going.
    !
    if ( all(gr%active(1,:,index) > (/ ixl, iyl, izl /) + 1 ) .and. &
         all(gr%active(2,:,index) < (/ ixu, iyu, izu /) - 1 ) ) then
      call TimerStop('reconcileOuterWall')
      return
    end if
    !
    call FieldAugmentBox(index,1,gr%active(1,1,0),gr%active(2,1,0), &
                                 gr%active(1,2,0),gr%active(2,2,0), &
                                 gr%active(1,3,0),gr%active(2,3,0))
    !
    !  Prepare for extrapolation, using full cube
    !
    call prepareExtrapolation(gr,ixl,iyl,izl)
    !
    !  Go over ther border points; for each point, calculate the interpolating
    !  polynomials, and extrapolate.
    !
    !$omp parallel do private(sz,sy,sx,cx,cy,cz,it,exf)
    do sz=zl,zu
      do sy=yl,yu
        do sx=xl,xu
          if (sx>=1 .and. sx<=xb .and. sy>=1 .and. sy<=yb .and. sz>=1 .and. sz<=zb) cycle
          !
          !  Figure the nearest inner point. To do this, it's enough to simply
          !  fold the point into the interior.
          !
          cx = min(max(ixl,sx),ixu)
          cy = min(max(iyl,sy),iyu)
          cz = min(max(izl,sz),izu)
          call extrapolateFast(gr,cx,cy,cz,index,it)
          exf = interp(it,gr%coords(1:3,sx,sy,sz))
          gr%fields(sx,sy,sz,index) = exf
          !
          if (verbose>=4) then
            write (out,"(' Extrapolate pt. ',3i6,' from cube at ',3i6,"// &
                       "' value = ',2f14.8,' cube centre = ',2f14.8)") &
                   sx, sy, sz, cx, cy, cz, exf, gr%fields(cx,cy,cz,index)
          end if
        end do
      end do
    end do
    !$omp end parallel do
    !
    !  Cleanup after the interpolation routines
    !
    call endInterpolation
    !
    call TimerStop('reconcileOuterWall')

  end subroutine reconcileOuterWall

  !
  !  Estimate charge flux from the box, by integrating particle flux at the
  !  box boundary
  !
  subroutine calculateChargeFlux(in)
    integer(ik), intent(in)    :: in         ! Field to work on

    type(SimpleGridT), pointer :: gr         ! Outermost grid
    real(rk)                   :: flux(2,3)  ! Flux for the low (1) or high (2) boundary
                                             ! of the dimension (1,2,3)
    integer(ik)                :: nx,ny,nz   ! Upper range of the grid
    !
    call TimerStart('calculateChargeFlux')
    if (verbose>=2) write (out,"(' calculateChargeFlux - ',i8)") in
    !
    gr  => grid%grids(1)
    chargeFlux = 0
    flux       = 0
    nx         = gr%npoints(1)
    ny         = gr%npoints(2)
    nz         = gr%npoints(3)
    !
    !  Possible shortcut - if the wavefunction had not reached the outer wall yet,
    !  we know that the flux is zero.
    !
    if ( all(gr%active(1,:,in) > (/ 2,  2,  2  /)     ) .and. &
         all(gr%active(2,:,in) < (/ nx, ny, nz /) - 1 ) ) then
      chargeFlux = 0
      call TimerStop('calculateChargeFlux')
      return
    end if
    !
    call FieldAugmentBox(in,1,gr%active(1,1,0),gr%active(2,1,0), &
                              gr%active(1,2,0),gr%active(2,2,0), &
                              gr%active(1,3,0),gr%active(2,3,0))
    !
    !  Particle current density is given by:
    !
    !   hbar       *
    !   ---- Im Psi  \Nabla Psi
    !    m
    !
    !  This current density must be integrated over each of the faces of the
    !  cube at the grid boundary. First, we'll take the integral without
    !  figuring in surface patch area, or the numerical differentiation step.
    !
    !  Fluxes are defined as positive if they flow towards outside of the box,
    !  so that we have to reverse signs for the low faces in each direction.
    !
    !                              Psi                            \Nabla Psi
    !$omp parallel
    !$omp sections
    !$omp section
      flux(1,1) = -sum(aimag( conjg(gr%fields(   1,1:ny,1:nz,in))*(gr%fields(   2,1:ny,1:nz,in)-gr%fields(   0,1:ny,1:nz,in)) ))
    !$omp section
      flux(2,1) =  sum(aimag( conjg(gr%fields(  nx,1:ny,1:nz,in))*(gr%fields(nx+1,1:ny,1:nz,in)-gr%fields(nx-1,1:ny,1:nz,in)) ))
    !$omp section
      flux(1,2) = -sum(aimag( conjg(gr%fields(1:nx,   1,1:nz,in))*(gr%fields(1:nx,   2,1:nz,in)-gr%fields(1:nx,   0,1:nz,in)) ))
    !$omp section
      flux(2,2) =  sum(aimag( conjg(gr%fields(1:nx,  ny,1:nz,in))*(gr%fields(1:nx,ny+1,1:nz,in)-gr%fields(1:nx,ny-1,1:nz,in)) ))
    !$omp section
      flux(1,3) = -sum(aimag( conjg(gr%fields(1:nx,1:ny,   1,in))*(gr%fields(1:nx,1:ny,   2,in)-gr%fields(1:nx,1:ny,   0,in)) ))
    !$omp section
      flux(2,3) =  sum(aimag( conjg(gr%fields(1:nx,1:ny,  nz,in))*(gr%fields(1:nx,1:ny,nz+1,in)-gr%fields(1:nx,1:ny,nz-1,in)) ))
    !$omp end sections
    !$omp end parallel
    !
    !  Fluxes must be normalized by multiplying them with the surface patch area,
    !  and dividing with the numerical differentiation step (which is twice the
    !  grid spacing in the that direction).
    !
    flux(:,:) = flux(:,:) * spread(gr%weight/(2.0_rk*gr%step(:)**2),1,2)
    !
    if (verbose>=2) then
      write (out,"('      raw particle fluxes X: ',2f12.8,' Y: ',2f12.8,' Z: ',2f12.8)") flux
    end if
    !
    !  Finally, charge flux is the sum of fluxes
    !
    chargeFlux = sum(flux)
    call TimerStop('calculateChargeFlux')
  end subroutine calculateChargeFlux

  subroutine prepareInterpolation(gr,sx,sy,sz)
    type(SimpleGridT), intent(in)     :: gr
    integer(ik), intent(in)           :: sx, sy, sz
    !
    integer(ik) :: ix, iy, iz, ipt
    integer(ik) :: px, py, pz
    real(rk)    :: coords (3,125)
    real(rk)    :: weights(  125)
    !
    if (sx<=1 .or. sx>=gr%npoints(1) .or. &
        sy<=1 .or. sy>=gr%npoints(2) .or. &
        sz<=1 .or. sz>=gr%npoints(3) ) then
      write (out,"(' Intepolation for Martian points ',3i4,' requested')") sx, sy, sz
      stop 'prepareInterpolation - Martian points'
    end if

    ipt = 0
    do iz=-2,2
      do iy=-2,2
        do ix=-2,2
          if (abs(ix)+abs(iy)+abs(iz)>3) cycle
          px = sx+ix ; py = sy+iy ; pz = sz+iz
          ipt = ipt + 1
          coords (:,ipt) = gr%coords(1:3,px,py,pz)
          weights(  ipt) = 1.0_rk/real(1+abs(ix)+abs(iy)+abs(iz),kind=rk)
        end do
      end do
    end do
    call beginInterpolation(interpolationOrder,ipt,coords,weights)
  end subroutine prepareInterpolation

  subroutine interpolateFast(gr,sx,sy,sz,index,it,polar)
    type(SimpleGridT), intent(in)     :: gr
    integer(ik), intent(in)           :: sx, sy, sz
    integer(ik), intent(in)           :: index
    type(InterpolationT), intent(out) :: it
    logical, intent(in)               :: polar

    integer(ik) :: ix, iy, iz, ipt
    integer(ik) :: px, py, pz
    complex(rk) :: values (  125)

    if (sx<=1 .or. sx>=gr%npoints(1) .or. &
        sy<=1 .or. sy>=gr%npoints(2) .or. &
        sz<=1 .or. sz>=gr%npoints(3) ) then
      write (out,"(' Intepolation for Martian points ',3i4,' requested')") sx, sy, sz
      stop 'interpolateFast - Martian points'
    end if

    ipt = 0
    do iz=-2,2
      do iy=-2,2
        do ix=-2,2
          if (abs(ix)+abs(iy)+abs(iz)>3) cycle
          px = sx+ix ; py = sy+iy ; pz = sz+iz
          ipt = ipt + 1
          !
          !  For the norm-conserving version, we cheat by packing (r,arg)
          !  pair into a standard complex type. This works, because our
          !  interpolation abscissas and weights are real, and don't mix
          !  the two components
          !
          if (polar) then
            values(ipt) = cmplx_c2p(gr%fields(px,py,pz,index))
          else
            values(ipt) = gr%fields(px,py,pz,index)
          end if
        end do
      end do
    end do
    call fastInterpolation(gr%coords(1:3,sx,sy,sz),values(1:ipt),it)
  end subroutine interpolateFast

  !
  !  Extrapolation does not have to be very accurate (first and second derivatives
  !  are sufficient). However, it must be be very, very stable! For this reason, we
  !  can't use the same routine we use for stitching grid boundaries.
  !
  subroutine prepareExtrapolation(gr,sx,sy,sz)
    type(SimpleGridT), intent(in)     :: gr
    integer(ik), intent(in)           :: sx, sy, sz

    integer(ik) :: ipt
    real(rk)    :: coords (3,125)
    real(rk)    :: weights(  125)

    if (sx<=1 .or. sx>=gr%npoints(1) .or. &
        sy<=1 .or. sy>=gr%npoints(2) .or. &
        sz<=1 .or. sz>=gr%npoints(3) ) then
      write (out,"(' Intepolation for Martian points ',3i4,' requested')") sx, sy, sz
      stop 'prepareInterpolation - Martian points'
    end if

    ipt = 3**3
    coords (:,1:ipt) = reshape( gr%coords(1:3,sx-1:sx+1,sy-1:sy+1,sz-1:sz+1), (/ 3, ipt /) )
    weights(1:ipt)   = 1.0_rk
    call beginInterpolation(extrapolationOrder,ipt,coords,weights)
  end subroutine prepareExtrapolation

  subroutine extrapolateFast(gr,sx,sy,sz,index,it)
    type(SimpleGridT), intent(in)     :: gr
    integer(ik), intent(in)           :: sx, sy, sz
    integer(ik), intent(in)           :: index
    type(InterpolationT), intent(out) :: it

    integer(ik) :: ipt
    complex(rk) :: values (  125)

    if (sx<=1 .or. sx>=gr%npoints(1) .or. &
        sy<=1 .or. sy>=gr%npoints(2) .or. &
        sz<=1 .or. sz>=gr%npoints(3) ) then
      write (out,"(' Intepolation for Martian points ',3i4,' requested')") sx, sy, sz
      stop 'interpolateFast - Martian points'
    end if

    ipt = 3**3
    values (1:ipt) = reshape( gr%fields(sx-1:sx+1,sy-1:sy+1,sz-1:sz+1,index), (/ ipt /) )
    call fastInterpolation(gr%coords(1:3,sx,sy,sz),values(1:ipt),it)
  end subroutine extrapolateFast

  !
  ! Calculates values of all border grid points within the parent cube (sx,sy,sz)
  !
  subroutine stitchInner(gr,sx,sy,sz,index,it,polar)
    type(SimpleGridT), intent(inout) :: gr         ! Grid containing the target points
    integer(ik), intent(in)          :: sx, sy, sz ! Cube of the -parent- grid
    integer(ik), intent(in)          :: index      ! Field to work with
    type(InterpolationT)             :: it         ! Interpolating object
    logical, intent(in)              :: polar      ! Requests norm-conserving version

    integer(ik) :: xl, xu, yl, yu, zl, zu       ! Border points, corresponding to cube (sx,sy,sz)
    integer(ik) :: ix, iy, iz                   ! Border point
    integer(ik) :: iix, iiy, iiz                ! Inner point, closest to the current border point
    complex(rk) :: ref_f                        ! Actual function value at nearest interior point
    complex(rk) :: int_f                        ! Interpolated function value at the same point
    complex(rk) :: val_f                        ! Interpolated value

    xl = gr%down_index(1,sx,1) ; xu = gr%down_index(2,sx,1) ;
    yl = gr%down_index(1,sy,2) ; yu = gr%down_index(2,sy,2) ;
    zl = gr%down_index(1,sz,3) ; zu = gr%down_index(2,sz,3) ;
    do ix=xl,xu
      do iy=yl,yu
        do iz=zl,zu
          if ( ix>=1 .and. ix<=gr%npoints(1) .and. &
               iy>=1 .and. iy<=gr%npoints(2) .and. &
               iz>=1 .and. iz<=gr%npoints(3) ) cycle ! Should not occur, actually?
          if (verbose>=4) write (out,"(' Interpolating at inner point: ',3i4)") ix, iy, iz
          !
          !  Figure the nearest inner point. To do this, it's enough to simply
          !  fold the point into the interior.
          !
          iix = min(max(1,ix),gr%npoints(1))
          iiy = min(max(1,iy),gr%npoints(2))
          iiz = min(max(1,iz),gr%npoints(3))
          if (verbose>=4) write (out,"(' Nearest internal point is: ',3i4)") iix, iiy, iiz
          int_f = interp(it,gr%coords(1:3,iix,iiy,iiz))
          ref_f = gr%fields(iix,iiy,iiz,index)
          if (polar) ref_f = cmplx_c2p(ref_f)
          if (verbose>=3) then
            write (out,"(3(1x,i4),' -> ',3(1x,i4))") ix, iy, iz, iix, iiy, iiz
            write (out,"(' At position ',3f8.3,' int = ',2f12.6,' ref = ',2f12.6,' diff = ',2f12.6)") &
                   gr%coords(1:3,iix,iiy,iiz), int_f, ref_f, int_f - ref_f
          end if
          val_f = interp(it,gr%coords(1:3,ix,iy,iz)) + ref_f - int_f
          if (polar) val_f = cmplx_p2c(val_f)
          gr%fields(ix,iy,iz,index) = val_f
        end do
      end do
    end do
  end subroutine stitchInner

end module multigrid
