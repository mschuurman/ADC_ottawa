!
!  Michael, some comments/suggections:
!   a) A more elegant handling of the field parameters would be very useful.
!      The goal should be to let the user to perform a sensible range
!      of calculations without having to recompile the program.
!   b) Document all global variables and function arguments. Since the s_*
!      notation does not make much sense anymore, it would be nice if you
!      could clean this up.
!   c) Where possible/necessary, also document local variables
!   d) Add explicit handling of error conditions (open, allocate, close, etc),
!      issuing a meanigful error message
!   e) Wherever possible, avoid passing parameters to functions through
!      global parameters - this prevents OpenMP parallelization.
!      Whenever a function/subroutine depends on non-constant data
!      passed through global parameters - document it in the routine 
!      itself. E.g. "s_fdt" is OK as a global variable (it is write-once), 
!      but s_fCx1, etc are probably not.
!   f) Give subroutines and functions meaningful names, which would be
!      unlikely to clash with other routines (e.g. "A" is not a good name
!      for anything but a subroutine-local variable).
!   f) All names to all non-trivial loops
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  eikpoly3D.cpp - Code to calculate the eikonal-Volkov states
!!                  using polynomial interpolation, in 3D.
!!
!!  Michael Spanner - 2008
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module eikonal_tools_eva
  use accuracy
  use phase_interpolation_eva
  use multigrid
  use fields
  use timer
  implicit none
  !
  private
  public EVAdataT
  public EVA_Initialize, EVA_Advance, EVA_Destroy
  !
  !  Michael, as a general rule: all non-constant module-global variables which must preserve
  !  their values across invocations of module procedures must have "save" attribute. The
  !  behavior of module variables without "save" is undefined once the the last active 
  !  module procedure returns.
  !
  integer(ik), parameter   :: io_unit = 46    ! Unit number to use for various files
  !
  !  We would like to allow several eikonal-Volkov solutions to be propagated 
  !  simultaneously, if necessary. To make this relatively easy, we collect
  !  all persistent state for the trajectorues in the EVAdataT data structure.
  !  Although the state can be examined from the calling modules, modifying
  !  it is NOT A GOOD IDEA.
  !
  !  When any of the subroutines from eikonal_tools_eva are active, the current
  !  trajectory descriptor is kept in the "tr" pointer. 
  !
  type EVAdataT
    !
    !  Input parameters: These values have to be initialized prior to EVA_Initialize call. 
    !                    They should not be modified after the call.
    !
    !  Program output and debugging
    !
    integer(ik)        :: verbose = 2                        ! Overall verbosity level
    logical            :: dump_trajectory = .true.           ! Save human-readable trajectory
    character(len=200) :: trajectory_file = 'trajectory.dat' ! Name of the file for the trajectory
    logical            :: dump_returns    = .true.           ! Save human-readable return time data
    character(len=200) :: returns_file    = 'returns.dat'    ! Name of the file for the returns table
    !
    !  General parameters for the simulation, defaults are likely OK
    !
    real(rk)           :: max_timestep       = 0.50_rk       ! Largest allowed time step for the trajectory
    integer(ik)        :: s_nScaleEdgePoints = 5             ! Factor to scale down the density of the edge points grid
    integer(ik)        :: time_subdivision   = 10            ! Additional subdivision of the time step used for
                                                             ! taking some of the line integrals
    real(rk)           :: derivative_step    = 1e-2_rk       ! Fraction of s_fdt time step used in finite-difference
                                                             ! calculation of the derivatives
    integer(ik)        :: io_edge            = 54            ! IO unit to use for the edge points cache
    !
    !  System-specific parameters, no sensible defaults
    !  Michael, perhaps other field parameters should be here, too?
    !  See evaluateExternalField() in cylindrical/h2p_recollide.f90 for an example
    !
    real(rk)           :: s_fK(3)                            ! Asymptotic field-free momentum of the electron
    real(rk)           :: s_fF                               ! Laser electric field
    real(rk)           :: s_fw                               ! Laser frequency
    real(rk)           :: s_fAn(3)                           ! Vector potential direction
    real(rk)           :: s_fMinT                            ! Starting time for propagation
    real(rk)           :: s_fMaxT                            ! Final time for propagation
    !
    !  ====================================================================
    !  All variables below this point are maintained by the module routines
    !
    integer(ik)        :: s_nNumT          ! Number of time steps
    integer(ik)        :: s_nt             ! Present time step. The valid range is 1 to s_nNumT
    !
    real(rk)           :: s_fAmax          ! Maximum value of Vector potential
    real(rk)           :: s_fOsc           ! Free electron oscillation radius
    !
    real(rk)           :: s_pfGrid(3,-1:1) ! Grid edges, 1:3=x,y,z, -1:1=min,undefined,max
                                           ! Note that the EVA phase grid does not include the
                                           ! first and the last layers along each direction -
                                           ! so it is not the same as the definition used by
                                           ! multigrid.f90
    integer(ik)        :: s_nNum(3)        ! Number of grid points
    real(rk)           :: s_fd(3)          ! Grid spacing
    real(rk)           :: s_fdt            ! Estimated minimum time step
    !
    !  Arrays
    !
    real(rk), pointer  :: s_pfX(:)         ! X coordinates of the grid points
    real(rk), pointer  :: s_pfY(:)         ! Y coordinates of the grid points
    real(rk), pointer  :: s_pfZ(:)         ! Z coordinates of the grid points
    !
    real(rk), pointer  :: s_pfRt(:,:)      ! Trajectory position for a given time step
    real(rk), pointer  :: s_pft(:)         ! Time value at a given time step
    !
    real(rk)           :: s_pfWG(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)  ! Stencil for interpolation of G
    real(rk)           :: s_pfWU(-1:1,-1:1,-1:1)                 ! Stencil for integration over U
    !
    !  Edge point for the returning trajectories
    !
    complex(rk), pointer :: edge_xy(:,:,:)       ! (nx,ny,6) - lower / upper edges of the grid
    complex(rk), pointer :: edge_xz(:,:,:)       ! (nx,6,nz) - ditto
    complex(rk), pointer :: edge_yz(:,:,:)       ! (6,ny,nz) - ...
    complex(rk), pointer :: interp_xy(:,:,:)     ! (nx,ny,1) - Integrpolated data for the edge points
    complex(rk), pointer :: interp_xz(:,:,:)     ! (nx,1,nz)
    complex(rk), pointer :: interp_yz(:,:,:)     ! (1,ny,nz)
    !
    !  Edge of grid interpolation stuff
    !
    integer(ik)          :: s_nNumE(3)           ! Number of points in the reduced edge grid
    real(rk)             :: s_fdE(3)             ! Spacing of the reduced edge grid
    !
    real(rk), pointer    :: s_pfXe(:)            ! 
    real(rk), pointer    :: s_pfYe(:)            ! Sparse edge points at which to calculate full integral
    real(rk), pointer    :: s_pfZe(:)            !
    !
    real(rk), pointer    :: s_pfGxy(:,:)         !
    real(rk), pointer    :: s_pfGyz(:,:)         ! Grid of 'sparse' full integrals
    real(rk), pointer    :: s_pfGxz(:,:)         !
    !
    integer(ik), pointer :: s_pnReturnToGrid(:)  ! Time step at which the trajectory returns to the grid
    integer(ik), pointer :: s_pbSaveThisStep(:)  ! Number of times a trajectory returns to the grid at this 
                                                 ! time step. If non-zero, save the edge wavefunction
    integer(ik)          :: next_edge_record     ! Next available edge record
    integer(ik), pointer :: edge_record(:)       ! Record number in the (io_edge) file containg the data
    !
  end type
  !
  type(EVAdataT), pointer :: tr  ! Persistent state for the currently-active trajectory
  !
  contains
  !
  !  Initialize persistent data for an eikonal-Volkov wavefunction
  !
  subroutine EVA_Initialize(trajectory)
    type(EVAdataT), target, intent(inout) :: trajectory ! Persistent trajectory data
    !
    !
    call TimerStart('EVA Initialize')
    !
    !  Set global state pointer
    !
    tr => trajectory
    ! 
    !  Field parameters
    !
    tr%s_fAmax = tr%s_fF/tr%s_fw                    ! Vector potential
    tr%s_fOsc  = tr%s_fAmax/tr%s_fw                 ! Free electron oscillations
    tr%s_fAn   = tr%s_fAn / sqrt(sum(tr%s_fAn**2))  ! Normalize direction to unity
    !
    !  Query relevant grid parameters - grids must be already initialized
    !  at this point. Eikonal-Volkov routines can't handle multigrids, so
    !  we'll be alwyas asking for the top-level grid.
    !
    !  We have two kinds of data on grid: the potential, and the phases.
    !  The potential is available for the entire grid, with the range of
    !  1:tr%s_nNum. The phases are propagated for grid points in the
    !  2:tr%s_nNum-1 range, while the the first and the last layers along
    !  each dimension remain undefined.
    !
    tr%s_nNum = FieldGridNPoints(1_ik)
    tr%s_fd   = FieldGridSpacing(1_ik)
    !
    !  Arrays
    !
    allocate(tr%s_pfX(tr%s_nNum(1)))                      ! X grid points
    allocate(tr%s_pfY(tr%s_nNum(2)))                      ! Y grid points
    allocate(tr%s_pfZ(tr%s_nNum(3)))                      ! Z grid points
    !  
    call FieldGridCoordinates(1_ik,tr%s_pfX,1_ik)
    call FieldGridCoordinates(2_ik,tr%s_pfY,1_ik)
    call FieldGridCoordinates(3_ik,tr%s_pfZ,1_ik)
    !
    !  As far as the definition of the EVA grid, the first and the last layers 
    !  do not exist. 
    !
    tr%s_pfGrid(:,-1) = (/ tr%s_pfX(2),              tr%s_pfY(2),            tr%s_pfZ(2)                /)
    tr%s_pfGrid(:, 1) = (/ tr%s_pfX(tr%s_nNum(1)-1), tr%s_pfY(tr%s_nNum(2)-1), tr%s_pfZ(tr%s_nNum(3)-1) /)
    !
    !  Buffers for the grid-edge points.
    !  Interpolation requires 3 layers at each extreme (2:4 and the nNum-3:nNum-1)
    !  Interpolation of the boundary trajectories needs just one layer on at most one side.
    !
    allocate (tr%edge_xy(tr%s_nNum(1),tr%s_nNum(2),           6))
    allocate (tr%edge_xz(tr%s_nNum(1),           6,tr%s_nNum(3)))
    allocate (tr%edge_yz(           6,tr%s_nNum(2),tr%s_nNum(3)))
    allocate (tr%interp_xy(tr%s_nNum(1),tr%s_nNum(2),           1))
    allocate (tr%interp_xz(tr%s_nNum(1),           1,tr%s_nNum(3)))
    allocate (tr%interp_yz(           1,tr%s_nNum(2),tr%s_nNum(3)))
    !
    !  Coarse grid for edge interpolation. We don't count the two boundary 
    !  points of the fine grid!
    !
    tr%s_nNumE = (tr%s_nNum-2) / tr%s_nScaleEdgePoints
    !
    tr%s_fdE(:) = (tr%s_pfGrid(:,1)-tr%s_pfGrid(:,-1))/(tr%s_nNumE(:)-1)
    !
    allocate(tr%s_pfXe(tr%s_nNumE(1)))
    allocate(tr%s_pfYe(tr%s_nNumE(2)))
    allocate(tr%s_pfZe(tr%s_nNumE(3)))
    !
    allocate(tr%s_pfGxy(tr%s_nNumE(1),tr%s_nNumE(2)))
    allocate(tr%s_pfGyz(tr%s_nNumE(2),tr%s_nNumE(3)))
    allocate(tr%s_pfGxz(tr%s_nNumE(1),tr%s_nNumE(3)))
    !
    call SetSpaceGrids
    !
    if (tr%verbose>0) then
      write (out,"()")
      write (out,"('            Grid parameters for the eikonal-Volkov integration:')")
      write (out,"('                          ',3a15)") '   X   ', '   Y   ', '   Z   '
      write (out,"('          Lower boundary: ',3f15.8)") tr%s_pfGrid(:,-1)
      write (out,"('          Upper boundary: ',3f15.8)") tr%s_pfGrid(:, 1)
      write (out,"('    Interior grid points: ',3i15)") tr%s_nNum
      write (out,"(' Interior grid step size: ',3f15.8)") tr%s_fd
      write (out,"('    Boundary grid points: ',3i15)") tr%s_nNumE
      write (out,"(' Boundary grid step size: ',3f15.8)") tr%s_fdE
      write (out,"()")
    end if
    !
    call SetTimeGrid
    call SetReturnToGrid
    !
    !  Finish initialization of the global data
    !
    tr%s_nt = 1
    call initialize_io
    call TimerStop('EVA Initialize')
  end subroutine EVA_Initialize
  !
  !  Release all resources associated with the trajectory
  !
  subroutine EVA_Destroy(trajectory)
    type(EVAdataT), target, intent(inout) :: trajectory ! Persistent data
    !
    call TimerStart('EVA Destroy')
    !
    tr => trajectory
    !
    close (tr%io_edge)
    !
    if (associated(tr%s_pfX           )) deallocate(tr%s_pfX)
    if (associated(tr%s_pfY           )) deallocate(tr%s_pfY)
    if (associated(tr%s_pfZ           )) deallocate(tr%s_pfZ)
    if (associated(tr%s_pfRt          )) deallocate(tr%s_pfRt)
    if (associated(tr%s_pft           )) deallocate(tr%s_pft)
    if (associated(tr%edge_xy         )) deallocate(tr%edge_xy)
    if (associated(tr%edge_xz         )) deallocate(tr%edge_xz)
    if (associated(tr%edge_yz         )) deallocate(tr%edge_yz)
    if (associated(tr%interp_xy       )) deallocate(tr%interp_xy)
    if (associated(tr%interp_xz       )) deallocate(tr%interp_xz)
    if (associated(tr%interp_yz       )) deallocate(tr%interp_yz)
    if (associated(tr%s_pfXe          )) deallocate(tr%s_pfXe)
    if (associated(tr%s_pfYe          )) deallocate(tr%s_pfYe)
    if (associated(tr%s_pfZe          )) deallocate(tr%s_pfZe)
    if (associated(tr%s_pfGxy         )) deallocate(tr%s_pfGxy)
    if (associated(tr%s_pfGyz         )) deallocate(tr%s_pfGyz)
    if (associated(tr%s_pfGxz         )) deallocate(tr%s_pfGxz)
    if (associated(tr%s_pnReturnToGrid)) deallocate(tr%s_pnReturnToGrid)
    if (associated(tr%s_pbSaveThisStep)) deallocate(tr%s_pbSaveThisStep)
    if (associated(tr%edge_record     )) deallocate(tr%edge_record)
    !
    call TimerStop('EVA Destroy')
  end subroutine EVA_Destroy
  !
  !  Advance eikonal-Volkov trajectory by one time step
  !
  !  Warning: EVA_Advance depends on FLlrePotential() returning a fair approximation
  !           of the long-range tail of the potential in f_core_pot. Don't forget to
  !           initialize FLlrePotential with a suitable call to FLsetLREmultipoles
  !           Furthermore, FLexactphase() should return a fair approximation to the
  !           long-range rail of the adiabatic eikonal phase. So don't forget calls
  !           to FLsetPlanewave() and FLsetMultipolesPhase(), either. Note that the
  !           FLsetPlanewave() should change the sign of the kvec argument: the adiabatic
  !           solutions are for the flat _outgoing_ wavefront, while here we expect the
  !           flat _incoming_ wavefront.
  !
  subroutine EVA_Advance(trajectory,potential,src,dst)
    type(EVAdataT), target, intent(inout) :: trajectory ! Persistent data
    integer(ik), intent(in)               :: potential  ! Input: Data field containing the potential
    integer(ik), intent(in)               :: src        ! Input: Phase at the previous time step
    integer(ik), intent(in)               :: dst        ! Output: Updated phase
    !
    call TimerStart('EVA Advance')
    !
    !  Set up global state pointer
    !
    tr => trajectory
    !
    !  Sanity check
    !
    if (tr%s_nt>=tr%s_nNumT) then
      write (out,"('EVA_Advance: time step ',i8,' is beyond the limit (',i8,')')") &
             tr%s_nt, tr%s_nNumT
      stop 'eikonal_tools_eva%EVA_Advance - trajectory terminated'
    end if
    !
    !  Save current grid if needed for 're-entry' edge trajectories
    !
    if ( tr%s_pbSaveThisStep(tr%s_nt)>0 ) then
      call SaveCurrentEdges(src)
    end if
    !
    !  Load a previous grid if needed for 're-entry' edge trajectories
    !
    if ( tr%s_pnReturnToGrid(tr%s_nt+1) > 0 ) then
      call LoadOldEdges( tr%s_pnReturnToGrid(tr%s_nt+1) )
    end if
    !
    !  Step eikonal function
    !
    call FieldZero(dst)
    call StepPhase(potential,src,dst)
    !
    tr%s_nt = tr%s_nt + 1
    !
    call TimerStop('EVA Advance')
  end subroutine EVA_Advance
  !
  ! ======================================================================================
  !   End of externally visible routines. Everything below is private to the module
  ! ======================================================================================
  !
  !  Open direct-access file for the edge data recall
  !
  subroutine initialize_io
    integer(ik) :: rec_len  ! Record length for the edge i/o file
    !
    inquire (iolength=rec_len) tr%edge_xy, tr%edge_xz, tr%edge_yz
    !
    if (tr%verbose>=1) then
      write (out,"('Record length for the edge I/O file = ',i8)") rec_len
    end if
    !
    open (tr%io_edge,access='DIRECT',status='SCRATCH',recl=rec_len)
    tr%next_edge_record = 1
    tr%edge_record      = 0
    !
  end subroutine initialize_io
  !
  function A( ft ) 
    !
    real(rk) :: A
    real(rk), intent(in) :: ft
    real(rk) :: sfFWHM
    real(rk) :: fsin
    !
    sfFWHM = tr%s_fMaxT*0.5*0.75
    fsin = sin(pi*ft/sfFWHM/2)
    !
    if (ft > sfFWHM*2.0) fsin = 0.0
    A = -tr%s_fAmax*sin(tr%s_fw*ft)*fsin*fsin
    !
  end function A
  !
  !  Calculate Cartesian displacement for the trajectory between two moments of time:
  !
  !    / t2               / t2                                / t2
  !    |   k_L(t') dt' =  |   ( k + A(t') ) dt' = k*(t2-t1) + |   A(t') dt'
  !    / t1               / t1                                / t1
  !
  ! Integrates using Euler method, since it is not important to be fast, with
  ! a step size 10 times smaller than the estimated 'smallest' step size...
  !
  !
  function trajectory_segment(ft1,ft2) result (val)
    real(rk), intent(in) :: ft1, ft2  ! The initial and final times
    real(rk)             :: val(3)    ! Displacement
    !
    real(rk)    :: fInterval
    real(rk)    :: fdtp
    integer(ik) :: ntp
    integer(ik) :: nNumtp
    !
    real(rk) :: fSum(3)
    !
    fInterval = ft2 - ft1
    !
    !  '+1' is to make sure there is at least _some_ steps
    !
    nNumtp = nint( (abs(fInterval) / tr%s_fdt + 1) * tr%time_subdivision ) 
    fdtp = fInterval / nNumtp
    !
    fSum = 0.0
    do ntp=1,nNumtp
      fSum = fSum + A( ft1 + (ntp-1)*fdtp ) * tr%s_fAn
    end do
    !
    val = tr%s_fK*fInterval + fSum*fdtp
    !
  end function trajectory_segment
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SetSpaceGrids
    !
    integer(ik) :: nx, ny, nz
    !
    do nx=1,tr%s_nNumE(1); tr%s_pfXe(nx) = tr%s_pfGrid(1,-1) + (nx-1)*tr%s_fdE(1); end do
    do ny=1,tr%s_nNumE(2); tr%s_pfYe(ny) = tr%s_pfGrid(2,-1) + (ny-1)*tr%s_fdE(1); end do
    do nz=1,tr%s_nNumE(3); tr%s_pfZe(nz) = tr%s_pfGrid(3,-1) + (nz-1)*tr%s_fdE(1); end do
    !
  end subroutine SetSpaceGrids

  subroutine SetTimeGrid
    real(rk)    :: fdt(3)       ! Estimated minimum time step for each direction
    real(rk)    :: fdts         ! 'smallest' time step used while finding best variable time steps 
    integer(ik) :: numTmax      ! Estimated maximum number of time steps
    !
    integer(ik) :: nt
    real(rk)    :: f1(3), ft1
    real(rk)    :: f2(3), ft2
    real(rk)    :: scale
    !
    call TimerStart('EVA SetTimeGrid')
    !
    ! Estimate the minimum time step using the known maximum slope
    ! of the trajectories and the grid spacings.  The goal is to 
    ! have the trajectory segments always smaller than the grid steps.
    ! 
    ! 'scale' determines the target size of the trajectory segments 
    ! as compared to the grid spacings. We have to be a little 
    ! conservative, to make sure that we can grow the trajectory by
    ! a small interval without exiting the cell.
    ! The overall factor of 0.5 is necessary to make sure that the 
    ! trajectory will always terminate within the edge layers, even
    ! if it reverses direction.
    !
    scale = 0.5_rk * ( 1._rk - 1._rk / tr%time_subdivision )
    ! 
    !  Adding spacing(tr%s_fd) to the canonical momentum prevents division by zero.
    !  The estimate relies on vector potential never ever exceeding tr%s_fAmax in
    !  magnitude.
    !
    fdt(:) = scale * tr%s_fd(:) / ( abs(tr%s_fK(:)) + abs(tr%s_fAmax*tr%s_fAn(:)) + spacing(tr%s_fd) )
    !
    tr%s_fdt = min( tr%max_timestep, minval(abs(fdt)) )
    !
    numTmax = int((tr%s_fMaxT - tr%s_fMinT) / tr%s_fdt + 1)
    !
    if (numTmax<=0) then
      write (out,"('eikonal_tools_eva%SetTimeGrid: Estimated size of the time grid is ',i8)") numTmax
      stop 'eikonal_tools_eva%SetTimeGrid - logic failure 1'
    end if
    !
    allocate(tr%s_pft(numTmax)) 
    allocate(tr%s_pfRt(3,numTmax)) 
    !
    ! Set variable time steps and the master trajectory
    !
    tr%s_pft(   1) = tr%s_fMinT
    tr%s_pfRt(:,1) = 0
    !
    nt   = 1
    fdts = tr%s_fdt / tr%time_subdivision
    !
    if (tr%dump_trajectory) then
      open (io_unit,file=tr%trajectory_file,status='replace')
      write (io_unit,"('# MinT = ',f15.7)") tr%s_fMinT
      write (io_unit,"('# MaxT = ',f15.7)") tr%s_fMaxT
      write (io_unit,"('#',4(1x,a15))") ' T ', ' X ', ' Y ', ' Z '
      write (io_unit,"(4(1x,f15.7))") tr%s_pft(1), tr%s_pfRt(:,1)
    end if
    !
    trajectory_step: do
      !
      f1  = 0                ! Position at the beginning of the interval
      ft1 = tr%s_pft(nt)     ! Time at the beginning of the interval
      ft2 = ft1 + tr%s_fdt   ! Time at the end of the interval (to be determined here)
      !
      !  Increase size of current step until either:
      !  a) trajectory segment hits the maximum allowed length
      !  b) maximum allowed time interval is reached
      !  c) the desired simulation time is reached
      !
      trajectory_microstep: do
        f2 = f1 + trajectory_segment(ft1,ft2)
        if ( any(abs(f2)>scale*tr%s_fd) .or. ft2+fdts-tr%s_pft(nt)>=tr%max_timestep .or. ft2>=tr%s_fMaxT ) then
          exit trajectory_microstep
        end if
        f1  = f2
        ft1 = ft2
        ft2 = ft2 + fdts
      end do trajectory_microstep
      !
      !  Save current trajectory and step information for future use
      !
      tr%s_pft   (nt+1) = ft2
      tr%s_pfRt(:,nt+1) = tr%s_pfRt(:,nt) + f2
      !
      !  Sanity check
      !
      if (tr%s_pft(nt+1)-tr%s_pft(nt)>tr%max_timestep .or. &
          any(abs(tr%s_pfRt(:,nt+1)-tr%s_pfRt(:,nt))>tr%s_fd) ) then
        write (out,"(' nt   = ',i6,' t = ',f16.8,' r = ',3f16.8)") nt,   tr%s_pft(nt),   tr%s_pfRt(:,nt)
        write (out,"(' nt+1 = ',i6,' t = ',f16.8,' r = ',3f16.8)") nt+1, tr%s_pft(nt+1), tr%s_pfRt(:,nt+1)
        stop 'eikonal_tools_eva%SetTimeGrid - Maximum allowed step size exceeded'
      end if
      !
      if (tr%dump_trajectory) then
        write (io_unit,"(4(1x,f15.7))") tr%s_pft(nt), tr%s_pfRt(:,nt)
      end if
      !
      nt = nt+1
      if ( tr%s_pft(nt) > tr%s_fMaxT ) exit trajectory_step
      !
      if (nt>=numTmax) then
        write (out,"('SetTimeGrid: Number of time steps (',i6,') exceeds the upper bound (',i6,')')") &
               nt, numTmax
        stop 'eikonal_tools_eva%SetTimeGrid - (impossible?) table overflow'
      end if
    end do trajectory_step
    !
    if (tr%dump_trajectory) then
      close(io_unit)
    end if
    !
    tr%s_nNumT = nt
    !
    write (out,"(/'Estimated maximum number of time points = ',i8)") numTmax
    write (out,"( '     Actual       number of time points = ',i8/)") tr%s_nNumT
    !
    call TimerStop('EVA SetTimeGrid')
  end subroutine SetTimeGrid
  !
  ! Checks if the point (x,y,z) is inside the main grid
  !
  logical function IsInGrid(xyz)
    real(rk), intent(in) :: xyz(3)  ! Coodinates of point to test
    !
    IsInGrid = all(xyz>=tr%s_pfGrid(:,-1)) .and. all(xyz<tr%s_pfGrid(:,1))
  end function IsInGrid
  !
  ! Check for edge trajectories that returns to grid so we know timesteps to save grid data.
  ! It is enough to check if all corner points return to grid...
  !
  subroutine SetReturnToGrid
    integer(ik) :: nt, nt2, ic
    integer(ik) :: pfDir(3)
    real(rk)    :: fc(3,3)  ! Positions where the trajectory has to start to end at
                            ! one of the corners at time step nt
    real(rk)    :: tret
    !
    call TimerStart('EVA SetReturnToGrid')
    !
    allocate(tr%s_pnReturnToGrid(tr%s_nNumT)) 
    allocate(tr%s_pbSaveThisStep(tr%s_nNumT)) 
    allocate(tr%edge_record     (tr%s_nNumT)) 
    !
    tr%s_pnReturnToGrid = 0
    tr%s_pbSaveThisStep = 0
    !
    if (tr%dump_returns) then
      open (io_unit,file=tr%returns_file,status='replace')
      write(io_unit,"('#',4(1x,a12))") ' Time step ', ' Time ', ' Return step ', ' Time '
    end if
    !
    reentry_times: do nt=2,tr%s_nNumT
      !
      ! Find initial (backward) direction of trajectory at current time
      !
      pfDir = nint(sign(1._rk,tr%s_pfRt(:,nt-1)-tr%s_pfRt(:,nt)))
      !
      ! Check relevant corner points
      !
      fc(:,1) = (/ tr%s_pfGrid(1, pfDir(1)), tr%s_pfGrid(2,-pfDir(2)), tr%s_pfGrid(3,-pfDir(3)) /)
      fc(:,2) = (/ tr%s_pfGrid(1,-pfDir(1)), tr%s_pfGrid(2, pfDir(2)), tr%s_pfGrid(3,-pfDir(3)) /)
      fc(:,3) = (/ tr%s_pfGrid(1,-pfDir(1)), tr%s_pfGrid(2,-pfDir(2)), tr%s_pfGrid(3, pfDir(3)) /)
      do ic=1,3 ; fc(:,ic) = fc(:,ic) - tr%s_pfRt(:,nt) ; end do
      departure_times: do nt2=(nt-1),1,-1
        do ic=1,3
          if ( IsInGrid(tr%s_pfRt(:,nt2) + fc(:,ic)) ) then
            !
            !  Remember time of the return, and count the number of times a particular
            !  grid gets returned to.
            !
            tr%s_pnReturnToGrid(nt ) = nt2 
            tr%s_pbSaveThisStep(nt2) = tr%s_pbSaveThisStep(nt2) + 1
            exit departure_times
          end if
        end do
      end do departure_times
      !
      if (tr%dump_returns) then
        nt2  = tr%s_pnReturnToGrid(nt)
        tret = -1e2_rk
        if (nt2>0) tret = tr%s_pft(nt2)
        write(io_unit,"(1x,i12,1x,f12.4,1x,i12,f12.4)") nt, tr%s_pft(nt), nt2, tret
      end if
      !
    end do reentry_times
    !
    if (tr%dump_returns) then
      close(io_unit)
    end if
    !
    if (tr%verbose>=0) then
      write (out,"(/' Trajectories re-enter the grid at ',i8,' distinct time steps (out of ',i8,')')") &
             count(tr%s_pnReturnToGrid/=0), tr%s_nNumT
      write (out,"( ' Re-entrant trajectories leave the grid at ',i8,' distinct time steps.')") &
             count(tr%s_pbSaveThisStep>0)
      write (out,"( ' ',i8,' time steps have more than one trajectory terminating there.')") &
             count(tr%s_pbSaveThisStep>1)
      write (out,"( ' Maximum number of distinct trajectories terminating at the same time step is ',i8/)") &
             maxval(tr%s_pbSaveThisStep)
    end if
    !
    call TimerStop('EVA SetReturnToGrid')
  end subroutine SetReturnToGrid
  !
  !  Sets to expansion coefficients for the trajectories in the form:
  !
  !     x(t) = Cx1*t + (Cx2/2)*t^2 + (Cx3/6)*t^3
  !     y(t) = Cy1*t + (Cy2/2)*t^2 + (Cy3/6)*t^3
  !     z(t) = Cz1*t + (Cz2/2)*t^2 + (Cz3/6)*t^3
  !
  !  with the convention that x(t=0) = y(t=0) = z(t=0) = 0.
  !
  !  The derivatives for the expansion are calculated using finite difference formulae.
  !
  subroutine SetTrajCoeffs(ftup,s_fC1,s_fC2,s_fC3)
    real(rk), intent(in)  :: ftup     ! Time point at which expansion is needed
    real(rk), intent(out) :: s_fC1(3)
    real(rk), intent(out) :: s_fC2(3)
    real(rk), intent(out) :: s_fC3(3)
    !
    real(rk) :: feps
    real(rk) :: fTmm(3), fTm(3), fT0(3), fTp(3), fTpp(3)
    !
    feps = tr%s_fdt*tr%derivative_step 
    !
    fTmm = trajectory_segment( ftup, ftup-feps*2.0_rk )
    fTm  = trajectory_segment( ftup, ftup-feps        )
    fT0  = trajectory_segment( ftup, ftup             )
    fTp  = trajectory_segment( ftup, ftup+feps        )
    fTpp = trajectory_segment( ftup, ftup+feps*2.0_rk )
    !
    s_fC1 = (fTp - fTm) / (2.0_rk*feps)
    s_fC2 = (fTp - 2.0_rk*fT0 + fTm) / (feps**2)
    s_fC3 = (fTpp - 2.0_rk*fTp + 2.0_rk*fTm - fTmm) / (2.0_rk*feps**3)
    !
  end subroutine SetTrajCoeffs

  subroutine SetAllGWeights(f0, pfWG)
    real(rk), intent(in)  :: f0(3)       ! Displacement of the initial point of the trajectory
                                         ! (the final point is at (0,0,0)
    real(rk), intent(out) :: pfWG(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)
    integer(ik)           :: ni, nj, nk
    !
    do nk=-1,1
    do nj=-1,1
    do ni=-1,1
      !
      call EVA_SetIntGWeights( tr%s_fd(1), tr%s_fd(2), tr%s_fd(3), &
                               f0(1)+ni*tr%s_fd(1), f0(2)+nj*tr%s_fd(2), f0(3)+nk*tr%s_fd(3), &
                               pfWG(:,:,:,ni,nj,nk) )
      !
    end do
    end do
    end do
    !
  end subroutine SetAllGWeights
  !
  !  Interpolates the starting phase for the 'fill edge' functions
  !
  function InterpStartingPhase(fxyz) result (v)
    real(rk), intent(in)    :: fxyz(3)               ! Point at which to interpolate the phase
                                                     ! The point MUST be at the grid boundary
    real(rk)                :: v                     ! Interpolated phase
    !
    integer(ik)             :: nx, ny, nz            ! Nearest grid point
    real(rk)                :: fD(3)                 ! Distance of interpolation point from grid point
    real(rk)                :: pfWG(-1:1,-1:1,-1:1)  ! Interpolation stencil
    !
    !  Find the nearest grid point
    !
    nx = nint( (fxyz(1)-tr%s_pfGrid(1,-1)) / tr%s_fd(1) ) + 1
    ny = nint( (fxyz(2)-tr%s_pfGrid(2,-1)) / tr%s_fd(2) ) + 1
    nz = nint( (fxyz(3)-tr%s_pfGrid(3,-1)) / tr%s_fd(3) ) + 1
    !
    !  Check to be away from last layer
    !
    nx = min(max(nx,3),tr%s_nNum(1)-2)
    ny = min(max(ny,3),tr%s_nNum(2)-2)
    nz = min(max(nz,3),tr%s_nNum(3)-2)
    !
    !  Get distance from nearest grid point
    !
    fD = fxyz - (/ tr%s_pfX(nx), tr%s_pfY(ny), tr%s_pfZ(nz) /)
    !
    !  Sanity check
    !
    if (any(abs(fD)>2._rk*tr%s_fd)) then
      write (out,"('Phase interpolation requested at xyz = ',3f14.5)") fxyz
      write (out,"('Displacement from the nearest grid point (',3i6,') is ',3f14.5)") nx, ny, nz, fD
      stop 'eikonal_tools_eva%InterpStartingPhase - unsafe extrapolation'
    end if
    !
    !  Set interpolation stencil
    !
    call EVA_SetIntGWeights( tr%s_fd(1), tr%s_fd(2), tr%s_fd(3), fD(1), fD(2), fD(3), pfWG )
    !
    !  Interpolate, choosing the appropriate edge
    !
         if (nx==3) then ;              v = sum(pfWG*real(tr%edge_yz(   1:3,   ny-1:ny+1,nz-1:nz+1),kind=rk))
    else if (ny==3) then ;              v = sum(pfWG*real(tr%edge_xz(nx-1:nx+1,   1:3,   nz-1:nz+1),kind=rk))
    else if (nz==3) then ;              v = sum(pfWG*real(tr%edge_xy(nx-1:nx+1,ny-1:ny+1,   1:3   ),kind=rk))
    else if (nx==tr%s_nNum(1)-2) then ; v = sum(pfWG*real(tr%edge_yz(   4:6,   ny-1:ny+1,nz-1:nz+1),kind=rk))
    else if (ny==tr%s_nNum(2)-2) then ; v = sum(pfWG*real(tr%edge_xz(nx-1:nx+1,   4:6,   nz-1:nz+1),kind=rk))
    else if (nz==tr%s_nNum(3)-2) then ; v = sum(pfWG*real(tr%edge_xy(nx-1:nx+1,ny-1:ny+1,   4:6   ),kind=rk))
    else
      write (out,"('Trajectory ends around point ',3i6)") nx, ny, nz
      write (out,"('                Grid ends at ',3i6)") tr%s_nNum
      stop 'eikonal_tools_eva%InterpStartingPhase - trajectory ends at non-edge point'
    end if
  end function InterpStartingPhase
  !
  !
  function GetFullEdgeInt(f0) result(v)
    real(rk), intent(in)  :: f0(3)  ! Desired final point at step s_nt+1
    real(rk)              :: v
    !
    real(rk)              :: f1(3)  ! Starting point of the trajectory at step 1
    integer(ik)           :: nt
    integer(ik)           :: ntstop
    !
    !  Shifts to move global trajectory to current time and position
    !
    f1 = f0 - tr%s_pfRt(:,tr%s_nt+1)   
    !
    !  Check for trajectory that returns to grid.  
    !  The point at which the grid is re-entered is 'ntstop'.
    !  If trajectory never re-enters grid, then ntstop = 1
    !  Note that we do not need to differentiate between trajectories
    !  re-entering at time step 1 and trajectories which do not re-enter
    !  at all: we can use the asymptotic phase at the boundary in either
    !  case.
    !
    find_return_point: do ntstop=tr%s_nt,2,-1
      if ( IsInGrid(tr%s_pfRt(:,ntstop)+f1) ) exit find_return_point
    end do find_return_point
    !
    !  Add up integral outside the main grid
    !
    v = 0.0_rk
    integrate_outside: do nt=ntstop,tr%s_nt
      v = v - real(FLlrePotential(tr%s_pfRt(:,nt)+f1),kind=rk) * ( tr%s_pft(nt+1) - tr%s_pft(nt) )
    end do integrate_outside
    !
    !  Add starting phase when trajectory hits grid (if applicable)
    !
    if (ntstop>1) then
      if (ntstop/=tr%s_pnReturnToGrid(tr%s_nt+1)) then
        write (out,"('Trajectory ending at ',3f12.5,' at time step ',i6,' returns at step ',i6,' but should at ',i6,'. Bummer.')") &
               f0, tr%s_nt+1, ntstop, tr%s_pnReturnToGrid(tr%s_nt+1)
        stop 'eikonal_tools_eva%GetFullEdgeInt - assumption on the return time failed!'
      end if
      v = v + InterpStartingPhase(tr%s_pfRt(:,ntstop)+f1)
    else
      !
      !  Note the change in sign for the value returned by FLexactphase. The sign change is
      !  necessary to reflect different boundary conditions: the adiabatic solution returned
      !  by FLexactphase() is for the flat _outgoing_ front with kvec=-s_fK. So here.
      !
      v = v - real(FLexactphase(tr%s_pfRt(:,1)+f1),kind=rk)
    end if
    if (isnan(v)) then
      write (out,"('GetFullEdgeInt: got Nan at timestep ',i8,' coord ',3f15.7)") tr%s_nt+1, f0
      stop 'eikonal_tools_eva%GetFullEdgeInt - Got a Nan'
    end if
    !
  end function GetFullEdgeInt
  !
  !  Fill tables of coarse-grained edge integrals
  !
  subroutine FillCoarseEdges( nxe, nye, nze )
    integer(ik), intent(in) :: nxe, nye, nze  ! location of the edges to be filled
    !
    integer(ik) :: nx, ny, nz
    !
    call TimerStart('EVA edge trajectories')
    !
    !  This routine implicitly relies on the edges of the fine and coarse grids
    !  coinciding. Let's make sure this is the case.
    !
    if ( .not. ( zero(tr%s_pfXe(1)-tr%s_pfX(nxe)).or.zero(tr%s_pfXe(tr%s_nNumE(1))-tr%s_pfX(nxe)) ) .or. &
         .not. ( zero(tr%s_pfYe(1)-tr%s_pfY(nye)).or.zero(tr%s_pfYe(tr%s_nNumE(2))-tr%s_pfY(nye)) ) .or. &
         .not. ( zero(tr%s_pfZe(1)-tr%s_pfZ(nze)).or.zero(tr%s_pfZe(tr%s_nNumE(3))-tr%s_pfZ(nze)) ) ) then
      write (out,"('    Coarse low: ',3f20.10)") tr%s_pfXe(1), tr%s_pfYe(1), tr%s_pfZe(1)
      write (out,"('Edge positions: ',3f20.10)") tr%s_pfX(nxe), tr%s_pfY(nye), tr%s_pfZ(nze)
      write (out,"('   Coarse high: ',3f20.10)") &
             tr%s_pfXe(tr%s_nNumE(1)), tr%s_pfYe(tr%s_nNumE(2)), tr%s_pfZe(tr%s_nNumE(3))
      stop 'eikonal_tools_eva%FillCoarseEdges - bad grid alignement'
    end if
    !
    !  Calculate coarse-grained edge integrals
    !
    !$omp parallel private(nx,ny,nz)
    !$omp do 
    do ny=1,tr%s_nNumE(2)
      do nx=1,tr%s_nNumE(1)
        tr%s_pfGxy(nx,ny) = GetFullEdgeInt( (/tr%s_pfXe(nx), tr%s_pfYe(ny), tr%s_pfZ(nze)/) )
      end do
    end do
    !$omp end do nowait
    !
    !$omp do
    do nz=1,tr%s_nNumE(3)
      do ny=1,tr%s_nNumE(2)
        tr%s_pfGyz(ny,nz) = GetFullEdgeInt( (/tr%s_pfX(nxe), tr%s_pfYe(ny), tr%s_pfZe(nz)/) )
      end do
    end do
    !$omp end do nowait
    !
    !$omp do
    do nz=1,tr%s_nNumE(3)
      do nx=1,tr%s_nNumE(1)
        tr%s_pfGxz(nx,nz) = GetFullEdgeInt( (/tr%s_pfXe(nx), tr%s_pfY(nye), tr%s_pfZe(nz)/) )
      end do
    end do
    !$omp end do
    !$omp end parallel
    !
    call TimerStop('EVA edge trajectories')
    !
    contains
    logical function zero(diff)
      real(rk), intent(in) :: diff
      zero = abs(diff)<=spacing(10._rk)
    end function zero
  end subroutine FillCoarseEdges
  !
  !  Interpolate coarse-grained exact integrals at the boundary
  !
  !  If necessary, this routine can be accelerated by caching the 
  !  interpolation matrices - they do not change between time steps.
  !
  subroutine InterpolateOutsideEdges
    integer(ik) :: nx, ny, nz, nIDx, nIDy, nIDz
    real(rk)    :: fXp, fYp, fZp
    real(rk)    :: pfW(-1:1,-1:1)
    !
    call TimerStart('EVA edge interpolation')
    !
    !  Fill extreme boundaries with zero, to avoid junk there
    !$omp parallel private(nx,ny,nz,nIDx,nIDy,nIDz,fXp,fYp,fZp,pfW)
    !$omp sections
    !$omp section
      tr%interp_xy(:,:,:) = 0
    !$omp section
      tr%interp_xz(:,:,:) = 0
    !$omp section
      tr%interp_yz(:,:,:) = 0
    !$omp end sections
    !
    !  Fill xy edge
    !
    !$omp do
    fill_xy: do ny=2,tr%s_nNum(2)-1
      nIDy = nint( (tr%s_pfY(ny) - tr%s_pfGrid(2,-1)) / tr%s_fdE(2) ) + 1
      nIDy = min(max(nIDy,2),tr%s_nNumE(2)-1)
      fYp  = tr%s_pfY(ny) - tr%s_pfYe(nIDy)
      do nx=2,tr%s_nNum(1)-1
        nIDx = nint( (tr%s_pfX(nx) - tr%s_pfGrid(1,-1)) / tr%s_fdE(1) ) + 1
        nIDx = min(max(nIDx,2),tr%s_nNumE(1)-1)
        fXp  = tr%s_pfX(nx) - tr%s_pfXe(nIDx)
        !
        call EVA_SetInterpWeights2D( tr%s_fdE(1), tr%s_fdE(2), fXp, fYp, pfW )
        tr%interp_xy(nx,ny,1) = sum(pfW*tr%s_pfGxy(nIDx-1:nIDx+1,nIDy-1:nIDy+1))
      end do
    end do fill_xy
    !$omp end do nowait
    !
    ! Fill yz edge
    !
    !$omp do
    fill_yz: do nz=2,tr%s_nNum(3)-1
      nIDz = nint( (tr%s_pfZ(nz) - tr%s_pfGrid(3,-1)) / tr%s_fdE(3) ) + 1
      nIDz = min(max(nIDz,2),tr%s_nNumE(3)-1)
      fZp  = tr%s_pfZ(nz) - tr%s_pfZe(nIDz)
      do ny=2,tr%s_nNum(2)-1
        nIDy = nint( (tr%s_pfY(ny) - tr%s_pfGrid(2,-1)) / tr%s_fdE(2) ) + 1
        nIDy = min(max(nIDy,2),tr%s_nNumE(2)-1)
        fYp  = tr%s_pfY(ny) - tr%s_pfYe(nIDy)
        !
        call EVA_SetInterpWeights2D( tr%s_fdE(2), tr%s_fdE(3), fYp, fZp, pfW )
        tr%interp_yz(1,ny,nz) = sum(pfW*tr%s_pfGyz(nIDy-1:nIDy+1,nIDz-1:nIDz+1))
      end do
    end do fill_yz
    !$omp end do nowait
    !
    ! Fill xz edge
    !
    !$omp do
    fill_xz: do nz=2,tr%s_nNum(3)-1
      nIDz = nint( (tr%s_pfZ(nz) - tr%s_pfGrid(3,-1)) / tr%s_fdE(3) ) + 1
      nIDz = min(max(nIDz,2),tr%s_nNumE(3)-1)
      fZp  = tr%s_pfZ(nz) - tr%s_pfZe(nIDz)
      do nx=2,tr%s_nNum(1)-1
        nIDx = nint( (tr%s_pfX(nx) - tr%s_pfGrid(1,-1)) / tr%s_fdE(1) ) + 1
        nIDx = min(max(nIDx,2),tr%s_nNumE(1)-1)
        fXp  = tr%s_pfX(nx) - tr%s_pfXe(nIDx)
        !
        call EVA_SetInterpWeights2D( tr%s_fdE(1), tr%s_fdE(3), fXp, fZp, pfW )
        tr%interp_xz(nx,1,nz) = sum(pfW*tr%s_pfGxz(nIDx-1:nIDx+1,nIDz-1:nIDz+1))
      end do
    end do fill_xz
    !$omp end do
    !$omp end parallel
    call TimerStop('EVA edge interpolation')
    !
  end subroutine InterpolateOutsideEdges
  !
  !
  !
  subroutine StepPhase(potential,src,dst)
    integer(ik), intent(in)  :: potential  ! Input: Data field containing the potential
    integer(ik), intent(in)  :: src        ! Input: Phase at the previous time step
    integer(ik), intent(in)  :: dst        ! Output: Updated phase
    !
    real(rk)    :: ft                               ! Current time
    real(rk)    :: fdt                              ! Time step
    integer(ik) :: nXedgeOut, nYedgeOut, nZedgeOut  ! edges where the trajectory falls outside the grid
    real(rk)    :: f0(3)
    real(rk)    :: s_fC1(3)                         ! Trajectory expansion coefficients
    real(rk)    :: s_fC2(3)
    real(rk)    :: s_fC3(3)
    integer(ik) :: top      ! Trajectory topology code - see FieldEVAExteriorStep in multigrid.f90
    !
    ft  = tr%s_pft(tr%s_nt)
    fdt = tr%s_pft(tr%s_nt+1)-tr%s_pft(tr%s_nt)
    !
    if (tr%verbose>=2) then
      write (out,"('Advancing to time step ',i6,'. Previous time ',f12.4,', time change ',f12.4)") &
             tr%s_nt+1, ft, fdt
    end if
    !
    !  Set trajectory info; f0(:) is the start of the trajectory segment
    !  leading to (0,0,0) position at the end of the time step.
    !
    f0 = tr%s_pfRt(:,tr%s_nt) - tr%s_pfRt(:,tr%s_nt+1)
    !
    call SetTrajCoeffs(ft+fdt, s_fC1, s_fC2, s_fC3 )
    !
    !  Set weights
    !
    call SetAllGWeights( f0, tr%s_pfWG )
    call EVA_SetIntUWeights( -fdt, tr%s_fd(1), tr%s_fd(2), tr%s_fd(3), &
                                  s_fC1(1), s_fC2(1), s_fC3(1), &
                                  s_fC1(2), s_fC2(2), s_fC3(2), &
                                  s_fC1(3), s_fC2(3), s_fC3(3), tr%s_pfWU )
    !
    !  Handle the interior points
    !
    if (tr%verbose>=2) then
      write (out,"('Advancing the interior points')") 
    end if
    !
    call FieldEVAInteriorStep(src_u=potential, src_g=src, dst_g=dst, &
                              wgt_u=tr%s_pfWU, wgt_g=tr%s_pfWG(-1:1,-1:1,-1:1,0,0,0))
    !
    !  Fill "inside" edges using a variant of the general interpolation scheme
    !
    top = 0
    if (f0(1) < 0._rk) top = top + 1
    if (f0(2) < 0._rk) top = top + 2
    if (f0(3) < 0._rk) top = top + 4
    !
    if (tr%verbose>=2) then
      write (out,"('Advancing the ""inside"" boundary, topology code ',i3)") top
    end if
    !
    call FieldEVAExteriorStep(src_u=potential, src_g=src, dst_g=dst, topology=top, &
                              wgt_u=tr%s_pfWU, wgt_g=tr%s_pfWG)
    !
    !  Fill "outside" edges, by interpolating coarse-grained exact integrals
    !
    nXedgeOut = 2 ; if (f0(1) >= 0.0) nXedgeOut = tr%s_nNum(1)-1
    nYedgeOut = 2 ; if (f0(2) >= 0.0) nYedgeOut = tr%s_nNum(2)-1
    nZedgeOut = 2 ; if (f0(3) >= 0.0) nZedgeOut = tr%s_nNum(3)-1
    !
    if (tr%verbose>=2) then
      write (out,"('Advancing ""outside"" points. Edge positions ',3i5)") &
             nXedgeOut, nYedgeOut, nZedgeOut
    end if
    !
    call FillCoarseEdges( nXedgeOut, nYedgeOut, nZedgeOut )
    call InterpolateOutsideEdges
    call FieldSetSlab(dst=dst,dir=1,n1=nXedgeOut,n2=nXedgeOut,data=tr%interp_yz)
    call FieldSetSlab(dst=dst,dir=2,n1=nYedgeOut,n2=nYedgeOut,data=tr%interp_xz)
    call FieldSetSlab(dst=dst,dir=3,n1=nZedgeOut,n2=nZedgeOut,data=tr%interp_xy)
    !
  end subroutine StepPhase
  !
  !  Save edge points of the current field
  !
  subroutine SaveCurrentEdges(src)
    integer(ik), intent(in) :: src ! Data field containing the phase
    !
    call TimerStart('EVA SaveCurrentEdges')
    !
    !  Extract edges of the current phase grid
    !
    call FieldFetchSlab(src,dir=3_ik,n1= 2_ik,n2= 4_ik,data=tr%edge_xy(:,:,1:3))
    call FieldFetchSlab(src,dir=3_ik,n1=-4_ik,n2=-2_ik,data=tr%edge_xy(:,:,4:6))
    call FieldFetchSlab(src,dir=2_ik,n1= 2_ik,n2= 4_ik,data=tr%edge_xz(:,1:3,:))
    call FieldFetchSlab(src,dir=2_ik,n1=-4_ik,n2=-2_ik,data=tr%edge_xz(:,4:6,:))
    call FieldFetchSlab(src,dir=1_ik,n1= 2_ik,n2= 4_ik,data=tr%edge_yz(1:3,:,:))
    call FieldFetchSlab(src,dir=1_ik,n1=-4_ik,n2=-2_ik,data=tr%edge_yz(4:6,:,:))
    !
    write(tr%io_edge,rec=tr%next_edge_record) tr%edge_xy, tr%edge_xz, tr%edge_yz
    !
    if (tr%verbose>=2) then
      write (out,"('Edges at time step ',i6,' saved to record ',i6)") &
             tr%s_nt, tr%next_edge_record
    end if
    !
    tr%edge_record(tr%s_nt) = tr%next_edge_record
    tr%next_edge_record = tr%next_edge_record + 1
    !
    call TimerStop('EVA SaveCurrentEdges')
  end subroutine SaveCurrentEdges

  subroutine LoadOldEdges( nt )
    integer(ik), intent(in) :: nt
    !
    integer(ik) :: rec
    !
    call TimerStart('EVA LoadOldEdges')
    rec = tr%edge_record(nt)
    if (rec==0) then
      write (out,"('eikonal_tools_eva%LoadOldEdges: asked for time step ',i6,' - but it aint here')") nt
      stop 'eikonal_tools_eva%LoadOldEdges - missing record'
    end if
    !
    if (tr%verbose>=2) then
      write (out,"('Loading edge data for time step ',i6,' from record ',i6)") nt, rec
    end if
    !
    read (tr%io_edge,rec=rec) tr%edge_xy, tr%edge_xz, tr%edge_yz
    !
    call TimerStop('EVA LoadOldEdges')
  end subroutine LoadOldEdges

end module eikonal_tools_eva
