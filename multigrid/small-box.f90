!
!  Stationary solutions in a relatively small box - direct
!  diagonalization. This example solves for the stationary 
!  states of hydrogen molecule adsorbed on a 2D periodic
!  graphene sheet.
!
  module small_box
    use accuracy
    use multigrid
    use lapack
    use fftw
    use timer
    implicit none
!
    integer(ik), parameter   :: verbose   = 3
    integer(ik), parameter   :: unit_save = 34
    !
    !  Particle mass, in electron mass units
    !
    real(rk), parameter      :: mass     = (1.0_rk+1838.7_rk) * 2.0_rk 
    !
    !  Bohr length
    !
    real(rk), parameter      :: a0       = 0.52918_rk
    !
    !  Carbon-carbon bond length, in Bohr
    !
    real(rk), parameter      :: bond_cc  = 1.421_rk / a0 ! Experimental value for graphite
    !
    !  Cell and grid parameters
    !
    ! Box 14, 18x21x20 - double layer, at 7.00A, cell starts at 1.5A and ends at 4.5A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  7.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.5 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  5.5 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 13, 18x21x20 - double layer, at 5.75A, cell starts at 1.5A and ends at 4.25A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  5.75/ a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.5 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  4.25/ a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 12, 18x21x20 - double layer, at 5.5A, cell starts at 1.0A and ends at 4.5A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  5.5 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  4.5 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 11, 18x21x20 - double layer, at 4A, cell starts at 1.0A and ends at 3.0A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  4.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  3.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 10, 18x21x20 - double layer, at 5A, cell starts at 1.0A and ends at 4.0A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  5.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  4.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 9, 18x21x20 - double layer, at 6A, cell starts at 1.5A and ends at 4.5A
    !                   (this is a special case - it's a very small box!)
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  6.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  1.5 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  4.5 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 8, 18x21x20 - double layer, at 8A, cell ends at 6A
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  8.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  6.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 7, 18x21x30 - double layer, at 14A, cell ends at 12A
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 14.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 12.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 30        ! 
    !
    ! Box 6, 18x21x20 - double layer, at 10A, cell ends at 8A
    !                   This spacing is slightly better than Box 3
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 10.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  8.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 6+, 18x21x20 - double layer, at 9A, cell ends at 7A
    !
      logical, parameter       :: layer_2    = .true.    ! Include second layer?
      real(rk), parameter      :: layer_2_z  =  9.0 / a0 ! Z-coordinate of the 2nd layer
      integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
      integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
      real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
      real(rk), parameter      :: max_z      =  7.0 / a0 ! Max. distance for the potential
      integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
      integer(ik)              :: n_points_b = 21        ! direction; 
      integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 5, 18x21x20 - double layer, at 12A, cell ends at 10A
    !                   This spacing is slightly better than Box 3
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 12.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 10.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 18        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 21        ! direction; 
    ! integer(ik)              :: n_points_c = 20        ! 
    !
    ! Box 4, 11x18x82 - extra convergence test
    !
    ! logical, parameter       :: layer_2    = .false.   ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 0.0       ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 1         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 1         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 30.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 11        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 18        ! direction; 14x16x72 seems to work;
    ! integer(ik)              :: n_points_c = 82        ! 15x17x68 seems to work as well;
    !
    !
    ! Box 3, 15x18x48 - pictures production run
    !
    ! logical, parameter       :: layer_2    = .false.   ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 0.0       ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 24.5 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a = 17        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 20        ! direction; 14x16x72 seems to work;
    ! integer(ik)              :: n_points_c = 48        ! 15x17x68 seems to work as well;
    !
    ! Box 2, 8x13x100 - thermodynamics production run
    !
    ! logical, parameter       :: layer_2    = .false.   ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 0.0       ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 1         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 1         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 45.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a =  8        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 13        ! direction; 14x16x72 seems to work;
    ! integer(ik)              :: n_points_c =100        ! 15x17x68 seems to work as well;
    !
    ! Box 1, 8x14x40 - used for debugging and testing
    !
    ! logical, parameter       :: layer_2    = .false.   ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  = 0.0       ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 1         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 1         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      = 15.0 / a0 ! Max. distance for the potential
    ! integer(ik)              :: n_points_a =  8        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 14        ! direction; 14x16x72 seems to work;
    ! integer(ik)              :: n_points_c = 40        ! 15x17x68 seems to work as well;
    !
    integer(ik)                :: sum_cells  = 10        ! Number of cells to include in potential
                                                         ! summation. 30 is a good number to reach
                                                         ! convergence, but we'll use 10 for testing
    integer(ik)                :: plot_cells = 5         ! Number of cells to include in background
                                                         ! image plot
    !
    !  Definitions of the interacting centres in the unit cell
    !
    real(rk), parameter      :: sqrt3        = 1.73205080756887729352_rk
    real(rk), parameter      :: centres(3,4) = reshape( &
      (/                                                &
         0.0_rk,            bond_cc/2.0,  0.0_rk,  &
         0.0_rk,           -bond_cc/2.0,  0.0_rk,  &
         bond_cc*sqrt3/2,   bond_cc,      0.0_rk,  &
         bond_cc*sqrt3/2,  -bond_cc,      0.0_rk   &
      /), (/ 3, 4 /) )
    !
    !  The rest of the stuff will be done automatically
    !
    integer(ik)              :: n_states        ! Total number of states
    complex(rk), allocatable :: h(:,:)          ! Hamiltonian matrix / eigenstates
    real(rk), allocatable    :: e(:)            ! Eigenvalues
    real(rk), allocatable    :: e_empty(:)      ! Eigenvalues - empty box
    complex(rk), allocatable :: heev_work(:)    ! Scratch space for LAPACK
    real(rk), allocatable    :: heev_rwork(:)
!
    integer(ik)              :: pot, fft
!
    contains
!
    function vdw_potential(r) result(v)
      real(rk), intent(in) :: r(3)
      real(rk)             :: v
      !
      real(rk)             :: rx
      !
      rx = sqrt(sum(r**2)) * a0  ! Bohr -> Angstroms
      !
      !  Potential A - fit to MP2/cc-pvtz in benzene
      !
      ! v = -617.58_rk/rx**6 + 14040.0_rk * exp(-3.2696_rk * rx)
      !
      !  Potential B - fit to CP-CCSD(T)/cc-pvtz in benzene
      !
      ! v = -481.94_rk/rx**6 + 12577.0_rk * exp(-3.2508_rk * rx)
      !
      !  Potential C - fit to MP2/cc-pvtz in coronene
      ! 
        v = -400.37_rk/rx**6 + 25352.0_rk * exp(-3.5763_rk * rx)
      !
      v = v / ( 627.51_rk ) ! kcal/mol -> Hartree
    end function vdw_potential
!
    function interactionLayer(x,z) result(v)
      real(rk), intent(in) :: x(3)  ! Coordinates of the point
      real(rk), intent(in) :: z     ! Z-coordinate of the layer
      complex(rk)          :: v
      !
      integer(ik) :: cell_a, cell_b, atom
      real(rk)    :: cell_d(3)
      !
      v = 0
      sum_a: do cell_a=-sum_cells,sum_cells
        sum_b: do cell_b=-sum_cells,sum_cells
          cell_d(1) = bond_cc * sqrt(3.0_rk) * cell_a
          cell_d(2) = bond_cc * 3.0_rk       * cell_b
          cell_d(3) = z
          sum_at: do atom=1,size(centres,dim=2)
            v = v + vdw_potential(x-centres(:,atom)-cell_d)
          end do sum_at
        end do sum_b
      end do sum_a
      return
    end function interactionLayer
    !
    function interactionPotential(x) result(v)
      real(rk), intent(in) :: x(3)
      complex(rk)          :: v
      !
      ! Debug - quadratic well
      !
      ! v = 0.20_rk * sum((x-(/0._rk,0._rk,3._rk/))**2)
      ! return
      v = interactionLayer(x,0.0_rk)
      if (layer_2) then
        v = v + interactionLayer(x,layer_2_z)
      end if
    end function interactionPotential
!
    subroutine nuc
      integer(ik)        :: root, info
      character(len=200) :: buf
!
      call TimerStart('nuc')
      call accuracyInitialize
!
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
!
      call plotGraphene
      call buildGrid
      call addFields
!
!    For a very small field we have here, we'll simply diagonalize
!    the Hamiltonian matrix. FFT transform of the potential will 
!    give us the potential part. The kinetic part is diagonal.
!
      call TimerStart('stationary problem')
      call FieldInit(pot,interactionPotential)
      call Visualize(pot,'Interaction potential')
!
      call buildHamiltonian
      !
      write (out,"(/t5'Empty box eigenvalues are:'/)")
      write (out,"((t3,8(f12.6,1x)))") e_empty
      !
      call TimerStart('Diagonalization')
      !
      !  Our interface to LAPACK is extremely inefficient (it's designed
      !  for small matrices). Therefore, we'll have to call LAPACK directly.
      !
      !  call lapack_heev(h,e)
      !
      allocate (heev_work(64*n_states),heev_rwork(3*n_states),stat=info)
      if (info/=0) then
        write (out,"('Memory allocation failed, code = ',i7)") info
        write (out,"('Array sizes were: ',5i10)") 64*n_states, 3*n_states
        pause 'smallbox - diagonalization scratch'
      end if
      call cheev('V','U',n_states,h,n_states,e,heev_work,size(heev_work),heev_rwork,info)
      deallocate (heev_work,heev_rwork)
      if (info/=0) then
        write (out,"('cheev failed: ',i6)") info
        stop 'cheev failed'
      end if
      !
      call TimerStop('Diagonalization')
      call TimerStop('stationary problem')
      call TimerReport
!
      write (out,"(/t5'Eigenvalues are:'/)")
      write (out,"((t3,8(f12.6,1x)))") e
!
      call save_all_roots
!
      display_roots: do root=1,n_states
        write(out,"('Choose state to visualize; 0 to exit',i5)")
!       read(5,*) root
!       if (root<=0) exit
        write (buf,"(' Root = ',i5,' energy = ',f12.6)") root, e(root)
        write (out,"(a)") trim(buf)
        call buildRoot(root)
        call Visualize(pot,trim(buf))
        call sleep(20)
      end do display_roots
!
      call TimerStop('nuc')
      call TimerReport
    end subroutine nuc

    subroutine save_all_roots
      integer(ik) :: root
      !
      open (unit=unit_save,file='small-box-eigenvalues.dat', &
            form='formatted',status='replace')
      write (unit_save,"('# ',a5,2x,a25,2x,a25)") &
             'Root', 'Empty box', 'Box with potential'
      do root=1,n_states
        write (unit_save,"(2x,i5,2x,f25.12,2x,f25.12)") &
               root, e_empty(root), e(root)
      end do 
      close (unit=unit_save)
    end subroutine save_all_roots

    subroutine buildHamiltonian
      integer(ik)           :: v0x, v0y, v0z  ! Zero-point of the potential
      integer(ik)           :: kRx, kRy, kRz  ! Right-hand side
      integer(ik)           :: kLx, kLy, kLz  ! Left-hand side
      integer(ik)           :: kDx, kDy, kDz  ! Difference
      integer(ik)           :: hL, hR         ! Left and right compount indices
      real(rk)              :: ekin
      complex(rk)           :: epot
      real(rk)              :: ekin_max
      integer(ik)           :: alloc
      integer(ik)           :: npts(3)        ! Number of points at this grid level
      real(rk)              :: step(3)        ! Grid spacing
      real(rk), pointer     :: coord(:,:,:,:) ! Coordinates and integration weights
      complex(rk), pointer  :: field(:,:,:)   ! Field values
      !
      complex(rk), pointer  :: oversample_fft(:,:,:)
      !
      call TimerStart('buildHamiltonian')
      call FieldCheckOut(typ='MOMENTUM',src=fft,grd=1,npts=npts, &
                         step=step,coord=coord,field=field,expand=.true.)
      !
      !  Build FFT of the potential on double-frequency grid
      !
      call build_potential_FFT(oversample_fft)
      !
      n_states = product(npts  ) ! Exclude states with zero momentum
      write (out,"('Total number of states = ',i8)") n_states
      allocate(e(n_states),e_empty(n_states),h(n_states,n_states),stat=alloc)
      if (alloc/=0) then
        write (out,"('Allocation failed, code = ',i5)") alloc
        write (out,"('Array sizes were: ',5i10)") n_states, n_states, n_states**2
        pause 'buildHamiltonian - no memory'
      end if
      !
      !  Find the zero-point of the transform
      !
      v0x = 1+(size(oversample_fft,dim=1)-1)/2
      v0y = 1+(size(oversample_fft,dim=2)-1)/2
      v0z = 1+(size(oversample_fft,dim=3)-1)/2
      !
      ekin_max = 0
      h  = 0 
      hR = 0
      right_z: do kRz=1,npts(3)
        right_y: do kRy=1,npts(2)
          right_x: do kRx=1,npts(1)
            hR = hR + 1
            !
            !  Diagonal kinetic part. The product(step) part is coming from
            !  the normalization of the wavefunction
            !
            ekin        = sum(coord(1:3,kRx,kRy,kRz)**2)/(2.0_rk*mass)
            ekin_max    = max(ekin,ekin_max)
            h(hR,hR)    = h(hR,hR) + ekin
            e_empty(hR) = ekin
            !
            hL = 0
            left_z: do kLz=1,npts(3)
              left_y: do kLy=1,npts(2)
                left_x: do kLx=1,npts(1)
                  hL = hL + 1
                  !
                  kDx = v0x + kRx - kLx
                  kDy = v0y + kRy - kLy
                  kDz = v0z + kRz - kLz
                  !
                  epot = oversample_fft(kDx,kDy,kDz)
                  h(hL,hR) = h(hL,hR) + epot
                  !
                end do left_x
              end do left_y
            end do left_z
            if (hL/=n_states) stop 'buildHamiltonian - count error 2'
            !
          end do right_x
        end do right_y
      end do right_z
      !
      if (hR/=n_states) stop 'buildHamiltonian - count error 1'
      !
      write (out,"(/'Max. kinetic energy in this box = ',f12.6,' Hartree')") ekin_max
      !
      ! if (verbose>=3) then
      !   do hR=1,n_states
      !     do hL=1,hR
      !       write (out,"(2x,2i4,2f14.6,2x,2f14.6)") hL, hR, h(hL,hR), h(hL,hR) - conjg(h(hR,hL))
      !     end do
      !   end do
      ! end if
      !
      deallocate (oversample_fft)
      !
      call TimerStop('buildHamiltonian')
      !
      call TimerStart('stupidSort')
      call stupid_sort(e_empty)
      call TimerStop('stupidSort')
      !
    end subroutine buildHamiltonian

    subroutine build_potential_FFT(over_fft)
      complex(rk), pointer     :: over_fft(:,:,:)
      complex(rk), allocatable :: over_pot(:,:,:)
      !
      integer(ik)           :: npts(3)        ! Number of points at this grid level
      real(rk)              :: step(3)        ! Grid spacing
      real(rk), pointer     :: coord(:,:,:,:) ! Coordinates and integration weights
      complex(rk), pointer  :: field(:,:,:)   ! Field values
      !
      integer(ik)           :: over_npts(3)
      integer(ik)           :: ox, oy, oz     ! Oversampled grid size
      integer(ik)           :: ix, iy, iz     ! Oversampled index
      integer(ik)           :: sx, sy, sz     ! Subsampled index
      real(rk)              :: x(3)           ! Oversampled position
      integer(ik)           :: alloc
      !
      call TimerStart('build_potential_FFT')
      !
      call FieldCheckOut(typ='COORDINATE',src=pot,grd=1,npts=npts, &
                         step=step,coord=coord,field=field,expand=.true.)
      !
      over_npts = 2*npts - 1
      ox = over_npts(1)
      oy = over_npts(2)
      oz = over_npts(3)
      !
      allocate (over_fft(ox,oy,oz),over_pot(ox,oy,oz),stat=alloc)
      if (alloc/=0) then
        write (out,"('Allocation failed, code = ',i6)") alloc
        write (out,"('Array sizes were: ',5i10)") ox*oy*oz, ox*oy*oz
        pause 'build_potential - no memory'
      end if
      !
      fill_z: do iz=1,oz
        sz   = (iz+1)/2
        fill_y: do iy=1,oy
          sy = (iy+1)/2
          fill_x: do ix=1,ox
            sx = (ix+1)/2
            x  = coord(1:3,sx,sy,sz) + 0.5_rk * (/ ix+1-2*sx, iy+1-2*sy, iz+1-2*sz /) * step
            over_pot(ix,iy,iz) = interactionPotential(x)
          end do fill_x
        end do fill_y
      end do fill_z
      !
      call fftw_3d(ox,oy,oz,over_pot,over_fft,.false.)
      over_fft = over_fft / product(over_npts)
      !
      deallocate (over_pot) 
      !
      call TimerStop('build_potential_FFT')
      !
    end subroutine build_potential_FFT

    subroutine stupid_sort(e)
      real(rk), intent(inout) :: e(:)
      !
      integer(ik) :: i
      logical     :: swap
      real(rk)    :: x
      !
      swap = .true.
      do while(swap)
        swap = .false.
        do i=2,size(e)
          if(e(i-1)>e(i)) then
            swap   = .true.
            x      = e(i-1)
            e(i-1) = e(i)
            e(i)   = x
          end if
        end do
      end do
    end subroutine stupid_sort

    subroutine buildRoot(root)
      integer(ik), intent(in) :: root           ! Root to plot
      integer(ik)             :: k0x, k0y, k0z  ! Zero-point
      integer(ik)             :: kRx, kRy, kRz  ! K-point indices
      integer(ik)             :: hR             ! Compound indices
      integer(ik)             :: npts(3)        ! Number of points at this grid level
      real(rk)                :: step(3)        ! Grid spacing
      real(rk), pointer       :: coord(:,:,:,:) ! Coordinates and integration weights
      complex(rk), pointer    :: field(:,:,:)   ! Field values
      !
      call TimerStart('buildRoot')
      call FieldCheckOut(typ='MOMENTUM',src=fft,grd=1,npts=npts, &
                         step=step,coord=coord,field=field,expand=.true.)
      !
      k0x = 1+(npts(1)-1)/2
      k0y = 1+(npts(2)-1)/2
      k0z = 1+(npts(3)-1)/2
      !
      write (out,"(/t4,'Root = ',i5,' eigenvector norm = ',f14.7)") &
             root, sqrt(sum(abs(h(:,root))**2))
      !
      field = 0
      hR    = 0
      field(:,:,:) = 0
      right_z: do kRz=1,npts(3)
        right_y: do kRy=1,npts(2)
          right_x: do kRx=1,npts(1)
            hR = hR + 1
            !
            field(kRx,kRy,kRz) = conjg(h(hR,root))
            !
          end do right_x
        end do right_y
      end do right_z
      !
      if (hR/=n_states) stop 'buildRoot - count error'
      !
      !  Normalization
      !
      field = field / sqrt(product(step))
      !
      call FieldFFT(fft,pot,inverse=.true.)
      call TimerStop('buildRoot')
      !
    end subroutine buildRoot

    subroutine Visualize(psi,title)
      integer(ik), intent(in)      :: psi       ! Wavefunction; FFT will always use Hpsi
      character(len=*), intent(in) :: title     ! Descriptive label
      real(rk)                     :: efield(3) ! Driving field
      real(rk)                     :: gnorm
      !
      call TimerStart('Visualization')
      !
      efield = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
      !
      !  Visualize density in the real and momentum space
      !
      call FieldVisualize(0,psi,trim(title),1,efield)
      call FieldFFT(psi,fft)
      call FieldVisualize(1,fft,trim(title),1,efield)
      call FieldShow
      !
      call TimerStop('Visualization')
    end subroutine Visualize

    subroutine buildGrid
      real(rk) :: box(2,3)
      !
      call TimerStart('Grid initialization')
      call MultiGridInit(max_grids=1,max_fields=4,nborder=1)
      !
      box(2,1) = bond_cc * sqrt(3.0_rk) * units_a / 2
      box(2,2) = bond_cc * 3.0_rk       * units_b / 2
      box(1,1:2) = -box(2,1:2)
      box(:,3) = (/ min_z, max_z /)
      !
      write (out,"(//t5,'Simulation box size (Bohrs)'/)") 
      write (out,"(t8,'X: ',2f14.6)") box(:,1)
      write (out,"(t8,'Y: ',2f14.6)") box(:,2)
      write (out,"(t8,'Z: ',2f14.6)") box(:,3)
      write (out,"()")
      !
      call SimpleGridNew('Rectangular box', (/  n_points_a, n_points_b, n_points_c /), box)
      !
      call TimerStop('Grid initialization')
    end subroutine buildGrid

    subroutine addFields
      call FieldNew('Potential (Real space)',     pot, scratch=.true. , wavefunction=.false.)
      call FieldNew('Potential (Momentum space)', fft, scratch=.true. , wavefunction=.false.)
    end subroutine addFields

    subroutine plotGraphene
      integer(ik)              :: cell_a, cell_b, cell_at
      integer(ik)              :: iatom, jatom, natoms, ibond, nbonds
      real(rk)                 :: cell_d(3), rij
      real(rk),allocatable     :: atoms (:,:)
      integer(ik), allocatable :: bonds(:,:)
      integer(ik)              :: mol
      !
      ! Count number of atoms in the replicated structure
      !
      natoms = size(centres,dim=2)*(2*plot_cells+1)**2
      if (layer_2) then
        natoms = natoms * 2
      end if
      allocate (atoms(3,natoms),bonds(2,(natoms*(natoms-1))/2))
      !
      ! Replicate atoms
      !
      iatom = 0
      scan_a: do cell_a=-plot_cells,plot_cells
        scan_b: do cell_b=-plot_cells,plot_cells
          cell_d(1) = bond_cc * sqrt(3.0_rk) * cell_a
          cell_d(2) = bond_cc * 3.0_rk       * cell_b
          cell_d(3) = 0
          scan_at: do cell_at=1,size(centres,dim=2)
            iatom = iatom + 1
            atoms(:,iatom) = centres(:,cell_at)+cell_d+(/0._rk,0._rk,2._rk/a0/)
            if (layer_2) then
              iatom = iatom + 1
              atoms(:,iatom) = centres(:,cell_at)+cell_d + (/0._rk,0._rk,layer_2_z-2._rk/a0/)
            end if
          end do scan_at
        end do scan_b
      end do scan_a
      if (iatom/=natoms) then
        write (out,"(' iatom = ',i5,' natoms = ',i5,'. Oops.')") iatom, natoms
        stop 'atom count error'
      endif
      !
      !  Count and record bonds
      !
      ibond = 0
      bond_scan_i: do iatom=1,natoms
        bond_scan_j: do jatom=iatom+1,natoms
          rij = sqrt(sum( (atoms(:,iatom)-atoms(:,jatom))**2 ))
          if (rij>bond_cc*1.10_rk) cycle bond_scan_j
          ibond = ibond + 1
          bonds(1,ibond) = iatom
          bonds(2,ibond) = jatom
        end do bond_scan_j
      end do bond_scan_i
      nbonds = ibond
      !
      !  Write OpenDX output file
      !
      open (mol,form='formatted',status='replace',file='graphene.dx')
      write (mol,"('object 1 class array type float rank 0 "// &
                 "items ',i6,' data follows')") natoms
      write (mol,"(8(1x,i5))") (6,iatom=1,natoms)
      write (mol,"('attribute ""dep"" string ""positions""')")
      !
      !  Our data fields go out to DX in the reverse order - so we'll
      !  need to swap the XYZ coordinates as well.
      !
      write (mol,"('object 2 class array type float rank 1 shape 3"// &
                 " items ',i6,' data follows')") natoms
      write (mol,"(3(1x,f15.8))") atoms(3:1:-1,1:natoms)
      write (mol,"('attribute ""dep"" string ""positions""')")
      !
      write (mol,"('object 3 class array type int rank 1 shape 2"// &
                 " items ',i6,' data follows')") nbonds
      write (mol,"(2(1x,i7))") bonds(:,1:nbonds)-1
      write (mol,"('attribute ""element type"" string ""lines""')")
      write (mol,"('attribute ""ref"" string ""positions""')")
      !
      write (mol,"('object ""molecule"" class field')")
      write (mol,"('component ""data"" value 1')")
      write (mol,"('component ""positions"" value 2')")
      write (mol,"('component ""connections"" value 3')")
      write (mol,"('attribute ""name"" string ""graphene""')")
      !
      close (mol)
      !
      deallocate (atoms,bonds)
      return
    end subroutine plotGraphene

  end module small_box
!
  subroutine driver
    use small_box

!   write (out,"(' This program needs an adjustment to unmapped_base ')") 
!   pause 'Please adjust unmapped_base now!'
!  This is no longer necessary - use maxram wrapper to run the code!
    call nuc
  end subroutine driver
