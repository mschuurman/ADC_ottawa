!
!  Stationary solutions in a relatively small box - direct diagonalization. This
!  example is for a very special case of C60-intercalated graphene, which contains
!  excluded volume inside C60. As a result, the potential evaluation is a bit
!  strange. One of the limitations of the code is that it supports only orthorhombic
!  cells. Sorry.
!
  module small_box
    use accuracy
    use multigrid
    use lapack
    use timer
    implicit none
    !
    !  ==== Fixed parameters, changing which does not make any sense ====
    !
    integer(ik), parameter :: unit_save   = 34
    integer(ik), parameter :: unit_struct = 35
    integer(ik), parameter :: max_exclude = 2                 ! Max number of excluded volume 
                                                              ! elements.
    real(rk), parameter    :: bond_cc     = 1.421_rk / abohr  ! Experimental C=C bond in graphite
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik) :: verbose       = 2                           ! Verbosity level
    real(rk)    :: mass          = (1.0_rk+1838.7_rk) * 2.0_rk ! Particle mass, in e.m.u
    integer(ik) :: n_points(3)   = (/ 25, 26, 26 /)            ! Number of sampling points along each
    integer(ik) :: sum_cells(3)  = (/ 5, 5, 5 /)               ! Number of cells to include in potential
    integer(ik) :: plot_cells(3) = (/ 1, 1, 1 /)               ! cells to include in background image
    integer(ik) :: plot_states   = -1                          ! Number of states to plot
    !
    !  ==== Excluded spherical volume elements ====
    !
    integer(ik) :: n_exclude     = 0            ! Number of excluded VEs (spherical)
    real(rk)    :: xyzr_exclude(4,max_exclude)  ! X/Y/Z/R (Bohr) of excluded VEs
    !
    !  ==== Structure data ====
    !
    character(len=80)     :: coord_file = ' '   ! Name of the file containing structure
                                                ! (blank if inline)
    real(rk)              :: lat_vec(3)         ! Orthorhombic periodic lattice parameters (Bohr)
                                                ! Note that the unit cell MUST be centered at
                                                ! the origin - a lot of things won't work otherwise.
    real(rk)              :: lat_clip(3) = (/ 0._rk, 0._rk, 0._rk /)
                                                ! Symmetrically reduce the box from a given direction 
                                                ! by lat_clip(:) amount. The reduction affects only 
                                                ! the wavefunction, and should be used when we -know-
                                                ! that the potential in the excluded area is so large
                                                ! that the wavefunction -will- be negligible.
    integer(ik)           :: n_centres          ! Number of centres within the unit cell
    real(rk), allocatable :: centres(:,:)       ! Coordinates of the carbon atoms (Bohr)
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)              :: n_states        ! Total number of states
    complex(rk), allocatable :: h(:,:)          ! Hamiltonian matrix / eigenstates
    real(rk), allocatable    :: e(:)            ! Eigenvalues
    real(rk), allocatable    :: e_empty(:)      ! Eigenvalues - empty box
    complex(rk), allocatable :: heev_work(:)    ! Scratch space for LAPACK
    real(rk), allocatable    :: heev_rwork(:)
    !
    integer(ik)              :: pot, fft        ! Field handles for multigrid
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /smallbox/ verbose, &
                        mass, &
                        n_points, sum_cells, plot_cells, &
                        lat_vec, lat_clip,                      &
                        coord_file,                             &
                        n_exclude, xyzr_exclude
    !
    !  ==== End of global data ====
    !
    contains
    !   
    !  Pairwise potential for C-H2 interaction, where H2 is taken as a
    !  feature-less particle at the centre of mass of H2.
    !
    function vdw_potential(r) result(v)
      real(rk), intent(in) :: r(3)
      real(rk)             :: v
      !
      real(rk)             :: rx
      !
      rx = sqrt(sum(r**2)) * abohr  ! Bohr -> Angstroms
      rx = max(1.8_rk,rx)           ! Cut off hard repulsive tail
      !
      !  Potential C - fit to MP2/cc-pvtz in coronene
      ! 
      v = -400.37_rk/rx**6 + 25352.0_rk * exp(-3.5763_rk * rx)
      !
      v = v / ( 627.51_rk ) ! kcal/mol -> Hartree
    end function vdw_potential
    !
    !  Total interaction potential. Input coordinate is assumed to be
    !  inside the unit cell (no checking!)
    !
    function pot_create(x) result(vr)
      real(rk), intent(in) :: x(3)       ! Coordinates of the point
      real(rk)             :: vr         ! Interaction potential
      integer(ik)          :: ixc        ! Exclusion zone counter
      real(rk)             :: push_x(3)  ! Adjusted coordinates, pushed "outside"
                                         ! of exclusion zones if necessary
      real(rk)             :: dx(3)      ! Relative position wrt exclusion zone centre
      real(rk)             :: lat_off(3) ! Replicated cell offset
      real(rk)             :: rx         ! Displacement from the exclusion centre
      real(rk)             :: rs         ! Exclusion radius
      integer(ik)          :: ia, ib, ic ! Replication counters
      integer(ik)          :: iat        ! Centres within the unit cell
      !
      !  Check for the exclusion zones. The check here works only for non-overlapping
      !  exclusion zones - otherwise, there is no guarantee we'll clear all exclusion
      !  zones in the system.
      !
      push_x = fold_vector(x)
      !
      exclusion: do ixc=1,n_exclude
        dx = fold_vector(push_x - xyzr_exclude(1:3,ixc))
        !
        rs = xyzr_exclude(4,ixc)
        rx = sqrt(sum(dx**2))
        if ( rx>=rs ) cycle exclusion
        !
        !  Point is within the exclusion zone, push it out. If point happens to
        !  be at the origin, arbitraliry push it along the third coordinate ("up")
        !
        if ( rx>=0.01_rk ) then
          push_x = xyzr_exclude(1:3,ixc) + rs*dx/rx
        else
          push_x = xyzr_exclude(1:3,ixc) + (/ 0._rk, 0._rk, rs /)
        end if
        !
        !  Pushing could lead to a vector outside the unit cell, renormalize
        !
        push_x = fold_vector(push_x)
      end do exclusion
      !
      !  push_x(:) coordinate should now be either outside or on the boundary of
      !            the exclusion zones. 
      !            
      !  Sum the potential over enough surrounding zones
      !
      vr = 0
      sum_a: do ia=-sum_cells(1),sum_cells(1)        ! Replicate the cell in x direction
         sum_b: do ib=-sum_cells(2), sum_cells(2)    ! Replicate the cell in y direction
            sum_c: do ic=-sum_cells(3),sum_cells(3)  ! Replicate the cell in z direction
               lat_off = lat_vec * (/ ia, ib, ic /)
               do iat=1,n_centres
                  dx = push_x - (centres(:,iat) + lat_off)
                  vr = vr + vdw_potential(dx)
               end do
            end do sum_c
         end do sum_b
      end do sum_a
      !
    end function pot_create
    !
    !  Fold a vector within a zero-centered unit cell 
    !
    function fold_vector(xin) result(x)
      real(rk), intent(in) :: xin(3)  ! Arbitrary vector
      real(rk)             :: x  (3)  ! Vector folded within the origin cell
      !
      x = xin
      normalize_x: do while (any(abs(x)>lat_vec/2))
        where (x>lat_vec/2) 
          x = x - lat_vec
        end where
        where (x<-lat_vec/2)
          x = x + lat_vec
        end where
      end do normalize_x
    end function fold_vector
    ! 
    !  External interface neede by multigrid routines - must return a complex field
    !
    function interactionPotential(x) result(v)
      real(rk), intent(in)  :: x(3)
      complex(rk)           :: v
      !
      ! Debug - quadratic well
      !
      ! v = 0.20_rk * sum((x-(/0._rk,0._rk,3._rk/))**2)
      ! return
      v = pot_create(x)
      !
    end function interactionPotential
    !
    !  Read unit cell structure, formatted as an XYZ file.
    !
    subroutine read_xyz_structure
      integer(ik)           :: inp     ! I/O channel used for reading the structure data
      integer(ik)           :: iat     ! Atom number
      integer(ik)           :: ic      ! Lattice component index
      character(len=2)      :: element ! Discarded later
      !
      !  If the input file name was provided in the namelist, use it.
      !  Otherwise, continue reading the in-line input.
      !
      inp = input
      if (coord_file/=' ') then
        inp = unit_struct
        open(inp,file=trim(coord_file),status='old',action='read')
      end if
      !
      !  Read the XYZ file. The first line contains the number of atoms.
      !
      n_centres = 0
      read(inp,*) n_centres
      allocate(centres(3,n_centres))
      !
      ! Read the atoms. Coordinates are expected to be in Angstroms,
      ! do that we'll do the conversion later.
      !
      read_atoms: do iat=1,n_centres
        read(inp,*) element, centres(:,iat)
        if (any(abs(centres(:,iat))>lat_vec/2)) then
          write (out,"('Atom at coordinates ',3f14.7,' is outside unit cell')") &
            centres(:,iat)
          stop 'bad geometry'
        end if
      end do read_atoms
      !
      !  If the structure was coming from an external file, close it
      !
      if (coord_file/=' ') then
        close(inp,status='keep')
      end if
    end subroutine read_xyz_structure

    !
    !  Problem driver
    !
    subroutine nuc
      integer(ik)        :: root, info, plot_roots
      character(len=200) :: buf
      !
      call TimerStart('nuc')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk 
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=smallbox,iostat=info)
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=smallbox)
      write (out,"()")
      !
      !  Additional structure input - comes as an XYZ file. Note that the
      !  structure must be consistent with the positions and radii of the
      !  exclusion centres
      !
      call read_xyz_structure
      !
      !  Units conversion: Angstrom -> Bohr
      !
      lat_vec      = lat_vec      / abohr
      lat_clip     = lat_clip     / abohr
      xyzr_exclude = xyzr_exclude / abohr
      centres      = centres      / abohr
      !
      call plotStructure
      call buildGrid
      call addFields
      !
      !  For a very small field we have here, we'll simply diagonalize
      !  the Hamiltonian matrix. FFT transform of the potential will 
      !  give us the potential part. The kinetic part is diagonal.
      !
      call TimerStart('stationary problem')
      call FieldInit(pot,interactionPotential)
      call Visualize(pot,'Interaction potential')
      !
      write (out,"(' Done potential visualization')")
      call flush(out)
      !
      call buildHamiltonian
      !
      write (out,"(' Done Hamiltonian construction')")
      call flush(out)
      !
      write (out,"(/t5'Empty box eigenvalues are:'/)")
      write (out,"((t3,8(f12.6,1x)))") e_empty
      !
      call flush(out)
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
        write (out,"('Array sizes were: ',5i15)") 64*n_states, 3*n_states
        pause 'smallbox - diagonalization scratch'
      end if
      call cheev('V','U',n_states,h,n_states,e,heev_work,size(heev_work),heev_rwork,info)
      deallocate (heev_work,heev_rwork)
      if (info/=0) then
        write (out,"('cheev failed: ',i6)") info
        stop 'cheev failed'
      end if
      call flush(out)
      !
      call TimerStop('Diagonalization')
      call TimerStop('stationary problem')
      call TimerReport
      !
      write (out,"(/t5'Eigenvalues are:'/)")
      write (out,"((t3,8(f12.6,1x)))") e
      !
      !  Store all eigenvalues with lots of significant digits, for use in
      !  partition function evaluation.
      !
      call save_all_roots
      !
      !  Visualize some or all of the states
      !
      plot_roots = n_states
      if (plot_states>=0) plot_roots = min(plot_roots,plot_states)
      display_roots: do root=1,plot_roots
        write (buf,"(' Root = ',i5,' energy = ',f12.6)") root, e(root)
        write (out,"(a)") trim(buf)
        call buildRoot(root)
        call Visualize(pot,trim(buf))
      end do display_roots
      !
      !  Halleluya!
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
        write (out,"('Array sizes were: ',5i15)") n_states, n_states, n_states**2
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
      logical               :: inverse
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
        write (out,"('Array sizes were: ',5i15)") ox*oy*oz, ox*oy*oz
        pause 'build_potential - no memory'
      end if
      !
      !$omp parallel do private(iz,sz,iy,sy,ix,sx,x)
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
      !$omp end parallel do
      !
      inverse = .false.
      call fftw_3d(ox,oy,oz,over_pot,over_fft,inverse)
      over_fft = over_fft / product(over_npts)
      !
      deallocate (over_pot) 
      !
      call TimerStop('build_potential_FFT')
      !
    end subroutine build_potential_FFT

    subroutine stupid_sort(esort)
      real(rk), intent(inout) :: esort(:)
      !
      integer(ik) :: i
      logical     :: swap
      real(rk)    :: x
      !
      swap = .true.
      do while(swap)
        swap = .false.
        do i=2,size(e)
          if(esort(i-1)>esort(i)) then
            swap   = .true.
            x      = esort(i-1)
            esort(i-1) = esort(i)
            esort(i)   = x
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
      ! Create the unit box 
      !
      box(2,:) = 0.5_rk*lat_vec - lat_clip
      !
      if (any(box(2,:)<=0)) then
        write (out,"(' Encountered negative upper dimension for the box: '/t8,3f14.6)") &
               box(2,:)
        stop 'bad box - too much clipping?'
      end if
      !
      box(1,:) = -box(2,:)
      !
      write (out,"(//t5,'Simulation box size (Bohrs)'/)") 
      write (out,"(t8,'X: ',2f14.6)") box(:,1)
      write (out,"(t8,'Y: ',2f14.6)") box(:,2)
      write (out,"(t8,'Z: ',2f14.6)") box(:,3)
      write (out,"()")
      !
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      call TimerStop('Grid initialization')
    end subroutine buildGrid

    subroutine addFields
      call FieldNew('Potential (Real space)',     pot, scratch=.true. , wavefunction=.false.)
      call FieldNew('Potential (Momentum space)', fft, scratch=.true. , wavefunction=.false.)
    end subroutine addFields
    !
    !  Generate image of the structure, suitable for plotting with OpenDX
    !
    subroutine plotStructure
      integer(ik)              :: cell_a, cell_b, cell_at
      integer(ik)              :: iatom, jatom, natoms, ibond, nbonds,k,j,i,n
      real(rk)                 :: cell_d(3), rij 
      real(rk),allocatable     :: atoms (:,:)
      integer(ik), allocatable :: bonds(:,:)
      integer(ik)              :: mol
      !
      ! Count number of atoms in the replicated structure
      !
      natoms = n_centres*product(2*plot_cells+1)
      !
      allocate (atoms(3,natoms),bonds(2,(natoms*(natoms-1))/2))
      !
      ! Explicitly replicate atoms
      !
      iatom = 0
      !
      sum_a: do k = -plot_cells(1),plot_cells(1)       ! Replicate the cell in x direction
         sum_b: do j = -plot_cells(2),plot_cells(2)    ! Replicate the cell in y direction
            sum_c: do n = -plot_cells(3),plot_cells(3) ! replicate the cell in z direction
               do i = 1,n_centres                      ! sum over all atoms
                  iatom = iatom+1
                  atoms(:,iatom) = centres(:,i) + lat_vec * (/k,j,n/)
               end do
            end do sum_c
         end do sum_b
      end do sum_a

      if (iatom/=natoms) then
         write (out,"(' iatom = ',i5,' natoms = ',i5,'. Oops.')") iatom, natoms
         stop 'atom count error'
      end if
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
    end subroutine plotStructure
  end module small_box
  !
  subroutine driver
    use small_box
    
    call nuc
  end subroutine driver
