!
!  Stationary solutions in a relatively small box - direct
!  diagonalization. This example solves for the stationary 
!  states of hydrogen molecule adsorbed on a 2D periodic
!  graphene sheet.
!
  module big_brother
    use accuracy
    use multigrid
    use lapack
    use fftw
    use timer
    use lanczos
    use liu
    use symmetry
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
    !  Math constants
    !
    real(rk), parameter      :: sqrt3        = 1.73205080756887729352_rk
    !
    !  Cell and grid parameters
    !
    !  ========================================================================
    !
    !  Example 1 - section of dual-layer graphene surface, zero boundary
    !              conditions for the wavefunction.
    !
    ! logical, parameter       :: layer_2    = .true.    ! Include second layer?
    ! real(rk), parameter      :: layer_2_z  =  9.0 / a0 ! Z-coordinate of the 2nd layer
    ! integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
    ! integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
    ! real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
    ! real(rk), parameter      :: max_z      =  7.0 / a0 ! Max. distance for the potential
    ! real(rk), parameter      :: shift_image_z = 2.0/a0 ! Distance to shift the image by
    ! integer(ik)              :: n_points_a = 36        ! Number of sampling points along each
    ! integer(ik)              :: n_points_b = 42        ! direction; 
    ! integer(ik)              :: n_points_c = 40        ! 
    ! !
    ! integer, parameter       :: roots_count = 1        ! Number of roots we need
    ! !
    ! integer(ik)              :: sum_cells   = 20       ! Number of cells to include in potential
    !                                                    ! summation. 30 is a good number to reach
    !                                                    ! convergence, but we'll use 10 for testing
    ! integer(ik)              :: plot_cells  = 5        ! Number of cells to include in background
    !                                                    ! image plot
    ! integer(ik), parameter   :: n_scratch   = 500      ! Number of extra scratch fields,
    !                                                    ! used for solving
    !                                                    ! stationary problem
    ! !
    ! !  Definitions of the interacting centres in the unit cell
    ! !
    ! real(rk), parameter      :: centres(3,4) = reshape( &
    !   (/                                                &
    !      0.0_rk,            bond_cc/2.0,  0.0_rk,  &
    !      0.0_rk,           -bond_cc/2.0,  0.0_rk,  &
    !      bond_cc*sqrt3/2,   bond_cc,      0.0_rk,  &
    !      bond_cc*sqrt3/2,  -bond_cc,      0.0_rk   &
    !   /), (/ 3, 4 /) )
    ! !
    ! real(rk), parameter      :: simulation_box(2,3) = reshape( (/ &
    !  -bond_cc*sqrt3*units_a/2,  bond_cc*sqrt3*units_a/2,          &
    !  -bond_cc*3.0_rk*units_b/2, bond_cc*3.0_rk*units_b/2,         &
    !   min_z, max_z /), (/ 2, 3 /) )
    ! !
    ! real(rk), parameter      :: pot_cell_d(3) = (/ bond_cc*sqrt3, bond_cc*3, 0.0_rk /)
    !
    !  End of Example 1
    !  ========================================================================
    !  ========================================================================
    !
    !  Example 2 - "benzene" in a box.
    !
      logical, parameter       :: layer_2    = .false.   ! Include second layer?
      real(rk), parameter      :: layer_2_z  =  0.0      ! Z-coordinate of the 2nd layer
      real(rk), parameter      :: max_z      =  8.0 / a0 ! Max. distance for the potential (Z)
      real(rk), parameter      :: max_xy     = 10.0 / a0 ! Max. distance for the potential (X/Y)
      real(rk), parameter      :: shift_image_z = 0.0/a0 ! Distance to shift the image by
      integer(ik)              :: n_points_a = 80        ! Number of sampling points along each
      integer(ik)              :: n_points_b = 80        ! direction; 
      integer(ik)              :: n_points_c = 64        ! 
      !
      integer(ik), parameter   :: roots_count = 10       ! Number of roots we need
      !
      integer(ik)              :: sum_cells   = 0        ! This is a non-periodic potenital
      integer(ik)              :: plot_cells  = 0        ! Number of cells to include in background
                                                         ! image plot
      integer(ik), parameter   :: n_scratch   = 200      ! Number of extra scratch fields,
                                                         ! used for solving stationary problem
      !
      !  Definitions of the interacting centres in the unit cell
      !
      real(rk), parameter      :: centres(3,6) = reshape( (/ &
           bond_cc,     0.0_rk,          0.0_rk,  &
          -bond_cc,     0.0_rk,          0.0_rk,  &
           bond_cc/2,   bond_cc*sqrt3/2, 0.0_rk,  &
          -bond_cc/2,   bond_cc*sqrt3/2, 0.0_rk,  &
           bond_cc/2,  -bond_cc*sqrt3/2, 0.0_rk,  &
          -bond_cc/2,  -bond_cc*sqrt3/2, 0.0_rk   &
        /), (/ 3, 6 /) )
      !
      real(rk), parameter      :: simulation_box(2,3) = reshape( (/ &
       -max_xy, max_xy, -max_xy, max_xy, -max_z, +max_z /), (/ 2, 3 /) )
      !
      real(rk), parameter      :: pot_cell_d(3) = 2 * (/ max_xy, max_xy, max_z /)
    !
    !  End of Example 2
    !  ========================================================================
    !
    !  The rest of the stuff will be done automatically
    !
    integer(ik)              :: n_states                 ! Total number of states
    !
    integer(ik)              :: pot1, Hpsi, psi1, psi2
    integer(ik)              :: extra_scr(n_scratch)     ! Extra scratch fields
    !
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
      rx = max(rx,1.7_rk)        ! Clip potential close to the nuclei
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
          cell_d    = pot_cell_d * (/ cell_a, cell_b, 0 /)
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
    function guessFunction(x) result(psi)
      real(rk), intent(in) :: x(3)
      complex(rk)          :: psi
      !
      real(rk) :: pot
      !
      pot = interactionPotential(x)
      if (pot>=0) then
        psi = 0
      else
        psi = -pot
      end if
    end function guessFunction
!
    subroutine nuc
      integer(ik)        :: root, info
      character(len=200) :: buf
      real(rk)           :: eval(roots_count)
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
      call FieldInit(psi1,guessFunction)
      call Visualize(psi1,'Initial guess')
!
      call TimerStart('stationary problem')
      call solveStationaryProblem(eval)
      call TimerStop('stationary problem')
      call TimerReport
!
!     display_roots: do root=1,n_states
!       call sleep(20)
      display_roots: do 
        write(out,"('Choose state to visualize; 0 to exit',i5)")
        read(5,*) root
        if (root<=0) exit
        write (buf,"(' Root = ',i5,' energy = ',f12.6)") root, eval(root)
        write (out,"(a)") trim(buf)
        call LUeigenvector(root,psi1)
        call Visualize(psi1,trim(buf))
      end do display_roots
!
!    Release the solution
!
      call LUeigenvector(0,0)
!
      call TimerStop('nuc')
      call TimerReport
    end subroutine nuc

    subroutine Visualize(psi,title)
      integer(ik), intent(in)      :: psi       ! Wavefunction
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
      call FieldFFT(psi,psi2)
      call FieldVisualize(1,psi2,trim(title),1,efield)
      call FieldShow
      !
      call TimerStop('Visualization')
    end subroutine Visualize

    subroutine solveStationaryProblem(eval)
      real(rk), intent(out) :: eval(:)
      real(rk) :: wf_norm
!
!    Solve for the eigenvalues
!
      call SMsetSymmetry('Diatomic',main_axis_order=10, &
                         main_axis=(/0.0_rk,0.0_rk,1.0_rk/), &
                         inversion_centre=sum(simulation_box,dim=1)/2)
      call FieldSetOuterWall('REFLECTING')
      call LUeigenvalues(mass,(/psi1,psi2,pot1,Hpsi,extra_scr/),eval, &
                         interactionPotential,guessFunction,interactionPotential)
!
      write (out,"('Begin all eigenvalues:')") 
      write (out,"((t8,f30.15))") eval
      write (out,"('End all eigenvalues')")
!
!    Get the ground-state eigenvector
!
      call LUeigenvector(1,psi1)
      wf_norm = FieldNorm(psi1)
      write (out,"(' Norm of the ground state solution wavefunction = ',f15.10)") wf_norm
!
!    Release scratch space
!
!     call LUeigenvector(0,0)
!
    end subroutine solveStationaryProblem

    subroutine buildGrid
      !
      call TimerStart('Grid initialization')
      call MultiGridInit(max_grids=1,max_fields=4+n_scratch,nborder=1)
      !
      write (out,"(//t5,'Simulation box size (Bohrs)'/)") 
      write (out,"(t8,'X: ',2f14.6)") simulation_box(:,1)
      write (out,"(t8,'Y: ',2f14.6)") simulation_box(:,2)
      write (out,"(t8,'Z: ',2f14.6)") simulation_box(:,3)
      write (out,"()")
      !
      call SimpleGridNew('Rectangular box', &
             (/  n_points_a, n_points_b, n_points_c /), simulation_box)
      !
      call TimerStop('Grid initialization')
    end subroutine buildGrid

    subroutine addFields
      integer(ik)       :: i
      character(len=20) :: tag
      !
      call FieldNew('Nuclear Hamiltonian', pot1, scratch=.true. , wavefunction=.false.)
      call FieldNew(' |psi-1>',            psi1, scratch=.false., wavefunction=.true.)
      call FieldNew(' |psi-2>',            psi2, scratch=.false., wavefunction=.true.)
      call FieldNew('H|psi>',              Hpsi, scratch=.true. , wavefunction=.false.)
      !
      do i=1,n_scratch
        write (tag,"('W.F. scratch ',i5)") i
        call FieldNew(trim(tag), extra_scr(i), scratch=.true., wavefunction=.false.)
      end do
    end subroutine addFields

    subroutine plotGraphene
      integer(ik)              :: cell_a, cell_b, cell_at
      integer(ik)              :: iatom, jatom, natoms, ibond, nbonds
      real(rk)                 :: cell_d(3), rij
      real(rk),allocatable     :: atoms (:,:)
      integer(ik), allocatable :: bonds(:,:)
      integer(ik), parameter   :: mol = 56
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
          cell_d = pot_cell_d * (/ cell_a, cell_b, 0 /)
          scan_at: do cell_at=1,size(centres,dim=2)
            iatom = iatom + 1
            atoms(:,iatom) = centres(:,cell_at)+cell_d+(/0._rk,0._rk,shift_image_z/)
            if (layer_2) then
              iatom = iatom + 1
              atoms(:,iatom) = centres(:,cell_at)+cell_d + (/0._rk,0._rk,layer_2_z-shift_image_z/)
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

  end module big_brother
!
  subroutine driver
    use big_brother

    call nuc
  end subroutine driver
