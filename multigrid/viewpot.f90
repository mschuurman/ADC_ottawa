!
!  Visualize potentials for the graphite system
!
  module v_small_box
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
    ! Box V - potential visualization
    !
      logical, parameter       :: layer_2    = .true.    ! Include second layer?
      real(rk), parameter      :: layer_2_z  = 10.0 / a0 ! Z-coordinate of the 2nd layer
      integer(ik), parameter   :: units_a    = 3         ! Replicated units in the X direction
      integer(ik), parameter   :: units_b    = 2         ! Replicated units in the Y direction
      real(rk), parameter      :: min_z      =  2.0 / a0 ! Min. distance for the potential
      real(rk), parameter      :: max_z      =  8.0 / a0 ! Max. distance for the potential
      integer(ik)              :: n_points_a = 17*8      ! Number of sampling points along each
      integer(ik)              :: n_points_b = 20*8      ! direction; 14x16x72 seems to work;
      integer(ik)              :: n_points_c = 24*8      ! 15x17x68 seems to work as well;
    !
      real(rk), parameter      :: test_kt    = 0.5_rk / 627.51_rk ! kt at 250K, in Hartree
    !
    integer(ik)                :: sum_cells  = 30        ! Number of cells to include in potential
                                                         ! summation
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
    integer(ik)              :: pot, fft, scr
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
      v = v + 1.0_rk
    end function interactionPotential
!
    complex(rk) function funex1(xyz,val)
      real(rk)    :: xyz(3)
      complex(rk) :: val
      !
      funex1 = exp(-(val-1.0_rk)/test_kt)
    end function funex1
!
    complex(rk) function funex2(xyz,val)
      real(rk)    :: xyz(3)
      complex(rk) :: val
      !
      funex2 = exp(-2*(val-1.0_rk)/test_kt)
    end function funex2
!
    subroutine nuc
      integer(ik)        :: root, info
      character(len=200) :: buf
      real(rk)           :: int1, int2
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
!
!    Compare integrals of exp(-v/kt) and exp(-2v/kt) - check for
!    non-ideal gas extrapolation correctness.
!
      call FieldCopy(pot,scr)
      call FieldProcess(scr,funex1)
      int1 = FieldNorm1(scr)
      call FieldCopy(pot,scr)
      call FieldProcess(scr,funex2)
      int2 = FieldNorm1(scr)
      write (out,"(' kt = ',g12.6,' exp(-v/kt) = ',g12.6,' sqrt(exp(-2v/kt)) = ',g12.6)") &
             test_kt, int1, sqrt(int2)
!
      call Visualize(pot,'')
!
      call TimerStop('nuc')
      call TimerReport
!
      wait_end: do 
        call sleep(20)
      end do wait_end
!
    end subroutine nuc

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
      call FieldNew('Scratch',                    scr, scratch=.false., wavefunction=.false.)
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

  end module v_small_box
!
  subroutine driver
    use v_small_box

!   write (out,"(' This program needs an adjustment to unmapped_base ')") 
!   pause 'Please adjust unmapped_base now!'
!  This is no longer necessary - use maxram wrapper to run the code!
    call nuc
  end subroutine driver
