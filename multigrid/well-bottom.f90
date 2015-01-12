!
!  Estimation of the potential well bottoms for exp-6 potential
!
  module well_bottom
    use accuracy
    implicit none
!
    !
    !  Bohr length
    !
    real(rk), parameter      :: a0       = 0.52918_rk
    !
    !  Carbon-carbon bond length, in Bohr
    !
    real(rk), parameter      :: bond_cc  = 1.421_rk / a0 ! Experimental value for graphite
    !
    integer(ik)              :: sum_cells  = 30          ! Number of cells to include in potential
                                                         ! summation
    logical, parameter       :: layer_2    = .true.      ! Include second layer?
    real(rk), parameter      :: layer_2_z  =  9.00/ a0   ! Z-coordinate of the 2nd layer
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
  end module well_bottom
!
  program driver
    use well_bottom

    integer(ik) :: i
    real(rk)    :: z, dz, minf, f
    !
    !  We know what the minimum is found above the ring
    !  centre - so this makes the search a lot easier.
    !
    z = 2.0_rk ; dz = 0.001_rk ; minf = 1000.0_rk
    scan: do i=1,25000
      z    = z + dz
      if (layer_2 .and. z>=layer_2_z) exit scan
      f    = interactionPotential( (/ bond_cc*sqrt3/2, 0.0_rk, z /) )
      minf = min(f,minf)
      if (mod(i,100)==1) then
        write (out,"(1x,f8.4,2x,f12.6)") z, f
      end if
    end do scan
    !
    write (out,"('Potential well bottom is at ',f13.8,' Hartree')") minf

  end program driver
