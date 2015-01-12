!
!  (A partial) rotational correction to the adsorption potential of H2 molecules.
!  If a hydrogen molecule is adsorbed at a site where different molecular orientations
!  possess different energy, the rotational states of H2 are perturbed, causing a 
!  change in the rotational partition function and therefore in the adsorption
!  potential. The routines included here provide an estimation of the effect.
!
!  The following (somewhat crude) assumptions are made:
!  1. The rotational states of the molecule are in thermal equilibrium at any
!     given adsorption site, and are not affected by the nearby sites.
!  2. The anisotropy of the potential is caused by:
!     a. Interaction of the electric field gradient at the site with the
!        permanent quadrupole moment of the H2. This interaction creates
!        L=2 term in the angular potential.
!     b. Interaction of the electric field with the induced dipole of H2.
!        This interaction creates L=0 and L=2 terms in the angular potential.
!  3. The molecule is ^1H_2, with the rotational constants not perturbed by
!     the environemnt.
!  4. H2 molecules in the uniform H2 fluid are assumed to undergo free rotation.
!
!  Under these assumptions, the problem reduces to calculation of:
!
!   dv = -kT log(Zv/Zo)
!
!  where Zv is the rotational partition function in the anisotropic potential
!  with L=0 and L=2 components included, Zo is the rotational partition function
!  of the free rotor, k is the Boltmann constant, and T is the rotational
!  temperature.
!
!  The quadrupole moment is the recommended value of D.E. Stogryn and A.P. Stogryn,
!  Mol. Phys. 11, 371 (1966).
!
module h2_rotation
  use accuracy
  use math
  use lapack
  use printing
  use timer
  implicit none
  private
  public h2_vrot
  !
  integer(ik), parameter :: verbose =-1
  integer(ik), parameter :: h2_jmax = 9        ! The largest J value to consider
  real(rk), parameter    :: h2_Q    = 0.492_rk ! Quadrupole moment of H2 molecule, au[q]-Bohr^2
  real(rk), parameter    :: h2_apar = 6.831_rk ! Parallel contribution to the H2 polarizability, Bohr^3
  real(rk), parameter    :: h2_aper = 4.743_rk ! Perperndicular contribution ...
  real(rk), parameter    :: h2_B0   = 60.853_rk+3.062_rk/2 ! Rotational constant for \nu=0, in cm^-1
  !
  contains
  !
  function h2_vrot(t,ef,efg) result(vrot)
    real(rk), intent(in) :: t         ! Rotational temperature, in Kelvin
    real(rk), intent(in) :: ef (3)    ! Electric field, in Hartree atomic units
    real(rk), intent(in) :: efg(3,3)  ! Electric field gradient, in Hartree atomic units
    real(rk)             :: vrot      ! Rotational correction to the interaction potential
    !
    real(rk)    :: vcart(6)   ! Energies of H2 molecule along X, Y, Z, XY, XZ, and YZ directions
    complex(rk) :: vsph (6)   ! Spherical components of the potential.
    real(rk)    :: viso       ! Isotropic part - treat separately
    complex(rk) :: c2s  (6,6) ! Transformation matrix: Cartesian to spherical
    real(rk)    :: z0, zv     ! Logarithms of partition functions - free-space and perturbed
    !
    ! Can't use timing routines here - we may be called from a parallel region.
    ! call TimerStart('h2_vrot')
    !
    !  Begin by calculating the energies of the H2 molecule for six orientations
    !
    vcart(1) = e_h2(     1.0_rk ,     0.0_rk ,     0.0_rk )
    vcart(2) = e_h2(     0.0_rk ,     1.0_rk ,     0.0_rk )
    vcart(3) = e_h2(     0.0_rk ,     0.0_rk ,     1.0_rk )
    vcart(4) = e_h2(sqrt(0.5_rk),sqrt(0.5_rk),     0.0_rk )
    vcart(5) = e_h2(sqrt(0.5_rk),     0.0_rk ,sqrt(0.5_rk))
    vcart(6) = e_h2(     0.0_rk ,sqrt(0.5_rk),sqrt(0.5_rk))
    !
    !  Convert potential to spherical components
    !
    call fill_c2s(c2s)
    if (verbose>=3) then
      write (out,"(' c2s = '/6(6(1x,'(',f8.4,',',f8.4,')')/))") c2s
    end if
    vsph = matmul(c2s,vcart)
    if (verbose>=1) then
      write (out,"(' vcart = '/6(1x,g12.6))") vcart
      write (out,"(' vsph = '/6(2(1x,'(',g12.6,',',g12.6,')')))") vsph
    end if
    !
    !  Separate out isotropic part of the interaction
    !
    if (abs(aimag(vsph(1)))>spacing(100._rk)) then
      write (out,"('V00 is not real?! V00 = ',2g14.7)") vsph(1)
      stop 'h2_rotation - bad isotropic part'
    end if
    viso    = real(vsph(1),kind=rk) * 0.5_rk / sqrtpi ;
    vsph(1) = 0
    !
    !  Calculate rotational partition function of H2 for this potential
    !
    zv = rotational_z(t,vsph)
    !
    !  Calculate rotational partition function of H2 for zero potential
    !  Using the general routine is an overkill, but as long as this does 
    !  become too expensive, we could handle it.
    !
    vsph = 0
    z0 = rotational_z(t,vsph)
    !
    vrot = viso - k_Boltzmann*t*(zv-z0)
    !
    if (verbose>=0) then
      write (out,"(' ef  = ',3(1x,f14.7))") ef
      write (out,"(' efg = ',(2(t8,3(1x,f14.7)/),t8,3(1x,f14.7)))") efg
      write (out,"(' isotropic polarizability gives = ',f14.7,' K')") viso/k_Boltzmann
      write (out,"('            reorientation gives = ',f14.7,' K')") -t*(zv-z0)
    end if
    !
    ! call TimerStop('h2_vrot')
    return
  contains
    real(rk) function e_h2(dx,dy,dz)
      real(rk), intent(in) :: dx, dy, dz ! Direction cosines defining H2 axis
      !
      real(rk) :: rxy, rx, ry    ! 
      real(rk) :: u(3,3)         ! Rotation matrix: R_lab = U R_mol
      real(rk) :: a(3,3), am(3,3)! Polarizability tensor - lab and molecular
      real(rk) :: q(3,3), qm(3,3)! Quadrupole tensor - lab and molecular
      !
      !  Rotational matrix. This is just one of many possible choices, which are
      !  all equivalent for a linear molecule. It is obtained by choosing Euler
      !  angle gamma to be zero, so that the molecular Y axis is chosen in the
      !  lab XY plane.
      !
      rxy = sqrt(dx**2+dy**2)
      if (abs(rxy)>=spacing(10._rk)) then
        rx = dx / rxy
        ry = dy / rxy
      else  ! The alpha Euler angle is arbitrary for this orientation
        rx = 1._rk
        ry = 0._rk
      end if
      u(1,:) = (/ rx*dz,   -ry, dx /)
      u(2,:) = (/ ry*dz,    rx, dy /)
      u(3,:) = (/ -rxy,  0._rk, dz /)
      !
      if (verbose>=3) then 
        write (out,"(' dr = ',3(2x,f16.10))") dx, dy, dz
        write (out,"(' u(m2l) = ',(3(/3(2x,f16.10))))") transpose(u)
        write (out,"(' u . u^t = ',(3(/3(2x,f16.10))))") matmul(u,transpose(u))
      end if
      !
      !  Polarizability tensor in this coordinate system:
      !
      am(1,:) = (/ h2_aper, 0._rk,   0._rk   /)
      am(2,:) = (/ 0._rk,   h2_aper, 0._rk   /)
      am(3,:) = (/ 0._rk,   0._rk,   h2_apar /)
      a = matmul(u,matmul(am,transpose(u)))
      !
      !  Quadrupole tensor
      !
      qm(1,:) = (/ -h2_Q/2, 0._rk,   0._rk /)
      qm(2,:) = (/ 0._rk,   -h2_Q/2, 0._rk /)
      qm(3,:) = (/ 0._rk,   0._rk,   h2_Q  /)
      q = matmul(u,matmul(qm,transpose(u)))
      !
      e_h2 = -0.5_rk * dot_product(ef,matmul(a,ef)) &
           + (1.0_rk/3.0_rk) * sum(efg*q)
      !
      if (verbose>=3) then
        write (out,"(' alab = ',(3(/3(2x,f16.10))))") transpose(a)
        write (out,"(' qlab = ',(3(/3(2x,f16.10))))") transpose(q)
        write (out,"(' e_h2 = ',g16.10)") e_h2
      end if
    end function e_h2
    !
    subroutine fill_c2s(c2s)
      complex(rk), intent(out) :: c2s(6,6)
      complex(rk)              :: c0, c1, c2r, c2i, c3
      !
      c0  = cmplx( 0._rk                      ,   0._rk                      , kind=rk )
      c1  = cmplx( sqrtpi * 2._rk / 3._rk     ,   0._rk                      , kind=rk )
      c2r = cmplx( sqrtpi * sqrt(2._rk/15._rk),   0._rk                      , kind=rk )
      c2i = cmplx( 0._rk                      ,   sqrtpi * sqrt(2._rk/15._rk), kind=rk )
      c3  = cmplx( sqrtpi * sqrt(4._rk/45._rk),   0._rk                      , kind=rk )
      !
      c2s(1,:) = (/ c1,       c1,       c1,       c0,     c0,    c0    /)
      c2s(2,:) = (/ c2r-c2i, -c2r-c2i,  c0,       2*c2i,  c0,    c0    /)
      c2s(3,:) = (/ -c2r,    -c2i,     -c2r-c2i,  c0,     2*c2r, 2*c2i /)
      c2s(4,:) = (/ -c3,     -c3,       2*c3,     c0,     c0,    c0    /)
      c2s(5,:) = (/  c2r,    -c2i,      c2r-c2i,  c0,    -2*c2r, 2*c2i /)
      c2s(6,:) = (/ c2r+c2i, -c2r+c2i,  c0,      -2*c2i,  c0,    c0    /)
    end subroutine fill_c2s
  end function h2_vrot
  !
  function rotational_z(t,vsph) result(z)
    real(rk), intent(in)    :: t       ! Rotational temperature
    complex(rk), intent(in) :: vsph(6) ! Components of the rotational potential, in Hartrees
    real(rk)                :: z       ! Logarithm of the partition function
    integer(ik)             :: nstates ! Number of rotational states (J,M)
    !
    real(rk) :: z_even, z_odd ! Logarithms of even and odd contributions
    !
    !  We'll split the rotational problem into two parts - even and odd J values.
    !  Due to the nuclear statistical factor, even and odd J values can't mix.
    !
    !  Even J - each solution is non-degenerate
    !
    nstates = (1+h2_jmax/2)*(1+2*(h2_jmax/2))
    z_even = rotational_z_block(0_ik,h2_jmax,2_ik,t,vsph,h2_B0,nstates)
    !
    !  Odd J - each solution is triply degenerate
    !
    nstates = (1+(h2_jmax-1)/2)*(3+2*((h2_jmax-1)/2))
    z_odd = log(3._rk) + rotational_z_block(1_ik,h2_jmax,2_ik,t,vsph,h2_B0,nstates)
    !
    z = MathLogSum(z_even,z_odd)
    if (verbose>=2) then
      write (out,"(' Final log(Z)= ',g14.7)") z
    end if
  end function rotational_z
  !
  function rotational_z_block(jmin,jmax,dj,t,v,b0,nstates) result(z)
    integer(ik), intent(in) :: jmin, jmax, dj ! Indices of the coupled J block
    real(rk), intent(in)    :: t              ! Temperature in Kelvin
    complex(rk), intent(in) :: v(6)           ! J=0 and J=2 components of the potential, in Hartree
    real(rk), intent(in)    :: b0             ! Rotational constant, in cm^-1
    integer(ik), intent(in) :: nstates        ! Number of states, must be consistent with jmin/jmax/dj
    real(rk)                :: z              ! Logarithm of the partition function, assuming no degeneracies.
    !
    complex(rk), allocatable :: h(:,:)      ! The rotational Hamiltonian matrix
    real(rk), allocatable    :: e(:)        ! Eigenvalues
    integer(ik)              :: jl, ml, sl  ! State indices - (j,m) and compound
    integer(ik)              :: jr, mr, sr
    integer(ik)              :: ir
    !
    allocate (h(nstates,nstates),e(nstates))
    !
    sr = 0
    j_right: do jr=jmin,jmax,dj
      m_right: do mr=-jr,jr
        sr = sr + 1
        sl = 0
        j_left: do jl=jmin,jmax,dj
          m_left: do ml=-jl,jl
            sl = sl + 1
            h(sl,sr) = rotational_matrix_element(jl,ml,jr,mr,v,b0)
          end do m_left
        end do j_left
        if (sl/=nstates) stop 'count error (left) in rotational_z_block'
      end do m_right
    end do j_right
    if (sr/=nstates) stop 'count error (right) in rotational_z_block'
    !
    if (verbose>=3) then
      write (out,"(' Rotational J block: ',3i3)") jmin, jmax, dj
      call print_matrix(h,12_ik,'g12.6')
      write (out,"(' Deviation from hermitian: ',g12.6)") maxval(abs(h-transpose(conjg(h))))
    end if
    h = 0.5_rk * (h + transpose(conjg(h)))
    !
    !  Diagonalize the perturbed solution
    !
    call lapack_heev(h,e)
    !
    if (verbose>=2) then
      write (out,"(' Eigenvalues in cm^-1')")
      write (out,"(10(1x,g12.6))") h2cm * e
    end if
    !
    !  Calculate partition function. This is a bit complicated since we
    !  need to work with logarithms.
    !
    !  z = sum( exp(-e/(k_Boltzmann*t)) )
    !
    z = -e(1)/(k_Boltzmann*t)
    add_z_terms: do ir=2,nstates
      z = MathLogSum(z,-e(ir)/(k_Boltzmann*t))
    end do add_z_terms
    !
    if (verbose>=2) then
      write (out,"(' log of partial partition function = ',g20.10)") z
    end if
    !
    deallocate (h,e)
  end function rotational_z_block
  !
  function rotational_matrix_element(jl,ml,jr,mr,v,b0) result(hlr)
    integer(ik), intent(in) :: jl,ml ! "Left" state indices
    integer(ik), intent(in) :: jr,mr ! "Right" state indices
    complex(rk), intent(in) :: v(6)  ! J=0 and J=2 potential, Hartree
    real(rk), intent(in)    :: b0    ! Rotational constant, in cm^-1
    complex(rk)             :: hlr   ! Matrix element
    complex(rk)             :: term  ! Current contribution
    integer(ik)             :: mv
    !
    !  Quick exit if matrix element is guaranteed to be zero
    !
    hlr = 0._rk
    if (abs(jl-jr)>2 .or. abs(mr-ml)>2) return
    !
    if (jl==jr .and. mr==ml) then
      term = b0*jl*(jl+1)/h2cm ! b0 is in cm^-1, convert it
      hlr  = hlr + term
      ! write (out,"('diag <',i2,',',i3,'|',i2,',',i3,'> += ',2g12.6)") jl, ml, jr, mr, term
    end if
    !
    !  L=0 term in the potential
    !
    term = v(1) * (-1)**ml                   &
         * sqrt((2*jl+1)*(2*jr+1)/fourpi)    &
         * Math3J(jl,0_ik,jr, -ml,0_ik,  mr) &
         * Math3J(jl,0_ik,jr,0_ik,0_ik,0_ik)
    hlr = hlr + term
    ! write (out,"('v00  <',i2,',',i3,'|',i2,',',i3,'> += ',2g12.6)") jl, ml, jr, mr, term
    !
    !  L=2 terms in the potential
    !
    vL2_terms: do mv=-2,2,1
      term = v(4+mv) * (-1)**ml                &
           * sqrt(5*(2*jl+1)*(2*jr+1)/fourpi)  &
           * Math3J(jl,2_ik,jr, -ml,  mv,  mr) &
           * Math3J(jl,2_ik,jr,0_ik,0_ik,0_ik)
      hlr = hlr + term
      ! write (out,*) ' v  = ', v(4+mv)
      ! write (out,*) ' sq = ', sqrt(5*(2*jl+1)*(2*jr+1)/fourpi)
      ! write (out,*) ' 3j1 = ',Math3J(jl,2_ik,jr, -ml,  mv,  mr)
      ! write (out,*) ' 3j2 = ',Math3J(jl,2_ik,jr,0_ik,0_ik,0_ik)
      ! write (out,"('v2',i2,' <',i2,',',i3,'|',i2,',',i3,'> += ',2g12.6)") mv, jl, ml, jr, mr, term
    end do vL2_terms
  end function rotational_matrix_element
end module h2_rotation
!
!* program test_rot
!*   use accuracy
!*   use math
!*   use timer
!*   use h2_rotation
!*   real(rk) :: t, ef(3), efg(3,3), vrot
!*   !
!*   call accuracyInitialize
!*   read(*,*) t, ef, efg
!*   vrot = h2_vrot(t,ef,efg)
!*   write (*,*) 'vrot = ',vrot
!*   !
!*   call TimerReport
!* end program test_rot
