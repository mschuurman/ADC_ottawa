!
!  Various useful external fields and functions
!
  module fields
    use accuracy
    use math
    use cubeint
    implicit none
    private
    public NucleiT, FLsetNuclei, FLsetGridParameters, FLnuclearPotential, FLpsiGuess
    public FLnuclearPotentialOld
    public FLnuclearMultipoles
    public FLsetField, FLelectricField
    public FLsetGauss, FLgauss, FLrandom
    public FLpotentialAttractionForces, FieldSetNuclearDimension
    public FLpotentialRepulsionForces
    public FLpsiGaussPacket
    public FLsetMaskZero, FLmaskZero
    public FLsetPlanewave, FLplanewave, FLplaneWaveWithPhase, FLplaneWaveWithMinusPhase
    public FLexactPhase, FLsetMultipolesPhase
    public FLsetLREmultipoles, FLlrePotential
    public FLtotalPotential, FLtotalPotentialReal
    public FLasymptoticSetParameters, FLasymptotic
    public FLharmonicsSetParameters, FLharmonics
!
!  Data fields
!
    type NucleiT
      real(rk) ::  xyz(3)       ! Coordinates of a nucleus
      real(rk) ::  velocity(3)  ! Velocity of a nucleus
      real(rk) ::  charge       ! Charge of a nucleus
      real(rk) ::  width        ! "Nucleus size"
      real(rk) ::  nelec        ! Number of electrons (for the starting guess)
      real(rk) ::  mass         ! Mass of a nucleus
    end type NucleiT
!
!  Parameters for the nuclear Hamiltonian
!
   type(NucleiT), allocatable, save :: nuc(:)
!
!  Coulomb Hamiltonian can also include a second-order correction to the
!  potential, arising from finite patch dimension. To activate the correction,
!  it is necessary to call FLsetGridParameters
!
   real(rk), save :: grid_dx(3) = (/ 0._rk, 0._rk, 0._rk /)
!
!  Parameters for the electric field Hamiltonian
!
   real(rk), save :: efield(3)
!
!  Parameters for the planewave
!
   real(rk), save :: planek(3)
!
!  Parameters for the multipoles used in eikonal boundary conditions
!  Zero-th element is the charge of ion, followed by dipole moment
!  and traceless quadrupole moment.
!
   real(rk), save :: multipoles_phase(0:9)
!
!  Parameters for the multipoles used in long-range potential expansion
!  See FLsetLREmultipoles and FLlrePotential below. Stone's conventions
!  are used for the multipoles
!
   complex(rk), save :: lre_c        ! Charge
   complex(rk), save :: lre_d(3)     ! Dipole
   complex(rk), save :: lre_q(3,3)   ! Quadrupole
   real(rk), save    :: lre_width    ! Distribution width, used to clamp the singularities
!
!  Parameters for the current nuclei (1) and dimension (2) at the field force calcs
!
   integer(ik), save :: NuclDim(2)
!
!  Parameters for Gaussian convolution filter
!
   real(rk), save :: gaussAlpha
!
!  Parameters for masking function
!
   real(rk), save :: mask_width, mask_step
   logical, save  :: mask_invert
!
!  Parameters for the asymptotic radial expansion used in MO-ADK
!
       real(rk), save    :: as_z, as_kappa, as_r0(3)
!
!  Parameters for the spherical harmonics expansion used in MO-ADK
!
    real(rk), save    :: harm_rmin, harm_rmax, harm_r0(3)
    integer(ik), save :: harm_l, harm_m
    complex(rk), save :: harm_fact
!
!  Routines
!
    contains
    !
    !  Sets parameters of a masking function
    !
    subroutine FLsetMaskZero(width,step,invert)
      real(rk), intent(in)          :: width  ! Characteristic size of the mask
      real(rk), intent(in)          :: step   ! Characteristic decay width of the mask
      logical, intent(in), optional :: invert ! Invert the mask if true
      !
      mask_width  = width
      mask_step   = step
      mask_invert = .false.
      if (present(invert)) mask_invert = invert
    end subroutine FLsetMaskZero
    !
    function FLmaskZero(coord) result(v)
      real(rk), intent(in) :: coord(3) ! Coordinates
      complex(rk)          :: v        ! Mask
      !
      real(rk) :: r, w, we
      !
      r   = sqrt(sum(coord**2))
      !
      we  = (r-mask_width)/mask_step
      if (we>max_exp) then
        w = 0._rk
      else if (we<-max_exp) then
        w = 1._rk
      else
        w = 1.0_rk/(1.0_rk + exp(we))
      end if
      !
      if (mask_invert) w = 1.0_rk - w
      v = w
    end function FLmaskZero
    !
    !  Sets positions, charges, and velocities of nuclei
    !
    subroutine FLsetNuclei(nuclei)
      type(NucleiT), intent(in) :: nuclei(:)

      integer(ik) :: nnuc, alloc

      nnuc = size(nuclei)
      if (allocated(nuc)) deallocate (nuc)
      allocate (nuc(nnuc),stat=alloc)
      if (alloc/=0) then
        write (out,"(' Error ',i8,' trying to allocate ',i8,' nuclei ')") alloc, nnuc
        stop 'FLsetNuclei - out of memory'
      end if
      nuc = nuclei
    end subroutine FLsetNuclei
    !
    !  Evaluates multipole moments of nuclear point charges.
    !  This is the point charge counterpart of FieldNormMultipoles
    !
    !  The order of the multipoles is:
    !    1 = X ; 2 = Y ; 3 = Z ;
    !    4 = XX ; 5 = YY ; 6 = ZZ ; 7 = XY ; 8 = XZ ; 9 = YZ
    !
    subroutine FLnuclearMultipoles(mult)
      real(rk), intent(out) :: mult(:)
      !
      integer(ik)           :: inuc
      !
      mult = 0
      !
      !  Dipole terms
      !
      if (size(mult)<3) return
      dip: do inuc=1,size(nuc)
        mult(1:3) = mult(1:3) + nuc(inuc)%charge * nuc(inuc)%xyz
      end do dip
      !
      !  Diagonal quadrupole terms
      !
      if (size(mult)<6) return
      quad_diag: do inuc=1,size(nuc)
        mult(4:6) = mult(4:6) + nuc(inuc)%charge * nuc(inuc)%xyz**2
      end do quad_diag
      !
      !  Off-diagonal quadrupole terms
      !
      quad_off: do inuc=1,size(nuc)
        mult(7:9) = mult(7:9) + nuc(inuc)%charge * &
                  (/ nuc(inuc)%xyz(1)*nuc(inuc)%xyz(2), &
                     nuc(inuc)%xyz(1)*nuc(inuc)%xyz(3), &
                     nuc(inuc)%xyz(2)*nuc(inuc)%xyz(3) /)
      end do quad_off
    end subroutine FLnuclearMultipoles
    !
    !  Set grid parameters for 2-order nuclear potential correction
    !
    subroutine FLsetGridParameters(dx)
      real(rk), intent(in) :: dx(3)  ! Grid spacing
      !
      grid_dx = dx
    end subroutine FLsetGridParameters
    !
    !  Nuclear potential, including the "exact" volume element integration.
    !  Hopefully, this should reduce the egg-box effect in our simulations.
    !  This routine _requires_ grid_dx to be set.
    !
    function FLnuclearPotential(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v
      !
      real(rk)    :: tv, d(3), q
      integer(ik) :: inuc
      !
      ! v = FLnuclearPotentialOld(r) ; return
      if (any(grid_dx<=0)) stop 'fields%FLnuclearPotential - called before FLsetGridParameters'
      !
      tv = 0
      sum_nuclei: do inuc=1,size(nuc)
        d   = r - nuc(inuc)%xyz
        q   = nuc(inuc)%charge
        tv  = tv - q * cube_potential(d,grid_dx)
      end do sum_nuclei
      v = tv
    ! write (out,"('At point: ',3f14.7,' new potential: ',g25.15,' old: ',g25.15,' diff: ',g25.15)") &
    !        r, tv, real(FLnuclearPotentialOld(r),kind=rk), tv - real(FLnuclearPotentialOld(r),kind=rk)
    end function FLnuclearPotential

    !
    !  The width parameter is not arbitrary - it is related to
    !  the smallest grid spacing. Our basis functions are, in fact,
    !  flat 3D patches. Therefore, all potentials should be given
    !  by integrals over these parches. This does not matter too
    !  much for the external electric or magnetic fields. However,
    !  the nuclear Coulomb potential is singular, and would behave
    !  qualitatively wrong unless it is capped in some way for patches
    !  close to the nucleus.
    !
    !  One possible capping is given by integrating a sphere with the
    !  same volume as the surface patch, centered at the origin. The
    !  average potential in this case is 1.5/R, where R is the sphere
    !  radius. For rectangular grids, taking R as the geometric mean
    !  of the smallest grid spacing should do the trick. With our
    !  choice of the capped potential, this means that:
    !
    !    w = (4/3) (ax*ay*az)**(1/3)
    !
    !  In practive, we get better energies for:
    !
    !    w = (2/3) (ax*ay*az)**(1/3)
    !
    !  but the differences are pretty minor.
    !
    function FLnuclearPotentialOld(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v

      real(rk)    :: tv, d(3), dr(3), r12, w, q, c2f
      integer(ik) :: inuc, ic

      tv = 0
      do inuc=1,size(nuc)
        d   = r - nuc(inuc)%xyz
        r12 = sqrt(sum(d**2))
        q   = nuc(inuc)%charge
        w   = nuc(inuc)%width
        if (r12 <= w) then
          tv  = tv - 2.0_rk * q/(r12+w)
        else
          tv  = tv - q/r12
          !
          !  Outside of the finite distribution, also consider the 2nd-order
          !  correction term to the potential
          !
          if (any(grid_dx>0)) then
            dr  = d/r12
            c2f = 0
            corr_loop: do ic=1,3
               c2f = c2f + grid_dx(ic)**2 * (3*dr(ic)**2-1._rk)
            end do corr_loop
            tv = tv - q*c2f/(24._rk*r12**3)
          end if
        end if
      end do
      v = tv
    end function FLnuclearPotentialOld

    function FLpsiGuess(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v

      real(rk)    :: tv, d(3), r12, ne, q
      integer(ik) :: inuc

      tv = 0
      do inuc=1,size(nuc)
        d   = r - nuc(inuc)%xyz
        r12 = sqrt(sum(d**2))
        q   = nuc(inuc)%charge
        ne  = nuc(inuc)%nelec
        if (q*r12<max_exp) then
          tv  = tv + ne*exp(-q*r12)/sqrt(3.14159265358_rk)
        end if
      end do
      v = tv
    end function FLpsiGuess
!
    subroutine FLsetField(ef)
      real(rk), intent(in) :: ef(3)

      efield = ef
    end subroutine FLsetField
!
!  For the moment, we'll ignore retardation effects - these could be
!  important for really large boxes?
!
    function FLelectricField(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v

      v =  dot_product(r,efield)  ! This plus comes from a product of -1 in E = - \Nabla v
                                  ! and the -1 charge of the electron.
    end function FLelectricField
!
    subroutine FLsetPlanewave(k)
      real(rk), intent(in) :: k(3)

      planek = k
    end subroutine FLsetPlanewave
!
!    This routien expects trace-less multipole moments
!
    subroutine FLsetMultipolesPhase(mult)
      real(rk), intent(in) :: mult(0:9)

      multipoles_phase = mult
    end subroutine FLsetMultipolesPhase
!
    function FLplanewave(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v
      real(rk)             :: arg

      arg = dot_product(planek,r)
      v = cmplx(cos(arg),sin(arg),kind=rk)
    end function FLplanewave
!
!  Asymptotic scattered wave solution for quadrupole potential
!
    function FLexactphase(r) result (v)
      real(rk), intent(in) :: r(3)

      complex(rk)          :: v
      real(rk)             :: arg
      real(rk)             :: Q(9)
      real(rk)             :: absr
      real(rk)             :: absk
      real(rk)             :: arglog
      real(rk)             :: G0                 ! The eikonal Coulomb phase, logarithm
      real(rk)             :: G1                 ! The eikonal Coulomb phase, dipole
      real(rk)             :: G2                 ! The eikonal Coulomb phase, quadrupole

      absk = sqrt(sum(planek**2))
      absr = sqrt(sum(r**2))
      arg  = dot_product(planek,r)
      arglog = arg+absk*absr
      if (arglog<=0) then
        ! write (out,"('ARGH: r = ',3f14.7,' arglog = ',g20.12)") r, arglog
        arglog = spacing(1._rk)
      end if
      G0   = multipoles_phase(0)*log(arglog)/absk
      !
      ! First derivatives of the phase G0
      !
      Q(1:3)  = (r/absr+planek/absk)/arglog
      !
      ! Second derivatives, diagonal
      !
      Q(4:6) = (1/absr-r(1:3)**2/absr**3-(planek(1:3)+absk*r(1:3)/absr)**2/(arglog*absk))/arglog
      !
      ! Off-diagonal terms
      ! G_xy
      Q(7) = (-r(1)*r(2)/absr**3-(planek(1)+absk*r(1)/absr)*(planek(2)+absk*r(2)/absr)/(arglog*absk))/arglog
      ! G_xz
      Q(8) = (-r(1)*r(3)/absr**3-(planek(1)+absk*r(1)/absr)*(planek(3)+absk*r(3)/absr)/(arglog*absk))/arglog
      ! G_yz
      Q(9) = (-r(2)*r(3)/absr**3-(planek(2)+absk*r(2)/absr)*(planek(3)+absk*r(3)/absr)/(arglog*absk))/arglog
      !
      G1 = dot_product(multipoles_phase(1:3),Q(1:3))
      G2 = sum(multipoles_phase(4:6)*Q(4:6))+2*sum(multipoles_phase(7:9)*Q(7:9))  ! simple way of calculating Tr[mult(4:9)*Q]
      v  = G0-G1+G2/3
    end function FLexactphase
!
!  Construct phase-corrected plane wave
!
    function FLplanewavewithphase(r,phase) result (v)
      real(rk), intent(in)    :: r(3)
      complex(rk), intent(in) :: phase
      complex(rk)             :: v
      real(rk)                :: arg

      arg = dot_product(planek,r)+real(phase,kind=rk)
      v   = cmplx(cos(arg),sin(arg),kind=rk)
    end function FLplanewavewithphase
!
!  Construct phase-corrected plane wave with a minus in the phase
!
    function FLplanewavewithminusphase(r,phase) result (v)
      real(rk), intent(in)    :: r(3)
      complex(rk), intent(in) :: phase
      complex(rk)             :: v
      real(rk)                :: arg

      arg = dot_product(planek,r)+real(phase,kind=rk)
      v   = cmplx(cos(-arg),sin(-arg),kind=rk)
    end function FLplanewavewithminusphase
!
!  Set half-width at half-height of the gaussian at the coordinate origin.
!
    subroutine FLsetGauss(width)
      real(rk), intent(in) :: width

      gaussAlpha = log(2.0_rk)/width**2
    end subroutine FLsetGauss
!
!  Evaluate Gaussian at coordinate origin, with specified half-width.
!
    function FLgauss(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v

      real(rk)             :: x2

      x2 = gaussAlpha * sum(r**2)
      if (x2<max_exp) then
        v = exp(-x2)
      else
        v = 0
      end if
    end function FLgauss
!
!   Nuclear forces:
!   a) electron-nuclei interaction
!           F_i_x  = - (  q_i/r_ie^3 )*(x_i-x_e) );
!   b) nuclei-nuclei interaction:
!           F_i_x  = - \sum_{j<>i} ( -  q_i*q_j/r_ij^3 )*(x_i-x_j) ).
!
      function FLpotentialAttractionForces(r) result (v)
      real(rk)   , intent(in)          :: r(3)
      integer(ik)                      :: inuc
      integer(ik)                      :: axis
      real(rk)                         :: v
      real(rk)                         :: tv_ne
      real(rk)                         :: d(3),r_ne, w, q

      stop 'FLpotentialAttractionForces - implementation inconsistent with potential. Fix.'
      inuc = NuclDim(1)
      axis = NuclDim(2)
      !
      ! a) Coulombic attraction
      !
      d    = nuc(inuc)%xyz - r
      r_ne = sqrt(sum(d**2))
      q    = nuc(inuc)%charge
      w    = nuc(inuc)%width
      if (r_ne <= w) then
        tv_ne  =  2.0_rk * q/(r_ne+w)**3 * d(axis)
      else
        tv_ne  =  q/r_ne**3 * d(axis)
      end if

      v = - tv_ne

    end function FLpotentialAttractionForces

    function FLpotentialRepulsionForces() result (v)
      integer(ik)                      :: inuc,in
      integer(ik)                      :: axis
      real(rk)                         :: v
      real(rk)                         :: tv12
      real(rk)                         :: d(3), r12, w, q12

      inuc = NuclDim(1)
      axis = NuclDim(2)

      !
      ! b) Coulombic repulsion
      !
      tv12 = 0.0_rk
      w    = 0.1_rk
      do in = 1,size(nuc)
        if (in/=inuc) then
          d   = nuc(inuc)%xyz - nuc(in)%xyz
          r12 = sqrt(sum(d**2))
          q12 = nuc(in)%charge*nuc(inuc)%charge
          if (r12 <= w) then
            write(out,"('Your nuclear dynamics just knocked two nuclei:"// &
                      " ',I2,'and',I2)") inuc, in
            stop 'Something wrong with it!'
          endif
          tv12 =tv12 + q12/(r12**3) * d(axis)
        endif
      end do

      v = tv12

    end function FLpotentialRepulsionForces
    !
    ! Set and save  the current nuclear and the current dimension
    ! to calculate nuclear forces for
    !
    subroutine FieldSetNuclearDimension(inuc,axis)
      integer(ik), intent(in) :: inuc,axis

      NuclDim(1) = inuc
      NuclDim(2) = axis
    end subroutine FieldSetNuclearDimension
    !
    ! Gaussian Wave Packet:  N*exp( i*p((r-r0) )*exp(-a*(r-r0)^2)
    !   p/m - the group velocity
    !   r0  - the center of the packet
    !   a   - (a measure of) the velocity dispersion
    ! There is one parameter for each dimension
    !
    ! The parameters of the packet are hard-coded!
    !
    function FLpsiGaussPacket(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v
      complex(rk)          :: cv
      !
      real(rk) :: tv, r0(3), p(3), a(3), Norm , pr
      !
!     p  = (/ 0.00_rk, 0.00_rk, -pi*0.50_rk /)
      p  = (/ 0.00_rk, 0.00_rk, -pi*0.25_rk /)
      a  = (/ 0.01_rk, 0.01_rk,     0.01_rk /)
      r0 = (/ 0.00_rk, 0.00_rk,    15.00_rk /)
      !
      ! Norm factor is from Patchkovsii's manuscript
      !
      Norm = sqrt(sqrt(2.0_rk/pi))**3*sqrt(sqrt(a(1)*a(2)*a(3)))

      pr = sum( p(:)*( r(:)-r0(:) ) )

      tv = exp( -( a(1)*( r(1)-r0(1) )**2+a(2)*( r(2)-r0(2) )**2+a(3)*( r(3)-r0(3) )**2 ) )
      cv = cmplx(cos(pr),sin(pr),kind=rk)

      v = Norm*tv*cv

    end function FLpsiGaussPacket
    !
    !  Process a field - scale the numbers by a random factors
    !
    function FLrandom(coord,val) result (res)
      real(rk), intent(in)    :: coord(3) ! Coordinate, not used
      complex(rk), intent(in) :: val
      complex(rk)             :: res
      !
      real(rk)                :: scale
      !
      call random_number(scale)
      res = val * (1.0_rk + 0.1_rk * (1 - scale)) ! Too much noise will hurt the convergence,
                                                  ! just break the symmetry.
!     res = val * 2.0_rk * (1 - scale)

    end function FLrandom
    !
    !  Return total potential (=molecular + electric field) at a grid point
    !
    function FLtotalPotentialReal(coord) result(res)
      real(rk), intent(in) :: coord(3)
      real(rk)             :: res
      !
      res = real(FLnuclearPotential(coord) + FLelectricField(coord),kind=rk)
    end function FLtotalPotentialReal
    !
    !DEC$ATTRIBUTES FORCEINLINE :: FLtotalPotential
    function FLtotalPotential(coord) result(res)
      real(rk), intent(in) :: coord(3)
      complex(rk)          :: res
      !
      res = FLnuclearPotential(coord) + FLelectricField(coord)
    end function FLtotalPotential
    !
    subroutine FLsetLREmultipoles(npoles,width)
      complex(rk), intent(in) :: npoles(10)  ! Raw multipole moments, in the order:
                                             ! 1      = monopole
                                             ! 2,3,4  = X, Y, Z
                                             ! 5,6,7  = XX, YY, ZZ
                                             ! 8,9,10 = XY, XZ, YZ
                                             ! The moments may be complex (exchange integrals!)
      real(rk), intent(in) :: width          ! Charactestic distribution width
      !
      complex(rk) :: r2 ! Second moment of r
      !
      lre_width = width
      lre_c     = npoles(1)
      lre_d     = npoles(2:4)
      !
      !  Quadrupoles need to be converted to Stone's conventions
      !
      r2 = sum(npoles(5:7))
      lre_q(1,1) = 1.5_rk * npoles(5) - 0.5_rk * r2
      lre_q(2,2) = 1.5_rk * npoles(6) - 0.5_rk * r2
      lre_q(3,3) = 1.5_rk * npoles(7) - 0.5_rk * r2
      lre_q(1,2) = 1.5_rk * npoles(8)
      lre_q(1,3) = 1.5_rk * npoles(9)
      lre_q(2,3) = 1.5_rk * npoles(10)
      lre_q(2,1) = lre_q(1,2)
      lre_q(3,1) = lre_q(1,3)
      lre_q(3,2) = lre_q(2,3)
    end subroutine FLsetLREmultipoles
    !
    !  All formulae are from: A.J. Stone, "The Theory of Intermolecular Forces",
    !  Clarendon Press, Oxford, 1996.
    !
    function FLlrePotential(r) result (v)
      real(rk), intent(in) :: r(3)
      complex(rk)          :: v
      !
      real(rk) :: ra, rc   ! Absolute distance, real and possibly clipped from below
      real(rk) :: t_c      ! Charge interaction
      real(rk) :: t_d(3)   ! Dipolar interaction
      real(rk) :: t_q(3,3) ! Quadrupolar interaction
      !
      ra = sqrt(sum(r**2))
      rc = max(ra,lre_width)
      !
      !  Compute interaction tensors for the potential
      !
      t_c = 1._rk/rc
      t_d = -r/rc**3
      t_q(1,1) = 3._rk*r(1)**2 - ra**2
      t_q(2,2) = 3._rk*r(2)**2 - ra**2
      t_q(3,3) = 3._rk*r(3)**2 - ra**2
      t_q(1,2) = 3._rk*r(1)*r(2)
      t_q(1,3) = 3._rk*r(1)*r(3)
      t_q(2,3) = 3._rk*r(2)*r(3)
      t_q(2,1) = t_q(1,2)
      t_q(3,1) = t_q(1,3)
      t_q(3,2) = t_q(2,3)
      t_q      = t_q/rc**5
      !
      !  Evaluate the interaction
      !
      v = lre_c*t_c - sum(lre_d*t_d) + (1._rk/3._rk)*sum(lre_q*t_q)
    end function FLlrePotential
    !
    subroutine FLasymptoticSetParameters(z,kappa,r0)
      real(rk), intent(in) :: z     ! Effective charge
      real(rk), intent(in) :: kappa ! Exponent parameter
      real(rk), intent(in) :: r0(3) ! Expansion centre
      !
      as_z     = z
      as_kappa = kappa
      as_r0    = r0
    end subroutine FLasymptoticSetParameters
    !
    function FLasymptotic(coord,val) result(res)
      real(rk), intent(in)    :: coord(3)
      complex(rk), intent(in) :: val
      complex(rk)             :: res
      !
      real(rk) :: r, ex
      !
      r   = sqrt(sum((coord-as_r0)**2))
      !
      ex  = as_kappa*r
      if (ex<max_exp) then
        res = r**(1._rk-as_z/as_kappa)*exp(ex)
      else
        res = 0._rk
      end if
    end function FLasymptotic
    !
    subroutine FLharmonicsSetParameters(l,m,r0,rmin,rmax)
      integer(ik), intent(in) :: l, m       ! Order of the harmonic
      real(rk), intent(in)    :: r0(3)      ! Expansion centre
      real(rk), intent(in)    :: rmin, rmax ! Radial range
      !
      harm_l    = l
      harm_m    = m
      harm_r0   = r0
      harm_rmin = rmin
      harm_rmax = rmax
      !
      !  Calculate and remember the overall constant factor for this (l,m)
      !  The prefactor is the Landau&Lifshitz phase convention
      !
      harm_fact = (-1)**((m+abs(m))/2)
      harm_fact = harm_fact * (0._rk,1._rk)**l
      harm_fact = harm_fact * sqrt((2*l+1)/fourpi)
      harm_fact = harm_fact * sqrt(MathFactorial(l-abs(m))/MathFactorial(l+abs(m)))
    end subroutine FLharmonicsSetParameters
    !
    function FLharmonics(coord,val) result(res)
      real(rk), intent(in)    :: coord(3)
      complex(rk), intent(in) :: val
      complex(rk)             :: res
      !
      real(rk)    :: xyz(3), r, ct, xymod
      complex(rk) :: xy
      !
      xyz = coord - harm_r0
      r   = sqrt(sum(xyz**2))
      !
      if (r>=harm_rmin .and. r<=harm_rmax) then
        ct = xyz(3)/r
        if (harm_m>0) then
          xy = cmplx(xyz(1), xyz(2),kind=rk)
        else
          xy = cmplx(xyz(1),-xyz(2),kind=rk)
        end if
        xymod = abs(xy)
        if (xymod>0._rk) then
          xy = xy / xymod
        else
          xy = 1._rk
        end if
        res = harm_fact * MathLegendrePnm(harm_l,abs(harm_m),ct) * xy**abs(harm_m)
      else
        res = 0
      end if
      !
    end function FLharmonics
    !
  end module fields
