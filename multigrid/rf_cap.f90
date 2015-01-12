!
!  Implementation of Moiseyev's RF-CAP for a general boundary geometry,
!  for use in our non-Hermitian static tunneling code.
!
!  See:
!    N. Moiseyev, "Derivations of universal exact complex absorption potentials
!       by the generalized complex coordinate method", J.Phys. B 31, 1431 (1998)
!    Y. Sajeev, M. Sidelka, and N. Moiseyev, "Reflection-free complex absorbing 
!       potential for electronic structure calculation: Feshbach type autoionization
!       resonance of Helium", Chem. Phys. 329, 307 (2006).
!
!  We implement the kinetic part of the potential completely; the only assumption
!  is that the real->complex mapping and its derivatives are available, and form
!  a non-singular Jacobian.
!
!  The potential part is implemented approximately, assuming a multipole expansion
!  valid at long range.
!
!  All integrals are taken numerically, using Becke's molecular integration scheme.
!
!  The implementation was tested by comparing CAP integrals against a naive
!  implementation in Mathematica for a Slater-type 1S function (compared to
!  and STO-6G) and a Gaussian 2Px function. The results agree to >=4 significant
!  digits.
!
  module rf_cap
    use accuracy
    use timer
    use math
    use gamess_internal
    use import_gamess
    use molecular_grid
    implicit none
    private
    public rfc_parameters
    public rfc_matrix_elements
    !
    !  Flavour of the RF-CAP 
    !
    type rfc_parameters
      character(len=20) :: cap_type      = 'moiseyev' ! Can be one of:
                                                      ! 'moiseyev'      - Use Moiseyev's RF-CAP with spherical boundary
                                                      ! 'atom moiseyev' - Use overlapping-atomic spheres boundary
      real(rk)          :: cap_centre(3) = 0._rk      ! Origin of the cap sphere (if applicable) and molecular potential
      real(rk)          :: cap_r0        = 0._rk      ! Characteristic size of the cap sphere
      real(rk)          :: cap_lambda    = 3._rk      ! Steepness of the cap switch-on region
      complex(rk)       :: cap_slope     = 0.1_rk     ! Scaling slope = exp(I*theta), where theta is the scaling angle
      real(rk)          :: cap_mpole(10) = 0._rk      ! Multipole moments of the long-range part of the potential
                                                      ! Multipoles are expected in the Cartesian-moment form
      real(rk)          :: cap_efield(3) = 0._rk      ! Static electric field in the long-range part of the potential
      real(rk)          :: diff_step     = 0.001_rk   ! Finite displacement step (Bohr, Cartesian) used for evaluating 
                                                      ! second derivatives
      !
      !  Multipole moments, using Stone's convention from:
      !  A.J. Stone, "The Theory of Intermolecular Forces", Clarendon Press, Oxford, 1996.
      !  These will be computed from cap_mpole values as needed
      !
      real(rk)          :: lre_c                      ! Monopole
      real(rk)          :: lre_d(3)                   ! Dipole
      real(rk)          :: lre_q(3,3)                 ! Quadrupole
    end type rfc_parameters
    !
    !  Numerical differentiation in the transformed Laplacian may require maintaining
    !  some state between different displaced geometries. This is where we keep it.
    !
    type jacobian_state
      logical           :: first         ! First call, initialize hidden state
      real(rk)          :: rc(3)         ! see scaling_jacobian() below
    end type jacobian_state
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose   = 1        ! Level of output
    real(rk), parameter :: inner_cut = 1e-10_rk ! When the RF-CAP weight drops below (inner_cut), set it to zero!
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    !
    contains
    !
    !  CAPs are evaluated using numerical integration on a grid
    !
    subroutine rfc_matrix_elements(gam,grid,par,fsf,ssf,lsf)
      type(gam_structure), intent(inout)  :: gam      ! Molecule and basis set descriptor
      type(mol_grid), intent(inout)       :: grid     ! Integration grid descriptor
      type(rfc_parameters), intent(inout) :: par      ! Parameters of the CAP
      complex(rk), intent(inout)          :: fsf(:,:) ! spin-free field-interaction Hamiltonian; we need to add CAP terms
      real(rk), intent(inout), optional   :: ssf(:,:) ! spin-free overlap integrals; done as a consistency check
      real(rk), intent(inout), optional   :: lsf(:,:) ! spin-free Laplacian integrals; done as a consistency check
      !
      integer(ik)              :: ib, ipt, npts
      integer(ik)              :: iaol, iaor
      integer(ik)              :: nao
      integer(ik)              :: nbatch
      real(rk), pointer        :: xyzw(:,:)       ! Coordinates and weights of the grid points
      real(rk), allocatable    :: basval(:,:,:,:) ! Values of the basis functions and gradients at the grid points
                                                  ! and displaced positions. Indices are:
                                                  ! [1] = f, df/dx, df/dy, df/dz
                                                  ! [2] = basis function index
                                                  ! [3] = displacement: c, -x, +x, -y, +y, -z, +z
                                                  ! [4] = grid point index
      complex(rk), allocatable :: lapval(:,:,:)   ! Laplacian in Cartesian and complex-scaled coordinates
                                                  ! [1] = cartesian, scaled
                                                  ! [2] = basis function index
                                                  ! [3] = grid point index
      complex(rk), allocatable :: potval(:,:)     ! Scalar part of the RF-CAP
                                                  ! [1] = cartesian, scaled
                                                  ! [2] = grid point index
      logical, allocatable     :: good_point(:)   ! True if point is affected by the scaling transformation
      type(jacobian_state)     :: state           ! Our boundary may be non-smooth; to get reasonable Laplacian,
                                                  ! we need to choose parameters consistently at all displaced points
      complex(rk), allocatable :: fsf_thread(:,:) ! Per-thread integration buffer
      real(rk), allocatable    :: ssf_thread(:,:) ! Per-thread integration buffer
      real(rk), allocatable    :: lsf_thread(:,:) ! Per-thread integration buffer
      !
      call TimerStart('Evaluate RF-CAP')
      !
      !  We assume that the numerical integration grid was prepared for us elsewhere
      !
      call GridPointsBatch(grid,'Batches count',count=nbatch)
      nao = gam%nbasis
      call convert_multipoles(par)
      !$omp parallel default(none) &
      !$omp& shared(nbatch,nao,grid,gam,fsf,ssf,lsf,par) &
      !$omp& private(ib,ipt,npts,iaol,iaor,xyzw,basval,lapval,potval,good_point) &
      !$omp& private(state,fsf_thread,ssf_thread,lsf_thread)
      nullify(xyzw)
      allocate (fsf_thread(nao,nao))
      fsf_thread = 0
      if (present(ssf)) then
        allocate (ssf_thread(nao,nao))
        ssf_thread = 0
      end if
      if (present(lsf)) then
        allocate (lsf_thread(nao,nao))
        lsf_thread = 0
      end if
      !$omp do schedule(dynamic)
      grid_batches: do ib=1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
        npts = size(xyzw,dim=2)
        if (allocated(basval)) then
          if (size(basval,dim=4)/=npts) deallocate (basval,lapval,potval,good_point)
        end if
        if (.not.allocated(basval)) allocate(basval(4,nao,7,npts),lapval(2,nao,npts),potval(2,npts),good_point(npts))
        !
        evaluate_basis_functions: do ipt=1,npts
          state%first = .true.
          call fill_potential(gam,par,state,xyzw(1:3,ipt),potval(:,ipt),good_point(ipt))
          call fill_functions(gam,par,state,xyzw(1:3,ipt),basval(:,:,:,ipt),lapval(:,:,ipt))
        end do evaluate_basis_functions
        !
        !  Accumulate the integrals
        !
        !  Accuracy check - overlap
        !
        if (present(ssf)) then
          overlap_right_aos: do iaor=1,nao
            overlap_left_aos: do iaol=1,nao
              ssf_thread(iaol,iaor) = ssf_thread(iaol,iaor) + sum(basval(1,iaol,1,:)*basval(1,iaor,1,:)*xyzw(4,:))
            end do overlap_left_aos
          end do overlap_right_aos
        end if
        !
        assemble_points: do ipt=1,npts
          !
          !  Accuracy check - Laplacian; integrate in the whole space
          !
          if (present(lsf)) then
            lap_right_aos: do iaor=1,nao
              lap_left_aos: do iaol=1,nao
                lsf_thread(iaol,iaor) = lsf_thread(iaol,iaor) &
                    + xyzw(4,ipt)*basval(1,iaol,1,ipt)*real(lapval(1,iaor,ipt),kind=rk)
              end do lap_left_aos
            end do lap_right_aos
          end if
          !
          !  CAP will only be integrated in the region where the laplacian is affected by complex scaling
          !
          if (.not.good_point(ipt)) cycle assemble_points
          cap_right_aos: do iaor=1,nao
            cap_left_aos: do iaol=1,nao
              fsf_thread(iaol,iaor) = fsf_thread(iaol,iaor) &
                  + xyzw(4,ipt)*basval(1,iaol,1,ipt)*basval(1,iaor,1,ipt)*(potval(2,ipt)-potval(1,ipt)) &
                  -0.5_rk*xyzw(4,ipt)*basval(1,iaol,1,ipt)*(lapval(2,iaor,ipt)-lapval(1,iaor,ipt))
            end do cap_left_aos
          end do cap_right_aos
        end do assemble_points
        !
      end do grid_batches
      !$omp end do nowait
      if (associated(xyzw)  ) deallocate (xyzw)
      if (allocated (basval)) deallocate (basval,lapval,potval,good_point)
      !$omp critical
      fsf = fsf + fsf_thread
      if (present(ssf)) ssf = ssf + ssf_thread
      if (present(lsf)) lsf = lsf + lsf_thread
      !$omp end critical
      deallocate (fsf_thread)
      if (present(ssf)) deallocate (ssf_thread)
      if (present(lsf)) deallocate (lsf_thread)
      !$omp end parallel
      call TimerStop('Evaluate RF-CAP')
    end subroutine rfc_matrix_elements
    ! 
    !  Convert multipole moments of charge distribution to Stone's
    !  (traceless) conventions
    !
    subroutine convert_multipoles(par)
      type(rfc_parameters), intent(inout) :: par        ! Parameters of the CAP
      !
      real(rk) :: r2 ! Second moment of r
      !
      par%lre_c      = par%cap_mpole(1)
      par%lre_d      = par%cap_mpole(2:4)
      !
      !  Quadrupoles need to be converted to Stone's conventions
      !
      r2 = sum(par%cap_mpole(5:7))
      par%lre_q(1,1) = 1.5_rk * par%cap_mpole(5) - 0.5_rk * r2
      par%lre_q(2,2) = 1.5_rk * par%cap_mpole(6) - 0.5_rk * r2
      par%lre_q(3,3) = 1.5_rk * par%cap_mpole(7) - 0.5_rk * r2
      par%lre_q(1,2) = 1.5_rk * par%cap_mpole(8)
      par%lre_q(1,3) = 1.5_rk * par%cap_mpole(9)
      par%lre_q(2,3) = 1.5_rk * par%cap_mpole(10)
      par%lre_q(2,1) = par%lre_q(1,2)
      par%lre_q(3,1) = par%lre_q(1,3)
      par%lre_q(3,2) = par%lre_q(2,3)
    end subroutine convert_multipoles
    !
    subroutine fill_potential(gam,par,state,xyz,potval,good_point)
      type(gam_structure), intent(inout)  :: gam        ! Molecule and basis set descriptor
      type(rfc_parameters), intent(in)    :: par        ! Parameters of the CAP
      type(jacobian_state), intent(inout) :: state      ! Choice of the reference point, if needed
      real(rk), intent(in)                :: xyz(:)     ! Grid position
      complex(rk), intent(out)            :: potval(:)  ! Potential in Cartesian (1) and complex-scaled (2)
                                                        ! coordinates at point xyz
      logical, intent(out)                :: good_point ! True if scaling transformation affects this point
      !
      complex(rk) :: c_xyz(3)   ! Complex-scaled coordinates
      !
      !  Even if we do not calculate the potential term here, we still must
      !  make a call scaling_map - this fixes the choice of the reference
      !  point, which will be needed in numerical differentiation later on
      !
      call scaling_map(gam,par,state,xyz,good_point,c_xyz)
      !
      !  Multipole potentials; All formulae are from: A.J. Stone, 
      !  "The Theory of Intermolecular Forces", Clarendon Press, Oxford, 1996.
      !
      potval(1) = evaluate_potential(cmplx(xyz,kind=kind(xyz)))
      potval(2) = evaluate_potential(c_xyz)
      !
      contains
      function evaluate_potential(r) result(v)
        complex(rk), intent(in) :: r(3)
        complex(rk)             :: v
        !
        complex(rk) :: ra       ! Absolute distance, possibly clipped from below
        complex(rk) :: t_c      ! Charge interaction
        complex(rk) :: t_d(3)   ! Dipolar interaction
        complex(rk) :: t_q(3,3) ! Quadrupolar interaction
        !
        ra = sqrt(sum((r-par%cap_centre)**2))
        if (abs(ra)<=1e-2_rk) ra = 1e-3_rk * (ra/abs(ra))
        !
        !  Compute interaction tensors for the potential
        !
        t_c = 1._rk/ra
        t_d = -r/ra**3
        t_q(1,1) = 3._rk*r(1)**2 - ra**2
        t_q(2,2) = 3._rk*r(2)**2 - ra**2
        t_q(3,3) = 3._rk*r(3)**2 - ra**2
        t_q(1,2) = 3._rk*r(1)*r(2)
        t_q(1,3) = 3._rk*r(1)*r(3)
        t_q(2,3) = 3._rk*r(2)*r(3)
        t_q(2,1) = t_q(1,2)
        t_q(3,1) = t_q(1,3)
        t_q(3,2) = t_q(2,3)
        t_q      = t_q/ra**5
        !
        !  Evaluate the interaction
        !
        v = par%lre_c*t_c - sum(par%lre_d*t_d) + (1._rk/3._rk)*sum(par%lre_q*t_q) - sum(par%cap_efield*r)
      end function evaluate_potential
    end subroutine fill_potential
    !
    subroutine fill_functions(gam,par,state,xyz,basval,lapval)
      type(gam_structure), intent(inout)  :: gam           ! Molecule and basis set descriptor
      type(rfc_parameters), intent(in)    :: par           ! Parameters of the CAP
      type(jacobian_state), intent(inout) :: state         ! Choice of the reference point, if needed
      real(rk), intent(in)                :: xyz(:)        ! Grid position
      real(rk), intent(out)               :: basval(:,:,:) ! Function values and gradients at this point
      complex(rk), intent(out)            :: lapval(:,:)   ! Laplacian in Cartesian (1,:) and complex-scaled (2,:)
                                                           ! coordinates at each point
      !
      real(rk)     :: tmp(3), scl
      integer(ik)  :: nbas, ibas, id
      complex(rk)  :: kmat(3,3,7)  ! Inverse of Jacobian matrices for the complex scaling transformation
      complex(rk)  :: vol (7)      ! Determinant of Jacobian matrices (volume element)
      complex(rk)  :: volg(3,7)    ! Gradient of the volume element
      complex(rk)  :: kmatg (3,7)  ! Transformed basis function gradients at displaced points
      complex(rk)  :: gvolpsi(3)   ! Scratch space for d/dr |jacobian|**-0.5 psi
      complex(rk)  :: outer(3,3)   ! 
      !
      !  Offset table for numerical differentiation
      !
      integer(ik), parameter :: diff_index(7) = (/ 0,      1,     1,      2,     2,      3,     3     /)
      real(rk), parameter    :: diff_scale(7) = (/ 0._rk, -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, 1._rk /)
      !
      nbas = size(basval,dim=2)
      if (size(lapval,dim=1)/=2 .or. size(lapval,dim=2)/=nbas .or. &
          size(basval,dim=1)/=4 .or. size(basval,dim=3)/=7) stop 'rf_cap%fill_functions - bad array dimensions'
      !
      numerical_differentiate: do id=1,7
        tmp = xyz
        if (diff_index(id)/=0) then
          tmp(diff_index(id)) = tmp(diff_index(id)) + diff_scale(id)*par%diff_step
        end if
        call gamess_evaluate_functions(tmp,basval(:,:,id),gam)
        call scaling_jacobian(gam,par,state,tmp,kmat(:,:,id),vol(id))
        call volume_gradient(tmp,volg(:,id))
      end do numerical_differentiate
      !
      !  Now, calculate the Laplacian in the Cartesian and complex-scaled coordinates
      !
      scl = 1/(2*par%diff_step) ! Prefactor due to numerical differentiation
      evaluate_laplacians: do ibas=1,nbas
        lapval(1,ibas) = scl * ( &
                         (basval(2,ibas,3)-basval(2,ibas,2))   & ! (2,:,3) = X gradient, +X displacement; (2,:,2) = X gradient, -X displacement
                       + (basval(3,ibas,5)-basval(3,ibas,4))   & ! (3,:,3) = Y gradient, +Y displacement; (3,:,4) = X gradient, -Y displacement
                       + (basval(4,ibas,7)-basval(4,ibas,6)) )   ! (4,:,3) = Z gradient, +Z displacement; (4,:,6) = X gradient, -Z displacement
        !
        !  Forming the transformed Laplacian is a bit ugly in the general case.
        !  The expression is:
        !
        !    |jacobian|**0.5 Tr( jacobian**-1 [d/dr (outer x) (jacobian**-1 d/dr (|jacobian|**-0.5 psi))] )
        !
        !  We'll assemble this quantity step-by-step, using chain rules and numerical
        !  differentiation where appropriate.
        !
        transform_gradients: do id=2,7
          !
          !  (d/dr) |j|**-0.5 psi = |j|**-0.5 (d/dr psi) - 0.5 |j|**-1.5 psi (d/dr |j|)
          !
          gvolpsi(:) = vol(id)**(-0.5_rk) * basval(2:4,ibas,id) &
                     - 0.5_rk*vol(id)**(-1.5_rk) * basval(1,  ibas,id) * volg(:,id)
          !
          !  j**(-1) [(d/dr) |j|**-0.5 psi]
          !
          kmatg(:,id) = matmul(kmat(:,:,id),gvolpsi)
        end do transform_gradients
        !
        !  d/dr (outer x) [j**(-1) [(d/dr) |j|**-0.5 psi]]
        !  using numerical differentiation to calculate the derivative from displaced values
        !
        outer(1,:) = scl * (kmatg(:,3)-kmatg(:,2)) ! d/dX
        outer(2,:) = scl * (kmatg(:,5)-kmatg(:,4)) ! d/dY
        outer(3,:) = scl * (kmatg(:,7)-kmatg(:,6)) ! d/dZ
        !
        !  Finally, the Laplacian; sum(kmat * outer) = Tr(matmul(kmat,outer))
        !
        lapval(2,ibas) = (vol(1)**0.5_rk) * sum(kmat(:,:,1)*outer)
      end do evaluate_laplacians
      !
      contains
        ! 
        !  This subroutine will recalculate several of displaced points;
        !  Because this part of the code is not expected to be on the
        !  critical path, I'd rather have simpler, marginally slower code.
        !
        subroutine volume_gradient(xyz2,volg)
          real(rk), intent(in)     :: xyz2(:) ! Coordinate where volume gradient is needed
          complex(rk), intent(out) :: volg(:) ! Gradient
          !
          real(rk)    :: tmp2(3)  ! Displacements from the displaced geometries
          integer(ik) :: id2      ! Displacement index for numerical differentiation
          complex(rk) :: dum(3,3) ! We do not care about the inverse jacobian here
          complex(rk) :: vol2(7)  ! The first element is not used
          !
          numdiff2: do id2=2,7 ! Central point is not needed, so no check for index==0
            tmp2 = xyz2
            tmp2(diff_index(id2)) = tmp2(diff_index(id2)) + diff_scale(id2)*par%diff_step
            call scaling_jacobian(gam,par,state,tmp2,dum,vol2(id2))
          end do numdiff2
          volg(1) = (vol2(3)-vol2(2))/(2*par%diff_step)
          volg(2) = (vol2(5)-vol2(4))/(2*par%diff_step)
          volg(3) = (vol2(7)-vol2(6))/(2*par%diff_step)
        end subroutine volume_gradient
    end subroutine fill_functions
    !
    subroutine radial_map(par,x,good,f,dfdx)
      type(rfc_parameters), intent(in)  :: par   ! Parameters of the CAP
      real(rk), intent(in)              :: x     ! Argment of the Moiseyev's switching function. 
                                                 ! x is assumed to be positive
      logical, intent(out)              :: good  ! .true. if this point is in the region affected by the cap
                                                 ! f and dfdx are evaluated regardless of the value of (good)
      complex(rk), intent(out)          :: f     ! f(x)
      complex(rk), intent(out)          :: dfdx  ! d f(x) / d x
      !
      real(rk) :: xplus, xminus, selector
      !
      xplus    = par%cap_lambda * (x+par%cap_r0)
      xminus   = par%cap_lambda * (x-par%cap_r0)
      !
      selector = 1._rk + 0.5_rk*(tanh(xminus) - tanh(xplus))
      if (selector<=inner_cut) then
        !
        !  We are in the inner asymptotic region, F(x)=x
        !
        good = .false.
        f    = x
        dfdx = 1
      else if (selector>=2.0_rk-inner_cut) then
        !
        !  We are in the outer asymptotic region, F(x)=r0 + slope*(x-r0)
        !
        good = .true.
        f    = par%cap_r0 + par%cap_slope * (x-par%cap_r0)
        dfdx = par%cap_slope
      else
        !
        !  We are in the transition region, do it the hard way.
        !  The expression for f(x) is not numerically stable as written, so we have to
        !  do be a little careful
        !
        good = .true.
        f    = x + (par%cap_slope-1)*(x + (logcosh(xminus)-logcosh(xplus))/(2*par%cap_lambda))
        dfdx = 1 + (par%cap_slope-1)*selector
      end if
      contains
      real(rk) function logcosh(x) 
        real(rk), intent(in) :: x ! Argument for which we need log(cosh(x))
        !
        if (abs(x)>=max_exp) then
          logcosh = abs(x) - log(2._rk)
        else
          logcosh = log(cosh(x))
        end if
      end function logcosh
    end subroutine radial_map
    !
    subroutine choose_reference_point(gam,par,xyz,state)
      type(gam_structure), intent(in)     :: gam     ! Molecule and basis set descriptor
      type(rfc_parameters), intent(in)    :: par     ! Parameters of the CAP
      real(rk), intent(in)                :: xyz(:)  ! Grid position
      type(jacobian_state), intent(inout) :: state   ! Internal state
      !
      integer(ik) :: iat
      real(rk)    :: at_xyz(3), rat, small_rat
      !
      if (.not.state%first) return
      state%first = .false.
      !
      select case (par%cap_type)
        case default
          write (out,"('rf_cap%choose_reference_point: CAP type ',a,' is not implemented')") trim(par%cap_type)
          stop 'rf_cap%choose_reference_point - bad CAP'
        case ('moiseyev')
          state%rc = par%cap_centre
        case ('atom moiseyev')
          small_rat = safe_max
          scan_atoms: do iat=1,gam%natoms
            at_xyz = real(gam%atoms(iat)%xyz,kind=rk)/abohr
            rat    = sqrt(sum((xyz-at_xyz)**2))
            if (rat<small_rat) then
              state%rc  = at_xyz
              small_rat = rat
            end if
          end do scan_atoms
          if (small_rat>=safe_max) then
            stop 'rf_cap%choose_reference_point - all atoms are impossibly far'
          end if
      end select
    end subroutine choose_reference_point
    !
    !  Our coordinate mapping is given by:
    !    F(x,y,z) = [x,y,z] (1/r) radial_map(r) + [Rx,Ry,Rz]
    !    r        = sqrt(x**2+y**2+z**2))
    !  where radial_map is Moiseyev's radial switching function
    !  and all distances are relative to the reference point
    !  [Rx,Ry,Rz]
    !
    subroutine scaling_map(gam,par,state,xyz,good,c_xyz)
      type(gam_structure), intent(in)     :: gam      ! Molecule and basis set descriptor
      type(rfc_parameters), intent(in)    :: par      ! Parameters of the CAP
      type(jacobian_state), intent(inout) :: state    ! Choice of the reference point, if needed
      real(rk), intent(in)                :: xyz(:)   ! Grid position
      logical, intent(out)                :: good     ! Grid position is in the region affected by scaling
      complex(rk), intent(out)            :: c_xyz(:) ! Complex-scaled grid position
      !
      real(rk)    :: xyz_rel(3), r_rel  ! Coordinates and distance relative to the reference
      complex(rk) :: map, dmap          ! Radial map and its derivative
      !
      call choose_reference_point(gam,par,xyz,state)
      !
      xyz_rel = xyz - state%rc
      r_rel   = sqrt(sum(xyz_rel**2))
      call radial_map(par,r_rel,good,map,dmap)
      c_xyz   = state%rc + map*xyz_rel/r_rel
    end subroutine scaling_map
    !
    subroutine scaling_jacobian(gam,par,state,xyz,kmat,vol)
      type(gam_structure), intent(in)     :: gam        ! Molecule and basis set descriptor
      type(rfc_parameters), intent(in)    :: par        ! Parameters of the CAP
      type(jacobian_state), intent(inout) :: state      ! Choice of the reference point, if needed
      real(rk), intent(in)                :: xyz(:)     ! Grid position
      complex(rk), intent(out)            :: kmat(:,:)  ! Inverse of the transformation's Jacobian
      complex(rk), intent(out)            :: vol        ! Determinant of the Jacobian (volume element)
      !
      complex(rk) :: jmat(3,3)          ! Jacobian matrix
      real(rk)    :: xyz_rel(3), r_rel  ! Coordinates and distance relative to the reference
      logical     :: good               ! Unused
      complex(rk) :: map, dmap          ! Radial map and its derivative
      integer(ik) :: mc, gc             ! Map and gradient components
      !
      call choose_reference_point(gam,par,xyz,state)
      !
      xyz_rel = xyz - state%rc
      r_rel   = sqrt(sum(xyz_rel**2))
      call radial_map(par,r_rel,good,map,dmap)
      !
      !  Evaluate the Jacobian using chain rule for differentiation:
      !   J(mc,gc) = d F_mc / d r_gc
      !  Note that for our transformation, J is symmetric, but not
      !  necessarily Hermitian.
      !
      jacobian_grad: do gc=1,3
        jacobian_map: do mc=1,3
          jmat(mc,gc) = xyz_rel(mc)*xyz_rel(gc)*(dmap-map/r_rel)/r_rel**2
        end do jacobian_map
        jmat(gc,gc) = jmat(gc,gc) + map/r_rel
      end do jacobian_grad
      !
      !  We are nearly done - calculate volume element and the inverse
      !
      vol  = MathDet3x3(jmat)
      if (vol==0._rk) stop 'rf_cap%scaling_jacobian - Jacobian is singular'
      kmat = MathInv3x3(jmat)
      !*ps debug
      ! write (out,"('Position: ',3f14.7,' r = ',f14.7,' volume element = ',f20.12,1x,f20.12)") xyz_rel, r_rel, vol
      ! write (out,"('Jacobian: ',3(1x,f20.12,1x,f20.12,2x))") jmat(1,:)
      ! write (out,"('          ',3(1x,f20.12,1x,f20.12,2x))") jmat(2,:)
      ! write (out,"('          ',3(1x,f20.12,1x,f20.12,2x))") jmat(3,:)
      ! write (out,"('inv. jac: ',3(1x,f20.12,1x,f20.12,2x))") kmat(1,:)
      ! write (out,"('          ',3(1x,f20.12,1x,f20.12,2x))") kmat(2,:)
      ! write (out,"('          ',3(1x,f20.12,1x,f20.12,2x))") kmat(3,:)
      !*ps debug
    end subroutine scaling_jacobian
  end module rf_cap
