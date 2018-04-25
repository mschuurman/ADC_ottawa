!######################################################################
!
!  Evaluation of 1-electron matrix elements of complex absorbing
!  potentials using numerical integration.
!
!  This module must not be called from a parallel region.
!  The interface design for the cap_ parameters is a bit messy, but
!  should do.
!
!  Taken from the ECS-MP2 code of Serguei Patchkovskii
!######################################################################

module basis_cap
  use accuracy
  use timer
  use math
  use rf_cap
  use gamess_internal
  use import_gamess
  use molecular_grid
  implicit none
  !
  private
  public cap_evaluate, cap_dump_effective
  public cap_type, cap_centre, cap_r0, cap_strength, cap_order, cap_lambda
  public cap_theta, cap_mpole, cap_efield, cap_diff_step
  public cap_kmin, cap_delta, cap_limit
  public cap_aocoeff
  !
  !  Adjustable parameters
  !
  character(len=clen), save :: cap_type   = 'monomial' ! Can be either of:
                                                       ! 'monomial'      = -I*eta*(r-r0)**n for r>r0
                                                       !                   Reading Riss and Meyer, JCP 105, 1409 (1996) is
                                                       !                   highly recommended before using this one
                                                       ! 'atom monomial' = Same as 'monomial', but r0 is chosen as the position
                                                       !                   of the nearest nucleus
                                                       ! 'moiseyev'      = Moiseyev's non-local perfect CAP, single sphere
                                                       ! 'atom moiseyev' = Moiseyev's non-local perfect CAP, atom-centered spheres
                                                       ! 'manolopoulous' = JWKB perfect absorber; DOT NOT USE
                                                       ! 'sigmoidal'     = Cavity-type sigmoidal CAP
                                                       ! The cap_centre and cap_r0 paraments apply to all CAPs
  real(rk), save      :: cap_centre(3)   = 0._rk       ! Centre of the spherical shell carrying the complex absorbing potential
  real(rk), save      :: cap_r0          = 12._rk      ! Radius at which CAP starts
                                                       ! === Parameters related to the monomial CAP
  real(rk), save      :: cap_strength    = 0.02_rk     ! eta parameter of Riss&Meyer
  integer(ik), save   :: cap_order       = 2_ik        ! Power of the cap 
                                                       ! The default cap_strength (0.02) and cap_order (2) correspond to the 1e-4 
                                                       ! R&M "error" in absorbing 1 Ry incoming wave over the effective
                                                       ! length of 10 Bohr.
                                                       ! === Parameters related to the Moiseyev's CAP
  real(rk), save      :: cap_lambda      = 3.0_rk      ! Steepness of the switch-off region; width of the region is
                                                       ! on the order of 1/cap_lambda
  real(rk), save      :: cap_theta       = 0.15_rk     ! Complex scaling angle, Radian
  real(rk), save      :: cap_mpole(10)   = 0._rk       ! Long-range multiplicative potential to use in the CAP; expressed as 
                                                       ! the moments of a charge distribution around the origin
  real(rk), save      :: cap_efield(3)   = 0._rk       ! Static electric field to use in the CAP
  real(rk), save      :: cap_diff_step   = 0.001_rk    ! Numerical differentiation step (Bohr) to use in constructing the Laplacian
                                                       ! === Parameters related to the Manolopolous CAP
  real(rk), save      :: cap_kmin        = 0.30_rk     ! Smallest momentum which CAP absorbs perfectly
  real(rk), save      :: cap_delta       = 0.2_rk      ! JWKB scaling constant; determines absorption efficiency for cap_kmin
                                                       ! cap_delta = 0.2 corresponds to 1% efficiency (?)
  real(rk), save      :: cap_limit       = 10._rk      ! Manolopoulos CAPs are singular; this plays havoc with Gussian basis functions,
                                                       ! so we'll cap the CAP (;o) at large, but still sensible value
                                                       ! === end of CAP-related parameters
  !
  !  Parameters related to the absorbing boundary
  !
  real(rk), parameter         :: ma_c = 2.622057554292119810464840_rk ! The position of the CAP singularity
  real(rk), parameter         :: ma_a = 1._rk - 16._rk*ma_c**(-3)
  real(rk), parameter         :: ma_b = (1._rk - 17._rk*ma_c**(-3)) * ma_c**(-2)
  real(rk)                    :: cap_width         ! Overall width of the CAP
  real(rk)                    :: cap_scale         ! Scaling parameter
  !

contains

!######################################################################

  function manolopoulosVcap(xyz) result(v)
    real(rk), intent(in) :: xyz(:) ! Coordinates of the point
    real(rk)             :: v      ! The potential
    !
    real(rk)             :: r0     ! Distance to the CAP origin
    real(rk)             :: x      ! Dimensionless penetration coordinate
    real(rk)             :: xt
    !
    r0 = sqrt(sum((xyz-cap_centre)**2))
    x  = (r0 - cap_r0)/cap_width
    !
    if (x<=0._rk) then
       v = 0._rk
    else
       xt = ma_c * min(max(0._rk,x),1._rk-spacing(100._rk))
       v  = ma_a * xt - ma_b * xt**3 + 4._rk * (ma_c-xt)**(-2) - 4._rk * (ma_c+xt)**(-2)
       v  = -min(cap_limit,cap_scale * v)
       !
       !  Before Dec 17, 2013, there was a sign error in our Manolopoulous CAP.
       !  Forttunately, this seems to just flip the sign of the imaginary part
       !  of the energy (and presumably flips the left and right eigenvectors)
       !
    end if
  end function manolopoulosVcap

!######################################################################

  function monomialVcap(xyz) result(v)
    real(rk), intent(in) :: xyz(:) ! Coordinates of the point
    real(rk)             :: v      ! The potential
    !
    real(rk)             :: rc     ! Penetration depth
    !
    rc = sqrt(sum((xyz-cap_centre)**2)) - cap_r0 
    !
    v = 0._rk
    if (rc>0._rk) then
       v = -cap_strength * rc**cap_order
    end if
  end function monomialVcap

!######################################################################
  
  function atomMonomialVcap(gam,xyz) result(v)
    type(gam_structure), intent(in) :: gam    ! Basis set and structure descriptor
    real(rk), intent(in)            :: xyz(:) ! Coordinates of the point
    real(rk)                        :: v      ! The potential
    !
    integer(ik)          :: iat
    real(rk)             :: at_xyz(3), rat
    real(rk)             :: rc     ! Penetration depth
    !
    rc = safe_max
    scan_atoms: do iat=1,gam%natoms
       at_xyz = real(gam%atoms(iat)%xyz,kind=rk)/abohr
       rat    = sqrt(sum((xyz-at_xyz)**2))
       rc     = min(rc,rat)
    end do scan_atoms
    !
    rc = rc - cap_r0 
    !
    v = 0._rk
    if (rc>0._rk) then
       v = -cap_strength * rc**cap_order
    end if
  end function atomMonomialVcap

!######################################################################

  function sigmoidalVcap(gam,xyz) result(v)
    type(gam_structure), intent(in) :: gam    ! Basis set and structure descriptor
    real(rk), intent(in)            :: xyz(:) ! Coordinates of the point
    real(rk)                        :: v      ! The potential

    integer(ik) :: i
    real(rk)    :: r,width
    real(rk)    :: val(gam%natoms)
    real(rk)    :: at_xyz(3)

    ! Temporary hardwiring of the cap width parameter
    width=10.0d0
    
    ! Loop over atoms
    do i=1,gam%natoms

       ! Skip dummy atoms
       if (gam%atoms(i)%name.eq.'x'&
            .or.gam%atoms(i)%name.eq.'X') cycle

       ! Cartesian coordinates of the current atom (in Bohr)
       at_xyz=real(gam%atoms(i)%xyz,kind=rk)/abohr

       ! Distance from the point xyz to the current atom
       r=sqrt(dot_product(xyz-at_xyz,xyz-at_xyz))
       
       ! Calculate the CAP value
       if (r.le.cap_r0) then
          val(i)=0.0d0
       else if (r.lt.cap_r0+width) then
          val(i)=-cap_strength*(sin(0.5d0*pi*(r-cap_r0)/width))**2
       else
          val(i)=-cap_strength
       endif

    enddo

    v=minval(val)
    
  end function sigmoidalVcap
    
!######################################################################
!
!  Multiplicative CAPs can all be treated alike
!
!######################################################################

  subroutine evaluate_caps_multiplicative(gam,grid,fsf,ssf)
    type(gam_structure), intent(inout) :: gam      ! Basis set and structure descriptor
    type(mol_grid), intent(inout)      :: grid
    complex(rk), intent(inout)         :: fsf(:,:) ! spin-free field-interaction Hamiltonian; we need to add CAP terms
    real(rk), intent(inout)            :: ssf(:,:) ! spin-free overlap integrals; done as a consistency check
    !
    integer(ik)           :: ib, ipt, npts
    integer(ik)           :: iaol, iaor
    integer(ik)           :: nbatch, nao
    real(rk), pointer     :: xyzw(:,:)       ! Coordinates and weights of the grid points
    real(rk), allocatable :: basval(:,:,:)   ! Values of the basis functions at the grid points
    real(rk), allocatable :: cap(:)          ! Values of the imaginary part of the potential at grid points
    real(rk), allocatable :: fsf_thread(:,:) ! Per-thread integration buffer
    real(rk), allocatable :: ssf_thread(:,:) ! Per-thread integration buffer
    !
    call GridPointsBatch(grid,'Batches count',count=nbatch)
    nao = gam%nbasis
    !$omp parallel default(none) &
    !$omp& shared(nbatch,nao,grid,gam,fsf,ssf,cap_type) &
    !$omp& private(ib,ipt,npts,iaol,iaor,xyzw,basval,cap,fsf_thread,ssf_thread)
    nullify(xyzw)
    allocate (fsf_thread(nao,nao),ssf_thread(nao,nao))
    fsf_thread = 0
    ssf_thread = 0
    !$omp do schedule(dynamic)
    grid_batches: do ib=1,nbatch
       !
       !  Get grid points
       !
       call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
       npts = size(xyzw,dim=2)
       if (allocated(basval)) then
          if (size(basval,dim=3)/=npts) deallocate (basval,cap)
       end if
       if (.not.allocated(basval)) allocate(basval(1,nao,npts),cap(npts))
       evaluate_basis_functions: do ipt=1,npts
          call gamess_evaluate_functions(xyzw(1:3,ipt),basval(:,:,ipt),gam)
       end do evaluate_basis_functions
       ! Premultiply integrand by the point weight; saves us an operation later
       evaluate_cap: do ipt=1,npts
          select case (cap_type)
          case default
             stop 'basis_cap%evaluate_caps_multiplicative - unrecognized cap_type'
          case ('manolopoulos')
             cap(ipt) = manolopoulosVcap(xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('monomial')
             cap(ipt) = monomialVcap(xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('atom monomial')
             cap(ipt) = atomMonomialVcap(gam,xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('sigmoidal')
             cap(ipt) = sigmoidalVcap(gam,xyzw(1:3,ipt)) * xyzw(4,ipt)
          end select
       end do evaluate_cap
       !
       !  Accumulate the integrals
       !
       right_aos: do iaor=1,nao
          left_aos: do iaol=1,nao
             fsf_thread(iaol,iaor) = fsf_thread(iaol,iaor) + sum(basval(1,iaol,:)*basval(1,iaor,:)*cap)
             ssf_thread(iaol,iaor) = ssf_thread(iaol,iaor) + sum(basval(1,iaol,:)*basval(1,iaor,:)*xyzw(4,:))
          end do left_aos
       end do right_aos
       !
    end do grid_batches
    !$omp end do nowait
    if (associated(xyzw)  ) deallocate (xyzw)
    if (allocated (basval)) deallocate (basval,cap)
    !$omp critical
    fsf = fsf + (0,1)*fsf_thread
    ssf = ssf + ssf_thread
    !$omp end critical
    deallocate (fsf_thread,ssf_thread)
    !$omp end parallel
  end subroutine evaluate_caps_multiplicative

!######################################################################
!
!  Calculate effective multiplicative CAP, and (if possible) compare
!  it to the reference CAP
!
!######################################################################

  subroutine dump_vcap(iu,gam,grid,vh)
    integer(ik)                        :: iu       ! Unit for the output; negative unit suppreses the dump
    type(gam_structure), intent(inout) :: gam      ! Basis set and structure descriptor
    type(mol_grid), intent(inout)      :: grid
    complex(rk), intent(in)            :: vh(:,:)  ! Matrix elements of the CAP
    !
    integer(ik)              :: ib, ipt, npts
    integer(ik)              :: iaol, iaor
    integer(ik)              :: nbatch, nao
    real(rk), pointer        :: xyzw(:,:)       ! Coordinates and weights of the grid points
    real(rk), allocatable    :: basval(:,:,:)   ! Values of the basis functions at the grid points
    real(rk), allocatable    :: cap(:)          ! Values of the imaginary part of the potential at grid points
    complex(rk), allocatable :: hcap(:)         ! Values of the effective multiplicative CAP
    !
    call GridPointsBatch(grid,'Batches count',count=nbatch)
    nao = gam%nbasis
    if (size(vh,dim=1)/=nao .or. size(vh,dim=2)/=nao) then
       stop 'basis_cap%dump_vcap - bad VH dimensions'
    end if
    if (iu>=0) then
       write (iu,"('#',7(1x,a23,1x))") ' X ', ' Y ', ' Z ', ' W ', ' Re(VH) ', ' Im(VH) ', ' Im(VC) '
    end if
    !
    !  This routine will mostl likely be I/O bound, and serialize on the $critical
    !  section; however, there is no obvious harm in keeping it parallel.
    !
    !$omp parallel default(none) &
    !$omp& shared(nbatch,nao,grid,gam,cap_type,iu,vh) &
    !$omp& private(ib,ipt,npts,iaol,iaor,xyzw,basval,cap,hcap)
    nullify(xyzw)
    !$omp do schedule(dynamic)
    grid_batches: do ib=1,nbatch
       !
       !  Get grid points
       !
       call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
       npts = size(xyzw,dim=2)
       if (allocated(basval)) then
          if (size(basval,dim=3)/=npts) deallocate (basval,cap,hcap)
       end if
       if (.not.allocated(basval)) allocate(basval(1,nao,npts),cap(npts),hcap(npts))
       evaluate_basis_functions: do ipt=1,npts
          call gamess_evaluate_functions(xyzw(1:3,ipt),basval(:,:,ipt),gam)
       end do evaluate_basis_functions
       ! 
       !  Evaluate CAP from the real-space definition
       !
       evaluate_cap: do ipt=1,npts
          select case (cap_type)
          case default
             cap(ipt) = 0
          case ('manolopoulos')
             cap(ipt) = manolopoulosVcap(xyzw(1:3,ipt))
          case ('monomial')
             cap(ipt) = monomialVcap(xyzw(1:3,ipt))
          case ('atom monomial')
             cap(ipt) = atomMonomialVcap(gam,xyzw(1:3,ipt))
          case ('sigmoidal')
             cap(ipt) = sigmoidalVcap(gam,xyzw(1:3,ipt))
          end select
       end do evaluate_cap
       !
       !  Evaluate the same thing from the Hilbert-space projection VH
       !
       evaluate_hcap: do ipt=1,npts
          hcap(ipt) = sum(basval(1,:,ipt)*matmul(vh,basval(1,:,ipt)))
       end do evaluate_hcap
       !
       !  Dump the result
       !
       if (iu>=0) then
          !$omp critical
          print_cap: do ipt=1,npts
             write (iu,"(7(1x,g24.16))") xyzw(:,ipt), hcap(ipt), cap(ipt)
          end do print_cap
          !$omp end critical
       end if
       !
    end do grid_batches
    !$omp end do nowait
    if (associated(xyzw)  ) deallocate (xyzw)
    if (allocated (basval)) deallocate (basval,cap)
    !$omp end parallel
  end subroutine dump_vcap

!######################################################################
  
  subroutine prepape_for_cap_integration(grid,gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang)
    type(mol_grid), intent(inout)      :: grid            ! Grid descriptor
    type(gam_structure), intent(inout) :: gam             ! Basis set and structure descriptor
    integer(ik), intent(in)            :: grid_nrad       ! Basic number of radial points in atomic spheres; actual number of points 
                                                          ! may depend on this value and atom types
    integer(ik), intent(in)            :: grid_nang       ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik), intent(in)            :: grid_outer_nrad ! Basic number of radial points in the outer sphere
    integer(ik), intent(in)            :: grid_outer_nang ! Angular points in the outer sphere
    !
    integer(ik)           :: iat
    real(rk)              :: mol_xyz(3,gam%natoms)
    character(len=20)     :: mol_labels(gam%natoms)
    real(rk)              :: rout
    !
    !  Scaling parameters for the CAP; the code is lifted from caps.f90; it makes no sense
    !  to provide an interface
    !
    cap_width = ma_c / ( 2._rk * cap_delta * cap_kmin )
    cap_scale = 0.5_rk * cap_kmin**2
    write (out,"('   CAP sphere origin = ',3g14.6,' Bohr')") cap_centre
    write (out,"('   CAP stating point = ',g14.6,' Bohr')") cap_r0
    write (out,"('            CAP type = ',a)") trim(cap_type)
    select case (cap_type)
    case default
       flush (out)
       stop 'basis_cap%prepare_for_cap_integration - unrecognized cap_type (1)'
    case ('manolopoulos')
       write (out,"('    Width of the CAP = ',g14.6,' Bohr')") cap_width
       write (out,"('        CAP strength = ',g14.6,' Hartree')") cap_scale
       write (out,"('     CAP upper limit = ',g14.6,' Hartree')") cap_limit
    case ('monomial','atom monomial')
       write (out,"('        CAP strength = ',g14.6,' Hartree Bohr**(-order)')") cap_strength
       write (out,"('  CAP monomial order = ',i0)") cap_order
    case ('moiseyev','atom moiseyev')
       write (out,"(' Switch-on steepness = ',g14.6,' Radian')") cap_lambda
       write (out,"('       Scaling angle = ',g14.6,' Radian')") cap_theta
       write (out,"('Numerical diff. step = ',g14.6,' Bohr')") cap_diff_step
       write (out,"('        Long-range Z = ',g14.6)") cap_mpole(1)
       write (out,"('        Long-range D = ',3(g14.6,1x))") cap_mpole(2:4)
       write (out,"('        Long-range Q = ',3(g14.6,1x))") cap_mpole(5:7)
       write (out,"('                       ',3(g14.6,1x))") cap_mpole(8:10)
       cap_width = 2.0_rk / cap_lambda
    case ('sigmoidal')
       ! We will sort this out later...
    end select
    rout = cap_r0 + 0.5_rk * cap_width
    !
    !  Numerical integration scheme also needs atoms positions and types
    !
    forall (iat=1:gam%natoms) mol_xyz(:,iat) = real(gam%atoms(iat)%xyz,kind=rk) / abohr
    mol_labels = gam%atoms(:gam%natoms)%name
    !
    !  The loop structure below is lifted from correlation_potential_v3.f90
    !
    call GridInitialize(grid,grid_nrad,grid_nang,mol_xyz,mol_labels, &
         outer=.true.,outer_rc=cap_centre,outer_rout=rout,outer_nrad=grid_outer_nrad,outer_nang=grid_outer_nang)
  end subroutine prepape_for_cap_integration

!######################################################################
!
!  CAPs are evaluated using numerical integration on a grid; this seems
!  easier than trying to figure out an analytical solution using RI.
!
!######################################################################

  subroutine cap_evaluate(gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang,fsf,ssf,lsf)
    type(gam_structure), intent(inout) :: gam             ! Basis set and structure descriptor
    integer(ik), intent(in)            :: grid_nrad       ! Basic number of radial points in atomic spheres; actual number of points 
                                                          ! may depend on this value and atom types
    integer(ik), intent(in)            :: grid_nang       ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik), intent(in)            :: grid_outer_nrad ! Basic number of radial points in the outer sphere
    integer(ik), intent(in)            :: grid_outer_nang ! Angular points in the outer sphere
    complex(rk), intent(out)           :: fsf(:,:)        ! spin-free field-interaction Hamiltonian; we need to add CAP terms
    real(rk), intent(out)              :: ssf(:,:)        ! spin-free overlap integrals; done as a consistency check
    real(rk), intent(out)              :: lsf(:,:)        ! spin-free Laplacian integrals; done as a consistency check,
                                                          ! but only when using a non-local CAP
    !
    type(mol_grid)        :: grid
    type(rfc_parameters)  :: rfc_par
    !
    call TimerStart('Evaluate CAPs')
    !
    !  Set up numerical integration grid and CAP parameters
    !
    call prepape_for_cap_integration(grid,gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang)
    !
    select case (cap_type)
    case default
       fsf = 0 ; ssf = 0
       call evaluate_caps_multiplicative(gam,grid,fsf,ssf)
    case ('moiseyev','atom moiseyev')
       rfc_par%cap_type   = trim(cap_type)
       rfc_par%cap_r0     = cap_r0
       rfc_par%cap_lambda = cap_lambda
       rfc_par%cap_slope  = exp((0,1)*cap_theta)
       rfc_par%cap_mpole  = cap_mpole
       rfc_par%cap_efield = cap_efield
       rfc_par%diff_step  = cap_diff_step
       fsf = 0 ; ssf = 0 ; lsf = 0
       call rfc_matrix_elements(gam,grid,rfc_par,fsf,ssf,lsf)
    end select
    call GridDestroy(grid)
    call TimerStop('Evaluate CAPs')
    !
  end subroutine cap_evaluate

!######################################################################
!
!  CAPs are evaluated using numerical integration on a grid; this
!  seems easier than trying to figure out an analytical solution
!  using RI.
!
!######################################################################

  subroutine cap_dump_effective(iu,gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang,vh)
    integer(ik), intent(in)            :: iu              ! Unit number to use for output. Negative unit number will suppress
                                                          ! the output, but still report basic quality diagnostics
    type(gam_structure), intent(inout) :: gam             ! Basis set and structure descriptor
    integer(ik), intent(in)            :: grid_nrad       ! Basic number of radial points in atomic spheres; actual number of points 
                                                          ! may depend on this value and atom types
    integer(ik), intent(in)            :: grid_nang       ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik), intent(in)            :: grid_outer_nrad ! Basic number of radial points in the outer sphere
    integer(ik), intent(in)            :: grid_outer_nang ! Angular points in the outer sphere
    complex(rk), intent(in)            :: vh(:,:)         ! Effective CAP matrix: VH = S^{-1} VC S^{-1}, where VC are the matrix 
                                                          ! elements computed by cap_evaluae
    !
    type(mol_grid)        :: grid
    !
    call TimerStart('Dump effective CAP')
    !
    !  Set up numerical integration grid and CAP parameters
    !
    call prepape_for_cap_integration(grid,gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang)
    !
    call dump_vcap(iu,gam,grid,vh)
    call GridDestroy(grid)
    call TimerStop('Dump effective CAP')
    !
  end subroutine cap_dump_effective

!######################################################################
!
! Evaluation of the coefficients entering into the expansion of the
! CAP potential in terms of the AO basis. Uses numerical integration
! on a grid
!
!######################################################################
  
  subroutine cap_aocoeff(gam,grid_nrad,grid_nang,grid_outer_nrad,&
       grid_outer_nang,coeff)
    type(gam_structure), intent(inout) :: gam             ! Basis set and structure descriptor
    integer(ik), intent(in)            :: grid_nrad       ! Basic number of radial points in atomic spheres; actual number of points 
                                                          ! may depend on this value and atom types
    integer(ik), intent(in)            :: grid_nang       ! Number of angular points; can be 110, 302, or 770 (ie Lebedev grids)
    integer(ik), intent(in)            :: grid_outer_nrad ! Basic number of radial points in the outer sphere
    integer(ik), intent(in)            :: grid_outer_nang ! Angular points in the outer sphere
    complex(rk), intent(out)           :: coeff(:)        ! Expansion coefficients
    !
    type(mol_grid)        :: grid
    !
    call TimerStart('Evaluate CAPs')
    !
    !  Set up numerical integration grid and CAP parameters
    !
    call prepape_for_cap_integration(grid,gam,grid_nrad,grid_nang,grid_outer_nrad,grid_outer_nang)
    ! Perform the numerical integration
    call evaluate_caps_aocoeff(gam,grid,coeff)
    call GridDestroy(grid)
    call TimerStop('Evaluate CAPs')
    !
  end subroutine cap_aocoeff

!######################################################################

  subroutine evaluate_caps_aocoeff(gam,grid,coeff)

    type(gam_structure), intent(inout) :: gam             ! Basis set and structure descriptor
    type(mol_grid), intent(inout)      :: grid
    complex(rk), intent(out)           :: coeff(:)        ! Expansion coefficients

    !
    integer(ik)           :: ib, ipt, npts
    integer(ik)           :: iaol, iaor
    integer(ik)           :: nbatch, nao
    real(rk), pointer     :: xyzw(:,:)       ! Coordinates and weights of the grid points
    real(rk), allocatable :: basval(:,:,:)   ! Values of the basis functions at the grid points
    real(rk), allocatable :: cap(:)          ! Values of the imaginary part of the potential at grid points
    real(rk), allocatable :: coeff_thread(:) ! Per-thread integration buffer
    !
    call GridPointsBatch(grid,'Batches count',count=nbatch)
    nao = gam%nbasis
    !$omp parallel default(none) &
    !$omp& shared(nbatch,nao,grid,gam,coeff,cap_type) &
    !$omp& private(ib,ipt,npts,iaol,iaor,xyzw,basval,cap,coeff_thread)
    nullify(xyzw)
    allocate (coeff_thread(nao))
    coeff_thread = 0
    !$omp do schedule(dynamic)
    grid_batches: do ib=1,nbatch
       !
       !  Get grid points
       !
       call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
       npts = size(xyzw,dim=2)
       if (allocated(basval)) then
          if (size(basval,dim=3)/=npts) deallocate (basval,cap)
       end if
       if (.not.allocated(basval)) allocate(basval(1,nao,npts),cap(npts))
       evaluate_basis_functions: do ipt=1,npts
          call gamess_evaluate_functions(xyzw(1:3,ipt),basval(:,:,ipt),gam)
       end do evaluate_basis_functions
       ! Premultiply integrand by the point weight; saves us an operation later
       evaluate_cap: do ipt=1,npts
          select case (cap_type)
          case default
             stop 'basis_cap%evaluate_caps_multiplicative - unrecognized cap_type'
          case ('manolopoulos')
             cap(ipt) = manolopoulosVcap(xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('monomial')
             cap(ipt) = monomialVcap(xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('atom monomial')
             cap(ipt) = atomMonomialVcap(gam,xyzw(1:3,ipt)) * xyzw(4,ipt)
          case ('sigmoidal')
             cap(ipt) = sigmoidalVcap(gam,xyzw(1:3,ipt)) * xyzw(4,ipt)
          end select
       end do evaluate_cap
       !
       !  Accumulate the integrals
       !
       right_aos: do iaor=1,nao
          coeff_thread(iaor) = coeff_thread(iaor) + sum(basval(1,iaor,:)*cap)
       end do right_aos
       !
    end do grid_batches
    !$omp end do nowait
    if (associated(xyzw)  ) deallocate (xyzw)
    if (allocated (basval)) deallocate (basval,cap)
    !$omp critical
    coeff = coeff + (0,1)*coeff_thread
    !$omp end critical
    deallocate (coeff_thread)
    !$omp end parallel
  end subroutine evaluate_caps_aocoeff
    
!######################################################################

end module basis_cap
