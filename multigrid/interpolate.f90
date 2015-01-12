module interpolate
!
!  Routines for least-squares interpolation on grid. 
!  In principle, arbitrary (real) interpolating functions are supported.
!  In practice, they must be mentioned in basisFunc() below.
!
!  There are two ways of building an iterpolate:
!
!  buildInterpolation will construct an interpolating object from scratch.
!      It is completely stateless.
!
!  beginInterpolation will precompute data, needed for determining interpolation
!      coefficients. 
!
!  fastInterpolation will then determine interpolating coefficients for any 
!      specific data set.
!
!  endInterpolation should be called to release data allocated by beginInterpolation
!
!  If there are many points, with the same mesh spacing and interpolation
!  order, it is more efficient to use the beginInterpolation/fastInterpolation/
!  endInterpolation sequence.
!
  use lapack
  implicit none
!
  integer(ik), parameter :: maxOrder        = 4  ! Max. supported interpolation order
  integer(ik), parameter :: maxCoefficients = 35 ! Max. number of coefficients in interpolating polynomial
  logical, parameter     :: debug           = .false.
  logical, parameter     :: trackNaNs       = .false.
!
  private
  public InterpolationT, buildInterpolation, interp
  public beginInterpolation, endInterpolation, fastInterpolation
!
!  Interpolating types
!
  type InterpolationT
    integer(ik)                :: ncoeff         ! Number of coefficients
    real(rk)                   :: r0(3)          ! Expansion centre (numerical stability reasons)
    complex(rk)                :: coeff(maxCoefficients)
                                                 ! See basisFunc below for the associated functions
  end type InterpolationT
!
!  Data structure for beginInterpolation/fastInterpolation/endInterpolation
!
!  Multiplying int_mat by function values on the right gives interpolation
!  coefficients wrt centre of mass.
!
  real(rk), allocatable, save :: int_mat(:,:)
!
  contains
!
!  buildInterpolation is simply a wrapper for LAPACK cgelss/zgelss
!
!  order can be one of:
!
!    0 - constant fit
!    1 - linear fit
!   -2 - quadratic fit with no off-diagonal terms
!    2 - full quadratic fit
!    3 - cubic fit
!    4 - quartic fit
!
  subroutine buildInterpolation(order,npts,args,weights,values,ip)
    integer(ik), intent(in)           :: order           ! Desired interpolation order
    integer(ik), intent(in)           :: npts            ! Number of data points
    real(rk), intent(in)              :: args   (3,npts) ! Coordinates of points
    real(rk), intent(in)              :: weights(  npts) ! Weights of points
    complex(rk), intent(in)           :: values (  npts) ! Values to interpolate
    type(InterpolationT), intent(out) :: ip              ! Interpolation object

    real(rk)      :: a(npts,maxCoefficients)             ! LSQ problem matrix for *gells
    real(rk)      :: b(max(npts,maxCoefficients),2)      ! LSQ right-hand side
!   integer(ik)   :: ipt, alloc, ifunc, j

    ip%ncoeff = parseOrder(order)
    call buildLSQmat(ip%r0,args,weights,a,ip%ncoeff,npts)
    !
    !  Right-hand side
    !
    b(1:npts,1)      = real (values(:)*weights(:),kind=rk)
    b(1:npts,2)      = aimag(values(:)*weights(:))
    !
    !  Solve it
    !
    call lapack_gelss(a(:,1:ip%ncoeff),b(:,:))
    ip%coeff(1:ip%ncoeff) = cmplx(b(1:ip%ncoeff,1),b(1:ip%ncoeff,2),kind=rk)
    !
    !  A bit of paranoia - try getting back the values we got in
    !
    if (debug) then
      write (out,"(' Interpolating coefficients (slow)')")
      write (out,"(1x,2f15.8/3(1x,2f15.8)/(6(1x,2f15.8)))") ip%coeff(1:ip%ncoeff)
!     write (out,"(1x,a4,3(2x,a30))")    'IPT', '  Original  ', '  Interpolation ', '  Error  '
!     do j=1,npts
!       write (out,"(1x,i4,3(1x,f8.4),3(1x,2f15.8))") j, args(:,j) - ip%r0, values(j), &
!              interp(ip,args(:,j)), values(j) - interp(ip,args(:,j))
!     end do
    end if
  end subroutine buildInterpolation
!
!  Prepare data structures for interpolation. After this call, and until the matchine
!  endInterpolation call, fastInterpolation can be used to calculate interpolating
!  polynomials on the -same- spatial grid (but possibly with different data values).
!
  subroutine beginInterpolation(order,npts,args,weights)
    integer(ik), intent(in)           :: order            ! Desired interpolation order
    integer(ik), intent(in)           :: npts             ! Number of data points
    real(rk), intent(in)              :: args   (3,npts)  ! Coordinates of points
    real(rk), intent(in)              :: weights(  npts)  ! Weights of points

    real(rk)      :: a(npts,maxCoefficients)              ! LSQ problem matrix for *gells
    real(ark)     :: lsq(maxCoefficients,maxCoefficients) ! correlation matrix for the LSQ 
                                                          ! basis functions
    real(rk)      :: r0(3)
    integer(ik)   :: ncoeff, alloc

    ncoeff = parseOrder(order)
    call buildLSQmat(r0,args,weights,a,ncoeff,npts)
    lsq(1:ncoeff,1:ncoeff) = matmul(transpose(real(a(:,1:ncoeff),kind=ark)), &
                                              real(a(:,1:ncoeff),kind=ark))
    !
!   if (trackNaNs) then
!     if (any(isnan(lsq(1:ncoeff,1:ncoeff)))) then
!       write (out,"('beginInterpolation: matmul returned NaNs at ',3f12.5)") r0
!       stop 'beginInterpolation - NaN in matmult'
!     end if
!     write (out,"()") 
!   end if
    !
    call lapack_ginverse(lsq(1:ncoeff,1:ncoeff))
    !
!   if (trackNaNs) then
!     if (any(isnan(lsq(1:ncoeff,1:ncoeff)))) then
!       write (out,"('beginInterpolation: lapack_ginverse returned NaNs at ',3f12.5)") r0
!       stop 'beginInterpolation - NaN in lapack_ginverse'
!     end if
!   end if
    !
    a(:,1:ncoeff) = a(:,1:ncoeff) * spread(weights,ncopies=ncoeff,dim=2)

    allocate (int_mat(ncoeff,npts),stat=alloc)
    if (alloc/=0) then
      write (out,"('beginInterpolation - error ',i8,' in allocate')") alloc
      stop 'beginInterpolation - no memory'
    end if
    int_mat = matmul(dble(lsq(1:ncoeff,1:ncoeff)),transpose(a(:,1:ncoeff)))
  end subroutine beginInterpolation
!
! endInterpolation deletes data structures created by beginInterpolation
!
  subroutine endInterpolation
    deallocate (int_mat)
  end subroutine endInterpolation
!
! fastInterpolation builds interpolation object from data precomputed by 
! beginInterpolation, and new function values at grid points.
!
  subroutine fastInterpolation(r0,values,ip)
    real(rk), intent(in)              :: r0(3)      ! Coordinates of the new "centre of mass"
    complex(rk), intent(in)           :: values (:) ! Values to interpolate
    type(InterpolationT), intent(out) :: ip         ! Interpolation object

    ip%ncoeff = size(int_mat,dim=1)
    ip%r0     = r0
    if (all(abs(values)<1.0e-10_rk)) then
      ip%coeff(1:ip%ncoeff) = 0.0_rk
    else 
      ip%coeff(1:ip%ncoeff) = matmul(int_mat,values)
    end if
    if (debug) then
      write (out,"(' Interpolating coefficients (fast)')")
      write (out,"(1x,2f15.8/3(1x,2f15.8)/(6(1x,2f15.8)))") ip%coeff(1:ip%ncoeff)
    end if
  end subroutine fastInterpolation
!
  function interp(in,arg) result (v)
    type (InterpolationT), intent(in) :: in        ! Interpolation object, set up by buildInterpolation
    real(rk), intent(in)              :: arg(3)    ! Coordinates for the function value
    complex(rk)                       :: v

    real(rk)         :: x(3)      ! arg, adjusted to centrepoint position
    integer(ik)      :: i

    x = arg - in%r0
    v = 0
    do i=1,in%ncoeff
      v = v + in%coeff(i) * basisFunc(i,x)
    end do
!   if (trackNaNs) then
!     if (isnan(abs(v))) then
!       write (out,"(' interp: interpolation at point ',3f12.5,' gave NaN')") arg
!       stop 'interp - NaN'
!     end if
!   end if
  end function interp
!
!  Useful service routines
!
!  Choose the number of coefficients, from the specifiec interpolation order
!    
  function parseOrder(order) result(nc)
    integer(ik), intent(in) :: order
    integer(ik)             :: nc

    select case (order)
      case default
        write (out,"(' Interpolation order ',i3,' not implemented')") order
        stop 'interpolation%parseOrder - bad order'
      case ( 0) ; nc =  1
      case ( 1) ; nc =  4
      case (-2) ; nc =  7
      case ( 2) ; nc = 10
      case (-3) ; nc = 13
      case ( 3) ; nc = 20
      case (-4) ; nc = 23
      case ( 4) ; nc = 35
    end select
  end function parseOrder
!
!  Construct basic matrix, required for the solutiuon of LSQ linear system.
!  It is in the form required by LAPACK *gelss
!
  subroutine buildLSQmat(r0,args,weights,a,ncoeff,npts)
    real(rk), intent(out)   :: r0     (3)                  ! "Centre of mass" position
    integer(ik), intent(in) :: npts                        ! Number of data points
    real(rk), intent(in)    :: args   (3,npts)             ! Coordinates of the interpolation points
    real(rk), intent(in)    :: weights(  npts)             ! Weights at the interpolation points
    real(rk), intent(out)   :: a(npts,maxCoefficients)     ! LSQ problem matrix for *gells
    integer(ik), intent(in) :: ncoeff                      ! Number of interpolating functions

    integer(ik) :: ifunc, j

    !
    !  Get the "centre of mass" position - we'll interpolate wrt this
    !
    r0 = sum(args(:,:)*spread(weights(:),ncopies=3,dim=1),dim=2)/sum(weights(:))
    if (debug) write (out,"('Interp. central point = ',3f15.7)") r0
    !
    !  Build up the LSQ problem matrix
    !
    a(:,1:ncoeff) = 0
    do ifunc=1,ncoeff
      do j=1,npts
        a(j,ifunc) = weights(j)*basisFunc(ifunc,args(:,j)-r0)
      end do
    end do
!   if (trackNaNs) then
!     if (any(isnan(a))) then
!       write (out,"(' buildingLSQmat:  at ',3f12.5,' got NaN')") r0
!       stop 'buildLSQmat - NaN'
!     end if
!   end if
  end subroutine buildLSQmat
!
!  basisFunc is Fortran's way of writing table of function pointers.
!  A good compiler will probably implement it that way, too.
!
  function basisFunc(code,arg) result(v)
    integer(ik), intent(in) :: code
    real(rk), intent(in)    :: arg(3)
    real(rk)                :: v
    real(rk)                :: x, y, z
    
    x = arg(1) ; y = arg(2) ; z = arg(3)
    select case (code)
      case default
        write(out,"(' Basis function ',i8,' is not implemented for interpolation')") code
        stop 'basisFunc - bad code'
      case (1)  ; v = 1.0_rk
      case (2)  ; v = x
      case (3)  ; v = y
      case (4)  ; v = z
      case (5)  ; v = x**2
      case (6)  ; v = y**2
      case (7)  ; v = z**2
      case (8)  ; v = x*y
      case (9)  ; v = x*z
      case (10) ; v = y*z
      case (11) ; v = x**3
      case (12) ; v = y**3
      case (13) ; v = z**3
      case (14) ; v = x**2 * y
      case (15) ; v = x**2 * z
      case (16) ; v = y**2 * x
      case (17) ; v = y**2 * z
      case (18) ; v = z**2 * x
      case (19) ; v = x**2 * y
      case (20) ; v = x * y * z
      case (21) ; v = x**4
      case (22) ; v = y**4
      case (23) ; v = z**4
      case (24) ; v = x**3 * y
      case (25) ; v = x**3 * z
      case (26) ; v = y**3 * x
      case (27) ; v = y**3 * z
      case (28) ; v = z**3 * x
      case (29) ; v = z**3 * y
      case (30) ; v = x**2 * y**2
      case (31) ; v = x**2 * z**2
      case (32) ; v = x**2 * y * z
      case (33) ; v = y**2 * z**2
      case (34) ; v = y**2 * x * z
      case (35) ; v = z**2 * x * y
    end select
  end function basisFunc

end module interpolate
