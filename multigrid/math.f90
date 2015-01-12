module math
!
!  Sundry mathematical definitions.
!
  use accuracy
  use lapack
  !$ use OMP_LIB
  implicit none
  private
  public MathFactorial, MathLogFactorial, MathLegendrePn, MathLegendrePnm
  public MathJacobiPn
  public MathDoubleFactorial, MathLogDoubleFactorial
  public MathBinomial, MathPochhammer
  public MathLogGamma, MathSimpleGamma, Math3J, MathLogSum
  public MathRotationMatrix, MathYJMRotationMatrix
  public MathSetUnitMatrix, MathIsUnitMatrix
  public MathGetQuadrature, MathErf
  public MathBoysF, MathDawsonF
  public MathYLM, MathAllYLM
  public MathDet3x3, MathInv3x3
!
!  Fortran does not allow choosing a specific function based on the type
!  of the result alone. Therefore, we are forced to use some trickery.
!
!  For functions which need such resolution, we will add an extra argument,
!  with the type matching the desired return type. In order to avoid breaking
!  the existing code, we also supply a function without the extra argument;
!  this one will make the "default" choice of the return type.
!
!  Unfortunately, this makes for a bit of clutter.
!
  interface MathBinomial
     module procedure MathBinomialIntegerReal
     module procedure MathBinomialReal
!*qd module procedure MathBinomialIntegerQuad
  end interface MathBinomial
!
  interface MathFactorial
     module procedure MathFactorialReal
!*qd module procedure MathFactorialQuad
  end interface MathFactorial
!
  interface MathDoubleFactorial
     module procedure MathDoubleFactorialReal
!*qd module procedure MathDoubleFactorialQuad
  end interface MathDoubleFactorial
!
  interface MathLogFactorial
     module procedure MathLogFactorialReal
!*qd module procedure MathLogFactorialQuad
  end interface MathLogFactorial
!
  interface MathDet3x3
     module procedure MathDet3x3Real
     module procedure MathDet3x3Complex
  end interface MathDet3x3
!
  interface MathInv3x3
     module procedure MathInv3x3Real
     module procedure MathInv3x3Complex
  end interface MathInv3x3
!
  interface MathIsUnitMatrix
     module procedure MathIsUnitMatrixReal
!*qd module procedure MathIsUnitMatrixQuad
  end interface MathIsUnitMatrix
!
  interface MathSetUnitMatrix
     module procedure MathSetUnitMatrixReal
!*qd module procedure MathSetUnitMatrixQuad
  end interface MathSetUnitMatrix
!
  interface MathBoysF
     module procedure MathBoysFReal
!*qd module procedure MathBoysFQuad
  end interface MathBoysF
!
  integer(ik), parameter       :: factorial_slack = 5         ! Extra factorials to produce while filling the cache
  integer(ik), save            :: factorial_max = -1          ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: factorial_table(:)
  integer(ik), save            :: factorial_quad_max = -1     ! Largest value of the factorials cached in the table
  real(xrk), allocatable, save :: factorial_quad_table(:)
  integer(ik), save            :: log_factorial_max = -1      ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: log_factorial_table(:)
  integer(ik), save            :: log_factorial_quad_max = -1 ! Largest value of the factorials cached in the table
  real(xrk), allocatable, save :: log_factorial_quad_table(:)
  integer(ik), save            :: dfactorial_max = -2         ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: dfactorial_table(:)
  integer(ik), save            :: dfactorial_quad_max = -2    ! Largest value of the double factorials (quad) cached in the table
  real(xrk), allocatable, save :: dfactorial_quad_table(:)
  integer(ik), save            :: log_dfactorial_max = -2     ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: log_dfactorial_table(:)
!
  contains
  !
  !  External interfaces
  !
  function MathFactorialReal(n,kind) result(v)
    integer(ik), intent(in)        :: n
    real(rk), intent(in), optional :: kind
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathFactorialReal - domain error'
    if (n>factorial_max) call fill_factorial_table(n+factorial_slack)
    v = factorial_table(n)
  end function MathFactorialReal
  !
  function MathFactorialQuad(n,kind) result(v)
    integer(ik), intent(in)        :: n
    real(xrk), intent(in)          :: kind
    real(xrk)                      :: v
    !
    if (n<0) stop 'math%MathFactorialQuad - domain error'
    if (n>factorial_quad_max) call fill_factorial_quad_table(n+factorial_slack)
    v = factorial_quad_table(n)
  end function MathFactorialQuad
  !
  function MathLogFactorialReal(n,kind) result(v)
    integer(ik), intent(in)        :: n
    real(rk), intent(in), optional :: kind
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathLogFactorial - domain error'
    if (n>log_factorial_max) call fill_log_factorial_table(n+factorial_slack)
    v = log_factorial_table(n)
  end function MathLogFactorialReal
  !
  function MathLogFactorialQuad(n,kind) result(v)
    integer(ik), intent(in)        :: n
    real(xrk), intent(in)          :: kind
    real(xrk)                      :: v
    !
    if (n<0) stop 'math%MathLogFactorialQuad - domain error'
    if (n>log_factorial_quad_max) call fill_log_factorial_quad_table(n+factorial_slack)
    v = log_factorial_quad_table(n)
  end function MathLogFactorialQuad
  !
  function MathDoubleFactorialReal(n,kind) result(v)
    integer(ik), intent(in)        :: n
    real(rk), intent(in), optional :: kind
    real(rk)                       :: v
    !
    if (n<-1) stop 'math%MathDoubleFactorial - domain error'
    if (n>dfactorial_max) call fill_dfactorial_table(n+factorial_slack)
    v = dfactorial_table(n)
  end function MathDoubleFactorialReal
  !
  function MathDoubleFactorialQuad(n,kind) result(v)
    integer(ik), intent(in)         :: n
    real(xrk), intent(in)           :: kind
    real(xrk)                       :: v
    !
    if (n<-1) stop 'math%MathDoubleFactorialQuad - domain error'
    if (n>dfactorial_quad_max) call fill_dfactorial_quad_table(n+factorial_slack)
    v = dfactorial_quad_table(n)
  end function MathDoubleFactorialQuad
  !
  function MathLogDoubleFactorial(n) result(v)
    integer(ik), intent(in) :: n
    real(rk)                :: v
    !
    if (n<-1) stop 'math%MathLogDoubleFactorial - domain error'
    if (n>log_dfactorial_max) call fill_log_dfactorial_table(n+factorial_slack)
    v = log_dfactorial_table(n)
  end function MathLogDoubleFactorial
  !
  !  Pochhammer function: product of (n) integers starting at (a)
  !
  function MathPochhammer(a,n) result(v)
    integer(ik), intent(in) :: a, n
    real(rk)                :: v
    !
    integer(ik) :: aa     ! Starting integer of the equivalent positive sequence
    logical     :: minus
    !
    if (n<0) stop 'math%MathPochhammer - domain error'
    !
    if (n==0) then
      v = 1._rk
      return
    end if
    !
    if (a<=0) then
      !
      !  Catch sequences containing zero factors: these are always zero
      !
      if (a+n-1>=0) then
        v = 0._rk
        return
      end if
      aa    = -(a+n-1)
      minus = mod(n,2)==1
    else
      aa    = a
      minus = .false.
    end if
    !
    v = MathLogFactorial(aa+n-1)-MathLogFactorial(aa-1)
    v = exp(v)
    if (minus) v = -v
  end function MathPochhammer
  !
  !  Binomial coefficients
  !
  function MathBinomialIntegerReal(n,m,kind) result(cnm)
    integer(ik), intent(in)        :: n, m
    real(rk), intent(in), optional :: kind
    real(rk)                       :: cnm
    !
    if (n<0 .or. m<0 .or. m>n) stop 'MathBinomialIntegerReal - domain error'
    cnm = exp(MathLogFactorial(n)-MathLogFactorial(n-m)-MathLogFactorial(m))
  end function MathBinomialIntegerReal
  !
  function MathBinomialIntegerQuad(n,m,kind) result(cnm)
    integer(ik), intent(in) :: n, m
    real(xrk), intent(in)   :: kind
    real(xrk)               :: cnm
    !
    if (n<0 .or. m<0 .or. m>n) stop 'MathBinomialIntegerQuad - domain error'
    cnm = exp(MathLogFactorial(n,kind)-MathLogFactorial(n-m,kind)-MathLogFactorial(m,kind))
  end function MathBinomialIntegerQuad
  !
  function MathBinomialReal(n,m) result(cnm)
    real(rk), intent(in)   :: n, m
    real(rk)               :: cnm
    complex(rk), parameter :: one = (1._rk,0._rk)
    !
    if (n<0._rk .or. m<0._rk .or. m>n) stop 'MathBinomialReal - domain error'
    cnm = exp(real(MathLogGamma(n+one)-MathLogGamma(n-m+one)-MathLogGamma(m+one),kind=rk))
  end function MathBinomialReal
  !
  !  Evaluate log(exp(a)+exp(b)), where a and b may be too large to fit 
  !  in the exponent range.
  !
  function MathLogSum(a,b) result(c)
    real(rk), intent(in) :: a, b
    real(rk)             :: c
    !
    if (a>=b) then
      c = a + log(1._rk + exp(b-a))
    else
      c = b + log(1._rk + exp(a-b))
    end if
  end function MathLogSum
  !
  !  If all values of Pn or Pnm up to a given order are required, it is much
  !  more efficient to compute them all at once!
  !
  !  Accuracy of the recursions used to calculate Ledendre functions
  !  deteriorates in the vicinity of the polynomial roots, especially
  !  for very high orders (n,m)
  !
  function MathLegendrePn(n,x) result(pn)
    integer(ik), intent(in) :: n    ! Order of the Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate
    real(rk)                :: pn   ! Value of the Legenre polynomial
    !
    real(rk) :: tab(0:n)
    !
    call legendrePn_table(n,x,tab)
    pn = tab(n)
  end function MathLegendrePn
  !
  function MathLegendrePnm(n,m,x) result(pnm)
    integer(ik), intent(in) :: n, m ! Order of the associated Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate, abs(x)<=1
    real(rk)                :: pnm  ! Value of the associated Legenre polynomial
    !
    real(rk) :: tab(0:n,0:m)
    !
    call legendrePnm_table(n,m,x,tab)
    pnm = tab(n,m)
  end function MathLegendrePnm
  !
  !  Use recurrence with respect to degree, see Abramowitz & Stegun 22.7.1
  !
  function MathJacobiPn(n,alp,bet,x) result(pn)
    integer(ik), intent(in) :: n        ! Order of the polynomial
    real(rk), intent(in)    :: alp, bet ! Powers of the weight function, must be > -1
    real(rk), intent(in)    :: x        ! Coordinate, abs(x)<=1
    real(rk)                :: pn
    !
    real(rk) :: tab(0:n)
    !
    call jacobiPn_table(n,alp,bet,x,tab)
    pn = tab(n)
  end function MathJacobiPn
  !
  !  See Abramowitz & Stegun, 22.18 for the evaluation procedure
  !  This routine works very well for low orders n, but becomes unstable
  !  for higher order. It loses 4 decimal digits of significance at n=6,
  !  deteriorating to complete loss of significance at n=20
  !
! function MathJacobiPnUnstable(n,alp,bet,x) result(pn)
!   integer(ik), intent(in) :: n        ! Order of the polynomial
!   real(rk), intent(in)    :: alp, bet ! Powers of the weight function, must be > -1
!   real(rk), intent(in)    :: x        ! Coordinate, abs(x)<=1
!   real(rk)                :: pn
!   !
!   integer(ik) :: m
!   real(rk)    :: am  ! Auxiliary function; see A&S 22.18
!   real(rk)    :: fx  ! f(x)
!   !
!   am = 1._rk     ! am(n)
!   fx = 1._rk - x ! See the second table in A&S 22.18
!   recurse_am: do m=n,1,-1
!     am = 1._rk - am*bm(m)*fx/cm(m)  ! Calculates a(m-1)
!   end do recurse_am
!   pn = dn()*am
!   !
!   contains
!   real(rk) function dn()
!     dn = MathBinomial(n+alp,real(n,kind=rk))
!   end function dn
!   real(rk) function bm(m)
!     integer(ik), intent(in) :: m
!     !
!     bm = (n-m+1._rk)*(alp+bet+n+m)
!   end function bm
!   real(rk) function cm(m)
!     integer(ik), intent(in) :: m
!     !
!     cm = 2._rk*m*(alp+m)
!   end function cm
! end function MathJacobiPnUnstable
  !
  !  Computes Wigner 3J symbols. The code below is a direct implementation
  !  of L&L 3J formulae. The accuracy of this routine is reduced relative to
  !  that is theoretically possible, due to the use of logarithms. The routine
  !  I had in MNDO99 is more accurate and can handle broader range of J values.
  !
  function Math3J(j1,j2,j3,m1,m2,m3) result(v)
    integer(ik), intent(in) :: j1, j2, j3  ! / J1 J2 J3 \ 3-j
    integer(ik), intent(in) :: m1, m2, m3  ! \ M1 M2 M3 / 
    real(rk)                :: v
    !
    integer(ik) :: ij0, ij1, ij2, ij3, im1a, im1b, im2a, im2b, im3a, im3b
    integer(ik) :: t1, t2, t3
    integer(ik) :: z, minz, maxz
    real(rk)    :: logscale, logterm
    !
    !  Before we do anything, check whether this 3J symbol satisfies the
    !  vector addition constraints
    !
    ij0  =   j1 + j2 + j3 + 1
    ij1  =   j1 + j2 - j3
    ij2  =   j1 - j2 + j3
    ij3  = - j1 + j2 + j3
    im1a =   j1 - m1 ; im1b = j1 + m1
    im2a =   j2 - m2 ; im2b = j2 + m2
    im3a =   j3 - m3 ; im3b = j3 + m3
    if (ij1<0 .or. ij2<0 .or. ij3<0 .or. im1a<0 .or. im1b<0 .or. im2a<0 .or. im2b<0 .or. im3a<0 .or. im3b<0 .or. m1+m2+m3/=0) then
      v = 0
      return
    end if
    !
    logscale = MathLogFactorial(ij1)  + MathLogFactorial(ij2)  + MathLogFactorial(ij3)  &
             + MathLogFactorial(im1a) + MathLogFactorial(im1b) + MathLogFactorial(im2a) &
             + MathLogFactorial(im2b) + MathLogFactorial(im3a) + MathLogFactorial(im3b) &
             - MathLogFactorial(ij0)
    logscale = 0.5_rk * logscale
    !
    t1   = j2 - j3 - m1
    t2   = j1 + m2 - j3
    t3   = j1 - j2 - m3
    minz = max(0_ik,t1,t2)
    maxz = min(ij1,im1a,im2b)
    v = 0
    sum_terms: do z=minz,maxz,1
      logterm = logscale - MathLogFactorial(z)      - MathLogFactorial(ij1-z)  - MathLogFactorial(im1a-z) &
                         - MathLogFactorial(im2b-z) - MathLogFactorial(z-t1)   - MathLogFactorial(z-t2)
      if (abs(logterm)>=max_exp) then
        write (out,"('Math3J: Intermediate logarithm ',g12.5,' exceeds the real(rk) dynamic range.')") logterm
        write (out,"('Math3J: The 3J arguments were: ',6i10)") j1, j2, j3, m1, m2, m3
        stop 'math%Math3J - exceeded dynamic range'
      end if
      if (mod(z+t3,2)==0) then
        v = v + exp(logterm)
      else
        v = v - exp(logterm)
      end if
    end do sum_terms
    !
  end function Math3J
  !
  !  Given Euler angles, construct rotation matrix for coordinate axes 
  !  in 3D space.  The Euler angles are defined as follows:
  !   1. Rotate coordinate axes by alpha around the Z axis
  !   2. Rotate axes by beta around the new Y axis
  !   3. Rotate axes by gamma around the new Z axis
  !  This definition of the Euler angles matches the definition from
  !  section 58 of L&L III.
  !
  !  Note that prior to June 10th, 2010 the definition of all three Euler 
  !  angles in this routine used a wrong sign, corresponding to anti-screw
  !  rotation sense. Thanks, Mike.
  !
  !  If you would like to rotate an object instead of the axes, take the
  !  transpose or (equivalently) replace (alpha,beta,gamma) with
  !  (-gamma,-beta,-alpha).
  !
  !  Some useful special cases are:
  !
  !    a     b     c  
  !   ---   ---   --- 
  !    x     0     0   Rotate coordinate system by x around the Z axis
  !    0     0     x   Rotate coordinate system by x around the Z axis
  !    0     x     0   Rotate coordinate system by x around the Y axis
  !  -pi/2   x   pi/2  Rotate coordinate system by x around the X axis
  !      
  subroutine MathRotationMatrix(euler_angles,mat)
    real(ark), intent(in)  :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
    real(ark), intent(out) :: mat(3,3)        ! Rotation matrix
    !
    real(ark) :: a, b, g, rma(3,3), rmb(3,3), rmg(3,3)
    !
    a = euler_angles(1)
    b = euler_angles(2)
    g = euler_angles(3)
    !
    rma(1,:) = (/  cos(a),  sin(a),   0._rk /)
    rma(2,:) = (/ -sin(a),  cos(a),   0._rk /)
    rma(3,:) = (/  0._rk,    0._rk,   1._rk /)
    !
    rmb(1,:) = (/  cos(b),   0._rk, -sin(b) /)
    rmb(2,:) = (/  0._rk,    1._rk,   0._rk /)
    rmb(3,:) = (/  sin(b),   0._rk,  cos(b) /)
    !
    rmg(1,:) = (/  cos(g),  sin(g),  0._rk  /)
    rmg(2,:) = (/ -sin(g),  cos(g),  0._rk  /)
    rmg(3,:) = (/  0._rk,    0._rk,  1._rk  /)
    !
    mat = matmul(rmg,matmul(rmb,rma))
  end subroutine MathRotationMatrix
  !
  !  Rotation matrix for angular momentum eigenfunctions, following L&L III Eq. 58.10
  !  Both integer and half-integer J values are OK.
  !
  !  The resulting rotation matrix is accurate to 2ulp for multiplicities up to 6,
  !  with error increasing to 4ulp for multiplicity 20. It loses about 11 decimal places
  !  of accuracy for multiplicity 81, and overflows IEEE double at higher multiplicities.
  !
  subroutine MathYJMRotationMatrix(euler_angles,mult,mat)
    real(ark), intent(in)     :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
                                                 ! See comments in MathRotationMatrix
    integer(ik), intent(in)   :: mult            ! Multipliplicity of the angular-momentum state,
                                                 ! mult = 2*j+1
    complex(ark), intent(out) :: mat(:,:)        ! Rotation matrix
    !
    real(ark)    :: a, b, g
    real(ark)    :: cosb2, sinb2
    complex(ark) :: expa2, expg2
    integer(ik)  :: j2, m2, mp2
    integer(ik)  :: im, imp
    !
    if (mult<1) then
      stop 'math%MathYJMRotationMatrix - multiplicity: domain error'
    end if
    if (size(mat,dim=1)/=mult .or. size(mat,dim=2)/=mult) then
      stop 'math%MathYJMRotationMatrix - rotation matrix dimensions do not match multiplicity'
    end if
    !
    a = euler_angles(1)
    b = euler_angles(2)
    g = euler_angles(3)
    !
    !  We need to take special care when angle beta approaches n*pi. For these angles,
    !  
    !
    sinb2 = sin(0.5_ark*b)
    cosb2 = cos(0.5_ark*b)
    !
    expa2 = exp(cmplx(0._ark,0.5_ark*a,kind=ark))
    expg2 = exp(cmplx(0._ark,0.5_ark*g,kind=ark))
    !
    j2  = mult - 1
    mp2 = -j2
    row_mp: do imp=1,mult
      m2 = -j2
      column_m: do im=1,mult
        mat(imp,im) = sqrt(hf(j2+mp2)*hf(j2-mp2)/(hf(j2+m2)*hf(j2-m2))) &
                    * modJacobiPn((j2-mp2)/2,(mp2-m2)/2,(mp2+m2)/2,sinb2,cosb2) &
                    * expa2**m2 * expg2**mp2
        m2 = m2 + 2
      end do column_m
      mp2 = mp2 + 2
    end do row_mp
    !
    contains
    !
    !  (n/2)! where n must be even
    !
    function hf(n) result(f)
      integer(ik), intent(in) :: n ! Must be even
      real(ark)               :: f
      !
      if (mod(n,2)/=0) stop 'math%MathYJMRotationMatrix%hf - domain error'
      f = MathFactorial(n/2)
    end function hf
    !
    !  Specialized derivative of Jacobi P polynomial:
    !
    !    y^a x^b JacobiP(n,a,b,x^2-y^2)
    !
    !  where x^2+y^2 is equal to 1. Care should be taken in evaluating this function
    !  when either a or b are negative: standard recursion with respect to degree 
    !  becomes undefined in this case.
    !
    !  As the result, we have to use series expansion around +1/-1 argument to
    !  evaluate this function.
    !
    function modJacobiPn(n,a,b,y,x) result(jp)
      integer(ik), intent(in) :: n, a, b ! Parameters of the Jacobi P polynomial
      real(ark), intent(in)   :: y, x    ! Coordinate
      real(ark)               :: jp
      !
      integer(ik) :: k
      !
      if (abs(y)<abs(x)) then
        !
        !  Small y, expand JacobiP around z=+1
        !
        jp = 0
        expand_plus1: do k=max(0,-a),n
          jp = jp + MathPochhammer(a+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * y**(2*k+a) / MathFactorial(k)
        end do expand_plus1
        jp = x**b * jp / MathFactorial(n)
!       write (out,"(a,3(1x,i3),3(1x,f35.25))") 'A: ',n,a,b,y,x,jp
      else 
        !
        !  Small x, expand JacobiP around z=-1
        !
        jp = 0
        expand_minus1: do k=max(0,-b),n
          jp = jp + MathPochhammer(b+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * x**(2*k+b) / MathFactorial(k)
        end do expand_minus1
        jp = y**a * jp / MathFactorial(n)
        if (mod(n,2)==1) jp = -jp
!       write (out,"(a,3(1x,i3),3(1x,f35.25))") 'B: ',n,a,b,y,x,jp
!       write (out,"()")
      endif
!       !
!       !  The general case; do not use
!       !
!       jp = y**a * x**b * MathJacobiPn(n,real(a,kind=rk),real(b,kind=rk),x**2-y**2)
    end function modJacobiPn
  end subroutine MathYJMRotationMatrix
  !
  !  Initialize unit matrix
  !
  subroutine MathSetUnitMatrixReal(m)
    real(ark), intent(out) :: m(:,:) ! Matrix to initialize
    !
    integer(ik) :: i
    !
    if (size(m,dim=1)/=size(m,dim=2)) then
      write (out,"('math%MathSetUnitMatrix - argument is not a square matrix?!')") 
      stop 'math%MathSetUnitMatrix - argument is not a square matrix?!'
    end if
    !
    m = 0
    do i=1,size(m,dim=1)
      m(i,i) = 1
    end do
  end subroutine MathSetUnitMatrixReal

  subroutine MathSetUnitMatrixQuad(m)
    real(xrk), intent(out) :: m(:,:) ! Matrix to initialize
    !
    integer(ik) :: i
    !
    if (size(m,dim=1)/=size(m,dim=2)) then
      write (out,"('math%MathSetUnitMatrix - argument is not a square matrix?!')") 
      stop 'math%MathSetUnitMatrix - argument is not a square matrix?!'
    end if
    !
    m = 0
    do i=1,size(m,dim=1)
      m(i,i) = 1
    end do
  end subroutine MathSetUnitMatrixQuad
  !
  !  Check matrix for unity
  !
  function MathIsUnitMatrixReal(m) result(isunit)
    real(ark), intent(in) :: m(:,:) ! Matrix to check
    logical               :: isunit
    !
    real(kind=kind(m))    :: eps 
    integer(ik)           :: i, j
    !
    eps = spacing(real(2,kind=kind(m)))
    isunit = .false.
    if (size(m,dim=1)/=size(m,dim=2)) return
    !
    scan_rows: do i=1,size(m,dim=2)
      if (abs(m(i,i)-1)>eps) return
      scan_upper_columns: do j=1,i-1
        if (abs(m(j,i))>eps) return
      end do scan_upper_columns
      scan_lower_columns: do j=i+1,size(m,dim=1)
        if (abs(m(j,i))>eps) return
      end do scan_lower_columns
    end do scan_rows
    isunit = .true.
  end function MathIsUnitMatrixReal

  function MathIsUnitMatrixQuad(m) result(isunit)
    real(xrk), intent(in) :: m(:,:) ! Matrix to check
    logical               :: isunit
    !
    real(kind=kind(m))    :: eps 
    integer(ik)           :: i, j
    !
    eps = spacing(real(2,kind=kind(m)))
    isunit = .false.
    if (size(m,dim=1)/=size(m,dim=2)) return
    !
    scan_rows: do i=1,size(m,dim=2)
      if (abs(m(i,i)-1)>eps) return
      scan_upper_columns: do j=1,i-1
        if (abs(m(j,i))>eps) return
      end do scan_upper_columns
      scan_lower_columns: do j=i+1,size(m,dim=1)
        if (abs(m(j,i))>eps) return
      end do scan_lower_columns
    end do scan_rows
    isunit = .true.
  end function MathIsUnitMatrixQuad
  !
  !  Calculate Gaussian quadrature rules. See: Golub and Welsch, "Calculation of Gaussian
  !  quadrature rules," mathematics of computation 23 (april, 1969), pp. 221-230
  !
  !  WARNING: Only the Gauss-Hermite rule has been verified against the tables
  !
  subroutine MathGetQuadrature(type,order,x,w,alpha)
    character(len=*), intent(in)   :: type  ! Desired quadrature type; see below for the allowed types
    integer(ik), intent(in)        :: order ! Desired quadrature order
    real(rk), intent(out)          :: x(:)  ! Integration points
    real(rk), intent(out)          :: w(:)  ! Integration weights
    real(rk), intent(in), optional :: alpha ! Laguerre alpha
    !
    integer(ik)           :: i, alloc
    real(rk), allocatable :: a(:), b(:)   ! Diagonal and sub-diagonal values for the eigensystem
    real(rk), allocatable :: z(:,:)       ! Space for the eigenvectors
    real(rk)              :: mu0          ! Zeroth momentum of the weight function
    !
    if (size(x)/=order .or. size(w)/=order) then
      stop 'math%MathGetQuadrature - inconsistent buffer sizes'
    end if
    !
    allocate (a(order),b(order-1),z(order,order),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i8,' allocating scratch for order-',i8,' quadrature')") alloc, order
      stop 'math%MathGetQuadrature - no memory'
    end if
    !
    !  The select case below is the only part of the routine which depends on the
    !  quadrature required; everything else is the same for all quadratures.
    !
    select case (type)
      case default
        write (out,"('MathGetQuadrature called for unknown quadrature type ',a)") trim(type)
        stop 'math%MathGetQuadrature - invalid quadrature type'
      case ('Legendre')
        !
        !  (-1,+1), 1.0
        !
        a = 0
        legendre: do i=1,order-1
          b(i) = i/sqrt(4*i**2 - 1.0_rk)
        end do legendre
        mu0 = 2._rk
      case ('Chebyshev T')  ! Chebychev of the 1st kind
        !
        !  (-1,+1), (1-X**2)**-0.5
        !
        a     = 0
        b(1)  = sqrt(0.5_rk)
        b(2:) = 0.5_rk
        mu0   = pi
      case ('Chebyshev U')  ! Chebychev of the 2nd kind
        !
        !  (-1,+1), (1-X**2)**0.5
        !
        a   = 0
        b   = 0.5_rk
        mu0 = 0.5_rk*pi
      case ('Hermite')
        !
        !  (-Infinity,+Infinity), EXP(-X**2)
        !
        a = 0
        hermite: do i=1,order-1
          b(i) = sqrt(0.5_rk*i)
        end do hermite
        mu0 = sqrt(pi)
      case ('Laguerre')
        !
        !  (0,+Infinity), EXP(-X) * X**ALPHA
        !
        if (.not.present(alpha)) stop 'MathGetQuadrature(Laguerre) - missing alpha'
        laguerre: do i=1,order-1
          a(i) = 2*i - 1 + alpha
          b(i) = sqrt(i*(i+alpha))
        end do laguerre
        a(order) = 2*order - 1 + alpha
        mu0 = real(MathSimpleGamma(cmplx(1._rk + alpha,0._rk,kind=rk)),kind=rk)
    end select
    !
    call lapack_stev(a,b,z)  
    !
    !  Copy integration abscissas and weights, and we are done
    !
    x = a
    w = mu0 * z(1,:)**2
    !
    deallocate (a,b,z)
  end subroutine MathGetQuadrature
  !
  !  Evaluate a given spherical harmonic. If you need to evaluate the same spherical
  !  harmonic at multiple points, it is preferable to use FLharmonics() in fields.f90
  !
  function MathYLM(l,m,dir) result(ylm)
    integer(ik), intent(in) :: l, m
    real(rk), intent(in)    :: dir(3)   ! (Unnormalized) direction vector
    complex(rk)             :: ylm
    !
    complex(rk) :: harm_fact
    real(rk)    :: r, ct, xymod
    complex(rk) :: xy
    !
    harm_fact = (-1)**((m+abs(m))/2)
    harm_fact = harm_fact * (0._rk,1._rk)**l
    harm_fact = harm_fact * sqrt((2*l+1)/fourpi)
    harm_fact = harm_fact * sqrt(MathFactorial(l-abs(m))/MathFactorial(l+abs(m)))
    !
    r  = sqrt(sum(dir**2))
    !
    ct = dir(3)/r
    if (m>0) then
      xy = cmplx(dir(1), dir(2),kind=rk)
    else
      xy = cmplx(dir(1),-dir(2),kind=rk)
    end if
    xymod = abs(xy)
    if (xymod>0._rk) then
      xy = xy / xymod
    else
      xy = 1._rk
    end if
    ylm = harm_fact * MathLegendrePnm(l,abs(m),ct) * xy**abs(m)
  end function MathYLM
  !
  !  Calculate all spherical harmonics up to angular momentum l_max
  !
  subroutine MathAllYLM(l_max,dir,ylm)
    integer(ik), intent(in)  :: l_max   ! Desired maximum angular momentum
    real(rk), intent(in)     :: dir(3)  ! (Unnormalized) direction vector
    complex(rk), intent(out) :: ylm(-l_max:l_max,0:l_max)
    !
    integer(ik) :: l, m
    complex(rk) :: harm_fact
    real(rk)    :: r, ct, xymod
    complex(rk) :: xy_plus, xy_minus
    real(rk)    :: legendre_tab(0:l_max,0:l_max)
    complex(rk) :: xy_powm(-l_max:l_max)
    !
    r  = sqrt(sum(dir**2))
    ct = dir(3)/r
    call legendrePnm_table(l_max,l_max,ct,legendre_tab)
    !
    xy_plus  = cmplx(dir(1), dir(2),kind=rk)
    xy_minus = cmplx(dir(1),-dir(2),kind=rk)
    xymod    = abs(xy_plus)
    if (xymod>0._rk) then
      xy_plus  = xy_plus  / xymod
      xy_minus = xy_minus / xymod
    else
      xy_plus  = 1._rk
      xy_minus = 1._rk
    end if
    !
    xy_pow1: do m=-l_max,-1
      xy_powm(m) = xy_minus ** abs(m)
    end do xy_pow1
    xy_powm(0) = 1._rk
    xy_pow2: do m=1,l_max
      xy_powm(m) = xy_plus ** abs(m)
    end do xy_pow2
    !
    tab_l: do l=0,l_max
      tab_m: do m=-l,l
        harm_fact = (-1)**((m+abs(m))/2) * (0._rk,1._rk)**l * sqrt((2*l+1)/fourpi) &
                  * sqrt(MathFactorial(l-abs(m))/MathFactorial(l+abs(m)))
        ylm(m,l) = harm_fact * legendre_tab(l,abs(m)) * xy_powm(m)
      end do tab_m
    end do tab_l
  end subroutine MathAllYLM
  !
  !  Auxiliary functions
  !
  subroutine legendrePn_table(nmax,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: invn
    !
    if (nmax<0) stop 'math%legendreP_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = x
    !
    n_recursion: do n=2,nmax
      invn = 1._rk / n
      pn(1+n) = (2._rk-invn)*x*pn(1+(n-1)) - (1._rk-invn)*pn(1+(n-2))
    end do n_recursion
  end subroutine legendrePn_table
  !
  subroutine legendrePnm_table(nmax,mmax,x,pnm)
    integer(ik), intent(in) :: nmax     ! Maximum order of the Legendre polynomials desired
    integer(ik), intent(in) :: mmax     ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x        ! Coordinate at which P values are needed, abs(x) must be <=1
    real(rk), intent(out)   :: pnm(:,:) ! Values of LegendreP from n,m=0 to n,m=nmax,mmax
                                        ! n is the first subscript; m is the second subscript
    !
    integer(ik) :: n, m
    real(rk)    :: sqfac ! sqrt(1-x**2)
    !
    sqfac = 1._rk - x**2
    if (sqfac<0) stop 'math%legendrePnm_table - domain error'
    sqfac = sqrt(sqfac)
    !
    call legendrePn_table(nmax,x,pnm(:,1))
    if (mmax<1) return
    !
    !  Special case for m=1: recursion is truncated for n=1, m=1
    !
    pnm(1+0,1+1) = 0._rk
    if (nmax>=1) pnm(1+1,1+1) = -sqfac
    m = 1
    n1_recursion: do n=2,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
    end do n1_recursion
    !
    m_recursion: do m=2,mmax
      pnm(1+0:1+min(nmax,(m-1)),1+m) = 0._rk
      nm_recursion: do n=m,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
      end do nm_recursion
    end do m_recursion
  end subroutine legendrePnm_table
  !
  !  Recurrence from A&S 22.7.1. These recurrences are not completely numerically
  !  stable, and lose up to 3 decimal digits for n>10.
  !  Unfortunately, these recurrences do not work for negative integer a and b:
  !  sooner or later, we always hit division by zero. Oops.
  !
  subroutine jacobiPn_table(nmax,a,b,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: a, b  ! alpha and beta parameters of the weight function
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: a1n, a2n, a3n, a4n
    !
    if (nmax<0) stop 'math%jacobiPn_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = 0.5_rk*((a-b) + (2._rk+a+b)*x)
    !
    !  A&S 22.7 recursions are written for a polynomial of degree (n+1).
    !  To avoid confusion, we'll do the same. Note that the a3n coefficient in A&S 22.7.1
    !  is written using non-standard series notation, which is likely to cause confusion ...
    !  Keep in mind that n-th degree polynomial is at pn(n+1)
    !
    n_recursion: do n=1,nmax-1
      a1n = 2 * (n+1)*(n+a+b+1) * (2*n+a+b)
      a2n = (2*n+a+b+1) * (a**2-b**2)
      a3n = (2*n+a+b) * (2*n+a+b+1) * (2*n+a+b+2)
      a4n = 2 * (n+a) * (n+b) * (2*n+a+b+2)
      !
      pn(n+2) = ( (a2n+a3n*x)*pn(n+1) - a4n*pn(n) ) / a1n
    end do n_recursion
  end subroutine jacobiPn_table
  !
  subroutine fill_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathFactorial'
    !$ end if
    !
    if (factorial_max>=0) then
      deallocate (factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 1._rk
    !
    allocate (factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element factorial table')") & 
             alloc, nmax
      stop 'math%fill_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      factorial_table(n) = fac
      n = n + 1
      !
      if (huge(fac)/n<=fac) then
        write (out,"(1x,i10,'! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_factorial_table - range exceeded'
      end if
      fac = fac * n
    end do fill_factorials
    factorial_table(n) = fac
    !
    factorial_max = nmax
  end subroutine fill_factorial_table
  !
  subroutine fill_factorial_quad_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(xrk)               :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_quad_table - unsafe call to MathFactorial'
    !$ end if
    !
    if (factorial_quad_max>=0) then
      deallocate (factorial_quad_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating factorial table')") alloc
        stop 'math%fill_factorial_quad_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 1._xrk
    !
    allocate (factorial_quad_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element factorial table')") & 
             alloc, nmax
      stop 'math%fill_factorial_quad_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      factorial_quad_table(n) = fac
      n = n + 1
      !
      if (huge(fac)/n<=fac) then
        write (out,"(1x,i10,'! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_factorial_quad_table - range exceeded'
      end if
      fac = fac * n
    end do fill_factorials
    factorial_quad_table(n) = fac
    !
    factorial_quad_max = nmax
  end subroutine fill_factorial_quad_table
  !
  subroutine fill_log_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathLogFactorial'
    !$ end if
    !
    if (log_factorial_max>=0) then
      deallocate (log_factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 0._rk
    !
    allocate (log_factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      log_factorial_table(n) = fac
      n = n + 1
      !
      fac = fac + log(real(n,kind=rk))
    end do fill_factorials
    log_factorial_table(n) = fac
    !
    log_factorial_max = nmax
  end subroutine fill_log_factorial_table
  !
  subroutine fill_log_factorial_quad_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(xrk)               :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table_quad - unsafe call to MathLogFactorial'
    !$ end if
    !
    if (log_factorial_quad_max>=0) then
      deallocate (log_factorial_quad_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-factorial table')") alloc
        stop 'math%fill_log_factorial_quad_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 0._xrk
    !
    allocate (log_factorial_quad_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_factorial_quad_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      log_factorial_quad_table(n) = fac
      n = n + 1
      !
      fac = fac + log(real(n,kind=xrk))
    end do fill_factorials
    log_factorial_quad_table(n) = fac
    !
    log_factorial_quad_max = nmax
  end subroutine fill_log_factorial_quad_table
  !
  subroutine fill_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathDoubleFactorial'
    !$ end if
    !
    if (dfactorial_max>=0) then
      deallocate (dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element double factorial table')") & 
             alloc, nmax
      stop 'math%fill_dfactorial_table - allocate'
    end if
    !
    dfactorial_table(-1:1) = 1._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      if (huge(1._rk)/n<=dfactorial_table(n-2)) then
        write (out,"(1x,i10,'!! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_dfactorial_table - range exceeded'
      end if
      dfactorial_table(n) = dfactorial_table(n-2) * n
      n = n + 1
    end do fill_factorials
    !
    dfactorial_max = nmax
  end subroutine fill_dfactorial_table
  !
  subroutine fill_dfactorial_quad_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_quad_table - unsafe call to MathDoubleFactorial'
    !$ end if
    !
    if (dfactorial_quad_max>=0) then
      deallocate (dfactorial_quad_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating double factorial table')") alloc
        stop 'math%fill_dfactorial_quad_table - deallocate'
      end if
    end if
    !
    allocate (dfactorial_quad_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element double factorial table')") & 
             alloc, nmax
      stop 'math%fill_dfactorial_quad_table - allocate'
    end if
    !
    dfactorial_quad_table(-1:1) = 1._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      if (huge(1._xrk)/n<=dfactorial_quad_table(n-2)) then
        write (out,"(1x,i10,'!! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_dfactorial_quad_table - range exceeded'
      end if
      dfactorial_quad_table(n) = dfactorial_quad_table(n-2) * n
      n = n + 1
    end do fill_factorials
    !
    dfactorial_quad_max = nmax
  end subroutine fill_dfactorial_quad_table
  !
  subroutine fill_log_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathLogDoubleFactorial'
    !$ end if
    !
    if (log_dfactorial_max>=0) then
      deallocate (log_dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (log_dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-double factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_dfactorial_table - allocate'
    end if
    !
    log_dfactorial_table(-1:1) = 0._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      log_dfactorial_table(n) = log_dfactorial_table(n-2) + log(real(n,kind=rk))
      n = n + 1
    end do fill_factorials
    !
    log_dfactorial_max = nmax
  end subroutine fill_log_dfactorial_table
  !
  function MathSimpleGamma(z) result(v)
    complex(rk), intent(in) :: z ! Argument of the gamma function
    complex(rk)             :: v ! Gamma function
    !
    v = exp(MathLogGamma(z))
  end function MathSimpleGamma
  !
  !  Evaluation of the gamma function of complex argument seems to be a bit
  !  complicated. Let's do something very straightforward, and hopefully accurate
  !  enough:
  !    1. Use recurrence relations to _increase_ the real part of the argument
  !       until modulus is large enough to let us use asymprotic expressions
  !    2. Then apply the Stirling formula to get the result.
  !
  !  The resulting routine is not fast, but is reasonably accurate: for
  !  arguments between 0 and 10, the relative error does not exceed 5e-14
  !
  function MathLogGamma(z) result(v)
    complex(rk), intent(in) :: z ! Argument of the gamma function
    complex(rk)             :: v ! Gamma function
    !
    !  Stirling formula coefficients
    !
    real(rk), parameter :: c01 =               1._rk/                12._rk
    real(rk), parameter :: c02 =               1._rk/               288._rk
    real(rk), parameter :: c03 =            -139._rk/             51840._rk
    real(rk), parameter :: c04 =            -571._rk/           2488320._rk
    real(rk), parameter :: c05 =          163879._rk/         209018880._rk
    real(rk), parameter :: c06 =         5246819._rk/       75246796800._rk
    real(rk), parameter :: c07 =      -534703531._rk/      902961561600._rk
    real(rk), parameter :: c08 =     -4483131259._rk/    86684309913600._rk
    real(rk), parameter :: c09 = 432261921612371._rk/514904800886784000._rk
    !
    complex(rk) :: zr    ! Reduced argument
    complex(rk) :: logs  ! Scaling coefficient needed to reduce the result of Stirling formula
    complex(rk) :: logv  ! Logarithm of the result
    real(rk)    :: zcut  ! Minimal safe argument for the Stirling formula
    complex(rk) :: vs    ! Series part of the Stirling formula
    !
    !  To get accurate results from Stirling formula, we must make sure that
    !  the modulus of the argument is large enough to make the last contribution
    !  to the expansion small enough.
    !
    zcut = (c09/spacing(1._rk))**(1._rk/9._rk)
    ! write (out,"(' zcut = ',f25.15)") zcut
    logs = 0._rk
    zr   = z
    inflate_z: do while(abs(zr)<zcut) 
      logs = logs + log(zr*(zr+1)*(zr+2))
      zr   = zr + 3
    end do inflate_z
    ! write (out,"(' zr = ',2(1x,g25.15),' logs = ',2(1x,g25.15))") zr, logs
    !
    !  It is safe to use Stirling formula now
    !
    vs = 1._rk + (c01 + (c02 + (c03 + (c04 + (c05 + (c06 + (c07 + (c08 + c09/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr)/zr
    ! write (out,"(' vs = ',2(1x,g25.15))") vs
    logv = log(sqrt2pi) - zr + (zr-0.5_rk)*log(zr) + log(vs) - logs
    ! write (out,"(' logv = ',2(1x,g25.15))") logv
    !
    ! v = exp(logv)
    v = logv
  end function MathLogGamma
  !
  !  If intrinsic erf() is not available (it is first mandatated by Fortran 2008), feel
  !  free to use the appropriate Slatec routine in erf.f/derf.f
  !
  function MathErf(x) result(v)
    real(rk), intent(in) :: x
    real(rk)             :: v
    intrinsic :: erf
    !
    v = erf (x)
  end function MathErf
  !
  !  MathBoysF implements the integral of x**(2*n)*Exp(-t*x**2) for x between 0 and 1. 
  !  Depending on the arguments, the implementation will switch between:
  !   a. (large t) Asymptotic form of the integral (with the upper bound replaced by infinity)
  !   b. (small t) Explicit power series
  !   c. Upward recursion in n.
  !  This implementation should produce the integral to machine accuracy, but does not have
  !  efficiency as its primary goal. If a large number of integrals is to be evaluated, it
  !  is probably better to precompute F(n,x) on a grid, then use series relationship between
  !  different F orders to interpolate. the values.
  !
  !  When compiled in double precision, this routine gives at least 14 significant digits
  !  for all inputs with n<=40, 1e-6<=t<=200. At larger t values, the integral becomes
  !  too small to be worth considering at all.
  !
  function MathBoysFReal(n,t) result(f)
    integer(ik), intent(in) :: n     ! Order of the Boys' function
    real(rk), intent(in)    :: t     ! Argument of the Boys' function
    real(rk)                :: f
    !
    integer(ik), parameter :: sum_terms = 30_ik ! Max. number of terms to include in the explicit sum
    real(rk)               :: loge              ! Desired precision
    real(rk)               :: asse              ! Estimated error in the asymptotic expression
    real(rk)               :: expt              ! Exp(-t)
    integer(ik)            :: safe_n            ! Min. order safe for explict summation
    real(rk)               :: safe_f            ! Boys' function for the safe order
    real(rk)               :: safe_f_term       ! Current term in the Boys' function
    integer(ik)            :: i                 ! 
    !
    !  Try to use the asymptotic formula first
    !
    loge = Log(spacing(1._rk))
    asse = 2*n*max(1._rk,Log(t))-T
    ! write (out,"('= loge = ',g25.16,' asse = ',g25.16)") loge, asse
    if (asse<loge) then
      !
      !  Asymptotic formula is expected to be accurate enough - hurray!
      !
      f = gamma_n_plus_half_real(n)/(2._rk*t**(n+0.5_rk))
      ! write (out,"('= assf = ',g25.16)") f
    else
      !
      !  Asymptotic formula does not work; find gamma order safe for the explicit sum.
      !  The estimate we use here is incredibly generous.
      !
      expt   = Exp(-t)                                ! We'll need this many times ...
      safe_n = int(t * Exp(-loge/sum_terms) + 0.5_rk) ! Actual estimate has -0.5; +0.5 gives safe rounding direction
      ! write (out,"('= expt = ',g25.16,' safe_n = ',i0)") expt, safe_n
      if (safe_n<n) safe_n = n                        ! We should only do downward recursions!
      !
      !  Evaluate F(safe_n,t) using explicit summation. Don't bother to truncate the sum.
      !  We'll trade a bit of accuracy for the dynamic range, and use logarithmic representation
      !
      safe_f_term = expt / (2*safe_n+1)
      safe_f      = safe_f_term
      boys_f_series: do i=1,sum_terms
        safe_f_term = safe_f_term * (2._rk*t) / (2*i+2*safe_n+1)
        safe_f      = safe_f + safe_f_term
      end do boys_f_series
      ! write (out,"('= safe_f = ',g25.16)") safe_f
      !
      !  Evaluate downward recursions until we have the desired order
      !
      f = safe_f
      boys_f_recurse: do i=safe_n-1,n,-1
        f = (2*t*f + expt)/(2*i+1)
      end do boys_f_recurse
      ! write (out,"('= rec_f = ',g25.16)") f
      !
      !  We should be OK now ...
      !
    endif
  end function MathBoysFReal
  !
  !  Gamma(n+1/2). What did you expect?
  !
  function gamma_n_plus_half_real(n) result(g)
    integer(ik), intent(in) :: n
    real(rk)                :: g
    !
    g = sqrtpi * MathDoubleFactorial(2*n-1,g) / 2._rk**n
  end function gamma_n_plus_half_real

  function MathBoysFQuad(n,t) result(f)
    integer(ik), intent(in) :: n     ! Order of the Boys' function
    real(xrk), intent(in)   :: t     ! Argument of the Boys' function
    real(xrk)               :: f
    !
    integer(ik), parameter :: sum_terms = 30_ik ! Max. number of terms to include in the explicit sum
    real(xrk)              :: loge              ! Desired precision
    real(xrk)              :: asse              ! Estimated error in the asymptotic expression
    real(xrk)              :: expt              ! Exp(-t)
    integer(ik)            :: safe_n            ! Min. order safe for explict summation
    real(xrk)              :: safe_f            ! Boys' function for the safe order
    real(xrk)              :: safe_f_term       ! Current term in the Boys' function
    integer(ik)            :: i                 ! 
    !
    !  Try to use the asymptotic formula first
    !
    loge = Log(spacing(1._xrk))
    asse = 2*n*max(1._xrk,Log(t))-T
    ! write (out,"('= loge = ',g25.16,' asse = ',g25.16)") loge, asse
    if (asse<loge) then
      !
      !  Asymptotic formula is expected to be accurate enough - hurray!
      !
      f = gamma_n_plus_half_quad(n)/(2._xrk*t**(n+0.5_xrk))
      ! write (out,"('= assf = ',g25.16)") f
    else
      !
      !  Asymptotic formula does not work; find gamma order safe for the explicit sum.
      !  The estimate we use here is incredibly generous.
      !
      expt   = Exp(-t)                                 ! We'll need this many times ...
      safe_n = int(t * Exp(-loge/sum_terms) + 0.5_xrk) ! Actual estimate has -0.5; +0.5 gives safe rounding direction
      ! write (out,"('= expt = ',g25.16,' safe_n = ',i0)") expt, safe_n
      if (safe_n<n) safe_n = n                         ! We should only do downward recursions!
      !
      !  Evaluate F(safe_n,t) using explicit summation. Don't bother to truncate the sum.
      !  We'll trade a bit of accuracy for the dynamic range, and use logarithmic representation
      !
      safe_f_term = expt / (2*safe_n+1)
      safe_f      = safe_f_term
      boys_f_series: do i=1,sum_terms
        safe_f_term = safe_f_term * (2._xrk*t) / (2*i+2*safe_n+1)
        safe_f      = safe_f + safe_f_term
      end do boys_f_series
      ! write (out,"('= safe_f = ',g25.16)") safe_f
      !
      !  Evaluate downward recursions until we have the desired order
      !
      f = safe_f
      boys_f_recurse: do i=safe_n-1,n,-1
        f = (2*t*f + expt)/(2*i+1)
      end do boys_f_recurse
      ! write (out,"('= rec_f = ',g25.16)") f
      !
      !  We should be OK now ...
      !
    endif
  end function MathBoysFQuad
  !
  !  Gamma(n+1/2). What did you expect?
  !
  function gamma_n_plus_half_quad(n) result(g)
    integer(ik), intent(in) :: n
    real(xrk)               :: g
    !
    g = sqrt(pi_xrk) * MathDoubleFactorial(2*n-1,g) / 2._xrk**n
  end function gamma_n_plus_half_quad
  !
  !  Determinant and inverse of a 3x3 matrix are so common that it makes 
  !  sense to code these explicitly. See lapack.f90 for the general case.
  !
  function MathDet3x3Real(m) result(d)
    real(rk), intent(in) :: m(3,3) ! Matrix for which we need the determinant; 
                                   ! it seems wasteful to call general routines in lapack.f90 for this.
    real(rk)             :: d
    !
    d = m(1,1)*m(2,2)*m(3,3) + m(1,3)*m(2,1)*m(3,2) + m(1,2)*m(2,3)*m(3,1) &
      - m(1,3)*m(2,2)*m(3,1) - m(1,2)*m(2,1)*m(3,3) - m(1,1)*m(2,3)*m(3,2)
  end function MathDet3x3Real
  !
  function MathInv3x3Real(m) result(inv)
    real(rk), intent(in) :: m(3,3)   ! Matrix which we need to invert
    real(rk)             :: inv(3,3) ! The inverse; again, lapack.f90 has a general routine ...
    !
    inv(1,1) = minor(2,2,3,3) ; inv(1,2) = minor(1,3,3,2) ; inv(1,3) = minor(1,2,2,3)
    inv(2,1) = minor(2,3,3,1) ; inv(2,2) = minor(1,1,3,3) ; inv(2,3) = minor(1,3,2,1)
    inv(3,1) = minor(2,1,3,2) ; inv(3,2) = minor(1,2,3,1) ; inv(3,3) = minor(1,1,2,2)
    inv = inv / MathDet3x3(m)
    !
    contains
    real(rk) function minor(i,j,k,l)
      integer(ik) :: i, j, k, l
      minor = m(i,j)*m(k,l) - m(i,l)*m(k,j)
    end function minor
  end function MathInv3x3Real
  !
  function MathDet3x3Complex(m) result(d)
    complex(rk), intent(in) :: m(3,3) ! Matrix for which we need the determinant; 
                                      ! it seems wasteful to call general routines in lapack.f90 for this.
    complex(rk)             :: d
    !
    d = m(1,1)*m(2,2)*m(3,3) + m(1,3)*m(2,1)*m(3,2) + m(1,2)*m(2,3)*m(3,1) &
      - m(1,3)*m(2,2)*m(3,1) - m(1,2)*m(2,1)*m(3,3) - m(1,1)*m(2,3)*m(3,2)
  end function MathDet3x3Complex
  !
  function MathInv3x3Complex(m) result(inv)
    complex(rk), intent(in) :: m(3,3)   ! Matrix which we need to invert
    complex(rk)             :: inv(3,3) ! The inverse; again, lapack.f90 has a general routine ...
    !
    inv(1,1) = minor(2,2,3,3) ; inv(1,2) = minor(1,3,3,2) ; inv(1,3) = minor(1,2,2,3)
    inv(2,1) = minor(2,3,3,1) ; inv(2,2) = minor(1,1,3,3) ; inv(2,3) = minor(1,3,2,1)
    inv(3,1) = minor(2,1,3,2) ; inv(3,2) = minor(1,2,3,1) ; inv(3,3) = minor(1,1,2,2)
    inv = inv / MathDet3x3(m)
    !
    contains
    complex(rk) function minor(i,j,k,l)
      integer(ik) :: i, j, k, l
      minor = m(i,j)*m(k,l) - m(i,l)*m(k,j)
    end function minor
  end function MathInv3x3Complex
  !
  !  Evaluation of Dawson's integral:
  !                       x
  !                       /
  !   F(x) = exp(-x**2) * | Exp(y**2) d y
  !                       /
  !                       0
  !  We use the method of George B. Rybicki, from Computets in Physics, 3, 85 (1989)
  !
  !   F(x) = (pi)**-0.5  Sum Exp(-(z-n*h)**2)/n
  !                     n odd
  !
  !  The error in this approximate integral is estimated to be less than Exp(-(pi/(2*h))**2)
  !
  !  We do not care about the speed too much, so we won't bother with trying to eliminate
  !  exponentiation inside the summation.
  !
  !  For double precision, this implementation produces a result with 14-15 significant
  !  digits.
  !
  function MathDawsonF(z) result(f)
    real(rk), intent(in)   :: z ! Argument to evaluate Dawson's integral for
    real(rk)               :: f 
    !
    real(rk), parameter    :: magic_h = 1.03517084701_rk      ! pi/(2*sqrt(log(10.)))
    real(rk), parameter    :: sqrt_nn = 6._rk                 ! sqrt(number of significant digits desired)
    real(rk), parameter    :: h       = magic_h / sqrt_nn     ! magic_h/sqrt(nn), where nn is the number of significant digits
                                                              ! desired in the integral [36 digits is definitly an overkill!]
    real(rk), parameter    :: magic_n = 1.51742712938_rk/2    ! sqrt(log(10.))/2
    integer(ik), parameter :: max_n   = 2*nint(magic_n * sqrt_nn / h) ! magic_n*sqrt(nn)/h
    real(rk), parameter    :: ser_c1  =   1._rk
    real(rk), parameter    :: ser_c3  =  -2._rk/3._rk
    real(rk), parameter    :: ser_c5  =   4._rk/15._rk
    real(rk), parameter    :: ser_c7  =  -8._rk/105._rk
    real(rk), parameter    :: ser_c9  =  16._rk/945._rk
    real(rk), parameter    :: ser_c11 = -32._rk/10395._rk
    real(rk), parameter    :: ser_c13 =  64._rk/135135._rk
    real(rk), parameter    :: ser_c15 =-128._rk/2027025._rk
    real(rk), parameter    :: ser_c17 = 256._rk/34459425._rk
    real(rk), parameter    :: ser_eps =  0.7_rk * (10._rk**(-(sqrt_nn**2)/17._rk))
    !
    real(rk)    :: z2    ! z**2
    real(rk)    :: zcorr ! z - n0*h
    real(rk)    :: n0    ! Central value of n, giving the largest term in the sum
                         ! This needs to be real to avoid an integer overflow for 
                         ! large arguments
    real(rk)    :: n
    integer(ik) :: dn
    !
    if (abs(z)<=ser_eps) then
      !
      !  Use series expansion
      !
      z2 = z**2
      f  = z*(ser_c1+z2*(ser_c3+z2*(ser_c5+z2*(ser_c7+z2*(ser_c9+z2*(ser_c11+z2*(ser_c13+z2*(ser_c15+z2*ser_c17))))))))
      return
    end if
    !
    !  Use Rybicki's expansion
    !
    f     = 0
    n0    = 2._rk*aint(0.5_rk*z/h)+1.0_rk  ! n0 has to be odd; it is OK if we don't hit the exact maximum
    zcorr = z-n0*h
    ! write (*,*) 'n0 = ',n0,' zcorr= ',zcorr
    if (abs(zcorr)>2*h) then
      stop 'math%MathDawsonF - domain error'
    end if
    sum_terms: do dn=-max_n,max_n,2
      n = n0 + dn
      f = f + exp(-(zcorr-dn*h)**2)/n
    end do sum_terms
    f = f / sqrt(pi)
    !
  end function MathDawsonF
  !
end module math
