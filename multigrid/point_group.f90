module point_group
!
!  Handling of point-group symmetry. Some useful externally-visible types and entry points:
!
!   pt_symop         - Rotation-inversion-translation operation
!   pt_group         - Point group. This is an opaque type
!
!   pg_symop         - Constructs a desired coordiate transformation operation
!   pg_symop_inverse - Constructs an inverse of a given operation
!   pg_symop_compose - Composes two operations
!
!   pg_initialize    - Initializes a desired point group and setting
!   pg_destroy       - Deletes point group table
!
  use accuracy
  use math
  implicit none
  private
  public pt_symop, pt_group
  public pg_symop, pg_symop_inverse, pg_symop_compose, pg_symop_print
  public pg_symop_apply, pg_symop_same
  public pg_initialize, pg_destroy
  !
  integer(ik), parameter :: verbose     =    1 ! How chatty do we want to be?
  integer(ik), parameter :: max_gen_ops =    4 ! Maximum number of generators within a group
  integer(ik), parameter :: max_grp_ops = 1000 ! Maximum number of operations within a group
  integer(ik), parameter :: max_name    =   50 ! Maximum length of group or operation name
  !
  interface pg_symop_apply
    module procedure pg_symop_apply_vector
  end interface pg_symop_apply
  !
  !  Coordinate transformation in 3D space
  !  Given spatial coordinate x, the transformed coordinate is:
  !   x' = r . x + t
  !
  type pt_symop
    character(len=max_name) :: name   ! Name of the operation
    real(rk)                :: r(3,3) ! Rotation matrix
    real(rk)                :: t(3)   ! Translation vector
  end type pt_symop
  !
  !  Point group descriptor
  !
  type pt_group
    character(len=max_name)     :: name     ! Name of the group
    integer(ik)                 :: nops     ! Number of symmetry operations
    type(pt_symop), allocatable :: ops(:)   ! Symmetry operations
  end type pt_group
  !
  !  Point group generators
  !  Generators are expected in the standard setting, namely:
  !  1. All operators must leave (0,0,0) invariant
  !  2. The main axis is along Z
  !
  type pt_generators
    character(len=max_name) :: name             ! Name of the group
    integer(ik)             :: nops             ! Number of symmetry operations
    type(pt_symop)          :: ops(max_gen_ops) ! Transformation matrices. The translation vector
                                                ! is expected to be zero.
  end type pt_generators
  !
  !
  !
  contains
  !
  !  Useful routines for handling symmetry operations
  !
  function pg_symop (code,angle,axis,origin,name) result(op)
    character(len=*), intent(in)           :: code      ! Code of the symmetry operation
    real(rk), intent(in), optional         :: angle     ! Rotation angle; must be provided for Cn and Sn operations
    real(rk), intent(in), optional         :: axis  (:) ! Direction of the axis; the default it (0,0,1)
    real(rk), intent(in), optional         :: origin(:) ! An invariant point of the operation; the default is (0,0,0)
    character(len=*), intent(in), optional :: name      ! A more specific name for this operation
    type(pt_symop)                         :: op        ! Return an identity operation
    !
    real(rk)    :: lang       ! Rotation angle
    real(rk)    :: cosa, sina ! Sine and cosine of the rotation angle
    real(rk)    :: cp         ! Coefficient for the "parallel" contribution
    real(rk)    :: lax(3)     ! Local copy of the axis parameter
    real(rk)    :: lor(3)     ! Local copy of the origin parameter
    integer(ik) :: i
    !
    op%name = code
    if (present(name)) op%name = name
    select case (code)
      case default
        write (out,"('point_group%pg_symop: symmetry operation code ',a,' is not recognized.')") trim(code)
      case ('E','T')  ! Identity or translation
        op%r = reshape( (/ 1._rk, 0._rk, 0._rk, &
                           0._rk, 1._rk, 0._rk, &
                           0._rk, 0._rk, 1._rk /), (/ 3, 3 /) )
      case ('I')  ! Inversion
        op%r = reshape( (/-1._rk, 0._rk, 0._rk, &
                           0._rk,-1._rk, 0._rk, &
                           0._rk, 0._rk,-1._rk /), (/ 3, 3 /) )
      case ('SG') ! Reflection
        call fill_axis
        do i=1,3
          op%r(:,i) = -2._rk*lax(:)*lax(i)
          op%r(i,i) = op%r(i,i) + 1._rk
        end do
      case ('Cn','Sn') ! Rotation and improper rotation
        call fill_axis
        call fill_angle
        !
        !  The expression is (for zero origin):
        !    r*cosa + n (n.r) [+/-1 - cosa] + (n x r) sina
        !
        cp = 1._rk - cosa
        if (code=='Sn') cp = -1._rk - cosa
        op%r(1,1) = cosa + cp*lax(1)*lax(1)
        op%r(2,1) =        cp*lax(2)*lax(1) + sina*lax(3)
        op%r(3,1) =        cp*lax(3)*lax(1) - sina*lax(2)
        op%r(1,2) =        cp*lax(1)*lax(2) - sina*lax(3)
        op%r(2,2) = cosa + cp*lax(2)*lax(2)
        op%r(3,2) =        cp*lax(3)*lax(2) + sina*lax(1)
        op%r(1,3) =        cp*lax(1)*lax(3) + sina*lax(2)
        op%r(2,3) =        cp*lax(2)*lax(3) - sina*lax(1)
        op%r(3,3) = cosa + cp*lax(3)*lax(3)
    end select
    !
    call fill_origin
    if (code=='T') then
      !
      !  Translation requires special handling
      !
      op%t = lor
    else
      !
      !  All remaining operations are the same as far as origin is concerned.
      !
      op%t = lor - matmul(op%r,lor)
    end if
    !
    if (verbose>=2) then
      call pg_symop_print(op,'pg_symop: Constructed a new operation:')
    end if
    !
    contains
    subroutine fill_angle
      if (.not.present(angle)) stop 'point_group%pg_symop - required angle parameter is missing'
      lang = angle
      cosa = cos(lang)
      sina = sin(lang)
    end subroutine fill_angle
    subroutine fill_axis
      lax = (/ 0._rk, 0._rk, 1._rk /)
      if (.not.present(axis)) return
      if (size(axis)/=3) stop 'point_group%pg_symop - parameter axis must have size = 3'
      lax = axis
      lax = lax / sqrt(sum(lax**2))
    end subroutine fill_axis
    subroutine fill_origin
      lor = (/ 0._rk, 0._rk, 0._rk /)
      if (.not.present(origin)) return
      if (size(origin)/=3) stop 'point_group%pg_symop - parameter origin must have size = 3'
      lor = origin
    end subroutine fill_origin
  end function pg_symop
  !
  function pg_symop_inverse (a) result(op)
    type(pt_symop), intent(in) :: a   ! Symmetry operation to invert
    type(pt_symop)             :: op  ! op * a must be an identity operation.
                                      ! a * op is necessarily also an identity operation
    !
    op%name = ' '
    if (a%name/=' ') op%name = '(' // trim(a%name) // ')^-1'
    op%r = MathInv3x3(a%r)
    op%t = -matmul(op%r,a%t)
    !
    if (verbose>=2) then
      call pg_symop_print(a, 'pg_symop_inverse: inverting transformation:')
      call pg_symop_print(op,'pg_symop_inverse: the inverse is:')
    end if
  end function pg_symop_inverse
  !
  function pg_symop_compose (a,b) result(op)
    type(pt_symop), intent(in) :: a, b ! Symmetry operations to compose
    type(pt_symop)             :: op   ! op = a * b; note that a and b do not have to commute
    !
    if(a%name==' ') then
      op%name = b%name
    else
      op%name = trim(a%name) // ' ' // trim(b%name)
    end if
    op%r = matmul(a%r,b%r)
    op%t = matmul(a%r,b%t) + a%t
    !
    if (verbose>=2) then
      call pg_symop_print(a, 'pg_symop_compose: multiplying transformation:')
      call pg_symop_print(b, 'pg_symop_compose: on the right by:')
      call pg_symop_print(op,'pg_symop_compose: giving:')
    end if
  end function pg_symop_compose
  !
  function pg_symop_same(a,b) result(s)
    type(pt_symop), intent(in) :: a, b ! Operators to compare
    logical                    :: s    ! Are they the same?
    !
    real(rk) :: eps 
    !
    eps = 0._rk
    eps = max(eps,maxval(abs(a%r)))
    eps = max(eps,maxval(abs(a%t)))
    eps = max(eps,maxval(abs(b%r)))
    eps = max(eps,maxval(abs(b%t)))
    eps = 100._rk * spacing(eps)
    !
    s = all(abs(a%r-b%r)<=eps) .and. all(abs(a%t-b%t)<=eps)
  end function pg_symop_same
  !
  subroutine pg_symop_print(op,label)
    type(pt_symop), intent(in)             :: op  ! Operation to report
    character(len=*), intent(in), optional :: label
    !
    if (present(label)) write (out,"(1x,a,':')",advance='no') trim(label)
    write (out,"(1x,a)") trim(op%name)
    write (out,"(4x,a,1x,3(1x,a12  ),2x,a12  )") '   ', ' Rotation, X ', ' Rotation, Y ', ' Rotation, Z ', ' Translation '
    write (out,"(4x,a,1x,3(1x,f12.9),2x,f12.5)") 'X: ', op%r(1,:), op%t(1)
    write (out,"(4x,a,1x,3(1x,f12.9),2x,f12.5)") 'Y: ', op%r(2,:), op%t(2)
    write (out,"(4x,a,1x,3(1x,f12.9),2x,f12.5)") 'Z: ', op%r(3,:), op%t(3)
    write (out,"(4x,'det = ',f18.14)") MathDet3x3(op%r)
  end subroutine pg_symop_print
  !
  function pg_symop_apply_vector(op,vec) result(vp)
    type(pt_symop), intent(in) :: op      ! Operation to apply
    real(rk), intent(in)       :: vec(3)  ! Vector to transform
    real(rk)                   :: vp (3)  ! The result
    !
    vp = matmul(op%r,vec) + op%t
  end function pg_symop_apply_vector
  !
  !  External interfaces
  !
  subroutine pg_initialize(name,pg,setting)
    character(len=*), intent(in)         :: name     ! Name of the point group
    type(pt_group), intent(out)          :: pg       ! Point group descriptor
    type(pt_symop), intent(in), optional :: setting  ! Possible non-default setting; acting with
                                                     ! the setting symop on a coordinate will bring it
                                                     ! into the standard setting.
    !
    type(pt_generators) :: gen              ! Generators for the group
    type(pt_symop)      :: set, invset      ! Setting for the group and the inverse
    integer(ik)         :: nops             ! Number of operators
    type(pt_symop)      :: ops(max_grp_ops) ! Current set of symmetry operators
    type(pt_symop)      :: test
    logical             :: grown            ! Did the group grow in the last pass?
    integer(ik)         :: alloc
    integer(ik)         :: i, j, k
    !
    !  Fill group generators in the standard setting
    !
    call fill_generators(name,gen)
    !
    !  Fill out the group, starting with the generators
    !
    nops = gen%nops
    ops(:nops) = gen%ops(:nops)
    grown = .true.
    expand_group: do while(grown)
      grown = .false.
      scan_operators_1: do i=1,nops
        scan_operators_2: do j=1,nops
          test = pg_symop_compose(ops(i),ops(j))
          test_operators: do k=1,nops
            if (pg_symop_same(test,ops(k))) cycle scan_operators_2
          end do test_operators
          !
          !  Operator is new; add it to the table
          !
          nops = nops + 1
          if (nops>max_grp_ops) stop 'point_group%pg_initialize - blown operator table'
          ops(nops) = test
          grown = .true.
        end do scan_operators_2
      end do scan_operators_1
    end do expand_group
    !
    !  Transform to the desired setting
    !
    set = pg_symop('E',name=' ')
    if (present(setting)) set = setting
    invset = pg_symop_inverse(set)
    convert_operators: do i=1,nops
      ops(i) = pg_symop_compose(invset,pg_symop_compose(ops(i),set))
    end do convert_operators
    !
    !  Copy final symops to the group table
    !
    allocate (pg%ops(nops),stat=alloc)
    if (alloc/=0) stop 'point_group%pg_initialize - allocation failed (1)'
    pg%nops = nops
    pg%ops  = ops(:nops)
    pg%name = gen%name // ' ' // trim(set%name)
    !
    if (verbose>=1) then
      write (out,"('Found ',i0,' operators in group ',a)") nops, trim(name)
      print_ops: do i=1,pg%nops
        call pg_symop_print(pg%ops(i))
      end do print_ops
    end if
  end subroutine pg_initialize
  !
  subroutine pg_destroy
  end subroutine pg_destroy
  !
  !  Tables of group generators
  !
  subroutine fill_generators(name,gen)
    character(len=*), intent(in)     :: name ! Name of the point group
    type(pt_generators), intent(out) :: gen  ! Point group generators
    !
    real(rk)       :: sih ! Icosahedral angle
    !
    gen%ops(1) = pg_symop('E')
    select case (name)
      case default
        write (out,"('point_group%fill_generators: Point group ',a,' is not recognized')") trim(name)
        stop 'point_group%fill_generators - bad point group'
      case ('C1') ! No operations beyond implied 'E'
        gen%nops = 1
      case ('Cs') 
        gen%nops = 2
        gen%ops(2) = pg_symop('SG')
      case ('Ci') 
        gen%nops = 2
        gen%ops(2) = pg_symop('I')
      case ('C2') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2")
      case ('C3') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
      case ('C4') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
      case ('C5') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
      case ('C6') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
      case ('S4') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Sn',angle=twopi/4,name="S4")
      case ('S6') 
        gen%nops = 2
        gen%ops(2) = pg_symop('Sn',angle=twopi/6,name="S6")
      case ('D2') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2z")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('D3') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('D4') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('D5') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('D6') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('D2h') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2z")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGh",axis=(/0.0_rk,0.0_rk,1.0_rk/))
      case ('D3h') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGh",axis=(/0.0_rk,0.0_rk,1.0_rk/))
      case ('D4h') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGh",axis=(/0.0_rk,0.0_rk,1.0_rk/))
      case ('D5h') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGh",axis=(/0.0_rk,0.0_rk,1.0_rk/))
      case ('D6h') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGh",axis=(/0.0_rk,0.0_rk,1.0_rk/))
      case ('D2d') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2z")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGd",axis=(/1.0_rk,1.0_rk,0.0_rk/))
      case ('D3d') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('I')
      case ('D4d') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGd",axis=(/cos(pi/8),sin(pi/8),0.0_rk/))
      case ('D5d') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('I')
      case ('D6d') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
        gen%ops(3) = pg_symop('Cn',angle=twopi/2,name="C2x",axis=(/1.0_rk,0.0_rk,0.0_rk/))
        gen%ops(4) = pg_symop('SG',name="SGd",axis=(/1.0_rk,1.0_rk,0.0_rk/))
      case ('C2v') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('C3v') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('C4v') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('C5v') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('C6v') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
        gen%ops(3) = pg_symop('SG',name="SGv",axis=(/1.0_rk,0.0_rk,0.0_rk/))
      case ('C2h') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2")
        gen%ops(3) = pg_symop('SG',name="SGh")
      case ('C3h') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/3,name="C3")
        gen%ops(3) = pg_symop('SG',name="SGh")
      case ('C4h') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('SG',name="SGh")
      case ('C5h') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5")
        gen%ops(3) = pg_symop('SG',name="SGh")
      case ('C6h') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/6,name="C6")
        gen%ops(3) = pg_symop('SG',name="SGh")
      case ('T') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2z")
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3", axis=(/1.0_rk,1.0_rk,1.0_rk/))
      case ('Th') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/2,name="C2z")
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3", axis=(/1.0_rk,1.0_rk,1.0_rk/))
        gen%ops(4) = pg_symop('I')
      case ('Td') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Sn',angle=twopi/4,name="S4")
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3", axis=(/1.0_rk,1.0_rk,1.0_rk/))
      case ('O') 
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3", axis=(/1.0_rk,1.0_rk,1.0_rk/))
      case ('Oh') 
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/4,name="C4")
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3", axis=(/1.0_rk,1.0_rk,1.0_rk/))
        gen%ops(4) = pg_symop('I')
      case ('I') 
        sih = atan(sqrt(1.25_rk)+0.5_rk)
        gen%nops = 3
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5",axis=(/cos(sih),0._rk,sin(sih)/))
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3",axis=(/1.0_rk,1.0_rk,1.0_rk/))
      case ('Ih') 
        sih = atan(sqrt(1.25_rk)+0.5_rk)
        gen%nops = 4
        gen%ops(2) = pg_symop('Cn',angle=twopi/5,name="C5",axis=(/cos(sih),0._rk,sin(sih)/))
        gen%ops(3) = pg_symop('Cn',angle=twopi/3,name="C3",axis=(/1.0_rk,1.0_rk,1.0_rk/))
        gen%ops(4) = pg_symop('I')
    end select
    !
  end subroutine fill_generators
  !
end module point_group
