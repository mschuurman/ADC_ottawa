module symmetry
!
!  Very primitive symmetry handling, mostly intended for 
!  generating starting guess for the stationary problem.
!  At best, this module qualifies as a quick hack.
!
!  The idea is to take an initial guess, which has A1
!  symmetry, and apply transformations, which will produce
!  functions of other symmetries. For degenerate irreps,
!  there must be enough transformations to span the irrep.
!
!  The key routines are:
!
!  SMsetSymmetry       - chooses symmetry group
!  SMgetTransformCount - quiries the number of irreps, supported 
!                        for the currently chosen group. Degenerate
!                        irreps enter enough times to span the irrep.
!  SMsetTransform      - selects an irrep, and returns it's symbolic
!                        code.
!  SMtransformField    - performs a transformation from an
!                        A1 field to the choosen irrep. Should be
!                        used with FieldProcess routine.
!
  use accuracy
  implicit none
  private
  public SMsetSymmetry, SMsetTransform, SMgetTransformCount
  public SMtransformField
!
!  Generalized symops
!
  integer(ik), parameter :: Symop_A1  =  1
  integer(ik), parameter :: Symop_CSn = -1
!
  type symop
    integer(ik)       :: code ! Generalzied subroutine (see Symop_* 
                              ! above) responsible for applying symop
    character(len=20) :: name ! Human-readable symop name
    !
    !  Symop-specific parameters from this point down
    !
    real(rk)          :: scalePhi   ! Parameters for the Symop_CSn
    real(rk)          :: scaleTheta ! Parameters for the Symop_CSn
  end type symop

!
  integer(ik), save              :: opsCount = 0 ! Number of operations in the
                                                 ! current transformation set
  type(symop), allocatable, save :: opsCodes(:)  ! List of the operations
  integer(ik), save              :: activeOp = 0 ! Currently selected transformation
!
  real(rk), save    :: invCentre (3)   ! Molecular inversion centre
  real(rk), save    :: mainAxis  (3)   ! Direction of the main axis
  real(rk), save    :: minorAxis1(3)   ! Arbitrary directions, perpendicular
  real(rk), save    :: minorAxis2(3)   ! to the main axis and each other
  integer(ik), save :: mainAxisOrder   ! Order of the main axis
!
  real(rk), save    :: scalePhi        ! Parameters for the Symop_CSn
  real(rk), save    :: scaleTheta      ! Parameters for the Symop_CSn
!
  contains
!
!  Public interface routines
!
!
!  SMsetSymmetry sets symmetry code. Currently supported codes
!  are:
!
!  'None'     - no symmetry at all. Leaves all functions
!               unchanged.
!  'Diatomic' - which is actually D_{\infty h}. This code
!               requires three additional parameters:
!          main_axis_order : Max. order of rotations
!          main_axis       : Orientation of the main symmtry
!                            axis
!          inversion_centre: Position of the inversion centre
!
!  Note that this symmetry code should also be used for non-
!  symmetric diatomic molecules as well - otherwise any states 
!  with a nodal plane between the two atoms will converge very
!  slowly.
!
    subroutine SMsetSymmetry(symm,main_axis_order,main_axis,inversion_centre)
      character(len=*), intent(in)      :: symm                ! Symmetry label
      integer(ik), intent(in), optional :: main_axis_order
      real(rk), intent(in), optional    :: main_axis       (3)
      real(rk), intent(in), optional    :: inversion_centre(3)
      !
      integer(ik) :: iord, iop
      !
      if (allocated(opsCodes)) deallocate(opsCodes)
      !
      select case (symm)
        case default
          write (out,"('SMsetSymmetry: symmetry code ',a,' not recognized')") &
                 trim(symm)
          stop 'SMsetSymmetry - bad symmetry'
        case ('None')
          opsCount    = 1
          allocate (opsCodes(1))
          opsCodes(1)%code = Symop_A1 
          opsCodes(1)%name = 'A1'
        case ('Diatomic','Dinfh')
          if (.not.present(main_axis_order) .or. &
              .not.present(main_axis)       .or. &
              .not.present(inversion_centre) ) then
            write (out,"('SMsetSymmetry(',a,'): required arguments missing')") &
                   trim(symm)
          end if
          !
          invCentre     = inversion_centre
          mainAxisOrder = main_axis_order
          call setupMainAxis(main_axis)
          !
          opsCount      = 4*main_axis_order - 2
          allocate (opsCodes(opsCount))
          !
          !  The first two symops break the pattern a little
          !
          opsCodes(1)%code       = Symop_A1
          opsCodes(1)%name       = 'Sigma G'
          !
          opsCodes(2)%code       = Symop_CSn
          opsCodes(2)%name       = 'Sigma U'
          opsCodes(2)%scalePhi   = 0
          opsCodes(2)%scaleTheta = 1
          !
          opsCodes(3:)%code = Symop_CSn 
          !
          iop = 3
          build_lin: do iord=2,main_axis_order
            if (mod(iord,2)==0) then
              write (opsCodes(iop)%name,"('L',i0,'U(1)')") iord-1
                     opsCodes(iop)%scalePhi   =  (iord-1)
                     opsCodes(iop)%scaleTheta =  0
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'U(2)')") iord-1
                     opsCodes(iop)%scalePhi   = -(iord-1)
                     opsCodes(iop)%scaleTheta =  0
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'G(1)')") iord-1
                     opsCodes(iop)%scalePhi   =  (iord-1)
                     opsCodes(iop)%scaleTheta =  1
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'G(2)')") iord-1
                     opsCodes(iop)%scalePhi   = -(iord-1)
                     opsCodes(iop)%scaleTheta =  1
                     iop = iop + 1
            else
              write (opsCodes(iop)%name,"('L',i0,'G(1)')") iord-1
                     opsCodes(iop)%scalePhi   =  (iord-1)
                     opsCodes(iop)%scaleTheta =  1
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'G(2)')") iord-1
                     opsCodes(iop)%scalePhi   = -(iord-1)
                     opsCodes(iop)%scaleTheta =  1
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'U(1)')") iord-1
                     opsCodes(iop)%scalePhi   =  (iord-1)
                     opsCodes(iop)%scaleTheta =  0
                     iop = iop + 1
              write (opsCodes(iop)%name,"('L',i0,'U(2)')") iord-1
                     opsCodes(iop)%scalePhi   = -(iord-1)
                     opsCodes(iop)%scaleTheta =  0
                     iop = iop + 1
            end if
          end do build_lin
      end select
    end subroutine SMsetSymmetry

    function SMsetTransform(it) result(name)
      integer(ik), intent(in) :: it   ! Ordinal number of the irrep
      character(len=10)       :: name ! name of the irrep

      if (it<=0 .or. it>opsCount) then
        write (out,"('SMsetTransform - request for symop ',i5,' but opsCount = ',i4)") &
               it, opsCount
        stop 'SMsetTransform'
      end if
      activeOp = opsCodes(it)%code
      name     = opsCodes(it)%name
      !
      ! Additional parameters
      !
      select case (activeOp)
        case (Symop_CSn)
          scalePhi   = opsCodes(it)%scalePhi 
          scaleTheta = opsCodes(it)%scaleTheta
      end select
    end function SMsetTransform

    function SMgetTransformCount() result(ni)
      integer(ik) :: ni

      ni = opsCount
    end function SMgetTransformCount

    function SMtransformField(coord,f) result (v)
      real(rk), intent(in)    :: coord(3)
      complex(rk), intent(in) :: f
      complex(rk)             :: v
      !
      real(rk)    :: r, theta, phi, ang
      complex(rk) :: scale
      !
      select case (activeOp)
        case default
          write (out,"('Symop ',i4,' not recognized')") activeOp
          stop 'SMtransformField - bad op'
        case (Symop_A1)
          v = f
        case (Symop_CSn)
          call polarMainAxis(coord,r,theta,phi)
          ang   = scalePhi*phi
          scale = r * cos(scaleTheta*theta)*cmplx(cos(ang),sin(ang),kind=rk)
          v     = f * scale
      end select
    end function SMtransformField

!
!  Private support routines
!
    subroutine setupMainAxis(main)
      real(rk), intent(in) :: main(3)
      !
      integer(ik) :: minor_direction
      real(rk)    :: proj
      !
      mainAxis = main
      mainAxis = mainAxis / sqrt(sum(mainAxis**2))
      !
      ! Choose the direction of the minor axis (needed to fix
      ! the rotational angle) such that it has largest component
      ! along the axis where mainAxis is smallest.
      !
      minor_direction = minloc(abs(mainAxis),dim=1)
      minorAxis1 = 0
      minorAxis1(minor_direction) = 1
      minor_direction = minloc(abs(mainAxis),dim=1,mask=(minorAxis1==0))
      minorAxis2 = 0
      minorAxis2(minor_direction) = 1
      !
      ! Make minor axis orthogonal to the major axis
      !
      proj = dot_product(mainAxis,minorAxis1)
      minorAxis1 = minorAxis1 - proj*mainAxis
      minorAxis1 = minorAxis1 / sqrt(sum(minorAxis1**2))
      !
      proj = dot_product(mainAxis,minorAxis2)
      minorAxis2 = minorAxis2 - proj*mainAxis
      proj = dot_product(minorAxis1,minorAxis2)
      minorAxis2 = minorAxis2 - proj*minorAxis1 
      minorAxis2 = minorAxis2 / sqrt(sum(minorAxis2**2))
    end subroutine setupMainAxis
    !
    !  Calculate polar coordinates of a given point, relative
    !  to the main axis.
    !
    subroutine polarMainAxis(coord,r,theta,phi)
      real(rk), intent(in)  :: coord(3)
      real(rk), intent(out) :: r, theta, phi
      !
      real(rk) :: rel(3)  ! Normalized coordinates relative to the 
                          ! inversion centre
      real(rk) :: pm (3)  ! Projection to the main axis' coordinate
                          ! system
      !
      rel = coord - invCentre
      r   = sqrt(sum(rel**2))
      !
      if (r==0) then
        theta = 0
        phi   = 0
        return
      end if
      !
      rel = rel / r
      pm(3) = dot_product(rel,mainAxis)
      pm(1) = dot_product(rel,minorAxis1)
      pm(2) = dot_product(rel,minorAxis2)
      !
      pm = max(-1.0_rk,min(1.0_rk,pm))
      !
      theta = acos(pm(3))
      phi   = atan2(pm(2),pm(1))
    end subroutine polarMainAxis

end module symmetry
