!######################################################################
! capmod: routines for the calculation of the MO representation of
!         the CAP operator. Employs numerical integration performed
!         using Becke's partitioning scheme (JCP 88, 2547 (1988)).
!         Makes use of the NUMGRID libraries of Radovan Bast.
!######################################################################

module capmod

  use constants
  use numgrid
  use, intrinsic :: iso_c_binding, only: c_ptr
  
  implicit none  
  
  save

  ! NUMGRID arrays and variables
  integer              :: min_num_angular_points
  integer              :: max_num_angular_points
  integer              :: num_points
  integer              :: num_centers
  integer              :: num_outer_centers
  integer              :: num_shells
  integer              :: num_primitives
  integer, allocatable :: center_elements(:)
  integer, allocatable :: shell_centers(:)
  integer, allocatable :: shell_l_quantum_numbers(:)
  integer, allocatable :: shell_num_primitives(:)
  integer, allocatable :: outer_center_elements(:)
  real(d), allocatable :: outer_center_coordinates(:)
  real(d), allocatable :: primitive_exponents(:)
  real(d)              :: radial_precision
  real(d), allocatable :: center_coordinates(:)
  real(d), pointer     :: grid(:)
  type(c_ptr)          :: context
  
contains
  
!######################################################################

  subroutine cap_mobas(gam)

    use channels
    use constants
    use parameters
    use import_gamess
    
    implicit none

    type(gam_structure) :: gam
    
!----------------------------------------------------------------------
! Set up the integration grid
!----------------------------------------------------------------------
    call initialise_intgrid(gam)

!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    call calc_mocapmat(gam)
    
!----------------------------------------------------------------------
! Deallocate arrays and destroy the grid context
!----------------------------------------------------------------------
    call finalise_intgrid
    
    return
    
  end subroutine cap_mobas

!######################################################################

  subroutine initialise_intgrid(gam)

    use channels
    use constants
    use parameters
    use import_gamess
    
    implicit none

    integer            :: i,j,k,l
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Grid parameters
!----------------------------------------------------------------------
    ! Grid parameters
    radial_precision=1.0d-12
    min_num_angular_points=86
    max_num_angular_points=302

!----------------------------------------------------------------------
! Atomic coordinates (in Bohr). Note that the gam_structure derived
! type holds the coordinates in Angstrom.
!----------------------------------------------------------------------
    ! Number of atoms
    num_centers=gam%natoms

    ! Atom coordinates
    allocate(center_coordinates(num_centers*3))

    do i=1,gam%natoms
       k=i*3-3
       do j=1,3
          k=k+1
          center_coordinates(k)=gam%atoms(i)%xyz(j)*0.529177249d0
       enddo
    enddo

!----------------------------------------------------------------------
! Atom types
!----------------------------------------------------------------------
    allocate(center_elements(num_centers))
    do i=1,gam%natoms
       center_elements(i)=int(gam%atoms(i)%znuc)
    enddo

!----------------------------------------------------------------------
! Basis information
!----------------------------------------------------------------------
    num_outer_centers=0

    ! Number of shells
    num_shells=0
    do i=1,gam%natoms
       num_shells=num_shells+gam%atoms(i)%nshell
    enddo

    ! Shell centers
    allocate(shell_centers(num_shells))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_centers(k)=i
       enddo
    enddo

    ! Angular momentum quantum numbers for each shell
    allocate(shell_l_quantum_numbers(num_shells))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_l_quantum_numbers(k)=gam%atoms(i)%sh_l(j)
       enddo
    enddo

    ! Number of primitives for each shell
    allocate(shell_num_primitives(num_shells))
    k=0
    do i=1,gam%natoms
       do j=2,gam%atoms(i)%nshell+1
          k=k+1
          shell_num_primitives(k)=gam%atoms(i)%sh_p(j)-gam%atoms(i)%sh_p(j-1)
       enddo
    enddo

    ! Total number of primitives
    num_primitives=sum(shell_num_primitives)

    ! Primitive exponents
    allocate(primitive_exponents(num_primitives))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
             k=k+1
             primitive_exponents(k)=gam%atoms(i)%p_zet(l)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Create the grid
!----------------------------------------------------------------------
    context=numgrid_new_context()

    call numgrid_generate_grid(context,&
         radial_precision,&
         min_num_angular_points,&
         max_num_angular_points,&
         num_centers,&
         center_coordinates,&
         center_elements,&
         num_outer_centers,&
         outer_center_coordinates,&
         outer_center_elements,&
         num_shells,&
         shell_centers,&
         shell_l_quantum_numbers,&
         shell_num_primitives,&
         primitive_exponents)

!----------------------------------------------------------------------
! Number of grid points
!----------------------------------------------------------------------
    num_points=numgrid_get_num_points(context)

    return
    
  end subroutine initialise_intgrid

!######################################################################

  subroutine calc_mocapmat(gam)

    use channels
    use constants
    use parameters
    use import_gamess
    
    implicit none

    type(gam_structure) :: gam

    print*,"FINISH WRITING THE CAP CONSTRUCTION CODE!"
    STOP
    
    return
    
  end subroutine calc_mocapmat
    
!######################################################################
  
  subroutine finalise_intgrid

    implicit none
    
    deallocate(center_coordinates)
    deallocate(center_elements)
    deallocate(shell_centers)
    deallocate(shell_l_quantum_numbers)
    deallocate(shell_num_primitives)
    deallocate(primitive_exponents)
    
    call numgrid_free_context(context)
    
  end subroutine finalise_intgrid
    
!######################################################################
  
end module capmod
