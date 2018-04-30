!######################################################################
! Routines for the automatic determination of a CAP box based on the
! analysis of the initial state density
!######################################################################

module autocapbox

  use constants
  
  implicit none

  ! Initial state density matrix
  real(dp), allocatable :: rho(:,:)
  
contains

!######################################################################

  subroutine autobox(gam,cstrt,gcent)

    use channels
    use parameters
    use iomod
    use import_gamess
    
    implicit none
    
    real(dp), dimension(3) :: cstrt,gcent
    type(gam_structure)    :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Initial state density matrix
    allocate(rho(nbas,nbas))
    rho=0.0d0
    
!----------------------------------------------------------------------
! Currently, we can only calculate the ground state density matrix,
! so exit here if the initial state is an excited state
!----------------------------------------------------------------------
    if (statenumber.ne.0) then
       errmsg='The density-based determination of the CAP box is &
            only supported for the ground state'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate the initial state density matrix
!----------------------------------------------------------------------
    call istate_density_matrix

!----------------------------------------------------------------------
! Determine the CAP box from an analysis of the initial state density
!----------------------------------------------------------------------
    call density_analysis(gam,cstrt,gcent)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rho)
    
    return
    
  end subroutine autobox

!######################################################################

  subroutine istate_density_matrix

    use parameters
    use mp2
    
    implicit none

    integer  :: i
    real(dp) :: trace
    
!----------------------------------------------------------------------
! ADC(1) (HF) ground state density matrix
!----------------------------------------------------------------------
    if (method.eq.1) then
       rho=0.0d0
       do i=1,nocc
          rho(i,i)=2.0d0
       enddo
    endif

!----------------------------------------------------------------------
! ADC(2) (MP2) ground state density matrix
!----------------------------------------------------------------------
    if (method.eq.2.or.method.eq.3) then
       call rho_mp2(rho)
    endif

    return
    
  end subroutine istate_density_matrix
    
!######################################################################

  subroutine density_analysis(gam,cstrt,gcent)

    use channels
    use parameters
    use iomod
    use electron_density
    use import_gamess

    implicit none

    integer                  :: i,k,dir,unit
    real(dp), parameter      :: dx=0.1d0
    real(dp), dimension(3)   :: cstrt,gcent,r
    real(dp)                 :: dens
    real(dp), dimension(3,2) :: rc
    type(gam_structure)      :: gam

!----------------------------------------------------------------------
! For each of the x-, y-, and z-directions, determine the distance at
! which the initial state density drops below the user-set threshold.
! We take these points to define the CAP box.
!----------------------------------------------------------------------
    ! Loop over negative and positive displacements
    do dir=1,2
    
       ! Loop over the x, y and z directions
       do i=1,3

          ! Loop over points until the density drops
          ! below threshold
          k=0
          do
             
             k=k+1

             ! Current coordinate (in Bohr)
             r=0.0d0
             r(i)=gcent(i)+(k-1)*dx*(-1)**dir
             
             ! Initial state density at the current coordinate
             dens=density_value(gam,rho,r)
             
             ! Exit the loop if the density has dropped below
             ! threshold
             if (dens.lt.densthrsh) exit

          enddo
          
          ! Coordinate value (relative to the geometric centre) at
          ! which the density dropped below threshold
          rc(i,dir)=(k-1)*dx*(-1)**dir

       enddo

    enddo

    ! Set the CAP box parameters
    do i=1,3
       cstrt(i)=maxval(abs(rc(i,:)))
    enddo

    return
    
  end subroutine density_analysis

!######################################################################
  
end module autocapbox
