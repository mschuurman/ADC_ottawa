!######################################################################
! Routines for the automatic determination of a CAP box based on the
! analysis of the initial state density
!######################################################################

module autocapbox
  
  implicit none

  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter    :: dp=selected_real_kind(8)
  
  ! Initial state density matrix
  real(dp), allocatable :: rho(:,:)
  
contains

!######################################################################

  subroutine autobox(gam,cstrt)

    use channels
    use parameters
    use iomod
    use import_gamess
    
    implicit none
    
    real(dp), dimension(3) :: cstrt
    type(gam_structure)    :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    call initialise
    
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
    call density_analysis(gam)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    call finalise
    
    return
    
  end subroutine autobox

!######################################################################

  subroutine initialise

    use parameters
    
    implicit none

    ! Initial state density matrix
    allocate(rho(nbas,nbas))
    rho=0.0d0
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine finalise

    implicit none

    deallocate(rho)
    
    return
    
  end subroutine finalise
    
!######################################################################

  subroutine istate_density_matrix

    use parameters
    use mp2
    
    implicit none

    integer  :: i
    real(dp) :: trace
    
!----------------------------------------------------------------------
! ADC(1) ground state density matrix
!----------------------------------------------------------------------
    if (method.eq.1) then
       rho=0.0d0
       do i=1,nocc
          rho(i,i)=2.0d0
       enddo
    endif

!----------------------------------------------------------------------
! ADC(2) ground state density matrix
!----------------------------------------------------------------------
    if (method.eq.2.or.method.eq.3) then
       call rho_mp2(rho)
    endif

    return
    
  end subroutine istate_density_matrix
    
!######################################################################

  subroutine density_analysis(gam)

    use channels
    use parameters
    use iomod
    use density
    use import_gamess

    implicit none

    integer                :: i,k
    integer, parameter     :: npnt=201
    real(dp), parameter    :: dx=0.1d0
    real(dp), dimension(3) :: r
    real(dp)               :: dens
    type(gam_structure)    :: gam


    
    ! REMEMBER TO SCAN IN THE NEGATIVE AS WELL AS
    ! POSITIVE DIRECTIONS

    
    
    ! Loop over the x, y and z directions
    do i=1,3

       print*,
       
       ! Loop over grid points
       do k=1,npnt

          ! Current coordinate (in Bohr)
          r=0.0d0
          r(i)=0.0d0+(k-1)*dx

          ! Initial state density at the current coordinate
          dens=density_value(gam,rho,r)

          print*,k,dens
          
       enddo
       
    enddo

    STOP
    
    return
    
  end subroutine density_analysis
    
!######################################################################
  
end module autocapbox
