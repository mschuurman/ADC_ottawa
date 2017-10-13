!######################################################################
! thetamod: calculation of the MO representation of the projector,
!           Theta, onto the CAP region (also known as the
!           characteristic function)
!######################################################################

module thetamod

  implicit none

  save
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter     :: dp=selected_real_kind(8)

contains
  
!######################################################################

  subroutine theta_mobas(gam,theta_mo)

    use channels
    use iomod
    use parameters
    use timingmod
    use monomial_analytic
    use import_gamess

    implicit none

    integer                               :: k
    real(dp)                              :: tw1,tw2,tc1,tc2
    real(dp), dimension(:,:), allocatable :: theta_mo
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the MO representation of the &
         projector onto the CAP region'
    write(ilog,'(72a,/)') ('-',k=1,72)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(theta_mo(nbas,nbas))
    theta_mo=0.0d0

!----------------------------------------------------------------------    
! Calculate the MO representation of the projector onto the CAP region
! At present, this is only supported for analytically calculated
! monomial-type caps
!----------------------------------------------------------------------
    if (icap.eq.1) then
       call monomial_ana(gam,theta_mo,0,1.0d0)
    else
       errmsg='Currently, flux analysis is only supported for &
            analytic monomial-type CAPs...'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"

    return
    
  end subroutine theta_mobas

!######################################################################
  
end module thetamod
