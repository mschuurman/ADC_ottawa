!#######################################################################
! Calculation of ionisation probabilities using the CAP-TD-ADC(2)
! method.
!#######################################################################

module adc2capmod

  use channels

contains

!#######################################################################

  subroutine adc2_cap(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    use capmod
    
    implicit none

    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimsf,&
                                            nout,ndimf,ndimd,noutf
    integer*8                            :: noffd,noffdf
    real(d)                              :: e0
    type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!-----------------------------------------------------------------------
    call mp2_master(e0)

!-----------------------------------------------------------------------
! Determine the 1h1p and 2h2p subspaces
!-----------------------------------------------------------------------
    call get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf,&
         ndims,ndimsf)

!-----------------------------------------------------------------------
! Set MO representation of the dipole operator
!-----------------------------------------------------------------------
    call set_dpl

!-----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!-----------------------------------------------------------------------
    call cap_mobas(gam)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    
    return
    
  end subroutine adc2_cap
    
!#######################################################################
  
end module adc2capmod
