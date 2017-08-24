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

    integer, dimension(:,:), allocatable  :: kpq,kpqd,kpqf
    integer                               :: i,ndim,ndims,ndimsf,&
                                             nout,ndimf,ndimd,noutf
    integer*8                             :: noffd,noffdf
    real(d)                               :: e0
    real(d), dimension(:,:), allocatable  :: cap_mo
    type(gam_structure)                   :: gam

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
    call cap_mobas(gam,cap_mo)

!-----------------------------------------------------------------------
! Calculate the IS representation of the CAP operator
!-----------------------------------------------------------------------
    call cap_isbas(cap_mo)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    deallocate(cap_mo)
    
    return
    
  end subroutine adc2_cap

!#######################################################################

  subroutine cap_isbas(cap_mo)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    
    implicit none

    real(d), dimension(nbas,nbas) :: cap_mo
    real(d), dimension(nbas,nbas) :: rho0
    
!----------------------------------------------------------------------
! Calculate the ground state density matrix
!----------------------------------------------------------------------
    call rhogs(rho0)

!----------------------------------------------------------------------
! Calculate the IS representation of the shifted CAP operator
!----------------------------------------------------------------------
!    call get_adc2_dipole_improved_omp(ndim,ndim,kpq,kpq,&
!         nbuf_vv(c),nel_vv(c),filename)
    
    return
    
  end subroutine cap_isbas
    
!#######################################################################
  
end module adc2capmod
