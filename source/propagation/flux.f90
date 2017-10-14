!######################################################################
! fluxmod: routines for the calculation of flux expectation values
!######################################################################

module fluxmod

  implicit none

contains

!######################################################################
! adc1_flux_cap: calculation of the flux expectation value for a
!                TD-ADC(1) calculation using a CAP operator
!                analysis
!######################################################################
  
  subroutine adc1_flux_cap(matdim,psi,dtpsi,flux)

    use constants
    use parameters
    
    implicit none

    integer                               :: matdim
    real(d)                               :: flux,val1,val2
    complex(d), dimension(matdim)         :: psi,dtpsi
    complex(d), dimension(:), allocatable :: oppsi

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(oppsi(matdim))
    oppsi=czero
    
!----------------------------------------------------------------------
! (1) 2 Re < d Psi/dt | Theta | Psi >
!----------------------------------------------------------------------
    oppsi=czero

    ! (a) IS-IS block
    !
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +matmul(thetaij,psi(1:matdim-1)) &
         +theta00*psi(1:matdim-1)

    ! (b) Ground state-ground state element
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim)+theta00*psi(matdim)
    endif

    ! (c) Ground state-IS block
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim) &
            +dot_product(theta0j(1:matdim-1),psi(1:matdim-1))
       oppsi(1:matdim-1)=oppsi(1:matdim-1) &
            +theta0j(1:matdim-1)*psi(matdim)
    endif

    val1=2.0d0*real(dot_product(dtpsi,oppsi))

!----------------------------------------------------------------------
! 2 < Psi | W | Psi >
!----------------------------------------------------------------------
    oppsi=czero

    ! (a) IS-IS block
    !
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +matmul(wij,psi(1:matdim-1)) &
         +w00*psi(1:matdim-1)

    ! (b) Ground state-ground state element
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim)+w00*psi(matdim)
    endif

    ! (c) Ground state-IS block
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim) &
            +dot_product(w0j(1:matdim-1),psi(1:matdim-1))
       oppsi(1:matdim-1)=oppsi(1:matdim-1) &
            +w0j(1:matdim-1)*psi(matdim)
    endif

    val2=2.0d0*real(dot_product(psi,oppsi))
    
!----------------------------------------------------------------------
! Flux expectation value
!----------------------------------------------------------------------
    flux=val1+val2
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(oppsi)
    
    return
    
  end subroutine adc1_flux_cap

!######################################################################
! adc2_flux_cap: calculation of the flux expectation value for a
!                TD-ADC(2) calculation using a CAP operator
!                analysis
!######################################################################
  
  subroutine adc2_flux_cap(matdim,psi,dtpsi,flux)

    use constants
    use parameters
    use tdsemod, only: opxvec_ext
    
    implicit none

    integer                               :: matdim
    real(d)                               :: flux,val1,val2
    complex(d), dimension(matdim)         :: psi,dtpsi
    complex(d), dimension(:), allocatable :: oppsi
    character(len=70)                     :: filename
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(oppsi(matdim))
    oppsi=czero
    
!----------------------------------------------------------------------
! (1) 2 Re < d Psi/dt | Theta | Psi >
!----------------------------------------------------------------------
    oppsi=czero

    ! (a) IS-IS block
    !
    ! Calculate W*v1
    filename='SCRATCH/theta'
    call opxvec_ext(matdim-1,psi(1:matdim-1),&
         oppsi(1:matdim-1),filename,nbuf_theta)
    
    oppsi(1:matdim-1)=oppsi(1:matdim-1)+theta00*oppsi(1:matdim-1)

    ! (b) Ground state-ground state element
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim)+theta00*psi(matdim)
    endif

    ! (c) Ground state-IS block
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim) &
            +dot_product(theta0j(1:matdim-1),psi(1:matdim-1))
       oppsi(1:matdim-1)=oppsi(1:matdim-1) &
            +theta0j(1:matdim-1)*psi(matdim)
    endif

    val1=2.0d0*real(dot_product(dtpsi,oppsi))    

!----------------------------------------------------------------------
! 2 < Psi | W | Psi >
!----------------------------------------------------------------------
    oppsi=czero

    ! (a) IS-IS block
    !
    filename='SCRATCH/cap'
    call opxvec_ext(matdim-1,psi(1:matdim-1),&
         oppsi(1:matdim-1),filename,nbuf_cap)
    
    oppsi(1:matdim-1)=oppsi(1:matdim-1)+w00*oppsi(1:matdim-1)

    ! (b) Ground state-ground state element
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim)+w00*psi(matdim)
    endif

    ! (c) Ground state-IS block
    !
    if (.not.lprojcap.or.statenumber.gt.0) then
       oppsi(matdim)=oppsi(matdim) &
            +dot_product(w0j(1:matdim-1),psi(1:matdim-1))
       oppsi(1:matdim-1)=oppsi(1:matdim-1) &
            +w0j(1:matdim-1)*psi(matdim)
    endif

    val2=2.0d0*real(dot_product(psi,oppsi))
    
!----------------------------------------------------------------------
! Flux expectation value
!----------------------------------------------------------------------
    flux=val1+val2
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(oppsi)
    
    return
    
  end subroutine adc2_flux_cap
    
!######################################################################
  
end module fluxmod
