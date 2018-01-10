!######################################################################
! nto: subroutines for the calculation of natural transition orbitals
!######################################################################

module nto
  
  use constants
  use parameters
  use iomod
  use channels
  
  implicit none
  
contains

!######################################################################
! adc2_nto: gateway subroutine for the calculation of ADC(2) NTOs
!######################################################################
  
  subroutine adc2_nto(ndimf,kpqf)
    
    implicit none

    integer                                   :: ndimf
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',i=1,70)
    write (ilog,'(2x,a)') 'Calculating natural transition orbitals'
    write(ilog,'(70a,/)') ('-',i=1,70)
    
!----------------------------------------------------------------------
! Calculate the NTOs
!----------------------------------------------------------------------
    if (statenumber.eq.0) then
       call adc2_nto_gs(ndimf,kpqf)
    else
       errmsg='WRITE THE EXCITED STATE-TO-EXCITED STATE NTO CODE!'
       call error_control
    endif
    
    return
    
  end subroutine adc2_nto

!######################################################################
! adc2_nto_gs: calculation of ADC(2) NTOs for excitation from the
!              ground state
!######################################################################
  
  subroutine adc2_nto_gs(ndimf,kpqf)

    use diagmod, only: readdavvc
    use density_matrix
    
    implicit none
    
    integer                                   :: ndimf
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    real(d), allocatable                      :: rvec(:,:)
    real(d), allocatable                      :: ener(:)
    real(d), allocatable                      :: trdens(:,:,:)
    character(len=36)                         :: filename

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rvec(ndimf,davstates_f))
    rvec=0.0d0

    allocate(ener(davstates_f))
    ener=0.0d0

    allocate(trdens(nbas,nbas,davstates_f))
    trdens=0.0d0
    
!----------------------------------------------------------------------
! Read the ADC(2) vectors
!----------------------------------------------------------------------
    call readdavvc(davstates_f,ener,rvec,'f',ndimf)

!----------------------------------------------------------------------
! Calculate the ADC(2) ground state-to-excited state transition
! density matrix
!----------------------------------------------------------------------
    call adc2_trden_gs(trdens,ndimf,kpqf,rvec,davstates_f)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rvec)
    deallocate(ener)
    deallocate(trdens)
    
    return
    
  end subroutine adc2_nto_gs
    
!######################################################################
  
end module nto
