!######################################################################
! nto: subroutines for the calculation of natural transition orbitals
!######################################################################

module nto
  
  use constants
  use parameters
  use iomod
  use channels
  
  implicit none

  save

  private :: dp
    
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter :: dp=selected_real_kind(8)
  
contains

!######################################################################
! adc2_nto: gateway subroutine for the calculation of ADC(2) NTOs
!######################################################################
  
  subroutine adc2_nto(gam,ndimf,kpqf)

    use gamess_internal
    
    implicit none

    integer                                   :: ndimf
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i
    type(gam_structure)                       :: gam
    
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
       call adc2_nto_gs(gam,ndimf,kpqf)
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
  
  subroutine adc2_nto_gs(gam,ndimf,kpqf)

    use gamess_internal
    use diagmod, only: readdavvc
    use density_matrix
    
    implicit none
    
    integer                                   :: ndimf
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i,j
    integer                                   :: nao,lwork,ierr
    real(dp), allocatable                     :: rvec(:,:)
    real(dp), allocatable                     :: ener(:)
    real(dp), allocatable                     :: trdens(:,:,:)
    real(dp), allocatable                     :: tmp(:,:)
    real(dp), allocatable                     :: sigma(:,:)
    real(dp), allocatable                     :: VT(:,:,:)
    real(dp), allocatable                     :: U(:,:,:)
    real(dp), allocatable                     :: work(:)
    real(dp), allocatable                     :: nto_p(:,:,:)
    real(dp), allocatable                     :: nto_h(:,:,:)
    character(len=36)                         :: filename
    type(gam_structure)                       :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rvec(ndimf,davstates_f))
    rvec=0.0d0

    allocate(ener(davstates_f))
    ener=0.0d0

    allocate(trdens(nbas,nbas,davstates_f))
    trdens=0.0d0

    allocate(tmp(nbas,nbas))
    tmp=0.0d0

    allocate(sigma(nbas,davstates_f))
    sigma=0.0d0

    allocate(VT(nbas,nbas,davstates_f))
    VT=0.0d0

    allocate(U(nbas,nbas,davstates_f))
    U=0.0d0

    lwork=5*nbas
    allocate(work(lwork))
    work=0.0d0

    nao=gam%nbasis
    allocate(nto_p(nao,nbas,davstates_f))
    allocate(nto_h(nao,nbas,davstates_f))
    nto_p=0.0d0
    nto_h=0.0d0
    
!----------------------------------------------------------------------
! Read the ADC(2) vectors
!----------------------------------------------------------------------
    call readdavvc(davstates_f,ener,rvec,'f',ndimf)

!----------------------------------------------------------------------
! Calculate the ADC(2) ground state-to-excited state transition
! density matrices
!----------------------------------------------------------------------
    call adc2_trden_gs(trdens,ndimf,kpqf,rvec,davstates_f)

!----------------------------------------------------------------------
! Calculate the SVDs of the ADC(2) ground state-to-excited state
! transition density matrices
!----------------------------------------------------------------------
    do i=1,davstates_f

       tmp=trdens(:,:,i)
       
       call dgesvd('A','A',nbas,nbas,tmp,nbas,sigma(:,i),U(:,:,i),&
            nbas,VT(:,:,i),nbas,work,lwork,ierr)

       if (ierr.ne.0) then
          errmsg='SVD of the transition density matrices failed in &
               subroutine adc2_nto_gs'
       endif
       
    enddo

!----------------------------------------------------------------------
! Calculate the NTOs
!----------------------------------------------------------------------
    do i=1,davstates_f

       print*,i,sigma(1,i)**2,sigma(nbas,i)**2
       
       ! Particle NTOs in terms of the AOs
              
       ! Hole NTOs in terms of the AOs
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rvec)
    deallocate(ener)
    deallocate(trdens)
    deallocate(tmp)
    deallocate(sigma)
    deallocate(VT)
    deallocate(U)
    deallocate(work)
    deallocate(nto_p)
    deallocate(nto_h)
    
    return
    
  end subroutine adc2_nto_gs
    
!######################################################################
  
end module nto
