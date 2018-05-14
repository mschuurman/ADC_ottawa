module adc1enermod

  use channels
  
contains

!#######################################################################

  subroutine adc1_ener

    use constants
    use parameters
    use fspace
    use misc
    use mp2
    use adc_common
    
    implicit none
    
    integer, dimension(:,:), allocatable  :: kpq
    integer                               :: ndim,ndimd
    integer                               :: i,itmp
    real(dp)                              :: time
    real(dp), dimension(:), allocatable   :: ener,mtm,tmvec,osc_str
    real(dp), dimension(:,:), allocatable :: arr
    
!-----------------------------------------------------------------------
! Determine the 1h1p subspace
!-----------------------------------------------------------------------
    allocate(kpq(7,0:nBas**2*4*nOcc**2))
    kpq(:,:)=-1

    if (lcvs) then
       if (lfakeip) then
          ! CVS-IP-ADC(1)
          call select_atom_is_cvs_fakeip(kpq(:,:))
       else
          ! CVS-ADC(1)
          call select_atom_is_cvs(kpq(:,:))
       endif
    else if (lfakeip) then
       ! IP-ADC(1)
       call select_atom_is_fakeip(kpq(:,:))
    else
       ! ADC(1)
       call select_atom_is(kpq(:,:))
    endif

!-----------------------------------------------------------------------
! Output supspace dimensions
!-----------------------------------------------------------------------
    ndim=kpq(1,0)
    
    write(ilog,'(/)')
    write(ilog,*) 'ADC(1) INITIAL Space dim',ndim

!-----------------------------------------------------------------------
! Set the dipole matrix
!-----------------------------------------------------------------------
    call set_dpl

!-----------------------------------------------------------------------
! Diagonalisation in the initial space
!-----------------------------------------------------------------------
    allocate(arr(ndim,ndim),mtm(ndim),tmvec(ndim),osc_str(ndim))
    allocate(ener(ndim))

    write(ilog,'(/,2x,a)') "Full diagonalisation of the ADC(1) &
         Hamiltonian matrix..."
    
    call get_fspace_tda_direct(ndim,kpq(:,:),arr,ener)

!-----------------------------------------------------------------------
! Transition dipole moments between the ground state and the 1h1p ISs
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Calculating the transition dipole &
         moments between the ground state and all excited states...'
    if (lcis) then
       call get_tm_cis(ndim,kpq(:,:),mtm(:))
    else
       call get_modifiedtm_tda(ndim,kpq(:,:),mtm(:))
    endif
       
    do i=1,ndim
       tmvec(i)=tm(ndim,arr(:,i),mtm(:))
       osc_str(i)=2.0d0/3.0d0*ener(i)*tmvec(i)**2
    enddo

    write(ilog,'(/,70a)') ('*',i=1,70)
    if (lcis) then
       write(ilog,'(2x,a)') &
            'Initial space CIS excitation energies'
    else
       write(ilog,'(2x,a)') &
            'Initial space ADC(1) excitation energies'
    endif
    write(ilog,'(70a)') ('*',i=1,70)
    itmp=1+nBas**2*4*nOcc**2
    call table2(ndim,ndim,ener(:),arr(:,:),tmvec(:),osc_str(:),kpq,&
         itmp,'i')

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(kpq)
    deallocate(arr,mtm,tmvec,osc_str)
    deallocate(ener)
    
    return
      
  end subroutine adc1_ener
      
!#######################################################################
  
end module adc1enermod
