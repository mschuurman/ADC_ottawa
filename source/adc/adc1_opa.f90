module adc1specmod

  use channels
    
contains

!#######################################################################

  subroutine adc1_spec()

    use constants
    use parameters
    use fspace
    use misc
    use mp2
    use adc2common
    
    implicit none        

    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimsf,&
                                            nout,ndimf,ndimd,&
                                            noutf,itmp
    real(d)                              :: time
    real(d), dimension(:), allocatable   :: ener,mtm,tmvec,osc_str
    real(d), dimension(:), allocatable   :: travec
    real(d)                              :: e_init,e0,s0
    real(d), dimension(:,:), allocatable :: rvec
    real(d), dimension(:), allocatable   :: vec_init
    real*8, dimension(:), allocatable    :: mtmf
    
    real(d), dimension(:,:), allocatable :: arr,arrd,arrf
    real(d), dimension(:), allocatable   :: autvec,tmvecf,osc_strf,&
                                            enerf

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic (if requested)
!-----------------------------------------------------------------------
    call mp2_master(e0)

!-----------------------------------------------------------------------
! Determine the 1h1p subspace
!-----------------------------------------------------------------------
    call get_subspaces_adc1(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf)

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
    
    itmp=1+nBas**2*4*nOcc**2
    call table2(ndim,ndim,ener(:),arr(:,:),tmvec(:),&
         osc_str(:),kpq,itmp,'i')
    
!-----------------------------------------------------------------------
! Output transition energies and oscillator strentghs to file
!-----------------------------------------------------------------------
    call get_sigma(ndim,ener(:),osc_str(:))

!-----------------------------------------------------------------------
! For checking purposes, calculate S(0)
!-----------------------------------------------------------------------
    s0=0.0d0
    do i=1,ndim
       s0=s0+osc_str(i)
    enddo

    write(6,'(/,2x,a,2x,F10.7,/)') "S(0):",s0    

    return

  end subroutine adc1_spec

!#######################################################################
      
end module adc1specmod
