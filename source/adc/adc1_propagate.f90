!#######################################################################
! TD-ADC(1) wavepacket propagation including the interaction with an
! applied laser pulse
!#######################################################################

module adc1propmod

  use channels

contains

!#######################################################################

  subroutine adc1_propagate(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    use capmod
    use propagate
    
    implicit none

    integer, dimension(:,:), allocatable  :: kpq,kpqd,kpqf
    integer                               :: i,ndim,ndims,ndimsf,&
                                             nout,ndimf,ndimd,noutf
    integer*8                             :: noffd,noffdf
    real(d)                               :: e0
    real(d), dimension(:,:), allocatable  :: cap_mo,hmat,eigvec
    type(gam_structure)                   :: gam

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!-----------------------------------------------------------------------
    call mp2_master(e0)

!-----------------------------------------------------------------------
! Determine the 1h1p subspace
!-----------------------------------------------------------------------
    call get_subspaces_adc1(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf)

!-----------------------------------------------------------------------
! Set MO representation of the dipole operator
!-----------------------------------------------------------------------
    call set_dpl

!-----------------------------------------------------------------------
! Calculate the final space Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian_incore(kpqf,ndimf,hmat,eigvec)

!-----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!-----------------------------------------------------------------------
    if (lcap) call cap_mobas(gam,cap_mo)
    
!-----------------------------------------------------------------------
! Calculate the matrix elements needed to represent the CAP operator
! in the the ground state + intermediate state basis
!-----------------------------------------------------------------------
    if (lcap) call cap_isbas_adc1(cap_mo,kpqf,ndimf)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    deallocate(hmat)
    deallocate(eigvec)
    if (allocated(cap_mo)) deallocate(cap_mo)
    if (allocated(w0j)) deallocate(w0j)
    
    return
    
  end subroutine adc1_propagate

!#######################################################################

  subroutine calc_hamiltonian_incore(kpqf,ndimf,hmat,eigvec)

    use parameters
    use constants
    use fspace
    
    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer                              :: ndimf
    real(d), dimension(:,:), allocatable :: hmat,eigvec
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! Hamiltonian matrix
    allocate(hmat(ndimf,ndimf))
    hmat=0.0d0

    ! ADC(1) eigenvectors
    allocate(eigvec(ndimf,ndimf))
    eigvec=0.0d0

!-----------------------------------------------------------------------
! Full calculation and in-core storage of the ADC(1) Hamiltonian
! Matrix
!-----------------------------------------------------------------------
    if (lcvsfinal) then
       write(ilog,'(/,2x,a)') 'Calculating the CVS-ADC(1) &
            Hamiltonian matrix'
       call get_fspace_tda_direct_cvs(ndimf,kpqf,hmat,eigvec)
    else
       write(ilog,'(/,2x,a)') 'Calculating the ADC(1) &
            Hamiltonian matrix'
       call get_fspace_tda_direct(ndimf,kpqf,hmat,eigvec)
    endif
       
    return
    
  end subroutine calc_hamiltonian_incore

!#######################################################################

  subroutine cap_isbas_adc1(cap_mo,kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: i,p,q,k
    real(d), dimension(nbas,nbas)             :: cap_mo
    real(d), dimension(nbas,nbas)             :: rho0
    real(d), dimension(nbas,nbas)             :: dpl_orig
    character(len=60)                         :: filename

!----------------------------------------------------------------------
! Ground state density matrix.
! Note that the 1st-order correction is zero.
!----------------------------------------------------------------------
    rho0=0.0d0

    ! Occupied-occupied block: 0th-order contribution
    do i=1,nocc
       rho0(i,i)=2.0d0
    enddo

!----------------------------------------------------------------------
! Calculate the CAP matrix element W_00 = < Psi_0 | W | Psi_0 >
!----------------------------------------------------------------------
    w00=0.0d0
    do p=1,nbas
       do q=1,nbas
          w00=w00+rho0(p,q)*cap_mo(p,q)
       enddo
    enddo

    print*,
    print*,'WRITE THE REST OF THE ADC(1) CAP MATRIX CODE!'
    print*,
    
    return

  end subroutine cap_isbas_adc1
    
!#######################################################################
  
end module adc1propmod
