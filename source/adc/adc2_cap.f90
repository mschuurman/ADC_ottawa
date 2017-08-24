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
! Calculate the matrix elements needed to represent the CAP operator
! in the the ground state + intermediate state basis
!-----------------------------------------------------------------------
    call cap_isbas(cap_mo,kpqf,ndimf)

!-----------------------------------------------------------------------
! Calculate the dipole matrices
!-----------------------------------------------------------------------
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    deallocate(cap_mo)
    deallocate(w0j)
    
    return
    
  end subroutine adc2_cap

!#######################################################################

  subroutine cap_isbas(cap_mo,kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: p,q,k
    real(d), dimension(nbas,nbas)             :: cap_mo
    real(d), dimension(nbas,nbas)             :: rho0
    real(d), dimension(nbas,nbas)             :: dpl_orig
    character(len=60)                         :: filename

!----------------------------------------------------------------------
! Calculate the ground state density matrix
!----------------------------------------------------------------------
    call rhogs(rho0)

!----------------------------------------------------------------------
! Calculate the CAP matrix element W_00 = < Psi_0 | W | Psi_0 >
!----------------------------------------------------------------------
    w00=0.0d0
    do p=1,nbas
       do q=1,nbas
          w00=w00+rho0(p,q)*cap_mo(p,q)
       enddo
    enddo

!----------------------------------------------------------------------
! In the following, we calculate CAP matrix elements using the shifted
! dipole code (D-matrix and f-vector code) by simply temporarily
! copying the MO CAP matrix into the dpl array.
!----------------------------------------------------------------------
    dpl_orig(1:nbas,1:nbas)=dpl(1:nbas,1:nbas)
    dpl(1:nbas,1:nbas)=cap_mo(1:nbas,1:nbas)
    
!----------------------------------------------------------------------
! Calculate the IS representation of the shifted CAP operator W-W_00
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         shifted CAP operator'
    write(ilog,'(72a,/)') ('-',k=1,72)
    
    filename='SCRATCH/cap_isr'
    
    call get_adc2_dipole_improved_omp(ndimf,ndimf,kpqf,kpqf,&
         nbuf_cap,nel_cap,filename)

!----------------------------------------------------------------------
! Calculate the vector W_0J = < Psi_0 | W | Psi_J >
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the vector &
         W_0J = < Psi_0 | W | Psi_J >'
    write(ilog,'(72a)') ('-',k=1,72) 

    allocate(w0j(ndimf))
    w0j=0.0d0
    
    call get_modifiedtm_adc2(ndimf,kpqf(:,:),w0j,1)
    
!----------------------------------------------------------------------
! Reset the dpl array
!----------------------------------------------------------------------
    dpl(1:nbas,1:nbas)=dpl_orig(1:nbas,1:nbas)
    
    return
    
  end subroutine cap_isbas
    
!#######################################################################
  
end module adc2capmod
