!#######################################################################
! Calculation of the electronic absorption spectrum from the Fourier
! transform of the autocorrelation function obtained by propagation
! of |Psi(t=0)> = D |Psi_0>
!#######################################################################

module adc2automod

  use channels

contains

!#######################################################################

  subroutine adc2_autospec(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    use fvecprop
    
    implicit none
    
    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimf,&
                                            ndimsf,nout,noutf,ndimd
    integer*8                            :: noffd,noffdf
    real(d)                              :: e0,time
    real(d), dimension(:), allocatable   :: fvec,dpsi
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
! Contraction of the dipole operator with the initial state
!-----------------------------------------------------------------------
    call contract_dipole_initial_state(dpsi,ndim,ndims,ndimf,kpq,kpqf,&
         noffd)
    
!-----------------------------------------------------------------------
! Calculate the final space Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian(ndimf,kpqf,noffdf)
    
!-----------------------------------------------------------------------
! Perform the wavepacket propagation and autocorrelation function
! calculation
!-----------------------------------------------------------------------
    hamflag='f'
    call propagate_fvec(dpsi,ndimf,noffdf)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(dpsi)

    return
      
  end subroutine adc2_autospec
    
!#######################################################################

  subroutine calc_hamiltonian(ndimf,kpqf,noffdf)

    use fspace
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer*8                                 :: noffdf
    
    write(ilog,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
    
    if (method.eq.2) then
       ! ADC(2)-s
       if (lcvsfinal) then
          call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
       else
          call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
       endif
    else if (method.eq.3) then
       ! ADC(2)-x
       if (lcvsfinal) then
          call write_fspace_adc2e_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
       else
          call write_fspace_adc2e_1(ndimf,kpqf(:,:),noffdf,'c')
       endif
    endif
    
    return
    
  end subroutine calc_hamiltonian

!#######################################################################

  subroutine contract_dipole_initial_state(dpsi,ndim,ndims,ndimf,&
       kpq,kpqf,noffd)

    use constants
    use parameters
    use get_moment
    use adc2common
    use get_matrix_dipole
    use guessvecs
    use misc
    use diagmod
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: i
    integer                                   :: ndim,ndims,ndimf,itmp
    integer*8                                 :: noffd
    real(d), dimension(:), allocatable        :: dpsi,ener,mtm,tmvec,&
                                                 osc_str
    real(d), dimension(:,:), allocatable      :: rvec
    real(d)                                   :: time
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(dpsi(ndimf))
    
!-----------------------------------------------------------------------
! Initial state = | Psi_0 >
! Calculate the f-vector, F_J = < Psi_J | D | Psi_0 >    
!-----------------------------------------------------------------------
    if (statenumber.eq.0) call get_modifiedtm_adc2(ndimf,kpqf(:,:),&
         dpsi(:),1)

!-----------------------------------------------------------------------
! Excited initial state: diagonalisation in the initial space and
! contraction of the D-matrix with the initial state vector
!-----------------------------------------------------------------------
    if (statenumber.gt.0) then

       ! (1) Initial space diagonalisation
       if (ladc1guess) call adc1_guessvecs
       call initial_space_diag(time,kpq,ndim,ndims,noffd)
       call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,kpq)
       
       ! (2) Output the initial space state information
       write(ilog,'(/,70a)') ('*',i=1,70)
       write(ilog,'(2x,a)') &
            'Initial space ADC(2)-s excitation energies'
       write(ilog,'(70a)') ('*',i=1,70)
       itmp=1+nBas**2*4*nOcc**2
       call table2(ndim,davstates,ener(1:davstates),&
            rvec(:,1:davstates),tmvec(1:davstates),&
            osc_str(1:davstates),kpq,itmp,'i')
       write(ilog,'(/,70a,/)') ('*',i=1,70)
       
       ! (3) Contraction of the D-matrix with the initial state vector
       call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,&
            rvec(:,statenumber),dpsi)

       ! Deallocate arrays that are no longer needed (N.B. these
       ! are allocated in initial_space_tdm)
       deallocate(ener,rvec,tmvec,osc_str)

    endif
    
    return
    
  end subroutine contract_dipole_initial_state
    
!#######################################################################
  
end module adc2automod
