!#######################################################################
! TD-ADC(2) wavepacket propagation including the interaction with an
! applied laser pulse
!#######################################################################

module adc2propmod

  use channels

contains

!#######################################################################

  subroutine adc2_propagate(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    use propagate_adc2
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
! Calculate the final space Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian(ndimf,kpqf,noffdf)
    
!-----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!-----------------------------------------------------------------------
    if (lcap) call cap_mobas(gam,cap_mo)
    
!-----------------------------------------------------------------------
! Calculate the matrix elements needed to represent the CAP operator
! in the the ground state + intermediate state basis
!-----------------------------------------------------------------------
    if (lcap) call cap_isbas_adc2(cap_mo,kpqf,ndimf)

!-----------------------------------------------------------------------
! Calculate the dipole matrices
!-----------------------------------------------------------------------
    call dipole_isbas_adc2(kpqf,ndimf)

!-----------------------------------------------------------------------
! Perform the wavepacket propagation
!-----------------------------------------------------------------------
    hamflag='f'
    call propagate_laser_adc2(ndimf,noffdf,kpqf)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    if (allocated(cap_mo)) deallocate(cap_mo)
    if (allocated(w0j)) deallocate(w0j)
    deallocate(d0j)
    deallocate(dpl_all)
    
    return
    
  end subroutine adc2_propagate

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

  subroutine cap_isbas_adc2(cap_mo,kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: p,q,k
    integer                                   :: error
    real(d), dimension(nbas,nbas)             :: cap_mo
    real(d), dimension(nbas,nbas)             :: rho0
    real(d), dimension(nbas,nbas)             :: dpl_orig
    character(len=60)                         :: filename

!----------------------------------------------------------------------
! Calculate the ground state density matrix
!----------------------------------------------------------------------
    call rho_mp2(rho0)

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
! Calculate the vector W_0J = < Psi_0 | W | Psi_J >
!
! Note that if the projected CAP is being used and the initial state is
! the ground state, then these matrix elements are zero
!----------------------------------------------------------------------
    allocate(w0j(ndimf))
    w0j=0.0d0

    if (.not.lprojcap.or.statenumber.gt.0) then

       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the vector &
            W_0J = < Psi_0 | W | Psi_J >'
       write(ilog,'(72a)') ('-',k=1,72)
       
       call get_modifiedtm_adc2(ndimf,kpqf(:,:),w0j,1)

    endif

!----------------------------------------------------------------------
! Calculate the IS representation of the shifted CAP operator W-W_00
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         shifted CAP operator'
    write(ilog,'(72a,/)') ('-',k=1,72)
    
    filename='SCRATCH/cap'
    
    call get_adc2_dipole_improved_omp(ndimf,ndimf,kpqf,kpqf,&
         nbuf_cap,nel_cap,filename)
    
!----------------------------------------------------------------------
! Reset the dpl array
!----------------------------------------------------------------------
    dpl(1:nbas,1:nbas)=dpl_orig(1:nbas,1:nbas)
    
    return
    
  end subroutine cap_isbas_adc2
    
!#######################################################################

  subroutine dipole_isbas_adc2(kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: p,q,k,c
    real(d), dimension(nbas,nbas)             :: rho0
    character(len=60)                         :: filename
    character(len=1), dimension(3)            :: acomp

    acomp=(/ 'x','y','z' /)

!----------------------------------------------------------------------
! Calculate the ground state density matrix
!----------------------------------------------------------------------
    call rho_mp2(rho0)

!----------------------------------------------------------------------
! Calculate the dipole matrix elements Dc_00 = < Psi_0 | Dc | Psi_0 >
! for c=x,y,z
!----------------------------------------------------------------------
    d00=0.0d0
    do p=1,nbas
       do q=1,nbas
          d00(1:3)=d00(1:3)+rho0(p,q)*dpl_all(1:3,p,q)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the vectors Dc_0J = < Psi_0 | D | Psi_J >, c=x,y,z
!----------------------------------------------------------------------
    allocate(d0j(3,ndimf))
    d0j=0.0d0

    ! Loop over the x, y, and z components
    do c=1,3

       ! Skip if the current component is not required
       if (pulse_vec(c).eq.0.0d0) cycle
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the vector D'//acomp(c)//&
            '_0J = < Psi_0 | D'//acomp(c)//' | Psi_J >'
       write(ilog,'(72a)') ('-',k=1,72)

       dpl(:,:)=dpl_all(c,:,:)
       
       call get_modifiedtm_adc2(ndimf,kpqf(:,:),d0j(c,:),1)

    enddo
        
!----------------------------------------------------------------------
! Calculate the IS representations of the shifted dipole operators
! Dc - Dc_0, c=x,y,z
!----------------------------------------------------------------------
    ! Loop over the x, y, and z components
    do c=1,3

       ! Skip if the current component is not required
       if (pulse_vec(c).eq.0.0d0) cycle
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the IS representation of &
            the shifted dipole operator D'//acomp(c)
       write(ilog,'(72a)') ('-',k=1,72)
       
       filename='SCRATCH/dipole_'//acomp(c)

       dpl(:,:)=dpl_all(c,:,:)
       
       call get_adc2_dipole_improved_omp(ndimf,ndimf,kpqf,kpqf,&
            nbuf_dip(c),nel_dip(c),filename)
       
    enddo

    return
    
  end subroutine dipole_isbas_adc2
    
!#######################################################################
  
end module adc2propmod
