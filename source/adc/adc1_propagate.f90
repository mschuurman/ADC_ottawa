!#######################################################################
! TD-ADC(1) and TD-CIS wavepacket propagation including the interaction
! with an applied laser pulse
!
! Note that we peform the full diagonalisation of the ADC(1)/CIS
! Hamiltonian in this module.
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
    use thetamod
    use propagate_adc1
    
    implicit none

    integer, dimension(:,:), allocatable  :: kpq,kpqd,kpqf
    integer                               :: i,ndim,ndims,ndimsf,&
                                             nout,ndimf,ndimd,noutf
    integer*8                             :: noffd,noffdf
    real(d)                               :: e0
    real(d), dimension(:,:), allocatable  :: cap_mo,theta_mo
    type(gam_structure)                   :: gam
    
!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!
! Do we really want to be doing this for an ADC(1) calculation? Note
! that the MP1 ground state is simply the HF ground state.
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
! If the initial state is an excited state, then diagonalise the
! initial state Hamiltonian
!-----------------------------------------------------------------------
    if (statenumber.gt.0.or.iprojcap.eq.2) &
         call get_initial_state_adc1(kpq,ndim)

!-----------------------------------------------------------------------
! Calculate the final space Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian_incore(kpqf,ndimf)

!-----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!-----------------------------------------------------------------------
    if (lcap) call cap_mobas(gam,cap_mo)

!-----------------------------------------------------------------------
! If a projected CAP is being used, determine which states are to be
! included in the projector
!-----------------------------------------------------------------------
    if (lprojcap) call get_proj_states_adc1(ndim)
    
!-----------------------------------------------------------------------
! If flux analysis is to be performed, then calculate the MO
! representation of the projector (Theta) onto the CAP region
!-----------------------------------------------------------------------
    if (lflux) call theta_mobas(gam,theta_mo)
    
!-----------------------------------------------------------------------
! Calculate the matrix elements needed to represent the CAP operator
! in the the ground state + intermediate state basis
!-----------------------------------------------------------------------
    if (lcap) call cap_isbas_adc1(cap_mo,kpqf,ndimf)

!-----------------------------------------------------------------------
! Calculate the dipole matrices
!-----------------------------------------------------------------------
    call dipole_isbas_adc1(kpqf,ndimf)

!-----------------------------------------------------------------------
! If flux analysis is to be performed, then calculate the matrix
! elements needed to represent the projector onto the CAP region in
! the ground state + intermediate state basis
!-----------------------------------------------------------------------
    if (lflux) call theta_isbas_adc1(theta_mo,kpqf,ndimf)
    
!-----------------------------------------------------------------------
! Perform the wavepacket propagation
!-----------------------------------------------------------------------
    call propagate_laser_adc1(ndimf,kpqf)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    deallocate(h1)
    if (allocated(cap_mo)) deallocate(cap_mo)
    if (allocated(w0j)) deallocate(w0j)
    if (allocated(wij)) deallocate(wij)
    deallocate(d0j)
    deallocate(dij)
    if (allocated(theta_mo)) deallocate(theta_mo)
    if (allocated(theta0j)) deallocate(theta0j)
    if (allocated(thetaij)) deallocate(thetaij)
    if (allocated(projmask)) deallocate(projmask)
    
    return
    
  end subroutine adc1_propagate

!#######################################################################

  subroutine calc_hamiltonian_incore(kpqf,ndimf)

    use parameters
    use constants
    use fspace
    use iomod
    
    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer                                             :: ndimf,i,j
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! ADC(1) Hamiltonian matrix
    allocate(h1(ndimf,ndimf))
    h1=0.0d0

!-----------------------------------------------------------------------
! Full calculation and in-core storage of the ADC(1) Hamiltonian
! Matrix
!-----------------------------------------------------------------------
    if (method.eq.1) then
       ! ADC(1)
       if (lcvsfinal) then
          write(ilog,'(/,2x,a)') 'Calculating the CVS-ADC(1) &
               Hamiltonian matrix'
          call get_fspace_tda_direct_nodiag_cvs(ndimf,kpqf,h1)
       else
          write(ilog,'(/,2x,a)') 'Calculating the ADC(1) &
               Hamiltonian matrix'
          call get_fspace_tda_direct_nodiag(ndimf,kpqf,h1)
       endif
    else if (method.eq.4) then
       ! ADC(1)-x
       ! (Note that the CVS-ADC(1) Hamiltonian can be
       ! calculated using the ADC(1) routines)
       call get_fspace_adc1ext_direct_nodiag(ndimf,kpqf,h1)
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
! Note that if the ground state is included in the CAP projector, Q,
! then these elements are zero
!----------------------------------------------------------------------
    allocate(w0j(ndimf))
    w0j=0.0d0

    if (.not.lprojcap &
         .or.(statenumber.gt.0.and.iprojcap.eq.1)) then
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the vector &
            W_0J = < Psi_0 | W | Psi_J >'
       write(ilog,'(72a)') ('-',k=1,72)

       if (lcis) then
          ! CIS
          call get_tm_cis(ndimf,kpqf,w0j)
       else if (method.eq.4) then
          ! ADC(1)-x
          call get_modifiedtm_adc1ext(ndimf,kpqf,w0j,1)
       else
          ! ADC(1)
          call get_modifiedtm_tda(ndimf,kpqf,w0j)
       endif
          
    endif

!----------------------------------------------------------------------
! Calculate the IS representation of the shifted CAP operator W-W_00
!
! Note that we are here assuming that the ADC(1) D-matrix is small
! enough to fit into memory
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         shifted CAP operator'
    write(ilog,'(72a,/)') ('-',k=1,72)
    
    allocate(wij(ndimf,ndimf))
    wij=0.0d0

    if (method.eq.4) then
       ! ADC(1)-x
       call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,wij)
    else
       ! ADC(1) and CIS
       call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,wij)
    endif

!----------------------------------------------------------------------
! Reset the dpl array
!----------------------------------------------------------------------
    dpl(1:nbas,1:nbas)=dpl_orig(1:nbas,1:nbas)
    
    return

  end subroutine cap_isbas_adc1

!#######################################################################

  subroutine dipole_isbas_adc1(kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: i,p,q,k,c
    real(d), dimension(nbas,nbas)             :: rho0
    character(len=60)                         :: filename
    character(len=1), dimension(3)            :: acomp

    acomp=(/ 'x','y','z' /)

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

       if (lcis) then
          ! CIS
          call get_tm_cis(ndimf,kpqf,d0j(c,:))
       else if (method.eq.4) then
          ! ADC(1)-x
          call get_modifiedtm_adc1ext(ndimf,kpqf,d0j(c,:),1)
       else
          ! ADC(1)
          call get_modifiedtm_tda(ndimf,kpqf,d0j(c,:))
       endif
          
    enddo
    
!----------------------------------------------------------------------
! Calculate the IS representations of the shifted dipole operators
! Dc - Dc_0, c=x,y,z
!----------------------------------------------------------------------
    allocate(dij(3,ndimf,ndimf))
    dij=0.0d0
    
    ! Loop over the x, y, and z components
    do c=1,3

       ! Skip if the current component is not required
       if (pulse_vec(c).eq.0.0d0) cycle
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the IS representation of &
            the shifted dipole operator D'//acomp(c)
       write(ilog,'(72a)') ('-',k=1,72)

       dpl(:,:)=dpl_all(c,:,:)

       if (method.eq.4) then
          ! ADC(1)-x
          call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,dij(c,:,:))
       else
          ! ADC(1) and CIS
          call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,&
               dij(c,:,:))
       endif
          
    enddo

    return
    
  end subroutine dipole_isbas_adc1

!#######################################################################

  subroutine theta_isbas_adc1(theta_mo,kpqf,ndimf)

    use constants
    use parameters
    use mp2
    use get_matrix_dipole
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer                                   :: i,p,q,k
    real(d), dimension(nbas,nbas)             :: theta_mo
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
! Calculate the ground state-ground state projector matrix element
! Theta_00 = < Psi_0 | Theta | Psi_0 >
!----------------------------------------------------------------------
    theta00=0.0d0
    do p=1,nbas
       do q=1,nbas
          theta00=theta00+rho0(p,q)*theta_mo(p,q)
       enddo
    enddo

!----------------------------------------------------------------------
! In the following, we calculate projector matrix elements using the
! shifted dipole code (D-matrix and f-vector code) by simply
! temporarily copying the MO Theta matrix into the dpl array.
!----------------------------------------------------------------------
    dpl_orig(1:nbas,1:nbas)=dpl(1:nbas,1:nbas)
    dpl(1:nbas,1:nbas)=theta_mo(1:nbas,1:nbas)

!----------------------------------------------------------------------
! Calculate the vector Theta_0J = < Psi_0 | Theta | Psi_J >
!----------------------------------------------------------------------
    allocate(theta0j(ndimf))
    theta0j=0.0d0

    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the vector &
         Theta_0J = < Psi_0 | Theta | Psi_J >'
    write(ilog,'(72a)') ('-',k=1,72)
    
    if (lcis) then
       ! CIS
       call get_tm_cis(ndimf,kpqf,theta0j)
    else if (method.eq.4) then
       ! ADC(1)-x
       call get_modifiedtm_adc1ext(ndimf,kpqf,theta0j,1)
    else
       ! ADC(1)
       call get_modifiedtm_tda(ndimf,kpqf,theta0j)
    endif

!----------------------------------------------------------------------
! Calculate the IS representation of the shifted projection operator
! Theta-Theta_00
!
! Note that we are here assuming that the ADC(1) D-matrix is small
! enough to fit into memory
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         shifted CAP-projector (Theta)'
    write(ilog,'(72a,/)') ('-',k=1,72)
    
    allocate(thetaij(ndimf,ndimf))
    thetaij=0.0d0

    if (method.eq.4) then
       ! ADC(1)-x
       call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,thetaij)
    else
       ! ADC(1) and CIS
       call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,&
            thetaij)
    endif

!----------------------------------------------------------------------
! Reset the dpl array
!----------------------------------------------------------------------
    dpl(1:nbas,1:nbas)=dpl_orig(1:nbas,1:nbas)
    
    return
    
  end subroutine theta_isbas_adc1

!#######################################################################

  subroutine get_initial_state_adc1(kpq,ndim)

    use constants
    use parameters
    use fspace
    use get_moment
    use misc
    use iomod
    
    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2) :: kpq
    integer                                 :: ndim,i,itmp,unit
    real(d), dimension(:,:), allocatable    :: eigvec
    real(d), dimension(:), allocatable      :: ener,mtm,tmvec,osc_str

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(eigvec(ndim,ndim))
    allocate(ener(ndim))
    allocate(mtm(ndim))
    allocate(tmvec(ndim))
    allocate(osc_str(ndim))
    
!-----------------------------------------------------------------------
! Diagonalisation in the initial space
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') "Full diagonalisation of the ADC(1)/CIS &
         Hamiltonian matrix..."
    call get_fspace_tda_direct(ndim,kpq,eigvec,ener)

!-----------------------------------------------------------------------
! Transition dipole moments between the ground state and the ADC(1)/CIS
! eigenstates
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Calculating the transition dipole &
         moments between the ground state and all excited states...'
    if (lcis) then
       call get_tm_cis(ndim,kpq(:,:),mtm(:))
    else
       call get_modifiedtm_tda(ndim,kpq(:,:),mtm(:))
    endif
       
    do i=1,ndim
       tmvec(i)=tm(ndim,eigvec(:,i),mtm(:))
       osc_str(i)=2.0d0/3.0d0*ener(i)*tmvec(i)**2
    enddo

    write(ilog,'(/,70a)') ('*',i=1,70)
    write(ilog,'(2x,a)') &
         'Initial space ADC(1)/CIS excitation energies'
    write(ilog,'(70a)') ('*',i=1,70)
    
    itmp=1+nBas**2*4*nOcc**2
    call table2(ndim,statenumber+5,ener(1:statenumber+5),&
         eigvec(:,1:statenumber+5),tmvec(1:statenumber+5),&
         osc_str(1:statenumber+5),kpq,itmp,'i')

!-----------------------------------------------------------------------
! Write the eigenpairs to file
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='SCRATCH/initvecs',status='unknown',&
         access='sequential',form='unformatted')

    do i=1,ndim
       write(unit) i,ener(i),eigvec(:,i)
    enddo

    close(unit)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(eigvec)
    deallocate(ener)
    deallocate(mtm)
    deallocate(tmvec)
    deallocate(osc_str)
    
    return
    
  end subroutine get_initial_state_adc1

!#######################################################################

  subroutine get_proj_states_adc1(ndim)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer              :: ndim,unit,itmp,i
    real(d), allocatable :: ener(:)
    real(d), allocatable :: vec(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(projmask(ndim))
    projmask=0
    
    allocate(vec(ndim))
    vec=0.0d0

    allocate(ener(ndim))
    ener=0.0d0
    
!----------------------------------------------------------------------
! Fill in the projmask array
!----------------------------------------------------------------------
    projmask=0

    if (iprojcap.eq.2) then
       call freeunit(unit)
       open(unit,file='SCRATCH/initvecs',status='old',&
            access='sequential',form='unformatted')
       do i=1,ndim
          read(unit) itmp,ener(i),vec
          if (ener(i).le.projlim) then
             projmask(i)=1
          else
             exit
          endif
       enddo
       close(unit)
    else if (iprojcap.eq.1.and.statenumber.gt.0) then
       projmask(statenumber)=1
    endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    deallocate(ener)
    
    return
    
  end subroutine get_proj_states_adc1
    
!#######################################################################
  
end module adc1propmod
