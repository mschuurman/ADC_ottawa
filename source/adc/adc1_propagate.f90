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
    use adc_common
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
    real(dp)                              :: e0
    real(dp), dimension(:,:), allocatable :: cap_mo,theta_mo
    type(gam_structure)                   :: gam

!-----------------------------------------------------------------------
! Determine the 1h1p subspace
!-----------------------------------------------------------------------
    call get_subspaces_adc1(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf)

!-----------------------------------------------------------------------
! For now, we will take the initial space to be equal to the final space
!-----------------------------------------------------------------------
    kpq=kpqf
    ndim=ndimf
    nout=noutf
    ndims=ndimsf
    
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
    if (npulse.gt.0) call dipole_isbas_adc1(kpqf,ndimf)

!-----------------------------------------------------------------------
! If flux analysis is to be performed, then calculate the matrix
! elements needed to represent the projector onto the CAP region in
! the ground state + intermediate state basis
!-----------------------------------------------------------------------
    if (lflux) call theta_isbas_adc1(theta_mo,kpqf,ndimf)

!-----------------------------------------------------------------------
! Transformation of operator matrices to the field-free Hamiltonian
! eigenstate basis
!-----------------------------------------------------------------------
    if (tdrep.eq.2) call transform_operators(ndimf)
    
!-----------------------------------------------------------------------
! Diagonalisation of the CAP-augmented Hamiltonian for analysis
! purposes.
!
! Note that this needs to be done before any CAP projection takes place.
!-----------------------------------------------------------------------
    if (lcapdiag) call diag_hcap_adc1(ndimf)

!-----------------------------------------------------------------------
! Projection of the CAP and Theta matrices
!-----------------------------------------------------------------------
    if (lprojcap) call cap_projection(ndimf)

!-----------------------------------------------------------------------
! Perform the wavepacket propagation
!-----------------------------------------------------------------------
    call propagate_laser_adc1(ndimf,kpqf)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------    
    deallocate(kpq,kpqf,kpqd)
    deallocate(h1)
    if (allocated(d1)) deallocate(d1)
    if (allocated(w1)) deallocate(w1)
    if (allocated(theta1)) deallocate(theta1)
    if (allocated(cap_mo)) deallocate(cap_mo)
    if (allocated(theta_mo)) deallocate(theta_mo)
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
    allocate(h1(ndimf+1,ndimf+1))
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
          call get_fspace_tda_direct_nodiag_cvs(ndimf,kpqf,&
               h1(1:ndimf,1:ndimf))
       else
          write(ilog,'(/,2x,a)') 'Calculating the ADC(1) &
               Hamiltonian matrix'
          call get_fspace_tda_direct_nodiag(ndimf,kpqf,&
               h1(1:ndimf,1:ndimf))
       endif
    else if (method.eq.4) then
       ! ADC(1)-x
       ! (Note that the CVS-ADC(1) Hamiltonian can be
       ! calculated using the ADC(1) routines)
       call get_fspace_adc1ext_direct_nodiag(ndimf,kpqf,&
            h1(1:ndimf,1:ndimf))
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
    real(dp), dimension(nbas,nbas)            :: cap_mo
    real(dp), dimension(nbas,nbas)            :: rho0
    real(dp), dimension(nbas,nbas)            :: dpl_orig
    character(len=60)                         :: filename

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(w1(ndimf+1,ndimf+1))
    w1=0.0d0
    
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
    do p=1,nbas
       do q=1,nbas
          w1(ndimf+1,ndimf+1)=w1(ndimf+1,ndimf+1)+rho0(p,q)*cap_mo(p,q)
       enddo
    enddo
 
!----------------------------------------------------------------------
! In the following, we calculate CAP matrix elements using the
! dipole code (D-matrix and f-vector code) by simply temporarily
! copying the MO CAP matrix into the dpl array.
!----------------------------------------------------------------------
    dpl_orig(1:nbas,1:nbas)=dpl(1:nbas,1:nbas)
    dpl(1:nbas,1:nbas)=cap_mo(1:nbas,1:nbas)

!----------------------------------------------------------------------
! Calculate the vector W_0J = < Psi_0 | W | Psi_J >
!----------------------------------------------------------------------
    if (lcis) then
       ! CIS
       call get_tm_cis(ndimf,kpqf,w1(1:ndimf,ndimf+1))
    else if (method.eq.4) then
       ! ADC(1)-x
       call get_modifiedtm_adc1ext(ndimf,kpqf,w1(1:ndimf,ndimf+1),1)
    else
       ! ADC(1)
       call get_modifiedtm_tda(ndimf,kpqf,w1(1:ndimf,ndimf+1))
    endif

    ! Hermitian conjugate
    w1(ndimf+1,1:ndimf)=w1(1:ndimf,ndimf+1)
    
!----------------------------------------------------------------------
! Calculate the IS representation of the CAP operator W
!
! Note that we are here assuming that the ADC(1) D-matrix is small
! enough to fit into memory
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         CAP operator'
    write(ilog,'(72a,/)') ('-',k=1,72)

    if (method.eq.4) then
       ! ADC(1)-x
       call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,&
            w1(1:ndimf,1:ndimf))
    else
       ! ADC(1) and CIS
       call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,&
            w1(1:ndimf,1:ndimf))
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
    real(dp), dimension(nbas,nbas)            :: rho0
    character(len=60)                         :: filename
    character(len=1), dimension(3)            :: acomp

    acomp=(/ 'x','y','z' /)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(d1(3,ndimf+1,ndimf+1))
    d1=0.0d0
    
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
    do p=1,nbas
       do q=1,nbas
          d1(1:3,ndimf+1,ndimf+1)=d1(1:3,ndimf+1,ndimf+1)&
               +rho0(p,q)*dpl_all(1:3,p,q)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Calculate the vectors Dc_0J = < Psi_0 | D | Psi_J >, c=x,y,z
!----------------------------------------------------------------------
    ! Loop over the x, y, and z components
    do c=1,3

       ! Skip if the current component is not required
       if (sum(abs(pulse_vec(c,1:npulse))).eq.0.0d0) cycle
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the vector D'//acomp(c)//&
            '_0J = < Psi_0 | D'//acomp(c)//' | Psi_J >'
       write(ilog,'(72a)') ('-',k=1,72)

       dpl(:,:)=dpl_all(c,:,:)

       if (lcis) then
          ! CIS
          call get_tm_cis(ndimf,kpqf,d1(c,1:ndimf,ndimf+1))
       else if (method.eq.4) then
          ! ADC(1)-x
          call get_modifiedtm_adc1ext(ndimf,kpqf,d1(c,1:ndimf,ndimf+1),1)
       else
          ! ADC(1)
          call get_modifiedtm_tda(ndimf,kpqf,d1(c,1:ndimf,ndimf+1))
       endif

       ! Hermitian conjugate
       d1(c,ndimf+1,1:ndimf)=d1(c,1:ndimf,ndimf+1)
       
    enddo
    
!----------------------------------------------------------------------
! Calculate the IS representations of the dipole operators Dc, c=x,y,z
!----------------------------------------------------------------------
    ! Loop over the x, y, and z components
    do c=1,3
    
       ! Skip if the current component is not required
       if (sum(abs(pulse_vec(c,1:npulse))).eq.0.0d0) cycle
       
       write(ilog,'(/,72a)') ('-',k=1,72)
       write(ilog,'(2x,a)') 'Calculating the IS representation of &
            the dipole operator D'//acomp(c)
       write(ilog,'(72a)') ('-',k=1,72)
    
       dpl(:,:)=dpl_all(c,:,:)
    
       if (method.eq.4) then
          ! ADC(1)-x
          call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,&
               d1(c,1:ndimf,1:ndimf))
       else
          ! ADC(1) and CIS
          call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,&
               d1(c,1:ndimf,1:ndimf))
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
    real(dp), dimension(nbas,nbas)            :: theta_mo
    real(dp), dimension(nbas,nbas)            :: rho0
    real(dp), dimension(nbas,nbas)            :: dpl_orig
    character(len=60)                         :: filename

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(theta1(ndimf+1,ndimf+1))
    theta1=0.0d0
    
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
    do p=1,nbas
       do q=1,nbas
          theta1(ndimf+1,ndimf+1)=theta1(ndimf+1,ndimf+1)&
               +rho0(p,q)*theta_mo(p,q)
       enddo
    enddo
    
!----------------------------------------------------------------------
! In the following, we calculate projector matrix elements using the
! dipole code (D-matrix and f-vector code) by simply
! temporarily copying the MO Theta matrix into the dpl array.
!----------------------------------------------------------------------
    dpl_orig(1:nbas,1:nbas)=dpl(1:nbas,1:nbas)
    dpl(1:nbas,1:nbas)=theta_mo(1:nbas,1:nbas)

!----------------------------------------------------------------------
! Calculate the vector Theta_0J = < Psi_0 | Theta | Psi_J >
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the vector &
         Theta_0J = < Psi_0 | Theta | Psi_J >'
    write(ilog,'(72a)') ('-',k=1,72)
    
    if (lcis) then
       ! CIS
       call get_tm_cis(ndimf,kpqf,theta1(1:ndimf,ndimf+1))
    else if (method.eq.4) then
       ! ADC(1)-x
       call get_modifiedtm_adc1ext(ndimf,kpqf,theta1(1:ndimf,ndimf+1),1)
    else
       ! ADC(1)
       call get_modifiedtm_tda(ndimf,kpqf,theta1(1:ndimf,ndimf+1))
    endif

    ! Hermitian conjugate
    theta1(ndimf+1,1:ndimf)=theta1(1:ndimf,ndimf+1)

!----------------------------------------------------------------------
! Calculate the IS representation of the projection operator Theta
!
! Note that we are here assuming that the ADC(1) D-matrix is small
! enough to fit into memory
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the IS representation of the &
         CAP-projector (Theta)'
    write(ilog,'(72a,/)') ('-',k=1,72)
    
    if (method.eq.4) then
       ! ADC(1)-x
       call get_adc1ext_dipole_omp(ndimf,ndimf,kpqf,kpqf,&
            theta1(1:ndimf,1:ndimf))
    else
       ! ADC(1) and CIS
       call get_offdiag_tda_dipole_direct_ok(ndimf,ndimf,kpqf,kpqf,&
            theta1(1:ndimf,1:ndimf))
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
    real(dp), dimension(:,:), allocatable   :: eigvec
    real(dp), dimension(:), allocatable     :: ener,mtm,tmvec,osc_str

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
    call table2(ndim,statenumber+15,ener(1:statenumber+15),&
         eigvec(:,1:statenumber+15),tmvec(1:statenumber+15),&
         osc_str(1:statenumber+15),kpq,itmp,'i')

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

    integer               :: ndim,unit,itmp,i
    real(dp), allocatable :: ener(:)
    real(dp), allocatable :: vec(:)
    
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

  subroutine transform_operators(ndimf)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer                               :: ndimf,workdim,error,i,j
    real(dp), dimension(:,:), allocatable :: eigvec,tmp
    real(dp), dimension(:), allocatable   :: ener,work

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(eigvec(ndimf+1,ndimf+1))
    eigvec=0.0d0

    allocate(ener(ndimf+1))
    ener=0.0d0

    workdim=3*(ndimf+1)
    allocate(work(workdim))
    work=0.0d0

    allocate(tmp(ndimf+1,ndimf+1))
    tmp=0.0d0
    
!-----------------------------------------------------------------------
! Diagonalise the field-free Hamiltonian matrix
!-----------------------------------------------------------------------
    ! Diagonalisation of the ADC(1)/CIS Hamiltonian matrix padded with
    ! zeros corresponding to the H_J0 elements
    eigvec=0.0d0
    eigvec(1:ndimf,1:ndimf)=h1
    call dsyev('V','U',ndimf+1,eigvec,ndimf+1,ener,work,workdim,error)

    ! Exit if the diagonalisation failed
    if (error.ne.0) then
       errmsg='Diagonalisation of the ADC(1)/CIS Hamiltonian failed &
            in subroutine transform_operators'
    endif

!-----------------------------------------------------------------------
! To be consistent with the TD-ADC(2) code, move the ground state vector
! to the last basis vector position
!-----------------------------------------------------------------------
    do i=1,ndimf
       eigvec(:,i)=eigvec(:,i+1)
    enddo
    eigvec(:,ndimf+1)=0.0d0
    eigvec(ndimf+1,ndimf+1)=1.0d0

!-----------------------------------------------------------------------
! Transform the operator matrices to the eigenstate representation
!-----------------------------------------------------------------------
    ! Field-free Hamiltonian
    tmp=matmul(transpose(eigvec),matmul(h1,eigvec))
    h1=tmp

    ! Dipole operator
    do i=1,3
       if (sum(abs(pulse_vec(i,1:npulse))).eq.0.0d0) cycle
       tmp=matmul(transpose(eigvec),matmul(d1(i,:,:),eigvec))
       d1(i,:,:)=tmp
    enddo
    
    ! CAP
    if (lcap) then
       tmp=matmul(transpose(eigvec),matmul(w1,eigvec))
       w1=tmp
    endif

    ! CAP projector
    if (lflux) then
       tmp=matmul(transpose(eigvec),matmul(theta1,eigvec))
       theta1=tmp
    endif
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(eigvec)
    deallocate(ener)
    deallocate(work)
    deallocate(tmp)
    
    return
    
  end subroutine transform_operators
  
!#######################################################################
  
  subroutine cap_projection(ndimf)

    use constants
    use parameters

    implicit none

    integer :: ndimf

    if (tdrep.eq.1) then
       call cap_projection_isr(ndimf)
    else if (tdrep.eq.2) then
       call cap_projection_eigen(ndimf)
    endif

    return
    
  end subroutine cap_projection

!#######################################################################

  subroutine cap_projection_isr(ndimf)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer                               :: ndimf,workdim,error,i,j
    real(dp), dimension(:,:), allocatable :: eigvec,tmp
    real(dp), dimension(:), allocatable   :: ener,work

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(eigvec(ndimf+1,ndimf+1))
    eigvec=0.0d0

    allocate(ener(ndimf+1))
    ener=0.0d0

    workdim=3*(ndimf+1)
    allocate(work(workdim))
    work=0.0d0

    allocate(tmp(ndimf+1,ndimf+1))
    tmp=0.0d0
    
!----------------------------------------------------------------------
! Diagonalise the field-free Hamiltonian matrix
!----------------------------------------------------------------------
    ! Diagonalisation of the ADC(1)/CIS Hamiltonian matrix padded with
    ! zeros corresponding to the H_J0 elements
    eigvec=0.0d0
    eigvec(1:ndimf,1:ndimf)=h1
    call dsyev('V','U',ndimf+1,eigvec,ndimf+1,ener,work,workdim,error)

    ! Exit if the diagonalisation failed
    if (error.ne.0) then
       errmsg='Diagonalisation of the ADC(1)/CIS Hamiltonian failed &
            in subroutine cap_projection_isr'
    endif

!-----------------------------------------------------------------------
! Move the ground state vector to the last basis vector position
!-----------------------------------------------------------------------
    do i=1,ndimf
       eigvec(:,i)=eigvec(:,i+1)
       ener(i)=ener(i+1)
    enddo
    eigvec(:,ndimf+1)=0.0d0
    eigvec(ndimf+1,ndimf+1)=1.0d0
    ener(ndimf+1)=0.0d0
    
!-----------------------------------------------------------------------
! Projection of the CAP matrix
!-----------------------------------------------------------------------
    ! Rotate the CAP matrix to the eigenvector representation
    tmp=matmul(transpose(eigvec),matmul(w1,eigvec))

    ! Projection
    if (iprojcap.eq.1) then
       ! Initial state projection
       if (statenumber.eq.0) then
          tmp(ndimf+1,:)=0.0d0
          tmp(:,ndimf+1)=0.0d0
       else
          tmp(statenumber,:)=0.0d0
          tmp(:,statenumber)=0.0d0
       endif
    else if (iprojcap.eq.2) then
       ! Bound state projection
       do i=1,ndimf+1
          if (ener(i).le.projlim) then
             tmp(i,:)=0.0d0
             tmp(:,i)=0.0d0
          endif
       enddo
    else if (iprojcap.eq.3) then
       ! Annihilation of bound-bound elements
       do i=1,ndimf+1
          if (ener(i).gt.projlim) cycle
          do j=1,ndimf+1
             if (ener(j).gt.projlim) cycle
             tmp(i,j)=0.0d0
          enddo
       enddo
    else if (iprojcap.eq.4) then
       ! Annihilation of bound-unbound and unbound-bound elements
       do i=1,ndimf+1
          do j=1,ndimf+1
             if (ener(i).le.projlim.and.ener(j).le.projlim) cycle
             if (ener(i).gt.projlim.and.ener(j).gt.projlim) cycle
             tmp(i,j)=0.0d0
          enddo
       enddo
    endif

    ! Rotate the CAP matrix back to the ISR representation
    w1=matmul(eigvec,matmul(tmp,transpose(eigvec)))

!-----------------------------------------------------------------------
! Projection of the CAP-projector matrix
!-----------------------------------------------------------------------
    ! Rotate the CAP-projector matrix to the eigenvector representation
    tmp=matmul(transpose(eigvec),matmul(theta1,eigvec))

    ! Projection
    if (iprojcap.eq.1) then
       ! Initial state projection
       if (statenumber.eq.0) then
          tmp(ndimf+1,:)=0.0d0
          tmp(:,ndimf+1)=0.0d0
       else
          tmp(statenumber,:)=0.0d0
          tmp(:,statenumber)=0.0d0
       endif
    else if (iprojcap.eq.2) then
       ! Bound state projection
       do i=1,ndimf+1
          if (ener(i).le.projlim) then
             tmp(i,:)=0.0d0
             tmp(:,i)=0.0d0
          endif
       enddo
    else if (iprojcap.eq.3) then
       ! Annihilation of bound-bound elements
       do i=1,ndimf+1
          if (ener(i).gt.projlim) cycle
          do j=1,ndimf+1
             if (ener(j).gt.projlim) cycle
             tmp(i,j)=0.0d0
          enddo
       enddo
    else if (iprojcap.eq.4) then
       ! Annihilation of bound-unbound and unbound-bound elements
       do i=1,ndimf+1
          do j=1,ndimf+1
             if (ener(i).le.projlim.and.ener(j).le.projlim) cycle
             if (ener(i).gt.projlim.and.ener(j).gt.projlim) cycle
             tmp(i,j)=0.0d0
          enddo
       enddo
    endif

    ! Rotate the CAP-projector matrix back to the ISR representation
    theta1=matmul(eigvec,matmul(tmp,transpose(eigvec)))
 
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    deallocate(eigvec)
    deallocate(ener)
    deallocate(work)
    deallocate(tmp)
    
    return
    
  end subroutine cap_projection_isr

!#######################################################################

  subroutine cap_projection_eigen(ndimf)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer               :: ndimf
    integer               :: i,j
    real(dp), allocatable :: ener(:)
    real(dp), allocatable :: vec(:)

!----------------------------------------------------------------------
! Initial state projection
!----------------------------------------------------------------------
    if (iprojcap.eq.1) then
       if (statenumber.eq.0) then
          w1(ndimf+1,:)=0.0d0
          w1(:,ndimf+1)=0.0d0
          theta1(ndimf+1,:)=0.0d0
          theta1(:,ndimf+1)=0.0d0
       else
          w1(statenumber,:)=0.0d0
          w1(:,statenumber)=0.0d0
          theta1(statenumber,:)=0.0d0
          theta1(:,statenumber)=0.0d0
       endif
    endif

!----------------------------------------------------------------------
! Bound state projection
!----------------------------------------------------------------------
    if (iprojcap.eq.2) then
       do i=1,ndimf+1
          if (h1(i,i).le.projlim) then
             w1(i,:)=0.0d0
             w1(:,i)=0.0d0
             theta1(i,:)=0.0d0
             theta1(:,i)=0.0d0
          endif
       enddo
    endif

!----------------------------------------------------------------------
! Annihilation of bound-bound elements
!----------------------------------------------------------------------
    if (iprojcap.eq.3) then
       do i=1,ndimf+1
          if (h1(i,i).gt.projlim) cycle
          do j=1,ndimf+1
             if (h1(j,j).gt.projlim) cycle
             w1(i,j)=0.0d0
             theta1(i,j)=0.0d0
          enddo
       enddo
    endif

!----------------------------------------------------------------------
! Annihilation of bound-unbound and unbound-bound elements
!----------------------------------------------------------------------
    if (iprojcap.eq.4) then
       do i=1,ndimf+1
          do j=1,ndimf+1
             if (h1(i,i).le.projlim.and.h1(j,j).le.projlim) cycle
             if (h1(i,i).gt.projlim.and.h1(j,j).gt.projlim) cycle
             w1(i,j)=0.0d0
             theta1(i,j)=0.0d0
          enddo
       enddo
    endif
       
    return
    
  end subroutine cap_projection_eigen

!#######################################################################

  subroutine diag_hcap_adc1(ndimf)

    use constants
    use parameters
    use iomod
    use misc, only: dsortindxa1
    
    implicit none

    integer                                  :: ndimf,lwork,ierr,i
    integer, dimension(:), allocatable       :: indx
    real(dp), dimension(:), allocatable      :: rwork
    complex(dp), dimension(:,:), allocatable :: capham,vecr,vecl
    complex(dp), dimension(:), allocatable   :: lambda,work
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(capham(ndimf+1,ndimf+1))
    capham=czero

    allocate(lambda(ndimf+1))
    lambda=czero

    allocate(vecr(ndimf+1,ndimf+1))
    vecr=czero

    allocate(vecl(ndimf+1,ndimf+1))
    vecl=czero

    lwork=5*(ndimf+1)
    allocate(work(lwork))
    work=czero

    allocate(rwork(2*(ndimf+1)))
    rwork=0.0d0

    allocate(indx(ndimf+1))
    indx=0
    
!----------------------------------------------------------------------
! Compute eigenvalues of H-iW
!----------------------------------------------------------------------
    capham=h1-ci*w1
    call zgeev('V','V',ndimf+1,capham,ndimf+1,lambda,vecl,ndimf+1,&
         vecr,ndimf+1,work,lwork,rwork,ierr)

    if (ierr.ne.0) then
       errmsg='Diagonalisation of H-iW failed in subroutine &
            diag_hcap_adc1'
       call error_control
    endif

!----------------------------------------------------------------------
! Print the eigenvalues of H-iW for checking purposes
!----------------------------------------------------------------------
    call dsortindxa1('A',ndimf+1,abs(lambda),indx)

    write(ilog,'(/,72a)') ('-',i=1,72)
    write(ilog,'(2x,a)') 'Eigenvalues of H-iW (eV)'
    write(ilog,'(72a)') ('-',i=1,72)
    do i=1,min(statenumber+50,ndimf+1)
       write(ilog,'(2x,i2,2x,F10.7,x,a,x,F11.8,a)') &
            i,real(lambda(indx(i)))*27.2113845d0,&
            '+',aimag(lambda(indx(i)))*27.2113845d0,'*i'
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(capham)
    deallocate(lambda)
    deallocate(vecr)
    deallocate(vecl)
    deallocate(work)
    deallocate(rwork)
    deallocate(indx)
    
    return
    
  end subroutine diag_hcap_adc1
    
!#######################################################################

end module adc1propmod
