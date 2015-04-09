  subroutine master_adc2_ener

    use constants
    use parameters
    use select_fano
    use davmod
    use band_lanczos
    use fspace
    use get_moment
    use misc
    use fspace2
    use get_matrix
    use get_matrix_DIPOLE
    use propagate_prepare
    use guessvecs

    implicit none

    integer, dimension(:,:), allocatable :: kpq
    integer                              :: ndim,ndims,i,itmp,ista
    integer*8                            :: noffd
    real(d), dimension(:), allocatable   :: ener,vec_init,mtm,tmvec,&
                                            osc_str
    real(d), dimension(:,:), allocatable :: rvec
    real(d)                              :: t1,t2

!-----------------------------------------------------------------------
! Allocate kpq and select configurations
!-----------------------------------------------------------------------
    allocate(kpq(7,0:nBas**2*4*nOcc**2))
    kpq(:,:)=-1

    if (lcvs) then
       ! CVS-ADC(2)
       call select_atom_is_cvs(kpq(:,:))
       call select_atom_d_cvs(kpq(:,:),-1)
    else
       ! ADC(2)-s
       call select_atom_is(kpq(:,:))
       call select_atom_d(kpq(:,:),-1)
    endif

!-----------------------------------------------------------------------
! Output supspace dimensions
!-----------------------------------------------------------------------
    ndims=kpq(1,0)
    ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)

    write(ilog,*) 'ADC(2) INITIAL Space dim',ndim

    write(ilog,*) 'dimension of various INITIAL configuration spaces'
    write(ilog,*) '      1p1h       2p2h_1      2p2h_2      2p2h_3      2p2h_4i      2p2h_4ii'
    write(ilog,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0),kpq(5,0)
    write(ilog,*)

!-----------------------------------------------------------------------
! Initialise dipole moment matrix
!-----------------------------------------------------------------------
    if (tranmom2 .eq. 'x') then
       dpl(:,:)=x_dipole(:,:)
    elseif (tranmom2 .eq. 'y') then
       dpl(:,:)=y_dipole(:,:)
    elseif (tranmom2 .eq. 'z') then
       dpl(:,:)=z_dipole(:,:)
    end if
    
    CHECK_dip=nirrep2

!-----------------------------------------------------------------------
! If we are performing a fake ip calculation, then determine the
! indices of the 1h1p configurations corresponding to excitation into
! the additional diffuse 'fake continuum' orbital
!-----------------------------------------------------------------------
    itmp=1+nBas**2*4*nOcc**2
    if (lfakeip) call get_fakeip_indices(kpq,itmp,ndims,ndim)

!-----------------------------------------------------------------------
! Calculate and save the Hamiltonian matrix to file
!-----------------------------------------------------------------------
    write(ilog,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
    if (lcvs) then
       call write_fspace_adc2_1_cvs(ndim,kpq(:,:),noffd,'i')
    else
       call write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i')
    endif

!-----------------------------------------------------------------------
! Block-Davidson diagonalisation of the Hamiltonian matrix
!-----------------------------------------------------------------------
    allocate(ener(davstates),rvec(ndim,davstates))
    allocate(vec_init(ndim))

    call master_dav(ndim,noffd,'i',ndims)
    
    call readdavvc(davstates,ener,rvec)

!-----------------------------------------------------------------------
! Calculate TDMs from the ground state
!-----------------------------------------------------------------------    
    allocate(mtm(ndim),tmvec(davstates),osc_str(davstates))
    
    tmvec=0.0d0
    osc_str=0.0d0

    if (.not.lfakeip) then
       call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:),1)
       do i=1,davstates
          tmvec(i)=tm(ndim,rvec(:,i),mtm(:))
          osc_str(i)=2._d/3._d*ener(i)*tmvec(i)**2
       end do
    endif

    itmp=1+nBas**2*4*nOcc**2
    call table2(ndim,davstates,ener(1:davstates),rvec(:,1:davstates),&
         tmvec(1:davstates),osc_str(1:davstates),kpq,itmp)
    
    return

  end subroutine master_adc2_ener
