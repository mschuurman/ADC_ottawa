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

        implicit none        

        integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
        integer                              :: i,ndim,ndims,ndimsf,&
                                                nout,ndimf,ndimd,&
                                                noutf,itmp
        integer*8                            :: noffd,noffdf
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
! Determine the 1h1p and 2h2p subspaces
!-----------------------------------------------------------------------
        allocate(kpq(7,0:nBas**2*4*nOcc**2))
        allocate(kpqf(7,0:nBas**2*4*nOcc**2))
        allocate(kpqd(7,0:nBas**2*4*nOcc**2))

        kpq(:,:)=-1
        call select_atom_is(kpq(:,:))

        kpqf(:,:)=-1
        call select_atom_isf(kpqf(:,:))
        
        kpqd(:,:)=-1
        call select_atom_ist(kpqd(:,:))

        ndim  = kpq(1,0)
        ndimf = kpqf(1,0)
        ndimd = kpqd(1,0)

        write(ilog,*) 'ADC(1) INITIAL Space dim',ndim
        write(ilog,*) 'ADC(1) FINAL Space dim',ndimf
        write(ilog,*) 'ADC(1) TOTAL Space dim WITHOUT GROUND STATE',ndimd
        write(ilog,*) 'dimension of various INITIAL configuration spaces'
        write(ilog,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0)
        write(ilog,*) 'dimension of various FINAL configuration spaces'
        write(ilog,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
        write(ilog,*) 'dimension of various TOTAL configuration spaces'
        write(ilog,*) kpqd(1,0),kpqd(2,0),kpqd(3,0),kpqd(4,0),kpqd(5,0)

!-----------------------------------------------------------------------
! Set the dipole matrix
!-----------------------------------------------------------------------
        if (tranmom2 .eq. 'x') then
           dpl(:,:)=x_dipole(:,:)
        elseif (tranmom2 .eq. 'y') then
           dpl(:,:)=y_dipole(:,:)
        elseif (tranmom2 .eq. 'z') then
           dpl(:,:)=z_dipole(:,:)
        endif

!-----------------------------------------------------------------------
! Set the irrep of the dipole operator
!-----------------------------------------------------------------------
        CHECK_dip = nirrep2

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
        call get_modifiedtm_tda(ndim,kpq(:,:),mtm(:))

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
