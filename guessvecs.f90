  module guessvecs
    
    implicit none

  contains

!#######################################################################

    subroutine adc1_guessvecs

      use constants
      use parameters
      use select_fano
      use fspace
      use fspace2

      implicit none

      integer                              :: ndim,i,iout
      integer, dimension(:,:), allocatable :: kpq
      real(d), dimension(:,:), allocatable :: eigvec
      real(d), dimension(:), allocatable   :: eigval

!-----------------------------------------------------------------------
! Select initial 1h1p space
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))
      kpq(:,:)=-1
      call select_atom_is(kpq(:,:))

!-----------------------------------------------------------------------
! Output initial 1h1p space information
!-----------------------------------------------------------------------
      ndim=kpq(1,0)
      write(6,*) 'ADC(1) INITIAL Space dim',ndim
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(eigvec(ndim,ndim))
      allocate(eigval(ndim))

!-----------------------------------------------------------------------
! Diagonalise the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------
      write(6,'(/,78a)') ('-',i=1,78)
      write(6,'(2x,a)') 'Generating guess Davidson vectors by diagonalising &
           the ADC(1) Hamiltonian'
      write(6,'(78a,/)') ('-',i=1,78)
      call get_fspace_tda_direct(ndim,kpq(:,:),eigvec,eigval)

!-----------------------------------------------------------------------
! Write the ADC(1) eigenvectors to file
!-----------------------------------------------------------------------
      iout=22
      open(iout,file='adc1_vecs',form='unformatted',status='unknown')
      write(iout) ndim,eigvec
      close(iout)

      return

    end subroutine adc1_guessvecs

!#######################################################################

    subroutine get_fakeip_indices(kpq,kpqdim2,dims)
      
      use parameters, only: ifakeorb,nocc,ifakeex,dmain,ncore,lcvs
      use misc, only: get_indices

      implicit none
      
      integer                            :: kpqdim2,dims,inda,indb,indj,&
                                             indk,spin,i,nex,k,norb
      integer, dimension(7,0:kpqdim2-1)  :: kpq
      integer, dimension(nocc,2)         :: exindx

!-----------------------------------------------------------------------
! Determine the indices of the 1h1p configurations corresponding to
! excitation into the diffuse orbital
!-----------------------------------------------------------------------
      exindx=0
      nex=0
      do i=1,dims
         call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
         if (inda.eq.ifakeorb) then
            nex=nex+1
            exindx(nex,1)=i
            exindx(nex,2)=indj
         endif         
      enddo

!-----------------------------------------------------------------------
! Fill in the ifakeex array
!
! Here we take the 1h1p configurations in the order:
!
! homo   -> inf.
! homo-1 -> inf.
! .
! .
! .
! homo-N+1 -> inf.
!
! Ideally we should change this such that we start using 2h2p
! configurations as guesses when dmain>N_elec
! (or dmain>ncore if the CVS approximation is being used)
!-----------------------------------------------------------------------
      allocate(ifakeex(dmain))

      if (lcvs) then
         norb=ncore
         if (dmain.gt.ncore) then
            write(6,'(/,2(2x,a,/))') 'Currently only 1h1p initial vectors &
                 are supported for a fake IP calculation.','Reduce &
                 dmain so that it is less than or equal to ncore and &
                 re-run'
            STOP
         endif
      else
         norb=nocc         
         if (dmain.gt.nocc) then
            write(6,'(/,2(2x,a,/))') 'Currently only 1h1p initial vectors &
                 are supported for a fake IP calculation.','Reduce &
                 dmain so that it is less than or equal to nocc and &
                 re-run'
            STOP
         endif
      endif

      do i=1,dmain
         k=norb-i+1
         ifakeex(i)=exindx(k,1)
      enddo

      return

    end subroutine get_fakeip_indices

!#######################################################################

  end module guessvecs
