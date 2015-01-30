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

  end module guessvecs
