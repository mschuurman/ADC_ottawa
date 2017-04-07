!#######################################################################
! Calculation of the electronic absorption spectrum from the Fourier
! transform of the autocorrelation function obtained by propagation
! of |Psi(t=0)> = mu |Psi_0>
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
    real(d), dimension(:), allocatable   :: fvec
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
! Calculate the f-vector, F_J = < Psi_J | D | Psi_0 >    
!-----------------------------------------------------------------------
    call fvector(fvec,ndimf,kpqf)

!-----------------------------------------------------------------------
! Calculate the Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian(ndimf,kpqf,noffdf)
    
!-----------------------------------------------------------------------
! Perform the wavepacket propagation and autocorrelation function
! calculation
!-----------------------------------------------------------------------
    hamflag='f'
    call propagate_fvec(fvec,ndimf)

!-----------------------------------------------------------------------
! Output the autocorrelation function
!-----------------------------------------------------------------------
    call wrauto

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(fvec)
    deallocate(auto)

    return
      
  end subroutine adc2_autospec

!#######################################################################

  subroutine fvector(fvec,ndimf,kpqf)

    use constants
    use parameters
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    real(d), dimension(:), allocatable        :: fvec
    
!-----------------------------------------------------------------------
! Allocate the F-vector array
!-----------------------------------------------------------------------
    allocate(fvec(ndimf))

!-----------------------------------------------------------------------
! Calculate the F-vector, F_J = < Psi_J | D | Psi_0 >
!-----------------------------------------------------------------------
    call get_modifiedtm_adc2(ndimf,kpqf(:,:),fvec(:),1)
    
    return
    
  end subroutine fvector
    
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
    
  subroutine wrauto

    use constants
    use parameters
    use iomod

    implicit none

    integer :: iauto,i
    real(d) :: t

!----------------------------------------------------------------------
! Open output file
!----------------------------------------------------------------------
    call freeunit(iauto)
    open(iauto,file='auto',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the autocorrelation function to file
!----------------------------------------------------------------------
    write(iauto,'(a)') '#    time[fs]         Re(autocorrel)     &
         Im(autocorrel)     Abs(autocorrel)'

    do i=1,nstep
       t=(i-1)*dt/41.341375d0
       write(iauto,'(F15.8,4x,3(2x,F17.14))') &
            t,real(auto(i)),aimag(auto(i)),abs(auto(i)) 
    enddo

!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(iauto)

    return

  end subroutine wrauto

!#######################################################################
  
end module adc2automod
