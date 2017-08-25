!#######################################################################
! molecule_laser: Routines to perform wavepacket propagations
!                 including the interaction of the molecule with a
!                 laser pulse
!#######################################################################

module propagate

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                               :: matdim
  integer*8                             :: noffdiag
  complex(d), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser: Wavepacket propagation including the molecule-laser
!                  interaction.
!                  The wavepacket is represented in a basis consisting
!                  of the intermediate state basis plus the MP2 ground
!                  state.
!########################################################################
! IMPORTANT:       For ease of implementation of the Hamiltonian and
!                  dipole matrix-vector products, the MP2 ground state
!                  is taken to be the last basis function in the set.
!#######################################################################
  
  subroutine propagate_laser(ndimf,noffdf)

    use tdsemod
    
    implicit none

    integer, intent(in)   :: ndimf
    integer*8, intent(in) :: noffdf
    integer               :: k
    real(d)               :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Output some information about the calculation to be performed
!----------------------------------------------------------------------
    call wrinfo

!----------------------------------------------------------------------
! Set up the initial wavefunction vector
!----------------------------------------------------------------------
    call initialise(ndimf,noffdf)

!----------------------------------------------------------------------
! Determine what can be held in memory
!----------------------------------------------------------------------
    call memory_managment

    print*,hincore

    STOP
    
!----------------------------------------------------------------------
! Peform the wavepacket propagation
!----------------------------------------------------------------------
    call propagate_wavepacket
    
!----------------------------------------------------------------------
! Finalise and deallocate arrays
!----------------------------------------------------------------------
    call finalise
    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine propagate_laser

!#######################################################################

  subroutine wrinfo

    implicit none

    integer :: k

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Wavepacket propagation'
    write(ilog,'(72a)') ('-',k=1,72)

!----------------------------------------------------------------------
! SIL parameters
!----------------------------------------------------------------------         
    if (lcap) then
       write(ilog,'(2x,a,/)') 'Wavepacket propagation performed &
            using the CSIL method'
    else
       write(ilog,'(2x,a,/)') 'Wavepacket propagation performed &
            using the SIL method'
    endif
       
    write(ilog,'(2x,a,x,i2,/)') 'Maximum Krylov subspace dimension:',&
         kdim

!    write(ilog,'(2x,a,x,ES15.8)') 'Error tolerance:',autotol

!!----------------------------------------------------------------------
!! Matrix-vector multiplication algorithm
!!----------------------------------------------------------------------
!    if (hincore) then
!       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
!            will proceed in-core'
!    else
!       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
!            will proceed out-of-core'
!    endif
    
    return
    
  end subroutine wrinfo

!#######################################################################

  subroutine initialise(ndimf,noffdf)

    implicit none

    integer   :: ndimf
    integer*8 :: noffdf

!----------------------------------------------------------------------
! No. non-zero off-diagonal matrix elements
!----------------------------------------------------------------------
    noffdiag=noffdf
    
!----------------------------------------------------------------------
! Wavepacket dimension: IS basis plus the ground state
!----------------------------------------------------------------------
    matdim=ndimf+1

!----------------------------------------------------------------------
! Allocate the wavepacket array
!----------------------------------------------------------------------
    allocate(psi(matdim))
    psi=0.0d0
    
!----------------------------------------------------------------------
! Set the initial wavepacket
!
! For now, we will only support excitation/ionisation from the ground
! state, corresponding to the vector (0,...,0,1)^T
!----------------------------------------------------------------------
    psi=0.0d0
    psi(matdim)=cone
    
    return
    
  end subroutine initialise

!#######################################################################

  subroutine memory_managment

    use tdsemod
    use omp_lib
    
    implicit none

    integer*8 :: maxrecl,reqmem
    integer   :: nthreads
    real(d)   :: memavail
    
!----------------------------------------------------------------------
! Available memory
!----------------------------------------------------------------------
    ! Maximum memory requested to be used by the user
    memavail=maxmem

    ! Two-electron integrals held in-core
    memavail=memavail-8.0d0*(nbas**4)/1024.0d0**2

    ! kpq
    memavail=memavail-8.0d0*7.0d0*(1+nbas**2*4*nocc**2)/1024.0d0**2
    
    ! Psi
    memavail=memavail-8.0d0*matdim/1024.0d0**2
    
    ! Lanczos vectors used in the SIL propagation method
    memavail=memavail-(kdim-1)*8.0d0*matdim/1024.0d0**2

    ! Be cautious and only use say 90% of the available memory
    memavail=memavail*0.9d0 

!----------------------------------------------------------------------
! Determine whether or not we can hold the non-zero Hamiltonian
! matrix elements in-core
!----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

    reqmem=0.0d0
    
    ! Parallelised matrix-vector multiplication
    reqmem=reqmem+8.0d0*nthreads*matdim/1024.0d0**2

    ! Non-zero off-diagonal Hamiltonian matrix elements and their
    ! indices
    reqmem=reqmem+8.0d0*2.0d0*noffdiag/1024.0d0**2

    ! On-diagonal Hamiltonian matrix elements
    reqmem=reqmem+8.0d0*matdim/1024.0d0**2

    ! Set the hincore flag controling whether the matrix-vector
    ! multiplication proceeds in-core
    if (reqmem.lt.memavail) then
       hincore=.true.
    else
       hincore=.false.
    endif
    
    return

  end subroutine memory_managment
    
!#######################################################################

  subroutine propagate_wavepacket

    use sillib
    use csillib
    
    implicit none

    
    
    return
    
  end subroutine propagate_wavepacket
    
!#######################################################################

  subroutine finalise

    implicit none

    deallocate(psi)
    
    return
    
  end subroutine finalise
    
!#######################################################################
  
end module propagate
