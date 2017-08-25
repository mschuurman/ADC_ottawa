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
! Set up the initial wavefunction vector
!----------------------------------------------------------------------
    call initialise(ndimf,noffdf)

!----------------------------------------------------------------------
! Determine what can be held in memory
!----------------------------------------------------------------------
    call memory_managment

!----------------------------------------------------------------------
! Output some information about the calculation to be performed
!----------------------------------------------------------------------
    call wrinfo

!----------------------------------------------------------------------
! Loading of the non-zero elements of the Hamiltonian matrix into
! memory
!----------------------------------------------------------------------
    if (hincore) call load_hamiltonian('SCRATCH/hmlt.diac',&
         'SCRATCH/hmlt.offc',ndimf,noffdf)
    
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

    use tdsemod
    
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

    write(ilog,'(2x,a,x,ES15.8)') 'Error tolerance:',proptol

!----------------------------------------------------------------------
! Matrix-vector multiplication algorithm
!----------------------------------------------------------------------
    if (hincore) then
       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
            will proceed in-core'
    else
       write(ilog,'(/,2x,a,/)') 'Matrix-vector multiplication &
            will proceed out-of-core'
    endif
    
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
    
    implicit none

    if (lcap) then
       ! Propagation using the short iterative Lanczos-Arnoldi
       ! algorithm
       errmsg='Propagation using a non-Hermitian Hamiltonian is &
            not yet supported'
       call error_control
       ! call propagate_wavepacket_csil
    else
       ! Propagation using the short iterative Lanczos algorithm
       call propagate_wavepacket_sil
    endif
       
    return
    
  end subroutine propagate_wavepacket

!#######################################################################

  subroutine propagate_wavepacket_sil

    use tdsemod
    use sillib
    use csillib
    
    implicit none

    integer                                 :: i
    real(d)                                 :: norm
    real(d), parameter                      :: tiny=1e-9_d
    complex(d), dimension(:), allocatable   :: dtpsi
    
    ! SIL arrays and variables
    integer                                 :: steps,trueorder,&
                                               errorcode
    real(d)                                 :: intperiod,stepsize,&
                                               truestepsize,time,&
                                               inttime
    real(d), dimension(:,:), allocatable    :: eigenvector
    real(d), dimension(:), allocatable      :: diagonal,eigenval
    real(d), dimension(:), allocatable      :: offdiag
    real(d), dimension(:), allocatable      :: offdg2    
    complex(d), dimension(:,:), allocatable :: krylov
    logical(kind=4)                         :: restart,relax,stdform

!----------------------------------------------------------------------
! sillib variables
!----------------------------------------------------------------------
    ! This is not a relaxation calculation
    relax=.false.

    ! func <-> -iH|Psi>
    stdform=.true.

    ! No. steps taken
    steps=0

    ! Time
    time=0.0d0

    ! Restart flag - if true, the Krylov space is built up
    ! before propagation, else the old Krylov vectors are used.       
    restart=.true.
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Wavefunction arrays
    allocate(dtpsi(matdim))
    dtpsi=czero

    ! sillib arrays
    allocate(krylov(matdim,kdim-1))
    krylov=czero

    allocate(eigenvector(kdim+1,kdim+3))
    eigenvector=0.0d0

    allocate(eigenval(kdim+1))
    eigenval=0.0d0

    allocate(diagonal(kdim+1))
    diagonal=0.0d0

    allocate(offdg2(kdim+1))
    offdg2=0.0d0

    allocate(offdiag(kdim))
    offdiag=0.0d0

!----------------------------------------------------------------------
! Propagate forwards in time with the light-matter interaction
! Hamiltonian, H - mu.E(t), using the Quantics sillib libraries
!----------------------------------------------------------------------
    ! Integration period
    intperiod=tout

    ! Loop over the timesteps
    do i=1,int(tfinal/intperiod)

       ! Time at the start of the current timestep
       time=(i-1)*intperiod

       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime

       ! dtpsi = -iH(t)|Psi>
       call matxvec_treal_laser(time,matdim,noffdiag,psi,dtpsi)
           
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dtpsi)
    deallocate(krylov)
    deallocate(eigenvector)
    deallocate(eigenval)
    deallocate(diagonal)
    deallocate(offdg2)
    deallocate(offdiag)
    
    return
    
  end subroutine propagate_wavepacket_sil
    
!#######################################################################

  subroutine finalise

    use tdsemod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(psi)
    if (hincore) call deallocate_hamiltonian
    
    return
    
  end subroutine finalise
    
!#######################################################################
  
end module propagate
