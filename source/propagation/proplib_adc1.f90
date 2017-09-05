!#######################################################################
! propagate_adc1: Routines to perform ADC(1) wavepacket propagations
!                 including the interaction of the molecule with a
!                 laser pulse
!#######################################################################

module propagate_adc1

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                               :: matdim
  complex(d), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser_adc1: ADC(1) wavepacket propagation including the
!                       molecule-laser interaction.
!                       The wavepacket is represented in a basis
!                       consisting of the intermediate state basis plus
!                       the HF ground state.
!########################################################################
! IMPORTANT:            For ease of implementation of the Hamiltonian
!                       and dipole matrix-vector products, the HF ground
!                       state is taken to be the last basis function in
!                       the set.
!#######################################################################
  
  subroutine propagate_laser_adc1(ndimf,kpqf)

    use tdsemod
    
    implicit none

    integer, intent(in)                       :: ndimf
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(d)                                   :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Set up the initial wavefunction vector
!----------------------------------------------------------------------
    call initialise(ndimf)

!----------------------------------------------------------------------
! Output some information about the calculation to be performed
!----------------------------------------------------------------------
    call wrinfo

!----------------------------------------------------------------------
! Peform the wavepacket propagation
!----------------------------------------------------------------------
    call propagate_wavepacket(kpqf)
    
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
    
  end subroutine propagate_laser_adc1

!#######################################################################

    subroutine initialise(ndimf)

    implicit none

    integer :: ndimf

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

    write(ilog,'(2x,a,x,ES15.8)') 'Error tolerance:',proptol

    return
    
  end subroutine wrinfo

!#######################################################################

    subroutine propagate_wavepacket(kpqf)
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf

    if (lcap) then
       ! Propagation using the short iterative Lanczos-Arnoldi
       ! algorithm
        call propagate_wavepacket_csil(kpqf)
    else
       ! Propagation using the short iterative Lanczos algorithm
       !call propagate_wavepacket_sil(kpqf)
       print*,"FINISH WRITING THE ADC(1) PROPAGATION CODE!"       
    endif

    return
    
  end subroutine propagate_wavepacket

!#######################################################################

  !######################################################################

  subroutine propagate_wavepacket_csil(kpqf)

    use tdsemod
    use csillib
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: i
    integer*8                                 :: dummy
    real(d)                                   :: norm
    real(d), parameter                        :: tiny=1e-9_d
    complex(d), dimension(:), allocatable     :: dtpsi,hpsi
    
    ! CSIL arrays and variables
    integer                                   :: steps,trueorder,&
                                                 errorcode
    real(d)                                   :: intperiod,stepsize,&
                                                 truestepsize,time,&
                                                 inttime,macheps
    real(d), dimension(:,:), allocatable      :: eigenvector
    real(d), dimension(:), allocatable        :: diagonal,eigenval
    real(d), dimension(:), allocatable        :: offdiag
    real(d), dimension(:), allocatable        :: offdg2
    complex(d), dimension(:,:), allocatable   :: krylov
    complex(d), dimension(0:kdim,0:kdim)      :: hessenberg,eigvec,&
                                                 auxmat
    logical(kind=4)                           :: restart,relax,&
                                                 stdform,olderrcri

!----------------------------------------------------------------------
! Machine epsilon
!----------------------------------------------------------------------
    macheps=epsilon(macheps)

!----------------------------------------------------------------------
! csillib variables
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

    ! Error criterion - if false, the improved error criterion is used
    olderrcri=.false.

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    ! Wavefunction arrays
    allocate(dtpsi(matdim))
    dtpsi=czero

    allocate(hpsi(matdim))
    hpsi=czero
    
    ! sillib arrays
    allocate(krylov(matdim,kdim-1))
    krylov=czero   

    hessenberg=czero

!----------------------------------------------------------------------
! Dummy variable: The sil/csil libraries need to be passed the number
! of non-zero off-diagonal Hamiltonian matrix elements for use in an
! ADC(2) propagation, for which only the non-zero elements are stored.
! For the ADC(1) propagation being performed, we are storing all
! matrices in full in-core, but still need to pass some dummy argument.
!----------------------------------------------------------------------
    dummy=0
    
!----------------------------------------------------------------------
! Propagate forwards in time with the CAP-augmented light-matter
! interaction Hamiltonian, H - mu.E(t) - i*W_CAP, using the Quantics
! csillib libraries
!----------------------------------------------------------------------
    ! Integration period
    intperiod=tout

    ! Loop over the timesteps
    do i=1,int(tfinal/intperiod)

       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! NOTE THAT IN THIS MODULE MATDIM = NO. ISs + 1
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       ! dtpsi = -iH(t)|Psi>
       call matxvec_treal_laser_adc1(time,matdim,dummy,psi,dtpsi)

       ! Take one step using the SIL algorithm
       call csilstep(psi,dtpsi,matdim,dummy,stepsize,kdim,&
            proptol,relax,restart,stdform,olderrcri,steps,&
            truestepsize,trueorder,errorcode,time,macheps,&
            matxvec_treal_laser_adc1,hessenberg,eigvec,krylov,&
            auxmat)
                     
       ! Exit if the CSIL integration failed
       if (errorcode.ne.0) then
          call csilerrormsg(errorcode,errmsg)
          call error_control
       endif

       ! Updtate the propagation time
       time=time+truestepsize

       ! Check whether the integration is complete
       inttime=inttime+truestepsize
       if (abs(intperiod-inttime).gt.abs(tiny*intperiod)) goto 100

       ! Output some information about our progress
       norm=real(sqrt(dot_product(psi,psi)))
       !call wrstepinfo(time,norm,kpqf)

       print*,time,norm
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dtpsi)
    deallocate(hpsi)
    deallocate(krylov)
    
    return
    
  end subroutine propagate_wavepacket_csil
    
!#######################################################################

    subroutine finalise

    use tdsemod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(psi)

    return

  end subroutine finalise
    
!#######################################################################
  
end module propagate_adc1
