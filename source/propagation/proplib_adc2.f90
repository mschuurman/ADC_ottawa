!#######################################################################
! propagate_adc2: Routines to perform ADC(2) wavepacket propagations
!                 including the interaction of the molecule with a
!                 laser pulse
!#######################################################################

module propagate_adc2

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                               :: matdim
  integer*8                             :: noffdiag
  integer                               :: iflux
  complex(d), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser_adc2: ADC(2) wavepacket propagation including the
!                       molecule-laser interaction.
!                       The wavepacket is represented in a basis
!                       consisting of the intermediate state basis plus
!                       the MP2 ground state.
!########################################################################
! IMPORTANT:       For ease of implementation of the Hamiltonian and
!                  dipole matrix-vector products, the MP2 ground state
!                  is taken to be the last basis function in the set.
!#######################################################################
  
  subroutine propagate_laser_adc2(ndimf,noffdf,kpqf)

    use tdsemod
    
    implicit none

    integer, intent(in)                       :: ndimf
    integer*8, intent(in)                     :: noffdf
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(d)                                   :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Open the output files
!----------------------------------------------------------------------
    call open_outfiles
    
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
    call propagate_wavepacket(kpqf)
    
!----------------------------------------------------------------------
! Finalise and deallocate arrays
!----------------------------------------------------------------------
    call finalise

!----------------------------------------------------------------------
! Close the output files
!----------------------------------------------------------------------
    call close_outfiles

!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine propagate_laser_adc2

!#######################################################################

  subroutine open_outfiles

    use iomod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Flux expectation values
!----------------------------------------------------------------------
    if (lflux) then
       call freeunit(iflux)
       open(iflux,file='flux.dat',form='formatted',status='unknown')
       write(iflux,'(72a)') ('#',i=1,72)
       write(iflux,'(a)') '# Time (au)    Flux'
       write(iflux,'(72a)') ('#',i=1,72)
    endif
       
    return
    
  end subroutine open_outfiles

!#######################################################################

  subroutine close_outfiles

    implicit none

!----------------------------------------------------------------------
! Flux expectation values
!----------------------------------------------------------------------
    if (lflux) close(iflux)
       
    return
    
  end subroutine close_outfiles
  
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

  subroutine propagate_wavepacket(kpqf)
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf

    if (lcap) then
       ! Propagation using the short iterative Lanczos-Arnoldi
       ! algorithm
        call propagate_wavepacket_csil(kpqf)
    else
       ! Propagation using the short iterative Lanczos algorithm
       call propagate_wavepacket_sil(kpqf)
    endif
       
    return
    
  end subroutine propagate_wavepacket

!#######################################################################

  subroutine propagate_wavepacket_sil(kpqf)

    use tdsemod
    use sillib
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: i
    real(d)                                   :: norm
    real(d), parameter                        :: tiny=1e-9_d
    real(d), parameter                        :: tinier=1e-10_d
    complex(d), dimension(:), allocatable     :: dtpsi,hpsi
    
    ! SIL arrays and variables
    integer                                   :: steps,trueorder,&
                                                 errorcode
    real(d)                                   :: intperiod,stepsize,&
                                                 truestepsize,time,&
                                                 inttime
    real(d), dimension(:,:), allocatable      :: eigenvector
    real(d), dimension(:), allocatable        :: diagonal,eigenval
    real(d), dimension(:), allocatable        :: offdiag
    real(d), dimension(:), allocatable        :: offdg2    
    complex(d), dimension(:,:), allocatable   :: krylov
    logical(kind=4)                           :: restart,relax,stdform

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

    allocate(hpsi(matdim))
    hpsi=czero
    
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

       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! NOTE THAT IN THIS MODULE MATDIM = NO. ISs + 1
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       ! dtpsi = -iH(t)|Psi>
       call matxvec_treal_laser(time,matdim,noffdiag,psi,dtpsi)

       ! Integrate the TDSE if || dtpsi || is non-zero
       if (real(sqrt(dot_product(dtpsi,dtpsi))).gt.tinier) then
       
          ! Take one step using the SIL algorithm
          call silstep(psi,dtpsi,matdim,noffdiag,stepsize,kdim,proptol,relax,&
               restart,stdform,steps,krylov,truestepsize,trueorder,&
               errorcode,time,matxvec_treal_laser,eigenvector,eigenval,&
               diagonal,offdg2,offdiag)
          
          ! Exit if the SIL integration failed
          if (errorcode.ne.0) then
             call silerrormsg(errorcode,errmsg)
             call error_control
          endif
          
          ! Updtate the propagation time
          time=time+truestepsize
          
          ! Check whether the integration is complete
          inttime=inttime+truestepsize
          if (abs(intperiod-inttime).gt.abs(tiny*intperiod)) goto 100

       else

          ! If the time derivative of the wavefunction is almost zero
          ! then skip the integration for this timestep
          time=time+intperiod
          
       endif
          
       ! Output some information about our progress
       norm=real(sqrt(dot_product(psi,psi)))
       call wrstepinfo(time,norm,kpqf)
       
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

!######################################################################

  subroutine propagate_wavepacket_csil(kpqf)

    use tdsemod
    use csillib
    use fluxmod
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: i
    real(d)                                   :: norm,flux
    real(d), parameter                        :: tiny=1e-9_d
    real(d), parameter                        :: tinier=1e-10_d
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
       call matxvec_treal_laser(time,matdim,noffdiag,psi,dtpsi)

       ! Flux analysis (we only need to do this at the start of
       ! each timestep)
       if (lflux.and.inttime.eq.0.0d0) then
          call adc2_flux_cap(matdim,psi,dtpsi,flux)
          write(iflux,'(F10.4,5x,ES15.8)') time,flux
       endif
       
       ! Integrate the TDSE if || dtpsi || is non-zero
       if (real(sqrt(dot_product(dtpsi,dtpsi))).gt.tinier) then
       
          ! Take one step using the SIL algorithm
          call csilstep(psi,dtpsi,matdim,noffdiag,stepsize,kdim,&
               proptol,relax,restart,stdform,olderrcri,steps,&
               truestepsize,trueorder,errorcode,time,macheps,&
               matxvec_treal_laser,hessenberg,eigvec,krylov,&
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

       else

          ! If the time derivative of the wavefunction is almost zero
          ! then skip the integration for this timestep
          time=time+intperiod
          
       endif
          
       ! Output some information about our progress
       norm=real(sqrt(dot_product(psi,psi)))
       call wrstepinfo(time,norm,kpqf)
       
    enddo

    ! Final flux expectation value
    if (lflux) then
       call matxvec_treal_laser(time,matdim,noffdiag,psi,dtpsi)
       call adc2_flux_cap(matdim,psi,dtpsi,flux)
       write(iflux,'(F10.4,5x,ES15.8)') time,flux
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dtpsi)
    deallocate(hpsi)
    deallocate(krylov)
    
    return
    
  end subroutine propagate_wavepacket_csil
    
!######################################################################

  subroutine wrstepinfo(t,norm,kpqf)

    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(d)                                   :: t,norm

    write(ilog,'(70a)') ('+',k=1,70)

    ! Propagation time
    write(ilog,'(/,2x,a,7x,F10.4)') 'Time:',t

    ! Wavefunction norm
    write(ilog,'(2x,a,11x,F6.4)') 'Norm:',norm

    ! 'Number of electrons'
    write(ilog,'(2x,a,3x,F8.4)') 'Norm x Nel:',norm*nocc*2.0d0
    
    ! Wavefunction analysis
    call wrpsi(kpqf)
    
    return
    
  end subroutine wrstepinfo

!#######################################################################

  subroutine wrpsi(kpqf)

    use misc, only: dsortindxa1,getspincase
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer, dimension(:), allocatable        :: indx
    integer                                   :: k,ilbl
    integer                                   :: kpqdim2
    real(d), dimension(:), allocatable        :: abscoeff
    real(d), parameter                        :: coefftol=0.01d0
    character(len=2)                          :: spincase

    kpqdim2=nbas**2*4*nocc**2+1
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(abscoeff(matdim))
    allocate(indx(matdim))

!-----------------------------------------------------------------------
! Sort the coefficients by magnitude
!-----------------------------------------------------------------------
    abscoeff=abs(psi)
    call dsortindxa1('D',matdim,abscoeff,indx)

!-----------------------------------------------------------------------
! Output the configurations contributing significantly to the
! wavepacket
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a,/)') 'Dominant Configurations:'
    write(ilog,'(2x,30a)') ('*',k=1,30)
    write(ilog,'(3x,a)') 'j   k -> a  b       |C_jkab|'
    write(ilog,'(2x,30a)') ('*',k=1,30)

    ! Ground state contribution
    write(ilog,'(3x,a,15x,F8.5)') 'Psi0',abs(psi(matdim))

    ! IS basis functions
    do k=1,50

       ilbl=indx(k)

       ! Skip the ground state
       if (ilbl.eq.matdim) cycle

       ! Skip if the coefficient is small
       if (abs(psi(ilbl)).lt.coefftol) cycle

       if (kpqf(4,ilbl).eq.-1) then
          ! Single excitations
          write(ilog,'(3x,i2,4x,a2,1x,i2,8x,F8.5)') &
               kpqf(3,ilbl),'->',kpqf(5,ilbl),abs(psi(ilbl))
       else
          ! Double excitations
          if (kpqf(3,ilbl).ne.kpqf(4,ilbl)&
               .and.kpqf(5,ilbl).ne.kpqf(6,ilbl)) then
             ! a|=b, i|=j
             spincase=getspincase(ilbl,kpqf,kpqdim2)
             write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),2x,a2,1x,F8.5)') &
                  kpqf(3,ilbl),kpqf(4,ilbl),'->',kpqf(5,ilbl),&
                  kpqf(6,ilbl),spincase,abs(psi(ilbl))
          else
             ! a=b,  i=j
             ! a|=b, i=j
             ! a=b,  i=|j
             write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),5x,F8.5)') &
                  kpqf(3,ilbl),kpqf(4,ilbl),'->',kpqf(5,ilbl),&
                  kpqf(6,ilbl),abs(psi(ilbl))
          endif
       endif

    enddo

    write(ilog,'(2x,30a,/)') ('*',k=1,30)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(abscoeff)
    deallocate(indx)
    
    return
    
  end subroutine wrpsi
    
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
  
end module propagate_adc2
