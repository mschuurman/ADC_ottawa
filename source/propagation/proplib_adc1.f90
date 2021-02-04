!#######################################################################
! propagate_adc1: Routines to perform ADC(1) and CIS wavepacket
!                 propagations including the interaction of the
!                 molecule with a laser pulse.
!#######################################################################

module propagate_adc1

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  save

  integer                                :: matdim
  integer                                :: iflux
  complex(dp), dimension(:), allocatable :: psi
  
contains

!#######################################################################
! propagate_laser_adc1: ADC(1) wavepacket propagation including the
!                       molecule-laser interaction.
!                       The wavepacket is represented in a basis
!                       consisting of the intermediate state basis plus
!                       the HF ground state.
!########################################################################
! IMPORTANT:            To be consistent with the TD-ADC(2) code, the
!                       HF ground state is taken to be the last basis
!                       function in the set.
!#######################################################################
  
  subroutine propagate_laser_adc1(ndimf,kpqf)

    use tdsemod
    
    implicit none

    integer, intent(in)                       :: ndimf
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(dp)                                  :: tw1,tw2,tc1,tc2

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
    call initialise(ndimf)

!----------------------------------------------------------------------
! Set up the initial wavefunction vector
!----------------------------------------------------------------------
    call initwf
    
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
    
  end subroutine propagate_laser_adc1

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
    psi=czero
    
    return
    
  end subroutine initialise

!#######################################################################

    subroutine initwf

    implicit none
    
!----------------------------------------------------------------------
! Fill in the initial wavefunction vector
!----------------------------------------------------------------------
    if (lmixistate) then
       ! Mixed initial state
       if (tdrep.eq.1) then
          ! ISR basis
          call initwf_mixed_isr
       else
          ! Eigenstate basis
          call initwf_mixed_eigen
       endif
    else if (statenumber.eq.0) then
       ! Ground state: the MP2 ground state is included in the basis
       ! as the ndimf'th + 1 IS basis function, i.e., the initial
       ! wavefunction vector is given by (0,...,0,1)^T
       psi=czero
       psi(matdim)=cone
    else
       ! Excited state
       if (tdrep.eq.1) then
          ! ISR basis: read the initial wavefunction vector from disk
          call initwf_exci_initstate
       else
          ! Eigenstate basis
          psi=czero
          psi(statenumber)=cone
       endif
    endif
       
    return
    
  end subroutine initwf

!#######################################################################

  subroutine initwf_mixed_isr

    use misc, only: dsortindxa1
    
    implicit none

    integer               :: i,j,klbl,unit,itmp
    integer, allocatable  :: indx(:)
    real(dp)              :: ftmp
    real(dp), allocatable :: vec(:)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(indx(nmix))
    allocate(vec(matdim-1))

!-----------------------------------------------------------------------
! Sort the state indices in the expansion in order of increasing value
!-----------------------------------------------------------------------
    call dsortindxa1('A',nmix,dble(imix),indx)

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
    psi=czero
    
!-----------------------------------------------------------------------
! Ground state contribution
!-----------------------------------------------------------------------
    if (imix(indx(1)).eq.0) then
       psi(matdim)=cmplx(cmix(indx(1)),0.0d0)
       klbl=2
    else
       klbl=1
    endif

!-----------------------------------------------------------------------    
! Excited state contributions
!-----------------------------------------------------------------------
    ! Open the file containing the initial space eigenpairs
    call freeunit(unit)
    open(unit,file='SCRATCH/initvecs',status='old',&
         access='sequential',form='unformatted')

    ! Read the initial state vectors from disk
    do i=1,imix(indx(nmix))
       read(unit) itmp,ftmp,vec(1:matdim-1)
       if (i.eq.imix(indx(klbl))) then
          do j=1,matdim-1
             psi(j)=psi(j)+cmix(indx(klbl))*cmplx(vec(j),0.0d0)
          enddo
          klbl=klbl+1
       endif
    enddo
    
    ! Close the file containing the initial space eigenpairs
    close(unit)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(indx)
    deallocate(vec)
    
    return
    
  end subroutine initwf_mixed_isr

!#######################################################################

  subroutine initwf_mixed_eigen

    use misc, only: dsortindxa1
    
    implicit none

    integer               :: klbl,i
    integer, allocatable  :: indx(:)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(indx(nmix))

!-----------------------------------------------------------------------
! Sort the state indices in the expansion in order of increasing value
!-----------------------------------------------------------------------
    call dsortindxa1('A',nmix,dble(imix),indx)

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
    psi=czero

!-----------------------------------------------------------------------
! Construct the initial wavepacket
!-----------------------------------------------------------------------
    do i=1,nmix
       if (imix(indx(i)).eq.0) then
          ! Ground state
          psi(matdim)=cmplx(cmix(indx(i)),0.0d0)
       else
          ! Excited state
          psi(imix(indx(i)))=cmplx(cmix(indx(i)),0.0d0)
       endif
    enddo
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(indx)
    
    return
    
  end subroutine initwf_mixed_eigen
  
!#######################################################################

  subroutine initwf_exci_initstate

    implicit none
    
    integer               :: unit,itmp,i
    real(dp)              :: ftmp
    real(dp), allocatable :: vec(:)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(vec(matdim-1))

!-----------------------------------------------------------------------
! Open the file containing the initial space eigenpairs
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='SCRATCH/initvecs',status='old',&
         access='sequential',form='unformatted')

!-----------------------------------------------------------------------
! Read the initial state vector from disk
!-----------------------------------------------------------------------
    do i=1,statenumber
       read(unit) itmp,ftmp,vec(1:matdim-1)
    enddo

!-----------------------------------------------------------------------
! Set up the initial wavefunction vector
!-----------------------------------------------------------------------
    do i=1,matdim-1
       psi(i)=cmplx(vec(i),0.0d0)
    enddo
    psi(matdim)=czero

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(unit)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(vec)

    return
    
  end subroutine initwf_exci_initstate
  
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

!----------------------------------------------------------------------
! N-electron basis
!----------------------------------------------------------------------
    if (tdrep.eq.1) then
       write(ilog,'(/,2x,a)') 'The ISR basis is used'
    else if (tdrep.eq.2) then
       write(ilog,'(/,2x,a,/)') &
            'The field-free eigenstate basis is used'
    endif
    
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
    integer*8                                 :: dummy
    real(dp)                                  :: norm,flux
    real(dp), parameter                       :: tiny=1e-9_dp
    real(dp), parameter                       :: tinier=1e-10_dp
    complex(dp), dimension(:), allocatable    :: dtpsi,hpsi
    
    ! SIL arrays and variables
    integer                                   :: steps,trueorder,&
                                                 errorcode
    real(dp)                                  :: intperiod,stepsize,&
                                                 truestepsize,time,&
                                                 inttime
    real(dp), dimension(:,:), allocatable     :: eigenvector
    real(dp), dimension(:), allocatable       :: diagonal,eigenval
    real(dp), dimension(:), allocatable       :: offdiag
    real(dp), dimension(:), allocatable       :: offdg2    
    complex(dp), dimension(:,:), allocatable  :: krylov
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
! Dummy variable: The sil/csil libraries need to be passed the number
! of non-zero off-diagonal Hamiltonian matrix elements for use in an
! ADC(2) propagation, for which only the non-zero elements are stored.
! For the ADC(1) propagation being performed, we are storing all
! matrices in full in-core, but still need to pass some dummy argument.
!----------------------------------------------------------------------
    dummy=0.5*(matdim-1)*matdim
    
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
       call matxvec_treal_laser_adc1(time,matdim,dummy,psi,dtpsi)

       ! Output some information about our progress
       ! (this only needs to be done at the start of
       ! each timestep)
       if (inttime.eq.0.0d0) then
          norm=real(sqrt(dot_product(psi,psi)))
          call wrstepinfo(time,norm,flux,kpqf)
       endif
       
       ! Integrate the TDSE if || dtpsi || is non-zero
       if (real(sqrt(dot_product(dtpsi,dtpsi))).gt.tinier) then

          ! Take one step using the SIL algorithm
          call silstep(psi,dtpsi,matdim,dummy,stepsize,kdim,proptol,relax,&
               restart,stdform,steps,krylov,truestepsize,trueorder,&
               errorcode,time,matxvec_treal_laser_adc1,eigenvector,eigenval,&
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

    enddo

!----------------------------------------------------------------------
! Final report
!----------------------------------------------------------------------
    ! Final timestep output
    norm=real(sqrt(dot_product(psi,psi)))
    call wrstepinfo(time,norm,flux,kpqf)
    
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

  subroutine propagate_wavepacket_csil(kpqf)

    use tdsemod
    use csillib
    use fluxmod
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: i
    integer*8                                 :: dummy
    real(dp)                                  :: norm,flux
    real(dp), parameter                       :: tiny=1e-9_dp
    real(dp), parameter                       :: tinier=1e-10_dp
    complex(dp), dimension(:), allocatable    :: dtpsi,hpsi
    
    ! CSIL arrays and variables
    integer                                   :: steps,trueorder,&
                                                 errorcode
    real(dp)                                  :: intperiod,stepsize,&
                                                 truestepsize,time,&
                                                 inttime,macheps
    real(dp), dimension(:,:), allocatable     :: eigenvector
    real(dp), dimension(:), allocatable       :: diagonal,eigenval
    real(dp), dimension(:), allocatable       :: offdiag
    real(dp), dimension(:), allocatable       :: offdg2
    complex(dp), dimension(:,:), allocatable  :: krylov
    complex(dp), dimension(0:kdim,0:kdim)     :: hessenberg,eigvec,&
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
    dummy=0.5*(matdim-1)*matdim
    
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

       ! Output some information about our progress
       ! (this only needs to be done at the start of
       ! each timestep)
       if (inttime.eq.0.0d0) then
          norm=real(sqrt(dot_product(psi,psi)))
          if (lflux) call adc1_flux_cap(matdim,psi,dtpsi,flux,time)
          call wrstepinfo(time,norm,flux,kpqf)
       endif
       
       ! Integrate the TDSE if || dtpsi || is non-zero
       if (real(sqrt(dot_product(dtpsi,dtpsi))).gt.tinier) then

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

       else
          
          ! If the time derivative of the wavefunction is almost zero
          ! then skip the integration for this timestep
          time=time+intperiod

       endif

    enddo

!----------------------------------------------------------------------
! Final report
!----------------------------------------------------------------------
    ! Final flux expectation value
    if (lflux) then
       call matxvec_treal_laser_adc1(time,matdim,dummy,psi,dtpsi)
       call adc1_flux_cap(matdim,psi,dtpsi,flux,time)
    endif

    ! Final timestep output
    norm=real(sqrt(dot_product(psi,psi)))
    call wrstepinfo(time,norm,flux,kpqf)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dtpsi)
    deallocate(hpsi)
    deallocate(krylov)
    
    return
    
  end subroutine propagate_wavepacket_csil

!#######################################################################

  subroutine wrstepinfo(t,norm,flux,kpqf)

    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: k
    real(dp)                                  :: t,norm,flux

    write(ilog,'(70a)') ('+',k=1,70)

    ! Propagation time
    write(ilog,'(/,2x,a,7x,F12.6)') 'Time:',t

    ! Wavefunction norm
    write(ilog,'(2x,a,11x,F8.6)') 'Norm:',norm

    ! 'Number of electrons'
    write(ilog,'(2x,a,3x,F10.6)') 'Norm x Nel:',norm*nocc*2.0d0

    ! Flux
    if (lflux) write(iflux,'(F10.4,5x,ES15.8)') t,flux
    
    ! Wavefunction analysis
    if (tdrep.eq.1) then
       call wrpsi_isr(kpqf)
    else
       call wrpsi_eigen
    endif
       
    return
    
  end subroutine wrstepinfo

!#######################################################################

    subroutine wrpsi_isr(kpqf)

    use misc, only: dsortindxa1,getspincase
    
    implicit none

    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer, dimension(:), allocatable        :: indx
    integer                                   :: k,ilbl
    integer                                   :: kpqdim2
    real(dp), dimension(:), allocatable       :: abscoeff
    real(dp), parameter                       :: coefftol=0.01d0
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
     if (abs(psi(matdim)).gt.coefftol) &
          write(ilog,'(3x,a,15x,F8.5)') 'Psi0',abs(psi(matdim))

    ! IS basis functions
    do k=1,min(50,matdim-1)

       ilbl=indx(k)

       ! Skip the ground state
       if (ilbl.eq.matdim) cycle

       ! Skip if the coefficient is small
       if (abs(psi(ilbl)).lt.coefftol) cycle

       ! Single excitations
       write(ilog,'(3x,i2,4x,a2,1x,i2,8x,F8.5)') &
            kpqf(3,ilbl),'->',kpqf(5,ilbl),abs(psi(ilbl))

    enddo

    write(ilog,'(2x,30a,/)') ('*',k=1,30)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(abscoeff)
    deallocate(indx)
    
    return
    
  end subroutine wrpsi_isr

!#######################################################################

  subroutine wrpsi_eigen

    use misc, only: dsortindxa1,getspincase

    implicit none

    integer                             :: k,ilbl
    integer, dimension(:), allocatable  :: indx
    real(dp), dimension(:), allocatable :: abscoeff
    real(dp), parameter                 :: coefftol=0.001d0
    character(len=3)                    :: as
    character(len=5)                    :: aket
    
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
! Output the field-free eigenstates contributing significantly to the
! wavepacket
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a,/)') 'Dominant States:'

    write(ilog,'(2x,30a)') ('*',k=1,30)
    write(ilog,'(3x,a)') '|J>                 |C_J|'
    write(ilog,'(2x,30a)') ('*',k=1,30)

    ! Ground state contribution
    if (abs(psi(matdim)).gt.coefftol) &
         write(ilog,'(3x,a,15x,F8.5)') '|HF>',abs(psi(matdim))

    ! Excited state contributions
    do k=1,min(50,matdim-1)

       ilbl=indx(k)

       ! Skip the ground state
       if (ilbl.eq.matdim) cycle

       ! Skip if the coefficient is small
       if (abs(psi(ilbl)).lt.coefftol) cycle

       ! State and absolute coefficient value
       write(as,'(i3)') ilbl
       write(aket,'(a)') '|'//trim(adjustl(as))//'>'       
       write(ilog,'(3x,a,14x,F8.5)') aket,abs(psi(ilbl))

    enddo
    
    write(ilog,'(2x,30a,/)') ('*',k=1,30)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(abscoeff)
    deallocate(indx)
    
    return
    
  end subroutine wrpsi_eigen
    
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
