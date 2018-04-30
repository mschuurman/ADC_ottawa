!######################################################################
! fvecprop: Routines for the calculation of wavepacket autocorrelation
!           functions using an initial state given by the f-vector
!           f_J = < Psi_J | D | Psi_0 >
!######################################################################

module fvecprop

  use constants
  use parameters
  use channels
  use iomod
  use timingmod
  
  implicit none

  integer                                :: matdim,iout0,iout1,iout2
  integer*8                              :: noffdiag
  real(dp), parameter                    :: au2fs=1.0d0/41.34137333656d0
  real(dp)                               :: norm0
  complex(dp), dimension(:), allocatable :: psi0

contains

!######################################################################

  subroutine propagate_fvec(fvec,ndimf,noffdf)

    use tdsemod
    
    implicit none

    integer, intent(in)                    :: ndimf
    integer*8, intent(in)                  :: noffdf
    integer                                :: k
    real(dp), dimension(ndimf), intent(in) :: fvec
    real(dp)                               :: tw1,tw2,tc1,tc2
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Initialisation and allocatation
!----------------------------------------------------------------------
    call initialise(ndimf,noffdf)

!----------------------------------------------------------------------
! Determine what can be held in memory
!----------------------------------------------------------------------
    call memory_managment

!----------------------------------------------------------------------
! Output where we are at and what we are doing
!----------------------------------------------------------------------
    call wrinfo
    
!----------------------------------------------------------------------
! Normalise the F-vector to form |Psi(t=0)>
!----------------------------------------------------------------------
    call init_wavepacket(fvec)

!----------------------------------------------------------------------
! Loading of the non-zero elements of the Hamiltonian matrix into
! memory
!----------------------------------------------------------------------
    if (hincore) call load_hamiltonian('SCRATCH/hmlt.diac',&
         'SCRATCH/hmlt.offc',matdim,noffdf)
    
!----------------------------------------------------------------------
! Open the autocorrelation function output files and write the file
! headers
!----------------------------------------------------------------------
    call open_autofiles

!----------------------------------------------------------------------
! Propagation and calculation of the autocorrelation functions
!----------------------------------------------------------------------
    call propagate_sillib
           
!----------------------------------------------------------------------
! Close the autocorrelation function output file
!----------------------------------------------------------------------
    call close_autofiles

!----------------------------------------------------------------------
! Output timings and the no. matrix-vector multiplications
!----------------------------------------------------------------------
    write(ilog,'(70a)') ('+',k=1,70)

    write(ilog,'(/,a,1x,i5)') 'No. matrix-vector multiplications:',&
           nmult

    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
!----------------------------------------------------------------------
! Finalisation and deallocation
!----------------------------------------------------------------------
    call finalise
    
    return
    
  end subroutine propagate_fvec

!######################################################################

  subroutine wrinfo

    use tdsemod
    
    implicit none

    integer :: k

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of the nth-order &
         autocorrelation functions:'
    write(ilog,'(70a)') ('-',k=1,70)
    write(ilog,'(4x,a)') 'a_n(t) = <Psi(0)| H^n |Psi(t)>, &
         |Psi(0)> = D|Psi_0>'
    write(ilog,'(70a,/)') ('-',k=1,70)

!----------------------------------------------------------------------
! Maximum autocorrelation function order
!----------------------------------------------------------------------
    write(ilog,'(2x,a,x,i1,x,a,/)') &
         'Autocorrelation functions up to order',autoord,&
         'will be calculated'

!----------------------------------------------------------------------
! Wavepacket propagation information
!----------------------------------------------------------------------         
    write(ilog,'(2x,a,/)') 'Wavepacket propagation performed &
         using the SIL method'
    write(ilog,'(2x,a,x,i2,/)') 'Maximum Krylov subspace &
         dimension:',kdim
    write(ilog,'(2x,a,x,ES15.8)') 'Error tolerance:',autotol
       
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

!----------------------------------------------------------------------
! Initialisation of the matrix-vector multiplication counter
!----------------------------------------------------------------------
    nmult=0
    
    return

  end subroutine wrinfo

!######################################################################

  subroutine initialise(ndimf,noffdf)

    implicit none

    integer, intent(in)   :: ndimf
    integer*8, intent(in) :: noffdf
    
    ! Hamiltonian matrix dimension
    matdim=ndimf

    ! No. non-zero off-diagonal matrix elements
    noffdiag=noffdf
    
    ! Psi(t=0)
    allocate(psi0(matdim))
    psi0=czero
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine memory_managment

    use tdsemod
    use omp_lib
    
    implicit none

    integer*8 :: maxrecl,reqmem
    integer   :: nthreads
    real(dp)  :: memavail

!----------------------------------------------------------------------
! Available memory
!----------------------------------------------------------------------
    ! Maximum memory requested to be used by the user
    memavail=maxmem

    ! Two-electron integrals held in-core
    memavail=memavail-8.0d0*(nbas**4)/1024.0d0**2

    ! kpq
    memavail=memavail-8.0d0*7.0d0*(1+nbas**2*4*nocc**2)/1024.0d0**2
    
    ! Psi(0) and Psi(t)
    memavail=memavail-2.0d0*16.0d0*matdim/1024.0d0**2
    
    ! Lanczos vectors used in the SIL propagation method
    memavail=memavail-(kdim-1)*16.0d0*matdim/1024.0d0**2

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
    
!######################################################################

  subroutine finalise

    use tdsemod
    
    implicit none
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(psi0)
    if (hincore) call deallocate_hamiltonian
    
    return
    
  end subroutine finalise

!######################################################################
    
  subroutine open_autofiles

    implicit none

!----------------------------------------------------------------------
! Zeroth-order autocorrelation function file
!----------------------------------------------------------------------
    call freeunit(iout0)
    open(iout0,file='auto',form='formatted',status='unknown')
    write(iout0,'(a)') '#    time[fs]         Re(autocorrel)         &
         Im(autocorrel)         Abs(autocorrel)'

!----------------------------------------------------------------------
! First-order autocorrelation function file
!----------------------------------------------------------------------
    if (autoord.ge.1) then
       call freeunit(iout1)
       open(iout1,file='auto1',form='formatted',status='unknown')
       write(iout1,'(a)') &
            '#    time[fs]         Re(autocorrel)         &
            Im(autocorrel)         Abs(autocorrel)'
    endif
    
!----------------------------------------------------------------------
! Second-order autocorrelation function file
!----------------------------------------------------------------------
    if (autoord.eq.2) then
       call freeunit(iout2)
       open(iout2,file='auto2',form='formatted',status='unknown')
       write(iout2,'(a)') &
            '#    time[fs]         Re(autocorrel)         &
            Im(autocorrel)         Abs(autocorrel)'
    endif

    return

  end subroutine open_autofiles

!######################################################################

  subroutine close_autofiles

    implicit none

!----------------------------------------------------------------------
! Zero-order autocorrelation function file
!----------------------------------------------------------------------
    close(iout0)

!----------------------------------------------------------------------
! First-order autocorrelation function file
!----------------------------------------------------------------------
    if (autoord.ge.1) close(iout1)

!----------------------------------------------------------------------
! Second-order autocorrelation function file
!----------------------------------------------------------------------
    if (autoord.eq.2) close(iout2)

    return

  end subroutine close_autofiles

!######################################################################

  subroutine init_wavepacket(fvec)

    implicit none

    integer                                 :: i
    real(dp), dimension(matdim), intent(in) :: fvec
    
!----------------------------------------------------------------------
! The initial wavepacket is taken as D|Psi_0>/||D|Psi_0>||
!----------------------------------------------------------------------
    do i=1,matdim
       psi0(i)=dcmplx(fvec(i),0.0d0)
    enddo

    norm0=sqrt(dot_product(psi0,psi0))
    psi0=psi0/norm0

    return
    
  end subroutine init_wavepacket

!######################################################################
! propagate_sillib: wavepacket propagation using the short iterative
!                   Lanczos method to yield the time-domain wavepacket
!                   autocorrelation functions a_n(t), n=0,1,2
!######################################################################
  
  subroutine propagate_sillib

    use sillib
    use tdsemod
    
    implicit none

    integer                                :: i
    real(dp)                               :: norm
    real(dp), parameter                    :: tiny=1e-9_dp
    complex(dp), dimension(:), allocatable :: psi,dtpsi
    complex(dp), dimension(:), allocatable :: hpsi,h2psi
    complex(dp)                            :: auto0,auto1,auto2
    
    ! SIL arrays and variables
    integer                                  :: steps,trueorder,&
                                                errorcode
    real(dp)                                 :: intperiod,stepsize,&
                                                truestepsize,time,&
                                                inttime
    real(dp), dimension(:,:), allocatable    :: eigenvector
    real(dp), dimension(:), allocatable      :: diagonal,eigenval
    real(dp), dimension(:), allocatable      :: offdiag
    real(dp), dimension(:), allocatable      :: offdg2    
    complex(dp), dimension(:,:), allocatable :: krylov
    logical(kind=4)                          :: restart,relax,stdform
    
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
    allocate(psi(matdim))
    psi=czero

    allocate(dtpsi(matdim))
    dtpsi=czero

    if (autoord.ge.1) then
       allocate(hpsi(matdim))
       hpsi=czero
    endif

    if (autoord.eq.2) then
       allocate(h2psi(matdim))
       h2psi=czero
    endif

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
! Output the autocorrelation functions at time t=0
!----------------------------------------------------------------------
    ! Zeroth-order autocorrelation function at time t=0
    auto0=dot_product(psi0,psi0)*norm0**2
    call wrauto(iout0,auto0,0.0d0)
    
    ! First-order autocorrelation functions at time t=0
    if (autoord.ge.1) then
       call matxvec_treal(time,matdim,noffdiag,psi0,hpsi)
       hpsi=hpsi*ci
       auto1=dot_product(psi0,hpsi)*norm0**2
       call wrauto(iout1,auto1,0.0d0)
    endif
    
    ! Second-order autocorrelation functions at time t=0
    if (autoord.eq.2) then
       call matxvec_treal(time,matdim,noffdiag,hpsi,h2psi)
       h2psi=h2psi*ci
       auto2=dot_product(psi0,h2psi)*norm0**2
       call wrauto(iout2,auto2,0.0d0)
    endif

!----------------------------------------------------------------------
! Propagate |Psi(t=0)> = D |Psi_0> forwards in time using the Quantics
! sillib libraries
!
! Note that we are using t/2 trick here in calculating the
! autocorrelation functions
!----------------------------------------------------------------------
    ! Integration period
    intperiod=0.5d0*tout

    ! Set the initial wavepacket
    psi=psi0

    ! Loop over the timesteps at which we want the autocorrelation
    ! function
    do i=1,int(tfinal/intperiod)
       
       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime
       
       ! dtpsi = -iH|Psi>
       call matxvec_treal(time,matdim,noffdiag,psi,dtpsi)
    
       ! Take one step using the SIL algorithm
       call silstep(psi,dtpsi,matdim,noffdiag,stepsize,kdim,autotol,relax,&
            restart,stdform,steps,krylov,truestepsize,trueorder,&
            errorcode,time,matxvec_treal,eigenvector,eigenval,&
            diagonal,offdg2,offdiag)

       ! Exit if the SIL integration failed
       if (errorcode.ne.0) then
          call silerrormsg (errorcode,errmsg)
          call error_control
       endif

       ! Updtate the propagation time
       time=time+truestepsize
       
       ! Check whether the integration is complete
       inttime=inttime+truestepsize
       if (abs(intperiod-inttime).gt.abs(tiny*intperiod)) goto 100

       ! Calculate and output the zeroth-order autocorrelation 
       ! function at the current timestep
       auto0=dot_product(conjg(psi),psi)*norm0**2
       call wrauto(iout0,auto0,i*intperiod*2.0d0)

       ! Calculate and output the first-order autocorrelation
       ! at the current timestep
       if (autoord.ge.1) then
          call matxvec_treal(time,matdim,noffdiag,psi,hpsi)
          hpsi=hpsi*ci
          auto1=dot_product(conjg(psi),hpsi)*norm0**2
          call wrauto(iout1,auto1,i*intperiod*2.0d0)
       endif

       ! Calculate and output the second-order autocorrelation
       ! at the current timestep
       if (autoord.eq.2) then
          call matxvec_treal(time,matdim,noffdiag,hpsi,h2psi)
          h2psi=h2psi*ci
          auto2=dot_product(conjg(psi),h2psi)*norm0**2
          call wrauto(iout2,auto2,i*intperiod*2.0d0)
       endif

       ! Output some information about our progress
       norm=real(sqrt(dot_product(psi,psi)))
       call wrprogress(i*intperiod,norm)

    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(psi)
    deallocate(dtpsi)
    deallocate(krylov)
    deallocate(eigenvector)
    deallocate(eigenval)
    deallocate(diagonal)
    deallocate(offdg2)
    deallocate(offdiag)
    if (allocated(hpsi)) deallocate(hpsi)
    if (allocated(h2psi)) deallocate(h2psi)

    return
    
  end subroutine propagate_sillib

!######################################################################
  
  subroutine wrprogress(t,norm)

    implicit none

    integer  :: k
    real(dp) :: t,norm

    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(a,x,F10.4)') 'Time:',t
    write(ilog,'(a,x,F8.6)') 'Norm',norm
    
    return
    
  end subroutine wrprogress

!######################################################################

  subroutine wrauto(iout,auto,t)

    implicit none

    integer     :: iout
    real(dp)    :: t
    complex(dp) :: auto

    write(iout,'(F15.8,4x,3(2x,E21.14))') &
         t*au2fs,real(auto),aimag(auto),abs(auto)

    return
    
  end subroutine wrauto

!######################################################################
    
end module fvecprop
