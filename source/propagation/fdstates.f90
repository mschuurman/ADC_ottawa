!######################################################################
! fdstates: Routines for the calculation of filter diagonalisation
!           states via the real-time propagation of the f-vector
!           f_J = < Psi_J | D | Psi_0 >
!######################################################################
module fdstates

  use constants
  use parameters
  use channels
  use iomod
  use timingmod

  implicit none

  integer                               :: matdim,fsunit
  integer*8                             :: noffdiag,buffsize,reclength
  real(d), dimension(:,:), allocatable  :: buffer
  real(d), dimension(:,:), allocatable  :: eigvec
  real(d), dimension(:), allocatable    :: eigval
  complex(d), dimension(:), allocatable :: psi0
  logical                               :: fsincore
  
contains

!######################################################################

  subroutine calc_fdstates(fvec,ndimf,noffdf)

    use tdsemod
    
    implicit none
    
    integer, intent(in)                   :: ndimf
    integer*8, intent(in)                 :: noffdf
    integer                               :: k
    real(d), dimension(ndimf), intent(in) :: fvec
    real(d)                               :: tw1,tw2,tc1,tc2

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
! Calculation of the filter states |Psi_Ej> via wavepacket propagation
!----------------------------------------------------------------------
    call calc_filterstates(noffdf)

!----------------------------------------------------------------------
! Calculation of the eigenstates of interest
!----------------------------------------------------------------------
    call calc_eigenstates
    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"

!----------------------------------------------------------------------
! Finalisation and deallocation
!----------------------------------------------------------------------
    call finalise
    
    return
    
  end subroutine calc_fdstates

!######################################################################

  subroutine wrinfo

    use tdsemod
    
    implicit none

    integer :: k

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of filter diagonalisation &
         eigenstates'
    write(ilog,'(70a,/)') ('-',k=1,70)

!----------------------------------------------------------------------
! SIL parameters
!----------------------------------------------------------------------         
    write(ilog,'(2x,a,/)') 'Wavepacket propagation performed using &
         the SIL method'
    write(ilog,'(2x,a,x,i2,/)') 'Maximum Krylov subspace dimension:',&
         kdim
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
    
    ! Psi(0) and Psi(t)
    memavail=memavail-2.0d0*8.0d0*matdim/1024.0d0**2
    
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

    ! If the non-zero Hamiltonian matrix elements are to be held
    ! in-core, then update memavail
    if (hincore) then
       memavail=memavail-8.0d0*nthreads*matdim/1024.0d0**2 &
            -8.0d0*2.0d0*noffdiag/1024.0d0**2 &
            -8.0d0*matdim/1024.0d0**2
    endif
    
!----------------------------------------------------------------------
! Determine whether or not we can store the filter states in-core
!----------------------------------------------------------------------
    ! Ammount of memory required for the in-core storage of the
    ! filter states
    reqmem=8.0d0*nfbas*matdim/1024.0d0**2

    ! Set the fsincore flag
    if (memavail.ge.reqmem) then
       fsincore=.true.
    else
       fsincore=.false.
    endif
    
!----------------------------------------------------------------------
! Determine the buffer size (equal to the no. filter states that can
! be held in-core) and the corresponding record length
!----------------------------------------------------------------------
    ! If all filter states can be held in-core, then set buffsize to
    ! the no. filter states, else set buffsize to the no. of filter
    ! states that can be held in-core
    if (fsincore) then
       buffsize=nfbas
    else
       buffsize=int(floor((memavail*1024.0d0**2)/(8.0d0*matdim)))
    endif
    
    ! Record length
    reclength=8*matdim*buffsize
    
    ! Make sure that the record length corresponding to the
    ! buffer size does not exceed the maximum value (2**31-1)
    maxrecl=2147483647
    if (reclength.gt.maxrecl) then
       buffsize=maxrecl/(8*matdim)
       reclength=8*matdim*buffsize
    endif

!----------------------------------------------------------------------
! Allocation of the buffer array and opening of the scratch file
! that will hold the filter states
!----------------------------------------------------------------------
    ! Allocate the buffer array
    allocate(buffer(matdim,buffsize))
    
    ! Open the scratch file
    if (.not.fsincore) then
       call freeunit(fsunit)
       open(fsunit,file='SCRATCH/fstates',form='unformatted',&
            status='unknown',access='direct',recl=reclength)
    endif
       
    return
    
  end subroutine memory_managment
    
!######################################################################

  subroutine init_wavepacket(fvec)

    implicit none

    integer                                :: i
    real(d), dimension(matdim), intent(in) :: fvec
    real(d)                                :: norm
    
!----------------------------------------------------------------------
! The initial wavepacket is taken as D|Psi_0>/||D|Psi_0>||
!----------------------------------------------------------------------
    do i=1,matdim
       psi0(i)=dcmplx(fvec(i),0.0d0)
    enddo

    norm=sqrt(dot_product(psi0,psi0))
    psi0=psi0/norm

    return
    
  end subroutine init_wavepacket

!######################################################################

  subroutine calc_filterstates(noffdf)

    use sillib
    use tdsemod
    
    implicit none

    integer*8, intent(in)                   :: noffdf
    integer                                 :: i
    real(d)                                 :: norm
    real(d), parameter                      :: tiny=1e-9_d
    complex(d), dimension(:), allocatable   :: psi,dtpsi
    
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
    allocate(psi(matdim))
    psi=czero

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
! Propagate |Psi(t=0)> = D |Psi_0> forwards in time using the Quantics
! sillib libraries and calculate the filter states on-the-fly
!----------------------------------------------------------------------
    ! Integration period
    intperiod=tout

    ! Set the initial wavepacket
    psi=psi0

    ! Initialise the buffer array that will hold the filter states
    buffer=0.0d0

    ! Contribution to the filter states from |Psi(t=0)>
    call update_filterstates(psi,intperiod,0)
    
    ! Loop over the timesteps at which we want the autocorrelation
    ! function
    do i=1,int(tfinal/intperiod)
       
       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime
       
       ! dtpsi = -iH|Psi>
       call matxvec_treal(matdim,noffdiag,psi,dtpsi)
    
       ! Take one step using the SIL algorithm
       call silstep(psi,dtpsi,matdim,noffdiag,stepsize,kdim,autotol,&
            relax,restart,stdform,steps,krylov,truestepsize,trueorder,&
            errorcode,time,matxvec_treal,eigenvector,eigenval,&
            diagonal,offdg2,offdiag)

       ! Exit if the SIL integration failed
       if (errorcode.ne.0) then
          call silerrormsg (errorcode,errmsg)
          call error_control
       endif

       ! Check whether the integration is complete
       inttime=inttime+truestepsize
       if (abs(intperiod-inttime).gt.abs(tiny*intperiod)) goto 100

       ! On-the-fly calculation of the filter states
       call update_filterstates(psi,intperiod,i)

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
    
    return
    
  end subroutine calc_filterstates

!######################################################################

  subroutine update_filterstates(psi,dt,istep)

    implicit none

    integer                       :: istep,j
    real(d)                       :: dt,fac,gk,t,ej,de
    complex(d), dimension(matdim) :: psi

!----------------------------------------------------------------------    
! We are performing the windowed Fourier transforms of the wavepacket
! using the trapezoidal rule.
! Accordingly, the contributions at t=0 or t=T are multiplied by
! 1/2 compared to those for other times.
!----------------------------------------------------------------------    
    if (istep.eq.0.or.istep.eq.int(tfinal/dt)) then
       fac=1.0d0
    else
       fac=2.0d0
    endif

!----------------------------------------------------------------------    
! Calculation of the contributions to the filter states from the
! current timestep
!----------------------------------------------------------------------    
    ! Current time
    t=istep*dt

    ! Window function value
    gk=windowfunc(t)

    ! Energy spacing
    de=(ebound(2)-ebound(1))/(nfbas-1)

    ! On-the-fly calculation of the filter states
    if (fsincore) then

       ! Loop over the filter states
       do j=1,nfbas
          
          ! Current energy
          ej=ebound(1)+(j-1)*de
          
          ! Contribution of the current timestep to the current filter
          ! state
          buffer(:,j)=buffer(:,j)+fac*dt*gk*real(exp(ci*ej*t)*psi)
          
       enddo

    else
       print*,"THE OUT-OF-CORE CODE NEEDS WRITING..."
       STOP
    endif
    
    return
    
  end subroutine update_filterstates

!######################################################################

  function windowfunc(t) result(gk)

    implicit none

    real(d) :: t,gk

    gk=cos((pi*t)/(2.0d0*tfinal))
    gk=gk**iwfunc
    
    return
    
  end function windowfunc
    
!######################################################################

  subroutine wrprogress(t,norm)

    implicit none
    
    integer :: istep,k
    real(d) :: t,norm

    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(a,x,F10.4)') 'Time:',t
    write(ilog,'(a,x,F8.6)') 'Norm',norm
    
    return
    
  end subroutine wrprogress

!######################################################################

  subroutine calc_eigenstates

    use tdsemod
    
    implicit none

    integer                            :: i,j,unit
    real(d)                            :: norm
    real(d), dimension(:), allocatable :: hpsi

!----------------------------------------------------------------------
! Calculation of the eigensates of interest
!----------------------------------------------------------------------
    if (fsincore) then
       !
       ! In-core calculation of the eigenstates of interest
       !

       ! Allocate arrays
       allocate(eigvec(matdim,nsel))
       eigvec=0.0d0
       allocate(eigval(nsel))
       eigval=0.0d0
       allocate(hpsi(matdim))
       hpsi=0.0d0
       
       ! Calculation of the eigenstates of interest
       do i=1,nsel
          do j=1,nfbas
             eigvec(:,i)=eigvec(:,i)+fbas2eig(isel(i),j)*buffer(:,j)
          enddo
       enddo
       
       ! Normalisation
       do i=1,nsel
          norm=sqrt(dot_product(eigvec(:,i),eigvec(:,i)))
          eigvec(:,i)=eigvec(:,i)/norm
       enddo
       
       ! Calculation of energies
       do i=1,nsel
          call matxvec(matdim,noffdiag,eigvec(:,i),hpsi)
          hpsi=-hpsi
          eigval(i)=dot_product(eigvec(:,i),hpsi)
       enddo
       
       ! Write the eigenstates to file
       call freeunit(unit)
       open(unit=unit,file='SCRATCH/fdstates',status='unknown',&
            access='sequential',form='unformatted')
       do i=1,nsel
          write(unit) i,eigval(i),eigvec(:,i)
       enddo
       close(unit)
       
       ! Deallocate arrays
       deallocate(eigvec)
       deallocate(eigval)
       deallocate(hpsi)
       
    else
       !
       ! Out-of-core calculation of the eigenstates
       !
       print*,"THE OUT-OF-CORE CODE NEEDS WRITING..."
       STOP

    endif
    
    return
    
  end subroutine calc_eigenstates
  
!######################################################################

  subroutine finalise

    use tdsemod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(psi0)
    if (allocated(buffer)) deallocate(buffer)

!----------------------------------------------------------------------
! Close scratch files
!----------------------------------------------------------------------
    if (.not.fsincore) close(fsunit)
    
!----------------------------------------------------------------------
! Deallocation of Hamiltonian arrays
!----------------------------------------------------------------------
    if (hincore) call deallocate_hamiltonian
    
    return
    
  end subroutine finalise
  
!######################################################################
  
end module fdstates
