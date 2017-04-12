module fvecprop

  use constants
  use parameters
  use channels
  use iomod
  use timingmod
  
  implicit none

  integer                               :: matdim,iout
  real(d), parameter                    :: au2fs=1.0d0/41.34137333656d0  
  complex(d), dimension(:), allocatable :: psi0

contains

!######################################################################

  subroutine propagate_fvec(fvec,ndimf)

    implicit none

    integer, intent(in)                   :: ndimf
    integer                               :: k
    real(d), dimension(ndimf), intent(in) :: fvec
    real(d)                               :: tw1,tw2,tc1,tc2
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Output where we are at
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a,/,3x,a)') &
         'Calculation of the autocorrelation function:',&
         'a(t) = < Psi_0| D exp(-iHt) D | Psi_0 >'
    write(ilog,'(70a,/)') ('-',k=1,70)

!----------------------------------------------------------------------
! Initialisation and allocatation
!----------------------------------------------------------------------
    call initialise(ndimf)
    
!----------------------------------------------------------------------
! Normalise the F-vector to form |Psi(t=0)>
!----------------------------------------------------------------------
    call init_wavepacket(fvec)

!----------------------------------------------------------------------
! Open the autocorrelation function output file and write the file
! header
!----------------------------------------------------------------------
    call freeunit(iout)
    open(iout,file='auto',form='formatted',status='unknown')
    write(iout,'(a)') '#    time[fs]         Re(autocorrel)     &
         Im(autocorrel)     Abs(autocorrel)'
    
!----------------------------------------------------------------------
! Propagation and calculation of the autocorrelation function
!----------------------------------------------------------------------
    call propagate_sillib

!----------------------------------------------------------------------
! Close the autocorrelation function output file
!----------------------------------------------------------------------
    close(iout)

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
    
  end subroutine propagate_fvec

!######################################################################

  subroutine initialise(ndimf)

    implicit none

    integer, intent(in) :: ndimf
    
    ! Hamiltonian matrix dimension
    matdim=ndimf
    
    ! Psi(t=0)
    allocate(psi0(matdim))
    psi0=czero
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine finalise

    implicit none

    deallocate(psi0)
    
    return
    
  end subroutine finalise
    
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

  subroutine propagate_sillib

    use sillib
    
    implicit none

    integer                               :: i,kdim
    real(d)                               :: norm,tol
    real(d), parameter                    :: tiny=1e-9_d
    complex(d), dimension(:), allocatable :: psi,dtpsi
    complex(d)                            :: auto
    
    ! SIL
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
    
    external matxvec_treal

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
! Temporary hard-wiring of the maximum Krylov subspace dimension and
! error tolerance
!----------------------------------------------------------------------
    kdim=15
    tol=1e-6_d
    
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
! Propagate |Psi(t=0)> = mu |Psi_0> forwards in time using the
! Quantics sillib libraries
!
! Note that we are using t/2 trick here in calculating the
! autocorrelation function
!----------------------------------------------------------------------
    ! Integration period
    intperiod=0.5d0*tout

    ! Set the initial wavepacket
    psi=psi0

    ! Output the autocorrelation function at time t=0
    auto=dot_product(psi0,psi0)
    call wrauto(auto,0.0d0)
    
    ! Loop over the timesteps at which we want the autocorrelation
    ! function
    do i=1,int(tfinal/intperiod)
       
       ! Propagate forwards one timestep
       inttime=0.0d0
100    continue

       ! Update the required stepsize
       stepsize=intperiod-inttime
       
       ! dtpsi = -iH|Psi>
       call matxvec_treal(matdim,psi,dtpsi)       
    
       ! Take one step using the SIL algorithm
       call silstep(psi,dtpsi,matdim,stepsize,kdim,tol,relax,restart,&
            stdform,steps,krylov,truestepsize,trueorder,errorcode,&
            time,matxvec_treal,eigenvector,eigenval,diagonal,&
            offdg2,offdiag)

       ! Exit if the SIL integration failed
       if (errorcode.ne.0) then
          call silerrormsg (errorcode,errmsg)
          call error_control
       endif

       ! Check whether the integration is complete
       inttime=inttime+truestepsize
       if (abs(intperiod-inttime).gt.abs(tiny*intperiod)) goto 100

       ! Calculate and output the autocorrelation function at the
       ! current timestep
       auto=dot_product(conjg(psi),psi)
       call wrauto(auto,i*intperiod*2.0d0)
       
       ! Output some information
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
    
  end subroutine propagate_sillib

!######################################################################

  subroutine wrprogress(t,norm)

    implicit none

    integer :: k
    real(d) :: t,norm

    write(ilog,'(70a)') ('+',k=1,70)
    write(ilog,'(a,x,F10.4)') 'Time:',t
    write(ilog,'(a,x,F8.6)') 'Norm',norm
    
    return
    
  end subroutine wrprogress

!######################################################################

  subroutine wrauto(auto,t)

    implicit none

    complex(d) :: auto
    real(d)    :: t

    write(iout,'(F15.8,4x,3(2x,F17.14))') &
         t*au2fs,real(auto),aimag(auto),abs(auto) 
    
    return
    
  end subroutine wrauto
    
!######################################################################
    
end module fvecprop
