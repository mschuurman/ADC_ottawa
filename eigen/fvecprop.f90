module fvecprop

  use constants
  use parameters
  use channels
  use iomod
  
  implicit none

  integer                               :: dim
  complex(d), dimension(:), allocatable :: psi0
  
contains

!######################################################################

  subroutine propagate_fvec(fvec,ndimf)

    implicit none

    integer, intent(in)                   :: ndimf
    real(d), dimension(ndimf), intent(in) :: fvec

!----------------------------------------------------------------------
! Initialisation and allocatation
!----------------------------------------------------------------------
    call initialise(ndimf)
    
!----------------------------------------------------------------------
! Normalise the F-vector to form |Psi(t=0)>
!----------------------------------------------------------------------
    call init_wavepacket(fvec)

!----------------------------------------------------------------------
! Propagation
!----------------------------------------------------------------------
    !call propagate
    
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
    
    ! Wavefunction dimension
    dim=ndimf
    
    ! Psi(t=0)
    allocate(psi0(dim))
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

    integer                             :: i
    real(d), dimension(dim), intent(in) :: fvec
    real(d)                             :: norm
    
!----------------------------------------------------------------------
! The initial wavepacket is taken as mu |Psi_0>/|| mu |Psi_0> ||
!----------------------------------------------------------------------
    do i=1,dim
       psi0(i)=complex(fvec(i),0.0d0)
    enddo

    norm=sqrt(dot_product(psi0,psi0))
    psi0=psi0/norm

    print*,sqrt(dot_product(psi0,psi0))
    
    STOP
    
    return
    
  end subroutine init_wavepacket

!######################################################################

  subroutine propagate

    implicit none

    integer :: i,nstep
    real(d) :: t

    ! Expokit
    integer                            :: lwsp,liwsp,itrace,iflag
    integer, dimension(:), allocatable :: iwsp
    real(d), dimension(:), allocatable :: wsp

    external matxvec

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    lwsp=2*(matdim*(krydim+1)+matdim+(krydim+2)**2+4*(krydim+2)**2+7)
    liwsp=2*(krydim+2)
    allocate(wsp(lwsp))
    allocate(iwsp(liwsp))

    ! Expokit flags: supress output
    itrace=0

!-----------------------------------------------------------------------
! Calculate the infinity norm of the Hamiltonian matrix
!-----------------------------------------------------------------------
    print*,"WRITE THE INFINITY NORM CODE!"
    ! See the relaxation code.
    STOP
    
    !call infnorm(hnorm,matdim)
    
!----------------------------------------------------------------------
! Propagate |Psi(t=0)> forwards in time using the expokit SIL library    
!----------------------------------------------------------------------
    ! Loop over timesteps
    do i=1,tf/dt

       t=i*dt
       print*,t

       
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(wsp)
    deallocate(iwsp)
    
    return
    
  end subroutine propagate
    
!######################################################################
  
end module fvecprop
