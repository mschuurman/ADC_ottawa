module fvecprop

  use constants
  use parameters
  use channels
  use iomod
  
  implicit none

  integer                               :: matdim
  complex*16, dimension(:), allocatable :: psi0

  integer                               :: maxbl,nrec
  integer, dimension(:), allocatable    :: indxi,indxj
  real(d), dimension(:), allocatable    :: hii,hij

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
    call propagate
    
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
    
    ! Autocorrelation function
    nstep=int(tf/dt)+1
    allocate(auto(nstep))
    auto=czero

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
! The initial wavepacket is taken as mu |Psi_0>/|| mu |Psi_0> ||
!----------------------------------------------------------------------
    do i=1,matdim
       psi0(i)=complex(fvec(i),0.0d0)
    enddo

    norm=sqrt(dot_product(psi0,psi0))
    psi0=psi0/norm

    return
    
  end subroutine init_wavepacket

!######################################################################

  subroutine propagate

    implicit none

    integer                               :: i,kdim
    real(d)                               :: t,hnorm,tol
    complex*16, dimension(:), allocatable :: psi1,psi2

    ! Expokit
    integer                               :: lwsp,liwsp,itrace,iflag
    integer, dimension(:), allocatable    :: iwsp
    complex*16, dimension(:), allocatable :: wsp

    external matxvec_treal

!----------------------------------------------------------------------
! Temporary hard-wiring of the maximum Krylov subspace dimension and
! error tolerance
!----------------------------------------------------------------------
    kdim=7
    tol=1e-6

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Wavefunction arrays
    allocate(psi1(matdim))
    allocate(psi2(matdim))
    psi1=czero
    psi2=czero

    ! Expokit work arrays
    lwsp=matdim*(kdim+1)+matdim+(kdim+2)**2+4*(kdim+2)**2+7
    liwsp=kdim+2
    allocate(wsp(lwsp))
    allocate(iwsp(liwsp))

    ! Expokit flags: supress output
    itrace=0

!-----------------------------------------------------------------------
! Calculate the infinity norm of the Hamiltonian matrix
!-----------------------------------------------------------------------
    call infnorm_outofcore(hnorm,matdim)
    
!----------------------------------------------------------------------
! Propagate |Psi(t=0)> = mu |Psi_0> forwards in time using the
! expokit libraries
!----------------------------------------------------------------------
    ! Set the initial wavepacket
    psi1=psi0

    ! Set the autocorrelation function at time t=0
    auto(1)=dot_product(psi0,psi0)

    ! Loop over timesteps
    do i=1,nstep-1

       ! Update the time
       t=i*dt

       ! Propagate the wavepacket forwards one timestep
       call zgexpv(matdim,kdim,dt,psi1,psi2,tol,hnorm,wsp,lwsp,iwsp,&
            liwsp,matxvec_treal,itrace,iflag)
       if (iflag.ne.0) then
          errmsg='Wavepacket propagation failed...'
          call error_control
       endif

       ! Calculate the autocorrelation function at the current
       ! timestep
       auto(i+1)=dot_product(psi0,psi2)

       ! Output some information
       print*,
       print*,"Time:",t
       print*,"Norm:",real(sqrt(dot_product(psi2,psi2)))

       ! Reset the initial wavepacket
       psi1=psi2

    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(wsp)
    deallocate(iwsp)
    deallocate(psi1)
    deallocate(psi2)
    
    return
    
  end subroutine propagate

!######################################################################
    
  subroutine infnorm_outofcore(norm,matdim)

    use iomod
    use misc, only: dsortindxa1

    implicit none

    integer                            :: k,l,iham,nlim
    integer, intent(in)                :: matdim
    integer, dimension(:), allocatable :: rindx
    real(d)                            :: norm
    real(d), dimension(:), allocatable :: rsum
    character(len=70)                  :: filename
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rsum(matdim))
    rsum=0.0d0

    allocate(rindx(matdim))
    rindx=0

!----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!----------------------------------------------------------------------
    allocate(hii(matdim))

    if (hamflag.eq.'i') then
       filename='SCRATCH/hmlt.diai'
    else if (hamflag.eq.'f') then
       filename='SCRATCH/hmlt.diac'
    endif

    call freeunit(iham)
    open(iham,file=filename,status='old',access='sequential',&
       form='unformatted')

    read(iham) maxbl,nrec
    read(iham) hii

    close(iham)

    rsum=abs(hii)

    deallocate(hii)
    
!----------------------------------------------------------------------
! Contribution from the off-diagonal elements
!----------------------------------------------------------------------
    ! Out-of-core storage of the Hamiltonian
    allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))

    if (hamflag.eq.'i') then
       filename='SCRATCH/hmlt.offi'
    else if (hamflag.eq.'f') then
       filename='SCRATCH/hmlt.offc'
    endif
    
    call freeunit(iham)
    open(iham,file=filename,status='old',access='sequential',&
         form='unformatted')
      
    do k=1,nrec
       read(iham) hij(:),indxi(:),indxj(:),nlim
       do l=1,nlim
          rsum(indxi(l))=rsum(indxi(l))+abs(hij(l))
       enddo
    enddo
      
    close(iham)
    
    deallocate(hij,indxi,indxj)

!----------------------------------------------------------------------
! Determination of the infinity norm of the Hamiltonian matrix.
! (Note that this is equal to the infinity norm of -iH)
!----------------------------------------------------------------------
    call dsortindxa1('D',matdim,rsum,rindx)
    norm=rsum(rindx(1))

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rsum)
    deallocate(rindx)
      
    return

  end subroutine infnorm_outofcore

!######################################################################
  
end module fvecprop
