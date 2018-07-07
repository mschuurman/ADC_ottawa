!######################################################################
! specbounds: Routines for the estimation of the spectral bounds of
!             the Hamiltonian matrix. These are required if we are
!             performing a wavepacket propagation using the
!             Chebyshev representation of the time-evolution
!             operator.
!######################################################################

module specbounds

  use constants
  use parameters
  use iomod
  
  implicit none

contains

!######################################################################
! spectral_bounds: Gateway routine for the estimation of the spectral
!                  bounds of the Hamiltonian
!######################################################################
  
  subroutine spectral_bounds(bounds,flag,estimation,dim,noff)

    implicit none

    integer, intent(in)    :: dim
    integer*8, intent(in)  :: noff
    real(dp), dimension(2) :: bounds
    character(len=1)       :: flag
    character(len=*)       :: estimation

    if (estimation.eq.'lanczos') then
       call spectral_bounds_lanczos(bounds,flag,dim,noff)
    else if (estimation.eq.'davlanc') then
       call spectral_bounds_davlanc(bounds,flag,dim,noff)
    endif
    
    return
    
  end subroutine spectral_bounds

!######################################################################
! spectral_bounds_lanczos: Estimation of the spectral bounds of the
!                          Hamiltonian matrix using a few-iteration
!                          Lanczos procedure.
!######################################################################

  subroutine spectral_bounds_lanczos(bounds,flag,dim,noff)

    use iomod
    use block_lanczos
    
    implicit none

    integer, intent(in)                 :: dim
    integer*8, intent(in)               :: noff
    integer                             :: i,ilanc,n
    real(dp), dimension(2)              :: bounds
    real(dp), dimension(:), allocatable :: vec
    real(dp)                            :: ener
    character(len=1)                    :: flag
    
!----------------------------------------------------------------------
! Set the block Lanczos parameters
!----------------------------------------------------------------------
    ! Initial vectors generated from a double subspace diagonalisation
    ! procedure
    lancguess=8

    ! Block size
    lmain=2

    ! No. iterations
    ncycles=50

!----------------------------------------------------------------------
! Perform the Lanczos iterations
!----------------------------------------------------------------------
    call lancdiag_block(dim,noff,flag)

!----------------------------------------------------------------------
! Read the spectral bounds from file
!----------------------------------------------------------------------
    open(ilanc,file=lancname,access='sequential',form='unformatted',&
         status='old')

    allocate(vec(dim))
    
    do i=1,lancstates
       read(ilanc) n,ener,vec
       if (i.eq.1) bounds(1)=ener
       if (i.eq.lancstates) bounds(2)=ener
    enddo

    deallocate(vec)

    close(ilanc)

    return
    
  end subroutine spectral_bounds_lanczos

!######################################################################
! spectral_bounds_dav_lanc: Estimation of the lower spectral bound
!                           using the block Davidson method, and the
!                           upper spectral bound using the Lanczos
!                           method.
!######################################################################
  
  subroutine spectral_bounds_davlanc(bounds,flag,dim,noff)

    use iomod
    use block_davidson
    use block_lanczos
    
    implicit none

    integer, intent(in)                 :: dim
    integer*8, intent(in)               :: noff
    integer                             :: i,n,idav,ilanc
    real(dp), dimension(2)              :: bounds
    real(dp), dimension(:), allocatable :: vec
    real(dp)                            :: ener
    character(len=1)                    :: flag

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(vec(dim))
    vec=0.0d0
    
!----------------------------------------------------------------------
! Lower spectral bound from a block Davidson calculation
!----------------------------------------------------------------------
    ! Hamiltonian matrix flag
    if (flag.eq.'i') then
       hamflag='i'
    else if (flag.eq.'c') then
       hamflag='f'
    endif

    ! No. states
    davstates_f=1

    ! Block size
    dmain_f=2

    ! Maximum subspace dimension
    maxsubdim_f=10
    
    ! Convergence threshold
    davtol_f=1e-8_dp

    ! Generate guess vectors from a subspace diagonalisation
    lsubdiag_f=.true.
    
    ! Block Davidson calculation
    call davdiag_block(dim,noff)

    ! Read the lower spectral bound from file
    call freeunit(idav)
    open(unit=idav,file=davname_f,status='unknown',&
         access='sequential',form='unformatted')
    read(idav) i,bounds(1),vec(:)
    close(idav)

!----------------------------------------------------------------------
! Upper spectral bound from a block Lanczos calculation
!----------------------------------------------------------------------
    ! Initial vectors generated from a high-energy subspace
    ! diagonalisation procedure
    lancguess=9

    ! Block size
    lmain=2

    ! No. iterations
    ncycles=50

    ! Perform the Lanczos iterations
    call lancdiag_block(dim,noff,flag)
    
    ! Read the upper spectral bound from file
    open(ilanc,file=lancname,access='sequential',form='unformatted',&
         status='old')
    do i=1,lancstates
       read(ilanc) n,ener,vec
    enddo
    bounds(2)=ener
    close(ilanc)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    
    return
    
  end subroutine spectral_bounds_davlanc
    
!######################################################################
  
end module specbounds
