!######################################################################
! chebyspec: routines for the calculation of the Chebyshev order-domain
!            autocorrelation function
!######################################################################

module chebyspec

  use constants
  use parameters
  use iomod
  
  implicit none

  integer               :: matdim
  integer*8             :: noffdiag
  real(d), dimension(2) :: bounds
  real(d), allocatable  :: auto(:)
  
contains

!######################################################################

  subroutine chebyshev_recursion(q0,ndimf,noffdf)

    use tdsemod
    use specbounds
    
    implicit none

    integer, intent(in)                   :: ndimf
    integer*8, intent(in)                 :: noffdf
    real(d), dimension(ndimf), intent(in) :: q0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    call chebyshev_initialise(ndimf,noffdf)

!----------------------------------------------------------------------
! Determine what can be held in memory
!----------------------------------------------------------------------
    call memory_managment

!----------------------------------------------------------------------
! Loading of the non-zero elements of the Hamiltonian matrix into
! memory
!----------------------------------------------------------------------
    if (hincore) call load_hamiltonian('SCRATCH/hmlt.diac',&
         'SCRATCH/hmlt.offc',matdim,noffdf)
    
!----------------------------------------------------------------------
! Estimation of the spectral bounds using a few-iteration Lanczos
! calculation
!----------------------------------------------------------------------
    call spectral_bounds(bounds,'c','lanczos',ndimf,noffdf)

!----------------------------------------------------------------------
! Calculate the order-domain autocorrelation function
!----------------------------------------------------------------------
   call chebyshev_auto(q0,ndimf,noffdf)

   return
    
 end subroutine chebyshev_recursion

!######################################################################

  subroutine chebyshev_initialise(ndimf,noffdf)

    implicit none

    integer, intent(in)                   :: ndimf
    integer*8, intent(in)                 :: noffdf
    
!----------------------------------------------------------------------
! Dimensions
!----------------------------------------------------------------------
    matdim=ndimf
    noffdiag=noffdf
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    allocate(auto(0:chebyord))
    auto=0.0d0
    
    return
    
  end subroutine chebyshev_initialise
    
!######################################################################

  subroutine chebyshev_finalise

    use tdsemod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(auto)
    if (hincore) call deallocate_hamiltonian
    
    return
    
  end subroutine chebyshev_finalise

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

    ! Chebyshev polynomial-vector products
    memavail=memavail-8.0d0*4*matdim/1024.0d0**2

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

    print*,
    print*,'hincore:',hincore
    print*,
    stop
    
    
    return
    
  end subroutine memory_managment
    
!######################################################################

  subroutine chebyshev_auto(q0,ndimf,noffdf)

    use tdsemod
    
    implicit none

    integer, intent(in)                   :: ndimf
    integer*8, intent(in)                 :: noffdf
    integer                               :: k
    real(d), dimension(ndimf), intent(in) :: q0
    real(d)                               :: dummy
    real(d), allocatable                  :: qk(:),q1(:),q2(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qk(matdim))
    qk=0.0d0

    allocate(q1(matdim))
    q1=0.0d0

    allocate(q2(matdim))
    q2=0.0d0
    
!----------------------------------------------------------------------
! Calculate the Chebyshev order-domain autocorrelation function
!----------------------------------------------------------------------
    ! Loop over Chebyshev polynomials
    do k=1,chebyord

       ! Calculate the kth Chebyhev polynomial-vector product
       !call matxvec_chebyshev(matdim,noffdiag,bounds,qk,q1,q2)

       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(qk)
    deallocate(q1)
    deallocate(q2)
    
    return
    
  end subroutine chebyshev_auto
    
!######################################################################
  
end module chebyspec
