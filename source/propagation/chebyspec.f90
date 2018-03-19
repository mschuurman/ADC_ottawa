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

  subroutine chebyshev_auto_order_domain(q0,ndimf,noffdf)

    use tdsemod
    use specbounds
    
    implicit none

    integer, intent(in)       :: ndimf
    integer*8, intent(in)     :: noffdf
    real(d), dimension(ndimf) :: q0
    
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
    
 end subroutine chebyshev_auto_order_domain

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

    return
    
  end subroutine memory_managment
    
!######################################################################

  subroutine chebyshev_auto(q0,ndimf,noffdf)

    use tdsemod
    
    implicit none

    integer, intent(in)       :: ndimf
    integer*8, intent(in)     :: noffdf
    integer                   :: k
    real(d), dimension(ndimf) :: q0
    real(d), allocatable      :: qk(:),qkm1(:),qkm2(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(qk(matdim))
    qk=0.0d0

    allocate(qkm1(matdim))
    qkm1=0.0d0

    allocate(qkm2(matdim))
    qkm2=0.0d0

!----------------------------------------------------------------------
! C_0
!----------------------------------------------------------------------
    auto(0)=dot_product(q0,q0)

    print*,0,auto(0)
    
!----------------------------------------------------------------------
! Calculate the Chebyshev order-domain autocorrelation function
!----------------------------------------------------------------------
    ! Initialisation
    qkm1=q0

    ! Loop over Chebyshev polynomials of order k >= 1
    do k=1,chebyord

       ! Calculate the kth Chebyhev polynomial-vector product
       call chebyshev_recursion(k,matdim,noffdiag,bounds,qk,qkm1,qkm2)

       ! Calculate C_k
       auto(k)=dot_product(q0,qk)

       print*,k,auto(k)
       
       ! Update qkm1 and qkm2
       qkm1=qk
       qkm2=qkm1
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(qk)
    deallocate(qkm1)
    deallocate(qkm2)
    
    return
    
  end subroutine chebyshev_auto
    
!######################################################################
  
end module chebyspec
