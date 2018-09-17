!######################################################################
! fdstates_cfd: Routines for the calculation of Chebyshev filter
!               diagonalisation eigenstates
!######################################################################
module fdstates_cfd

  use constants
  use parameters
  use iomod
  use channels
  
  implicit none

  integer                               :: matdim
  integer*8                             :: noffdiag
  real(dp), dimension(:,:), allocatable :: eigvec
  real(dp), dimension(:), allocatable   :: eigval
  real(dp), dimension(2)                :: bounds
  
contains

!######################################################################
  
  subroutine calc_fdstates_cfd(q0,ndimf,noffdf)

    use tdsemod
    use timingmod
  
    implicit none
    
    integer, intent(in)        :: ndimf
    integer*8, intent(in)      :: noffdf
    integer                    :: k
    real(dp), dimension(ndimf) :: q0
    real(dp)                   :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of Chebyshev filter &
         diagonalisation eigenstates'
    write(ilog,'(70a,/)') ('-',k=1,70)

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    call fdstates_cfd_initialise(ndimf,noffdf)

!----------------------------------------------------------------------
! Determine what can be held in memory
!----------------------------------------------------------------------
    call memory_managment

!----------------------------------------------------------------------
! Loading of the non-zero elements of the Hamiltonian matrix into
! memory
!----------------------------------------------------------------------
    if (hincore) call load_hamiltonian('SCRATCH/hmlt.diac',&
         'SCRATCH/hmlt.offc',matdim,noffdiag)

!----------------------------------------------------------------------
! Spectral bounds
!----------------------------------------------------------------------
    bounds=ebound

!----------------------------------------------------------------------
! Calculate the Chebyshev filter diagonalisation eigenstates
!----------------------------------------------------------------------
    call calc_eigenstates(q0)


    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
   call times(tw2,tc2)
   write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
   return
    
 end subroutine calc_fdstates_cfd

!######################################################################

 subroutine fdstates_cfd_initialise(ndimf,noffdf)
   
   implicit none
   
   integer, intent(in)   :: ndimf
   integer*8, intent(in) :: noffdf

!----------------------------------------------------------------------
! Dimensions
!----------------------------------------------------------------------
   matdim=ndimf
   noffdiag=noffdf
   
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
   ! Filter diagonalisation eigenstates
   allocate(eigvec(matdim,nsel))
   eigvec=0.0d0

   ! Filter diagonalisation eigenvalues
   allocate(eigval(nsel))
   eigval=0.0d0
   
   return
   
 end subroutine fdstates_cfd_initialise

!######################################################################

   subroutine fdstates_cfd_finalise
    
    use tdsemod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(eigvec)
    deallocate(eigval)
    
    if (hincore) call deallocate_hamiltonian
    
    return
    
  end subroutine fdstates_cfd_finalise
 
!######################################################################

  subroutine memory_managment

    use tdsemod
    use omp_lib
    
    implicit none

    integer*8  :: maxrecl,reqmem
    integer    :: nthreads
    real(dp)   :: memavail

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

    ! Filter diagonalisation eigenstates
    memavail=memavail+8.0d0*neig/1024.0d0**2
    
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

  subroutine calc_eigenstates(q0)

    use tdsemod
        
    implicit none

    integer                     :: k,i,unit
    real(dp), dimension(matdim) :: q0
    real(dp), allocatable       :: qk(:),qkm1(:),qkm2(:)
    real(dp)                    :: norm

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
! On-the-fly calculation of the Chebyshev filter diagonalisation
! eigenstates
!----------------------------------------------------------------------
    ! Zeroth-order contribution to the eigenstates
    do i=1,nsel
       eigvec(:,i)=k2eig(0,isel(i))*q0(:)
    enddo

    ! Initialisation
    qkm1=q0

    ! Loop over Chebyshev polynomials of order k >= 1
    do k=1,chebyord

       ! Output our progress
       if (mod(k,10).eq.0) then
          write(ilog,'(70a)') ('+',i=1,70)
          write(ilog,'(a,x,i6)') 'Order:',k
       endif

       ! Calculate the kth Chebyhev polynomial-vector product
       call chebyshev_recursion(k,matdim,noffdiag,bounds,qk,qkm1,qkm2)

       ! kth-order contribution to the eigenstates
       do i=1,nsel
          eigvec(:,i)=eigvec(:,i)+k2eig(k,isel(i))*qk(:)
       enddo

       ! Update qkm1 and qkm2
       qkm2=qkm1       
       qkm1=qk
       qk=0.0d0
       
    enddo

    ! Last report
    if (mod(k,10).ne.0) then
       write(ilog,'(70a)') ('+',i=1,70)
       write(ilog,'(a,x,i6)') 'Order:',k
    endif
    write(ilog,'(70a)') ('+',i=1,70)

!----------------------------------------------------------------------
! Normalisation
!----------------------------------------------------------------------
    do i=1,nsel
       norm=sqrt(dot_product(eigvec(:,i),eigvec(:,i)))
       eigvec(:,i)=eigvec(:,i)/norm
    enddo
       
!----------------------------------------------------------------------
! Calculation of energies
!----------------------------------------------------------------------
    do i=1,nsel
       call matxvec(matdim,noffdiag,eigvec(:,i),qk)
       qk=-qk
       eigval(i)=dot_product(eigvec(:,i),qk)
    enddo

!----------------------------------------------------------------------
! Write the eigenstates to file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit=unit,file='SCRATCH/fdstates',status='unknown',&
         access='sequential',form='unformatted')
    do i=1,nsel
       write(unit) i,eigval(i),eigvec(:,i)
    enddo
    close(unit)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(qk)
    deallocate(qkm1)
    deallocate(qkm2)
    
    return
    
  end subroutine calc_eigenstates
    
!######################################################################

end module fdstates_cfd
