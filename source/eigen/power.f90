!#######################################################################
! power: Routine for the calculation of the highest-lying eigenpair
!        of the ADC Hamiltonian matrix using the power method.
!#######################################################################

module power

  use constants
  use parameters
  use channels
  use iomod
  use timingmod
  
  implicit none

  integer               :: nmult,nrec,maxbl
  integer, allocatable  :: indxi(:),indxj(:)
  real(dp), allocatable :: hii(:),hij(:)
  real(dp), allocatable :: v1(:),v2(:),vres(:)
  real(dp)              :: eval
  logical               :: lincore

contains

!#######################################################################

  subroutine power_method(matdim,noffd,flag)

    implicit none

    integer, intent(in)          :: matdim
    integer*8, intent(in)        :: noffd
    integer                      :: k
    real(dp)                     :: tw1,tw2,tc1,tc2
    character(len=1), intent(in) :: flag
    character(len=120)           :: atmp

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)

!-----------------------------------------------------------------------
! Write to the log file
!-----------------------------------------------------------------------
    atmp='Power iterations in the'
    if (flag.eq.'i') then
       atmp=trim(atmp)//' initial space'
    else if (flag.eq.'c') then
       atmp=trim(atmp)//' final space'
    endif
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') trim(atmp)
    write(ilog,'(70a,/)') ('-',k=1,70)
    
!-----------------------------------------------------------------------    
! Initialisation
!-----------------------------------------------------------------------    
    call power_initialise(matdim)

!-----------------------------------------------------------------------
! Determine whether or not we can perform the matrix-vector
! multiplications in-core
!-----------------------------------------------------------------------
    call power_isincore(matdim,noffd)
    
!-----------------------------------------------------------------------
! Read the on-diagonal elements of the Hamiltonian matrix from disk
!-----------------------------------------------------------------------
    call rdham_on(matdim,flag)

!-----------------------------------------------------------------------
! If the matrix-vector multiplication is to proceed in-core, then
! read the off-diagonal Hamiltonian matrix from disk
!-----------------------------------------------------------------------
    if (lincore) call rdham_off(noffd,flag)
    
!-----------------------------------------------------------------------
! Set the guess vector
!-----------------------------------------------------------------------
    call power_guess(matdim,noffd)

!-----------------------------------------------------------------------    
! Perform the power iterations
!-----------------------------------------------------------------------    
    call power_iterations(matdim,noffd,flag)

!-----------------------------------------------------------------------    
! Write the eigenpair to file
!-----------------------------------------------------------------------    
    call wreigenpair(flag)
    
!-----------------------------------------------------------------------    
! Finalisation
!-----------------------------------------------------------------------    
    call power_finalise
    
!-----------------------------------------------------------------------    
! Output timings and the no. matrix-vector multiplications
!-----------------------------------------------------------------------    
    write(ilog,'(/,a,1x,i4)') 'No. matrix-vector multiplications:',&
         nmult

    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine power_method
    
!#######################################################################

  subroutine power_initialise(matdim)

    implicit none

    integer, intent(in) :: matdim

!-----------------------------------------------------------------------
! Initialisation of the matrix-vector multiplication counter
!-----------------------------------------------------------------------
    nmult=0
    
!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
    ! Vectors
    allocate(v1(matdim))
    allocate(v2(matdim))
    allocate(vres(matdim))
    v1=0.0d0
    v2=0.0d0
    vres=0.0d0
    
    ! On-diagonal elements of the Hamiltonian matrix
    allocate(hii(matdim))
    hii=0.0d0    
    
    return
    
  end subroutine power_initialise

!#######################################################################

  subroutine power_finalise

    implicit none

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(v1)
    deallocate(v2)
    deallocate(vres)
    deallocate(hii)
    if (allocated(hij)) deallocate(hij)
    if (allocated(indxi)) deallocate(indxi)
    if (allocated(indxj)) deallocate(indxj)
    
    return
    
  end subroutine power_finalise

!#######################################################################

  subroutine power_isincore(matdim,noffd)

    use constants
    use parameters, only: maxmem
    use omp_lib
    
    implicit none
    
    integer, intent(in)   :: matdim
    integer*8, intent(in) :: noffd
    integer               :: nthreads
    real(dp)              :: mem

    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    
    mem=0.0d0

    ! On-diagonal Hamiltonian matrix elements
    mem=mem+8.0d0*matdim/1024.0d0**2

    ! Non-zero off-diagonal Hamiltonian matrix element values
    mem=mem+8.0d0*noffd/1024.0d0**2

    ! Indices of the non-zero off-diagonal Hamiltonian matrix elements
    mem=mem+8.0d0*noffd/1024.0d0**2
    
    ! Vectors
    mem=mem+8.0d0*(nthreads+3)*matdim/1024.0d0**2

    ! Two-electron integrals held in-core
    mem=mem+8.0d0*(nbas**4)/1024.0d0**2

    ! kpq
    mem=mem+8.0d0*7.0d0*(1+nbas**2*4*nocc**2)/1024.0d0**2

    if (mem.lt.maxmem) then
       lincore=.true.
    else
       lincore=.false.
    endif

    return
    
  end subroutine power_isincore
    
!#######################################################################

  subroutine rdham_on(matdim,flag)

    implicit none

    integer, intent(in)          :: matdim
    integer                      :: iham
    character(len=1), intent(in) :: flag
    character(len=70)            :: filename
    
!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    call freeunit(iham)

    if (flag.eq.'i') then
       filename='SCRATCH/hmlt.diai'
    else if (flag.eq.'c') then
       filename='SCRATCH/hmlt.diac'
    endif

    open(iham,file=filename,status='old',access='sequential',&
           form='unformatted')

    read(iham) maxbl,nrec
    read(iham) hii

    close(iham)

    return

  end subroutine rdham_on

!#######################################################################

  subroutine rdham_off(noffd,flag)

    use iomod, only: freeunit
    use constants
    use parameters

    implicit none

    integer*8, intent(in)               :: noffd
    integer                             :: iham,count,k,nlim
    integer, dimension(:), allocatable  :: itmp1,itmp2
    real(dp), dimension(:), allocatable :: ftmp
    character(len=1), intent(in)        :: flag
    character(len=70)                   :: filename

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(hij(noffd))
    allocate(indxi(noffd))
    allocate(indxj(noffd))

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
    allocate(ftmp(maxbl))
    allocate(itmp1(maxbl))
    allocate(itmp2(maxbl))
    
    if (flag.eq.'i') then
       filename='SCRATCH/hmlt.offi'
    else if (flag.eq.'c') then
       filename='SCRATCH/hmlt.offc'
    endif

    call freeunit(iham)
    open(iham,file=filename,status='old',access='sequential',&
         form='unformatted')

    count=0
    do k=1,nrec
       read(iham) ftmp(:),itmp1(:),itmp2(:),nlim
       hij(count+1:count+nlim)=ftmp(1:nlim)
       indxi(count+1:count+nlim)=itmp1(1:nlim)
       indxj(count+1:count+nlim)=itmp2(1:nlim)
       count=count+nlim
    enddo

    deallocate(ftmp,itmp1,itmp2)

    close(iham)
    
    return
    
  end subroutine rdham_off
    
!#######################################################################

  subroutine power_guess(matdim,noffd)

    use misc, only: dsortindxa1
    
    implicit none

    integer, intent(in)                   :: matdim
    integer*8, intent(in)                 :: noffd
    integer, dimension(:), allocatable    :: full2sub,sub2full,indxhii
    integer                               :: subdim
    integer*8                             :: k,i,j
    integer                               :: i1,j1
    integer                               :: e2,error,iham,nlim,l
    real(dp), dimension(:,:), allocatable :: hsub
    real(dp), dimension(:), allocatable   :: subeig,work
    character(len=70)                     :: filename
    
!-----------------------------------------------------------------------
! Temporary hardwiring of the subspace dimension
!-----------------------------------------------------------------------
    subdim=800

!-----------------------------------------------------------------------
! Subspace dimension check
!-----------------------------------------------------------------------
    if (subdim.gt.matdim) subdim=matdim

!-----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix elements in order of
! descending value
!-----------------------------------------------------------------------
    allocate(indxhii(matdim))
    call dsortindxa1('D',matdim,hii,indxhii)
    
!-----------------------------------------------------------------------
! Ensure that the subdim'th IS is not degenerate with subdim+1'th IS,
! and if it is, increase subdim accordingly
!-----------------------------------------------------------------------
5   continue
    if (subdim.lt.matdim) then
       if (abs(hii(indxhii(subdim))-hii(indxhii(subdim+1))).lt.1e-6_dp) then
          subdim=subdim+1
          goto 5
       endif
    endif

!-----------------------------------------------------------------------
! Allocate the subspace-associated arrays
!-----------------------------------------------------------------------
    allocate(full2sub(matdim))
    allocate(sub2full(subdim))
    allocate(hsub(subdim,subdim))
    allocate(subeig(subdim))
    allocate(work(3*subdim))

!-----------------------------------------------------------------------
! Set the full space-to-subsace mappings
!-----------------------------------------------------------------------
    full2sub=0
    do i=1,subdim
       k=indxhii(i)
       sub2full(i)=k
       full2sub(k)=i
    enddo

!-----------------------------------------------------------------------
! Construct the Hamiltonian matrix in the subspace
!-----------------------------------------------------------------------
    hsub=0.0d0

    ! On-diagonal elements
    do i=1,subdim
       k=sub2full(i)
       hsub(i,i)=hii(k)
    enddo

    ! Off-diagonal elements
    if (lincore) then
       ! Loop over all off-diagonal elements of the full space
       ! Hamiltonian
       do k=1,noffd
          ! Indices of the current off-diagonal element of the
          ! full space Hamiltonian
          i=indxi(k)
          j=indxj(k)
          ! If both indices correspond to subspace ISs, then
          ! add the element to subspace Hamiltonian
          if (full2sub(i).ne.0.and.full2sub(j).ne.0) then
             i1=full2sub(i)
             j1=full2sub(j)
             hsub(i1,j1)=hij(k)               
             hsub(j1,i1)=hsub(i1,j1)
          endif
       enddo
    else
       ! Open the off-diagonal element file
       call freeunit(iham)
       if (hamflag.eq.'i') then
          filename='SCRATCH/hmlt.offi'
       else if (hamflag.eq.'f') then
          filename='SCRATCH/hmlt.offc'
       endif
       open(iham,file=filename,status='old',access='sequential',&
            form='unformatted')
       ! Allocate arrays
       allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
       ! Loop over records
       do k=1,nrec
          read(iham) hij(:),indxi(:),indxj(:),nlim
          ! Loop over the non-zero elements of the full space
          ! Hamiltonian in the current record
          do l=1,nlim
             ! Indices of the current off-diagonal element of the
             ! full space Hamiltonian
             i=indxi(l)
             j=indxj(l)
             ! If both indices correspond to subspace ISs, then
             ! add the element to subspace Hamiltonian
             if (full2sub(i).ne.0.and.full2sub(j).ne.0) then
                i1=full2sub(i)
                j1=full2sub(j)
                hsub(i1,j1)=hij(l)
                hsub(j1,i1)=hsub(i1,j1)
             endif
          enddo
       enddo
       ! Close the off-diagonal element file
       close(iham)
       ! Deallocate arrays
       deallocate(hij,indxi,indxj)
    endif

!-----------------------------------------------------------------------
! Diagonalise the subspace Hamiltonian
!-----------------------------------------------------------------------
    e2=3*subdim
    call dsyev('V','U',subdim,hsub,subdim,subeig,work,e2,error)
    
    if (error.ne.0) then
       errmsg='The diagonalisation of the subspace Hamiltonian failed.'
       call error_control
    endif

!-----------------------------------------------------------------------
! Set the guess vector as the highest-lying eigenvector of the subspace
! Hamiltonian.
! Note that after calling dsyev, hsub now holds the eigenvectors of the
! subspace Hamiltonian.
!-----------------------------------------------------------------------
    do j=1,subdim
       k=sub2full(j)
       v1(k)=hsub(j,subdim)
    enddo

    return
    
  end subroutine power_guess

!#######################################################################

  subroutine power_iterations(matdim,noffd,flag)

    implicit none

    integer, intent(in)          :: matdim
    integer*8, intent(in)        :: noffd
    integer                      :: i,k
    integer, parameter           :: maxit=350
    real(dp), parameter          :: thrsh=1e-8_dp
    real(dp)                     :: resnorm
    character(len=1), intent(in) :: flag

!-----------------------------------------------------------------------
! Perform the power iterations
!-----------------------------------------------------------------------
    do i=1,maxit
       
       ! v2 = H*v1
       call power_hxvec(matdim,noffd,flag)

       ! Ritz value for v1
       eval=dot_product(v1,v2)

       ! Residual for v1
       vres=eval*v1-v2
       resnorm=sqrt(dot_product(vres,vres))

       ! Normalise v2
       v2=v2/sqrt(dot_product(v2,v2))

       ! Output our progress
       write(ilog,'(70a)') ('+',k=1,70)
       write(ilog,'(a,x,i5)') 'Iteration:',i
       write(ilog,'(a,x,ES21.14)') 'Residual: ',resnorm
       
       ! Terminate if we have converged...
       if (resnorm.le.thrsh) then
          exit
       endif
          
       ! ...else v2 <- v1 and we keep on iterating
       v1=v2
       
    enddo

    write(ilog,'(70a)') ('+',k=1,70)
    
!-----------------------------------------------------------------------
! Exit here if we did not reach convergence
!-----------------------------------------------------------------------
    if (resnorm.gt.thrsh) then
       errmsg='Convergence not reached...'
       call error_control
    endif
    
    return
    
  end subroutine power_iterations

!#######################################################################

  subroutine power_hxvec(matdim,noffd,flag)

    implicit none

    integer, intent(in)          :: matdim
    integer*8, intent(in)        :: noffd
    character(len=1), intent(in) :: flag

    if (lincore) then
       call power_hxvec_incore(matdim,noffd)
    else
       print*,"WRITE THE OUT-OF-CORE MATRIX-VECTOR MULTIPLICATION CODE!"
       stop
       !call power_hxvec_ext(matdim,noffd,flag)
    endif

    nmult=nmult+1
    
    return

  end subroutine power_hxvec

!#######################################################################

  subroutine power_hxvec_incore(matdim,noffd)

    use omp_lib
    
    implicit none
    
    integer, intent(in)          :: matdim
    integer*8, intent(in)        :: noffd
    integer*8                    :: npt
    integer                      :: i,k,nthreads,tid
    integer*8, allocatable       :: irange(:,:)
    real(dp), allocatable        :: tmpvec(:,:)

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(irange(nthreads,2))
    allocate(tmpvec(matdim,nthreads))

!-----------------------------------------------------------------------
! Partitioning of the off-diagonal elements: one chunk per thread
!-----------------------------------------------------------------------
    npt=int(floor(real(noffd)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffd

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=hii(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=0.0d0

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(indxi(k),tid)=tmpvec(indxi(k),tid)&
               +hij(k)*v1(indxj(k))
          tmpvec(indxj(k),tid)=tmpvec(indxj(k),tid)&
               +hij(k)*v1(indxi(k))
       enddo
          
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo
       
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(irange)
    deallocate(tmpvec)
    
    return
    
  end subroutine power_hxvec_incore

!#######################################################################

  subroutine wreigenpair(flag)

    implicit none

    integer                      :: ipower
    character(len=60)            :: filename
    character(len=1), intent(in) :: flag

    call freeunit(ipower)
    
    filename='SCRATCH/'//'powerstate'
    if (flag.eq.'c') filename=trim(filename)//'_final'

    open(unit=ipower,file=filename,status='unknown',access='sequential',&
         form='unformatted')

    write(ipower) eval,v1(:)
    
    close(ipower)
        
    return
    
  end subroutine wreigenpair
    
!#######################################################################
  
end module power

