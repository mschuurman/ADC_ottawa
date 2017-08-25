!#######################################################################
! tsdemod: common routines used in electronic wavepacket propagation
!          calculations.
!#######################################################################

module tdsemod

  use constants
  
  save

  integer                            :: mxbl,nrecord
  integer, dimension(:), allocatable :: iindx,jindx
  real(d), dimension(:), allocatable :: hon,hoff
  logical                            :: hincore
  
contains

!#######################################################################
! load_hamiltonian: loads the non-zero elements of the Hamiltonian
!                   matrix into memory
!#######################################################################

  subroutine load_hamiltonian(fileon,fileoff,matdim,noffdiag)

    use iomod
    
    implicit none

    integer, intent(in)                :: matdim
    integer*8, intent(in)              :: noffdiag
    integer                            :: unit,count,k,nlim
    integer, dimension(:), allocatable :: itmp,jtmp
    real(d), dimension(:), allocatable :: ftmp
    character(len=*), intent(in)       :: fileon,fileoff

!-----------------------------------------------------------------------
! Get the next free file unit
!-----------------------------------------------------------------------
    call freeunit(unit)
    
!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    if (allocated(hon)) deallocate(hon)
    allocate(hon(matdim))

    open(unit,file=fileon,status='old',access='sequential',&
         form='unformatted')

    read(unit) mxbl,nrecord
    read(unit) hon
    
    close(unit)

!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    if (allocated(hoff)) deallocate(hoff)
    if (allocated(iindx)) deallocate(iindx)
    if (allocated(jindx)) deallocate(jindx)
    allocate(hoff(noffdiag))
    allocate(iindx(noffdiag))
    allocate(jindx(noffdiag))

    allocate(ftmp(mxbl))
    allocate(itmp(mxbl))
    allocate(jtmp(mxbl))

    open(unit,file=fileoff,status='old',access='sequential',&
           form='unformatted')

    count=0
    do k=1,nrecord
       read(unit) ftmp(:),itmp(:),jtmp(:),nlim
       hoff(count+1:count+nlim)=ftmp(1:nlim)
       iindx(count+1:count+nlim)=itmp(1:nlim)
       jindx(count+1:count+nlim)=jtmp(1:nlim)
       count=count+nlim
    enddo
      
    close(unit)
    
    deallocate(ftmp,itmp,jtmp)

    return
    
  end subroutine load_hamiltonian

!#######################################################################
! deallocate_hamiltonian: deallocates the hamiltonian arrays
!#######################################################################

  subroutine deallocate_hamiltonian

    implicit none
    
    if (allocated(hon)) deallocate(hon)
    if (allocated(hoff)) deallocate(hoff)
    if (allocated(iindx)) deallocate(iindx)
    if (allocated(jindx)) deallocate(jindx)
    
    return
    
  end subroutine deallocate_hamiltonian
  
!#######################################################################
! matxvec: calculates v2 = -H.v1 (Note the minus sign!)
!          (used in imaginary time wavepacket propagations)
!#######################################################################

  subroutine matxvec(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
  
    implicit none

    integer                    :: matdim
    integer*8                  :: noffdiag
    real(d), dimension(matdim) :: v1,v2

    if (hincore) then
       call matxvec_incore(matdim,noffdiag,v1,v2)
    else
       call matxvec_ext(matdim,noffdiag,v1,v2)
    endif
    
    return
    
  end subroutine matxvec

!#######################################################################
! matxvec_incore: calculates v2 = -H.v1 using a Hamiltonian stored
!                 in-core
!#######################################################################

  subroutine matxvec_incore(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
    
    implicit none

    integer                                :: matdim,maxbl,nrec,nlim,&
                                              ion,i,j,k,l,nthreads,tid
    integer*8                              :: noffdiag,npt
    integer*8, dimension(:,:), allocatable :: irange
    real(d), dimension(matdim)             :: v1,v2
    real(d), dimension(:,:), allocatable   :: tmpvec
    character(len=70)                      :: fileon,fileoff

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
    npt=int(floor(real(noffdiag)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffdiag

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=-hon(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=czero

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(iindx(k),tid)=tmpvec(iindx(k),tid)&
               -hoff(k)*v1(jindx(k))
          tmpvec(jindx(k),tid)=tmpvec(jindx(k),tid)&
               -hoff(k)*v1(iindx(k))
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
    
  end subroutine matxvec_incore
  
!#######################################################################
! matxvec_ext: calculates v2 = -H.v1 using an externally stored
!              Hamiltonian
!#######################################################################
  
  subroutine matxvec_ext(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
  
    implicit none

    integer                            :: matdim,maxbl,nrec,nlim,&
                                          ion,ioff,k,l
    integer*8                          :: noffdiag
    integer, dimension(:), allocatable :: indxi,indxj
    real(d), dimension(matdim)         :: v1,v2
    real(d), dimension(:), allocatable :: hii,hij
    character(len=70)                  :: fileon,fileoff

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
    ! Filenames
    if (hamflag.eq.'i') then
       fileon='SCRATCH/hmlt.diai'
       fileoff='SCRATCH/hmlt.offi'
    else if (hamflag.eq.'f') then
       fileon='SCRATCH/hmlt.diac'
       fileoff='SCRATCH/hmlt.offc'
    endif

    ! On-diagonal elements
    call freeunit(ion)
    open(ion,file=fileon,status='old',access='sequential',&
         form='unformatted')
    
    ! Off-diagonal elements
    call freeunit(ioff)
    open(ioff,file=fileoff,status='old',access='sequential',&
         form='unformatted')
  
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    rewind(ion)
    read(ion) maxbl,nrec

    allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
    allocate(hii(matdim))

!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    read(ion) hii
  
    do k=1,matdim
       v2(k)=-hii(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
    rewind(ioff)

    do k=1,nrec
       read(ioff) hij(:),indxi(:),indxj(:),nlim
       do l=1,nlim
          v2(indxi(l))=v2(indxi(l))-hij(l)*v1(indxj(l))
          v2(indxj(l))=v2(indxj(l))-hij(l)*v1(indxi(l))
       enddo
    enddo
  
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(hii,hij,indxi,indxj)
  
!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(ion)
    close(ioff)
    
    return
    
  end subroutine matxvec_ext

!#######################################################################
! matxvec_treal: calculates v2 = -iH.v1 
!                (used in real time wavepacket propagations)
!#######################################################################

  subroutine matxvec_treal(matdim,noffdiag,v1,v2)

    use constants

    implicit none
    
    integer                       :: matdim
    integer*8                     :: noffdiag
    complex(d), dimension(matdim) :: v1,v2

    if (hincore) then
       call matxvec_treal_incore(matdim,noffdiag,v1,v2)
    else
       call matxvec_treal_ext(matdim,noffdiag,v1,v2)
    endif
       
    return

  end subroutine matxvec_treal

!#######################################################################
! matxvec_treal_incore: calculates v2 = -iH.v1 using a Hamiltonian
!                       stored in-core
!#######################################################################

  subroutine matxvec_treal_incore(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
  
    implicit none

    integer                                 :: matdim,maxbl,nrec,nlim,&
                                               ion,i,j,k,l,nthreads,tid
    integer*8                               :: noffdiag,npt
    integer*8, dimension(:,:), allocatable  :: irange
    complex(d), dimension(matdim)           :: v1,v2
    complex(d), dimension(:,:), allocatable :: tmpvec
    character(len=70)                       :: fileon,fileoff

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
    npt=int(floor(real(noffdiag)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffdiag

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=-ci*hon(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=czero

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(iindx(k),tid)=tmpvec(iindx(k),tid)&
               -ci*hoff(k)*v1(jindx(k))
          tmpvec(jindx(k),tid)=tmpvec(jindx(k),tid)&
               -ci*hoff(k)*v1(iindx(k))
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

  end subroutine matxvec_treal_incore
    
!#######################################################################
! matxvec_treal_ext: calculates v2 = -iH.v1 using an externally
!                    stored Hamiltonian
!#######################################################################
  
  subroutine matxvec_treal_ext(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
  
    implicit none

    integer                                 :: matdim,maxbl,nrec,nlim,&
                                               ion,i,j,k,l,nthreads,tid
    integer*8                               :: noffdiag
    integer, dimension(:), allocatable      :: indxi,indxj,ioff
    real(d), dimension(:), allocatable      :: hii,hij
    complex(d), dimension(matdim)           :: v1,v2
    complex(d), dimension(:,:), allocatable :: tmpvec
    character(len=70)                       :: fileon,fileoff

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
    ! On-diagonal file
    if (hamflag.eq.'i') then
       fileon='SCRATCH/hmlt.diai'
    else
       fileon='SCRATCH/hmlt.diac'
    endif
    call freeunit(ion)
    open(ion,file=fileon,status='old',access='sequential',&
         form='unformatted')
    
    ! Off-diagonal files
    allocate(ioff(nthreads))
    do i=1,nthreads
       call freeunit(ioff(i))
       if (hamflag.eq.'i') then
          fileoff='SCRATCH/hmlt.offi.'
       else
          fileoff='SCRATCH/hmlt.offc.'
       endif
       k=len_trim(fileoff)+1
       if (i.lt.10) then
          write(fileoff(k:k),'(i1)') i
       else
          write(fileoff(k:k+1),'(i2)') i
       endif
       open(unit=ioff(i),file=fileoff,status='unknown',&
            access='sequential',form='unformatted')
    enddo

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    rewind(ion)
    read(ion) maxbl,nrec

    allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
    allocate(hii(matdim))
    allocate(tmpvec(matdim,nthreads))

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    read(ion) hii
  
    do k=1,matdim
       v2(k)=-ci*hii(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=czero

    !$omp parallel do &
    !$omp& private(i,k,l,tid,hij,indxi,indxj,nlim) &
    !$omp& shared(tmpvec,v1,nrec_omp)
    do i=1,nthreads
       tid=1+omp_get_thread_num()
       do k=1,nrec_omp(i)
          read(ioff(i)) hij(:),indxi(:),indxj(:),nlim
          do l=1,nlim
             tmpvec(indxi(l),tid)=tmpvec(indxi(l),tid) &
                  -ci*hij(l)*v1(indxj(l))
             tmpvec(indxj(l),tid)=tmpvec(indxj(l),tid) &
                  -ci*hij(l)*v1(indxi(l))
          enddo
       enddo
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(ion)

    do i=1,nthreads
       close(ioff(i))
    enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(hii,hij,indxi,indxj)
    deallocate(tmpvec)
    deallocate(ioff)
    
    return

  end subroutine matxvec_treal_ext

!#######################################################################

  subroutine matxvec_treal_laser(time,matdim,noffdiag,v1,v2)

    use constants

    implicit none
    
    integer                       :: matdim
    integer*8                     :: noffdiag
    real(d)                       :: time
    complex(d), dimension(matdim) :: v1,v2

    print*,"WRITE THE MATXVEC_TREAL_LASER CODE!"
    stop

!----------------------------------------------------------------------
! (1) -i * H * v1
!----------------------------------------------------------------------
!    if (hincore) then
!       call matxvec_treal_incore(matdim,noffdiag,v1,v2)
!    else
!       call matxvec_treal_ext(matdim,noffdiag,v1,v2)
!    endif

!----------------------------------------------------------------------
! (2) -i * -mu.E(t) * v1
!----------------------------------------------------------------------
    
    return
    
  end subroutine matxvec_treal_laser
    
!#######################################################################

end module tdsemod
