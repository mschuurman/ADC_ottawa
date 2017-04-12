!#######################################################################
! matxvec: calculates v2 = -H.v1 (Note the minus sign!)
!          (used in imaginary time wavepacket propagations)
!#######################################################################

subroutine matxvec(matdim,v1,v2)

  use constants
  use parameters
  use iomod
  
  implicit none

  integer                            :: matdim,maxbl,nrec,nlim,&
                                        ion,ioff,k,l
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
  
end subroutine matxvec

!#######################################################################
! matxvec_treal: calculates v2 = -iH.v1 
!                (used in real time wavepacket propagations)
!#######################################################################

!subroutine matxvec_treal(matdim,v1,v2)
!
!  use constants
!  use parameters
!  use iomod
!
!  implicit none
!
!  integer                            :: matdim,maxbl,nrec,nlim,&
!                                        ion,ioff,k,l
!  integer, dimension(:), allocatable :: indxi,indxj
!  complex*16, dimension(matdim)      :: v1,v2
!  real(d), dimension(:), allocatable :: hii,hij
!  character(len=70)                  :: fileon,fileoff
!
!  ! Filenames
!  if (hamflag.eq.'i') then
!     fileon='SCRATCH/hmlt.diai'
!     fileoff='SCRATCH/hmlt.offi'
!  else if (hamflag.eq.'f') then
!     fileon='SCRATCH/hmlt.diac'
!     fileoff='SCRATCH/hmlt.offc'
!  endif
!
!  ! On-diagonal elements
!  call freeunit(ion)
!  open(ion,file=fileon,status='old',access='sequential',&
!       form='unformatted')
!
!  ! Off-diagonal elements
!  call freeunit(ioff)
!  open(ioff,file=fileoff,status='old',access='sequential',&
!       form='unformatted')
!
!!-----------------------------------------------------------------------
!! Allocate arrays
!!-----------------------------------------------------------------------
!  rewind(ion)
!  read(ion) maxbl,nrec
!
!  allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
!  allocate(hii(matdim))
!
!!-----------------------------------------------------------------------
!! On-diagonal elements
!!-----------------------------------------------------------------------
!  read(ion) hii
!  
!  do k=1,matdim
!     v2(k)=-ci*hii(k)*v1(k)
!  enddo
!
!!-----------------------------------------------------------------------
!! Off-diagonal elements
!!-----------------------------------------------------------------------
!  rewind(ioff)
!
!  do k=1,nrec
!     read(ioff) hij(:),indxi(:),indxj(:),nlim
!     do l=1,nlim
!        v2(indxi(l))=v2(indxi(l))-ci*hij(l)*v1(indxj(l))
!        v2(indxj(l))=v2(indxj(l))-ci*hij(l)*v1(indxi(l))
!     enddo
!  enddo
!
!!-----------------------------------------------------------------------
!! Deallocate arrays
!!-----------------------------------------------------------------------
!  deallocate(hii,hij,indxi,indxj)
!
!!-----------------------------------------------------------------------
!! Close files
!!-----------------------------------------------------------------------
!  close(ion)
!  close(ioff)
!
!  return
!
!end subroutine matxvec_treal

subroutine matxvec_treal(matdim,v1,v2)

  use constants
  use parameters
  use iomod
  use channels
  use omp_lib
  
  implicit none

  integer                                 :: matdim,maxbl,nrec,nlim,&
                                             ion,i,j,k,l,nthreads,tid
  integer, dimension(:), allocatable      :: indxi,indxj,ioff
  real(d), dimension(:), allocatable      :: hii,hij
  complex*16, dimension(matdim)           :: v1,v2
  complex*16, dimension(:,:), allocatable :: tmpvec
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

end subroutine matxvec_treal
