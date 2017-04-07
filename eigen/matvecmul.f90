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

subroutine matxvec_treal(matdim,v1,v2)

  use constants
  use parameters
  use iomod
  
  implicit none

  integer                            :: matdim,maxbl,nrec,nlim,&
                                        ion,ioff,k,l
  integer, dimension(:), allocatable :: indxi,indxj
  complex*16, dimension(matdim)      :: v1,v2
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
     v2(k)=-ci*hii(k)*v1(k)
  enddo

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
  rewind(ioff)

  do k=1,nrec
     read(ioff) hij(:),indxi(:),indxj(:),nlim
     do l=1,nlim
        v2(indxi(l))=v2(indxi(l))-ci*hij(l)*v1(indxj(l))
        v2(indxj(l))=v2(indxj(l))-ci*hij(l)*v1(indxi(l))
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

end subroutine matxvec_treal
