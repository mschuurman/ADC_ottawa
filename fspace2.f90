module fspace2
  
  use constants
  use parameters
  use get_matrix
  use get_moment
  use select_fano
  use filetools
  use misc
  use get_matrix_DIPOLE  
  use fspace

  implicit none
  
contains






!  subroutine get_tranmom_2(ndim,ndimf,kpq,kpqf,negvc,name,nstate,fen,travec,ndims)
    
!    integer, intent(in) :: ndim,negvc,ndims
!    integer, intent(out) :: nstate
!    real(d), dimension(:,:), allocatable :: arrd
!    real(d), dimension(negvc), intent(out) :: travec
!    real(d), intent(out) :: fen   
!    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
!    integer, intent(in) :: ndimf
!    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf

    
!    character(36), intent(in) :: name

!    integer :: i,j,num
!    real(d) :: enr,cntr
!    real(d), dimension(:), allocatable :: vec 
!    logical :: log1

!    INQUIRE(file=name,exist=log1)
    
!    if(.not. log1) then
!       write(6,*) 'The file ', name,' does not exist'
!       stop
!    end if
    
!    allocate(arrd(ndim,ndim))
!    allocate(vec(ndim))
!    nstate=0
!    cntr=0.0

!    if(dlim.lt.0.0) then

!    write(*,*) 'Take all states with 2h2p part greater than',-dlim

!    OPEN(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc

!       read(78,END=78) num,enr,vec(:)

!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do

!      if(i.eq.statenumber) then

!       if(cntr.lt.-dlim) then
!        fen=enr
!        call   get_fspace_adc2_DIPOLE_direct(ndim,ndimf,kpq,kpqf,vec(:),arrd,travec(:))
!!        travec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if
!      end if

!       cntr=0.0
!    end do
!78  CLOSE(78)
!    write(6,*) 'There is the ', statenumber,' energy and tran. moment'


!    else

!    write(*,*) 'Take all states with 1h1p part greater than',dlim
    
!    OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc
!
!       read(77,END=77) num,enr,vec(:)
!
!!       write(6,*) 'Read vector',num,enr
!
!
!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do
!

!
!      if(i.eq.statenumber) then
!
!
!       if(cntr.gt.dlim) then
!!        write(6,*) 'State',num,'taken'
!!        nstate=nstate+1   
!        fen=enr
!!        call   get_fspace_adc2_DIPOLE_direct(ndim,kpq,vec(:),arrd,travec(:))
!        call   get_fspace_adc2_DIPOLE_direct(ndim,ndimf,kpq,kpqf,vec(:),arrd,travec(:))
!!        travec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if
!      end if
!
!       cntr=0.0
!    end do
!77  CLOSE(77)
!    write(6,*) 'There is the ',statenumber,' energy and tran. moment'
!   
!   deallocate(arrd)
! 
!    end if 
    
!       
!  end subroutine get_tranmom_2
!





  subroutine get_tranmom_3(ndim,negvc,name,travec,nstate,fen,tmvec,ndims)
    
    integer, intent(in) :: ndim,negvc,ndims
    integer, intent(out) :: nstate
    real(d), dimension(ndim), intent(in) :: travec 
    real(d), dimension(negvc), intent(out) :: fen,tmvec
    
    character(36), intent(in) :: name

    integer :: i,j,num
    real(d) :: enr,cntr
    real(d), dimension(:), allocatable :: vec 
    logical :: log1

    INQUIRE(file=name,exist=log1)
    
    if(.not. log1) then
       write(6,*) 'The file ', name,' does not exist'
       stop
    end if
    
    allocate(vec(ndim))
    nstate=0
    cntr=0.0

    if(dlim.lt.0.0) then

    write(*,*) 'Take all states with 2h2p part greater than',-dlim

    OPEN(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
    do i=1,negvc

       read(78,END=78) num,enr,vec(:)

!       write(6,*) 'Read vector',num,enr
!      if (i.gt.statenumber) then
     

        do j=1,ndims
          cntr=cntr+vec(j)**2
        end do

       if(cntr.lt.-dlim) then
!        write(6,*) 'State',num,'taken'
        nstate=nstate+1
        fen(nstate)=enr
        tmvec(nstate)=tm(ndim,vec(:),travec(:))
       end if

!    end if

       cntr=0.0
    end do
78  CLOSE(78)
    write(6,*) 'There are ',nstate,' energies and tran. moments'


    else

    write(*,*) 'Take all states with 1h1p part greater than',dlim
    
    OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
    do i=1,negvc

       read(77,END=77) num,enr,vec(:)

!       write(6,*) 'Read vector',num,enr
!      if (i.gt.statenumber) then
 
        do j=1,ndims
          cntr=cntr+vec(j)**2
        end do

       if(cntr.gt.dlim) then
!        write(6,*) 'State',num,'taken'
        nstate=nstate+1   
        fen(nstate)=enr
        tmvec(nstate)=tm(ndim,vec(:),travec(:))
       end if


!     end if


       cntr=0.0
    end do
77  CLOSE(77)
    write(6,*) 'There are ',nstate,' energies and tran. moments'
    
    end if 
    
       
  end subroutine get_tranmom_3





!!$---------------------------------------------------------


end module fspace2
