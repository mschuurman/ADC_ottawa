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
  use channels

  implicit none
  
contains

!#######################################################################

  subroutine get_tranmom_3(ndim,negvc,name,travec,nstate,fen,tmvec,ndims)
    
    integer, intent(in)                     :: ndim,negvc,ndims
    integer, intent(out)                    :: nstate
    real(dp), dimension(ndim), intent(in)   :: travec 
    real(dp), dimension(negvc), intent(out) :: fen,tmvec
    character(36), intent(in)               :: name

    integer                                 :: i,j,num
    real(dp)                                :: enr,cntr
    real(dp), dimension(:), allocatable     :: vec 
    logical                                 :: log1

    inquire(file=name,exist=log1)
    
    if (.not.log1) then
       write(ilog,*) 'The file ', name,' does not exist'
       stop
    endif
    
    allocate(vec(ndim))
    nstate=0
    cntr=0.0d0

    open(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')

    do i=1,negvc
       read(78,end=78) num,enr,vec(:)
       nstate=nstate+1
       fen(nstate)=enr
       tmvec(nstate)=tm(ndim,vec(:),travec(:))
    enddo

78  close(78)
    write(ilog,*) 'There are ',nstate,' energies and tran. moments'
    
    return
       
  end subroutine get_tranmom_3

!#######################################################################

end module fspace2
