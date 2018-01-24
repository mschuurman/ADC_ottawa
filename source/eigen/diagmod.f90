module diagmod
  
  use constants
  use parameters
  use channels
  
  implicit none 

contains

!#######################################################################
! master_eig: Calls the main eigensoler routines
!#######################################################################

  subroutine master_eig(ndim,noffd,flag)
    
    use block_davidson
    use relaxmod
    
    implicit none
    
    integer, intent(in)      :: ndim
    integer*8, intent(in)    :: noffd
    integer                  :: es
    character(1), intent(in) :: flag

!-----------------------------------------------------------------------
! Set the Hamiltonian flag: 'i' <-> initial space
!                           'f' <-> final space
!-----------------------------------------------------------------------
    hamflag=flag

!-----------------------------------------------------------------------
! Set the eigensolver
!-----------------------------------------------------------------------
    if (hamflag.eq.'i') then
       es=solver
    else if (hamflag.eq.'f') then
       es=solver_f
    endif
      
!-----------------------------------------------------------------------
! Call to the main eigensolver routine
!-----------------------------------------------------------------------
    if (es.eq.1) then
       call davdiag_block(ndim,noffd)
    else if (es.eq.2) then
       call relaxation(ndim,noffd)
    endif
    
    return

  end subroutine master_eig

!#######################################################################
! readdavvc: Reads the Davidson vectors and energies from disk
!#######################################################################

  subroutine readdavvc(nvec,evec,rvec,flag,dim)
      
    use iomod
      
    implicit none

    integer, intent(in)                                :: nvec,dim
    double precision, dimension(dim,nvec), intent(out) :: rvec
    double precision, dimension(nvec), intent(out)     :: evec      
    integer                                            :: i,nr,nfl,count
    character(len=1)                                   :: flag
    character(len=36)                                  :: filename
      
!-----------------------------------------------------------------------
! Open the file containing the Davidson eigenpairs
!-----------------------------------------------------------------------
    call freeunit(nfl)
      
    if (flag.eq.'i') then
       filename=davname
    else if (flag.eq.'f') then
       filename=davname_f
    endif
    
    open(unit=nfl,file=filename,status='old',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Read the Davidson vectors and energies
!-----------------------------------------------------------------------
    count=0
    do i=1,nvec
       read(nfl,end=77) nr,evec(i),rvec(:,i)
       count=count+1
    end do
    
77  close(nfl)

!-----------------------------------------------------------------------
! Die if we haven't found the correct no. vectors
!-----------------------------------------------------------------------
    if (count.lt.nvec) then
       errmsg=''
       write(errmsg,'(2(x,a,x,i3),x,a)') 'Error:',nvec,&
            'vectors requested,',count,'vectors found on disk'
       call error_control
    endif

    return
      
  end subroutine readdavvc

!#######################################################################

end module diagmod
