!######################################################################
! slepian: A program to calculate discrete prolate spheroidal
!          sequences. Makes use of code lifted from the Multitaper
!          Spectrum Estimation library.
!######################################################################

module global

  use constants
  
  implicit none

  ! Number of points in the series
  integer                               :: npts

  ! Time-bandwidth product
  real(dp)                              :: fw

  ! Desired number of Slepian functions
  integer                               :: nev

  ! Slepian functions
  real(dp), dimension(:,:), allocatable :: v

  ! Eigenvalues of the Slepians
  real(dp), dimension(:), allocatable   :: lambda

  ! 1-eigenvalues
  real(dp), dimension(:), allocatable   :: theta
  
  save
  
end module global

!######################################################################

program slepian

  use global
  use dpssmt
  
  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinp

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  call alloc_arr

!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------  
  call dpss(npts,fw,nev,v,lambda)

  theta=1.0d0-lambda
  
!----------------------------------------------------------------------
! Output the DPSSs
!----------------------------------------------------------------------
  call wrout

contains

!######################################################################
  
  subroutine rdinp

    use global
    
    implicit none

    character(len=60) :: string

!----------------------------------------------------------------------
! Die here if the incorrect no. arguments have been given
!----------------------------------------------------------------------
    if (iargc().ne.3) then
       write(6,'(/,2(2x,a,/))') 'Incorrect no. command line &
            arguments','Correct input: npts fw nev'
       stop
    endif
    
!----------------------------------------------------------------------
! Argument 1: no. points
!----------------------------------------------------------------------
    call getarg(1,string)
    read(string,*) npts

!----------------------------------------------------------------------    
! Argument 2: time-bandwidth product
!----------------------------------------------------------------------    
    call getarg(2,string)
    read(string,*) fw

!----------------------------------------------------------------------    
! Argument 3: number of Slepians
!----------------------------------------------------------------------    
    call getarg(3,string)
    read(string,*) nev
    
    return
    
  end subroutine rdinp

!######################################################################
  
  subroutine alloc_arr

    use global

    implicit none

    ! Slepians
    allocate(v(npts,nev))
    v=0.0d0

    ! Eigenvalues
    allocate(lambda(nev))
    lambda=0.0d0

    ! 1-eigenvalues
    allocate(theta(nev))
    theta=0.0d0
    
    return
    
  end subroutine alloc_arr
  
!######################################################################

  subroutine wrout

    use iomod
    use global
    
    implicit none

    integer           :: unit,i,j
    real(dp)          :: ovrlp
    character(len=60) :: filename
    character(len=5)  :: ai

!----------------------------------------------------------------------
! Output the Slepians to file
!----------------------------------------------------------------------
    call freeunit(unit)

    ! Loop over the Slepians
    do i=1,nev

       ! Open the output file
       write(ai,'(i5)') i
       filename='dpss.'//trim(adjustl(ai))//'.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       
       ! Write the output file
       do j=1,npts
          write(unit,*) j,v(j,i)
       enddo
       
       ! Close the output file
       close(unit)
       
    enddo

!----------------------------------------------------------------------
! Output some information about the Slepians to file
!----------------------------------------------------------------------
    ! Open file
    open(unit,file='dpss.info',form='formatted',status='unknown')

    ! fw
    write(unit,'(a,2x,F10.7)') '# Time half-bandwidth product:',fw
    
    ! Eigenvalues
    write(unit,'(/,41a)') ('#',i=1,41)
    write(unit,'(a)') '  DPSS       lambda         1 - lambda'
    write(unit,'(41a)') ('#',i=1,41)
    write(6,'(41a)') ('#',i=1,41)
    write(6,'(a)') '  DPSS       lambda         1 - lambda'
    write(6,'(41a)') ('#',i=1,41)
    do i=1,nev
       write(unit,'(2x,i3,2(2x,ES15.6))') i,lambda(i),theta(i)
       write(6,'(2x,i3,2(2x,ES15.6))') i,lambda(i),theta(i)
    enddo

    ! Overlaps
    write(unit,'(/,41a)') ('#',i=1,41)
    write(unit,'(a)') '  Overlaps'
    write(unit,'(41a)') ('#',i=1,41)
    do i=1,nev-1
       do j=i,nev
          ovrlp=dpss_overlap(i,j)
          write(unit,'(2(2x,i3),2x,ES15.6)') i,j,ovrlp
       enddo
    enddo
          
    ! Close file
    close(unit)
    
    return

  end subroutine wrout

!######################################################################

  function dpss_overlap(i,j) result(ovrlp)

    use global
    
    implicit none

    integer  :: i,j,n
    real(dp) :: ovrlp

    ovrlp=v(1,i)*v(1,j)/2.0d0
    do n=2,npts-1
       ovrlp=ovrlp+v(n,i)*v(n,j)
    enddo
    ovrlp=ovrlp+v(npts,i)*v(npts,j)/2.0d0
    
    return
    
  end function dpss_overlap
  
!######################################################################
  
end program slepian

