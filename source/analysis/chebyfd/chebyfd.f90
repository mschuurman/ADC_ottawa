!######################################################################
! chebyfd: a low-storage Chebyshev filter diagonalisation program
!          based on the use of Slepian filter functions
!######################################################################

module chebyfdmod

  use constants

  implicit none
  
  save

  ! Order domain autocorrelation function
  integer                :: order,kfinal
  real(dp), allocatable  :: auto(:)

  ! Spectral bounds
  real(dp), dimension(2) :: bounds  
  
  ! DPSS information
  integer                :: nslepian
  integer                :: fw

  ! Energy interval
  real(dp), dimension(2) :: Ea,Eb
  
  ! Output
  integer                :: idat
  character(len=120)     :: adat
  
end module chebyfdmod

!######################################################################

program chebyfd

  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call openchebyfdfiles
  
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdchebyfdinp

contains
  
!######################################################################

  subroutine openchebyfdfiles

    use constants
    use channels
    use iomod
    use chebyfdmod
    
    implicit none

    integer :: k

!----------------------------------------------------------------------
! Get the name of the input file
!----------------------------------------------------------------------
    ain=''

    call getarg(1,ain)

    if (ain.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the input file has not &
            been given!'
       stop
    endif

    if (index(ain,'.inp').eq.0) ain=trim(ain)//'.inp'

!----------------------------------------------------------------------
! Log file name
!----------------------------------------------------------------------
    k=index(ain,'.inp')-1
    alog=ain(1:k)//'.log'

!----------------------------------------------------------------------
! Data file name
!----------------------------------------------------------------------
    k=index(ain,'.inp')-1
    adat=ain(1:k)//'.dat'

!----------------------------------------------------------------------
! Open the input file
!----------------------------------------------------------------------
    call freeunit(iin)
    open(iin,file=ain,form='formatted',status='old')

!----------------------------------------------------------------------
! Open the log file
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file=alog,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Open the data file
!----------------------------------------------------------------------
    call freeunit(idat)
    open(idat,file=adat,form='unformatted',status='unknown')
    
    return
    
  end subroutine openchebyfdfiles
    
!######################################################################

  subroutine rdchebyfdinp

    use constants
    use iomod
    use parsemod
    use channels
    use chebyfdmod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Number of Slepian filter functions
    nslepian=0

    ! Time half-bandwidth parameter
    fw=0.0d0

    ! Energy window
    Ea=-999.0d0
    Eb=-999.0d0

!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

5   call rdinp(iin)

    print*,keyword(1)
    stop
    
    return
    
  end subroutine rdchebyfdinp
    
!######################################################################
  
end program chebyfd
