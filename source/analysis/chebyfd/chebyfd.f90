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
  integer                :: npts
  real(dp)               :: fw
  
  ! Slepian functions
  real(dp), dimension(:,:), allocatable :: v

  ! Eigenvalues of the Slepians
  real(dp), dimension(:), allocatable   :: lambda

  ! 1-eigenvalues
  real(dp), dimension(:), allocatable   :: theta
  
  ! Energy interval
  real(dp)               :: Ea,Eb

  ! Output
  integer                :: idat
  character(len=120)     :: adat

  ! Unit conversion factors
  real(dp), parameter    :: eh2ev=27.2113845d0

  
end module chebyfdmod

!######################################################################

program chebyfd

  use chebyfdmod
  
  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call openchebyfdfiles
  
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdchebyfdinp

!----------------------------------------------------------------------
! Read the Chebyshev order domain autocorrelation function file
!----------------------------------------------------------------------
  call rdautofile

!----------------------------------------------------------------------
! If the lower bound of the energy window is less than the lower
! spectral bound, then reset it
!----------------------------------------------------------------------
  if (Ea.lt.bounds(1)) Ea=bounds(1)

!----------------------------------------------------------------------
! Calculate the Slepian filter functions
!----------------------------------------------------------------------
  call calc_slepians

!----------------------------------------------------------------------
! Output some information about the Slepians to the log file
!----------------------------------------------------------------------
  call wrslepinfo

!----------------------------------------------------------------------
! Calculate the coefficients entering into the expansion of the
! Slepian filter functions with respect to the Chebyshev polynomials
!----------------------------------------------------------------------
  call calc_expansion_coeffs
  
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

    ! Maximum order
    kfinal=1e+6
    
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

5   call rdinp(iin)
    
    i=0
    if (.not.lend.and.keyword(i).ne.'end-input') then
10     continue
       i=i+1
       
       if (keyword(i).eq.'window') then
          if (keyword(i+1).eq.'=') then
             ! Lower bound
             i=i+2
             read(keyword(i),*) Ea
             ! Upper bound
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) Eb
             else
                errmsg='The upper energy bound has not been given &
                     with the window keyword'
                call error_control
             endif
             ! Conversion to a.u.
             Ea=Ea/eh2ev
             Eb=Eb/eh2ev
          else
             goto 100
          endif

       else if (keyword(i).eq.'slepians') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             ! No. slepian filter functions
             read(keyword(i),*) nslepian
             ! Time half-bandwidth product factor
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) fw
             else
                errmsg='The ime half-bandwidth product factor has &
                     not been given with the slepians keyword'
                call error_control
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'kf') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) kfinal
             if (mod(kfinal,2).ne.0) kfinal=kfinal-1
          else
             goto 100
          endif
          
       else
          ! Exit if the keyword is not recognised
          errmsg='Unknown keyword: '//trim(keyword(i))
          call error_control
       endif
       
       ! If there are more keywords to be read on the current line,
       ! then read them, else read the next line
       if (i.lt.inkw) then
          goto 10
       else
          goto 5
       endif
       
       ! Exit if a required argument has not been given with a keyword
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       call error_control
       
    endif

!----------------------------------------------------------------------
! Make sure that all required information has been given
!----------------------------------------------------------------------
    ! Energy window
    if (Ea.eq.-999.0d0) then
       errmsg='The energy window has not been given'
       call error_control
    endif

    ! Slepian filter functions
    if (nslepian.eq.0) then
       errmsg='The slepian filter functions have not been specified'
       call error_control
    endif
    
    return
    
  end subroutine rdchebyfdinp
    
!######################################################################

  subroutine rdautofile
  
    use constants
    use iomod
    use parsemod
    use chebyfdmod
    
    implicit none

    integer :: unit,k,itmp
    logical :: exists

!----------------------------------------------------------------------
! Make sure that the chebyauto file exists
!----------------------------------------------------------------------
    inquire(file='chebyauto',exist=exists)

    if (.not.exists) then
       errmsg='The chebyauto file could not be found'
       call error_control
    endif

!----------------------------------------------------------------------
! Open the chebyauto file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='chebyauto',form='formatted',status='old')

!----------------------------------------------------------------------
! Read the spectral bounds
!----------------------------------------------------------------------
    read(unit,'(21x,2(2x,E21.14))') bounds(1),bounds(2)

!----------------------------------------------------------------------
! Determine the order of the Chebyshev expansion and allocate the auto
! array
!----------------------------------------------------------------------
    read(unit,*)
    read(unit,*)    
    
    order=-1
5   read(unit,*,end=10)
    order=order+1
    goto 5

10  continue

    if (order.gt.kfinal) order=kfinal
    
    allocate(auto(0:order))
    auto=0.0d0

!----------------------------------------------------------------------
! Read the Chebyshev order domain autocorrelation function
!----------------------------------------------------------------------
    rewind(unit)

    do k=1,3
       read(unit,*)
    enddo

    do k=0,order
       read(unit,*) itmp,auto(k)
    enddo
    
!----------------------------------------------------------------------
! Close the chebyauto file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rdautofile

!######################################################################

  subroutine calc_slepians

    use chebyfdmod
    use dpssmt

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the DPSSs...'
    
!----------------------------------------------------------------------
! Set the number of points at which the DPSSs will be evaluated.
!
! ??? This has to be greater that the no. Chebyshev polynomials entering
! into the expansion of the Slepian filter functions in order to use
! Gauss-Chebyshev quadrature in the calculation of the expansion
! coefficients ???
!----------------------------------------------------------------------
    npts=order+1

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Slepians
    allocate(v(npts,nslepian))
    v=0.0d0

    ! Eigenvalues
    allocate(lambda(nslepian))
    lambda=0.0d0

    ! 1-eigenvalues
    allocate(theta(nslepian))
    theta=0.0d0

!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
    call dpss(npts,fw,nslepian,v,lambda,theta)
    
    return
    
  end subroutine calc_slepians
    
!######################################################################

  subroutine wrslepinfo

    use channels
    use chebyfdmod
    
    implicit none

    integer  :: i,j
    real(dp) :: ovrlp
    
!----------------------------------------------------------------------
! Output some information about the Slepians to file
!----------------------------------------------------------------------
    ! fw
    write(ilog,'(a,2x,F10.7)') '# Time half-bandwidth product:',fw
    
    ! Eigenvalues
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(a)') '  DPSS       lambda         1 - lambda'
    write(ilog,'(41a)') ('#',i=1,41)
    do i=1,nslepian
       write(ilog,'(2x,i3,2(2x,ES15.6))') i,lambda(i),theta(i)
    enddo

    ! Overlaps
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(a)') '  Overlaps'
    write(ilog,'(41a)') ('#',i=1,41)
    do i=1,nslepian-1
       do j=i,nslepian
          ovrlp=dpss_overlap(i,j)
          write(ilog,'(2(2x,i3),2x,ES15.6)') i,j,ovrlp
       enddo
    enddo
    
    return
    
  end subroutine wrslepinfo

!######################################################################

  function dpss_overlap(i,j) result(ovrlp)

    use chebyfdmod
    
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

  subroutine calc_expansion_coeffs

    use chebyfdmod

    implicit none

    integer :: n,k
    
!----------------------------------------------------------------------
! Calculation of the expansion coefficients using Gauss-Chebyshev
! quadrature
!----------------------------------------------------------------------
    ! Loop over Slepians
    do n=1,nslepian

       ! Loop over Chebyshev polynomials
       do k=0,order

          
          
       enddo
       
    enddo

    
    return
    
  end subroutine calc_expansion_coeffs
    
!######################################################################
  
end program chebyfd
