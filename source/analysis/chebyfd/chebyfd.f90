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
  real(dp), dimension(2)  :: bounds
  
  ! DPSS information
  integer  :: nslepian
  integer  :: npts
  real(dp) :: fw
  
  ! Slepian functions
  real(dp), dimension(:,:), allocatable :: v

  ! Eigenvalues of the Slepians
  real(dp), dimension(:), allocatable :: lambda
  
  ! Energy interval
  real(dp) :: Ea,Eb

  ! Scaled energy interval
  real(dp) :: Eabar,Ebbar

  ! Expansion coefficients
  real(dp), dimension(:,:), allocatable :: fkn
  
  ! Output
  integer            :: idat
  character(len=120) :: adat

  ! Unit conversion factors
  real(dp), parameter :: eh2ev=27.2113845d0

  
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
! spectral bound, then reset it.
!
! Note that we cannot set Ea to bounds(1) as this will end up
! introducing a singularity into the sqrt(1-Ebar^2) term appearing
! in the equation for the coefficients...
!----------------------------------------------------------------------
  if (Ea.lt.bounds(1)) Ea=bounds(1)*1.001d0

!----------------------------------------------------------------------
! Calculate the scaled energy window bounds
!----------------------------------------------------------------------
  Eabar=scalefunc(Ea)
  Ebbar=scalefunc(Eb)

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

  function scalefunc(e) result(escale)

    use constants
    use chebyfdmod

    real(dp) :: e,escale

    escale=e-(0.5d0*(bounds(2)-bounds(1))+bounds(1))
    escale=escale/(bounds(2)-bounds(1))
    escale=2.0d0*escale
    
    return
    
  end function scalefunc
  
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
    npts=max(1000,order+1)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Slepians
    allocate(v(npts,nslepian))
    v=0.0d0

    ! Eigenvalues
    allocate(lambda(nslepian))
    lambda=0.0d0

!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
    call dpss(npts,fw,nslepian,v,lambda)

!----------------------------------------------------------------------
! Renormalisation on the interval [Eabar,Ebbar]
!----------------------------------------------------------------------
    v=v/((Ebbar-Eabar)/npts)
    
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
       write(ilog,'(2x,i3,2(2x,ES15.6))') i,lambda(i),1.0d0-lambda(i)
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
    use iomod
    
    implicit none

    integer                               :: n,j,k,unit
    real(dp)                              :: debar,ebar,theta
    real(dp), dimension(:,:), allocatable :: Tk,Tkw,val
    character(len=3)                      :: an
    character(len=20)                     :: filename
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the expansion coefficients...'
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(fkn(0:order,nslepian))
    fkn=0.0d0

    allocate(Tk(0:order,npts))
    Tkw=0.0d0
    
    allocate(Tkw(0:order,npts))
    Tkw=0.0d0

    allocate(val(npts,nslepian))
    val=0.0d0
    
!----------------------------------------------------------------------
! Calculation of the expansion coefficients using the trapezoidal
! rule
!----------------------------------------------------------------------
    ! Precalculate the values of the Chebyshev polynomials at the
    ! quadrature points weighted by 1/sqrt(1-Ebar)
    debar=((Ebbar-Eabar)/(npts-1))
    do j=1,npts
       ebar=Eabar+debar*(j-1)
       theta=acos(ebar)
       do k=0,order
          Tk(k,j)=cos(k*theta)
       enddo
       Tkw(:,j)=Tk(:,j)/sqrt(1.0d0-ebar**2)
    enddo
    
    ! Scale the first and last values of the weighted Chebyshev
    ! polynomials at the quadrature points s.t. we can take dot
    ! products to calculate the overlaps
    Tkw(:,1)=Tkw(:,1)/2.0d0
    Tkw(:,npts)=Tkw(:,npts)/2.0d0

    ! Calculate the expansion coefficients
    fkn=matmul(Tkw,v)
    
    ! Multiplication by DeltaE
    fkn=fkn*(Ebbar-Eabar)/npts

    ! Prefactors
    fkn(0,:)=fkn(0,:)/pi
    fkn(1:order,:)=fkn(1:order,:)*2.0d0/pi
    
!----------------------------------------------------------------------
! TEST: Do the calculated coefficients reproduce the DPSSs?
!----------------------------------------------------------------------
    ! Calculate the Chebyshev expansions of the DPSSs
    val=matmul(transpose(Tk),fkn)

    ! Output the Chebyshev expansions of the DPSSs
    call freeunit(unit)
    do n=1,nslepian
       write(an,'(i3)') n
       filename='dpss.'//trim(adjustl(an))//'.cheby.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       do j=1,npts
          write(unit,'(i5,2x,ES15.6)') j,val(j,n)*((Ebbar-Eabar)/npts)
       enddo
       close(unit)
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Tk)
    deallocate(Tkw)

    return
    
  end subroutine calc_expansion_coeffs
    
!######################################################################
  
end program chebyfd
