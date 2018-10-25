module cheby2specmod

  use constants

  save

  integer                :: order,epoints,kfinal
  real(dp), dimension(2) :: bounds
  real(dp)               :: emin,emax
  real(dp), allocatable  :: auto(:)
  real(dp), parameter    :: eh2ev=27.2113845d0
  real(dp)               :: convfac
  real(dp), allocatable  :: a(:),b(:),c(:),d(:)
  real(dp)               :: tau
  logical                :: lau
  logical                :: lpade
  
end module cheby2specmod

!######################################################################

program cheby2spec

  use cheby2specmod

  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdcheby2specinp

!----------------------------------------------------------------------
! Read the Chebyshev order domain autocorrelation function file
!----------------------------------------------------------------------
  call rdautofile

!----------------------------------------------------------------------
! Calculate the spectrum from the cosine Fourier transform of the
! Chebyshev order-domain autocorrelation function
!----------------------------------------------------------------------
  call calc_spectrum
  
contains

!######################################################################

  subroutine rdcheby2specinp

    use constants
    use iomod
    use cheby2specmod
    
    implicit none

    integer           :: i
    character(len=30) :: string1,string2

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
    ! Energy range
    emin=-1.0d0
    emax=-1.0d0

    ! No. energy points
    epoints=1000

    ! Maximum order
    kfinal=1e+6

    ! Energies in a.u.
    lau=.false.

    ! Pade approximant of the cosine transform
    lpade=.false.

    ! Damping factor
    tau=0.0d0
    
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    if (iargc().ne.0) then

       i=1
5      call getarg(i,string1)

       if (string1.eq.'-e') then
          ! Energy range in eV
          i=i+1
          call getarg(i,string2)
          read(string2,*) emin
          i=i+1
          call getarg(i,string2)
          read(string2,*) emax         

       else if (string1.eq.'-np') then
          ! No. energy points
          i=i+1
          call getarg(i,string2)
          read(string2,*) epoints

       else if (string1.eq.'-kf') then
          ! Maximum order
          i=i+1
          call getarg(i,string2)
          read(string2,*) kfinal
          if (mod(kfinal,2).ne.0) kfinal=kfinal-1

       else if (string1.eq.'-au') then
          ! Energies are given in a.u.
          lau=.true.

       else if (string1.eq.'-pade') then
          ! Pade approximant of the Fourier transform
          lpade=.true.
          ! Damping factor
          i=i+1
          call getarg(i,string2)
          read(string2,*) tau
          
       else
          errmsg='Unknown keyword: '//trim(string1)
          call error_control
       endif
       
       if (i.lt.iargc()) then
          i=i+1
          goto 5
       endif

    endif
       
!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (emin.eq.-1.0d0) then
       errmsg='The energy range has not been given'
       call error_control
    endif

!----------------------------------------------------------------------
! Set the energy conversion factor
!----------------------------------------------------------------------
    if (lau) then
       convfac=1.0d0
    else
       convfac=eh2ev
    endif

!----------------------------------------------------------------------
! Conversion of the spectrum bounds
!----------------------------------------------------------------------    
    emin=emin/convfac
    emax=emax/convfac
    
    return
    
  end subroutine rdcheby2specinp
    
!######################################################################

  subroutine rdautofile

    use constants
    use iomod
    use parsemod
    use cheby2specmod
    
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

    ! If we are calculating the pade approximant of the cosine
    ! transform, then make sure that order is even.
    if (lpade.and.mod(order,2).ne.0) order=order-1
    
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

  subroutine calc_spectrum

    use constants
    use iomod
    use cheby2specmod

    implicit none

    integer                        :: i,k,unit
    real(dp)                       :: e,escaled,theta
    real(dp), dimension(0:4)       :: spec
    real(dp), dimension(0:order,3) :: window
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='chebyspec.dat',form='formatted',status='unknown')

!----------------------------------------------------------------------
! File header
!----------------------------------------------------------------------
    call write_header(unit)

!----------------------------------------------------------------------
! Precompute window function values
!----------------------------------------------------------------------
    do k=0,order
       ! Gaussian window
       window(k,1)=gausswindow(k)
       ! cos^2 window
       window(k,2)=cos2window(k)
       ! Jackson window
       window(k,3)=jacksonwindow(k)
    enddo

!----------------------------------------------------------------------
! Calculate the Pade approximant coefficients if needed
!----------------------------------------------------------------------
    if (lpade) call pade_coeff_cheby
    
!----------------------------------------------------------------------
! Calculate and output the spectrum
!----------------------------------------------------------------------
    ! Loop over energies
    do i=1,epoints

       ! Unscaled energy
       e=(i-1)*(emax-emin)/(epoints-1)+emin

       ! Scaled energy
       escaled=scalefunc(e)

       ! Skip if the scaled energy in not in the interval [-1,1]
       if (escaled.lt.-1.0d0.or.escaled.gt.1.0d0) cycle

       ! Calculate the angle theta
       theta=acos(escaled)

       ! Calculate the spectrum in the angle domain
       spec(0)=auto(0)
       spec(1)=spec(0)*window(0,1)
       spec(2)=spec(0)*window(0,2)
       spec(3)=spec(0)*window(0,3)
       do k=1,order
          spec(0)=spec(0)+2.0d0*cos(k*theta)*auto(k)
          spec(1)=spec(1)+2.0d0*cos(k*theta)*auto(k)*window(k,1)
          spec(2)=spec(2)+2.0d0*cos(k*theta)*auto(k)*window(k,2)
          spec(3)=spec(3)+2.0d0*cos(k*theta)*auto(k)*window(k,3)
       enddo

       ! Pade approximant
       if (lpade) spec(4)=padespec_cheby(theta)
       
       ! Spectrum in the energy domain
       spec=spec/sin(theta)

       ! Prefactor
       spec=spec*e*2.0d0/3.0d0/pi/(order/2)
       
       ! Output the energy and spectrum values
       if (lpade) then
          write(unit,'(ES15.6,4(2x,ES15.6))') e*convfac,spec(4),&
               spec(1),spec(3),spec(0)
       else
          write(unit,'(ES15.6,3(2x,ES15.6))') e*convfac,spec(1),spec(3),&
               spec(0)
       endif
          
    enddo

!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine calc_spectrum

!######################################################################

  subroutine write_header(unit)

    use constants
    use cheby2specmod
    
    implicit none

    integer :: unit,i

    if (lpade) then

       write(unit,'(84a)') ('#',i=1,84)

       if (lau) then
          write(unit,'(a)') '#  Energy (au)      Pade             &
               Gaussian         Jackson          No'
       else
          write(unit,'(a)') '#  Energy (eV)      Pade             &
               Gaussian         Jackson          No'
       endif
          
          write(unit,'(a)') '#                   Approximant      &
               Window           Window           Window'
          
       write(unit,'(84a)') ('#',i=1,84)

    else

       write(unit,'(67a)') ('#',i=1,67)

       if (lau) then
          write(unit,'(a)') '#  Energy (au)      Gaussian         &
               Jackson          No'
       else
          write(unit,'(a)') '#  Energy (eV)      Gaussian         &
               Jackson          No'
       endif

       write(unit,'(a)') '#                   Window           &
            Window           Window'

       write(unit,'(67a)') ('#',i=1,67)

    endif
    
    return
    
  end subroutine write_header
    
!######################################################################
  
  function scalefunc(e) result(escale)

    use constants
    use cheby2specmod

    real(dp) :: e,escale

    escale=e-(0.5d0*(bounds(2)-bounds(1))+bounds(1))
    escale=escale/(bounds(2)-bounds(1))
    escale=2.0d0*escale
    
    return
    
  end function scalefunc

!######################################################################

  function cos2window(k) result(func)

    use constants
    use cheby2specmod
    
    implicit none

    integer  :: k
    real(dp) :: func

    func=cos(k*pi/(2.0d0*order))
    func=func**2

    return
    
  end function cos2window

!######################################################################

  function gausswindow(k) result(func)

    use constants
    use cheby2specmod
    
    implicit none

    integer  :: k
    real(dp) :: func

    func=exp(-(2.24d0*k/order)**2)

    return
    
  end function gausswindow

!######################################################################

  function jacksonwindow(k) result(func)

    use constants
    use cheby2specmod
    
    implicit none

    integer  :: k
    real(dp) :: func,alpha

    alpha=pi/(order+2)

    func=(order-k+2)*cos(k*alpha) &
         -cos((order+1)*alpha)*sin((order-k+2)*alpha)/sin(alpha)
    func=func/(order+2)

    return
    
  end function jacksonwindow

!######################################################################

  subroutine pade_coeff_cheby

    use constants
    use cheby2specmod
    
    implicit none

    integer               :: k,m,n
    integer               :: info
    integer, allocatable  :: ipiv(:)
    real(dp), allocatable :: gmat(:,:)
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    n=order/2
    
    allocate(a(0:n))
    a=0.0d0
    allocate(b(0:n))
    b=0.0d0
    allocate(c(0:order))
    c=0.0d0
    allocate(d(n))
    d=0.0d0
    allocate(gmat(n,n))
    gmat=0.0d0
    allocate(ipiv(n))
    ipiv=0
    
!----------------------------------------------------------------------
! c-vector
!----------------------------------------------------------------------
    c(0)=auto(0)
    do k=1,order
       c(k)=2.0d0*auto(k)*exp(-k/tau)
    enddo

!----------------------------------------------------------------------
! d-vector
!----------------------------------------------------------------------
    do k=1,n
       d(k)=-c(n+k)
    enddo
    
!----------------------------------------------------------------------
! G-matrix
!----------------------------------------------------------------------
    do k=1,n
       do m=1,n
          gmat(k,m)=c(n-m+k)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the b-coefficients by colving the system of linear
! equations G.b = d
! N.B., we use the usual convention for diagonal Pade approximant
! schemes and set b0=1
!----------------------------------------------------------------------
    call dgetrf(n,n,gmat,n,ipiv,info)

    if (info.ne.0) then
       write(6,'(/,2x,a,/)') 'LU decomposition of the G-matrix &
            failed in subroutine pade_coeff_cheby'
       stop
    endif

!----------------------------------------------------------------------
! Calculate the b-coefficients
!----------------------------------------------------------------------
    b(0)=1.0d0

    b(1:n)=d(1:n)
    call dgetrs('N',n,1,gmat,n,ipiv,b(1:n),n,info)

    if (info.ne.0) then
       write(6,'(/,2x,a,/)') 'Failed call to dgetrs in subroutine &
            pade_coeff_cheby'
       stop
    endif

!----------------------------------------------------------------------
! Calculate the a-coefficients
!----------------------------------------------------------------------
    a=0.0d0

    a(0)=c(0)
    
    do k=1,n
       do m=0,k
          a(k)=a(k)+b(m)*c(k-m)
       enddo
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(gmat)
    deallocate(c)
    deallocate(d)
    deallocate(ipiv)
    
    return
    
  end subroutine pade_coeff_cheby

!######################################################################

  function padespec_cheby(theta) result(func)

    use constants
    use cheby2specmod

    implicit none

    integer     :: k
    real(dp)    :: theta,func
    complex(dp) :: z,numer,denom,val1,val2

    numer=czero
    denom=czero

    z=exp(ci*theta)
    
    do k=0,order/2
       numer=numer+a(k)*(z**k)
       denom=denom+b(k)*(z**k)
    enddo

    val1=numer/denom

    func=real(val1)

    return
    
  end function padespec_cheby
  
!######################################################################
  
end program cheby2spec
