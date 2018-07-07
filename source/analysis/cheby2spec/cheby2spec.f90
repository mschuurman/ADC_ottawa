module cheby2specmod

  use constants

  save

  integer                :: order,epoints
  real(dp), dimension(2) :: bounds
  real(dp)               :: emin,emax
  real(dp)               :: gam
  real(dp), allocatable  :: auto(:)
  real(dp), parameter    :: eh2ev=27.2113845d0
  
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

    ! Gaussian window function width parameter
    gam=1e-3_dp
    
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
          ! Convert to au
          emin=emin/eh2ev
          emax=emax/eh2ev

       else if (string1.eq.'-np') then
          ! No. energy points
          i=i+1
          call getarg(i,string2)
          read(string2,*) epoints

       else if (string1.eq.'-gamma') then
          ! Gaussian window function width parameter
          i=i+1
          call getarg(i,string2)
          read(string2,*) gam
          
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

    integer  :: i,k,unit
    real(dp) :: e,escaled,theta,spec
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='chebyspec.dat',form='formatted',status='unknown')
    
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
       spec=auto(0)*gam/sqrt(pi)
       do k=1,order
          spec=spec+2.0d0*cos(k*theta)*auto(k)*gam/sqrt(pi)*exp(-gam**2*k**2/2.0d0)
       enddo
       spec=spec/pi
       
       ! Spectrum in the energy domain
       spec=spec/sin(theta)

       ! Prefactor
       spec=spec*e*2.0d0/3.0d0
       
       ! Output the energy and spectrum value
       write(unit,'(ES15.6,2x,ES15.6)') e*eh2ev,spec
       
    enddo

!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine calc_spectrum

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
  
end program cheby2spec
