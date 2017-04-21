module auto2specmod

  use constants

  save

  integer                               :: maxtp,iauto,epoints
  real(d)                               :: dt,t0,emin,emax,tau,sigma,&
                                           dele,a0,b0
  real(d), dimension(:,:), allocatable  :: sp
  real(d), parameter                    :: eh2ev=27.2113845d0
  real(d), parameter                    :: fs2au=41.3413745758d0
  complex(d), dimension(:), allocatable :: auto
  character(len=70)                     :: outfile
  
end module auto2specmod

!######################################################################

program auto2spec

  use auto2specmod
  
  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdauto2specinp
  
!----------------------------------------------------------------------
! Read the autocorrelation function from file
!----------------------------------------------------------------------
  call rdauto
  
!----------------------------------------------------------------------
! Calculate the spectrum
!----------------------------------------------------------------------
  call calc_spectrum

!----------------------------------------------------------------------
! Write the spectrum to file
!----------------------------------------------------------------------
  call wrspectrum
  
contains

!######################################################################

  subroutine rdauto2specinp

    use auto2specmod
    use constants
    use iomod
    
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

    ! Damping time
    tau=-1.0d0
    
    ! Name of the output file
    outfile=''

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    i=1
5   call getarg(i,string1)

    if (string1.eq.'-sigma') then
       ! FWHM in eV
       i=i+1
       call getarg(i,string2)
       read(string2,*) sigma
       ! Convert to tau in au
       sigma=sigma/eh2ev
       tau=2.0d0/sigma
    else if (string1.eq.'-e') then
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
    else if (string1.eq.'-o') then
       ! Name of the output file
       i=i+1
       call getarg(i,outfile)
    else
       errmsg='Unknown keyword: '//trim(string1)
       call error_control
    endif
       
    if (i.lt.iargc()) then
       i=i+1
       goto 5
    endif

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (emin.eq.-1.0d0) then
       errmsg='The energy range has not been given'
       call error_control
    endif

    if (tau.eq.-1.0d0) then
       errmsg='The spectral broadening has not been given'
       call error_control
    endif
    
    if (outfile.eq.'') then
       outfile='spectrum.dat'
    endif

    return
    
  end subroutine rdauto2specinp
    
!######################################################################

  subroutine rdauto

    use auto2specmod
    use constants
    use iomod
    
    implicit none

    integer :: iin,i
    real(d) :: t,re,im
    
!----------------------------------------------------------------------
! Open the autocorrelation function file
!----------------------------------------------------------------------
    call freeunit(iin)
    open(iin,file='auto',form='formatted',status='old')

!----------------------------------------------------------------------
! Determine the no. timesteps and the time interval, and allocate
! arrays
!----------------------------------------------------------------------
    ! Read past the comment line
    read(iin,*)

    ! Determine the timestep
    read(iin,*) t0
    read(iin,*) t
    t0=t0*fs2au
    t=t*fs2au
    dt=t-t0

    ! Determine the number of timesteps
    maxtp=2
5   read(iin,*,end=10)
    maxtp=maxtp+1
    goto 5
10 continue
    iauto=maxtp-1

    ! Allocate arrays
    allocate(auto(maxtp))
    auto=czero

!----------------------------------------------------------------------
! Read the autocorrelation function
!----------------------------------------------------------------------
    rewind(iin)
    read(iin,*)
    read(iin,*) t,a0,b0

    do i=1,iauto
       read(iin,*) t,re,im
       auto(i)=dcmplx(re,im)
    enddo
       
!----------------------------------------------------------------------
! Close the autocorrelation function file
!----------------------------------------------------------------------
    close(iin)
    
    return

  end subroutine rdauto

!######################################################################

  subroutine calc_spectrum

    use auto2specmod
    
    implicit none

    integer :: i,j
    real(d) :: eau,t,cc,sum0,sum1,sum2,sum3,pia,gfac

!----------------------------------------------------------------------
! Allocate the spectrum arrays
!----------------------------------------------------------------------
    allocate(sp(0:epoints,0:3))
    
!----------------------------------------------------------------------
! Calculate the spectrum
!----------------------------------------------------------------------
    dele=(emax-emin)/dble(epoints)
    pia=pi/dble(2*maxtp+2)    
    gfac=dble(2.24d0/(maxtp+1))

    ! Loop over energy points
    do j=0,epoints

       ! Current energy
       eau=emin+j*dele

       sum0=0.5d0*a0
       sum1=0.5d0*a0
       sum2=0.5d0*a0
       sum3=0.5d0*a0
       
       ! Loop over timesteps
       do i=1,iauto

          t=i*dt

          cc=dble(exp(dcmplx(0.0d0,eau*t))*auto(i))*exp(-(t/tau))

          sum0=sum0+cc
          sum1=sum1+cc*cos(pia*i)
          sum2=sum2+cc*cos(pia*i)**2
          sum3=sum3+cc*exp(-(gfac*i)**2)

       enddo

       sp(j,0)=eau*sum0*dt/pi
       sp(j,1)=eau*sum1*dt/pi
       sp(j,2)=eau*sum2*dt/pi          
       sp(j,3)=eau*sum3*dt/pi
          
    enddo

    return
    
  end subroutine calc_spectrum
    
!######################################################################

  subroutine wrspectrum

    use auto2specmod
    use constants
    use iomod
    
    implicit none

    integer :: i,iout
    real(d) :: e

!----------------------------------------------------------------------
! Open the spectrum file
!----------------------------------------------------------------------
    call freeunit(iout)
    open(iout,file=outfile,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the spectra to file
!----------------------------------------------------------------------
    write(iout,'(a,x,F7.4)') '# FHWM (eV):',sigma*eh2ev
    write(iout,'(/,84a)') ('#',i=1,84)
    write(iout,'(a)') '#  Energy (eV)      cos^2 filter     &
         cos filter       Gaussian filter  no filter'
    write(iout,'(84a)') ('#',i=1,84)

    do i=0,epoints
       e=(emin+i*dele)*eh2ev
       write(iout,'(5(ES15.6,2x))') e,sp(i,2),sp(i,1),sp(i,3),sp(i,0)
    enddo
    
!----------------------------------------------------------------------
! Close the spectrum file
!----------------------------------------------------------------------
    close(iout)
    
    return
    
  end subroutine wrspectrum

!######################################################################
  
end program auto2spec
