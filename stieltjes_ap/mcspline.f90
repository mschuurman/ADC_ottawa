  program mcspline
    
    implicit none

    integer                             :: ndat,nint
    real*8, dimension(:), allocatable   :: dat,x,deriv,s,dx
    real*8, dimension(:,:), allocatable :: coeff
    character(len=80)                   :: datfile

!-----------------------------------------------------------------------
! Determine the name of the data file
!-----------------------------------------------------------------------
    datfile=''
    call getarg(1,datfile)

    if (datfile.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the data file has not been &
            given'
       STOP
    endif

!-----------------------------------------------------------------------
! Read the data file
!-----------------------------------------------------------------------
    call rddat(datfile,ndat,dat,x)

!-----------------------------------------------------------------------
! Perform the monotonicity-constrained cubic Hermite interpolation
!-----------------------------------------------------------------------
    nint=ndat-1
    allocate(deriv(ndat))
    allocate(s(nint))
    allocate(dx(nint))
    allocate(coeff(4,nint))
    call interpolate(ndat,nint,dat,x,deriv,coeff,s,dx)

    STOP

  contains

!#######################################################################

    subroutine rddat(datfile,ndat,dat,x)

      use parsemod

      implicit none

      integer                           :: unit,ndat,n
      real*8, dimension(:), allocatable :: dat,x
      character(len=80)                 :: datfile,atmp

!-----------------------------------------------------------------------
! Open the data file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=datfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the data file
!-----------------------------------------------------------------------
! (1) Determine the no. points and allocate arrays 
      ndat=0
5     continue
      read(unit,*,end=10) atmp
      if (atmp.ne.'') ndat=ndat+1
      goto 5

10    continue
      allocate(dat(ndat))
      allocate(x(ndat))
      dat=0.0d0
      x=0.0d0

! (2) Read the data
      rewind(unit)
      n=0
15    continue
      read(unit,'(a)',end=20) atmp
      if (atmp.ne.'') then
         n=n+1
         read(atmp,*) x(n),dat(n)
      endif
      goto 15

20    continue

!-----------------------------------------------------------------------
! Close the data file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rddat

!#######################################################################

    subroutine interpolate(ndat,nint,dat,x,deriv,coeff,s,dx)

      implicit none

      integer                   :: ndat,nint
      real*8, dimension(ndat)   :: dat,x,deriv
      real*8, dimension(nint)   :: s,dx
      real*8, dimension(4,nint) :: coeff

!-----------------------------------------------------------------------
! Compute the approximate derivatives as well as the slopes
!-----------------------------------------------------------------------
      call getderiv(ndat,nint,dat,x,deriv,s,dx)

!-----------------------------------------------------------------------
! Constrain the derivatives to the Boor-Swartz monotonicity limit
!-----------------------------------------------------------------------
      call constrain_deriv(ndat,nint,deriv,s)

!-----------------------------------------------------------------------
! Compute the coefficients
!-----------------------------------------------------------------------
      call getcoeff(ndat,nint,dat,x,deriv,coeff,s,dx)

!-----------------------------------------------------------------------
! Perform the interpolation
!-----------------------------------------------------------------------
      call calcp(ndat,nint,x,coeff,dx)

      return

    end subroutine interpolate

!#######################################################################

    subroutine getderiv(ndat,nint,dat,x,deriv,s,dx)

      implicit none

      integer                 :: ndat,nint,i
      real*8, dimension(ndat) :: dat,x,deriv
      real*8, dimension(nint) :: s,dx

      ! Interval lengths
      do i=1,nint
         dx(i)=x(i+1)-x(i)
      enddo

      ! Derivatives
      deriv=0.0d0
      do i=1,ndat-1
         deriv(i)=(dat(i+1)-dat(i))/dx(i)
      enddo
      deriv(ndat)=deriv(ndat-1)

      ! Slopes within the intervals
      do i=1,nint
         s(i)=(dat(i+1)-dat(i))/dx(i)
      enddo

      return

    end subroutine getderiv

!#######################################################################

    subroutine constrain_deriv(ndat,nint,deriv,s)

      implicit none

      integer                 :: ndat,nint,i
      real*8, dimension(ndat) :: deriv,smin,smax
      real*8, dimension(nint) :: s
      real*8                  :: ftmp

      do i=2,ndat-1
         smin(i)=min(s(i-1),s(i+1))
         smax(i)=max(s(i-1),s(i+1))
      enddo

      smin(1)=min(s(1),s(2))
      smin(ndat)=min(s(ndat-1),s(ndat))

      smax(1)=max(s(1),s(2))
      smax(ndat)=max(s(ndat-1),s(ndat))

      do i=2,ndat-1
         if (smin(i).gt.0) then
            ftmp=max(0.0d0,deriv(i))
            deriv(i)=min(ftmp,3.0d0*smin(i))
         else if (smax(i).lt.0) then
            ftmp=min(0.0d0,deriv(i))
            deriv(i)=max(ftmp,3.0d0*smax(i))
         else if (s(i-1)*s(i).le.0) then
            deriv(i)=0.0d0
         endif
      enddo

      return

    end subroutine constrain_deriv

!#######################################################################

    subroutine getcoeff(ndat,nint,dat,x,deriv,coeff,s,dx)
      
      implicit none
      
      integer                   :: ndat,nint,i
      real*8, dimension(ndat)   :: dat,x,deriv
      real*8, dimension(nint)   :: s,dx
      real*8, dimension(4,nint) :: coeff

      ! Loop over intervals
      do i=1,nint

         ! C1
         coeff(1,i)=dat(i)
         
         ! C2
         coeff(2,i)=deriv(i)

         ! C3
         coeff(3,i)=3.0d0*s(i)-deriv(i+1)-2.0d0*deriv(i)
         coeff(3,i)=coeff(3,i)/dx(i)

         ! C4
         coeff(4,i)=-(2.0d0*s(i)-deriv(i+1)-deriv(i))
         coeff(4,i)=coeff(4,i)/dx(i)**2

      enddo

      return

    end subroutine getcoeff

!#######################################################################

    subroutine calcp(ndat,nint,x,coeff,dx)

      implicit none

      integer                   :: ndat,nint,k,iint,unit
      real*8, dimension(ndat)   :: x
      real*8, dimension(nint)   :: dx
      real*8, dimension(4,nint) :: coeff
      real*8                    :: diff,p,xcurr

!-----------------------------------------------------------------------
! Determine the step length
!-----------------------------------------------------------------------
      diff=(x(ndat)-x(1))/1000.0d0

!-----------------------------------------------------------------------
! (1) Cumulative oscillator stength distribution
!-----------------------------------------------------------------------
      unit=20      
      open(unit,file='cosc_int.dat',form='formatted',status='unknown')

      do k=1,1001
         
         ! Set the current value of x
         xcurr=x(1)+(k-1)*diff

         ! Determine which interval we are in
         call getiint(iint,x,xcurr,ndat,nint)

         ! Compute and output the value of the interpolant at the 
         ! current value of x
         p=coeff(1,iint) &
              + coeff(2,iint)*(xcurr-x(iint)) &
              + coeff(3,iint)*(xcurr-x(iint))**2 &
              + coeff(4,iint)*(xcurr-x(iint))**3

        write(unit,*) xcurr,p        

      enddo

      close(unit)

!-----------------------------------------------------------------------
! (2) Oscillator strengths
!-----------------------------------------------------------------------
      unit=20      
      open(unit,file='osc_int.dat',form='formatted',status='unknown')

      do k=1,1001
         
         ! Set the current value of x
         xcurr=x(1)+(k-1)*diff

         ! Determine which interval we are in
         call getiint(iint,x,xcurr,ndat,nint)

         ! Compute and output the value of the interpolant at the 
         ! current value of x
         p=coeff(2,iint)&
              + 2.0d0*coeff(3,iint)*(xcurr-x(iint)) &
              + 3.0d0*coeff(4,iint)*(xcurr-x(iint))**2

        write(unit,*) xcurr,p        

      enddo

      close(unit)

      return

    end subroutine calcp

!#######################################################################

    subroutine getiint(iint,x,xcurr,ndat,nint)

      implicit none

      integer                 :: iint,ndat,nint,i
      real*8, dimension(ndat) :: x
      real*8                  :: xcurr

      do i=1,nint
         if (xcurr.ge.x(i).and.xcurr.lt.x(i+1)) then
            iint=i
            exit
         endif
      enddo

      return

    end subroutine getiint

!#######################################################################

  end program mcspline
