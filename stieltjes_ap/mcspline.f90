  program mcspline
    
    use mcspmod

    implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
    call rdmcspinp

!-----------------------------------------------------------------------
! Read the data files
!-----------------------------------------------------------------------
    call rddat

!-----------------------------------------------------------------------
! Perform the monotonicity-constrained cubic Hermite interpolation
!-----------------------------------------------------------------------
    allocate(deriv(nfiles,maxdat))
    allocate(s(nfiles,maxintvl))
    allocate(dx(nfiles,maxintvl))
    allocate(coeff(nfiles,4,maxintvl))

    call interpolate

    STOP

  contains

!#######################################################################

    subroutine rdmcspinp

      use mcspmod
      use parsemod

      implicit none

      integer                           :: iin,i,k
      character(len=120)                :: inpfile,errmsg
      character(len=120), dimension(50) :: atmp

!-----------------------------------------------------------------------
! Determine the input file name
!-----------------------------------------------------------------------      
      inpfile=''
      call getarg(1,inpfile)

      if (inpfile.eq.'') then
         write(6,'(/,2x,a,/)') 'The input file name has not been given'
         STOP
      endif

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      nfiles=0
      erange=-999.0d0
      npoints=0

!-----------------------------------------------------------------------
! Open the input file
!-----------------------------------------------------------------------
      iin=20
      open(iin,file=inpfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(iin)
        
      i=0
      if (keyword(1).ne.'end-input') then
10       continue
         i=i+1
         
         if (keyword(i).eq.'files') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Determine the no. files
15             continue
               nfiles=nfiles+1
               atmp(nfiles)=keyword(i)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  goto 15
               endif
               ! Allocate and fill in the datfile array
               allocate(datfile(nfiles))
               datfile(:)=atmp(1:nfiles)               
            else
               goto 100
            endif

         else if (keyword(i).eq.'interval') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) erange(1)               
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) erange(2)
               endif
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'points') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) npoints
            else
               goto 100
            endif

         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            write(6,'(/,a,/)') trim(errmsg)
            STOP
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         write(6,'(/,a,/)') trim(errmsg)
         STOP

      endif

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
      if (erange(1).eq.-999.0d0) erange(1)=0.0d0
      if (erange(2).eq.-999.0d0) erange(2)=1000000.0d0

      if (nfiles.eq.0) then
         write(6,'(/,2x,a,/)') 'No data filenames have been given'
         STOP
      endif

      if (npoints.eq.0) npoints=1001

!-----------------------------------------------------------------------
! Close the input file
!-----------------------------------------------------------------------
      close(iin)

      return

    end subroutine rdmcspinp
      
!#######################################################################

    subroutine rddat

      use mcspmod
      use parsemod

      implicit none

      integer           :: unit,i,np
      character(len=80) :: atmp

      unit=20

!-----------------------------------------------------------------------
! Determine the maximum no. points and allocate arrays accordingly
!-----------------------------------------------------------------------
      maxdat=0
      
      ! Loop over data files
      do i=1,nfiles

         ! Open the current data file
         open(unit,file=datfile(i),form='formatted',status='old')

         ! Determine the no. points
         np=0
5        call rdinp(unit)
         if (.not.lend) then
            np=np+1
            goto 5
         endif

         ! If np>maxdat, then update maxdat
         if (np.gt.maxdat) maxdat=np

         ! Close the current data file
         close(unit)

      enddo

      maxintvl=maxdat-1

      allocate(dat(nfiles,maxdat))
      allocate(x(nfiles,maxdat))

      dat=0.0d0
      x=0.0d0

!-----------------------------------------------------------------------
! Read the data
!-----------------------------------------------------------------------
      allocate(ndat(nfiles))
      allocate(nintvl(nfiles))

      ! Loop over data files
      do i=1,nfiles

         ! Open the current data file
         open(unit,file=datfile(i),form='formatted',status='old')

         ! Read the current data file
         ndat(i)=0
10       call rdinp(unit)
         if (.not.lend) then
            ndat(i)=ndat(i)+1
            read(keyword(1),*) x(i,ndat(i))
            read(keyword(2),*) dat(i,ndat(i))
            goto 10
         endif
         nintvl(i)=ndat(i)-1

         ! Close the current data file
         close(unit)
         
      enddo

      return

    end subroutine rddat

!#######################################################################

    subroutine interpolate

      use mcspmod

      implicit none

!-----------------------------------------------------------------------
! Compute the approximate derivatives and the slopes (which in the 
! current, crude(ish) incarnation are equal)
!-----------------------------------------------------------------------
      call getderiv

!-----------------------------------------------------------------------
! Constrain the derivatives to the Boor-Swartz monotonicity limit
!-----------------------------------------------------------------------
      call constrain_deriv

!-----------------------------------------------------------------------
! Compute the coefficients
!-----------------------------------------------------------------------
      call getcoeff

!-----------------------------------------------------------------------
! Perform the interpolation
!-----------------------------------------------------------------------
      call calcp

      return

    end subroutine interpolate

!#######################################################################

    subroutine getderiv

      implicit none

      integer :: i,j

      ! Loop over data sets
      do i=1,nfiles

         ! Interval lengths for the current set
         do j=1,nintvl(i)
            dx(i,j)=x(i,j+1)-x(i,j)
         enddo

         ! Derivatives
         deriv(i,:)=0.0d0
         do j=1,nintvl(i)
            deriv(i,j)=(dat(i,j+1)-dat(i,j))/dx(i,j)
         enddo
         deriv(i,ndat(i))=deriv(i,ndat(i-1))
         
         ! Slopes within the intervals
         do j=1,nintvl(i)
            s(i,j)=(dat(i,j+1)-dat(i,j))/dx(i,j)
         enddo

      enddo

      return

    end subroutine getderiv

!#######################################################################

    subroutine constrain_deriv

      implicit none

      integer                   :: i,k
      real*8, dimension(maxdat) :: smin,smax
      real*8                    :: ftmp

      ! Loop over data sets
      do k=1,nfiles
         
         ! Calculate the smin and smax values for the current set
         do i=2,ndat(k)-1
            smin(i)=min(s(k,i-1),s(k,i+1))
            smax(i)=max(s(k,i-1),s(k,i+1))
         enddo

         smin(1)=min(s(k,1),s(k,2))
         smin(ndat)=min(s(k,ndat-1),s(k,ndat))

         smax(1)=max(s(k,1),s(k,2))
         smax(ndat)=max(s(k,ndat-1),s(k,ndat))

         ! Constrain the derivatives for the current set according
         ! to the Boor-Swartz monotonicity limit
         do i=2,ndat(k)-1
            if (smin(i).gt.0) then
               ftmp=max(0.0d0,deriv(k,i))
               deriv(k,i)=min(ftmp,3.0d0*smin(i))
            else if (smax(i).lt.0) then
               ftmp=min(0.0d0,deriv(k,i))
               deriv(k,i)=max(ftmp,3.0d0*smax(i))
            else if (s(k,i-1)*s(k,i).le.0) then
               deriv(k,i)=0.0d0
            endif
         enddo
         
      enddo

      return

    end subroutine constrain_deriv

!#######################################################################

    subroutine getcoeff
      
      implicit none
      
      integer :: k,i

      ! Loop over data sets
      do k=1,nfiles

         ! Loop over the intervals for the current set
         do i=1,nintvl(k)

            ! C1
            coeff(k,1,i)=dat(k,i)
            
            ! C2
            coeff(k,2,i)=deriv(k,i)
            
            ! C3
            coeff(k,3,i)=3.0d0*s(k,i)-deriv(k,i+1)-2.0d0*deriv(k,i)
            coeff(k,3,i)=coeff(k,3,i)/dx(k,i)
            
            ! C4
            coeff(k,4,i)=-(2.0d0*s(k,i)-deriv(k,i+1)-deriv(k,i))
            coeff(k,4,i)=coeff(k,4,i)/dx(k,i)**2
            
         enddo

      enddo

      return

    end subroutine getcoeff

!#######################################################################

    subroutine calcp

      use mcspmod

      implicit none

      integer              :: k,i,iint,unit
      real*8               :: diff,p,xcurr,lower,ftmp
      real*8, dimension(2) :: bound

!      npoints=1001

!-----------------------------------------------------------------------
! Determine the bounds of the data sets
!-----------------------------------------------------------------------
      bound(1)=9999.0d0
      bound(2)=-9999.0d0
      do k=1,nfiles
         do i=1,ndat(k)            
            if (x(k,i).lt.bound(1)) bound(1)=x(k,i)
            if (x(k,i).gt.bound(2)) bound(2)=x(k,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! Set the step length and initial energy
!-----------------------------------------------------------------------
      if (bound(1).lt.erange(1)) bound(1)=erange(1)
      if (bound(2).gt.erange(2)) bound(2)=erange(2)

      diff=(bound(2)-bound(1))/(real(npoints-1))

!-----------------------------------------------------------------------
! Calculate the interpolant at the grid points
!-----------------------------------------------------------------------
      ! Loop over points
      do i=1,npoints
         
         ! Set the current value of x
         xcurr=bound(1)+(i-1)*diff
         
         p=0.0d0

         ! Loop over data sets
         do k=1,nfiles
            
            if (xcurr.lt.x(k,1)) cycle
            if (xcurr.gt.x(k,ndat(k))) cycle

            ! Determine which interval we are in
            call getiint(iint,xcurr,k)

            ! Compute the value of the current interpolant at the 
            ! current value of x
            p=p+coeff(k,1,iint) &
                 + coeff(k,2,iint)*(xcurr-x(k,iint)) &
                 + coeff(k,3,iint)*(xcurr-x(k,iint))**2 &
                 + coeff(k,4,iint)*(xcurr-x(k,iint))**3

         enddo

         write(6,'(2(F10.4,2x))') xcurr,p

      enddo

      return

    end subroutine calcp

!#######################################################################

    subroutine getiint(iint,xcurr,k)

      implicit none

      integer :: iint,k,i
      real*8  :: xcurr

      do i=1,nintvl(k)
         if (xcurr.ge.x(k,i).and.xcurr.lt.x(k,i+1)) then
            iint=i
            exit
         endif
      enddo

      return

    end subroutine getiint

!#######################################################################

  end program mcspline
