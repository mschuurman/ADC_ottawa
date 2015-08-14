  program knit

    implicit none

!-----------------------------------------------------------------------
! Read the command line arguments
!-----------------------------------------------------------------------
    call rdcommands

!-----------------------------------------------------------------------
! Read the bound and continuum oscillator strengths
!-----------------------------------------------------------------------
    call rdoscfiles

!-----------------------------------------------------------------------
! Form the total spectrum
!-----------------------------------------------------------------------
    call mktotspec

    STOP
    
  contains

!#######################################################################

    subroutine rdcommands

      use iomod
      use knitmod

      implicit none

      integer            :: iarg
      character(len=120) :: string,atmp

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      abound=''
      acont=''
      et=-999d0

!-----------------------------------------------------------------------
! Read the command line arguments
!-----------------------------------------------------------------------
      iarg=0

10    iarg=iarg+1
      call getarg(iarg,string)
      
      if (string.eq.'-b') then
         iarg=iarg+1
         call getarg(iarg,abound)

      else if (string.eq.'-c') then
         iarg=iarg+1
         call getarg(iarg,acont)

      else if (string.eq.'-et') then
         iarg=iarg+1
         call getarg(iarg,atmp)
         read(atmp,*) et

      else
         errmsg='Unknown argument: '//trim(string)
         call error_control
      endif

      if (iarg.lt.iargc()) goto 10

!-----------------------------------------------------------------------
! Check that all the required information has been given
!-----------------------------------------------------------------------
      if (abound.eq.'') then
         errmsg='The bound file name has not been given'
         call error_control
      endif

      if (acont.eq.'') then
         errmsg='The continuum file name has not been given'
         call error_control
      endif

      if (et.eq.-999d0) then
         errmsg='The ionisation threshold has not been given'
         call error_control
      endif

      return

    end subroutine rdcommands

!#######################################################################

    subroutine rdoscfiles

      use iomod
      use parsemod
      use knitmod

      implicit none
      
      integer :: unit,k

!-----------------------------------------------------------------------
! Bound part
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file=abound,form='formatted',status='old')

      ! (1) Determine the no. points and allocate arrays
      npntb=0
10    call rdinp(unit)
      if (.not.lend) then
         npntb=npntb+1
         goto 10
      endif
     
      allocate(eb(npntb))
      allocate(fb(npntb))

      ! (2) Read the bound-state energies and oscillator strengths
      rewind(unit)
      k=0
20    call rdinp(unit)
      if (.not.lend) then
         k=k+1
         read(keyword(1),*) eb(k)
         read(keyword(2),*) fb(k)
         goto 20
      endif

      close(unit)

!-----------------------------------------------------------------------
! Continuum part
!-----------------------------------------------------------------------
      open(unit,file=acont,form='formatted',status='old')

      ! (1) Determine the no. points and allocate arrays
      npntc=0
30    call rdinp(unit)
      if (.not.lend) then
         npntc=npntc+1
         goto 30
      endif
     
      allocate(ec(npntc))
      allocate(fc(npntc))

      ! (2) Read the bound-state energies and oscillator strengths
      rewind(unit)
      k=0
40    call rdinp(unit)
      if (.not.lend) then
         k=k+1
         read(keyword(1),*) ec(k)
         read(keyword(2),*) fc(k)
         goto 40
      endif

      close(unit)

      return

    end subroutine rdoscfiles

!#######################################################################

    subroutine mktotspec

      use constants
      use knitmod

      implicit none

      real(d) :: fac
      integer :: i,k,lubb,glbc

!-----------------------------------------------------------------------
! Determine where to truncate the bound and continuum parts
!-----------------------------------------------------------------------
      lubb=npntb
      do i=1,npntb
         if (eb(i).ge.et) then
            lubb=i
            exit
         endif
      enddo

      glbc=1
      do i=1,npntc
         if (ec(i).ge.et) then
            glbc=i-1
            exit
         endif
      enddo

!-----------------------------------------------------------------------
! Scale the continnum part s.t. it agrees with the bound part at the
! cutoff
!-----------------------------------------------------------------------
      fac=fc(glbc)/fb(lubb)
      fc=fc/fac

!-----------------------------------------------------------------------
! Output the combined spectrum
!-----------------------------------------------------------------------
      do i=1,lubb
         print*,eb(i),fb(i)
      enddo

      do i=glbc,npntc
         print*,ec(i),fc(i)
      enddo

      return

    end subroutine mktotspec

!#######################################################################

  end program knit
