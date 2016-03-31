  program numhess

    use hessmod
    use prepmod
    use calcmod

    implicit none

!-----------------------------------------------------------------------
! Read the command line arguments
!-----------------------------------------------------------------------
    call rdargs

!-----------------------------------------------------------------------
! Read the reference file
!-----------------------------------------------------------------------
      call rdreffile

!-----------------------------------------------------------------------
! Perform the requested job:
!
! ijob = 1 <-> ADC input preparation
!        2 <-> Hessian calculation
!-----------------------------------------------------------------------
    if (ijob.eq.1) then
       call hessprep
    else
       call hesscalc
    endif

    stop

  contains

!#######################################################################

    subroutine rdargs

      use iomod
      use channels
      use hessmod

      implicit none

      integer            :: k,narg
      character(len=120) :: string,string2

!-----------------------------------------------------------------------
! Read the command line arguments
!-----------------------------------------------------------------------
      infile=''
      listfile=''
      ijob=-1

      narg=iargc()

      k=0
5     k=k+1
      call getarg(k,string)

      if (string.eq.'-r') then
         k=k+1
         call getarg(k,infile)
      else if (string.eq.'-j') then         
         k=k+1
         call getarg(k,string2)
         if (string2.eq.'prep') then
            ijob=1
         else if (string2.eq.'calc') then
            ijob=2
         else
            errmsg='Unknown job type: '//trim(string2)
            call error_control
         endif
      else if (string.eq.'-l') then
         k=k+1
         call getarg(k,listfile)
      else
         errmsg='Unknown keyword: '//trim(string)
         call error_control
      endif

      if (k.lt.narg) goto 5

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
      if (ijob.eq.-1) then
         errmsg='The job type has not been given'
         call error_control
      endif

      if (infile.eq.'') then
         errmsg='The name of the reference file has not been given'
         call error_control
      endif

      if (ijob.eq.2.and.listfile.eq.'') then
         errmsg='The name of the list file has not been given'
         call error_control
      endif

      return

    end subroutine rdargs

!#######################################################################

    subroutine rdreffile

      use iomod
      use channels
      use hessmod

      implicit none

      integer            :: i,j,k
      character(len=120) :: string

!-----------------------------------------------------------------------
! Open the reference file
!-----------------------------------------------------------------------
      call freeunit(iin)
      open(iin,file=infile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the no. lines in the reference file and allocate arrays
!-----------------------------------------------------------------------
      nlines=0
10    read(iin,'(a)',end=15) string
      nlines=nlines+1
      goto 10

15    continue
      allocate(aline(nlines))

!-----------------------------------------------------------------------
! Read the reference file
!-----------------------------------------------------------------------
      rewind(iin)
      do i=1,nlines
         read(iin,'(a)') aline(i)
      enddo

!-----------------------------------------------------------------------
! Determine the start of the geometry section
!-----------------------------------------------------------------------
      geomline=-1
      do i=1,nlines
         if (index(aline(i),'geometry').ne.0) then
            geomline=i
            exit
         endif
      enddo

      if (geomline.eq.-1) then
         errmsg='No geometry section could be found in the reference &
              file'
         call error_control
      endif

!-----------------------------------------------------------------------
! Determine the no. atoms and allocate arrays
!-----------------------------------------------------------------------
      natm=0
      k=geomline+1
20    continue
      if (index(aline(k),'end-geometry').eq.0) then
         natm=natm+1
         k=k+1
         goto 20
      endif

      ncoo=3*natm

      allocate(xcoo(ncoo))
      allocate(xcoo0(ncoo))
      allocate(aatm(natm))

!-----------------------------------------------------------------------
! Read the reference Cartesian coordinates
!-----------------------------------------------------------------------
      k=0
      do i=geomline+1,geomline+natm
         k=k+1
         read(aline(i),*) aatm(k),(xcoo0(j),j=k*3-2,k*3)
      enddo

!-----------------------------------------------------------------------
! Close the reference file
!-----------------------------------------------------------------------
      close(iin)

      return

    end subroutine rdreffile

!#######################################################################

  end program numhess
