  module iomod

    save

    character(len=120) :: errmsg

  contains

!#######################################################################

    subroutine open_files

      use constants
      use channels

      implicit none

      integer :: k
      logical :: lexists

!-----------------------------------------------------------------------
! Determine the input and log file names
!-----------------------------------------------------------------------
      ain=''

      call getarg(1,ain)

      if (ain.eq.'') then
         write(6,'(/,a,/)') 'The input file name has not been given'
         STOP
      endif

      k=index(ain,'.inp')

      if (k.eq.0) then
         alog=trim(ain)//'.log'
         ain=trim(ain)//'.inp'
      else
         alog=ain(1:k-1)//'.log'
      endif

!-----------------------------------------------------------------------
! Open the input and log files
!-----------------------------------------------------------------------
      open(iin,file=ain,form='formatted',status='old',err=999)

      open(ilog,file=alog,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Create the scratch directory
!-----------------------------------------------------------------------
      inquire(file='SCRATCH/.',exist=lexists)
      if (lexists) call system('rm -rf SCRATCH')
      call system('mkdir SCRATCH')

      return

999   continue
      write(6,'(/,3(1x,a))') 'The file',trim(ain),'does not exist'
      STOP

    end subroutine open_files

!#######################################################################

    subroutine freeunit(unit)

      implicit none

      integer         :: unit,i
      logical(kind=4) :: lopen

      ! N.B. Save the first 20 io units for standard files

      do i=20,1000
         inquire(unit=i,opened=lopen)
         if (.not.lopen) then
            unit=i
            exit
         endif
      enddo

      return

    end subroutine freeunit

!#######################################################################
!
! error_control: writes the passed string to the screen and, if open,
!                the log file, then terminates the program
!
!#######################################################################

    subroutine error_control

      use channels, only: ilog

      implicit none
      
      logical :: lopen

      ! Write error message to the screen
      write(6,'(/,2x,a,/)') trim(errmsg)

      ! If a log file is open, write the error message to the log file
      inquire(unit=ilog,opened=lopen)
      if (lopen) write(ilog,'(/,2x,a,/)') trim(errmsg)

      ! Terminate the program
      STOP

    end subroutine error_control

!#######################################################################

end module iomod
