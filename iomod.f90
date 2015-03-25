  module iomod

  contains

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

  end module iomod
