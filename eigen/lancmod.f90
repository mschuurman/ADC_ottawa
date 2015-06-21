  module lancmod

  contains

!#######################################################################
    
    subroutine master_lancdiag(ndim,noff,flag)

      use parameters
      use band_lanczos
      use block_lanczos
      
      implicit none

      integer, intent(in)     :: ndim
      integer*8, intent(in)   :: noff
      character(1),intent(in) :: flag

      if (lanctype.eq.1) then
         call lancdiag_band(ndim,noff,flag)
      else if (lanctype.eq.2) then
         call lancdiag_block(ndim,noff,flag)
      endif
      
      return
      
    end subroutine master_lancdiag

!#######################################################################
    
  end module lancmod
