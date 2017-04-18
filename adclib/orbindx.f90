  module orbindx

  contains

!#######################################################################

    subroutine get_hcentre

      use constants
      use parameters

      implicit none

      integer :: i
      
!-----------------------------------------------------------------------
! Determine the indices of the orbitals that can carry a hole, i.e.,
! the occupied orbitals
!
! N.B. hcenter(0) is the number of 'hole-carrying' orbitals
!-----------------------------------------------------------------------
      hcentre(0)=nocc
      do i=1,nocc
         hcentre(i)=i
      enddo
      
      return

    end subroutine get_hcentre

!#######################################################################
    
    subroutine coreindx

      use constants
      use parameters

      implicit none

      integer            :: i
      real(d), parameter :: tol=-5.0d0

!-----------------------------------------------------------------------
! Determine the number of core orbitals and their indices
!-----------------------------------------------------------------------
      ncore=0
      do i=1,nbas
         if (e(i).lt.tol) then
            ncore=ncore+1
            icore(ncore)=i
         endif
      enddo



      return

    end subroutine coreindx

!#######################################################################

    subroutine contindx

      use constants
      use parameters
      use channels

      implicit none

      integer            :: i,ncont
      real(d), parameter :: tol=1e-15_d

!-----------------------------------------------------------------------
! Determine the index of the 'continuum' orbital
!-----------------------------------------------------------------------
      ncont=0
      do i=1,nbas
         if (abs(e(i)).le.tol) then
            ncont=ncont+1
            ifakeorb=i
         endif
      enddo

!-----------------------------------------------------------------------
! Exit if the number of 'continuum' orbitals is not one
!-----------------------------------------------------------------------
      if (ncont.eq.0) then
         write(ilog,'(/,2x,a,/)') 'No continuum orbital has been found'
         STOP
      endif

      if (ncont.gt.1) then
         write(ilog,'(/,2x,a,/)') 'More than one continuum orbital has &
              been found'
         STOP
      endif

      return

    end subroutine contindx

!#######################################################################

  end module orbindx
