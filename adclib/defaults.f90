  module defaults

  contains

!#######################################################################

    subroutine set_defaults

      use constants
      use parameters

      implicit none

!-----------------------------------------------------------------------
! General ADC parameters
!-----------------------------------------------------------------------
      method=0
      nirrep=0
      nirrep2=0
      statenumber=-1
      tranmom=''
      tranmom2=''
      norder=2
      motype='incore'
      dlim=0.0d0
      ltdm_gs2i=.true.
      dmatmem=250.0d0

!-----------------------------------------------------------------------
! CVS-ADC parameters
!-----------------------------------------------------------------------      
      lcvs=.false.
      lcvsfinal=.false.
      icore=0

!-----------------------------------------------------------------------
! Ionisation potential calculation parameters
!-----------------------------------------------------------------------      
      lfakeip=.false.
      ifakeorb=0

!-----------------------------------------------------------------------
! Davidson parameters
!-----------------------------------------------------------------------
      davstates=0
      maxiter=0
      dmain=0
      davtol=1d-7
      ladc1guess=.false.
      davname='SCRATCH/davstates'

!-----------------------------------------------------------------------
! Lanczos parameters
!-----------------------------------------------------------------------      
      lmain=0
      ncycles=0
      lancname='SCRATCH/lancstates'

!-----------------------------------------------------------------------
! I/O channels
!-----------------------------------------------------------------------
      iin=1
      ilog=2

      return

    end subroutine set_defaults

!#######################################################################

  end module defaults
