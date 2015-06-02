  module defaults

  contains

!#######################################################################

    subroutine set_defaults

      use constants
      use parameters
      use channels

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
      lifrzcore=.false.
      lffrzcore=.false. 
      ldavfinal=.false.

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
      ! Initial space
      davstates=0
      maxiter=0
      dmain=0
      davtol=1d-7
      ladc1guess=.false.
      davname='SCRATCH/davstates'

      ! Final space
      davstates_f=0
      maxiter_f=0
      dmain_f=0
      davtol_f=1d-7
      ladc1guess_f=.false.
      davname_f='SCRATCH/davstates_final'
      
      ! Common
      ndavcalls=0

!-----------------------------------------------------------------------
! Lanczos parameters
!-----------------------------------------------------------------------      
      lmain=0
      ncycles=0
      lancguess=1
      lancname='SCRATCH/lancstates'

!-----------------------------------------------------------------------
! I/O channels
!-----------------------------------------------------------------------
      iin=1
      ilog=2

!-----------------------------------------------------------------------
! SCF parameters
!-----------------------------------------------------------------------
      scfiter=10

      return

    end subroutine set_defaults

!#######################################################################

  end module defaults
