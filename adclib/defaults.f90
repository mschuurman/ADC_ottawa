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
      method_f=0
      denord=2
      nirrep=0
      nirrep2=0
      statenumber=-1
      tranmom=''
      tranmom2=''
      norder=2
      motype='incore'
      dlim=0.0d0
      ltdm_gs2i=.true.
      lifrzcore=.false.
      lffrzcore=.false. 
      ldiagfinal=.false.
      hinit=1
      maxmem=250.0d0

!-----------------------------------------------------------------------
! CVS-ADC parameters
!-----------------------------------------------------------------------      
      lcvs=.false.
      lcvsfinal=.false.
      icore=0
      iexpfrz=0

!-----------------------------------------------------------------------
! Ionisation potential calculation parameters
!-----------------------------------------------------------------------      
      lfakeip=.false.
      ifakeorb=0

!-----------------------------------------------------------------------
! Diagonalisation parameters
!-----------------------------------------------------------------------
      ! Initial space
      davstates=0
      maxiter=0
      dmain=0
      davtol=1d-7
      ladc1guess=.false.
      davname='SCRATCH/davstates'
      precon=1
      maxsubdim=-1
      
      ! Final space
      davstates_f=0
      maxiter_f=0
      dmain_f=0
      davtol_f=1d-7
      ladc1guess_f=.false.
      davname_f='SCRATCH/davstates_final'
      precon_f=1
      maxsubdim_f=-1
      
      ! Common
      ndavcalls=0
      eigentype=1
      solver=1

!-----------------------------------------------------------------------
! Lanczos parameters
!-----------------------------------------------------------------------      
      lmain=0
      ncycles=0
      lancguess=1
      lancname='SCRATCH/lancstates'      
      ldynblock=.false.
      tdtol=0.0005d0
      orthotype=0

!-----------------------------------------------------------------------
! I/O channels
!-----------------------------------------------------------------------
      iin=1
      ilog=2

!-----------------------------------------------------------------------
! SCF parameters
!-----------------------------------------------------------------------
      scfiter=10

!-----------------------------------------------------------------------
! GAMESS parameters
!-----------------------------------------------------------------------
      lrungamess=.false.
      basname=''
      natm=0
      difftype=0
      ndiff=0
      pntgroup=''

!-----------------------------------------------------------------------
! Dyson orbital calculation parameters
!-----------------------------------------------------------------------
      ldyson=.false.
      dysirrep=0
      dyslim=9999d0
      ldysfulldiag=.false.
      
      return

    end subroutine set_defaults

!#######################################################################

  end module defaults
