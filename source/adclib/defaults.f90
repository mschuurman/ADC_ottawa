  module defaults

  contains

!#######################################################################

    subroutine set_defaults

      use constants
      use parameters
      use channels

      implicit none

!-----------------------------------------------------------------------
! I/O, verbosity and debugging
!-----------------------------------------------------------------------
      debug=.false.

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
      ldipole=.false.

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
      davtol=1e-7_d
      ladc1guess=.false.
      lsubdiag=.false.
      davname='SCRATCH/davstates'
      precon=2
      maxsubdim=-1
      ldfl=.true.
      solver=1
      guessdim=800
      
      ! Final space
      davstates_f=0
      maxiter_f=0
      dmain_f=0
      davtol_f=1e-7_d
      ladc1guess_f=.false.
      lsubdiag_f=.false.
      davname_f='SCRATCH/davstates_final'
      precon_f=2
      maxsubdim_f=-1
      ldfl_f=.true.
      solver_f=1
      guessdim_f=800

      ! Common
      ndavcalls=0
      eigentype=1

!-----------------------------------------------------------------------
! Wavepacket propagation parameters common to multiple calculation
! types
!-----------------------------------------------------------------------
      tfinal=0.0d0
      tout=0.0d0

      ! Initial space
      kdim=10
      stepsize=10.0d0
      siltol=1e-5_d
      
      ! Final space
      kdim_f=10
      stepsize_f=10.0d0
      siltol_f=1e-5_d
      
!-----------------------------------------------------------------------
! Relaxation parameters
!-----------------------------------------------------------------------
      ! Initial space
      lnoise=.false.
      rlxortho=2
      rlxtype=1
      
      ! Final space
      lnoise_f=.false.
      rlxortho_f=2
      rlxtype_f=1
      
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
! RIXS parameters
!-----------------------------------------------------------------------
      lrixs=.false.

!-----------------------------------------------------------------------
! TPA parameters
!-----------------------------------------------------------------------
      ltpa=.false.

!-----------------------------------------------------------------------
! Autospec parameters
!-----------------------------------------------------------------------
      lautospec=.false.
      autotol=1e-5_d
      autoord=0

!-----------------------------------------------------------------------
! Filter diagonalisation state calculation parameters
!-----------------------------------------------------------------------
      fdiagdat=''
      fdiagsel=''
      lfdstates=.false.
      
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
      ! Main Dyson orbital calculation parameters
      ldyson=.false.
      dysirrep=0
      dyslim=9999d0
      dysdiag=0
      dysout=0

      ! ezdyson input parameters
      lmax=4
      zcore=1.0d0
      nelen=10
      eleni=0.1d0
      elenf=10.0d0
      ngrdpnts=201
      grdi=-10.0d0
      grdf=10.0d0

!-----------------------------------------------------------------------
! Target state matching
!-----------------------------------------------------------------------
      ltarg=.false.
      detfile=''
      mofile=''
      detthrsh=-1.0d0
      ovrthrsh=-1.0d0

      return

    end subroutine set_defaults

!#######################################################################

  end module defaults
