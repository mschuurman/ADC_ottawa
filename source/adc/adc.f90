program adc

  use constants
  use parameters
  use read_param
  use guessvecs
  use adc1specmod
  use adc2specmod
  use adc2enermod
  use adc2dysonmod
  use adc2rixsmod
  use adc2tpamod
  use adc2automod
  use adc2fdstatesmod
  use adc2propmod
  use rdinput
  use orbindx
  use defaults
  use iomod
  use channels
  use vpqrsmod
  use rungamess
  use import_gamess
  use timingmod
  
  implicit none
  
  integer               :: i
  integer, dimension(2) :: shp
  real(d)               :: tw1,tw2,tc1,tc2
  character(len=72)     :: gam_chkpt,gam_log
  type(gam_structure)   :: gamess_info

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
  call times(tw1,tc1)

!-----------------------------------------------------------------------
! Setting up the multiplication table
!-----------------------------------------------------------------------
  shp(:)=(/ 8,8 /)
  MT=reshape(mtrow, shp)

!-----------------------------------------------------------------------
! Initialise parameter values: this must be done before read_input
! is called
!
! Note that common I/O channel units are also set here
!-----------------------------------------------------------------------
  call set_defaults

!-----------------------------------------------------------------------
! Open the input and log files and create the scratch directory
!-----------------------------------------------------------------------
  call open_files

!-----------------------------------------------------------------------
! Read the input file: new code based on the use of keywords
!-----------------------------------------------------------------------
  call read_input

  close(iin)

!-----------------------------------------------------------------------
! Set the pointer to the vpqrs function
!-----------------------------------------------------------------------
!  if (motype.eq.'disk') then
!     vpqrs => vpqrs_ext
!  else if (motype.eq.'incore') then
!     vpqrs => vpqrs_incore
!  endif

!-----------------------------------------------------------------------
! GAMESS interface
!-----------------------------------------------------------------------
  if (lrungamess) call rungamess_main

  gam_chkpt  = 'gamess.dat'
  gam_log = 'gamess.log'

  call rdgeom(gam_log)
  call load_gamess(gam_chkpt,gam_log,gamess_info)

!-----------------------------------------------------------------------  
! Rearranging orbitals such that occ. orbs preceed unoccupied orbs.
!-----------------------------------------------------------------------  
  call rearrange_occ()

!-----------------------------------------------------------------------
! Determine various orbital indices
!-----------------------------------------------------------------------
  call get_hcentre

  call coreindx

  if (lfakeip) call contindx

!-----------------------------------------------------------------------  
! Set dipole array
!-----------------------------------------------------------------------  
  if (tranmom2 .eq. 'x') then
     dpl(:,:)=x_dipole(:,:)
  elseif (tranmom2 .eq. 'y') then
     dpl(:,:)=y_dipole(:,:)
  elseif (tranmom2 .eq. 'z') then
     dpl(:,:)=z_dipole(:,:)
  endif

!-----------------------------------------------------------------------
! Perform the ADC calculation
!-----------------------------------------------------------------------
  if (ldyson) then
     
     ! Approximate ADC(2) Dyson orbital calculation
     call adc2_dyson(gamess_info)

  else 
     
     select case(method)
        
     case(-3:-2) ! ADC(2)-s or ADC(2)-x, energy calculation

        call adc2_ener()

     case(1) ! ADC(1) OPA spectrum, full diagonalisation
           
        call adc1_spec()
        
     case (2:3) ! ADC(2)-s or ADC(2)-x, spectrum calculation

        if (lrixs) then

           ! ADC(2) RIXS spectrum
           call adc2_rixs(gamess_info)

        else if (ltpa) then

           ! ADC(2) TPA spectrum
           call adc2_tpa(gamess_info)

        else if (lautospec) then

           ! ADC(2) autocorrelation function calculation
           call adc2_autospec(gamess_info)

        else if (lfdstates) then

           ! ADC(2) filter diagonalisation state calculation
           call adc2_fdstates(gamess_info)

        else if (lpropagation) then

           ! TD-ADC(2) wavepacket propagation
           call adc2_propagate(gamess_info)
           
        else

           ! ADC(2) OPA spectrum
           call adc2_spec(gamess_info)

        endif

     end select
     
  endif

!-----------------------------------------------------------------------    
! Clean up the scratch files and exit
!-----------------------------------------------------------------------    
  call system('rm -rf SCRATCH')

!-----------------------------------------------------------------------    
! Output timings and stop
!-----------------------------------------------------------------------    
  call times(tw2,tc2)
  write(ilog,'(/,70a)') ('*',i=1,70)
  write(ilog,'(a,1x,F9.2,1x,a)') 'Final wall time:',tw2-tw1," s"
  write(ilog,'(a,1x,F9.2,1x,a)') 'Final cpu time: ',tc2-tc1," s"
  write(ilog,'(70a)') ('*',i=1,70)

end program adc

