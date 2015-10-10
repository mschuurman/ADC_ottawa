program adc

  use constants
  use parameters
  use read_param
  use guessvecs
  use adc1specmod
  use adc2specmod
  use adc2enermod
  use adc2dysonmod
  use rdinput
  use orbindx
  use defaults
  use iomod
  use channels
  use vpqrsmod
  use rungamess

  implicit none
  
  integer, dimension(2) :: shp
  real(d)               :: time
  character(len=72)     :: gam_chkpt,gam_log

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
  call load_gamess(gam_chkpt,gam_log)

!-----------------------------------------------------------------------
! Read orbital symmetries
!-----------------------------------------------------------------------
!  call rdorbsym

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
  ! Approximate Dyson orbital calculation
  if (ldyson) then
     call adc2_dyson()
  else
  ! Spectrum calculation
     select case(method)
     case(1) ! ADC(1), full diagonalisation
        call adc1_spec()
     case(2:3) ! ADC(2)-s or ADC(2)-x, spectrum calculation
        call adc2_spec()
     case(-3:-2) ! ADC(2)-s or ADC(2)-x, energy calculation
        call adc2_ener()
     end select
  endif

!-----------------------------------------------------------------------    
! Clean up the scratch files and exit
!-----------------------------------------------------------------------    
  call system('rm -rf SCRATCH')

  call cpu_time(time)

  write(ilog,'(/,a,1x,F9.2,1x,a)') 'Final Time:',time," s"

  STOP
  
end program adc

