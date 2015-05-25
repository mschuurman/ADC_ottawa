program main
  
  use constants
  use parameters
  use read_param
  use guessvecs
  use adc1mod
  use adc2mod
  use adc2extmod
  use rdinput
  use orbindx
  use defaults
  use iomod
  use channels

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
! GAMESS interface
!-----------------------------------------------------------------------  
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
     tranflag='y'
  elseif (tranmom2 .eq. 'y') then
     dpl(:,:)=y_dipole(:,:)
     tranflag='y'
  elseif (tranmom2 .eq. 'z') then
     dpl(:,:)=z_dipole(:,:)
     tranflag='y'
  end if

!-----------------------------------------------------------------------  
! Perform the ADC calculation
!-----------------------------------------------------------------------    
  select case(method)

  case(1) ! ADC(1), full diagonalisation
     call master_adc1_prop()

  case(2) ! ADC(2)-s, Lanczos pseudo spectrum
     call master_adc2_prop()

  case(3) ! ADC(2)-x, Lanczos pseudo spectrum
     call master_adc2ext_prop()

  case(-2) ! ADC(2)-s, standard calculation
     call master_adc2_ener()

  case(-3) ! ADC(2)-x, standard calculation
     call master_adc2ext_ener()

  end select

!-----------------------------------------------------------------------    
! Clean up the scratch files and exit
!-----------------------------------------------------------------------    
  call system('rm -rf SCRATCH')

  call cpu_time(time)

  write(ilog,'(/,a,1x,F9.2,1x,a)') 'Final Time:',time," s"

  STOP
  
end program main
  
  
