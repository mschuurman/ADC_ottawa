program main
  
  use constants
  use parameters
  use read_param
  use guessvecs
  
  implicit none
  
  
  integer               :: i,j
  integer               :: ipr,jpr
  integer, dimension(2) :: shp
  real(d)               :: time
  character(len=72)     :: gam_chkpt,gam_log

!-----------------------------------------------------------------------
! Setting up the multiplication table
!-----------------------------------------------------------------------
  shp(:)=(/ 8,8 /)
  MT=reshape(mtrow, shp)

!-----------------------------------------------------------------------
! GAMESS interface
!-----------------------------------------------------------------------  
  gam_chkpt  = 'gamess.dat'
  gam_log = 'gamess.log'
  call read_user()
  
  call read_gamess(gam_chkpt,gam_log)

!-----------------------------------------------------------------------  
! Reading user's data
!-----------------------------------------------------------------------  
!  call read_user()
  call check_user()

!-----------------------------------------------------------------------  
! Rearranging orbitals such that occ. orbs preceed unoccupied orbs.
!-----------------------------------------------------------------------  
  call rearrange_occ()

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
! Calculate guess initial space vectors from an ADC(1) calculation if 
! requested
!-----------------------------------------------------------------------  
  if (method.ne.1.and.ladc1guess) call adc1_guessvecs

!-----------------------------------------------------------------------  
! Perform the ADC Stieltjes imaging calculation
!-----------------------------------------------------------------------  
  if (method .eq. 0) then
     write(6,*) "Sorry not included!"
  
  ! Calculation of spectral moments
  elseif (method .eq. 1) then
     write(6,*) "Activating ADC1"
     call master_adc1_prop()
     
  elseif (method .eq. 2) then
     write(6,*) "Activating ADC2 "
     call master_adc2_prop()
     
  elseif (method .eq. 3) then
     write(6,*) "Activating ADC2 EXT"
     call master_adc2ext_prop()
  
  ! Calculations of VEEs only
  elseif (method.eq.-2) then
     call master_adc2_ener()

  end if
  
  call cpu_time(time)
  write(6,*) 'Final Time=',time," s"
  
end program main
  
  
