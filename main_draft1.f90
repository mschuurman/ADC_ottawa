program main
  
  use constants
  use parameters
  use read_param
  
  implicit none
  
  
  integer :: i,j
  real(d) :: time
  integer :: ipr,jpr

  integer, dimension(2) :: shp
  character(len=72)     :: gam_chkpt,gam_log

!!$ Setting up the multiplication table

  shp(:)=(/ 8,8 /)
  MT=reshape(mtrow, shp)

!!$ Reading basic information from Molcas generated files

!!!!!!! MOLCAS INTERFACE !!!!!!!!!!!!!!!!
!  call read_molcas0()

!  allocate(e(nBas),occNum(nBas),orbSym(nBas),roccnum(nBas))
!  allocate(x_dipole(nBas,nBas),y_dipole(nBas,nBas),z_dipole(nBas,nBas),dpl(nBas,nBas))

!  call read_user()
!  call read_molcas1(info)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! GAMESS INTERFACE !!!!!!!!!!!!!!!!
  
   gam_chkpt  = 'gamess.dat'
   gam_log = 'gamess.log'
   call read_user()

   call read_gamess(gam_chkpt,gam_log)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$ Reading user's data

!  call read_user()
  call check_user()

!!$ Rearranging orbitals such that occ. orbs preceed unoccupied orbs.

  call rearrange_occ()

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

!!!  write(6,*) "basic dipole elements"
!!!   do i=1,nBas
!!!      ipr=roccnum(i)
!!!   do j=1,i
!!!      jpr=roccnum(j)
!!!    write(8,*) i,j, dpl(ipr,jpr)
!!!   end do
!!!   end do

 
  if (method .eq. 0) then
     write(6,*) "Sorry not included!"

  elseif (method .eq. 1) then
     write(6,*) "Activating ADC1"
     call master_adc1_prop()

  elseif (method .eq. 2) then
     write(6,*) "Activating ADC2 "
     call master_adc2_prop()

  elseif (method .eq. 3) then
     write(6,*) "Activating ADC2 EXT"
     call master_adc2ext_prop()

  end if

  call cpu_time(time)
  write(6,*) 'Final Time=',time," s"
  
end program main
  
  
  
