 ! 
 ! A little hack-y, but temporary for debugging purposes. Simply return the
 ! value of the integral buffer using the native spin-indices.
 !
 function vpqrs(r,s,u,v)
   use parameters
   implicit none
   integer,intent(in) :: r,s,u,v
   integer            :: r2
   real(d)            :: vpqrs

   if(moType == 'disk') then
     r2 = u
     if(moIntegrals%mo_l /= v .and. moIntegrals%mo_l /= u)then
       call fetch_moint2e(moIntegrals,v)
     else
       if(moIntegrals%mo_l == u)r2 = v
     endif
     vpqrs = real(moIntegrals%buffer_real(r,s,r2,1),kind=d)
   else
     vpqrs = real(moIntegrals%buffer_real(r,s,u,v),kind=d)
   endif

   return
 end function vpqrs



 subroutine errmsg(message)
  implicit none
  character*144         :: message

  write(6,*)message
  stop '-- Program excited abnormally -- '
  return 
 end subroutine

 !
 ! read a gamess output file
 !
 subroutine read_gamess_output(nbasis,nelec,ncen,nirr,orbsym,symlab,&
      ehf,earr,occnum)
  use constants
  integer,intent(in)            :: nbasis
  integer*4,intent(out)         :: nelec,nirr,ncen
  integer*4,intent(out)         :: orbsym(nbasis)
  real(d),intent(out)           :: ehf      ! The Hf energy
  real(d),intent(out)           :: earr(nbasis)  ! The orbital energies
  character*2,intent(out)       :: symlab(1024)
  real(d), dimension(nbasis)    :: occnum

  logical                       :: gexist
  integer*4                     :: cnt,i
  integer                       :: off1,off2
  integer                       :: ios
  integer                       :: orbfnd
  integer                       :: gamess=11
  real(d)                       :: two = 2.
  character(len=144)            :: line,scr

  inquire(file='gamess.log',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.log')
  else
   scr = 'gamess.log does not exist.'
   call errmsg(scr)
  endif

  nirr   = -1
  ncen   = -1
  occnum = 0.0d0
  ehf    = 0.0d0
  earr   = 0.0d0
  orbfnd = 0

  scan_lines: do while (orbfnd==0)
   call read_gamess_line(gamess,line)

   ! read number of electrons
   if(index(line,'NUMBER OF ELECTRONS').ne.0)then
     read(line,'(a47,i5)')scr,nelec
   endif

   ! read number of atoms
   if(index(line,'TOTAL NUMBER OF ATOMS').ne.0)then
     read(line,'(a47,i5)')scr,ncen
   endif

   ! read in orbital symmetry information
   if(index(line,'DIMENSIONS OF THE SYMMETRY SUBSPACES ARE').ne.0)then
     nirr=0
     scan_olabels: do
      read(gamess,'(a144)')line
      if(index(line,'=').eq.0)exit
      off1=1
      do
       if(index(line,'=').eq.0)exit
       nirr = nirr + 1
       off1 = index(line,'=')+1
       symlab(nirr) = trim(adjustl(line(off1-5:off1-2)))
       line = trim(adjustl(line(off1:len_trim(line))))
      enddo
     enddo scan_olabels
   endif

   ! read in HF energy
   if(index(line,'FINAL RHF ENERGY IS').ne.0 .and. ehf==0.) then
    read(line,'(a20,f20.10,a20)')scr,ehf,scr
   endif

   ! read in orbital occupations and energies
   if(index(line,'EIGENVECTORS').ne.0 .and. orbfnd==0) then
    call read_gamess_line(gamess,line) ! ---- format line
    call read_gamess_line(gamess,line) ! blank line
    off1=1
    off2=0
    scan_orbs: do while(index(line,'END OF RHF')==0)
      call read_gamess_line(gamess,line) ! orbital indices
      read(gamess,'(a15,5(f11.4))')scr,earr(off1:min(nbasis,off1+4))
      off1 = off1 + 5
      call read_gamess_line(gamess,line) ! irreps
      scan_irreps: do while(len_trim(line).gt.0)
       line = adjustl(line)
       cnt = 1
       do while(line(1:2) .ne. symlab(cnt))
        cnt = cnt + 1
       enddo
       off2 = off2 + 1
       orbsym(off2) = cnt
       line = line(3:len_trim(line))
      enddo scan_irreps

      do i = 1,nbasis
       call read_gamess_line(gamess,line) ! mos
      enddo
      call read_gamess_line(gamess,line) ! blank line/termination string
    enddo scan_orbs
    orbfnd=1
   endif

   ! BODGE
!   earr(20)=1d-10
   ! BODGE

  enddo scan_lines
  close(gamess)

 end subroutine read_gamess_output

 !
 ! subroutine to read a line of gamess output
 !
 subroutine read_gamess_line(unitn,line)
  implicit none
  integer,intent(in)                    :: unitn
  character(len=144),intent(out)        :: line
  integer                               :: ios

  read(unitn,'(a144)',iostat=ios)line
  if(ios < 0) then
   write(6,*)'ERROR in reading gamess log file: '//trim(line)
   stop
  endif
  if(ios > 0) then
   write(6,*)'EOF reached in gamess log file.'
   stop
  endif

  return
 end subroutine read_gamess_line

 ! 
 ! load MO integrals
 ! Jan. 30, 2015 -- assume we are working with RHF reference in a basis of
 !                  spatial orbitals.  This really only effects vpqrs...
 !
 subroutine load_mo_integrals(gam)
  use accuracy
  use math
  use parameters
  use import_gamess
  use integral_tools
  use integrals_mo2e
  use printing
  implicit none
  type(gam_structure)               :: gam ! gamess info (orbitals, geom, etc.) 
  type(int2e_cache)                 :: int2e    ! Currently active integrals context
  integer, parameter                :: iu_2e_ao   = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(hik)                      :: iosize_2e  = 220000000   ! Integral I/O

  real(xrk), dimension(:,:),allocatable                       :: mos_real,hao,hmo,fmo
  complex(xrk),dimension(:,:),allocatable                     :: mos_cmplx,hao_cmplx
  character(len=clen)                                         :: ao_mode,int_type
  integer(ik)                                                 :: a,i,j,nmo,nao,nao_spin,nocc_spin
  real(d)                                                     :: vpqrs
  real(xrk)                                                   :: xyz(3), d_cf, q, ov, eps, refval, e1,e2
  real(xrk)                                                   :: nuc_repulsion

  ! make a call to MathLogFactorial before we get into the parallel parts, since
  ! it isn't openmp safe...
  q = MathFactorial(nelec)
  q = MathLogFactorial(nelec)
  q = MathDoubleFactorial(nelec)

  eps      = 1.d-5
  nao      = gam%nbasis
  nao_spin = 2*nao
  nmo      = gam%nvectors

  ! allocate arrays
  allocate(mos_cmplx(nao_spin,nmo),mos_real(nao_spin,nmo),hao(nao,nao),hao_cmplx(nao_spin,nao_spin),hmo(nmo,nmo),fmo(nmo,nmo))
  allocate(x_dipole(nmo,nmo),y_dipole(nmo,nmo),z_dipole(nmo,nmo),dpl(nmo,nmo))

  ! assume for the time being ao integrals stored in core
  ao_mode = 'incore'
  call prepare_2e(int2e,gam,ao_mode,iu_2e_ao,iosize_2e,ints_math='real')

  ! form 1e Hamiltonian (in AO basis)
  call core_hamiltonian(gam,hao)
  hao_cmplx                                = 0 
  hao_cmplx(1:nao,1:nao)                   = hao(1:nao,1:nao)
  hao_cmplx(nao+1:nao_spin,nao+1:nao_spin) = hao(1:nao,1:nao)

  call scf_loop(int2e,gam,hao_cmplx,mos_cmplx)
  mos_real = realpart(mos_cmplx)

  ! load up variables depending on RHF/UHF case 
  if(nmo == nao) then         ! RHF case
      nocc_spin = nelec / 2
      d_cf = 2._xrk
   else if(nmo == nao_spin) then ! UHF case
      nocc_spin = nelec
      d_cf = 1._xrk
   else
      stop 'load_mo_integrals - confusing number of mos'
   endif 

  ! actual call to general integrals in MO basis 
  call transform_moint2e_real(int2e,moType,mos_cmplx,mos_cmplx,mos_cmplx,mos_cmplx,moIntegrals,io_unit=99,l_block=10)

  ! form 1e Hamiltonian in MO basis
  hmo = realpart(matmul(matmul(transpose(mos_cmplx),hao_cmplx),mos_cmplx))

  ! form Fock matrix in spin-MO basis (i.e. diagonal)
  fmo = hmo
  sum_orb: do a = 1,nocc_spin
   sum_row: do i = 1,nmo
    sum_col: do j = 1,nmo
     fmo(i,j) = fmo(i,j) + d_cf*vpqrs(i,j,a,a) - vpqrs(i,a,j,a)
    enddo sum_col
   enddo sum_row
  enddo sum_orb

  if(debug) then
   write(6,"(/t5,'FOCK MATRIX IN MO BASIS')")
   call print_matrix(fmo,15,'f15.10')
  endif

  ! ensure Fock matrix is diagonal
  refval = 0._rk
  scan_row: do i = 1,nmo
   scan_col: do j = 1,nmo
    if(i /= j .and. abs(fmo(i,j)-refval)>eps)write(6,"('warning: <',i3,'|f|',i3,'> = ',f15.10,',expecting ',f15.10)")i,j,fmo(i,j),refval
   enddo scan_col
  enddo scan_row

  ! Compute HF energy and print, using spatial orbitals
  e1 = 0._xrk
  e2 = 0._xrk
  sum_i: do j = 1,nocc_spin
   e1 = e1 + d_cf*hmo(j,j)
   sum_j: do i = 1,nocc_spin
    e2 = e2 + d_cf*vpqrs(i,i,j,j) - vpqrs(i,j,i,j)
   enddo sum_j
  enddo sum_i
  
  ! in the future, we will fill orbitals using aufbau principle. However, I get
  ! the impression this code works primarily with spatial orbitals (i.e. under
  ! RHF). So, we trivially fill up the lowest nelec/2 orbitals, but this is
  ! easily changed in the future.
  do i = 1,nocc_spin
   occNum(i) = int(2,kind=4)
  enddo

  ! load up the dipole moment integrals
  int_type = 'AO DIPOLE X'
  call dipole_integrals(gam,int_type,nao_spin,nmo,x_dipole,mos_real)
  int_type = 'AO DIPOLE Y'
  call dipole_integrals(gam,int_type,nao_spin,nmo,y_dipole,mos_real)
  int_type = 'AO DIPOLE Z'
  call dipole_integrals(gam,int_type,nao_spin,nmo,z_dipole,mos_real)

  ! assuming dpl is the norm of the dipole moment vector
  do i = 1,nmo
   do j = 1,nmo
     dpl(i,j) = sqrt(x_dipole(i,j)**2 + y_dipole(i,j)**2 + z_dipole(i,j))
   enddo
  enddo

  ! print out summary (compare to GAMESS output)
  write(6,"(60('-'))")
  write(6,1000)
  write(6,1001)e1
  write(6,1002)e2
  write(6,1003)nuc_repulsion(gam)
  write(6,1004)e1+e2+nuc_repulsion(gam)
  write(6,"(60('-'))")
  call flush

  ! clear integral cache
!  call clear_2e(int2e)

  ! clean up all the allocated arrays
  deallocate(mos_real,mos_cmplx,hao,hmo,fmo)

  return
1000 format(' Computation of HF energy in mo basis:')
1001 format('  one electron energy:      ',f17.10)
1002 format('  two electron energy:      ',f17.10)
1003 format('  nuclear repulsion energy: ',f17.10)
1004 format('  total energy:             ',f17.10)
 end subroutine load_mo_integrals

 ! 
 ! 
 !
 subroutine phis_init
  implicit none

  return
 end subroutine

 ! 
 ! Determine the number of irreps, number of basis functions and number of atoms
 !
 subroutine phis_get_info(nirr,symlab,nbas,ncen)
  implicit none
  integer       :: nirr
  character(2)  :: symlab 
  integer       :: nbas
  integer       :: ncen
 
  return
 end subroutine

 ! 
 ! Get Hartree-Fock energy information: total energy, and energy of the orbitals
 !
 subroutine phis_get_epsi(ehf,earr,nbas)
  use constants
  implicit none
  integer*4,intent(in)   :: nbas     ! The number of basis functions
  real(d),intent(out)    :: ehf      ! The Hf energy
  real(d),intent(out)    :: earr(nbas)  ! The orbital energies

  return
 end subroutine

 !
 ! Get symmetries of each of the orbitals
 !
 subroutine phis_get_sym(orbsym,symlab,nbas)
  implicit none
  integer*4,intent(in)    :: nbas      ! The number of  basis functions
  integer*4,intent(inout) :: orbsym(nbas) ! Symmetry each orbital
  character*2,intent(in):: symlab(1024)

  return
 end subroutine

 !
 ! Get the occupation number for each orbital
 !
 subroutine phis_get_occ(occnum,nbas)
  use constants
  implicit none
  integer*4,intent(in)    :: nbas      ! The total number of basis functions
  real(d),intent(inout) :: occnum(nbas) ! The occupation number for each orbital 

  return
 end subroutine

 !
 ! Load up the 
 !
 subroutine phis_mc_dip(xdip,ydip,zdip,nbas)
  implicit none
  integer,intent(in)  :: nbas
  real,intent(inout)  :: xdip(nbas,nbas),ydip(nbas,nbas),zdip(nbas,nbas)

  return
 end subroutine

 !
 ! Load the MO integrals into memory
 !
 subroutine phis_load_vpqrs
  implicit none

  return
 end subroutine

