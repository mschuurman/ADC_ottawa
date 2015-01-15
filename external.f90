
 real function vpqrs(r,s,u,v)
   use parameters
   implicit none
   integer,intent(in) :: r,s,u,v

   vpqrs = real(moIntegrals%buffer_real(r,s,u,v),kind=d)      

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
 subroutine read_gamess_output(nbasis,ncen,nirr,orbsym,symlab,occnum,ehf,earr)
  use constants
  integer,intent(in)            :: nbasis
  integer*4,intent(out)         :: nirr,ncen
  integer*4,intent(out)         :: orbsym(nbasis)
  real(d),intent(out)           :: occnum(nbasis) ! The occupation number for each orbital 
  real(d),intent(out)           :: ehf      ! The Hf energy
  real(d),intent(out)           :: earr(nbasis)  ! The orbital energies
  character*2,intent(out)       :: symlab(1024)

  logical                       :: gexist
  integer*4                     :: cnt
  integer                       :: off1,off2
  integer                       :: nelec
  integer                       :: ios
  integer                       :: orbfnd
  integer                       :: gamess=11
  real(d)                       :: two = 2.
  character(len=144)            :: line,scr

  print *,'nbasis=',nbasis

  inquire(file='gamess.log',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.log')
  else
   scr = 'gamess.log does not exist.'
   call errmsg(scr)
  endif

  nirr   = -1
  ncen   = -1
  occnum = 0
  ehf    = 0.
  earr   = 0.
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
      read(gamess,'(a15,5(f11.4))')scr,earr(off1:off1+4)
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

  enddo scan_lines
  close(gamess)


  cnt = 1
  do while(nelec.gt.0)
   occnum(cnt) = min(1.*nelec,2.)
   nelec = nelec - occnum(cnt)
   cnt = cnt + 1
  enddo


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
 ! 
 subroutine load_mo_integrals(gam_info)
  use parameters
  use accuracy
  use import_gamess
  use integral_tools
  use integrals_mo2e
  use printing
  implicit none
  type(gam_structure)               :: gam_info ! gamess info (orbitals, geom, etc.) 
  type(int2e_cache)                 :: int2e    ! Currently active integrals context
  integer, parameter                :: iu_2e_ao   = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(8)                       :: iosize_2e  = 220000000   ! Integral I/O
  complex(8),dimension(2*gam_info%nbasis,gam_info%nbasis)  :: mos
  character(len=clen)               :: mo_mode,ao_mode
  integer                           :: iat,i,j,k,l
  real                              :: vpqrs
  real(xrk),dimension(gam_info%nbasis,gam_info%nbasis) :: hao,sao,hmo,tmp_xk
  real(xrk)                            :: xyz(3), q


  ao_mode = 'conventional'
  mo_mode = 'incore'
  mos = transpose(cmplx(gam_info%vectors))
  call prepare_2e(int2e,gam_info,ao_mode,iu_2e_ao,iosize_2e,ints_math='real')
  call transform_moint2e_real(int2e,mo_mode,mos,mos,mos,mos,moIntegrals) 
 
  call gamess_1e_integrals('AO OVERLAP',sao,bra=gam_info,ket=gam_info  )
  write (6,"(/t5,'AO OVERLAP INTEGRALS'/)")
  call gamess_print_1e_integrals(sao,bra=gam_info,ket=gam_info)

  call gamess_1e_integrals('AO KINETIC',hao,bra=gam_info,ket=gam_info)
  write (6,"(/t5,'KINETIC ENERGY INTEGRALS'/)")
  call gamess_print_1e_integrals(hao,bra=gam_info,ket=gam_info)
  
  nuclear_attraction: do iat=1,nCen
    xyz = real(gam_info%atoms(iat)%xyz,kind=kind(xyz)) / 0.5291772083 
    q   = real(gam_info%atoms(iat)%znuc,kind=kind(q))
    call gamess_1e_integrals('AO 3C 1/R',tmp_xk,bra=gam_info,ket=gam_info,op_xyz=xyz)
    tmp_xk = -q * tmp_xk
    write (6,"(/t5,'NUCLEAR ATTRACTION TO ATOM ',i3,' Z= ',f12.5,' XYZ=',3f12.5/)") iat, q, xyz
    call gamess_print_1e_integrals(tmp_xk,bra=gam_info,ket=gam_info)
    hao = hao + tmp_xk
  end do nuclear_attraction

  hao = matmul(sao,hao)
  hao = matmul(hao,transpose(sao))

  ! transform to MO basis
  hao = matmul(mos(1:gam_info%nbasis,:),hao)
  hmo = matmul(hao,transpose(mos(1:gam_info%nbasis,:)))

  print *,'one electron hamiltionian, mo basis: '
  call print_matrix(hmo,15,'f12.6')

  return
100 format(i3,i3,i3,i3,f16.10)
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

