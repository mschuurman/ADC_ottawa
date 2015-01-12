
 real function vpqrs(r,s,u,v)
   use parameters
   implicit none
   integer,intent(in) :: r,s,u,v
   integer            :: moIndex

   vpqrs = moIntegrals(moIndex(r,s,u,v))
         
   return
 end function vpqrs

 integer function moIndex(r,s,u,v)
  use parameters 
  implicit none
  integer,intent(in)    :: r,s,u,v

  moIndex = 0
  moIndex = (r-1)*nBas**3 + (s-1)*nBas**2 + (u-1)*nBas + v  
  return
 end function moIndex


 subroutine errmsg(message)
  implicit none
  character*144         :: message

  write(6,*)message
  stop '-- Program excited abnormally -- '
  return 
 end subroutine

 subroutine phis_init
  implicit none

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
  use accuracy
  use import_gamess
  use integral_tools
  use integrals_mo2e
  implicit none
  type(gam_structure)               :: gam_info ! gamess info (orbitals, geom, etc.) 
  type(int2e_cache)                 :: int2e    ! Currently active integrals context
  type(moint2e_cache)               :: moint2e  ! Currently active MO integrals context
  integer, parameter                :: iu_2e_ao   = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(8)                       :: iosize_2e  = 220000000   ! Integral I/O
  complex(8),dimension(gam_info%nbasis,gam_info%nbasis)  :: mos
  character(len=clen)               :: mo_mode,ao_mode

  ao_mode = 'conventional'
  mo_mode = 'incore'
  mos = cmplx(gam_info%vectors)
  call prepare_2e(int2e,gam_info,ao_mode,iu_2e_ao,iosize_2e,ints_math='real')
  call transform_moint2e_real(int2e,mo_mode,mos,mos,mos,mos,moint2e) 
 
  return
 end subroutine load_mo_integrals

 ! 
 ! Determine the number of irreps, number of basis functions and number of atoms
 !
 subroutine phis_get_info(nirr,symlab,nbas,ncen)
  implicit none
  integer*4,intent(out)    :: nirr,nbas,ncen
  character*2,intent(out):: symlab(1024)
  integer                :: gamess=11,offset
  logical                :: gexist
  character*144          :: line,scr

  nirr = -1
  ncen = -1

  inquire(file='gamess.out',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.out')
  else 
   scr = 'gamess.out does not exist.'
   call errmsg(scr)
  endif  

  scan_lines: do while (nbas==-1 .or. ncen==-1 .or. nirr==-1)
   read(gamess,'(a144)')line
   if(index(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS').ne.0)then
     read(line,'(a47,i5)')scr,nbas
   endif
   if(index(line,'TOTAL NUMBER OF ATOMS').ne.0)then
     read(line,'(a47,i5)')scr,ncen
   endif
   if(index(line,'DIMENSIONS OF THE SYMMETRY SUBSPACES ARE').ne.0)then
     nirr=0
     scan_irreps: do 
      read(gamess,'(a144)')line
      if(index(line,'=').eq.0)exit
      offset=1
      do
       if(index(line,'=').eq.0)exit
       nirr = nirr + 1
       offset = index(line,'=')+1
       symlab(nirr) = trim(adjustl(line(offset-5:offset-2)))
       line = trim(adjustl(line(offset:len_trim(line))))
      enddo
     enddo scan_irreps
   endif
  enddo scan_lines
  close(gamess)
  
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
  integer                :: gamess=11
  integer                :: i,offset,orbfnd
  logical                :: gexist
  character*144          :: line,scr

  inquire(file='gamess.out',exist=gexist)
  if(gexist)then
   open(unit=gamess,file='gamess.out')
  else 
   scr = 'gamess.out does not exist.'
   call errmsg(scr)
  endif

  ehf    = 0.
  earr   = 0.
  orbfnd = 0
  scan_lines: do while(ehf==0. .or. orbfnd==0) 
   read(gamess,'(a144)')line 
   if(index(line,'FINAL RHF ENERGY IS').ne.0 .and. ehf==0.) then
    read(line,'(a20,f20.10,a20)')scr,ehf,scr
   endif
   if(index(line,'EIGENVECTORS').ne.0 .and. orbfnd==0) then
    read(gamess,'(a144)')line ! ---- format line
    read(gamess,'(a144)')line ! blank line
    offset=1
    scan_orbs: do while(index(line,'END OF RHF')==0)
      read(gamess,'(a144)')line  ! orbital indices
      read(gamess,'(a15,5(f11.4))')scr,earr(offset:offset+4)
      offset = offset + 5
      read(gamess,'(a144)')line  ! irreps
      do i = 1,nbas
       read(gamess,'(a144)')line  ! mos
      enddo
      read(gamess,'(a144)')line ! blank line/termination string
    enddo scan_orbs
    orbfnd=1
   endif  
  enddo scan_lines
  close(gamess)

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
  integer               :: gamess=11
  integer               :: i,orbfnd,offset
  character*144         :: line
  logical               :: gexist

  inquire(file='gamess.out',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.out')
  else 
   line = 'gamess.out does not exist.'
   call errmsg(line)
  endif

  orbfnd = 0
  scan_lines: do while(orbfnd==0)
   read(gamess,'(a144)')line
   if(index(line,'EIGENVECTORS').ne.0) then
    read(gamess,'(a144)')line ! ---- format line
    read(gamess,'(a144)')line ! blank line
    offset=0
    scan_orbs: do while(index(line,'END OF RHF')==0)
      read(gamess,'(a144)')line  ! orbital indices
      read(gamess,'(a144)')line
      read(gamess,'(a144)')line  ! irrep labels
      scan_irreps: do while(len_trim(line).gt.0)
       line = adjustl(line)
       i = 1
       do while(line(1:2) .ne. symlab(i))
        i = i + 1
       enddo
       offset = offset + 1
       orbsym(offset) = i
       line = line(3:len_trim(line))
      enddo scan_irreps
      do i = 1,nbas
       read(gamess,'(a144)')line  ! mos
      enddo
      read(gamess,'(a144)')line ! blank line/termination string
    enddo scan_orbs
    orbfnd=1
   endif
  enddo scan_lines
  close(gamess)

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
  integer               :: gamess=11
  integer               :: cnt
  real(d)               :: nelec
  logical               :: gexist
  character*144         :: scr,line

  inquire(file='gamess.out',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.out')
  else 
   scr = 'gamess.out does not exist.'
   call errmsg(scr)
  endif

  scan_lines: do
   read(gamess,'(a144)')line
   if(index(line,'NUMBER OF ELECTRONS').ne.0)then
     read(line,'(a47,f5.0)')scr,nelec
     exit
   endif
  enddo scan_lines
  close(gamess)

  cnt = 1
  do while(nelec.gt.0)
   occnum(cnt) = min(nelec,2.)
   nelec = nelec - occnum(cnt)
   cnt = cnt + 1
  enddo

  return
 end subroutine

 !
 ! Load up the 
 !
 subroutine phis_mc_dip(xdip,ydip,zdip,nbas)
  implicit none
  integer,intent(in)  :: nbas
  real,intent(inout)  :: xdip(nbas,nbas),ydip(nbas,nbas),zdip(nbas,nbas)

  xdip = 0.
  ydip = 0.
  zdip = 0.

  return
 end subroutine

 !
 ! Load the MO integrals into memory
 !
 subroutine phis_load_vpqrs
  implicit none

  return
 end subroutine

