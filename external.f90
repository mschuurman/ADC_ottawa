
 function vpqrs(r,s,u,v)
   use parameters
   implicit none
   integer,intent(in) :: r,s,u,v
   real(kind=8)            :: vpqrs

   ! we're going to hack this a bit. I think this code assumes spatial orbitals.
   ! Given our orbitals ordering, and assuming alpha = beta spatial orbitals
   ! (i.e. rhf), we should pull out the odd indices corresponding to the
   ! alpha spin orbitals. This is easily changed if need be
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
 subroutine read_gamess_output(nbasis,nelec,ncen,nirr,orbsym,symlab,ehf,earr)
  use constants
  integer,intent(in)            :: nbasis
  integer*4,intent(out)         :: nelec,nirr,ncen
  integer*4,intent(out)         :: orbsym(nbasis)
  real(d),intent(out)           :: ehf      ! The Hf energy
  real(d),intent(out)           :: earr(nbasis)  ! The orbital energies
  character*2,intent(out)       :: symlab(1024)

  logical                       :: gexist
  integer*4                     :: cnt
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
 ! 
 subroutine load_mo_integrals(gam_info)
  use accuracy
  use parameters
  use import_gamess
  use integral_tools
  use integrals_mo2e
  use printing
  implicit none
  type(gam_structure)               :: gam_info ! gamess info (orbitals, geom, etc.) 
  type(int2e_cache)                 :: int2e    ! Currently active integrals context
  integer, parameter                :: iu_2e_ao   = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(8)                        :: iosize_2e  = 220000000   ! Integral I/O

  complex(8),dimension(2*gam_info%nbasis,2*gam_info%nbasis)   :: mo_cmplx
  real(xrk), dimension(2*gam_info%nbasis,2*gam_info%nbasis)   :: mo_spin,hao_spin,hmo_spin,fmo_spin
  real(xrk), dimension(  gam_info%nbasis,  gam_info%nbasis)   :: sao,hao,tmp_xk
  real(xrk), dimension(  gam_info%nbasis,  gam_info%nvectors) :: mos
  character(len=clen)                                         :: mo_mode,ao_mode
  integer                                                     :: a,i,j,nvec,nmo,nao
  real(d)                                                     :: vpqrs
  real(xrk)                                                   :: xyz(3), q, ov, eps, refval, e1,e2, nuc_repulsion

  eps  = 1.d-5
  nao  = gam_info%nbasis
  nmo  = 2*nao
  nvec = gam_info%nvectors
  mos  = gam_info%vectors(1:nao, 1:nvec)

  ! Note: most other parts of multigrid assume cartesian basis functions, so,
  ! for rhf, nmo = nao, for uhf nmo = 2*nao.  While this is not the case in
  ! general, we will work under that assumption for the time being.

  ! If rhf, number of orbitals = nao
  if(nvec == nao) then
   ! alpha and beta spin orbitals are the same
   mo_spin(1:nao,    1:nmo-1:2) = mos(1:nao,1:nao)
   mo_spin(nao+1:nmo,2:nmo  :2) = mos(1:nao,1:nao)
  else if(nvec == 2*nao) then
   mo_spin(1:nao,    1:nmo-1:2) = mos(1:nao,1:nao)
   mo_spin(nao+1:nmo,2:nmo  :2) = mos(1:nao,nao+1:nmo)
  else
   stop 'number of mos != nao or 2*nao'
  endif
  mo_cmplx = cmplx(mo_spin,kind=xrk)

  ! Compute ao overlap integrals
  call gamess_1e_integrals('AO OVERLAP',sao,bra=gam_info,ket=gam_info  )
  if(debug) then
   write (6,"(/t5,'AO OVERLAP INTEGRALS'/)")
   call gamess_print_1e_integrals(sao,bra=gam_info,ket=gam_info)
  endif

  ! check orthonormality of GAMESS orbitals
  scan_left: do i = 1,nvec
   scan_right: do j = 1,nvec
    ov = dot_product(mos(:,i),matmul(sao,mos(:,j)))
    refval = 0._xrk
    if(i==j)refval = 1._xrk
    if(abs(ov-refval)>eps)write(6,"('warning: <',i3,'|',i3,'> = ',f15.10,', expecting ',f15.10)")i,j,ov,refval
   enddo scan_right
  enddo scan_left

  ! Compute kinetic energy integrals
  call gamess_1e_integrals('AO KINETIC',hao,bra=gam_info,ket=gam_info)
  if(debug) then
   write (6,"(/t5,'KINETIC ENERGY INTEGRALS'/)")
   call gamess_print_1e_integrals(hao,bra=gam_info,ket=gam_info)
  endif  

  ! compute nuclear attraction integrals
  nuclear_attraction: do i=1,nCen
    xyz = real(gam_info%atoms(i)%xyz,kind=kind(xyz)) / abohr 
    q   = real(gam_info%atoms(i)%znuc,kind=kind(q))
    call gamess_1e_integrals('AO 3C 1/R',tmp_xk,bra=gam_info,ket=gam_info,op_xyz=xyz)
    tmp_xk = -q * tmp_xk
    if(debug) then
     write (6,"(/t5,'NUCLEAR ATTRACTION TO ATOM ',i3,' Z= ',f12.5,' XYZ=',3f12.5/)") i, q, xyz
     call gamess_print_1e_integrals(tmp_xk,bra=gam_info,ket=gam_info)
    endif
    hao = hao + tmp_xk
  end do nuclear_attraction
 
  hao_spin                       = 0._xrk
  hao_spin(1:nao    , 1:nao)     = hao
  hao_spin(nao+1:nmo, nao+1:nmo) = hao

  ! form 1e Hamiltonian in spin-MO basis
  hmo_spin = matmul(matmul(transpose(mo_spin),hao_spin),mo_spin)

  ! this is for debugging only! for production -- move this to disk/conventional
  ao_mode = 'incore'
  mo_mode = 'incore'
  call prepare_2e(int2e,gam_info,ao_mode,iu_2e_ao,iosize_2e,ints_math='real')
  call transform_moint2e_real(int2e,mo_mode,mo_cmplx,mo_cmplx,mo_cmplx,mo_cmplx,moIntegrals)

  ! form Fock matrix in spin-MO basis (i.e. diagonal)
  fmo_spin = hmo_spin
  sum_row: do i = 1,nmo
   sum_col: do j = 1,nmo
    sum_orb: do a = 1,nelec
     fmo_spin(i,j) = fmo_spin(i,j) + vpqrs(i,j,a,a) - vpqrs(i,a,j,a)
    enddo sum_orb
   enddo sum_col
  enddo sum_row

  if(debug) then
   write(6,"(/t5,'FOCK MATRIX IN MO BAIS')")
   call print_matrix(fmo_spin,15,'f15.10')
  endif

  ! ensure Fock matrix is diagonal
  refval = 0._rk
  scan_row: do i = 1,nmo
   scan_col: do j = 1,nmo
    if(i /= j .and. abs(fmo_spin(i,j)-refval)>eps)write(6,"('warning: <',i3,'|f|',i3,'> = ',f15.10,',expecting ',f15.10)")i,j,fmo_spin(i,j),refval
   enddo scan_col
  enddo scan_row

  ! Compute HF energy and print
  e1 = 0._xrk
  e2 = 0._xrk
  sum_i: do i = 1,nelec
   e1 = e1 + hmo_spin(i,i)
   sum_j: do j = 1,nelec
    e2 = e2 + 0.5*(vpqrs(i,i,j,j) - vpqrs(i,j,i,j))
   enddo sum_j
  enddo sum_i
  
  ! in the future, we will fill orbitals using aufbau principle. However, I get
  ! the impression this code works primarily with spatial orbitals (i.e. under
  ! RHF). So, we trivially fill up the lowest nelec/2 orbitals, but this is
  ! easily changed in the future.
  do i = 1,nelec/2
   occNum(i) = int(2,kind=4)
  enddo

  ! print out summary (compare to GAMESS output)
  write(6,"(60('-'))")
  write(6,1000)
  write(6,1001)e1
  write(6,1002)e2
  write(6,1003)nuc_repulsion(gam_info)
  write(6,1004)e1+e2+nuc_repulsion(gam_info)
  write(6,"(60('-'))")
  call flush

  return
1000 format(' Computation of HF energy in mo basis:')
1001 format('  one electron energy:      ',f17.10)
1002 format('  two electron energy:      ',f17.10)
1003 format('  nuclear repulsion energy: ',f17.10)
1004 format('  total energy:             ',f17.10)
 end subroutine load_mo_integrals

 ! 
 ! compute nuclear repulsion
 !
 function nuc_repulsion(gam_info)
  use accuracy
  use import_gamess
  type(gam_structure)               :: gam_info ! gamess info (orbitals, geom,etc.) 
  real(xrk)                         :: nuc_repulsion
  integer                           :: i,j
  real(xrk)                         :: xyz(3),q,r

  loop_atom1: do i=1,gam_info%natoms-1
    xyz = real(gam_info%atoms(i)%xyz, kind=kind(xyz))
      q = real(gam_info%atoms(i)%znuc,kind=kind(q))
      loop_atom2: do j=i+1,gam_info%natoms
        r = sqrt(sum((real(gam_info%atoms(j)%xyz,kind=kind(xyz))-xyz)**2)) / abohr
        nuc_repulsion = nuc_repulsion + q*real(gam_info%atoms(j)%znuc,kind=kind(q)) / r
      enddo loop_atom2
   enddo loop_atom1

   return
 end function nuc_repulsion

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

