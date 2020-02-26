
 subroutine errmsg(message)
  use channels, only: ilog
  implicit none
  character*144         :: message

  write(ilog,*)message
  stop '-- Program excited abnormally -- '
  return 
 end subroutine

 !
 ! read a gamess output file
 !
 subroutine read_gamess_output(naos,nbasis,nelec,ncen,nirr,eorb,orbsym,symlab,&
      pntgroup)
  use constants
  integer,intent(in)            :: nbasis,naos
  integer*4,intent(out)         :: nelec,ncen,nirr
  real(dp),intent(inout)        :: eorb(nbasis)
  integer*4,intent(out)         :: orbsym(nbasis)
  character*3,intent(out)       :: symlab(8)

  logical                       :: gexist
  integer*4                     :: cnt,i
  integer                       :: off1,off2
  integer                       :: ios
  integer                       :: orbfnd
  integer                       :: gamess=11
  real(dp)                      :: two = 2.
  real(dp)                      :: orbener(5)
  character(len=144)            :: line,scr
  character(len=3)              :: pntgroup

  inquire(file='gamess.log',exist=gexist)
  if(gexist) then
   open(unit=gamess,file='gamess.log')
  else
   scr = 'gamess.log does not exist.'
   call errmsg(scr)
  endif

  nirr   = -1
  ncen   = -1
  orbfnd = 0
  orbsym = 0
  symlab=''

  scan_lines: do while (orbfnd==0)
   call read_gamess_line(gamess,line)

   ! read the point group label
   if (index(line,'POINT GROUP').ne.0) then
      read(line,'(36x,a3)') pntgroup
   endif

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
       symlab(nirr) = trim(adjustl(line(off1-5:off1-3)))
       line = trim(adjustl(line(off1:len_trim(line))))
      enddo
     enddo scan_olabels
   endif

!   print *,'nbasis: ',nbasis
   ! read in orbital occupations and energies
   if (index(line,'------------').ne.0) then
      call read_gamess_line(gamess,line)
      if(index(line,'EIGENVECTORS').ne.0 .and. orbfnd==0) then
         call read_gamess_line(gamess,line) ! ---- format line
         call read_gamess_line(gamess,line) ! blank line
         off1=0
         off2=0
         scan_orbs: do while(index(line,'END OF RHF')==0)
            call read_gamess_line(gamess,line) ! orbital indices
            read(gamess,'(a15,5(f11.4))')scr,eorb(off1+1:off1+min(nbasis-off1,5))
            off1 = off1 + min(nbasis-off1,5)
            call read_gamess_line(gamess,line) ! irreps
            scan_irreps: do while(len_trim(line).gt.0)
               line = adjustl(line)
               cnt = 1
               do while(line(1:3) .ne. symlab(cnt))
                  cnt = cnt + 1
               enddo
               off2 = off2 + 1
               orbsym(off2) = cnt
               line = line(4:len_trim(line))
            enddo scan_irreps
!            do i = 1,nbasis
            do i = 1,naos
               call read_gamess_line(gamess,line) ! mos
            enddo
            call read_gamess_line(gamess,line) ! blank line/termination string
         enddo scan_orbs
         orbfnd=1
      endif      
   endif

  enddo scan_lines

  close(gamess)

 end subroutine read_gamess_output

 !
 ! subroutine to read a line of gamess output
 !
 subroutine read_gamess_line(unitn,line)
  use channels, only: ilog
  implicit none
  integer,intent(in)                    :: unitn
  character(len=144),intent(out)        :: line
  integer                               :: ios

  read(unitn,'(a144)',iostat=ios)line
  if(ios < 0) then
   write(ilog,*)'ERROR in reading gamess log file: '//trim(line)
   stop
  endif
  if(ios > 0) then
   write(ilog,*)'EOF reached in gamess log file.'
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
  use channels, only: ilog
  use vpqrsmod
  implicit none
  type(gam_structure)               :: gam ! gamess info (orbitals, geom, etc.) 
  type(int2e_cache)                 :: int2e    ! Currently active integrals context
  integer, parameter                :: iu_2e_ao   = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(hik)                      :: iosize_2e  = 120000000   ! Integral I/O

  real(xrk), dimension(:,:),allocatable                       :: mos_real,hao,hmo,fmo,smat
  complex(xrk),dimension(:,:),allocatable                     :: mos_cmplx,hao_cmplx
  character(len=clen)                                         :: ao_mode,int_type
  integer(ik)                                                 :: a,i,j
  integer(ik)                                                 :: nmo,nao,nao_spin,nocc_spin
  integer(ik)                                                 :: factable
  real(xrk)                                                   :: xyz(3), d_cf, q, ov, eps, refval, e1, e2
  real(xrk)                                                   :: nuc_repulsion
!  real(dp)                                                    :: vpqrs
  complex(xrk)                                                :: cz = (0._xrk,0._xrk),testnum

  ! make a call to MathLogFactorial before we get into the parallel parts, since
  ! it isn't openmp safe...
  factable = 30 ! function of the maximum angular momentum in atomic basis set
  q = MathFactorial(factable)
  q = MathLogFactorial(factable)
  q = MathDoubleFactorial(factable)

  eps      = 1.d-5
  nao      = gam%nbasis
  nao_spin = 2*nao
  nmo      = gam%nvectors

  ! allocate arrays
  allocate(occnum(nmo),roccnum(nmo))
  allocate(mos_cmplx(nao_spin,nmo),mos_real(nao_spin,nmo),hao_cmplx(nao_spin,nao_spin),hao(nao,nao),hmo(nmo,nmo),fmo(nmo,nmo),smat(nmo,nmo))
  allocate(x_dipole(nmo,nmo),y_dipole(nmo,nmo),z_dipole(nmo,nmo),dpl(nmo,nmo))

  occnum   = 0

  ! assume for the time being ao integrals stored in core
  ao_mode = 'incore'
  call prepare_2e(int2e,gam,ao_mode,iu_2e_ao,iosize_2e,ints_math='real')

  ! form 1e Hamiltonian (in AO basis)
  call core_hamiltonian(gam,hao)
  hao_cmplx                                = cz 
  hao_cmplx(1:nao,1:nao)                   = cmplx(hao(1:nao,1:nao),kind=xrk)
  hao_cmplx(nao+1:nao_spin,nao+1:nao_spin) = cmplx(hao(1:nao,1:nao),kind=xrk)

  ! read orbitals from gamess output
  mos_cmplx                          = cz
  if(nmo <= nao) then     ! Cartesian RHF case
   mos_cmplx(    1:nao     ,1:nmo)   = cmplx(gam%vectors(1:nao,1:nmo),kind=xrk)
  else if(nmo > nao) then ! Cartesian UHF case
   mos_cmplx(1:nao         ,1:nmo:2) = cmplx(gam%vectors(1:nao,1:nao),kind=xrk)
   mos_cmplx(nao+1:nao_spin,2:nmo:2) = cmplx(gam%vectors(1:nao,nao+1:nmo),kind=xrk)
  endif

  ! tighten up orbitals with a couple scf runs, if requested
  if (scfiter.gt.0) call scf_loop(int2e,gam,hao_cmplx,mos_cmplx)

  ! Reset the MOs in gam as these are used later in some parts of the code
  gam%vectors(1:nao,1:nmo)=real(real(mos_cmplx(1:nao,1:nmo)))
  
  ! ensure orbitals have not developed any complex character
  mos_real = real(real(mos_cmplx))
  refval = sum(real(aimag(mos_cmplx)))
  if(refval.gt.eps)then
   write(ilog,"(/'******* WARNING: imag part of mos_cmplx .gt. ',f14.10,': ',f14.10,' *********')")eps,refval
  endif

  ! Save the spatial orbitals for use later
  allocate(ao2mo(nao,nmo))
  ao2mo(1:nao,1:nmo)=real(mos_cmplx(1:nao,1:nmo))
  nbas_ao=nao

  ! load up variables depending on RHF/UHF case 
  if(nmo <= nao) then         ! RHF case, allowing for dropped orbitals
      nocc_spin = nelec / 2
      d_cf = 2._xrk
   else if(nmo == nao_spin) then ! UHF case
      nocc_spin = nelec
      d_cf = 1._xrk
   else
      stop 'load_mo_integrals - confusing number of mos'
   endif 

  ! actual call to generate integrals in MO basis 
  call transform_moint2e_real(int2e,moType,mos_cmplx,mos_cmplx,mos_cmplx,mos_cmplx,moIntegrals,io_unit=99,l_block=10)

  ! form 1e Hamiltonian in MO basis  (only take the real part)
  hmo = real(real((matmul(matmul(transpose(mos_cmplx),hao_cmplx),mos_cmplx))))
  refval = sum(real(aimag((matmul(matmul(transpose(mos_cmplx),hao_cmplx),mos_cmplx)))))
  if(refval.gt.eps)then
   write(ilog,"(/'******* WARNING: imag part of hmo .gt. ',f14.10,': ',f14.10,' *********')")eps,refval
  endif

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
   write(ilog,"(/t5,'FOCK MATRIX IN MO BASIS')")
   call print_matrix(fmo,15,'f15.10')
  endif

  if (debug.and.nmo <= nao) then
     call mo_energy_decomposition(nmo,nocc_spin,hmo)
  endif

  ! ensure Fock matrix is diagonal
  refval = 0._rk
  scan_row: do i = 1,nmo
   ! this check just looks for a mapping (i.e. major change due to SCF
   ! procedure) error. The number of figures in e(i) from reading the
   ! log is small enough that an error of 10*^-4 is reasonable.
   if(abs(e(i)-fmo(i,i)) > 0.001)write(ilog,"('warning: orbital ',i3,' energy mismatch - gamess.log=',f15.10,' internal=',f15.10)")i,e(i),fmo(i,i)
   e(i) = fmo(i,i)
   scan_col: do j = 1,nmo
    if(i /= j .and. abs(fmo(i,j)-refval)>eps)write(ilog,"('warning: <',i3,'|f|',i3,'> = ',f15.10,',expecting ',f15.10)")i,j,fmo(i,j),refval
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
  Ehf = e1 + e2 + nuc_repulsion(gam)

  ! in the future, we will fill orbitals using aufbau principle. However, I get
  ! the impression this code works primarily with spatial orbitals (i.e. under
  ! RHF). So, we trivially fill up the lowest nelec/2 orbitals, but this is
  ! easily changed in the future.
  do i = 1,nocc_spin
   occNum(i) = 2
  enddo

  ! load up the dipole moment integrals
  int_type = 'AO DIPOLE X'
  call dipole_integrals(gam,int_type,nao_spin,nmo,x_dipole,mos_real)
  int_type = 'AO DIPOLE Y'
  call dipole_integrals(gam,int_type,nao_spin,nmo,y_dipole,mos_real)
  int_type = 'AO DIPOLE Z'
  call dipole_integrals(gam,int_type,nao_spin,nmo,z_dipole,mos_real)

  ! print out summary (compare to GAMESS output)
  write(ilog,"(60('-'))")
  write(ilog,1000)
  write(ilog,1001)e1
  write(ilog,1002)e2
  write(ilog,1003)nuc_repulsion(gam)
  write(ilog,1004)Ehf
  write(ilog,"(60('-'))")

  ! clear integral cache
  call clear_2e(int2e)
 
  ! clean up all the allocated arrays
  deallocate(mos_real)
  deallocate(mos_cmplx)
  deallocate(hao_cmplx)
  deallocate(hao)
  deallocate(hmo)
  deallocate(fmo)

  return
1000 format(' Computation of HF energy in mo basis:')
1001 format('  one electron energy:      ',f17.10)
1002 format('  two electron energy:      ',f17.10)
1003 format('  nuclear repulsion energy: ',f17.10)
1004 format('  total energy:             ',f17.10)
 end subroutine load_mo_integrals

!#######################################################################
! rdgeom: reads the nuclear geometry from the GAMESS log file
!#######################################################################
 subroutine rdgeom(filename)

   use parameters
   use iomod, only: freeunit
   use parsemod, only: lowercase
   use channels

   implicit none
   
   integer                           :: igms,i,j
   real*8, dimension(:), allocatable :: x
   character(len=72)                 :: filename
   character(len=80)                 :: string

!-----------------------------------------------------------------------
! Open GAMESS output file
!-----------------------------------------------------------------------
   call freeunit(igms)
   open(igms,file=filename,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the Cartesian coordinates
!-----------------------------------------------------------------------
   ! Read to the coordinate section
5  read(igms,'(a)') string
   if (index(string,'COORDINATES').eq.0) goto 5
   read(igms,*)

   ! Determine the no. atoms and allocate arrays
   natm=0
10 read(igms,'(a)') string
   if (string.ne.'') then
      natm=natm+1
      goto 10
   endif
   
   allocate(x(natm*3))
   allocate(aatm(natm))

   ! Read the Cartesian coordinates and atom labels
   do i=1,natm+1
      backspace(igms)
   enddo
   do i=1,natm
      read(igms,'(1x,a2,10x,3(6x,F14.10))') aatm(i),&
           (x(j),j=i*3-2,i*3)
   enddo

!-----------------------------------------------------------------------
! Close the GAMESS output file
!-----------------------------------------------------------------------
   close(igms)

!-----------------------------------------------------------------------
! If the Cartesian coordinates have not been read from the input file
! then save the coordinates read from the gamess log file
!-----------------------------------------------------------------------
   ! Cartesian coordinates
   if (.not.allocated(xcoo)) then
      allocate(xcoo(natm*3))
      xcoo=x
   endif

   ! Atom labels.
   ! N.B. these have to be lowercase for use in other parts of the
   ! code.
   if (.not.allocated(atlbl)) then
      allocate(atlbl(natm))
      do i=1,natm
         call lowercase(aatm(i))
         atlbl(i)=aatm(i)
      enddo
   endif

!-----------------------------------------------------------------------
! Write the coordinates to the log file
!-----------------------------------------------------------------------
   write(ilog,'(/,82a)') ('-',i=1,82)
   write(ilog,'(27x,a)') 'Atomic Coordinates (Angstrom)'
   write(ilog,'(82a)') ('-',i=1,82)
   write(ilog,'(a4,22x,a1,18x,a1,20x,a1)') 'Atom','X','Y','Z'
   write(ilog,'(82a)') ('-',i=1,82)
   do i=1,natm
      write(ilog,'(1x,a2,10x,3(6x,F14.10))') aatm(i),&
           (x(j)*0.529177249d0,j=i*3-2,i*3)
   enddo
   write(ilog,'(82a,/)') ('-',i=1,82)

   return

 end subroutine rdgeom

!#######################################################################
 
 subroutine mo_energy_decomposition(nmo,nocc_spin,hmo)

   use accuracy
   use vpqrsmod
   use channels, only: ilog

   implicit none

   integer(ik)                  :: nmo,nocc_spin,i,j,a
   real(xrk),dimension(nmo,nmo) :: hmo,fmo
   real(xrk),dimension(nmo,nmo) :: jterm,kterm

!-----------------------------------------------------------------------
! Calculate the Fock matrix and its various contributions
!-----------------------------------------------------------------------
   jterm=0.0d0
   kterm=0.0d0

   fmo = hmo
   sum_orb: do a = 1,nocc_spin
      sum_row: do i = 1,nmo
         sum_col: do j = 1,nmo
            fmo(i,j) = fmo(i,j) + 2.0d0*vpqrs(i,j,a,a) - vpqrs(i,a,j,a)
            jterm(i,j)=jterm(i,j) + 2.0d0*vpqrs(i,j,a,a)
            kterm(i,j)=kterm(i,j) - vpqrs(i,a,j,a)
         enddo sum_col
      enddo sum_row
   enddo sum_orb

!-----------------------------------------------------------------------
! Output the orbital energies and their various contributions
!-----------------------------------------------------------------------
   write(ilog,'(/,75a)') ('-', i=1,75)
   write(ilog,'(26x,a)') 'MO Energy Decomposition'
   write(ilog,'(75a)') ('-', i=1,75)
   write(ilog,'(a)')'  i   | Total Energy  | Core Energy    |  &
        Coulomb Energy | Exchange Energy'
   write(ilog,'(75a)') ('-', i=1,75)
   do i=1,nmo
      write(ilog,'(2x,i3,4(2x,F15.10))') i,fmo(i,i),hmo(i,i),&
           jterm(i,i),kterm(i,i)
   enddo
   write(ilog,'(75a,/)') ('-', i=1,75)

   return

 end subroutine mo_energy_decomposition
