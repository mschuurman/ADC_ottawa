 !
 ! compute the 1e hamiltonian
 !
 subroutine core_hamiltonian(gam,hmat)
  use accuracy
  use math
  use parameters
  use integral_tools
  use import_gamess
  implicit none
  type(gam_structure),intent(in)         :: gam ! gamess info (orbitals, geom,etc.) 
  real(rk),intent(out)                   :: hmat(gam%nbasis,gam%nbasis)  ! Current 1-electron Hamiltonian matrix

  integer(ik)                            :: i,nao
  real(rk)                               :: xyz(3),q
  real(rk), allocatable                  :: tmp(:,:)         ! Overlap matrix (AObasis),null-space projected out

  nao = gam%nbasis
  allocate(tmp(nao,nao))
  hmat = 0.

  ! Compute kinetic energy integrals
  call gamess_1e_integrals('AO KINETIC',tmp,bra=gam,ket=gam)
  if(debug) then
   write (ilog,"(/t5,'KINETIC ENERGY INTEGRALS'/)")
   call gamess_print_1e_integrals(tmp,bra=gam,ket=gam)
  endif
  hmat = tmp

  ! compute nuclear attraction integrals
  nuclear_attraction: do i=1,nCen
    xyz = real(gam%atoms(i)%xyz,kind=kind(xyz)) / abohr
    q   = real(gam%atoms(i)%znuc,kind=kind(q))
    call gamess_1e_integrals('AO 3C 1/R',tmp,bra=gam,ket=gam,op_xyz=xyz)
    tmp = -q * tmp
    if(debug) then
     write (ilog,"(/t5,'NUCLEAR ATTRACTION TO ATOM ',i3,' Z=',f12.5,'XYZ=',3f12.5/)") i, q, xyz
     call gamess_print_1e_integrals(tmp,bra=gam,ket=gam)
    endif
    hmat = hmat + tmp
  end do nuclear_attraction

  deallocate(tmp)

  return
 end subroutine core_hamiltonian

 ! 
 ! Use the GAMESS orbitals that we've read in as a starting guess for computing
 ! HF orbitals.  Generally this means an iteration or two just to "tighten"
 ! things up.
 !
 subroutine scf_loop(int2e,gam,hmat,mos_conv)
  use accuracy
  use math
  use parameters
  use integral_tools
  use import_gamess
  use diis
  use biorthogonal_tools
  use fock_tools
  use scf_tools
  use printing
  implicit none
  type(int2e_cache),intent(inout)    :: int2e    ! Currently active integrals context
  type(gam_structure),intent(in)     :: gam ! gamess info (orbitals, geom,etc.) 
  complex(rk),intent(in)             :: hmat(2*gam%nbasis,2*gam%nbasis) ! 1e Hamiltonian
  complex(rk),intent(out)            :: mos_conv(2*gam%nbasis,gam%nvectors) ! converged mos

  type(diis_state)        :: diis_st           ! DIIS state, see diis.f90
  integer(ik)             :: itermax         = 10
  integer(ik)             :: max_diis_nvec   = 50      ! Maximum number of DIIS vectors allowed
  real(rk)                :: max_diis_coeff  = 20._rk  ! Restart DIIS if any of the coefficients exceed this threshold
  real(rk)                :: eps_geev        = 1e-7_rk        ! Threshold for declaring ZGEEV eigenvalues degenerate;
  integer(ik), parameter  :: iu_2e_ao        = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(hik)            :: iosize_2e       = 220000000   ! Integral I/O
  real(rk)                :: energy_toler    = 1e-12_rk        ! Desired SCF convergence for the total energy
  real(rk)                :: rho_toler       = 1e-8_rk        ! Desired SCF convergence for the density matrix

  logical                 :: converged
  real(rk)                :: drho
  real(xrk)                :: nuc_repulsion
  integer(ik)             :: i,iter,nao_spin,nmo,nao,nvec,nmo_null,nmo_act
  complex(rk)             :: efock,eg,escf, escf_old, enuc, de

  complex(rk), allocatable    :: mosg(:,:,:)     ! guess molecular orbitals
  complex(rk), allocatable    :: mos(:,:,:)      ! molecular orbitals
  complex(rk), allocatable    :: rho (:,:)       ! Electronic density matrix (AObasis)
  complex(rk), allocatable    :: rho_old(:,:)    ! Electronic density matrixfrom previous SCF iteration
  complex(rk), allocatable    :: fmat(:,:)       ! Fock matrix
  complex(rk), allocatable    :: fmat_old(:,:)   ! Fock matrix from the previous iteration
  complex(rk), allocatable    :: gmat(:,:)       ! 2-electron contribution tothe Fock matrix
  complex(rk), allocatable    :: mo_energy(:)    ! MO orbital energies

  real(rk), allocatable       :: tmp(:,:)        ! temporary matrix
  real(rk), allocatable       :: smat(:,:)       ! Overlap matrix (AObasis),null-space projected out
  real(rk), allocatable       :: sphalf(:,:)     ! S^{+1/2}, null-space isprojected out
  real(rk), allocatable       :: smhalf(:,:)     ! S^{-1/2}, null-space isprojected out
  real(rk), allocatable       :: mo_occ(:)       ! MO occupation vector


  nao = gam%nbasis
  nvec = gam%nvectors
  nao_spin = 2*nao
  enuc = nuc_repulsion(gam)
  nmo = nao_spin

  allocate(mo_occ(nmo),mo_energy(nmo),mos(nao_spin,nmo,2),mosg(nao_spin,nmo,2),tmp(nao,nao), &
           rho(nao_spin,nao_spin),rho_old(nao_spin,nao_spin),fmat(nao_spin,nao_spin),fmat_old(nao_spin,nao_spin), &
           gmat(nao_spin,nao_spin),smat(nao_spin,nao_spin),sphalf(nao_spin,nao_spin),smhalf(nao_spin,nao_spin))

  ! initialize the converger
  call diis_initialize(diis_st,nao_spin,max_diis_nvec,max_diis_coeff,'real')

  ! load the overlap matrix into memory
  call gamess_1e_integrals('AO OVERLAP',tmp,bra=gam,ket=gam )
  smat                                = 0.
  smat(1:nao,1:nao)                   = tmp(1:nao,1:nao)
  smat(nao+1:nao_spin,nao+1:nao_spin) = tmp(1:nao,1:nao)
  ! construct S^(1/2) and S^(-1/2)
  call st_invert_smat(nmo_null,smat,smhalf,sphalf)

  ! load gamess orbitals into starting guess as well as orbital occupations
  mosg(:,:,1) = 0
  mo_occ = 0
  if(nvec == nao) then    ! Cartesian RHF case
   mosg(    1:nao     ,1:nao_spin-1:2,1) = gam%vectors(1:nao,1:nao)
   mosg(nao+1:nao_spin,2:nao_spin  :2,1) = gam%vectors(1:nao,1:nao)
   do i = 1,nelec
    mo_occ(i) = 1.
   enddo
  else if(nvec == nao_spin) then ! Cartesian UHF case
   mosg(1:nao         ,1:nmo:2,1) = gam%vectors(1:nao,1:nao)
   mosg(nao+1:nao_spin,2:nmo:2,1) = gam%vectors(1:nao,nao+1:nmo)
   do i = 1,nelec
    mo_occ(i) = 1.
   enddo
  else                   ! Anything else..
   stop 'confusing number of MOs read from gamess output'
  endif
  mosg(:,:,2) = mosg(:,:,1)

  converged = .false.
  itermax = 10
  escf = 0
  rho = 0
  fmat = 0
  mos = mosg
  iterate_scf: do iter=1,itermax
     ! reset quantities from previous interation
     rho_old  = rho
     fmat_old = fmat
     escf_old = escf

     ! construct the density matrix (ao basis)
     call st_density_matrix(mo_occ,mos,rho)

     ! construct the gmatrix (2e contribution to fock matrix)
     gmat = 0
     call fock_g_matrix(int2e,rho,gmat)

     ! form new fock matrix
     fmat = hmat + gmat

     fock_diagonal: do i=1,nao_spin
       mo_energy(i) = dot_product(mos(:,i,1),matmul(fmat,mos(:,i,2)))
     end do fock_diagonal
     efock = sum(rho * fmat)
     eg    = sum(rho * gmat)
     escf  = efock - 0.5_rk*eg + enuc

     ! extrapolate fock matrix and diagonalize
     call diis_extrapolate(diis_st,iter,smat,rho,fmat)
     call st_diagonalize_fmat(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy)

     ! update the (non-null space) MOs
     nmo_act = nmo - nmo_null
     call bt_follow_mos(nmo_act,mo_occ,eps_geev,sphalf,mosg,mo_energy,mos)

     if (iter<=1) cycle iterate_scf
     de   = escf - escf_old
     drho = maxval(abs(rho-rho_old))
     write (out,"('SCF Iteration ',i4,' escf=',2(1x,g20.12),'de=',2(1x,g12.5),'drho= ',g12.5)") iter, escf, de, drho
     if (abs(de)<=energy_toler .and. drho<=rho_toler) then
       converged = .true.
       write (out,"('Final SCF energy = ',g24.14,1x,g24.14)") escf
       call flush(out)
       exit iterate_scf
     end if

  end do iterate_scf

  if(.not.converged)stop 'unable to determine converged orbitals'

  if(nvec /= nao_spin) then
   mos_conv = mos(1:nao_spin,1:nao_spin:2,1) ! for RHF case simply pull out the alpha orbitals
  else
   mos_conv = mos(:,:,1)                     ! for UHF case, take all orbitals
  endif

  deallocate(mo_occ,mo_energy,mos,mosg,tmp,rho,rho_old,fmat,fmat_old,gmat,smat,sphalf,smhalf)

  return
 end subroutine scf_loop

 ! 
 ! compute nuclear repulsion
 !
 function nuc_repulsion(gam)
  use accuracy
  use import_gamess
  type(gam_structure)               :: gam ! gamess info (orbitals, geom,etc.) 
  real(xrk)                         :: nuc_repulsion
  integer                           :: i,j
  real(xrk)                         :: xyz(3),q,r

  loop_atom1: do i=1,gam%natoms-1
    xyz = real(gam%atoms(i)%xyz, kind=kind(xyz))
      q = real(gam%atoms(i)%znuc,kind=kind(q))
      loop_atom2: do j=i+1,gam%natoms
        r = sqrt(sum((real(gam%atoms(j)%xyz,kind=kind(xyz))-xyz)**2)) / abohr
        nuc_repulsion = nuc_repulsion + q*real(gam%atoms(j)%znuc,kind=kind(q)) / r
      enddo loop_atom2
  enddo loop_atom1

   return
 end function nuc_repulsion

 !
 ! compute dipole moment integrals
 !
 subroutine dipole_integrals(gam,int_type,nbas,nmo,int_dipole,mos)
   use accuracy
   use integral_tools
   use import_gamess
   implicit none
   type(gam_structure),intent(in)  :: gam ! gamess info (orbitals, geom, etc.) 
   character(clen),intent(in)      :: int_type
   integer(ik),intent(in)          :: nbas,nmo
   real(rk),intent(out)            :: int_dipole(nmo,nmo)
   real(rk),intent(in)             :: mos(nbas,nmo)
   real(rk),allocatable            :: dao(:,:),dao_spin(:,:)
   integer(ik)                     :: nao,nao_spin

   nao      = gam%nbasis
   nao_spin = 2*nao
   allocate(dao(nao,nao),dao_spin(nao_spin,nao_spin))

   select case(int_type)
    case default; stop 'dipole moment integral not recognized.'
    case ('AO DIPOLE X')
      call gamess_1e_integrals('AO DIPOLE X',dao,bra=gam,ket=gam  )
    case ('AO DIPOLE Y')
      call gamess_1e_integrals('AO DIPOLE Y',dao,bra=gam,ket=gam  )
    case ('AO DIPOLE Z')
      call gamess_1e_integrals('AO DIPOLE Z',dao,bra=gam,ket=gam  )
   end select

   if(nbas == nao_spin) then
     dao_spin                       = 0._xrk
     dao_spin(1:nao    , 1:nao)     = dao
     dao_spin(nao+1:nmo, nao+1:nmo) = dao
     int_dipole = matmul(matmul(transpose(mos),dao_spin),mos)
     return
   else
     int_dipole = matmul(matmul(transpose(mos),dao),mos)
     return
   endif

 end subroutine dipole_integrals


