 !
 ! compute the 1e hamiltonian
 !
 subroutine core_hamiltonian(gam,hmat)
  use accuracy
  use math
  use parameters
  use integral_tools
  use import_gamess
  use channels, only: ilog
  implicit none
  type(gam_structure),intent(in)         :: gam ! gamess info (orbitals, geom,etc.) 
  real(xrk),intent(out)                  :: hmat(gam%nbasis,gam%nbasis)  ! Current 1-electron Hamiltonian matrix

  integer(ik)                            :: i,nao
  real(xrk)                              :: xyz(3),q
  real(xrk), allocatable                 :: tmp(:,:)         ! Overlap matrix (AObasis),null-space projected out

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
  use channels, only : ilog
  implicit none
  type(int2e_cache),intent(inout)    :: int2e    ! Currently active integrals context
  type(gam_structure),intent(in)     :: gam ! gamess info (orbitals, geom,etc.) 
  complex(xrk),intent(in)             :: hmat(2*gam%nbasis,2*gam%nbasis) ! 1e Hamiltonian
  complex(xrk),intent(inout)          :: mos_conv(2*gam%nbasis,gam%nvectors) ! initially contain guess orbitals converged mos

  type(diis_state)        :: diis_st           ! DIIS state, see diis.f90
  integer(ik)             :: itermax         = 30
  integer(ik)             :: max_diis_nvec   = 50      ! Maximum number of DIIS vectors allowed
  real(xrk)                :: max_diis_coeff  = 20._xrk  ! Restart DIIS if any of the coefficients exceed this threshold
  real(xrk)                :: eps_geev        = 1e-7_xrk        ! Threshold for declaring ZGEEV eigenvalues degenerate;
  integer(ik), parameter  :: iu_2e_ao        = 12 ! I/O unit used for storing 2e integrals over the atomic bfs
  integer(hik)            :: iosize_2e       = 220000000   ! Integral I/O
  real(xrk)                :: energy_toler    = 1e-12_xrk        ! Desired SCF convergence for the total energy
  real(xrk)                :: rho_toler       = 1e-8_xrk        ! Desired SCF convergence for the density matrix
  complex(xrk)             :: cz=(0._xrk,0._xrk)
 
  logical                 :: converged
  real(xrk)                :: drho
  real(xrk)                :: nuc_repulsion
  real(xrk)                :: null_cutoff
  integer(ik)             :: i,iter,nao_spin,nmo,nao,nvec,nmo_null,nmo_act
  complex(xrk)             :: efock,eg,escf, escf_old, enuc, de

  complex(xrk), allocatable    :: mos_init(:,:,:)     ! guess molecular orbitals
  complex(xrk), allocatable    :: mos(:,:,:)      ! molecular orbitals
  complex(xrk), allocatable    :: rho (:,:)       ! Electronic density matrix (AObasis)
  complex(xrk), allocatable    :: rho_old(:,:)    ! Electronic density matrixfrom previous SCF iteration
  complex(xrk), allocatable    :: fmat(:,:)       ! Fock matrix
  complex(xrk), allocatable    :: fmat_old(:,:)   ! Fock matrix from the previous iteration
  complex(xrk), allocatable    :: gmat(:,:)       ! 2-electron contribution tothe Fock matrix
  complex(xrk), allocatable    :: mo_energy(:)    ! MO orbital energies

  real(xrk), allocatable       :: tmp(:,:)        ! temporary matrix
  real(xrk), allocatable       :: smat(:,:)       ! Overlap matrix (AObasis),null-space projected out
  real(xrk), allocatable       :: sphalf(:,:)     ! S^{+1/2}, null-space isprojected out
  real(xrk), allocatable       :: smhalf(:,:)     ! S^{-1/2}, null-space isprojected out
  real(xrk), allocatable       :: mo_occ(:)       ! MO occupation vector

  null_cutoff = 1.e-6
  nao = gam%nbasis
  nvec = gam%nvectors
  nao_spin = 2*nao
  enuc = nuc_repulsion(gam)

  nmo = nvec                  ! number of spin orbitals, UHF case
  if(nvec <= nao)nmo = 2*nvec ! RHF case

  allocate(mo_occ(nao_spin),mo_energy(nao_spin),mos(nao_spin,nao_spin,2),mos_init(nao_spin,nao_spin,2),tmp(nao,nao), &
           rho(nao_spin,nao_spin),rho_old(nao_spin,nao_spin),fmat(nao_spin,nao_spin),fmat_old(nao_spin,nao_spin), &
           gmat(nao_spin,nao_spin),smat(nao_spin,nao_spin),sphalf(nao_spin,nao_spin),smhalf(nao_spin,nao_spin))

  ! initialize the converger
  call diis_initialize(diis_st,nao_spin,max_diis_nvec,max_diis_coeff,'real')

  ! load the overlap matrix into memory
  call gamess_1e_integrals('AO OVERLAP',tmp,bra=gam,ket=gam )
  smat                                = 0._xrk
  smat(1:nao,1:nao)                   = tmp(1:nao,1:nao)
  smat(nao+1:nao_spin,nao+1:nao_spin) = tmp(1:nao,1:nao)
  ! construct S^(1/2) and S^(-1/2)
  call st_invert_smat(nmo_null,smat,smhalf,sphalf,eps_smat=null_cutoff)

  ! load starting guess orbitals
  mos_init = cz
  if(nvec==nmo) then
    ! UHF case -- alpha and beta spin orbitals present
    mos_init(:,:,1) = mos_conv
  else
    ! RHF case -- only alpha orbitals present in mos_conv
    mos_init(1:nao         ,1:nmo:2,1) = mos_conv(1:nao,1:nvec)
    mos_init(nao+1:nao_spin,2:nmo:2,1) = mos_conv(1:nao,1:nvec)
  endif
  mos_init(:,:,2) = mos_init(:,:,1)

  ! load up orbital occupations.  Assume ordered by energy for now
  mo_occ = 0
  do i=1,nelec
    mo_occ(i) = 1._xrk
  enddo

  converged = .false.
  itermax = scfiter
  escf = cz
  rho = cz 
  fmat = cz
  mos = mos_init
  
  iterate_scf: do iter=1,itermax
     ! reset quantities from previous interation
     rho_old  = rho
     fmat_old = fmat
     escf_old = escf

     ! construct the density matrix (ao basis)
     call st_density_matrix(mo_occ,mos,rho)

     ! construct the gmatrix (2e contribution to fock matrix)
     gmat = cz
     call fock_g_matrix(int2e,rho,gmat)

     ! form new fock matrix
     fmat = hmat + gmat

     fock_diagonal: do i=1,nmo
       mo_energy(i) = dot_product(mos(:,i,1),matmul(fmat,mos(:,i,2)))
     end do fock_diagonal
     efock = sum(rho * fmat)
     eg    = sum(rho * gmat)
     escf  = efock - 0.5_xrk*eg + enuc

     ! extrapolate fock matrix and diagonalize
     call diis_extrapolate(diis_st,iter,smat,rho,fmat)

     call st_diagonalize_fmat(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy)

     ! update the (non-null space) MOs
     nmo_act = nao_spin - nmo_null
     call bt_follow_mos(nmo_act,mo_occ,eps_geev,sphalf,mos_init,mo_energy,mos)

!     if (iter<=1) cycle iterate_scf
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

 !  call print_matrix(realpart(mos(1:nao_spin,1:nao_spin,1)),11,'f10.5')

  if(nvec <= nao) then
   mos_conv = mos(1:nao_spin,1:nmo:2,1) ! for RHF case simply pull out the alpha orbitals
  else
   mos_conv = mos(:,:,1)                     ! for UHF case, take all orbitals
  endif

  deallocate(mo_occ,mo_energy,mos,mos_init,tmp,rho,rho_old,fmat,fmat_old,gmat,smat,sphalf,smhalf)

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

  nuc_repulsion = 0._xrk
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
   type(gam_structure),intent(in) :: gam ! gamess info (orbitals, geom, etc.) 
   character(clen),intent(in)     :: int_type
   integer(ik),intent(in)         :: nbas,nmo
   real(xrk),intent(out)          :: int_dipole(nmo,nmo)
   real(xrk),intent(in)           :: mos(nbas,nmo)
   real(xrk),allocatable          :: dao(:,:),dao_spin(:,:)
   integer(ik)                    :: nao,nao_spin

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


