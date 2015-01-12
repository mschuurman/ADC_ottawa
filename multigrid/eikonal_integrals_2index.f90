!
!  Multigrid test - calculation of one- and two-electron integrals involving one eikonal
!                   wavefunction and three MOs given by GTO expansion.
!
  module eikonal_integrals_2index
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use eikonal_tools
    implicit none
    private
    public run_integrals
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: min_fields    = 4                   ! The rock bottom in the number of fields
    integer(ik), parameter :: max_fields    = 50                  ! Let's hope we never try to allocate these many ...
    integer(ik), parameter :: max_naturals  = 50                  ! Max number of natural orbitals allowed
    integer(ik), parameter :: max_orbitals  = 50                  ! Doing 4-index orbitals for too many would already
                                                                  ! kill us - so it is not much of a restriction ...
    integer(ik), parameter :: unit_integral = 34                  ! A more or less random unit
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_points(3)     = (/ 200, 200, 200 /)   ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    character(len=100) :: natural_file    = 'natural.dat'         ! Name of the file containing MO coefficients
                                                                  ! for the natural orbitals of the cation
                                                                  ! (used for calculating Hartree potential)
    integer(ik)        :: natural_count   = 0                     ! Number of natural orbitals in the density matrix
    real(rk)           :: natural_occ(max_naturals)               ! Natural orbital occupation numbers
    real(rk)           :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                  ! and exchange potential(s)
    real(rk)           :: sor_rate        = 1.95_rk               ! Over-relaxation rate in solving Poisson equations
    logical            :: eikonal_pref    = .true.                ! Include eikonal pre-factor
    character(len=100) :: orbitals_file   = 'orbitals.dat'        ! Name of the file containing the MOs
    integer(ik)        :: orbitals_count  = 0                     ! Number of molecular orbitals to build integrals for
    integer(ik)        :: block_size      = 3                     ! Number of orbitals to load simultaneously. The cost
                                                                  ! of loading orbitals will decrease with the square
                                                                  ! of block_size, at the cost of 2*block_size extra
                                                                  ! memory fields
    real(rk)           :: eikonal_k(3)    = 0                     ! Desired momentum of the eikonal orbital
    integer(ik)        :: hole_orbital    = 1                     ! Hole orbital, fixed for a given process
    character(len=100) :: integrals_file  = 'integrals.dat'       ! Output file used to dump (unsorted) integrals
    character(len=100) :: v_xc            = ' '                   ! Exchange-correlation potential to use in construction of the
                                                                  ! eikonal functions. Possible choices are:
                                                                  ! ' ' or 'none' - no exchange potential
                                                                  ! 'Slater' - Dirac-Slater exchange
                                                                  ! 'SVWN'   - Slater exchange + VWN local correlation
    logical            :: plot_eikonal    = .false.               ! Plot eikonal solution before and after orthogonalization
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_core_pot                              ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_e_pot                                 ! Coulomb potential of the <i|k> cloud (4-index 2-e integrals)
    integer(ik)        :: f_eikonal                               ! Eikonal orbital (orthogonalized)
    integer(ik)        :: f_hole                                  ! Core hole orbital
    integer(ik)        :: f_rho                                   ! eikonal * core product density
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /intdata/ verbose, &
                       n_points, box_extent, &
                       eps_hartree, sor_rate, &
                       natural_count, natural_file, natural_occ, &
                       eikonal_pref, &
                       orbitals_count, orbitals_file, &
                       block_size, &
                       eikonal_k, hole_orbital, v_xc, &
                       plot_eikonal, &
                       integrals_file
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Prepare simulation box and allocate a sufficient number of data fields
    !
    subroutine initialize_grid
      real(rk)    :: box(2,3)
      integer(ik) :: field
      !
      call TimerStart('Grid initialization')
      !
      box(1,:) = -0.5_rk*box_extent
      box(2,:) =  0.5_rk*box_extent
      !
      write (out,"(//t5,'Simulation box size (Bohrs)'/)")
      write (out,"(t8,i4,' pts for X: ',2f14.6)") n_points(1), box(:,1)
      write (out,"(t8,i4,' pts for Y: ',2f14.6)") n_points(2), box(:,2)
      write (out,"(t8,i4,' pts for Z: ',2f14.6)") n_points(3), box(:,3)
      write (out,"()")
      !
      f_free = min_fields + 2*block_size
      write (out,"(t5,'Allocating ',i4,' data fields.')") f_free
      write (out,"(t5,'if allocation fails, try decreasing block_size and/or n_points'/)")
      call MultiGridInit(max_grids=1,max_fields=f_free,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate_fields: do field=1,f_free
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Be verbose about the molecule we load - but just once is enough
    !
    subroutine describe_molecule
      logical, save         :: have_run = .false.
      integer(ik)           :: nnuc, iat, alloc
      real(rk), allocatable :: xyzq(:,:)
      !
      if (have_run) return
      have_run = .true.
      !
      call gamess_report_nuclei(nnuc)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(orbitals_file), nnuc
      !
      allocate (xyzq(4,nnuc),stat=alloc)
      if (alloc/=0) then
        write (out,"('describe_molecule: Error ',i8,' allocating space for ',i6,' atoms')") alloc, nnuc
        stop 'eikonal_integrals%describe_molecule - out of memory'
      end if
      !
      call gamess_report_nuclei(nnuc,xyzq)
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,nnuc
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
      deallocate (xyzq)
    end subroutine describe_molecule
    !
    !  Load molecular orbitals defined GAMESS data file
    !
    subroutine load_mos(is,ie,block)
      integer(ik), intent(in) :: is       ! Index of the first orbital
      integer(ik), intent(in) :: ie       ! Index of the last orbital
      integer(ik), intent(in) :: block(:) ! List of field handles to fill with data
      !
      logical, save         :: norm_reported(max_orbitals) = .false.
      integer(ik)           :: mo_list(ie-is+1) ! List of the MOs to load
      integer(ik)           :: imo
      real(rk)              :: norm
      !
      call TimerStart('Load MOs')
      !
      if ( (ie-is+1)/=size(block) ) then
        write (out,"('load_mos: asked to fit orbitals ',i5,' - ',i5,' into '"// &
                   ",i3,'-element buffer')") is, ie, size(block)
        stop 'eikonal_integrals%load_mos - field set size mismatch'
      end if
      !
      fill_indices: do imo=is,ie
        mo_list(imo-is+1) = imo
      end do fill_indices
      !
      call FieldImport('GAMESS',orbitals_file,block,mo_list)
      call describe_molecule
      !
      !  Renormalize the MOs, to compensate for grid-related noise
      !
      renormalize_mos: do imo=is,ie
        call QMNormalize(block(imo-is+1),1.0_rk,norm)
        if (norm_reported(imo)) cycle renormalize_mos
        norm_reported(imo) = .true.
        write (out,"(1x,a,' orbital ',i5,' was normalized to ',g14.8)") &
               trim(orbitals_file), imo, norm
      end do renormalize_mos
      !
      call TimerStop('Load MOs')
    end subroutine load_mos
    !
    !  Compute the requested block of eikonal orbitals
    !
    subroutine load_eikonal
      integer(ik)   :: scr
      !
      call TimerStart('Load Eikonal')
      !
      if (f_free<1) stop 'eikonal_integrals%load_eikonal - out of fields'
      scr = f_table(f_free) ; f_free = f_free - 1
      !
      call eikonal_build_function(dst=f_eikonal,kvec=eikonal_k,potential=f_core_pot, &
                                  prefactor=eikonal_pref,scratch=scr,norm='energy')
      !
      f_free = f_free + 1
      ! 
      call TimerStop('Load Eikonal')
    end subroutine load_eikonal
    !
    !  Fill data fields with orbitals - either from the data file, or eikonal
    !  The first (orbitals_count) orbitals are from file; the rest are eikonal
    !
    subroutine fill_orbitals(is,ie,block)
      integer(ik), intent(in) :: is       ! Index of the first orbital
      integer(ik), intent(in) :: ie       ! Index of the last orbital
      integer(ik), intent(in) :: block(:) ! List of field handles to fill with data
      !
      integer(ik) :: n_mos     ! Number of MOs in this block
      !
      if ( (ie-is+1)>size(block) ) then
        write (out,"('fill_orbitals: asked to fill orbitals ',i5,' - ',i5,' into '"// &
                   ",i3,'-element buffer')") is, ie, size(block)
        stop 'eikonal_integrals%fill_orbitals - out of fields'
      end if
      !
      n_mos     = max(0,ie-is+1)
      !
      call load_mos(is,ie,block(1:n_mos))
    end subroutine fill_orbitals
    !
    !  Construct Coulomb potential of the <i|k> electron cloud
    !  The potential is implicitly stored in f_e_pot
    !
    subroutine build_e_e_potential(i,k)
      integer(ik), intent(in) :: i, k ! Orbital indices
      !
      integer(ik) :: bra, ket, rho
      !
      call TimerStart('4-i 2-e potential')
      if (f_free<3) stop 'eikonal_integrals%build_e_e_potential - no scratch fields'
      bra = f_table(f_free) ; f_free = f_free - 1
      ket = f_table(f_free) ; f_free = f_free - 1
      rho = f_table(f_free) ; f_free = f_free - 1
      !
      !  The two separate calls to fill_orbitals are potentially reading the
      !  external GAMESS file twice - but these calls are not on the critical
      !  inner loop, so we are unlikely to suffer for it.
      !
      call fill_orbitals(i,i,(/bra/))
      call fill_orbitals(k,k,(/ket/))
      call fock_exchange_potential(psic=bra,psi=ket,v=f_e_pot,product=rho)
      !
      f_free = f_free + 3
      call TimerStop('4-i 2-e potential')
    end subroutine build_e_e_potential
    !
    !  Build a <j|l> batch of 4-i 2-e integrals
    !
    subroutine build_ijkl_integrals(i,k)
      integer(ik), intent(in) :: i, k ! Integral labels for output
      !
      integer(ik) :: block_j(block_size)  ! Fields used for left-hand-side orbitals
      integer(ik) :: block_l(block_size)  ! Fields used for right-hand-side orbitals
      integer(ik) :: j_start, j_end, j    ! Indices of the current orbital set
      integer(ik) :: l_start, l_end, l
      complex(rk) :: int
      integer(ik) :: bra, ket             ! Field indices
      !
      call TimerStart('4-i 2-e batch')
      !
      !  Allocate data fields
      !
      if (f_free<2*block_size) stop 'eikonal_integrals%build_ijkl_integrals - out of scratch'
      block_j = f_table(f_free-block_size+1:f_free) ; f_free = f_free - block_size
      block_l = f_table(f_free-block_size+1:f_free) ; f_free = f_free - block_size
      !
      block_bra: do j_start=1,orbitals_count,block_size
        j_end = min(orbitals_count,j_start+block_size-1)
        call fill_orbitals(j_start,j_end,block_j)
        block_ket: do l_start=1,orbitals_count,block_size
          l_end = min(orbitals_count,l_start+block_size-1)
          call fill_orbitals(l_start,l_end,block_l)
          !
          !  Block of orbitals is ready for processing
          !
          loop_bra: do j=j_start,j_end
            bra = block_j(j-j_start+1)
            loop_ket: do l=l_start,l_end
              ket = block_l(l-l_start+1)
              int = FieldBraVKet(bra=bra,v=f_e_pot,ket=ket)
              !
              !  Print results for later processing
              !
              write (unit_integral,"(' %4',4(1x,i5),1x,e14.7,1x,e14.7)")  i, j, k, l, int
              !
            end do loop_ket
          end do loop_bra
        end do block_ket
      end do block_bra
      f_free = f_free + 2*block_size
      call TimerStop('4-i 2-e batch')
    end subroutine build_ijkl_integrals
    !
    !  Construct eikonal orbital and orthogonalize it to all occupied orbitals
    !
    subroutine orthogonalize_eikonal_orbital
      integer(ik) :: block_i(block_size)  ! Fields used for left-hand-side orbitals
      integer(ik) :: i_start, i_end, i    ! Indices of the current orbital set
      integer(ik) :: orb
      complex(rk) :: int
      !
      call TimerStart('Orthogonalized Eikonal')
      !
      !  Prepare the eikonbal orbital - no guarantees on orthogonality to
      !  occupied orbitals (yet).
      !
      call load_eikonal
      if (plot_eikonal) then
        call visualize_wavefunction('Eikonal function before projection',f_eikonal)
      end if
      !
      !  Allocate data fields
      !
      if (f_free<block_size) stop 'eikonal_integrals%orthogonalize_eikonal_orbital - out of scratch'
      block_i = f_table(f_free-block_size+1:f_free) ; f_free = f_free - block_size
      !
      ortho_block: do i_start=1,orbitals_count,block_size
        i_end = min(orbitals_count,i_start+block_size-1)
        call fill_orbitals(i_start,i_end,block_i)
        ortho: do i=i_start,i_end
          orb = block_i(i-i_start+1)
          int = FieldConjgIntegrate(left=orb,right=f_eikonal)
          call FieldAXPY(alpha=int,src=orb,dst=f_eikonal)
          if (verbose>=1) then
            write (out,"('Overlap <',i4,'|eikonal k> = ',2g15.8)") i, int
          endif
        end do ortho
      end do ortho_block
      f_free = f_free + block_size
      !
      if (plot_eikonal) then
        call visualize_wavefunction('Eikonal function after projection',f_eikonal)
      end if
      call TimerStop('Orthogonalized Eikonal')
    end subroutine orthogonalize_eikonal_orbital
    !
    !  Very simple visualization routine
    !
    subroutine visualize_wavefunction(text,psi)
      character(len=*), intent(in) :: text
      integer(ik), intent(in)      :: psi
      !
      call FieldVisualize(slot=0,src=psi,text=trim(text))
    end subroutine visualize_wavefunction
    !
    !  Problem driver
    !
    subroutine run_integrals
      integer(ik)        :: info
      !
      call TimerStart('Integrals')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=intdata,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=intdata)
      write (out,"()")
      !
      !  Error stop for parameter overflows
      !
      if (natural_count>max_naturals .or. orbitals_count>max_orbitals) then
        write (out,"('  natural_count = ',i8,' max. ',i8)") natural_count, max_naturals
        write (out,"(' orbitals_count = ',i8,' max. ',i8)") orbitals_count, max_orbitals
        stop 'eikonal_integrals%run_integrals - compile-time limits exceeded'
      end if
      !
      call initialize_grid
      !
      !  Construct the density and the Hartree potential - eikonal functions need them
      !
      if (natural_count>max_naturals) then
        write (out,"('Too many natural orbitals. Increase ''max_naturals''" &
                   // " to at least ',i0,' and recompile')") natural_count
        stop 'eikonal_integrals - too many natural orbitals'
      end if
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      if (f_free<2) stop 'eikonal_integrals - not enough fields for core potential'
      f_core_pot = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_potential(natural_file,natural_occ(:natural_count),f_table(:f_free),f_core_pot,v_xc=v_xc)
      !
      write (out,"(/'Done with preliminaries'/)")
      !
      open (unit_integral,form='formatted',status='replace',file=integrals_file)
      !
      call TimerReport
      !
      if (f_free<1) stop 'eikonal_integrals - not enough fields for the exchange potential'
      f_e_pot = f_table(f_free) ; f_free = f_free - 1
      !
      !  Temporarily allocate fields for the eikonal and core hole orbitals, as well as
      !  for the product density.
      !
      if (f_free<3) stop 'eikonal_integrals - not enough fields for the eikonal orbital'
      f_eikonal = f_table(f_free) ; f_free = f_free - 1
      f_hole    = f_table(f_free) ; f_free = f_free - 1
      f_rho     = f_table(f_free) ; f_free = f_free - 1
      !
      !  Construct the eikonal orbital
      !
      call orthogonalize_eikonal_orbital
      !
      !  Construct the core hole
      !
      call fill_orbitals(hole_orbital,hole_orbital,(/f_hole/))
      !
      !  Construct the Poisson solution for the (hole,eikonal) density
      !
      call fock_exchange_potential(psic=f_hole,psi=f_eikonal,v=f_e_pot,product=f_rho)
      !
      !  Release space for the eikonal & hole orbitals - they are no 
      !  longer needed.
      !
      f_free = f_free + 3
      !
      !  Construct the 2e integrals with the remaining occupied orbitals
      !
      call build_ijkl_integrals(hole_orbital,orbitals_count+1)
      !
      close (unit_integral)
      !
      call TimerStop('Integrals')
      call TimerReport
    end subroutine run_integrals

  end module eikonal_integrals_2index
!
  subroutine driver
    use eikonal_integrals_2index

    call run_integrals
  end subroutine driver
