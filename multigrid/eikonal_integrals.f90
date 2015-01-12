!
!  Multigrid test - calculation of one- and two-electron integrals involving one eikonal
!                   wavefunction and three MOs given by GTO expansion.
!
  module eikonal_integrals
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
    integer(ik), parameter :: max_eikonals  = 50                  ! Max number of eikonal orbitals allowed
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
    logical            :: all_integrals   = .true.                ! Calculate all integrals, even if the eikonal function
                                                                  ! does not enter.
    integer(ik)        :: block_size      = 3                     ! Number of orbitals to load simultaneously. The cost
                                                                  ! of loading orbitals will decrease with the square
                                                                  ! of block_size, at the cost of 2*block_size extra
                                                                  ! memory fields
    integer(ik)        :: eikonal_count   = 0                     ! Number of eikonal orbitals to include
    real(rk)           :: eikonal_orbitals(4,max_eikonals)        ! Eikonal orbitals parameters: energy in eV and
                                                                  ! unnormalized direction vectors
    character(len=100) :: integrals_file  = 'integrals.dat'       ! Output file used to dump (unsorted) integrals
                  
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_core_pot                              ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_bare_pot                              ! Total potential of the molecular core (electrons+nuclei)
    integer(ik)        :: f_e_pot                                 ! Coulomb potential of the <i|k> cloud (4-index 2-e integrals)
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /intdata/ verbose, &
                       n_points, box_extent, &
                       eps_hartree, sor_rate, &
                       natural_count, natural_file, natural_occ, &
                       eikonal_pref, &
                       orbitals_count, orbitals_file, &
                       all_integrals, block_size, &
                       eikonal_count, eikonal_orbitals, &
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
    subroutine load_eikonal(is,ie,block)
      integer(ik), intent(in) :: is       ! Index of the first eikonal orbital
      integer(ik), intent(in) :: ie       ! Index of the last eikonal orbital
      integer(ik), intent(in) :: block(:) ! List of field handles to fill with data
      !
      integer(ik)   :: imo
      integer(ik)   :: scr
      real(rk)      :: kvec(3) ! Asymptotic momentum of the eikonal function
      real(rk)      :: lambda  ! Wavelength of the eikonal function
      logical, save :: eikonal_reported(max_eikonals) = .false.
      !
      call TimerStart('Load Eikonal')
      !
      if ( (ie-is+1)/=size(block) ) then
        write (out,"('load_eikonal: asked to fit orbitals ',i5,' - ',i5,' into '"// &
                   ",i3,'-element buffer')") is, ie, size(block)
        stop 'eikonal_integrals%load_eikonal - field set size mismatch'
      end if
      !
      if (f_free<1) stop 'eikonal_integrals%load_eikonal - out of fields'
      scr = f_table(f_free) ; f_free = f_free - 1
      !
      build_eikonal: do imo=is,ie
        kvec   = eikonal_orbitals(2:4,imo)
        kvec   = kvec / sqrt(sum(kvec**2))
        kvec   = kvec * sqrt(2*eikonal_orbitals(1,imo)/h2ev)
        lambda = twopi / sqrt(sum(kvec**2))
        call eikonal_build_function(dst=block(imo-is+1),kvec=kvec,potential=f_core_pot, &
                                    prefactor=eikonal_pref,scratch=scr,norm='energy')
        if (eikonal_reported(imo)) cycle build_eikonal
        eikonal_reported(imo) = .true.
        write (out,"('Eikonal orbital ',i3,' E= ',f10.3,' eV; k= ',3(g12.6,1x),' au; l= ',f12.6,' Bohr')") &
               imo, eikonal_orbitals(1,imo), kvec, lambda
      end do build_eikonal
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
      integer(ik) :: ie_mo     ! Last MO in this block
      integer(ik) :: n_mos     ! Number of MOs in this block
      integer(ik) :: n_eikonal ! Number of eikonal orbitals in this block
      !
      if ( (ie-is+1)>size(block) ) then
        write (out,"('fill_orbitals: asked to fill orbitals ',i5,' - ',i5,' into '"// &
                   ",i3,'-element buffer')") is, ie, size(block)
        stop 'eikonal_integrals%fill_orbitals - out of fields'
      end if
      !
      ie_mo     = min(ie,orbitals_count)
      n_mos     = max(0,ie_mo-is+1)
      n_eikonal = (ie-is+1) - n_mos
      !
      if (n_mos>0) then
        call load_mos(is,ie_mo,block(1:n_mos))
      end if
      !
      if (n_eikonal>0) then
        call load_eikonal(is+n_mos-orbitals_count,ie-orbitals_count,block(n_mos+1:n_mos+n_eikonal))
      end if
    end subroutine fill_orbitals
    !
    !  Calculate two-index integrals: kinetic energy, overlaps, and bare-nuclei attraction
    !
    subroutine integrals_two_index
      integer(ik) :: orbitals_total       ! Total number of orbitals: MOs + eikonal
      integer(ik) :: block_i(block_size)  ! Fields used for left-hand-side orbitals
      integer(ik) :: block_j(block_size)  ! Fields used for right-hand-side orbitals
      integer(ik) :: i_start, i_end, i    ! Indices of the current orbital set
      integer(ik) :: j_start, j_end, j
      complex(rk) :: sint, vint, tint
      integer(ik) :: ic, nc               ! Gradient components
      integer(ik) :: bra, ket             ! Field indices
      integer(ik) :: gr_bra, gr_ket       ! Scratch fields for gradient components
      !
      call TimerStart('2-index integrals')
      write (out,"(/' Calculating 2-index integrals (S,T,V)'/)")
      !
      !  Allocate data fields
      !
      if (f_free<2*block_size) stop 'eikonal_integrals%integrals_two_index - out of scratch (A)'
      block_i = f_table(f_free-block_size+1:f_free) ; f_free = f_free - block_size
      block_j = f_table(f_free-block_size+1:f_free) ; f_free = f_free - block_size
      !
      orbitals_total = orbitals_count+eikonal_count
      nc = FieldComponentCount()
      write (unit_integral,"(('    ',a5,1x,a5,3(1x,a14,1x,a14,1x)))")  &
             ' i ', ' j ', '<i|j>', ' ', '<i|Vnuc|j>', ' ', '<i|Tkin|j>', ' ', &
             '---', '---', '-----', ' ', '----------', ' ', '----------', ' '
      block_bra: do i_start=1,orbitals_total,block_size
        write (out,"(' processing bra vector ',i5,' (out of ',i5,')')") i_start, orbitals_total
        i_end = min(orbitals_total,i_start+block_size-1)
        call fill_orbitals(i_start,i_end,block_i)
        block_ket: do j_start=1,orbitals_total,block_size
          j_end = min(orbitals_total,j_start+block_size-1)
          call fill_orbitals(j_start,j_end,block_j)
          !
          !  Block of orbitals is ready for processing
          !
          loop_bra: do i=i_start,i_end
            bra = block_i(i-i_start+1)
            loop_ket: do j=j_start,j_end
              ket = block_j(j-j_start+1)
              if (f_free<2) stop 'eikonal_integrals%integrals_two_index - out of scratch (B)'
              gr_bra = f_table(f_free) ; f_free = f_free - 1
              gr_ket = f_table(f_free) ; f_free = f_free - 1
              sint = FieldConjgIntegrate(left=bra,right=ket)
              vint = FieldBraVKet(bra=bra,v=f_bare_pot,ket=ket)
              tint = 0
              laplacian: do ic=1,nc
                call FieldGradientComponentRight(dir=ic,src=bra,dst=gr_bra)
                call FieldGradientComponentRight(dir=ic,src=ket,dst=gr_ket)
                tint = tint + FieldConjgIntegrate(left=gr_bra,right=gr_ket)
              end do laplacian
              f_free = f_free + 2
              tint = 0.5_rk * tint
              !
              !  Print results for later processing
              !
              write (unit_integral,"(' %2 ',i5,1x,i5,3(1x,e14.7,1x,e14.7,1x))")  i, j, sint, vint, tint
              !
            end do loop_ket
          end do loop_bra
        end do block_ket
      end do block_bra
      f_free = f_free + 2*block_size
      call TimerStop('2-index integrals')
    end subroutine integrals_two_index
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
      integer(ik) :: orbitals_total       ! Total number of orbitals: MOs + eikonal
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
      orbitals_total = orbitals_count+eikonal_count
      block_bra: do j_start=1,orbitals_total,block_size
        j_end = min(orbitals_total,j_start+block_size-1)
        call fill_orbitals(j_start,j_end,block_j)
        block_ket: do l_start=1,orbitals_total,block_size
          l_end = min(orbitals_total,l_start+block_size-1)
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
    !  Calculate four-index electron repulsion integrals
    !
    subroutine integrals_four_index
      integer(ik) :: orbitals_total  ! Total number of orbitals: MOs + eikonal
      integer(ik) :: i, k            ! Indices of the r1 orbitals
      !
      call TimerStart('4-index integrals')
      write (out,"(/' Calculating 4-index 2-electron integrals'/)")
      orbitals_total = orbitals_count+eikonal_count
      !
      if (f_free<1) stop 'eikonal_integrals%integrals_four_index - no scratch'
      f_e_pot = f_table(f_free) ; f_free = f_free - 1
      !
      write (unit_integral,"(('   ',4(1x,a5),1x,a14,1x,a14))")  &
          ' i ', ' j ', ' k ', ' l ', '<ij||kl>', ' ', &
          '---', '---', '---', '---', '--------', ' '
      loop_i: do i=1,orbitals_total
        loop_k: do k=1,orbitals_total
          write (out,"(/' processing orbital pair ',2i5,' (',f6.2,'%)')") &
                 i, k, 100._ark*(real(i-1,ark)*orbitals_total+k-1)/real(orbitals_total,ark)**2
          !
          !  Start by calculating the (possibly complex) product density, and the
          !  corresponding interaction potential.
          !
          call build_e_e_potential(i,k)
          !
          !  Now, the integrals involving the remaining two indices can be calculated
          !  in the same way as the 2-electron nuclear attraction integrals
          !
          call build_ijkl_integrals(i,k)
        end do loop_k
      end do loop_i
      !
      f_free = f_free + 1
      call TimerStop('4-index integrals')
    end subroutine integrals_four_index
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
      f_bare_pot = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_potential(natural_file,natural_occ(:natural_count), &
                                   f_table(:f_free),f_core_pot,f_bare_pot)
      !
      write (out,"(/'Done with preliminaries'/)")
      !
      open (unit_integral,form='formatted',status='replace',file=integrals_file)
      !
      !  Do the easy two-index integrals first: the overlaps, kinetic,
      !  and potential energy terms over the MOs (and the eikonal continuum)
      !
      call integrals_two_index
      call TimerReport
      !
      !  After the warm-up, do the 4-index 2-electron integrals. This one
      !  is going to cost A LOT :(
      !
      call integrals_four_index
      !
      close (unit_integral)
      !
      call TimerStop('Integrals')
      call TimerReport
    end subroutine run_integrals

  end module eikonal_integrals
!
  subroutine driver
    use eikonal_integrals

    call run_integrals
  end subroutine driver
