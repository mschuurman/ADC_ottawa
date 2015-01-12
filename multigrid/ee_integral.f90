!
!  Calculation of continuum-continuum matrix elements.
!  Note that ECPs are not supported for the scattering states.
!
  module ee_integral
    use accuracy
    use multigrid
    use qmech
    use fields
    use fock
    use timer
    use import_gamess
    use eikonal_tools
    use ecp_gamess
    use caps
    use math

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_naturals = 200                   ! Max number of natural orbitals allowed
    integer(ik), parameter :: max_rdm      = 200                   ! Max number of RDM orbitals allowed in a single rdm file
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)         :: verbose         = 2                     ! Verbosity level
    integer(ik)         :: n_points(3)     = (/ 100, 100, 100 /)   ! Number of sampling points along each direction
    real(rk)            :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    !
    real(rk)            :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                   ! and exchange potential(s)
    real(rk)            :: sor_rate        = 1.95_rk               ! Over-relaxation rate in solving Poisson equations
    character(len=clen) :: v_xc            = ' '                   ! Exchange-correlation potential to use in construction of the
                                                                   ! eikonal functions. Possible choices are:
                                                                   ! ' ' or 'none' - no exchange potential
                                                                   ! 'Slater' - Dirac-Slater exchange
                                                                   ! 'SVWN'   - Slater exchange + VWN local correlation
    character(len=clen) :: file_rdm        = 'rdm.dat'             ! File containing the RDM used for the integrals
    character(len=clen) :: file_vbra       = 'vbra.dat'            ! Scattering potential for the bra state
    character(len=clen) :: file_vket       = 'vket.dat'            ! Scattering potential for the ket state
                                                                   ! Blank file name will force plane wave scattering states
                                                                   ! Otherwise, eikonal scattering states will be generated.
    character(len=clen) :: eikonal_norm    = 'momentum'            ! Normalization of the scattering solution. Can be one of:
                                                                   ! 'natural', 'momentum', or 'energy'.
                                                                   ! See eikonal_build_function in eikonal_tools.f90 for more
    logical             :: eikonal_pref    = .true.                ! Include eikonal prefactor
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)              :: total_fields = 10_ik               ! Total number of actual fields required for the input value of num_channels
    integer(ik)              :: f_free                             ! Index of the last unused position in
                                                                   ! f_table. All fields in f_table(1:f_free) are
                                                                   ! available for scratch use.
    integer(ik), allocatable :: f_table (:)                        ! List of all fields allocated
    integer(ik)              :: f_vbra                             ! Field with bra potential. Not present if file_vbra==' '
    integer(ik)              :: f_vket                             ! Field with ket potential. Not present if file_vket==' '
    integer(ik)              :: f_vrdm                             ! Field with the transition potential.
    integer(ik)              :: f_bra                              ! Scattering solution
    integer(ik)              :: f_ket                              ! Scattering solution
    integer(ik)              :: f_pref                             ! Scratch space for eikonal prefactor
    complex(rk)              :: multipoles_bra(0:9)                ! Total multipole of the scattering potential
    complex(rk)              :: multipoles_ket(0:9)                ! Total multipole of the scattering potential
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /ee_data/ verbose, &
                       total_fields, n_points, box_extent, &
                       eps_hartree, sor_rate, v_xc, &
                       file_rdm, file_vbra, file_vket, &
                       eikonal_norm, eikonal_pref
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
      write (out,"('Total number of grid fields required = ',i10)") total_fields
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
      call MultiGridInit(max_grids=1,max_fields=total_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate (f_table(total_fields))
      allocate_fields: do field=1,total_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      f_free = total_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Report molecular geometry
    !
    subroutine report_structure(level,name)
      integer(ik), intent(in)      :: level ! Minimal verbosity level to print the structure;
                                            ! Report just the number of atoms otherwise
      character(len=*), intent(in) :: name  ! Name of the structure/file
      !
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call gamess_report_nuclei(nnuc)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      if (level>verbose) return
      !
      !  Tell a bit more!
      !
      allocate (xyzq(4,nnuc))
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
    end subroutine report_structure
    !
    !  Prepare ion-specific potentials: core+Hartree and transition potentials
    !
    subroutine load_ion_potentials
      integer(ik) :: natural_count
      real(rk)    :: natural_occ(max_naturals)
      !
      call TimerStart('Load ion potentials')
      !     
      !  Do we have natural orbitals to construct the density and the Hartree potential?
      !
      if (file_vbra/=' ') then
        write (out,"(/'Loading ',a)") trim(file_vbra)
        call gamess_load_natocc(file_vbra,natural_occ,natural_count)
        if (f_free<1) stop 'amplitudes - not enough fields for core potentials (1)'
        call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
        f_vbra = f_table(f_free) ; f_free = f_free - 1
        call eikonal_build_potential(file_vbra,natural_occ(:natural_count), &
                                     f_table(:f_free),f_vbra,v_xc=v_xc,multipoles=multipoles_bra)
        multipoles_bra(4:6) = 3*multipoles_bra(4:6) - sum(multipoles_bra(4:6))
        call report_structure(1_ik,trim(file_vbra))
      else
        f_vbra = -1
      end if
      !
      if (file_vket/=' ') then
        write (out,"(/'Loading ',a)") trim(file_vket)
        call gamess_load_natocc(file_vket,natural_occ,natural_count)
        if (f_free<1) stop 'amplitudes - not enough fields for core potentials (2)'
        call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
        f_vket = f_table(f_free) ; f_free = f_free - 1
        call eikonal_build_potential(file_vket,natural_occ(:natural_count), &
                                     f_table(:f_free),f_vket,v_xc=v_xc,multipoles=multipoles_ket)
        multipoles_ket(4:6) = 3*multipoles_ket(4:6) - sum(multipoles_ket(4:6))
        call report_structure(1_ik,trim(file_vket))
      else
        f_vket = -1
      end if
      !
      call TimerStop('Load ion potentials')
      !
      !  Load transition potential - this is something we didn't have in single-channel
      !
      call TimerStart('Load transition potential')
      write (out,"(/'Loading ',a)") trim(file_rdm)
      call gamess_load_rdmsv(file_rdm,natural_occ,natural_count)
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      if (f_free<1) stop 'setup_fields - not enough fields for transition potential'
      f_vrdm = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_transition_potential(file_rdm,natural_occ(:natural_count), &
                                              f_table(:f_free),f_vrdm)
      call report_structure(1_ik,trim(file_rdm))
      call TimerStop('Load transition potential')
    end subroutine load_ion_potentials
    !
    !  Initializes the fields
    !
    subroutine setup_fields
      call TimerStart('setup_fields')
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call initialize_grid
      !
      call load_ion_potentials
      !
      if (f_free<3) stop 'setup_fields - not enough fields for bra and ket'
      f_bra  = f_table(f_free) ; f_free = f_free - 1
      f_ket  = f_table(f_free) ; f_free = f_free - 1
      f_pref = f_table(f_free) ; f_free = f_free - 1
      !
      write (out,"()")
      write (out,"('[Grid fields free at end of setup_fields: ',i3,']')") f_free
      write (out,"()")
      !
      call TimerStop('setup_fields')
    end subroutine setup_fields
    !
    !
    !
    subroutine start
      integer(ik) :: info
      real(rk)    :: k_bra(3), k_ket(3)
      real(rk)    :: e_bra, e_ket
      complex(rk) :: braVket
      
      call TimerStart('start')

      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=ee_data,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=ee_data)
      write (out,"()")
      !
      call setup_fields

      call TimerReport

      write (out,"((1x,2(1x,3(a12,1x),1x,a10,1x),1x,2(a14,1x)))") &
             ' k_x(bra) ', ' k_y(bra) ', ' k_z(bra) ', ' e(bra), eV ', &
             ' k_x(ket) ', ' k_y(ket) ', ' k_z(ket) ', ' e(ket), eV ', &
             ' Re(V) ', ' Im(V) ', &
             '----------', '----------', '----------', '------------', &
             '----------', '----------', '----------', '------------', &
             '-------', '-------'
      read_kvectors: do
        read (input,*,iostat=info) k_bra, k_ket
        if (info/=0) exit read_kvectors
        !
        if (file_vbra/=' ') then
          call FLsetMultipolesPhase(mult=real(multipoles_bra,kind=rk))
          call eikonal_build_function(dst=f_bra,kvec=k_bra,potential=f_vbra, &
                                      prefactor=eikonal_pref,scratch=f_pref,norm=eikonal_norm)
        else
          call FLsetPlanewave(k=k_bra)
          call FieldInit(f_bra,FLplanewave)
        end if
        !
        if (file_vket/=' ') then
          call FLsetMultipolesPhase(mult=real(multipoles_ket,kind=rk))
          call eikonal_build_function(dst=f_ket,kvec=k_ket,potential=f_vket, &
                                      prefactor=eikonal_pref,scratch=f_pref,norm=eikonal_norm)
        else
          call FLsetPlanewave(k=k_ket)
          call FieldInit(f_ket,FLplanewave)
        end if
        !
        braVket = FieldBraVKet(bra=f_bra,v=f_vrdm,ket=f_ket)
        e_bra   = h2ev * 0.5_rk * sum(k_bra**2)
        e_ket   = h2ev * 0.5_rk * sum(k_ket**2)
        write (out,"('@',2(1x,3(f12.8,1x),1x,f10.4,1x),1x,2(g14.7,1x))") &
               k_bra, e_bra, k_ket, e_ket, braVket
      end do read_kvectors

      call TimerStop('start')
      
      call TimerReport
    end subroutine start
    
  end module ee_integral
  !
  !
  !
  subroutine driver
    use ee_integral
    use accuracy

    call accuracyInitialize

    call start

  end subroutine driver
