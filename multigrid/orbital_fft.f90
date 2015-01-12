!
!  Multigrid test - calculation of matrix elements of a planewave
!
  module orbital_fft
    use accuracy
    use multigrid
    use qmech
    use fields
    use timer
    use import_gamess
    implicit none
    private
    public run_integrals
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: min_fields    = 1                   ! The rock bottom in the number of fields
    integer(ik), parameter :: max_orbitals  = 50
    integer(ik), parameter :: max_fields    = min_fields+max_orbitals
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_points(3)     = (/ 200, 200, 200 /)   ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    character(len=100) :: orbitals_file   = 'orbitals.dat'        ! Name of the file containing the MOs
    integer(ik)        :: orbital_count  = 0                      ! Number of molecular orbitals to build integrals for
    integer(ik)        :: orbital_list(max_orbitals)              ! List of orbital indices from the input file
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_oper                                  ! Space for the multiplicative operator
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_orbital(max_orbitals)                 ! List of orbtials, what else?
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /intdata/ verbose, &
                       n_points, box_extent, &
                       orbitals_file, &
                       orbital_count, orbital_list
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
      f_free = min_fields + orbital_count
      write (out,"(t5,'Allocating ',i4,' data fields.')") f_free
      write (out,"(t5,'if allocation fails, try decreasing orbital_count and/or n_points'/)")
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
    subroutine load_mos
      integer(ik)           :: imo
      real(rk)              :: norm
      !
      call TimerStart('Load MOs')
      !
      if (f_free<orbital_count) stop 'load_mos - out of fields'
      f_free = f_free - orbital_count
      f_orbital(:orbital_count) = f_table(f_free+1:f_free+orbital_count)
      !
      call FieldImport('GAMESS',orbitals_file,f_orbital(:orbital_count),orbital_list(:orbital_count))
      call describe_molecule
      !
      !  Renormalize the MOs, to compensate for grid-related noise
      !
      renormalize_mos: do imo=1,orbital_count
        call QMNormalize(f_orbital(imo),1.0_rk,norm)
        write (out,"(1x,a,' orbital ',i5,' was normalized to ',g15.8)") trim(orbitals_file), imo, norm
      end do renormalize_mos
      !
      call TimerStop('Load MOs')
    end subroutine load_mos
    !
    !  Problem driver
    !
    subroutine run_integrals
      integer(ik)        :: info
      real(rk)           :: kvec(3)
      integer(ik)        :: i_bra, i_ket
      complex(rk)        :: val, ovr
      !
      call TimerStart('Integrals')
      call accuracyInitialize
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
      call initialize_grid
      !
      call load_mos
      !
      write (out,"(/'Done with preliminaries'/)")
      !
      call TimerReport
      !
      !  We are ready to take matrix elements!
      !
      if (f_free<1) stop 'orbital_fft - no field for the operator!'
      f_oper = f_table(f_free) ; f_free = f_free - 1
      write (out,"(('#',3(1x,a14),2(1x,a5),4(1x,a20)))") &
      '   kx   ', '   ky   ', '   kz   ', ' bra ', ' ket ', '  Re[<|e^ikr|>]  ', '  Im[<|e^ikr|>]  ', '  Re[<|>]  ', ' Im[<|>]  ', &
      '  ----  ', '  ----  ', '  ----  ', '-----', '-----', ' --------------- ', ' --------------- ', ' --------- ', '--------- '
      !
      read_kvecs: do 
        read(input,*,iostat=info) kvec
        if (info/=0) exit read_kvecs
        !
        call FLsetPlanewave(k=kvec)
        call FieldInit(f_oper,FLplanewave)
        loop_bra: do i_bra=1,orbital_count
          loop_ket: do i_ket=1,i_bra
            val = FieldBraVKet(bra=f_orbital(i_bra),v=f_oper,ket=f_orbital(i_ket))
            ovr = FieldConjgIntegrate(left=f_orbital(i_bra),right=f_orbital(i_ket))
            write (out,"('@',3(1x,g14.6),2(1x,i5),4(1x,g20.12))") kvec, i_bra, i_ket, val, ovr
          end do loop_ket
        end do loop_bra
      end do read_kvecs
      write (out,"()")
      !
      call TimerStop('Integrals')
      call TimerReport
    end subroutine run_integrals

  end module orbital_fft
!
  subroutine driver
    use orbital_fft

    call run_integrals
  end subroutine driver
