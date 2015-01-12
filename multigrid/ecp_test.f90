!
!  Simple sanity check for ECP projectors - evaluate ECP integrals for an atom.
!
  module ecp_test
    use accuracy
    use timer
    use multigrid
    use import_gamess
    use ecp_gamess
    implicit none
    private
    public start
    !
    !  User-controlled stuff
    !
    integer(ik)         :: verbose         = 2                     ! Verbosity level
    integer(ik)         :: n_points(3)     = (/ 100, 100, 100 /)   ! Number of sampling points along each direction
    real(rk)            :: box_extent(3)   = (/ 10., 10., 10. /)   ! Total size of the box
    integer(ik)         :: total_fields    = 4                     ! Total number of actual fields
    character(len=clen) :: ecp_file        = ' '                   ! Name of a file containing ECP
    real(rk)            :: ecp_eps_min     = 1e-3_rk               ! Small shift value cut-off (absolute) in ECP
    real(rk)            :: ecp_eps_max     = 1e2_rk                ! Large shift value cut-off (positive) in ECP
    real(rk)            :: ecp_eps_grid    = 0.1_rk                ! Characteristic grid spacing 
    character(len=clen) :: ecp_report      = ' '                   ! File to report ECP level-shift projectors to
    character(len=clen) :: mos_file        = ' '                   ! File containing MOs for which ECP integrals are needed
    integer(ik)         :: mos_count       = 0                     ! Number of MOs to evaluate integrals for
    !
    !  Internal stuff
    !
    integer(ik)              :: f_free                             ! Index of the last unused position in
                                                                   ! f_table. All fields in f_table(1:f_free) are
    integer(ik), allocatable :: f_table (:)                        ! List of all fields allocated
    type(ecp_molecule)       :: ecp                                ! The effective core potential of this molecule
    !
    namelist /ecp_test_data/ &
             verbose, &
             total_fields, n_points, box_extent, &
             ecp_file, ecp_eps_min, ecp_eps_max, ecp_eps_grid, ecp_report, &
             mos_file, mos_count
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
    !  Load and convert the ECP
    !
    subroutine load_ion_ecp
      type(gam_structure) :: gam
      !
      if (ecp_file==' ') return
      !
      call TimerStart('Prepare ECP')
      call gamess_load_orbitals(file=ecp_file,structure=gam)
      call ecp_convert(gam,ecp,ecp_eps_min,ecp_eps_max,ecp_eps_grid,ecp_report)
      call gamess_destroy(gam)
      call TimerStop('Prepare ECP')
      call TimerStart('Instantiate ECP')
      call ecp_build_cartesian_projector(ecp%ecps(1))
      call TimerStop('Instantiate ECP')
    end subroutine load_ion_ecp
    !
    !  This subroutine is a bit of a kludge. For each MOs of the test wavefunction,
    !  we'll evaluate the ECP projectors. The resulting overlaps are enough to
    !  figure out the ECP integrals.
    !
    subroutine do_ecp_integrals
      integer(ik)              :: mo, mo1, mon, nmo, mof
      integer(ik)              :: mol, mor
      complex(rk), allocatable :: ovr(:,:)
      integer(ik)              :: mo_ind(f_free)
      complex(rk)              :: ecp_int
      !
      allocate (ovr(ecp%ecps(1)%projectors%nvectors,mos_count))
      !
      load_mo_block: do mo1=1,mos_count,f_free
        mon = min(mo1+f_free-1,mos_count)
        nmo = mon - mo1 + 1
        fill_mos: do mo=mo1,mon
          mo_ind(mo-mo1+1) = mo
        end do fill_mos
        !
        call FieldImport('GAMESS',mos_file,f_table(:nmo),mo_ind(:nmo))
        !
        project_mos: do mo=mo1,mon
          mof = f_table(mo-mo1+1)
          call FieldECPProject(mof,ecp%ecps(1)%projectors,ecp%ecps(1)%grid,ecp%ecps(1)%rmax,ovr(:,mo))
          write (out,"(/'For the MO ',i5,' projectors are:')") mo
          write (out,"((t3,5(f16.10,1x,e12.5,2x)))") ovr(:,mo)
        end do project_mos
      end do load_mo_block
      !
      !  Do the integrals
      !
      write (out,"()")
      mor_loop: do mor=1,mos_count
        mol_loop: do mol=1,mor
          ecp_int = sum(conjg(ovr(:,mol))*ovr(:,mor)*ecp%ecps(1)%vshift)
          write (out,"('<',i4,'|ECP|',i4,'> = ',f16.8,1x,e12.5)") mol, mor, ecp_int
        end do mol_loop
      end do mor_loop
      !
      deallocate (ovr)
    end subroutine do_ecp_integrals
    !
    subroutine start
      call TimerStart('start')
      !
      read (input,nml=ecp_test_data)
      write (out,nml=ecp_test_data)
      write (out,"()")
      !
      call initialize_grid
      !
      call load_ion_ecp
      !
      if (ecp%necps/=1) stop 'need exactly one ECP-carrying centre. Sorry'
      !
      call do_ecp_integrals
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
  end module ecp_test
  !
  subroutine driver
    use ecp_test
    use accuracy
    !
    call accuracyInitialize
    call start
  end subroutine driver
