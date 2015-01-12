!
!  Evaluation of "correlation pulse" potential. See: Z.B. Walters and O. Smirnova,
!  J Phys B 43, 161002 (2010) for the idea behind this stuff.
!
!
  module correlation_potential
    use dyson_tools

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !       Keep in mind that some of the relevant constants are in dyson_tools.f90
    !
    integer(ik), parameter :: max_slices = 100     ! Max number of slices allowed
    integer(ik), parameter :: max_rdm_sv = 200     ! Max number of singular values in an RDM
    !
    !  ==== User-adjustable parameters =====
    !
    character(len=clen) :: rdm_file = ' '          ! Name of the file containing rdm orbital coefficients
    character(len=clen) :: pot_file = ' '          ! Name of the file for the output potential slice
    integer(ik)         :: n_slices = 1            !
    integer(ik)         :: slice_axis (max_slices) ! An axis perpendicuar to the desired slice
    integer(ik)         :: slice_shift(max_slices) ! Displacement of this slice from grid mid-point
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)         :: f_trans_pot             ! Field for the transition potential
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /correlation/ verbose, &
                       n_points, box_extent, euler_angles, &
                       eps_hartree, sor_rate, &
                       rdm_file, pot_file, &
                       n_slices, slice_axis, slice_shift, &
                       total_fields
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Prepare ion-specific potentials: core+Hartree and transition potentials
    !
    subroutine load_transition_potential
      real(rk)    :: rdm_sv(max_rdm_sv)
      integer(ik) :: rdm_count
      !
      call TimerStart('Load transition potential')
      write (out,"(/'Loading ',a)") trim(rdm_file)
      call gamess_load_rdmsv(trim(rdm_file),rdm_sv,rdm_count)
      write (out,"( 'Found ',i4,' singular values')") rdm_count
      write (out,"( 'Values are: '/)")
      write (out,"(10(1x,f12.8))") rdm_sv(:rdm_count)
      !
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      if (f_free<1) stop 'not enough fields for transition potential'
      f_trans_pot = f_table(f_free) ; f_free = f_free - 1
      call eikonal_build_transition_potential(trim(rdm_file),rdm_sv(:rdm_count), &
                                              f_table(:f_free),f_trans_pot,rot=rotmat)
      call report_structure(0,trim(rdm_file))
      call TimerStop('Load transition potential')
    end subroutine load_transition_potential
    !
    subroutine dump_slices
      integer(ik)               :: is, dir, ix, iy, iz
      integer(ik)               :: npts(3), nslb(3)
      integer(ik)               :: li(3), ui(3)
      integer(ik)               :: pos
      complex(rk), allocatable  :: slice(:,:,:)
      real(rk), allocatable     :: x(:), y(:), z(:)
      !
      call TimerStart('dump slices')
      npts = FieldGridNPoints()
      !
      allocate(x(npts(1)),y(npts(2)),z(npts(3)))
      call FieldGridCoordinates(1,x)
      call FieldGridCoordinates(2,y)
      call FieldGridCoordinates(3,z)
      !
      open(unit=unit_dump,form='formatted',status='replace',file=trim(pot_file))
      loop_over_slices: do is=1,n_slices
        dir       = slice_axis(is)
        nslb      = npts
        nslb(dir) = 1    ! Only one point along defining direction
        !
        allocate (slice(nslb(1),nslb(2),nslb(3)))
        pos       = npts(dir)/2_ik + slice_shift(is)
        call FieldFetchSlab(f_trans_pot,dir=dir,n1=pos,n2=pos,data=slice)
        !
        li = 1 ; 
        ui = npts
        li(dir) = pos 
        ui(dir) = pos 
        !
        write (out,"('Writing slice ',i4,' axis ',i1,' displacement ',i4)") is, dir, slice_shift(is)
        write (out,"(' nslb = ',3i8)") nslb
        write (out,"('   li = ',3i8)") li
        write (out,"('   ui = ',3i8)") ui
        !
        write (unit_dump,"('# ',' Slice ',i3,' axis ',i1,' displacement ',i4)") is, dir, slice_shift(is)
        write (unit_dump,"('# ',a13,2x,a13,2x,a13,3x,a20,1x,a20)") &
               ' X ', ' Y ', ' Z ', ' Re(V) ', ' Im(V) '
        dump_z: do iz=li(3),ui(3)
          dump_y: do iy=li(2),ui(2)
            dump_x: do ix=li(1),ui(1)
              write (unit_dump,"(3(1x,f14.6),2x,2(1x,g20.14))") &
                     x(ix), y(iy), z(iz), slice(1+ix-li(1),1+iy-li(2),1+iz-li(3))
            end do dump_x
          end do dump_y
        end do dump_z
        write(unit_dump,"()")
        !
        deallocate(slice)
      end do loop_over_slices
      close(unit_dump)
      deallocate (x,y,z)
      call TimerStop('dump slices')
    end subroutine dump_slices
    !
    subroutine start
      integer(ik) :: info
      
      call TimerStart('start')

      total_fields = 8
      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=correlation,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=correlation)
      write (out,"()")
      !
      call TimerStart('setup_fields')
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Set number of required fields
      !
      tMin = 0 ; tMax = 0 ; dt = 1 ;
      call initialize_grid
      !
      call initialize_rotation
      !
      call load_transition_potential
      !
      if (verbose>=1) then
        call visualize_wavefunction('Transition potential '//trim(rdm_file),f_trans_pot)
      end if
      !
      call dump_slices
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
    
  end module correlation_potential
  !
  !
  !
  subroutine driver
    use correlation_potential
    use accuracy

    call accuracyInitialize

    call start

  end subroutine driver
