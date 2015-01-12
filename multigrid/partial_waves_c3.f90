!
!  Multigrid test - expansion of a 3D orbital in spherical partial waves
!                   (this is an essential step for an MO-ADK calculation)
!                   Cartesian version: insanely sensitive to basis set 
!                                      quality in the diffuse region.
!
  module partial_waves_c3
    use accuracy
    use multigrid
    use qmech
    use fields
    use timer
    use import_gamess
    implicit none
    private
    public run_partial
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_fields   = 2                   ! Max. number of fields used in the code
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_points(3)     = (/ 200, 200, 200 /)   ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    character(len=100) :: orbital_file    = 'orbital.dat'         ! Name of the file containing MO coefficients
    integer(ik)        :: orbital_index   = 1                     ! Orbital components to load from "orbital_file"
    real(rk)           :: orbital_ip      = 0.5_rk                ! Ionization potential - for the long-range tail
    real(rk)           :: ion_charge      = 1.0_rk                ! Charge of the ion core - for the long-range tail
    integer(ik)        :: wave_l_max      = 4                     ! Max. L value to scan in the expansion
    real(rk)           :: wave_r_min      = 5._rk                 ! Min. distance cut-off for the long-range tail
    real(rk)           :: wave_r_max      = 10._rk                ! Max. distance cut-off for the long-range tail
    character(len=20)  :: wave_origin     = 'orbital'             ! Origin of the long-range expansion centre. Can be:
                                                                  ! 'orbital' - centre of charge for the given orbital
                                                                  ! 'zero'    - the coordinate origin
                                                                  ! 'specify' - value entered in wave_r0
    real(rk)           :: wave_r0(3)      = (/ 0._rk, 0._rk, 0._rk /)
    logical            :: visualize       = .true.                ! General OpenDX visualization files
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: f_free                                  ! Index of the last unused position in
                                                                  ! f_table. All fields in f_table(1:f_free) are
                                                                  ! available for scratch use.
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: f_orbital                               ! Orbital to be expanded in partial waves
    integer(ik)        :: f_mask                                  ! Various masking functions
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /wavedata/ verbose, &
                       n_points, box_extent, &
                       orbital_file, orbital_index, &
                       orbital_ip, ion_charge, &
                       wave_origin, wave_r0, &
                       wave_l_max, wave_r_min, wave_r_max, &
                       visualize
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
      call MultiGridInit(max_grids=1,max_fields=max_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate all data fields
      !
      allocate_fields: do field=1,max_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      f_free = max_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Allocate fields for and load target Gamess orbital
    !
    subroutine load_gamess_mos
      real(rk)              :: norm
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call FieldImport('GAMESS',orbital_file,(/f_orbital/),(/orbital_index/))
      write (out,"('Loaded orbital ',i4,' from ',a)") orbital_index, trim(orbital_file)
      call QMNormalize(f_orbital,1.0_rk,norm)
      write (out,"('Squared norm before rescaling = ',f20.10)") norm**2
      !
      call gamess_report_nuclei(nnuc)
      write (out,"('Data file contained ',i5,' nuclei')") nnuc
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
    end subroutine load_gamess_mos
    !
    !  Determine the origin of the long-range expansion
    !
    subroutine calculate_wave_origin
      select case (wave_origin)
        case default
          write (out,"('Wave origin specification ""',a,'"" is not recognized')") trim(wave_origin)
          stop 'partial_waves%calculate_wave_origin - bad origin specification'
        case ('orbital','Orbital','ORBITAL')
          call FieldNormMultipoles(src=f_orbital,mult=wave_r0)
        case ('zero','Zero','ZERO')
          wave_r0 = 0._rk
        case ('specify','Specify','SPECIFY')
      end select
      write (out,"(/'Centre of the long-range expansion is at ',3f18.10/)") wave_r0
    end subroutine calculate_wave_origin
    !
    !  Divide the orbital by the expected long-range radial function
    !
    subroutine scale_long_range_tail
      real(rk) :: kappa
      !
      kappa = sqrt(2._rk*orbital_ip)
      write (out,"('Long-range exponent (kappa) = ',f18.12)") kappa
      call FLasymptoticSetParameters(z=ion_charge,kappa=kappa,r0=wave_r0)
      call FieldInit(f_mask,FLasymptotic)
      call FieldMul(src=f_mask,dst=f_orbital)
    end subroutine scale_long_range_tail
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
    !  Descriptive message
    !
    subroutine simulation_banner
      write (out,"()")
      write (out,"('  ===================================================================')")
      write (out,"('  One-centre expansion of the long-range tail of a molecular orbital ')")
      write (out,"('  in terms of spherical harmonics:')")
      write (out,"('  ')")
      write (out,"('     Psi(r,t,p) = Rad(r) Sum Sum C(L,M) Y(L,M,t,p)')")
      write (out,"('                          L   M')")
      write (out,"('  ')")
      write (out,"('  where the asymptotic radial dependence is in the form:')")
      write (out,"('  ')")
      write (out,"('     Rad(r) = r**(Z/kappa-1) exp(-kappa r)')")
      write (out,"('  ')")
      write (out,"('  where Z is the charge of the remaining ion, kappa in (2*Ip)**0.5,')")
      write (out,"('  and the phase of spherical harmonics is follows Landau&Lifshitz')")
      write (out,"('  conventions. The coefficients C(L,M) are calculated as:')")
      write (out,"('  ')")
      write (out,"('              /      *                      /  *')")
      write (out,"('     C(L,M) = | Psi Y (L,M) / Rad(r) d r  / | Y (L,M) Y(L,M) d r')")
      write (out,"('              /                             /')")
      write (out,"('  ')")
      write (out,"('  and both integrals are taken over a spherical shell of the inner')")
      write (out,"('  (outer) radius Rmin (Rmax).')")
      write (out,"('  ===================================================================')")
      write (out,"()")
    end subroutine simulation_banner
    !
    !  Problem driver
    !
    subroutine run_partial
      integer(ik)        :: info, vl, vm
      character(len=200) :: comment
      complex(rk)        :: ovr, nrm
      !
      call TimerStart('Partial Waves')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=wavedata,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=wavedata)
      write (out,"()")
      !
      call TimerStart('Preliminaries')
      call initialize_grid
      !
      if (f_free<2) stop 'partial_waves - not enough fields'
      f_orbital = f_table(f_free) ; f_free = f_free - 1
      f_mask    = f_table(f_free) ; f_free = f_free - 1
      !
      call load_gamess_mos
      if (visualize) call visualize_wavefunction('Initially loaded MO',f_orbital)
      !
      call calculate_wave_origin
      !
      call scale_long_range_tail
      if (visualize) call visualize_wavefunction('MO scaled by long-range tail',f_orbital)
      call TimerStop('Preliminaries')
      !
      !  Expand it!
      !
      call simulation_banner
      write (out,"(/' Expected norm of spherical harmonics is ',g18.10/)") &
             (1._rk/3._rk) * (wave_r_max**3 - wave_r_min**3)
      call TimerStart('Expansion')
      scan_spherical_l: do vl=0,wave_l_max
        scan_spherical_m: do vm=-vl,vl
          call FLharmonicsSetParameters(l=vl,m=vm,r0=wave_r0,rmin=wave_r_min,rmax=wave_r_max)
          call FieldInit(f_mask,FLharmonics)
          ovr = FieldConjgIntegrate(f_mask,f_orbital)
          nrm = FieldConjgIntegrate(f_mask,f_mask)
          write (out,"(' L= ',i2,' M= ',i3,' C= ',2f18.10,' (ovr= ',2g12.5,' nrm= ',2g12.5,')')") &
                 vl, vm, ovr/nrm, ovr, nrm
          if (visualize) then
            !
            !  Visualization sequence below is correct only for real orbitals!
            !
            call FieldMul(src=f_orbital,dst=f_mask)
            write (comment,"('Conjugate product: L=',i2,' M=',i3,' weight= ',2f18.10)") vl, vm, ovr/nrm
            call visualize_wavefunction(trim(comment),f_mask)
          end if
        end do scan_spherical_m
      end do scan_spherical_l
      call TimerStop('Expansion')
      !
      call TimerStop('Partial Waves')
      call TimerReport
    end subroutine run_partial

  end module partial_waves_c3
!
  subroutine driver
    use partial_waves_c3

    call run_partial
  end subroutine driver
