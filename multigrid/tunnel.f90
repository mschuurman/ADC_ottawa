!
!  Multigrid test - tunnel ionization rates of small many-electron molecules
!
!  The tunnel driver needs converged Hartree-Fock MOs and corresponding
!  eigenvalues on input (effectively fixing the Fock operator). The orbitals
!  must be calculated in the presence of the -same- electric field which
!  will be used in ionization rate calculation.
!
  module tunnel
    use accuracy
    use multigrid
    use qmech
    use fields
    use timer
    use import_gamess
    use fock
    use caps
    implicit none
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_mos     = 20                  ! Max number of MOs we support
    integer(ik), parameter :: max_fields  = 10+max_mos          ! Max number of data fields we support
    integer(ik), parameter :: max_nuclei  = 10                  ! Max number of nuclei we accept
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)        :: verbose         = 2                     ! Verbosity level
    integer(ik)        :: n_mos           = 1                     ! Actual number of the MOs
    integer(ik)        :: n_points(3)     = (/ 50, 50, 50 /)      ! Number of sampling points along each direction
    real(rk)           :: box_extent(3)   = (/ 20., 20., 20. /)   ! Total size of the box
    real(rk)           :: eps_hartree     = 1e-8_rk               ! Desired convergence in the Hartree potential
                                                                  ! and exchange potential(s)
    real(rk)           :: sor_rate        = 1.9_rk                ! Over-relaxation rate in solving Poisson equations
    integer(ik)        :: e_axis          = 3                     ! Axis along which electric field is applied
    real(rk)           :: e_field         = 0.                    ! Constant electric field strendth
    integer(ik)        :: sink_order      = 6                     ! Order of the absorbing potential at the boundary
    real(rk)           :: sink_en         = 2._rk                 ! Characteristic sink energy
    real(rk)           :: sink_width      = 5._rk                 ! Width of the absorbing potential
    real(rk)           :: src_rate        = 0.0001                ! Initial value of the source rate
    character(len=100) :: gamess_file = 'mos.dat'                 ! Name of the file containing MO coefficients
    integer(ik)        :: mo_ind(max_mos)                         ! Indices of the MOs in the data file
    real(rk)           :: mo_occ(max_mos)                         ! Occupation of the MOs in the original molecule
                                                                  ! Active electron is assumed to originate from 
                                                                  ! the last MO in the list. Note that only 0 or 2
                                                                  ! are in fact expected here!
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)        :: n_fields                                ! Number of fields allocated
    integer(ik)        :: n_free                                  ! Number of fields still free 
                                                                  ! Free fields are at the beginning of f_table(:)
    integer(ik)        :: f_table (max_fields)                    ! List of all fields allocated
    integer(ik)        :: mo_table(max_mos)                       ! List of fields used by the MOs
    integer(ik)        :: v_external                              ! External potential in the Fock operator
    integer(ik)        :: v_hartree                               ! Hartree potential in the Fock operator
    integer(ik)        :: v_sink                                  ! Multiplicative sink potential
    integer(ik)        :: rho_total                               ! Total electron density of the molecule
    integer(ik)        :: n_nuclei                                ! Number of nuclei
    real(rk)           :: xyzq(4,max_nuclei)                      ! Coordinates and charges of the nuclei
    real(rk)           :: e_input                                 ! Hartree-Fock energy of the input wavefunction
    real(rk)           :: eps_input(max_mos,max_mos)              ! Lagrangian multipliers for the input MOs
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /tunneldata/ verbose, &
                          n_mos, &
                          n_points, box_extent, &
                          e_axis, e_field, &
                          sink_order, sink_en, sink_width, &
                          src_rate, &
                          eps_hartree, sor_rate, &
                          gamess_file, &
                          mo_ind, mo_occ
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
      n_fields = n_mos + 13
      if (n_fields>max_fields) then
        write (out,"('Need ',i5,' fields, but compile-time max_fields is ',i5)") &
               n_fields, max_fields
        stop 'initialize_grid'
      end if
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
      call MultiGridInit(max_grids=1,max_fields=n_fields,nborder=1)
      call SimpleGridNew('Rectangular box', n_points, box)
      !
      !  Allocate data fields
      !
      allocate_fields: do field=1,n_fields
        call FieldNew(' ',f_table(field),scratch=.true.,wavefunction=.true.)
      end do allocate_fields
      n_free = n_fields
      !
      call TimerStop('Grid initialization')
    end subroutine initialize_grid
    !
    !  Allocate fields for and load molecular orbitals defining the Fock operator
    !  We'll also grab the nuclear positions and charges, and do a bit of
    !  output to give the user a chance to spot stupid input problems
    !
    subroutine load_gamess_mos
      integer(ik) :: iat, imo
      real(rk)    :: norm
      !
      if (n_mos>max_mos) then
        write (out,"('Input requested ',i4,' MOs, but compile-time limit is ',i4)") &
               n_mos, max_mos
        stop 'load_gamess_mos - too many MOs'
      end if
      !
      mo_table(1:n_mos) = f_table(n_free-n_mos+1:n_free)
      n_free = n_free - n_mos
      if (n_free<0) then
        stop 'load_gamess_mos - out of fields!'
      end if
      !
      call FieldImport('GAMESS',gamess_file,mo_table(1:n_mos),mo_ind(1:n_mos))
      !
      call gamess_report_nuclei(n_nuclei,xyzq)
      !
      write (out,"()")
      write (out,"('            Number of nuclei is ',i5)") n_nuclei
      write (out,"('Number of molecular orbitals is ',i5)") n_mos
      write (out,"('        Total nuclear charge is ',f10.3)") sum(xyzq(4,1:n_nuclei))
      write (out,"('              Electron count is ',f10.3)") sum(mo_occ(1:n_mos))
      write (out,"('  Total charge of the system is ',f10.3)") sum(xyzq(4,1:n_nuclei))-sum(mo_occ(1:n_mos))
      write (out,"()")
      !
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,n_nuclei
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
      !
      !  Enforce orbital normalization
      !
      write (out,"(t2,a3,t7,a7,t16,a10)") 'MO', 'Occup.', 'Input norm'
      write (out,"(t2,a3,t7,a7,t16,a10)") '--', '------', '----------'
      orbital_norm: do imo=1,n_mos
        call QMNormalize(mo_table(imo),1.0_rk,norm)
        write (out,"(t2,i3,t7,f7.3,t16,f10.7)") imo, mo_occ(imo), norm
      end do orbital_norm
      write (out,"()")
    end subroutine load_gamess_mos
    !
    !  Add absorbing potential far downfield
    !
    subroutine add_absorbing_potential
      real(rk)    :: box_max
      !
      v_sink = f_table(n_free) ; n_free = n_free - 1
      if (n_free<0) stop 'add_core_absorbing_potential - out of fields!'
      !
      call TimerStart('Absorbing potential')
      box_max = box_extent(e_axis)/2
      call CAPsetUpPolynomial1D(sink_order,e_axis,sink_en,box_max-sink_width,box_max)
      call FieldInit(v_sink,CAPpolynomial1D)
      call TimerStop('Absorbing potential')
      !
    end subroutine add_absorbing_potential
    !
    !  A wrapper for applying the Fock operator of the input w.f. to a field
    !
    subroutine do_fock(psi,fpsi)
      integer(ik), intent(in) :: psi  ! Input:  Field containing the orbital. The orbital
                                      !         is assumed to be of spin alpha
      integer(ik), intent(in) :: fpsi ! Output: Field containing the result of application 
                                      !         of the Fock operator
      !
      call fock_operator(psi=psi,fpsi=fpsi,v_ext=v_external,v_hartree=v_hartree, &
                         mos=mo_table(:n_mos),f_scr=f_table(:n_free))
    end subroutine do_fock
    !
    !  Calculate total energy and Lagrangian multipliers for the input orbitals
    !  These are needed later, and serve as a good test of correctness.
    !
    subroutine calculate_input_hartree_fock
      integer(ik) :: imo, jmo  ! Orbital indices
      integer(ik) :: jfock
      !
      !  First, the total energy.
      !
      e_input = fock_energy(v_ext=v_external,v_hartree=v_hartree,rho=rho_total, &
                            mos=mo_table(:n_mos),occ=mo_occ(:n_mos),f_scr=f_table(:n_free))
      if (verbose>=1) then
        write(out,"('Hartree-Fock energy of the input w.f.: ',f16.8,' Hartree')") e_input
      end if
      !
      !  Now, the Lagrangian multipliers.
      !
      jfock = f_table(n_free) ; n_free = n_free - 1
      right_orbitals: do jmo=1,n_mos
        call do_fock(mo_table(jmo),jfock)
        left_orbitals: do imo=1,n_mos
          eps_input(imo,jmo) = FieldConjgIntegrate(left=mo_table(imo),right=jfock)
        end do left_orbitals
      end do right_orbitals
      !
      if (verbose>=2) then
        write (out,"(/t10,'Lagrangian multipliers before symmetrization')") 
        call pretty_print
      end if
      !
      eps_input(:n_mos,:n_mos) = 0.5_rk*(eps_input(:n_mos,:n_mos) + &
                               transpose(eps_input(:n_mos,:n_mos)))
      !
      if (verbose>=1) then
        write (out,"(/t10,'Symmetrized Lagrangian multipliers for the input orbitals')") 
        call pretty_print
      end if
      return
      !
      contains
        subroutine pretty_print
          integer(ik) :: cap_1, cap_2
          paginate: do cap_1=1,n_mos,5
            cap_2 = min(cap_1+4,n_mos)
            !
            write (out,"()")
            write (out,"(t20,5(5x,i5,  5x))") (       jmo, jmo=cap_1,cap_2)
            write (out,"(t20,5(8x,f6.3,1x))") (mo_occ(jmo),jmo=cap_1,cap_2)
            !
            print_page: do imo=1,n_mos
              write (out,"(1x,i5,2x,f6.3,2x,t20,5(1x,g14.7))") &
                     imo, mo_occ(imo), (eps_input(imo,jmo),jmo=cap_1,cap_2)
            end do print_page
            write (out,"()")
          end do paginate
        end subroutine pretty_print
    end subroutine calculate_input_hartree_fock
    !
    !  Very simple visualization routine
    !
    subroutine visualize_wavefunction(pass,en,psi)
      integer(ik), intent(in) :: pass
      complex(rk), intent(in) :: en 
      integer(ik), intent(in) :: psi
      !
      character(len=200) :: buf
      !
      write (buf,"('Source rate = ',g12.5,' pass ',i4,' energy ',2g14.7)") &
             src_rate, pass, en
      !
      call FieldVisualize(slot=0,src=psi,text=trim(buf))
    end subroutine visualize_wavefunction
    !
    !  Problem driver
    !
    subroutine run_tunnel
      integer(ik) :: info
      real(rk)    :: efield(3)
      !
      call TimerStart('Tunnel')
      call accuracyInitialize
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=tunneldata,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=tunneldata)
      write (out,"()")
      !
      call initialize_grid
      !
      call load_gamess_mos
      !
      call fock_set_options(sor_rate=sor_rate,eps=eps_hartree)
      !
      !  Construct multiplicative part of the potential in the Fock operator.
      !
      efield = 0
      efield(e_axis) = e_field
      !
      if (n_free<3) stop 'run_tunnel - out of data fields'
      v_external = f_table(n_free) ; n_free = n_free - 1
      v_hartree  = f_table(n_free) ; n_free = n_free - 1
      rho_total  = f_table(n_free) ; n_free = n_free - 1
      !
      !  External potential in the Fock operator
      !
      call fock_external_potential(v=v_external,   &
                       f_scr=f_table(:n_free),     &
                       xyzq=xyzq(:,1:n_nuclei),    &
                       efield=efield)
      !
      !  Hartree potential in the Fock operator
      !
      call fock_hartree_potential(v=v_hartree,     &
                       rho=rho_total,              &
                       mos=mo_table(1:n_mos),occ=mo_occ(1:n_mos))
      !
      !  Prepare imaginary absorbing potential
      !
      call add_absorbing_potential
      !
      !  Calculate Hartree-Fock total energy and Lagrangian multipliers
      !  corresponding to the orbitals supplied on input.
      !
      call calculate_input_hartree_fock
      !

      !
      !
      write (out,"(/' Done with preliminaries'/)")
      call TimerReport
!     !
!     !  We are ready to start iterations
!     !
!     psi  = f_table(n_free) ; n_free = n_free - 1
!     if (n_free<0) stop 'run_tunnel - out of data fields'
!     !
!     !  Initial guess - the field-perturbed MO
!     !
!     call FieldCopy(src=mo_table(n_mos),dst=psi)
!     !
!     !  Final test - calculate at report the energy expectation value
!     !
!     call starting_guess_energy
!     !
!     !  Get stationary solution in the presence of the source and the sink
!     !
!     call solve_stationary_problem
!     !
!     !  Eeek.
!     !
      call TimerStop('Tunnel')
      call TimerReport
    end subroutine run_tunnel

  end module tunnel
!
  subroutine driver
    use tunnel

    call run_tunnel
  end subroutine driver
