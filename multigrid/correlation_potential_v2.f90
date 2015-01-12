!
!  Evaluation of "correlation pulse" potential. See: Z.B. Walters and O. Smirnova,
!  J Phys B 43, 161002 (2010) for the idea behind this stuff.
!
!  This is the second version of the code, using analytical evaluation of
!  all molecular integrals. As the result, it should be able to handle complex
!  trajectories.
!
  module correlation_potential_v2
    use accuracy
    use timer
    use import_gamess
    !$ use OMP_LIB

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)           :: verbose    = 0          ! Verbosity level
    integer(ik)           :: max_rdm_sv = 1000       ! Should be good enough ...
    integer(ik)           :: max_points = 100000     ! Should be good enough ...
    character(len=100)    :: rdm_file   = ' '        ! Name of the file containing rdm orbital coefficients
    character(len=100)    :: algorithm  = 'AO'       ! Either 'MO' or 'AO'. 'AO' should always be faster,
                                                     ! so there should be no reason to use 'MO'
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)           :: rdm_count               ! Number of singular values
    real(rk), allocatable :: rdm_sv(:)               ! Singular values of the transition 1-RDM
    real(rk), allocatable :: rdm_ao(:,:)             ! Reduced density matrix in the AO basis
    type(gam_structure)   :: mol                     ! Gamess basis set and orbital data
    integer(ik)           :: npoints                 ! Actual number of points
    real(rk), allocatable :: points(:,:)             ! Buffer for data points; makes parallel programming a bit simpler
    real(rk), allocatable :: vc(:,:)                 ! Coulomb potential at the nearest real point and at the actual point
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /correlation/ verbose, rdm_file, max_rdm_sv, max_points, algorithm
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Report molecular geometry - lifted from dyson_tools.f90
    !
    subroutine report_structure(level,name)
      integer(ik), intent(in)      :: level ! Minimal verbosity level to print the structure;
                                            ! Report just the number of atoms otherwise
      character(len=*), intent(in) :: name  ! Name of the structure/file
      !
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call gamess_report_nuclei(nnuc,structure=mol)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      if (level>verbose) return
      !
      !  Tell a bit more!
      !
      allocate (xyzq(4,nnuc))
      call gamess_report_nuclei(nnuc,xyzq,structure=mol)
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
    subroutine load_molecular_data
      real(rk)    :: tmp_rdm_sv(max_rdm_sv)
      !
      call TimerStart('Load molecular data')
      write (out,"(/'Loading ',a)") trim(rdm_file)
      call gamess_load_rdmsv(trim(rdm_file),tmp_rdm_sv,rdm_count)
      write (out,"( 'Found ',i4,' singular values')") rdm_count
      write (out,"( 'Values are: '/)")
      write (out,"(10(1x,f12.8))") tmp_rdm_sv(:rdm_count)
      !
      allocate (rdm_sv(rdm_count))
      rdm_sv = tmp_rdm_sv(:rdm_count)
      !
      call gamess_load_orbitals(file=trim(rdm_file),structure=mol)
      !
      call report_structure(0,trim(rdm_file))
      call TimerStop('Load molecular data')
    end subroutine load_molecular_data
    !
    !  Transform 1-RDM to the AO basis
    !
    subroutine transform_1rdm
      integer(ik) :: nao, nmo, ird, imo, mu, nu
      real(rk)    :: tmp
      !
      call TimerStart('1-RDM: MO->AO')
      !
      nao = mol%nbasis
      nmo = mol%nvectors
      allocate (rdm_ao(nao,nao))
      !
      !$omp parallel do default(none) private(nu,mu,imo,tmp) &
      !$omp& shared(nao,nmo,mol,rdm_ao,rdm_count,rdm_sv)
      rdm_ao_nu: do nu=1,nao
        rdm_ao_mu: do mu=1,nao
          tmp = 0._rk
          rdm_ao_sv: do ird=1,rdm_count
            imo = 2*ird-1
            tmp = tmp + rdm_sv(ird) * mol%vectors(mu,imo) * mol%vectors(nu,imo+1)
          end do rdm_ao_sv
          rdm_ao(mu,nu) = tmp
        end do rdm_ao_mu
      end do rdm_ao_nu
      !$omp end parallel do
      !
      call TimerStop('1-RDM: MO->AO')
    end subroutine transform_1rdm
    !
    function evaluate_property(what,coord) result(v)
      character(len=*), intent(in) :: what     ! Operator to evaluate
      real(rk), intent(in)         :: coord(:) ! X,Y,Z,A of the operator
      real(rk)                     :: v        ! Desired integral
      !
      integer(ik)           :: nao, nmo, ird, imo, nu
      real(rk), allocatable :: ao_ints(:,:)
      real(rk), allocatable :: mo_ints(:,:)
      !
      call TimerStart('Property '//trim(what))
      nao = mol%nbasis
      nmo = mol%nvectors
      allocate (ao_ints(nao,nao))
      call gamess_1e_integrals('AO 3C '//what,ao_ints,bra=mol,ket=mol,op_xyz=coord(1:3),op_param=coord(4:4))
      if (algorithm=='AO') then
        v = 0
        sum_property_ao: do nu=1,nao
          v = v + dot_product(ao_ints(:,nu),rdm_ao(:,nu))
        end do sum_property_ao
      else
        allocate (mo_ints(nmo,nmo))
        mo_ints = matmul(transpose(mol%vectors),matmul(ao_ints,mol%vectors))
        v = 0
        sum_property_mo: do ird=1,rdm_count
          imo = 2*ird-1
          v = v + rdm_sv(ird) * mo_ints(imo,imo+1)
        end do sum_property_mo
        deallocate (mo_ints)
      end if
      deallocate (ao_ints)
      call TimerStop('Property '//trim(what))
    end function evaluate_property
    !
    subroutine start
      integer(ik) :: info, ipt
      !
      call TimerStart('start')
      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=correlation,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=correlation)
      write (out,"()")
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call load_molecular_data
      !
      if (algorithm=='AO') call transform_1rdm
      !
      call TimerStart('Load trajectory')
      allocate (points(4,max_points),vc(3,max_points))
      npoints = 0
      load_trajectory_loop: do while(npoints+1<=max_points)
        read(input,*,iostat=info) points(:,npoints+1)
        if (info/=0) exit load_trajectory_loop
        npoints = npoints + 1
      end do load_trajectory_loop
      call TimerStop('Load trajectory')
      !
      write (out,"('Loaded ',i8,' trajectory points.')") npoints
      if (npoints+1>=max_points) then
        write (out,"(/'WARNING: Reached the end of buffer loading trajectory. Some points may be left behind.'/)")
      end if
      !
      call TimerStart('Run trajectory')
      !$omp parallel default(none) private(ipt) shared(vc,points,npoints)
      !$omp do schedule(guided)
      run_trajectory: do ipt=1,npoints
        if (mod(ipt,100)==0) write (out,"('Processing trajectory point ',i8)") ipt
        vc(1,ipt) = evaluate_property('1/R',          points(:,ipt))
        !
        !  Zero imaginary coordinate requires special handling
        !
        if (abs(points(4,ipt))<=spacing(1._rk)) then
          vc(2,ipt) = vc(1,ipt)
          vc(3,ipt) = 0._rk
        else
          vc(2,ipt) = evaluate_property('R/(R**2+A**2)',points(:,ipt))
          vc(3,ipt) = evaluate_property('A/(R**2+A**2)',points(:,ipt))
        end if
      end do run_trajectory
      !$omp end do
      !$omp end parallel
      call TimerStop('Run trajectory')
      !
      call TimerStart('Print trajectory')
      write (out,"(('#',4(1x,a12),2x,a18,1x,2(1x,a18)))") &
         '   X   ', '   Y   ', '   Z   ', '   A   ', '   V(1/R)   ', '  V(R/(R^2+A^2)) ', ' -V(A/(R^2+A^2) ', &
         '-------', '-------', '-------', '-------', '------------', '-----------------', '----------------'
      print_trajectory: do ipt=1,npoints
        !
        !  Note that the imaginary part of the potential has to change sign.
        !
        write (out,"('@',4(1x,f12.6),2x,g18.11,1x,2(1x,g18.11))") points(:,ipt), vc(:,ipt)
      end do print_trajectory
      call TimerStop('Print trajectory')
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
    
  end module correlation_potential_v2
  !
  !
  !
  subroutine driver
    use accuracy
    use math
    use correlation_potential_v2
    !
    real(rk) :: dummy
    !
    call accuracyInitialize
    dummy = MathDoubleFactorial(100_ik)
    dummy = MathFactorial(100_ik)
    dummy = MathLogDoubleFactorial(1000_ik)
    dummy = MathLogFactorial(1000_ik)

    call start

  end subroutine driver
