!
!  Evaluation of a quadratic forcefield for a fixed transition 1-RDM,
!  corresponding to a "crude adiabatic" approximation for the electronic
!  structure.
!
!  WARNING: "Hessians" calculates here are just plain wrong; gradients
!  WARNING: are only as good as the basis quality in the core region;
!  WARNING: the reasons are the same as to why Hellman-Feynman forces
!  WARNING: tend to be pretty bad.
!
!  This program is a companion to franck-condon/multisurface_ac.f90
!
!  2014 Jul 04: Intital version, derived from correlation_potential_v2.f90
!
  module crude_adiabatic_ff
    use accuracy
    use timer
    use import_gamess
    use atoms
    !$ use OMP_LIB
    implicit none
    private
    public start
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: iu_scr     = 44_ik      ! I/O unit for all our writing and temp files
    integer(ik), parameter :: name_len   = 128        ! Max. length of a file name
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)            :: verbose    = 0          ! Verbosity level
    character(len=130)     :: comment    = ' '        ! Optional comment, to be added to all output files
    integer(ik)            :: max_rdm_sv = 1000       ! Should be good enough ...
    character(len=name_len):: rdm_file   = ' '        ! Name of the file containing rdm orbital coefficients
    real(rk)               :: step       = 0.0001_rk  ! Finite-differentiation step, in Bohr
    character(len=name_len):: xyz_file   = 'xyz.dat'  ! Name of the file for the structure, .xyz format
    character(len=name_len):: dip_file   = 'dip.dat'  ! Name of the file for the transition dipole; 
                                                      ! blank suppresses writing (but not evaluation) of the dipole
    character(len=name_len):: grad_file  = 'grad.dat' ! Name of the file for the gradient
                                                      ! blank suppresses writing (but not evaluation) of the gradient
    character(len=name_len):: hess_file  = ' '        ! Name of the file for the gradient
                                                      ! blank suppresses writing (and calculation) of hessian
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    integer(ik)            :: rdm_count               ! Number of singular values
    real(rk), allocatable  :: rdm_sv(:)               ! Singular values of the transition 1-RDM
    real(rk), allocatable  :: rdm_ao(:,:)             ! Reduced density matrix in the AO basis
    real(rk)               :: rdm_trace               ! Trace of the density matrix; determines whether
                                                      ! we need to include nuclear contributions or not
    logical                :: add_nuclear             ! True if the nuclear term is called for
    type(gam_structure)    :: mol                     ! GAMESS basis set and orbital data
    integer(ik)            :: nnuc                    ! Number of nuclei in the structure
    real(rk), allocatable  :: xyzq(:,:)               ! Coordinates and charges of the nuclei (Bohr)
    real(rk), allocatable  :: grad(:)                 ! Gradient of the "energy"
    real(rk), allocatable  :: hess(:,:)               ! Hessian of the same
    real(rk)               :: dipole(3)               ! Transition dipole in a.u., including electron charge
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /cra_ff/ verbose,                          &
                      comment,                          &
                      rdm_file, max_rdm_sv, step,       &
                      xyz_file,                         &
                      dip_file, grad_file, hess_file
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
      integer(ik) :: iat
      !
      call gamess_report_nuclei(nnuc,structure=mol)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      !
      if (allocated(xyzq)) deallocate (xyzq)
      allocate (xyzq(4,nnuc))
      call gamess_report_nuclei(nnuc,xyzq,structure=mol)
      if (level>verbose) return
      !
      !  Tell a bit more!
      !
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,nnuc
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
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
    subroutine add_nuclear_dipole
      integer(ik) :: iat
      real(rk)    :: dip(3)
      !
      dip = 0
      nuc_dip: do iat=1,nnuc
        dip = dip + xyzq(4,iat)*xyzq(1:3,iat)
      end do nuc_dip
      dipole = dipole + dip
    end subroutine add_nuclear_dipole
    !
    !  Transform 1-RDM to the AO basis
    !  Also calculate transition dipole between the two components
    !
    subroutine transform_1rdm
      integer(ik)           :: nao, nmo, ird, imo, mu, nu
      real(rk)              :: tmp
      real(rk), allocatable :: tmp_ao(:,:) ! AO integrals
      !
      call TimerStart('1-RDM: MO->AO')
      !
      nao = mol%nbasis
      nmo = mol%nvectors
      if (allocated(rdm_ao) ) deallocate (rdm_ao)
      allocate (rdm_ao(nao,nao),tmp_ao(nao,nao))
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
      call TimerStop('1-RDM: MO->AO')
      call TimerStart('1-RDM trace&overlap')
      !
      !  It is much easier to deal with the 1-RDM trace in the AO basis
      !  (our represenation of 1-RDM uses SVD, so we don't know the
      !  eigenvalues per se). However, we'll need to have overlap
      !  matrix before we can proceed.
      !
      call gamess_1e_integrals('AO OVERLAP',tmp_ao,bra=mol,ket=mol)
      rdm_trace   = sum(rdm_ao*tmp_ao)
      add_nuclear = abs(rdm_trace)>1e-5_rk
      !  
      call gamess_1e_integrals('AO DIPOLE X',tmp_ao,bra=mol,ket=mol)
      dipole(1) = -sum(rdm_ao*tmp_ao)
      call gamess_1e_integrals('AO DIPOLE Y',tmp_ao,bra=mol,ket=mol)
      dipole(2) = -sum(rdm_ao*tmp_ao)
      call gamess_1e_integrals('AO DIPOLE Z',tmp_ao,bra=mol,ket=mol)
      dipole(3) = -sum(rdm_ao*tmp_ao)
      !
      if (add_nuclear) call add_nuclear_dipole
      !
      write (out,"(/'                          1-RDM trace is: ',f17.12)") rdm_trace
      write (out,"( '                   Include nuclear terms: ',l4)") add_nuclear
      write (out,"( 'Transition dipole is [au, incl e charge]: ',3(f12.7,1x)/)") dipole
      !
      call TimerStop('1-RDM trace&overlap')
    end subroutine transform_1rdm
    !
    function electron_nuclear(xyzq) result(v)
      real(rk), intent(in) :: xyzq(:,:) ! X,Y,Z,Q of the nuclei
      real(rk)             :: v         ! e-N interaction energy, Hartree
      !
      integer(ik)           :: nao, nmo, iat
      real(rk)              :: v_nuc
      real(rk), allocatable :: ao_ints(:,:)
      !
      call TimerStart('e-n interaction')
      nao = mol%nbasis
      nmo = mol%nvectors
      allocate (ao_ints(nao,nao))
      v = 0
      atom_loop: do iat=1,nnuc
        call gamess_1e_integrals('AO 3C 1/R',ao_ints,bra=mol,ket=mol,op_xyz=xyzq(1:3,iat))
        v_nuc = sum(ao_ints*rdm_ao)
        v     = v - xyzq(4,iat)*v_nuc
      end do atom_loop
      deallocate (ao_ints)
      call TimerStop('e-n interaction')
    end function electron_nuclear
    !
    function nuclear_repulsion(xyzq) result(v)
      real(rk), intent(in) :: xyzq(:,:) ! X,Y,Z,Q of the nuclei
      real(rk)             :: v         ! e-N interaction energy, Hartree
      !
      integer(ik) :: iat, jat
      real(rk)    :: r
      !
      v = 0
      atom_i: do iat=1,nnuc-1
        atom_j: do jat=iat+1,nnuc
          r = sqrt(sum((xyzq(1:3,iat)-xyzq(1:3,jat))**2))
          if (r <= 1e-2_rk) stop 'nuclear_repulsion - nuclear collision!'
          v = v + xyzq(4,iat)*xyzq(4,jat)/r
        end do atom_j
      end do atom_i
    end function nuclear_repulsion
    !
    function total_energy(ic,ic_step) result(v)
      integer(ik), intent(in) :: ic     (:) ! Coordinate(s) to perturb; 0 if no change
      real(rk), intent(in)    :: ic_step(:) ! Desired change (Bohr)
      real(rk)                :: v          ! e-N interaction energy, Hartree
      !
      integer(ik) :: ip, icv, iat, ixyz
      real(rk)    :: temp_xyzq(4,nnuc) ! X,Y,Z,Q of the nuclei after displacement
      !
      if (size(ic)/=size(ic_step)) stop 'total_energy - displacement oops'
      !
      temp_xyzq = xyzq
      displacements: do ip=1,size(ic)
        icv = ic(ip)
        if (icv<0 .or. icv>3*nnuc) stop 'total_energy - range oops'
        if (icv==0) cycle displacements
        iat  = (icv+2)/3
        ixyz = icv - 3*(iat-1)
        ! write (out,"('Displacing atom ',i0,' index ',i0,' by ',f12.6)") iat, ixyz, ic_step(ip)
        temp_xyzq(ixyz,iat) = temp_xyzq(ixyz,iat) + ic_step(ip)
      end do displacements
      !
      v = electron_nuclear(temp_xyzq) 
      if (add_nuclear) v = v + nuclear_repulsion(temp_xyzq)
    end function total_energy
    !
    subroutine energy_test
      real(rk) :: e_en, e_nn
      !
      e_en = electron_nuclear (xyzq)
      e_nn = nuclear_repulsion(xyzq)
      write (out,"(/' Electron-nuclear interaction energy: ',g25.15)") e_en
      write (out,"( '            Nuclear repulsion energy: ',g25.15/)") e_nn
    end subroutine energy_test
    !
    !  Trying to construct an optimal displacement list is an overkill
    !  for our problem; however, I find any other solution to be
    !  aesthetically grating, so here it goes ...
    !
    subroutine calculate_hessian
      integer(ik)              :: ic, jc, ij_max, ij2_max
      integer(ik)              :: disp_max, ndisp, idisp
      integer(ik), allocatable :: disp_directory(:,:,:,:)
      integer(ik), allocatable :: disp_list_ic  (:,:)
      real(rk), allocatable    :: disp_list_step(:,:)
      real(rk), allocatable    :: disp_energy   (  :)
      real(rk)                 :: fp, fm, fpp, fpm, fmp, fmm
      !
      ij_max   = 3*nnuc
      ij2_max  = ij_max
      if (hess_file==' ') ij2_max = 0  ! Hessian is not requested; do not do these displacements
      disp_max = (2*(1+ij_max))**2
      !
      allocate (disp_directory(2,0:ij_max,2,0:ij_max))
      allocate (disp_list_ic(2,disp_max),disp_list_step(2,disp_max),disp_energy(disp_max))
      disp_directory = -1
      !
      call TimerStart('Displacement directory')
      ndisp = 0
      coordinate_ic: do ic=0,ij_max
        coordinate_jc: do jc=0,ij2_max
          disp_directory(1,ic,1,jc) = add_displacement(-step,-step)
          disp_directory(1,ic,2,jc) = add_displacement(-step, step)
          disp_directory(2,ic,1,jc) = add_displacement( step,-step)
          disp_directory(2,ic,2,jc) = add_displacement( step, step)
        end do coordinate_jc
      end do coordinate_ic
      write (out,"('Found ',i0,' unique displacements.')") ndisp
      call TimerStop('Displacement directory')
      !
      call TimerStart('Evaluate displacements')
      !$omp parallel do default(none) shared(ndisp,disp_energy,disp_list_ic,disp_list_step) private(idisp)
      evaluate_displacements: do idisp=1,ndisp
        disp_energy(idisp) = total_energy(disp_list_ic(:,idisp),disp_list_step(:,idisp))
      end do evaluate_displacements
      !$omp end parallel do
      call TimerStop('Evaluate displacements')
      !
      !  Construct the gradient and the Hessian using numerical differentiation
      !
      call TimerStart('Numerical diffentiation')
      if (allocated(grad)) deallocate (grad)
      if (allocated(hess)) deallocate (hess)
      allocate (grad(ij_max),hess(ij_max,ij_max))
      gradient_loop: do ic=1,ij_max
        fp = get_energy(2,ic,1,0) 
        fm = get_energy(1,ic,1,0) 
        grad(ic) = (fp-fm)/(2*step)
      end do gradient_loop
      if (hess_file/=' ') then
        hess_loop_ic: do ic=1,ij_max
          hess_loop_jc: do jc=1,ij_max
            fmm = get_energy(1,ic,1,jc)
            fmp = get_energy(1,ic,2,jc)
            fpm = get_energy(2,ic,1,jc)
            fpp = get_energy(2,ic,2,jc)
            !
            !  This expression also works for the diagonal term in the Hessian,
            !  since fpm = fmp = f00 in this case, and the effective step 
            !  is doubled
            !
              write (out,"(' ic, jc = ',i5,1x,i5)") ic, jc
              write (out,"('    fmm = ',f25.15)") fmm
              write (out,"('    fmp = ',f25.15)") fmp
	      write (out,"('    fpm = ',f25.15)") fpm
              write (out,"('    fpp = ',f25.15)") fpp
            hess(ic,jc) = (fpp + fmm - fmp - fpm)/(4*(step**2))
          end do hess_loop_jc
        end do hess_loop_ic
      end if
      call TimerStop('Numerical diffentiation')
      !
      deallocate(disp_directory,disp_list_ic,disp_list_step,disp_energy)
      !
      contains
      function get_energy(id1,ic1,id2,ic2) result (en)
        integer(ik) :: id1, ic1, id2, ic2 ! Coordinates and displacements
        real(rk)    :: en
        integer(ik) :: pos
        !
        if (id1<1 .or. id1>2 .or. ic1<0 .or. ic1>ij_max .or. &
            id2<1 .or. id2>2 .or. ic2<0 .or. ic2>ij_max ) stop 'calculate_hessian%get_energy - bad indices'
        pos = disp_directory(id1,ic1,id2,ic2)
        if (pos<1 .or. pos>ndisp) stop 'calculate_hessian%get_energy - unexpected displacement'
        en = disp_energy(pos)
      end function get_energy
      function add_displacement(dic,djc) result(pos)
        real(rk), intent(in) :: dic, djc ! Displacements
        integer(ik)          :: pos
        !
        integer(ik) :: can_ic(2)  ! "Canonical" displacement order; displacements commute!
        real(rk)    :: can_dc(2)
        !
        if (ic==jc) then
          can_ic = (/ic,      0_ik /)
          can_dc = (/dic+djc, 0._rk/)
        else if (ic>=jc) then
          can_ic = (/ic, jc /)
          can_dc = (/dic,djc/)
        else
          can_ic = (/jc, ic /)
          can_dc = (/djc,dic/)
        endif
        where (can_ic==0) can_dc = 0._rk
        !
        pos = 1 ! Needed in case ndisp was zero, and the loop did not execute
                ! This is a well-defined case in Fortran-2005, but Intel compiler
                ! (incorrectly) leaves loop counter unchanged for zero-trip DO loops.
        scan_duplicates: do pos=1,ndisp
          if (any(    disp_list_ic  (:,pos)/=can_ic)) cycle scan_duplicates
          if (any(abs(disp_list_step(:,pos)- can_dc)>spacing(100._rk))) cycle scan_duplicates
          exit scan_duplicates
        end do scan_duplicates
        if (pos>ndisp) then
          ndisp = max(ndisp,pos)
          if (ndisp>disp_max) stop 'calculate_hessian%add_displacement - blown table'
          disp_list_ic  (:,ndisp) = can_ic
          disp_list_step(:,ndisp) = can_dc
        end if
      end function add_displacement
    end subroutine calculate_hessian
    !
    subroutine report_xyz
      integer(ik) :: ios, iat
      !
      open (iu_scr,file=trim(xyz_file),form='formatted',action='write',status='new',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new file ',a,' for the .xyz coordinate')") ios, trim(xyz_file)
        stop 'report_xyz - open failed'
      end if
      write (iu_scr,"(i0/a)") nnuc, trim(comment)
      print_atoms: do iat=1,nnuc
        write (iu_scr,"(a2,3(1x,f16.8))") AtomElementSymbol(xyzq(4,iat)), xyzq(1:3,iat)*abohr
      end do print_atoms
      !
      close (iu_scr)
    end subroutine report_xyz
    !
    subroutine report_dipole
      integer(ik) :: ios
      integer(ik) :: irow, icol, ichk, icont, ic2
      !
      open (iu_scr,file=trim(dip_file),form='formatted',action='write',status='new',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new file ',a,' for the dipole')") ios, trim(dip_file)
        stop 'report_dipole - open failed'
      end if
      write (iu_scr,"('# ',a)") trim(comment)
      write (iu_scr,"('# dipole, atomic units, incl. electron charge')")
      write (iu_scr,"('#')")
      !
      write (iu_scr,"(3(1x,e25.15))") dipole
      !
      close (iu_scr)
    end subroutine report_dipole
    !
    subroutine report_gradient
      integer(ik) :: ic
      integer(ik) :: ios
      !
      open (iu_scr,file=trim(grad_file),form='formatted',action='write',status='new',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new file ',a,' for the gradient')") ios, trim(grad_file)
        stop 'report_gradient - open failed'
      end if
      write (iu_scr,"('# ',a)") trim(comment)
      write (iu_scr,"('# Gradients, atomic units')")
      write (iu_scr,"('#')")
      !
      print_gradient: do ic=1,3*nnuc,3
        write (iu_scr,"(3(1x,f25.15))") grad(ic:ic+2)
      end do print_gradient
      !
      close (iu_scr)
    end subroutine report_gradient
    !
    subroutine report_hessian
      integer(ik) :: ios
      integer(ik) :: irow, icol, ichk, icont, ic2
      !
      open (iu_scr,file=trim(hess_file),form='formatted',action='write',status='new',iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new file ',a,' for the Hessian')") ios, trim(hess_file)
        stop 'report_hessian - open failed'
      end if
      write (iu_scr,"('# ',a)") trim(comment)
      write (iu_scr,"('#')")
      write (iu_scr,"('# This hessian is known to be of extremely poor quality')")
      write (iu_scr,"('# Most likely, it will lead to completely unphysical results')")
      write (iu_scr,"('# in the nuclear dynamics simulation. You have been warned')")
      write (iu_scr,"('#$HESS   ')")
      !
      irow = 1 ; icol = 1 ; 
      write_rows: do irow=1,3*nnuc
        ichk = mod(irow,100)
        icont = 1
        write_columns: do icol=1,3*nnuc,5
          ic2  = min(icol+4,3*nnuc)
          write (iu_scr,"(i2,i3,5e15.8)") ichk, icont, hess(irow,icol:ic2)
          icont = icont + 1
        end do write_columns
      end do write_rows
      !
      write (iu_scr,"('#$END    ')")
      close (iu_scr)
    end subroutine report_hessian
    !
    subroutine start
      integer(ik) :: info, ipt
      !
      call TimerStart('start')
      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=cra_ff,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=cra_ff)
      write (out,"()")
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call load_molecular_data
      !
      call transform_1rdm

      !
      call energy_test
      !
      call flush(out)
      call calculate_hessian
      !
      if (xyz_file /=' ') call report_xyz
      !
      if (dip_file /=' ') call report_dipole
      !
      if (grad_file/=' ') call report_gradient
      !
      if (hess_file/=' ') call report_hessian
      !
      call TimerStop('start')
      call TimerReport
      call flush(out)
    end subroutine start
    
  end module crude_adiabatic_ff
  !
  !
  !
  program driver
    use accuracy
    use math
    use crude_adiabatic_ff
    !
    real(rk) :: dummy
    !
    write (out,"('Version: ',a/)") __BUILD_ID__
    !
    call accuracyInitialize
    dummy = MathDoubleFactorial(100_ik)
    dummy = MathFactorial(100_ik)
    dummy = MathLogDoubleFactorial(1000_ik)
    dummy = MathLogFactorial(1000_ik)

    call start

  end program driver
