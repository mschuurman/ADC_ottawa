!
!  Import simple GAMESS checkpoint files.
!
!  Restrictions: 
!    Only C1 symmetry (no symmetry) is allowed.
!    Basis sets must be in-line.
!    Only simple, single-SCF files are supported.
!    Only one $VEC section is supported.
!
!  In most cases, these restrictions mean that the checkpoint file
!  will have to be edited by hand to be loadable by this import filter.
!
!  The module is derived from my (much more capable) gamess-dx.c filter.
!
!  Added Feb 2010 (Michael Spanner): Support for loading natural orbital occupation numbers
!                                    and rdm_sv values from the GAMESS checkpoint files.
!  Added March 2010 (Michael Spanner): Support for rotating the gamess data relative to multigrid 
!
!  Added Aug 2011 (SP): Support for 3-centre 1e integral evaluation. The 3-centre integrals
!                       implementation is essentiall Obara-Saika, with the extensions to
!                       the supported operator classes based on Ahlrichs PCCP 8, 3072 (2006).
!                       Note that OS and A use slightly different definitions of the primitive
!                       integrals, which adds a tiny bit of complexity.
!
!        WARNING: Non-unit rotation matrix is only implemented for
!                 evaluation of orbitals on grid. Obara-Saika integrals
!                 or ECP integrals will prodice incorrect results.
!
!        WARNING: R/(R**2+A**2) and A/(R**2+A**2) integrals have limited range and accuracy.
!                 Please make sure you understand the limits when you use those!
!
!  Added May 2012 (SP): Support for 2e integral evaluation, following Ahlrichs' PCCP.
!                       Support for kinetic-energy integrals, following Obara & Saika.
!
!       WARNING: 2e integral implementation is extremely inefficient.
!
!  Added Dec 21, 2012 (SP): Support for matrix elements of the planewave
!
!  Feb 06, 2013 (SP): Changed batching of the integrals, to improve 2e efficiency
!                     for heavily-contracted basis sets (think correlation-consistent
!                     basis sets).
!
!  Jul 17, 2013 (SP): Added support for r_i (d/d r_j) integrals.
!
!  Apr 23, 2014 (SP): Added support for (most) integrals in quadruple precision.
!
!  In order to remain compatible with historical usage in multigrid, we 
!  maintain a default structure context; all calls which do not specify 
!  (one of the) context(s) default to that one.
!
  module import_gamess
    use accuracy
    use math
    use lapack
    use timer
    use gamess_internal
    use os_integral_operators
    !$ use OMP_LIB
    implicit none
    private
!   public os_basic_integral     ! Debugging!
    public gamess_load_orbitals
    public gamess_evaluate_functions
    public gamess_load_natocc
    public gamess_load_rdmsv
    public gamess_report_nuclei
    public gamess_destroy
    public gamess_1e_integrals
    public gamess_print_1e_integrals
    public gam_atom, gam_structure
    public gamess_oversample
    public gamess_2e_info
    public gamess_2e_integrals
    !
    !  Conditional compilation is needed for the interfaces:
    !  the _quad versions differ only in their argument kinds; if we
    !  compile with real and quad kinds being the same, the interface 
    !  part will fail!
    !
    interface gamess_print_1e_integrals
      module procedure gamess_print_1e_integrals_real
      module procedure gamess_print_1e_integrals_complex
!*qd  module procedure gamess_print_1e_integrals_quad
    end interface gamess_print_1e_integrals
    !
    interface gamess_1e_integrals
      module procedure gamess_1e_integrals_real
      module procedure gamess_1e_integrals_complex
!*qd  module procedure gamess_1e_integrals_quad
    end interface gamess_1e_integrals
    !
    interface gamess_2e_integrals
      module procedure gamess_2e_integrals_real
!*qd  module procedure gamess_2e_integrals_quad
    end interface gamess_2e_integrals
    !
    !  Interfaces below are for the internal routines; they exist to make it possible to
    !  share code between different-precision versions of the routines.
    !
    interface os_1e_matrix
      module procedure os_1e_matrix_real
      module procedure os_1e_matrix_complex
!*qd  module procedure os_1e_matrix_quad
    end interface os_1e_matrix
    !
    interface os_1e_contraction
      module procedure os_1e_contraction_real
      module procedure os_1e_contraction_complex
!*qd  module procedure os_1e_contraction_quad
    end interface os_1e_contraction
    !
    interface os_1e_overlap
      module procedure os_1e_overlap_real
!*qd  module procedure os_1e_overlap_quad
    end interface os_1e_overlap
    !
    interface os_1e_dipole
      module procedure os_1e_dipole_real
!*qd  module procedure os_1e_dipole_quad
    end interface os_1e_dipole
    !
    interface os_1e_grad
      module procedure os_1e_grad_real
!*qd  module procedure os_1e_grad_quad
    end interface os_1e_grad
    !
    interface os_1e_rddr
      module procedure os_1e_rddr_real
!*qd  module procedure os_1e_rddr_quad
    end interface os_1e_rddr
    !
    interface os_1e_kinetic
      module procedure os_1e_kinetic_real
!*qd  module procedure os_1e_kinetic_quad
    end interface os_1e_kinetic
    !
    interface os_1e_3c
      module procedure os_1e_3c_real
!*qd  module procedure os_1e_3c_quad
    end interface os_1e_3c
    !
    interface os_2e_batch
      module procedure os_2e_batch_real
!*qd  module procedure os_2e_batch_quad
    end interface os_2e_batch
    !
    interface os_2e_primitives
      module procedure os_2e_primitives_real
!*qd  module procedure os_2e_primitives_quad
    end interface os_2e_primitives
    !
    interface os_basic_integral
      module procedure os_basic_integral_real
!*qd  module procedure os_basic_integral_quad
    end interface os_basic_integral
    !
    interface os_boys_table
      module procedure os_boys_table_real
!*qd  module procedure os_boys_table_quad
    end interface os_boys_table
    !
    interface os_common_primitives
      module procedure os_common_primitives_real
!*qd  module procedure os_common_primitives_quad
    end interface os_common_primitives
    !
    !  Internal data structures and magic constants
    !
    integer(ik), parameter            :: verbose    =  0      ! 0 is too verbose for reading files repeatedly
    !
    logical, save                     :: oversample = .true.  ! Activate oversampling
    !
    type(gam_structure), save, target :: gam_def              ! Default GAMESS structure description
    !
    contains
    !
    !  External interfaces
    !
    subroutine gamess_oversample(sample)
      logical, intent(in) :: sample
      !
      oversample = sample
    end subroutine gamess_oversample
    !
    !  Report positions and charges of the nuclei loaded alongside the MOs
    !
    subroutine gamess_report_nuclei(natom,xyzq,structure,true_charge)
      integer(ik), intent(out)        :: natom       ! Number of atoms in the molecule loaded by the
                                                     ! last call to gamess_load_orbitals, or 0 if none
      real(rk), intent(out), optional :: xyzq(:,:)   ! Coordinates and charges of the atoms
      type(gam_structure), intent(in), &
                     target, optional :: structure   ! Context to use
      logical, intent(in), optional   :: true_charge ! Report true nuclear charge, including the ECP charge
      !
      type(gam_structure), pointer :: gam          
      integer(ik)                  :: na
      logical                      :: add_ecp
      !
      if (present(structure)) then
        gam => structure
      else
        gam => gam_def
      end if
      natom = gam%natoms
      !
      add_ecp = .false.
      if (present(true_charge)) add_ecp = true_charge
      !
      if (present(xyzq)) then
        if ( size(xyzq,dim=1)<4 .or. size(xyzq,dim=2)<gam%natoms ) then
          stop 'import_gamess%gamess_report_nuclei - buffer too small'
        end if
        xyzq(1,1:gam%natoms) = real(gam%atoms(1:gam%natoms)%xyz(1),kind=kind(xyzq)) / abohr
        xyzq(2,1:gam%natoms) = real(gam%atoms(1:gam%natoms)%xyz(2),kind=kind(xyzq)) / abohr
        xyzq(3,1:gam%natoms) = real(gam%atoms(1:gam%natoms)%xyz(3),kind=kind(xyzq)) / abohr
        xyzq(4,1:gam%natoms) = real(gam%atoms(1:gam%natoms)%znuc,kind=kind(xyzq))
        if (add_ecp) then
          xyzq(4,1:gam%natoms) = xyzq(4,1:gam%natoms) + nint(gam%atoms(1:gam%natoms)%ecp_zcore)
        end if
        !
        !  Rotate coordinates
        !
        rotate_atoms: do na=1,gam%natoms
          xyzq(1:3,na) = real(matmul(gam%rotmat,xyzq(1:3,na)),kind=kind(xyzq))
        end do rotate_atoms
      end if
    end subroutine gamess_report_nuclei
    !
    !  Destroy data context
    !
    subroutine gamess_destroy(structure)
      type(gam_structure), intent(inout), target, optional :: structure ! Context to use
      !
      type(gam_structure), pointer :: gam
      integer(ik)                  :: alloc
      !
      if (present(structure)) then
        gam => structure
      else
        gam => gam_def
      end if
      if (allocated(gam%atoms)) then
        deallocate(gam%atoms,stat=alloc)
        if (alloc/=0) then
          write (out,"('gamess_import%gamess_destroy: Error ',i6,' deallocating atoms table')") alloc
          stop 'gamess_import%gamess_destroy: Error deallocating atoms table'
        end if
      end if
      if (allocated(gam%bas_labels)) then
        deallocate(gam%bas_labels,stat=alloc)
        if (alloc/=0) then
          write (out,"('gamess_import%gamess_destroy: Error ',i6,' deallocating basis labels table')") alloc
          stop 'gamess_import%gamess_destroy: Error deallocating basis labels table'
        end if
      end if
      if (allocated(gam%vectors)) then
        deallocate(gam%vectors,stat=alloc)
        if (alloc/=0) then
          write (out,"('gamess_import%gamess_destroy: Error ',i6,' deallocating vectors table')") alloc
          stop 'gamess_import%gamess_destroy: Error deallocating vectors table'
        end if
      end if
      if (allocated(gam%batches)) then
        deallocate(gam%batches,stat=alloc)
        if (alloc/=0) then
          write (out,"('gamess_import%gamess_destroy: Error ',i6,' deallocating batches table')") alloc
          stop 'gamess_import%gamess_destroy: Error deallocating batches table'
        end if
      end if
      gam%natoms   = 0
      gam%nbasis   = 0
      gam%nvectors = 0
      gam%nbatches = 0
    end subroutine gamess_destroy
    !
    !  Load structure file and evaluate the MOs. This subroutine is a little 
    !  complicated, since it can perform three separate functions:
    !
    !  1. Parse GAMESS checkpoint file, and remember the result. In this case,
    !     it should be called with the argument "file" (and possibly "structure")
    !  2. Evaluate orbitals on grid. In this case, the arguments "mos", "dst",
    !     "coord", and "grid" (and possibly "structure") should all be present,
    !     but "file" should not be.
    !  3. Perform (1), then (2), then deallocate the orbitals. In this case, all
    !     arguments (except possibly "structure") should be present.
    !
    !  Historically, only the #3 semantics was used in multigrid.
    !
    !  See also gamess_evaluate_functions() for a less integrated interface
    !
    subroutine gamess_load_orbitals(file,mos,dst,coord,grid,structure,rot,dx)
      character(len=*), intent(in), optional :: file           ! Name of the GAMESS checkpoint file. WARNING: Presence of the
                                                               ! file= argument will affect semantics of the rot= argument
                                                               ! below
      integer(ik), intent(in), optional      :: mos(:)         ! Indices of the orbitals to be computed
      integer(ik), intent(in), optional      :: dst(:)         ! Last indices of the grid array, receiving
                                                               ! the MOs.
      real(rk), intent(in), optional         :: coord(:,:,:,:) ! Coordinates of the points, at which MOs
                                                               ! must be evaluated. First index is X,Y,Z.
      complex(rk), intent(out), optional     :: grid(:,:,:,:)  ! Data field for the MOs. First three
                                                               ! indices match last three indices of coord().
                                                               ! The MO mos(i) should end up in grid(:,:,:,dst(i))
      type(gam_structure), intent(inout), target, optional &
                                             :: structure      ! Context to use
      real(rk), intent(in), optional        :: rot(:,:)        ! Rotation matrix to apply to the structure. It is
                                                               ! perfectly OK to use different rotation matrix
                                                               ! in different calls using the same context. 
                                                               ! If the rotation matrix is not given for a call where
                                                               ! file= is specified, it will be reset to unit matrix.
                                                               ! If the rotation matrix is not given for a pure-evaluation
                                                               ! call (file= not specified), rotation matrix from the
                                                               ! previous call will be used.
      real(rk), intent(in), optional         :: dx(:)          ! Grid spacing, used for oversampling. Supplying all-negative
                                                               ! input, or omitting this argument, suppresses oversampling.
      !
      type(gam_structure), pointer :: gam
      integer(ik)                  :: ios
      logical                      :: have_geo, have_vec, have_ecp
      logical                      :: load_mos, evaluate_mos
      logical                      :: eof
      real(xrk)                    :: det_rot
      !
      if (present(structure)) then
        gam => structure
      else
        gam => gam_def
      end if
      !
      !  Decide on the call semantics
      !
      load_mos     = present(file)
      evaluate_mos = present(mos)
      !
      !  Deal with the rotation matrix
      !
      if (present(rot)) then
        gam%rotmat    = real(rot,kind=kind(gam%rotmat))
        gam%rotmat_rk = real(gam%rotmat,kind=kind(gam%rotmat_rk))
      else if (load_mos) then
        gam%rotmat      = 0.0_xrk
        gam%rotmat(1,1) = 1.0_xrk
        gam%rotmat(2,2) = 1.0_xrk
        gam%rotmat(3,3) = 1.0_xrk
        gam%rotmat_rk = real(gam%rotmat,kind=kind(gam%rotmat_rk))
      end if
      !
      det_rot = lapack_determinant(gam%rotmat)
      if (abs(abs(det_rot)-1)>spacing(100._rk)) then
         write (out,"('Rotation matrix has non-unit determinant: ',g20.13)") det_rot
         stop 'import_gamess%gamess_load_orbitals: Bad rotation matrix'
      end if
      !
      !  Grid spacing; if not supplied upon orbital load, disable oversampling.
      !
      if (present(dx)) then
        gam%dx = dx
      else if (load_mos) then
        gam%dx = -1._rk
      end if
      !
      !  Make sure present arguments are consistent with the semantics
      !
      if (evaluate_mos) then
        if (.not.present(dst) .or. .not.present(coord) .or. .not.present(grid)) then
          write (out,"('import_gamess%gamess_load_orbitals: Arguments needed for orbital evaluation are missing')")
          stop 'import_gamess%gamess_load_orbitals: Arguments needed for orbital evaluation are missing'
        end if
        !
        !  A little bit of sanity checking on dimensions
        !
        if ( (size(mos)/=size(dst))                .or. &
             (any(dst<1))                          .or. &
             (any(dst>size(grid,dim=4)))           .or. &
             (size(coord,dim=1)/=3)                .or. &
             (size(coord,dim=2)/=size(grid,dim=1)) .or. &
             (size(coord,dim=3)/=size(grid,dim=2)) .or. &
             (size(coord,dim=4)/=size(grid,dim=3)) ) then
          stop 'import_gamess%gamess_load_orbitals - Argument array sizes are not compatible'
        end if
      end if
      !
      !  If no implicit load is performed, orbitals must already be present
      !
      if (evaluate_mos .and. .not.load_mos) then
        if (.not.allocated(gam%atoms) .or. .not.allocated(gam%vectors)) then
          write (out,"('import_gamess%gamess_load_orbitals: basis set and/or orbitals are not present')") 
          stop 'import_gamess%gamess_load_orbitals: basis set and/or orbitals are not present'
        end if
      end if
      !
      !  After this point, we should not have major surprises accessing out data.
      !
      !
      !  Open GAMESS .DAT file, and parse the first $DATA and $VEC sections we see.
      !
      if (load_mos) then
        call gamess_destroy(gam)
        call TimerStart('GAMESS import - parse')
        open (gam_file,file=file,action='read',position='rewind',status='old',iostat=ios)
        if (ios/=0) then
          write (out,"('import_gamess%gamess_load_orbitals - error ',i8,' opening ',a)") ios, trim(file)
          stop 'import_gamess%gamess_load_orbitals - nofile'
        end if
        have_geo = .false.
        have_vec = .false.
        have_ecp = .false.
        scan_lines: do while (.not.have_geo .or. .not.have_vec .or. .not.have_ecp)
          call gam_readline(eof)
          if (eof) exit scan_lines
          !
          if (gam_line_buf==' $DATA  ') then
            if (have_geo) stop 'import_gamess%gamess_load_orbitals - multiple $DATA ?!'
            call parse_geometry_data(gam)
            have_geo = .true.
          end if
          if (gam_line_buf(1:5)==' $ECP') then
            if (have_ecp) stop 'import_gamess%gamess_load_orbitals - multiple $ECP ?!'
            if (.not.have_geo) stop 'import_gamess%gamess_load_orbitals - not $DATA before $ECP ?!'
            call parse_ecp_data(gam)
            have_ecp = .true.
          end if
          if (gam_line_buf(1:5)==' $VEC') then
            if (have_vec) stop 'import_gamess%gamess_load_orbitals - multiple $VEC ?!'
            if (.not.have_geo) stop 'import_gamess%gamess_load_orbitals - no $DATA before $VEC ?!'
            call parse_orbital_data(gam)
            have_vec = .true.
          end if
        end do scan_lines
        close(gam_file,iostat=ios)
        call TimerStop('GAMESS import - parse')
      end if ! load_mos
      !
      !  MO vectors are loaded; proceed to evaluate MOs on grid.
      !
      if (evaluate_mos) then
        call TimerStart('GAMESS import - orbitals')
        call build_orbitals(gam,mos,dst,coord,grid)
        call TimerStop('GAMESS import - orbitals')
      end if ! evaluate_mos
      !
      !  Special logic for the legacy semantics - MO matrix may be quite large, 
      !  and it won't be needed after this point.
      !
      if (load_mos .and. evaluate_mos) then
        deallocate (gam%vectors)
      end if
    end subroutine gamess_load_orbitals
    !
    !  Evaluate basis functions and their gradients on grid
    !
    subroutine gamess_evaluate_functions(xyz,basval,structure,rot)
      real(rk), intent(in)                :: xyz(:)      ! Coordinates of the point, at which AOs must be evaluated.
      real(rk), intent(out)              :: basval(:,:) ! Basis functions and gradients at a grid point. The first
                                                         ! index is (function,d/dx,d/dy,d/dz) or just (function)
      type(gam_structure), intent(inout), target, optional &
                                          :: structure   ! Context to use
      real(rk), intent(in), optional     :: rot(:,:)    ! Rotation matrix to apply to the structure. It is
                                                         ! perfectly OK to use different rotation matrix
                                                         ! in different calls using the same context. 
                                                         ! If the rotation matrix is not given, rotation matrix 
                                                         ! from the previous call will be used.
      !
      type(gam_structure), pointer :: gam
      real(xrk)                    :: det_rot
      !
      call TimerStart('GAMESS Evaluate functions')
      !
      if (present(structure)) then
        gam => structure
      else
        gam => gam_def
      end if
      !
      !  Deal with the rotation matrix
      !
      if (present(rot)) then
        gam%rotmat = real(rot,kind=kind(gam%rotmat))
        det_rot = lapack_determinant(gam%rotmat)
        if (abs(det_rot-1)>spacing(100._rk)) then
           write (out,"('Rotation matrix has non-unit determinant: ',g20.13)") det_rot
           stop 'import_gamess%gamess_evaluate_functions: Bad rotation matrix'
        end if
        gam%rotmat_rk = real(gam%rotmat,kind=kind(gam%rotmat_rk))
      end if
      !
      if (size(xyz)/=3 .or. all(size(basval,dim=1)/=(/1,4/)) .or. size(basval,dim=2)/=gam%nbasis) then
        stop 'import_gamess%gamess_evaluate_functions - Argument array sizes are not valid'
      end if
      !
      if (.not.allocated(gam%atoms)) then
        write (out,"('import_gamess%gamess_evaluate_functions: basis set is not present')") 
        stop 'import_gamess%gamess_evaluate_functions: basis set is not present'
      end if
      !
      !  After this point, we should not have major surprises accessing out data.
      !
      if (size(basval,dim=1)==1) then
        call evaluate_basis_functions(gam,xyz,basval(1,:))
      else
        call evaluate_basis_functions_and_gradients(gam,xyz,basval)
      end if
      call TimerStop('GAMESS Evaluate functions')
      !
    end subroutine gamess_evaluate_functions
    !
    !  Load the natural orbital occupation numbers
    !
    subroutine gamess_load_natocc(file,natural_occ,natural_count)
      character(len=*), intent(in) :: file            ! Name of the GAMESS checkpoint file
      real(rk), intent(out)        :: natural_occ(:)  ! Natural orbital occupation numbers
      integer(ik), intent(out)     :: natural_count   ! Number of natural orbitals in the density matrix
      !
      integer(ik) :: ios
      logical     :: have_natocc
      !
      !  Open GAMESS .DAT file, and parse the first $OCCNO section we see
      !
      call TimerStart('GAMESS import natocc')
      open (gam_file,file=file,action='read',position='rewind',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"('import_gamess%gamess_load_orbitals - error ',i8,' opening ',a)") ios, trim(file)
        stop 'import_gamess%gamess_load_natocc - nofile'
      end if
      have_natocc = .false.
      scan_lines: do while (.not.have_natocc)
        call gam_readline 
        if (gam_line_buf==' $OCCNO ') then
          call parse_natocc_data(natural_occ,natural_count)
          have_natocc = .true.
        end if
      end do scan_lines
      close(gam_file,iostat=ios)
      call TimerStop('GAMESS import natocc')
    end subroutine gamess_load_natocc
    !
    !  Load the transition density decomposition singular values 
    !
    subroutine gamess_load_rdmsv(file,rdm_sv,rdm_count)
      character(len=*), intent(in) :: file            ! Name of the GAMESS checkpoint file
      real(rk), intent(out)        :: rdm_sv(:)       ! Density matrix singular values
      integer(ik), intent(out)     :: rdm_count       ! Number of rdm orbitals
      !
      integer(ik) :: ios
      logical     :: have_rdm
      !
      !  Open GAMESS .DAT file, and parse the first $RDMPCE section we see
      !
      call TimerStart('GAMESS import rdmsv')
      open (gam_file,file=file,action='read',position='rewind',status='old',iostat=ios)
      if (ios/=0) then
        write (out,"('import_gamess%gamess_load_orbitals - error ',i8,' opening ',a)") ios, trim(file)
        stop 'import_gamess%gamess_load_rdm_sv - nofile'
      end if
      have_rdm = .false.
      scan_lines: do while (.not.have_rdm)
        call gam_readline 
        if (gam_line_buf==' $RDMPCE') then
          call parse_rdmsv_data(rdm_sv,rdm_count)
          have_rdm = .true.
        end if
      end do scan_lines
      close(gam_file,iostat=ios)
      call TimerStop('GAMESS import rdmsv')
    end subroutine gamess_load_rdmsv
    !
    !  Evaluate 1-electron integrals. Not much is implemented at the moment ...
    !  Numerical efficiency is NOT an issue - in every conceivable usage case,
    !  the overall cost will be dominated by something else.
    !
    !  Note that complex operators are handled in gamess_1e_integrals_complex below
    !
    subroutine gamess_1e_integrals_real(what,v,bra,ket,op_xyz,op_param,op_index)
      character(len=*), intent(in)      :: what        ! What to evaluate
      real(rk), intent(out)             :: v(:,:)      ! Buffer for the results
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      real(rk), intent(in), optional    :: op_xyz(:)   ! Position of the operator (3c integrals)
      real(rk), intent(in), optional    :: op_param(:) ! Additional parameters for the operator
                                                       ! Required size and allowed values depend on
                                                       ! the operator; see below.
      integer(ik), intent(in), optional :: op_index(:) ! Additional (integer) parameters for the operator.
      !
      include 'import_gamess_1e_integrals_common.f90'
    end subroutine gamess_1e_integrals_real
    !
    subroutine gamess_1e_integrals_quad(what,v,bra,ket,op_xyz,op_param,op_index)
      character(len=*), intent(in)      :: what        ! What to evaluate
      real(xrk), intent(out)            :: v(:,:)      ! Buffer for the results
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      real(xrk), intent(in), optional   :: op_xyz(:)   ! Position of the operator (3c integrals)
      real(xrk), intent(in), optional   :: op_param(:) ! Additional parameters for the operator
                                                       ! Required size and allowed values depend on
                                                       ! the operator; see below.
      integer(ik), intent(in), optional :: op_index(:) ! Additional (integer) parameters for the operator.
      !
      include 'import_gamess_1e_integrals_common.f90'
    end subroutine gamess_1e_integrals_quad
    !
    !  1-electron integrals for complex operators.
    !
    subroutine gamess_1e_integrals_complex(what,v,bra,ket,op_xyz,op_param,op_index)
      character(len=*), intent(in)   :: what         ! What to evaluate
      complex(rk), intent(out)       :: v(:,:)       ! Buffer for the results
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      real(rk), intent(in), optional  :: op_xyz(:)   ! Position of the operator (3c integrals)
      real(rk), intent(in), optional  :: op_param(:) ! Additional parameters for the operator
                                                     ! Required size and allowed values depend on
                                                     ! the operator; see below.
      integer(ik), intent(in), optional :: op_index(:) ! Additional (integer) parameters for the operator.
      !
      character(len=40)            :: lwhat          ! Local copy of "what"
      type(gam_structure), pointer :: l, r
      type(gam_operator_data)      :: op_data        ! Operator data
      !
      call TimerStart('GAMESS '//trim(what))
      l => gam_def ; r => gam_def
      if (present(bra)) l => bra
      if (present(ket)) r => ket
      !
      !  What is it we want to evaluate?
      !
      lwhat = what
      select case (lwhat(1:3))
        case ('AO ')
          if ( (size(v,dim=1)<l%nbasis) .or. (size(v,dim=2)<r%nbasis) ) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: Output buffer is not large enough')")
            stop 'import_gamess%gamess_1e_integrals_complex: Output buffer is not large enough'
          end if
      end select
      !
      !  3-centre integrals need position of the operator; make sure it was supplied
      !  Operator-specific parameters will be handled later.
      !
      if (lwhat(4:5)=='3C') then
        if (.not.present(op_xyz)) then
          write (out,"('import_gamess%gamess_1e_integrals_complex: operator location missing for ',a)") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals_complex: Missing operator position'
        end if
        op_data%op_xyz = real(op_xyz,kind=kind(op_data%op_xyz))
      end if
      !
      op_data%op_name = lwhat(4:)
      !
      select case (lwhat)
        case default
          write (out,"('import_gamess%gamess_1e_integrals_complex: ',a,' are not implemented')") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals_complex - integral not implemented'
        !
        !  First the integrals where no additional parmaters are needed:
        !
        !  'AO PLANEWAVE'  - Complex planewave operator [Exp(I K R)]; needs K vector
        !
        case ('AO PLANEWAVE')
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Missing operator parameter'
          end if
          if (size(op_param)/=3) then
            write (out,"('import_gamess%gamess_1e_integrals_complex: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals_complex: Wrong parameter count'
          end if
          op_data%kvec = real(op_param(1:3),kind=kind(op_data%kvec))
          call os_1e_matrix(op_data,l,r,v)
      end select
      call TimerStop('GAMESS '//trim(what))
    end subroutine gamess_1e_integrals_complex
    !
    !  Print 1-electron integrals in a somewhat reasonable form
    !
    subroutine gamess_print_1e_integrals_real(v,bra,ket,symmetry,heading)
      real(rk), intent(in)                              :: v(:,:)    ! Integrals to print
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      character(len=*), intent(in), optional            :: symmetry  ! Expected symmetry of matrix elements, either
                                                                     ! 'HERMITIAN' (default), 'ANTI-HERMITIAN', or 'NONE'
      character(len=*), intent(in), optional            :: heading   ! Headings to include, either 'BOTH', 'LEFT', or 'NONE'
                                                                     ! (the latter options are useful for MO printing)
      !
      include 'import_gamess_print_1e_integrals_common.f90'
    end subroutine gamess_print_1e_integrals_real
    !
    !  Print 1-electron integrals in a somewhat reasonable form (quad precision)
    !
    subroutine gamess_print_1e_integrals_quad(v,bra,ket,symmetry,heading)
      real(xrk), intent(in)                             :: v(:,:)    ! Integrals to print
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      character(len=*), intent(in), optional            :: symmetry  ! Expected symmetry of matrix elements, either
                                                                     ! 'HERMITIAN' (default), 'ANTI-HERMITIAN', or 'NONE'
      character(len=*), intent(in), optional            :: heading   ! Headings to include, either 'BOTH', 'LEFT', or 'NONE'
                                                                     ! (the latter options are useful for MO printing)
      !
      include 'import_gamess_print_1e_integrals_common.f90'
    end subroutine gamess_print_1e_integrals_quad
    !
    subroutine gamess_print_1e_integrals_complex(v,bra,ket,symmetry,heading)
      complex(rk), intent(in)                           :: v(:,:) ! Integrals to print
      type(gam_structure), intent(in), target, optional :: bra
      type(gam_structure), intent(in), target, optional :: ket
      character(len=*), intent(in), optional            :: symmetry  ! Expected symmetry of matrix elements, either
                                                                     ! 'HERMITIAN' (default), 'ANTI-HERMITIAN', 'SYMMETRIC', or 'NONE'
      character(len=*), intent(in), optional            :: heading   ! Headings to include, either 'BOTH', 'LEFT', or 'NONE'
                                                                     ! (the latter options are useful for MO printing)
      !
      type(gam_structure), pointer :: l, r
      integer(ik)                  :: ir, ic, ic_p1, ic_pn
      integer(ik), parameter       :: cpp = 5               ! Columns per page
      logical                      :: usegfmt
      character(len=20)            :: sym, head
      character(len=20)            :: lab
      !
      call TimerStart('GAMESS print 1e')
      l => gam_def ; r => gam_def
      if (present(bra)) l => bra
      if (present(ket)) r => ket
      !
      !  Symmetry check
      !
      sym = 'HERMITIAN'
      if (present(symmetry)) sym = symmetry
      if (size(v,dim=1)==size(v,dim=2)) then
        write (out,"('Expected integral symmetry: ',a)") trim(sym)
        select case (sym)
          case default
            write (out,"('Symmetry code ',a,' is not recognized')") trim(sym)
            stop 'import_gamess%gamess_print_1e_integrals_real - bad argument'
          case ('NONE')
          case ('HERMITIAN')
            write (out,"('Maximum deviation from symmetry = ',g12.5)") maxval(abs(v-conjg(transpose(v))))
            write (out,"('     Maximum deviation at index = ',2i6  )") maxloc(abs(v-conjg(transpose(v))))
          case ('ANTI-HERMITIAN')
            write (out,"('Maximum deviation from symmetry = ',g12.5)") maxval(abs(v+conjg(transpose(v))))
            write (out,"('     Maximum deviation at index = ',2i6  )") maxloc(abs(v+conjg(transpose(v))))
          case ('SYMMETRIC')
            write (out,"('Maximum deviation from symmetry = ',g12.5)") maxval(abs(v-transpose(v)))
            write (out,"('     Maximum deviation at index = ',2i6  )") maxloc(abs(v-transpose(v)))
        end select
            write (out,"('             Largest matrix element = ',g12.5)") maxval(abs(v))
            write (out,"(' Largest matrix element is at index = ',2i6)") maxloc(abs(v))
      end if
      !
      head = 'BOTH'
      if (present(heading)) head = heading
      usegfmt = any(abs(v)>=1e5) .or. all(abs(v)<=1e-3)
      print_pages: do ic_p1=1,r%nbasis,cpp
        ic_pn = min(r%nbasis,ic_p1+cpp-1)
        write (out,"(t17,5(2x,i12,1x,12x))") (ic,               ic=ic_p1,ic_pn)
        if (head=='BOTH') write (out,"(t23,5(2x,a12,1x,12x))") (r%bas_labels(ic), ic=ic_p1,ic_pn)
        print_rows: do ir=1,l%nbasis
          lab = ' '
          if (head/='NONE') lab = l%bas_labels(ir)
          if (usegfmt) then
            write (out,"(1x,i5,1x,a11,1x,5(2x,g12.5,1x,g12.5))") ir, lab, v(ir,ic_p1:ic_pn)
          else
            write (out,"(1x,i5,1x,a11,1x,5(2x,f12.5,1x,g12.5))") ir, lab, v(ir,ic_p1:ic_pn)
          endif
        end do print_rows
        write (out,"()")
      end do print_pages
      !
      call TimerStop('GAMESS print 1e')
    end subroutine gamess_print_1e_integrals_complex
    !
    !  Return key information on 2e integrals and integral batches
    !
    subroutine gamess_2e_info(what,gam,op_iparam)
      character(len=*), intent(in)                         :: what         ! What to do, see below
      type(gam_structure), intent(inout), target, optional :: gam          ! Context to work on
      integer(ik), intent(inout), optional                 :: op_iparam(:) ! Additional parameters
      !
      type(gam_structure), pointer :: g
      integer(ik)                  :: batch
      !
      call TimerStart('GAMESS 2e '//trim(what))
      g => gam_def
      if (present(gam)) g => gam
      !
      select case (what)
        case default
          write (out,"('import_gamess%gamess_2e_info: ',a,' are not implemented')") trim(what)
          stop 'import_gamess%gamess_2e_info - operation not implemented'
        case ('GET INFO')
          call need_iparam(3_ik)
          op_iparam(1) = g%nbatches   ! Output: number of batches
          op_iparam(2) = g%max_batch  ! Output: number of orbitals in the largest batch
          op_iparam(3) = g%nbasis     ! Output: total number of atomic orbitals
        case ('GET BATCH INFO')
          call need_iparam(4_ik)
          batch = op_iparam(1)        ! Input: desired batch
          if (batch<1 .or. batch>g%nbatches) then
            write (out,"('import_gamess%gamess_2e_info: batch ',i10,' is outside of the valid range')") batch
            stop 'import_gamess%gamess_2e_info - bad batch'
          end if
          op_iparam(2) = g%batches(4,batch) ! Output: size of this batch
          op_iparam(3) = g%batches(3,batch) ! Output: first orbital of this batch
          op_iparam(4) = g%batches(1,batch) ! Output: atom this batch belongs to
      end select
      !
      call TimerStop('GAMESS 2e '//trim(what))
      !
      contains
      !
      subroutine need_iparam(count)
        integer(ik), intent(in) :: count ! Required number of entries in op_iparam
        !
        if (.not.present(op_iparam)) stop 'import_gamess%gamess_2e_info - required parameter op_iparam missing'
        if (size(op_iparam)<count) stop 'import_gamess%gamess_2e_info - op_iparam is too small'
      end subroutine need_iparam
    end subroutine gamess_2e_info
    !
    !  Evaluate a batch of 2e integrals. Due to the number of 2e integrals in
    !  a typical calculation, it is impractical to evaluate them all in-memory,
    !  as we do for the 1e integrals above.
    !
    subroutine gamess_2e_integrals_real(what,v,batch,a,b,c,d,op_param,accuracy)
      character(len=*), intent(in)                      :: what        ! What to evaluate
      real(rk), intent(out)                             :: v(:,:,:,:)  ! Buffer for the results
                                                                       ! Integrals are in the charge-cloud order: that is,
                                                                       ! first two indices correspond to the first variable;
                                                                       ! the last two indices correspond to the second variable
      integer(ik), intent(in)                           :: batch(:)    ! Batch indices
      type(gam_structure), intent(in), target, optional :: a, b, c, d  ! Contexts to work on. In principle, we can handle
                                                                       ! the situation with every index being on a separate
                                                                       ! molecule - not that this is hard for an AO integral
      real(rk), intent(in), optional                    :: op_param(:) ! Additional parameters for the operator
      real(rk), intent(in), optional                    :: accuracy    ! Desired accuracy of the integrals; contributions
                                                                       ! smaller than this can be neglected. 
                                                                       ! The default is to be as accurate as possible.
                                                                       ! Negative accuracy will restore the default behaviour
      !
      include 'import_gamess_gamess_2e_integrals_common.f90'
    end subroutine gamess_2e_integrals_real
    !
    subroutine gamess_2e_integrals_quad(what,v,batch,a,b,c,d,op_param,accuracy)
      character(len=*), intent(in)                      :: what        ! What to evaluate
      real(xrk), intent(out)                            :: v(:,:,:,:)  ! Buffer for the results
                                                                       ! Integrals are in the charge-cloud order: that is,
                                                                       ! first two indices correspond to the first variable;
                                                                       ! the last two indices correspond to the second variable
      integer(ik), intent(in)                           :: batch(:)    ! Batch indices
      type(gam_structure), intent(in), target, optional :: a, b, c, d  ! Contexts to work on. In principle, we can handle
                                                                       ! the situation with every index being on a separate
                                                                       ! molecule - not that this is hard for an AO integral
      real(xrk), intent(in), optional                   :: op_param(:) ! Additional parameters for the operator
      real(xrk), intent(in), optional                   :: accuracy    ! Desired accuracy of the integrals; contributions
                                                                       ! smaller than this can be neglected. 
                                                                       ! The default is to be as accurate as possible.
                                                                       ! Negative accuracy will restore the default behaviour
      !
      include 'import_gamess_gamess_2e_integrals_common.f90'
    end subroutine gamess_2e_integrals_quad
    !
    !  Internal routines, only callable from module external interface parts
    !
    subroutine parse_geometry_data(gam)
      type(gam_structure), intent(inout) :: gam      ! Context to work on
      integer(ik)                        :: iat, ish
      !
      allocate (gam%atoms(gam_max_atoms))
      ! 
      !  Comment line - ignore
      !  
      call gam_readline ; call echo(1_ik)
      !
      !  Symmetry line - ignore as we now support symmetry in the ADC code
      !
      call gam_readline ; call echo(1_ik)
      !
      !  Begin loading geometry, COORD=UNIQUE format
      !
      gam%natoms = 0
      load_atoms: do 
        call gam_readline ; call echo(1_ik)
        if (gam_line_buf==' $END') exit load_atoms
        gam%natoms = gam%natoms + 1
        if (gam%natoms>gam_max_atoms) then
          write (out,"('import_gamess%parse_geometry_data: Exceeded the max. number of atoms: ',i6)") gam_max_atoms
          stop 'import_gamess%parse_geometry_data - too many atoms'
        end if
        call parse_one_atom(gam%atoms(gam%natoms))
      end do load_atoms
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' atoms')") gam%natoms
      end if
      if (gam%natoms==0) then
        stop 'import_gamess%parse_geometry_data - Found no atoms!'
      end if
      ! 
      !  Count basis functions and contractions
      !  
      gam%nbasis = 0
      gam%max_contractions = 0
      count_basis: do iat=1,gam%natoms
        count_shells: do ish=1,gam%atoms(iat)%nshell
          gam%nbasis = gam%nbasis + gam_orbcnt(gam%atoms(iat)%sh_l(ish))
          gam%max_contractions = max(gam%max_contractions,gam%atoms(iat)%sh_p(ish+1)-gam%atoms(iat)%sh_p(ish))
        end do count_shells
      end do count_basis
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' Cartesian basis functions')") gam%nbasis
        write (out,"('import_gamess: Longest contraction has ',i8,' primitives')") gam%max_contractions
      end if
      !
      !  Count batches, may be needed for 2e integrals
      !
      call count_2e_batches(gam)
      !
      !  Label basis functions, in case anybody cares
      !
      call label_basis_set(gam)
    end subroutine parse_geometry_data
    !
    !  Initialize integral batch tables for this structure
    !
    subroutine count_2e_batches(gam)
      type(gam_structure), intent(inout) :: gam ! Context to work on
      !
      integer(ik)  :: iat, ish           ! Atom and shell counters
      integer(ik)  :: batch_shells       ! Number of shells in the current batch
      integer(ik)  :: alloc
      integer(ik)  :: ibatch             ! Current shell batch
      integer(ik)  :: next_batch_index   ! Orbital index for the next shell batch
      integer(ik)  :: batch_size         ! Number of orbitals in the current batch.
                                         ! Should not exceed gam_max_batch_size
      integer(ik)  :: shell_size         ! Number of orbitals in the current shell
      integer(ik)  :: sh_l               ! Angular momentum of the current shell
      !
      !  Count batches 
      !
      gam%nbatches   = 0
      gam%max_shells = 0
      !
      count_2e_batches_atoms: do iat=1,gam%natoms
        if (gam%atoms(iat)%nshell<=0) cycle count_2e_batches_atoms
        gam%nbatches   = gam%nbatches + 1 
        batch_shells   = 1
        gam%max_shells = max(gam%max_shells,batch_shells)
        sh_l           = gam%atoms(iat)%sh_l(1)
        batch_size     = gam_orbcnt(sh_l)
        count_2e_batches_shells: do ish=2,gam%atoms(iat)%nshell
          sh_l       = gam%atoms(iat)%sh_l(ish)
          shell_size = gam_orbcnt(sh_l)
          if (sh_l==gam%atoms(iat)%sh_l(ish-1) .and. batch_size+shell_size<=gam_max_batch_size) then
            batch_shells   = batch_shells + 1
            gam%max_shells = max(gam%max_shells,batch_shells)
            batch_size     = batch_size + shell_size
          else
            !
            !  L changed or shell became too large, advance shell counter
            !
            batch_size   = shell_size
            gam%nbatches = gam%nbatches + 1
            batch_shells = 1
          end if
        end do count_2e_batches_shells
      end do count_2e_batches_atoms
      !
      !  Allocate batch table
      !
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' shell batches')") gam%nbatches
        write (out,"('import_gamess: Largest batch contains ',i8,' shells')") gam%max_shells
      end if
      allocate (gam%batches(4+gam%max_shells,gam%nbatches),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating batches table of size ',i8,' by ',i8)") &
               alloc, 4+gam%max_shells, gam%nbatches
        stop 'import_gamess%parse_geometry_data - batch table allocation failed'
      end if
      gam%batches = 0 ! This initialization is important for stuff_shell below!
      !
      ibatch = 0
      next_batch_index = 1
      fill_2e_batches_atoms: do iat=1,gam%natoms
        if (gam%atoms(iat)%nshell<=0) cycle fill_2e_batches_atoms
        ibatch     = ibatch + 1
        sh_l       = gam%atoms(iat)%sh_l(1)
        batch_size = gam_orbcnt(sh_l)
        call stuff_shell(iat,1_ik)
        fill_2e_batches_shells: do ish=2,gam%atoms(iat)%nshell
          sh_l       = gam%atoms(iat)%sh_l(ish)
          shell_size = gam_orbcnt(sh_l)
          if (sh_l/=gam%atoms(iat)%sh_l(ish-1) .or. batch_size+shell_size>gam_max_batch_size) then
             ibatch     = ibatch + 1
             batch_size = shell_size
          else
             batch_size = batch_size + shell_size
          end if
          call stuff_shell(iat,ish)
        end do fill_2e_batches_shells
      end do fill_2e_batches_atoms
      !
      if (ibatch/=gam%nbatches) stop 'import_gamess%count_2e_batches - batch count error (2)'
      !
      !  Fill in the max. number of orbitals in a batch
      !
      gam%max_batch = maxval(gam%batches(4,:))
      !
      if (verbose>=0) then
        write (out,"('import_gamess: Largest batch contains ',i8,' orbitals')") gam%max_batch
      end if
      if (verbose>=1) then
        write (out,"(/'Integral batch table:')")
        write (out,"(1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a10)") 'Batch', 'Atom', '# shells', '1st orb', '# orbs', 'shells ...'
        print_batches: do ibatch=1,gam%nbatches
          write (out,"(1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,20(1x,i4))") ibatch, gam%batches(:,ibatch)
        end do print_batches
        write (out,"()")
      end if
      !
      contains
      subroutine stuff_shell(iat,ish)
        integer(ik), intent(in) :: iat, ish
        integer(ik) :: ns ! Number of shells in this batch
        integer(ik) :: no ! Number of orbitals in this shell
        !
        ns = gam%batches(2,ibatch) + 1
        no = gam_orbcnt(gam%atoms(iat)%sh_l(ish))
        !
        if (ns    >gam%max_shells) stop 'import_gamess%count_2e_batches - shell count error'
        if (ibatch>gam%nbatches  ) stop 'import_gamess%count_2e_batches - batch count error (1)'
        !
                   gam%batches(1   ,ibatch) = iat
                   gam%batches(2   ,ibatch) = ns
        if (ns==1) gam%batches(3   ,ibatch) = next_batch_index
                   gam%batches(4   ,ibatch) = gam%batches(4,ibatch) + no
                   gam%batches(4+ns,ibatch) = ish
        !
        next_batch_index = next_batch_index + no
      end subroutine stuff_shell
    end subroutine count_2e_batches
    !
    subroutine parse_one_atom(atom)
      type(gam_atom), intent(out) :: atom
      !
      integer(ik)       :: ios
      character(len=20) :: lcode
      integer(ik)       :: idum
      integer(ik)       :: nprim, ip, pprim, sh, pl, pu, ip1, ip2
      real(xrk)         :: norm
      !
      !  Coordinates line
      !
      read (gam_line_buf,*,iostat=ios) atom%name, atom%znuc, atom%xyz
      if (ios/=0) then
        write (out,"('import_gamess%parse_one_atom: Error ',i8,' parsing line:'/1x,a)") ios, trim(gam_line_buf)
        stop 'import_gamess%parse_one_atom - coord read'
      end if
      atom%xyz_rk = real(atom%xyz,kind=kind(atom%xyz_rk))
      !
      !  Dummy ECP entry
      !
      atom%ecp_name   = 'NONE'
      atom%ecp_zcore  = 0._xrk
      atom%ecp_nterms = 0
      !
      !  Basis set lines
      !
      atom%nshell  = 0 ;
      atom%sh_p(1) = 1 ;
      read_shells: do 
        call gam_readline ; call echo(1_ik)
        if (gam_line_buf==' ') exit read_shells ! End of atom data
        read (gam_line_buf,*,iostat=ios) lcode, nprim
        if (ios/=0) then
          write (out,"('import_gamess%parse_one_atom: Error ',i8,' parsing line:'/1x,a)") ios, trim(gam_line_buf)
          stop 'import_gamess%parse_one_atom - lcode read'
        end if
        !
        if (nprim>GAM_MAX_CONTRACTION) then
          write (out,"('import_gamess%parse_one_atom: too many primitives at line ',a)") trim(gam_line_buf)
          write (out,"('import_gamess%parse_one_atom: increase GAM_MAX_CONTRACTION and recompile')")
          stop 'import_gamess%parse_one_atom - too many contractions'
        end if
        !
        atom%nshell = atom%nshell + 1
        if (atom%nshell>=GAM_MAX_SHELLS) then
          write (out,"('import_gamess%parse_one_atom: too many shells at line ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_one_atom - too many shells'
        end if
        !
        select case (lcode)
          case default
            write (out,"('import_gamess%parse_one_atom: lcode ',a,' is not recognized')") trim(gam_line_buf)
            stop 'import_gamess%parse_one_atom - bad lcode'
          case ('S') ; atom%sh_l(atom%nshell) = 0
          case ('P') ; atom%sh_l(atom%nshell) = 1
          case ('D') ; atom%sh_l(atom%nshell) = 2
          case ('F') ; atom%sh_l(atom%nshell) = 3
          case ('G') ; atom%sh_l(atom%nshell) = 4
        end select
        !
        pprim = atom%sh_p(atom%nshell)
        read_contractions: do ip=1,nprim
          call gam_readline ; call echo(1_ik)
          if (pprim>=GAM_MAX_PRIMITIVE) then
          write (out,"('import_gamess%parse_one_atom: too many primitives at line ',a)") trim(gam_line_buf)
            stop 'import_gamess%parse_one_atom - too many primitives'
          end if
          read (gam_line_buf,*,iostat=ios) idum, atom%p_zet(pprim),atom%p_c(pprim)
          if (ios/=0) then
            write (out,"('import_gamess%parse_one_atom: Error ',i8,' parsing line:'/1x,a)") &
                   ios, trim(gam_line_buf)
            stop 'import_gamess%parse_one_atom - primitive read'
          end if
          !
          !  Save a copy of the original contraction coefficients
          ! 
          atom%p_c_orig(pprim)=atom%p_c(pprim)
          !
          !  Rescale contraction coefficient to include primitive normalization factor
          !
          atom%p_c(pprim) = atom%p_c(pprim) * primitive_norm(atom%sh_l(atom%nshell),atom%p_zet(pprim))
          ! write (out,*) 'ip = ', ip, ' pprim = ', pprim, ' zeta = ', atom%p_zet(pprim), ' c = ', atom%p_c(pprim)
          pprim = pprim + 1
        end do read_contractions
        !
        atom%sh_p(atom%nshell+1) = pprim ;
      end do read_shells
      !
      !  Check normalization of each shell, and renormalize if needed.
      !
      renormalize: do sh=1,atom%nshell
        pl   = atom%sh_p(sh)
        pu   = atom%sh_p(sh+1) - 1
        norm = 0._xrk
        left: do ip1=pl,pu
          right: do ip2=pl,pu
            norm = norm + atom%p_c(ip1) * atom%p_c(ip2) &
                    / (atom%p_zet(ip1)+atom%p_zet(ip2))**(1.5_xrk+atom%sh_l(sh))
          end do right
        end do left
        norm = norm * (pi_xrk**1.5_xrk)*gam_normc(atom%sh_l(sh))
        !
        if ( (verbose>=0) .and. (abs(norm-1.0_xrk)>1e-6_xrk) ) then
          write (out,"('import_gamess: Shell ',i4,' lval = ',i2,' with ',i4,' primitives, norm = ',g20.12)") &
                 sh, atom%sh_l(sh), pu-pl+1, norm
        end if
        !
        !  Renormalize
        !
        atom%p_c(pl:pu) = atom%p_c(pl:pu) / sqrt(norm)
        !
        !  Make reduced-precision copy
        !
        atom%p_c_rk  (pl:pu) = real(atom%p_c  (pl:pu),kind=kind(atom%p_c_rk  ))
        atom%p_zet_rk(pl:pu) = real(atom%p_zet(pl:pu),kind=kind(atom%p_zet_rk))
      end do renormalize
    end subroutine parse_one_atom
    !
    subroutine parse_ecp_data(gam)
      type(gam_structure), intent(inout) :: gam      ! Context to work on
      !
      integer(ik)                        :: iat      ! Atom we are currently working on
      integer(ik)                        :: iax      ! Another atom
      character(len=10)                  :: ptype    ! Dummy variable for reading
      integer(ik)                        :: ios      ! I/O status
      !
      !  Can only parse a subset of possible PP=READ $ECP inputs; see GAMESS manual
      !  for the detailed description of what we are missing here
      !
      iat = 0
      load_atoms: do 
        call gam_readline ; call echo(1_ik)
        if (gam_line_buf==' $END') exit load_atoms
        !
        iat = iat + 1
        if (iat>gam%natoms) then
          write (out,"('import_gamess%parse_ecp_data: More ECPs than there are atoms.')") 
          write (out,"('Last processed line: ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_ecp_data - too many ECPs'
        end if
        !
        !  First line: PNAME, PTYPE, IZCORE, LMAX+1
        !
        ptype = 'GEN'
        read(gam_line_buf,*,iostat=ios) gam%atoms(iat)%ecp_name, ptype, gam%atoms(iat)%ecp_zcore, gam%atoms(iat)%ecp_lmaxp1
        if (ios>0) then  ! It's OK to have an end-of-record condition (ios<0)
          write (out,"('import_gamess%parse_ecp_data: Error ',i8,' parsing ECP for atom ',i8)") ios, iat
          write (out,"('Last processed line: ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_ecp_data - bad first line of an ECP'
        end if
        if (ptype=='NONE') cycle load_atoms ! No ECP on this atom, continue to the next
        if (ptype/='GEN') then
          write (out,"('import_gamess%parse_ecp_data: Atom ',i8,' has an unsupported ECP type ',a)") iat, trim(ptype)
          write (out,"('Last processed line: ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_ecp_data - unsupported ECP type'
        end if
        !
        !  Have we seen this ECP name before? If yes, copy it and proceed to the next atom
        !
        scan_known_ecps: do iax=1,iat-1
          if (gam%atoms(iat)%ecp_name==gam%atoms(iax)%ecp_name) then
            call copy_ecp(gam%atoms(iat),gam%atoms(iax))
            cycle load_atoms
          end if
        end do scan_known_ecps
        call parse_one_ecp(gam%atoms(iat))
      end do load_atoms
      if (iat/=gam%natoms) then
        write (out,"('import_gamess%parse_ecp_data: Found ',i6,' ECPs and ',i6,' atoms. Oops.')") iat, gam%natoms
        stop 'import_gamess%parse_ecp_data - $ECP does not match $DATA'
      end if
      !
      !  Adjust nuclear charges to reflect the long-range part of the potential, and do some printing
      !
      adjust_znuc: do iat=1,gam%natoms
        gam%atoms(iat)%znuc = gam%atoms(iat)%znuc - gam%atoms(iat)%ecp_zcore
        if (gam%atoms(iat)%ecp_name=='NONE' .or. verbose<=0) cycle adjust_znuc
        write (out,"('In atom ',i6,', ',f8.2,' core electrons are treated by an ECP ',a)") &
               iat, gam%atoms(iat)%ecp_zcore, gam%atoms(iat)%ecp_name
      end do adjust_znuc
    end subroutine parse_ecp_data
    !
    subroutine copy_ecp(dst,src)
      type(gam_atom), intent(inout) :: dst
      type(gam_atom), intent(in)    :: src
      !
      dst%ecp_name   = src%ecp_name
      dst%ecp_zcore  = src%ecp_zcore
      dst%ecp_lmaxp1 = src%ecp_lmaxp1
      dst%ecp_nterms = src%ecp_nterms
      dst%ecp_l      = src%ecp_l
      dst%ecp_n      = src%ecp_n
      dst%ecp_c      = src%ecp_c
      dst%ecp_d      = src%ecp_d
    end subroutine copy_ecp
    !
    subroutine parse_one_ecp(at)
      type(gam_atom), intent(inout) :: at
      !
      integer(ik)  :: lv, ig, it, ng, ios
      !
      !  The first entry is the "Lmax+1" term, which appears for
      !  all channels
      !
      at%ecp_nterms = 0
      read_channels: do lv=-1,at%ecp_lmaxp1-1
        call gam_readline ; call echo(2_ik)
        read(gam_line_buf,*,iostat=ios) ng
        if (ios/=0) then
          write (out,"('Error ',i6,' parsing channel header for L=',i3,' in ECP ',a)") ios, lv, trim(at%ecp_name)
          write (out,"('Last processed line: ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_one_ecp - error parsing channel header'
        end if
        if (at%ecp_nterms+ng>gam_max_ecp_terms) then
          write (out,"('Exceeded maximum number of ECP terms while reading L=',i3,' in ECP ',a,'. The maximum is ',i8)") &
                 lv, trim(at%ecp_name), gam_max_ecp_terms
          write (out,"('Last processed line: ',a)") trim(gam_line_buf)
          stop 'import_gamess%parse_one_ecp - error parsing channel header'
        end if
        read_gaussians: do ig=1,ng
          call gam_readline ; call echo(2_ik)
          at%ecp_nterms = at%ecp_nterms + 1
          it            = at%ecp_nterms
          read(gam_line_buf,*,iostat=ios) at%ecp_c(it), at%ecp_n(it), at%ecp_d(it)
          if (ios/=0) then
            write (out,"('Error ',i6,' parsing gaussian ',i3,' for L=',i3,' in ECP ',a)") ios, ig, lv, trim(at%ecp_name)
            write (out,"('Last processed line: ',a)") trim(gam_line_buf)
            stop 'import_gamess%parse_one_ecp - error parsing channel data'
          end if
          at%ecp_l(it) = lv
        end do read_gaussians
      end do read_channels
      !
      if (verbose>=1) then
        write (out,"(/'Parsed ECP ',a,' Zcore = ',f12.4,' Lmax+1 = ',i3)") trim(at%ecp_name), at%ecp_zcore, at%ecp_lmaxp1
        write (out,"(1x,a3,2x,a3,2x,a18,2x,a18)") ' L ', ' N ', '     C     ', '     D     ', &
                                                  '---', '---', '-----------', '-----------'
        print_terms: do ig=1,at%ecp_nterms
          write (out,"(1x,i3,2x,i3,2x,f18.10,2x,f18.10)") at%ecp_l(ig), at%ecp_n(ig), at%ecp_c(ig), at%ecp_d(ig)
        end do print_terms
        write (out,"()")
      end if
    end subroutine parse_one_ecp
    !
    subroutine parse_orbital_data(gam)
      type(gam_structure), intent(inout) :: gam
      !
      integer(ik) :: iao, line, ios, ifield
      integer(ik) :: chk_mo, chk_line
      integer(ik) :: wrap_mo
      real(rk)   :: loc_val(5)
      !
      allocate (gam%vectors(gam%nbasis,2*gam%nbasis))
      !
      gam%nvectors = 0
      iao          = 1
      line         = 0
      read_vectors: do
        call gam_readline ; call echo(2_ik)
        if (gam_line_buf==' $END') exit read_vectors
        line = line + 1
        if (gam%nvectors+1>2*gam%nbasis) then
          write (out,"('import_gamess%parse_orbital_data: Too many lines in $VEC at line:'/1x,a)") &
                 trim(gam_line_buf)
          stop 'import_gamess%parse_orbital_data - too many lines'
        end if
        !
        read (gam_line_buf,"(i2,i3,5g15.10)",iostat=ios) chk_mo, chk_line, loc_val
        if (ios/=0) then
          write (out,"('import_gamess%parse_orbital_data: format error ',i8,' at line:'/1x,a)") &
                 ios, trim(gam_line_buf)
          stop 'import_gamess%parse_orbital_data - format error'
        end if
        ! 
        !  Line parsed, check it for validity
        !
        wrap_mo = gam%nvectors+1
        if (wrap_mo>gam%nbasis) wrap_mo = wrap_mo - gam%nbasis
        if (modulo(wrap_mo,100)/=chk_mo) then
          write (out,"('import_gamess%parse_orbital_data: Incorrect data: expecting MO ',i8,', got ',i8,'. Line: '/1x,a)") &
                   wrap_mo, chk_mo, trim(gam_line_buf)
          stop 'import_gamess%parse_orbital_data - MO mismatch'
        end if
        if (line/=chk_line) then
          write (out,"('import_gamess%parse_orbital_data: Incorrect data: expecting line ',i8,', got ',i8,'. Line: '/1x,a)") &
                   line, chk_line, trim(gam_line_buf)
          stop 'import_gamess%parse_orbital_data - line mismatch'
        end if
        !
        !  Stuff vector in
        !
        stuff_c: do ifield=1,5
          gam%vectors(iao,gam%nvectors+1) = loc_val(ifield)
          iao = iao + 1
          if (iao>gam%nbasis) then
            gam%nvectors = gam%nvectors + 1
            iao          = 1
            line         = 0
            exit stuff_c
          end if
        end do stuff_c
      end do read_vectors
      !
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' vectors in the set')") gam%nvectors
      end if
    end subroutine parse_orbital_data
    !
    subroutine parse_natocc_data(natural_occ,natural_count)
      real(rk), intent(out)        :: natural_occ(:)  ! Natural orbital occupation numbers
      integer(ik), intent(out)     :: natural_count   ! Number of natural orbitals in the density matrix
      !
      integer(ik) :: ios
      integer(ik) :: line_numocc
      integer(ik) :: max_numocc
      real(rk)   :: loc_val(5)
      !
      max_numocc = size(natural_occ,dim=1)
      !
      natural_count = 0
      !
      read_natocc: do
        call gam_readline ; call echo(2_ik)
        if (gam_line_buf==' $END') exit read_natocc
        !
        read (gam_line_buf,"(5f16.10)",iostat=ios) loc_val
        if (ios/=0) then
          write (out,"('import_gamess%parse_natocc_data: format error ',i8,' at line:'/1x,a)") &
                 ios, trim(gam_line_buf)
          stop 'import_gamess%parse_natocc_data - format error (1)'
        end if
        ! 
        line_numocc = len_trim(gam_line_buf)/16
        !
        if (line_numocc>5) then
          write (out,"('import_gamess%parse_natocc_data: unexpected data at the end of line:'/1x,a)") &
                 trim(gam_line_buf)
          stop 'import_gamess%parse_natocc_data - format error (2)'
        end if
        !
        if ( (natural_count+line_numocc)>max_numocc) then
          write (out,"('import_gamess%parse_natocc_data: too many natural orbitals in file. Increase max_naturals and recompile.')") 
          stop 'import_gamess%parse_natocc_data - too many natural orbitals - increase max_naturals and recompile.'
        end if
        !
        natural_occ(natural_count+1:natural_count+line_numocc) = loc_val(1:line_numocc)
        natural_count = natural_count + line_numocc
        ! 
      end do read_natocc
      !
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' natural orbitals in the set')") natural_count
      end if
    end subroutine parse_natocc_data
    !
    subroutine parse_rdmsv_data(rdm_sv,rdm_count)
      real(rk), intent(out)        :: rdm_sv(:)   ! Density maxtrix singular values
      integer(ik), intent(out)     :: rdm_count   ! Number of density matrix orbitals
      !
      integer(ik) :: ios
      integer(ik) :: line_numsv
      integer(ik) :: max_numsv
      real(rk)   :: loc_val(5)
      !
      max_numsv = size(rdm_sv,dim=1)
      !
      rdm_count = 0
      !
      read_rdmsv: do
        call gam_readline ; call echo(2_ik)
        if (gam_line_buf==' $END') exit read_rdmsv
        !
        read (gam_line_buf,"(5f15.12)",iostat=ios) loc_val
        if (ios/=0) then
          write (out,"('import_gamess%parse_rdmsv_data: format error ',i8,' at line:'/1x,a)") &
                 ios, trim(gam_line_buf)
          stop 'import_gamess%parse_rdmsv_data - format error (1)'
        end if
        ! 
        line_numsv = len_trim(gam_line_buf)/15
        !
        if (line_numsv>5) then
          write (out,"('import_gamess%parse_rdmsv_data: unexpected data at the end of line:'/1x,a)") &
                 trim(gam_line_buf)
          stop 'import_gamess%parse_rdmsv_data - format error (2)'
        end if
        !
        if ( (rdm_count+line_numsv)>max_numsv) then
          write (out,"('import_gamess%parse_rdmsv_data: too many rdm orbitals in file. Increase max_rdm and recompile.')") 
          stop 'import_gamess%parse_rdmsv_data - too many rdm orbitals - increase max_rdm and recompile.'
        end if
        !
        rdm_sv(rdm_count+1:rdm_count+line_numsv) = loc_val(1:line_numsv)
        rdm_count = rdm_count + line_numsv
        ! 
      end do read_rdmsv
      !
      if (verbose>=0) then
        write (out,"('import_gamess: Found ',i8,' rdm orbitals in the set')") rdm_count
      end if
    end subroutine parse_rdmsv_data
    !
    subroutine echo(level)
      integer(ik), intent(in) :: level
      !
      if (level>verbose) return
      write (out,"('import_gamess: ',a)") trim(gam_line_buf)
    end subroutine echo
    !
    !  Actual evaluation of MOs on grid - everything above was just the preparatory
    !  work. The parameters are:
    !
    !    grid  - Locations of grid points, (X,Y,Z) triples
    !    nmos  - Number of MOs to evaluate
    !    mos   - Indices of MOs to evaluate (these may be non-sensical - watch out!)
    !    moval - Buffer for MO values, with values at each grid point grouped together.
    !
    subroutine build_orbitals(gam,mos,dst,coord,grid)
      type(gam_structure), intent(in) :: gam             ! Context to use
      integer(ik), intent(in)         :: mos(:)          ! Indices of the orbitals to be computed
      integer(ik), intent(in)         :: dst(:)          ! Indices of grid to receive the orbitals
      real(rk), intent(in)            :: coord(:,:,:,:)  ! Coordinates of the points, at which MOs
                                                         ! must be evaluated. First index is X,Y,Z.
      complex(rk), intent(out)        :: grid(:,:,:,:)   ! Data field for the MOs. First three
                                                         ! indices match last three indices of coord().
                                                         ! The last index is for the MOs.
      !
      integer(ik), allocatable  :: special_pt(:,:)    ! Grid points, where oversampling is deemed
                                                      ! necessary. The final 3 entries are the expected 
                                                      ! number or samples per point.
      real(rk)                  :: r_cut              ! Max distance from the nearest nucleus to 
                                                      ! trigger oversampling
      integer(ik)               :: num_special        ! Number of grid points requiring special treatment
      integer(ik)               :: num_special_max    ! (Estimate of the maximum possible # of the special points
      integer(ik)               :: ipx, ipy, ipz      ! Grid coordinates
      integer(ik)               :: isp                ! Special-point index
!
!    For some reason, GNU Fortran has trouble with private automatic arrays - they seem to
!    cause segmentation faults in parallel runs. We'll simply have to allicate basval 
!    inside the parallel section ....
!
!     real(rk)                 :: basval(gam%nbasis) ! local buffer for basis function values
      real(rk), allocatable    :: basval(:)          ! local buffer for basis function values
      integer(ik)              :: alloc_             ! For some reason, Intel Fortran 14.0 objects to "alloc" in !$omp directive
                                                     ! These guys introduce new bugs much faster than they can fix them ...
      !
      !  We would like to test for special points in parallel; this requires a
      !  reasonable estimate of how many such points we are goins to see.
      !
      call prepare_for_oversampling(gam,num_special_max,r_cut)
      if (num_special_max>0) then
        allocate (special_pt(6,num_special_max),stat=alloc_)
        if (alloc_/=0) stop 'import_gamess%build_orbitals - no memory for special points!'
      end if
      num_special = 0
      !
      !  Loop over grid points. We can run in parallel over Z or, if there is only one value
      !  of Z, over Y, or (if there is only one Y as well) over X.
      !
      !$omp parallel default(none) private(ipz,ipy,ipx,basval,alloc_) &
      !$omp&         shared(grid,gam,coord,dst,mos,num_special,special_pt,r_cut)
      allocate (basval(gam%nbasis),stat=alloc_)
      if (alloc_/=0) stop 'import_gamess%build_orbitals - no (per-thread) memory for basis function values!'
      if (size(coord,dim=4)/=1) then
        !$omp do
        grid_z1: do ipz=1,size(coord,dim=4)
          grid_y1: do ipy=1,size(coord,dim=3)
            grid_x1: do ipx=1,size(coord,dim=2)
              call check_for_oversampling(gam,num_special,special_pt,r_cut,ipx,ipy,ipz,coord(1:3,ipx,ipy,ipz))
              call evaluate_basis_functions(gam,coord(1:3,ipx,ipy,ipz),basval)
              !
              !  Loop over MOs
              ! 
              grid(ipx,ipy,ipz,dst) = matmul(basval,gam%vectors(:,mos))
            end do grid_x1
          end do grid_y1
        end do grid_z1
        !$omp end do
      else if (size(coord,dim=3)/=1) then
        !$omp do
        grid_y2: do ipy=1,size(coord,dim=3)
          grid_x2: do ipx=1,size(coord,dim=2)
            call check_for_oversampling(gam,num_special,special_pt,r_cut,ipx,ipy,1_ik,coord(1:3,ipx,ipy,1))
            call evaluate_basis_functions(gam,coord(1:3,ipx,ipy,1),basval)
            !
            !  Loop over MOs
            ! 
            grid(ipx,ipy,1,dst) = matmul(basval,gam%vectors(:,mos))
          end do grid_x2
        end do grid_y2
        !$omp end do
      else
        !$omp do
        grid_x3: do ipx=1,size(coord,dim=2)
          call check_for_oversampling(gam,num_special,special_pt,r_cut,ipx,1_ik,1_ik,coord(1:3,ipx,1,1))
          call evaluate_basis_functions(gam,coord(1:3,ipx,1,1),basval)
          !
          !  Loop over MOs
          ! 
          grid(ipx,1,1,dst) = matmul(basval,gam%vectors(:,mos))
        end do grid_x3
        !$omp end do
      end if
      if (verbose>=0 .and. num_special>0) then
        !$omp single
        write (out,"(/'Number of grid points requiring oversampling: ',i8/)") num_special
        if (num_special>0 .and. verbose>3) then 
          report_special: do isp=1,num_special
            ipx = special_pt(1,isp)
            ipy = special_pt(2,isp)
            ipz = special_pt(3,isp)
            write (out,"(t10,3(i6,1x),2x,3(1x,f14.7),2x,3(1x,i4))") &
                   special_pt(1:3,isp), coord(1:3,ipx,ipy,ipz), special_pt(4:6,isp)
          end do report_special
        end if
        !$omp end single
      end if
      !
      !  The required degree of oversampling is likely to be (very) different for each volume
      !  element. In this situation, the only reasonable schedule is "dynamic". An alternative 
      !  would be to parallelize individual volume elements; if load balancing becomes an issue,
      !  that's the natural thing to try.
      !
      !$omp do schedule(dynamic,1)
      oversampling_loop: do isp=1,num_special
        ipx = special_pt(1,isp)
        ipy = special_pt(2,isp)
        ipz = special_pt(3,isp)
        call oversample_volume_element(gam,mos,dst,coord(1:3,ipx,ipy,ipz),grid(ipx,ipy,ipz,:),basval,special_pt(4:6,isp))
      end do oversampling_loop
      !$omp end do
      !
      deallocate (basval)
      !$omp end parallel
      if (allocated(special_pt)) deallocate (special_pt)
    end subroutine build_orbitals
    !
    subroutine prepare_for_oversampling(gam,num_special_max,r_cut)
      type(gam_structure), intent(in) :: gam             ! Context to use
      integer(ik), intent(out)        :: num_special_max ! Expected max. number of oversampled points
      real(rk), intent(out)           :: r_cut           ! Cut-off distance for oversampling
      !
      real(rk)    :: diag, vol
      integer(ik) :: per_nucleus
      !
      if (any(gam%dx<=0._rk) .or. .not.oversample) then
        per_nucleus     = 0
        num_special_max = 0
        r_cut           = -1._rk
      else
        diag  = sqrt(sum(gam%dx**2))  ! Diagonal of the grid volume
        vol   = product(gam%dx)       ! Volume
        r_cut = 10._rk * diag         ! Make sure a reasonable, but small number of grid points is included
        per_nucleus     = int(( (4._rk*pi/3._rk)*(r_cut+diag)**3 + vol ) / vol)
        num_special_max = per_nucleus * gam%natoms
      end if
      if (verbose>=2) then
        write (out,"('Expecting up to ',i8,' special points (',i7,' per nucleus)')") num_special_max, per_nucleus
      end if
    end subroutine prepare_for_oversampling
    !
    subroutine check_for_oversampling(gam,num_special,special_pt,r_cut,ipx,ipy,ipz,abs_xyz)
      type(gam_structure), intent(in) :: gam             ! Context to use
      integer(ik), intent(inout)      :: num_special     ! Current number of special points
      integer(ik), intent(inout)      :: special_pt(:,:) ! List of special points
      real(rk), intent(in)            :: r_cut           ! Distance cut-off
      integer(ik), intent(in)         :: ipx, ipy, ipz   ! Grid point indices
      real(rk), intent(in)            :: abs_xyz(:)      ! Coordinates of the grid point to check
      !
      integer(ik)  :: iat
      real(rk)     :: rot_xyz(3), r12
      integer(ik)  :: ord(3)          ! Number of samples
      !
      rot_xyz = matmul(transpose(real(gam%rotmat,kind=kind(abs_xyz))),abs_xyz)
      scan_atoms: do iat=1,gam%natoms
        r12 = sqrt(sum((rot_xyz-real(gam%atoms(iat)%xyz,kind=kind(abs_xyz))/abohr)**2))
        if (r12>r_cut) cycle scan_atoms
        !
        !  Volume element is within the sphere of interest; worth checking sampling
        !
        call get_oversampling_order(ord,gam,abs_xyz)
        !
        !  This volume element needs to be sub-sampled - Update the hit list. 
        !  Since this routine is executing concurrently, we have to excercise 
        !  a bit of care! 
        !
        if (product(ord)>1) then
        !$omp critical
           num_special = num_special + 1
           if (num_special>size(special_pt,dim=2)) then
             stop 'import_gamess%check_for_oversampling - blown the estimate'
           end if
           special_pt(:,num_special) = (/ ipx, ipy, ipz, ord /)
        !$omp end critical
        end if
        !
        !  Nothing more to check for this volume element, bail out
        !
        return
      end do scan_atoms
    end subroutine check_for_oversampling
    !
    !  Evaluate a more accurate representation of the MOs, by oversampling within the volume
    !  element. For each MO, the goal of the exercise is to calculate a number x, which:
    !  a) has the same sign as the integral of the MO over the volume element
    !  b) has a modulus equal to the square root of the square of the MO integrated over the
    !     volume element.
    !  We'll use Gauss-Legendre product quadrature, although this is by no means the perfect
    !  choice.
    !
    subroutine oversample_volume_element(gam,mos,dst,coord,grid,basval,ord)
      type(gam_structure), intent(in) :: gam             ! Context to use
      integer(ik), intent(in)         :: mos(:)          ! Indices of the orbitals to be computed
      integer(ik), intent(in)         :: dst(:)          ! Indices of grid to receive the orbitals
      real(rk), intent(in)            :: coord(:)        ! Coordinates of the centre of the volume element
      complex(rk), intent(inout)      :: grid(:)         ! Data field for the MOs. Index by (dst) to get
                                                         ! the MOs in the order specified by (mos)
                                                         ! On input, contains non-oversampled values
      real(rk)                       :: basval(:)       ! Scratch; enough to hold (gam%nbasis) worth of data
      integer(ik), intent(in)         :: ord(3)          ! Integration order needed to hanlde each Cartesian dimension
                                                         !
      !
      real(rk), allocatable  :: ptx(:), wgx(:)  ! Integration abscissas and weights, X direction
      real(rk), allocatable  :: pty(:), wgy(:)  ! Integration abscissas and weights, Y direction
      real(rk), allocatable  :: ptz(:), wgz(:)  ! Integration abscissas and weights, Z direction
      real(rk), allocatable :: pv(:)           ! MO value at sub-sample grid point
      real(rk), allocatable :: spv(:), s2pv(:) ! Accumulated average MO value and its square
      real(rk)              :: wgt, swgt
      integer(ik)            :: alloc
      integer(ik)            :: ix, iy, iz
      real(rk)               :: cpt(3)
      integer(ik)            :: nmos, imo
      logical, save          :: warn_symmetry = .true.
      !
      nmos = size(mos)
      allocate (ptx(ord(1)),wgx(ord(1)),pty(ord(2)),wgy(ord(2)),ptz(ord(3)),wgz(ord(3)), &
               pv(nmos),spv(nmos),s2pv(nmos),stat=alloc)
      if (alloc/=0) then
        write (out,"('oversample_volume_element: Allocation failed. code = ',i8,' ord = ',3i8,' nmos = ',i8)") &
               alloc, ord, nmos
        stop 'import_gamess%oversample_volume_element - no memory'
      end if
      call MathGetQuadrature('Legendre',order=ord(1),x=ptx,w=wgx)
      call MathGetQuadrature('Legendre',order=ord(2),x=pty,w=wgy)
      call MathGetQuadrature('Legendre',order=ord(3),x=ptz,w=wgz)
      !
      swgt = 0
      spv  = 0
      s2pv = 0
      samples_z: do iz=1,ord(3)
        cpt(3) = coord(3) + 0.5_rk*gam%dx(3) * ptz(iz)
        samples_y: do iy=1,ord(2)
          cpt(2) = coord(2) + 0.5_rk*gam%dx(2) * pty(iy)
          samples_x: do ix=1,ord(1)
            cpt(1) = coord(1) + 0.5_rk*gam%dx(1) * ptx(ix)
            wgt    = wgx(ix) * wgy(iy) * wgz(iz)
            if (verbose>=3) then
              write (out,"(t4,'sub-sample ',3i4,' is at ',3f14.7,' wgt = ',g14.7)") ix, iy, iz, cpt, wgt
            end if
            call evaluate_basis_functions(gam,cpt,basval)
            pv   = matmul(basval,gam%vectors(:,mos))
            spv  = spv  + wgt * pv
            s2pv = s2pv + wgt * pv**2
            swgt = swgt + wgt
          end do samples_x
        end do samples_y
      end do samples_z
      !
      spv  =      spv  / swgt
      s2pv = sqrt(s2pv / swgt)
      if (verbose>=3) then
        write (out,"('Total weight is: ',g25.15,'; points: ',i6)") swgt, product(ord)
        report_mos: do imo=1,nmos
          write (out,"(4x,' mo= ',i4,' old= ',g14.7,1x,g14.7,' ave= ',g14.7,' ave2= ',g14.7)") &
                 imo, grid(dst(imo)), spv(imo), s2pv(imo)
        end do report_mos
      end if
      !
      !  There is a potential problem here: if the grid point lies on a high-symmetry
      !  plane, the average of the orbital may vanish, even though the average of the
      !  square will not. There is nothing we could really do in this situation, so
      !  if we detect a potential inconsistency, we'll issue a warning (just once!)
      !  and use the non-oversampled value at the high-symmetry position. Bummer, but ...
      !
      assign_mos: do imo=1,nmos
        if (abs(abs(spv(imo))-s2pv(imo))>=0.9_rk*max(abs(spv(imo)),s2pv(imo),1e-5_rk)) then
          !$omp flush(warn_symmetry)
          !
          !  The odd structure with an apparently redundant double test is to avoid producing
          !  untidy duplicate output, while at the same time to avoid an unnecessary serialization
          !  point at the critical section. Concurrency could be quite odd ...
          !
          if (verbose>=2 .or. warn_symmetry) then
            !$omp critical
          
            if (verbose>=2 .or. warn_symmetry) then
              write (out,"(/'Unable to consistently oversample orbital at grid point ',3(1x,f14.7))") coord
              write (out,"('Central value: ',g14.7,' average: ',g14.7,' root of average square: ',g14.7/)") &
                     real(grid(dst(imo)),kind=rk), spv(imo), s2pv(imo)
              if (warn_symmetry) then
                write (out,"(/' The most likely reason for this problem is the presence of a nodal plane, which')")
                write (out,"( ' accidentally, or by symmetry, coincides with the exact position of the centre')")
                write (out,"( ' of the volume element. Changing the number of grid points along the affected ')")
                write (out,"( ' direction by 1 may help to resolve this warning. Or possibly not.'/)")
                warn_symmetry = .false.
                !$omp flush(warn_symmetry)
              end if
            end if
            !$omp end critical
          end if
        else
          grid(dst(imo)) = sign(s2pv(imo),spv(imo))
        end if
      end do assign_mos
      !
      deallocate(ptx,wgx,pty,wgy,ptz,wgz)
    end subroutine oversample_volume_element
    !
    subroutine get_oversampling_order(ord,gam,abs_xyz)
      integer(ik), intent(out)                :: ord(:)      ! Required oversampling order along each axis
      type(gam_structure), intent(in), target :: gam         ! Context to use
      real(rk), intent(in)                    :: abs_xyz(:)  ! Coordinates of the grid point to check
      !
      integer(ik)             :: iat
      real(rk)               :: rot_xyz(3)      ! Rotated coordinates of the centre of the volume
      real(rk)               :: r12             ! Distance from current atom to the centre
      real(rk)               :: rmin, rmax      ! Min and max possible distance
      real(rk)               :: diag            ! Diagonal of the grid volume
      real(rk)               :: r0, fwhm        ! Estimates of the maximum position and width
                                                 ! for the primitive
      real(rk)               :: pmax            ! Estimate of the maximum value of the primitive
      type(gam_atom), pointer :: at              ! Atom, what else?
      integer(ik)             :: ish             ! Shell index
      integer(ik)             :: ip, ip_l, ip_h  ! Primitive index and the limits
      integer(ik)             :: lv              ! Power of r associated with this primitive
      real(rk)               :: zet             ! Exponent
      integer(ik)             :: npts_ang        ! Number of angular points we'd need to get current primitive right
      integer(ik)             :: npts_rad        ! Number of radial points we'd need to get current primitive right
      integer(ik)             :: npts            ! Number of points across the diagonal we need for all primitives
      !
      rot_xyz = matmul(transpose(real(gam%rotmat,kind=kind(abs_xyz))),abs_xyz)
      diag    = sqrt(sum(gam%dx**2))
      !
      npts = 0
      scan_atoms: do iat=1,gam%natoms
        at   => gam%atoms(iat)
        r12  = sqrt(sum((rot_xyz-real(at%xyz,kind=kind(abs_xyz))/abohr)**2))
        rmin = max(0.0_rk,r12 - 0.5_rk * diag)
        rmax =            r12 + 0.5_rk * diag
        !
        scan_shells: do ish=1,at%nshell
          lv   = at%sh_l(ish)
          ip_l = at%sh_p(ish)
          ip_h = at%sh_p(ish+1)-1
          scan_primitives: do ip=ip_l,ip_h
            zet = real(at%p_zet(ip),kind=kind(abs_xyz))
            !
            !  For each primitive, estimate the position of the maximum, full width at
            !  half height, and the value at the maximum. We'll use a quadratic approximation
            !  for FWHH - there is no point in being too accurate here!
            !
            r0   = sqrt(0.5_rk*lv/zet)
            fwhm = 1._rk/sqrt(zet)
            pmax = real(at%p_c(ip),kind=kind(abs_xyz)) * exp(-zet*r0**2)
            if (lv>0) pmax = pmax * r0**lv
            !
            !  If primitive is too far away from our volume, we can ignore it for oversampling.
            !  Ditto if it can't contribute anything to the total integral. Note that the
            !  factors "2._rk" and "1e-5_rk" are completely empirical!
            !
            if (r0+2._rk*fwhm<rmin) cycle scan_primitives
            if (pmax<1e-5_rk) cycle scan_primitives
            !
            !  We probably have to worry about describing this guy. Again, the number "1.5_rk" is
            !  purely empirical - adjust to taste.
            !
            npts_rad = ceiling(3.0_rk*diag/fwhm)
            npts_ang = ceiling((2*lv+1)*diag/(twopi*max(r0,rmax)))
            npts = max(npts,npts_rad,npts_ang)
          end do scan_primitives
        end do scan_shells
      end do scan_atoms
      !
      !  Now scale the number of primitives across the diagonal according to the 
      !  linear extent of each dimension
      !
      ord = ceiling( (gam%dx * npts) / diag)
    end subroutine get_oversampling_order
    !
    !  Evaluates values of all basis functions at a grid point.
    !  We do not need to be terribly efficient - just reasonably
    !  careful should do.
    !
    subroutine evaluate_basis_functions(gam,abs_xyz,basval)
      type(gam_structure), target, intent(in) :: gam        ! Context to use
      real(rk), intent(in)                    :: abs_xyz(:) ! Coordinates of the point
      real(rk), intent(out)                   :: basval(:)  ! Values of the basis functions
      !
      integer(ik)             :: iat, ish
      integer(ik)             :: ip1, ip2, ip
      integer(ik)             :: ic1, ic2, ic
      integer(ik)             :: p_bas
      type(gam_atom), pointer :: at
      real(rk)                :: rot_xyz(3)
      real(rk)                :: xyz(3,0:5), r2
      real(rk)                :: ang(0:34)
      real(rk)                :: radial, exparg
!     real(rk)                :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
!     !
!     !  Prepare angular factors of the right kind
!     !
!     ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      !
      !  Rotate grid position into molecular frame
      !
      rot_xyz = matmul(transpose(gam%rotmat_rk),abs_xyz)
      !
      !  Process all centres containing basis functions
      ! 
      p_bas = 1
      atom_loop: do iat=1,gam%natoms
        at => gam%atoms(iat)
        !
        !  Calculate relative position and their powers. We need
        !  those for doing angular parts of the orbitals. Atomic
        !  coordinates are stored in Angstroms - convert to Bohr
        !
        !  Furthermore, we need to rotate the current position into molecular
        !  orientation before calculating anything.
        !
        xyz(:,0) = 1._rk
        xyz(:,1) = rot_xyz - at%xyz_rk/abohr
        xyz(:,2) = xyz(:,1)**2
        xyz(:,3) = xyz(:,1)*xyz(:,2)
        xyz(:,4) = xyz(:,2)**2
        xyz(:,5) = xyz(:,2)*xyz(:,3)
        !
        !  Calculate square distance
        !
        r2 = sum(xyz(:,2))
        !
        !  Precalculate angular parts for all possible momenta
        !
        ang = ang_c_rk*xyz(1,ang_nx)*xyz(2,ang_ny)*xyz(3,ang_nz)
        !
        !  Process all shells on the currently selected centre
        !
        shell_loop: do ish=1,at%nshell
          !
          !  Calculate contracted radial part
          !
          ip1 = at%sh_p(ish)
          ip2 = at%sh_p(ish+1)-1
          !
          radial = 0 ;
          contract_radial: do ip=ip1,ip2
            exparg = -at%p_zet_rk(ip)*r2
            if (exparg<=-max_exp) cycle contract_radial
            radial = radial + at%p_c_rk(ip)*exp(exparg)
          end do contract_radial
          !
          !  Combine angular and radial parts in the output buffer
          !
          ic1 = ang_loc(at%sh_l(ish))
          ic2 = ang_loc(at%sh_l(ish)+1)-1
          combine_angular: do ic=ic1,ic2
            basval(p_bas) = radial * ang(ic)
            p_bas = p_bas + 1
          end do combine_angular
        end do shell_loop
      end do atom_loop
      !
      if (p_bas-1/=gam%nbasis) then
        stop 'import_gamess%evaluate_basis_functions - count error'
      end if
    end subroutine evaluate_basis_functions
    !
    !  Evaluates values of all basis functions and their gradients at a grid point.
    !  This is a modification of evaluate_basis_functions above.
    !
    subroutine evaluate_basis_functions_and_gradients(gam,abs_xyz,basval)
      type(gam_structure), target, intent(in) :: gam          ! Context to use
      real(rk), intent(in)                    :: abs_xyz(:)   ! Coordinates of the point
      real(rk), intent(out)                  :: basval(:,:)  ! First index:
                                                              ! 1 = Values of the basis functions
                                                              ! 2 = X gradient of the basis function, etc
                                                               
      !
      integer(ik)             :: iat, ish
      integer(ik)             :: ip1, ip2, ip
      integer(ik)             :: ic1, ic2, ic
      integer(ik)             :: p_bas
      type(gam_atom), pointer :: at
      real(rk)                :: rot_xyz(3)
      real(rk)                :: xyz(3,-1:6), r2
      real(rk)                :: ang(0:34), angd(3,0:34), angu(3,0:34)
      real(rk)                :: radial, radial2, exparg
!     real(rk)                :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
!     !
!     !  Prepare angular factors of the right kind
!     !
!     ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      !
      !  Rotate grid position into molecular frame
      !
      rot_xyz = matmul(transpose(gam%rotmat_rk),abs_xyz)
      !
      !  Process all centres containing basis functions
      ! 
      p_bas = 1
      atom_loop: do iat=1,gam%natoms
        at => gam%atoms(iat)
        !
        !  Calculate relative position and their powers. We need
        !  those for doing angular parts of the orbitals. Atomic
        !  coordinates are stored in Angstroms - convert to Bohr
        !
        !  Furthermore, we need to rotate the current position into molecular
        !  orientation before calculating anything.
        !
        xyz(:,-1) = 0._rk    ! The value is immaterial, as long as its finite
        xyz(:, 0) = 1._rk
        xyz(:, 1) = rot_xyz - at%xyz_rk/abohr
        xyz(:, 2) = xyz(:,1)**2
        xyz(:, 3) = xyz(:,1)*xyz(:,2)
        xyz(:, 4) = xyz(:,2)**2
        xyz(:, 5) = xyz(:,2)*xyz(:,3)
        xyz(:, 6) = xyz(:,3)*xyz(:,3)
        !
        !  Calculate square distance
        !
        r2 = sum(xyz(:,2))
        !
        !  Precalculate angular parts for all possible momenta
        !
        ang       = ang_c_rk*xyz(1,ang_nx)  *xyz(2,ang_ny)  *xyz(3,ang_nz)
        !  Step-down in power along each direction; include the corresponding power
        angd(1,:) = ang_c_rk*xyz(1,ang_nx-1)*xyz(2,ang_ny)  *xyz(3,ang_nz)  *ang_nx
        angd(2,:) = ang_c_rk*xyz(1,ang_nx)  *xyz(2,ang_ny-1)*xyz(3,ang_nz)  *ang_ny
        angd(3,:) = ang_c_rk*xyz(1,ang_nx)  *xyz(2,ang_ny)  *xyz(3,ang_nz-1)*ang_nz
        !  Step-up in power along each direction; include the overall -2 factor
        angu(1,:) = ang_c_rk*xyz(1,ang_nx+1)*xyz(2,ang_ny)  *xyz(3,ang_nz)  *(-2)
        angu(2,:) = ang_c_rk*xyz(1,ang_nx)  *xyz(2,ang_ny+1)*xyz(3,ang_nz)  *(-2)
        angu(3,:) = ang_c_rk*xyz(1,ang_nx)  *xyz(2,ang_ny)  *xyz(3,ang_nz+1)*(-2)
        !
        !  Process all shells on the currently selected centre
        !
        shell_loop: do ish=1,at%nshell
          !
          !  Calculate contracted radial part
          !
          ip1 = at%sh_p(ish)
          ip2 = at%sh_p(ish+1)-1
          !
          radial = 0 ; radial2 = 0
          contract_radial: do ip=ip1,ip2
            exparg = -at%p_zet_rk(ip)*r2
            if (exparg<=-max_exp) cycle contract_radial
            exparg = at%p_c_rk(ip)*exp(exparg)
            radial  = radial  + exparg
            radial2 = radial2 + exparg * at%p_zet_rk(ip)
          end do contract_radial
          !
          !  Combine angular and radial parts in the output buffer
          !
          ic1 = ang_loc(at%sh_l(ish))
          ic2 = ang_loc(at%sh_l(ish)+1)-1
          combine_angular: do ic=ic1,ic2
            basval(1,p_bas) = radial * ang(ic)
            basval(2,p_bas) = radial * angd(1,ic) + radial2*angu(1,ic)
            basval(3,p_bas) = radial * angd(2,ic) + radial2*angu(2,ic)
            basval(4,p_bas) = radial * angd(3,ic) + radial2*angu(3,ic)
            p_bas = p_bas + 1
          end do combine_angular
        end do shell_loop
      end do atom_loop
      !
      if (p_bas-1/=gam%nbasis) then
        stop 'import_gamess%evaluate_basis_functions_and_gradients - count error'
      end if
    end subroutine evaluate_basis_functions_and_gradients
    !
    !  A very simple Obara-Saika integral package (real operators)
    !
    subroutine os_1e_matrix_real(what,l,r,v)
      type (gam_operator_data), intent(in) :: what   ! Operator to evaluate
      type (gam_structure), intent(in)     :: l      ! Basis set on the left
      type (gam_structure), intent(in)     :: r      ! Basis set on the right
      real(rk), intent(out)                :: v(:,:) ! Buffer for the integrals
      !
      include 'import_gamess_os_1e_matrix_common.f90'
    end subroutine os_1e_matrix_real
    !
    !  Real operators, quad accuracy. 
    !
    subroutine os_1e_matrix_quad(what,l,r,v)
      type (gam_operator_data), intent(in) :: what   ! Operator to evaluate
      type (gam_structure), intent(in)     :: l      ! Basis set on the left
      type (gam_structure), intent(in)     :: r      ! Basis set on the right
      real(xrk), intent(out)               :: v(:,:) ! Buffer for the integrals
      !
      include 'import_gamess_os_1e_matrix_common.f90'
    end subroutine os_1e_matrix_quad
    !
    !  A very simple Obara-Saika integral package (complex operators)
    !
    subroutine os_1e_matrix_complex(what,l,r,v)
      type (gam_operator_data), intent(in) :: what   ! Operator to evaluate
      type (gam_structure), intent(in)     :: l      ! Basis set on the left
      type (gam_structure), intent(in)     :: r      ! Basis set on the right
      complex(rk), intent(out)             :: v(:,:) ! Buffer for the integrals
      !
      include 'import_gamess_os_1e_matrix_common.f90'
    end subroutine os_1e_matrix_complex
    !
    !  Calculate contraction of a primitive integral block (real version)
    !
    subroutine os_1e_contraction_real(what,vb,r_l,l_l,z_l,c_l,r_r,l_r,z_r,c_r)
      type (gam_operator_data), intent(in) :: what    ! Integral to evaluate
      real(rk), intent(out)                :: vb(:,:) ! Place for the integrals
      real(rk), intent(in)                :: r_l(:)  ! Coordinate on the left, in Bohr
      integer(ik), intent(in)              :: l_l     ! "Angular momentum" on the left
      real(rk), intent(in)                :: z_l(:)  ! Primitive exponents on the left
      real(rk), intent(in)                :: c_l(:)  ! Contractions on the left
      real(rk), intent(in)                :: r_r(:)  ! Ditto on the right
      integer(ik), intent(in)              :: l_r     ! 
      real(rk), intent(in)                :: z_r(:)  ! 
      real(rk), intent(in)                :: c_r(:)  ! 
      !
      include 'import_gamess_os_1e_contraction_common.f90'
    end subroutine os_1e_contraction_real
    !
    !  Calculate contraction of a primitive integral block (high-precision version)
    !
    subroutine os_1e_contraction_quad(what,vb,r_l,l_l,z_l,c_l,r_r,l_r,z_r,c_r)
      type (gam_operator_data), intent(in) :: what    ! Integral to evaluate
      real(xrk), intent(out)               :: vb(:,:) ! Place for the integrals
      real(xrk), intent(in)                :: r_l(:)  ! Coordinate on the left, in Bohr
      integer(ik), intent(in)              :: l_l     ! "Angular momentum" on the left
      real(xrk), intent(in)                :: z_l(:)  ! Primitive exponents on the left
      real(xrk), intent(in)                :: c_l(:)  ! Contractions on the left
      real(xrk), intent(in)                :: r_r(:)  ! Ditto on the right
      integer(ik), intent(in)              :: l_r     ! 
      real(xrk), intent(in)                :: z_r(:)  ! 
      real(xrk), intent(in)                :: c_r(:)  ! 
      !
      include 'import_gamess_os_1e_contraction_common.f90'
    end subroutine os_1e_contraction_quad
    !
    !  Calculate contraction of a primitive integral block (complex version)
    !
    subroutine os_1e_contraction_complex(what,vb,r_l,l_l,z_l,c_l,r_r,l_r,z_r,c_r)
      type (gam_operator_data), intent(in) :: what    ! Integral to evaluate
      complex(rk), intent(out)             :: vb(:,:) ! Place for the integrals
      real(rk), intent(in)                :: r_l(:)  ! Coordinate on the left, in Bohr
      integer(ik), intent(in)              :: l_l     ! "Angular momentum" on the left
      real(rk), intent(in)                :: z_l(:)  ! Primitive exponents on the left
      real(rk), intent(in)                :: c_l(:)  ! Contractions on the left
      real(rk), intent(in)                :: r_r(:)  ! Ditto on the right
      integer(ik), intent(in)              :: l_r     ! 
      real(rk), intent(in)                :: z_r(:)  ! 
      real(rk), intent(in)                :: c_r(:)  ! 
      !
      integer(ik)              :: ic_l, p1_l, p2_l, p_l, int_l
      integer(ik)              :: ic_r, p1_r, p2_r, p_r, int_r
      complex(rk), allocatable :: vb_r(:,:,:)  ! Buffer for recursions
      real(rk)                 :: wgt
      integer(ik)              :: n_rec        ! Number of additional recursion intermediated
      integer(ik)              :: alloc
      real(rk)                :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
      !
      !  Prepare angular factors of the right kind
      !
      ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      !
      p1_l = ang_loc(l_l)
      p2_l = ang_loc(l_l+1)-1
      p1_r = ang_loc(l_r)
      p2_r = ang_loc(l_r+1)-1
      !
      select case (what%op_name)
        case default
          !
          !  Assume general 3-centre integral; each unit of angular momentum requires
          !  an extra order in recursion.
          !
          n_rec = l_l + l_r
        case ('PLANEWAVE')
          n_rec = 0
      end select
      allocate (vb_r(0:p2_l,0:p2_r,0:n_rec),stat=alloc)
      if (alloc/=0) then
        write (out,"('import_gamess%os_1e_contraction: allocation failed. Parameters are:',3(1x,i8))") p2_l, p2_r, n_rec
        stop 'import_gamess%os_1e_contraction: Buffer allocation failed'
      end if
      !
      vb = 0._rk
      r_contractions_loop: do ic_r=1,size(z_r)
        l_contractions_loop: do ic_l=1,size(z_l)
          !
          !  The only part of the 1e integral code which is property-specific
          !  is the primitive-integral call below. Everything else is identical
          !  for all properties.
          !
          select case (what%op_name)
            case default
              write (out,"('import_gamess%os_1e_contraction_complex: Property ',a,' is not implemented')") trim(what%op_name)
              stop 'import_gamess%os_1e_contraction_complex: Unimplemented property'
            case ('PLANEWAVE')
              call os_1e_planewave(p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r),real(what%kvec,kind=kind(vb_r)))
          end select
          !
          wgt = c_l(ic_l) * c_r(ic_r)
          r_accumulate: do p_r=p1_r,p2_r
            int_r = p_r - p1_r + 1
            l_accumulate: do p_l=p1_l,p2_l
              int_l = p_l - p1_l + 1
              vb(int_l,int_r) = vb(int_l,int_r) + wgt*ang_c_kind(p_l)*ang_c_kind(p_r)*vb_r(p_l,p_r,0)
            end do l_accumulate
          end do r_accumulate
          !
        end do l_contractions_loop
      end do r_contractions_loop
      !
      deallocate (vb_r,stat=alloc)
      if (alloc/=0) then
        stop 'import_gamess%os_1e_contraction_complex: Buffer deallocation failed'
      end if
    end subroutine os_1e_contraction_complex
    !
    !  Calculate primitive overlaps (real version)
    !
    subroutine os_1e_overlap_real(p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)    :: r_l(:)             ! Centre of the left b.f.
      real(rk), intent(in)    :: z_l                ! Orbital exponent of the left b.f.
      real(rk), intent(in)    :: r_r(:)             ! ditto, for the right b.f.
      real(rk), intent(in)    :: z_r                ! 
      !
      include 'import_gamess_os_1e_overlap_common.f90'
    end subroutine os_1e_overlap_real
    !
    !  Calculate primitive overlaps (high-accuracy version)
    !
    subroutine os_1e_overlap_quad(p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(xrk), intent(out)  :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(xrk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(xrk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(xrk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_overlap_common.f90'
    end subroutine os_1e_overlap_quad
    !
    !  Calculate primitive planewave integrals. This is a simple modification of the
    !  os_1e_overlap routine
    !
    subroutine os_1e_planewave(p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r,kvec)
      integer(ik), intent(in)  :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      complex(rk), intent(out) :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)    :: r_l(:)              ! Centre of the left b.f.
      real(rk), intent(in)    :: z_l                 ! Orbital exponent of the left b.f.
      real(rk), intent(in)    :: r_r(:)              ! ditto, for the right b.f.
      real(rk), intent(in)    :: z_r                 ! 
      real(rk), intent(in)    :: kvec(3)             ! Planewave momentum
      !
      real(rk)    :: zeta, xi, p(3), s00
      complex(rk) :: s
      integer(ik)  :: m_l, m_r, ic, id1, id2, id3
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,xi=xi,zeta=zeta,p=p,s00=s00)
      !
      !  The primitive overlap has an extra term in the case of planevawe
      !
      vb_p(0,0) = s00 * exp(-sum(kvec**2)/(4*zeta)) * exp((0,1)*sum(kvec*p))
      !
      !  Bootstrap: calculate all overlaps on the right, keeping left at zero. 
      !
      right_bootstrap: do m_r=1,p2_r
        find_right: do ic=1,3
          id1 = drop_xyz(m_r,ic)
          if (id1<0) cycle find_right
          s = (p(ic)-r_r(ic)+(0,1)*kvec(ic)/(2.0_rk*zeta)) * vb_p(0,id1)
          id2 = drop_xyz(id1,ic)
          if (id2>=0) s = s + vb_p(0,id2)*ang_nxyz(id1,ic)/(2.0_rk*zeta)
          vb_p(0,m_r) = s
          exit find_right
        end do find_right
      end do right_bootstrap
      !
      !  Now the general case: always recurse on the left. 
      !
      right_loop: do m_r=0,p2_r
        left_recurrence: do m_l=1,p2_l
          find_left: do ic=1,3
            id1 = drop_xyz(m_l,ic)
            if (id1<0) cycle find_left
            s = (p(ic)-r_l(ic)+(0,1)*kvec(ic)/(2.0_rk*zeta)) * vb_p(id1,m_r)
            id2 = drop_xyz(id1,ic)
            if (id2>=0) s = s + vb_p(id2,m_r)*ang_nxyz(id1,ic)/(2.0_rk*zeta)
            id3 = drop_xyz(m_r,ic)
            if (id3>=0) s = s + vb_p(id1,id3)*ang_nxyz(m_r,ic)/(2.0_rk*zeta)
            vb_p(m_l,m_r) = s
            exit find_left
          end do find_left
        end do left_recurrence
      end do right_loop
    end subroutine os_1e_planewave
    !
    !  Calculate primitive dipoles
    !
    subroutine os_1e_dipole_real(id,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: id                  ! Dipole component to evaluate
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)    :: r_l(:)             ! Centre of the left b.f.
      real(rk), intent(in)    :: z_l                ! Orbital exponent of the left b.f.
      real(rk), intent(in)    :: r_r(:)             ! ditto, for the right b.f.
      real(rk), intent(in)    :: z_r                ! 
      !
      include 'import_gamess_os_1e_dipole_common.f90'
    end subroutine os_1e_dipole_real
    !
    !  Calculate primitive dipoles (quad-precision version)
    !
    subroutine os_1e_dipole_quad(id,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: id                  ! Dipole component to evaluate
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(xrk), intent(in)   :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(xrk), intent(out)  :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(xrk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(xrk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(xrk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_dipole_common.f90'
    end subroutine os_1e_dipole_quad
    !
    !  Calculate primitive matrix elements of (d/dx) operator.
    !  These reduce to linear combinations of overlap integrals, so it's extra easy.
    !  See 'Electron current implementation notes.PDF' for the derivation of the expression.
    !
    subroutine os_1e_grad_real(ic,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: ic                  ! Gradient component to evaluate
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(rk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(rk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(rk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_grad_common.f90'
    end subroutine os_1e_grad_real
    !
    subroutine os_1e_grad_quad(ic,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: ic                  ! Gradient component to evaluate
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(xrk), intent(in)   :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(xrk), intent(out)  :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(xrk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(xrk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(xrk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_grad_common.f90'
    end subroutine os_1e_grad_quad
    !
    !  Matrix elements of r_i (d/d r_j) operator
    !  See "integral-notes-rddr.pdf" for the derivation
    !
    subroutine os_1e_rddr_real(ic,jc,vb_s,vb_d,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: ic                  ! Component of the dipole operator
      integer(ik), intent(in) :: jc                  ! Component of the derivative operator
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(rk), intent(in)    :: vb_d(0:p2_l,0:p2_r) ! Buffer containing primitive dipole matrix elements
      real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(rk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(rk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(rk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_rddr_common.f90'
    end subroutine os_1e_rddr_real
    !
    subroutine os_1e_rddr_quad(ic,jc,vb_s,vb_d,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: ic                  ! Component of the dipole operator
      integer(ik), intent(in) :: jc                  ! Component of the derivative operator
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(xrk), intent(in)   :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(xrk), intent(in)   :: vb_d(0:p2_l,0:p2_r) ! Buffer containing primitive dipole matrix elements
      real(xrk), intent(out)  :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(xrk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(xrk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(xrk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_rddr_common.f90'
    end subroutine os_1e_rddr_quad
    !
    !  Calculate primitive kinetic energy integrals (-0.5 \Delta)
    !  We follow eq. A12, A13 of Obara & Saika
    !  This routine is a modification of os_1e_dipole (in case this wasn't clear ;-)
    !
    subroutine os_1e_kinetic_real(vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(rk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(rk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(rk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_kinetic_common.f90'
    end subroutine os_1e_kinetic_real
    !
    !  Quad-precision version of the kinetic integrals
    !
    subroutine os_1e_kinetic_quad(vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
      integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
      real(xrk), intent(in)   :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
      real(xrk), intent(out)  :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)   :: r_l(:)              ! Centre of the left b.f.
      real(xrk), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
      real(xrk), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
      real(xrk), intent(in)   :: z_r                 ! 
      !
      include 'import_gamess_os_1e_kinetic_common.f90'
    end subroutine os_1e_kinetic_quad
    !
    !  Calculate primitive 1e 3c integrals. We essentially follow the treatment
    !  of R. Ahlrichs (PCCP 8, 3072-3077 (2006)), suitably modified for three-
    !  centre integrals. The actual recursion is (of course) identical to OS eq. A19
    !
    !  We are not particularly interested in making the code efficient, so the
    !  order of recursions may be less than optimal.
    !
    subroutine os_1e_3c_real(what,p2_l,p2_r,n_rec,vb_p,r_l,z_l,r_r,z_r)
      type(gam_operator_data), intent(in) :: what             ! Operator parameters
      integer(ik), intent(in)  :: p2_l, p2_r, n_rec           ! Upper dimensions of the vb_p array
      real(rk), intent(out)    :: vb_p(0:p2_l,0:p2_r,0:n_rec) ! Buffer used for recurrences, and the final result
      real(rk), intent(in)    :: r_l(:)                      ! Centre of the left b.f.
      real(rk), intent(in)    :: z_l                         ! Orbital exponent of the left b.f.
      real(rk), intent(in)    :: r_r(:)                      ! ditto, for the right b.f.
      real(rk), intent(in)    :: z_r                         ! 
      !
      include 'import_gamess_os_1e_3c_common.f90'
    end subroutine os_1e_3c_real
    !
    !  Quad-precision version of 3-centre integral recursions
    !
    subroutine os_1e_3c_quad(what,p2_l,p2_r,n_rec,vb_p,r_l,z_l,r_r,z_r)
      type(gam_operator_data), intent(in) :: what             ! Operator parameters
      integer(ik), intent(in)  :: p2_l, p2_r, n_rec           ! Upper dimensions of the vb_p array
      real(xrk), intent(out)   :: vb_p(0:p2_l,0:p2_r,0:n_rec) ! Buffer used for recurrences, and the final result
      real(xrk), intent(in)    :: r_l(:)                      ! Centre of the left b.f.
      real(xrk), intent(in)    :: z_l                         ! Orbital exponent of the left b.f.
      real(xrk), intent(in)    :: r_r(:)                      ! ditto, for the right b.f.
      real(xrk), intent(in)    :: z_r                         ! 
      !
      include 'import_gamess_os_1e_3c_common.f90'
    end subroutine os_1e_3c_quad
    !
    !  4-centre 2-electron integrals
    !  We follow the Ahlrichs' PCCP paper. However, instead of using the transfer relations, 
    !  advocated by that paper, we directly use the general recursion formula (31).
    !  This is less computationally efficient than the Ahlrichs' preferred route, but
    !  avoids the need for extending our data tables in gamess_internal.f90.
    !  If efficiency ever becomes a problem, this can be changed easily.
    !
    subroutine os_2e_batch_real(what,v,xyz,sh_l,sh_ns,sh_np,p_zet,p_c,p_c_max,p_master,p_z_same,cutoff_2e)
      type(gam_operator_data), intent(in) :: what              ! Operator parameters
      real(rk), intent(out)               :: v(:,:,:,:)        ! Buffer for the results
      real(rk), intent(in)               :: xyz(:,:)          ! Coordinates of the basis functions within each batch
      integer(ik), intent(in)             :: sh_l (:)          ! Angular momentum of each shell within a batch is the same
      integer(ik), intent(in)             :: sh_ns(:)          ! Number of shells in each batch
      integer(ik), intent(in)             :: sh_np(:,:)        ! Number of primitives in each shell
      real(rk), intent(in)               :: p_zet   (:,:,:)   ! Primitive exponents in each shell
      real(rk), intent(in)               :: p_c     (:,:,:)   ! Primitive contraction coefficients in each shell
      real(rk), intent(in)               :: p_c_max (:,:,:)   ! Absolute maximum of a contraction coefficient for all instances
                                                               ! of a primitive. This value is only defined for the "master" entry
      logical, intent(in)                 :: p_master(:,:,:)   ! "Master" copy of the exponent
      integer(ik), intent(in)             :: p_z_same(:,:,:,:) ! Next contraction with an identical exponent
                                                               ! First index: 1 = shell index; 2 = contraction within shell
                                                               ! 2nd, 3rd, and 4th indices are as the 1st etc in p_zet/p_c/p_master
      real(rk), intent(in)                :: cutoff_2e         ! Absolute accuracy required in the final integrals;
                                                               ! Negative means retain full accuracy
      !
      include 'import_gamess_os_2e_batch_common.f90'
    end subroutine os_2e_batch_real
    !
    subroutine os_2e_batch_quad(what,v,xyz,sh_l,sh_ns,sh_np,p_zet,p_c,p_c_max,p_master,p_z_same,cutoff_2e)
      type(gam_operator_data), intent(in) :: what              ! Operator parameters
      real(xrk), intent(out)              :: v(:,:,:,:)        ! Buffer for the results
      real(xrk), intent(in)               :: xyz(:,:)          ! Coordinates of the basis functions within each batch
      integer(ik), intent(in)             :: sh_l (:)          ! Angular momentum of each shell within a batch is the same
      integer(ik), intent(in)             :: sh_ns(:)          ! Number of shells in each batch
      integer(ik), intent(in)             :: sh_np(:,:)        ! Number of primitives in each shell
      real(xrk), intent(in)               :: p_zet   (:,:,:)   ! Primitive exponents in each shell
      real(xrk), intent(in)               :: p_c     (:,:,:)   ! Primitive contraction coefficients in each shell
      real(xrk), intent(in)               :: p_c_max (:,:,:)   ! Absolute maximum of a contraction coefficient for all instances
                                                               ! of a primitive. This value is only defined for the "master" entry
      logical, intent(in)                 :: p_master(:,:,:)   ! "Master" copy of the exponent
      integer(ik), intent(in)             :: p_z_same(:,:,:,:) ! Next contraction with an identical exponent
                                                               ! First index: 1 = shell index; 2 = contraction within shell
                                                               ! 2nd, 3rd, and 4th indices are as the 1st etc in p_zet/p_c/p_master
      real(xrk), intent(in)               :: cutoff_2e         ! Absolute accuracy required in the final integrals;
                                                               ! Negative means retain full accuracy
      !
      include 'import_gamess_os_2e_batch_common.f90'
    end subroutine os_2e_batch_quad
    !
    !  2-electron 4-centre primitive integrals
    !
    subroutine os_2e_primitives_real(what,xyz,zet,p2,nrec,vb,w_cut,zero)
      type(gam_operator_data), intent(in) :: what        ! Operator parameters
      real(rk), intent(in)               :: xyz(:,:)    ! Coordinates of the basis functions within each batch
      real(rk), intent(in)                :: zet  (:)    ! Primitive exponents
      integer(ik), intent(in)             :: p2   (:)    ! Ending position of the desired block of integrals
      integer(ik), intent(in)             :: nrec        ! Maximum recursion order
      real(rk), intent(out)               :: vb(0:p2(1),0:p2(2),0:p2(3),0:p2(4),0:nrec)
      real(rk), intent(in)                :: w_cut       ! Absolute cutoff for the integrals
      logical, intent(out)                :: zero        ! True if all integrals are smaller than w_cut
      !
      include 'import_gamess_os_2e_primitives_common.f90'
    end subroutine os_2e_primitives_real
    !
    subroutine os_2e_primitives_quad(what,xyz,zet,p2,nrec,vb,w_cut,zero)
      type(gam_operator_data), intent(in) :: what        ! Operator parameters
      real(xrk), intent(in)               :: xyz(:,:)    ! Coordinates of the basis functions within each batch
      real(xrk), intent(in)               :: zet  (:)    ! Primitive exponents
      integer(ik), intent(in)             :: p2   (:)    ! Ending position of the desired block of integrals
      integer(ik), intent(in)             :: nrec        ! Maximum recursion order
      real(xrk), intent(out)              :: vb(0:p2(1),0:p2(2),0:p2(3),0:p2(4),0:nrec)
      real(xrk), intent(in)               :: w_cut       ! Absolute cutoff for the integrals
      logical, intent(out)                :: zero        ! True if all integrals are smaller than w_cut
      !
      include 'import_gamess_os_2e_primitives_common.f90'
    end subroutine os_2e_primitives_quad
    !
    !  Calculate basic integrals Gn(rho,T) for OS/Ahlrichs recursions.
    !  Unless qualified, references below at to the equation numbers in the Ahlrichs' 2006 PCCP.
    !  Since some of the more complicated Gn(rho,T) versions are defined in terms of simpler Gns,
    !  we may enter this subroutine more than once.
    !
    recursive subroutine os_basic_integral_real(what,rho,t,n_max,gn)
      type(gam_operator_data), intent(in) :: what        ! Operator parameters
      real(rk), intent(in)               :: rho         ! Effective basis exponent
      real(rk), intent(in)               :: t           ! Argument of the basic integral
      integer(ik), intent(in)             :: n_max       ! Maximum order of basic integral needed
      real(rk), intent(out)               :: gn(0:n_max) ! Basic integrals
      !
      include 'import_gamess_os_basic_integral_common.f90'
    end subroutine os_basic_integral_real
    !
    recursive subroutine os_basic_integral_quad(what,rho,t,n_max,gn)
      type(gam_operator_data), intent(in) :: what        ! Operator parameters
      real(xrk), intent(in)               :: rho         ! Effective basis exponent
      real(xrk), intent(in)               :: t           ! Argument of the basic integral
      integer(ik), intent(in)             :: n_max       ! Maximum order of basic integral needed
      real(xrk), intent(out)              :: gn(0:n_max) ! Basic integrals
      !
      include 'import_gamess_os_basic_integral_common.f90'
    end subroutine os_basic_integral_quad
    !
    subroutine os_boys_table_real(n_max,t,fn)
      integer(ik), intent(in) :: n_max       ! Max. desired order of the Boys' function
      real(rk), intent(in)    :: t           ! Argument of the Boys' function
      real(rk), intent(out)   :: fn(0:n_max) ! Table for the function values
      !
      include 'import_gamess_os_boys_table_common.f90'
    end subroutine os_boys_table_real
    !
    subroutine os_boys_table_quad(n_max,t,fn)
      integer(ik), intent(in) :: n_max       ! Max. desired order of the Boys' function
      real(xrk), intent(in)   :: t           ! Argument of the Boys' function
      real(xrk), intent(out)  :: fn(0:n_max) ! Table for the function values
      !
      include 'import_gamess_os_boys_table_common.f90'
    end subroutine os_boys_table_quad
    !
    !  Calculate quantities which appear in (nearly) all Obara-Saika recurrences
    !
    subroutine os_common_primitives_real(r_l,z_l,r_r,z_r,xi,zeta,p,s00,s00a)
      real(rk), intent(in)            :: r_l(:)  ! Centre of the left b.f.
      real(rk), intent(in)            :: z_l     ! Orbital exponent of the left b.f.
      real(rk), intent(in)            :: r_r(:)  ! ditto, for the right b.f.
      real(rk), intent(in)            :: z_r     ! 
      real(rk), intent(out), optional :: xi      ! Eq. 13 of Obara-Saika
      real(rk), intent(out), optional :: zeta    ! Eq. 14
      real(rk), intent(out), optional :: p(:)    ! Eq. 15
      real(rk), intent(out), optional :: s00     ! Eq. 22
      real(rk), intent(out), optional :: s00a    ! Ahlrichs' version of s00 - no prefactor
      !
      include 'import_gamess_os_common_primitives_common.f90'
    end subroutine os_common_primitives_real
    !
    !  Calculate quantities which appear in (nearly) all Obara-Saika recurrences, high-accuracy version
    !
    subroutine os_common_primitives_quad(r_l,z_l,r_r,z_r,xi,zeta,p,s00,s00a)
      real(xrk), intent(in)            :: r_l(:)  ! Centre of the left b.f.
      real(xrk), intent(in)            :: z_l     ! Orbital exponent of the left b.f.
      real(xrk), intent(in)            :: r_r(:)  ! ditto, for the right b.f.
      real(xrk), intent(in)            :: z_r     ! 
      real(xrk), intent(out), optional :: xi      ! Eq. 13 of Obara-Saika
      real(xrk), intent(out), optional :: zeta    ! Eq. 14
      real(xrk), intent(out), optional :: p(:)    ! Eq. 15
      real(xrk), intent(out), optional :: s00     ! Eq. 22
      real(xrk), intent(out), optional :: s00a    ! Ahlrichs' version of s00 - no prefactor
      !
      include 'import_gamess_os_common_primitives_common.f90'
    end subroutine os_common_primitives_quad
  end module import_gamess
!  program xxx
!    use accuracy
!    use import_gamess
!    !
!    type(gam_structure)   :: l, r
!    real(rk), allocatable :: ov(:,:), o3(:,:)
!    integer(ik)           :: il, ir
!    real(rk)              :: op_xyz(3)
!    !
!    call accuracyInitialize
!    call gamess_load_orbitals(file='xxx_left.dat',structure=l)
!    call gamess_load_orbitals(file='xxx_right.dat',structure=r)
!    allocate (ov(l%nbasis,r%nbasis),o3(l%nbasis,r%nbasis))
!    op_xyz = 0
!!   call gamess_1e_integrals('AO OVERLAP',ov,l,r)
!!   call gamess_1e_integrals('AO 3C ONE',o3,l,r,op_xyz=op_xyz)
!!   write (out,"(' RESULTS ')")
!!   lr: do ir=1,r%nbasis
!!     ll: do il=1,l%nbasis
!!       if (abs(o3(il,ir)-ov(il,ir))>4*spacing(ov(il,ir))) then
!!         write (out,"(1x,i4,1x,i4,2x,4g28.18)") il, ir, o3(il,ir)-ov(il,ir), o3(il,ir), ov(il,ir)
!!       end if
!!     end do ll
!!   end do lr
!!   call gamess_print_1e_integrals(o3-ov,l,r)
!!   op_xyz(3) = 1.0371572111_rk
!!   call gamess_1e_integrals('AO 3C 1/R',o3,l,r,op_xyz=op_xyz)
!!   ov = 7*o3
!!   op_xyz(3) =-1.0371572111_rk
!!   call gamess_1e_integrals('AO 3C 1/R',o3,l,r,op_xyz=op_xyz)
!!   ov = ov + 7*o3
!!   call gamess_print_1e_integrals(ov,l,r)
!    call gamess_1e_integrals('AO 3C 1/R',ov,l,r,op_xyz=op_xyz)
!    call gamess_1e_integrals('AO 3C R/(R**2+A**2)',o3,l,r,op_xyz=op_xyz,op_param=(/0.1_rk/))
!    write (out,"('RESULT: 1/R')")
!    call gamess_print_1e_integrals(ov,l,r)
!    write (out,"('RESULT: R/(R**2+0.1**2)')")
!    call gamess_print_1e_integrals(o3,l,r)
!    write (out,"('RESULT: DIFFERENCE')")
!    call gamess_print_1e_integrals(o3-ov,l,r)
!    deallocate (ov,o3)
!    call gamess_destroy(l)
!    call gamess_destroy(r)
!  end program xxx
