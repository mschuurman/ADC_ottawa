!
!  A very simple example of point group symmetry, relevant for strong-field
!  ionization of molecules. This code takes a few bits of information:
!
!  1. Overall point group to work with. All point groups with axes up to 
!     order 8 are implemented; the limitation to 8-th order axes is not
!     fundamental; any finite point group should be trivial to add - see point_group.f90
!  2. Molecular structure (either an .xyz file or a GAMESS checkpoint file)
!  3. Set of laser field polarizations
!  4. (Optionally) list of Dyson orbitals
!
!  See "task" variable below for the description of the actual things which can be done
!
!  It produces list of equivalent polarizations for this point group and
!  molecule.
!
!  Alternatively, we can take a list of structures (XYZ format), and decide
!  which ones are symmetry-equivalent.
!
module symgen
  use accuracy
  use math
  use timer
  use point_group
  use import_gamess
  use atoms
  use molecular_grid
  use printing
  use sort_tools
  !
  implicit none
  private
  public run_symgen
  !
  integer(ik), parameter :: iu_temp       = 22    ! I/O unit for scratch use inside subroutines
  !
  !  User-controller input
  !
  integer(ik)        :: verbose           = 2     ! Verbosity level
  character(len=20)  :: pg_name           = 'C1'  ! Point group name
  real(rk)           :: setting_euler (3) = 0._rk ! Euler angles for setting transformation
  real(rk)           :: setting_origin(3) = 0._rk ! Origin for setting transofrmation
                                                  ! The defauls yield standard setting
  character(len=20)  :: task     = 'generate'     ! Task to perform; could be one of:
                                                  ! 'generate'  - takes a single structure and field orientations
                                                  !               as input; produces the list of equivalent field
                                                  !               orientations as the output. This option can take
                                                  !               either GAMESS checkpoint of XYZ inputs
                                                  ! 'compare'   - takes a list of structures as input; produces the
                                                  !               list of equivalent structures. The structures must
                                                  !               be in an XYZ format, and must be in the same file
                                                  !               (possibly the standard input)
                                                  ! 'transform' - takes a structure, a list of Dyson orbitals, and
                                                  !               two field orientations as input. Finds transformation
                                                  !               matrix mapping the second field direction onto the
                                                  !               first, plus the transformation among the Dyson orbitals.
  character(len=200) :: xyz_file = ' '            ! Name of the .xyz structure to read. If both xyz_file and
                                                  ! gam_file are blank, the structre in .xyz format must follow
                                                  ! the namelist.
  character(len=200) :: gam_file = ' '            ! Name of the GAMESS structure to read
  integer(ik)        :: n_molecule = 1            ! Number of structures; only relevant for task='compare'
  integer(ik)        :: n_field = 0               ! Number of field orientations to consider
                                                  ! Is reset to 2 if task == 'transform'
  character(len=200) :: field_file = ' '          ! Name of file containing field polarization directions;
                                                  ! the first three values on a line define orientation; the
                                                  ! rest of the line will be ignored. Blank value means read
                                                  ! field directions from the standard input
  integer(ik)        :: n_dysons   = 0            ! Number of Dyson orbitals to consider in task='transform'
  character(len=200) :: dysons_file = ' '         ! File containing the list of file names, which in turn contain Dyson
                                                  ! orbitals to examine. Blank mens read from standard input.
                                                  ! Each line contains a single file name, followed by a single orbital index to use
  real(rk)           :: match_eps  = -1_rk        ! Desired accuracy in matching; negative means go for 
                                                  ! nearly-machine accuracy
  integer(ik)        :: grid_nrad  = 400_ik       ! Number of radial points in the numerical grid
  integer(ik)        :: grid_nang  = 302_ik       ! Number of angular points in the numerical grid; 110, 302, or 770
  real(rk)           :: complete_eps = 0.001_rk   ! Thereshold for deciding whether degenerate rotation is complete;
                                                  ! only relevant for task=='transform'
  !
  !  Rest is automatically generated
  !
  type(pt_symop)           :: setting             ! Group setting
  type(pt_group)           :: pg                  ! Point group
  integer(ik)              :: nnuc                ! Number of atoms
  real(rk), allocatable    :: xyzq(:,:,:)         ! Coordinates and nuclear charges of atoms. The last index 
                                                  ! is the structure number.
  real(rk), allocatable    :: fields(:,:)         ! Field polarization directions
  integer(ik), allocatable :: ref_field(:,:)      ! The second index is the polarization direction.
                                                  ! The first index is:
                                                  !  1 = 0 if this field direction is unique
                                                  !      otherwise, the index of the first equivalent field orientation
                                                  !  2 = index of the symmetry operation which transforms this orientation
                                                  !      to the equivalent reference field orientation.
  character(len=200), allocatable  :: dyson_files(:)   ! File names containing Dyson orbitals
  integer(ik), allocatable         :: dyson_indices(:) ! Indices of the orbitals to take from the Dyson files
  type(gam_structure), allocatable :: dysons(:,:)      ! Gamess orbital files containing Dyson orbitals.
                                                       ! The second index corresponds to original (=1) and rotated (=2) structures
  real(rk), allocatable            :: dysons_ovr(:,:,:)! Overlaps between Dyson orbitals;
                                                       ! First index: left Dyson
                                                       ! Second index: right Dyson
                                                       ! Third index:
                                                       !   1 = left and right for the reference orientation
                                                       !   2 = left and right after symmetry transformation
                                                       !   3 = left for the reference; right after the transformation
  !
  namelist /symgen_input/ &
           verbose, &
           pg_name, setting_euler, setting_origin, &
           task, &
           n_molecule, xyz_file, gam_file, &
           n_field, field_file, &
           n_dysons, dysons_file, &
           match_eps, grid_nrad, grid_nang, &
           complete_eps
  !
  contains
  !
  ! Read simple .xyz structure input
  !
  subroutine read_xyz_file
    integer(ik)       :: inp     ! I/O channel used for reading the structure data
    integer(ik)       :: iat     ! Atom number
    integer(ik)       :: imol    ! Molecule number
    integer(ik)       :: tnuc    ! Temp for the number of atoms
    character(len=20) :: label
    !
    !  If the input file name was provided in the namelist, use it.
    !  Otherwise, continue reading the in-line input.
    !
    if (xyz_file/=' ') then
      inp = iu_temp
      open(inp,file=trim(xyz_file),status='old',action='read',form='formatted')
    else
      inp = input
    end if
    !
    !  Read the XYZ file. The first line contains the number of atoms.
    !
    nnuc = 0
    read_molecules: do imol=1,n_molecule
      read(inp,*) tnuc
      if (imol==1) then
        !
        !  This is the first structure in the list; allocate memory
        !
        nnuc = tnuc
        if (nnuc<=0) stop 'symgen%read_xyz_file - structure must have _some_ atoms!'
        allocate(xyzq(4,nnuc,n_molecule))
      else
        if (tnuc/=nnuc) then
          write (out,"('Expecting ',i5,' atoms per structure; however structure ',i5,' has ',i5)") nnuc, imol, tnuc
          stop 'symgen%read_xyz_file - inconsistent atom counts!'
        end if
      end if
      !
      !  Skip the XYZ comment line
      !
      read(inp,"(a)")
      !
      ! Read the atoms. Coordinates are expected to be in Angstroms,
      ! do that we'll do the conversion later.
      !
      read_atoms: do iat=1,nnuc
        read(inp,*) label, xyzq(1:3,iat,imol)
        xyzq(1:3,iat,imol) = xyzq(1:3,iat,imol) / abohr
        xyzq(  4,iat,imol) = AtomNuclearCharge(label)
      end do read_atoms
    end do read_molecules
    !
    !  If the structure was coming from an external file, close it
    !
    if (xyz_file/=' ') then
      close(inp,status='keep')
    end if
  end subroutine read_xyz_file
  !
  subroutine read_gam_file
    type(gam_structure) :: gam
    !
    call gamess_load_orbitals(file=gam_file,structure=gam)
    call gamess_report_nuclei(nnuc,structure=gam)
    allocate (xyzq(4,nnuc,1))
    call gamess_report_nuclei(nnuc,xyzq(:,:,1),structure=gam,true_charge=.true.)
    call gamess_destroy(gam)
  end subroutine read_gam_file
  !
  subroutine report_structures
    integer(ik)      :: iat, imol
    character(len=7) :: lab
    !
    print_molecule: do imol=1,n_molecule
      write (out,"()")
      write (out,"(t10,'structure = ',i5)") imol
      write (out,"(      t16,a36,t56,a36)")            'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t16,a36,t56,a36)")            '------------------    ', '----------------------    '
      write (out,"(t10,a5,t16,3a12,t56,3a12)")  'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,nnuc
        lab = AtomElementSymbol(xyzq(4,iat,imol))
        write (out,"(t2,a6,t10,f5.2,t16,3f12.5,t56,3f12.5)") lab, xyzq(4,iat,imol), xyzq(1:3,iat,imol), xyzq(1:3,iat,imol)*abohr
      end do print_atoms
      write (out,"()")
    end do print_molecule
  end subroutine report_structures
  !
  subroutine read_field_directions
    integer(ik) :: inp     ! I/O channel used for reading the structure data
    integer(ik) :: idir
    !
    !  If the input file name was provided in the namelist, use it.
    !  Otherwise, continue reading the in-line input.
    !
    if (field_file/=' ') then
      inp = iu_temp
      open(inp,file=trim(field_file),status='old',action='read',form='formatted')
    else
      inp = input
    end if
    !
    if (n_field<=0) stop 'symgen%read_field_directions - we do need some!'
    allocate(fields(3,n_field),ref_field(2,n_field))
    ref_field = 0
    !
    !  Read field polarizations. We'll normalize them as we go.
    !
    read_directions: do idir=1,n_field
      read(inp,*) fields(:,idir)
      fields(:,idir) = fields(:,idir) / sqrt(sum(fields(:,idir)**2))
      write (out,"('Field polarization ',i0,' is ',3(1x,f14.10))") idir, fields(:,idir)
    end do read_directions
    !
    !  If the structure was coming from an external file, close it
    !
    if (field_file/=' ') then
      close(inp,status='keep')
    end if
  end subroutine read_field_directions
  !
  subroutine transform_molecule(op,from,to)
    type(pt_symop), intent(in) :: op        ! Operation to apply
    real(rk), intent(in)       :: from(:,:) ! XYZQ molecular coordinates to transform
    real(rk), intent(out)      :: to  (:,:) ! Result of the transformation
    !
    integer(ik) :: iat
    !
    if (size(from,dim=1)/=4 .or. size(to,dim=1)/=4 .or. size(from,dim=2)/=size(to,dim=2)) then
      stop 'symgen%transform_molecule - bad from/to arrays'
    end if
    !
    transform: do iat=1,size(from,dim=2)
      to(1:3,iat) = pg_symop_apply(op,from(1:3,iat))
      to(  4,iat) = from(4,iat)
    end do transform
  end subroutine transform_molecule
  !
  function same_molecule(a,b) result(same)
    real(rk), intent(in) :: a(:,:), b(:,:) ! Molecules to compare
    logical              :: same           ! Molecules are the same, up to rearranging the order of atoms
    !
    real(rk)    :: eps
    integer(ik) :: ia, ib               ! Atom indices
    logical     :: taken(size(b,dim=2)) ! True if atom has been matched in "b"
    !
    if (size(a,dim=1)/=4 .or. size(b,dim=1)/=4 .or. size(a,dim=2)/=size(b,dim=2)) then
      stop 'symgen%same_molecule - bad a/b arrays'
    end if
    !
    same = .false.
    taken = .false.
    if (match_eps>0) then
      eps = match_eps
    else
      eps = 100._rk * spacing(maxval(abs(a(1:3,:))))
    end if
    scan_a: do ia=1,size(a,dim=2)
      scan_b: do ib=1,size(b,dim=2)
        if (taken(ib)) cycle scan_b
        if (a(4,ia)/=b(4,ib)) cycle scan_b
        if (any(abs(a(1:3,ia)-b(1:3,ib))>eps)) cycle scan_b
        taken(ib) = .true.
        cycle scan_a
      end do scan_b
      !
      ! If we ended up here, then no match was found for atom iat; declare failure
      !
      return
    end do scan_a
    !
    ! If we end up here, all atoms were successfully matched.
    !
    same = .true.
  end function same_molecule
  !
  function same_field(a,b) result(same)
    real(rk), intent(in) :: a(:), b(:) ! Field polarizations to compare
    logical              :: same
    real(rk)             :: eps
    !
    if (match_eps>0) then
      eps = match_eps
    else 
      eps = spacing(100._rk)
    end if
    same = all( abs(a-b)<spacing(100._rk) )
  end function same_field
  !
  subroutine scan_field_directions
    integer(ik) :: idir              ! Field polarization direction we are lookig at
    integer(ik) :: iref              ! Field polarization direction we are comparing to
    integer(ik) :: iop, iox, iot     ! Symmetry operation
    real(rk)    :: op_xyzq(4,nnuc)   ! Symmetry-transformed molecule
    real(rk)    :: op_field(3)       ! Symmetry-transformed field direction
    logical     :: good_ops(pg%nops) ! Operations which transform molecule into itself
    !
    find_good_ops: do iop=1,pg%nops
      call transform_molecule(op=pg%ops(iop),from=xyzq(:,:,1),to=op_xyzq)
      good_ops(iop) = same_molecule(xyzq(:,:,1),op_xyzq)
    end do find_good_ops
    if (verbose>=1) then
      print_good_ops: do iop=1,pg%nops,5
        iox = min(iop+4,pg%nops)
        write (out,"('    operation:',10(1x,i16))") (iot,                   iot=iop,iox)
        write (out,"('        label:',10(1x,a16))") (trim(pg%ops(iot)%name),iot=iop,iox)
        write (out,"(' molecular op:',10(1x,l16))") (good_ops(iot),         iot=iop,iox)
        write (out,"()")
      end do print_good_ops
      write (out,"('Number of operations preserving molecular symmetry: ',i4)") count (good_ops)
      write (out,"('  Number of operations breaking molecular symmetry: ',i4)") count (.not.good_ops)
    end if
    !
    scan_directions: do idir=1,n_field
      ref_field(1,idir) = 0 ! Provisionally unique
      scan_previous_directions: do iref=1,idir-1
        scan_operations: do iop=1,pg%nops
          if (.not.good_ops(iop)) cycle scan_operations
          !
          op_field = pg_symop_apply(pg%ops(iop),fields(:,idir))
          if (.not.same_field(op_field,fields(:,iref))) cycle scan_operations
          !
          !  Remember reference and the symop and the transforming operation
          !
          ref_field(:,idir) = (/ iref, iop /)
          !
          !  Field polarization (idir) is equivalent to polarization (iref)
          !
          write (out,"('[SUM] Field polarization ',i4,' is equivalent to ',i4)") idir, iref
          if (verbose>0) then
            call pg_symop_print(label='     Transformation symmetry element',op=pg%ops(iop))
          end if
          cycle scan_directions
        end do scan_operations
      end do scan_previous_directions
      !
      !  We only get here if this field polarization direction is unique
      !
      write (out,"('[SUM] Field polarization ',i4,' is unique.')") idir
    end do scan_directions
    
  end subroutine scan_field_directions
  !
  subroutine compare_structures
    integer(ik) :: imol              ! Structure we are looking at
    integer(ik) :: iref              ! Structure we are comparing to
    integer(ik) :: iop               ! Symmetry operation
    real(rk)    :: op_xyzq(4,nnuc)   ! Symmetry-transformed molecule
    !
    scan_structures: do imol=1,n_molecule
      scan_refs: do iref=1,imol-1
        scan_operations: do iop=1,pg%nops
          call transform_molecule(op=pg%ops(iop),from=xyzq(:,:,imol),to=op_xyzq)
          if (.not.same_molecule(xyzq(:,:,iref),op_xyzq)) cycle scan_operations
          !
          !  Structure (imol) is equivalent to structure (iref)
          !
          write (out,"('[SUM] Structure ',i4,' is equivalent to ',i4)") imol, iref
          if (verbose>0) then
            call pg_symop_print(label='     Transformation symmetry element',op=pg%ops(iop))
          end if
          cycle scan_structures
        end do scan_operations
      end do scan_refs
      !
      !  We only get here if this structure is unique
      !
      write (out,"('[SUM] Structure ',i4,' is unique.')") imol
    end do scan_structures
  end subroutine compare_structures
  !
  !  Load Dyson orbitals for which we want to build the transformation matrix
  !
  subroutine load_dysons
    integer(ik) :: idys, inp, ios
    real(ark)   :: rot(3,3)   ! Rotation matrix
    !
    !  Figure out rotation matrices. Note that our symops rotate the molecule;
    !  in gamess_load_orbitals, we rotate the world, but keep the molecule fixed
    !
    !  Rotation matrix -must- be applied at the time structure is loaded: it is
    !  ultimately stored in the structure descriptor, and can't be changed inside
    !  a parallel loop.
    !
    if (ref_field(1,2)/=1) stop 'symgen%load_dysons - no symmetry relation between fields!'
    if (any(pg%ops(ref_field(2,2))%t/=0._rk)) stop 'symgen%load_dysons - translations are not allowed!'
    rot = transpose(pg%ops(ref_field(2,2))%r)
    !
    allocate (dyson_files(n_dysons),dyson_indices(n_dysons),dysons(n_dysons,2))
    !
    inp = input
    if (field_file/=' ') then
      inp = iu_temp
      open(inp,file=trim(dysons_file),status='old',action='read',form='formatted')
    end if
    !
    load_dyson_files: do idys=1,n_dysons
      read (inp,*,iostat=ios) dyson_files(idys), dyson_indices(idys)
      if (ios/=0) then
        write (out,"('Error ',i0,' reading name of Dyson orbital #',i0)") ios, idys
        stop 'symgen%load_dysons - bad Dyson name input'
      end if
      write (out,"(/'Loading orbital ',i0,' from ',a,' (no rotation)')") idys, trim(dyson_files(idys))
      call gamess_load_orbitals(file=trim(dyson_files(idys)),structure=dysons(idys,1))
      !
      write (out,"(/'Loading orbital ',i0,' from ',a,' (rotated)')") idys, trim(dyson_files(idys))
      call gamess_load_orbitals(file=trim(dyson_files(idys)),structure=dysons(idys,2),rot=rot)
    end do load_dyson_files
    !
    if (dysons_file/=' ') then
      close(inp,status='keep')
    end if
  end subroutine load_dysons
  !
  subroutine calculate_dysons_overlap_matrix
    integer(ik)           :: iat, nbatch, ib, npts, idys, il, ir
    character (len=10)    :: atom_labels(nnuc)
    type(mol_grid)        :: grid
    real(rk), allocatable :: lcl(:,:,:)               ! Per-thread overlap buffers; matches dysons_ovr
                                                      ! 1 = rot_0-rot_0; 2 = rot_1-rot_1; 3 = rot_0-rot_1
    real(rk), pointer     :: xyzw(:,:)                ! Integration grid points and weights
    complex(rk), allocatable :: orbs(:,:,:,:,:)       ! Orbitals at grid points. The last index is: 1 = rot_0; 2 = rot_1
                                                      ! The second and the third indices are dummies, and are equal to 1
    !
    !  Prepare numerical integration grid
    !
    fill_labels: do iat=1,nnuc
      atom_labels(iat) = AtomElementSymbol(xyzq(4,iat,1))
    end do fill_labels
    call GridInitialize(grid,grid_nrad,grid_nang,xyzq(1:3,:,1),atom_labels)
    call GridPointsBatch(grid,'Batches count',count=nbatch)
    write (out,"(/'Numerical integration grid contains ',i0,' point batches'/)") nbatch
    !
    allocate (dysons_ovr(n_dysons,n_dysons,3))
    dysons_ovr = 0
    !
    !$omp parallel default(none) &
    !$omp& shared(n_dysons,nbatch,grid,dyson_indices,dysons,dysons_ovr) &
    !$omp& private(xyzw,lcl,npts,ib,orbs,idys)
    nullify(xyzw)
    allocate (lcl(n_dysons,n_dysons,3))
    lcl = 0
    !
    !$omp do schedule(dynamic,1)
    point_batches: do ib=1,nbatch
      call GridPointsBatch(grid,'Next batch',xyzw=xyzw)
      npts = size(xyzw,dim=2)
      if (allocated(orbs)) then
        if (size(orbs,dim=1)/=npts) deallocate (orbs)
      end if
      if (.not.allocated(orbs)) allocate (orbs(npts,1,1,n_dysons,2))
      !
      fill_dysons: do idys=1,n_dysons
        !
        !  We need to evaluate each orbital twice: once before the rotation, then after the rotation
        !
        call gamess_load_orbitals(structure=dysons(idys,1),mos=dyson_indices(idys:idys),dst=(/1_ik/), &
                                  coord=reshape(xyzw(1:3,:),(/3_ik,npts,1_ik,1_ik/)),grid=orbs(:,:,:,idys,1:1))
        call gamess_load_orbitals(structure=dysons(idys,2),mos=dyson_indices(idys:idys),dst=(/1_ik/), &
                                  coord=reshape(xyzw(1:3,:),(/3_ik,npts,1_ik,1_ik/)),grid=orbs(:,:,:,idys,2:2))
        !
        !  Accumulate the integrals. Note that Dysons are guaranteed to be real, so that even
        !  though orbs() is a complex array, it is safe to ignore the imaginary part.
        !
      end do fill_dysons
      update_overlaps_r: do ir=1,n_dysons
        update_overlaps_l: do il=1,n_dysons
          lcl (il,ir,1) = lcl(il,ir,1) + sum(xyzw(4,:)*real(orbs(:,1,1,il,1),kind=rk)*real(orbs(:,1,1,ir,1),kind=rk))
          lcl (il,ir,2) = lcl(il,ir,2) + sum(xyzw(4,:)*real(orbs(:,1,1,il,2),kind=rk)*real(orbs(:,1,1,ir,2),kind=rk))
          lcl (il,ir,3) = lcl(il,ir,3) + sum(xyzw(4,:)*real(orbs(:,1,1,il,1),kind=rk)*real(orbs(:,1,1,ir,2),kind=rk))
        end do update_overlaps_l
      end do update_overlaps_r
    end do point_batches
    !$omp end do nowait
    !
    !$omp critical
    dysons_ovr = dysons_ovr + lcl
    !$omp end critical
    !
    if (associated(xyzw)) deallocate (xyzw)
    if (allocated (orbs)) deallocate (orbs)
    if (allocated (lcl))  deallocate (lcl)
    !$omp end parallel
    !
    call GridDestroy(grid)
  end subroutine calculate_dysons_overlap_matrix
  !
  subroutine report_dysons_transformation_matrix
    integer(ik) :: ir, il, nl
    real(rk)    :: dys_nl, dys_nr   ! Dyson norm before and after the transformation
    real(rk)    :: wgt_cum          ! Cumulative weight of the transformation
    integer(ik) :: weight(n_dysons) ! Ordering of transformation coefficients magnitude
    !
    if (verbose>=1) then
      write (out,"(/t5,'Dyson orbitals overlap before transformation:'/)") 
      call print_matrix(dysons_ovr(:,:,1),9_ik,"f9.6")
      write (out,"(/t5,'Dyson orbitals overlap after transformation:'/)") 
      call print_matrix(dysons_ovr(:,:,2),9_ik,"f9.6")
      write (out,"(/t5,'Dyson orbitals overlap original-transformed:'/)") 
      call print_matrix(dysons_ovr(:,:,3),9_ik,"f9.6")
    end if
    !
    !  Renormalize cross-orbital overlaps
    !
    renormalize_right: do ir=1,n_dysons
      dys_nr = dysons_ovr(ir,ir,2)
      renormalize_left: do il=1,n_dysons
        dys_nl = dysons_ovr(il,il,1)
        dysons_ovr(il,ir,3) = dysons_ovr(il,ir,3) / sqrt(dys_nr*dys_nl)
      end do renormalize_left
    end do renormalize_right
    !
    if (verbose>=1) then
      write (out,"(/t5,'Original-transformed overlap, normalized to unity:'/)") 
      call print_matrix(dysons_ovr(:,:,3),9_ik,"f9.6")
    end if
    !
    !  Try to build per-orbital transformation. 
    !
    transformed_orbital: do ir=1,n_dysons
      dys_nl = dysons_ovr(ir,ir,1)
      dys_nr = dysons_ovr(ir,ir,2)
      write (out,"(/'Dyson orbital ',i0,' norm = ',f12.8,' ; Change upon rotation = ',g12.5)") ir, dys_nl, dys_nr-dys_nl
      if (abs(dys_nr-dys_nl)>complete_eps) stop 'symgen%report_dysons_transformation_matrix - bad norm'
      !
      call order_keys(-abs(dysons_ovr(:,ir,3)),weight)  ! Give me the decreasing order
      !
      !  Choose reference orbitals which give me the desired Dyson
      !
      wgt_cum = 0 ; nl = 0
      choose_references: do il=1,n_dysons
        wgt_cum = wgt_cum + dysons_ovr(weight(il),ir,3)**2
        if (wgt_cum>=1.0_rk-complete_eps) then
          nl = il
          exit choose_references
        end if
      end do choose_references
      if (nl<=0) then
        write (out,"('Can''t find complete representation of Dyson ',i0,'. Only got to ',f12.8)") ir, wgt_cum
        stop 'symgen%report_dysons_transformation_matrix - bad decomposition'
      end if
      write (out,advance='no',fmt="('[DEC] Dyson ',i0,' = ')") ir
      print_references: do il=1,nl
        write (out,advance='no',fmt="(1x,sp,f12.8,' x ',ss,i0)") dysons_ovr(weight(il),ir,3), weight(il)
      end do print_references
      write (out,"()")
    end do transformed_orbital
  end subroutine report_dysons_transformation_matrix
  !
! subroutine clone_vector
!   integer(ik) :: iop
!   !
!   write (out,"('++++++++++++++')")
!   write (out,"(i0)") 4*pg%nops
!   do iop=1,pg%nops
!     write (out,"(i2,3(1x,f16.12))") 1, pg_symop_apply(pg%ops(iop),(/ 2.1_rk, 3.2_rk, 4.3_rk /))
!     write (out,"(i2,3(1x,f16.12))") 2, pg_symop_apply(pg%ops(iop),(/ 5.0_rk,-6.0_rk, 7.0_rk /))
!     write (out,"(i2,3(1x,f16.12))") 3, pg_symop_apply(pg%ops(iop),(/ 0.7_rk, 0.9_rk, 1.1_rk /))
!     write (out,"(i2,3(1x,f16.12))") 4, pg_symop_apply(pg%ops(iop),(/ 3.2_rk,-3.6_rk, 4.0_rk /))
!   end do
!   write (out,"('==============')")
! end subroutine clone_vector
  !
  subroutine run_symgen
    real(rk)          :: dummy
    !
    call TimerStart('Symgen')
    call accuracyInitialize
    dummy = MathFactorial(80_ik)
    !
    write (out,"('Default integer kind = ',i5)") ik
    write (out,"('Default real kind    = ',i5)") rk
    !
    read (input,nml=symgen_input)
    write (out,"(' ========= ')")
    write (out,nml=symgen_input)
    write (out,"(' ========= ')")
    !
    call MathRotationMatrix(euler_angles=setting_euler,mat=setting%r)
    setting%r    = transpose(setting%r) ! We rotate the "object" (in this case, the coordinate system)
    setting%t    = setting_origin
    setting%name = ' '
    !
    call pg_initialize(pg_name,pg,setting=setting)
    !
    select case (task)
      case default
        stop 'symgen%run_symgen - bad task'
      case ('generate')
        n_molecule = 1 ! 1 molecule is the only sensible choice for this task
        if (gam_file/=' ') then
          call read_gam_file
        else
          call read_xyz_file
        end if
        call report_structures
        call read_field_directions
        call scan_field_directions
      case ('compare')
        call read_xyz_file
        call report_structures
        call compare_structures
      case ('transform')
        n_molecule = 1 ! 1 molecule is the only sensible choice for this task
        n_field    = 2 ! Must have exactly two field orientations for this task
        if (gam_file/=' ') then
          call read_gam_file
        else
          call read_xyz_file
        end if
        call report_structures
        call read_field_directions
        call scan_field_directions  ! Fills out ref_field, which we need to process Dysons
        call load_dysons
        call calculate_dysons_overlap_matrix
        call report_dysons_transformation_matrix
    end select
    !
    call TimerStop('Symgen')
    call TimerReport
  end subroutine run_symgen
end module symgen
!
program main
  use symgen

  call run_symgen
end program main
