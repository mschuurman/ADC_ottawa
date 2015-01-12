!
!  Revision history: 2011 Sep 15 - The first implementation using correct formulae
!                    2013 Feb 04 - Separate out conversion routine, which does not
!                                  depend on a particular grid implementation
!                    2013 Jul 30 - Fixed a stupid bug in evaluating <R> expectations,
!                                  which caused the choice of projectors to be origin-
!                                  dependent
!
!  A very simple implementation of effective core potentials. We implement short-range
!  ECPs as a sum of level-shift operators, centered on each of the ECP-carrying atoms:
!
!   v_{ECP} = Sum |phi_a> V_a <phi_a|
!              a
!
!  where orbitals phi_a and shift constants V_a are determined from the "standard" 
!  GAMESS ECPs in the angular-projector form:
!
!   v_{ECP} = Sum Sum |Y_{lm}> (v_l(r) - v_{L+1)(r)) <Y_{lm}|
!              l   m
!
!  The conversion is carried out by solving the standard generalized eigenvalue
!  problem for the eigenvalues and eigenvectors of the shift operator:
!
!     G = S^{-1/2} V S^{-1/2}
!     G D = D E
!     C = S^{-1/2} D
!
!  where V are the matrix elements of v_{ECP}, S is the overlap matrix of the
!  auxiliary basis, G are the matrix elements of v_{ECP} in the Lowdin-orthogonalized
!  basis, D and E are the eigenvalues and eigenvectors of v_{ECP} in the same basis,
!  and C are the eigenvectors in the original aux. basis.
!
!  The auxiliary basis used to represent the Hilbert space for V_{ECP} is constructed 
!  by decontracting the AO basis on the ECP-carrying centre.
!
!  For historical reasons, this module is implemented as an "after-market" addition
!  to import_gamess: it will convert the projector-form ECP loaded by import_gamess
!  into the level-shift form.
!
!  The shift-operator projectors are intended for use with FieldECPProject and
!  FieldECPApply in multigrid.f90
!
  module ecp_convert_gamess
    use accuracy
    use timer
    use gamess_internal
    use import_gamess
    use lebedev
    use math
    use fields
    use block_diag
    use matrix_tools
    use printing

    implicit none

    private
    public ecp_molecule, ecp_atom
    public ecp_convert, ecp_destroy
    public ecp_evaluate_matrix_elements
    !
    type ecp_atom
      real(rk)              :: rmax            ! Projectors vanish this far away from the atom
      real(rk), allocatable :: vshift(:)       ! Per-projector portential shifts
                                               ! The number of projectors is projectors%nvectors
      real(rk), allocatable :: rv(:)           ! Expectation value of <R> for each projector; measures
                                               ! how tight this projector is
      type(gam_structure)   :: projectors      ! Projectors, basis set representation
      real(rk), allocatable :: grid(:,:,:,:)   ! Projectors, instantiated on grid. In memory-constrained
                                               ! cases, it is OK to deallocate grid at any time - it will
                                               ! be re-created on the next call to ecp_apply.
    end type ecp_atom
    !
    type ecp_molecule
      logical                     :: active = .false.
      integer(ik)                 :: necps    ! Number of ECP-carrying centres
      type(ecp_atom), allocatable :: ecps(:)  ! Per-atom ECPs in the shift-operator form
    end type ecp_molecule
    !
    !  Some global constants ...
    !
    integer(ik), parameter :: verbose     = 0
    integer(ik), parameter :: unit_report = 47       ! A semi-random unit used for reporting
    real(rk), parameter    :: ecp_vcut    = 1e-5_rk  ! Ignore smaller potential contributions
                                                     ! Prior to 2013 Jan 02, the cutoff was 1e-6, which is
                                                     ! unreasonably tight considering quality of our grids.
    real(rk), parameter    :: zero_frac   = 1e-9_rk  ! Set elements smaller than zero_frac*max to zero
                                                     ! for all matrix elements
    !
    interface clean_matrix
      module procedure clean_matrix_real
      module procedure clean_matrix_complex
    end interface clean_matrix

    contains

    subroutine ecp_convert(gam,ecp,eps_min,eps_max,grid_delta,report_file,rot)
      type(gam_structure), intent(in)        :: gam         ! Gamess data structure, containing basis and an ECP
      type(ecp_molecule), intent(out)        :: ecp         ! ECP converted into a level-shift form
      real(rk), intent(in), optional         :: eps_min     ! Threshold for eliminating small absolute shift eigenvalues.
                                                            ! The default is to keep all eigenvalues, no matter how small
                                                            ! The default is numerically harmless, but wasteful
      real(rk), intent(in), optional         :: eps_max     ! Threshold for eliminating very large positive shift eigenvalues.
                                                            ! The default is to keep all eigenvalues, no matter how large.
                                                            ! Keeping eigenvalues much greater than d**-2 where d is grid 
                                                            ! spacing is likely to cause numerical trouble.
                                                            ! One could also cut-off steep projectors by evaluating their
                                                            ! <r^2>, and eliminating those which can't be resolved on the
                                                            ! grid. I'll leave this one as an excercise for the future -
                                                            ! visual inspection is likely good enough for the time being.
      real(rk), intent(in), optional         :: grid_delta  ! Characteristic grid spacing of the target grid; projectors
                                                            ! which can't be resolved on this grid should be dropped.
                                                            ! About twice the grid spacing is a good rule of thumb for this
                                                            ! parameter.
      character(len=*), intent(in), optional :: report_file ! Produce a pseudo-GAMESS checkpoint file, which could be useful
                                                            ! for visual inspection of the level-shift projectors
      real(rk), intent(in), optional         :: rot(:,:)    ! Rotation matrix
      !
      integer(ik) :: ia, ie, alloc
      real(rk)    :: eps_low, eps_high, eps_grid
      logical     :: do_report
      real(rk)    :: rotmat(3,3)
      !
      !  Deal with optional parameters
      !
      eps_low   = -1._rk
      eps_high  = huge(1._rk)
      eps_grid  = 0._rk
      do_report = .false.
      if (present(eps_min)) eps_low  = eps_min
      if (present(eps_max)) eps_high = eps_max
      if (present(grid_delta)) eps_grid = grid_delta
      if (present(report_file)) then
        if (report_file/=' ') then
          do_report = .true.
          open (unit=unit_report,form='formatted',status='replace',file=trim(report_file))
        end if
      end if
      if (present(rot)) then
        rotmat = rot
      else
        call MathSetUnitMatrix(rotmat)
      end if
      !
      ecp%necps = count(gam%atoms(1:gam%natoms)%ecp_nterms>0)
      allocate (ecp%ecps(ecp%necps),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating table of ',i6,' ECPs')") alloc, ecp%necps
        stop 'ecp_gamess%ecp_convert - allocate failed'
      end if
      !
      ie = 0
      scan_atoms_for_ecps: do ia=1,gam%natoms
        if (gam%atoms(ia)%ecp_nterms<=0) cycle scan_atoms_for_ecps
        ie = ie + 1
        !
        !  Since the cost of the ECP conversion is trivial compared to the run-time costs,
        !  we won't bother trying to detect multiple instances of the same ECP.
        !
        call convert_one_atom(gam%atoms(ia),ecp%ecps(ie),eps_low,eps_high,eps_grid,do_report)
        call ecp_sanity_check(ecp%ecps(ie))
        ecp%ecps(ie)%projectors%rotmat    = rotmat
        ecp%ecps(ie)%projectors%rotmat_rk = real(rotmat,kind=kind(ecp%ecps(ie)%projectors%rotmat_rk))
      end do scan_atoms_for_ecps
      !
      if (do_report) then
        close (unit_report)
        write (out,"(/'Level-shift ECP projectors dumped to file ',a/)") trim(report_file)
      end if
      !
      ecp%active = .true.
    end subroutine ecp_convert
    !
    subroutine ecp_destroy(ecp)
      type(ecp_molecule), intent(inout) :: ecp
      !
      integer(ik) :: iprj
      !
      ecp%active = .false.
      destroy_prjs: do iprj=1,ecp%necps
        call gamess_destroy(ecp%ecps(iprj)%projectors)
        deallocate (ecp%ecps(iprj)%vshift,ecp%ecps(iprj)%rv)
        if (allocated(ecp%ecps(iprj)%grid)) deallocate(ecp%ecps(iprj)%grid)
      end do destroy_prjs
      deallocate (ecp%ecps)
    end subroutine ecp_destroy
    !
    subroutine ecp_evaluate_matrix_elements(gam,ecp,hecp)
      type(gam_structure), intent(inout) :: gam       ! Structure descriptor; for use with both the field-free and static-field calculations
      type(ecp_molecule), target         :: ecp       ! ECP descriptor, previosly prepared by ecp_convert
      real(rk), intent(out)              :: hecp(:,:) ! ECP contribution to the 1-electron Hamiltonian
      !
      integer(ik)                  :: iecp      ! Current ECP index
      integer(ik)                  :: nprj      ! Projector count in a ECP
      integer(ik)                  :: prj_nbas  ! Number of basis functions in a projector
      integer(ik)                  :: alloc
      integer(ik)                  :: nao
      type(ecp_atom), pointer      :: aecp      ! Current ECP
      type(gam_structure), pointer :: prj       ! Projector orbitals for the current ECP
      integer(ik)                  :: iprj, ia, ib
      !
      real(rk), allocatable :: ecp_prj_ao(:,:), ecp_prj(:,:)
      !
      if (ecp%necps<=0) return
      !
      call TimerStart('ECP Matrix elements')
      nao  = gam%nbasis
      hecp = 0
      apply_ecps: do iecp=1,ecp%necps
        aecp     => ecp%ecps(iecp)
        nprj     =  size(aecp%vshift)
        prj      => aecp%projectors
        prj_nbas =  prj%nbasis
        allocate (ecp_prj_ao(nao,prj_nbas),ecp_prj(nao,nprj),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i0,' allocating buffers for ECP projectors')") alloc
          stop 'static_tunnel%evaluate_ecps - memory allocation failed'
        end if
        !
        !  Evaluate projectors
        !
        call gamess_1e_integrals('AO OVERLAP',ecp_prj_ao,bra=gam,ket=prj)
        ecp_prj = mt_matmul(ecp_prj_ao,prj%vectors(:,1:nprj))
        !
        !  Accumulate the ECP contribution; this can be done more efficiently
        !  using DSYR2K, but we won't bother - this routine is not on a critical path (?)
        !
        apply_projectors: do iprj=1,nprj
          aos_right: do ib=1,nao
            aos_left: do ia=1,nao
              hecp(ia,ib) = hecp(ia,ib) + aecp%vshift(iprj)*ecp_prj(ia,iprj)*ecp_prj(ib,iprj)
            end do aos_left
          end do aos_right
        end do apply_projectors
        !
        deallocate (ecp_prj_ao,ecp_prj)
      end do apply_ecps
      if (verbose>=3) then
        write (out,"(/t5,'ECP INTEGRALS'/)")
        call gamess_print_1e_integrals(hecp,bra=gam,ket=gam)
      end if
      call TimerStop('ECP Matrix elements')
    end subroutine ecp_evaluate_matrix_elements
    !
    !  Internal routines
    !
    !  Make sure ECP projector orbitals we constructed are orthonormal
    !
    subroutine ecp_sanity_check(ec)
      type(ecp_atom), intent(in) :: ec  ! Single-atom ECP 
      !
      integer(ik)           :: nb          ! Number of basis functions
      integer(ik)           :: nv          ! Number of orbitals
      real(rk), allocatable :: ov_ao(:,:)  ! AO overlap integrals
      real(rk), allocatable :: ov_mo(:,:)  ! MO overlap integrals
      integer(ik)           :: ia, ib
      integer(ik)           :: alloc
      real(rk)              :: maxerr, err
      real(rk)              :: critical
      !
      nb = ec%projectors%nbasis
      nv = ec%projectors%nvectors
      allocate (ov_mo(nv,nv),ov_ao(nb,nb),stat=alloc)
      if (alloc/=0) then
        stop 'ecp_gamess%ecp_sanity_check - can''t allocate overlap matrix'
      end if
      !
      call gamess_1e_integrals('AO OVERLAP',ov_ao,ec%projectors,ec%projectors)
      ov_mo = matmul(transpose(ec%projectors%vectors),matmul(ov_ao,ec%projectors%vectors))
      !
      maxerr = 0.0_rk
      unity_ket: do ib=1,nv
        err = abs(ov_mo(ib,ib)-1.0_rk)
        maxerr = max(err,maxerr)
        unity_bra: do ia=1,nv
          if (ib==ia) cycle unity_bra
          err = abs(ov_mo(ia,ib))
          maxerr = max(err,maxerr)
        end do unity_bra
      end do unity_ket
      critical = max(spacing(1e5_rk),10._rk*zero_frac)
      if (verbose>=1 .or. maxerr>=critical) then
        write (out,"('ECP projector matrix for atom ',a,' deviates from unity by ',g12.5)") &
               trim(ec%projectors%atoms(1)%name), maxerr
      end if
      if (verbose>=2 .or. maxerr>=critical) then
        call print_matrix(ov_mo,14,'g14.7')
      end if
      deallocate (ov_ao,ov_mo)
      if (maxerr>=critical) stop 'ecp_gamess%ecp_sanity_check - bad ECP conversion'
    end subroutine ecp_sanity_check
    !
    !  Since the decontracted basis is very close to being linearly dependent,
    !  conversion of the ECP is extremely numerically delicate. In order to
    !  minimize chances for symmetry-breaking numerical errors to accumulate,
    !  we'll work with the spectral representation of all matrices.
    !
    subroutine convert_one_atom(ga,ec,eps_min,eps_max,eps_grid,do_report)
      type(gam_atom), intent(in)  :: ga        ! An ECP-carrying atom
      type(ecp_atom), intent(out) :: ec        ! Corresponding ECP in the shift-operator form
      real(rk), intent(in)        :: eps_min   ! Small (absolute) cut-off
      real(rk), intent(in)        :: eps_max   ! Large (positive) cut-off
      real(rk), intent(in)        :: eps_grid  ! Characteristic grid spacing where the projector will
                                               ! be used.
      logical, intent(in)         :: do_report ! Report projectors to unit_report
      !
      integer(ik)           :: nbas
      integer(ik)           :: alloc
      integer(ik)           :: ip, jp, np
      character(len=80)     :: comment
      real(rk), allocatable :: s (:,:)   ! Overlap matrix / Eigenvectors of same
      real(rk), allocatable :: es(:)     ! Eigenvalues of the overlap matrix or s^{-1/2}
      real(rk), allocatable :: u (:,:)   ! Matrix elements of the ECP / Eigenvectors of same
      real(rk), allocatable :: eu(:)     ! Eigenvalues of the ECP matrix
      real(rk), allocatable :: t (:,:)   ! Projector matrix / Eigenvectors of same
      real(rk), allocatable :: et(:)     ! Eigenvalues of the projector matrix
      real(rk), allocatable :: rv(:)     ! Expectation values of R for the projectors
      real(rk)              :: eps_s     ! Smallest eigenvalue of S to be included in inversion
      !
      if (verbose>=0) then
        write (out,"('Converting ECP on atom ',a,' at ',3(1x,f12.5),' Angstrom')") trim(ga%name), ga%xyz
      end if
      !
      !  Each projector is a one-atom "molecule"; do the basic initialization first
      !
      call initialize_blank_projector(ga,ec%projectors)
      !
      !  Construction of the shift operators requires a basis for
      !  resolution of the identity. We get that by decontracting 
      !  the orbital basis - at the very least, this should be enough
      !  for working with the occupied orbitals in that basis.
      !
      call decontract_basis(ec%projectors%atoms(1),ec%projectors%nbasis,eps_grid)
      call label_basis_set(ec%projectors)
      if (verbose>=2) then
        write (out,"(/'Decontracted basis set is:'/)")
        call report_gamess_data(out,ec%projectors,'Decontracted basis')
      end if
      !
      !  We require the inverse of the overlap matrix and ECP matrix elements
      !  to construct the RI expression
      !
      nbas = ec%projectors%nbasis
      allocate (s(nbas,nbas),es(nbas),u(nbas,nbas),eu(nbas),t(nbas,nbas),et(nbas),rv(nbas),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating temporary matrices of dimension ',i6)") alloc, nbas
        stop 'ecp_gamess%convert_one_atom - allocate failure (1)'
      end if
      !
      !  Overlap matrix. Once we have the spectral form, the inverse is trivial
      !
      call gamess_1e_integrals('AO OVERLAP',s,bra=ec%projectors,ket=ec%projectors)
      call symmetrize_matrix('RI basis overlap',s)
      call clean_matrix(s)
      if (verbose>=3) then
        write (out,"(/t10,'ECP RI basis overlap'/)") 
        call gamess_print_1e_integrals(s,bra=ec%projectors,ket=ec%projectors)
      end if
      ! call lapack_syev(s,es)
      call block_syev(s,es)
      if (verbose>=1) then
        write (out,"('Smallest and largest eigenvalues of the RI overlap: ',2g16.8)") minval(es), maxval(es)
        if (verbose>=2) then
          write (out,"(/'Eigenvalues of the RI overlap are:')")
          write (out,"(2x,10(1x,g12.5))") es
          write (out,"()")
        end if
      end if
      !
      !  Matrix elements of the ECP in the uncontracted basis
      !
      call one_centre_ecp_matrix(u,ec%projectors)
      call symmetrize_matrix('ECP matrix elements',u)
      call clean_matrix(u)
      if (verbose>=3) then
        write (out,"(/t10,'ECP matrix elements'/)") 
        call gamess_print_1e_integrals(u,bra=ec%projectors,ket=ec%projectors)
      end if
      ! call lapack_syev(u,eu)
      call block_syev(u,eu)
      if (verbose>=1) then
        write (out,"('Smallest and largest eigenvalues of the ECP matrix: ',2g16.8)") minval(eu), maxval(eu)
        if (verbose>=2) then
          write (out,"(/'Eigenvalues of the ECP matrix are:')")
          write (out,"(2x,10(1x,g12.5))") eu
          write (out,"()")
        end if
      end if
      !
      !  Lowdin-orthogonalized ECP matrix. We'd like to calculate S^{-1/2} U S^{-1/2}. 
      !  Using the spectral representation it becomes: 
      !     Cs Es^{-1/2} Cs^T Cu Eu Cu^T Cs Es^{-1/2} Cs^T, where
      !  Cs and Es are eigenvectors and eigenvalues of the overlap matrix, while
      !  Cu and Eu are eigenvectors and eigenvalues of the ECP matrix.
      !
      eps_s = 1e4_rk*sqrt(spacing(max(maxval(eu)-minval(eu),maxval(es))))
      if (verbose>=1) then
        write (out,"(1x,'S(RI) components smaller than ',g12.5,' will be eliminated')") eps_s
        write (out,"(1x,i4,' linearly-dependent components of S(RI) (out of ',i6,') are excluded')") &
               count(es<eps_s), nbas
      end if
      !
      !  Invert the overlap matrix
      !
      where (es>=eps_s)
        es = 1.0_rk/sqrt(es)
      elsewhere
        es = 0.0_rk
      end where
      !
      !  Construct the overall projector, taking care to preserve symmetries
      !
      t = matmul(matmul(transpose(u),s) * spread(es,dim=1,ncopies=nbas),transpose(s))
      !
      !  After this point, we no longer need U, but we'll still need S^{-1/2}
      !
      u = t
      call ut_ev_u(t,eu,u)
      call symmetrize_matrix('Lowdin-orthogonalized ECP',t)
      call clean_matrix(t)
      if (verbose>=3) then
        write (out,"(/t10,'Lowdin-orthogonalized ECP'/)") 
        call gamess_print_1e_integrals(t,bra=ec%projectors,ket=ec%projectors)
      end if
      ! call lapack_syev(t,et)
      call block_syev(t,et)
      if (verbose>=1) then
        write (out,"(/'Eigenvalues of the Lowdin-orthogonalized ECP projector matrix are:')")
        write (out,"(2x,10(1x,g12.5))") et
        write (out,"()")
      end if
      !
      !  Explicitly construct S^{-1/2} for eigenvector transformation
      !
      u = transpose(s)
      call ut_ev_u(s,es,u)
      call symmetrize_matrix('S^{-1/2}',s)
      call clean_matrix(s)
      if (verbose>=3) then
        write (out,"(/t10,'S^{-1/2}'/)") 
        call gamess_print_1e_integrals(s,bra=ec%projectors,ket=ec%projectors)
      end if
      !
      !  Transform the eigenvectors to the AO basis, avoiding temporaries
      !
      u = matmul(s,t) ; t = u
      !
      !  After this point, we do not need S^{-1/2}
      !  Construct matrix elements of <R> operator for the projector eigenvectors
      !
      call gamess_1e_integrals('AO 3C R',s,bra=ec%projectors,ket=ec%projectors, &
                               op_xyz=real(ec%projectors%atoms(1)%xyz,kind=kind(s))/abohr)
      call symmetrize_matrix('<R>',s)
      if (verbose>=3) then
        write (out,"(/t10,'<R>'/)") 
        call gamess_print_1e_integrals(s,bra=ec%projectors,ket=ec%projectors)
      end if
      u = matmul(s,t)
      convert_rv: do ip=1,nbas
        rv(ip) = dot_product(t(:,ip),u(:,ip))
      end do convert_rv
      if (verbose>=1) then
        write (out,"(/'<R> expectation values for the projectors are:')")
        write (out,"(2x,10(1x,g12.5))") rv
        write (out,"()")
      end if
      !
      !  Choose projectors to store
      !
      np = count((abs(et)>=eps_min).and.(et<=eps_max).and.(rv>=eps_grid))
      if (verbose>=0) then
        write (out,"('ECP eigenvalue cut-offs: low = ',g13.6,' high = ',g13.6)") eps_min, eps_max
        write (out,"('ECP small distance cut-off: ',g13.6)") eps_grid
        write (out,"('Keeping ',i5,' ECP level-shift operators out of ',i5)") np, nbas
      end if
      ec%projectors%nvectors = np
      allocate (ec%vshift(np),ec%rv(np),ec%projectors%vectors(nbas,np),stat=alloc)
      if (alloc/=0) then
        write (out,"('ecp_gamess%convert_one_atom: Error ',i8,' allocating space for ',i5,' projectors')") alloc, np
        stop 'ecp_gamess%convert_one_atom - allocate failure (2)'
      end if
      !
      ip = 0
      fill_projectors: do jp=1,nbas
        if (abs(et(jp))<eps_min .or. et(jp)>eps_max .or. rv(jp)<eps_grid) cycle fill_projectors
        ip = ip + 1
        ec%vshift(ip) = et(jp)
        ec%rv    (ip) = rv(jp)
        ec%projectors%vectors(:,ip) = t(:,jp)
      end do fill_projectors
      if (ip/=np) stop 'ecp_gamess%convert_one_atom - count error'
      call clean_matrix(ec%projectors%vectors)
      !
      call guess_cutoff_distance(ga,ec%rmax)
      if (verbose>=0) then
        write (out,"('ECP large distance cut-off is ',f14.7)") ec%rmax
      end if
      if (verbose>=0) then
        write (out,"(/t10,'ECP projection operators:'/)")
        write (out,"(3x,a5,2x,a14,2x,a14)") ' prj ', ' Shift, H ', '  <R>, Bohr  ', &
                                            '-----', '----------', '-------------'
        print_projectors: do ip=1,np
          write (out,"(3x,i5,2x,g14.7,2x,g14.7)") ip, ec%vshift(ip), ec%rv(ip)
        end do print_projectors
        write (out,"()")
      end if
      !
      deallocate (s,es,u,eu,t,et,rv,stat=alloc)
      !
      if (do_report) then
        write (comment,"('Level-shift representation of ECP ',a,' on ',a)") trim(ga%ecp_name), trim(ga%name)
        call report_gamess_data(unit_report,ec%projectors,trim(comment),'$VSHIFT',ec%vshift)
      end if
      !
    end subroutine convert_one_atom
    !
    subroutine initialize_blank_projector(ga,pr)
      type(gam_atom), intent(in)         :: ga ! An ECP-carrying atom
      type(gam_structure), intent(inout) :: pr ! The projector to initialize
      !
      integer(ik) :: alloc
      !
      pr%natoms   = 1
      pr%nbasis   = 0
      pr%nvectors = 0
      allocate (pr%atoms(1),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating atom table for a projector')") alloc
        stop 'ecp_gamess%initialize_blank_projector - allocate failed'
      end if
      pr%atoms = ga
      call MathSetUnitMatrix(pr%rotmat)
      pr%rotmat_rk = real(pr%rotmat,kind=kind(pr%rotmat_rk))
    end subroutine initialize_blank_projector
    !
    subroutine decontract_basis(pa,nbas,eps_grid)
      type(gam_atom), intent(inout) :: pa       ! Atom to be decontracted
      integer(ik), intent(out)      :: nbas     ! Total number of basis functions
      real(rk), intent(in)          :: eps_grid ! Grid spacing; do not generate functions which have
                                                ! <r> smaller than eps_grid, and thus are not representable.
      !
      integer(ik) :: c_ns                   ! Copy of the contracted nshell
      integer(ik) :: c_l(gam_max_shells)    ! Copy of the contracted sh_l
      integer(ik) :: c_p(gam_max_shells+1)  ! Copy of the contracted sh_p
      real(xrk)   :: c_z(gam_max_primitive) ! Copy of the contracted p_zet
      real(xrk)   :: c_c(gam_max_primitive) ! Copy of the contracted p_c
      integer(ik) :: is_c, is_u             ! Current contracted/uncontracted shell index
      integer(ik) :: ip_c, ip_u             ! Current primitive
!     real(rk)    :: rv                     ! Expectation of r for this shell
      !
      c_ns = pa%nshell
      c_l  = pa%sh_l
      c_p  = pa%sh_p
      c_z  = pa%p_zet
      c_c  = pa%p_c
      !
      pa%nshell = 0
      nbas      = 0
      contracted_shells: do is_c=1,c_ns
        contracted_primitives: do ip_c=c_p(is_c),c_p(is_c+1)-1
          ! write (out,*) 'cbf = ', is_c, ' l= ', c_l(is_c), ' z= ', c_z(ip_c), ' c= ', c_c(ip_c)
          !
          !  Let's see whether this primitive is already in the list
          !
          uncontracted_shells: do is_u=1,pa%nshell
            if (pa%sh_l(is_u)/=c_l(is_c)) cycle uncontracted_shells
            uncontracted_primitives: do ip_u=pa%sh_p(is_u),pa%sh_p(is_u+1)-1
              if (abs(pa%p_zet(ip_u)-c_z(ip_c))<=10*spacing(max(pa%p_zet(ip_u),c_z(ip_c)))) cycle contracted_primitives
            end do uncontracted_primitives
          end do uncontracted_shells
          ! !
          ! !  Primitive is new; check whether it is representable on the grid
          ! !
          ! rv = primitive_r_expectation(c_l(is_c),c_z(ip_c))
          ! if (rv<=eps_grid) cycle contracted_primitives
          !
          !  Pruning primitives according to their representability actually causes the results
          !  to detiriorate on coarse grids; do not do this. Pruning of non-representable projectors
          !  after the diagonalization should be sufficient.
          !
          !  Primitive is representable; add it. Since the new basis is uncontracted, this
          !  is a really easy job.
          !
          pa%nshell = pa%nshell + 1
          if (pa%nshell>gam_max_shells) then
            write (out,"('Generated too many shells decontracting atom ',a,'. The limit is ',i8)") trim(pa%name), gam_max_shells
            stop 'ecp_gamess%decontract_basis - too many shells'
          end if
          pa%sh_l (pa%nshell)   = c_l(is_c)
          pa%sh_p (pa%nshell)   = pa%nshell
          pa%sh_p (pa%nshell+1) = pa%nshell + 1
          pa%p_zet(pa%nshell)   = c_z(ip_c)
          pa%p_c  (pa%nshell)   = primitive_norm(c_l(is_c),c_z(ip_c))
          nbas                  = nbas + gam_orbcnt(c_l(is_c))
          ! write (out,"(' + ubf = ',i5,' l= ',i3,' z= ',g18.10,' c= ',g18.10)") &
          !        pa%nshell, pa%sh_l(pa%nshell), pa%p_zet(pa%nshell), pa%p_c(pa%nshell)
          !
          !  Keep reduced-precision copy, too
          !
          pa%p_zet_rk(pa%nshell) = real(pa%p_zet(pa%nshell),kind=kind(pa%p_zet_rk))
          pa%p_c_rk  (pa%nshell) = real(pa%p_c  (pa%nshell),kind=kind(pa%p_c_rk))
        end do contracted_primitives
      end do contracted_shells
      !
      if (verbose>0) then
        write (out,"('Decontraction of atom ',a,' gave ',i6,' shells and ',i8,' basis functions ')") &
               trim(pa%name), pa%nshell, nbas
      end if
    end subroutine decontract_basis
    !
    !  Although it would be perfectly possible to do these intergals analytically,
    !  one would have to be careful in handling the cartesian product functions
    !  in the basis set together with the pure spherical functions in the projectors.
    !  It is easier and hopefully less error-prone to evaluate the angular part 
    !  numerically instead. The radial part has a simple, closed-form analytical
    !  expression, so no worries there.
    !
    !  For the time being, we do not care over numerical efficiency overmuch - this
    !  is a one-shot operation.
    !
    subroutine one_centre_ecp_matrix(u,gam)
      real(rk), intent(out)           :: u(:,:)  ! ECP matrix elements
      type(gam_structure), intent(in) :: gam     ! One-atom "molecule"
      !
      integer(ik)              :: nshell            ! Number of shells/primitives in basis
      integer(ik)              :: ecp_nterms        ! Number of terms in the ECP
      integer(ik)              :: max_lbas          ! Maximum angular momentum in the basis
      integer(ik)              :: max_lecp          ! Maximum angular momentum in the ECP
      integer(ik)              :: max_necp          ! Maximum power of R in the ECP
      integer(ik)              :: nang              ! Number of angular terms in the basis
      integer(ik)              :: s_mu, s_nu, s_ecp ! Shell indices for the orbitals and ECP
      integer(ik)              :: l_mu, l_nu, n_ecp ! Powers or r appearing in the radial integral
      real(rk)                 :: z_mu, z_nu, d_ecp ! Exponents in the radial integral
      real(rk)                 :: c_mu, c_nu, c_ecp ! Weights in the radial integral
      integer(ik)              :: l_ecp             ! Angular momentum of the projector
      integer(ik)              :: p_mu, p_nu        ! Current position of the orbital in the u() buffer
      integer(ik)              :: alloc
      complex(rk), allocatable :: mu_ml(:,:,:)      ! Overlaps of angular parts of the basis functions
                                                    ! with spherical harmonics
      real(rk), allocatable    :: mu_nu(:,:)        ! Overlaps of angular parts of the basis functions
      real(rk)                 :: vrad              ! Radial integral
      !
      nshell     = gam%atoms(1)%nshell
      ecp_nterms = gam%atoms(1)%ecp_nterms
      max_lbas   = maxval(gam%atoms(1)%sh_l(1:nshell))
      max_lecp   = maxval(gam%atoms(1)%ecp_l(1:ecp_nterms))
      max_necp   = maxval(gam%atoms(1)%ecp_n(1:ecp_nterms))
      !
      !  Begin by precomputing the angular factors
      !
      nang = ang_loc(max_lbas+1)
      if (verbose>=3) then
        write (out,"('Angular buffer: nang = ',i3,' max lecp = ',i4)") nang, max_lecp
      end if
      allocate (mu_ml(0:nang-1,-max_lecp:max_lecp,0:max_lecp),mu_nu(0:nang-1,0:nang-1),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating angular factors buffer. nang = ',i5,' max_lecp = ',i5)") &
               alloc, nang, max_lecp
        stop 'ecp_gamess%one_centre_ecp_matrix - allocating angular factors'
      end if
      call fill_angular_factors(max_lbas,max_lecp,nang,mu_ml,mu_nu)
      !
      u    = 0
      p_nu = 1
      basis_shell_nu: do s_nu=1,nshell
        !
        !  Since we know that the basis is uncontracted, shell == primitive.
        !
        l_nu =      gam%atoms(1)%sh_l (s_nu)
        z_nu = real(gam%atoms(1)%p_zet(s_nu),kind=kind(vrad))
        c_nu = real(gam%atoms(1)%p_c  (s_nu),kind=kind(vrad))
        p_mu = 1
        basis_shell_mu: do s_mu=1,nshell
          l_mu =      gam%atoms(1)%sh_l (s_mu)
          z_mu = real(gam%atoms(1)%p_zet(s_mu),kind=kind(vrad))
          c_mu = real(gam%atoms(1)%p_c  (s_mu),kind=kind(vrad))
          ecp_term: do s_ecp=1,ecp_nterms
            n_ecp =      gam%atoms(1)%ecp_n(s_ecp)
            d_ecp = real(gam%atoms(1)%ecp_d(s_ecp),kind=kind(vrad))
            c_ecp = real(gam%atoms(1)%ecp_c(s_ecp),kind=kind(vrad))
            l_ecp =      gam%atoms(1)%ecp_l(s_ecp)
            ! In the basis functions, the polynomial terms are r**l_nu and r**l_mu
            ! In the ECP term, it is r**(n_ecp-2). 
            ! Finally, we have 4*pi*r**2 from the Jacobian.
            vrad = radial_integral(l_nu+l_mu+n_ecp,z_nu+z_mu+d_ecp)
            ! write (out,"(' l_mu = ',i3,' l_nu = ',i3,' n_ecp = ',i3,' exp = ',g14.7,' vrad = ',g16.8)") &
            !        l_nu, l_mu, n_ecp, z_nu+z_mu+d_ecp, vrad
            vrad = 4*real(pi_xrk,kind=kind(vrad))*c_mu*c_nu*c_ecp*vrad
            if (verbose>=4) then
              write (out,"(' s_mu = ',i4,' s_nu = ',i4,' s_ecp = ',i4,' l_ecp = ',i4,' rad = ',g16.8)") &
                     s_mu, s_nu, s_ecp, l_ecp, vrad
            end if
            !
            if (l_ecp< 0) call stuff_universal_integral
            if (l_ecp>=0) call stuff_projected_integral
          end do ecp_term
          p_mu = p_mu + gam_orbcnt(l_mu)
        end do basis_shell_mu
        if (p_mu/=size(u,dim=1)+1) stop 'ecp_gamess%one_centre_ecp_matrix - count error (mu)'
        p_nu = p_nu + gam_orbcnt(l_nu)
      end do basis_shell_nu
      if (p_nu/=size(u,dim=2)+1) stop 'ecp_gamess%one_centre_ecp_matrix - count error (nu)'
      !
      deallocate (mu_ml,mu_nu)
      !
      !  An internal subroutine would not be my first choice - but since I require
      !  access to pretty much all the local state, the alternative is too ugly
      !
      contains
      !
      !  Multiplicative contribution - weight by the angular overlap on the basis
      !
      subroutine stuff_universal_integral
        integer(ik) :: m_mu, m_nu 
        integer(ik) :: a_mu, a_nu
        real(rk)    :: tmp
        !
        stuff_nu: do m_nu=0,gam_orbcnt(l_nu)-1
          stuff_mu: do m_mu=0,gam_orbcnt(l_mu)-1
            a_mu = ang_loc(l_mu)+m_mu
            a_nu = ang_loc(l_nu)+m_nu
            tmp = vrad * mu_nu(a_mu,a_nu)
            u(p_mu+m_mu,p_nu+m_nu) = u(p_mu+m_mu,p_nu+m_nu) + tmp
            if (verbose>=4) then
              write (out,"('   element ',i6,' ',i6,' gets ',g16.8,' (LOC)')") p_mu+m_mu, p_nu+m_nu, + tmp
            end if
          end do stuff_mu
        end do stuff_nu
      end subroutine stuff_universal_integral
      !
      !  Non-local integral - weight by angular projectors.
      !  We'll do a dirty trick here: although the angular projectors are complex,
      !  we know that the final result must be purely real - so we'll accumulate
      !  just the real contribution.
      !
      subroutine stuff_projected_integral
        integer(ik) :: m_ecp
        integer(ik) :: m_mu, m_nu
        integer(ik) :: a_mu, a_nu
        complex(rk) :: tmp
        !
        stuff_nu: do m_nu=0,gam_orbcnt(l_nu)-1
          stuff_mu: do m_mu=0,gam_orbcnt(l_mu)-1
            a_mu = ang_loc(l_mu)+m_mu
            a_nu = ang_loc(l_nu)+m_nu
            tmp  = 0
            stuff_ecp_m: do m_ecp=-l_ecp,l_ecp
              ! In mu_ml, spherical harmonics are on the right-hand-side of the braket.
              tmp = tmp + mu_ml(a_mu,m_ecp,l_ecp) * conjg(mu_ml(a_nu,m_ecp,l_ecp))
            end do stuff_ecp_m
            ! Sanity check - just in case we got projectors wrong
            if (abs(aimag(tmp))>10._rk*spacing(1._rk+abs(real(tmp,kind=rk)))) then
              write (out,"('Result must be real, but isn''t: ',2g23.16)") tmp
              stop 'ecp_gamess%one_centre_ecp_matrix - spurious imaginary part'
            end if
            !
            !  An extra factor of 4*pi is because there are two angular integrals here,
            !  and we've only accounted for one so far
            !
            tmp = tmp * fourpi * vrad
            u(p_mu+m_mu,p_nu+m_nu) = u(p_mu+m_mu,p_nu+m_nu) + real(tmp,kind=rk)
            if (verbose>=4) then
              write (out,"('   element ',i6,' ',i6,' gets ',g16.8,' (NL)')") p_mu+m_mu, p_nu+m_nu, real(tmp,kind=rk)
            end if
          end do stuff_mu
        end do stuff_nu
      end subroutine stuff_projected_integral
    end subroutine one_centre_ecp_matrix
    !
    !  The somewhat awkward calling sequence is because the assumed-share calling convention 
    !  always resets the lower bound to 1. Bloody annoying.
    !
    subroutine fill_angular_factors(max_lbas,max_lecp,nang,mu_ml,mu_nu)
      integer(ik), intent(in)  :: max_lbas ! Maximum angular momentum in the basis
      integer(ik), intent(in)  :: max_lecp ! Maximum angular momentum in the ECP
      integer(ik), intent(in)  :: nang     ! Number of the angular factors in the basis; same as ang_loc(max_lbas+1)
      complex(rk), intent(out) :: mu_ml(0:nang-1,-max_lecp:max_lecp,0:max_lecp)
      real(rk), intent(out)    :: mu_nu(0:nang-1,0:nang-1)
      !
      real(rk), pointer        :: lg(:,:)     ! Lebedev grid of the chosen order
      integer(ik)              :: nleb        ! Number of grid points in the chosen Lebedev grid
      real(rk), allocatable    :: xyz(:,:,:)  ! Powers of distance at Lebedev points
      real(rk), allocatable    :: ang(:,:)    ! Basis angular factors at Lebedev points
      complex(rk), allocatable :: harm(:)     ! Values of the current harmonic at Lebedev points
      integer(ik)              :: alloc   
      integer(ik)              :: lecp, mecp  ! Angular momenta
      integer(ik)              :: ipt         ! Current Lebedev point
      integer(ik)              :: mu, nu      ! Current angular factor
      real(rk)                 :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
      !
      !  Prepare angular factors of the right kind
      !
      ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      !
      !  We'll just use 17-th order Lebedev grid. This should be more than enough to
      !  cover all angular momenta we'll ever see.
      ! 
      if (max_lbas+max(max_lbas,max_lecp)>17) stop 'ecp_gamess%fill_angular_factors - Blown LG17!'
      lg   => lebedev_gr17
      nleb = size(lg,dim=2)
      !
      allocate (xyz(nleb,3,0:5),ang(nleb,0:34),harm(nleb),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i6,' allocating angular intermediates. nleb = ',i7)") alloc, nleb
        stop 'ecp_gamess%fill_angular_factors - allocation error'
      end if
      !
      !  Calculate all possible angular factors for the basis functions at Lebedev points
      !  The code below is lifted directly from evaluate_basis_functions() in import_gamess.f90
      !  We rely on Lebedev grid points being on the surface of a unit sphere
      !
      xyz(:,:,0) = 1._ark
      xyz(:,:,1) = transpose(lg(1:3,:))
      xyz(:,:,2) = xyz(:,:,1)**2
      xyz(:,:,3) = xyz(:,:,1)*xyz(:,:,2)
      xyz(:,:,4) = xyz(:,:,2)**2
      xyz(:,:,5) = xyz(:,:,2)*xyz(:,:,3)
      ang = spread(ang_c_kind,1,nleb)*xyz(:,1,ang_nx)*xyz(:,2,ang_ny)*xyz(:,3,ang_nz)
      !
      !  Calculate spherical harmonics up to order max_lecp at all grid points.
      !  For the time being, we'll be calling FLharmonics; if this becomes too 
      !  expensive, a significant speedup can be obtained by calling legendrePnm_table
      !  in math.f90 directly - which gives all harmonics at once.
      !
      harmonics_l: do lecp=0,max_lecp
        harmonics_m: do mecp=-lecp,lecp
          call FLharmonicsSetParameters(lecp,mecp,(/0._rk,0._rk,0._rk/),0._rk,2._rk)
          harmonics_pt: do ipt=1,nleb
            harm(ipt) = FLharmonics(lg(1:3,ipt),(0._rk,0._rk))
          end do harmonics_pt
          !
          !  Ready to do integrals of this harmonic
          !
          basis_mu: do mu=0,nang-1
            mu_ml(mu,mecp,lecp) = sum(ang(:,mu)*harm*lg(4,:))
          end do basis_mu
        end do harmonics_m
      end do harmonics_l
      !
      !  Calculate overlaps of angular factors
      !
      overlap_nu: do nu=0,nang-1
        overlap_mu: do mu=0,nu
          mu_nu(mu,nu) = sum(ang(:,mu)*ang(:,nu)*lg(4,:))
          mu_nu(nu,mu) = mu_nu(mu,nu)  ! nop if mu==nu, so no worries
        end do overlap_mu
      end do overlap_nu
      !
      deallocate (xyz,ang,harm)
      !
      !  A bit of debugging output. Note that Lebedev grids are normalized to give 1 if the the
      !  integrand is unity; therefore, spherical harmonics' normalization will appear to be off
      !  by a factor (4*pi)**(-0.5)
      !
      if (verbose>=3) then
        write (out,"(/'Non-zero angular integrals: <bas_ang|bas_ang>')")
        write (out,"(1x,a3,1x,a4,2x,a3,1x,a4,2x,a12)") ' I ', ' SI ', ' J ', ' SJ ', ' <I|J> '
        print_mu_1: do mu=0,nang-1
          print_nu: do nu=0,nang-1
            if (abs(mu_nu(mu,nu))<spacing(10._rk)) cycle print_nu
            write (out,"(1x,i3,1x,a4,2x,i3,1x,a4,2x,f18.14)") mu, ang_name(mu), nu, ang_name(nu), mu_nu(mu,nu)
          end do print_nu
        end do print_mu_1
        write (out,"(/'Non-zero angular integrals: <bas_ang|ylm>')")
        write (out,"(1x,a3,1x,a4,2x,a3,1x,a3,3x,2a18)") ' I ', ' SI ', ' L ', ' M ', ' Re(<I|YLM>) ', ' Im(<I|YLM>) '
        print_mu_2: do mu=0,nang-1
          print_l: do lecp=0,max_lecp
            print_m: do mecp=-lecp,lecp
              if (abs(mu_ml(mu,mecp,lecp))<spacing(10._rk)) cycle print_m
              write (out,"(1x,i3,1x,a4,2x,i3,1x,i3,3x,2f18.14)") mu, ang_name(mu), lecp, mecp, mu_ml(mu,mecp,lecp)
            end do print_m
          end do print_l
        end do print_mu_2
        write (out,"()")
      end if
    end subroutine fill_angular_factors
    !
    !           inf
    !  Evaluate  |  r**n exp(-a*r**2) dr using Gauss-Hermite or Gauss-Laguerre quadrature
    !            0
    !  The answer is simply Gamma(0.5*(n+1))/(2*a**((n+1)/2)).
    !  Fortran 2008 also has gamma() intrisic, but who has a working F2008 compiler?
    !
    function radial_integral(n,a) result(v)
      integer(ik), intent(in) :: n       ! Power of r in the integral
      real(rk), intent(in)    :: a       ! Exponent
      real(rk)                :: v       ! The integral
      !
      real(rk) :: pow
      !
      pow = 0.5_rk * (n+1)
      v = real(MathSimpleGamma(cmplx(pow,0._rk,kind=rk))/(2._rk*(a**pow)),kind=rk)
   end function radial_integral
   !
   subroutine symmetrize_matrix(name,t)
     character(len=*), intent(in) :: name
     real(rk), intent(inout)      :: t(:,:)
     !
     integer(ik) :: ir, ic
     integer(ik) :: ir_max, ic_max
     real(rk)    :: ave, dif, max_dif
     !
     if (size(t,dim=1)/=size(t,dim=2)) stop 'ecp_gamess%symmetrize_matrix - can''t symmetrize rectangular matrix'
     max_dif = 0._rk
     ir_max  = 0
     ic_max  = 0
     scan_columns: do ic=2,size(t,dim=2)
       scan_rows: do ir=1,ic-1
         ave = 0.5_rk * (t(ir,ic) + t(ic,ir))
         dif =          (t(ir,ic) - t(ic,ir))
         t(ir,ic) = ave
         t(ic,ir) = ave
         if (abs(dif)>max_dif) then
           max_dif  = abs(dif)
           ir_max   = ir
           ic_max   = ic
         end if
       end do scan_rows
     end do scan_columns
     if (verbose>=1 .or. (max_dif>spacing(100._rk))) then
       write (out,"('Maximum deviation from symmetry in ',a,' = ',g12.5,' found at indices ',2i5)") &
              trim(name), max_dif, ir_max, ic_max
       if (max_dif>1e-3) stop 'ecp_gamess%symmetrize_matrix - error too large'
     end if
   end subroutine symmetrize_matrix
   !
   !  Inverse of diagonalization: calculate matrix product X = U^T e U, where e is diagonal
   !  We make no assumptions on U - in particular, it does not have to be unitary. If U is
   !  a matrix of left eigenvectors (transpose of the usual right eigenvectors), this operation
   !  will undo diagonalization.
   !
   subroutine ut_ev_u(x,eu,u)
     real(rk), intent(out) :: x (:,:) ! Reconstructed matrix
     real(rk), intent(in)  :: eu(:)   ! Diagonal matrix
     real(rk), intent(in)  :: u (:,:) ! Square matrix
     !
     integer(ik) :: i, j
     !
     !  This is not the most efficient way of calculating the result, but
     !  it does guarantee that the symmetry is maintained.
     !
     index_j: do j=1,size(x,dim=2)
       index_i: do i=1,j
         x(i,j) = sum(u(:,i)*u(:,j)*eu)
         x(j,i) = x(i,j)   ! a nop if i==j, so no worries
       end do index_i
     end do index_j
   end subroutine ut_ev_u
   !
   subroutine clean_matrix_real(mat)
     real(rk), intent(inout) :: mat(:,:) ! Matrix to cleanup
     !
     real(rk) :: thr
     !
     thr = zero_frac * maxval(abs(mat))
     !
     where (abs(mat)<=thr)
       mat = 0._rk
     end where
   end subroutine clean_matrix_real
   !
   subroutine clean_matrix_complex(mat)
     complex(rk), intent(inout) :: mat(:,:) ! Matrix to cleanup
     !
     real(rk) :: thr
     !
     thr = zero_frac * maxval(abs(mat))
     !
     where (abs(real(mat,kind=rk))<=thr)
       mat = cmplx(0._rk,aimag(mat),kind=rk)
     end where
     where (abs(aimag(mat))<=thr)
       mat = cmplx(real(mat,kind=rk),0._rk,kind=rk)
     end where
   end subroutine clean_matrix_complex
   !
   subroutine guess_cutoff_distance(ga,rmax)
     type(gam_atom), intent(in) :: ga    ! Atom with an ECP
     real(rk), intent(out)      :: rmax  ! Distance at which ECP can be neglected
     !
     integer(ik) :: ie
     real(rk)    :: rmax_term
     !
     rmax = -1._rk
     scan_ecp_terms: do ie=1,ga%ecp_nterms
       rmax_term = drop_dead_point(ga%ecp_n(ie),real(ga%ecp_c(ie),kind=kind(rmax)),real(ga%ecp_d(ie),kind=kind(rmax)))
       rmax      = max(rmax,rmax_term)
     end do scan_ecp_terms
     !
     contains
     function drop_dead_point(n,c,d) result(r)
       integer(ik), intent(in) :: n ! Power of r (including the Jacobian)
       real(rk), intent(in)    :: c ! Linear coefficient
       real(rk), intent(in)    :: d ! Exponential term
       real(rk)                :: r
       !
       real(rk) :: ac, rmax, vmax, r_inner, r_outer, v
       !
       ac = abs(c)
       !
       !  Find location of the maximum first; we'll monotonically decrease 
       !  past that point.
       !
       rmax = sqrt(0.5_rk*n/d)
       vmax = rmax**n * exp(-d*rmax**2)
       if (ac*vmax<=ecp_vcut) then
         !
         !  This term never rises above the threshold; ignore it
         !
         r = 0._rk
         return
       end if
       !
       !  Locate the outer point which satisfies the cut-off
       !
       r_inner = rmax
       scan_outwards: do 
         r_outer = r_inner * 2._rk
         v       = ac * r_outer**n * exp(-d*r_outer**2)
         if (v<=ecp_vcut) exit scan_outwards
         r_inner = r_outer
       end do scan_outwards
       !
       !  Bisect until we have a reasonably tight bound
       !
       bisect_interval: do while(r_outer-r_inner>100._rk*spacing(r_outer+r_inner))
         r = 0.5_rk * (r_inner + r_outer)
         v = ac * r**n * exp(-d*r**2)
         if (v<=ecp_vcut) then
           r_outer = r
         else
           r_inner = r
         end if
       end do bisect_interval
       r = r_outer
     end function drop_dead_point 
 
   end subroutine guess_cutoff_distance
  end module ecp_convert_gamess
