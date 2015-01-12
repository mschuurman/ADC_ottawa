!
!  A very naive implementation of complex spin-orbit CI using FON-SCF
!  references. This code is derived from (the currently stagnated)
!  DFT-CI project with Jochen Autschbach.
!
!  This implementation assumes that:
!
!  1. The list of determinants can be held in memory
!  2. Two-electron integrals can be held in memory
!  3. Sparse Hamiltonian matrix can be held in memory
!
!  Furthermore, we assume that the MOs are come from variationally-optimized
!  single-determinantal wavefunction (possibly with fractionally-occupied MOs),
!  so that the all one-electron matrix elements are subsumed by the orbital
!  eigenvalues. The requirement could be easity removed if necessary.
!
module ci_tools
  use accuracy
  use timer
  use integrals_mo2e
  use block_diag
  use sort_tools
  use arpack
  implicit none
  private
  public ci_data
  public ci_initialize_determinants
  public ci_direct_diagonalization  ! This is mostly a debugging routine!
  public ci_build_hamiltonian
  public ci_apply_hamiltonian
  public ci_diagonalization_arpack
  public ci_destroy
  !
  type ci_formula                    
    integer(hik)                     :: bra                  ! Determinant on the left
    integer(hik)                     :: ket                  ! Determinant on the right
    integer(sik)                     :: ijkl(4)              ! Orbital indices; j=0, k=0, l=1 means 1-electron formula
    integer(hik)                     :: ind_h                ! Index of the corresponding hij entry in ci%hamiltonian
    real(rk)                         :: wgt                  ! Coefficient 
  end type ci_formula   
  !
  type ci_helement   
    integer(hik)                     :: bra = 0              ! Determinant on the left
    integer(hik)                     :: ket = 0              ! Determinant on the right
    complex(rk)                      :: hbk = 0              ! Matrix element
  end type ci_helement   
  !                                  
  type ci_data
    integer(ik)                      :: norbs                ! Number of active orbitals in current CI
    integer(ik)                      :: nelectrons           ! Number of electrons in current CI
    character(len=20)                :: mode                 ! Either 'sz ormas' or 'ormas'
    real(rk)                         :: fix_sz               ! Non-negative value restricts CI to this Sz.
    real(rk), allocatable            :: orbital_sz(:)        ! Per-orbital Sz expectation: <L|Sz|R>
    real(rk), allocatable            :: occ_ref(:)           ! Occupation numbers for the reference determinant
    integer(ik)                      :: ormas_npart          ! Number of ORMAS partitions
    integer(ik), allocatable         :: ormas_orbitals(:)    ! Starting orbital index of each orbital in an ORMAS partitions
                                                             ! The first entry must be 1
                                                             ! The last entry [ormas_orbitals(ormas_npart)] marks the orbital
                                                             ! one past the last active partition
    integer(ik), allocatable         :: ormas_mine(:)        ! Minimum number of electrons allowed in a partition
    integer(ik), allocatable         :: ormas_maxe(:)        ! Maximum number ... etc
    integer(hik)                     :: ndets                ! Number of determinants in current CI
    integer(sik), allocatable        :: dets(:,:)            ! List of determinants. The first index (1:nelectrons)
                                                             ! lists the occupied spin-orbitals in this determinant.
                                                             ! The second index goes from 1 to ndets
    integer(hik)                     :: nform                ! Number of Hamiltonian formulae
    integer(hik)                     :: nham                 ! Number of (potentially) non-zero matrix elements in the Hamiltonian
    type(ci_formula), allocatable    :: formulae(:)          ! Hamiltonian matrix element formulae
                                                             ! Note that all formula entries for the same determinant pair is
                                                             ! guaranteed to be grouped together; this is important for hamiltonian
                                                             ! construction
    type(ci_helement), allocatable   :: hamiltonian(:)       ! CI hamiltonian, in general sparse-matrix format
  end type ci_data
  !
  !  Fixed parameters
  !
  integer(ik), parameter   :: verbose       = 1
  !
  contains
  !
  !  Externally-visible entry points
  !
  subroutine ci_initialize_determinants(ci,mode,norbs,nelectrons,occ_ref, &
                                        fix_sz,orbital_sz, &
                                        ormas_orbitals,ormas_mine,ormas_maxe)
    type(ci_data), intent(inout)      :: ci                ! CI data descriptor
    character(len=*), intent(in)      :: mode              ! Currently this must be 'ormas' or 'sz ormas'
                                                           ! If mode=='sz *', fix_sz and orbital_sz must be supplied
                                                           ! If mode=='*ormas', ormas_orbitals, ormas_mine, and ormas_maxe must be
                                                           ! given
    integer(ik), intent(in)           :: norbs             ! Number of orbitals in current CI
    integer(ik), intent(in)           :: nelectrons        ! Number of electrons in current CI
    real(rk), intent(in)              :: occ_ref(:)        ! Orbital occupation numbers for the reference determinant
    real(rk), intent(in), optional    :: fix_sz            ! Using mode='sz *' will restrict CI to the specified Sz value
                                                           ! orbital_sz(:) must be supplied as well.
    real(rk), intent(in), optional    :: orbital_sz(:)     ! Per-orbital Sz value
    integer(ik), intent(in), optional :: ormas_orbitals(:) ! Starting orbital index of each orbital in an ORMAS partitions
                                                           ! The last entry marks the orbital one past the last active partition
    integer(ik), intent(in), optional :: ormas_mine(:)     ! Minimum number of electrons allowed in a partition
    integer(ik), intent(in), optional :: ormas_maxe(:)     ! Maximum number ... etc
    !
    integer(ik) :: alloc
    !
    call TimerStart('CI Determinant list')
    ci%mode       = mode
    ci%norbs      = norbs
    ci%nelectrons = nelectrons
    if (size(occ_ref)/=ci%norbs) then
      stop 'ci_tools%ci_initialize_determinants - size(occ_ref)/=norbs'
    end if
    allocate (ci%occ_ref(norbs),stat=alloc)
    if (alloc/=0) stop 'ci_tools%ci_initialize_determinants - allocation failure (AZ)'
    ci%occ_ref = occ_ref
    !
    if (ci%mode(1:3)=='sz ') then
      if (.not.present(fix_sz) .or. .not.present(orbital_sz)) then
        stop 'ci_tools%ci_initialize_determinants - fix_sz and/or orbital_sz arguments are missing'
      end if
      if (size(orbital_sz)/=ci%norbs) then
        stop 'ci_tools%ci_initialize_determinants - size(orbitals_sz)/=norbs'
      end if
      ci%fix_sz = fix_sz
      allocate (ci%orbital_sz(norbs),stat=alloc)
      if (alloc/=0) stop 'ci_tools%ci_initialize_determinants - allocation failure (A)'
      ci%orbital_sz = orbital_sz
    end if
    !
    select case (ci%mode)
      case default
        write (out,"('CI Active space choice ',a,' is not recognized. Abort.')") trim(ci%mode)
        stop 'ci_tools%ci_initialize_determinants - unknown active space'
      case ('ormas','sz ormas')
        if (.not.present(ormas_orbitals) .or. .not.present(ormas_mine) .or. .not.present(ormas_maxe)) then
          stop 'ci_tools%ci_initialize_determinants - some of ormas_* arguments are missing'
        end if
        ci%ormas_npart = size(ormas_maxe)
        if (size(ormas_orbitals)/=ci%ormas_npart+1 .or. size(ormas_mine)/=ci%ormas_npart) then
          stop 'ci_tools%ci_initialize_determinants - ORMAS array sizes are inconsistent'
        end if
        allocate (ci%ormas_orbitals(ci%ormas_npart+1),ci%ormas_mine(ci%ormas_npart),ci%ormas_maxe(ci%ormas_npart),stat=alloc)
        if (alloc/=0) stop 'ci_tools%ci_initialize_determinants - allocation failure (B)'
        ci%ormas_orbitals = ormas_orbitals
        ci%ormas_mine     = ormas_mine
        ci%ormas_maxe     = ormas_maxe
        !
        call count_determinants_ormas(ci,record=.false.)
        write (out,"('Found ',i0,' determinants')") ci%ndets
        if (ci%ndets==0) stop 'ci_tools%ci_initialize_determinants - can''t do CI with no determinants'
        !
        allocate (ci%dets(ci%nelectrons,ci%ndets),stat=alloc)
        if (alloc/=0) stop 'ci_tools%ci_initialize_determinants - allocation failure (C)'
        call count_determinants_ormas(ci,record=.true.)
    end select
    !
    !  Print determinant list
    !
    if (verbose>=2) then
      write (out,"(/'Determinants are:'/)")
      call print_determinants(ci)
      write (out,"()")
    end if
    call TimerStop('CI Determinant list')
    call TimerStart('CI Formula tape')
    !
    !  Prepare formula tapes for evaluating the Hamiltonian
    !
    call count_formulae(ci,record=.false.)
    write (out,"('Number of (potentially) non-zero CI Hamiltonian matrix elements: ',i0)") ci%nham
    write (out,"('                              Number of entries in formula tape: ',i0)") ci%nform
    write (out,"('Estimated memory requirements for the formula tape: ',f0.3,' Gbytes')") &
           1e-9_rk*((2*hik_bytes+4*sik_bytes+rk_bytes+0._rk)*ci%nform)
    call flush(out)
    allocate (ci%formulae(ci%nform),stat=alloc)
    if (alloc/=0) stop 'ci_tools%ci_initialize_determinants - allocation failure (D)'
    call count_formulae(ci,record=.true.)
    !
    call TimerStop('CI Formula tape')
  end subroutine ci_initialize_determinants
  !
  subroutine ci_destroy(ci)
    type(ci_data), intent(inout) :: ci ! CI data descriptor
    !
    call TimerStart('CI Destroy')
    if (allocated(ci%orbital_sz)    ) deallocate (ci%orbital_sz)
    if (allocated(ci%occ_ref)       ) deallocate (ci%occ_ref)
    if (allocated(ci%ormas_orbitals)) deallocate (ci%ormas_orbitals)
    if (allocated(ci%ormas_mine)    ) deallocate (ci%ormas_mine)
    if (allocated(ci%ormas_maxe)    ) deallocate (ci%ormas_maxe)
    if (allocated(ci%dets)          ) deallocate (ci%dets)
    if (allocated(ci%formulae)      ) deallocate (ci%formulae)
    if (allocated(ci%hamiltonian)   ) deallocate (ci%hamiltonian)
    call TimerStop('CI Destroy')
  end subroutine ci_destroy
  !
  !  Debugging aid: direct diagonalization of the CI Hamiltonian
  !
  subroutine ci_direct_diagonalization(ci,eref,eps,moint2e)
    type(ci_data), intent(inout)       :: ci      ! CI data descriptor
    complex(rk), intent(in)            :: eref    ! Energy of the reference state
    complex(rk), intent(in)            :: eps(:)  ! Orbital eigenvalues of the active orbitals
    type(moint2e_cache), intent(inout) :: moint2e ! 2e integrals over the active orbitals
    !
    complex(rk), allocatable :: hmat(:,:)   ! Full Hamiltonian matrix
    complex(rk), allocatable :: evec(:,:,:) ! Left and right eigenvectors
    complex(rk), allocatable :: eval(:)     ! Eigenvalues
    integer(ik), allocatable :: order(:)    ! Ordering table for the eigenvalues
    integer(ik)              :: alloc
    integer(hik)             :: iev, idet
    integer(ik)              :: iorb
    real(rk)                 :: maxc
    !
    call TimerStart('CI Diagonalization')
    allocate (hmat(ci%ndets,ci%ndets),evec(ci%ndets,ci%ndets,2),eval(ci%ndets),order(ci%ndets),stat=alloc)
    if (alloc/=0) then
      stop 'ci_tools%ci_direct_diagonalization - allocation failed'
    end if
    ! call stupid_build_ci_hamiltonian(hmat,ci,eps,moint2e)
    call ci_build_hamiltonian(ci,eps,moint2e)
    call ci_expand_hamiltonian(ci,hmat)
    !
    evec(:,:,1) = hmat
    call block_geev(evec,eval)
    call order_keys(real(eval,kind=rk),order)
    eval = eval(order)
    evec = evec(:,order,:)
    !
    write (out,"(/t5,'CI eigenvalues')")
    print_eigenvalues: do iev=1,ci%ndets
      write (out,"('CI state ',i5,' Energy = ',g20.13,1x,g20.13)") iev, eref + eval(iev)
      if (verbose>=1) then
        maxc = maxval(abs(evec(:,iev,:)))
        print_dets: do idet=1,ci%ndets
          if (maxval(abs(evec(idet,iev,:)))<0.05_rk*maxc) cycle print_dets
          print_orbitals: do iorb=1,ci%nelectrons
            write (out,"(1x,i3)",advance='no') ci%dets(iorb,idet)
          end do print_orbitals
          write (out,"(' | ',i8,1x,2(1x,f12.8,1x,f12.8,3x))") idet, evec(idet,iev,:)
        end do print_dets
      end if
    end do print_eigenvalues
    !
    deallocate (hmat,evec,eval,order)
    call TimerStop('CI Diagonalization')
    call TimerReport
    call flush(out)
  end subroutine ci_direct_diagonalization
  !
  subroutine ci_diagonalization_arpack(ci,eref,guess_e,guess_c,eps_root,max_iter,n_root,success,final_e,final_c)
    type(ci_data), intent(in)    :: ci          ! CI data descriptor
    complex(rk), intent(in)      :: eref        ! Reference energy; all CI matrix elements are relative to this value
    complex(rk), intent(in)      :: guess_e     ! Guess for the energy
    complex(rk), intent(in)      :: guess_c(:)  ! Guess for the right eigenvector. Passing zero vector
                                                ! will start diagonalization from a random vector
    real(rk), intent(in)         :: eps_root    ! Desired convergence
    integer(ik), intent(in)      :: max_iter    ! Maximum allowed number of iterations
    integer(ik), intent(in)      :: n_root      ! Number of eigenvalues to compute
    logical, intent(out)         :: success     ! Diagonalization converged?
    complex(rk), intent(out)     :: final_e     ! Final energy
    complex(rk), intent(out)     :: final_c(:)  ! Final right eigenvector; chosen to be most similar in 
                                                ! Cartesian space to the guess
    !
    !  Parameters passed to ARPACK; use ARPACK type conventions instead of ours
    !
    integer                   :: ido, info, ldv, lworkl, n, ncv, nev, ierr, ldz
    real(kind=drk)            :: tol
    character(len=1)          :: bmat
    character(len=2)          :: which
    integer                   :: iparam(11), ipntr(14)
    complex(drk), allocatable :: resid(:), v(:,:), workd(:), workl(:), d(:), workev(:), z(:,:)
    complex(drk)              :: sigma
    real(drk), allocatable    :: rwork(:)
    logical                   :: rvec
    logical, allocatable      :: select(:)
    !
    integer(ik)              :: alloc, product_count, iev, good_iev, iorb
    logical                  :: random_guess
    real(rk)                 :: lg, lv, maxc
    complex(rk), allocatable :: lgv(:)
    character(len=1)         :: good_lab
    integer(hik)             :: idet
    !
    if (ci%ndets>huge(n)) stop 'ci_tools%ci_diagonalization_arpack - number of determinants too large for an integer'
    !
    ido       = 0                     ! Indicates first call to znaupd
    bmat      = 'I'                   ! Solving standard eigenvalue problem
    n         = int(ci%ndets)         ! Dimension of the eigenproblem
    which     = 'SM'                  ! Compute eigenvalues smallest in magnitude (ie closest to guess_e)
    nev       = n_root                ! How many eigenvalues to compute
    tol       = eps_root              ! Desired convergence
    ncv       = min(n,max(10,3*nev))  ! Max. number of expansion vectors
    ldv       = n                     ! Size of the expansion vector again
    ldz       = n                     
    iparam    = 0                     
    iparam(1) = 1                     ! 1 means use exact shifts
    iparam(3) = max_iter              ! Max number of Arnoldi update iterations allowed
    iparam(7) = 1                     ! Solve eigenproblem
    lworkl    = 3*ncv**2 + 5*ncv      ! Size of the work area
    info      = count(abs(guess_c)>1e-5_rk)
    random_guess = info==0
    !
    if (size(guess_c)/=n) stop 'ci_tools%ci_diagonalization_arpack - inconsistent eigenvector sizes'
    if (size(final_c)/=n) stop 'ci_tools%ci_diagonalization_arpack - inconsistent eigenvector sizes'
    !
    write (out,"()") 
    allocate (resid(n), v(ldv,ncv), z(ldz,nev+1), workd(3*n), workl(lworkl), &
              rwork(ncv), d(nev+1), select(ncv), workev(3*ncv), lgv(nev), stat=alloc)
    if (alloc/=0) stop 'ci_tools%ci_diagonalization_arpack - allocation failed'
    !
    product_count = 0
    resid = guess_c
    arpack_iterations: do
      call arpack_znaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)
      select case (ido)
        case default
          write (out,"('znaupd returned ido = ',i0,' which we can''t handle')") ido
          stop 'ci_tools%ci_diagonalization_arpack - unexpected ido'
        case (-1,1) 
          if (kind(1d0)/=rk) stop 'ci_tools%ci_diagonalization_arpack - assumption on real kinds failed'
          call ci_apply_hamiltonian(ci,'right',shift=-(guess_e-eref),vec=workd(ipntr(1):ipntr(1)+n-1), &
                                                             Hvec=workd(ipntr(2):ipntr(2)+n-1))
          product_count = product_count + 1
          if (verbose>=1) then
            if (mod(product_count,50)==0) then
              write (out,"('CI: evaluated ',i0,' H.vec products so far')") product_count
              call flush(out)
            end if
          end if
        case (99)
         exit arpack_iterations ! Converged
      end select
    end do arpack_iterations
    !
    if (info<0) then
      success = .false.
      write (out,"('znaupd returned error ',i0)") info
      stop 'ci_tools%ci_diagonalization_arpack - znaupd failed'
    else
      write (out,"('CI diagonalization converged for ',i0,' eigenstates')") iparam(5)
      write (out,"('Number of Hamiltonian-vector products: ',i0)") product_count
      write (out,"('Number of Arnoldi iterations: ',i0)") iparam(3)
      call flush(out)
      rvec = .true.
      call arpack_zneupd  (rvec, 'A', select, d, z, ldz, sigma,  &
           workev, bmat, n, which, nev, tol, resid, ncv,         &
           v, ldv, iparam, ipntr, workd, workl, lworkl,          &
           rwork, ierr)
      if (ierr/=0) then
        write (out,"('Error ',i0,' occured in arpack_zneupd')") ierr
        stop 'ci_tools%ci_diagonalization_arpack - znaupd failed'
      end if
      !
      !  We now have to decide which eigenpair is the correct solution.
      !  There are two options: If we had a guess for the starting vector
      !  (random_guess=.false.), we will use Cartesian similarity between
      !  the old and the new vector. If we started with a random vector,
      !  go for the smallest distance from the guess energy.
      !
      lgv = 0
      if (.not.random_guess) then
        lg = sqrt(abs(dot_product(guess_c,guess_c)))
        calculate_distance: do iev=1,iparam(5)
          lv = sqrt(abs(dot_product(z(:,iev),z(:,iev))))
          lgv(iev) = dot_product(z(:,iev),guess_c)/(lg*lv)
        end do calculate_distance
        good_iev = maxloc(abs(lgv(:iparam(5))),dim=1)
      else
        good_iev = minloc(abs(d(:iparam(5))),dim=1)
      end if
      !
      print_eigenvalues: do iev=1,iparam(5)
        good_lab = ' '
        if (iev==good_iev) good_lab = '*'
        write (out,"('CI state ',i5,' Energy = ',g20.13,1x,g20.13,' <final|guess> =',2(1x,f14.9),1x,a)") &
               iev, eref + d(iev), lgv(iev), good_lab
        if (verbose>=1) then
          maxc = maxval(abs(z(:,iev)))
          print_dets: do idet=1,ci%ndets
            if (abs(z(idet,iev))<0.05_rk*maxc) cycle print_dets
            print_orbitals: do iorb=1,ci%nelectrons
              write (out,"(1x,i3)",advance='no') ci%dets(iorb,idet)
            end do print_orbitals
            write (out,"(' | ',i8,1x,2(1x,f12.8,1x,f12.8,3x))") idet, z(idet,iev)
          end do print_dets
        end if
      end do print_eigenvalues
      success = .true.
      final_e = eref + d(good_iev)
      final_c = z(:,good_iev)
    end if
    !
  end subroutine ci_diagonalization_arpack
  !
  !  Apply CI Hamiltonian at a given right vector
  !
  subroutine ci_build_hamiltonian(ci,eps,moint2e)
    type(ci_data), intent(inout)       :: ci      ! CI data descriptor
    complex(rk), intent(in)            :: eps(:)  ! Orbital eigenvalues of the active orbitals
    type(moint2e_cache), intent(inout) :: moint2e ! 2e integrals over the active orbitals
    !
    integer(ik) :: alloc
    integer(ik) :: lv
    !
    call TimerStart('CI Build Hamiltonian')
    if (allocated(ci%hamiltonian)) deallocate (ci%hamiltonian)
    allocate (ci%hamiltonian(ci%nham),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i0,' allocating ',i0,'-element sparse Hamiltonian')") alloc, ci%nham
      stop 'ci_tools%ci_build_hamiltonian - allocation failure'
    end if
    select case (moint2e%mode)
      case default
        write (out,"('Integral storage mode ',a,' is not supported in ci_build_hamiltonian')") trim(moint2e%mode)
        stop 'ci_tools%ci_build_hamiltonian - bad integral storage mode'
      case ('incore')
        call ci_hamiltonian_pass(ci,eps,moint2e,0)
      case ('disk')
        integral_blocks: do lv=1,moint2e%nmo(4)
          if (verbose>=1) then
            write (out,"('CI Hamiltonian: integral block ',i0,' of ',i0)") lv, moint2e%nmo(4)
          end if
          call fetch_moint2e(moint2e,lv)
          call ci_hamiltonian_pass(ci,eps,moint2e,lv)
        end do integral_blocks
    end select
    call TimerStop('CI Build Hamiltonian')
  end subroutine ci_build_hamiltonian
  !
  subroutine ci_apply_hamiltonian(ci,side,shift,vec,Hvec)
    type(ci_data), intent(in)          :: ci      ! CI data descriptor
    character(len=*), intent(in)       :: side    ! 'right' = Calculate H vec
                                                  ! 'left'  = Calculate H^T vec
    complex(rk), intent(in)            :: shift   ! Overall energy shift
    complex(rk), intent(in)            :: vec(:)  ! Amplitudes of the determinants, in the order of ci%dets
    complex(rk), intent(out)           :: Hvec(:) ! The product.
    !
    integer(hik) :: ind_h, bra, ket
    complex(rk)  :: term
    !
    call TimerStart('CI Apply Hamiltonian')
    !
    !$omp parallel do default(none) private(bra) shared(ci,Hvec)
    zero_result: do bra=1,ci%ndets
      Hvec(bra) = 0
    end do zero_result
    !$omp end parallel do
    select case (side)
      case default
        write (out,"('ci_apply_hamiltonian: side ',a,' is not recognized')") trim(side)
        stop 'ci_tools%ci_apply_hamiltonian - bad side'
      case ('right')
        !$omp parallel do default(none) shared(ci,vec,Hvec,shift) private(ind_h,bra,ket,term)
        apply_to_the_right: do ind_h=1,ci%nham
          bra = ci%hamiltonian(ind_h)%bra
          ket = ci%hamiltonian(ind_h)%ket
          if (bra==0 .or. ket==0) cycle apply_to_the_right ! It turned out to be zero after all
          term = (ci%hamiltonian(ind_h)%hbk+shift) * vec(ket)
          !$omp atomic
          Hvec(bra) = Hvec(bra) + term
        end do apply_to_the_right
        !$omp end parallel do
      case ('left')
        !$omp parallel do default(none) shared(ci,vec,Hvec,shift) private(ind_h,bra,ket,term)
        apply_to_the_left: do ind_h=1,ci%nham
          bra = ci%hamiltonian(ind_h)%bra
          ket = ci%hamiltonian(ind_h)%ket
          if (bra==0 .or. ket==0) cycle apply_to_the_left ! It turned out to be zero after all
          term = (ci%hamiltonian(ind_h)%hbk+shift) * vec(bra)
          !$omp atomic
          Hvec(ket) = Hvec(ket) + term
        end do apply_to_the_left
        !$omp end parallel do
    end select
    call TimerStop('CI Apply Hamiltonian')
  end subroutine ci_apply_hamiltonian
  !
  !  Internal subroutines
  !
  !  Enumerate all determinants within a given ORMAS space.
  !  We have two versions of this routine: one counts the determinants,
  !  the second stores the list of the determinants it found.
  !
  !  The algorithm is simple:
  !    a. In the outermost loop, we'll go through all distributions of
  !       electron occupation numbers within the ORMAS spaces consistent 
  !       with the ORMAS ranges and the total electron count.
  !    b. For fixed population counts, each ORMAS partition is simply
  !       a CAS, so we'll go through each one of those, too
  !    c. If necessary, we'll discard determinants with wrong Sz.
  !
  subroutine count_determinants_ormas(ci,record)
    type(ci_data), intent(inout) :: ci     ! CI data descriptor
    logical, intent(in)          :: record ! Set to .true. if determinants should be recorded, rather than just counted
    !
    integer(ik)               :: part_size(ci%ormas_npart) ! Number of orbitals in each partition
    integer(ik)               :: part_nume(ci%ormas_npart) ! Number of electrons in each partition
    integer(ik)               :: max_maxe                  ! Largest possible number of orbitals in any of the partitions
    integer(sik), allocatable :: part_dets(:,:)            ! Orbital lists within each subspace; note that spin-orbitals 
                                                           ! within each subspace are numbered starting from 1.
    integer(sik)              :: det(ci%nelectrons)        ! Overall determinant
    integer(hik)              :: idet                      ! Current determinant index
    integer(ik)               :: alloc
    integer(ik)               :: ip, ip2
    integer(ik)               :: iorb
    real(rk)                  :: sz
    logical                   :: good_det
    !
    max_maxe = maxval(ci%ormas_maxe)
    allocate (part_dets(max_maxe,ci%ormas_npart),stat=alloc)
    if (alloc/=0) stop 'ci_tools%count_determinants_ormas - allocation failed'
    !
    idet = 0
    call first_occupations_set(ci,part_size,part_nume)
    scan_occupation_sets: do 
      ! write (out,"('next occ: ',10i5)") part_nume
      if (sum(part_nume)==ci%nelectrons) then
        !
        !  This is a valid set of occupations, generate determinants
        !
        ! write (out,"('good occ: ',10i5)") part_nume
        initialize_partitions: do ip=1,ci%ormas_npart
          call first_cas_determinant(part_nume(ip),part_dets(:part_nume(ip),ip))
        end do initialize_partitions
        scan_partition_determinants: do
          !
          !  Construct the overall determinant from partial orbital lists in the subspaces
          !
          iorb = 1
          merge_partitions: do ip2=1,ci%ormas_npart
            det(iorb:iorb+part_nume(ip2)-1) = int(ci%ormas_orbitals(ip2),kind=sik) + part_dets(:part_nume(ip2),ip2) - 1_sik
            iorb = iorb + part_nume(ip2)
          end do merge_partitions
          if (iorb-1/=ci%nelectrons) stop 'ci_tools%count_determinants_ormas - electron count error'
          !
          good_det = .true.
          if (ci%mode(1:3)=='sz ') then
            sz = sum(ci%orbital_sz(det))
            good_det = (abs(sz-ci%fix_sz)<=1e-3_rk)
          end if
          !
          if (good_det) then
            ! write (out,"('det: ',10i5)") det
            idet = idet + 1
            if (record) then
              if (idet>ci%ndets) stop 'ci_tools%count_determinants_ormas - determinants count error'
              ci%dets(:,idet) = det
            end if
          end if
          !
          if (.not.next_partition_determinant(ci%ormas_npart,part_size,part_nume,part_dets)) exit scan_partition_determinants
        end do scan_partition_determinants
      end if
      if (.not.next_occupation_set(ci,part_nume)) exit scan_occupation_sets
    end do scan_occupation_sets
    deallocate (part_dets)
    !
    ! If we were counting determinants, remember the total
    !
    if (.not.record) ci%ndets = idet
  end subroutine count_determinants_ormas
  !
  !  We'll start from maximum possible occupation numbers in all partitions,
  !  then count down, decrementing highest-numbered partitions first.
  !
  subroutine first_occupations_set(ci,part_size,part_nume)
    type(ci_data), intent(in) :: ci 
    integer(ik), intent(out)  :: part_size(:)
    integer(ik), intent(out)  :: part_nume(:)
    !
    part_size = ci%ormas_orbitals(2:) - ci%ormas_orbitals(1:ci%ormas_npart)
    part_nume = ci%ormas_maxe
  end subroutine first_occupations_set
  !
  function next_occupation_set(ci,part_nume) result(success)
    type(ci_data), intent(in)  :: ci 
    integer(ik), intent(inout) :: part_nume(:)
    logical                    :: success
    !
    integer(ik) :: ipart, ip2
    !
    success = .false.
    advance_positions: do ipart=ci%ormas_npart,1,-1
      if (part_nume(ipart)<=ci%ormas_mine(ipart)) cycle advance_positions
      part_nume(ipart) = part_nume(ipart) - 1
      ! Once an occupation number has been decremented, we need to reset all
      ! higher-numbered partitions to the maximum value
      advance_rest: do ip2=ipart+1,ci%ormas_npart
        part_nume(ip2) = ci%ormas_maxe(ip2)
      end do advance_rest
      success = .true.
      exit advance_positions
    end do advance_positions
  end function next_occupation_set

  function next_partition_determinant(npart,part_size,part_nume,part_dets) result(success)
    integer(ik), intent(in)     :: npart          ! Number of partitions
    integer(ik), intent(in)     :: part_size(:)   ! Number of orbitals in each partition
    integer(ik), intent(in)     :: part_nume(:)   ! Number of electrons in each partition
    integer(sik), intent(inout) :: part_dets(:,:) ! Determinants in each partition
    logical                     :: success
    !
    integer(ik) :: ip, ip2
    !   
    success = .false.
    advance_partitions: do ip=1,npart
      if (.not.next_cas_determinant(part_size(ip),part_nume(ip),part_dets(:part_nume(ip),ip))) cycle advance_partitions
      ! Reset orbital counters for all preceeding partitions
      advance_rest: do ip2=1,ip-1
        call first_cas_determinant(part_nume(ip2),part_dets(:part_nume(ip2),ip2))
      end do advance_rest
      success = .true.
      exit advance_partitions
    end do advance_partitions
  end function next_partition_determinant

  subroutine first_cas_determinant(nelectrons,det)
    integer(ik), intent(in)   :: nelectrons
    integer(sik), intent(out) :: det(nelectrons)
    integer(sik)              :: ipos
    !
    fill_first_determinant: do ipos=1,int(nelectrons,kind=sik)
      det(ipos) = ipos
    end do fill_first_determinant
  end subroutine first_cas_determinant
  !
  !  Go to the next determinant within the complete active space
  !
  function next_cas_determinant(norbs,nelectrons,det) result(success)
    integer(ik), intent(in)     :: norbs
    integer(ik), intent(in)     :: nelectrons
    integer(sik), intent(inout) :: det(nelectrons) 
    logical                     :: success         ! Will be .false. if the input was 
                                                  ! the last determinant in sequence
    integer(sik) :: ipos, ip2
    !
    success = .false.
    advance_positions: do ipos=0,int(nelectrons-1,kind=sik)
      if (det(nelectrons-ipos)==norbs-ipos) cycle advance_positions
      det(nelectrons-ipos) = det(nelectrons-ipos) + 1_sik
      advance_rest: do ip2=int(nelectrons-ipos+1,kind=sik),int(nelectrons,kind=sik)
        det(ip2) = det(ip2-1) + 1_sik
      end do advance_rest
      success = .true.
      exit advance_positions
    end do advance_positions
  end function next_cas_determinant
  !
  !  Print determinants
  !
  subroutine print_determinants(ci)
    type(ci_data), intent(inout) :: ci ! CI data descriptor
    !
    integer(hik) :: idet
    integer(ik)  :: iorb
    !
    print_dets: do idet=1,ci%ndets
      print_orbitals: do iorb=1,ci%nelectrons
        write (out,"(1x,i3)",advance='no') ci%dets(iorb,idet)
      end do print_orbitals
      write (out,"(' | ',i8)") idet
    end do print_dets
  end subroutine print_determinants
  !
  subroutine count_formulae(ci,record)
    type(ci_data), intent(inout) :: ci      ! CI data descriptor
    logical, intent(in)          :: record  ! Record the coefficients; don't just count them
    ! 
    integer(hik) :: bra, ket
    integer(ik)  :: far
    !
    !  Reset counters
    !
    ci%nham  = 0
    ci%nform = 0
    !
    !  First, construct the diagonal of the Hamiltonian matrix.
    !
    h_diagonal: do bra=1,ci%ndets
      ci%nham = ci%nham + 1
      call field_free_h_diagonal(ci,record,bra,ci%dets(:,bra))
    end do h_diagonal
    !
    !  The off-diagonal terms. There will be a lot of zeros; it is therefore
    !  OK to check distance in parallel, then compute formulae inside the
    !  critical section.
    !
    !$omp parallel do default(none) private(bra,ket,far) shared(ci,record)
    h_off_diagonal_ket: do ket=1,ci%ndets
      h_off_diagonal_bra: do bra=1,ci%ndets
        if (bra == ket) cycle h_off_diagonal_bra
        !
        ! far = determinant_distance(ci%dets(:,bra),ci%dets(:,ket))
        far = determinant_distance_short(ci%dets(:,bra),ci%dets(:,ket))
        ! far = determinant_distance_long(ci%dets(:,bra),ci%dets(:,ket),ci%norbs)
        if (far<=0 .or. far>2) cycle h_off_diagonal_bra
        !$omp critical
        ci%nham = ci%nham + 1
        if (far==1) then
          call field_free_h_1far(ci,record,bra,ket,ci%dets(:,bra),ci%dets(:,ket))
        else
          call field_free_h_2far(ci,record,bra,ket,ci%dets(:,bra),ci%dets(:,ket))
        endif
        !$omp end critical
      end do h_off_diagonal_bra
    end do h_off_diagonal_ket
    !$omp end parallel do
  end subroutine count_formulae
  !
  subroutine add_form_1e(ci,record,bra,i,wgt)
    type(ci_data), intent(inout) :: ci      ! CI data descriptor
    logical, intent(in)          :: record  ! Record the coefficients; don't just count them
    integer(hik), intent(in)     :: bra     ! Determinant index
    integer(ik), intent(in)      :: i       ! Orbital index
    real(rk), intent(in)         :: wgt     ! Coefficient to record
    !
    if (abs(wgt)<=spacing(1e2_rk)) return ! Ignore zero entries
    ci%nform = ci%nform + 1
    if (.not.record) return
    !
    if (.not.allocated(ci%formulae)) stop 'ci_tools%add_form_1e - record==.true., but formulae is not allocated'
    if (ci%nform > size(ci%formulae)) stop 'ci_tools%add_form_1e - counting error'
    !
    ci%formulae(ci%nform)%bra   = bra
    ci%formulae(ci%nform)%ket   = bra
    ! Placing 1 in the L position makes it easier to do multipass hamiltonian construction
    ci%formulae(ci%nform)%ijkl  = int((/i,0_ik,0_ik,1_ik/),kind=sik)
    ci%formulae(ci%nform)%ind_h = ci%nham
    ci%formulae(ci%nform)%wgt   = wgt
  end subroutine add_form_1e
  !
  subroutine add_form_2e(ci,record,bra,ket,ijkl,wgt)
    type(ci_data), intent(inout) :: ci       ! CI data descriptor
    logical, intent(in)          :: record   ! Record the coefficients; don't just count them
    integer(hik), intent(in)     :: bra, ket ! Determinant indices
    integer(ik), intent(in)      :: ijkl(:)  ! Orbital indices
    real(rk), intent(in)         :: wgt      ! Coefficient to record
    !
    if (abs(wgt)<=spacing(1e2_rk)) return ! Ignore zero entries
    ci%nform = ci%nform + 1
    if (.not.record) return
    !
    if (.not.allocated(ci%formulae)) stop 'ci_tools%add_form_2e - record==.true., but formulae is not allocated'
    if (ci%nform > size(ci%formulae)) stop 'ci_tools%add_form_2e - counting error'
    !
    ci%formulae(ci%nform)%bra   = bra
    ci%formulae(ci%nform)%ket   = ket
    ci%formulae(ci%nform)%ind_h = ci%nham
    ci%formulae(ci%nform)%wgt   = wgt
    !
    !  Our 2e integrals are guaranteed to have (ij,kl) = (kl,ij) symmetry
    !  Since it is considerably more expensive to transform over the fourth index
    !  in our implementation, we'll try to choose the copy where the last index is
    !  as small as possible.
    !
    if (ijkl(4)<=ijkl(2)) then
      ci%formulae(ci%nform)%ijkl = int(ijkl,kind=sik)
    else
      ci%formulae(ci%nform)%ijkl = int(ijkl((/3,4,1,2/)),kind=sik)
    end if
  end subroutine add_form_2e
  !
  !  Calculate energy of a given determinant, relative to the reference energy.
  !
  subroutine field_free_h_diagonal(ci,record,bra,det)
    type(ci_data), intent(inout) :: ci      ! CI data descriptor
    logical, intent(in)          :: record  ! Record the coefficients; don't just count them
    integer(hik), intent(in)     :: bra     ! Determinant index
    integer(sik), intent(in)     :: det(:)  ! Orbitals occupied in this determinant
    !
    integer(ik) :: i, j           ! Orbital indices
    real(rk)    :: wgt         
    real(rk)    :: docc(ci%norbs) ! Orbital occupations in the current determinant
    !
    !  Fill orbital occupations for each active orbital - this makes loops below simpler
    !
    docc      = 0
    docc(det) = 1
    !
    !  ... one-particle part of the energy
    !
    hf_1e: do i=1,ci%norbs
      ! vo = vo + eps(i)*(docc(i)-nocc_ref(i))
      call add_form_1e(ci,record,bra,i,docc(i)-ci%occ_ref(i))
    end do hf_1e
    !
    !  Two-electron energy
    !
    hf_2e_j: do j=1,ci%norbs
      hf_2e_i: do i=1,ci%norbs
        if (i==j) cycle hf_2e_i
        wgt = 0.5_rk * (docc(i)*(docc(j)-ci%occ_ref(j)) - (docc(i)-ci%occ_ref(i))*ci%occ_ref(j))
        call add_form_2e(ci,record,bra,bra,(/i,i,j,j/), wgt)  ! Integrals are in charge-cloud notation
        call add_form_2e(ci,record,bra,bra,(/i,j,j,i/),-wgt)
      end do hf_2e_i
    end do hf_2e_j
  end subroutine field_free_h_diagonal
! !
! !  Off-diagonal elements of the field-free Hamiltonian. Again we could have 
! !  multiple versions, depending on the mode_hij1
! !
  subroutine field_free_h_1far(ci,record,p_bra,p_ket,bra,ket)
    type(ci_data), intent(inout) :: ci             ! CI data descriptor
    logical, intent(in)          :: record         ! Record the coefficients; don't just count them
    integer(hik), intent(in)     :: p_bra, p_ket   ! Determinant indices
    integer(sik), intent(in)     :: bra(:), ket(:) ! Orbitals occupied in the two determinants
    !
    integer(ik)             :: ub, uk         ! Orbital indices which are not the same in bra and ket
    integer(ik)             :: par            ! Parity factor: +1/-1
    integer(ik)             :: k              ! Orbital index
    real(rk)                :: docc(ci%norbs) ! Orbital occupations
    real(rk)                :: wgt
    !
    !  Figure out which orbitals are not the same, and get the parity factor
    !
    call parity_1far(bra,ket,ub,uk,par)
    !
    docc(:)   = 0._rk
    docc(bra) = 1._rk
    !
    hf_ref: do k=1,ci%norbs
      if (k==ub .or. k==uk) cycle hf_ref ! This term must be zero due to integrals symmetry
      wgt = par*(docc(k)-ci%occ_ref(k))
      call add_form_2e(ci,record,p_bra,p_ket,(/ub,uk,k,k /), wgt)
      call add_form_2e(ci,record,p_bra,p_ket,(/ub,k, k,uk/),-wgt)
    end do hf_ref
  end subroutine field_free_h_1far
  !
  subroutine field_free_h_2far(ci,record,p_bra,p_ket,bra,ket)
    type(ci_data), intent(inout) :: ci             ! CI data descriptor
    logical, intent(in)          :: record         ! Record the coefficients; don't just count them
    integer(hik), intent(in)     :: p_bra, p_ket   ! Determinant indices
    integer(sik), intent(in)     :: bra(:), ket(:) ! Orbitals occupied in the two determinants
    !
    integer(ik) :: ub1, uk1    ! Orbital indices which are not the same in bra and ket
    integer(ik) :: ub2, uk2    ! Orbital indices which are not the same in bra and ket
    integer(ik) :: par         ! Parity factor: +1/-1
    real(rk)    :: wgt
    !
    !  Figure out which orbitals are not the same, and get the parity factor
    !
    call parity_2far(bra,ket,ub1,ub2,uk1,uk2,par)
    !
    wgt = par
    call add_form_2e(ci,record,p_bra,p_ket,(/ub1,uk1,ub2,uk2/), wgt)
    call add_form_2e(ci,record,p_bra,p_ket,(/ub1,uk2,ub2,uk1/),-wgt)
  end subroutine field_free_h_2far
  !
  !  Useful routines for dealing with determinants: counting minimal
  !  sets of substitutions and parities.
  !
  function determinant_distance(bra,ket) result(lbk)
    integer(sik), intent(in) :: bra(:), ket(:)  ! Determinants, represented by sorted orbital lists
    integer(ik)              :: lbk             ! Number of substitutions separating the determinants
    !
    integer(ik) :: nelectrons
    integer(ik) :: ob, ok
    !
    if (size(bra)/=size(ket)) stop 'ci_tools%determinant_distance - bra and ket differ in size'
    nelectrons = size(bra)
    lbk = 0
    if (nelectrons>0) then ! It's unlikely we get called with no electrons, but it'a a valid input
      ob = 1 ; ok = 1
      scan_bra: do
        !
        !  Compare the current orbitals, and advance the pointers
        !
        if (bra(ob)==ket(ok)) then
          ob = ob + 1 ; ok = ok + 1
        else if (bra(ob)<ket(ok)) then
          ob = ob + 1 ; lbk = lbk + 1
        else !   bra(ob)>ket(ok)
          ok = ok + 1 ; lbk = lbk + 1
        end if
        !
        !  Exit tests. If any orbitals are left uncounted at the end, they are
        !  obviously different
        !
        if (ob>nelectrons .and. ok>nelectrons) exit scan_bra
        if (ob>nelectrons) then
          lbk = lbk + (nelectrons-ok+1)
          exit scan_bra
        end if
        if (ok>nelectrons) then
          lbk = lbk + (nelectrons-ob+1)
          exit scan_bra
        end if
      end do scan_bra
    end if
    !
    !  The distance we calculate here should always be an even number
    !
    if (mod(lbk,2)/=0) stop 'determinant_distance - impossible distance!'
    lbk = lbk / 2
    if (verbose>=3) then
      write (out,"('bra: ',20i3)") bra
      write (out,"('ket: ',20i3)") ket
      write (out,"('distance = ',i3/)") lbk
    end if
  end function determinant_distance
  !
  !  Same as determinant_distance, but stop counting once distance exceeds 2.
  !  For all integrals we care, we only need to distinguish cases of the distance
  !  being 0, 1, 2, or more than 2.
  !
  function determinant_distance_short(bra,ket) result(lbk)
    integer(sik), intent(in) :: bra(:), ket(:)  ! Determinants, represented by sorted orbital lists
    integer(ik)              :: lbk             ! Number of substitutions separating the determinants
    !
    integer(ik) :: nelectrons
    integer(ik) :: ob, ok
    !
    if (size(bra)/=size(ket)) stop 'ci_tools%determinant_distance_short - bra and ket differ in size'
    nelectrons = size(bra)
    lbk = 0
    if (nelectrons>0) then ! It's unlikely we get called with no electrons, but it'a a valid input
      ob = 1 ; ok = 1
      scan_bra: do
        !
        !  Compare the current orbitals, and advance the pointers
        !
        if (bra(ob)==ket(ok)) then
          ob = ob + 1 ; ok = ok + 1
        else if (bra(ob)<ket(ok)) then
          ob = ob + 1 ; lbk = lbk + 1
        else !   bra(ob)>ket(ok)
          ok = ok + 1 ; lbk = lbk + 1
        end if
        !
        !  Quick exit if too many differences already
        !
        if (lbk>=3*2) exit scan_bra
        !
        !  Exit tests. If any orbitals are left uncounted at the end, they are
        !  obviously different
        !
        if (ob>nelectrons .and. ok>nelectrons) exit scan_bra
        if (ob>nelectrons) then
          lbk = lbk + (nelectrons-ok+1)
          exit scan_bra
        end if
        if (ok>nelectrons) then
          lbk = lbk + (nelectrons-ob+1)
          exit scan_bra
        end if
      end do scan_bra
    end if
    !
    !  The distance we calculate here should always be an even number
    !
    if (mod(lbk,2)/=0) stop 'determinant_distance_short - impossible distance!'
    lbk = lbk / 2
  end function determinant_distance_short
  !
  function determinant_distance_long(bra,ket,maxorb) result(lbk)
    integer(sik), intent(in) :: bra(:), ket(:)  ! Determinants, represented by sorted orbital lists
    integer(sik), intent(in) :: maxorb          ! Maximum number of orbitals
    integer(ik)              :: lbk             ! Number of substitutions separating the determinants
    !
    logical :: orbmask(maxorb)
    !
    orbmask      = .false.
    orbmask(bra) = .true.
    orbmask(ket) = .not.orbmask(ket)
    lbk = count(orbmask)
    if (mod(lbk,2)/=0) stop 'determinant_distance_long - impossible distance!'
    lbk = lbk / 2
  end function determinant_distance_long
  !
  !  For two determinants, figure out which orbitals are not the same, 
  !  and get the parity factor coming from swapping the different orbitals
  !  between the maximally aligned and canonical orbital order.
  !
  subroutine parity_nfar(bra,ket,ub,uk,par)
    integer(sik), intent(in) :: bra(:), ket(:) ! Orbitals occupied in the two determinants
    integer(ik), intent(out) :: ub(:)          ! Unique orbitals in bra/ket. Must be long enough to list all
    integer(ik), intent(out) :: uk(:)          ! unique orbitals; there is no checking. The remainder of the
                                               ! array after the unique orbitals will be initialized to zero.
                                               ! ib/ik orbitals are always in pairs; For each pair, replacing
                                               ! ib(i) by ik(i) will yield the parity factor par
    integer(ik), intent(out) :: par            ! Parity factor: +1/-1
    !
    integer(ik) :: nelectrons ! Number of electrons
    integer(ik) :: ex         ! Number of excess unique orbitals in bra
    integer(ik) :: ob, ok     ! Current electron in bra or ket
    integer(ik) :: pb, pk     ! Unique orbital pointer in ub/uk
    !
    if (size(bra)/=size(ket)) stop 'ci_tools%parity_nfar - bra and ket differ in length'
    nelectrons = size(bra)
    ex = 0 ; par = 0 ;
    ub = 0 ; uk = 0 ;
    pb = 1 ; pk = 1 ;
    ob = 1 ; ok = 1 ;
    scan_bra: do
      if (bra(ob)==ket(ok)) then
        ob = ob + 1 ; ok = ok + 1 ;
        par = par + ex                 ! We need to exchange this orbital with abs(ex) others
      else if (bra(ob)<ket(ok)) then
        ub(pb) = bra(ob) ; pb = pb + 1 ! Add one more unique orbital in bra
        ob = ob + 1 ; 
        ex = ex + 1                    ! One more excess bra orbital to swap through
      else !   bra(ob)>ket(ok)
        uk(pk) = ket(ok) ; pk = pk + 1 ! Add one more unique orbital in ket
        ok = ok + 1 ; 
        ex = ex - 1                    ! One less excess bra orbital, or an excess ket
      end if
      if (ob<=nelectrons .and. ok<=nelectrons) cycle scan_bra 
      !
      !  Exit. If both strings run out simultaneously, there is nothing to do.
      !  Otherwise, the remainder of the longer string has unique orbitals
      !
      if (ob<=nelectrons) then ! All orbitals remaining in the bra are unique
        ub(pb:pb+nelectrons-ob) = bra(ob:nelectrons)
        pb = pb + 1 + nelectrons-ob
      end if
      if (ok<=nelectrons) then ! All orbitals remaining in the ket are unique
        uk(pk:pk+nelectrons-ok) = ket(ok:nelectrons)
        pk = pk + 1 + nelectrons-ok
      end if
      exit scan_bra
    end do scan_bra
    !
    !  Sanity check
    !
    if (pk/=pb) stop 'parity_nfar - logic error'
    !
    if (mod(par,2)/=0) then
      par = -1
    else
      par = 1
    end if
    !
    if (verbose>=3) then
      write (out,"('        bra: ',20i3)") bra
      write (out,"('        ket: ',20i3)") ket
      write (out,"(' unique bra: ',20i3)") ub
      write (out,"(' unique ket: ',20i3)") uk
      write (out,"('     parity: ',i3/)") par
    end if
  end subroutine parity_nfar
  !
  !  Wrapper around parity_nfar - determinants should be 1 orbital apart
  !
  subroutine parity_1far(bra,ket,ub,uk,par)
    integer(sik), intent(in) :: bra(:), ket(:) ! Orbitals occupied in the two determinants
    integer(ik), intent(out) :: ub, uk         ! Unique orbitals pair in bra/ket
    integer(ik), intent(out) :: par            ! Parity factor: +1/-1
    !
    integer(ik) :: uba(size(bra)+1), uka(size(ket)+1)
    !
    call parity_nfar(bra,ket,uba,uka,par)
    if (uba(1)==0 .or. uba(2)/=0 .or. uka(1)==0 .or. uka(2)/=0) stop 'parity_1far - determinants not 1 apart'
    ub = uba(1)
    uk = uka(1)
  end subroutine parity_1far
  !
  !  Wrapper around parity_nfar - determinants should be 2 orbital apart
  !
  subroutine parity_2far(bra,ket,ub1,ub2,uk1,uk2,par)
    integer(sik), intent(in) :: bra(:), ket(:) ! Orbitals occupied in the two determinants
    integer(ik), intent(out) :: ub1, ub2       ! Unique orbital pairs: ub1/uk1 and  ub2/uk2
    integer(ik), intent(out) :: uk1, uk2       ! 
    integer(ik), intent(out) :: par            ! Parity factor: +1/-1
    !
    integer(ik) :: uba(size(bra)+1), uka(size(ket)+1)
    !
    call parity_nfar(bra,ket,uba,uka,par)
    if (uba(2)==0 .or. uba(3)/=0 .or. uka(2)==0 .or. uka(3)/=0) stop 'parity_2far - determinants not 2 apart'
    ub1 = uba(1) ; ub2 = uba(2)
    uk1 = uka(1) ; uk2 = uka(2)
  end subroutine parity_2far
  !
  !  Do one pass over integrals for CI Hamiltonian construction
  !
  subroutine ci_hamiltonian_pass(ci,eps,moint2e,lv)
    type(ci_data), intent(inout)       :: ci        ! CI data descriptor
    complex(rk), intent(in)            :: eps(:)    ! Orbital eigenvalues of the active orbitals
    type(moint2e_cache), intent(inout) :: moint2e   ! 2e integrals over the active orbitals
    integer(ik)                        :: lv        ! L index to process; L=0 means process everything
    !
    integer(hik) :: it, bra, ket, ind_h, old_bra, old_ket
    integer(sik) :: i, j, k, l
    complex(rk)  :: hbk
    real(rk)     :: wgt
    !A
    if (.not.allocated(ci%hamiltonian)) stop 'ci_tools%build_ci_hamiltonian - hamiltonian is not allocated'
    evaluate_formulae: do it=1,ci%nform
      l     = ci%formulae(it)%ijkl(4)
      if (l/=lv .and. lv/=0) cycle evaluate_formulae
      !
      bra     = ci%formulae(it)%bra
      ket     = ci%formulae(it)%ket
      i       = ci%formulae(it)%ijkl(1)
      j       = ci%formulae(it)%ijkl(2)
      k       = ci%formulae(it)%ijkl(3)
      ind_h   = ci%formulae(it)%ind_h
      wgt     = ci%formulae(it)%wgt
      if (ind_h<=0 .or. ind_h>ci%nham) stop 'ci_tools%build_ci_hamiltonian - hamiltonian entry out of range'
      old_bra = ci%hamiltonian(ind_h)%bra
      old_ket = ci%hamiltonian(ind_h)%ket
      hbk     = ci%hamiltonian(ind_h)%hbk
      if (old_bra==0 .and. old_ket==0) then
        ci%hamiltonian(ind_h)%bra = bra
        ci%hamiltonian(ind_h)%ket = ket
      else if (old_bra/=bra .or. old_ket/=ket) then
        stop 'ci_tools%build_ci_hamiltonian - inconstent formulae encountered'
      end if
      if (j==0) then
        !
        ! Special case: 1e formula
        !
        if (ket/=bra .or. k/=0 .or. l/=1) stop 'ci_tools%build_ci_hamiltonian - bad 1e formula'
        hbk = hbk + wgt*eps(i)
      else
        if (lv==0) then
          hbk = hbk + wgt*moint2e%buffer_real(i,j,k,l)
        else
          hbk = hbk + wgt*moint2e%buffer_real(i,j,k,1)
        end if
      end if
      ci%hamiltonian(ind_h)%hbk = hbk
    end do evaluate_formulae
  end subroutine ci_hamiltonian_pass
  !
  subroutine ci_expand_hamiltonian(ci,hmat)
    type(ci_data), intent(in)          :: ci        ! CI data descriptor
    complex(rk), intent(out)           :: hmat(:,:) ! Full Hamiltonian matrix
    !
    integer(hik) :: ind_h, bra, ket
    !
    hmat = 0
    nonzero_elements: do ind_h=1,ci%nham
      bra = ci%hamiltonian(ind_h)%bra
      ket = ci%hamiltonian(ind_h)%ket
      if (bra==0 .or. ket==0) cycle nonzero_elements ! It turned out to be zero after all
      hmat(bra,ket) = ci%hamiltonian(ind_h)%hbk
    end do nonzero_elements
  end subroutine ci_expand_hamiltonian
  !
  !  Explicitly construct CI Hamiltonian using formula tape and integrals
  !  This is a debugging aid, not a production routine
  !
  subroutine stupid_build_ci_hamiltonian(hmat,ci,eps,moint2e)
    complex(rk), intent(out)           :: hmat(:,:) ! Full Hamiltonian matrix
    type(ci_data), intent(in)          :: ci        ! CI data descriptor
    complex(rk), intent(in)            :: eps(:)    ! Orbital eigenvalues of the active orbitals
    type(moint2e_cache), intent(inout) :: moint2e   ! 2e integrals over the active orbitals
    !
    integer(hik) :: it, bra, ket
    integer(sik) :: i, j, k, l
    real(rk)     :: wgt
    !
    if (moint2e%mode/='incore') stop 'ci_tools%build_ci_hamiltonian - only incore integrals are supported'
    !
    hmat = 0
    evaluate_formulae: do it=1,ci%nform
      bra = ci%formulae(it)%bra
      ket = ci%formulae(it)%ket
      i   = ci%formulae(it)%ijkl(1)
      j   = ci%formulae(it)%ijkl(2)
      k   = ci%formulae(it)%ijkl(3)
      l   = ci%formulae(it)%ijkl(4)
      wgt = ci%formulae(it)%wgt
      if (j==0) then
        !
        ! Special case: 1e formula
        !
        if (ket/=bra .or. k/=0 .or. l/=1) stop 'ci_tools%build_ci_hamiltonian - bad 1e formula'
        hmat(bra,bra) = hmat(bra,bra) + wgt*eps(i)
      else
        hmat(bra,ket) = hmat(bra,ket) + wgt*moint2e%buffer_real(i,j,k,l)
      end if
    end do evaluate_formulae
  end subroutine stupid_build_ci_hamiltonian
end module ci_tools
