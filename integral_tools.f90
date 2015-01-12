!
!  Useful routines for managing complete sets of molecular 2e integrals
!  Everything here is really a bunch of wrappers around gamess_2e_integrals
!
!  The standard call sequence is expected to be:
!
!    call prepare_2e ! At the very beginning of the calculation
!
!       iobatch=1 ... 
!         call fetch_ao2e_batch
!
!         [ parallel loop over integral batches]
!                             call process_block
!            if (swap_jikl()) call process_block
!            ...
!
!    call clear_2e   ! Once integrals are no longer needed
!
!  See fock_tools.f90 for a usage example
!
!  Note that all integrals produced bt this module are spin-less.
!
!  Integrals are in the charge-cloud order: the first two indices
!  correspond to the variable x1; the last two indices correspond
!  to the variable x2.
!
!  prepare_2e, fetch_ao2e_batch, and clear_2e are not thread-safe, but will
!  run in parallel (OpenMP) where appropriate
!  
!  swap_* can be called from parallel regions
!
!  swap_*_s will only swap the indices, but will not copy the integrals.
!
!  Multiple integral contexts can be active at the same time, provided
!  that the corresponding I/O units are different.
!
  module integral_tools
    use accuracy
    use timer
    use math
    use import_gamess
    use gamess_internal
    implicit none
    private
    public int2e_cache
    public prepare_2e, clear_2e, fetch_ao2e_batch
    public swap_jikl, swap_ijlk, swap_jilk, swap_lkij, swap_klij, swap_klji, swap_lkji
    public swap_jikl_s, swap_ijlk_s, swap_jilk_s, swap_lkij_s, swap_klij_s, swap_klji_s, swap_lkji_s
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter  :: verbose  = 1              ! Level of output
    !
    !  ==== User-adjustable parameters =====
    !
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    !  Variables related to 2-electron integral handling
    !
    type int2e_cache
      integer(ik)                    :: nbatch_ao        ! Number of single-index AO batches
      integer(ik)                    :: maxbatch_ao      ! Largest AO batch
      integer(ik), allocatable       :: indices_ao(:,:)  ! Index 1: first orbital in batch; # of orbitals in batch; atom batch belongs to
                                                         ! Index 2: the batch 
      integer(hik)                   :: nbatch_2e        ! Number of 2e integral batches
      integer(ik), allocatable       :: indices_2e(:,:)  ! Index 1: (i, j, k, l)
                                                         ! Index 2: the batch
      integer(hik), allocatable      :: offset_2e(:)     ! Starting position of each 2e integral batch within its I/O buffer
      integer(hik)                   :: nbatch_io        ! Number of I/O batches; may contain multiple integral batches
      integer(hik)                   :: current_batch_io ! I/O batch currently in the buffer
      integer(hik), allocatable      :: indices_io(:,:)  ! Index 1: First position in indices_2e; # of entries in indices_2e; total size in real(rk)
                                                         ! Index 2: the I/O batch
      integer(ik)                    :: io_unit          ! Only meaningful for scf_type='conventional'
      real(rk), allocatable          :: buffer_real(:)   ! Integral buffer (conventional) or all integrals (incore), active if ints_math=='real'
      real(xrk), allocatable         :: buffer_quad(:)   ! ditto, for ints_math=='quad'
      type(gam_structure), pointer   :: gam              ! Molecule descriptor associated with this integral cache
      character(len=20)              :: scf_type         ! Can be one of:
                                                         ! 'direct'        - compute all 2e AO integrals on the fly
                                                         ! 'incore'        - keep 2e integrals on disk
                                                         ! 'conventional'  - store integrals on disk
      character(len=20)              :: ints_math        ! Can be one of:
                                                         ! 'real'          - use standard precision
                                                         ! 'quad'          - use extended precision
    end type int2e_cache
    !
    interface swap_jikl
      module procedure swap_jikl_real
!*qd  module procedure swap_jikl_quad
    end interface swap_jikl
    interface swap_ijlk
      module procedure swap_ijlk_real
!*qd  module procedure swap_ijlk_quad
    end interface swap_ijlk
    interface swap_jilk
      module procedure swap_jilk_real
!*qd  module procedure swap_jilk_quad
    end interface swap_jilk
    interface swap_lkij
      module procedure swap_lkij_real
!*qd  module procedure swap_lkij_quad
    end interface swap_lkij
    interface swap_klij
      module procedure swap_klij_real
!*qd  module procedure swap_klij_quad
    end interface swap_klij
    interface swap_klji
      module procedure swap_klji_real
!*qd  module procedure swap_klji_quad
    end interface swap_klji
    interface swap_lkji
      module procedure swap_lkji_real
!*qd  module procedure swap_lkji_quad
    end interface swap_lkji
    !
    contains
    !
    !  Initialize basic data on 2e integrals.
    !  For scf_type='conventional' or 'incore', we'll also precalculate and store the integrals
    !
    subroutine prepare_2e(int2e,gam,scf_type,io_unit,iosize_2e,ints_math)
      type(int2e_cache), intent(inout)       :: int2e     ! Integral cache to be initialized
      type(gam_structure), target            :: gam       ! Molecule descriptor to be associated with the integral cache
      character(len=*), intent(in)           :: scf_type  ! Desired SCF type; Can be 'direct', 'incore', or 'conventional'
      integer(ik), intent(in)                :: io_unit   ! I/O unit used for creating scratch file; only relevant for 'conventional'
      integer(hik), intent(in)               :: iosize_2e ! Maximum size of the 2e integral buffer
      character(len=*), intent(in), optional :: ints_math ! Integral precision; either 'real' or 'quad'. The default if
                                                          ! the argument is missing is 'real'
      !
      integer(ik)               :: info_2e(4), alloc, ib, ios
      real(ark)                 :: mem_size, nr
      integer(ik)               :: i, j, k, l, cimax
      integer(hik)              :: nhi, intbatch, iobatch, iosize
      integer(hik), allocatable :: indices_io(:,:)   ! To avoid making two passes over integral batches
      !
      call TimerStart('Prepare 2E Buffer')
      !
      int2e%gam       => gam
      int2e%scf_type  =  scf_type
      int2e%io_unit   = io_unit
      int2e%ints_math = 'real'
      if (present(ints_math)) int2e%ints_math = ints_math
      !
      !  AO index batches
      ! 
      call gamess_2e_info('GET INFO',gam=int2e%gam,op_iparam=info_2e)
      int2e%nbatch_ao   = info_2e(1)
      int2e%maxbatch_ao = info_2e(2)
      !
      if (verbose>=0) then
        write (out,"('         Number of AO blocks: ',i0)") int2e%nbatch_ao
        write (out,"('Orbitals in largest AO block: ',i0)") int2e%maxbatch_ao
      end if
      !
      allocate (int2e%indices_ao(3,int2e%nbatch_ao),stat=alloc)
      if (alloc/=0) stop 'integral_tools%prepare_2e - no memory for batch table (A)'
      !
      remember_ao_batches: do ib=1,int2e%nbatch_ao
        info_2e(1) = ib
        call gamess_2e_info('GET BATCH INFO',gam=int2e%gam,op_iparam=info_2e)
        int2e%indices_ao(1,ib) = info_2e(3)  ! Initial position
        int2e%indices_ao(2,ib) = info_2e(2)  ! Size
        int2e%indices_ao(3,ib) = info_2e(4)  ! Atom batch belongs to
      end do remember_ao_batches
      !
      !  2e integral batches. We'll calculate the number of batches 
      !  in real arithmetics first, to make sure we won't blow an
      !  integer ...
      !
      nr       = real(int2e%nbatch_ao,kind=ark)
      mem_size = 0.125_ark * nr * (nr+1) * (nr**2+nr+2)
      if (mem_size>huge(1_hik)/8_hik) then
        write (out,"('Number of 2e AO integral batches too large for a long integer: ',f20.0)") mem_size
        stop 'integral_tools%prepare_2e - blown a long integer'
      end if
      nhi = int2e%nbatch_ao
      int2e%nbatch_2e = (nhi*(nhi+1)*(nhi**2+nhi+2))/8
      !
      if (verbose>=0) then
        write (out,"('     Number of AO 2e batches: ',i0)") int2e%nbatch_2e
      end if
      !
      allocate (int2e%indices_2e(4,int2e%nbatch_2e), &
                int2e%offset_2e(int2e%nbatch_2e), &
                indices_io(3,int2e%nbatch_2e),stat=alloc)
      if (alloc/=0) stop 'integral_tools%prepare_2e - no memory for batch table (B)'
      !
      !  Initialize batch indices and count the total number of integrals
      !
      intbatch = 0
      mem_size = 0
      iosize   = 0
      iobatch  = 1
      indices_io(1,1) = 1
      batch_l: do l=1,int2e%nbatch_ao
        batch_k: do k=1,l
          batch_j: do j=1,l
            cimax = j
            if (j==l) cimax = k
            batch_i: do i=1,cimax
              intbatch = intbatch + 1
              int2e%indices_2e(1:4,intbatch) = (/ i, j, k, l /)
              nhi = product(int(int2e%indices_ao(2,(/i,j,k,l/)),kind=hik))
              if (nhi>iosize_2e) stop 'integral_tools%prepare_2e - must increase iosize_2e!'
              if (iosize+nhi>iosize_2e .and. int2e%scf_type/='incore') then
                !
                !  I/O buffer is full; advance to the next buffer
                !
                iosize  = 0
                iobatch = iobatch + 1
                indices_io(1,iobatch) = intbatch
              endif
              indices_io(2,iobatch)     = intbatch - indices_io(1,iobatch) + 1
              int2e%offset_2e(intbatch) = iosize + 1
              iosize                    = iosize + nhi
              indices_io(3,iobatch)     = iosize
              !
              mem_size = mem_size + nhi
            end do batch_i
          end do batch_j
        end do batch_k
      end do batch_l
      if (intbatch/=int2e%nbatch_2e) stop 'integral_tools%prepare_2e - count error'
      !
      int2e%nbatch_io = iobatch
      int2e%current_batch_io = 0 ! Invalid; buffer is not filled
      !
      select case (int2e%scf_type)
        case default
          write (out,"('integral_tools%prepare_2e: SCF type ',a,' is not recognized')") trim(int2e%scf_type)
          stop 'integral_tools%prepare_2e - unknown scf_type'
        case ('incore')
          if (iosize/=mem_size) stop 'integral_tools%prepare_2e - integral count error (incore)'
        case ('conventional','direct')
          iosize = iosize_2e
      end select
      !
      if (verbose>=0) then
        write (out,"('   Number of AO 2e integrals: ',f0.0)") mem_size
        write (out,"('       Number of I/O batches: ',i0)") iobatch
        write (out,"('     Size of integral buffer: ',i0,' words')") iosize
        write (out,"('           Integral accuracy: ',a)") trim(int2e%ints_math)
        call flush(out)
      end if
      !
      allocate (int2e%indices_io(3,int2e%nbatch_io),stat=alloc)
      if (alloc/=0) stop 'integral_tools%prepare_2e - no memory for batch table (C)'
      int2e%indices_io(:,1:iobatch) = indices_io(:,1:iobatch)
      deallocate (indices_io)
      !
      select case (int2e%ints_math)
        case default ; stop 'integral_tools%prepare_2e - integral accuracy not recognized'
        case ('real') ; allocate (int2e%buffer_real(iosize),stat=alloc) ; mem_size = iosize * rk_bytes
        case ('quad') ; allocate (int2e%buffer_quad(iosize),stat=alloc) ; mem_size = iosize * xrk_bytes
      end select
      if (verbose >=0) then
        write (out,"('    Size of integrals buffer: ',f0.6,' GBytes')") mem_size / (1024._rk**3)
        call flush(out)
      end if
      if (alloc/=0) stop 'integral_tools%prepare_2e - no memory for integral buffer (D)'
      !
      !  For conventional or incore, we'll also need to calculate the integrals
      !
      if (int2e%scf_type=='incore' .or. int2e%scf_type=='conventional') then
        if (int2e%scf_type=='conventional') then
          open (int2e%io_unit,form='unformatted',action='readwrite',status='scratch',iostat=ios)
          if (ios/=0) then
            write (out,"('Error ',i0,' creating scratch file for 2E AO integrals')") ios
            stop 'integral_tools%prepare_2e - error creating AO 2E integral file'
          end if
          rewind (int2e%io_unit)
        end if
        !
        !  For in-core, there will be exactly one batch, so this also works
        !
        integral_buffers: do iobatch=1,int2e%nbatch_io
          call recalculate_ao2e(int2e,iobatch)
          if (int2e%scf_type=='conventional') then
            select case (int2e%ints_math)
              case ('real') ; write (int2e%io_unit) int2e%buffer_real(:int2e%indices_io(3,iobatch))
              case ('quad') ; write (int2e%io_unit) int2e%buffer_quad(:int2e%indices_io(3,iobatch))
            end select
          end if
        end do integral_buffers
      end if
      !
      call TimerStop('Prepare 2E Buffer')
      call TimerReport
    end subroutine prepare_2e
    !
    subroutine clear_2e(int2e)
      type (int2e_cache), intent(inout) :: int2e ! Integral cache to cleanup
      logical :: opened
      !
      call TimerStart('Clear 2E buffer')
      if (int2e%scf_type=='conventional') then
        inquire (int2e%io_unit,opened=opened)
        if (opened) close (int2e%io_unit,status='delete')
      end if
      if (allocated(int2e%indices_ao )) deallocate (int2e%indices_ao )
      if (allocated(int2e%indices_2e )) deallocate (int2e%indices_2e )
      if (allocated(int2e%offset_2e  )) deallocate (int2e%offset_2e  )
      if (allocated(int2e%indices_io )) deallocate (int2e%indices_io )
      if (allocated(int2e%buffer_real)) deallocate (int2e%buffer_real)
      if (allocated(int2e%buffer_quad)) deallocate (int2e%buffer_quad)
      nullify (int2e%gam)
      call TimerStop('Clear 2E buffer')
    end subroutine clear_2e
    !
    subroutine recalculate_ao2e(int2e,iobatch)
      type (int2e_cache), intent(inout) :: int2e   ! Integral cache to operate on
      integer(hik), intent(in)          :: iobatch
      !
      integer(hik)           :: ijkl_first, ijkl_count, ijkl_last, ijkl_size
      integer(hik)           :: ijkl, ijkl_position
      integer(ik)            :: bi(4), p0(4), sz(4)
      integer(ik)            :: mxsz, blk_sz
      integer(ik)            :: alloc_             ! Work around for Intel Fortran 14.0 bug
      real(rk), allocatable  :: ints_real(:,:,:,:)
      real(xrk), allocatable :: ints_quad(:,:,:,:)
      !
      call TimerStart('Calculate AO 2E integrals')
      ijkl_first = int2e%indices_io(1,iobatch)
      ijkl_count = int2e%indices_io(2,iobatch)
      ijkl_last  = ijkl_first + ijkl_count - 1
      ijkl_size  = int2e%indices_io(3,iobatch)
      mxsz       = int2e%maxbatch_ao
      ! write (out,"('+BW ',4(1x,i10))") iobatch, ijkl_first, ijkl_last, ijkl_size
      !
      !$omp parallel default(none) &
      !$omp& shared(int2e,mxsz,ijkl_first,ijkl_last) &
      !$omp& private(ints_real,ints_quad,alloc_,ijkl,ijkl_position,bi,p0,sz,blk_sz)
      select case (int2e%ints_math)
        case default; stop 'integral_tools%recalculate_ao2e - ints_math not recognized'
        case ('real') ; allocate (ints_real(mxsz,mxsz,mxsz,mxsz),stat=alloc_)
        case ('quad') ; allocate (ints_quad(mxsz,mxsz,mxsz,mxsz),stat=alloc_)
      end select
      if (alloc_/=0) stop 'integral_tools%fill_integral_buffer - integral buffer allocation failed'
      !$omp do schedule(dynamic)
      integral_batches: do ijkl=ijkl_first,ijkl_last
        ijkl_position = int2e%offset_2e(ijkl)
        bi            = int2e%indices_2e(:,ijkl)
        p0            = int2e%indices_ao(1,bi)
        sz            = int2e%indices_ao(2,bi)
        blk_sz        = product(sz)
        !
        !  Calculate this batch of integrals and put in into the integral buffer.
        !
        select case (int2e%ints_math)
          case ('real')
            call gamess_2e_integrals('AO 4C 1/R',ints_real,bi,a=int2e%gam,b=int2e%gam,c=int2e%gam,d=int2e%gam)
            int2e%buffer_real(ijkl_position:ijkl_position+blk_sz-1) = reshape(ints_real(:sz(1),:sz(2),:sz(3),:sz(4)),(/blk_sz/))
          case ('quad')
            call gamess_2e_integrals('AO 4C 1/R',ints_quad,bi,a=int2e%gam,b=int2e%gam,c=int2e%gam,d=int2e%gam)
            int2e%buffer_quad(ijkl_position:ijkl_position+blk_sz-1) = reshape(ints_quad(:sz(1),:sz(2),:sz(3),:sz(4)),(/blk_sz/))
        end select
      end do integral_batches
      !$omp end do
      if (allocated(ints_real)) deallocate (ints_real)
      if (allocated(ints_quad)) deallocate (ints_quad)
      !$omp end parallel
      int2e%current_batch_io = iobatch
      call TimerStop('Calculate AO 2E integrals')
    end subroutine recalculate_ao2e
    !
    !  Load (or recalculate) the next batch of the 2e integrals over the AOs
    !
    subroutine fetch_ao2e_batch(int2e,iobatch)
      type (int2e_cache), intent(inout) :: int2e  ! Integral cache to operate on
      integer(hik), intent(in) :: iobatch         ! Integral batch to process. This parameter is
                                                  ! largerly for consistency: We are expecting to process
                                                  ! all integrals sequentially, in the order of the I/O batches.
      !
      if (int2e%current_batch_io==iobatch) return ! Nothing to see here
      if (iobatch<0 .or. iobatch>int2e%nbatch_io) stop 'integral_tools%fetch_ao2e_batch - bad batch index'
      call TimerStart('Fetch AO 2E integrals')
      !
      select case (int2e%scf_type)
        case default
          stop 'integral_tools%fetch_ao2e_batch - impossible scf_type'
        case ('conventional')
          if (iobatch==1) then 
            !
            !  We are starting from the beginning
            !
            rewind (int2e%io_unit)
            int2e%current_batch_io = 0
          end if
          if (iobatch/=int2e%current_batch_io+1) stop 'integral_tools%fetch_ao2e_batch - nonsequential integral batch?!'
          int2e%current_batch_io = iobatch
          select case (int2e%ints_math)
            case ('real') ; read (int2e%io_unit) int2e%buffer_real(:int2e%indices_io(3,iobatch))
            case ('quad') ; read (int2e%io_unit) int2e%buffer_quad(:int2e%indices_io(3,iobatch))
          end select
        case ('direct')
          call recalculate_ao2e(int2e,iobatch)
      end select
      !
      call TimerStop('Fetch AO 2E integrals')
    end subroutine fetch_ao2e_batch
    !
    !  The routines below apply index permutations to integral blocks, so
    !  that we do not need to compute any of the blocks more than once.
    !  The only tricky bit are the conditional expressions used to cull
    !  redundant permutations.
    !
    !  The swap_* functions are a little verbose and repetitive; Unfortunately, 
    !  I do not see a more elegant way of writing them :(
    !
    !  Serial versions of the swapping routines (standard real kind)
    !
    function swap_jikl_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_jikl_common.f90'
    end function swap_jikl_real
    !
    function swap_ijlk_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_ijlk_common.f90'
    end function swap_ijlk_real
    !
    function swap_jilk_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_jilk_common.f90'
    end function swap_jilk_real
    !
    function swap_klij_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_klij_common.f90'
    end function swap_klij_real
    !
    function swap_klji_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_klji_common.f90'
    end function swap_klji_real
    !
    function swap_lkij_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_lkij_common.f90'
    end function swap_lkij_real
    !
    function swap_lkji_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_lkji_common.f90'
    end function swap_lkji_real
    !
    !  Serial versions of the swapping routines (quad real kind)
    !
    function swap_jikl_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_jikl_common.f90'
    end function swap_jikl_quad
    !
    function swap_ijlk_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_ijlk_common.f90'
    end function swap_ijlk_quad
    !
    function swap_jilk_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_jilk_common.f90'
    end function swap_jilk_quad
    !
    function swap_klij_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_klij_common.f90'
    end function swap_klij_quad
    !
    function swap_klji_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_klji_common.f90'
    end function swap_klji_quad
    !
    function swap_lkij_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_lkij_common.f90'
    end function swap_lkij_quad
    !
    function swap_lkji_quad(bi,p0,sz,p2,s2,a2e,a22) result(go)
      real(xrk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
      real(xrk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      !
      include 'integral_tools_swap_lkji_common.f90'
    end function swap_lkji_quad
    !
    !  Stripped-down versions; no integral swapping
    !
    function swap_jikl_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter  :: code(4) = (/2,1,3,4/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ci/=cj)
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_jikl_s
    !
    function swap_ijlk_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/1,2,4,3/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ck/=cl)
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_ijlk_s
    !
    function swap_jilk_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/2,1,4,3/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ci/=cj) .and. (ck/=cl)
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_jilk_s
    !
    function swap_klij_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/3,4,1,2/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_klij_s
    !
    function swap_klji_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/3,4,2,1/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ci/=cj) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_klji_s
    !
    function swap_lkij_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/4,3,1,2/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ck/=cl) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_lkij_s
    !
    function swap_lkji_s(bi,p0,sz,p2,s2) result(go)
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter :: code(4) = (/4,3,2,1/)
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ci/=cj) .and. (ck/=cl) .and. ((ci/=ck).or.(cj/=cl)) .and. ((ci/=cl).or.(cj/=ck))
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
    end function swap_lkji_s
  end module integral_tools
