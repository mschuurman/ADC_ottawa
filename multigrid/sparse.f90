!
!  02 Jan 2007 - Modified to reduce number of allocate/deallocate calls,
!                at the cost of higher memory usage. The old behavior
!                can be recovered by setting allocate_slack and 
!                allocate_extra to zero.
!
!  Handling of sparse matrices. This module implements only the functions
!  needed for matrix exponentiation in evaluation of partitition function,
!  and is likely not appropriate for any other purpose
!
!  Matrices are restricted to be real and symmetric
!
module sparse
  use accuracy
  use timer
  private
  public SparseMatrixT
  public SparseCreateUnitMatrix, SparseDestroyMatrix, SparseExtractDiagonal
  public SparseScale, SparseScreen, SparseHamiltonianMultiply, SparseStatistics
  !
  integer(ik), parameter :: verbose        = 0
  integer(ik), parameter :: allocate_slack = 20 ! Max. number of unused memory
                                                ! entries before reallocation
  integer(ik), parameter :: allocate_extra = 10 ! Number of additional words to
                                                ! allocate, should be around half
                                                ! of allocate_slack
  !
  type SparseRowT
    integer(ik)          :: nel           ! Number of non-zero elements
    integer(ik), pointer :: pos(:)        ! Indices of non-zero elements
    real(rk), pointer    :: x  (:)        ! Values of non-zero elements
  end type SparseRowT
  !
  type SparseMatrixT
    private                                    ! It is an opaque type
    real(ark)                 :: max_allocated ! Maximum number of elements allocated
    real(ark)                 :: max_used      ! Maximum number of elements used
    integer(ik)               :: nrows         ! Number of rows in the matrix
    type(SparseRowT), pointer :: rows(:)       ! Sparse rows of the matrix. All rows
                                               ! have to be present, even if they happen 
                                               ! to be empty
  end type SparseMatrixT
  !
  contains
  !
  !  Externally visible subroutines
  !
  function SparseCreateUnitMatrix(nrows) result(m)
    integer(ik), intent(in) :: nrows
    type(SparseMatrixT)     :: m
    !
    integer(ik) :: ir
    !
    call TimerStart('SparseCreateUnitMatrix')
    !
    !  Allocate space for the matrix
    !
    call preallocateSparse(m,nrows,1+allocate_extra)
    !
    !  Fill the elements
    !
    !$omp parallel do private(ir)
    fill_rows: do ir=1,nrows
      m%rows(ir)%nel    = 1
      m%rows(ir)%pos(1) = ir
      m%rows(ir)%x  (1) = 1._rk
    end do fill_rows
    !$omp end parallel do
    m%max_allocated = nrows * real(1+allocate_extra,kind=ark)
    m%max_used      = nrows
    call TimerStop('SparseCreateUnitMatrix')
  end function SparseCreateUnitMatrix
  !
  subroutine SparseDestroyMatrix(m)
    type(SparseMatrixT), intent(inout) :: m
    !
    integer(ik) :: ir, alloc
    !
    call TimerStart('SparseDestroyMatrix')
    !$omp parallel do private (ir,alloc)
    destroy_rows: do ir=1,m%nrows
      deallocate (m%rows(ir)%pos,m%rows(ir)%x,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating row ',i10,'(of ',i10,')')") &
               alloc, ir, m%nrows
        stop 'sparse%SparseDestroyMatrix - deallocate (1)'
      end if
    end do destroy_rows
    !$omp end parallel do
    !
    deallocate (m%rows,stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' deallocating row table')") alloc
      stop 'sparse%SparseDestroyMatrix - deallocate (2)'
    end if
    m%rows => NULL()
    call TimerStop('SparseDestroyMatrix')
  end subroutine SparseDestroyMatrix
  !
  subroutine SparseExtractDiagonal(m,d)
    type(SparseMatrixT), intent(in) :: m
    real(rk), intent(out)           :: d(:)
    !
    integer(ik) :: ir, ind
    !
    call TimerStart('SparseExtractDiagonal')
    if (size(d)<m%nrows) then
      write (out,"('Extract diagonal of ',i10,'-row matrix to ',i10,'-element array?!')") &
             m%nrows, size(d)
      stop 'sparse%SparseExtractDiagonal - too small'
    end if
    !
    !$omp parallel do private(ir,ind) schedule(guided)
    extract_row: do ir=1,m%nrows
      ind = locateElement(ir,m%rows(ir))
      if (ind<=0) then
        d(ir) = 0._rk
      else
        d(ir) = m%rows(ir)%x(ind)
      end if
    end do extract_row
    !$omp end parallel do
    call TimerStop('SparseExtractDiagonal')
  end subroutine SparseExtractDiagonal
  !
  subroutine SparseScale(m,scl)
    type(SparseMatrixT), intent(inout) :: m
    real(rk), intent(in)               :: scl
    !
    integer(ik) :: ir, nel
    !
    call TimerStart('SparseScale')
    !$omp parallel do private(ir,nel) schedule(guided)
    scale_row: do ir=1,m%nrows
      nel = m%rows(ir)%nel
      m%rows(ir)%x(:nel) = scl*m%rows(ir)%x(:nel)
    end do scale_row
    !$omp end parallel do
    call TimerStop('SparseScale')
  end subroutine SparseScale
  !
  subroutine SparseScreen(m,threshold)
    type(SparseMatrixT), intent(inout) :: m
    real(rk), intent(in)               :: threshold
    !
    integer(ik)          :: ir       ! Row index
    real(ark)            :: total    ! Non-zero elements upon entry
    real(ark)            :: left     ! Non-zero elements after screening
    real(ark)            :: allocated! Non-zero elements after screening
    integer(ik)          :: nel_src  ! Non-zero elements before screening
    integer(ik)          :: nel      ! Screened non-zero elements in each row
    integer(ik), pointer :: pos(:)   ! Placeholder for position array
    real(rk), pointer    :: x  (:)   ! Placeholder for value array
    integer(ik)          :: alloc
    integer(ik)          :: ip, op   ! Input/Output positions in compaction
    logical              :: realloc  ! Reallocation is needed
    !
    call TimerStart('SparseScreen')
    total     = 0
    left      = 0
    allocated = 0
    !$omp parallel do schedule(guided) reduction(+:total,left,allocated) &
    !$omp&            private(ir,nel_src,nel,pos,x,ip,op,alloc,realloc)
    screen_row: do ir=1,m%nrows
      nel_src   = m%rows(ir)%nel
      total     = total + nel_src
      nel       = count(abs(m%rows(ir)%x(:nel_src))>=threshold)
      left      = left + nel
      if (nel==nel_src) then
        allocated = allocated + size(m%rows)
        cycle screen_row
      end if
      !
      !  Row needs compaction, decide between reallocation and in-place
      !
      pos => m%rows(ir)%pos
      x   => m%rows(ir)%x
      realloc = (nel+allocate_slack) < size(x)
      !
      if (realloc) then 
        allocate (m%rows(ir)%pos(nel+allocate_extra),m%rows(ir)%x(nel+allocate_extra),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i10,' reallocating row to size ',i10)") alloc, nel
          stop 'sparse%SparseScreen - allocate'
        end if
      end if
      !
      op = 0
      copy_row: do ip=1,nel_src
        if (abs(x(ip))>=threshold) then
          op = op + 1
          !
          !  Without this test, we trigger a compiler bug in Intel's 64-bit
          !  compiler, where the sense of the comparison is inverted.
          !
          ! if (op>nel) then
          !   stop 'sparse%SparseScreen - compiler bug: inconsistent >= operator'
          ! end if
          m%rows(ir)%pos(op) = pos(ip)
          m%rows(ir)%x  (op) = x  (ip)
        end if
      end do copy_row
      if (op/=nel) then
        write (out,"('Count error - expected ',i10,' elements, but got ',i10)") nel, op
        stop 'sparse%SparseScreen - count'
      end if
      !
      m%rows(ir)%nel = nel
      allocated = allocated + size(m%rows)
      !
      if (realloc) then
        deallocate (pos,x,stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i10,' releasing old row')") alloc
          stop 'sparse%SparseScreen - deallocate'
        end if
      end if
    end do screen_row
    !$omp end parallel do 
    if (verbose>=1) then
      write (out,"(5x,'Screen: full = ',f15.0,' input = ',f15.0,' output = ',f15.0,' allocated = ',f15.0)") &
             real(m%nrows,kind=rk)**2, total, left, allocated
    end if
    m%max_allocated = max(allocated,m%max_allocated)
    m%max_used      = max(left,m%max_used)
    call TimerStop('SparseScreen')
  end subroutine SparseScreen
  !
  !  Multiply a sparse matrix in-place from the right by an implicitly 
  !  defined sparse Hamiltonian matrix. The kinetic energy part is 
  !  given by the neighbour list, and the list of distinct matrix elements. 
  !  The format is:
  !
  !  neigbour(1:6,ipt) - gives indices of the neighbouring points,
  !                      displaced in -/+ directions along each coordinate
  !  kmat(1:6)         - gives matrix elements for each of the displacement
  !                      directions.
  !  diag(ipt)         - gives the diagonal value of the Hamiltonian matrix,
  !                      which is always non-zero (except by accident)
  !
  subroutine SparseHamiltonianMultiply(mat,neighbour,kmat,diag,screen)
 !$ use OMP_LIB
    type(SparseMatrixT), intent(inout) :: mat            ! Input/output matrix
    integer(ik), intent(in)            :: neighbour(:,:) ! Neighbour list 
    real(rk), intent(in)               :: kmat(1:6)      ! Off-diagonal matrix elements
    real(rk), intent(in)               :: diag(:)        ! Diagonal matrix elements
    real(rk), intent(in)               :: screen         ! Magnitude of non-zero elements
                                                         ! to preserve in the output
    !
    integer(ik)       :: q_max           ! Maximum possible length of the queue
    integer(ik)       :: max_threads     ! Maximum possible number of execution threads
    integer(ik)       :: irow, erow      ! Initial and final rows in a chunk
    integer(ik)       :: row_step        ! Chunk size to use when allocating rows
    real(ark)         :: allocated, used ! Allocation statistics
    !
    call TimerStart('SparseHamiltonianMultiply')
    !
    !   Scan the input matrix, to determine the size of the longest row. 
    !   Due to the sparse-ness of the Hamiltonian matrix, the accumulation queue 
    !   is limited to 7*number of non-zero elements in a row.
    !
    q_max = 7*maxval(mat%rows(:)%nel)
    !
    if (verbose>1) then
      write (out,"('Multiply: Max. queue length = ',i10)") q_max
      call flush(out)
    end if
    !
    !  Decide on the number of threads
    !
    max_threads = 1
 !$ max_threads = omp_get_max_threads()
    row_step    = (mat%nrows+8*max_threads-1)/(8*max_threads)
    !
    !  Parallel loop. We do it in a slightly ugly way to avoid allocating
    !  and deallocating queue memory for each row
    !
    !$omp parallel do private(irow,erow) schedule(guided)
    parallel_multiply: do irow=1,mat%nrows,row_step
      erow = min(mat%nrows,irow+row_step-1)
      call doHamiltonianMultiply(irow,erow,mat,neighbour,kmat,diag,q_max,screen)
    end do parallel_multiply
    !$omp end parallel do
    !
    !  Update memory allocation statistics
    !
    allocated = 0 ; used = 0
    !$omp parallel do private(irow) reduction(+:allocated,used)
    memory_stats: do irow=1,mat%nrows
      allocated = allocated + size(mat%rows(irow)%x)
      used      = used      + mat%rows(irow)%nel
    end do memory_stats
    !$omp end parallel do
    mat%max_allocated = max(mat%max_allocated,allocated)
    mat%max_used      = max(mat%max_used,used)
    !
    call TimerStop('SparseHamiltonianMultiply')
  end subroutine SparseHamiltonianMultiply
  !
  !  Return memory allocation statistics
  !
  subroutine SparseStatistics(mat,max_allocated,max_used)
    type(SparseMatrixT), intent(in)  :: mat           ! Matrix to return statistics for
    real(ark), intent(out), optional :: max_allocated ! Maximum amount of memory allocated
    real(ark), intent(out), optional :: max_used      ! Maximum amount of memory used
    !
    if (present(max_allocated)) max_allocated = mat%max_allocated*(ik_bytes+rk_bytes)
    if (present(max_used))      max_used      = mat%max_used     *(ik_bytes+rk_bytes)
  end subroutine SparseStatistics
  !
  !  Auxiliary routines, used only inside the module
  !
  !  Multiply by Hamiltonian matrix from the right using accumulation queue
  !
  recursive subroutine doHamiltonianMultiply(start_row,end_row,mat,neighbour,kmat,diag,q_max,screen)
    integer(ik)                        :: start_row      ! First row to process in this block
    integer(ik)                        :: end_row        ! Last row to process in this block
    type(SparseMatrixT), intent(inout) :: mat            ! Input/output matrix
    integer(ik), intent(in)            :: neighbour(:,:) ! Neighbour list
    real(rk), intent(in)               :: kmat(1:6)      ! Non-zero off-diagonal matrix elements
    real(rk), intent(in)               :: diag(:)        ! Diagonal matrix elements
    integer(ik), intent(in)            :: q_max          ! Maximum possible queue length
    real(rk), intent(in)               :: screen         ! Magnitude of non-zero elements
                                                         ! to preserve in the output
    !
    integer(ik), allocatable  :: q_key(:)      ! Space for the keys queue
    real(rk), allocatable     :: q_val(:)      ! Space for the values queue
                                               ! q_key & q_val were originally on the stack;
                                               ! this causes all kinds of memory problems
                                               ! in large OpenMP runs.
    integer(ik)               :: q_p           ! Next position in the accumulation queue
    integer(ik)               :: nel           ! Number of non-zero elements
    type(SparseRowT), pointer :: row           ! Source row element
    integer(ik)               :: isr           ! Source matrix row
    integer(ik)               :: isc_p         ! Source matrix column pointer
    integer(ik)               :: isc           ! Source matrix column (position)
    integer(ik)               :: alloc, i
    !
    allocate (q_key(q_max),q_val(q_max),stat=alloc)
    if (alloc/=0) then
      write (out,"('Allocate ',i10,' reals+integers: ',i6)") q_max, alloc
      stop 'sparse%doHamiltonianMultiply - out of memory for the sort queue'
    end if
    !
    scan_source_rows: do isr=start_row,end_row
      if (verbose>2) then
        write (out,"('Multiply: Processing row ',i10)") isr
        call flush(out)
      end if
      row => mat%rows(isr)
      q_p = 1
      scan_source_columns: do isc_p=1,row%nel
        isc = row%pos(isc_p)
        !
        !  We have a non-zero element at (isr,isc). The potential columns are 
        !  (isr,isc), plus the neighbours of the element (isc,isc) in the 
        !  kinetic matrix
        !
        q_key(q_p)         = isc
        q_val(q_p)         = diag(isc)*row%x(isc_p)
        q_key(q_p+1:q_p+6) = neighbour(1:6,isc)
        q_val(q_p+1:q_p+6) = kmat(1:6)*row%x(isc_p)
        q_p                = q_p + 7
      end do scan_source_columns
      !
      if (verbose>2) then
        write (out,"(/'Before sorting:')")
        write (out,"(' qpos= ',i10,' key= ',i10,' val= ',g14.7)") &
               (i, q_key(i), q_val(i),i=1,q_p-1)
        call flush(out)
      end if
      !
      !  Collapse the queue, and copy the result back to the matrix
      !
      nel = mergeSort2(q_p-1,q_key,q_val)
      !
      !  Screen the output queue
      !
      nel = screenQueue(nel,q_key,q_val,screen)
      !
      if (verbose>2) then
        write (out,"(/'After sorting:')")
        write (out,"(' qpos= ',i10,' key= ',i10,' val= ',g14.7)") &
               (i, q_key(i), q_val(i),i=1,nel)
        call flush(out)
      end if
      !
      !  Decide on whether reallocation is needed
      !
      if ( ((nel+allocate_slack) < size(row%pos)) .or. &
           ( nel                 > size(row%pos)) ) then
        deallocate (row%pos,row%x,stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i10,' deallocating old row arrays')") alloc
          stop 'sparse%doHamiltonianMultiply - deallocate'
        end if
        allocate (row%pos(nel+allocate_extra),row%x(nel+allocate_extra),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i10,' allocating ',i10,'-element row arrays')") &
                 alloc, nel
          stop 'sparse%doHamiltonianMultiply - allocate'
        end if
      end if
      !
      row%nel       = nel
      row%pos(:nel) = q_key(1:nel)
      row%x  (:nel) = q_val(1:nel)
    end do scan_source_rows
    !
    deallocate (q_key,q_val,stat=alloc)
    if (alloc/=0) then
      write (out,"('Deallocate: ',i6)") alloc
      stop 'sparse%doHamiltonianMultiply - deallocate failed'
    end if
  end subroutine doHamiltonianMultiply
  !
  !  Preallocate space for a sparse matrix
  !
  subroutine preallocateSparse(m,nrows,nel)
    type(SparseMatrixT), intent(out) :: m
    integer(ik), intent(in)          :: nrows ! Number of rows/columns
    integer(ik), intent(in)          :: nel   ! Number of non-zero elements in each row
    !
    integer(ik) :: ir
    integer(ik) :: alloc
    !
    m%nrows = nrows
    allocate(m%rows(nrows),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating sparse matrix with ',i10,' rows')") alloc, nrows
      stop 'sparse%preallocateSparse - no memory (1)'
    end if
    !
    !$omp parallel do private (ir,alloc)
    fill_rows: do ir=1,nrows
      m%rows(ir)%nel  = nel
      allocate(m%rows(ir)%pos(nel),m%rows(ir)%x(nel),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' allocating ',i8,'-element row ',i10,' (total ',i10,')')") &
               alloc, nel, ir, nrows
      stop 'sparse%preallocateSparse - no memory (2)'
      end if
    end do fill_rows
    !$imp end parallel do
  end subroutine preallocateSparse
  !
  recursive function screenQueue(n_in,key,val,screen) result(n_out)
    integer(ik), intent(in)    :: n_in   ! Number of elements to sort
    integer(ik), intent(inout) :: key(:) ! Keys
    real(rk), intent(inout)    :: val(:) ! Values; elements with identical keys are summed up
    real(rk), intent(in)       :: screen ! Screening threshold
    integer(ik)                :: n_out  ! Number of surviving keys
    !
    integer(ik)                :: in
    !
    n_out = 0
    screen_loop: do in=1,n_in
      if (abs(val(in))<screen) cycle screen_loop
      n_out = n_out + 1
      key(n_out) = key(in)
      val(n_out) = val(in)
    end do screen_loop
  end function screenQueue
  !
  !  Recursive 2-way merge sort with element reduction
  !
  recursive function mergeSort2(n_in,key,val) result(n_out)
    integer(ik), intent(in)    :: n_in   ! Number of elements to sort
    integer(ik), intent(inout) :: key(:) ! Keys
    real(rk), intent(inout)    :: val(:) ! Values; elements with identical keys are summed up
    integer(ik)                :: n_out  ! Number of surviving keys
    !
    integer(ik) :: k_left (2+n_in/2)  ! Temporary buffers for keys
    integer(ik) :: k_right(2+n_in/2)  
    real(rk)    :: v_left (2+n_in/2)  ! Temporary buffers for values
    real(rk)    :: v_right(2+n_in/2)
    integer(ik) :: c_left, c_right    ! Initial number of left/right keys
    integer(ik) :: n_left, n_right    ! Number of surviving keys in left/right arrays
    integer(ik) :: p_left, p_right    ! Positions of the left/right sections
    !
    if (n_in<=1) then
      n_out = n_in
      return ! Already sorted
    end if
    !
    !  Partition the input array
    !
    p_left  = 1
    c_left  = n_in/2
    p_right = p_left + c_left
    c_right = n_in - c_left
    !
    if (verbose>3) then
      write (out,"('mergeSort2: Asked to sort ',i10,' elements')") n_in
      write (out,"((t12,8(i10,1x)))") key(1:n_in)
      write (out,"('mergeSort2:  Left: ',i10,' elements from ',i10)") c_left, p_left
      write (out,"('mergeSort2: Right: ',i10,' elements from ',i10)") c_right, p_right
      call flush(out)
    end if
    !
    !  Copy out the keys and values
    !
    k_left (1:c_left ) = key( p_left:p_left+c_left-1  )
    v_left (1:c_left ) = val( p_left:p_left+c_left-1  )
    k_right(1:c_right) = key(p_right:p_right+c_right-1)
    v_right(1:c_right) = val(p_right:p_right+c_right-1)
    !
    !  Merge sort left and right arrays
    !
    n_left  = mergeSort2(c_left, k_left, v_left )
    n_right = mergeSort2(c_right,k_right,v_right)
    !
    if (verbose>3) then
      write (out,"('mergeSort2: Survivors left: ',i10,' right: ',i10)") n_left, n_right
      call flush(out)
    end if
    !
    !  Add impossibly large, key values as sentinels at the right
    !
    k_left (n_left+1)  = huge(1_ik)
    k_right(n_right+1) = huge(1_ik)
    !
    !  Merge pre-sorted arrays. At this point, left and right sections are already
    !  sorted, and contain no duplicate keys.
    !
    n_out   = 0
    p_left  = 1
    p_right = 1
    merge_left_right: do 
      if (k_left(p_left)==k_right(p_right)) then
        !
        !  Have we reached the end of both arrays? Thanks to sentinels, only
        !  one index need to be tested.
        !
        if (p_left>n_left) exit merge_left_right
        !
        !  Key are identical - merge and consume both values
        !
        n_out      = n_out + 1
        key(n_out) = k_left(p_left)
        val(n_out) = v_left(p_left) + v_right(p_right)
        p_left     = p_left + 1
        p_right    = p_right + 1
        cycle merge_left_right
      end if
      if (k_left(p_left)<k_right(p_right)) then
        !
        !  Left key is smaller - copy the left value
        !
        n_out      = n_out + 1
        key(n_out) = k_left(p_left)
        val(n_out) = v_left(p_left)
        p_left     = p_left + 1
      else
        !
        !  Right key is smaller - copy the right value
        !
        n_out      = n_out + 1
        key(n_out) = k_right(p_right)
        val(n_out) = v_right(p_right)
        p_right    = p_right + 1
      end if
    end do merge_left_right
    !
    if (verbose>3) then
      write (out,"('mergeSort2: Final merge gave ',i10)") n_out
    end if
  end function mergeSort2
  !
  !  Locate column index within a row
  !
  function locateElement(ip,row) result(mp)
    integer(ik), intent(in)      :: ip   ! Desired index
    type(SparseRowT), intent(in) :: row  ! Row to search in
    integer(ik)                  :: mp   ! Index within the table, or 0 if not present
    !
    integer(ik) :: lp, rp  ! Left and right search intervals
    integer(ik) :: tp      ! Index at the midpoint
    !
    if (verbose>3) then
      write (out,"('locateElement: Looking for key ',i10)") ip
      write (out,"((t5,'keys: ',8(i10,1x)))") row%pos
    end if
    !
    mp = 0 
    lp = 1 
    rp = row%nel
    !
    !  Range checks for the search
    !
    if (row%nel==0) return
    if (row%pos(lp)>ip) return
    if (row%pos(rp)<ip) return
    !
    !  Desired position -may- be within the interval; do the search
    !
    bisection_search: do while(rp>=lp)
      mp = (lp+rp)/2
      tp = row%pos(mp)
      if (verbose>3) then
        write (out,"(t10,'left= ',i10,' right= ',i10,' centre= ',i10,' kval= ',i10)") &
               lp, rp, mp, tp
      end if
      if (tp==ip) return ! Position found
      if (tp<ip) then
        lp = mp + 1   ! Desired index -may- be in the right interval
      else
        rp = mp - 1   ! Desired index -may- be in the left interval
      end if
    end do bisection_search
    !
    !  We reach here only if the desired element is not present
    !
    mp = 0
  end function locateElement
end module sparse
