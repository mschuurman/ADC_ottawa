module block_diag
!
!  Simplistic block-diagonalization routines
!
!  Note that this set of routines is not complete: only the data types 
!  We encountered so far are implemented. Sorry.
!
!  2013 Nov 02 - Added block_geev routine
!  2013 Jul 16 - Added sorting of eigenvalues / singular values
!  2014 Jun ?? - Added quad-precision
!
  use accuracy
  use timer
  use lapack
  use printing
  use sort_tools
  implicit none
  private
  public block_syev, block_gesvd, block_geev
  !
  interface block_syev
     module procedure block_syev_real
!*qd module procedure block_syev_quad
  end interface block_syev
  !
  interface block_geev
     module procedure block_zgeev_real
!*qd module procedure block_zgeev_quad
  end interface block_geev
  !
  integer(ik), parameter :: verbose = 0

  contains

  subroutine block_syev_real(h,e,eps)
    real(rk), intent(inout)        :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                              ! Out: Eigenvectors
    real(rk), intent(out)          :: e(:)    ! Out: Eigenvalues
    real(rk), intent(in), optional :: eps     ! Threshold for treating elements as zero;
                                              ! Missing or negative value will be replaces by the
                                              ! default threshold of 1000*spacing(maxval(abs(h))
    !
    include 'block_diag_syev_common.f90'
  end subroutine block_syev_real
!
  subroutine block_syev_quad(h,e,eps)
    real(xrk), intent(inout)        :: h(:,:) ! In:  symmetric matrix to be diagonalized
                                              ! Out: Eigenvectors
    real(xrk), intent(out)          :: e(:)   ! Out: Eigenvalues
    real(xrk), intent(in), optional :: eps    ! Threshold for treating elements as zero;
                                              ! Missing or negative value will be replaces by the
                                              ! default threshold of 1000*spacing(maxval(abs(h))
    !
    include 'block_diag_syev_common.f90'
  end subroutine block_syev_quad
  !
  !  This routine is modeled on lapack_zgeev2
  !
  subroutine block_zgeev_real(h,e,eps)
    complex(rk), intent(inout)     :: h(:,:,:)! In:  h(:,:,1) = general matrix to be diagonalized
                                              ! Out: h(:,:,1) = Left eigenvectors
                                              !      h(:,:,2) = Right eigenvectors
    complex(rk), intent(out)       :: e(:)    ! Out: Eigenvalues. Note that eigenvalues are NOT sorted
    real(rk), intent(in), optional :: eps     ! Threshold for treating elements as zero;
                                              ! Missing or negative value will be replaces by the
                                              ! default threshold of 1000*spacing(maxval(abs(h))
    !
    include "block_diag_zgeev_common.f90"
  end subroutine block_zgeev_real

  subroutine block_zgeev_quad(h,e,eps)
    complex(xrk), intent(inout)     :: h(:,:,:)
                                              ! In:  h(:,:,1) = general matrix to be diagonalized
                                              ! Out: h(:,:,1) = Left eigenvectors
                                              !      h(:,:,2) = Right eigenvectors
    complex(xrk), intent(out)       :: e(:)   ! Out: Eigenvalues. Note that eigenvalues are NOT sorted
    real(xrk), intent(in), optional :: eps    ! Threshold for treating elements as zero;
                                              ! Missing or negative value will be replaces by the
                                              ! default threshold of 1000*spacing(maxval(abs(h))
    !
    include "block_diag_zgeev_common.f90"
  end subroutine block_zgeev_quad

  subroutine block_gesvd(a,s,u,vth,eps)
    real(rk), intent(inout)        :: a  (:,:) ! In:  Matrix to be decomposed
                                               ! Out: Matrix destroyed
    real(rk), intent(out)          :: s  (:)   ! Out: Singular values
    real(rk), intent(out)          :: u  (:,:) ! Out: Left singular vectors
    real(rk), intent(out)          :: vth(:,:) ! Out: Right singular vectors, transposed & conjugated
                                               ! The overall result is A = U S VTH
    real(rk), intent(in), optional :: eps      ! Threshold for treating elements as zero;
                                               ! Missing or negative value will be replaces by the
                                               ! default threshold of 1000*spacing(maxval(abs(a))
    !
    integer(ik)              :: nrow, ncol     ! Number of rows and columns in matrix a
    integer(ik)              :: nsing          ! Number of singular values in matrix a (=min(nrow,ncol))
    integer(ik)              :: ns_out         ! Number of singular values produced so far
    integer(ik)              :: blk_nrow       ! Number of rows in a block
    integer(ik)              :: blk_ncol       ! ditto, columns
    integer(ik)              :: blk_nsing      ! smaller of blk_nrow and blk_ncol
    integer(ik)              :: ir, ic         ! Row/column indices
    integer(ik)              :: alloc
    real(rk)                 :: thr            ! zero threshold
    logical, allocatable     :: row_avail(:)   ! Rows still not processed
    logical, allocatable     :: col_avail(:)   ! ditto, columns
    integer(ik), allocatable :: blk_rows(:)    ! Rows within the current block
    integer(ik), allocatable :: blk_cols(:)    ! ditto, columns
    integer(ik), allocatable :: order(:)       ! Singular-value ordering table
    real(rk), allocatable    :: tmp_a(:,:)     ! Temporary used for decomposing the connected block
    real(rk), allocatable    :: tmp_u(:,:)     ! ditto
    real(rk), allocatable    :: tmp_vth(:,:)   ! ditto
    real(rk), allocatable    :: tmp_s(:)       ! ditto
    logical                  :: block_grown    ! Indicates that the block has grown, so that an extra
                                               ! connections pass may be required
    !
    call TimerStart('block_gesvd')
    nrow  = size(a,dim=1)
    ncol  = size(a,dim=2)
    nsing = min(nrow,ncol)
    if (nrow/=size(u,dim=1) .or. nsing/=size(u,dim=2) .or. &
        nsing/=size(s) .or. ncol/=size(vth,dim=1) .or. nsing/=size(vth,dim=2) ) then
      stop 'block_diag%block_gesvd - dimensions mismatch'
    end if
    !
    thr = 1000._rk * spacing(maxval(abs(a)))
    if (present(eps)) then
      if (eps>=0) thr = eps
    end if
    if (verbose>=1) then
      write (out,"('Block-SVD zero threshold is ',g20.10)") thr
    end if
    !
    allocate (row_avail(nrow),col_avail(ncol),blk_rows(nrow),blk_cols(ncol),order(nsing),stat=alloc)
    if (alloc/=0) stop 'block_diag%block_gesvd - allocation failed (1)'
    !
    !  In order to simplify block extraction, mark all zero rows and columns as
    !  unavailable. These rows/columns are part of the null-space; we'll deal
    !  with them later, once non-zero singular values have been determined.
    !
    row_avail = any(abs(a(:,:))>thr,dim=2)
    col_avail = any(abs(a(:,:))>thr,dim=1)
    ns_out    = 0
    !
    block_loop: do while(any(row_avail) .and. any(col_avail))
      blk_nrow = 0 
      blk_ncol = 0
      !
      !  The number of singular values is the same as the smaller dimension,
      !  so this is what we'll need to scan along. The logic of block scanning
      !  below is a little complicated: since we are working with a general
      !  matrix, adding a row and collected columns may introduce connection
      !  to a row which was previously unconnected to the block. To catch this
      !  situation, we have to restart the scan each time we've added some rows
      !
      if (nrow<=ncol) then
        block_grown = .true.
        restart_scan_rows: do while (block_grown)
          block_grown = .false.
          scan_unused_rows: do ir=1,nrow
            if (.not.row_avail(ir)) cycle scan_unused_rows
            if (row_connected(ir)) then
              call add_row(ir)
              block_grown = .true.
            end if
          end do scan_unused_rows
        end do restart_scan_rows
      else
        block_grown = .true.
        restart_scan_cols: do while (block_grown)
          block_grown = .false.
          scan_unused_cols: do ic=1,ncol
            if (.not.col_avail(ic)) cycle scan_unused_cols
            if (col_connected(ic)) then
              call add_col(ic)
              block_grown = .true.
            end if
          end do scan_unused_cols
        end do restart_scan_cols
      end if
      if (blk_nrow==0 .or. blk_ncol==0) stop 'block_diag%block_gesvd - logic error: no block found'
      !
      if (verbose>=1) then
        write (out,"('Found connected block of dimension ',i8,' x ',i8)") blk_nrow, blk_ncol
        write (out,"(('rows: ',10(1x,i8)))") blk_rows(1:blk_nrow)
        write (out,"(('cols: ',10(1x,i8)))") blk_cols(1:blk_ncol)
        write (out,"()")
      end if
      !
      !  Extract the connected block and decompose it
      !
      blk_nsing = min(blk_nrow,blk_ncol)
      allocate (tmp_a(blk_nrow,blk_ncol),tmp_u(blk_nrow,blk_nsing),tmp_s(blk_nsing),tmp_vth(blk_nsing,blk_ncol),stat=alloc)
      if (alloc/=0) stop 'block_diag%block_gesvd - allocation failed (2)'
      tmp_a(:,:) = a(blk_rows(1:blk_nrow),blk_cols(1:blk_ncol))
      !
      if (verbose>=2) then
        write (out,"('The connected block is: ')")
        call print_matrix(tmp_a,12_ik,'g12.4')
      end if
      !
      !  A bit of extra debugging here: zero out the connected block in a(:,:),
      !  then make sure that there are no additional connections for the columns
      !  in the current block. Once this is done, zero the remainder of the
      !  rows/columns
      !
      a(blk_rows(1:blk_nrow),blk_cols(1:blk_ncol)) = 0
      if (any(abs(a(blk_rows(1:blk_nrow),:))>1.1_rk*thr) .or. any(abs(a(:,blk_cols(1:blk_ncol)))>1.1_rk*thr)) then
        write (out,"(/'Logic error in block_gesvd')") 
        write (out,"('threshold = ',g24.15)") thr
        write (out,"('blk_nrow = ',i8,' blk_ncol = ',i8)") blk_nrow, blk_ncol
        write (out,"('blk_rows = '/(12(1x,i8)))") blk_rows(1:blk_nrow)
        write (out,"('blk_cols = '/(12(1x,i8)))") blk_cols(1:blk_ncol)
        write (out,"('coupling row-block is:')")
        call print_matrix(a(blk_rows(1:blk_nrow),:),16,'e16.8')
        write (out,"('coupling column-block is:')")
        call print_matrix(a(:,blk_cols(1:blk_ncol)),16,'e16.8')
        stop 'block_diag%block_gesvd - logic error: found non-zeros in an uncoupled block'
      end if
      a(blk_rows(1:blk_nrow),:) = 0
      a(:,blk_cols(1:blk_ncol)) = 0
      !
      !  Decompose the sub-block
      !
      call lapack_svd(tmp_a,tmp_s,tmp_u,tmp_vth)
      if (verbose>=1) then
        write (out,"('Block singular values are:')")
        write (out,"(5(1x,g20.12))") tmp_s
        write (out,"()")
      end if
      if (verbose>=2) then
        write (out,"('Block left singular vectors are:')")
        call print_matrix(tmp_u,20,'g20.12')
        write (out,"('Block transposed right singular vectors are:')")
        call print_matrix(tmp_vth,20,'g20.12')
      end if
      !
      !  Store singular values and vectors into the output matrices
      !
      u(:,ns_out+1:ns_out+blk_nsing) = 0
      vth(ns_out+1:ns_out+blk_nsing,:) = 0
      s(  ns_out+1:ns_out+blk_nsing) = tmp_s
      u(blk_rows(1:blk_nrow),ns_out+1:ns_out+blk_nsing) = tmp_u
      vth(ns_out+1:ns_out+blk_nsing,blk_cols(1:blk_ncol)) = tmp_vth
      !
      deallocate (tmp_a,tmp_s,tmp_u,tmp_vth)
      !
      ns_out = ns_out + blk_nsing
    end do block_loop
    !
    !  At this point, all singular values not accounted for must be zero
    !  We'll pair them with zero singular vectors. Strictly speaking, this 
    !  is not the correct thing to do, since zeros are not proper singular 
    !  vectors. As long as we actually discard the null-space later on, 
    !  our choice won't do any harm.
    !
    if (ns_out<nsing) then
      if (verbose>=1) write (out,"('Adding ',i5,' zero singular values')") nsing-ns_out
      s  (  ns_out+1:nsing)   = 0._rk
      u  (:,ns_out+1:nsing)   = 0._rk
      vth(  ns_out+1:nsing,:) = 0._rk
    end if
    !
    !  Order singular values. order_keys() returns -ascending- order, while we require
    !  the descending order in this routine.
    !
    call order_keys(s,order)
    order = order(nsing:1:-1)
    s     = s  (  order)
    u     = u  (:,order)
    vth   = vth(  order,:)
    !
    deallocate (row_avail,col_avail,blk_rows,blk_cols,order)
    !
    call TimerStop('block_gesvd')
    !
    contains
    !
    !  Return true if row (ir) is connected to the rows already in the block
    !
    logical function row_connected(ir)
      integer(ik), intent(in) :: ir
      !
      if (blk_nrow==0) then  ! Any row is connected to an empty block!
        row_connected = .true.
      else
        row_connected = any(abs(a(ir,blk_cols(:blk_ncol)))>thr)
      end if
    end function row_connected
    !
    !  Return true if column (ic) is connected to the columns already in the block
    !
    logical function col_connected(ic)
      integer(ik), intent(in) :: ic
      !
      if (blk_ncol==0) then  ! Any column is connected to an empty block!
        col_connected = .true.
      else
        col_connected = any(abs(a(blk_rows(:blk_nrow),ic))>thr)
      end if
    end function col_connected
    !
    !  Add row to the block; we may also need to add a few columns
    !
    subroutine add_row(ir)
      integer(ik), intent(in) :: ir
      integer(ik)             :: ic
      !
      if (.not.row_avail(ir)) stop 'block_diag%block_gesvd%add_row - logical impossibility trap (1)'
      !
      blk_nrow           = blk_nrow + 1
      blk_rows(blk_nrow) = ir
      row_avail(ir)      = .false.
      !
      add_connected_columns: do ic=1,ncol
        if (abs(a(ir,ic))<=thr) cycle add_connected_columns
        if (any(blk_cols(:blk_ncol)==ic)) cycle add_connected_columns ! Already in the block
        if (.not.col_avail(ic)) stop 'block_diag%block_gesvd%add_row - logical impossibility trap (2)'
        !
        blk_ncol           = blk_ncol + 1
        blk_cols(blk_ncol) = ic
        col_avail(ic)      = .false.
      end do add_connected_columns
    end subroutine add_row
    !
    !  Add column to the block; we may also need to add a few rows
    !
    subroutine add_col(ic)
      integer(ik), intent(in) :: ic
      integer(ik)             :: ir
      !
      if (.not.col_avail(ic)) stop 'block_diag%block_gesvd%add_col - logical impossibility trap (1)'
      !
      blk_ncol           = blk_ncol + 1
      blk_cols(blk_ncol) = ic
      col_avail(ic)      = .false.
      !
      add_connected_rows: do ir=1,nrow
        if (abs(a(ir,ic))<=thr) cycle add_connected_rows
        if (any(blk_rows(:blk_nrow)==ir)) cycle add_connected_rows ! Already in the block
        if (.not.row_avail(ir)) stop 'block_diag%block_gesvd%add_col - logical impossibility trap (2)'
        !
        blk_nrow           = blk_nrow + 1
        blk_rows(blk_nrow) = ir
        row_avail(ir)      = .false.
      end do add_connected_rows
    end subroutine add_col

  end subroutine block_gesvd

end module block_diag
