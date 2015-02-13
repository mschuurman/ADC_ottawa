! subroutine block_zgeev(h,e,eps)
!   complex(rk), intent(inout)     :: h(:,:,:)! In:  h(:,:,1) = general matrix to be diagonalized
!                                             ! Out: h(:,:,1) = Left eigenvectors
!                                             !      h(:,:,2) = Right eigenvectors
!   complex(rk), intent(out)       :: e(:)    ! Out: Eigenvalues. Note that eigenvalues are NOT sorted
!   real(rk), intent(in), optional :: eps     ! Threshold for treating elements as zero;
!                                             ! Missing or negative value will be replaces by the
!                                             ! default threshold of 1000*spacing(maxval(abs(h))
    !
    integer(ik)                   :: sz            ! Number of rows in matrix h
    integer(ik)                   :: n_avail       ! Number of available rows in matrix h
    integer(ik)                   :: ir            ! Row index
    integer(ik)                   :: alloc
    logical, allocatable          :: available(:)  ! True for each row of h which has not been processed yet
    real(kind(h))                 :: thr           ! zero threshold
    integer(ik), allocatable      :: blk_rows(:)   ! List of rows in the current connected block
    integer(ik)                   :: blk_sz        ! Size of the current connected block
    complex(kind(h)), allocatable :: tmp_h(:,:,:)  ! Temporary used for diagonalizing the connected block
    complex(kind(h)), allocatable :: tmp_e(:)      ! ditto
    logical                       :: block_grown   ! .true. if block has grown on the current pass; do another pass
    !
    call TimerStart('block_geev (complex)')
    sz = size(h,dim=1)
    if (sz/=size(h,dim=2) .or. sz/=size(e)) stop 'block_diag%block_geev - dimensions mismatch'
    if (size(h,dim=3)/=2) stop 'block_diag%block_geev - third dimension of h is not 2?!'
    !
    thr = real(1e3_xrk,kind(thr)) * spacing(maxval(abs(h(:,:,1))))
    if (present(eps)) then
      if (eps>=0) thr = eps
    end if
    if (verbose>=1) then
      write (out,"('Block-diagonalization zero threshold is ',g20.10)") thr
    end if
    !
    allocate (available(sz),blk_rows(sz),stat=alloc)
    if (alloc/=0) stop 'block_diag%block_geev - allocation failed (1)'
    !
    available(:) = .true.
    n_avail = sz
    block_loop: do while(n_avail>0)
      blk_sz = 0
      grow_block: do
        block_grown = .false.
        find_first_unused: do ir=1,sz
          if (.not.available(ir)) cycle find_first_unused
          if (blk_sz==0) then
            !
            !  Find the first unused row
            !
            call add_row
          else
            !
            !  If this row is connected to any of the rows already in the list,
            !  add it as well. Make sure to check both row and column connections
            !
            if (any(abs(h(   blk_rows(1:blk_sz),ir,1))>thr) .or. &
                any(abs(h(ir,blk_rows(1:blk_sz)   ,1))>thr)) call add_row
          end if
        end do find_first_unused
        if (blk_sz==0) stop 'block_diag%block_geev - logic error: no block found' 
        if (.not.block_grown) exit grow_block
        ! Otherwise, we need another pass to make sure we have not missed any connected rows
      end do grow_block
      !
      if (verbose>=1) then
        write (out,"('Found connected block of dimension ',i8,' containing rows/columns:')") blk_sz
        write (out,"(10(1x,i8))") blk_rows(1:blk_sz)
        write (out,"()")
      end if
      !
      !  Extract the connected block
      !
      allocate (tmp_h(blk_sz,blk_sz,2),tmp_e(blk_sz),stat=alloc)
      if (alloc/=0) stop 'block_diag%block_geev - allocation failed (2)'
      tmp_h(:,:,1) = h(blk_rows(1:blk_sz),blk_rows(1:blk_sz),1)
      !
      if (verbose>=2) then
        write (out,"('The connected block is: ')")
        call print_matrix(tmp_h(:,:,1),12_ik,'g12.4')
      end if
      !
      !  A bit of extra debugging here: zero out the connected block in h(:,:),
      !  then make sure that there are no additional connections for the columns
      !  in the current block. Once this is done, zero the remainder of the
      !  rows/columns
      !
      h(blk_rows(1:blk_sz),blk_rows(1:blk_sz),1) = 0
      if (any(abs(h(  blk_rows(1:blk_sz),:,1))>real(1.1_xrk,kind(thr))*thr) .or. &
          any(abs(h(:,blk_rows(1:blk_sz),  1))>real(1.1_xrk,kind(thr))*thr)) then
        stop 'block_diag%block_geev - logic error: found non-zeros in an uncoupled block'
      end if
      h(  blk_rows(1:blk_sz),:,:) = 0
      h(:,blk_rows(1:blk_sz),  :) = 0
      !
      !  Diagonalize the sub-block
      !
      call lapack_geev(tmp_h,tmp_e)
      if (verbose>=1) then
        write (out,"('Block eigenvalues are:')")
        write (out,"(5(1x,g20.12))") tmp_e
        write (out,"()")
      end if
      if (verbose>=2) then
        write (out,"('Block left eigenvectors are:')")
        call print_matrix(tmp_h(:,:,1),20,'g20.12')
        write (out,"('Block right eigenvectors are:')")
        call print_matrix(tmp_h(:,:,2),20,'g20.12')
      end if
      !
      !  Store eigenvalues and eigenvectors back into the input/output matrix
      !
      e(blk_rows(1:blk_sz)) = tmp_e
      h(blk_rows(1:blk_sz),blk_rows(1:blk_sz),:) = tmp_h
      !
      deallocate (tmp_h,tmp_e)
    end do block_loop
    !
    !  There is no universal sensible order for complex eigenvalues, so leave them alone
    !
    deallocate (available,blk_rows)
    call TimerStop('block_geev (complex)')
    !
    contains
    !
    subroutine add_row
      blk_sz           = blk_sz + 1
      blk_rows(blk_sz) = ir
      available(ir)    = .false.
      n_avail          = n_avail - 1
      block_grown      = .true.
    end subroutine add_row

! end subroutine block_zgeev
