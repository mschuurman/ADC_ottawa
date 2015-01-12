!   subroutine g_matrix_mixed(int2e,rho,gmat)
      type(int2e_cache), intent(inout)  :: int2e     ! 4-centre 2-electron integral descriptor; caller is responsible for initialization
!     complex(xrk), intent(in)          :: rho (:,:) ! Density matrix
!     complex(xrk), intent(out)         :: gmat(:,:) ! 2-electron contribution to the Fock matrix
      !
      !  Code below this line should not use explicit real kinds; it must refer to the kinds of 
      !  the input arguments, so that it can be reused in different-precision contexts.
      !
      integer(hik)                          :: iobatch
      integer(hik)                          :: ijkl_first, ijkl_count, ijkl_last
      integer(hik)                          :: ijkl, ijkl_position
      !
      integer(ik)                           :: bi(4), p0(4), sz(4)
      integer(ik)                           :: nao_spin, nao
      integer(ik)                           :: mxsz, blk_sz
      integer(ik)                           :: alloc_               ! Work-around for Intel Fortran 14.0 bug
      complex(kind(gmat)), allocatable      :: glocal(:,:)
      real(rk), allocatable                 :: ints_real(:,:,:,:)
      real(xrk), allocatable                :: ints_quad(:,:,:,:)
      !
      call TimerStart('2e Fock contribution')
      !
      nao_spin = size(rho,dim=1)
      nao      = nao_spin / 2
      mxsz     = int2e%maxbatch_ao
      !
      !  A bit of sanity/consistency checking
      !
      if (2*nao/=nao_spin .or. size(rho,dim=2)/=nao_spin .or. size(gmat,dim=1)/=nao_spin .or. size(gmat,dim=2)/=nao_spin) then
        stop 'fock_tools%g_matrix - inconsistent input matrices'
      end if
      !
      gmat = 0
      !
      !  The outermost loop is over the I/O batches. Since we use the same buffer
      !  for all I/O, this part can't run in parallel
      !
      iobatch_loop: do iobatch=1,int2e%nbatch_io
        call fetch_ao2e_batch(int2e,iobatch)      ! This could do a number of things depending on scf_type; this is the only bit
                                                  ! which actually cares about the specific scf_type.
        ijkl_first = int2e%indices_io(1,iobatch)  ! First 2e integral batch in this buffer
        ijkl_count = int2e%indices_io(2,iobatch)  ! Number of 2e integral batches in this buffer
        ijkl_last  = ijkl_first + ijkl_count - 1  ! Last 2e integral batch in the buffer
        ! write (out,"('+BR ',4(1x,i10))") iobatch, ijkl_first, ijkl_last, int2e%indices_io(3,iobatch)
        !
        !  Accumulate contributions to the G matrix, using integrals currently in the I/O buffer.
        !  This bit is embarassingly parallel.
        !
        !$omp parallel default(none) &
        !$omp& shared(mxsz,nao,nao_spin,ijkl_first,ijkl_last,int2e,rho,gmat) &
        !$omp& private(ints_real,ints_quad,glocal,alloc_,ijkl,ijkl_position,bi,p0,sz,blk_sz)
        select case (int2e%ints_math)
          case default; stop 'fock_tools%g_matrix - ints_type is not recognized'
          case ('real') ; allocate (ints_real(mxsz,mxsz,mxsz,mxsz),glocal(nao_spin,nao_spin),stat=alloc_)
          case ('quad') ; allocate (ints_quad(mxsz,mxsz,mxsz,mxsz),glocal(nao_spin,nao_spin),stat=alloc_)
        end select
        if (alloc_/=0) stop 'fock_tools%g_matrix - out of memory for parallel buffers'
        glocal = 0
        !$omp do schedule(dynamic)
        integral_batches: do ijkl=ijkl_first,ijkl_last
          ijkl_position = int2e%offset_2e(ijkl)       ! Position of this integral batch in the I/O buffer
          bi            = int2e%indices_2e(:,ijkl)    ! AO batch indices (spinless)
          p0            = int2e%indices_ao(1,bi)      ! Initial AO indices (spinless; 1-nao)
          sz            = int2e%indices_ao(2,bi)      ! Number of AOs in the block
          blk_sz        = product(sz)                 ! Total number of 2e integrals in the block
          !
          !  Accumulate contributions to the 2e-part of the Fock matrix
          !
          select case (int2e%ints_math)
            case ('real')
              ints_real(:sz(1),:sz(2),:sz(3),:sz(4)) = reshape(int2e%buffer_real(ijkl_position:ijkl_position+blk_sz-1),sz)
              call g_matrix_block(nao,bi,p0,sz,mxsz,ints_real,rho,glocal)
            case ('quad')
              ints_quad(:sz(1),:sz(2),:sz(3),:sz(4)) = reshape(int2e%buffer_quad(ijkl_position:ijkl_position+blk_sz-1),sz)
              call g_matrix_block(nao,bi,p0,sz,mxsz,ints_quad,rho,glocal)
          end select
        end do integral_batches
        !$omp end do nowait
        !$omp critical
        !
        gmat = gmat + glocal
        !$omp end critical
        deallocate (glocal)
        if (allocated(ints_real)) deallocate(ints_real)
        if (allocated(ints_quad)) deallocate(ints_quad)
        !$omp barrier
        !$omp end parallel
      end do iobatch_loop
      !
      call TimerStop('2e Fock contribution')
!   end subroutine g_matrix_mixed
