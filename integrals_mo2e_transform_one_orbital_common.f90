!   subroutine transform_one_orbital_rr(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
!     type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
!     integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
!     integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
!     complex(rk), intent(in)          :: mo_i(:,:)          ! First-index MOs
!     complex(rk), intent(in)          :: mo_j(:,:)          ! Second-index MOs
!     complex(rk), intent(in)          :: mo_k(:,:)          ! Third-index MOs
!     complex(rk), intent(in)          :: mo_l(:,:)          ! Current fourth-index MO(s)
!     complex(rk), intent(out)         :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L); 
!     complex(rk), intent(out)         :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals, single L
!     complex(rk), intent(out)         :: buf_kl (:,:,:)     ! ditto
!     complex(rk), intent(out)         :: buf_l  (:,:,:,:,:) ! fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
!     real(rk), intent(in)             :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      integer(hik)                       :: iobatch
      integer(hik)                       :: ijkl_first, ijkl_count, ijkl_last
      integer(hik)                       :: ijkl, ijkl_position
      integer(ik)                        :: bi(4), p0(4), sz(4)
      integer(ik)                        :: mxsz, alloc, blk_sz
      integer(ik)                        :: l, l_count
      real(kind(ints_kind)), allocatable :: ints(:,:,:,:)
      !
      !  The first pass is over I/O batches of 2e AO integrals
      !  Very little processing occurs in this pass; we'll try to run it serially for simplicity.
      !
      call TimerStart('AO>MO 2e pass L')
      l_count = size(mo_l,dim=2)
      if (l_count/=size(buffer,dim=4) .or. l_count>size(buf_l,dim=5)) then
        stop 'integrals_mo2e%transform_one_orbital - inconsistent buffer dimensions'
      end if
      !
      !$omp parallel do shared(buf_l,l_count) private(l)
      zero_buf_l: do l=1,l_count
        buf_l(:,:,:,:,l) = 0
      end do zero_buf_l
      !$omp end parallel do
      !
      !  The integral pass potentially involves I/O. We would like to reduce its cost,
      !  by sharing the integrals between multiple one-orbital transforms.
      !
      mxsz = int2e%maxbatch_ao
      allocate (ints(mxsz,mxsz,mxsz,mxsz),stat=alloc)
      if (alloc/=0) stop 'integrals_mo2e%transform_one_orbital - out of memory for integral buffers'
      iobatch_loop: do iobatch=1,int2e%nbatch_io
        call fetch_ao2e_batch(int2e,iobatch)      ! This could do a number of things depending on scf_type; this is the only bit
                                                  ! which actually cares about the specific scf_type.
        ijkl_first = int2e%indices_io(1,iobatch)  ! First 2e integral batch in this buffer
        ijkl_count = int2e%indices_io(2,iobatch)  ! Number of 2e integral batches in this buffer
        ijkl_last  = ijkl_first + ijkl_count - 1  ! Last 2e integral batch in the buffer
        !
        integral_batches: do ijkl=ijkl_first,ijkl_last
          ijkl_position = int2e%offset_2e(ijkl)       ! Position of this integral batch in the I/O buffer
          bi            = int2e%indices_2e(:,ijkl)    ! AO batch indices (spinless)
          p0            = int2e%indices_ao(1,bi)      ! Initial AO indices (spinless; 1-nao)
          sz            = int2e%indices_ao(2,bi)      ! Number of AOs in the block
          blk_sz        = product(sz)                 ! Total number of 2e integrals in the block
          !
          !  In the select case below, the type-changing version will never execute; however, it
          !  has to be present to make the code kind-agnostic. A good compiler should be able to
          !  delete the dead branch ...
          !
          select case (int2e%ints_math)
            case default; stop 'integrals_mo2e%transform_one_orbital - bad AO ints_math'
            case ('real') 
              ints(:sz(1),:sz(2),:sz(3),:sz(4)) = real(reshape(int2e%buffer_real(ijkl_position:ijkl_position+blk_sz-1),sz),kind(ints))
            case ('quad')
              ints(:sz(1),:sz(2),:sz(3),:sz(4)) = real(reshape(int2e%buffer_quad(ijkl_position:ijkl_position+blk_sz-1),sz),kind(ints))
          end select
          !
          !  Accumulate contributions to the partially transformed integrals
          !  The fourth-index MO is fixed: it's the L'th right eigenorbital
          !  Transformation is parallelized, so no parallel loop here
          !
          last_index_transform: do l=1,l_count
            call transform_ao_integral_block(bi,p0,sz,ints,nao,mo_l(:,l),buf_l(:,:,:,:,l))
          end do last_index_transform
        end do integral_batches
      end do iobatch_loop
      deallocate (ints)
      call TimerStop('AO>MO 2e pass L')
      !
      three_index_transform: do l=1,l_count
        !
        ! The remaining passes are trivial to run in parallel ...
        ! Third index transform runs over left MOs
        !
        call TimerStart('AO>MO 2e pass K')
        call transform_index_k(nao,nmo,mo_k,buf_l(:,:,:,:,l),buf_kl)
        call TimerStop('AO>MO 2e pass K')
        !
        !
        ! Second index transform runs over right MOs
        !
        call TimerStart('AO>MO 2e pass J')
        call transform_index_j(nao,nmo,mo_j,buf_kl,buf_jkl)
        call TimerStop('AO>MO 2e pass J')
        !
        !
        ! First index transform runs over left MOs
        !
        call TimerStart('AO>MO 2e pass I')
        call transform_index_i(nao,nmo,mo_i,buf_jkl,buffer(:,:,:,l))
        call TimerStop('AO>MO 2e pass I')
      end do three_index_transform
!   end subroutine transform_one_orbital_rr
