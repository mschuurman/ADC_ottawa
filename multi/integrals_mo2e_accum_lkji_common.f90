!   subroutine accum_lkji(p0,sz,a2e,mol,buf_l)
      integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
      integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
!     real(rk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
!     complex(rk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
!     complex(rk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
      !
      integer(ik) ::  i,  j,  k ! Indices in the local integral block
      integer(ik) :: gi, gj, gk ! Matching spin-less integral indices
      integer(ik) :: l0, le     ! First and last spinless AO index for the fourth index
      !
      l0 = p0(4)
      le = p0(4) + sz(4) - 1
      ! This inner loop is not the best place to parallelize; however, it is the easiest one
      !$omp parallel do default(none) private(k,gk,j,gj,i,gi) shared(sz,p0,buf_l,a2e,mol,l0,le)
      scan_k: do k=1,sz(3)
        gk = p0(3) + k - 1
        scan_j: do j=1,sz(2)
          gj = p0(2) + j - 1
          scan_i: do i=1,sz(1)
            gi = p0(1) + i - 1
            buf_l(gi,gj,gk) = buf_l(gi,gj,gk) + cmplx(sum(a2e(:sz(4),k,j,i)*mol(l0:le)),kind=kind(buf_l))
          end do scan_i
        end do scan_j
      end do scan_k
      !$omp end parallel do
!   end subroutine accum_lkji
