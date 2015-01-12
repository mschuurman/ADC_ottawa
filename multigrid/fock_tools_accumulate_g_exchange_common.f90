!   subroutine accumulate_g_exchange_real(p0,sz,mxsz,a2e,rho,glocal)
      integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
      integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
      integer(ik), intent(in)    :: mxsz          ! Largest possible AO block
!     real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
!     complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
!     complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      integer(ik)           :: pe(4)          ! Last global index
      integer(ik)           :: i, j, k, l     ! Index in the local integral buffer
      integer(ik)           :: gj, gk         ! Index in the global array
      complex(kind(glocal)) :: s
      complex(kind(glocal)) :: bex(mxsz,mxsz) ! Exchange accumulation buffer
      !
      !  Exchange contribution
      !
      xchg_l: do l=1,sz(4)
        xchg_i: do i=1,sz(1)
          s = 0
          xchg_j: do j=1,sz(2)
            gj = p0(2)+j-1
            xchg_k: do k=1,sz(3)
              gk = p0(3)+k-1
              s = s + rho(gk,gj) * a2e(i,j,k,l)
            end do xchg_k
          end do xchg_j
          bex(i,l) = s
        end do xchg_i
      end do xchg_l
      !
      !  Update local copy of G matrix
      !
      pe = p0 + sz - 1
      glocal(p0(1):pe(1),p0(4):pe(4)) = glocal(p0(1):pe(1),p0(4):pe(4)) - bex(:sz(1),:sz(4))
!   end subroutine accumulate_g_exchange_real
