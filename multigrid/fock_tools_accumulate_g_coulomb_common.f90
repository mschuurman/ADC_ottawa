!   subroutine accumulate_g_coulomb_real(p0,sz,mxsz,a2e,rho,glocal)
      integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
      integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
      integer(ik), intent(in)    :: mxsz          ! Largest possible AO block
!     real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
!     complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
!     complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      integer(ik)           :: pe(4)          ! Last global index
      integer(ik)           :: i, j, k, l     ! Index in the local integral buffer
      integer(ik)           :: gk, gl         ! Index in the global array (density matrix)
      complex(kind(glocal)) :: s
      complex(kind(glocal)) :: bcl(mxsz,mxsz) ! Coulomb accumulation buffer
      !
      !  Coulomb contribution
      !
      coul_j: do j=1,sz(2)
        coul_i: do i=1,sz(1)
          s = 0
          coul_l: do l=1,sz(4)
            gl = p0(4)+l-1
            coul_k: do k=1,sz(3)
              gk = p0(3)+k-1
              s = s + rho(gk,gl) * a2e(i,j,k,l)
            end do coul_k
          end do coul_l
          bcl(i,j) = s
        end do coul_i
      end do coul_j
      !
      !  Update local copy of G matrix
      !
      pe = p0 + sz - 1
      glocal(p0(1):pe(1),p0(2):pe(2)) = glocal(p0(1):pe(1),p0(2):pe(2)) + bcl(:sz(1),:sz(2))
!   end subroutine accumulate_g_coulomb_real
