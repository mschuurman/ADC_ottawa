!   subroutine accumulate_g_real(nao,p0,sz,mxsz,a2e,rho,glocal)
      integer(ik), intent(in)    :: nao           ! Number of spin-less AOs
      integer(ik), intent(in)    :: p0(:)         ! First AO in the block; SPINLESS
      integer(ik), intent(in)    :: sz(:)         ! Number of AOs
      integer(ik), intent(in)    :: mxsz          ! Largest possible AO block
!     real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
!     complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
!     complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
!     !
      integer(ik) :: spin_ij, spin_kl  ! Spin labels associated with integral indices; either 0 (alpha) or 1 (beta)
      integer(ik) :: shift(4)          ! Spin-related shift in orbital indices
      !
      scan_spin_k: do spin_kl=0,1
        scan_spin_i: do spin_ij=0,1
          shift = nao * (/ spin_ij, spin_ij, spin_kl, spin_kl /)
          call accumulate_g_coulomb (p0+shift,sz,mxsz,a2e,rho,glocal)
          call accumulate_g_exchange(p0+shift,sz,mxsz,a2e,rho,glocal)
        end do scan_spin_i
      end do scan_spin_k
!   end subroutine accumulate_g_real
