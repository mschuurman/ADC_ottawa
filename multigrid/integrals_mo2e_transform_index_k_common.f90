!   subroutine transform_index_k(nao,nmo,mo_k,buf_l,buf_kl)
      integer(ik), intent(in)  :: nao            ! Number of spin-less atomic orbitals
      integer(ik), intent(in)  :: nmo(:)         ! Number of MOs for each index
!     complex(rk), intent(in)  :: mo_k(:,:)      ! Third-index MOs coefficients
!     complex(rk), intent(in)  :: buf_l(:,:,:,:) ! Integrals transformed over the last index; L index is fixed
!                                                ! The fourth index is spin of the L component.
!     complex(rk), intent(out) :: buf_kl(:,:,:)  ! Integrals transformed over the last two indices; L index is fixed
      !
      integer(ik) :: j
      !
      if (size(buf_kl,dim=1)/=size(buf_l,dim=1)) stop 'integrals_mo2e%transform_index_k - dimension mismatch (A)'
      if (size(buf_kl,dim=2)/=size(buf_l,dim=2)) stop 'integrals_mo2e%transform_index_k - dimension mismatch (B)'
      if (size(buf_kl,dim=3)/=nmo(3))            stop 'integrals_mo2e%transform_index_k - dimension mismatch (C)'
      if (size(buf_l, dim=3)/=nao)               stop 'integrals_mo2e%transform_index_k - dimension mismatch (D)'
      if (size(buf_l, dim=4)/=2)                 stop 'integrals_mo2e%transform_index_k - dimension mismatch (E)'
      if (size(mo_k,  dim=1)/=2*nao)             stop 'integrals_mo2e%transform_index_k - dimension mismatch (F)'
      if (size(mo_k,  dim=2)/=nmo(3))            stop 'integrals_mo2e%transform_index_k - dimension mismatch (G)'
      !
      !$omp parallel do default(none) private(j) shared(nao,buf_kl,buf_l,mo_k)
      transform_j: do j=1,nao
        buf_kl(:,j,:) = matmul(buf_l(:,j,:,1),mo_k(  :nao,:)) &
                      + matmul(buf_l(:,j,:,2),mo_k(nao+1:,:))
      end do transform_j
      !$omp end parallel do
!   end subroutine transform_index_k
