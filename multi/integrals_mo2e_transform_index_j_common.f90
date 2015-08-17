!   subroutine transform_index_j(nao,nmo,mo_j,buf_kl,buf_jkl)
      integer(ik), intent(in)  :: nao              ! Number of spin-less atomic orbitals
      integer(ik), intent(in)  :: nmo(:)           ! Number of MOs for each index
!     complex(rk), intent(in)  :: mo_j(:,:)        ! Second-index MO coefficients to transform over
!     complex(rk), intent(in)  :: buf_kl(:,:,:)    ! Integrals transformed over the last two indices; L index is fixed
!     complex(rk), intent(out) :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
!                                                  ! The fourth index is spin of the J component.
      !
      integer(ik) :: k
      !
      if (size(buf_jkl,dim=1)/=size(buf_kl,dim=1)) stop 'integrals_mo2e%transform_index_j - dimension mismatch (A)'
      if (size(buf_jkl,dim=3)/=size(buf_kl,dim=3)) stop 'integrals_mo2e%transform_index_j - dimension mismatch (B)'
      if (size(buf_jkl,dim=2)/=nmo(2))             stop 'integrals_mo2e%transform_index_j - dimension mismatch (C)'
      if (size(buf_kl, dim=2)/=nao)                stop 'integrals_mo2e%transform_index_j - dimension mismatch (D)'
      if (size(buf_jkl,dim=4)/=2)                  stop 'integrals_mo2e%transform_index_j - dimension mismatch (E)'
      if (size(mo_j,   dim=1)/=2*nao)              stop 'integrals_mo2e%transform_index_j - dimension mismatch (F)'
      if (size(mo_j,   dim=2)/=nmo(2))             stop 'integrals_mo2e%transform_index_j - dimension mismatch (G)'
      !
      !$omp parallel do default(none) private(k) shared(nao,nmo,buf_kl,buf_jkl,mo_j)
      transform_k: do k=1,nmo(3)
        buf_jkl(:,:,k,1) = matmul(buf_kl(:,:,k),mo_j(  :nao,:))
        buf_jkl(:,:,k,2) = matmul(buf_kl(:,:,k),mo_j(nao+1:,:))
      end do transform_k
      !$omp end parallel do
!   end subroutine transform_index_j
