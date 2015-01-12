!   subroutine transform_index_i(nao,nmo,mo_i,buf_jkl,buffer)
      integer(ik), intent(in)  :: nao              ! Number of spin-less atomic orbitals
      integer(ik), intent(in)  :: nmo(:)           ! Number of MOs for each index
!     complex(rk), intent(in)  :: mo_i(:,:)        ! First-index MO coefficient to transform over
!     complex(rk), intent(in)  :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
!                                                  ! The fourth index is spin of the J component.
!     complex(rk), intent(out) :: buffer(:,:,:)    ! Fully transformed MOs; L index is fixed
!     !
      integer(ik) :: k
      !
      if (size(buf_jkl,dim=2)/=size(buffer,dim=2)) stop 'integrals_mo2e%transform_index_i - dimension mismatch (A)'
      if (size(buf_jkl,dim=3)/=size(buffer,dim=3)) stop 'integrals_mo2e%transform_index_i - dimension mismatch (B)'
      if (size(buf_jkl,dim=1)/=nao)                stop 'integrals_mo2e%transform_index_i - dimension mismatch (C)'
      if (size(buffer, dim=1)/=nmo(1))             stop 'integrals_mo2e%transform_index_i - dimension mismatch (D)'
      if (size(buf_jkl,dim=4)/=2)                  stop 'integrals_mo2e%transform_index_i - dimension mismatch (E)'
      if (size(mo_i,   dim=1)/=2*nao)              stop 'integrals_mo2e%transform_index_i - dimension mismatch (F)'
      if (size(mo_i,   dim=2)/=nmo(1))             stop 'integrals_mo2e%transform_index_i - dimension mismatch (G)'
      !
      !$omp parallel do default(none) private(k) shared(nao,nmo,buf_jkl,buffer,mo_i)
      transform_k: do k=1,nmo(3)
        buffer(:,:,k) = cmplx(matmul(transpose(mo_i(  :nao,:)),buf_jkl(:,:,k,1)) &
                            + matmul(transpose(mo_i(nao+1:,:)),buf_jkl(:,:,k,2)),kind=kind(buffer))
      end do transform_k
      !$omp end parallel do
!   end subroutine transform_index_i
