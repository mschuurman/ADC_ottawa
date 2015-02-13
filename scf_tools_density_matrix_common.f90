!   !
!   !  This is not the most efficient implementation. The "correct" way
!   !  is to pre-multiply by square root of occupation numbers, and call 
!   !  ZHERK BLAS3 routine.
!   !
!   subroutine density_matrix_r(mo_occ,mos,rho)
!     real(rk), intent(in)     :: mo_occ(:)   ! MO occupation numbers
!     complex(rk), intent(in)  :: mos(:,:,:)  ! Left and right orbitals
!     complex(rk), intent(out) :: rho(:,:)
      !
      integer(ik) :: nmo, nao_spin
      integer(ik) :: mu, nu, i
      !
      call TimerStart('AO Density matrix')
      nmo      = size(mo_occ)
      nao_spin = size(mos,dim=1)
      if (size(mos,dim=2)/=nmo .or. size(mos,dim=3)/=2 .or. &
          size(rho,dim=1)/=nao_spin .or. size(rho,dim=2)/=nao_spin) then
        stop 'scf_tools%density_matrix - inconsistent input arrays'
      end if
      !
      rho = 0
      rho_mos: do i=1,nmo
        if (mo_occ(i)<=0) cycle rho_mos
        !$omp parallel do default(none) private(mu,nu) shared(mos,mo_occ,rho,i,nao_spin)
        rho_nu: do nu=1,nao_spin
          rho_mu: do mu=1,nao_spin
            !
            !  Since we dealing with non-self-adjoint operators, we half
            !  left eigenvectors [mos(:,:,1)] and right eigenvectors [mos(:,:,2)].
            !  In the more familiar self-adjoint case, we would have had
            !  mos(:,:,1) = conjg(mos(:,:,2))
            !
            rho(mu,nu) = rho(mu,nu) + mo_occ(i)*mos(mu,i,1)*mos(nu,i,2)
          end do rho_mu
        end do rho_nu
        !$omp end parallel do
      end do rho_mos
      call TimerStop('AO Density matrix')
!   end subroutine density_matrix_r
