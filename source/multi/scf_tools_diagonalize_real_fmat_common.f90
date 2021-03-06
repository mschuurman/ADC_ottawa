!   !
!   !  Solve generalized eigenvalue problem, real symmetric matrix
!   !
!   subroutine diagonalize_fmat_rr(smhalf,sphalf,fmat,nmo_null,eps_geev,use_block_diag,mos,mo_energy,use_block_diag)
!     real(rk), intent(in)          :: smhalf   (:,:)   ! S**(-0.5)
!     real(rk), intent(in)          :: sphalf   (:,:)   ! S**(+0.5)
!     real(rk), intent(in)          :: fmat     (:,:)   ! Fock matrix
!     integer(ik), intent(in)       :: nmo_null         ! Expected size of the null-space
!     real(rk), intent(in)          :: eps_geev         ! Null-space threshold
!     real(rk), intent(out)         :: mos      (:,:)   ! (Right) molecular orbitals
!     real(rk), intent(out)         :: mo_energy(:)     ! Molecular orbital energy
!     logical, intent(in), optional :: use_block_diag   ! The default is .true.
      !
      !  Keep things kind-generic below this line
      !
      integer(ik) :: nmo
      logical     :: block
      !
      call TimerStart('Eigenvalue problem (real)')
      !
      nmo = size(fmat,dim=1)
      block = .true.
      if (present(use_block_diag)) block = use_block_diag
      !
      !  Construct S^{-1/2} F S^{-1/2}
      !
      mos(:,:) = mt_matmul(smhalf,mt_matmul(fmat,smhalf))
      mos = real(0.5_xrk,kind(mos)) * (mos + transpose(mos))  ! Symmetrize it!
      if (block) then
        call block_syev(mos,mo_energy)
      else
        call lapack_syev(mos,mo_energy)
      end if
      ! call bt_eigenvectors(nmo,mo_energy,mos,eps_geev) ! LAPACK already arranges Hermitian eigenvectors the way we want
      !
      !  Construct eigenvectors
      !
      mos(:,:) = mt_matmul(smhalf,mos(:,:))
      !
      !  Move null orbitals to the end of the MO table; this is the last time 
      !  we are looking at them.
      !
      call shift_null_mos
      !
      if (verbose>=1) then
        write (out,"(/t5,'MO energies:')")
        write (out,"(5(1x,f12.6,2x,12x))") mo_energy 
        write (out,"()")
      end if
      !
      call TimerStop('Eigenvalue problem (real)')
      !
      contains
      !
      !  Move null-space MOs to the end of the MO table
      !
      subroutine shift_null_mos
        real(kind(mos))            :: mo_norms  (nmo)
        integer(ik)                :: mo_reorder(nmo)
        integer(ik)                :: imo, iout
        real(kind(mos)), parameter :: eps = real(0.1_xrk,kind(mos))
        !
        call TimerStart('Shift null MOs')
        mo_norms  = sum(mt_matmul(sphalf,mos(:,:))**2,dim=1)  ! Expression valid for real MOs only
        !
        iout = 0
        pick_nonzero: do imo=1,nmo
          if (abs(mo_norms(imo))<=eps) cycle pick_nonzero
          iout = iout + 1
          mo_reorder(iout) = imo
        end do pick_nonzero
        if (iout/=nmo-nmo_null) then
          write (out,"('scf_tools%shift_null_mos: Expected ',i0,' non-null MOs, but got ',i0,'. Ooops.')") nmo-nmo_null, iout
          stop 'scf_tools%shift_null_mos - count error (1)'
        end if
        !
        pick_zero: do imo=1,nmo
          if (abs(mo_norms(imo))>eps) cycle pick_zero
          iout = iout + 1
          mo_reorder(iout) = imo
        end do pick_zero
        if (iout/=nmo) then
          write (out,"('scf_tools%shift_null_mos: Expected ',i0,' total MOs, but got ',i0,'. Ooops.')") nmo, iout
          stop 'scf_tools%shift_null_mos - count error (2)'
        end if
        !
        mos       = mos(:,mo_reorder)
        mo_energy = mo_energy(mo_reorder)
        !
        call TimerStop('Shift null MOs')
      end subroutine shift_null_mos
!   end subroutine diagonalize_fmat_rr
