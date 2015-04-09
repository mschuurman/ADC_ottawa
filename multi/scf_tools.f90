!
!  Useful routines which occur again and again in SCF implementation
!
  module scf_tools
    use accuracy
    use timer
    use lapack
    use block_diag
    use biorthogonal_tools
    use matrix_tools
    implicit none
    private
    public st_invert_smat, st_diagonalize_fmat, st_density_matrix
    !
    interface st_invert_smat
      module procedure invert_smat_r
!*qd  module procedure invert_smat_q
    end interface st_invert_smat
    !
    interface st_diagonalize_fmat
      module procedure diagonalize_fmat_cr
!*qd  module procedure diagonalize_fmat_cq
      module procedure diagonalize_fmat_rr
!*qd  module procedure diagonalize_fmat_rq
    end interface st_diagonalize_fmat
    !
    interface st_density_matrix
      module procedure density_matrix_r
!*qd  module procedure density_matrix_q
    end interface st_density_matrix
    !
    integer(ik), parameter :: verbose    = 0           ! Level of output
    !
    contains
    !
    subroutine invert_smat_r(nmo_null,smat,smhalf,sphalf,use_block_diag,eps_smat)
      integer(ik), intent(out)       :: nmo_null        ! Number of null-space vectors
      real(rk), intent(inout)        :: smat(:,:)       ! In: Raw overlap matrix
                                                        ! Out: Overlap matrix with null-space vectors removed
      real(rk), intent(out)          :: smhalf(:,:)     ! S**(-0.5), with null-space removed
      real(rk), intent(out)          :: sphalf(:,:)     ! S**(+0.5), with null-space removed
      logical, intent(in), optional  :: use_block_diag  ! Default is .true.
      real(rk), intent(in), optional :: eps_smat        ! Desired null-space cutoff. Default is 0 (near-machice accuracy)

      include 'scf_tools_invert_smat_common.f90'
    end subroutine invert_smat_r
    !
    subroutine invert_smat_q(nmo_null,smat,smhalf,sphalf,use_block_diag,eps_smat)
      integer(ik), intent(out)        :: nmo_null        ! Number of null-space vectors
      real(xrk), intent(inout)        :: smat(:,:)       ! In: Raw overlap matrix
                                                         ! Out: Overlap matrix with null-space vectors removed
      real(xrk), intent(out)          :: smhalf(:,:)     ! S**(-0.5), with null-space removed
      real(xrk), intent(out)          :: sphalf(:,:)     ! S**(+0.5), with null-space removed
      logical, intent(in), optional   :: use_block_diag  ! Default is .true.
      real(xrk), intent(in), optional :: eps_smat        ! Desired null-space cutoff. Default is 0 (near-machice accuracy)

      include 'scf_tools_invert_smat_common.f90'
    end subroutine invert_smat_q
    !
    !  Solve generalized eigenvalue problem
    !
    subroutine diagonalize_fmat_cr(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy,use_block_diag)
      real(rk), intent(in)          :: smhalf   (:,:)   ! S**(-0.5)
      real(rk), intent(in)          :: sphalf   (:,:)   ! S**(+0.5)
      complex(rk), intent(in)       :: fmat     (:,:)   ! Fock matrix
      integer(ik), intent(in)       :: nmo_null         ! Expected size of the null-space
      real(rk), intent(in)          :: eps_geev         ! Null-space threshold
      complex(rk), intent(out)      :: mos      (:,:,:) ! Molecular orbitals
      complex(rk), intent(out)      :: mo_energy(:)     ! Molecular orbital energy
      logical, intent(in), optional :: use_block_diag   ! The default is .true.
      !
      include 'scf_tools_diagonalize_fmat_common.f90'
    end subroutine diagonalize_fmat_cr
    !
    subroutine diagonalize_fmat_cq(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy,use_block_diag)
      real(xrk), intent(in)         :: smhalf   (:,:)   ! S**(-0.5)
      real(xrk), intent(in)         :: sphalf   (:,:)   ! S**(+0.5)
      complex(xrk), intent(in)      :: fmat     (:,:)   ! Fock matrix
      integer(ik), intent(in)       :: nmo_null         ! Expected size of the null-space
      real(xrk), intent(in)         :: eps_geev         ! Null-space threshold
      complex(xrk), intent(out)     :: mos      (:,:,:) ! Molecular orbitals
      complex(xrk), intent(out)     :: mo_energy(:)     ! Molecular orbital energy
      logical, intent(in), optional :: use_block_diag   ! The default is .true.
      !
      include 'scf_tools_diagonalize_fmat_common.f90'
    end subroutine diagonalize_fmat_cq
    !
    subroutine diagonalize_fmat_rr(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy,use_block_diag)
      real(rk), intent(in)          :: smhalf   (:,:)   ! S**(-0.5)
      real(rk), intent(in)          :: sphalf   (:,:)   ! S**(+0.5)
      real(rk), intent(in)          :: fmat     (:,:)   ! Fock matrix
      integer(ik), intent(in)       :: nmo_null         ! Expected size of the null-space
      real(rk), intent(in)          :: eps_geev         ! Null-space threshold
      real(rk), intent(out)         :: mos      (:,:)   ! (Right) molecular orbitals
      real(rk), intent(out)         :: mo_energy(:)     ! Molecular orbital energy
      logical, intent(in), optional :: use_block_diag   ! The default is .true.
      !
      include 'scf_tools_diagonalize_real_fmat_common.f90'
    end subroutine diagonalize_fmat_rr
    !
    subroutine diagonalize_fmat_rq(smhalf,sphalf,fmat,nmo_null,eps_geev,mos,mo_energy,use_block_diag)
      real(xrk), intent(in)         :: smhalf   (:,:)   ! S**(-0.5)
      real(xrk), intent(in)         :: sphalf   (:,:)   ! S**(+0.5)
      real(xrk), intent(in)         :: fmat     (:,:)   ! Fock matrix
      integer(ik), intent(in)       :: nmo_null         ! Expected size of the null-space
      real(xrk), intent(in)         :: eps_geev         ! Null-space threshold
      real(xrk), intent(out)        :: mos      (:,:)   ! (Right) molecular orbitals
      real(xrk), intent(out)        :: mo_energy(:)     ! Molecular orbital energy
      logical, intent(in), optional :: use_block_diag   ! The default is .true.
      !
      include 'scf_tools_diagonalize_real_fmat_common.f90'
    end subroutine diagonalize_fmat_rq
    !
    subroutine density_matrix_r(mo_occ,mos,rho)
      real(rk), intent(in)     :: mo_occ(:)   ! MO occupation numbers
      complex(rk), intent(in)  :: mos(:,:,:)  ! Left and right orbitals
      complex(rk), intent(out) :: rho(:,:)
      !
      include 'scf_tools_density_matrix_common.f90'
    end subroutine density_matrix_r
    !
    subroutine density_matrix_q(mo_occ,mos,rho)
      real(xrk), intent(in)     :: mo_occ(:)   ! MO occupation numbers
      complex(xrk), intent(in)  :: mos(:,:,:)  ! Left and right orbitals
      complex(xrk), intent(out) :: rho(:,:)
      !
      include 'scf_tools_density_matrix_common.f90'
    end subroutine density_matrix_q
    !
  end module scf_tools
