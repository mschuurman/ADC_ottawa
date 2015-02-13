!
!  A very naive implementation of Pulay's DIIS for convergence acceleration,
!  minimizing norm of (FPS - SPF). The routines here are a little more
!  complicated than they could have been since we want to support multiple
!  precisions.
!
  module diis
    use accuracy
    use timer
    use math
    use lapack
    use matrix_tools
    implicit none
    private
    public diis_state
    public diis_initialize
    public diis_destroy
    public diis_extrapolate
    !
    interface diis_extrapolate
      module procedure diis_extrapolate_wrapper_rk
!*qd  module procedure diis_extrapolate_wrapper_xk
    end interface diis_extrapolate
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: verbose = 0              ! Level of output
    !
    !  ==== Internal DIIS state ====
    !
    type diis_state
      character(len=20)          :: math = 'real'        ! Can be 'real' or 'quad'
      integer(ik)                :: max_nvec   = 20      ! Maximum number of DIIS vectors allowed
      real(rk)                   :: max_coeff  = 20._rk  ! Restart DIIS if any of the coefficients exceed this threshold
      integer(ik)                :: nvec                 ! Current number of DIIS expansion vectors
      !
      !  These arrays will be active if diis_math=='real'
      !
      complex(rk), allocatable   :: fmat_rk(:,:,:)       ! Table of Fock operators for DIIS extrapolation
      complex(rk), allocatable   :: error_rk(:,:,:)      ! Table of DIIS residuals [F P^T S - S P^T F]
      complex(rk), allocatable   :: eiej_rk(:,:)         ! Overlaps of DIIS error vectors
      !
      !  These arrays will be active if diis_math=='quad'
      !
      complex(xrk), allocatable  :: fmat_xrk(:,:,:)      ! Table of Fock operators for DIIS extrapolation
      complex(xrk), allocatable  :: error_xrk(:,:,:)     ! Table of DIIS residuals [F P^T S - S P^T F]
      complex(xrk), allocatable  :: eiej_xrk(:,:)        ! Overlaps of DIIS error vectors
    end type diis_state
    !
    contains
    !
    subroutine diis_initialize(diis,nao,max_diis_nvec,max_diis_coeff,diis_math)
      type(diis_state), intent(out)          :: diis           ! DIIS state
      integer(ik), intent(in)                :: nao            ! Size of the matrices in DIIS
      integer(ik), intent(in), optional      :: max_diis_nvec  ! Max number of DIIS vectors
      real(rk), intent(in), optional         :: max_diis_coeff ! Largest allowed DIIS extratpolation coefficient
      character(len=*), intent(in), optional :: diis_math      ! Numerical accuracy to use in DIIS
      !
      integer(ik) :: alloc
      real(rk)    :: rnao, ram
      !
      !  Initialize the limits and dimensions
      !
      diis%math      = 'real'
      diis%max_nvec  = 20
      diis%max_coeff = 20._rk
      diis%nvec      = 0
      if (present(diis_math)     ) diis%math      = diis_math
      if (present(max_diis_nvec) ) diis%max_nvec  = max_diis_nvec
      if (present(max_diis_coeff)) diis%max_coeff = max_diis_coeff
      !
      rnao = nao ! To force calculation to floating point!
      ram  = (2*2*rnao**2*diis%max_nvec)
      !
      select case (diis%math)
        case default
          stop 'diis%diis_initialize - bad diis_math'
        case ('real')
          ram = ram * rk_bytes
          allocate (diis%fmat_rk (nao,nao,diis%max_nvec), &
                    diis%error_rk(nao,nao,diis%max_nvec), &
                    diis%eiej_rk (diis%max_nvec,diis%max_nvec),stat=alloc)
        case ('quad')
          ram = ram * xrk_bytes
          allocate (diis%fmat_xrk (nao,nao,diis%max_nvec), &
                    diis%error_xrk(nao,nao,diis%max_nvec), &
                    diis%eiej_xrk (diis%max_nvec,diis%max_nvec),stat=alloc)
      end select
      !
      if (verbose>=0) then
        write (out,"(/'Allocating about ',f0.6,' GBytes in ',a,'-precision DIIS-related arrays.'/)") &
               (1.0_rk/1024.0_rk**3) * ram, trim(diis%math)
      end if
      if (alloc/=0) stop 'diis%diis_initialize - allocation failed'
    end subroutine diis_initialize
    !
    subroutine diis_destroy(diis)
      type(diis_state), intent(inout) :: diis ! DIIS state
      !
      if (allocated(diis%fmat_rk )) deallocate (diis%fmat_rk)
      if (allocated(diis%error_rk)) deallocate (diis%error_rk)
      if (allocated(diis%eiej_rk )) deallocate (diis%eiej_rk)
      !
      if (allocated(diis%fmat_xrk )) deallocate (diis%fmat_xrk)
      if (allocated(diis%error_xrk)) deallocate (diis%error_xrk)
      if (allocated(diis%eiej_xrk )) deallocate (diis%eiej_xrk)
    end subroutine diis_destroy
    !
    !  A very naive implementation of DIIS
    !
    subroutine diis_extrapolate_wrapper_rk(diis,iter,smat,rho,fmat)
      type(diis_state), intent(inout) :: diis      ! DIIS state
      integer(ik), intent(in)         :: iter      ! Current iteration; starts at 1
      real(rk), intent(in)            :: smat(:,:) ! Overlap matrix
      complex(rk), intent(in)         :: rho (:,:) ! Current density matrix
      complex(rk), intent(inout)      :: fmat(:,:) ! Input: Fock matrix corresponding to rho(:,:)
                                                   ! Output: extrapolated Fock matrix
      !
      if (diis%math/='real') stop 'diis%diis_extrapolate_wrapper_rk - inconsistent diis%math'
      call diis_extrapolate_rk(diis,iter,smat,rho,fmat,diis%fmat_rk,diis%error_rk,diis%eiej_rk)
    end subroutine diis_extrapolate_wrapper_rk
    !
    subroutine diis_extrapolate_wrapper_xk(diis,iter,smat,rho,fmat)
      type(diis_state), intent(inout) :: diis      ! DIIS state
      integer(ik), intent(in)         :: iter      ! Current iteration; starts at 1
      real(xrk), intent(in)           :: smat(:,:) ! Overlap matrix
      complex(xrk), intent(in)        :: rho (:,:) ! Current density matrix
      complex(xrk), intent(inout)     :: fmat(:,:) ! Input: Fock matrix corresponding to rho(:,:)
                                                   ! Output: extrapolated Fock matrix
      !
      if (diis%math/='quad') stop 'diis%diis_extrapolate_wrapper_rk - inconsistent diis%math'
      call diis_extrapolate_xk(diis,iter,smat,rho,fmat,diis%fmat_xrk,diis%error_xrk,diis%eiej_xrk)
    end subroutine diis_extrapolate_wrapper_xk
    !
    subroutine diis_extrapolate_rk(diis,iter,smat,rho,fmat,diis_fmat,diis_error,diis_eiej)
      type(diis_state), intent(inout)  :: diis              ! DIIS state
      integer(ik), intent(in)          :: iter              ! Current iteration; starts at 1
      real(rk), intent(in)             :: smat(:,:)         ! Overlap matrix
      complex(rk), intent(in)          :: rho (:,:)         ! Current density matrix
      complex(rk), intent(inout)       :: fmat(:,:)         ! Input: Fock matrix corresponding to rho(:,:)
                                                            ! Output: extrapolated Fock matrix
      complex(rk), intent(inout)       :: diis_fmat (:,:,:) ! The appropriate-precision field in diis%
      complex(rk), intent(inout)       :: diis_error(:,:,:) ! The appropriate-precision field in diis%
      complex(rk), intent(inout)       :: diis_eiej (:,:)   ! The appropriate-precision field in diis%
      !
      include 'diis_extrapolate_common.f90'
    end subroutine diis_extrapolate_rk
    !
    subroutine diis_extrapolate_xk(diis,iter,smat,rho,fmat,diis_fmat,diis_error,diis_eiej)
      type(diis_state), intent(inout)  :: diis              ! DIIS state
      integer(ik), intent(in)          :: iter              ! Current iteration; starts at 1
      real(xrk), intent(in)            :: smat(:,:)         ! Overlap matrix
      complex(xrk), intent(in)         :: rho (:,:)         ! Current density matrix
      complex(xrk), intent(inout)      :: fmat(:,:)         ! Input: Fock matrix corresponding to rho(:,:)
                                                            ! Output: extrapolated Fock matrix
      complex(xrk), intent(inout)      :: diis_fmat (:,:,:) ! The appropriate-precision field in diis%
      complex(xrk), intent(inout)      :: diis_error(:,:,:) ! The appropriate-precision field in diis%
      complex(xrk), intent(inout)      :: diis_eiej (:,:)   ! The appropriate-precision field in diis%
      !
      include 'diis_extrapolate_common.f90'
    end subroutine diis_extrapolate_xk
  end module diis
