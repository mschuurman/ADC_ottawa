!
!  Evaluation of the 2-electron part of the Fock matrix, using 4-centre 2-electron integrals
!
!  Depending on the available real kinds, this module has up to three distinct entry points:
!
!  a. Standard-kind (rk) for both the density matrix and the integrals
!  b. Standard-kind (rk) for the integrals, accurate kind (xrk) for the density
!  c. Accurate-kind (xrk) for both density and the integrals
!
!  The routines in this moduel will run in parallel. Because of the dependence on
!  the integrals descriptor in int2e_cache (see integral_tools.f90) it may be tricky
!  (although possible?) to call this module while already inside a parallel region.
!
  module fock_tools
    use accuracy
    use timer
    use import_gamess
    use integral_tools
    implicit none
    private
    public fock_g_matrix
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: verbose    = 1           ! Level of output
    !
    !  We are reusing a bunch of code using includes to implement several versions of the
    !  Fock matrix construction routines. In order to dispatch the calls correctly, we
    !  have to use interface blocks here; everything else uses the type-generic name.
    !
    interface accumulate_g_coulomb
      module procedure accumulate_g_coulomb_rr
!*qd  module procedure accumulate_g_coulomb_rq
!*qd  module procedure accumulate_g_coulomb_qq
    end interface accumulate_g_coulomb
    !
    interface accumulate_g_exchange
      module procedure accumulate_g_exchange_rr
!*qd  module procedure accumulate_g_exchange_rq
!*qd  module procedure accumulate_g_exchange_qq
    end interface accumulate_g_exchange
    !
    interface accumulate_g
      module procedure accumulate_g_rr
!*qd  module procedure accumulate_g_rq
!*qd  module procedure accumulate_g_qq
    end interface accumulate_g
    !
    interface g_matrix_block
      module procedure g_matrix_block_rr
!*qd  module procedure g_matrix_block_rq
!*qd  module procedure g_matrix_block_qq
!*qd  module procedure g_matrix_block_qr
    end interface g_matrix_block
    !
    interface fock_g_matrix
      module procedure g_matrix_real
!*qd  module procedure g_matrix_quad
    end interface fock_g_matrix
    !
    contains
    !
    subroutine accumulate_g_coulomb_rr(p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_coulomb_common.f90'
    end subroutine accumulate_g_coulomb_rr
    !
    subroutine accumulate_g_coulomb_rq(p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)        :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_coulomb_common.f90'
    end subroutine accumulate_g_coulomb_rq
    !
    subroutine accumulate_g_coulomb_qq(p0,sz,mxsz,a2e,rho,glocal)
      real(xrk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_coulomb_common.f90'
    end subroutine accumulate_g_coulomb_qq
    !
    subroutine accumulate_g_exchange_rr(p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_exchange_common.f90'
    end subroutine accumulate_g_exchange_rr
    !
    subroutine accumulate_g_exchange_rq(p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)        :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_exchange_common.f90'
    end subroutine accumulate_g_exchange_rq
    !
    subroutine accumulate_g_exchange_qq(p0,sz,mxsz,a2e,rho,glocal)
      real(xrk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_exchange_common.f90'
    end subroutine accumulate_g_exchange_qq
    !
    subroutine accumulate_g_rr(nao,p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(rk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(rk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_common.f90'
    end subroutine accumulate_g_rr
    !
    subroutine accumulate_g_rq(nao,p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)        :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_common.f90'
    end subroutine accumulate_g_rq
    !
    subroutine accumulate_g_qq(nao,p0,sz,mxsz,a2e,rho,glocal)
      real(xrk), intent(in)       :: a2e(:,:,:,:)  ! 2-e integrals over AOs
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Per-thread copy of gmat
      !
      include 'fock_tools_accumulate_g_common.f90'
    end subroutine accumulate_g_qq
    !
    subroutine g_matrix_block_rr(nao,bi,p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)       :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      complex(rk), intent(in)    :: rho   (:,:)    ! Shared density matrix
      complex(rk), intent(inout) :: glocal(:,:)    ! Pre-thread copy of gmat; safe to update 
      !
      include 'fock_tools_g_matrix_block_common.f90'
    end subroutine g_matrix_block_rr
    !
    subroutine g_matrix_block_rq(nao,bi,p0,sz,mxsz,a2e,rho,glocal)
      real(rk), intent(in)        :: a2e(:,:,:,:)  ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Pre-thread copy of gmat; safe to update 
      !
      include 'fock_tools_g_matrix_block_common.f90'
    end subroutine g_matrix_block_rq
    !
    subroutine g_matrix_block_qr(nao,bi,p0,sz,mxsz,a2e,rho,glocal)
      integer(ik), intent(in)    :: nao            ! Number of spin-less AOs
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      integer(ik), intent(in)    :: mxsz           ! Largest size of an integral block
      real(xrk), intent(in)        :: a2e(:,:,:,:)  ! Integrals block, for the canonical order of block indices
                                                    ! (i<j) < (k<l)
      complex(rk), intent(in)      :: rho   (:,:)   ! Shared density matrix
      complex(rk), intent(inout)   :: glocal(:,:)   ! Pre-thread copy of gmat; safe to update 
      !
      write (out,"('Combining quad-precision integrals with standard-precision density does not make sense.')")
      call flush(out)
      stop 'fock_tools%g_matrix_block_qr - meaningless case not implemented'
    end subroutine g_matrix_block_qr
    !
    subroutine g_matrix_block_qq(nao,bi,p0,sz,mxsz,a2e,rho,glocal)
      real(xrk), intent(in)       :: a2e(:,:,:,:)  ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      complex(xrk), intent(in)    :: rho   (:,:)   ! Shared density matrix
      complex(xrk), intent(inout) :: glocal(:,:)   ! Pre-thread copy of gmat; safe to update 
      !
      include 'fock_tools_g_matrix_block_common.f90'
    end subroutine g_matrix_block_qq
    !
    subroutine g_matrix_real(int2e,rho,gmat)
      complex(rk), intent(in)               :: rho (:,:) ! Density matrix
      complex(rk), intent(out)              :: gmat(:,:) ! 2-electron contribution to the Fock matrix
      !
      include 'fock_tools_g_matrix_common.f90'
    end subroutine g_matrix_real
    !
    subroutine g_matrix_quad(int2e,rho,gmat)
      complex(xrk), intent(in)          :: rho (:,:) ! Density matrix
      complex(xrk), intent(out)         :: gmat(:,:) ! 2-electron contribution to the Fock matrix
      !
      include 'fock_tools_g_matrix_common.f90'
    end subroutine g_matrix_quad
  end module fock_tools
