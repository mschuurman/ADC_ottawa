!
!  Transformation of 2e integrals to the bi-orthogonal spin-MO basis.
!  The MO format is as described in static_tunnel.f90, namely:
!
!  All alpha-spin AOs are coming first, at positions (1:nao). 
!  The beta-spin AOs follow, at positions (nao+1:nao_spin). 
!
!  The spin-less AO 2e integrals are supplied by integrals_tools.f90
!
!  The resulting 2e integrals generally have only one symmetry operation:
!
!  g(i,k,j,l) = g(j,l,i,k)
!
!  where indices i and j refer to the left MOs; indices k and l are for the
!  right MOs. At least in the initial implementation, we will not be using
!  this symmetry, and explicitly compute all integrals.
!
!  Note that we use the charge-cloud order for storing the integrals;
!  the first two indices correspond to the variable x1; the last two
!  indices correspond to the variable x2. Furhermode, the first and 
!  the third indices correspond to the left eigenvectors; the second
!  and the fourth indices are the right eigenvectors. This all can get
!  a little confusing, so once more:
!
!    g(i,   j,    k,   l)
!      left right left right
!      x1   x1    x2   x2
!
!  The actual integral routines here are a little more general: we do not
!  make any assumptions on the relationships between the MOs used for 
!  three transformation indices. Although -usually- we'll see the MOs
!  coming from a bi-orthogonal set, this is neither expected nor enforced.
!
  module integrals_mo2e
    use accuracy
    use timer
    use math
    use import_gamess
    use gamess_internal
    use integral_tools
    use matrix_tools
    implicit none
    private
    public moint2e_cache
    public transform_moint2e, destroy_moint2e
    public transform_moint2e_real
    public fetch_moint2e
    !
    interface transform_moint2e
      module procedure transform_moint2e_real
!*qd  module procedure transform_moint2e_quad
    end interface transform_moint2e
    !
    !  Internal interfaces; not for public consumption
    !
    interface transform_one_orbital
      module procedure transform_one_orbital_rrr  ! real AO integrals, real MO integrals, real orbitals
!*qd  module procedure transform_one_orbital_rqr  ! real AO integrals, quad MO integrals, real orbitals
!*qd  module procedure transform_one_orbital_rrq  ! real AO integrals, real MO integrals, quad orbitals
!*qd  module procedure transform_one_orbital_rqq  ! real AO integrals, quad MO integrals, quad orbitals
!*qd  module procedure transform_one_orbital_qrr  ! quad AO integrals, real MO integrals, real orbitals
!*qd  module procedure transform_one_orbital_qqr  ! quad AO integrals, quad MO integrals, real orbitals
!*qd  module procedure transform_one_orbital_qrq  ! quad AO integrals, real MO integrals, quad orbitals
!*qd  module procedure transform_one_orbital_qqq  ! quad AO integrals, quad MO integrals, quad orbitals
    end interface transform_one_orbital
    !
    interface transform_ao_integral_block
       module procedure transform_ao_integral_block_rr
!*qd   module procedure transform_ao_integral_block_rq
!*qd   module procedure transform_ao_integral_block_qr
!*qd   module procedure transform_ao_integral_block_qq
    end interface transform_ao_integral_block
    !
    interface transform_index_k
      module procedure transform_index_k_real
!*qd  module procedure transform_index_k_quad
    end interface transform_index_k
    !
    interface transform_index_j
      module procedure transform_index_j_real
!*qd  module procedure transform_index_j_quad
    end interface transform_index_j
    !
    interface transform_index_i
      module procedure transform_index_i_rr   ! Real partial transform, real MO integrals
!*qd  module procedure transform_index_i_rq   ! Real partial transform, quad MO integrals
!*qd  module procedure transform_index_i_qr   ! Quad partial transform, real MO integrals
!*qd  module procedure transform_index_i_qq   ! Quad partial transform, quad MO integrals
    end interface transform_index_i
    !
    interface accumulate_transform
      module procedure accumulate_transform_rr
!*qd  module procedure accumulate_transform_rq
!*qd  module procedure accumulate_transform_qr
!*qd  module procedure accumulate_transform_qq
    end interface accumulate_transform
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: verbose    = 1           ! Level of output
    !
    !  ==== User-adjustable parameters =====
    !
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type moint2e_cache
      integer(ik)           :: nmo(4)                ! Number of MOs for each one of the indices
      character(len=20)     :: mode                  ! Mode of operation; can be one of:
                                                     ! 'incore' - All integrals are kept in memory
                                                     ! 'disk'   - In-memory buffer holds a specific first index
      character(len=20)     :: ints_math             ! Can be either 'real' or 'quad'. Depending on the accuracy
                                                     ! of the AO integrals used in constructing the MO integrals,
                                                     ! actual accuracy may be lower.
      integer(ik)           :: io_unit               ! Direct I/O file containing transformed integrals integrals
      integer(ik)           :: mo_l                  ! Current value for the l index; only meaningful if mode='disk'
      complex(rk), pointer  :: buffer_real(:,:,:,:)  ! Integral buffer for the (i,j,k,l) indices, for ints_math=='real'
                                                     ! The meaning of the last index depends on the mode:
                                                     ! For mode='incore', it is l;
                                                     ! For mode='disk', the last index must be 1.
      complex(xrk), pointer :: buffer_quad(:,:,:,:)  ! Integral buffer for the (i,j,k,l) indices, for ints_math=='quad'
    end type moint2e_cache
    !
    contains
    !
    !  There are thirty-two (8x4) copies of the accum_???? routines below; these are needed to avoid making
    !  integral copies during the transformation and to maintain the possible type mixtures. Oops.
    !
    ! @@ real integrals real mos @@
    subroutine accum_ijkl_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_ijkl_common.f90'
    end subroutine accum_ijkl_rr
    !
    subroutine accum_jikl_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_jikl_common.f90'
    end subroutine accum_jikl_rr
    !
    subroutine accum_ijlk_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_ijlk_common.f90'
    end subroutine accum_ijlk_rr
    !
    subroutine accum_jilk_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_jilk_common.f90'
    end subroutine accum_jilk_rr
    !
    subroutine accum_klij_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_klij_common.f90'
    end subroutine accum_klij_rr
    !
    subroutine accum_klji_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_klji_common.f90'
    end subroutine accum_klji_rr
    !
    subroutine accum_lkij_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_lkij_common.f90'
    end subroutine accum_lkij_rr
    !
    subroutine accum_lkji_rr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rr.f90'
      include 'integrals_mo2e_accum_lkji_common.f90'
    end subroutine accum_lkji_rr
    !
    ! @@ real integrals quad mos @@
    subroutine accum_ijkl_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_ijkl_common.f90'
    end subroutine accum_ijkl_rq
    !
    subroutine accum_jikl_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_jikl_common.f90'
    end subroutine accum_jikl_rq
    !
    subroutine accum_ijlk_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_ijlk_common.f90'
    end subroutine accum_ijlk_rq
    !
    subroutine accum_jilk_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_jilk_common.f90'
    end subroutine accum_jilk_rq
    !
    subroutine accum_klij_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_klij_common.f90'
    end subroutine accum_klij_rq
    !
    subroutine accum_klji_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_klji_common.f90'
    end subroutine accum_klji_rq
    !
    subroutine accum_lkij_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_lkij_common.f90'
    end subroutine accum_lkij_rq
    !
    subroutine accum_lkji_rq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_rq.f90'
      include 'integrals_mo2e_accum_lkji_common.f90'
    end subroutine accum_lkji_rq
    !
    ! @@ quad integrals real mos @@
    subroutine accum_ijkl_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_ijkl_common.f90'
    end subroutine accum_ijkl_qr
    !
    subroutine accum_jikl_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_jikl_common.f90'
    end subroutine accum_jikl_qr
    !
    subroutine accum_ijlk_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_ijlk_common.f90'
    end subroutine accum_ijlk_qr
    !
    subroutine accum_jilk_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_jilk_common.f90'
    end subroutine accum_jilk_qr
    !
    subroutine accum_klij_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_klij_common.f90'
    end subroutine accum_klij_qr
    !
    subroutine accum_klji_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_klji_common.f90'
    end subroutine accum_klji_qr
    !
    subroutine accum_lkij_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_lkij_common.f90'
    end subroutine accum_lkij_qr
    !
    subroutine accum_lkji_qr(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qr.f90'
      include 'integrals_mo2e_accum_lkji_common.f90'
    end subroutine accum_lkji_qr
    !
    ! @@ quad integrals quad mos @@
    subroutine accum_ijkl_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_ijkl_common.f90'
    end subroutine accum_ijkl_qq
    !
    subroutine accum_jikl_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_jikl_common.f90'
    end subroutine accum_jikl_qq
    !
    subroutine accum_ijlk_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_ijlk_common.f90'
    end subroutine accum_ijlk_qq
    !
    subroutine accum_jilk_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_jilk_common.f90'
    end subroutine accum_jilk_qq
    !
    subroutine accum_klij_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_klij_common.f90'
    end subroutine accum_klij_qq
    !
    subroutine accum_klji_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_klji_common.f90'
    end subroutine accum_klji_qq
    !
    subroutine accum_lkij_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_lkij_common.f90'
    end subroutine accum_lkij_qq
    !
    subroutine accum_lkji_qq(p0,sz,a2e,mol,buf_l)
      include 'integrals_mo2e_accum_head_qq.f90'
      include 'integrals_mo2e_accum_lkji_common.f90'
    end subroutine accum_lkji_qq
    !
    !  Unfortunately, include-based cloning accumulate_transform_* functions is more trouble than it's worth
    !
    subroutine accumulate_transform_rr(func,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: p0(:)          ! First AO in the block; SPINLESS
      integer(ik), intent(in)    :: sz(:)          ! Number of AOs
      real(rk), intent(in)       :: a2e(:,:,:,:)   ! 2-e integrals over AOs
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(rk), intent(in)    :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(rk), intent(inout) :: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      interface 
        subroutine func(p0,sz,a2e,mol,buf_l)
          use accuracy
          integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
          integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
          real(rk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
          complex(rk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
          complex(rk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
        end subroutine func
      end interface
      !
      integer(ik) :: spin_l   ! Spin labels associated with integral indices; either 0 (alpha) or 1 (beta)
      !
      scan_spin_l: do spin_l=1,2
        call func(p0,sz,a2e,mol(1+(spin_l-1)*nao:),buf_l(:,:,:,spin_l))
      end do scan_spin_l
    end subroutine accumulate_transform_rr
    !
    subroutine accumulate_transform_rq(func,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: p0(:)          ! First AO in the block; SPINLESS
      integer(ik), intent(in)    :: sz(:)          ! Number of AOs
      real(rk), intent(in)       :: a2e(:,:,:,:)   ! 2-e integrals over AOs
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(xrk), intent(in)   :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(xrk), intent(inout):: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      interface 
        subroutine func(p0,sz,a2e,mol,buf_l)
          use accuracy
          integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
          integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
          real(rk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
          complex(xrk), intent(in)   :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
          complex(xrk), intent(inout):: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
        end subroutine func
      end interface
      !
      integer(ik) :: spin_l   ! Spin labels associated with integral indices; either 0 (alpha) or 1 (beta)
      !
      scan_spin_l: do spin_l=1,2
        call func(p0,sz,a2e,mol(1+(spin_l-1)*nao:),buf_l(:,:,:,spin_l))
      end do scan_spin_l
    end subroutine accumulate_transform_rq
    !
    subroutine accumulate_transform_qr(func,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: p0(:)          ! First AO in the block; SPINLESS
      integer(ik), intent(in)    :: sz(:)          ! Number of AOs
      real(xrk), intent(in)      :: a2e(:,:,:,:)   ! 2-e integrals over AOs
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(rk), intent(in)    :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(rk), intent(inout) :: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      interface 
        subroutine func(p0,sz,a2e,mol,buf_l)
          use accuracy
          integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
          integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
          real(xrk), intent(in)      :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
          complex(rk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
          complex(rk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
        end subroutine func
      end interface
      !
      integer(ik) :: spin_l   ! Spin labels associated with integral indices; either 0 (alpha) or 1 (beta)
      !
      scan_spin_l: do spin_l=1,2
        call func(p0,sz,a2e,mol(1+(spin_l-1)*nao:),buf_l(:,:,:,spin_l))
      end do scan_spin_l
    end subroutine accumulate_transform_qr
    !
    subroutine accumulate_transform_qq(func,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: p0(:)          ! First AO in the block; SPINLESS
      integer(ik), intent(in)    :: sz(:)          ! Number of AOs
      real(xrk), intent(in)      :: a2e(:,:,:,:)   ! 2-e integrals over AOs
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(xrk), intent(in)   :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(xrk), intent(inout):: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      interface 
        subroutine func(p0,sz,a2e,mol,buf_l)
          use accuracy
          integer(ik), intent(in)    :: p0(:)         ! First orbital in the block, including the spin
          integer(ik), intent(in)    :: sz(:)         ! Number of orbitals
          real(xrk), intent(in)      :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
          complex(xrk), intent(in)   :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
          complex(xrk), intent(inout):: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
        end subroutine func
      end interface
      !
      integer(ik) :: spin_l   ! Spin labels associated with integral indices; either 0 (alpha) or 1 (beta)
      !
      scan_spin_l: do spin_l=1,2
        call func(p0,sz,a2e,mol(1+(spin_l-1)*nao:),buf_l(:,:,:,spin_l))
      end do scan_spin_l
    end subroutine accumulate_transform_qq
    !
    !  A very naive routine for last-index transformation of the 2e integrals
    !  Unfortunately, this routine can't be made type-agnostic and cloned because
    !  of the function argument to accumulate_transform()
    !
    subroutine transform_ao_integral_block_rr(bi,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      real(rk), intent(in)       :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(rk), intent(in)    :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(rk), intent(inout) :: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      integer(ik) :: p2(4)                         ! Same as p0, after index swapping 
      integer(ik) :: s2(4)                         ! Same as sz, after index swapping
      !
                                       call accumulate_transform(accum_ijkl_rr,p0,sz,a2e,nao,mol,buf_l)
      if (swap_jikl_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jikl_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_ijlk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_ijlk_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_jilk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jilk_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkij_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klij_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klji_rr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkji_rr,p2,s2,a2e,nao,mol,buf_l) 
    end subroutine transform_ao_integral_block_rr
    !
    subroutine transform_ao_integral_block_rq(bi,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      real(rk), intent(in)       :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(xrk), intent(in)   :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(xrk), intent(inout):: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      integer(ik) :: p2(4)                         ! Same as p0, after index swapping 
      integer(ik) :: s2(4)                         ! Same as sz, after index swapping
      !
                                       call accumulate_transform(accum_ijkl_rq,p0,sz,a2e,nao,mol,buf_l)
      if (swap_jikl_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jikl_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_ijlk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_ijlk_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_jilk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jilk_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkij_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klij_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klji_rq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkji_rq,p2,s2,a2e,nao,mol,buf_l) 
    end subroutine transform_ao_integral_block_rq
    !
    subroutine transform_ao_integral_block_qr(bi,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      real(xrk), intent(in)      :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(rk), intent(in)    :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(rk), intent(inout):: buf_l(:,:,:,:)  ! Partially-transformed integrals
      !
      integer(ik) :: p2(4)                         ! Same as p0, after index swapping 
      integer(ik) :: s2(4)                         ! Same as sz, after index swapping
      !
                                       call accumulate_transform(accum_ijkl_qr,p0,sz,a2e,nao,mol,buf_l)
      if (swap_jikl_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jikl_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_ijlk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_ijlk_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_jilk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jilk_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkij_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klij_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klji_qr,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkji_qr,p2,s2,a2e,nao,mol,buf_l) 
    end subroutine transform_ao_integral_block_qr
    !
    subroutine transform_ao_integral_block_qq(bi,p0,sz,a2e,nao,mol,buf_l)
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      real(xrk), intent(in)      :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
                                                   ! (i<j) < (k<l)
      integer(ik), intent(in)    :: nao            ! Number of spin-less atomic orbitals
      complex(xrk), intent(in)   :: mol(:)         ! AO coefficients for the fourth-index orbital
      complex(xrk), intent(inout):: buf_l(:,:,:,:) ! Partially-transformed integrals
      !
      integer(ik) :: p2(4)                         ! Same as p0, after index swapping 
      integer(ik) :: s2(4)                         ! Same as sz, after index swapping
      !
                                       call accumulate_transform(accum_ijkl_qq,p0,sz,a2e,nao,mol,buf_l)
      if (swap_jikl_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jikl_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_ijlk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_ijlk_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_jilk_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_jilk_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkij_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klij_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klij_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_klji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_klji_qq,p2,s2,a2e,nao,mol,buf_l) 
      if (swap_lkji_s(bi,p0,sz,p2,s2)) call accumulate_transform(accum_lkji_qq,p2,s2,a2e,nao,mol,buf_l) 
    end subroutine transform_ao_integral_block_qq
    !
    subroutine transform_index_k_real(nao,nmo,mo_k,buf_l,buf_kl)
      complex(rk), intent(in)  :: mo_k(:,:)      ! Third-index MOs coefficients
      complex(rk), intent(in)  :: buf_l(:,:,:,:) ! Integrals transformed over the last index; L index is fixed
                                                 ! The fourth index is spin of the L component.
      complex(rk), intent(out) :: buf_kl(:,:,:)  ! Integrals transformed over the last two indices; L index is fixed
      !
      include 'integrals_mo2e_transform_index_k_common.f90'
    end subroutine transform_index_k_real
    !
    subroutine transform_index_k_quad(nao,nmo,mo_k,buf_l,buf_kl)
      complex(xrk), intent(in)  :: mo_k(:,:)      ! Third-index MOs coefficients
      complex(xrk), intent(in)  :: buf_l(:,:,:,:) ! Integrals transformed over the last index; L index is fixed
                                                  ! The fourth index is spin of the L component.
      complex(xrk), intent(out) :: buf_kl(:,:,:)  ! Integrals transformed over the last two indices; L index is fixed
      !
      include 'integrals_mo2e_transform_index_k_common.f90'
    end subroutine transform_index_k_quad
    !
    subroutine transform_index_j_real(nao,nmo,mo_j,buf_kl,buf_jkl)
      complex(rk), intent(in)  :: mo_j(:,:)        ! Second-index MO coefficients to transform over
      complex(rk), intent(in)  :: buf_kl(:,:,:)    ! Integrals transformed over the last two indices; L index is fixed
      complex(rk), intent(out) :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                   ! The fourth index is spin of the J component.
      !
      include 'integrals_mo2e_transform_index_j_common.f90'
    end subroutine transform_index_j_real
    !
    subroutine transform_index_j_quad(nao,nmo,mo_j,buf_kl,buf_jkl)
      complex(xrk), intent(in)  :: mo_j(:,:)        ! Second-index MO coefficients to transform over
      complex(xrk), intent(in)  :: buf_kl(:,:,:)    ! Integrals transformed over the last two indices; L index is fixed
      complex(xrk), intent(out) :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                   ! The fourth index is spin of the J component.
      !
      include 'integrals_mo2e_transform_index_j_common.f90'
    end subroutine transform_index_j_quad
    !
    subroutine transform_index_i_rr(nao,nmo,mo_i,buf_jkl,buffer)
      complex(rk), intent(in)   :: mo_i(:,:)        ! First-index MO coefficient to transform over
      complex(rk), intent(in)   :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                    ! The fourth index is spin of the J component.
      complex(rk), intent(out)  :: buffer(:,:,:)    ! Fully transformed MOs; L index is fixed
      !
      include 'integrals_mo2e_transform_index_i_common.f90'
    end subroutine transform_index_i_rr
    !
    subroutine transform_index_i_rq(nao,nmo,mo_i,buf_jkl,buffer)
      complex(rk), intent(in)   :: mo_i(:,:)        ! First-index MO coefficient to transform over
      complex(rk), intent(in)   :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                    ! The fourth index is spin of the J component.
      complex(xrk), intent(out) :: buffer(:,:,:)    ! Fully transformed MOs; L index is fixed
      !
      include 'integrals_mo2e_transform_index_i_common.f90'
    end subroutine transform_index_i_rq
    !
    subroutine transform_index_i_qr(nao,nmo,mo_i,buf_jkl,buffer)
      complex(xrk), intent(in)  :: mo_i(:,:)        ! First-index MO coefficient to transform over
      complex(xrk), intent(in)  :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                    ! The fourth index is spin of the J component.
      complex(rk), intent(out)  :: buffer(:,:,:)    ! Fully transformed MOs; L index is fixed
      !
      include 'integrals_mo2e_transform_index_i_common.f90'
    end subroutine transform_index_i_qr
    !
    subroutine transform_index_i_qq(nao,nmo,mo_i,buf_jkl,buffer)
      complex(xrk), intent(in)  :: mo_i(:,:)        ! First-index MO coefficient to transform over
      complex(xrk), intent(in)  :: buf_jkl(:,:,:,:) ! Integrals transformed over the all but first index; L index is fixed
                                                    ! The fourth index is spin of the J component.
      complex(xrk), intent(out) :: buffer(:,:,:)    ! Fully transformed MOs; L index is fixed
      !
      include 'integrals_mo2e_transform_index_i_common.f90'
    end subroutine transform_index_i_qq
    !
    subroutine transform_one_orbital_rrr(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(rk), intent(in)          :: mo_i(:,:)          ! First-index MOs
      complex(rk), intent(in)          :: mo_j(:,:)          ! Second-index MOs
      complex(rk), intent(in)          :: mo_k(:,:)          ! Third-index MOs
      complex(rk), intent(in)          :: mo_l(:,:)          ! Current fourth-index MO block
      complex(rk), intent(out)         :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(rk), intent(out)         :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(rk), intent(out)         :: buf_kl (:,:,:)     ! ditto
      complex(rk), intent(out)         :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(rk), intent(in)             :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_rrr
    !
    subroutine transform_one_orbital_rqr(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(rk), intent(in)          :: mo_i(:,:)          ! First-index MOs
      complex(rk), intent(in)          :: mo_j(:,:)          ! Second-index MOs
      complex(rk), intent(in)          :: mo_k(:,:)          ! Third-index MOs
      complex(rk), intent(in)          :: mo_l(:,:)          ! Current fourth-index MO block
      complex(xrk), intent(out)        :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(rk), intent(out)         :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(rk), intent(out)         :: buf_kl (:,:,:)     ! ditto
      complex(rk), intent(out)         :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(rk), intent(in)             :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_rqr
    !
    subroutine transform_one_orbital_qrr(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(rk), intent(in)          :: mo_i(:,:)          ! First-index MOs
      complex(rk), intent(in)          :: mo_j(:,:)          ! Second-index MOs
      complex(rk), intent(in)          :: mo_k(:,:)          ! Third-index MOs
      complex(rk), intent(in)          :: mo_l(:,:)          ! Current fourth-index MO block
      complex(rk), intent(out)         :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L); 
      complex(rk), intent(out)         :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(rk), intent(out)         :: buf_kl (:,:,:)     ! ditto
      complex(rk), intent(out)         :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(xrk), intent(in)            :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_qrr
    !
    subroutine transform_one_orbital_qqr(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(rk), intent(in)          :: mo_i(:,:)          ! First-index MOs
      complex(rk), intent(in)          :: mo_j(:,:)          ! Second-index MOs
      complex(rk), intent(in)          :: mo_k(:,:)          ! Third-index MOs
      complex(rk), intent(in)          :: mo_l(:,:)          ! Current fourth-index MO block
      complex(xrk), intent(out)        :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L); 
      complex(rk), intent(out)         :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(rk), intent(out)         :: buf_kl (:,:,:)     ! ditto
      complex(rk), intent(out)         :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(xrk), intent(in)            :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_qqr
    !
    subroutine transform_one_orbital_rrq(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(xrk), intent(in)         :: mo_i(:,:)          ! First-index MOs
      complex(xrk), intent(in)         :: mo_j(:,:)          ! Second-index MOs
      complex(xrk), intent(in)         :: mo_k(:,:)          ! Third-index MOs
      complex(xrk), intent(in)         :: mo_l(:,:)          ! Current fourth-index MO block
      complex(rk), intent(out)         :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(xrk), intent(out)        :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(xrk), intent(out)        :: buf_kl (:,:,:)     ! ditto
      complex(xrk), intent(out)        :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(rk), intent(in)             :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_rrq
    !
    subroutine transform_one_orbital_rqq(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(xrk), intent(in)         :: mo_i(:,:)          ! First-index MOs
      complex(xrk), intent(in)         :: mo_j(:,:)          ! Second-index MOs
      complex(xrk), intent(in)         :: mo_k(:,:)          ! Third-index MOs
      complex(xrk), intent(in)         :: mo_l(:,:)          ! Current fourth-index MO block
      complex(xrk), intent(out)        :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(xrk), intent(out)        :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(xrk), intent(out)        :: buf_kl (:,:,:)     ! ditto
      complex(xrk), intent(out)        :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(rk), intent(in)             :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_rqq
    !
    subroutine transform_one_orbital_qrq(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(xrk), intent(in)         :: mo_i(:,:)          ! First-index MOs
      complex(xrk), intent(in)         :: mo_j(:,:)          ! Second-index MOs
      complex(xrk), intent(in)         :: mo_k(:,:)          ! Third-index MOs
      complex(xrk), intent(in)         :: mo_l(:,:)          ! Current fourth-index MO block
      complex(rk), intent(out)         :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(xrk), intent(out)        :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(xrk), intent(out)        :: buf_kl (:,:,:)     ! ditto
      complex(xrk), intent(out)        :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(xrk), intent(in)            :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_qrq
    !
    subroutine transform_one_orbital_qqq(int2e,nao,nmo,mo_i,mo_j,mo_k,mo_l,buffer,buf_jkl,buf_kl,buf_l,ints_kind)
      type(int2e_cache), intent(inout) :: int2e              ! 2E integrals over spin-less AOs
      integer(ik), intent(in)          :: nao                ! Number of spin-less atomic orbitals
      integer(ik), intent(in)          :: nmo(:)             ! Number of MOs for each index
      complex(xrk), intent(in)         :: mo_i(:,:)          ! First-index MOs
      complex(xrk), intent(in)         :: mo_j(:,:)          ! Second-index MOs
      complex(xrk), intent(in)         :: mo_k(:,:)          ! Third-index MOs
      complex(xrk), intent(in)         :: mo_l(:,:)          ! Current fourth-index MO block
      complex(xrk), intent(out)        :: buffer (:,:,:,:)   ! Integrals over (i,j,k,L);
      complex(xrk), intent(out)        :: buf_jkl(:,:,:,:)   ! Scratch for partially transformed integrals
      complex(xrk), intent(out)        :: buf_kl (:,:,:)     ! ditto
      complex(xrk), intent(out)        :: buf_l  (:,:,:,:,:) ! ditto; fourth index is spin - samee as in buf_jkl
                                                             ! last index is L
      real(xrk), intent(in)            :: ints_kind          ! Serves only to fix the kind of the AO integrals
      !
      include 'integrals_mo2e_transform_one_orbital_common.f90'
    end subroutine transform_one_orbital_qqq
    !
    subroutine transform_moint2e_real(int2e,mode,mo_i,mo_j,mo_k,mo_l,moint2e,io_unit,l_block,storage_mode)
      type(int2e_cache), intent(inout)       :: int2e        ! 2E integrals over spin-less AOs
      character(len=*), intent(in)           :: mode         ! Mode of operation; either 'incore' or 'disk'
      complex(rk), intent(in)                :: mo_i(:,:)    ! Spin-MOs to be used for the first integral index
                                                             ! Index 1: all spin-alpha AOs, followed by all spin-beta AOs
                                                             ! Index 2: MO index
      complex(rk), intent(in)                :: mo_j(:,:)    ! Spin-MOs to be used for the second integral index
      complex(rk), intent(in)                :: mo_k(:,:)    ! Spin-MOs to be used for the third integral index
      complex(rk), intent(in)                :: mo_l(:,:)    ! Spin-MOs to be used for the fourth integral index
                                                             ! Note that all mo_* arrays must have the same first dimension;
                                                             ! the second dimension however does not need to be identical
      type(moint2e_cache), intent(inout)     :: moint2e      ! MO integral descriptor
      integer(ik), intent(in), optional      :: io_unit      ! I/O unit for storing MO integrals; not needed for mode='incore'
      integer(ik), intent(in), optional      :: l_block      ! Number of the first-index MOs to transform at the same time 
                                                             ! (this is an I/O optionization)
      character(len=*), intent(in), optional :: storage_mode ! Override precision of the final transformed integrals;
                                                             ! 'real', 'quad', or 'as-is'. The default is the same kind
                                                             ! as used for the MOs.
      !
      include 'integrals_mo2e_transform_moint2e_common.f90'
    end subroutine transform_moint2e_real
    !
    subroutine transform_moint2e_quad(int2e,mode,mo_i,mo_j,mo_k,mo_l,moint2e,io_unit,l_block,storage_mode)
      type(int2e_cache), intent(inout)       :: int2e        ! 2E integrals over spin-less AOs
      character(len=*), intent(in)           :: mode         ! Mode of operation; either 'incore' or 'disk'
      complex(xrk), intent(in)               :: mo_i(:,:)    ! Spin-MOs to be used for the first integral index
                                                             ! Index 1: all spin-alpha AOs, followed by all spin-beta AOs
                                                             ! Index 2: MO index
      complex(xrk), intent(in)               :: mo_j(:,:)    ! Spin-MOs to be used for the second integral index
      complex(xrk), intent(in)               :: mo_k(:,:)    ! Spin-MOs to be used for the third integral index
      complex(xrk), intent(in)               :: mo_l(:,:)    ! Spin-MOs to be used for the fourth integral index
                                                             ! Note that all mo_* arrays must have the same first dimension;
                                                             ! the second dimension however does not need to be identical
      type(moint2e_cache), intent(inout)     :: moint2e      ! MO integral descriptor
      integer(ik), intent(in), optional      :: io_unit      ! I/O unit for storing MO integrals; not needed for mode='incore'
      integer(ik), intent(in), optional      :: l_block      ! Number of the first-index MOs to transform at the same time 
                                                             ! (this is an I/O optionization)
      character(len=*), intent(in), optional :: storage_mode ! Override precision of the final transformed integrals;
                                                             ! 'real', 'quad', or 'as-is'. The default is the same kind
                                                             ! as used for the MOs.
      !
      include 'integrals_mo2e_transform_moint2e_common.f90'
    end subroutine transform_moint2e_quad
    !
    subroutine destroy_moint2e(moint2e)
      type(moint2e_cache), intent(inout) :: moint2e
      logical :: opened
      !
      call TimerStart('Destroy MO 2E cache')
      !
      if (associated(moint2e%buffer_real)) deallocate (moint2e%buffer_real)
      if (associated(moint2e%buffer_quad)) deallocate (moint2e%buffer_quad)
      nullify (moint2e%buffer_real)
      nullify (moint2e%buffer_quad)
      if (moint2e%mode=='disk') then
        inquire (moint2e%io_unit,opened=opened)
        if (opened) close (moint2e%io_unit,status='delete')
      end if
      !
      call TimerStop('Destroy MO 2E cache')
    end subroutine destroy_moint2e
    !
    subroutine fetch_moint2e(moint2e,l_index)
      type(moint2e_cache), intent(inout) :: moint2e   ! MO integral descriptor
      integer(ik), intent(in)            :: l_index   ! Desired fourth MO index
      !
      integer(ik) :: ios
      !
      if (moint2e%mode/='disk') return
      if (moint2e%mo_l==l_index) return
      call TimerStart('Fetch MO 2e integrals')
      !
      if (l_index<=0 .or. l_index>moint2e%nmo(4)) then
        write (out,"('fetch_moint2e: Orbital index ',i0,' is out of bounds.')") l_index
        stop 'integrals_mo2e%fetch_moint2e - bad MO index'
      end if
      !
      select case (moint2e%ints_math)
        case default; stop 'integrals_mo2e%fetch_moint2e - bad ints_math value'
        case ('real') ; read (moint2e%io_unit,rec=l_index,iostat=ios) moint2e%buffer_real
        case ('quad') ; read (moint2e%io_unit,rec=l_index,iostat=ios) moint2e%buffer_quad
      end select
      if (ios/=0) then
        write (out,"('Error ',i0,' reading transformed integrals for L= ',i0)") ios, l_index
        stop 'integrals_mo2e%fetch_moint2e - I/O error reading transformed integrals'
      end if
      moint2e%mo_l = l_index
      !
      call TimerStop('Fetch MO 2e integrals')
      ! call TimerReport
    end subroutine fetch_moint2e
  end module integrals_mo2e
