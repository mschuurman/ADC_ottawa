!
!  Various tools for dealing with bio-orthogonal orbital sets.
!  The notable entry points are:
!
!    bt_follow_mos   - Try to maximally align set of MOs to a given reference
!    bt_sanity_check - Confirm orbital bi-orthogonality
!    bt_eigenvectors - Biorthogonalize eigenvectrors returned by LAPACK zgeev
!
  module biorthogonal_tools
    use accuracy
    use timer
    use sort_tools
    use lapack
    use matrix_tools
    implicit none
    private
    public bt_follow_mos
    public bt_sanity_check
    public bt_eigenvectors
    !
    interface bt_follow_mos
      module procedure bt_follow_mos_rk
!*qd  module procedure bt_follow_mos_xk
    end interface bt_follow_mos
    !
    interface bt_sanity_check
      module procedure bt_sanity_check_rk
!*qd  module procedure bt_sanity_check_xk
    end interface bt_sanity_check
    !
    interface bt_eigenvectors
      module procedure bt_eigenvectors_rk
!*qd  module procedure bt_eigenvectors_xk
    end interface bt_eigenvectors
    !
    !  Internal-use interfaces
    !
    interface find_optimal_rotation
      module procedure find_optimal_rotation_rk
!*qd  module procedure find_optimal_rotation_xk
    end interface find_optimal_rotation
    !
    interface biorthogonalize_block_of_eigenvectors
      module procedure biorthogonalize_block_of_eigenvectors_rk
!*qd  module procedure biorthogonalize_block_of_eigenvectors_xk
    end interface biorthogonalize_block_of_eigenvectors
    !
    interface eigenvalue_degeneracy_epsilon
      module procedure eigenvalue_degeneracy_epsilon_rk
!*qd  module procedure eigenvalue_degeneracy_epsilon_xk
    end interface eigenvalue_degeneracy_epsilon
    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: verbose         = 0              ! Level of output
    !
    contains
    !
    subroutine bt_sanity_check_rk(evec,msg,sphalf,null_count)
      complex(rk), intent(in)           :: evec(:,:,:)  ! Eigenvectors
      character(len=*), intent(in)      :: msg          ! Name of the caller for error messages
      real(rk), intent(in), optional    :: sphalf(:,:)  ! Metric tensor; if missing assume unity.
      integer(ik), intent(in), optional :: null_count   ! Number of null orbitals expected
      !
      include 'biorthogonal_tools_bt_sanity_check_common.f90'
    end subroutine bt_sanity_check_rk
    !
    subroutine bt_sanity_check_xk(evec,msg,sphalf,null_count)
      complex(xrk), intent(in)          :: evec(:,:,:)  ! Eigenvectors
      character(len=*), intent(in)      :: msg          ! Name of the caller for error messages
      real(xrk), intent(in), optional   :: sphalf(:,:)  ! Metric tensor; if missing assume unity.
      integer(ik), intent(in), optional :: null_count   ! Number of null orbitals expected
      !
      include 'biorthogonal_tools_bt_sanity_check_common.f90'
    end subroutine bt_sanity_check_xk
    !
    !  Bi-orthonormalize a block of degenerate eigenvectors
    !  We call them "MOs" here, but they do not have to be ...
    !
    subroutine biorthogonalize_block_of_eigenvectors_rk(mo_block)
      complex(rk), intent(inout) :: mo_block(:,:,:)
      !
      include 'biorthogonal_tools_block_of_eigenvectors_common.f90'
    end subroutine biorthogonalize_block_of_eigenvectors_rk
    !
    subroutine biorthogonalize_block_of_eigenvectors_xk(mo_block)
      complex(xrk), intent(inout) :: mo_block(:,:,:)
      !
      include 'biorthogonal_tools_block_of_eigenvectors_common.f90'
    end subroutine biorthogonalize_block_of_eigenvectors_xk
    !
    subroutine bt_eigenvectors_rk(nev,eval,evec,epsx)
      integer(ik), intent(in)     :: nev        ! Number of eigenvalues
      complex(rk), intent(inout) :: eval(:)     ! In:  Eigenvalues, in no particular order
                                                ! Out: Eigenvalues, in order of increasing real part
      complex(rk), intent(inout) :: evec(:,:,:) ! In/out: Left [(:,:,1)] and right [(:,:,2)] eigenvectors
                                                ! On output, the eigenvectors are biorthogonal.
      real(rk), intent(in)       :: epsx        ! The threshold for detecting degenerate eigenvalues
                                                ! Zero or negative uses tighest sensible value given the accuracy
      !
      include 'biorthogonal_tools_bt_eigenvectors_common.f90'
    end subroutine bt_eigenvectors_rk
    !
    subroutine bt_eigenvectors_xk(nev,eval,evec,epsx)
      integer(ik), intent(in)     :: nev         ! Number of eigenvalues
      complex(xrk), intent(inout) :: eval(:)     ! In:  Eigenvalues, in no particular order
                                                 ! Out: Eigenvalues, in order of increasing real part
      complex(xrk), intent(inout) :: evec(:,:,:) ! In/out: Left [(:,:,1)] and right [(:,:,2)] eigenvectors
                                                 ! On output, the eigenvectors are biorthogonal.
      real(xrk), intent(in)       :: epsx        ! The threshold for detecting degenerate eigenvalues
                                                 ! Zero or negative uses tighest sensible value given the accuracy
      !
      include 'biorthogonal_tools_bt_eigenvectors_common.f90'
    end subroutine bt_eigenvectors_xk
    !
    !  Find the optimal rotation within a degenerate sub-block, aligning
    !  degenerate MOs to the reference. 
    !
    !  We use Newton-Rafson iterations to locate the stationary point in the sum
    !  of the diagonal occupation numbers, using exponential represenation of the
    !  rotation matrix.
    !
    subroutine find_optimal_rotation_rk(tmat,umat,vcur)
      complex(rk), intent(in)  :: tmat(:,:)  ! Subset overlap of the current left MOs with reference right MOs.
      complex(rk), intent(in)  :: umat(:,:)  ! Subset overlap of the reference left MOs with current right MOs.
      complex(rk), intent(out) :: vcur(:,:)  ! Unitary transformation of the right MOs.
                                             ! Left MOs are to be transformed by conjg(vcur)
      !
      include 'biorthogonal_tools_find_optimal_rotation_common.f90'
    end subroutine find_optimal_rotation_rk
    !
    subroutine find_optimal_rotation_xk(tmat,umat,vcur)
      complex(xrk), intent(in)  :: tmat(:,:)  ! Subset overlap of the current left MOs with reference right MOs.
      complex(xrk), intent(in)  :: umat(:,:)  ! Subset overlap of the reference left MOs with current right MOs.
      complex(xrk), intent(out) :: vcur(:,:)  ! Unitary transformation of the right MOs.
                                              ! Left MOs are to be transformed by conjg(vcur)
      !
      include 'biorthogonal_tools_find_optimal_rotation_common.f90'
    end subroutine find_optimal_rotation_xk
    !
    !  Establish correspondence between current and reference MOs, 
    !  rearranging and rotating degenate orbitals as needed.
    !  The procedure is in principle simple, but having left and right
    !  MOs causes some complications. 
    !  See: biorthogonal-orbital-similarity-good.pdf for details.
    !
    !  At the entry to this routine, MOs should be sorted in the order of
    !  increasing real part of the eigenvalue.
    !
    !  Upon exit, they will be sorted in the same order as the corresponding
    !  reference MOs, with the eigenvalues reordered to match.
    !
    subroutine bt_follow_mos_rk(nmo_act,mo_occ,eps_eval,sphalf,mosg,mo_energy,mos)
      integer(ik), intent(in)    :: nmo_act       ! Number of orbitals to include in the matching process
      real(rk), intent(in)       :: eps_eval      ! Degeneracy threshold for eigenvalues. Zero or negative means machine accuracy
      real(rk), intent(in)       :: sphalf(:,:)   ! Overlap matrix to the power +1/2
      complex(rk), intent(in)    :: mosg(:,:,:)   ! Reference MOs; the last index is 1 of the left vectors; 2 for the right vectors
      real(rk), intent(in)       :: mo_occ(:)     ! Orbital occupations, only used for printing
      complex(rk), intent(inout) :: mo_energy(:)  ! In/Out: Orbital energies, same order as mos()
      complex(rk), intent(inout) :: mos(:,:,:)    ! In: MOs to be aligned
                                                  ! Out: MOs where degenerate subblocks are in maximal conincidence with the guess
      !
      include 'biorthogonal_tools_bt_follow_mos_common.f90'
    end subroutine bt_follow_mos_rk
    !
    subroutine bt_follow_mos_xk(nmo_act,mo_occ,eps_eval,sphalf,mosg,mo_energy,mos)
      integer(ik), intent(in)     :: nmo_act       ! Number of orbitals to include in the matching process
      real(xrk), intent(in)       :: eps_eval      ! Degeneracy threshold for eigenvalues. Zero or negative means machine accuracy
      real(xrk), intent(in)       :: sphalf(:,:)   ! Overlap matrix to the power +1/2
      complex(xrk), intent(in)    :: mosg(:,:,:)   ! Reference MOs; the last index is 1 of the left vectors; 2 for the right vectors
      real(xrk), intent(in)       :: mo_occ(:)     ! Orbital occupations, only used for printing
      complex(xrk), intent(inout) :: mo_energy(:)  ! In/Out: Orbital energies, same order as mos()
      complex(xrk), intent(inout) :: mos(:,:,:)    ! In: MOs to be aligned
                                                   ! Out: MOs where degenerate subblocks are in maximal conincidence with the guess
      !
      include 'biorthogonal_tools_bt_follow_mos_common.f90'
    end subroutine bt_follow_mos_xk
    !
    function eigenvalue_degeneracy_epsilon_rk(epsx,eval,msg) result(eps)
      real(rk), intent(in)         :: epsx    ! Desired threshold; negative or zero means default
      complex(rk), intent(in)      :: eval(:) ! Eigenvalues
      character(len=*), intent(in) :: msg     ! Name of the calling routine for warnings
      real(rk)                     :: eps     ! Final threshold
      !
      eps = 1e5*spacing(maxval(abs(eval))) ! ????
      if (epsx>0) then
        if (epsx < eps) then
          write (out,"(/'WARNING: degeneracy threshold (',g12.5,') is tigher than accuracy warrants (',g12.5,')')") &
                 epsx, eps
          write (out,"( 'WARNING: Expect numerical trouble in ',a/)") trim(msg)
        end if
        eps = epsx
      end if
    end function eigenvalue_degeneracy_epsilon_rk
    !
    function eigenvalue_degeneracy_epsilon_xk(epsx,eval,msg) result(eps)
      real(xrk), intent(in)        :: epsx    ! Desired threshold; negative or zero means default
      complex(xrk), intent(in)     :: eval(:) ! Eigenvalues
      character(len=*), intent(in) :: msg     ! Name of the calling routine for warnings
      real(xrk)                    :: eps     ! Final threshold
      !
      eps = 1e5*spacing(maxval(abs(eval))) ! ????
      if (epsx>0) then
        if (epsx < eps) then
          write (out,"(/'WARNING: degeneracy threshold (',g12.5,') is tigher than accuracy warrants (',g12.5,')')") &
                 epsx, eps
          write (out,"( 'WARNING: Expect numerical trouble in ',a/)") trim(msg)
        end if
        eps = epsx
      end if
    end function eigenvalue_degeneracy_epsilon_xk
    !
    subroutine stop(message)
      character(len=*), intent(in) :: message
      !
      write (out,"('STOP: ',a)") trim(message)
      write (0,"('STOP: ',a)") trim(message)
      call flush(out)
      call flush(0)
      stop 'error in biorthogonal_tools.f90'
    end subroutine stop
  end module biorthogonal_tools
