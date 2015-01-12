!
!  We also have a very naive implementation of single-determinant MP2 energies
!  for a FON-SCF orbital set. This implementation requires the complete set of
!  transformed 2-electron integrals to be stored in memory. 
!  Most of the code is copied/adapted from ci_tools.f90
!
!  We also have a more reasonable in-core version, which only stores two-virtual
!  index integrals; unfortunately, since integral symmetry is much reduced in the 
!  non-Hermitian case compared to the familiar real, symmetric case, integral 
!  transformations get a little awkward.
!
!  Core here is definitely not designed to be efficient.
!
module mp2_tools
  use accuracy
  use timer
  use import_gamess
  use integral_tools
  use integrals_mo2e
  implicit none
  private
  public mp2_energy_incore
  !
  interface mp2_energy_incore
     module procedure mp2_energy_incore_real
!*qd module procedure mp2_energy_incore_quad
  end interface mp2_energy_incore
  !
  interface incore_mp2_single_shift
     module procedure incore_mp2_single_shift_real
!*qd module procedure incore_mp2_single_shift_quad
  end interface incore_mp2_single_shift
  !
  !  Fixed parameters
  !
  integer(ik), parameter   :: verbose       = 1
  !
  !  Real kind employed in MP2 calculation can be adjusted at compile-time.
  !  It does not have to match the real kind used in SCF.
  !
  type isa_statistics 
    integer(hik) :: terms       = 0 ! Total number of terms
    integer(hik) :: big_terms   = 0 ! Total number of "large" contributions to the imaginary part of the energy
    integer(hik) :: hard_zeros  = 0 ! Number of hard zeros
    integer(hik) :: count_p001  = 0 ! Number of denominators with positive real part and magnitude less than 0.1*de
    integer(hik) :: count_m001  = 0 ! ditto, negative real part
    integer(hik) :: count_p010  = 0 ! ditto, magnitude less than de
    integer(hik) :: count_m010  = 0 ! 
    integer(hik) :: count_p100  = 0 ! ditto, magnitude less than 10*de
    integer(hik) :: count_m100  = 0 ! 
  end type isa_statistics 
  !
  contains
  !
  !  Externally-visible entry points
  !
  subroutine mp2_energy_incore_real(occ_det,occ_ref,eref,eps,int2e,mos,isa_de,eps_big,storage_mode)
    integer(sik), intent(in)               :: occ_det(:)    ! Orbital occupation numbers for the MP2 reference determinant
    real(rk), intent(in)                   :: occ_ref(:)    ! Orbital occupation numbers for the FON-SCF reference determinant
    complex(rk), intent(in)                :: eref          ! Energy of the FON-SCF reference determinant
    complex(rk), intent(in)                :: eps(:)        ! Orbital eigenvalues of the active orbitals
    type(int2e_cache), intent(inout)       :: int2e         ! 2e AO integrals
    complex(rk), intent(in)                :: mos(:,:,:)    ! Bi-orthogonal spin-MOs.
    real(rk), intent(in)                   :: isa_de(:)     ! Intruder-state avoidance parameter(s) list
    real(rk), intent(in), optional         :: eps_big       ! Threshold for picking out important orbitals; default is 1e-2
    character(len=*), intent(in), optional :: storage_mode  ! Accuracy used for storing transformed integrals; can be
                                                            ! 'real', 'quad', or 'as-is'
    !
    include 'mp2_tools_mp2_energy_incore_common.f90'
  end subroutine mp2_energy_incore_real
  !
  subroutine mp2_energy_incore_quad(occ_det,occ_ref,eref,eps,int2e,mos,isa_de,eps_big,storage_mode)
    integer(sik), intent(in)               :: occ_det(:)    ! Orbital occupation numbers for the MP2 reference determinant
    real(xrk), intent(in)                  :: occ_ref(:)    ! Orbital occupation numbers for the FON-SCF reference determinant
    complex(xrk), intent(in)               :: eref          ! Energy of the FON-SCF reference determinant
    complex(xrk), intent(in)               :: eps(:)        ! Orbital eigenvalues of the active orbitals
    type(int2e_cache), intent(inout)       :: int2e         ! 2e AO integrals
    complex(xrk), intent(in)               :: mos(:,:,:)    ! Bi-orthogonal spin-MOs.
    real(xrk), intent(in)                  :: isa_de(:)     ! Intruder-state avoidance parameter(s) list
    real(xrk), intent(in), optional        :: eps_big       ! Threshold for picking out important orbitals; default is 1e-2
    character(len=*), intent(in), optional :: storage_mode  ! Accuracy used for storing transformed integrals; can be
                                                            ! 'real', 'quad', or 'as-is'
    !
    include 'mp2_tools_mp2_energy_incore_common.f90'
  end subroutine mp2_energy_incore_quad
  !
  !  Internal subroutines
  !
  subroutine incore_mp2_single_shift_real(nocc,nvir,ndlt,occ,vir,dlt,occ_det,occ_ref,eps,isa_de,m2e_dd,m2e_vo,m2e_ov, &
                                          eps_bigplayer,e_mp2s,e_mp2d,isa_stat,bigplayer)
    real(rk), intent(in)                :: occ_ref(:)     ! Orbital occupation numbers for the FON-SCF reference determinant
    complex(rk), intent(in)             :: eps(:)         ! Orbital eigenvalues of the active orbitals
    real(rk), intent(in)                :: isa_de         ! Value of the intruder-state avoidance parameter
    real(rk), intent(in)                :: eps_bigplayer  ! Actual numerical threshold for eps_bigplayer
    complex(rk), intent(out)            :: e_mp2s         ! MP2(singles) energy correction; this must be zero if there are no fractional occupations
    complex(rk), intent(out)            :: e_mp2d         ! MP2(doubles) energy correction
    !
    include 'mp2_tools_incore_mp2_single_shift_common.f90'
  end subroutine incore_mp2_single_shift_real
  !
  subroutine incore_mp2_single_shift_quad(nocc,nvir,ndlt,occ,vir,dlt,occ_det,occ_ref,eps,isa_de,m2e_dd,m2e_vo,m2e_ov, &
                                          eps_bigplayer,e_mp2s,e_mp2d,isa_stat,bigplayer)
    real(xrk), intent(in)               :: occ_ref(:)     ! Orbital occupation numbers for the FON-SCF reference determinant
    complex(xrk), intent(in)            :: eps(:)         ! Orbital eigenvalues of the active orbitals
    real(xrk), intent(in)               :: isa_de         ! Value of the intruder-state avoidance parameter
    real(xrk), intent(in)               :: eps_bigplayer  ! Actual numerical threshold for eps_bigplayer
    complex(xrk), intent(out)           :: e_mp2s         ! MP2(singles) energy correction; this must be zero if there are no fractional occupations
    complex(xrk), intent(out)           :: e_mp2d         ! MP2(doubles) energy correction
    !
    include 'mp2_tools_incore_mp2_single_shift_common.f90'
  end subroutine incore_mp2_single_shift_quad
  !
  subroutine isa_report(isa_de,isa_stat)
    real(rk), intent(in)             :: isa_de      ! ISA parameter
    type(isa_statistics), intent(in) :: isa_stat    ! ISA statistics
    !
    if (isa_de<=0) then
      return
    end if
    write (out,"()")
    write (out,"('                  MP2 intruder-state avoidance parameter: ',f12.8,' Hartree')") isa_de
    write (out,"('                                   Total number of terms: ',i0)") isa_stat%terms
    write (out,"('                     Number of hard zeros in denominator: ',i0)") isa_stat%hard_zeros
    write (out,"('Number of denominators -10.0*isa_de <= de <  -1.0*isa_de: ',i0)") isa_stat%count_m100
    write (out,"('Number of denominators  -1.0*isa_de <= de <  -0.1*isa_de: ',i0)") isa_stat%count_m010
    write (out,"('Number of denominators  -0.1*isa_de <= de <   0         : ',i0)") isa_stat%count_m001
    write (out,"('Number of denominators            0 <  de <=  0.1*isa_de: ',i0)") isa_stat%count_p001
    write (out,"('Number of denominators   0.1*isa_de <  de <=  1.0*isa_de: ',i0)") isa_stat%count_p010
    write (out,"('Number of denominators   1.0*isa_de <  de <= 10.0*isa_de: ',i0)") isa_stat%count_p100
    write (out,"()")
  end subroutine isa_report
end module mp2_tools
