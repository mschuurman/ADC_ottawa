!
!  Orientation-independent part of the H2-H2 interaction potential.
!  The potential is based on the 0-0-0 term of Diep and Johnson,
!  JCP 113, 3480 (2000). The tabulated potential was fit to a
!  sum of Morse potential and the capped 6-12 potential [Zhechkov
!  et al, JCTC 1, 841 (2005)]. The fit is accurate to better than
!  10 cm^-1 within the entire attractive region, deteriorating to
!  50 cm^-1 in the strongly repulsive region. It also behaves 
!  sensibly in the entire range of internuclear disrances.
!
module h2_potential
  use accuracy
  implicit none
  private
  public h2_v12
  !
  contains
  !
  function h2_v12(rb) result(v12)
    real(rk), intent(in) :: rb  ! Distance between centres of mass, in bohr
    real(rk)             :: v12 ! Calculated interaction potential, in hartree
    !
    real(rk), parameter :: wexp = 5.13653454390905_rk
    real(rk), parameter :: aexp = 2.105248788432718_rk
    real(rk), parameter :: rexp = 3.5409282325241667_rk
    real(rk), parameter :: dvdw = 34.97278714201808_rk
    real(rk), parameter :: rvdw = 3.3229678275098564_rk
    !
    real(rk), parameter :: r0   = rvdw * 2._rk**(-1._rk/6._rk)
    real(rk), parameter :: u0   = 396._rk/25._rk
    real(rk), parameter :: u1   =-(672._rk/25._rk)*2._rk**(5._rk/6._rk)
    real(rk), parameter :: u2   = (552._rk/25._rk)*2._rk**(2._rk/3._rk)
    !
    real(rk) :: r, earg, rr
    !
    !  The fit is in Kelvin, with the argument in Angstrom - convert units
    !
    r = max(rb*abohr,2.0_rk)  ! The tabulated potential stops at 2.0 Angstrom
    !
    !  Morse potential first
    !
    earg = -aexp * (r-rexp)
    if (earg<-max_exp) then
      v12 = 0.0_rk
    else 
      v12 = wexp * (1._rk - exp(earg))**2 - wexp
    end if  
    !
    !  Now add the truncated 6-12
    !
    if (r>=r0) then
      rr = (rvdw/r)**6 
      v12 = v12 + dvdw * (-2._rk*rr + rr**2)
    else
      rr = (r/rvdw)**5
      v12 = v12 + dvdw * (u0 + u1*rr + u2*rr**2)
    end if
    !
    !  Switch back to atomic units
    !
    v12 = v12 * k_Boltzmann
  end function h2_v12
end module h2_potential
!
!program test_v12
!  use accuracy
!  use h2_potential
!  real(rk) :: r, v, v_ref
!  !
!  call accuracyInitialize
!  do 
!    read(5,*,end=100) r, v_ref
!    v = h2_v12(r)
!    call compare('V',v_ref,v)
!  end do
!  100 continue
!  !
!  contains
!    subroutine compare(lbl,ref,val)
!      character(len=*), intent(in) :: lbl
!      real(rk), intent(in)         :: ref, val
!      !
!      if (abs(ref-val)<=max(1e-14_rk*abs(ref),100*spacing(abs(ref)))) return
!      write (6,"('ERROR AT R=',g14.7,1x,a,': REF=',g20.12,' VAL=',g20.12,' ERR=',2g20.12)") &
!             r, lbl, ref, val, val-ref, abs((val/ref-1.0_rk))
!    end subroutine compare
!end program test_v12
