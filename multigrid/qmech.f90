module qmech
!
!  Common higher-level routines, required for solving stationary and
!  time-dependent Schroedinger equations on grid
!
  use accuracy
  use multigrid
  use timer
  implicit none
  private
  public QMNormalize, QMHpsi, QMExpectation, QMProject

  contains

  !
  !  QMNormalize renormalizes the specified field, to a specified values.
  !
  subroutine QMNormalize(wf,norm_desired,norm_original)
    integer(ik), intent(in) :: wf             ! In: Field, containing wavefunction to normalize
                                              ! Out: this field is updated, with the normalized w.f.
    real(rk), intent(in)    :: norm_desired   ! Desired norm
    real(rk), intent(out)   :: norm_original  ! Actiual norm of the input w.f.

    norm_original = FieldNorm(wf)
    call FieldScale(wf,cmplx(norm_desired/norm_original,kind=rk))
  end subroutine QMNormalize
  !
  !  QMHpsi applies non-relativistic one-particle Hamiltonian operator to a
  !         specified field. The external potential must be supplied as a
  !         scalar field as well.
  !
  !  If energy expectation value is needed, use QMExpectation with the 
  !  wavefunction and output of QMHpsi.
  !
  subroutine QMHpsi(mass,pot,psi,Hpsi)
    real(rk), intent(in)    :: mass           ! Mass of the particle
    integer(ik), intent(in) :: pot            ! Field, containing the potential part of the Hamiltonian
    integer(ik), intent(in) :: psi            ! Field, containing the wavefunction
    integer(ik), intent(in) :: Hpsi           ! Field, containing calculated Hpsi
    !
    !  Kinetic energy part. We assume an electron in atomic units, so the prefactor
    !                       is simply 1/2
    !
    call FieldSetWavefunction(psi, .true. )
    call FieldSetWavefunction(pot, .false.)
    call FieldSetWavefunction(Hpsi,.false.)
    call FieldLaplacian(psi,Hpsi)
    call FieldScale(Hpsi,cmplx(-0.5_rk/mass,0.0_rk,kind=rk))
    !
    !  Potential energy part
    !
    call FieldMulAdd(pot,psi,Hpsi)
    call FieldSanitize(Hpsi)
  end subroutine QMHpsi
  !
  !  QMExpectation calculates the expectation value of an operator, given the
  !  wavefunction, and the result of applying the operator on the -same- 
  !  wavefunction. It is simply a syntactic sugar pill for FieldConjgIntegrate.
  !
  function QMExpectation(psi,Opsi) result(v)
    integer(ik), intent(in) :: psi            ! Field, containing |psi>
    integer(ik), intent(in) :: Opsi           ! Field, containing O|psi>
    real(rk)                :: v
    !
    complex(rk) :: val
    real(rk)    :: im_fraction

    val = FieldConjgIntegrate(psi,Opsi)
    im_fraction = abs(aimag(val)) / abs(val)
    if (im_fraction>10.0d0/spacing(1.0_rk)) then
      write(out,"('Calculated expectation value contains large imaginary part: ',2f15.8)") im_fraction
    end if
    if (im_fraction>1000.0d0/spacing(1.0_rk)) then
      write(out,"('Imaginary part of an expectation value is too large - aborting')")
      stop 'QMHpsi - nonsense expectation value'
    end if
    v = real(val,kind=rk)

  end function QMExpectation
  !
  !  Project out specified states
  !
  subroutine QMProject(cores,psi)
    integer(ik), intent(in) :: cores(:) ! List of core states to project out
    integer(ik), intent(in) :: psi      ! State to project cores from
    !
    integer(ik) :: ic
    complex(rk) :: wgt
    !
    call TimerStart('QMProject')
    project_cores: do ic=1,size(cores)
      wgt = FieldConjgIntegrate(left=cores(ic),right=psi)
    ! write (out,"('core ',i2,' weight ',2g14.7)") ic, wgt
      call FieldAXPY(alpha=-wgt,src=cores(ic),dst=psi)
    end do project_cores
    call TimerStop('QMProject')
  end subroutine QMProject
end module qmech
