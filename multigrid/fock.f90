!
!  Construction of the Fock operator and related quantities
!
!  All routines in this module implicitly assume particle mass of 1 and
!  charge of -1.
!
  module fock
    use accuracy
    use multigrid
    use fields
    use qmech
    use timer
    implicit none
    private
    public fock_set_options
    public fock_external_potential, fock_hartree_potential
    public fock_exchange_potential, fock_operator, fock_energy
    !
    !  ==== Global data =====
    !
    integer(ik), parameter :: verbose  = 0       ! Verbosity level 
    real(rk), save         :: fock_eps = 1e-8_rk ! Desired convergence of Poisson equations
    real(rk), save         :: fock_sor = 1.9_rk  ! Over-relaxation rate
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  ==== External interface routines ====
    !
    subroutine fock_set_options(eps,sor_rate)
      real(rk), intent(in), optional :: eps      ! Desired convergence for the potential
      real(rk), intent(in), optional :: sor_rate ! Over-relaxation rate
      !
      if (present(eps)     ) fock_eps = eps
      if (present(sor_rate)) fock_sor = sor_rate
    end subroutine fock_set_options
    !
    subroutine fock_external_potential(v,f_scr,xyzq,efield)
      integer(ik), intent(in)        :: v          ! Field for the external potential
      integer(ik), intent(in)        :: f_scr(:)   ! Fields which can be used for scratch
      real(rk), intent(in)           :: xyzq(:,:)  ! Coordinates and charges of the nuclei
      real(rk), intent(in), optional :: efield(:)  ! Uniform electric field
      !
      integer(ik)  :: n_scr  ! Number of available scratch fields
      integer(ik)  :: scr    ! Field for the total density
      !
      call TimerStart('External potential')
      !
      n_scr = size(f_scr)
      scr   = f_scr(n_scr) ; n_scr = n_scr - 1
      if (n_scr<0) stop 'fock_external_potential - no scratch field for the total density!'
      !
      call FieldZero(v)
      !
      !  Add Coulomb potential due to nuclei
      !
      call add_nuclear_potential(xyzq,v,scr=scr)
      !
      !  Add uniform electric field (if given)
      !
      if (present(efield)) then
        call add_uniform_electric_field(efield,v,scr=scr)
      end if
      !
      call TimerStop('External potential')
    end subroutine fock_external_potential
    !
    subroutine fock_hartree_potential(v,rho,mos,occ)
      integer(ik), intent(in)           :: v       ! Field for the Hartree potential
      integer(ik), intent(in)           :: rho     ! Field for the total density. If mos()
                                                   ! and occ() arguments are not present,
                                                   ! total electron density must be supplied
                                                   ! on input
      integer(ik), intent(in), optional :: mos (:) ! List of fields containing the MOs
      real(rk), intent(in), optional    :: occ (:) ! Occupation numbers of the MOs
      !
      !  Construct the total electron density
      !
      if (present(mos) .and. present(occ)) then
        call build_total_density(mos,occ,rho)
      end if
      !
      !  Calculate the long-range part of the Hartree potential, which will
      !  serve both as the initial guess, and as the boundary conditions.
      !
      call build_long_range_tail(rho,v)
      !
      !  Solve Poisson equations for the density, giving us Hartree potential
      !
      call solve_Poisson_equations(rho,v)
      !
    end subroutine fock_hartree_potential
    !
    !  Evaluate the Fock operator.
    !
    subroutine fock_operator(psi,fpsi,v_ext,v_hartree,mos,f_scr)
      integer(ik), intent(in)        :: psi        ! Input: Right-hand-side orbital
      integer(ik), intent(in)        :: fpsi       ! Output: The result of applying the Fock operator
      integer(ik), intent(in)        :: v_ext      ! Input: External potential
      integer(ik), intent(in)        :: v_hartree  ! Input: Hartree potential
      integer(ik), intent(in)        :: mos (:)    ! List of fields containing the MOs
      integer(ik), intent(in)        :: f_scr(:)   ! Fields which can be used for scratch
      !
      integer(ik) :: n_scr    ! Number of scratch fields available/remaining
      integer(ik) :: v_scr    ! "Scratch potential"
      integer(ik) :: f_xchg   ! Temporary for accumulated exchange contribution
      integer(ik) :: imo, mo
      !
      call TimerStart('Fock operator')
      !
      n_scr = size(f_scr)
      if (n_scr<2) stop 'fock%fock_operator - out of scratch fields!'
      v_scr  = f_scr(n_scr) ; n_scr = n_scr - 1 ;
      f_xchg = f_scr(n_scr) ; n_scr = n_scr - 1 ;
      !
      !  Deal with the (cheapish) local part of the potential first
      !
      call FieldCopy(src=v_ext,dst=v_scr)
      call FieldAXPY(alpha=(1._rk,0._rk),src=v_hartree,dst=v_scr)
      call QMHpsi(mass=1._rk,pot=v_scr,psi=psi,Hpsi=fpsi)
      !
      !  The expensive part involves non-local exchange integrals,
      !  which requires solving many Poisson equations. Ouch.
      !
      call FieldZero(f_xchg)
      exchange_potential: do imo=1,size(mos)
        mo = mos(imo)
        call fock_exchange_potential(psic=mo,psi=psi,v=v_scr,product=f_scr(n_scr))
        call FieldMulAdd(src_a=mo,src_b=v_scr,dst=f_xchg)
      end do exchange_potential
      !
      !  Add the exchange part to the total, and we are done.
      !
      call FieldAXPY(alpha=(-1._rk,0._rk),src=f_xchg,dst=fpsi)
      !
      call TimerStop('Fock operator')
    end subroutine fock_operator
    !
    !  Evaluate total energy of a Hartree-Fock wavefunction
    !
    function fock_energy(v_ext,v_hartree,rho,mos,occ,f_scr) result(etot)
      integer(ik), intent(in)        :: v_ext      ! Input: External potential
      integer(ik), intent(in)        :: v_hartree  ! Input: Hartree potential
      integer(ik), intent(in)        :: rho        ! Input: Total electron density
      integer(ik), intent(in)        :: mos (:)    ! List of fields containing the MOs
      real(rk), intent(in)           :: occ (:)    ! Occupation numbers of the molecular orbitals
      integer(ik), intent(in)        :: f_scr(:)   ! Fields which can be used for scratch
      real(rk)                       :: etot       ! Total energy
      !
      integer(ik)   :: n_scr    ! Number of scratch fields available/remaining
      integer(ik)   :: v_exc    ! Scratch area for the exchange potentian
      integer(ik)   :: rho_prod ! Scratch area for orbital products
      integer(ik)   :: n_mos    ! Number of molecular orbitals
      integer(ik)   :: imo, jmo
      complex(rk)   :: e_term
      real(rk)      :: e_kin, e_pot, e_hartree, e_self, e_exchange
      !
      call TimerStart('Hartree-Fock total energy')
      !
      n_mos = size(mos)
      !
      n_scr = size(f_scr)
      if (n_scr<2) stop 'fock%fock_energy - out of scratch fields!'
      v_exc    = f_scr(n_scr) ; n_scr = n_scr - 1 ;
      rho_prod = f_scr(n_scr) ; n_scr = n_scr - 1 ;
      !
      !  Kinetic energy contribution. 
      !
      e_kin = 0
      kinetic_energy: do imo=1,n_mos
        call FieldLaplacian(src=mos(imo),dst=v_exc)
        e_term = -0.5_rk*FieldConjgIntegrate(left=mos(imo),right=v_exc)
        e_kin  = e_kin + occ(imo)*real(e_term,kind=rk)
      end do kinetic_energy
      !
      !  Potential and Hartree energy contributions.
      !
      e_pot     = real(FieldProductIntegrate(v_ext,rho),kind=rk)
      e_hartree = 0.5_rk*real(FieldProductIntegrate(v_hartree,rho),kind=rk)
      !
      !  Now the expensive part - the exchange energy
      !
      e_self     = 0
      e_exchange = 0
      exchange_i: do imo=1,n_mos
        !
        !  First the self-interaction term. 
        !
        call fock_exchange_potential(mos(imo),mos(imo),v_exc,rho_prod)
        e_term = FieldConjgIntegrate(rho_prod,v_exc)
        e_self = e_self - 0.5_rk*occ(imo)*real(e_term,kind=rk)
        if (verbose>=1) then
          write (out,"(t3,'<',i2,',',i2,'|',i2,',',i2,'>=',2g20.10)") imo, imo, imo, imo, e_term
        end if
        !
        !  Now the "true" exchange terms. The weighting is only OK for
        !  high-spin cases, so beware!
        !
        exchange_j: do jmo=1,imo-1
          call fock_exchange_potential(mos(imo),mos(jmo),v_exc,rho_prod)
          e_term = FieldConjgIntegrate(rho_prod,v_exc)
          e_exchange = e_exchange - min(occ(imo),occ(jmo))*real(e_term,kind=rk)
          if (verbose>=1) then
            write (out,"(t3,'<',i2,',',i2,'|',i2,',',i2,'>=',2g20.10)") imo, jmo, jmo, imo, e_term
          end if
        end do exchange_j
      end do exchange_i
      !
      !  Combine the results, and we are done, hurra!
      !
      etot = e_kin + e_pot + e_hartree + e_self + e_exchange
      !
      if (verbose>=1) then
        write (out,"(/t3,'Total Hartree-Fock energy = ',f16.8,' Hartree')") etot
        write (out,"( t3,'           Kinetic energy = ',f16.8 )") e_kin
        write (out,"( t3,'     E-N potential energy = ',f16.8 )") e_pot
        write (out,"( t3,'     E-E potential energy = ',f16.8 )") e_hartree+e_self+e_exchange
        write (out,"( t3,'           Hartree energy = ',f16.8 )") e_hartree
        write (out,"( t3,'         Self-interaction = ',f16.8 )") e_self
        write (out,"( t3,'          Exchange energy = ',f16.8/)") e_exchange
      end if
      !
      call TimerStop('Hartree-Fock total energy')
    end function fock_energy
    !
    !  Evaluate exchange potential corresponding to a single orbital product
    !
    subroutine fock_exchange_potential(psic,psi,v,product)
      integer(ik), intent(in) :: psic    ! Input: Orbital to be conjugated in the product
      integer(ik), intent(in) :: psi     ! Input: Orbital to be taken as is in the product
      integer(ik), intent(in) :: v       ! Output: Exchange potential
      integer(ik), intent(in) :: product ! Output: Conj(psic)*psi orbital product
      !
      call TimerStart('Exchange potential')
      !
      !  Construct orbital product. We'll use v as a temporary scratch, too.
      !
      call FieldConjugate(src=psic,dst=v)
      call FieldZero(product)
      call FieldMulAdd(src_a=v,src_b=psi,dst=product)
      !
      call build_long_range_tail(product,v)
      !
      !  Solve Poisson equations for the product, giving us the exchange potential
      !  (or Hartree potential, of psic and psi are the same).
      !
      call solve_Poisson_equations(product,v)
      !
      call TimerStop('Exchange potential')
    end subroutine fock_exchange_potential
    !
    !  ==== End of external interface routines ====
    !
    !  Accumulate total electron density
    !
    subroutine build_total_density(mos,occ,rho)
      integer(ik), intent(in)        :: mos (:)    ! List of fields containing the MOs
      real(rk), intent(in)           :: occ (:)    ! Occupation numbers of the MOs
      integer(ik), intent(in)        :: rho        ! Field for the total density
      !
      integer(ik) :: imo
      !
      call TimerStart('Total density')
      call FieldZero(rho)
      sum_mos: do imo=1,size(mos)
        call FieldRhoAccumulate(occ(imo),src=mos(imo),dst=rho)
      end do sum_mos
      call TimerStop('Total density')
    end subroutine build_total_density
    !
    !  Calculate the long-range tail of the potential given electron density
    !  (or orbital product, if this is an exchange integral). We'll use
    !  multipole explansion including terms through quadrupole - so the box
    !  must not be too small!
    !
    subroutine build_long_range_tail(rho,vfar)
      integer(ik), intent(in) :: rho  ! Right-hand side of the Poisson equations,
                                      ! expected to be localized around the origin
      integer(ik), intent(in) :: vfar ! Long-range part of the potential, far away
                                      ! from the origin
      !
      complex(rk)  :: npoles(10)      ! rho field may be complex
      real(rk)     :: width
      !
      !  Boundary conditions for the Poisson solver are a little tricky:
      !  our boxes will never be large enough to use zero boundary. As
      !  a reasonable compromise, use multipole expansion up to quadrupole
      !  - this converges as r**-3, which is probably fast enough.
      !
      call TimerStart('Poisson solver guess')
      npoles(1) = FieldNorm1(rho) ! Monopole
      call FieldNorm1Multipoles(rho,npoles(2:)) ! Dipole and quadripole (tracefull)
      !
      !  We'll need an estimate of the distribution extent, to avoid nasty
      !  singularities at the origin by damping the inverse terms. The
      !  choice of the extent is quite arbitrary, as long as it a) does
      !  not touch the walls; and b) contains all charge distribution
      !  inside. 
      !
      width = 0.0_rk
      if (abs(npoles(1))>1e-2) then
        width = sqrt(abs(sum(npoles(5:7))/npoles(1)))
      end if
      width = min(width,7._rk)  ! Just some estimate of a reasonable size of the molecule/atom
      width = max(width,1._rk)
      !
      !  Prepare for multipole expansion
      !
      call FLsetLREmultipoles(npoles,width)
      call FieldInit(vfar,FLlrePotential)
      !
      call TimerStop('Poisson solver guess')
    end subroutine build_long_range_tail
    !
    !  Construct Hartree potential of the ion core.
    !
    subroutine solve_Poisson_equations(rho,v)
      integer(ik), intent(in) :: rho   ! Input: Field containing right-hand side of
                                       !        the Poisson equations
      integer(ik), intent(in) :: v     ! Input: Initial guess and boundary conditions
                                       ! Output: Converged potential
      integer(ik)   :: iter, iter_p
      real(rk)      :: delta
      !
      !  Run SOR iterations on the potential until it converges ... hopefully
      !
      iter   = 1
      iter_p = 200
      if (verbose>=2) iter_p = 1
      call TimerStart('Poisson solver')
      potential_iterations: do
        call FieldIterationSOR(rho=rho,pot=v,delta=delta,omega=fock_sor,boundary=1)
        !
        !  Produce output each time the number of iterations increases by 20%
        !
        ! write (out,"(1x,'iter = ',i5,' delta = ',e20.10)") iter, delta
        if ((5*(iter-iter_p))/iter>0) then
          iter_p = iter
          if (verbose>=1) then
            write (out,"(t5,'poisson: iter = ',i6,' delta = ',g12.5,' ratio = ',g12.4)") &
                   iter, delta, delta/fock_eps
          end if
        end if
        if (delta<fock_eps) exit potential_iterations
        iter = iter + 1
      end do potential_iterations
      call TimerStop('Poisson solver')
      if (verbose>=0) then
        write (out,"(t3,'Poisson solution converged to ',g12.5,' after ',i8,' iterations')") &
               delta, iter
      end if
      !
    end subroutine solve_Poisson_equations
    !
    !  Add electrostatic potential created by the nuclei to the multiplicative
    !  potential
    !
    subroutine add_nuclear_potential(xyzq,v,scr)
      real(rk), intent(in)    :: xyzq(:,:)  ! Coordinates and charges of the nuclei
      integer(ik), intent(in) :: v          ! Input/Output: multiplicative potential
      integer(ik), intent(in) :: scr        ! Scratch field
      !
      integer(ik)   :: n_nuclei
      type(NucleiT) :: nuc(size(xyzq,dim=2))
      integer(ik)   :: inuc
      !
      call TimerStart('Nuclear potential')
      n_nuclei = size(xyzq,dim=2)
      fill_nuclei: do inuc=1,n_nuclei
        nuc(inuc)%xyz    = xyzq(1:3,inuc)
        nuc(inuc)%charge = xyzq(  4,inuc)
        nuc(inuc)%width  = (1._rk/3._rk) * product(FieldGridSpacing())**(1._rk/3._rk)
      end do fill_nuclei
      !
      call FLsetNuclei(nuc)
      call FLsetGridParameters(FieldGridSpacing())
      call FieldInit(scr,FLnuclearPotential)
      !
      !  Add nuclear potential to the total
      !
      call FieldAXPY((1.0_rk,0.0_rk),src=scr,dst=v)
      call TimerStop('Nuclear potential')
    end subroutine add_nuclear_potential
    !
    !  Add uniform electric field to the multiplicative potential
    !
    subroutine add_uniform_electric_field(e,v,scr)
      real(rk), intent(in)    :: e(:)   ! Applied uniform electric field
      integer(ik), intent(in) :: v      ! Input/Output: multiplicative potential
      integer(ik), intent(in) :: scr    ! Scratch field
      !
      call TimerStart('Uniform electric field')
      call FLsetField(e)
      call FieldInit(scr,FLelectricField)
      call FieldAXPY((1._rk,0._rk),src=scr,dst=v)
      call TimerStop('Uniform electric field')
    end subroutine add_uniform_electric_field
!   !
!   !  Apply Hamiltonian operator to our wavefunction. Because we have
!   !  a partitioned Hamiltonian, orthogonalization constraints, and
!   !  and additional electron source, this routine is a a long sequence 
!   !  of nested projection operators.
!   !
!   subroutine apply_hamiltonian(psi,hermit,not_hermit)
!     integer(ik), intent(in) :: psi        ! Field containing input wavefunction
!                                           ! psi get modified here - contributions from inner
!                                           ! shells are projected out.
!     integer(ik), intent(in) :: hermit     ! Hermitian parts of the result
!     integer(ik), intent(in) :: not_hermit ! Possibly non-Hermitian parts
!     !
!     integer(ik) :: psi_cur, psi_inner, f_scr, f_herm, f_noth
!     integer(ik) :: imo
!     complex(rk) :: ovl
!     complex(rk) :: w_source
!     !
!     if (psi==hermit .or. psi==not_hermit .or. hermit==not_hermit) then
!       stop 'apply_hamiltonian - oops with fields'
!     end if
!     !
!     psi_cur   = f_table(n_free) ; n_free = n_free - 1 
!     psi_inner = f_table(n_free) ; n_free = n_free - 1 
!     f_herm    = f_table(n_free) ; n_free = n_free - 1 
!     f_noth    = f_table(n_free) ; n_free = n_free - 1 
!     f_scr     = f_table(n_free) ; n_free = n_free - 1 
!     if (n_free<0) stop 'apply_hamiltonian - out of scratch fields'
!     !
!     !  Both Hermitian and non-Hermitian contributions to the 
!     !  Hamiltonian are initially zero
!     !
!     call FieldZero(hermit)      ! These are the final contributions
!     call FieldZero(not_hermit)
!     call FieldZero(f_herm)      ! These are contributions which require
!     call FieldZero(f_noth)      ! projection to unoccupied fock-space
!     !
!     !  Part I: Any overlap with occupied orbitals gets penalized 
!     !          by FOCK_VSHIFT. This part is Hermitian. We'll also 
!     !          collect the residual in psi_cur.
!     !
!     call FieldCopy(src=psi,dst=psi_cur)
!     project_right: do imo=1,n_mos-1
!       ovl = FieldConjgIntegrate(left=mo_table(imo),right=psi_cur)
!       call FieldAXPY(-ovl,src=mo_table(imo),dst=psi_cur)
!       call FieldAXPY(fock_vshift*ovl,src=mo_table(imo),dst=hermit)
!       if (verbose>=2) then
!         write (out,"(' Right overlap with MO ',i3,': ',2g14.7,' abs: ',g14.7)") &
!                imo, ovl, abs(ovl)
!       end if
!     end do project_right
!     !
!     !  All contributions below are derived from psi_cur, and
!     !  must be projected to unoccupied Fock-space again.
!     !
!     !  Part II: The source function applies to all w.f. satisfying
!     !           orthogonalization constraints. The source is not
!     !           Hermitian.
!     !
!     src_ovl = FieldConjgIntegrate(left=mo_table(n_mos),right=psi_cur)
!     call FieldAXPY(src_ovl*cmplx(0._rk,src_rate,kind=rk), &
!                    src=mo_table(n_mos),dst=f_noth)
!     if (verbose>=2) then
!       write (out,"(' HOMO overlap for the source: ',2g14.7,' abs: ',g14.7)") &
!              src_ovl, abs(src_ovl)
!     end if
!     !
!     !  Part III: Separate remaining wavefunction into inner 
!     !            and outer regions.
!     !
!   ! call FieldZero(psi_inner)
!   ! call FieldMulAdd(src_a=mask_inner,src_b=psi_cur,dst=psi_inner)
!   ! call FieldAXPY((-1._rk,0._rk),src=psi_inner,dst=psi_cur)
!     !
!     !  Part IV: Inner region matching the HOMO gets the corresponding
!     !           orbital energy mo_eps(n_mos). This part is Hermitian.
!     !           Because it was masked before, we need to mask it again
!     !           on the left to keep things Hermitian.
!     !
!   ! call FieldZero(f_scr)
!   ! ovl = FieldConjgIntegrate(left=mo_table(n_mos),right=psi_inner)
!   ! call FieldAXPY(ovl*mo_eps(n_mos),src=mo_table(n_mos),dst=f_scr)
!   ! call FieldMulAdd(src_a=mask_inner,src_b=f_scr,dst=f_herm)
!   ! if (verbose>=2) then
!   !   write (out,"(' Inner HOMO overlap: ',2g14.7,' abs: ',g14.7)") &
!   !          ovl, abs(ovl)
!   ! end if
!     !
!     !  Part V: The remainder of the inner region gets penalized by
!     !          MASK_VSHIFT. Because we use a constant shift, there
!     !          is no need to project the HOMO out on the left.
!     !          However, we still need to mask to the inner region.
!     !
!   ! call FieldAXPY(-ovl,src=mo_table(n_mos),dst=psi_inner)
!   ! call FieldScale(con=cmplx(mask_vshift,0._rk,kind=rk),dst=psi_inner)
!   ! call FieldMulAdd(src_a=mask_inner,src_b=psi_inner,dst=f_herm)
!     !
!     !  Part VI: The outer region sees the long-range Hamiltonian
!     !           The regular potential gives Hermitian part; the
!     !           sink contributes a non-Hermitian term.
!     !
!   ! if (verbose>=2) then
!   !   write (out,"(' Outer region norm: ',g14.7)") FieldNorm(psi_cur)
!   ! end if
!     call FieldLaplacian(src=psi_cur,dst=f_scr)
!     call FieldScale(con=cmplx(-0.5_rk/mass,0.0_rk,kind=rk),dst=f_scr)
!     call FieldMulAdd(src_a=v_core,src_b=psi_cur,dst=f_scr)
!     call FieldAXPY((1._rk,0._rk),src=f_scr,dst=f_herm)
!   ! call FieldMulAdd(src_a=mask_outer,src_b=f_scr,dst=f_herm)
!     call FieldZero(f_scr)
!     call FieldMulAdd(src_a=v_sink,src_b=psi_cur,dst=f_noth)
!   ! call FieldMulAdd(src_a=v_sink,src_b=psi_cur,dst=f_scr)
!   ! call FieldMulAdd(src_a=mask_outer,src_b=f_scr,dst=f_noth)
!     !
!     !  Part VII: Project contributions in f_herm and f_noth to 
!     !            unoccupied Fock-space again. To be consistent,
!     !            we'll have to run the MO loop in reverse here.
!     !
!     project_left: do imo=n_mos-1,1,-1
!       ovl = FieldConjgIntegrate(left=mo_table(imo),right=f_herm)
!       call FieldAXPY(-ovl,src=mo_table(imo),dst=f_herm)
!       if (verbose>=2) then
!         write (out,"(' Left Hermitian overlap with MO ',i3,': ',2g14.7,' abs: ',g14.7)") &
!                imo, ovl, abs(ovl)
!       end if
!       ovl = FieldConjgIntegrate(left=mo_table(imo),right=f_noth)
!       call FieldAXPY(-ovl,src=mo_table(imo),dst=f_noth)
!       if (verbose>=2) then
!         write (out,"(' Left non-Herm. overlap with MO ',i3,': ',2g14.7,' abs: ',g14.7)") &
!                imo, ovl, abs(ovl)
!       end if
!     end do project_left
!     call FieldAXPY((1._rk,0._rk),src=f_herm,dst=hermit)
!     call FieldAXPY((1._rk,0._rk),src=f_noth,dst=not_hermit)
!     !
!     !  Release scratch
!     !
!     n_free = n_free + 5
!     !
!   end subroutine apply_hamiltonian
  end module fock
