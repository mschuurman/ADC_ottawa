!   recursive subroutine os_basic_integral_real(what,rho,t,n_max,gn)
!     type(gam_operator_data), intent(in) :: what        ! Operator parameters
!     real(ark), intent(in)               :: rho         ! Effective basis exponent
!     real(ark), intent(in)               :: t           ! Argument of the basic integral
!     integer(ik), intent(in)             :: n_max       ! Maximum order of basic integral needed
!     real(rk), intent(out)               :: gn(0:n_max) ! Basic integrals
!     !
      integer(ik)             :: n, i
      real(kind(t))           :: gx(0:n_max+2)     ! Extra room, in case recursion needs a higher-order term
      type(gam_operator_data) :: wx                ! Operator parameters for recursion
      real(kind(t))           :: rx                ! rho for recursion
      real(kind(t))           :: tx                ! t for recursion
      real(rk), pointer       :: gp(:,:)           ! Gaussian expansion of an operator, only available in
                                                   ! standard precision
      !
      ! write (out,"('os_basic: ',a,' rho = ',g14.7,' t = ',g14.7,' n = ',i4)") trim(what%op_name), rho, t, n_max
      ! write (out,"('os_basic: omega = ',g14.7,' a = ',g14.7)") what%omega, what%imag_rc
      !
      select case (what%op_name(4:))
        case default
          write (out,"('import_gamess%os_basic_integral: Operator ',a,' is not implemented.')") trim(what%op_name(4:))
          stop 'import_gamess%os_basic_integral - unimplemented'
        !
        !  "Primitive" cases first.
        !
        case ('DELTA') ! Delta function - A37
          gn(:) = exp(-t)
        case ('ONE')   ! Unity: an elaborate way to compute 2c overlaps - A38
          gn(0)  = (real(pi_xrk,kind=kind(rho))/rho)**real(1.5_xrk,kind=kind(rho))
          gn(1:) = 0
        case ('1/R')   ! Nuclear attraction - A39
          call os_boys_table(n_max,t,gn)
          gn = (2*real(pi_xrk,kind=kind(rho))/rho) * gn
        case ('R')     ! Linear attraction - A44
          !
          !  Linear attraction operator requires Boys' function table up to order n_max+1
          !
          call os_boys_table(n_max+1,t,gx(0:n_max+1))
          gn(0) = (t+1)*gx(0) - t*gx(1)
          gn_r_sum: do n=1,n_max
            gn(n) = -n*gx(n-1) + (t+n+1)*gx(n) - t*gx(n+1)
          end do gn_r_sum
          gn = (2*real(pi_xrk,kind=kind(rho))/rho**2) * gn
        case ('R**2')  ! Harmonic attraction - A45
          gn(0) = (t + real(1.5_xrk,kind=kind(t)))
          if (n_max>=1) gn(1)  = -1
          if (n_max>=2) gn(2:) =  0
          gn = sqrt(real(pi_xrk,kind=kind(rho))**3/rho**5) * gn
        case ('GAUSS') ! Gaussian - same as 'GAUSS ONE', but likely slightly faster
          rx = real(what%omega,kind=kind(rho))/(rho+real(what%omega,kind=kind(rho)))
          gn_gauss: do n=0,n_max
            gn(n) = rx**n
          end do gn_gauss
          gn = gn * exp(-rx*t) * real(pi_xrk,kind=kind(rho))/ &
               (rho+real(what%omega,kind=kind(rho)))**real(1.5_xrk,kind=kind(rho))
        !
        !  Recursive cases - reduction to simpler Gn(t)
        !
        case ('GAUSS ONE','GAUSS R')
          !
          !  g(r) = exp(-omega*r**2) * gx(r), for an arbitrary gx(r) - A41
          !  There is no increase in the Gn order. Note that the full reduction
          !  below is an overkill if Gx does not have t dependence. However, any
          !  savings are likely trivial.
          !
          wx = what
          wx%op_name(4:) = what%op_name(4+6:)
          rx = rho + real(what%omega,kind=kind(rho))
          tx = t * rho / rx
          call os_basic_integral(wx,rx,tx,n_max,gx(0:n_max))
          !
          !  This is not the most elegant way to organize this loop, but it probably
          !  does not matter ...
          !
          gn_gaussian_recursion: do n=0,n_max
            gn(n) = 0
            gn_gaussian_sum: do i=0,n
              gn(n) = gn(n) + MathBinomial(n,i,rx) * (real(what%omega,kind=kind(rx))/rx)**(n-i) * (rho/rx)**i * gx(i)
            end do gn_gaussian_sum
          end do gn_gaussian_recursion
          gn = gn * exp(-real(what%omega,kind=kind(t))*t/rx)
        !
        !  The next two interactions - 'A/(R**2+A**2)' and 'R/(R**2+A**2)' are
        !  needed to calculate Coulomb potential at a complex position. I can't
        !  find an analytical expression for these kernels, so the code gets a 
        !  little more complicated, and relies on expansion of the 1/(R**2+A**2)
        !  term in Gaussians. Using the 200-term expansion in os_integral_operators,
        !  A values between 0.01 and 100 Bohr, rho values between 0.001 and 100,
        !  and T values from 10^(-3) to 120 yield errors (the smallest of the absolute
        !  and relative error) of better than 1e-7. This is not great, but OK for 
        !  our intended use.
        !
        case ('A/(R**2+A**2)')
          if (kind(t)/=rk) stop 'import_gamess%os_basic_integral - requested integral is not implemented in quad precision'
          call warn_integral_accuracy
          wx = what
          wx%op_name = '3C GAUSS'
          gp => operator_1_1pr2
          gn = 0
          gn_complex_integrate_a: do i=1,size(gp,dim=2)
            wx%omega = gp(1,i)/what%imag_rc**2
            call os_basic_integral(wx,rho,t,n_max,gx(0:n_max))
            gn = gn + gx(0:n_max) * gp(2,i)/real(what%imag_rc,kind=kind(gp))
          end do gn_complex_integrate_a
        case ('R/(R**2+A**2)')
          if (kind(t)/=rk) stop 'import_gamess%os_basic_integral - requested integral is not implemented in quad precision'
          call warn_integral_accuracy
          wx = what
          wx%op_name = '3C GAUSS R'
          gp => operator_1_1pr2
          gn = 0
          gn_complex_integrate_r: do i=1,size(gp,dim=2)
            wx%omega = gp(1,i)/what%imag_rc**2
            call os_basic_integral(wx,rho,t,n_max,gx(0:n_max))
            gn = gn + gx(0:n_max) * gp(2,i)/real(what%imag_rc,kind=kind(gp))**2
          end do gn_complex_integrate_r
      end select
      ! write (out,"('os_basic: gn =',10(1x,g14.7))") gn
      contains
      subroutine warn_integral_accuracy
        integer(ik) :: this_thread
        !
        logical, save :: warned = .false.
        this_thread = 0
        !$ this_thread = omp_get_thread_num()
        if (this_thread/=0 .or. warned) return
        warned = .true.
        write (out,"(/'WARNING: This program utilizes ',a,'-type integrals.')") trim(what%op_name)
        write (out,"( 'WARNING: The numerical accuracy domain for these integrals is restricted.')")
        write (out,"( 'WARNING: Please make sure you understand code in import_gamess%os_basic_integral')")
        write (out,"( 'WARNING: before using these integrals for production calculations.'/)")
      end subroutine warn_integral_accuracy
!   end subroutine os_basic_integral_real
