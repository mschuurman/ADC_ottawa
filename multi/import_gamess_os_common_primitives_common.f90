!   subroutine os_common_primitives_real(r_l,z_l,r_r,z_r,xi,zeta,p,s00,s00a)
!     real(ark), intent(in)            :: r_l(:)  ! Centre of the left b.f.
!     real(ark), intent(in)            :: z_l     ! Orbital exponent of the left b.f.
!     real(ark), intent(in)            :: r_r(:)  ! ditto, for the right b.f.
!     real(ark), intent(in)            :: z_r     ! 
!     real(ark), intent(out), optional :: xi      ! Eq. 13 of Obara-Saika
!     real(ark), intent(out), optional :: zeta    ! Eq. 14
!     real(ark), intent(out), optional :: p(:)    ! Eq. 15
!     real(ark), intent(out), optional :: s00     ! Eq. 22
!     real(ark), intent(out), optional :: s00a    ! Ahlrichs' version of s00 - no prefactor
!     !
      if (present(xi)  ) xi   = z_l*z_r/(z_l+z_r)
      if (present(zeta)) zeta = z_l + z_r
      if (present(p)   ) p    = (z_l*r_l + z_r*r_r)/(z_l+z_r)
      if (present(s00) ) s00  = (real(pi,kind=kind(z_l))/(z_l+z_r))**real(1.5_xrk,kind=kind(z_l)) &
                                                       * exp(-(z_l*z_r/(z_l+z_r))*sum((r_l-r_r)**2))
      if (present(s00a)) s00a =                          exp(-(z_l*z_r/(z_l+z_r))*sum((r_l-r_r)**2))
!   end subroutine os_common_primitives_real
