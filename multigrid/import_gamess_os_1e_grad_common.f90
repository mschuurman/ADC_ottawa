!   subroutine os_1e_grad_real(ic,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
!     integer(ik), intent(in) :: ic                  ! Gradient component to evaluate
!     integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
!     real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
!     real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
!     real(ark), intent(in)   :: r_l(:)              ! Centre of the left b.f.
!     real(ark), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
!     real(ark), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
!     real(ark), intent(in)   :: z_r                 ! 
      !
      real(kind(z_l)) :: zeta, p(3), s
      integer(ik)     :: m_l, m_r, id1, id2
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,zeta=zeta,p=p)
      !
      right: do m_r=0,p2_r
        left: do m_l=0,p2_l
          s = -2*z_r*(p(ic)-r_r(ic)) * vb_s(m_l,m_r)
          id1 = drop_xyz(m_r,ic)
          if (id1>=0) s = s + (ang_nxyz(m_r,ic)*z_l/zeta) * vb_s(m_l,id1)
          id2 = drop_xyz(m_l,ic)
          if (id2>=0) s = s - (ang_nxyz(m_l,ic)*z_r/zeta) * vb_s(id2,m_r)
          vb_p(m_l,m_r) = s
        end do left
      end do right
!   end subroutine os_1e_grad_real
