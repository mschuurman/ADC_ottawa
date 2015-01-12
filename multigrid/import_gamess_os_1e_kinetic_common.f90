!   subroutine os_1e_kinetic_real(vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
!     integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
!     real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
!     real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
!     real(ark), intent(in)   :: r_l(:)              ! Centre of the left b.f.
!     real(ark), intent(in)   :: z_l                 ! Orbital exponent of the left b.f.
!     real(ark), intent(in)   :: r_r(:)              ! ditto, for the right b.f.
!     real(ark), intent(in)   :: z_r                 ! 
!     !
      real(kind(z_l)) :: zeta, xi, p(3), s00, s
      integer(ik)     :: m_l, m_r, ic, id1, id2, id3
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,xi=xi,zeta=zeta,p=p,s00=s00)
      !
      !  S-S kinetic energy; O&S eq. A13
      !
      vb_p(0,0) = xi*(3-2*xi*sum((r_l-r_r)**2))*s00
      !
      !  Bootstrap: calculate all overlaps on the right, keeping left
      !             at zero. This is a special case of eq. A12
      !
      right_bootstrap: do m_r=1,p2_r
        find_right: do ic=1,3
          id1 = drop_xyz(m_r,ic)
          if (id1<0) cycle find_right
          id2 = drop_xyz(id1,ic)
          s = (p(ic)-r_r(ic)) * vb_p(0,id1)
          s = s + 2*xi * vb_s(0,m_r)
          if (id2>=0) then
            s = s + (ang_nxyz(id1,ic)/(2*zeta)) * vb_p(0,id2)
            s = s - (xi/z_r)*ang_nxyz(id1,ic) * vb_s(0,id2)
          end if
          vb_p(0,m_r) = s
          exit find_right
        end do find_right
      end do right_bootstrap
      !
      !  Now the general case: always recurse on the left. This is the 
      !                        full eq. A12 of Obara-Saika
      !
      right_loop: do m_r=0,p2_r
        left_recurrence: do m_l=1,p2_l  ! m_l=0 is already done
          find_left: do ic=1,3
            id1 = drop_xyz(m_l,ic)
            if (id1<0) cycle find_left
            id2 = drop_xyz(id1,ic)
            id3 = drop_xyz(m_r,ic)
            s = (p(ic)-r_l(ic)) * vb_p(id1,m_r)
            s = s + 2*xi * vb_s(m_l,m_r)
            if (id2>=0) then
              s = s + (ang_nxyz(id1,ic)/(2*zeta)) * vb_p(id2,m_r)
              s = s - (xi/z_l)*ang_nxyz(id1,ic) * vb_s(id2,m_r)
            end if
            if (id3>=0) then
              s = s + (ang_nxyz(m_r,ic)/(2*zeta)) * vb_p(id1,id3)
            end if
            vb_p(m_l,m_r) = s
            exit find_left
          end do find_left
        end do left_recurrence
      end do right_loop
!   end subroutine os_1e_kinetic_real
