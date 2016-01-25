!   subroutine os_1e_dipole_real(id,vb_s,p2_l,p2_r,vb_p,r_l,z_l,r_r,z_r)
!     integer(ik), intent(in) :: id                  ! Dipole component to evaluate
!     integer(ik), intent(in) :: p2_l, p2_r          ! Upper dimensions of the vb_p array
!     real(rk), intent(in)    :: vb_s(0:p2_l,0:p2_r) ! Buffer containing primitive overlaps
!     real(rk), intent(out)   :: vb_p(0:p2_l,0:p2_r) ! Buffer used for recurrences, and the final result
!     real(ark), intent(in)    :: r_l(:)             ! Centre of the left b.f.
!     real(ark), intent(in)    :: z_l                ! Orbital exponent of the left b.f.
!     real(ark), intent(in)    :: r_r(:)             ! ditto, for the right b.f.
!     real(ark), intent(in)    :: z_r                ! 
      !
      real(kind(z_l)) :: zeta, xi, p(3), s00, s
      integer(ik)     :: m_l, m_r, ic, id1, id2, id3
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,xi=xi,zeta=zeta,p=p,s00=s00)
      !
      !*ps
      ! write (out,"('   xi = ',g12.6,' zeta = ',g12.6,' p = ',3f13.6,' s00 = ',g12.6)") xi, zeta, p, s00
      !
      !  Bootstrap: calculate all overlaps on the right, keeping left
      !             at zero. This is a special case of eq. A7
      !
      vb_p(0,0) = p(id)*s00
      right_bootstrap: do m_r=1,p2_r
        find_right: do ic=1,3
          id1 = drop_xyz(m_r,ic)
          if (id1<0) cycle find_right
          s = (p(ic)-r_r(ic)) * vb_p(0,id1)
          id2 = drop_xyz(id1,ic)
          if (id2>=0) s = s + vb_p(0,id2)*ang_nxyz(id1,ic)/(2*zeta)
          if (ic==id) s = s + vb_s(0,id1)/(2*zeta)
          vb_p(0,m_r) = s
          exit find_right
        end do find_right
      end do right_bootstrap
      !
      !  Now the general case: always recurse on the left. This is the 
      !                        full eq. A7 of Obara-Saika
      !
      right_loop: do m_r=0,p2_r
        left_recurrence: do m_l=1,p2_l
          find_left: do ic=1,3
            id1 = drop_xyz(m_l,ic)
            if (id1<0) cycle find_left
            s = (p(ic)-r_l(ic)) * vb_p(id1,m_r)
            id2 = drop_xyz(id1,ic)
            if (id2>=0) s = s + vb_p(id2,m_r)*ang_nxyz(id1,ic)/(2*zeta)
            id3 = drop_xyz(m_r,ic)
            if (id3>=0) s = s + vb_p(id1,id3)*ang_nxyz(m_r,ic)/(2*zeta)
            if (ic==id) s = s + vb_s(id1,m_r)/(2*zeta)
            vb_p(m_l,m_r) = s
            exit find_left
          end do find_left
        end do left_recurrence
      end do right_loop
!   end subroutine os_1e_dipole_real
