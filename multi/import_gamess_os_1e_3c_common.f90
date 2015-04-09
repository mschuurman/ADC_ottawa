!   subroutine os_1e_3c_real(what,p2_l,p2_r,n_rec,vb_p,r_l,z_l,r_r,z_r)
!     type(gam_operator_data), intent(in) :: what             ! Operator parameters
!     integer(ik), intent(in)  :: p2_l, p2_r, n_rec           ! Upper dimensions of the vb_p array
!     real(rk), intent(out)    :: vb_p(0:p2_l,0:p2_r,0:n_rec) ! Buffer used for recurrences, and the final result
!     real(ark), intent(in)    :: r_l(:)                      ! Centre of the left b.f.
!     real(ark), intent(in)    :: z_l                         ! Orbital exponent of the left b.f.
!     real(ark), intent(in)    :: r_r(:)                      ! ditto, for the right b.f.
!     real(ark), intent(in)    :: z_r                         ! 
!     !
      real(kind(z_l)) :: zeta, p(3), c(3), s00, s, t
      integer(ik)     :: l_l, l_r                            ! "Angular momentum" on the the left and right
      integer(ik)     :: m_l, m_r, ic, id1, id2, id3, m_i
      !
      !  Sanity check: if the recursion order enough to reach the desided angular momentum?
      !
      l_l = sum(ang_nxyz(p2_l,:))
      l_r = sum(ang_nxyz(p2_r,:))
      if (n_rec/=l_l+l_r) stop 'import_gamess%os_1e_3c - called with an unreasonable order of recursion'
      !
      call os_common_primitives(r_l,z_l,r_r,z_r,zeta=zeta,p=p,s00a=s00)
      c = real(what%op_xyz,kind=kind(c))
      !
      !  Initialize primitive integrals. This is the only part which depends on the
      !  specific operator; the recursions are identical for all kernels.
      !
      t = zeta * sum((p-c)**2)
      call os_basic_integral(what,zeta,t,n_rec,vb_p(0,0,:))
      !
      !  (0|f(r)|0) integrals
      !
      vb_p(0,0,:) = s00 * vb_p(0,0,:)
      ! write (out,"('Basic integrals = ',10g16.9)") vb_p(0,0,:)
      !
      !  (0|f(r)|n) integrals; special-case recursion.
      !
      right_bootstrap: do m_r=1,p2_r  ! S-functions (m_r=0) are already taken care of
        l_r = sum(ang_nxyz(m_r,:))    ! Total "angular momentum" in the right function
        find_right: do ic=1,3
          id1 = drop_xyz(m_r,ic)
          if (id1<0) cycle find_right
          order_bootstrap: do m_i=0,n_rec-l_r
            s = (p(ic)-r_r(ic)) * vb_p(0,id1,m_i) - (p(ic)-c(ic)) * vb_p(0,id1,m_i+1)
            id2 = drop_xyz(id1,ic)
            if (id2>=0) s = s + (ang_nxyz(id1,ic)/(2*zeta)) * (vb_p(0,id2,m_i) - vb_p(0,id2,m_i+1))
            vb_p(0,m_r,m_i) = s
            ! write (out,"(' set ',3i5,' to ',g16.9)") 0, m_r, m_i, s
          end do order_bootstrap
          exit find_right
        end do find_right
      end do right_bootstrap
      !
      !  Now the general case: always recurse on the left.
      !
      right_loop: do m_r=0,p2_r
        l_r = sum(ang_nxyz(m_r,:))    ! Total "angular momentum" in the right function
        left_recurrence: do m_l=1,p2_l
          l_l = sum(ang_nxyz(m_l,:))    ! Total "angular momentum" in the left function
          find_left: do ic=1,3
            id1 = drop_xyz(m_l,ic)
            if (id1<0) cycle find_left
            order_loop: do m_i=0,n_rec-l_l-l_r
              s = (p(ic)-r_l(ic)) * vb_p(id1,m_r,m_i) - (p(ic)-c(ic)) * vb_p(id1,m_r,m_i+1)
              id2 = drop_xyz(id1,ic)
              if (id2>=0) s = s + (ang_nxyz(id1,ic)/(2*zeta)) * (vb_p(id2,m_r,m_i) - vb_p(id2,m_r,m_i+1))
              id3 = drop_xyz(m_r,ic)
              if (id3>=0) s = s + (ang_nxyz(m_r,ic)/(2*zeta)) * (vb_p(id1,id3,m_i) - vb_p(id1,id3,m_i+1))
              vb_p(m_l,m_r,m_i) = s
              ! write (out,"(' set ',3i5,' to ',g16.9)") 0, m_r, m_i, s
            end do order_loop
            exit find_left
          end do find_left
        end do left_recurrence
      end do right_loop
!   end subroutine os_1e_3c_real
