!   subroutine os_1e_contraction_real(what,vb,r_l,l_l,z_l,c_l,r_r,l_r,z_r,c_r)
!     type (gam_operator_data), intent(in) :: what    ! Integral to evaluate
!     real(rk), intent(out)                :: vb(:,:) ! Place for the integrals
!     real(ark), intent(in)                :: r_l(:)  ! Coordinate on the left, in Bohr
!     integer(ik), intent(in)              :: l_l     ! "Angular momentum" on the left
!     real(ark), intent(in)                :: z_l(:)  ! Primitive exponents on the left
!     real(ark), intent(in)                :: c_l(:)  ! Contractions on the left
!     real(ark), intent(in)                :: r_r(:)  ! Ditto on the right
!     integer(ik), intent(in)              :: l_r     ! 
!     real(ark), intent(in)                :: z_r(:)  ! 
!     real(ark), intent(in)                :: c_r(:)  ! 
      !
      integer(ik)                 :: ic_l, p1_l, p2_l, p_l, int_l
      integer(ik)                 :: ic_r, p1_r, p2_r, p_r, int_r
      real(kind(vb)), allocatable :: vb_r(:,:,:)  ! Buffer for recursions
      real(kind(vb))              :: wgt
      integer(ik)                 :: n_rec        ! Number of additional recursion intermediated
      integer(ik)                 :: alloc
      real(kind(vb))              :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
      integer(ik)                 :: out=6
      !
      p1_l = ang_loc(l_l)
      p2_l = ang_loc(l_l+1)-1
      p1_r = ang_loc(l_r)
      p2_r = ang_loc(l_r+1)-1
      !
      select case (what%op_name)
        case default
          !
          !  Assume general 3-centre integral; each unit of angular momentum requires
          !  an extra order in recursion.
          !
          n_rec = l_l + l_r
        case ('OVERLAP')
          n_rec = 0
        case ('DIPOLE X','DIPOLE Y','DIPOLE Z','KINETIC','D/DX','D/DY','D/DZ')
          n_rec = 1
        case ('R-D/DR')
          n_rec = 2
      end select
      allocate (vb_r(0:p2_l,0:p2_r,0:n_rec),stat=alloc)
      if (alloc/=0) then
        write (out,"('import_gamess%os_1e_contraction: allocation failed. Parameters are:',3(1x,i8))") p2_l, p2_r, n_rec
        stop 'import_gamess%os_1e_contraction: Buffer allocation failed'
      end if
      !
      vb = 0
      ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      r_contractions_loop: do ic_r=1,size(z_r)
        l_contractions_loop: do ic_l=1,size(z_l)
          !*ps
          ! write (out,"('doing contraction ',2i3,' c = ',2g13.6)") ic_l, ic_r, c_l(ic_l), c_r(ic_r)
          ! write (out,"('p2 = '2i4,' r_l = ',3f12.5,' l_l = ',i2,' z_l = ',g12.6,' r_r = ',3f12.5,' l_r = ',i2,' z_r = ',g12.6)") &
          !        p2_l,p2_r,r_l,l_l,z_l(ic_l),r_r,l_r,z_r(ic_r)
          !
          !  The only part of the 1e integral code which is property-specific
          !  is the primitive-integral call below. Everything else is identical
          !  for all properties.
          !
          select case (what%op_name)
            case default
              write (out,"('import_gamess%os_1e_contraction: Property ',a,' is not implemented')") trim(what%op_name)
              stop 'import_gamess%os_1e_contraction: Unimplemented property'
            case ('OVERLAP')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('DIPOLE X')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_dipole (1_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('DIPOLE Y')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_dipole (2_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('DIPOLE Z')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_dipole (3_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('D/DX')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_grad   (1_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('D/DY')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_grad   (2_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('D/DZ')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_grad   (3_ik,vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('R-D/DR')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,2),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_dipole (what%op_i, &
                                      vb_r(:,:,2),p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_rddr   (what%op_i,what%op_j, &
                          vb_r(:,:,2),vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('KINETIC')
              call os_1e_overlap(                 p2_l,p2_r,vb_r(:,:,1),r_l,z_l(ic_l),r_r,z_r(ic_r))
              call os_1e_kinetic(     vb_r(:,:,1),p2_l,p2_r,vb_r(:,:,0),r_l,z_l(ic_l),r_r,z_r(ic_r))
            case ('3C 1/R','3C DELTA','3C ONE','3C R','3C R**2','3C GAUSS', &
                  '3C GAUSS ONE', '3C GAUSS R', '3C A/(R**2+A**2)', '3C R/(R**2+A**2)')
              call os_1e_3c(what,p2_l,p2_r,n_rec,vb_r,r_l,z_l(ic_l),r_r,z_r(ic_r))
          end select
          !
          wgt = c_l(ic_l) * c_r(ic_r)
          r_accumulate: do p_r=p1_r,p2_r
            int_r = p_r - p1_r + 1
            l_accumulate: do p_l=p1_l,p2_l
              int_l = p_l - p1_l + 1
              vb(int_l,int_r) = vb(int_l,int_r) + wgt*ang_c_kind(p_l)*ang_c_kind(p_r)*vb_r(p_l,p_r,0)
              !*ps
              ! write (out,"(' int ',2i4,' increased by ',f12.6)") int_l, int_r, wgt*ang_c(p_l)*ang_c(p_r)*vb_p(p_l,p_r)
            end do l_accumulate
          end do r_accumulate
          !
        end do l_contractions_loop
      end do r_contractions_loop
      !
      deallocate (vb_r,stat=alloc)
      if (alloc/=0) then
        stop 'import_gamess%os_1e_contraction: Buffer deallocation failed'
      end if
!   end subroutine os_1e_contraction_real
