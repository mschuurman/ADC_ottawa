!   subroutine os_1e_matrix_real(what,l,r,v)
!     type (gam_operator_data), intent(in) :: what   ! Operator to evaluate
!     type (gam_structure), intent(in)     :: l      ! Basis set on the left
!     type (gam_structure), intent(in)     :: r      ! Basis set on the right
!     real(rk), intent(out)                :: v(:,:) ! Buffer for the integrals
      !
      integer(ik)  :: at_l, sh_l, c1_l, c2_l, l_l, p1_l, p2_l
      integer(ik)  :: at_r, sh_r, c1_r, c2_r, l_r, p1_r, p2_r
      !
      !  Make sure neither bra nor ket include a rotation matrix
      !
      if (.not. MathIsUnitMatrix(l%rotmat) .or. .not. MathIsUnitMatrix(r%rotmat) ) then
        stop 'import_gamess%os_1e_matrix - rotation is not implemented'
      end if
      !
      !  Iterate over all shell blocks. Since l and r are not necessarily the same,
      !  matrix elements do not have to be symmetric, and we have to process all
      !  shell blocks explicitly.
      ! 
      p1_r = 1
      r_atom_loop: do at_r=1,r%natoms
        r_shell_loop: do sh_r=1,r%atoms(at_r)%nshell
          !
          !  Get "angular momentum" (this is a Cartesian shell, so not really)
          !  and the block length.
          !
          l_r  = r%atoms(at_r)%sh_l(sh_r)
          p2_r = p1_r + ang_loc(l_r+1) - ang_loc(l_r) - 1
          !
          !  Figure out the range of contractions for this shell
          !
          c1_r = r%atoms(at_r)%sh_p(sh_r  )
          c2_r = r%atoms(at_r)%sh_p(sh_r+1)-1
          p1_l = 1
          l_atom_loop: do at_l=1,l%natoms
            l_shell_loop: do sh_l=1,l%atoms(at_l)%nshell
              l_l  = l%atoms(at_l)%sh_l(sh_l)
              p2_l = p1_l + ang_loc(l_l+1) - ang_loc(l_l) - 1
              c1_l = l%atoms(at_l)%sh_p(sh_l  )
              c2_l = l%atoms(at_l)%sh_p(sh_l+1)-1
              !*ps
              ! write (out,"('doing pair: ',2i4,' shells: ',2i4,' block: ',4i6)") at_l, at_r, sh_l, sh_r, p1_l, p2_l, p1_r, p2_r
              ! write (out,"('atoms at:  '6f12.6)") l%atoms(at_l)%xyz/abohr, r%atoms(at_r)%xyz/abohr
              ! write (out,"('l: ',2i2)") l_l, l_r
              ! write (out,"('zeta_l: ',12g12.6)") l%atoms(at_l)%p_zet(c1_l:c2_l)
              ! write (out,"('zeta_r: ',12g12.6)") r%atoms(at_r)%p_zet(c1_r:c2_r)
              ! write (out,"('   c_l: ',12g12.6)") l%atoms(at_l)%p_c(c1_l:c2_l)
              ! write (out,"('   c_r: ',12g12.6)") r%atoms(at_r)%p_c(c1_r:c2_r)
              call os_1e_contraction(what,v(p1_l:p2_l,p1_r:p2_r),        &
                   real(l%atoms(at_l)%xyz,kind=kind(v))/abohr,l_l,       &
                     real(l%atoms(at_l)%p_zet(c1_l:c2_l),kind=kind(v)),  &
                       real(l%atoms(at_l)%p_c(c1_l:c2_l),kind=kind(v)),  &
                   real(r%atoms(at_r)%xyz,kind=kind(v))/abohr,l_r,       &
                     real(r%atoms(at_r)%p_zet(c1_r:c2_r),kind=kind(v)),  &
                       real(r%atoms(at_r)%p_c(c1_r:c2_r),kind=kind(v)))
              p1_l = p2_l + 1
            end do l_shell_loop
          end do l_atom_loop
          if (p1_l-1/=l%nbasis) then
            stop 'import_gamess%os_1e_matrix - count error (l)'
          end if
          p1_r = p2_r + 1
        end do r_shell_loop
      end do r_atom_loop
      !
      if (p1_r-1/=r%nbasis) then
        stop 'import_gamess%os_1e_matrix - count error (r)'
      end if
!   end subroutine os_1e_matrix_real
