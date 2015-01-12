!   subroutine gamess_1e_integrals_real(what,v,bra,ket,op_xyz,op_param,op_index)
!     character(len=*), intent(in)      :: what        ! What to evaluate
!     real(rk), intent(out)             :: v(:,:)      ! Buffer for the results
!     type(gam_structure), intent(in), target, optional :: bra
!     type(gam_structure), intent(in), target, optional :: ket
!     real(rk), intent(in), optional    :: op_xyz(:)   ! Position of the operator (3c integrals)
!     real(rk), intent(in), optional    :: op_param(:) ! Additional parameters for the operator
!                                                      ! Required size and allowed values depend on
!                                                      ! the operator; see below.
!     integer(ik), intent(in), optional :: op_index(:) ! Additional (integer) parameters for the operator.
      !
      character(len=40)            :: lwhat          ! Local copy of "what"
      type(gam_structure), pointer :: l, r
      type(gam_operator_data)      :: op_data        ! Operator data
      !
      call TimerStart('GAMESS '//trim(what))
      l => gam_def ; r => gam_def
      if (present(bra)) l => bra
      if (present(ket)) r => ket
      !
      !  What is it we want to evaluate?
      !
      lwhat = what
      select case (lwhat(1:3))
        case ('AO ')
          if ( (size(v,dim=1)<l%nbasis) .or. (size(v,dim=2)<r%nbasis) ) then
            write (out,"('import_gamess%gamess_1e_integrals: Output buffer is not large enough')")
            stop 'import_gamess%gamess_1e_integrals: Output buffer is not large enough'
          end if
      end select
      !
      !  3-centre integrals need position of the operator; make sure it was supplied
      !  Operator-specific parameters will be handled later.
      !
      if (lwhat(4:5)=='3C') then
        if (.not.present(op_xyz)) then
          write (out,"('import_gamess%gamess_1e_integrals: operator location missing for ',a)") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals: Missing operator position'
        end if
        op_data%op_xyz = real(op_xyz,kind=kind(op_data%op_xyz))
        ! write (out,"('Operator at ',3f20.12)") op_data%op_xyz
      end if
      !
      op_data%op_name = lwhat(4:)
      !
      select case (lwhat)
        case default
          write (out,"('import_gamess%gamess_1e_integrals_real: ',a,' are not implemented')") trim(lwhat)
          stop 'import_gamess%gamess_1e_integrals - integral not implemented'
        !
        !  2-centre integrals. There is no extra checking and no parameters
        !
        case ('AO OVERLAP','AO DIPOLE X','AO DIPOLE Y','AO DIPOLE Z','AO KINETIC','AO D/DX','AO D/DY','AO D/DZ')
          !
          !  Each evaluation of the dipole integrals will also evaluate (and discard)
          !  the overlap integrals. This is not the most efficient way of doing things,
          !  but do we really care?
          !
          call os_1e_matrix(op_data,l,r,v)
        case ('AO R-D/DR')
          !
          !  Two-centre integarls of operator r_i (d/d r_j).
          !  These integrals are reduced to linear combinations of overlap and dipole integrals
          !  We need two integer indices, both in the range 1-3
          !
          if (.not.present(op_index)) then
            write (out,"('import_gamess%gamess_1e_integrals: indices missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Missing operator parameter'
          end if
          if (size(op_index)/=2) then
            write (out,"('import_gamess%gamess_1e_integrals: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Wrong parameter count'
          end if
          op_data%op_i = op_index(1)
          op_data%op_j = op_index(2)
          if (op_data%op_i<1 .or. op_data%op_i>3 .or. op_data%op_j<1 .or. op_data%op_j>3) then
            write (out,"('import_gamess%gamess_1e_integrals: index parameters out of range for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Parameters out of range'
          end if
          call os_1e_matrix(op_data,l,r,v)
        !
        !  3-centre integrals. We may need to do extra work here.
        !
        !  First the integrals where no additional parmaters are needed:
        !
        !  'AO 3C DELTA'     - Delta-function at C (I know it's silly)
        !  'AO 3C ONE'       - Very complicated 2-centre overlap implementation
        !  'AO 3C 1/R'       - Nuclear attraction
        !  'AO 3C R'         - Linear attraction
        !  'AO 3C R**2'      - Harmonic attraction
        !
        case ('AO 3C DELTA', 'AO 3C ONE', 'AO 3C 1/R', 'AO 3C R', 'AO 3C R**2')
          call os_1e_matrix(op_data,l,r,v)
        !
        !  Operators containing a Gaussian at C; need the exponent.
        !
        !  'AO 3C GAUSS'     - Three-centre overlap
        !  'AO 3C GAUSS ONE' - Ditto, but a more complicated implementation
        !  'AO 3C GAUSS R'   - Gaussian-damped linear attraction
        !
        case ('AO 3C GAUSS', 'AO 3C GAUSS ONE', 'AO 3C GAUSS R')
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_1e_integrals: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Missing operator parameter'
          end if
          if (size(op_param)/=1) then
            write (out,"('import_gamess%gamess_1e_integrals: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Wrong parameter count'
          end if
          op_data%omega = real(op_param(1),kind=kind(op_data%omega))
          call os_1e_matrix(op_data,l,r,v)
        !
        !  Operators with C at a complex point; need the imaginary component
        !
        !  'AO 3C R/(R**2+A**2)' - Real part of the nuclear attraction
        !  'AO 3C A/(R**2+A**2)' - Imaginary part of the nuclear attraction (with negative sign)
        !
        case ('AO 3C R/(R**2+A**2)','AO 3C A/(R**2+A**2)')
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_1e_integrals: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Missing operator parameter'
          end if
          if (size(op_param)/=1) then
            write (out,"('import_gamess%gamess_1e_integrals: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_1e_integrals: Wrong parameter count'
          end if
          op_data%imag_rc = real(op_param(1),kind=kind(op_data%imag_rc))
          call os_1e_matrix(op_data,l,r,v)
      end select
      call TimerStop('GAMESS '//trim(what))
!   end subroutine gamess_1e_integrals_real
