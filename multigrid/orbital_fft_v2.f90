!
!  Evaluation of matrix elements of a planewave.
!
!  This is the second version of the code, using analytical evaluation of
!  molecular integrals. 
!
  module orbital_fft_v2
    use accuracy
    use timer
    use import_gamess
    !$ use OMP_LIB

    implicit none

    private

    public start

    !
    !  ==== Fixed parameters. These can't be changed without recompilation ====
    !
    integer(ik), parameter :: max_orbitals   = 1000
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)           :: verbose         = 2                 ! Verbosity level
    character(len=100)    :: orbitals_file   = 'orbitals.dat'    ! Name of the file containing the MOs
    integer(ik)           :: orbital_count  = -1                 ! Number of molecular orbitals to build integrals for
                                                                 ! -1 means all
    integer(ik)           :: orbital_list(max_orbitals)          ! List of orbital indices from the input file
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    type(gam_structure)      :: mol                              ! Gamess basis set and orbital data
    complex(rk), allocatable :: ao_ints(:,:)                     ! Integrals in the AO basis
    complex(rk), allocatable :: mo_ints(:,:)                     ! Integrals in the MO basis
    real(rk), allocatable    :: ao_over(:,:)                     ! Overlaps in the AO basis
    real(rk), allocatable    :: mo_over(:,:)                     ! Overlaps in the MO basis
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /intdata2/ verbose, &
                       orbitals_file, &
                       orbital_count, orbital_list
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Report molecular geometry - lifted from dyson_tools.f90
    !
    subroutine report_structure(level,name)
      integer(ik), intent(in)      :: level ! Minimal verbosity level to print the structure;
                                            ! Report just the number of atoms otherwise
      character(len=*), intent(in) :: name  ! Name of the structure/file
      !
      integer(ik)           :: nnuc, iat
      real(rk), allocatable :: xyzq(:,:)
      !
      call gamess_report_nuclei(nnuc,structure=mol)
      write (out,"('Data file ',a,' contained ',i5,' nuclei')") trim(name), nnuc
      if (level>verbose) return
      !
      !  Tell a bit more!
      !
      allocate (xyzq(4,nnuc))
      call gamess_report_nuclei(nnuc,xyzq,structure=mol)
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do iat=1,nnuc
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
      end do print_atoms
      write (out,"()")
      deallocate (xyzq)
    end subroutine report_structure
    !
    !  Prepare ion-specific potentials: core+Hartree and transition potentials
    !
    subroutine load_molecular_data
      integer(ik) :: iorb
      !
      call TimerStart('Load molecular data')
      write (out,"(/'Loading ',a)") trim(orbitals_file)
      !
      call gamess_load_orbitals(file=trim(orbitals_file),structure=mol)
      !
      write (out,"(/'Number of AOs = ',i6 )") mol%nbasis
      write (out,"( 'Number of MOs = ',i6/)") mol%nvectors
      !
      if (orbital_count<=0) then
        orbital_count = mol%nvectors
        fill_orbitals: do iorb=1,orbital_count
          orbital_list(iorb) = iorb
        end do fill_orbitals
      end if
      !
      allocate (ao_ints(mol%nbasis,  mol%nbasis))
      allocate (mo_ints(mol%nvectors,mol%nvectors))
      allocate (ao_over(mol%nbasis,  mol%nbasis))
      allocate (mo_over(mol%nvectors,mol%nvectors))
      !
      call report_structure(0,trim(orbitals_file))
      !
      call TimerStop('Load molecular data')
    end subroutine load_molecular_data
    !
    subroutine check_overlaps
      integer(ik) :: imo, jmo
      real(rk)    :: val, eval
      !
      call TimerStart('MO overlap check')
      call gamess_1e_integrals('AO OVERLAP',ao_over,bra=mol,ket=mol)
      mo_over = matmul(transpose(mol%vectors(:,1:mol%nvectors)),matmul(ao_over,mol%vectors(:,1:mol%nvectors)))
      test_jmo: do jmo=1,mol%nvectors
        test_imo: do imo=1,mol%nvectors
          val  = mo_over(imo,jmo)
          eval = 0.0_rk ; if (imo==jmo) eval = 1.0_rk
          if (abs(val-eval)>=1e-6_rk) then
            write (out,"('Overlap integral between MOs ',i5,' and ',i5,' is in error by ',g22.14)") imo, jmo, val-eval
          end if
        end do test_imo
      end do test_jmo
      call TimerStop('MO overlap check')
    end subroutine check_overlaps
    !
    subroutine evaluate_integrals(kvec)
      real(rk), intent(in) :: kvec(3) ! K-vector of the planewave
      !
      call TimerStart('MO planewave integrals')
      call gamess_1e_integrals('AO PLANEWAVE',ao_ints,bra=mol,ket=mol,op_param=kvec)
      mo_ints = matmul(transpose(mol%vectors(:,1:mol%nvectors)),matmul(ao_ints,mol%vectors(:,1:mol%nvectors)))
      call TimerStop('MO planewave integrals')
    end subroutine evaluate_integrals
    !
    subroutine start
      integer(ik) :: info
      integer(ik) :: bra, i_bra, ket, i_ket
      real(rk)    :: kvec(3)
      complex(rk) :: val, ovr
      !
      call TimerStart('start')
      !  Read and echo input parameters. Don't you love namelists?
      read (input,nml=intdata2,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=intdata2)
      write (out,"()")
      !
      write (out,"('Default integer kind = ',i5)") ik
      write (out,"('Default real kind    = ',i5)") rk
      !
      call load_molecular_data
      !
      call check_overlaps
      !
      write (out,"(/'Done with preliminaries'/)")
      !
      call TimerReport
      !
      !  We are ready to take matrix elements!
      !
      write (out,"(('#',3(1x,a14),2(1x,a5),4(1x,a20)))") &
      '   kx   ', '   ky   ', '   kz   ', ' bra ', ' ket ', '  Re[<|e^ikr|>]  ', '  Im[<|e^ikr|>]  ', '  Re[<|>]  ', ' Im[<|>]  ', &
      '  ----  ', '  ----  ', '  ----  ', '-----', '-----', ' --------------- ', ' --------------- ', ' --------- ', '--------- '
      !
      read_kvecs: do
        read(input,*,iostat=info) kvec
        if (info/=0) exit read_kvecs
        !
        call evaluate_integrals(kvec)
        loop_bra: do i_bra=1,orbital_count
          bra = orbital_list(i_bra)
          loop_ket: do i_ket=1,i_bra
            ket = orbital_list(i_ket)
            val = mo_ints(bra,ket)
            ovr = mo_over(bra,ket)
            write (out,"('@',3(1x,g14.6),2(1x,i5),4(1x,g20.12))") kvec, bra, ket, val, ovr
          end do loop_ket
        end do loop_bra
      end do read_kvecs
      write (out,"()")
      !
      call TimerStop('start')
      call TimerReport
    end subroutine start
    
  end module orbital_fft_v2
  !
  !
  !
  subroutine driver
    use accuracy
    use math
    use orbital_fft_v2
    !
    real(rk) :: dummy
    !
    call accuracyInitialize
    dummy = MathDoubleFactorial(100_ik)
    dummy = MathFactorial(100_ik)
    dummy = MathLogDoubleFactorial(1000_ik)
    dummy = MathLogFactorial(1000_ik)

    call start

  end subroutine driver
