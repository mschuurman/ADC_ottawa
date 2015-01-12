!
!  Multigrid test - Sample orbital on a specified surface; for the exclusive
!  use of M.Yu.Ivanov. Since this routine uses a trivial amount of memory,
!  there is no detailed error reporting.
!
  module sample_orbital
    use accuracy
    use fields
    use timer
    use import_gamess
    implicit none
    private
    public run_sample
    !
    !  ==== Compile-time parameters ====
    !
    integer(ik), parameter :: buffer_size     = 10000         ! Grid-point buffer
    !
    !  ==== User-adjustable parameters =====
    !
    integer(ik)            :: verbose         = 0             ! Verbosity level
    character(len=100)     :: orbital_file    = 'orbital.dat' ! Name of the file containing MO coefficients
    integer(ik)            :: orbital_index   = 2             ! Orbital components to load from "orbital_file"
    !
    !  ==== The rest of the parameters will be calculated internally ====
    !
    !  ==== All user-controllable input parameters are in the namelist below ====
    !
    namelist /sample/ verbose, orbital_file, orbital_index
    !
    !  ==== End of global data ====
    !
    contains
    !
    !  Load target Gamess orbital
    !
    subroutine load_gamess_mos(coord,wf)
      real(rk), intent(in)     :: coord(:,:)  ! Coordinates at which we want the orbitals
      real(rk), intent(out)    :: wf   (:)    ! Data field for the orbital values
      !
      integer(ik)              :: npts
      integer(ik)              :: nnuc, iat
      real(rk), allocatable    :: xyzq(:,:)
      logical, save            :: first_run = .true.
      real(rk), allocatable    :: tmp_coord(:,:,:,:)
      complex(rk), allocatable :: tmp_data (:,:,:,:)
      !
      !  The interface to GAMESS import routine is for a 3D grid, so this
      !  will get a little messy. 
      !
      npts = size(coord,dim=2)
      if (npts/=size(wf)) stop 'load_gamess_mos - mismatched arrays'
      allocate (tmp_coord(3,npts,1,1),tmp_data(npts,1,1,1))
      !
      tmp_coord(:,:,1,1) = coord(1:3,:)
      call gamess_load_orbitals(orbital_file,(/orbital_index/),(/1_ik/),tmp_coord,tmp_data)
      wf = real(tmp_data(:,1,1,1),kind=rk)
      !
      deallocate (tmp_coord,tmp_data)
      !
      if (first_run .or. verbose>=1) then
        first_run = .false.
        write (out,"('Loaded orbital ',i4,' from ',a)") orbital_index, trim(orbital_file)
        !
        call gamess_report_nuclei(nnuc)
        write (out,"('Data file contained ',i5,' nuclei')") nnuc
        allocate (xyzq(4,nnuc))
        call gamess_report_nuclei(nnuc,xyzq)
        write (out,"()")
        write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
        write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
        write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
        print_atoms: do iat=1,nnuc
          write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") xyzq(4,iat), xyzq(1:3,iat), xyzq(1:3,iat)*abohr
        end do print_atoms
        write (out,"()")
        deallocate (xyzq)
      end if
    end subroutine load_gamess_mos
    !
    !  Problem driver
    !
    subroutine run_sample
      integer(ik)        :: info
      character(len=30)  :: label(  buffer_size)
      real(rk)           :: coord(3,buffer_size)
      real(rk)           :: wf(buffer_size)
      integer(ik)        :: ipt, npt
      !
      call TimerStart('Orbital Sample')
      call accuracyInitialize
      !
      !  Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=sample,iostat=info)
      if (info/=0) then
        write (out,"(/'**** ENCOUNTERED ERROR ',i4,' READING INPUT NAMELIST ****'/)") info
      end if
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=sample)
      write (out,"()")
      !
      write (out,"((' # ',a16,2x,3(1x,a16),2x,a24))") &
         ' Label ', ' X (Bohr) ', ' Y (Bohr) ', ' Z (Bohr) ', ' Orbital [Bohr^(-3/2)] ', &
         '-------', '----------', '----------', '----------', '-----------------------'
      npt = 0
      read_points: do 
        read(input,*,iostat=info) label(npt+1), coord(1:3,npt+1)
        if (info==0) npt = npt+1
        if (info==0 .and. npt<buffer_size-1) cycle read_points
        call load_gamess_mos(coord(:,1:npt),wf(1:npt))
        print_points: do ipt=1,npt
          write (out,"(' @ ',a16,2x,3(1x,f16.9),2x,g24.14)") trim(label(ipt)), coord(1:3,ipt), wf(ipt)
        end do print_points
        if (info/=0) exit read_points
      end do read_points
      !
      call TimerStop('Orbital Sample')
      call TimerReport
    end subroutine run_sample

  end module sample_orbital
!
  subroutine driver
    use sample_orbital

    call run_sample
  end subroutine driver
