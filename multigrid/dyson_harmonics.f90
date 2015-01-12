!
!  This is a trivial companion tool for post-processing recombination dipoles
!  from the "dyson" runs.
!
  module dyson_harmonics
    use accuracy
    use math
    use fftw
    implicit none
    !
    private
    public calculate_harmonics
    !
    integer(ik), parameter   :: max_components    = 3        ! Max number of dipole components
    integer(ik), parameter   :: max_reclen        = 256      ! Max record length in the data file
    integer(ik), parameter   :: iu_data           = 27       ! Unit for the data file
    integer(ik), parameter   :: iu_report         = 28       ! Unit for the report
    !
    integer(ik)              :: verbose           = 0_ik     ! Verbosity level
    integer(ik)              :: skip_timesteps    = 0_ik     ! Number of time steps to skip at the beginning - may be reset
    integer(ik)              :: timesteps         = 0_ik     ! Number of time steps - this value may be reset while
                                                             ! loading the data
    integer(ik)              :: skip_timesteps_input = 0_ik  ! Number of time steps to skip at the beginning
    integer(ik)              :: timesteps_input   = 0_ik     ! Number of time steps in the input file
    character(len=100)       :: harmonic_operator = 'dipole' ! Can be 'dipole' or 'acceleration'
    logical                  :: subtract_baseline = .true.   ! Subtract linear correction, to make the signal "periodic"
    character(len=100)       :: window_function   = 'Welch'  ! Windowing function to use for spectral estimates. 
                                                             ! Can be either "Welch" or "square".
    character(len=100)       :: detailed_output   = ' '      ! Name of the .table detailed output
    integer(ik)              :: dipole_components = 3        ! Number of spatial components
    real(rk)                 :: pad_to_time       = 0._rk    ! If necessary, pad dipoles array with 
                                                             ! zeros until the specified time.
    integer(ik)              :: n_replicas        = 1        ! Replicate the padded array before FFW
    real(rk)                 :: frequency_cut     = 10._rk   ! Discard Fourier components at higher frequencies
    real(rk)                 :: frequency_tolerance = 1e-3_rk! Maximum tolerance in matching frequencies
    character(len=100)       :: in_prefix         = ' '      ! Required prefix on the input line; skip the rest
    integer(ik)              :: in_time_c         = 12       ! Starting column for the time field
    integer(ik)              :: in_time_e         = 27       ! Final column for the time field
                                                             ! Ditto for the dipole components
    integer(ik)              :: in_dipole_c(max_components) = (/ 131, 148, 165 /) 
    integer(ik)              :: in_dipole_e(max_components) = (/ 146, 163, 180 /)
    !
    real(rk)                 :: timestep                     ! Timestep of the dipole data
    real(rk), allocatable    :: times (:)                    ! Buffer for the times of measurement
    real(rk), allocatable    :: frequencies(:)               ! Buffer for the frequencies
    complex(rk), allocatable :: dipole(:,:)                  ! Buffer for the dipoles/accelerations
    integer(ik)              :: n_harmonics                  ! Number of harmonics in current dataset
    real(rk), allocatable    :: h_frequencies(:)             ! Harmonics frequencies for the current dataset
    complex(rk), allocatable :: h_harmonics(:,:)             ! Complex harmonics amplitudes for the current dataset. 
                                                             ! The first index is the dipole component; 
                                                             ! the second index is the sample index
    integer(ik)              :: ns_harmonics                 ! Number of harmonics in the total spectrum
    real(rk), allocatable    :: s_frequencies(:)             ! Harmonics frequencies for the total spectrum
    complex(rk), allocatable :: s_harmonics(:,:)             ! Complex harmonics amplitudes for the total spectrul
    real(ark)                :: euler_angles(3)              ! Euler angles, what else?
    character(len=256)       :: report_destination = ' '     ! File for the report; stdout if blank
    !
    contains
    !
    !  Load simulation results. There is some flexibility in the format - see in_time_* and in_dipole_* above
    !
    subroutine read_data
      character(len=max_reclen) :: line
      integer(ik)               :: ios, irec, itic, ic
      real(rk)                  :: temp
      !
      if (verbose>0) then
        write (out,"(/'           Loading file: ',a)") trim(detailed_output)
        write (out,"( 'Time steps to skip over: ',i10)") skip_timesteps
        write (out,"( '  Time steps to process: ',i10)") timesteps
      end if
      !
      open (iu_data,form='formatted',status='old',position='rewind', &
                    action='read',recl=max_reclen,file=trim(detailed_output),iostat=ios)
      if (ios/=0) then
        write (out,"('Error ',i8,' trying to open ',a,' for reading')") ios, trim(detailed_output)
        stop 'dyson_harmonics%read_data - missing file'
      end if
      !
      irec = 1
      itic = 0
      load_records: do
        read (iu_data,"(a)",iostat=ios) line
        if (ios<0) exit load_records
        if (ios>0) then
          write (out,"('Error ',i8,' reading line ',i10,' of ',a)") ios, irec, trim(detailed_output)
          stop 'dyson_harmonics%read_data - I/O problem'
        end if
        if (line(1:1)=='#') cycle load_records
        if (in_prefix/=' ') then
          if (in_prefix/=line(:len_trim(in_prefix)))  cycle load_records
        end if
        !
        !  Do we need to skip any timesteps?
        !
        if (skip_timesteps>0) then
          skip_timesteps = skip_timesteps - 1
          cycle load_records
        end if
        !
        !  This is not a comment - parse it. The code below looks more complicated than it is:
        !  we simply assume fixed-form Fortran files. For each field, we cut out the substring,
        !  and let the run-time library to decide on the correct conversion routine.
        !
        itic = itic + 1
        read (line(in_time_c:in_time_e),*,iostat=ios) times(itic)
        if (ios/=0) then
          write (out,"('Error ',i8,' converting ""',a,'"" to a time stamp in line ',i10,' of file ',a)") &
                 ios, line(in_time_c:in_time_e), irec, trim(detailed_output)
          stop 'dyson_harmonics%read_data - bad time'
        end if
        components: do ic=1,dipole_components
          read (line(in_dipole_c(ic):in_dipole_e(ic)),*,iostat=ios) temp
          if (ios/=0) then
            write (out,"('Error ',i8,' converting ""',a,'"" to a dipole component ',i2,' in line ',i10,' of file ',a)") &
                   ios, line(in_dipole_c(ic):in_dipole_e(ic)), ic, irec, trim(detailed_output)
            stop 'dyson_harmonics%read_data - bad time'
          end if
          dipole(itic,ic) = temp ! dipole is defined as complex, but we only read in the real part.
        end do components
        if (itic>=timesteps) exit load_records ! Had enough
      end do load_records
      !
      if (itic/=timesteps) then
        write (out,"('Expected ',i10,' data record, but found only ',i10)") timesteps, itic
        stop 'dyson_harmonics%read_data - still hungry'
      end if
      !
      close (iu_data)
    end subroutine read_data
    !
    subroutine check_time_grid
      integer(ik) :: ic
      real(rk)    :: time
      !
      !  Figure out the time step
      !
      timestep = (times(timesteps)-times(1))/(timesteps-1)
      if (verbose>0) then
        write (out,"('Start time = ',g14.7)") times(1)
        write (out,"('  End time = ',g14.7)") times(timesteps)
        write (out,"(' Time step = ',g14.7)") timestep
      end if
      !
      !  Make sure the time grid is really uniform(ish)
      !
      time_grid: do ic=2,timesteps-1
        time = times(1) + (ic-1)*timestep
        if (abs(time-times(ic))>max(1e-5_rk,spacing(100._rk))) then
          write (out,"('For the sample ',i10,' in ',a,' expected time ',g14.7,', but got ',g14.7)") &
                 ic, trim(detailed_output), time, times(ic)
          stop 'dyson_harmonics%check_time_grid - non-uniform time grid'
        end if
      end do time_grid
    end subroutine check_time_grid
    !
    subroutine baseline
      complex(rk) :: base_rate
      integer(ik) :: ic, its
      !
      if (.not.subtract_baseline) return
      !
      component_loop: do ic=1,dipole_components
        base_rate = (dipole(timesteps,ic)-dipole(1,ic))/(timesteps-1)
        if (verbose>0) then
          write (out,"('For the dipole component ',i2,' baseline drift is ',2(g14.7,1x),' per timestep')") ic, base_rate
        end if
        sample_loop: do its=2,timesteps
          dipole(its,ic) = dipole(its,ic) - (its-1)*base_rate
        end do sample_loop
      end do component_loop
    end subroutine baseline
    !
    subroutine pad_data
      integer(ik)              :: padding_samples  ! Number of padding samples in one replica
      integer(ik)              :: it, ir, lb
      integer(ik)              :: total_samples    ! Total number of smaples in all replicas
      real(rk), allocatable    :: t_times (:)      ! There is no need to worry about frequencies
      complex(rk), allocatable :: t_dipole(:,:)    ! - that array is not allocated yet
      integer(ik)              :: alloc
      !
      !  Figure out the size of the padded data array
      !
      if (pad_to_time<=times(timesteps)) then
        padding_samples = 0
      else
        padding_samples = nint((pad_to_time-times(timesteps))/timestep)
      end if
      if (n_replicas<=0) n_replicas = 1
      total_samples = n_replicas * (timesteps + padding_samples)
      if (total_samples<=timesteps) return ! Nothing to be done
      !
      !  Data reallocation needed
      !
      if (verbose>0) then
        write (out,"('Need ',i10,' padding elements, and ',i3,' replicas.')") &
               padding_samples, n_replicas
        write (out,"('New number of samples is ',i10)") total_samples
      end if
      !
      allocate (t_times(timesteps),t_dipole(timesteps,dipole_components),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating memory in pad_data (1)')") alloc
        write (out,"(' timesteps = ',i10,' dipole_components = ',i3)") timesteps, dipole_components
        stop 'dyson_harmonics%pad_data - out of memory (1)'
      end if
      t_times  = times
      t_dipole = dipole
      deallocate (times,dipole)
      allocate (times(total_samples),dipole(total_samples,dipole_components),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating memory in pad_data (2)')") alloc
        write (out,"(' total_samples = ',i10,' dipole_components = ',i3)") total_samples, dipole_components
        stop 'dyson_harmonics%pad_data - out of memory (2)'
      end if
      dipole = 0
      dipole(:timesteps,:) = t_dipole
      lb = timesteps + padding_samples
      augment_dipoles: do ir=2,n_replicas
        dipole(1+(ir-1)*lb:ir*lb,:) = dipole(1:lb,:)
      end do augment_dipoles
      !
      times (:timesteps) = t_times
      augment_times: do it=timesteps+1,total_samples
        times(it) = times(timesteps) + (it-timesteps) * timestep
      end do augment_times
      !
      deallocate (t_times,t_dipole)
      timesteps = total_samples
    end subroutine pad_data
    !
    subroutine fill_frequency_grid
      real(rk)    :: omega_step, omega_min, omega_max
      integer(ik) :: ic, alloc
      !
      if (allocated(frequencies)) deallocate(frequencies)
      allocate (frequencies(timesteps),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating memory in fill_frequency_grid')") alloc
        write (out,"(' timesteps = ',i10)") timesteps
        stop 'dyson_harmonics%fill_frequency_grid - out of memory'
      end if
      !
      omega_step = twopi / (times(timesteps)-times(1))
      omega_min  = -omega_step * ( 0 + (timesteps-1)/2 )
      omega_max  =  omega_min + omega_step * timesteps
      write (out,"('Frequency range is ',g14.7,' to ',g14.7)") omega_min, omega_max
      write (out,"('Frequency resolution is ',g14.7)") omega_step
      frequency_grid: do ic=1,timesteps
        frequencies(ic) = omega_min + (ic-1)*omega_step
      end do frequency_grid
    end subroutine fill_frequency_grid
    !
    subroutine apply_window_function
      integer(ik) :: is
      real(rk)    :: a, wgt
      !
      if (verbose>0) then
        write (out,"('Applying window function: ',a)") trim(window_function)
      end if
      !
      select case (window_function)
        case default
          write (out,"('Window function ',a,' is not implemented')") trim(window_function)
          stop 'dyson_harmonics%apply_window_function'
        case ('Square','square','SQUARE','None','none','NONE')
          ! Nothing to be done for this one
        case ('Welch','welch','WELCH')
          a = 2._rk/(timesteps-1)
          welch_window: do is=1,timesteps
            wgt = max(0._rk,1._rk - (a*(is-1)-1._rk)**2)
            dipole(is,:) = wgt * dipole(is,:)
          end do welch_window
      end select
    end subroutine apply_window_function
    !
    subroutine fft_data
      call fftw_1d(timesteps,dipole_components,dipole,.false.)
      !
      !  Normalization factor
      !
      dipole = dipole * (timestep/sqrt(twopi))
    end subroutine fft_data
    !
    !  Multiply each component by its squared frequency
    !
    subroutine postprocess_dipole
       integer(ik) :: ic
       !
       do ic=1,dipole_components
         dipole(:,ic) = dipole(:,ic) * frequencies**2
       end do
    end subroutine postprocess_dipole
    !
    !  We do not have to do anything for the acceleration form
    !
    subroutine postprocess_acceleration
    end subroutine postprocess_acceleration
    !
    subroutine process_harmonics
      integer(ik) :: zero, ic, alloc
      !
      zero        = 1+(timesteps-1)/2
      n_harmonics = timesteps-zero
      if (allocated(h_frequencies)) deallocate(h_frequencies,h_harmonics)
      allocate (h_frequencies(n_harmonics),h_harmonics(dipole_components,n_harmonics),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating memory in process_harmonics')") alloc
        write (out,"(' n_harmonics = ',i10,' dipole_components = ',i3)") n_harmonics, dipole_components
        stop 'dyson_harmonics%process_harmonics'
      end if
      !
      !  Special case: zero frequency. 
      !
      h_frequencies(1) = 0._rk
      h_harmonics(:,1) = dipole(zero,:) ! / sqrt(real(n_replicas,kind=rk))
      convert_frequencies: do ic=zero+1,timesteps-1
        h_frequencies(1+ic-zero) = frequencies(ic)
        h_harmonics(:,1+ic-zero) = 2*dipole(ic,:) ! / sqrt(real(n_replicas,kind=rk))
      end do convert_frequencies
    end subroutine process_harmonics
    !
    !  Rotate harmonics in 3D
    !
    subroutine rotate_harmonics_3D
      real(ark)   :: rmat(3,3)   ! Rotation matrix
      integer(ik) :: is
      !
      call MathRotationMatrix(euler_angles,rmat)
      !
      if (verbose>0) then
        write (out,"('Euler angles for dipole rotation:',3(1x,g14.5),' degree')") euler_angles * 180.0_rk/pi
        write (out,"('Euler angles for dipole rotation:',3(1x,g14.7),' rad')") euler_angles
        write (out,"('Rotation matrix: ',3(t20,3(1x,f14.10)/))") transpose(rmat)
      end if
      !
      rotate_samples: do is=1,n_harmonics
        h_harmonics(:,is) = matmul(rmat,h_harmonics(:,is))
      end do rotate_samples
    end subroutine rotate_harmonics_3D
    !
    !  The overall prefactor is from L&L II eq. 67.10, divided by 2 (we keep both
    !  positive and negative frequency components explicitly) is given by:
    !
    !   2/(3*vlight**3*2pi)
    !
    subroutine report_harmonics(nh,frequencies,harmonics)
      integer(ik), intent(in) :: nh             ! Number of harmonics in the dataset
      real(rk), intent(in)    :: frequencies(:) ! Frequency
      complex(rk), intent(in) :: harmonics(:,:) ! Complex harmonic amplitudes
      !
      integer(ik) :: iu          ! I/O unit - could be either out or iu_report
      integer(ik) :: ic, i2, ios
      real(rk)    :: rd(dipole_components), id(dipole_components)
      real(rk)    :: amp(dipole_components), phase(dipole_components)
      !
      if (report_destination/=' ') then
        iu = iu_report
        open (iu,iostat=ios,action='write',position='rewind',status='replace',file=trim(report_destination))
        if (ios/=0) then
          write (out,"('Error ',i6,' opening file ',a,' for writing')") ios, trim(report_destination)
          stop 'dyson_harmonics%report_harmonics - creation failed'
        end if
      else
        iu = out
      end if
      !
      write (iu,"('#',8x,a10,2x,a14,1x,a14,2x,a14)") 'Frequency', 'Amplitude', 'Phase', '...'
      list_frequencies: do ic=1,nh
        if (frequencies(ic)>frequency_cut) exit list_frequencies
        rd    = real (harmonics(:,ic),kind=rk)
        id    = aimag(harmonics(:,ic))
        amp   = sqrt(rd**2+id**2)
        phase = atan2(id,rd)
        write (iu,"(' [DATA] ',1x,f10.5,3(2x,g14.7,1x,g14.7))") &
               frequencies(ic), (amp(i2), phase(i2), i2=1,dipole_components)
      end do list_frequencies
      !
      if (report_destination/=' ') then
        close (iu)
        write (out,"('Saved harmonic spectrum to ',a)") trim(report_destination)
      end if
      !
    end subroutine report_harmonics
    !
    subroutine load_dipole_file
      integer(ik) :: alloc
      !
      allocate (times(timesteps),dipole(timesteps,dipole_components),stat=alloc)
      if (alloc/=0) then
        write (out,"('Memory allocation failed with error ',i8)") alloc
        write (out,"('  timesteps         = ',i10)") timesteps
        write (out,"('  dipole_components = ',i3)" ) dipole_components
        stop 'dyson_harmonics%load_dipole_file - out of memory'
      end if
      !
      call read_data
      !
      call check_time_grid
      !
      call baseline
      !
      call pad_data
      !
      call fill_frequency_grid
      !
      call apply_window_function
      !
      call fft_data
      !
      select case (harmonic_operator)
        case default
          write (out,"()") trim(harmonic_operator)
          stop 'dyson_harmonics - bad harmonic_operator'
        case ('dipole')
          call postprocess_dipole
        case ('acceleration')
          call postprocess_acceleration
      end select
      !
      call process_harmonics
      !
      deallocate (times,dipole)
    end subroutine load_dipole_file
    !
    !  Accumulate harmonics data
    !
    subroutine accumulate_harmonics
      integer(ik) :: alloc
      integer(ik) :: is, ib
      integer(ik) :: matched, missed_high, skipped   ! Counters
      !
      if (.not.allocated(s_frequencies)) then
        !
        !  The easy part - this is the first batch of updates
        !
        ns_harmonics = n_harmonics
        allocate (s_frequencies(ns_harmonics),s_harmonics(dipole_components,ns_harmonics),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i6,' allocating memory in accumulate_harmonics')") alloc
          write (out,"(' ns_harmonics = ',i10,' dipole_components = ',i3)") ns_harmonics, dipole_components
          stop 'dyson_harmonics%accumulate_harmonics - out of memory'
        end if
        s_frequencies = h_frequencies
        s_harmonics   = h_harmonics
      else
        !
        !  The hard part. Since the dipole profiles may come from different simulations, their
        !                 frequency content may be different. We do not want to do a full 
        !                 resampling of the frequencies; instead, we'll try to match frequencies
        !                 already in the spectrum to those in the incoming buffer. We'll assume
        !                 that the spectrum with the smallest frequency content is loaded first
        !
        matched     = 0
        missed_high = 0
        skipped     = 0
        ib          = 1 ! Next unused position in the incoming buffer
        scan_total: do is=1,ns_harmonics
          !
          !  Find matching frequency in the incoming buffer
          !
          find_buffer_match: do while(abs(h_frequencies(ib)-s_frequencies(is))>frequency_tolerance)
            if (h_frequencies(ib)>s_frequencies(is)+frequency_tolerance) then
              call missing_match
              cycle scan_total
            end if
            if (ib+1>n_harmonics) then
              call missing_match
              cycle scan_total
            end if
            ib      = ib + 1
            skipped = skipped + 1
          end do find_buffer_match
          !
          !  This is a match
          !
          matched           = matched + 1
          s_harmonics(:,is) = s_harmonics(:,is) + h_harmonics(:,ib)
          if (verbose>1) then
            write (out,"('Matched sample ',i8,' frequency ',g14.7,' to buffer sample ',i7,' frequency ',g14.7)") &
                   is, s_frequencies(is), ib, h_frequencies(ib)
          end if
        end do scan_total
        !
        if (verbose>0) then
          !
          !  Note that all matched values are also skipped over - so there is a correction below
          !
          write (out,"('Matched ',i8,' frequencies, discarded ',i8,', and missed ',i8)") matched, skipped-matched, missed_high
        end if
      end if
      !
      contains
        subroutine missing_match
          !
          !  There is no match. This is not a big deal of frequency is above the report
          !  cut-off, but a fatal error below the cutoff.
          !
          if (s_frequencies(is)>frequency_cut) then
            missed_high = missed_high + 1
            return
          end if
          write (out,"('Can''t find frequency match for sample ',i8,' at frequency ',g14.7)") is, s_frequencies(is)
          write (out,"('Total spectrum can not be calculated. Either loosen the frequency_tolerance,')")
          write (out,"('or rearrange the input to process the shortest simulation first.')")
          stop 'dyson_harmonics%accumulate_harmonics'
        end subroutine missing_match
    end subroutine accumulate_harmonics
    !
    !  An extensible driver for harmonics evaluation - includes coherent summation
    !  for directional averaging and/or ensemble averaging. The input is handled
    !  by a (very stupid) interpreter.
    !
    subroutine calculate_harmonics
      character(len=256) :: cmd          ! Command line buffer
      integer(ik)        :: line         ! Input line position
      integer(ik)        :: ios          ! I/O status variable
      integer(ik)        :: ip, ip1, ip2 ! Buffer pointers. cmd(ip1:ip2) is the current command
      integer(ik)        :: ichr
      logical            :: lock_in      ! Data processing started - lock in essential dimensionalities
      complex(rk)        :: scale        ! Scaling factor for the harmonics.
      !
      call accuracyInitialize
      !
      lock_in = .false.
      line    = 0
      command_loop: do
        line = line + 1
        read (input,"(a)",iostat=ios) cmd
        if (ios/=0) exit command_loop
        !
        if (verbose>=0) then
          write (out,"('@ ',a)") trim(cmd)
        end if
        !
        skip_leading_blanks: do ip1=1,len(cmd)
          if ( cmd(ip1:ip1)/=' ' ) exit skip_leading_blanks
        end do skip_leading_blanks
        if (ip1>=len(cmd)) cycle command_loop     ! Empty line, skip over it
        if (cmd(ip1:ip1)=='#') cycle command_loop ! A comment, skip over it
        locate_keyword: do ip2=ip1,len(cmd)-1
          if (cmd(ip2+1:ip2+1)==' ') exit locate_keyword
        end do locate_keyword
        uppercase_keyword: do ip=ip1,ip2
          ichr = iachar(cmd(ip:ip)) ! ASCII collating sequence
          if (ichr>=97 .and. ichr<=122) then  ! lower-case letters
            cmd(ip:ip) = achar(ichr+65-97)    ! convert to upper-case
          end if
        end do uppercase_keyword
        !
        select case (cmd(ip1:ip2))
          case default
            write (out,"('Unrecognized directive in line ',i8,': ',a)") line, cmd(ip1:ip2)
            stop 'dyson_harmonics - bad command'
          !
          !  Declarations - this is just a very complicated way of writing the namelist, 
          !                 with each keyword corresponding to a variable in the static
          !                 data section.
          !
          case ('VERBOSE')
            read (cmd(ip2+1:),*,iostat=ios) verbose
            call parse_error_check
          case ('SKIP_TIMESTEPS')
            read (cmd(ip2+1:),*,iostat=ios) skip_timesteps_input
            call parse_error_check
          case ('TIMESTEPS')
            read (cmd(ip2+1:),*,iostat=ios) timesteps_input
            call parse_error_check
          case ('HARMONIC_OPERATOR')
            read (cmd(ip2+1:),*,iostat=ios) harmonic_operator
            call parse_error_check
          case ('DIPOLE_COMPONENTS')
            call lockin_check
            read (cmd(ip2+1:),*,iostat=ios) dipole_components
            call parse_error_check
            if (dipole_components>max_components) then
              write (out,"('Too many dipole components in input line ',i8)") line
              stop 'dyson_harmonics - parameter out of range'
            end if
          case ('SUBTRACT_BASELINE')
            read (cmd(ip2+1:),*,iostat=ios) subtract_baseline
            call parse_error_check
          case ('WINDOW_FUNCTION')
            read (cmd(ip2+1:),*,iostat=ios) window_function
            call parse_error_check
          case ('PAD_TO_TIME')
            read (cmd(ip2+1:),*,iostat=ios) pad_to_time
            call parse_error_check
          case ('N_REPLICAS')
            read (cmd(ip2+1:),*,iostat=ios) n_replicas
            call parse_error_check
          case ('FREQUENCY_CUT')
            read (cmd(ip2+1:),*,iostat=ios) frequency_cut
            call parse_error_check
          case ('FREQUENCY_TOLERANCE')
            read (cmd(ip2+1:),*,iostat=ios) frequency_tolerance
            call parse_error_check
          case ('IN_PREFIX')
            read (cmd(ip2+1:),*,iostat=ios) in_prefix
            call parse_error_check
          case ('IN_TIME_C')
            read (cmd(ip2+1:),*,iostat=ios) in_time_c
            call parse_error_check
          case ('IN_TIME_E')
            read (cmd(ip2+1:),*,iostat=ios) in_time_e
            call parse_error_check
          case ('IN_DIPOLE_C')
            read (cmd(ip2+1:),*,iostat=ios) in_dipole_c(1:dipole_components)
            call parse_error_check
          case ('IN_DIPOLE_E')
            read (cmd(ip2+1:),*,iostat=ios) in_dipole_e(1:dipole_components)
            call parse_error_check
          !
          !  Operators - these might actually do something
          !
          case ('CLEAR')
            lock_in = .false.
            if (allocated(s_frequencies)) deallocate (s_frequencies,s_harmonics)
          case ('LOAD')
            read (cmd(ip2+1:),*,iostat=ios) detailed_output
            call parse_error_check
            lock_in        = .true.
            timesteps      = timesteps_input
            skip_timesteps = skip_timesteps_input
            call load_dipole_file
          case ('SCALE')
            read (cmd(ip2+1:),*,iostat=ios) scale
            call parse_error_check
            call dataset_present_check
            if (verbose>0) then
              write (out,"('Scaling harmonic amplitudes in the current dataset by ',2g14.7)") scale
            end if
            h_harmonics(:,:) = scale * h_harmonics(:,:)
          case ('ROTATE','ROTATED','ROTATER')  ! Angles in degree or radians (default)
            if (dipole_components/=3) then
              write (out,"('Coordinate rotation encontered in line ',i8,', but the world in not 3D. Oops.')") line
              stop 'dyson_harmonics%command_loop - unimplemented dimension for rotation'
            end if
            read (cmd(ip2+1:),*,iostat=ios) euler_angles
            if (cmd(ip1:ip2)=='ROTATED') then
              euler_angles = euler_angles * pi / 180.0_rk
            end if
            call parse_error_check
            call dataset_present_check
            call rotate_harmonics_3D
          case ('ACCUMULATE')
            call accumulate_harmonics
          case ('REPORT','DATASET')
            report_destination = ' '
            if (len_trim(cmd(ip2+1:))>0) then
              read (cmd(ip2+1:),*,iostat=ios) report_destination
              call parse_error_check
            end if
            if (cmd(ip1:ip2)=='DATASET') then
              call dataset_present_check
              call report_harmonics(n_harmonics,h_frequencies,h_harmonics)
            else ! 'REPORT' - total accumulated spectrum
              call total_present_check
              call report_harmonics(ns_harmonics,s_frequencies,s_harmonics)
            end if
        end select
      end do command_loop
      !
      !  Exceptions and error handling
      !
      contains
        subroutine parse_error_check
          if (ios==0) return ! All is OK
          write (out,"(/'Input line ',i8,': error ',i6,' parsing arguments of ',a)") line, ios, cmd(ip1:ip2)
          write (out,"('The arguments were: ',a)") trim(cmd(ip2+1:))
          stop 'dyson_harmonics - command argument parse error'
        end subroutine parse_error_check
        !
        subroutine lockin_check
          if (.not.lock_in) return ! All is OK
          write (out,"(/'Input line ',i8,' changes a locked-in parameter ',a)") line, cmd(ip1:ip2)
          stop 'dyson_harmonics - invalid input sequence (A)'
        end subroutine lockin_check
        !
        subroutine dataset_present_check
          if (allocated(h_harmonics)) return ! All is OK
          write (out,"(/'Input line ',i8,' operates on a non-existent dataset. Command was: ',a)") line, cmd(ip1:ip2)
          stop 'dyson_harmonics - invalid input sequence (B)'
        end subroutine dataset_present_check
        !
        subroutine total_present_check
          if (allocated(s_harmonics)) return ! All is OK
          write (out,"(/'Input line ',i8,' operates on a non-existent total dataset. Command was: ',a)") line, cmd(ip1:ip2)
          stop 'dyson_harmonics - invalid input sequence (C)'
        end subroutine total_present_check
    end subroutine calculate_harmonics
  end module dyson_harmonics
  !
  subroutine driver
    use dyson_harmonics
 
    call calculate_harmonics
  end subroutine driver
