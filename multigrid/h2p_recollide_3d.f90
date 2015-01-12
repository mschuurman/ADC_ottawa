!
!  3D multigrid sample problem - H2(1+) ion in laser field, using Cartesian
!  coordinates. This program is a back-port of the cylindrical h2p_recollide
!  example, and is trying to remain compatible with it as far as possible
!
  module cart_h2p
    use accuracy
    use multigrid
    use fields
    use qmech
    use timer
    use liu
    use lanczos
    use symmetry
    implicit none
    private
    public h2p
!
    integer(ik), parameter :: n_total_max =  6        ! Number of extra scratch fields,
                                                      ! used for solving stationary problem
    integer(ik), parameter :: max_checkpoints = 10    ! Max. number of checkpoints
    integer(ik), parameter :: max_nuclei      = 5     ! Max. number of nuclei allowed
    integer(ik)            :: f_list   (n_total_max)  ! Fields directory
    integer(ik)            :: n_scratch               ! Number of fields available for scratch
    integer(ik)            :: f_psi                   ! Current wavefunction
    integer(ik)            :: f_psi2                  ! Previous wavefunction
    integer(ik)            :: f_start                 ! Initial wavefunction (ionization fraction)
                                                      ! If continous_reference is not blank, this
                                                      ! function will be loaded from disk instead
!
    integer(ik)            :: visible_box = 1         ! Which of the simulation boxes to plot?
!
!  ===== Simulation parameters, which allow meaningful changes            =====
!  ===== Values below are just the defaults - use namelist to change them =====
!
!  Plotting and checkpointing options
!
    logical                :: plot_difference  =   .false.      ! Subtract the reference first
    integer(ik)            :: plot_each        =   100          ! Reduce plotting frequency by
    integer(ik)            :: n_checkpoints    = 0              ! Number of checkpoints
    real(ark)              :: checkpoint_times(max_checkpoints) ! Times at which full checkpoints will be
                                                                ! created
!                                                 
!  Field, pulse and propagation parameters        
!                                                 
    integer(ik)            :: n_nuclei         = 0            ! Number of nuclei
    real(rk)               :: nuc_xyzq(4,max_nuclei)          ! Coordinates and charges of the nuclei
!
    real(rk)               :: dt               =   0.0020_rk  ! Time step, in atomic units
    real(ark)              :: tMin             =   0.0000_ark ! Time at which simulation starts
    real(ark)              :: tMax             = 110.0000_ark ! Time at which simulation ends
    real(ark)              :: tProject         =  55.12_ark   ! Zero electric field at omega=0.057
                                                              ! tProject < tMin means do not project
!
!  Pulse parameters. We have two choices of pulse functions
!
    character(len=clen)    :: pulse_form        = 'sine'               ! Choice of the pulse shape: 'sine'
    real(ark)              :: e0_IR(3)          = (/0,0,1/)* 0.10_ark  ! Amplitude of IR electric field
    real(ark)              :: omega_IR          =   0.057_ark          ! Frequency of the IR field
    real(ark)              :: aField0(3)        = (/0,0,1/)*0.0_ark    ! Initial value of the gauge potential
!
!  Grid parameters                                
!                                                 
    real(rk)               :: box_l(3) = (/ 30._rk, 30._rk, 30._rk /)  ! Box extent 
    integer(ik)            :: box_n(3) = (/ 200, 200, 200 /)           ! Number of points  in the grid
!                                                 
!  Parameters for the initial solver
!
    integer(ik)            :: root_count       = 1           ! Number of roots to solve for
    integer(ik)            :: root_use         = 1           ! Root to use for the initial w.f.
    real(rk)               :: cnv_e            =   1e-10_rk  ! Required energy change
    real(rk)               :: cnv_psi          =   1e-5_rk   ! Required w.f. change
    integer(ik)            :: max_lanczos      =   100       ! Max number of Lanczos restarts
!
!  Quantum particle mass, in electron mass units
!
    real(rk)               :: mass             =   1.0_rk    ! This is just an electron!
!
!  W.f. checkpointing parameters. Two fields will be written/loaded - one
!  at the current time, the second at the previous time. This is the minimum
!  needed for checkpointing with leap-frog. The wavefunctions will be identified
!  by suffixes "_tzero.ckpt" and "_tminus.ckpt", respectively.
!
    character(len=clen)    :: initial_wf_prefix   = ' '       ! Blank means use internal guess
                                                              ! 'gauss' means use gaussian packet
                                                              ! 'plane' means use "planewave" packet
    character(len=clen)    :: final_wf_prefix     = ' '       ! Blank means do not checkpoint
    character(len=clen)    :: final_reference     = ' '       ! Blank means do not re-analyze with
                                                              ! reference subtracted.
    character(len=clen)    :: continuous_reference= ' '       ! Blank means do not subtract reference
                                                              ! during the analysis.
!
!  All parameters - collected in a namelist
!
    namelist /h2p_irxuv/                                         &
                     tMin, tMax, tProject,                       &
                     aField0, e0_IR, omega_IR,                   &
                     pulse_form,                                 &
                     box_l, box_n,                               &
                     plot_each,                                  &
                     n_nuclei, nuc_xyzq,                         &
                     dt,                                         &
                     root_count, root_use,                       &
                     cnv_e, cnv_psi, max_lanczos,                &
                     initial_wf_prefix, final_wf_prefix,         &
                     final_reference,                            &
                     continuous_reference,                       &
                     plot_difference,                            &
                     n_checkpoints, checkpoint_times
!
!  ===== Everything below here should NOT be changed casually =====
!
!  Time propagation parameters
!
    integer(ik)            :: tstep                          ! Current time step
    real(ark)              :: time                           ! Current time
    real(rk)               :: totalNorm          = 1.0_rk    ! Current normalization
    real(rk)               :: eField(3)          = 0.0_rk    ! Current electric field
    real(ark)              :: aField(3)          = 0.0_rk    ! Current gauge field
    logical                :: projected          = .false.   ! Projection of the initial state done
!
!  Nuclear charge of the atom
!
!
    contains

    subroutine h2p
      character(len=clen) :: title, ckpt_name
      real(rk)            :: eref, energy, wf_norm, wf_norm2, w_start, w1, w2
      real(rk)            :: dipole(3)
      complex(rk)         :: c_start
      integer(ik)         :: f_pot, ios
      integer(ik)         :: icheck
      !
      ! Read and echo input parameters. Don't you love namelists?
      !
      read (input,nml=h2p_irxuv,iostat=ios)
      write (out,"(/'==== Simulation parameters ====')")
      write (out,nml=h2p_irxuv)
      write (out,"()")
      if (n_checkpoints>max_checkpoints) stop 'Too many checkpoints - recompile'
      !
      call TimerStart('h2p')
      call accuracyInitialize
      !
      if (tProject<tMin) projected = .true.
      !
      call buildGridFields
      call FieldSetOuterWall('REFLECTING')
      !
      !  Fix Hamiltonian parameters
      !
      call set_nuclear_potential                ! Cannibalized from fock.f90
      call FLsetField((/0._rk,0._rk,0._rk/))
      ! call SMsetSymmetry('None') ! Very bad convergence!
      call SMsetSymmetry('Diatomic',main_axis_order=3, &
                         main_axis=(/1._rk,0._rk,0._rk/), &
                         inversion_centre=(/0._rk,0._rk,0._rk/))
      !
      ! Do we have a checkpoint available?
      !
      if (initial_wf_prefix==' ') then
        !
        !  Nope, no checkpoint - solve stationary problem first. The wavefunction is 
        !  expected in f_psi
        !
        call solveStationaryProblem(eref)
        call report_etv
        call FieldSetOuterWall('REFLECTING')
        !
        !  We know that f_psi is a stationary solution in the absence of
        !  the field. Therefore, we can use analytical solution to determine
        !  its phase at the previous time step - this should (slightly)
        !  improve numerical stability of propagation.
        !
        call FieldCopy(src=f_psi,dst=f_psi2)
        call FieldScale(dst=f_psi2,con=exp(-cmplx(0._rk,eref*(-dt),kind=rk)))
      else
        !
        !  User thinks he has a checkpoint (hur hur), so let's try to
        !  load it, and see how this fails.
        !
        call FieldImport('binary',trim(initial_wf_prefix)//'_tzero.ckpt' ,(/f_psi /),(/1/))
        call FieldImport('binary',trim(initial_wf_prefix)//'_tminus.ckpt',(/f_psi2/),(/1/))
        write (out,"('Loaded ',a)") trim(initial_wf_prefix)
        call report_etv
      end if
      !
      f_pot     = f_list(n_scratch) 
      n_scratch = n_scratch - 1
      f_start   = f_list(n_scratch) 
      n_scratch = n_scratch - 1
      aField    = aField0
      time      = tMin 
      tstep     = 0
      wf_norm   = totalNorm
      !
      if (continuous_reference==' ') then
        call FieldCopy(src=f_psi,dst=f_start)
      else
        call FieldImport('binary',trim(continuous_reference),(/f_start/),(/1/))
        call QMNormalize(f_start,1.0_rk,w1)
        write (out,"(/'Continuous reference wavefunction ',a,' had norm ',g14.7/)") &
               trim(continuous_reference), w1
      end if
      !
      !
      write (out,"(/'Beginning time evolution'/)")
      call TimerStart('Time evolution')
      !
      write (out,"(2x,a10,1x,a12,1x,a18,3(1x,a12),3(1x,a12),1x,a14,1x,a14)") &
             'timestep', 'time [au]', ' energy[H] ', ' E-X[au] ', ' E-Y[au] ', ' E-Z[au] ', &
              'dip-X[e-Bohr]', 'dip-Y[e-Bohr]', 'dip-Z[e-Bohr]', &
              ' <psi|psi> ', '<psi|psi(t=0)>'
      time_loop: do while(time<tMax)
        !
        !  Propagate wavefunction
        !
        eField = evaluateExternalField(time)
        call FLsetField(ef=eField)
        call FieldInit(f_pot,FLtotalPotential,mask=f_psi,grace=2)
        !
        !  We have to be careful with times to which various parameters
        !  refer. performTimeStep returns total energy for the wavefunction
        !  -before- the time step, but the norm for the wavefunction -after-
        !  the time step. This makes reporting and checkpointing a little
        !  tricky.
        !
        wf_norm2 = wf_norm
        call performTimeStep(f_pot,dt,eval=energy,norm=wf_norm)
        !
        !  Special case - cancel out remaining population in the initial state
        !
        if (.not.projected .and. time>=tProject) then
          call cancelInitialState
          projected = .true.
        end if
        !
        !  Weight of the initial wavefunction
        !
        c_start = FieldConjgIntegrate(f_start,f_psi2)
        w_start = abs(c_start)**2
        !
        !  Dipole moment - we need this for harmonics
        !
        call FieldNormMultipoles(f_psi2,mult=dipole)
        !
        !  Report a few results, and possibly save the wavefunction for
        !  later analysis. We need to save f_psi2 - f_psi now refers to
        !  next moment in time.
        !
        write (title,"(i10,1x,f12.4,1x,g18.10,3(1x,g12.6),3(1x,g12.6),1x,g14.8,1x,g14.8)") &
               tstep, time, energy, eField, dipole, wf_norm2, w_start
        write (out,"(' * ',a)") trim(title)
        if (mod(tstep,plot_each)==0) then
          call Visualize(f_psi2,trim(title),eField,plot_difference)
        end if
        !
        !  Advance time, and repeat
        !
        time  = time + dt 
        tstep = tstep + 1
        aField = aField - eField*dt
        !
        !  Test for possible checkpoints
        !
        test_for_checkpoint: do icheck=1,n_checkpoints
          if ( .not.((time-dt)<checkpoint_times(icheck) .and. time>=checkpoint_times(icheck)) ) &
            cycle test_for_checkpoint
          write (ckpt_name,"('checkpoint_',i0,'_step_',i0)") icheck, tstep
          call checkpointWavefunction(ckpt_name)
        end do test_for_checkpoint
      end do time_loop
      call Visualize(f_psi2,"Final wavefunction",(/ eField, 0._rk, 0._rk /),plot_difference)
      !
      !  Checkpoint the wavefunction at the end of propagation
      !
      if (final_wf_prefix/=' ') then
        call checkpointWavefunction(final_wf_prefix)
      end if
      !
      !  If you are going to use checkpoint from previous run as a reference here, 
      !  you should probably choose the '_tminus' one!
      !
      if (final_reference/=' ') then
        call FieldImport('binary',trim(final_reference),(/f_start/),(/1/))
        w1      = FieldNorm(f_psi2)
        w2      = FieldNorm(f_start)
        c_start = FieldConjgIntegrate(f_start,f_psi2)
        w_start = (abs(c_start)/(w1*w2))**2
        write (out,"(//'Comparing final wavefunction at timestep ',i10,' to ',a)") tstep-1, trim(final_reference)
        write (out,"('N(fin) = ',g16.10,' N(Ref) = ',g16.10,' <ref|fin> = ',2g16.10,' weight = ',g16.10)") &
               w1, w2, c_start, w_start
        call FieldAXPY(alpha=-c_start/w2**2,src=f_start,dst=f_psi2)
        call Visualize(f_psi2,"Final wavefunction, with reference subtracted",(/ eField, 0._rk, 0._rk /),.false.)
      end if
      !
      !  Release field for the potential and initial w.f. mask
      !
      n_scratch = n_scratch + 2
      !
      call TimerStop('Time evolution')
      call TimerStop('h2p')
      call TimerReport
    end subroutine h2p
    !
    subroutine set_nuclear_potential
      type(NucleiT) :: nuc(n_nuclei)
      integer(ik)   :: inuc
      !
      fill_nuclei: do inuc=1,n_nuclei
        nuc(inuc)%xyz    = nuc_xyzq(1:3,inuc)
        nuc(inuc)%charge = nuc_xyzq(  4,inuc)
        nuc(inuc)%nelec  = 1.0_rk / n_nuclei
        nuc(inuc)%width  = (1._rk/3._rk) * product(FieldGridSpacing())**(1._rk/3._rk)
      end do fill_nuclei
      !
      call FLsetNuclei(nuc)
      !
      write (out,"()")
      write (out,"('     Number of nuclei is ',i5)") n_nuclei
      write (out,"(' Total nuclear charge is ',f10.3)") sum(nuc_xyzq(4,1:n_nuclei))
      write (out,"()")
      !
      write (out,"()")
      write (out,"(      t8,a36,t48,a36)") 'Coordinates (Bohr)    ', 'Coordinates (Angstrom)    '
      write (out,"(      t8,a36,t48,a36)") '------------------    ', '----------------------    '
      write (out,"(t2,a5,t8,3a12,t48,3a12)") 'ZNUC', '  X  ', '  Y  ', '  Z  ', '  X  ', '  Y  ', '  Z  '
      print_atoms: do inuc=1,n_nuclei
        write (out,"(t2,f5.2,t8,3f12.5,t48,3f12.5)") nuc_xyzq(4,inuc), nuc_xyzq(1:3,inuc), nuc_xyzq(1:3,inuc)*abohr
      end do print_atoms
      write (out,"()")
    end subroutine set_nuclear_potential
    !
    subroutine checkpointWavefunction(name)
      character(len=clen), intent(in) :: name ! Base name of the checkpoint
      !
      write (out,"(//'Checkpointing wavefunction for time step ',i12,' time ',f20.10)") tstep, time
      write (out,"(  'Base name for the checkpoint is ',a//)") trim(name)
      call FieldExport('binary',f_psi, trim(name)//'_tzero.ckpt' )
      call FieldExport('binary',f_psi2,trim(name)//'_tminus.ckpt')
    end subroutine checkpointWavefunction
    !
    function evaluateExternalField(time) result(e)
      real(ark), intent(in) :: time
      real(ark)             :: e(3)
      !
      select case(pulse_form)
        case default
          write (out,"()") trim(pulse_form)
          stop 'evaluateExternalField - bad pulse'
        case ('sine')
          e = evaluateExternalFieldSine(time)
      end select
    end function evaluateExternalField
    !
    function evaluateExternalFieldSine(time) result(e)
      real(ark), intent(in) :: time
      real(ark)             :: e(3)
      !
      e = e0_IR*sin(omega_IR*time)
    end function evaluateExternalFieldSine
    !
    subroutine solveStationaryProblem(e)
      real(rk), intent(out) :: e
      !
      !  In 2D, we had to use Lanczos, which has lower requirements
      !  for quantities needed, but converges relatively slowly. In
      !  3D, we can use Liu/Davidson. Provided that it actually works ...
      !
      call TimerStart('Stationary problem')
      call solveStationaryProblemDavidson(e)
      ! call solveStationaryProblemLanczos(e)
      !
      call Visualize(f_psi,"Initial wavefunction",(/ 0.0_rk, 0.0_rk, 0.0_rk /),.false.)
      !
      call TimerStop('Stationary problem')
      call TimerReport
    end subroutine solveStationaryProblem
    !
    subroutine solveStationaryProblemDavidson(e)
      real(rk), intent(out) :: e
      !
      real(rk)              :: wf_norm
      real(rk)              :: eval(root_count)
      !
      !  Solve for the eigenvalues
      !
      call FieldSetOuterWall('REFLECTING')
      call LUeigenvalues(mass,(/f_psi,f_psi2,f_list(:n_scratch)/),eval, & 
                         FLnuclearPotential,FLpsiGuess,FLnuclearPotential)
      !
      !  Get the desired eigenvector
      !
      call LUeigenvector(root_use,f_psi)
      e       = eval(root_use)
      wf_norm = FieldNorm(f_psi)
      !
      write (out,"(' Using eigenstate ',i4,' E= ',g16.10)") root_use, e
      write (out,"(' Norm of the solution wavefunction = ',f15.10)") wf_norm
      !
      !  Release scratch space
      !
      call LUeigenvector(0,0)
      !
    end subroutine solveStationaryProblemDavidson
    !
    subroutine solveStationaryProblemLanczos(e)
      real(rk), intent(out) :: e
      integer(ik)           :: pass
      integer(ik)           :: f_pot             ! Field index for the potential
      real(rk)              :: eval(root_count)  ! Current and previous energies
      logical               :: c1, c2            ! Converged flags
      !
      !
      !  Prepare external potential and scratch
      !
      f_pot = f_list(n_scratch) ; n_scratch = n_scratch - 2
      call FieldInit(f_pot,FLnuclearPotential)
      call FieldInit(f_psi,FLpsiGuess)
      !
      if (n_scratch<0) then
        stop 'solveStationaryProblemLanczos - not enough scratch'
      end if
      !
      c1 = .false.
      passes: do pass=1,max_lanczos
        c2 = c1
        call LCeigenvalues(mass,f_pot,f_psi,f_psi2,f_list(n_scratch+1),eval,cnv_e,c1, &
                           f_scr_=f_list(:n_scratch))
        call LCeigenvector(f_psi,f_list(n_scratch+1),a_root=root_use)
        write (out,"('**** Lanczos iteration ',i2,' returned ',6g20.10)") pass, eval
        call report_etv
        write (out,"()")
        if (c2 .and. c1) exit passes
      end do passes
      e = eval(root_use)
      !
      !  Release scratch
      !
      n_scratch = n_scratch + 2
      !
    end subroutine solveStationaryProblemLanczos
    !
    !  performTimeStep implicitly operates on f_psi and f_psi2
    !                  Both fields are assumed to be properly normalized
    !                  upon entry (to totalNorm)
    !
    subroutine performTimeStep(f_pot,dt,eval,norm)
      integer(ik), intent(in) :: f_pot  ! Field index containing the potential
      real(rk), intent(in)    :: dt     ! Desired time step
      real(rk), intent(out)   :: eval   ! Total energy of f_psi wavefunction
      real(rk), intent(out)   :: norm   ! Initial norm of the propagated W.F.
      !
      integer(ik)             :: f_HPsi ! Field index for the H|psi> result
      !
      if (n_scratch<1) then
         stop 'performTimeStep - not enough scratch fields'
      end if
      !
      !  Evaluate H|psi> for the wavefunction at T=t
      !
      f_HPsi = f_list(n_scratch) 
      n_scratch = n_scratch - 1 
      if (n_scratch<0) then
        stop 'performTimeStep - not enough scratch fields!'
      end if
      call QMHpsi(mass,f_pot,f_psi,f_HPsi)
      eval = QMExpectation(f_psi,f_HPsi)
      !
      !  Propagate the wavefunction at T=t-dt
      !
      call FieldAXPY(cmplx(0._rk,-2._rk,kind=rk)*dt,f_HPsi,f_psi2)
      !
      !  Renormalize the new wavefunction at T=t+dt
      !
      call QMNormalize(f_psi2,totalNorm,norm)
      !
      !  Release scratch
      !
      n_scratch = n_scratch + 1
      !
      !  Swap field for the wavefunctions
      !
      f_HPsi = f_psi ; f_psi = f_psi2 ; f_psi2 = f_HPsi
    end subroutine performTimeStep
    !
    !  cancelInitialState projects overlap with the initial state (in f_start)
    !                     from wavefunctions in f_psi and f_psi2
    !
    subroutine cancelInitialState
      complex(rk) :: cp     ! Projection factor
      real(rk)    :: norm   ! Initial norm of the residual
      !
      write (out,"(/'Removing initial state overlap')")
      cp = FieldConjgIntegrate(f_start,f_psi)
      call FieldAXPY(alpha=-cp,src=f_start,dst=f_psi)
      call QMNormalize(f_psi,totalNorm,norm)
      write (out,"('f_psi  overlap was:',2(1x,g14.7),' residual: ',g14.7)") cp, norm
      !
      cp = FieldConjgIntegrate(f_start,f_psi2)
      call FieldAXPY(alpha=-cp,src=f_start,dst=f_psi2)
      call QMNormalize(f_psi2,totalNorm,norm)
      write (out,"('f_psi2 overlap was:',2(1x,g14.7),' residual: ',g14.7)") cp, norm
    end subroutine cancelInitialState
    !
    subroutine Visualize(psi_in,title,efield,difference)
      integer(ik), intent(in)      :: psi_in     ! Wavefunction
      character(len=*), intent(in) :: title      ! Descriptive label
      real(rk), intent(in)         :: efield(3)  ! Driving field
      logical, intent(in)          :: difference ! Subtract reference (in f_start) before plotting
      !
      integer(ik)                  :: psi
      complex(rk)                  :: cp
      !
      integer(ik) :: f_scr
      !
      call TimerStart('Visualization')
      f_scr = f_list(n_scratch) ; n_scratch = n_scratch - 1
      if (n_scratch<0) stop 'Visualize - not enough scratch fields (1)'
      !
      if (.not.difference) then
        psi = psi_in
      else
        psi = f_list(n_scratch) ; n_scratch = n_scratch - 1
        if (n_scratch<0) stop 'Visualize - not enough scratch fields (2)'
        cp      = FieldConjgIntegrate(f_start,psi_in)
        call FieldCopy(src=psi_in,dst=psi)
        call FieldAXPY(alpha=-cp,src=f_start,dst=psi)
      end if
      !
      !  Visualize wavefunction, real space
      !
      call FieldVisualize(0,psi,trim(title),visible_box,efield)
      !
      !  Visualize wavefunction, momentum space. Unfortunately, FFTW is broken
      !  for arrays larger than 2GB. Ooopsie.
      !
      ! call FieldFFT(src=psi,dst=f_scr)
      ! call FieldVisualize(1,f_scr,trim(title),1,efield)
      !
      !  Finalize output
      !
      call FieldShow
      !
      n_scratch = n_scratch + 1
      if (difference) n_scratch = n_scratch + 1
      call TimerStop('Visualization')
    end subroutine Visualize

    subroutine buildGridFields
      integer(ik)       :: i, n_total
      character(len=20) :: tag
      real(rk)          :: box(2,3), spacing(3)
      !
      call TimerStart('Grid initialization')
      n_total = 5
      if (plot_difference) n_total = n_total + 1
      if (n_total>n_total_max) stop 'buildGridFields - catastrophe with fields'
      call MultiGridInit(max_grids=1,max_fields=n_total,nborder=1)
      !
      box(1,:) = -box_l
      box(2,:) =  box_l
      call SimpleGridNew('Total',box_n,box)
      spacing = FieldGridSpacing()
      call FLsetGridParameters(spacing)
      !
      do i=1,n_total
        write (tag,"('W.F. scratch ',i5)") i
        call FieldNew(trim(tag),f_list(i),scratch=.true.,wavefunction=.true.)
      end do
      n_scratch = n_total
      f_psi  = f_list(n_scratch) ; n_scratch = n_scratch - 1
      f_psi2 = f_list(n_scratch) ; n_scratch = n_scratch - 1
      if (n_scratch<2) then
        stop 'cyl - disaster with scratch fields!'
      end if
      !
      call TimerStop('Grid initialization')
    end subroutine buildGridFields
    !
    subroutine report_etv
      real(rk) :: epot, ekin
      !
      call FieldLaplacian(f_psi,f_list(1))
      ekin = -(0.5_rk/mass)*QMExpectation(f_psi,f_list(1))
      epot = FieldScalarIntegrate(f_psi,FLtotalPotentialReal)
      write (out,"('Inital guess E= ',g20.10,' T=',g20.10,' V=',g20.10)") ekin+epot, ekin, epot
    end subroutine report_etv

  end module cart_h2p
!
  subroutine driver
    use cart_h2p

    call h2p
  end subroutine driver
