!
!  Multigrid test - hydrogen molecular ion
!
  module h_module
    use accuracy
    use multigrid
    use qmech
    use fields
    use lanczos
    use liu
    use pulse
    use vandalyse
    use symmetry
    use timer
    implicit none
!
    integer(ik), parameter :: n_scratch   = 2         ! Number of extra scratch fields,
                                                      ! used for solving stationary problem
    integer(ik)            :: pot1, Hpsi, psi1, psi2
    integer(ik)            :: extra_scr(n_scratch)    ! Extra scratch fields
    type(NucleiT)          :: nucxyzqw(2) = (/   &
!
!  X,Y,Z are coordinates of the nucleus, with nuclear charge Q
!  Vx,Vy,Vz are velocity of the nucleus
!  W     is the minimum distance for the regularized Coulomb potential
!  NE    is the number of electrons, assigned to this nucleus
!  mass  is the mass of the nucleaus in a.u. (1.007825)
!
!  R0(H2+) is 1.05713 Angstrom, or 1.99768 Bohr (aug-ccPVQZ basis)
!  This corresponds to non-relativistic total energy of -1.1031117
!
!  The first few excited states (CIS/aug-ccPVQZ)
!       -0.667295
!       -0.428183 (2x)
!       -0.359242
!       -0.254442 (2x)
!
!  Note that w is not arbitrary, and should be adjusted according to
!  the grid spacing close to the nucleus. The magic expression is
!  (see comments in fields.f90):
!
!     w = (2/3) (ax*ay*az)**(1/3)
!
!  where ax, ay, and az are grid spacings along each coordinate
!  direction.
!
!                   X         Y         Z             Vx       Vy      Vz          Q        W       NE    mass
! Parallel to Z
!   NucleiT( (/ 0.00000,  0.00000, -0.998844 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.178,  1.00, 1836.15 ), &
!   NucleiT( (/ 0.00000,  0.00000,  0.998844 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.178,  1.00, 1836.15 )  &
! Perpendicular to Z
    NucleiT( (/-0.998844, 0.00000, -0.000000 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.222,  1.00, 1836.15 ), &
    NucleiT( (/ 0.998844, 0.00000,  0.000000 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.222,  1.00, 1836.15 )  &
!   NucleiT( (/-2.000000, 0.00000, -0.000000 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.167,  1.00, 1836.15 ), &
!   NucleiT( (/ 2.000000, 0.00000, -0.000000 /), (/ 0.0000,  0.0000, 0.0000 /),  1.00000,  0.167,  1.00, 1836.15 )  &
! Perpendicular to Z
!   NucleiT( (/-0.998844, 0.00000,  0.000000 /),  1.00000,  0.140,  1.00),     &
!   NucleiT( (/ 0.998844, 0.00000,  0.000000 /),  1.00000,  0.140,  1.00)      &
! 45 degree from Z
!   NucleiT( (/-0.706289, 0.00000, -0.706289 /),  1.00000,  0.140,  1.00),     &
!   NucleiT( (/ 0.706289, 0.00000,  0.706289 /),  1.00000,  0.140,  1.00)      &
    /)
!
!  Plotting and checkpointing options
!
    integer(ik), parameter :: visible_box        =   1       ! Which of the simulation boxes to plot?
    integer(ik), parameter :: plot_each          =  10       ! Reduce plotting frequency by
    integer(ik), parameter :: checkpoint_each    = 5000      ! Checkpoint after this many steps
    integer(ik)            :: checkpoint_restart = 0         ! Cycle to restart at. Must match the
                                                             ! checkpoint index -exactly-
    integer(ik), parameter :: chkptIO            = 67        ! I/O unit for checkpoint
    
!
!  Time propagation parameters
!
    integer(ik), parameter :: timesteps          = 100000    ! Number of time steps
    real(rk), parameter    :: dt                 = 0.0020_rk ! Time step, in atomic units
    real(rk), parameter    :: waitingTime        = 0.0000_rk ! Time before the real action goes, the el. field is off 
    integer(ik)            :: tstep                          ! Current time step
    real(rk)               :: time                           ! Current time
    real(rk)               :: totalNorm          = 1.0_rk    ! Current normalization
    real(rk)               :: groundGauss        = 7.0_rk    ! Ground-state separation half-width at
                                                             ! half-height
    logical, parameter     :: nucleardyn         = .false.   ! Do we want to do nuclear dynamics?
!
!  Choose the solver and number of roots (only for davidson)
!
    character(len=20), parameter :: solver       = 'davidson' ! Ground state solver - 'lanczos' or 'davidson'
!   integer(ik), parameter       :: roots_count  = 20         ! is 20 enough for a decent projector?
    integer(ik), parameter       :: roots_count  = 1          ! 
!
!  Quantum particle mass, in electron mass units
!
    real(rk)               :: mass               = 1.0_rk    ! This is just the electron!
!
!  Pulse definition
!
    real(rk), parameter    :: sqrt1_2 = 0.7071067811865_rk
    real(rk), parameter    :: pi_2    = 1.5707963267948_rk
!
!  Wavelength is usually given in nanometers. Conversion coefficient between nanometers and photon
!  energies in Hartrees is:
!
!     E[H] = 190.6/lambda[nm]
!
!  For (circular) frequencies in inverse atomic time units, divide by extra 2 pi
!
!  Field strength is usually given in watts per square meter. Conversion factor to
!
!  "Param" (pulse width) is given with respect to pulse frequency, so that the pulse
!   keeps the same shape for any frequency.
!
!                                    Name     Delay     Freq     Phase    Strength   Param         Direction
!
!   type(PulseT) :: aPulse = PulseT('Const',   0._rk, .04000_rk, 0._rk, .0500_rk, .000000_rk, (/ .000_rk, sqrt1_2, sqrt1_2  /) )  ! X
!   type(PulseT) :: aPulse = PulseT('Gauss', 200._rk, .10000_rk, 0._rk, .0500_rk, .016000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! A
!   type(PulseT) :: aPulse = PulseT('Gauss', 380._rk, .05000_rk, 0._rk, .0500_rk, .016000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! B
!   type(PulseT) :: aPulse = PulseT('Gauss', 600._rk, .05000_rk, 0._rk, .0500_rk, .005000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! C
!   type(PulseT) :: aPulse = PulseT('Gauss', 800._rk, .05695_rk, 0._rk, .1310_rk, .002501_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D1 (B1)
!   type(PulseT) :: aPulse = PulseT('Gauss', 400._rk, .05695_rk, 0._rk, .0924_rk, .010000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D2 (B1)
!   type(PulseT) :: aPulse = PulseT('Gauss', 400._rk, .03468_rk, 0._rk, .0924_rk, .002698_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D3 (B1)
!   type(PulseT) :: aPulse = PulseT('Gauss', 800._rk, .05695_rk, 0._rk, .1310_rk, .002501_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D4 (B2)
!   type(PulseT) :: aPulse = PulseT('Gauss', 400._rk, .05695_rk, 0._rk, .1310_rk, .010000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D7
!   type(PulseT) :: aPulse = PulseT('Const',1000._rk, .05695_rk, 0._rk, .0000_rk, .000000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! D10/D11/D12
!   type(PulseT) :: aPulse = PulseT('Gauss',123.6_rk, .05695_rk, 0._rk, .1310_rk, .0039_rk,   (/ .000_rk,    .000_rk, 1.000_rk /) )  ! F01  8 fs 800 nm
!   type(PulseT) :: aPulse = PulseT('Gauss',132.0_rk, .05695_rk, 0._rk, .1848_rk, .0100_rk,   (/ .000_rk,    .000_rk, 1.000_rk /) )  ! GD01  5 fs 800 nm I = 1.2x10^15 W/cm^2  
!   type(PulseT) :: aPulse = PulseT('Gauss',132.0_rk, .05695_rk, 0._rk, .1310_rk, .0100_rk,   (/ .000_rk,    .000_rk, 1.000_rk /) )  ! GD01  5 fs 800 nm I = 6.0x10^14 W/cm^2  
    type(PulseT) :: aPulse = PulseT('Const',   0._rk, .05695_rk, 0._rk, .1848_rk, .000000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! HD01
!   type(PulseT) :: aPulse = PulseT('Const',   0._rk, .05695_rk, 0._rk, .0000_rk, .000000_rk, (/ .000_rk, .000_rk, 1.000_rk /) )  ! HD01
!
!   Pulse A does not cause ionization
!   Pulse B does not cause ionization, either
!
! If you add anything here, change MainCheckpoint as well
!
    real(rk) :: norm    (  timesteps+1)     ! Wavefunction norm at each time step
    real(rk) :: etot    (  timesteps+1)     ! Total energy at each time step
    real(rk) :: field   (  timesteps+1)     ! Electric field strength
    real(rk) :: mult    (3,timesteps+1)     ! Multipole moments of the electron density
    real(rk) :: wallLeak(  timesteps+1)     ! Amount of charge leaked through the wall
                                            ! on this step
    real(rk) :: mult_nuc(3,timesteps+1)     ! Multipole moments of the nuclei. There is 
                                            ! no need to checkpoint them - they are
                                            ! trivial to calculate
    !
    !  Dynamics information. Only forces need to be included in the
    !  checkpoint - everything else can be recalculated at no cost.
    !
    real(rk) :: forces     (3,size(nucxyzqw),  timesteps+1) ! Forces on the nuclei
    real(rk) :: coordinates(3,size(nucxyzqw),0:timesteps+1) ! Coordinates of the nuclei
    real(rk) :: velocities (3,size(nucxyzqw),0:timesteps+1) ! Velocities of the nuclei
!
    contains

    subroutine h
      character(len=200) :: buf
      integer(ik)        :: inuc, axis

      call TimerStart('prog-h')
      call accuracyInitialize
      call buildGrid
      call addFields
      call FLsetNuclei(nucxyzqw)
!     call SMsetSymmetry('None')
      call SMsetSymmetry('Diatomic',main_axis_order=3, &
                         main_axis=(/1._rk,0._rk,0._rk/), &
                         inversion_centre=(/0._rk,0._rk,0._rk/))
!
!     Set initial values for the nucleus coordinates and velocity  
!
      if (nucleardyn) then
        forall (inuc=1:size(nucxyzqw))
          coordinates(:,inuc,0) = nucxyzqw(inuc)%xyz
          velocities (:,inuc,0) = nucxyzqw(inuc)%velocity
        end forall
        coordinates(:,:,1) = coordinates(:,:,0)
        velocities (:,:,1) = velocities (:,:,0)
      end if
!
!     Adjust the Pulse delay taking into account the waiting Time 
!
      aPulse%delay = aPulse%delay+waitingTime

!
!    We start from the stationary solution, provided that we don't restart
!
      if (checkpoint_restart==0) then
!
!    Choose the initial wavefunction here. At the moment,
!    there are two choices:
!
!      1. Solve stationary problem in the absence of the field
!
!        call solveStationaryProblem
!
!      2. Start with a Gaussian wavepacket
!
        call FieldInit (psi1,FLpsiGaussPacket)
        norm(1) = FieldNorm(psi1)
        write(out,"('norm of the Gauss Packet = ',f16.8)") norm(1)
        call QMNormalize(psi1,totalNorm,norm(1))
!
!    End of the initial guess
!
      else
        write(buf,"('checkpoint_main.',i6.6)") checkpoint_restart
        call MainCheckpoint ('RESTORE',trim(buf),checkpoint_restart)
      end if
!
!    Time propagation. We do time propagation with simple Verlet-like
!    algorithm - wavefunction at time t is used to propagate wavefunction
!    at time t-dt, to get wavefunction at t+dt.
!
!    1. Initialize: psi1 gets wavefunction at t=0; psi2 gets wavefunction
!                   at t=-dt (which is exactly the same)
      time = 0
!
!     call FieldSetOuterWall('ADSORBING') 
      call FieldSetOuterWall('REFLECTING')
      if (checkpoint_restart==0) then
        call FieldCopy(psi1,psi2)
        norm(1) = FieldNorm(psi1)
        call FieldNormMultipoles(psi1,mult    (:,1))
        call FLnuclearMultipoles(     mult_nuc(:,1))
      end if
!
      write (out,"('# ',a5,1x,a9,1x,a12,1x,a20,1x,a20,3(1x,a6))") &
             ' TS ', ' Time ', ' Efield ', ' W.F. norm ', ' Etot ', ' DipX ', ' DipY ', ' DipZ '
!
      call TimerStart('Time evolution')
      time_loop: do tstep=1,timesteps
!
!    2. Calculate the potential part of the time-dependent Hamiltonian
!       Here, we'll use electric field Hamiltonian, on top of the
!       static nuclear Hamiltonian. We'll use sine wave for the field,
!       so that it starts up gradually at t=0
!
        field(tstep) = FElectric_Pulse(aPulse,time)
        if (checkpoint_restart==0) then
          call FieldInit (pot1,FLnuclearPotential,mask=psi1,grace=1)
          ! 
          !  We start the pulse after an initial delay, which should
          !  allow any residual disturbances in the w.f. to settle
          !
          if (time>waitingTime) then 
            call FLsetField(aPulse%direction*field(tstep))
            call FieldInit (Hpsi,FLelectricField,mask=psi1,grace=1)
            call FieldSum  (Hpsi,pot1)
          endif
!
!     3. Electronic time step: calculate H|psi> at t=time. Wavefunction is in 
!        psi1. Change in the wavefunction is given by -2*i dt H|psi>, which 
!        should be added to wavefunction at t=time-dt (in psi2)
!
          call QMHpsi(mass,pot1,psi1,Hpsi)
          wallLeak(tstep) = dt*FieldGetLeakage()
          etot(tstep) = QMExpectation(psi1,Hpsi)
          call FieldAXPY(cmplx(0.0_rk,-(2.0_rk)*dt,kind=rk),Hpsi,psi2)
          call QMNormalize(psi2,totalNorm,norm(tstep+1))
          call FieldNormMultipoles(psi2,mult(:,tstep+1))
!
!     4. Calculate nuclear forces, provided that we have waited out our 
!        initial "shakedown" time. The forces are calculated from the 
!        Helmann-Feynman theorem. Because our basis is fixed in space,
!        this is not an approximation. 
!
          if (nucleardyn .and. time>waitingTime) then
            !
            !  Forces
            !
            force_nuc: do inuc = 1,size(nucxyzqw)
              force_axis: do axis = 1,3
                !
                !  The two terms in the force are nuclear repulsion (cheap), and
                !  electron-nuclear attraction (expensive)
                !
                call FieldSetNuclearDimension(inuc,axis)
                forces(axis,inuc,tstep) = FLpotentialRepulsionForces() + &
                                          FieldScalarIntegrate(psi1,FLpotentialAttractionForces)
              end do force_axis
            end do force_nuc
          else
            forces(:,:,tstep) = 0
          end if ! nucleardyn .and. time>waitingTime
          !
        end if ! checkpoint_restart == 0
!
!     5. Propagate nuclei, using velocity Verlet algorithm. There is no
!        need to checkpoint this part - once forces are known, the cost
!        is trivial.
!
        if (nucleardyn) then
          !
          coordinates(:,:,tstep+1) = coordinates(:,:,tstep-1) + &
                                 2*dt*velocities(:,:,tstep) 
          velocities (:,:,tstep+1) = velocities (:,:,tstep-1) + &
                                 2*dt*forces    (:,:,tstep) / spread(nucxyzqw(:)%mass,1,3)
          !
          print_dyn: do inuc=1,size(nucxyzqw)
            write(out,"('d',1x,f12.4,3(1x,f12.4),3(1x,f12.4),3(1x,f12.6))") &
              time, inuc, coordinates(:,inuc,tstep), velocities(:,inuc,tstep), &
                          forces     (:,inuc,tstep)
          end do print_dyn
          !
          !  Transfer nuclear coordinates to the fields module, for the
          !  next step.
          !
          forall (inuc=1:size(nucxyzqw))
            nucxyzqw(inuc)%xyz      = coordinates(:,inuc,tstep+1)
            nucxyzqw(inuc)%velocity = velocities (:,inuc,tstep+1)
          end forall
          call FLsetNuclei(nucxyzqw)
        end if ! nucleardyn
!
!        *** Important *** Nuclear multipoles should be calculated only 
!        *** Important *** after we moved the nuclei - or we'll have
!        *** Important *** artefactual higher harmonics generation!
!
        call FLnuclearMultipoles(mult_nuc(:,tstep+1))
!
!     5. Do some reporting and checkpointing
!
        write (buf,"('* ',i7,1x,f9.3,1x,f12.5,1x,f20.9,1x,f20.9,3(1x,f7.3))") &
               tstep, time, field(tstep), norm(tstep), etot(tstep), &
               mult(1:3,tstep)+mult_nuc(1:3,tstep)
        write (out,"(a)") trim(buf)
        if (checkpoint_restart==0) then
          !
          !  This is the normal operation sequence
          !
          if ( mod(tstep-1,plot_each)==0 ) then
            write (buf,"(f9.3,1x,f9.5,1x,f12.5,1x,f10.8,3(1x,f7.3))") &
                   time, field(tstep), etot(tstep), norm(tstep), &
                   mult(1:3,tstep) + mult_nuc(1:3,tstep)
            !
            !  Visualize is allowed to destroy Hpsi and pot1 - neither field
            !  is needed at this point.
            !
            call Visualize(psi2,buf,aPulse%direction*field(tstep))
            call TimerReport
          end if
          if ( mod(tstep-1,checkpoint_each)==0 ) then
            write(buf,"('checkpoint_field.',i6.6)") tstep
            call FieldCheckpoint('SAVE',trim(buf))
            write(buf,"('checkpoint_main.',i6.6)") tstep
            call MainCheckpoint ('SAVE',trim(buf),tstep)
            write(out,"(' Checkpoint file ',a,' written')") trim(buf)
          end if
        else
          !
          !  This is the restart. Once it is done, reset the restart indicator, and
          !  proceed with the normal run.
          !
          if (tstep==checkpoint_restart) then
            write(buf,"('checkpoint_field.',i6.6)") tstep
            call FieldCheckpoint('RESTORE',trim(buf))
            checkpoint_restart = 0
          end if
        end if
!
!     6. Now, psi1 is at t=time, while psi2 is at t=time+dt. Swap them,
!        increment time, and repeat.
!
        call swap(psi1,psi2)
        time = time + dt
        totalNorm = totalNorm - wallLeak(tstep)
        !
      end do time_loop
      !
      call TimerStop('Time evolution')
      call TimerStop('prog-h')
      call TimerReport
    end subroutine h

    subroutine swap(i,j)
      integer(ik), intent(inout) :: i, j
      integer(ik)                :: t
      t = i ; i = j ; j = t
    end subroutine swap

    subroutine solveStationaryProblem
      call TimerStart('Stationary problem')
      select case (solver)
        case default
          write (out,"('Eigenvalue solver not recognized: ',a)") trim(solver)
        case ('lanczos')
          call solveStationaryProblemLanczos
        case ('davidson')
          call solveStationaryProblemDavidson
      end select
      call Visualize(psi1,"Initial wavefunction",(/ 0.0_rk, 0.0_rk, 0.0_rk /))
      call TimerStop('Stationary problem')
    end subroutine solveStationaryProblem

    subroutine solveStationaryProblemDavidson
      real(rk) :: wf_norm
      real(rk) :: eval(roots_count)
!
!    Solve for the eigenvalues
!
      call FieldSetOuterWall('REFLECTING')
      call LUeigenvalues(mass,(/psi1,psi2,pot1,Hpsi,extra_scr/),eval, &
                         FLnuclearPotential,FLpsiGuess,FLnuclearPotential)
!
!    Get the ground-state eigenvector
!
      call LUeigenvector(1,psi1)
      wf_norm = FieldNorm(psi1)
      write (out,"(' Norm of the solution wavefunction = ',f15.10)") wf_norm
!
!    Release scratch space
!
      call LUeigenvector(0,0)
!
    end subroutine solveStationaryProblemDavidson

    subroutine solveStationaryProblemLanczos
      real(rk) :: wf_norm
      real(rk) :: eval(1)
!
!    Time-independent Hamiltonian (nuclei)
!
      call FieldInit(pot1,FLnuclearPotential)
!
!    Initial guess for the wavefunction
!
      call FieldInit(psi1,FLpsiGuess)
      call QMNormalize(psi1,1.0_rk,wf_norm)
      write (out,"(' Guess wavefunction norm = ',f15.8)") wf_norm
!
      call FieldSetOuterWall('REFLECTING')
!
!    Find a rougth solution first, using our initial guess
!
      call LCeigenvalues(mass,pot1,psi1,psi2,Hpsi,eval)
      write (out,"(' Rougth eigenvalues = ',5(1x,f15.8))") eval
      call LCeigenvector(psi1,psi2)
      wf_norm = FieldNorm(psi1)
      write (out,"(' Norm of the rougth solution wavefunction = ',f15.10)") wf_norm
!
!    Find a refined solution - this works around numerical instabilities
!    in our Davidson/Lanczos solver
!
      call LCeigenvalues(mass,pot1,psi1,psi2,Hpsi,eval)
      write (out,"(' Refined eigenvalues = ',5(1x,f15.8))") eval
      call LCeigenvector(psi1,psi2)
      wf_norm = FieldNorm(psi1)
      write (out,"(' Norm of the refined solution wavefunction = ',f15.10)") wf_norm
    end subroutine solveStationaryProblemLanczos

    subroutine Visualize(psi,title,efield)
      integer(ik), intent(in)      :: psi       ! Wavefunction; FFT will always use Hpsi
      character(len=*), intent(in) :: title     ! Descriptive label
      real(rk), intent(in)         :: efield(3) ! Driving field
      real(rk)                     :: gnorm
      !
      call TimerStart('Visualization')
      !
      !  Try to project out the contribution due to the bound states.
      !  This is rather difficult to define cleanly; therefore, we
      !  simply assign wavefunction component spatially close to the
      !  nucleus to bound state, and assume that everything else is
      !  in the continuum.
      !
      call FLsetGauss (groundGauss)
      call FieldInit  (pot1,FLgauss,mask=psi)
      call FieldZero  (Hpsi)
      call FieldMulAdd(psi,pot1,Hpsi)
      gnorm = FieldNorm(Hpsi)
      write (out,"(' Total tensity inside the central region = ',f12.9)") gnorm
      !
      !  FFT needs properly reconciled grids - let's (temporarily mark Hpsi 
      !  as a wavefunction field.
      !
      call FieldSetWavefunction(Hpsi,.true.)
      call FieldFFT(Hpsi,pot1)
      call FieldSetWavefunction(Hpsi,.false.)
!     !
!     !  Visualize density in the real and momentum space
!     !
      call FieldVisualize(0,psi,trim(title),visible_box,efield)
      call FieldFFT(psi,Hpsi)
      call FieldVisualize(1,Hpsi,trim(title),1,efield)
      call FieldShow
!     !
!     !  Whole space analysis
!     !
      write (out,"(' 1D Momentum projection - whole space ')")
      call AnalyzeSlice1D('MOMENTUM',Hpsi,1)
      write (out,"(' 1D Coordinate projection - whole space ')")
      call AnalyzeSlice1D('COORDINATE',psi,1)
!     write (out,"(' Spherical momentum distribution - whole space ')")
!     call AnalyzePolar('MOMENTUM',Hpsi,1,Nphi=64,Ntheta=32)
      !
      !  Whole space minus central region
      !
      call FieldAXPY((-1.0_rk,0.0_rk),pot1,Hpsi)
      write (out,"(' 1D Momentum projection - outer region ')")
      call AnalyzeSlice1D('MOMENTUM',Hpsi,1)
!     write (out,"(' Spherical momentum distribution - outer region ')")
!     call AnalyzePolar('MOMENTUM',Hpsi,1,Nphi=64,Ntheta=32)
      !
      call TimerStop('Visualization')
    end subroutine Visualize

    subroutine MainCheckpoint(action,name,tstep)
      character(len=*), intent(in)  :: action    ! 'SAVE' or 'RESTORE'
      character(len=*), intent(in)  :: name      ! File for save/restore
      integer(ik), intent(in)       :: tstep     ! tstep to restore at

      select case (action)
        case default
          stop 'MainCheckpoint - bad action'
        case ('SAVE')
          open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=name)
          write(chkptIO) norm    (    1:tstep+1)
          write(chkptIO) etot    (    1:tstep+1)
          write(chkptIO) field   (    1:tstep+1)
          write(chkptIO) mult    (  :,1:tstep+1)
          write(chkptIO) wallLeak(    1:tstep+1)
          if (nucleardyn) &
          write(chkptIO) forces  (:,:,1:tstep+1)
          close(chkptIO,status='keep')
        case ('RESTORE')
          open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=name)
          read(chkptIO) norm    (    1:tstep+1)
          read(chkptIO) etot    (    1:tstep+1)
          read(chkptIO) field   (    1:tstep+1)
          read(chkptIO) mult    (  :,1:tstep+1)
          read(chkptIO) wallLeak(    1:tstep+1)
          if (nucleardyn) &
          read(chkptIO) forces  (:,:,1:tstep+1)
          close(chkptIO,status='keep')
      end select
    end subroutine MainCheckpoint

    subroutine buildGrid
      real(rk) :: box(2,3)

!
      call TimerStart('Grid initialization')
      call MultiGridInit(max_grids=6,max_fields=4+n_scratch,nborder=1)
!
!  Box 0 (B0)
!
      box(:,1) = (/  -20.00d0, 20.00d0 /)
      box(:,2) = (/  -20.00d0, 20.00d0 /)
      box(:,3) = (/  -40.00d0, 40.00d0 /)
      call SimpleGridNew('The Atom (20)',   (/  45, 45, 90 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 30 /), box)
!
!
!  Box 1 (B1)
!
!     box(:,1) = (/  -20.00d0, 20.00d0 /)
!     box(:,2) = (/  -20.00d0, 20.00d0 /)
!     box(:,3) = (/  -40.00d0, 40.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/  30, 30, 60 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 30 /), box)
!
!  Box 2 (B2)
!
!     box(:,1) = (/  -30.00d0, 30.00d0 /)
!     box(:,2) = (/  -30.00d0, 30.00d0 /)
!     box(:,3) = (/  -60.00d0, 60.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/  45, 45, 90 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 30 /), box)
!
!  Box 3 (B3)
!
!     box(:,1) = (/  -45.00d0, 45.00d0 /)
!     box(:,2) = (/  -45.00d0, 45.00d0 /)
!     box(:,3) = (/  -80.00d0, 80.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/  68, 68,120 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 30 /), box)
!
!  Box 4 (B4)
!
!     box(:,1) = (/  -70.00d0, 70.00d0 /)
!     box(:,2) = (/  -70.00d0, 70.00d0 /)
!     box(:,3) = (/ -120.00d0,120.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/ 140,140,240 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 30 /), box)
!
!  Box G1 - aiming at accurate results for the ground state, limited
!           ionization.
!
!     box(:,1) = (/  -80.00d0, 80.00d0 /)
!     box(:,2) = (/  -80.00d0, 80.00d0 /)
!     box(:,3) = (/ -100.00d0,100.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/  80, 80,300 /), box)
!
!     box(:,1) = (/   -6.00d0,  6.00d0 /)
!     box(:,2) = (/   -6.00d0,  6.00d0 /)
!     box(:,3) = (/   -8.00d0,  8.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  30, 30, 60 /), box)

!  Box H1 - aiming at accurate results for the ground state, limited
!           ionization.
!
!     box(:,1) = (/  -40.00d0, 40.00d0 /)
!     box(:,2) = (/  -40.00d0, 40.00d0 /)
!     box(:,3) = (/  -40.00d0, 40.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/ 200,200,200 /), box)
 
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  40, 40, 40 /), box) 


!  Box XXX- aiming tom check difraction on the grid 
!
!     box(:,1) = (/  -60.00d0, 60.00d0 /)
!     box(:,2) = (/  -60.00d0, 60.00d0 /)
!     box(:,3) = (/  -60.00d0, 60.00d0 /)
!     call SimpleGridNew('The Atom (20)',   (/ 120,120,120 /), box)
!
!     box(:,1) = (/   -5.00d0,  5.00d0 /)
!     box(:,2) = (/   -5.00d0,  5.00d0 /)
!     box(:,3) = (/   -5.00d0,  5.00d0 /)
!     call SimpleGridNew('The Atom (5)',    (/  40, 40, 40 /), box) 

      call TimerStop('Grid initialization')
    end subroutine buildGrid

    subroutine addFields
      integer(ik)       :: i
      character(len=20) :: tag
      !
      call FieldNew('Nuclear Hamiltonian', pot1, scratch=.true. , wavefunction=.false.)
      call FieldNew(' |psi-1>',            psi1, scratch=.false., wavefunction=.true.)
      call FieldNew(' |psi-2>',            psi2, scratch=.false., wavefunction=.true.)
      call FieldNew('H|psi>',              Hpsi, scratch=.true. , wavefunction=.false.)
      !
      do i=1,n_scratch
        write (tag,"('W.F. scratch ',i5)") i
        call FieldNew(trim(tag), extra_scr(i), scratch=.true., wavefunction=.false.)
      end do
    end subroutine addFields

  end module h_module
!
  subroutine driver
    use h_module

    call h
  end subroutine driver
