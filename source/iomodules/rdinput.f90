  module rdinput

  contains

!#######################################################################

    subroutine read_input
      
      use parameters
      use parsemod
      use iomod
      use channels
      
      implicit none
      
      integer            :: i,k,n,l
      character(len=120) :: atmp1,atmp2
      logical            :: iscvs,energyonly,ldiag,llanc,ldiffsection,&
                            llasersection
           
!-----------------------------------------------------------------------
! Set 'traps'
!-----------------------------------------------------------------------
      energyonly=.false.
      ldiag=.false.
      llanc=.false.
      iscvs=.false.
      ldiffsection=.false.
      lmatvec=.false.
      llasersection=.false.
      llci=.false.

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(iin)
        
      i=0
      if (keyword(1).ne.'end-input') then
10       continue
         i=i+1
         
         if (keyword(i).eq.'method') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'adc1') then
                  method=1
               else if (keyword(i).eq.'cis') then
                  method=1
                  lcis=.true.
               else if (keyword(i).eq.'adc2' &
                    .or.keyword(i).eq.'adc2-s' &
                    .or.keyword(i).eq.'adc2s') then
                  method=2
               else if (keyword(i).eq.'adc2-x' &
                    .or.keyword(i).eq.'adc2x') then
                  method=3
               else if (keyword(i).eq.'adc1-x' &
                    .or.keyword(i).eq.'adc1x') then
                  method=4
               endif
            else
               goto 100
            endif
            ! Final space method, if present
            if (keyword(i+1).eq.',') then
               i=i+2
               if (keyword(i).eq.'adc2'.or.keyword(i).eq.'adc2-s') then
                  method_f=2
               else if (keyword(i).eq.'adc2-x') then
                  method_f=3
               else
                  goto 100
               endif
            else
               method_f=method
            endif

         else if (keyword(i).eq.'istate_symm') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nirrep
            else
               goto 100
            endif

         else if (keyword(i).eq.'dipole_symm') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nirrep2
            else
               goto 100
            endif

         else if (keyword(i).eq.'initial_state') then

            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'target') then
                  ! Initial state determined by comparison to a
                  ! user-supplied target state.
                  !
                  ! Assign any integer i>0: this is 
                  ! necessary to get around checks on the 
                  ! statenumber
                  statenumber=999
               else if (keyword(i).eq.'mix') then
                  ! Mixed initial state
                  lmixistate=.true.
                  ! Set statenumber to something positive to
                  ! get around checks later on
                  statenumber=999
                  ! Read the state indices and coefficientes
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) nmix
                     allocate(imix(nmix),cmix(nmix))
                     do k=1,nmix
                        i=i+2
                        read(keyword(i),*) imix(k)
                     enddo
                     do k=1,nmix
                        i=i+2
                        read(keyword(i),*) cmix(k)
                     enddo
                     ! Normalisation of the coefficient vector
                     cmix=cmix/sqrt(dot_product(cmix,cmix))
                  else
                     errmsg='The initial mixed state numbers and &
                          coefficients have not been given'
                     call error_control
                  endif
               else
                  ! ADC eigenstate number
                  read(keyword(i),*) statenumber
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'dipole_component') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') tranmom
               tranmom2=tranmom
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'cvs') then
            iscvs=.true.
            
         else if (keyword(i).eq.'rixs') then
            lrixs=.true.

         else if (keyword(i).eq.'tpa') then
            ltpa=.true.

         else if (keyword(i).eq.'matvec') then
            lmatvec=.true.

         else if (keyword(i).eq.'energy_only') then
            energyonly=.true.

         else if (keyword(i).eq.'fakeip') then
            lfakeip=.true.
            energyonly=.true.

         else if (keyword(i).eq.'motype') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') motype
            else
               goto 100
            endif

         else if (keyword(i).eq.'lci') then
            llci=.true.

         else if (keyword(i).eq.'init_energy') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') init_energy
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'diag_section') then
            ldiag=.true.
15          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-diag_section') goto 15
            i=inkw

         else if (keyword(i).eq.'lanczos_section') then
            llanc=.true.
20          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-lanczos_section') goto 20
            i=inkw

         else if (keyword(i).eq.'diag_final_section') then
            ldiagfinal=.true.
25          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-diag_final_section') goto 25
            i=inkw            
            
         else if (keyword(i).eq.'no_tdm') then
            ltdm_gs2i=.false.

         else if (keyword(i).eq.'istate_frozen_core') then
            lifrzcore=.true.

         else if (keyword(i).eq.'fstate_frozen_core') then
            lffrzcore=.true.

         else if (keyword(i).eq.'frozen_core') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               atmp1=keyword(i)               
               if (keyword(i+1).eq.',') then
                  i=i+2
                  atmp2=keyword(i)
               endif
               if (atmp1.eq.'initial'.or.atmp2.eq.'initial') &
                    lifrzcore=.true.
               if (atmp1.eq.'final'.or.atmp2.eq.'final') &
                    lffrzcore=.true.
            else
               goto 100
            endif

         else if (keyword(i).eq.'scf_iter') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) scfiter
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'freeze') then
            if (keyword(i+1).eq.'=') then
               i=i+1
30             continue
               i=i+1
               read(keyword(i),*) k
               iexpfrz(k)=1
               if (keyword(i+1).eq.',') then
                  i=i+1
                  goto 30
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'basis') then             
            lrungamess=.true.
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) basname
            else
               goto 100
            endif

         else if (keyword(i).eq.'diffuse') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Diffuse function type
               ! Note that difftype=4 is reserved for explicit
               ! lists of diffuse functions
               if (keyword(i).eq.'kbj'.or.keyword(i).eq.'kbj_cont') then
                  difftype=1
               else if (keyword(i).eq.'even') then
                  difftype=2
                  i=i+2
                  read(keyword(i),*) diffratio
               else if (keyword(i).eq.'kbj_ryd') then
                  difftype=3
               else
                  goto 100
               endif
               ! Numbers of diffuse functions
35             continue
               if (keyword(i+1).eq.',') then
                  i=i+2
                  call getdiffinfo(n,l,keyword(i))
                  ndiff(l)=n
                  goto 35
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'geometry') then
            lrungamess=.true.
40          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-geometry') then
               natm=natm+1
               goto 40
            endif
            ncoo=natm*3
            
         else if (keyword(i).eq.'point_group') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') pntgroup
            else
               goto 100
            endif

         else if (keyword(i).eq.'denord') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) denord
            else
               goto 100
            endif

         else if (keyword(i).eq.'mem') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Maximum memory value
               read(keyword(i),*) maxmem
               ! Units
               if (keyword(i+1).eq.',') then
                  i=i+2
                  if (keyword(i).eq.'mb') then
                     ! Do nothing: this is the default
                  else if (keyword(i).eq.'gb') then
                     maxmem=maxmem*1024.0d0
                  else
                     errmsg='Unknown memory unit: '//trim(keyword(i))
                     call error_control
                  endif
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'dyson_section') then
            ldyson=.true.
            lfakeip=.true.
45          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-dyson_section') goto 45
            i=inkw
            
         else if (keyword(i).eq.'target_section') then
            ltarg=.true.
50          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-target_section') goto 50
            i=inkw

         else if (keyword(i).eq.'dipole_moment') then
            ldipole=.true.

         else if (keyword(i).eq.'debug') then
            debug=.true.

         else if (keyword(i).eq.'autospec_section') then
            lautospec=.true.
55          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-autospec_section') goto 55
            i=inkw

         else if (keyword(i).eq.'fdstates_section') then
            lfdstates=.true.
60          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-fdstates_section') goto 60
            i=inkw

         else if (keyword(i).eq.'cap_section') then
            lcap=.true.
65          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-cap_section') goto 65
            i=inkw

         else if (keyword(i).eq.'propagation_section') then
            lpropagation=.true.
70          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-propagation_section') goto 70
            i=inkw
            
         else if (keyword(i).eq.'diffuse_section') then
            ldiffsection=.true.
            difftype=4
75          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-diffuse_section') goto 75
            i=inkw

         else if (keyword(i).eq.'nto') then
            lnto=.true.

         else if (keyword(i).eq.'laser_section') then
            llasersection=.true.
            npulse=npulse+1
80          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-laser_section') goto 80
            i=inkw
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
            STOP
         endif
      
         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
           
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         STOP
         
      endif

!-----------------------------------------------------------------------
! Set CVS flags
!-----------------------------------------------------------------------
      if (iscvs) then
         if (energyonly) then
            lcvs=.true.
         else
            lcvsfinal=.true.
         endif
      endif
      
      if (lrixs) lcvsfinal=.true.

!-----------------------------------------------------------------------
! If an energy-only calculation has been requested, reset method
! accordingly
!-----------------------------------------------------------------------
      if (energyonly) method=-method

!-----------------------------------------------------------------------
! Final space method.
! Default: same as the initial space method.
!-----------------------------------------------------------------------
      if (method_f.eq.0) method_f=method

!-----------------------------------------------------------------------
! Read the Davidson section(s)
!-----------------------------------------------------------------------
      if (ldiag) call rddiaginp

      if (ldiagfinal)  call rddiagfinalinp

!-----------------------------------------------------------------------
! Read the Lanczos section
!-----------------------------------------------------------------------
      if (llanc) call rdlancinp

!-----------------------------------------------------------------------
! Read the Dyson section
!-----------------------------------------------------------------------
      if (ldyson) call rddysoninp

!-----------------------------------------------------------------------
! Read the target state section
!-----------------------------------------------------------------------
      if (ltarg) call rdtargetinp

!-----------------------------------------------------------------------
! Read the autospec section
!-----------------------------------------------------------------------
      if (lautospec) call rdautospecinp

!-----------------------------------------------------------------------
! Read the filter diagonalisation states section
!-----------------------------------------------------------------------
      if (lfdstates) call rdfdstatesinp

!-----------------------------------------------------------------------
! Read the CAP section
!-----------------------------------------------------------------------
      if (lcap) call rdcapinp

!-----------------------------------------------------------------------
! Read the propagation section
!-----------------------------------------------------------------------
      if (lpropagation) call rdpropagationinp

!-----------------------------------------------------------------------
! Read the laser sections (after allocating arrays now that the number
! of pulses in known)
!-----------------------------------------------------------------------
      if (llasersection) then

         ! Allocate arrays
         call alloclaserpar

         ! Read the laser sections
         do i=1,npulse
            call rdlaserinp(i)
         enddo

      endif
      
!-----------------------------------------------------------------------
! Read the diffuse function section
!-----------------------------------------------------------------------
      if (ldiffsection) call rddiffbaslist

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
      call checkinp(ldiag,llanc,energyonly)
      
!-----------------------------------------------------------------------
! Read the geometry section if necessary
!-----------------------------------------------------------------------
      if (lrungamess) call rdgeometry

      return

    end subroutine read_input

!#######################################################################

    subroutine checkinp(ldiag,llanc,energyonly)

      use parameters
      use iomod

      implicit none

      integer            :: i
      character(len=170) :: msg
      character(len=2)   :: ai
      logical            :: ldiag,llanc,energyonly

!-----------------------------------------------------------------------
! ADC level
!-----------------------------------------------------------------------
      if (method.eq.0) then
         msg='The method has not been been given'
         goto 999
      endif

!-----------------------------------------------------------------------
! Density order
!-----------------------------------------------------------------------
      if (denord.lt.2.or.denord.gt.3) then
         msg='Only values of 2 or 3 are allowed for denord'
         goto 999
      endif

!-----------------------------------------------------------------------
! Initial state symmetry
!-----------------------------------------------------------------------
      if (nirrep.eq.0) then
         msg='The initial state symmetry has not been given'
         goto 999
      endif

!-----------------------------------------------------------------------
! Dipole operator symmetry
!-----------------------------------------------------------------------
      if (nirrep2.eq.0.and..not.ldyson) then
         msg='The dipole operator symmetry has not been given'
         goto 999
      endif

!-----------------------------------------------------------------------
! Initial state number
!-----------------------------------------------------------------------
      if (statenumber.eq.-1.and..not.energyonly) then
         msg='The initial state number has not been given'
         goto 999
      endif

      if (energyonly) statenumber=0

      ! Mixed initial states: only supported for wavepacket
      ! propagations
      if (lmixistate.and..not.lpropagation) then
         msg='Mixed initial states are only supported in TD-ADC and &
              TD-CIS calculations'
         goto 999
      endif
      
!-----------------------------------------------------------------------
! Dipole operator component
!-----------------------------------------------------------------------
      if (tranmom2.eq.''.and..not.ldyson) then
         msg='The dipole operator component has not been given'
         goto 999
      endif

!-----------------------------------------------------------------------
! MO storage
!-----------------------------------------------------------------------
      if (motype.ne.'incore'.and.motype.ne.'disk') then
         msg='Unknown MO storage flag - '//trim(motype)
         goto 999
      endif

!-----------------------------------------------------------------------
! Initial space diagonalisation section: only required if either:
!
! (1) We are condisdering excitation from an excited state, or;
! (2) We are performing an energy-only calculation.
!-----------------------------------------------------------------------
      if (.not.lautospec.and.abs(method).ne.1) then

         if (.not.ldiag) then
            if (statenumber.gt.0.or.energyonly) then
               msg='The diagonalisation section has not been found.'
               goto 999
            endif
         endif

         if (statenumber.gt.0.or.energyonly) then
            if (davstates.eq.0) then
               msg='The number of initial space states has not been given'
               goto 999
            endif
         endif

      endif
         
!-----------------------------------------------------------------------
! Final space diagonalisation section
!-----------------------------------------------------------------------
      if (ldiagfinal) then
         if (davstates_f.eq.0) then
            msg='The number of final space states has not been given'
            goto 999
         endif
      endif

!-----------------------------------------------------------------------
! Lanczos section: only required if we are considering ionization,
! a block-Lanczos-RIXS or block-Lanczos-TPA calculation
!-----------------------------------------------------------------------
      ! Photoionisation cross-section calculation
      if (.not.energyonly &
           .and..not.ldiagfinal &
           .and..not.ldyson &
           .and..not.lrixs &
           .and..not.ltpa &
           .and..not.lautospec &
           .and..not.lfdstates &
           .and..not.lpropagation &
           .and..not.lcap) then

         if (.not.llanc) then
            msg='No Lanczos section has been found'
            goto 999
         endif
         
         if (lmain.eq.0.and..not.ldynblock) then
            msg='The Lanczos block size has not been given'
            goto 999
         endif

         if (ncycles.eq.0) then
            msg='The no. Lanczos iterations has not been given'
         endif

      endif

      ! Block Lanczos-RIXS calculation
      if (lrixs.and..not.ldiagfinal) then

         if (.not.llanc) then
            msg='No Lanczos section has been found'
            goto 999
         endif
      
         if (ncycles.eq.0) then
            msg='The no. Lanczos iterations has not been given'
         endif
      
         ! Set lancguess=5: this is required so that the appropriate
         ! initial vectors are read from fie
         lancguess=5
         
      endif

      ! Block Lanczos-TPA calculation
      if (ltpa) then

         if (.not.llanc) then
            msg='No Lanczos section has been found'
            goto 999
         endif

         if (ncycles.eq.0) then
            msg='The no. Lanczos iterations has not been given'
         endif

         ! Set lancguess=6: this is required so that the appropriate
         ! initial vectors are used
         lancguess=6

         ! If a TPXAS spectrum is to be calculated, then make sure
         ! that a diag_final section is present
         if (lcvsfinal.and..not.ldiagfinal) then
            msg='No diag_final section has been found'
            goto 999
         endif

      endif

!-----------------------------------------------------------------------
! GAMESS
!-----------------------------------------------------------------------
      if (lrungamess) then

         if (basname.eq.'') then
            msg='The basis set has not been given'
            goto 999
         endif
         
         if (natm.eq.0) then
            msg='The geometry has not been given'
            goto 999
         endif

      endif

!-----------------------------------------------------------------------
! Dyson section
!-----------------------------------------------------------------------
      if (ldyson) then
         if (dysirrep.eq.0) then
            msg='The symmetry of the final states in the Dyson &
                 orbital calculation has not been given'
            goto 999
         endif
         if (dysdiag.eq.2.and..not.ldiagfinal) then
            msg='Iterative diagonalisation of the IP-ADC(2) &
                 Hamiltonian has been requested, but the diag_final &
                 section is missing...'
            goto 999
         endif
      endif

!-----------------------------------------------------------------------
! Target state section
!-----------------------------------------------------------------------
      if (ltarg) then
         ! Make sure that statenumber>0 so that the initial space
         ! diagonalisation is performed
         statenumber=1
         if (detthrsh.eq.-1.0d0) then
            msg='The Slater determinant threshold has not been given'
            goto 999
         endif
         if (ovrthrsh.eq.-1.0d0) then
            msg='The initial state overlap threshold has not been given'
            goto 999
         endif
         if (detfile.eq.'') then
            msg='The target Slater determinant file name has not been &
                 given'
            goto 999
         endif
         if (mofile.eq.'') then
            msg='The target MO file has not been given'
            goto 999
         endif
      endif

!-----------------------------------------------------------------------
! Autospec section
!-----------------------------------------------------------------------
      ! Wavepacket propagation
      if (lautospec.and.autoprop.eq.1) then

         ! Final propagation time
         if (tfinal.eq.0.0d0) then
            msg='The final propagation time has not been given'
            goto 999
         endif

         ! Timestep
         if (tout.eq.0.0d0) then
            msg='The timestep has not been given'
            goto 999
         endif
         
         ! Autocorrelation function order
         if (autoord.lt.0.or.autoord.gt.2) then
            msg='Invalid autocorrelation order. Valid values are &
                 0, 1 or 2 only'
            goto 999
         endif

         ! Saving of the 1h1p parts of the wavepacket to file is not
         ! yet supported
         if (save1h1p) then
            msg='Saving of the 1h1p parts of the wavepacket to file is &
                 not yet supported'
            goto 999
         endif
         
      endif

      ! Chebyshev recursion
      if (lautospec.and.autoprop.eq.2) then

         ! Order of the Chebyshev expansion
         if (chebyord.eq.0) then
            msg='The number of terms in the Chebyshev expansion &
                 of delta(E-H) has not been given'
            goto 999
         endif

      endif
      
!-----------------------------------------------------------------------
! Filter diagonalisation states section
!-----------------------------------------------------------------------
      ! Wavepacket propagation
      if (lfdstates.and.autoprop.eq.1) then

         ! Final propagation time
         if (tfinal.eq.0.0d0) then
            msg='The final propagation time has not been given'
            goto 999
         endif

         ! Timestep
         if (tout.eq.0.0d0) then
            msg='The timestep has not been given'
            goto 999
         endif

         ! fdiag data file
         if (fdiagdat.eq.'') then
            msg='The name of the FDIAG data filehas not been given'
            goto 999
         endif
         
         ! State selection file
         if (fdiagsel.eq.'') then
            msg='The name of the FDIAG state selection file has &
                 not been given'
            goto 999
         endif

         ! Calculation of NTOs using the 1h1p parts of the wavepacket
         ! read from file is not yet supported
         if (save1h1p) then
            msg='The calculation of NTOs using the 1h1p parts of the &
                 wavepacket read from file is not yet supported'
            goto 999
         endif
         
      endif

      ! Chebyshev recursion
      if (lfdstates.and.autoprop.eq.2) then

         ! Data file
         if (fdiagdat.eq.'') then
            msg='The name of the CHEBYFD data filehas not been given'
            goto 999
         endif
         
         ! State selection file
         if (fdiagsel.eq.'') then
            msg='The name of the CHEBYFD state selection file has &
                 not been given'
            goto 999
         endif
         
      endif
         
!-----------------------------------------------------------------------
! CAP section
!-----------------------------------------------------------------------
      if (lcap) then
         
         ! CAP type
         if (icap.eq.0) then
            msg='The CAP type has not been given'
            goto 999
         endif

         ! CAP strength
         if (capstr.eq.0.0d0) then
            msg='The CAP strength has not been given'
            goto 999
         endif

         ! CAP order
         if (icap.eq.1.or.icap.eq.2.or.icap.eq.3.or.icap.eq.7) then
            if (capord.eq.-1) then
               msg='The monomial CAP order has not been set'
               goto 999
            endif
         endif

         ! Integration grid
         if (icap.eq.2 &
              .or.icap.eq.3 &
              .or.icap.eq.4 &
              .or.icap.eq.5 &
              .or.icap.eq.6) then
            if (nang(1).ne.110.and.nang(1).ne.302&
                 .and.nang(1).ne.770) then
               msg='The number of angular grid points can only be &
                    one of: 110, 302 or 770'
               goto 999
            endif
            if (nang(2).ne.110.and.nang(2).ne.302&
                 .and.nang(2).ne.770) then
               msg='The number of angular grid points can only be &
                    one of: 110, 302 or 770'
               goto 999
            endif
         endif
         
         ! Propagation section
         if (.not.lpropagation) then
            msg='The propagation sesction has not been given'
            goto 999
         endif

         ! If the CAP is to projected out of the space spanned by a
         ! set of excited states, make sure that a diag_section has
         ! been given
         if (iprojcap.eq.2.and..not.ldiag.and.method.ne.1) then
            msg='For ADC(2) caclulations and projection=all, a &
                 diag_section is required'
            goto 999
         endif

         ! Selected bound-unbound and bound-bound projection is
         ! only available for TD-ADC(1) and TD-CIS calculations
         if (iprojcap.eq.3.or.iprojcap.eq.4) then
            if (method.ne.1) then
               errmsg='Bound-unbound and bound-bound projection is &
                    not available for TD-ADC(2) calculations'
            endif
         endif
         
         ! Diagonalisation of H-iW is only currently supported at
         ! the ADC(1)/CIS level
         if (lcapdiag.and.method.ne.1) then
            msg='Diagonalisation of H-iW is not yet supported for &
                 ADC(2) calculations'
         endif
         
      endif

!-----------------------------------------------------------------------
! Propagation section
!-----------------------------------------------------------------------
      if (lpropagation) then

         ! Final propagation time
         if (tfinal.eq.0.0d0) then
            msg='The final propagation time has not been given'
            goto 999
         endif

         ! Timestep
         if (tout.eq.0.0d0) then
            msg='The timestep has not been given'
            goto 999
         endif

         ! Represenation of the wavefunction: using the Hamiltonian
         ! eigenstate basis is only supported for TD-ADC(1) and TD-CIS
         ! calculations
         if (tdrep.eq.2) then
            if (method.eq.2.or.method.eq.3) then
               msg='The eigenstate representation cannot be used with &
                    ADC(2)'
               goto 999
            endif
         endif
         
      endif

!-----------------------------------------------------------------------
! Laser sections
!-----------------------------------------------------------------------
      do i=1,npulse
         write(ai,'(i2)') i

         ! Carrier wave
         if (ipulse(i).eq.0) then
            msg='The pulse type (carrier wave) has not been given &
                 for laser_section '//trim(adjustl(ai))
            goto 999
         endif

         ! Pulse envelope type
         if (ienvelope(i).eq.0) then
            msg='The pulse envelope type has not been given &
                 for laser_section '//trim(adjustl(ai))
            goto 999
         endif
         
         ! Laser-molecule orientation
         if (sum(abs(pulse_vec(:,i))).eq.0.0d0) then
            msg='The laser pulse orientation has not been given &
                 for laser_section '//trim(adjustl(ai))
            goto 999
         endif

         ! Laser frequency
         if (freq(i).eq.0.0d0) then
            msg='The laser frequency has not been given &
                 for laser_section '//trim(adjustl(ai))
            goto 999
         endif

         ! Laser strength
         if (strength(i).eq.0.0d0) then
            msg='The laser strength has not been given &
                 for laser_section '//trim(adjustl(ai))
            goto 999
         endif
         
      enddo
      
!-----------------------------------------------------------------------
! ADC(1)-x calculations - at present, this is only supported for a
! wavepacket propagation calculation
!-----------------------------------------------------------------------
      if (method.eq.4.and..not.lpropagation) then
         msg='ADC(1)-x calculations are only currently supported for &
              wavepacket propagation calculations'
         goto 999
      endif

!-----------------------------------------------------------------------
! Natural transition orbitals: currently these are only available for
! excitation from the ground state or for wavepacket propagations
!-----------------------------------------------------------------------
      if (lnto.and.statenumber.ne.0.and..not.lpropagation) then
         msg='NTOs are currently only available for excitation from &
              the ground state'
         goto 999
      endif
      
      return

999   continue
      errmsg='Problem with the input file: '//trim(msg)
      call error_control

    end subroutine checkinp

!#######################################################################

    subroutine rddiaginp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the diagonalisation section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'diag_section') goto 1

!-----------------------------------------------------------------------
! Read the diagonalisation parameters
!-----------------------------------------------------------------------
5    call rdinp(iin)
      
      i=0

      if (keyword(1).ne.'end-diag_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'nstates') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) davstates
            else
               goto 100
            endif

         else if (keyword(i).eq.'block_size') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dmain
            else
               goto 100
            endif

         else if (keyword(i).eq.'tol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) davtol
            else
               goto 100
            endif

         else if (keyword(i).eq.'guess') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'adc1') then
                  ladc1guess=.true.
               else if (keyword(i).eq.'noise') then
                  lnoise=.true.
               else if (keyword(i).eq.'subdiag') then
                  lsubdiag=.true.
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) guessdim
                  endif
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'maxit') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxiter
            else
               goto 100
            endif

         else if (keyword(i).eq.'target') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               eigentype=2
               read(keyword(i),*) davtarg
               davtarg=davtarg/eh2ev
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'method') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'davidson') then
                  solver=1
               else if (keyword(i).eq.'relaxation') then
                  solver=2
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     if (keyword(i).eq.'xsil') then
                        integrator=1
                     else if (keyword(i).eq.'bs') then
                        integrator=2
                     else if (keyword(i).eq.'rkf45') then
                        integrator=3
                     else
                        errmsg='Unknown keyword: '//trim(keyword(i))
                        call error_control
                     endif
                  endif
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif
         
         ! Deprecated: included for backwards compatibility
         else if (keyword(i).eq.'mem') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxmem
            else
               goto 100
            endif

         else if (keyword(i).eq.'preconditioner') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'dpr') then
                  precon=1
               else if (keyword(i).eq.'olsen') then
                  precon=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'maxvec') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxsubdim
            else
               goto 100
            endif

         else if (keyword(i).eq.'deflate') then
            ldfl=.true.
            
         else if (keyword(i).eq.'nodfl') then
            ldfl=.false.

         else if (keyword(i).eq.'krydim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) kdim
            else
               goto 100
            endif

         else if (keyword(i).eq.'timestep') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) stepsize
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'siltol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) siltol
            else
               goto 100
            endif

         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif

         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control

      endif

      return

    end subroutine rddiaginp

!#######################################################################

    subroutine rddiagfinalinp

      use parameters
      use parsemod
      use iomod
      use channels
      
      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the final space diagonalisation section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'diag_final_section') goto 1

!-----------------------------------------------------------------------
! Read the final space diagonalisation parameters
!-----------------------------------------------------------------------
5    call rdinp(iin)
      
      i=0

      if (keyword(1).ne.'end-diag_final_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'nstates') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) davstates_f
            else
               goto 100
            endif

         else if (keyword(i).eq.'block_size') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dmain_f
            else
               goto 100
            endif

         else if (keyword(i).eq.'tol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) davtol_f
            else
               goto 100
            endif

         else if (keyword(i).eq.'guess') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'adc1') then
                  ladc1guess_f=.true.
               else if (keyword(i).eq.'noise') then
                  lnoise_f=.true.
               else if (keyword(i).eq.'subdiag') then
                  lsubdiag_f=.true.
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) guessdim_f
                  endif
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'maxit') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxiter_f
            else
               goto 100
            endif

         else if (keyword(i).eq.'target') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               eigentype=2
               read(keyword(i),*) davtarg
               davtarg=davtarg/eh2ev
            else
               goto 100
            endif

         else if (keyword(i).eq.'method') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'davidson') then
                  solver_f=1
               else if (keyword(i).eq.'relaxation') then
                  solver_f=2
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     if (keyword(i).eq.'xsil') then
                        integrator=1
                     else if (keyword(i).eq.'bs') then
                        integrator=2
                     else if (keyword(i).eq.'rkf45') then
                        integrator=3
                     else
                        errmsg='Unknown keyword: '//trim(keyword(i))
                        call error_control
                     endif
                  endif
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif
            
         ! Deprecated: included for backwards compatibility
         else if (keyword(i).eq.'mem') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxmem
            else
               goto 100
            endif 

         else if (keyword(i).eq.'preconditioner') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'dpr') then
                  precon_f=1
               else if (keyword(i).eq.'olsen') then
                  precon_f=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'maxvec') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxsubdim_f
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'deflate') then
            ldfl_f=.true.
                        
         else if (keyword(i).eq.'nodfl') then
            ldfl_f=.false.

         else if (keyword(i).eq.'krydim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) kdim_f
            else
               goto 100
            endif

         else if (keyword(i).eq.'timestep') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) stepsize_f
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'siltol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) siltol_f
            else
               goto 100
            endif
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
           ! then read them, else read the next line
           if (i.lt.inkw) then
              goto 10
           else
              goto 5
           endif

         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control

        endif

      return

    end subroutine rddiagfinalinp

!#######################################################################

    subroutine rdlancinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the Lanczos section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'lanczos_section') goto 1

!-----------------------------------------------------------------------
! Read the Lanczos parameters
!-----------------------------------------------------------------------
5    call rdinp(iin)
      
      i=0

      if (keyword(1).ne.'end-lanczos_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'block_size') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'dynamic') then
                  ldynblock=.true.
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) maxblock
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) tdtol
                  endif
               else
                  read(keyword(i),*) lmain
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'iter') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) ncycles
            else
               goto 100
            endif

         else if (keyword(i).eq.'guess') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'adc1') then
                  lancguess=2
               else if (keyword(i).eq.'is_mix') then
                  lancguess=3
               else if (keyword(i).eq.'adc1_mix') then
                  lancguess=4
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif               
            else
               goto 100
            endif
            ! Taking linear combinations of 1h1p and 2h2p unit vectors
            ! doesn't make much sense for ADC(2)-x, so reset lancguess
            ! accordingly if necessary
            if (method.eq.3.and.lancguess.eq.3) lancguess=1
            if (method.eq.3.and.lancguess.eq.4) lancguess=2

         ! Deprecated: included for backwards compatibility
         else if (keyword(i).eq.'mem') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) maxmem
            else
               goto 100
            endif

         else if (keyword(i).eq.'ortho') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'pro') then
                  orthotype=1
               else if (keyword(i).eq.'mpro') then
                  orthotype=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control

        endif

      return

    end subroutine rdlancinp

!#######################################################################

    subroutine rddysoninp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the Dyson section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'dyson_section') goto 1

!-----------------------------------------------------------------------
! Read the Dyson orbital parameters
!-----------------------------------------------------------------------
5    call rdinp(iin)
      
      i=0

      if (keyword(1).ne.'end-dyson_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'symm') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dysirrep
            else
               goto 100
            endif

         else if (keyword(i).eq.'elim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dyslim
               if (keyword(i+1).eq.',') then
                  i=i+2
                  if (keyword(i).eq.'ev') then
                     dyslim=dyslim/eh2ev
                  else
                    errmsg='Unknown unit'//trim(keyword(i))
                    call error_control
                  endif
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'diag') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'full') then
                  dysdiag=1
               else if (keyword(i).eq.'iter') then
                  dysdiag=2
               else
                  errmsg='Unknown diagonalisation type:'//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif
               
         else if (keyword(i).eq.'out') then
            if (keyword(i+1).eq.'=') then
               i=i+2
15             continue
               if (keyword(i).eq.'molden') then
                  dysout(1)=1
               else if (keyword(i).eq.'ezdyson') then
                  dysout(2)=1
               else
                  errmsg='Unknown Dyson orbital output type: '&
                       //trim(keyword(i))
                  call error_control
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  goto 15
               endif               
            else
               goto 100
            endif

         else if (keyword(i).eq.'pe_wf') then
            if (keyword(i+1).eq.'=') then
               i=i+2              
               read(keyword(i),*) lmax
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) zcore
               else
                  errmsg='The core charge has not been given with &
                       the pe_wf keyword'
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'pe_en') then
            if (keyword(i+1).eq.'=') then
               i=i+2              
               read(keyword(i),*) nelen
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) eleni
                  eleni=max(eleni,0.1d0)
               else
                  errmsg='The lower photoelectron K.E. has not been &
                       given with the pe_en keyword'
                  call error_control
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) elenf
               else
                  errmsg='The upper photoelectron K.E. has not been &
                       given with the pe_en keyword'
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'grid') then
            if (keyword(i+1).eq.'=') then
               i=i+2              
               read(keyword(i),*) ngrdpnts
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) grdi
               else
                  errmsg='The lower bound of the grid has not been &
                       given with the grid keyword'
                  call error_control
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) grdf
               else
                  errmsg='The upper bound of the grid has not been &
                       given with the grid keyword'
                  call error_control
               endif
            else
               goto 100
            endif

         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control

      endif
            
      return
      
    end subroutine rddysoninp

!#######################################################################

    subroutine rdtargetinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the target state section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'target_section') goto 1

!-----------------------------------------------------------------------
! Read the target state matching parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-target_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'detthrsh') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) detthrsh
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'ovrthrsh') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) ovrthrsh
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'detfile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') detfile
            else
               goto 100
            endif

         else if (keyword(i).eq.'mofile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') mofile
            else
               goto 100
            endif
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif
         
      return

    end subroutine rdtargetinp

!#######################################################################

    subroutine rdautospecinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the autospec section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'autospec_section') goto 1

!-----------------------------------------------------------------------
! Read the autospec parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-autospec_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'tfinal') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tfinal
            else
               goto 100
            endif

         else if (keyword(i).eq.'tout') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tout
            else
               goto 100
            endif

         else if (keyword(i).eq.'tol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) autotol
            else
               goto 100
            endif

         else if (keyword(i).eq.'krydim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) kdim
            else
               goto 100
            endif

         else if (keyword(i).eq.'order') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) autoord
            else
               goto 100
            endif

         else if (keyword(i).eq.'method') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'propagation') then
                  autoprop=1
               else if (keyword(i).eq.'chebyshev') then
                  autoprop=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'nterms') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) chebyord
            else
               goto 100
            endif

         else if (keyword(i).eq.'filter') then
            lprojpsi0=.true.
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) projen
            else
               goto 100
            endif
            ! Optional energy unit
            if (keyword(i+1).eq.',') then
               i=i+2
               call convert_energy(keyword(i),projen)
            endif

         else if (keyword(i).eq.'save1h1p') then
            save1h1p=.true.
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif
         
      return
      
    end subroutine rdautospecinp

!#######################################################################

    subroutine rdfdstatesinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Read to the fdstates section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'fdstates_section') goto 1

!-----------------------------------------------------------------------
! Read the fdstates parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-fdstates_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'tfinal') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tfinal
            else
               goto 100
            endif

         else if (keyword(i).eq.'tout') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tout
            else
               goto 100
            endif         

         else if (keyword(i).eq.'krydim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) kdim
            else
               goto 100
            endif

         else if (keyword(i).eq.'data') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               fdiagdat=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'selection') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               fdiagsel=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'method') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'propagation') then
                  autoprop=1
               else if (keyword(i).eq.'chebyshev') then
                  autoprop=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'nterms') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) chebyord
            else
               goto 100
            endif

         else if (keyword(i).eq.'read1h1p') then
            read1h1p=.true.
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif

      return
      
    end subroutine rdfdstatesinp

!#######################################################################

    subroutine rdcapinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer          :: i,n
      character(len=2) :: ai
      
!-----------------------------------------------------------------------
! Read to the CAP section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'cap_section') goto 1

!-----------------------------------------------------------------------
! Read to the CAP parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-cap_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'cap_type') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'monomial') then
                  icap=1
               else if (keyword(i).eq.'monomial_grid') then
                  icap=2
               else if (keyword(i).eq.'atom_monomial') then
                  icap=3
               else if (keyword(i).eq.'moiseyev') then
                  icap=4
               else if (keyword(i).eq.'atom_moiseyev') then
                  icap=5
               else if (keyword(i).eq.'sigmoidal') then
                  icap=6
               else
                  errmsg='Unknown CAP type: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'cap_strength') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) capstr
            else
               goto 100
            endif

         else if (keyword(i).eq.'projection') then
            lprojcap=.true.
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'initial') then
                  iprojcap=1
               else if (keyword(i).eq.'all') then
                  iprojcap=2
               else if (keyword(i).eq.'all_bb') then
                  iprojcap=3
               else if (keyword(i).eq.'all_bu'.or.&
                    keyword(i).eq.'all_ub') then
                  iprojcap=4
               else
                  errmsg='Unkown projection type: '//trim(keyword(i))
                  call error_control
               endif
            else
               iprojcap=1
            endif
            ! Energy limit for states entering into the projector
            if (iprojcap.eq.2.or.iprojcap.eq.3.or.iprojcap.eq.4) then
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) projlim
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     call convert_energy(keyword(i),projlim)
                  endif
               endif
            endif
            
         else if (keyword(i).eq.'cap_order') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) capord
            else
               goto 100
            endif

         else if (keyword(i).eq.'grid_rad') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nrad(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) nrad(2)
               else
                  errmsg='The no. outer radial points has not been given'
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'grid_ang') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nang(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) nang(2)
               else
                  errmsg='The no. outer angular points has not been given'
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'cap_box') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'auto') then
                  ! Automatic CAP box parameterisation
                  lautobox=.true.
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) densthrsh
                  else
                     errmsg='The density threshold for the automatic CAP &
                          box parameterisation has not been given'
                     call error_control
                  endif
               else
                  ! User specified CAP box
                  read(keyword(i),*) boxpar(1)
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) boxpar(2)
                  else
                     errmsg='Only 1 out of 3 CAP box dimensions were given'
                     call error_control
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) boxpar(3)
                  else
                     errmsg='Only 3 out of 3 CAP box dimensions were given'
                     call error_control
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     do n=1,3
                        call convert_length(keyword(i),boxpar(n))
                     enddo
                  endif
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'flux') then
            lflux=.true.

         else if (keyword(i).eq.'cap_diag') then
            lcapdiag=.true.
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif

      return
      
    end subroutine rdcapinp

!#######################################################################

    subroutine rdpropagationinp

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer          :: i,n
      real(dp)         :: theta,phi
      character(len=2) :: ai
      
!-----------------------------------------------------------------------
! Read to the propagation section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'propagation_section') goto 1

!-----------------------------------------------------------------------
! Read to the popagation parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-propagation_section') then

10       continue
         i=i+1
         
         if (keyword(i).eq.'tfinal') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tfinal
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'tout') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tout
            else
               goto 100
            endif

         else if (keyword(i).eq.'tol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) proptol
            else
               goto 100
            endif

         else if (keyword(i).eq.'krydim') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) kdim
            else
               goto 100
            endif

         else if (keyword(i).eq.'representation') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'isr') then
                  tdrep=1
               else if (keyword(i).eq.'eigenstate') then
                  tdrep=2
               else
                  errmsg='Unknown representation: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif

      return
      
    end subroutine rdpropagationinp

!#######################################################################

    subroutine rdlaserinp(ip)

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer          :: ip,count,i,n
      real(dp)         :: theta,phi
      character(len=2) :: ai

!-----------------------------------------------------------------------
! Read to the current laser section
!-----------------------------------------------------------------------
      rewind(iin)

      count=0
      
1     call rdinp(iin)
      if (keyword(1).ne.'laser_section') goto 1
      count=count+1
      if (count.ne.ip) goto 1
      
!-----------------------------------------------------------------------
! Read the laser parameters
!-----------------------------------------------------------------------
5     call rdinp(iin)
      
      i=0
      
      if (keyword(1).ne.'end-laser_section') then

10       continue
         i=i+1

         if (keyword(i).eq.'pulse_vec') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'polar') then
                  ! Spherical polar specification
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) theta
                  else
                     errmsg='The polar angle was not given'
                     call error_control
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) phi
                  else
                     errmsg='The azimuthal angle was not given'
                     call error_control
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     call convert_angle(keyword(i),theta)
                     call convert_angle(keyword(i),phi)
                  endif
                  pulse_vec(1,ip)=sin(theta)*cos(phi)
                  pulse_vec(2,ip)=sin(theta)*sin(phi)
                  pulse_vec(3,ip)=cos(theta)
               else
                  ! Cartesian specification
                  read(keyword(i),*) pulse_vec(1,ip)
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) pulse_vec(2,ip)
                  else
                     errmsg='Only one argument out of three was given &
                          with the pulse_vec keyword'
                  endif
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) pulse_vec(3,ip)
                  else
                     errmsg='Only two arguments out of three was given &
                          with the pulse_vec keyword'
                  endif
               endif
               ! Normalisation of the vector
               pulse_vec(:,ip)=pulse_vec(:,ip)&
                    /sqrt(dot_product(pulse_vec(:,ip),pulse_vec(:,ip)))
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse_freq') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) freq(ip)
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse_strength') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) strength(ip)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  call convert_intensity(keyword(i),strength(ip))
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse_t0') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) t0(ip)
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse_phase') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) phase(ip)
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse'.or.keyword(i).eq.'carrier') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'cos') then
                  ipulse(ip)=1
               else if (keyword(i).eq.'sin') then
                  ipulse(ip)=2
               else
                  errmsg='Unknown pulse type: '//trim(keyword(i))
                  call error_control
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'envelope') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Read the envelope type
               if (keyword(i).eq.'cos2') then
                  ienvelope(ip)=1
                  nenvpar(ip)=2
               else if (keyword(i).eq.'sin2-ramp') then
                  ienvelope(ip)=2
                  nenvpar(ip)=1
               else if (keyword(i).eq.'box') then
                  ienvelope(ip)=3
                  nenvpar(ip)=2
               else if (keyword(i).eq.'sin2') then
                  ienvelope(ip)=4
                  nenvpar(ip)=2
               else
                  errmsg='Unknown evelope type: '//trim(keyword(i))
                  call error_control
               endif
               ! Read the envelope parameters
               do n=1,nenvpar(ip)
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) envpar(n,ip)
                  else
                     write(ai,'(i2)') nenvpar(ip)
                     errmsg='Not enough evelope parameters were given: '&
                          //trim(adjustl(ai))//' parameters expected.'
                     call error_control
                  endif
               enddo
            else
               goto 100
            endif
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         call error_control
         
      endif
      
      return
      
    end subroutine rdlaserinp
      
!#######################################################################

    subroutine rddiffbaslist

      use parameters
      use parsemod
      use iomod
      use channels

      implicit none

      integer :: i,l

!-----------------------------------------------------------------------
! Read to the diffuse basis section
!-----------------------------------------------------------------------
      rewind(iin)

1     call rdinp(iin)
      if (keyword(1).ne.'diffuse_section') goto 1

!-----------------------------------------------------------------------
! Read the diffuse basis information
!-----------------------------------------------------------------------
      5    call rdinp(iin)
      
      i=0

      if (keyword(1).ne.'end-diffuse_section') then
         
10       continue
         i=i+1
         
         if (keyword(i).eq.'s') then
            ndiff(1)=ndiff(1)+1
            i=i+1
            read(keyword(i),*) difflist(1,ndiff(1))
         else if (keyword(i).eq.'p') then
            ndiff(2)=ndiff(2)+1
            i=i+1
            read(keyword(i),*) difflist(2,ndiff(2))
         else if (keyword(i).eq.'d') then
            ndiff(3)=ndiff(3)+1
            i=i+1
            read(keyword(i),*) difflist(3,ndiff(3))
         else if (keyword(i).eq.'f') then
            ndiff(4)=ndiff(4)+1
            i=i+1
            read(keyword(i),*) difflist(4,ndiff(4))
         else if (keyword(i).eq.'g') then
            ndiff(5)=ndiff(5)+1
            i=i+1
            read(keyword(i),*) difflist(5,ndiff(5))
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            call error_control
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
      endif

      return

    end subroutine rddiffbaslist

!#######################################################################

    subroutine getdiffinfo(n,l,string)
      
      implicit none
      
      integer          :: n,l,i,k
      character(len=*) :: string
      character(len=5) :: labels
      
      labels='spdfg'
      
      do i=1,len_trim(string)
         k=index(labels,string(i:i))
         if (k.ne.0) then
            read(string(1:i-1),*) n
            l=k
            exit
         endif
      enddo

      return

    end subroutine getdiffinfo

!#######################################################################

    subroutine rdgeometry

      use parameters
      use parsemod
      use channels

      implicit none

      integer :: i,j

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(xcoo(ncoo))
      allocate(atlbl(natm))

!-----------------------------------------------------------------------
! Read the Cartesian coordinates from the input file
!-----------------------------------------------------------------------
      rewind(iin)

5     continue
      call rdinp(iin)
      if (keyword(1).ne.'geometry') goto 5

      do i=1,natm
         call rdinp(iin)
         atlbl(i)=keyword(1)
         do j=1,3
            read(keyword(j+1),*) xcoo(i*3-3+j)
         enddo
      enddo

      return 

    end subroutine rdgeometry

!#######################################################################

    subroutine alloclaserpar

      use parameters

      implicit none

      allocate(pulse_vec(3,npulse))
      pulse_vec=0.0d0

      allocate(freq(npulse))
      freq=0.0d0

      allocate(strength(npulse))
      strength=0.0d0

      allocate(t0(npulse))
      t0=0.0d0

      allocate(phase(npulse))
      phase=0.0d0

      allocate(ipulse(npulse))
      ipulse=0

      allocate(ienvelope(npulse))
      ienvelope=0

      allocate(envpar(mxenvpar,npulse))
      envpar=0.0d0

      allocate(nenvpar(npulse))
      nenvpar=0
      
      return
      
    end subroutine alloclaserpar
      
!#######################################################################

    subroutine convert_intensity(unit,val)

      use constants
      use iomod
      
      implicit none
      
      real(dp)         :: val
      character(len=*) :: unit

      if (unit.eq.'au') then
         ! Do nothing, atomic units are the default         
      else if (unit.eq.'winvcm2') then
         ! W cm^-2 -> au
         val=sqrt(val/3.5e+16_dp)
      else
         ! Unrecognised unit
         errmsg='Unrecognised unit: '//trim(unit)
         call error_control
      endif
      
      return
      
    end subroutine convert_intensity

!#######################################################################

    subroutine convert_length(unit,val)

      use constants
      use iomod
      
      implicit none
      
      real(dp)         :: val
      character(len=*) :: unit

      if (unit.eq.'au') then
         ! Do nothing, atomic units are the default         
      else if (unit.eq.'angstrom') then
         ! Angstrom -> au
         val=val*1.889725989d0
      else
         ! Unrecognised unit
         errmsg='Unrecognised unit: '//trim(unit)
         call error_control
      endif
      
      return
            
    end subroutine convert_length

!#######################################################################

    subroutine convert_angle(unit,val)

      use constants
      use iomod
      
      implicit none

      real(dp)         :: val
      character(len=*) :: unit

      if (unit.eq.'rad'.or.unit.eq.'radians') then
         ! Do nothing, radians are the default         
      else if (unit.eq.'deg'.or.unit.eq.'degrees') then
         ! radians -> degrees
         val=val*pi/180.0d0
      else
         ! Unrecognised unit
         errmsg='Unrecognised unit: '//trim(unit)
         call error_control
      endif
      
      return
      
    end subroutine convert_angle

!#######################################################################

    subroutine convert_energy(unit,val)

      use constants
      use parameters
      use iomod
      
      implicit none

      real(dp)         :: val
      character(len=*) :: unit

      if (unit.eq.'au') then
         ! Do nothing, atomic units are the default
      else if (unit.eq.'ev') then
         val=val/eh2ev
      else
         ! Unrecognised unit
         errmsg='Unrecognised unit: '//trim(unit)
         call error_control
      endif
      
      return
      
    end subroutine convert_energy
    
!#######################################################################

  end module rdinput
