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
      logical            :: iscvs,energyonly,ldiag,llanc

!-----------------------------------------------------------------------
! Set 'traps'
!-----------------------------------------------------------------------
      energyonly=.false.
      ldiag=.false.
      llanc=.false.
      iscvs=.false.

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
               else if (keyword(i).eq.'adc2'.or.keyword(i).eq.'adc2-s') then
                  method=2
               else if (keyword(i).eq.'adc2-x') then
                  method=3
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
                  ! Assign any integer i>0: this is 
                  ! necessary to get around checks on the 
                  ! statenumber
                  statenumber=999
               else
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
               if (keyword(i).eq.'kbj'.or.keyword(i).eq.'kbj_cont') then
                  difftype=1
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

      character(len=120) :: msg
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
! (1) We are ionizing from an excited state, or;
! (2) We are performing an energy-only calculation.
!-----------------------------------------------------------------------
      if (.not.lautospec) then

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
            else if (maxiter.eq.0) then
               msg='The maximum no. iterations for the initial space &
                    diagonalisation has not been &
                    given'
               goto 999
            else if (dmain.eq.0) then
               msg='The block size for the initial space diagonalisation &
                    has not been given'
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
         else if (maxiter_f.eq.0) then
            msg='The maximum no. iterations for the final space &
                 diagonalisation has not been given'
            goto 999
         else if (dmain_f.eq.0) then
            msg='The block size for the final space diagonalisation &
                 has not been given'
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
      if (lautospec) then

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

      endif

!-----------------------------------------------------------------------
! Filter diagonalisation states section
!-----------------------------------------------------------------------
      if (lfdstates) then

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
            msg='The name of the fdiag data filehas not been given'
            goto 999
         endif
         
         ! State selection file
         if (fdiagsel.eq.'') then
            msg='The name of the fdiag state selection file has &
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

         ! CAP width
         if (capwid.eq.0.0d0) then
            msg='The CAP width has not been given'
            goto 999
         endif         

         ! Propagation section
         if (.not.lpropagation) then
            msg='The propagation sesction has not been given'
            goto 999
         endif
         
      endif

!-----------------------------------------------------------------------
! Propagation section
!-----------------------------------------------------------------------
      if (lpropagation) then

         ! Pulse type
         if (ipulse.eq.0) then
            msg='The pulse type has not been given'
            goto 999
         endif

         ! Pulse envelope type
         if (ienvelope.eq.0) then
            msg='The pulse envelope type has not been given'
            goto 999
         endif
         
         ! Laser-molecule orientation
         if (sum(abs(pulse_vec)).eq.0.0d0) then
            msg='The laser pulse orientation has not been given'
            goto 999
         endif

         ! Laser frequency
         if (freq.eq.0.0d0) then
            msg='The laser frequency has not been given'
            goto 999
         endif

         ! Laser strength
         if (strength.eq.0.0d0) then
            msg='The laser strength has not been given'
            goto 999
         endif
         
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

      integer :: i
      
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
               if (keyword(i).eq.'sigmoidal') then
                  icap=1
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

         else if (keyword(i).eq.'cap_width') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) capwid
               if (keyword(i+1).eq.',') then
                  i=i+2
                  call convert_length(keyword(i),capwid)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'projection') then
            lprojcap=.true.
            
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

         if (keyword(i).eq.'pulse_vec') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) pulse_vec(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) pulse_vec(2)
               else
                  errmsg='Only one argument out of three was given &
                       with the pulse_vec keyword'
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) pulse_vec(3)
               else
                  errmsg='Only two arguments out of three was given &
                       with the pulse_vec keyword'
               endif
               ! Normalisation of the vector
               pulse_vec=pulse_vec/sqrt(dot_product(pulse_vec,pulse_vec))
            else
               goto 100
            endif

         else if (keyword(i).eq.'tfinal') then
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

         else if (keyword(i).eq.'pulse_freq') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) freq
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse_strength') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) strength
               if (keyword(i+1).eq.',') then
                  i=i+2
                  call convert_intensity(keyword(i),strength)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'pulse') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'cos') then
                  ipulse=1
               else if (keyword(i).eq.'sin') then
                  ipulse=2
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
                  ienvelope=1
                  nenvpar=2
               else if (keyword(i).eq.'sin2-ramp') then
                  ienvelope=2
                  nenvpar=1
               else
                  errmsg='Unknown evelope type: '//trim(keyword(i))
                  call error_control
               endif
               ! Read the envelope parameters
               do n=1,nenvpar
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     read(keyword(i),*) envpar(n)
                  else
                     write(ai,'(i2)') nenvpar
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
      
    end subroutine rdpropagationinp
      
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

    subroutine convert_intensity(unit,val)

      use constants
      use iomod
      
      implicit none
      
      real(d)          :: val
      character(len=*) :: unit

      if (unit.eq.'au') then
         ! Do nothing, atomic units are the default         
      else if (unit.eq.'winvcm2') then
         ! W cm^-2 -> au
         val=sqrt(val/3.5e+16_d)
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
      
      real(d)          :: val
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

  end module rdinput
