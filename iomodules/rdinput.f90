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
                     errmsg='Unknown unit: '//trim(keyword(i))
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
      if (.not.energyonly.and..not.ldiagfinal&
           .and..not.ldyson.and..not.lrixs.and..not.ltpa) then

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
                  rlxtype=1
               else if (keyword(i).eq.'relaxation_liu') then
                  solver=2
                  rlxtype=2
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
            
         else if (keyword(i).eq.'ortho') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'unmodified') then
                  rlxortho=1
               else if (keyword(i).eq.'modified') then
                  rlxortho=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
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
                  rlxtype_f=1
               else if (keyword(i).eq.'relaxation_liu') then
                  solver_f=2
                  rlxtype_f=2
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
            
         else if (keyword(i).eq.'ortho') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'unmodified') then
                  rlxortho_f=1
               else if (keyword(i).eq.'modified') then
                  rlxortho_f=2
               else
                  errmsg='Unknown keyword: '//trim(keyword(i))
                  call error_control
               endif
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

  end module rdinput
