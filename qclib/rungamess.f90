  module rungamess

  contains

!#######################################################################

    subroutine rungamess_main
      
      use channels
      use iomod
      use parameters
      
      implicit none
      
      integer            :: igms,k
      character(len=120) :: filename

!-----------------------------------------------------------------------
! Read the AO basis from file
!-----------------------------------------------------------------------
      call rdbas

!-----------------------------------------------------------------------
! Calculate the exponents of any additional diffuse functions
!-----------------------------------------------------------------------
      if (difftype.gt.0) call calc_diffexp

!-----------------------------------------------------------------------
! Determine where to place the additional diffuse functions
!-----------------------------------------------------------------------
      ldiffcom=.false.
      contcent=0
      if (difftype.gt.0.or.lfakeip) call place_diff

!-----------------------------------------------------------------------
! Write the gamess input file
!-----------------------------------------------------------------------      
      k=(index(ain,'.inp'))
      filename=ain(1:k-1)//'_gamess.inp'
      call freeunit(igms)
      open(igms,file=filename,form='formatted',status='unknown')

      call wrgamessinp('inp',igms)

      close(igms)

!-----------------------------------------------------------------------
! Run GAMESS
!-----------------------------------------------------------------------
      call rungms

!-----------------------------------------------------------------------
! Edit the .dat file
!-----------------------------------------------------------------------
      call rewrite_dat

      return
      
    end subroutine rungamess_main

!#######################################################################

    subroutine rdbas
      
      use iomod
      use parameters
      use parsemod

      implicit none

      integer            :: i,maxfunc,ibasfile
      character(len=120) :: string
      character(len=250) :: adcdir
      logical            :: lbas

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Number of AOs per atom
      allocate(naogms(natm))
      naogms=0

      ! Number of primitive functions per AO for each atom
      allocate(nprim(natm,maxao))
      nprim=0

      ! Coefficients
      allocate(aocoeff(natm,maxao,maxprim))
      aocoeff=0.0d0

      ! Exponents
      allocate(aoexp(natm,maxao,maxprim))
      aoexp=0.0d0

      ! Angular momentum quantum numbers
      allocate(ilquant(natm,maxao))
      ilquant=-1

!-----------------------------------------------------------------------
! Open the basis set file
!-----------------------------------------------------------------------
      call get_environment_variable("ADC_DIR",adcdir)

      call freeunit(ibasfile)

      open(ibasfile,file=trim(adcdir)//'/qclib/basis.dat',&
           form='formatted',status='old')      

!-----------------------------------------------------------------------
! Read the basis file
!-----------------------------------------------------------------------
5     continue      
      call rdinp(ibasfile)

      if (.not.lend) then

         ! Determine whether the current line is the start of a
         ! block of AOs of interest
         if (keyword(2).eq.basname) then
            lbas=.true.
         else
            lbas=.false.
         endif
 
         ! If we are at the start of a block of interest, then
         ! read in the coefficients
         if (lbas) call rdbas_1block(ibasfile)

         goto 5

      endif

!-----------------------------------------------------------------------
! Close the basis set file
!-----------------------------------------------------------------------
      close(ibasfile)

!-----------------------------------------------------------------------
! Check that we have found the requested basis for all atoms
!-----------------------------------------------------------------------
      do i=1,natm
         if (naogms(i).eq.0) then
            errmsg='The basis '//trim(basname)&
                 //' could not be found for atom '//trim(atlbl(i))
            call error_control
         endif
      enddo

      return

    end subroutine rdbas

!#######################################################################

    subroutine rdbas_1block(ibasfile)      

      use parameters
      use parsemod

      implicit none

      integer            :: ibasfile,i,j,norb
      character(len=2)   :: curratm
      character(len=120) :: string
      logical            :: latm

      real*8, dimension(maxao,maxprim) :: coeff,expnt
      integer, dimension(maxao)        :: npbas,iang

!-----------------------------------------------------------------------
! Go back to the start of the block
!-----------------------------------------------------------------------
      backspace(ibasfile)

!-----------------------------------------------------------------------
! Read the atom label
!-----------------------------------------------------------------------
      call rdinp(ibasfile)

      curratm=keyword(1)

      ! If the atom does not appear in the molecule then skip the
      ! current block...
      latm=.false.
      do i=1,natm
         if (curratm.eq.atlbl(i)) latm=.true.
      enddo
      if (.not.latm) goto 100

      ! ...else read the current block
      norb=0
10    continue
      norb=norb+1
      call rdinp(ibasfile)
      if (keyword(1).eq.'s') then
         iang(norb)=1
      else if (keyword(1).eq.'p') then
         iang(norb)=2
      else if (keyword(1).eq.'d') then
         iang(norb)=3
      else if (keyword(1).eq.'f') then
         iang(norb)=4
      else if (keyword(1).eq.'g') then
         iang(norb)=5
      endif

      read(keyword(2),*) npbas(norb)
      
      do i=1,npbas(norb)
         call rdinp(ibasfile)
         read(keyword(2),*) expnt(norb,i)
         read(keyword(3),*) coeff(norb,i)
      enddo

      ! If we are not at the end of the current block, then read
      ! the next AO
      read(ibasfile,'(a)',end=888) string
      if (string.ne.'') then
         backspace(ibasfile)
         goto 10
      endif
888   continue
      ! Go back one line so that we don't die in rdbas
      backspace(ibasfile)

      ! Fill in the coefficient, exponent and angular momentum quantum
      ! number arrays for all atoms of the current type
      do i=1,natm
         if (atlbl(i).eq.curratm) then
            naogms(i)=norb
            aocoeff(i,:,:)=coeff(:,:)
            aoexp(i,:,:)=expnt(:,:)
            nprim(i,:)=npbas(:)
            ilquant(i,:)=iang(:)-1
         endif
      enddo

      return

100   continue
      ! Skip the current block
110   read(ibasfile,'(a)',end=999) string
      if (string.ne.'') goto 110

999   continue

      return

    end subroutine rdbas_1block

!#######################################################################

    subroutine calc_diffexp

      use parameters

      implicit none

      integer              :: i,j,k,maxfunc
      real*8, dimension(5) :: minexp

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Exponents for the uncontracted diffuse functions
      maxfunc=0
      do i=1,5
         if (ndiff(i).gt.maxfunc) maxfunc=ndiff(i)
      enddo

      allocate(diffexp(5,maxfunc))
      diffexp=0.0d0

!-----------------------------------------------------------------------
! Determine the exponents of the most diffuse non-contracted AOs in the
! parent basis
!-----------------------------------------------------------------------
      minexp=99999.9d0
      do i=1,natm
         do j=1,naogms(i)
            if (nprim(i,j).eq.1) then
               k=ilquant(i,j)+1
               if (aoexp(i,j,1).lt.minexp(k)) then
                  minexp(k)=aoexp(i,j,1)
               endif
            endif
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the exponents of the additional diffuse functions
!-----------------------------------------------------------------------
      do i=1,5
         do j=1,ndiff(i)
            if (difftype.eq.1) then               
               ! KBJ-type diffuse functions
               diffexp(i,j)=kbjexp(i,minexp(i))
               minexp(i)=diffexp(i,j)
            else if (difftype.eq.2) then
               ! Even-tempered diffuse functions
               write(6,'(a)') 'You need to write the even-tempered basis &
                    function code...'
               STOP
            endif
         enddo
      enddo

      return

    end subroutine calc_diffexp

!#######################################################################

    subroutine place_diff

      use channels
      use parameters

      implicit none

      integer                           :: i,j,ilbl,count,currlquant
      integer, dimension(maxao)         :: tmpprim,tmplquant
      real(d)                           :: dist,mindist
      real(d), dimension(3)             :: xcom
      real(d), dimension(maxao,maxprim) :: tmpcoeff,tmpexp
 
!-----------------------------------------------------------------------
! Determine whether any atoms are close (<0.15 Angstrom) to the centre
! of mass
!-----------------------------------------------------------------------
      call getcom(xcom)

      ldiffcom=.true.
      mindist=0.15d0
      do i=1,natm
         dist=0.0d0
         do j=1,3
            dist=dist+(xcoo(i*3-3+j)-xcom(j))**2
         enddo
         dist=sqrt(dist)
         if (dist.lt.mindist) then
            ldiffcom=.false.
            ilbl=i
            mindist=dist
         endif
      enddo

      if (ldiffcom) then
         write(ilog,'(/,2x,a)') 'Placing diffuse functions at the &
              centre of mass'
      else
         write(ilog,'(/,2x,a,x,i2)') 'Placing diffuse functions on &
              atom number',ilbl
      endif

!-----------------------------------------------------------------------
! If the diffuse functions are to be placed at an atomic centre, then
! re-arrange the corresponding AO parameter arrays
!-----------------------------------------------------------------------
      if (.not.ldiffcom) then

         if (lfakeip) contcent=ilbl

         count=0
         currlquant=0         

         ! Parent functions and diffuse functions up to 
         ! l=iquant(ilbl,naogms(ilbl))-1
         do i=1,naogms(ilbl)

            ! AOs from the set of additional diffuse functions
            if (ilquant(ilbl,i).gt.currlquant) then
               ! Add any diffuse functions of angular momentum quantum
               ! number currlquant-1
               do j=1,ndiff(currlquant+1)
                  count=count+1
                  tmpprim(count)=1
                  tmpcoeff(count,1)=1.0d0
                  tmpexp(count,1)=diffexp(currlquant+1,j)
                  tmplquant(count)=currlquant
               enddo
               currlquant=ilquant(ilbl,i)
            endif

            ! AOs from the parent basis
            count=count+1
            tmpprim(count)=nprim(ilbl,i)
            tmpcoeff(count,:)=aocoeff(ilbl,i,:)
            tmpexp(count,:)=aoexp(ilbl,i,:)
            tmplquant(count)=ilquant(ilbl,i)
            
         enddo

         ! Remaining sets of diffuse functions for
         ! l=iquant(ilbl,naogms(ilbl))-1,...,4
         do i=currlquant,4
            do j=1,ndiff(i+1)
               count=count+1
               tmpprim(count)=1
               tmpcoeff(count,1)=1.0d0
               tmpexp(count,1)=diffexp(i+1,j)
               tmplquant(count)=i
            enddo
         enddo

         ! Update the no. AOs
         naogms(ilbl)=count

         ! Update the AO basis arrays
         nprim(ilbl,:)=tmpprim(:)
         aocoeff(ilbl,:,:)=tmpcoeff(:,:)
         aoexp(ilbl,:,:)=tmpexp(:,:)
         ilquant(ilbl,:)=tmplquant(:)

      endif

      return
      
    end subroutine place_diff

!#######################################################################

    function kbjexp(i,thrsh)

      implicit none

      integer :: i,l,k
      real*8  :: kbjexp,thrsh,al,bl

!----------------------------------------------------------------------- 
! Set the coefficients a_l and b_l according to the principal quantum
! number
!----------------------------------------------------------------------- 
      l=i-1

      select case(l)

      case(0) ! s-functions
         al=0.584342d0
         bl=0.424483d0
         
      case(1) ! p-functions
         al=0.452615d0
         bl=0.309805d0

      case(2) ! d-functions
         al=0.382362d0
         bl=0.251333d0

      case(3) ! f-functions
         al=0.337027d0
         bl=0.215013d0

      case(4) ! g-functions
         al=0.304679d0
         bl=0.189944d0

      end select

!----------------------------------------------------------------------- 
! Calculate the next KBJ exponent
!-----------------------------------------------------------------------       
      do k=1,200
         kbjexp=0.25d0*(1.0d0/(al*k+bl)**2)
         if (kbjexp.lt.thrsh) exit
      enddo

      return

    end function kbjexp

!#######################################################################

    subroutine wrgamessinp(type,iout)

      use parameters
      use iomod

      implicit none

      integer              :: iout,i,j,k
      real*8               :: atnum
      real*8, dimension(3) :: xcom
      character(len=2)     :: atupper
      character(len=1)     :: lquantlbl
      character(len=120)   :: aout
      character(len=3)     :: type

!-----------------------------------------------------------------------
! If diffuse functions are to be placed at the centre of mass, then
! determine the centre of mass
!-----------------------------------------------------------------------
      call getcom(xcom)

!-----------------------------------------------------------------------
! CONTRL section
!-----------------------------------------------------------------------
      if (pntgroup.eq.'c1'.or.pntgroup.eq.'') then
         if (type.eq.'inp') write(iout,'(a,/)') &
              ' $CONTRL SCFTYP=RHF RUNTYP=ENERGY COORD=UNIQUE NZVAR=0 $END'
      else
         if (type.eq.'inp') write(iout,'(a,/)') &
              ' $CONTRL SCFTYP=RHF RUNTYP=ENERGY COORD=PRINAXIS NZVAR=0 $END'
      endif

!-----------------------------------------------------------------------
! Start of the DATA section
!-----------------------------------------------------------------------
      write(iout,'(a,/)') ' $DATA'
      if (pntgroup.eq.'c1'.or.pntgroup.eq.'') then
         write(iout,'(a)') 'C1       0'
      else if (pntgroup.eq.'c2v') then
         write(iout,'(a,/)') 'Cnv       2'
      else
         errmsg='The point group '//trim(pntgroup)//' is not supported'
         call error_control
      endif

!-----------------------------------------------------------------------
! Centre-of-mass centred basis functions
!-----------------------------------------------------------------------      
      if (ldiffcom) then
         write(iout,'(a1,11x,a3,3(6x,F12.10))') 'X','0.0',&
              (xcom(i),i=1,3)
         ! Fake continuum orbital
         if (lfakeip) then
            write(iout,'(a)') '   S       1'
            write(iout,'(5x,a)') &
                 '1            1d-30  1.000000'
         endif
         ! s functions
         if (ndiff(1).gt.0) then
            do i=1,ndiff(1)
               write(iout,'(a)') '   S       1'
               write(iout,'(5x,a1,12x,2(2x,F8.6))') &
                    '1',diffexp(1,i),1.0d0
            enddo
         endif
         ! p functions
         if (ndiff(2).gt.0) then
            do i=1,ndiff(2)
               write(iout,'(a)') '   P       1'
               write(iout,'(5x,a1,12x,2(2x,F8.6))') &
                    '1',diffexp(2,i),1.0d0
            enddo
         endif
         ! d functions
         if (ndiff(3).gt.0) then
            do i=1,ndiff(3)
               write(iout,'(a)') '   D       1'
               write(iout,'(5x,a1,12x,2(2x,F8.6))') &
                    '1',diffexp(3,i),1.0d0
            enddo
         endif
         ! f functions
         if (ndiff(4).gt.0) then
            do i=1,ndiff(4)
               write(iout,'(a)') '   F       1'
               write(iout,'(5x,a1,12x,2(2x,F8.6))') &
                    '1',diffexp(4,i),1.0d0
            enddo
         endif
         ! g functions
         if (ndiff(5).gt.0) then
            do i=1,ndiff(5)
               write(iout,'(a)') '   G       1'
               write(iout,'(5x,a1,12x,2(2x,F8.6))') &
                    '1',diffexp(5,i),1.0d0
            enddo
         endif
         ! Blank line
         write(iout,*)
      endif

!-----------------------------------------------------------------------      
! Atom centred basis functions
!-----------------------------------------------------------------------      
      do i=1,natm
         ! Atom label, atomic number and coordinates
         atupper=uppercase(atlbl(i))
         atnum=atomic_number(atlbl(i))
         write(iout,'(a2,9x,F4.1,3(4x,F14.10))') atupper,atnum,&
              (xcoo(j),j=i*3-2,i*3)
         ! 'Continuum' orbital for IP-ADC calculations
         if (lfakeip.and.i.eq.contcent) then
            write(iout,'(a)') '   S       1'
            write(iout,'(5x,a)') &
                 '1            1d-30  1.000000'            
         endif
         ! AO exponents and contraction coefficients
         do j=1,naogms(i)
            lquantlbl=get_lquantlbl(ilquant(i,j))
            write(iout,'(3x,a1,6x,i2)') lquantlbl,nprim(i,j)
            do k=1,nprim(i,j)
               write(iout,'(4x,i2,8x,F18.10,1x,F11.8)') &
                    k,aoexp(i,j,k),aocoeff(i,j,k)
            enddo
         enddo
         ! Blank line
         write(iout,*)
      enddo

!-----------------------------------------------------------------------
! End of the DATA section
!-----------------------------------------------------------------------
      write(iout,'(a,/)') ' $END'

      return

    end subroutine wrgamessinp

!#######################################################################

       function uppercase(string1) result(string2)

      implicit none

      integer                     :: dim,i,j
      character(len=*)            :: string1
      character(len=len(string1)) :: string2

      do i = 1, len(string1)
         j = iachar(string1(i:i))
         if (j>= iachar("a") .and. j<=iachar("z") ) then
            string2(i:i) = achar(iachar(string1(i:i))-32)
         else
            string2(i:i) = string1(i:i)
         endif
     enddo

      return

    end function uppercase

!#######################################################################

    function atomic_number(lbl) result(num)

      real*8           :: num
      character(len=*) :: lbl

      if (lbl.eq.'h') then
         num=1.0d0
      else if (lbl.eq.'c') then
         num=6.0d0
      else if (lbl.eq.'n') then
         num=7.0d0
      else if (lbl.eq.'s') then
         num=16.0d0
      else
         write(6,'(/,2(2x,a),/)') 'Atomic no. not known for atom:',lbl
      endif

      return

    end function atomic_number

!#######################################################################

    function get_lquantlbl(l) result(lbl)

      implicit none

      integer          :: l
      character(len=1) :: lbl

      if (l.eq.0) then
         lbl='S'
      else if (l.eq.1) then
         lbl='P'
      else if (l.eq.2) then
         lbl='D'
      else if (l.eq.3) then
         lbl='F'
      else if (l.eq.4) then
         lbl='G'
      endif

      return

    end function get_lquantlbl

!#######################################################################

    subroutine getcom(xcom)

      use parameters

      implicit none

      integer              :: i,j
      real*8, dimension(3) :: xcom
      real*8               :: m,tmass

      xcom=0.0d0
      tmass=0.0d0

      do i=1,natm
         m=mass(atlbl(i))
         tmass=tmass+m
         do j=1,3
            xcom(j)=xcom(j)+xcoo(i*3-3+j)*m
         enddo
      enddo

      xcom=xcom/tmass

      return

    end subroutine getcom 

!#######################################################################

    function mass(label)

      implicit none

      real*8           :: mass
      character(len=*) :: label

      if (label.eq.'h') then
         mass=1.00794d0
      else if (label.eq.'c') then
         mass=12.0107d0
      else if (label.eq.'n') then
         mass=14.0067d0
      else if (label.eq.'s') then
         mass=32.065d0
      else if (label.eq.'o') then
         mass=15.9994d0
      else
         write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
         STOP
      endif

      return

    end function mass

!#######################################################################

    subroutine rungms
      
      use channels
      use iomod
      
      implicit none

      integer            :: unit,k,nthreads
      character(len=120) :: string,targ,stem,filename
      character(len=350) :: workdir,adcdir,command
      logical            :: found

      integer            :: omp_get_num_threads

      write(ilog,'(/,2x,a)') 'Running GAMESS...'

!-----------------------------------------------------------------------
! Run GAMESS
!-----------------------------------------------------------------------
      ! Determine the no. threads
      !$omp parallel
      nthreads=omp_get_num_threads()
      !$omp end parallel

      ! BODGE
      nthreads=1

      ! Determine the working and ADC directories
      call getcwd(workdir)

      call get_environment_variable("ADC_DIR",adcdir)
      if (adcdir.eq.'') then
         errmsg='ADC_DIR has not been set. Quitting.'
         call error_control
      endif

      ! Write the call to the run_gamess script
      k=(index(ain,'.inp'))
      stem=ain(1:k-1)

      command=trim(adcdir)//'/run_gamess '//trim(stem)//'_gamess '&
           //trim(workdir)

      k=len_trim(command)
      if (nthreads.lt.10) then
         write(command(k+1:k+2),'(1x,i1)') nthreads
      else
         write(command(k+1:k+3),'(1x,i2)') nthreads
      endif

      ! Execute the run_gamess script
      call system(command)

!-----------------------------------------------------------------------
! Check that the calculation didn't fail
!-----------------------------------------------------------------------
      call freeunit(unit)

      filename=trim(stem)//'_gamess.log'
      open(unit,file=filename,form='formatted',status='old')

      targ='EXECUTION OF GAMESS TERMINATED NORMALLY'
      found=.false.
10    read(unit,'(a)',end=20) string
      print*,string
      if (index(string,trim(targ)).ne.0) found=.true.      
      goto 10
      
20    continue
      if (found) then
         write(ilog,'(2x,a,/)') 'Successful termination of GAMESS'
      else
         errmsg='Execution of GAMESS failed...'
         call error_control
      endif

      close(unit)

!-----------------------------------------------------------------------
! Rename the log and dat files to be consitent with the rest of the
! program
!-----------------------------------------------------------------------
      filename=trim(stem)//'_gamess.log'
      command='mv '//trim(filename)//' gamess.log'
      call system(command)

      filename=trim(stem)//'_gamess.dat'
      command='mv '//trim(filename)//' gamess.dat'
      call system(command)

      return

    end subroutine rungms

!#######################################################################

    subroutine rewrite_dat

      use iomod

      implicit none

      integer            :: itmp,idat
      character(len=120) :: string,adat

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
      call freeunit(itmp)
      open(itmp,file='tmp.dat',form='formatted',status='unknown')

      call freeunit(idat)
      
      open(idat,file='gamess.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Read to the end of the DATA section of the gamess.dat file 
!-----------------------------------------------------------------------
10    continue
      read(idat,'(a)') string
      if (index(string,'$END').eq.0) goto 10

!-----------------------------------------------------------------------
! Write the new data section
!-----------------------------------------------------------------------
      call wrgamessinp('dat',itmp)

!-----------------------------------------------------------------------
! Add in the rest of the .dat file
!-----------------------------------------------------------------------
      backspace(itmp)

15    read(idat,'(a)',end=20) string
      write(itmp,'(a)') trim(string)
      goto 15

20    continue

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
      close(itmp)
      close(idat)

!-----------------------------------------------------------------------
! Clean up
!-----------------------------------------------------------------------
      call system('cp tmp.dat gamess.dat')
      call system('rm tmp.dat')

      return

    end subroutine rewrite_dat

!#######################################################################

  end module rungamess
