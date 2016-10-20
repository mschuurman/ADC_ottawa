!#######################################################################
! rixsplt: a simple program to plot the RIXS spectrum using transition
!          matrix elements and energies from a block-Lanczos ADC
!          calculation
!#######################################################################
  program rixsplt

    use constants
    use iomod

    implicit none

!-----------------------------------------------------------------------
! Read the input
!-----------------------------------------------------------------------
    call rdrixsinp

!-----------------------------------------------------------------------
! Read the RIXS data file
!-----------------------------------------------------------------------
    call rddatfile

!-----------------------------------------------------------------------
! Calculation of the transition matrix elements
! <f| D (Ebar_i(E_incident) - H)^-1 D |i>
!-----------------------------------------------------------------------
    call calc_transmat

!-----------------------------------------------------------------------
! Calculate the RIXS spectrum
!-----------------------------------------------------------------------
    call calc_rixs

  contains

!#######################################################################

    subroutine rdrixsinp

      use channels
      use constants
      use iomod
      use parsemod
      use rixsmod

      implicit none

      integer :: i,k

!-----------------------------------------------------------------------
! Initialise parameters
!-----------------------------------------------------------------------
      gammaint=0.0d0
      gammaf=0.0d0
      nener1=0
      nener2=0
      ener1=0.0d0
      ener2=0.0d0
      de1=0.0d0
      de2=0.0d0
      datfile='rixs.dat'
      istate=0

!-----------------------------------------------------------------------
! Read the name of the input file from the command line
!-----------------------------------------------------------------------
      ain=''

      call getarg(1,ain)

      if (ain.eq.'') then
         write(6,'(/,2x,a,/)') 'The name of the input file has not been &
              given'
         STOP
      endif

!-----------------------------------------------------------------------
! Open the input and log files
!-----------------------------------------------------------------------
      ! File names
      k=index(ain,'.inp')
      if (k.eq.0) then
         alog=trim(ain)//'.log'
         ain=trim(ain)//'.inp'
      else
         alog=ain(1:k-1)//'.log'
      endif
      
      ! Open the input file
      call freeunit(iin)
      open(iin,file=ain,form='formatted',status='old')

      ! Open the log file      
      call freeunit(ilog)
      open(ilog,file=alog,form='formatted',status='unknown')     

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(iin)

      i=0
      if (.not.lend) then
10       continue
         i=i+1

         if (keyword(i).eq.'gamma_int') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) gammaint
            else
               goto 100
            endif

         else if (keyword(i).eq.'gamma_final') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) gammaf
            else
               goto 100
            endif

         else if (keyword(i).eq.'energy1') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Lower bound
               read(keyword(i),*) ener1(1)
               ! Upper bound
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) ener1(2)
               else
                  errmsg='The incident photon energy upper bound &
                       has not been given'
                  call error_control
               endif
               ! Number of grid points
                if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) nener1
               else
                  errmsg='The number of incident photon energies &
                       has not been given'
                  call error_control
               endif
               ! Grid spacing
               de1=(ener1(2)-ener1(1))/nener1
            else
               goto 100
            endif

         else if (keyword(i).eq.'energy2') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Lower bound
               read(keyword(i),*) ener2(1)
               ! Upper bound
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) ener2(2)
               else
                  errmsg='The emitted photon energy upper bound &
                       has not been given'
                  call error_control
               endif
               ! Number of grid points
                if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) nener2
               else
                  errmsg='The number of emitted photon energies &
                       has not been given'
                  call error_control
               endif
               ! Grid spacing
               de2=(ener2(2)-ener2(1))/nener1
            else
               goto 100
            endif

         else if (keyword(i).eq.'datfile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               datfile=keyword(i)
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'istate') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) istate
               ! Note that, somewhat confusingly, in the input 
               ! file, state indices start from 0 but in the
               ! code, they start from 1
               istate=istate+1
            else
               goto 100
            endif

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
! Close the input file
!-----------------------------------------------------------------------
      close(iin)

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
      if (gammaint.eq.0.0d0) then
         errmsg='The intermediate state broadening has not been given'
         call error_control
      endif

      if (gammaf.eq.0.0d0) then
         errmsg='The final state broadening has not been given'
         call error_control
      endif

      if (de1.eq.0.0d0) then
         errmsg='The incident energy grid has not been given'
         call error_control
      endif
      
      if (de2.eq.0.0d0) then
         errmsg='The emission energy grid has not been given'
         call error_control
      endif

      if (istate.eq.0) then
         errmsg='The initial state has not been given'
         call error_control
      endif

!-----------------------------------------------------------------------
! Write the spectrum parameters to the log file
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,F7.4,/)') 'Intermediate state broadening (eV):',&
           gammaint
      
      write(ilog,'(2x,a,7x,F7.4,/)') 'Final state broadening (eV):',&
           gammaf
           
      write(ilog,'(2x,a,9x,3(x,F8.4),/)') 'Incident energy grid (eV):',&
           ener1(1),ener1(2),de1

      write(ilog,'(2x,a,9x,3(x,F8.4),/)') 'Emission energy grid (eV):',&
           ener2(1),ener2(2),de2

      write(ilog,'(2x,a,21x,a,/)') 'RIXS data file:',&
           trim(datfile)

      return

    end subroutine rdrixsinp

!#######################################################################

    subroutine rddatfile

      use constants
      use iomod
      use rixsmod

      implicit none

      integer :: unit,i,j

!-----------------------------------------------------------------------
! Read the RIXS data file
!-----------------------------------------------------------------------
      ! Open the RIXS data file
      call freeunit(unit)
      open(unit,file=datfile,status='old',form='unformatted')

      ! Dimensions
      read(unit) nval
      read(unit) nlanc

      ! Allocate arrays
      allocate(tdm(nval,nlanc))
      allocate(enerval(nval))
      allocate(enerlanc(nlanc))

      ! Valence state energies
      do i=1,nval
         read(unit) enerval(i)
      enddo
      
      ! Lanczos pseudo-state energies
      do i=1,nlanc
         read(unit) enerlanc(i)
      enddo

      ! Matrix elements <Psi_i | D | chi_j> between the valence
      ! states |Psi_i> (ground + excited) and the Lanczos pseudo-states
      ! |chi_j> from file
      do i=1,nval
         do j=1,nlanc
            read(unit) tdm(i,j)
         enddo
      enddo

      ! Close the RIXS data file
      close(unit)

      !! TEST
      !do i=1,nlanc
      !   print*,(enerlanc(i)-enerval(1))*27.2113845d0,&
      !        (2.0/3.0)*(enerlanc(i)-enerval(1))*tdm(1,i)**2
      !enddo
      !STOP

      return

    end subroutine rddatfile

!#######################################################################

     subroutine calc_transmat

      use constants
      use iomod
      use rixsmod

      implicit none

      integer                            :: i,j,k,f,e,l
      real(d)                            :: einc,gamma
      real(d), parameter                 :: eh2ev=27.2113845d0
      complex*16, dimension(nval,nener1) :: trans
      complex*16                         :: ctmp

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(transsq(nval,nener1))
      trans=0.0d0

!-----------------------------------------------------------------------
! Calculate the terms
!
! |<f| D (Ebar_i(E_incident) - H)^-1 D |i>|^2
!
! for all final states |f> and all incident photon energies
!
! N.B. All quantities here are in a.u.
!-----------------------------------------------------------------------
      ! Intermediate state broadening in a.u.
      gamma=gammaint/eh2ev

      ! Loop over incident photon energies
      do k=1,nener1

         ! Current incident photon energy (a.u.)
         einc=(ener1(1)+(i-1)*de1)/eh2ev

         ! Loop over final valence states (ground + excited)
         do f=1,nval

            !
            ! <f| D (Ebar_i(E_incident) - H)^-1 D |i>
            !

            ! Loop over Lanczos pseudostates
            do l=1,nlanc
               ctmp=(tdm(f,l)*tdm(istate,l))
               ctmp=ctmp/(enerval(istate) - enerlanc(l) + einc - ci*gamma/2)
               trans(f,k)=trans(f,k)+ctmp
            enddo

            !
            ! |<f| D (Ebar_i(E_incident) - H)^-1 D |i>|^2
            !
            transsq(f,k)=real(trans(f,k)*conjg(trans(f,k)))

         enddo
      enddo

      return

    end subroutine calc_transmat

!#######################################################################

    subroutine calc_rixs

      use iomod
      use constants
      use rixsmod

      implicit none

      integer            :: i,j,f,unit
      real(d)            :: einc,eemit,lineshape,func
      real(d), parameter :: eh2ev=27.2113845d0

!-----------------------------------------------------------------------
! Open the RIXS spectrum file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file='rixsspec.dat',form='formatted',status='unknown')      

!-----------------------------------------------------------------------
! Calculate and output the RIXS spectrum
!-----------------------------------------------------------------------
      ! Loop over incident photon energies
      do i=1,nener1
         
         write(unit,'(a)') ' '

         ! Current incident photon energy (eV)
         einc=ener1(1)+(i-1)*de1
         
         ! Loop over emitted photon energies
         do j=1,nener2

            ! Current emitted photon energy (eV)
            eemit=ener2(1)+(j-1)*de2

            ! Loop over final states
            func=0.0d0
            do f=1,nval
               ! Lineshape
               lineshape=lorentzian(enerval(istate)*eh2ev,&
                    enerval(f)*eh2ev,einc,eemit)
               
               ! Contribution to the function value
               func=func+transsq(f,i)*lineshape
               
            enddo

            func=func*(eemit/einc)

            ! Write the current point to file
            write(unit,*) einc,eemit,func

         enddo

      enddo

!-----------------------------------------------------------------------
! Close the RIXS spectrum file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine calc_rixs

!#######################################################################

    function lorentzian(ei,ef,einc,eemit) result(func)

      use constants
      use rixsmod

      implicit none

      real(d) :: ei,ef,einc,eemit,func,numer,denom
      
      numer=gammaf/(2.0d0*pi)

      denom=(ei-ef+einc-eemit)**2 + 0.25d0*gammaf**2

      func=numer/denom

      return

    end function lorentzian

!#######################################################################

  end program rixsplt
