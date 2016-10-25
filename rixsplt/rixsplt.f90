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
! Calculate the zeta tensor
!-----------------------------------------------------------------------
    call calc_zeta
    
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
      theta=-999.0d0

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

         else if (keyword(i).eq.'theta') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) theta
               ! Convert to radians
               theta=theta*pi/180.0d0
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

      if (theta.eq.-999.0d0) then
         errmsg='The detection angle has not been given'
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

      write(ilog,'(2x,a,13x,F6.2,/)') 'Detection angle (deg):',theta/pi*180.0d0
      
      write(ilog,'(2x,a,21x,a,/)') 'RIXS data file:',trim(datfile)

      return

    end subroutine rdrixsinp

!#######################################################################

    subroutine rddatfile

      use channels
      use constants
      use iomod
      use rixsmod

      implicit none

      integer :: unit,i,j,c

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
      allocate(tdm(3,nval,nlanc))
      tdm=0.0d0
      allocate(enerval(nval))
      enerval=0.0d0
      allocate(enerlanc(nlanc))
      enerlanc=0.0d0
      
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
         do c=1,3
            do j=1,nlanc
               read(unit) tdm(c,i,j)
            enddo
         enddo
      enddo
         
      ! Close the RIXS data file
      close(unit)

      ! Write the number of states to the log file
      write(ilog,'(/,2x,a,x,i2)') 'Number of valence states:',nval
      write(ilog,'(/,2x,a,x,i5)') 'Number of Lanczos states:',nlanc
      
      return

    end subroutine rddatfile

!#######################################################################

    subroutine calc_zeta

      use constants
      use rixsmod

      implicit none

      integer :: alpha,beta,f
      real(d) :: term1,term2
      
!-----------------------------------------------------------------------
! Calculation of the zeta tensor defined in Equation 13 of
! Phys. Rev. A, 49, 4378 (1994)      
!-----------------------------------------------------------------------
      ! Allocate arrays
      allocate(zeta(nval,nlanc,nlanc))
      zeta=0.0d0

      ! Loop over valence states
      do f=1,nval

         ! Loop over pairs of Lanczos states
         do alpha=1,nlanc
            do beta=1,nlanc

               term1 = dot_product(tdm(:,istate,alpha),tdm(:,istate,beta)) &
                    * dot_product(tdm(:,f,alpha),tdm(:,f,beta)) &
                    * (2.0d0 - cos(theta)**2)

               term2 = dot_product(tdm(:,istate,alpha),tdm(:,f,alpha)) &
                    * dot_product(tdm(:,istate,beta),tdm(:,f,beta))
               term2 = term2 + dot_product(tdm(:,istate,alpha),tdm(:,f,beta)) &
                    * dot_product(tdm(:,istate,beta),tdm(:,f,alpha))
               term2 = term2 * (3.0d0*cos(theta)**2 -1)/2.0d0

               zeta(f,alpha,beta) = term1 + term2
               
            enddo
         enddo

      enddo

      ! Prefactor
      zeta=zeta/15.0d0
      
      return
      
    end subroutine calc_zeta
    
!#######################################################################

    subroutine calc_rixs

      use iomod
      use constants
      use rixsmod

      implicit none

      integer            :: i,j,f,unit,alpha,beta
      real(d)            :: einc,eemit,lineshape,func,gamma
      real(d), parameter :: eh2ev=27.2113845d0
      complex*16         :: denom

!-----------------------------------------------------------------------
! Open the RIXS spectrum file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file='rixsspec.dat',form='formatted',status='unknown')      

!-----------------------------------------------------------------------
! Calculate and output the RIXS spectrum
!
! WE NEED TO CHANGE THIS SO THAT WE FIRST CONTRACT OVER THE LANCZOS
! STATES TO FORM AN INTERMEDIATE CORRESPONDING TO THE FINAL STATE
! INDEX, AND THEN CONTRACT THIS INTERMEDIATE WITH THE LINESHAPE
!-----------------------------------------------------------------------
      ! Intermediate state lifetime broadening in a.u.
      gamma=gammaint/eh2ev

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

               ! Loop over pairs of Lanczos states
               do alpha=1,nlanc
                  do beta=1,nlanc

                     denom=enerval(istate) - enerlanc(alpha) &
                          + einc/eh2ev - ci*gamma/2.0d0
                     denom=denom * (enerval(istate) - enerlanc(beta) &
                          + einc/eh2ev + ci*gamma/2.0d0)

                     ! WHY IS THE DENOMONATOR COMPLEX?
                     print*,imag(denom)
                     
                     func=func+lineshape*zeta(f,alpha,beta)/real(denom)
                     
                  enddo
               enddo

            enddo

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
