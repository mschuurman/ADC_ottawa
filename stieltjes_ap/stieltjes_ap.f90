  program stieltjes_ap

    use simod

    implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
    call rdstieltjesinp

!-----------------------------------------------------------------------
! Read the pseudostate energies and oscillator strengths
!-----------------------------------------------------------------------
    call rdosc

!-----------------------------------------------------------------------
! Perform the Stieltjes imaging calculation
!-----------------------------------------------------------------------
    call get_xsec

    STOP

  contains

!#######################################################################

    subroutine rdstieltjesinp

      use simod
      use parsemod

      implicit none

      integer            :: iin,i
      character(len=120) :: errmsg

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      asiinp=''
      aosc=''
      erange=-999.9d0
      ntrial=100

!-----------------------------------------------------------------------
! Determine the input file name
!-----------------------------------------------------------------------
      call getarg(1,asiinp)

      if (asiinp.eq.'') then
         write(6,'(/,a,/)') 'The input file name has not been given'
         STOP
      endif

!-----------------------------------------------------------------------
! Open the input file
!-----------------------------------------------------------------------
      iin=20
      open(iin,file=asiinp,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(iin)
        
      i=0
      if (keyword(1).ne.'end-input') then
10       continue
         i=i+1
         
         if (keyword(i).eq.'osc_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') aosc
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'interval') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) erange(1)
               i=i+1
               if (keyword(i).eq.',') then
                  i=i+1
                  read(keyword(i),*) erange(2)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'maxpol') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) ntrial
            else
               goto 100
            endif

         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            write(6,'(/,a,/)') trim(errmsg)
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
         write(6,'(/,a,/)') trim(errmsg)
         STOP

        endif

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
        if (aosc.eq.'') then
           errmsg='The oscillator strength filename has not been given'
           write(6,'(/,a,/)') trim(errmsg)
           STOP
        endif

        if (erange(1).eq.-999.9d0) then
           errmsg='The energy bounds have not been given'
           write(6,'(/,a,/)') trim(errmsg)
           STOP
        endif

        if (erange(2).eq.-999.9d0) then
           errmsg='The upper energy bound has not been given'
           write(6,'(/,a,/)') trim(errmsg)
           STOP
        endif

!-----------------------------------------------------------------------
! Close the input file
!-----------------------------------------------------------------------
        close(iin)

      return

    end subroutine rdstieltjesinp

!#######################################################################

    subroutine rdosc

      use simod
      use parsemod

      implicit none

      integer :: iosc,n
      real*8  :: ftmp

!-----------------------------------------------------------------------
! Open the oscillator strength file
!-----------------------------------------------------------------------
      iosc=20
      open(iosc,file=aosc,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the number of points and allocate the transition energy
! and oscillator strength arrays
!-----------------------------------------------------------------------
      npoints=0
5     call rdinp(iosc)      
      if (keyword(1).ne.'end-file') then
         read(keyword(1),*) ftmp
         if (ftmp.ge.erange(1).and.ftmp.le.erange(2)) npoints=npoints+1
         goto 5
      endif

      allocate(ener(npoints),osc(npoints))

!-----------------------------------------------------------------------
! Read the transition energies and oscillator strengths
!-----------------------------------------------------------------------
      rewind(iosc)

      n=0
10    call rdinp(iosc)      
      if (keyword(1).ne.'end-file') then
         read(keyword(1),*) ftmp
         if (ftmp.ge.erange(1).and.ftmp.le.erange(2)) then
            n=n+1
            ener(n)=ftmp
            read(keyword(2),*) osc(n)
         endif
         goto 10
      endif

!-----------------------------------------------------------------------
! Open the oscillator strength file
!-----------------------------------------------------------------------
      close(iosc)

      return

    end subroutine rdosc

!#######################################################################

    subroutine get_xsec

      use mpmodule
      use simod

      implicit none

      integer                     :: i
      real*8, dimension(0:ntrial) :: adbl,bdbl
      real*16                     :: tmp

!-----------------------------------------------------------------------
! MPFUN variables
!-----------------------------------------------------------------------
      ! Transition energies in arb. precision
      type(mp_real) e(npoints)

      ! Oscillator strengths in arb. precision
      type(mp_real) f(npoints)
      
      ! Orthogonal polynomials in arb. precision
      type(mp_real) Q(0:ntrial,npoints)
      
      ! Expansion coefficients in arb. precision
      type(mp_real) a(0:ntrial),b(0:ntrial)

!-----------------------------------------------------------------------
! Transform the transition energies and oscillator strengths to
! arbitrary precision
!-----------------------------------------------------------------------
      do i=1,npoints
         e(i)=mpreald(ener(i))
         f(i)=mpreald(osc(i))
      enddo

!-----------------------------------------------------------------------
! Calculate the coefficients entering into the three-term recursion
! relations for the polynomials orthogonal on the interval
! [erange(1),erange(2)] wrt to a weight function equal to the
! oscillator strength distribution
!-----------------------------------------------------------------------
      call Qrecursion(a,b,e,f,Q)

!-----------------------------------------------------------------------
! Calculate the cross-sections for successive orders up to nord
!
! From here on in, we convert the coefficients a and b to double
! precision
!-----------------------------------------------------------------------
      do i=0,ntrial
         tmp=qreal(a(i))
         adbl(i)=tmp
         tmp=qreal(b(i))
         bdbl(i)=tmp
      enddo

      call image_osc(adbl,bdbl)

      return

    end subroutine get_xsec

!#######################################################################
! Qrecursion: generates the polynomials orthogonal on the interval of 
!             interest wrt to the oscillator strength distribution.
!
! ntrial is the number of orthogonal polynomials Q_n that will be
! initially recursively generated.
!
! After checking the orthogonality of the Q_n, we will then use a
! smaller number nord < ntrial in the SI calculation.
!
! nord will be determined by the requirement that the first nord
! polynomials Q_n are orthogonal amongst themselves.
!
! Note that using a numerical recursion relation to generate the 
! Q_n will mean that there will exist an integer M s.t. for
! m or n greater than M <Q_m|Q_n> .ne. N_n delta_mn.
!#######################################################################

    subroutine Qrecursion(a,b,e,f,Q)

      use mpmodule
      use simod

      implicit none

      integer                    :: i,j,aunit,bunit
      real*8                     :: tmp
      character*1, dimension(30) :: string

      type(mp_real) asum,bprod,qnorm,qoverlap,thrsh,upper
      type(mp_real) a(0:ntrial),b(0:ntrial)
      type(mp_real) e(npoints),f(npoints)
      type(mp_real) Q(0:ntrial,npoints)
      type(mp_real) mpone,mpi,ainf,binf,mptmp

      write(6,'(/,2x,a)') 'Calculating the polynomials orthogonal &
           wrt the oscillator strength distribution...'

!-----------------------------------------------------------------------
! (1) Special cases
!-----------------------------------------------------------------------
      mpone='1.0'

      do i=0,ntrial
         a(i)='0.0'
         b(i)='0.0'
      enddo

      do i=1,npoints
         b(0)=b(0)+f(i)
         a(1)=a(1)+f(i)/e(i)
      enddo
      a(1)=a(1)/b(0)

      do i=1,npoints
         Q(0,i)='1.0'
         Q(1,i)=mpone/e(i)-a(1)
      enddo

      b(1)='0.0'
      a(2)='0.0'
      do i=1,npoints
         b(1)=b(1)+Q(1,i)*f(i)/e(i)
         a(2)=a(2)+Q(1,i)*f(i)/(e(i)**2)
      enddo
      b(1)=b(1)/b(0)
      a(2)=a(2)/(b(0)*b(1))-a(1)

!-----------------------------------------------------------------------
! (2) Remaining coefficients and polynomials
!-----------------------------------------------------------------------
      asum=a(1)
      do i=3,ntrial

!         print*,i,ntrial

         asum=asum+a(i-1)

         do j=1,npoints
            Q(i-1,j)=(mpone/e(j)-a(i-1))*Q(i-2,j)-b(i-2)*Q(i-3,j)
         enddo

         bprod=b(0)
         do j=1,i-2
            bprod=bprod*b(j)
         enddo

         b(i-1)='0.0'
         tmp=real(i-1)
         mpi=mpreald(tmp)
         do j=1,npoints
            b(i-1)=b(i-1)+Q(i-1,j)*f(j)/(e(j)**mpi)
         enddo
         b(i-1)=b(i-1)/bprod
         
         bprod=bprod*b(i-1)

         a(i)='0.0'
         tmp=real(i)
         mpi=mpreald(tmp)
         do j=1,npoints
            a(i)=a(i)+Q(i-1,j)*f(j)/(e(j)**mpi)
         enddo
         a(i)=a(i)/bprod-asum

      enddo

!-----------------------------------------------------------------------
! For the purposes of checking orthogonality, calculate the ntrial-th 
! order polynomial
!-----------------------------------------------------------------------
      do j=1,npoints
         Q(ntrial,j)=(mpone/e(j)-a(ntrial))*Q(ntrial-1,j)&
              -b(ntrial-1)*Q(ntrial-2,j)
      enddo

!-----------------------------------------------------------------------
! For checking purposes, output the differences between the
! coefficients and their limiting values
!
! We may monitor the a_n and b_n for divergence with increasing n in
! order to determine whether the maximum order used is acceptable
!
! lim n->inf a_n = 1/(2*E_thrsh)
!
! lim n->inf b_n = 1/(16*E_thrsh**2)
!
! where E_thrsh is the ionisation threshold.
!-----------------------------------------------------------------------
      aunit=20
      bunit=30
      open(aunit,file='a_ap.dat',form='formatted',status='unknown')
      open(bunit,file='b_ap.dat',form='formatted',status='unknown')

      ! Limiting values of a and b
      mptmp='2.0'
      ainf=mpone/(mptmp*e(1))
      mptmp='16.0'
      binf=mpone/(mptmp*e(1)*e(1))

      print*,
      call mpwrite(6,20,10,ainf)
      call mpwrite(6,20,10,binf)
      print*,

      do i=0,ntrial
         call mpwrite(aunit,20,10,a(i))
         call mpwrite(bunit,20,10,b(i))
      enddo

      close(aunit)
      close(bunit)

!-----------------------------------------------------------------------
! Determine the maximum approximation order by checking orthogonality
! of the polynomials Q_n
!-----------------------------------------------------------------------
      write(6,'(/,2x,a)') 'Checking orthogonality...'

      thrsh=mpreald(1.d-300)
      upper='10.0'

      nord=0
      ! Loop over polynomial orders
      do i=1,ntrial
         qnorm='0.0'
         qoverlap='0.0'
         ! Loop over data points
         do j=1,npoints
            ! Calculate the norm of the ith-order polynomial
            qnorm=qnorm+Q(i,j)*Q(i,j)*f(j)
            ! Calculate the overlap of the ith- and (i-1)th-order 
            ! polynomials
            qoverlap=qoverlap+Q(i,j)*Q(i-1,j)*f(j)
         enddo
         ! If we have lost orthogonality then set the maximum
         ! Stieltjes order (nord) to i-1 and exit
         if (abs(qoverlap).lt.thrsh) qoverlap=thrsh
         if (qnorm/abs(qoverlap).le.upper) then
            nord=i-1
            exit
         endif         
      enddo

      ! If nord is still zero, then all orders up to ntrial-1 are OK,
      ! which may well be the case when using arbitrary precision
      ! Note that we do not calculate the a and b coefficients, and thus
      ! can only go up to order ntrial-1
      if (nord.eq.0) nord=ntrial-1

      write(6,'(/,2x,a,1x,i3)') 'Maximum order:',nord

      return

    end subroutine Qrecursion

!#######################################################################
! image_osc: Performs the imaging of the oscillator strengths
!
! Note that in this subroutine the a and b coefficients are in double
! precision
!#######################################################################
    subroutine image_osc(a,b)

      use simod

      implicit none

      integer                      :: min,max,iord,i,j,iout,ierr,ioutF
      real*8, dimension(0:ntrial)  :: a,b
      real*8, dimension(nord,nord) :: abvec
      real*8, dimension(nord)      :: diag,offdiag
      real*8, dimension(nord)      :: ecent
      character(len=120)           :: outfile

      write(6,'(/,2x,a)') 'Calculating the cross-sections...'

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Gaussian quadrature points (Stieltjes energies)
      allocate(si_e(ntrial))

      ! Gaussian quadrature weights
      allocate(si_f(ntrial))

      ! Stieltjes oscillator strengths
      allocate(si_osc(ntrial))

      ! Stieltjes cumulative oscillator strengths:
      !  si_cosc1: histogram (see Eq. 3.5.1 of the Muller-Plathe and
      !            Diercksen paper)
      !  si_cosc2: LHS-RHS average at the midpoints (see Eq. 3.5.2 of
      !            the Muller-Plathe and Diercksen paper)
      allocate(si_cosc1(ntrial))
      allocate(si_cosc2(ntrial))

!-----------------------------------------------------------------------
! Set the minimum and maximum orders
!-----------------------------------------------------------------------
      if (nord.lt.5) then
         write(6,'(/,2x,a)') 'Only a very low-order approximation &
              is available'
         write(6,'(2x,a,1x,i2)')'Maximum order = ',nord
         min=nord
         max=nord
      else
         min=5
         max=nord
      endif

!-----------------------------------------------------------------------
! Calculate the photoionisation cross-sections using successive
! approximation orders from min to max
!-----------------------------------------------------------------------
        iout=20
        ioutF=30

        do iord=min,max

           ! Open the output files
           if (iord.lt.10) then
              write(outfile,'(a12,i1)') 'xsec_order00',iord
           else if (iord.lt.100) then
              write(outfile,'(a11,i2)') 'xsec_order0',iord
           else
              write(outfile,'(a10,i3)') 'xsec_order',iord
           endif
           open(iout,file=outfile,form='formatted',status='unknown')
           if (iord.lt.10) then
              write(outfile,'(a9,i1)') 'F_order00',iord
           else if (iord.lt.100) then
              write(outfile,'(a8,i2)') 'F_order0',iord
           else
              write(outfile,'(a7,i3)') 'F_order',iord
           endif
           open(ioutF,file=outfile,form='formatted',status='unknown')

           ! Construct the tridiagonal recursion coefficient matrix
           !
           ! (1) On-diagonal elements
           do i=1,iord
              diag(i)=a(i)
           enddo
           ! (2) Off-diagonal elements
           do i=2,iord
              offdiag(i)=-sqrt(b(i-1))
           enddo

           ! Diagonalise the tridiagonal recursion coefficient matrix
           abvec=0.0d0
           do i=1,nord
              abvec(i,i)=1.0d0
           enddo
           call tql2(nord,iord,diag,offdiag,abvec,ierr)
           if (ierr.ne.0) then              
              write(6,'(/,2x,a,/)') 'Diagonalisation of the &
                   recurrence coefficient matrix failed'
              STOP
           endif

           ! Calculate the Stieltjes energies and cumulative
           ! oscilator strength values
           ! 
           ! Note that the eigenvalues of the recurrence coefficient 
           ! matrix are the inverse energies in ascending order
           do i=1,iord
              si_e(i)=1.0d0/diag(iord+1-i)
              si_f(i)=b(0)*abvec(1,iord+1-i)**2
          enddo
          
          ! Calculate and output the oscillator strengths using 
          ! numerical differentiation of the Stieltjes distribution 
          ! function F
          do i=1,iord-1
             ecent(i)=(si_e(i)+si_e(i+1))/2.0d0
             si_osc(i)=(si_f(i+1)+si_f(i))/(2.d0*(si_e(i+1)-si_e(i)))
             write(iout,*) ecent(i)*27.211d0,si_osc(i)
          enddo

          ! Calculate the Stieltjes distribution function, i.e., the
          ! cumulative oscillator strength distribution
          si_cosc1=0.0d0
          do i=1,iord
             do j=1,i
                si_cosc1(i)=si_cosc1(i)+si_f(j)
             enddo
          enddo          
          si_cosc2=0.0d0
          do i=1,iord-1
             si_cosc2(i)=(si_cosc1(i)+si_cosc1(i+1))/2.0d0
             write(ioutF,*) ecent(i)*27.211d0,si_cosc2(i)
          enddo

          ! Close the output files
          close(iout)
          close (ioutF)

       enddo

      return

    end subroutine image_osc

!#######################################################################

  end program stieltjes_ap
