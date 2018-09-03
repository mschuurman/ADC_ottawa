!######################################################################
! slepianmod: routines for the calculation of the coefficients in the
!             expansion of Slepian filter functions in terms of
!             Chebyshev polynomials of the first kind
!######################################################################

module slepianmod

contains
  
!######################################################################
  
  subroutine get_coeffs_slepians

    use cfdmod
    
    implicit none

!----------------------------------------------------------------------
! If the time-bandwidth product has not been supplied by the user,
! then look up the optimal value
!----------------------------------------------------------------------
  if (fw.eq.0.0d0.and..not.varfw) call get_optimal_fw
    
!----------------------------------------------------------------------
! Get the Slepian filter functions
!----------------------------------------------------------------------
    call get_slepians

!----------------------------------------------------------------------
! Output some information about the Slepians to the log file
!----------------------------------------------------------------------
    call wrslepinfo

!----------------------------------------------------------------------
! Calculate the coefficients entering into the expansion of the
! Slepian filter functions with respect to the Chebyshev polynomials
!----------------------------------------------------------------------
    call calc_expansion_coeffs
    
    return
    
  end subroutine get_coeffs_slepians

!######################################################################
  
  subroutine get_optimal_fw

    use channels
    use iomod
    use cfdmod
    
    implicit none

    character(len=4) :: ai

!----------------------------------------------------------------------
! Exit if the optimal time-bandwidth product is not available for the
! value of nfsbas
!----------------------------------------------------------------------
    if (nfsbas.gt.nopt) then
       write(ai,'(i4)') nopt
       errmsg='Optimal time-bandwidth products are not currently &
            available for numbers of Slepians greater than '&
            //trim(adjustl(ai))
       call error_control
    endif
    
!----------------------------------------------------------------------
! Set the optimal time-bandwidth product for the current no. Slepians
!----------------------------------------------------------------------
    call fill_optfw
    fw=optfw(nfsbas)
    
    return
    
  end subroutine get_optimal_fw
    
!######################################################################
  
  subroutine fill_optfw

    use channels
    use cfdmod
    
    implicit none

!----------------------------------------------------------------------
! Fill in the optfw array with the time-bandwidth products that yield
! 1-lambda_N values less than 5e-7
!----------------------------------------------------------------------
    optfw(1)=2.8d0
    optfw(2)=3.6d0
    optfw(3)=4.2d0
    optfw(4)=4.8d0
    optfw(5)=5.4d0
    optfw(6)=6.0d0
    optfw(7)=6.6d0
    optfw(8)=7.2d0
    optfw(9)=7.7d0
    optfw(10)=8.3d0
    optfw(11)=8.8d0
    optfw(12)=9.4d0
    optfw(13)=9.9d0
    optfw(14)=10.5d0
    optfw(15)=11.0d0
    optfw(16)=11.5d0
    optfw(17)=12.1d0
    optfw(18)=12.6d0
    optfw(19)=13.1d0
    optfw(20)=13.7d0
    optfw(21)=14.2d0
    optfw(22)=14.7d0
    optfw(23)=15.2d0
    optfw(24)=15.8d0
    optfw(25)=16.3d0
    optfw(26)=16.8d0
    optfw(27)=17.3d0
    optfw(28)=17.9d0
    optfw(29)=18.4d0
    optfw(30)=18.9d0
    optfw(31)=19.4d0
    optfw(32)=20.0d0
    optfw(33)=20.5d0
    optfw(34)=21.0d0
    optfw(35)=21.5d0
    optfw(36)=22.0d0
    optfw(37)=22.5d0
    optfw(38)=23.1d0
    optfw(39)=23.6d0
    optfw(40)=24.1d0
    optfw(41)=24.6d0
    optfw(42)=25.1d0
    optfw(43)=25.6d0
    optfw(44)=26.2d0
    optfw(45)=26.7d0
    optfw(46)=27.2d0
    optfw(47)=27.7d0
    optfw(48)=28.2d0
    optfw(49)=28.7d0
    optfw(50)=29.2d0
    
    return
    
  end subroutine fill_optfw
  
!######################################################################

  subroutine get_slepians

    use cfdmod
    
    implicit none

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Slepians
    allocate(v(npts,nfsbas))
    v=0.0d0

    ! Eigenvalues
    allocate(lambda(nfsbas))
    lambda=0.0d0

!----------------------------------------------------------------------
! Get the DPSSs
!----------------------------------------------------------------------
    if (varfw) then
       ! Variable time-bandwidth product for each Slepian filter: read
       ! in the pre-calculated DPSSs
       call read_slepians
    else
       ! Constant time-bandwidth product for all Slepian filters:
       ! calculate the DPSSs
       call calc_slepians
    endif

!----------------------------------------------------------------------
! Renormalisation on the interval [Eabar,Ebbar]
!----------------------------------------------------------------------
    v=v/sqrt(((Ebbar-Eabar)/(npts-1)))

    return
    
  end subroutine get_slepians

!######################################################################

  subroutine read_slepians

    use iomod
    use cfdmod
    
    implicit none

    integer            :: n,i,unit,itmp
    character(len=350) :: adcdir
    character(len=500) :: path
    character(len=550) :: filename
    character(len=4)   :: an
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Reading the DPSSs from file...'

!----------------------------------------------------------------------
! Exit if the no. DPSSs is greater than the no. pre-calculated
!----------------------------------------------------------------------
    if (nfsbas.gt.nprecalc) then
       errmsg='The no. Slepians requested is greater than the no. &
            that has been pre-calculated'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Read in the pre-calculated DPSSs
!----------------------------------------------------------------------
    ! Path to the DPSS files
    call get_environment_variable("ADC_DIR",adcdir)
    path=trim(adcdir)//'/source/analysis/chebyfd/dpss.var/'

    ! Next free unit
    call freeunit(unit)
    
    ! Loop over the Slepians
    do n=1,nfsbas

       ! Current filename
       write(an,'(i4)') n
       filename=trim(path)//'dpss.var.'//trim(adjustl(an))//'.dat'

       ! Open the DPSS file
       open(unit,file=filename,form='formatted',status='old')
       
       ! Read the DPSS file
       do i=1,npts
          read(unit,*) itmp,v(i,n)
       enddo
       
       ! Close the DPSS file
       close(unit)
       
    enddo

    return
    
  end subroutine read_slepians
    
!######################################################################

  subroutine calc_slepians

    use cfdmod
    use dpssmt

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the DPSSs...'

!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
    call dpss(npts,fw,nfsbas,v,lambda)
    
    return
    
  end subroutine calc_slepians

!######################################################################

  subroutine wrslepinfo

    use channels
    use cfdmod
    
    implicit none

    integer  :: i,j
    real(dp) :: ovrlp
    
!----------------------------------------------------------------------
! Output some information about the Slepians to file
!----------------------------------------------------------------------
    ! fw
    if (.not.varfw) &
         write(ilog,'(a,2x,F10.7)') '# Time half-bandwidth product:',fw
    
    ! Eigenvalues
    if (.not.varfw) then
       write(ilog,'(/,41a)') ('#',i=1,41)
       write(ilog,'(a)') '  DPSS       lambda         1 - lambda'
       write(ilog,'(41a)') ('#',i=1,41)
       do i=1,nfsbas
          write(ilog,'(2x,i3,2(2x,ES15.6))') i,lambda(i),1.0d0-lambda(i)
       enddo
    endif
       
    ! Overlaps on the interval [Eabar,Ebbar]
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(a)') '  DPSS Overlaps'
    write(ilog,'(41a)') ('#',i=1,41)
    do i=1,nfsbas
       do j=i,nfsbas
          ovrlp=dpss_overlap(i,j)
          write(ilog,'(2(2x,i3),2x,ES15.6)') i,j,ovrlp
       enddo
    enddo
    
    return
    
  end subroutine wrslepinfo

!######################################################################

  function dpss_overlap(i,j) result(ovrlp)

    use cfdmod
    
    implicit none

    integer  :: i,j,n
    real(dp) :: ovrlp

    ovrlp=v(1,i)*v(1,j)/2.0d0
    do n=2,npts-1
       ovrlp=ovrlp+v(n,i)*v(n,j)
    enddo
    ovrlp=ovrlp+v(npts,i)*v(npts,j)/2.0d0

    ovrlp=ovrlp*((Ebbar-Eabar)/(npts-1))
    
    return
    
  end function dpss_overlap

!######################################################################

  subroutine calc_expansion_coeffs

    use cfdmod
    use iomod
    
    implicit none

    integer                               :: n,j,k,unit
    real(dp)                              :: debar,ebar,theta
    real(dp), dimension(:,:), allocatable :: Tk,Tkw,val
    character(len=3)                      :: an
    character(len=60)                     :: filename
    logical                               :: exists
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Tk(0:Kdim,npts))
    Tk=0.0d0
    
    allocate(Tkw(0:Kdim,npts))
    Tkw=0.0d0

    allocate(val(npts,nfsbas))
    val=0.0d0
    
!----------------------------------------------------------------------
! Calculation of the expansion coefficients using the trapezoidal
! rule
!----------------------------------------------------------------------
    ! Precalculate the values of the Chebyshev polynomials at the
    ! quadrature points weighted by 1/sqrt(1-Ebar)
    debar=((Ebbar-Eabar)/(npts-1))
    do j=1,npts
       ebar=Eabar+debar*(j-1)
       theta=acos(ebar)
       do k=0,Kdim
          Tk(k,j)=cos(k*theta)
       enddo
       Tkw(:,j)=Tk(:,j)/sqrt(1.0d0-ebar**2)
    enddo
    
    ! Scale the first and last values of the weighted Chebyshev
    ! polynomials at the quadrature points s.t. we can take dot
    ! products to calculate the overlaps
    Tkw(:,1)=Tkw(:,1)/2.0d0
    Tkw(:,npts)=Tkw(:,npts)/2.0d0

    ! Calculate the expansion coefficients
    fkn=matmul(Tkw,v)
    
    ! Multiplication by DeltaE
    fkn=fkn*(Ebbar-Eabar)/npts

    ! Prefactors
    fkn(0,:)=fkn(0,:)/pi
    fkn(1:Kdim,:)=fkn(1:Kdim,:)*2.0d0/pi
    
!----------------------------------------------------------------------
! For checking purposes, output the Chebyshev expansions of the DPSSs
!----------------------------------------------------------------------
    ! Calculate the Chebyshev expansions of the DPSSs
    val=matmul(transpose(Tk),fkn)

    ! Output the Chebyshev expansions of the DPSSs
    inquire(file='dpss.cheby/.',exist=exists)
    if (exists) call system('rm dpss.cheby/*')
    if (.not.exists) call system('mkdir dpss.cheby')
    call freeunit(unit)
    do n=1,nfsbas
       write(an,'(i3)') n
       filename='dpss.cheby/dpss.cheby.'//trim(adjustl(an))//'.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       write(unit,'(40a)') ('#',k=1,40)
       write(unit,'(a)') '# Point   Approximate      Actual'
       write(unit,'(40a)') ('#',k=1,40)
       do j=1,npts
          write(unit,'(i5,2(2x,ES15.6))') j,val(j,n),v(j,n)
       enddo
       close(unit)
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Tk)
    deallocate(Tkw)
    deallocate(val)
    
    return
    
  end subroutine calc_expansion_coeffs
  
!######################################################################
  
end module slepianmod
