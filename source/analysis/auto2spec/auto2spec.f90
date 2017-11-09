module auto2specmod

  use constants

  save

  integer                               :: maxtp,iauto,epoints
  integer                               :: padesolver
  real(d)                               :: dt,t0,emin,emax,tau,sigma,&
                                           dele,a0,b0,tcutoff
  real(d), dimension(:,:), allocatable  :: sp
  real(d), parameter                    :: eh2ev=27.2113845d0
  real(d), parameter                    :: fs2au=41.3413745758d0
  complex(d), dimension(:), allocatable :: auto,avec,bvec
  character(len=70)                     :: outfile
  logical                               :: lpade,lnormalise
  
end module auto2specmod

!######################################################################

program auto2spec

  use auto2specmod
  
  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdauto2specinp
  
!----------------------------------------------------------------------
! Read the autocorrelation function from file
!----------------------------------------------------------------------
  call rdauto
  
!----------------------------------------------------------------------
! Calculate the spectrum
!----------------------------------------------------------------------
  call calc_spectrum

!----------------------------------------------------------------------
! Write the spectrum to file
!----------------------------------------------------------------------
  call wrspectrum
  
contains

!######################################################################

  subroutine rdauto2specinp

    use auto2specmod
    use constants
    use iomod
    
    implicit none

    integer           :: i
    character(len=30) :: string1,string2
    
!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
    ! Energy range
    emin=-1.0d0
    emax=-1.0d0

    ! No. energy points
    epoints=1000

    ! Damping time
    tau=-1.0d0
    
    ! Name of the output file
    outfile=''

    ! Pade approximant of the Fourier transform
    lpade=.false.

    ! Method for calculating the Pade approximant coefficients
    padesolver=1
    
    ! Timestep cutofff (in fs)
    tcutoff=1e+10_d

    ! Normalisation of the spectra
    lnormalise=.false.
    
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    i=1
5   call getarg(i,string1)

    if (string1.eq.'-sigma') then
       ! FWHM in eV
       i=i+1
       call getarg(i,string2)
       read(string2,*) sigma
       ! Convert to tau in au
       sigma=sigma/eh2ev
       tau=2.0d0/sigma

    else if (string1.eq.'-e') then
       ! Energy range in eV
       i=i+1
       call getarg(i,string2)
       read(string2,*) emin
       i=i+1
       call getarg(i,string2)
       read(string2,*) emax
       ! Convert to au
       emin=emin/eh2ev
       emax=emax/eh2ev

    else if (string1.eq.'-np') then
       ! No. energy points
       i=i+1
       call getarg(i,string2)
       read(string2,*) epoints

    else if (string1.eq.'-o') then
       ! Name of the output file
       i=i+1
       call getarg(i,outfile)

    else if (string1.eq.'-pade') then
       ! Pade approximant of the Fourier transform
       lpade=.true.
       ! Calculation of the Pade approximant coefficients using
       ! matrix inversion
       padesolver=1
       
    else if (string1.eq.'-pade_lu') then
       ! Pade approximant of the Fourier transform
       lpade=.true.
       ! Calculation of the Pade approximant coefficients using
       ! LU factorisation
       padesolver=2
       
    else if (string1.eq.'-tcut') then
       ! Time cutoff
       i=i+1
       call getarg(i,string2)
       read(string2,*) tcutoff
!       tcutoff=tcutoff*fs2au
       
    else if (string1.eq.'-norm') then
       ! Normalisation of the spectra
       lnormalise=.true.

    else
       errmsg='Unknown keyword: '//trim(string1)
       call error_control
    endif
       
    if (i.lt.iargc()) then
       i=i+1
       goto 5
    endif

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (emin.eq.-1.0d0) then
       errmsg='The energy range has not been given'
       call error_control
    endif

    if (tau.eq.-1.0d0) then
       errmsg='The spectral broadening has not been given'
       call error_control
    endif
    
    if (outfile.eq.'') then
       outfile='spectrum.dat'
    endif

    return
    
  end subroutine rdauto2specinp
    
!######################################################################

  subroutine rdauto

    use auto2specmod
    use constants
    use iomod
    
    implicit none

    integer :: iin,i
    real(d) :: t,re,im
    
!----------------------------------------------------------------------
! Open the autocorrelation function file
!----------------------------------------------------------------------
    call freeunit(iin)
    open(iin,file='auto',form='formatted',status='old')

!----------------------------------------------------------------------
! Determine the no. timesteps and the time interval, and allocate
! arrays
!----------------------------------------------------------------------
    ! Read past the comment line
    read(iin,*)

    ! Determine the timestep
    read(iin,*) t0
    read(iin,*) t
    t0=t0*fs2au
    t=t*fs2au
    dt=t-t0

    ! Determine the number of timesteps
    maxtp=2
5   read(iin,*,end=10) t
    if (t.lt.tcutoff) then
       maxtp=maxtp+1
       goto 5
    endif
10  continue
    iauto=maxtp-1

    ! Allocate arrays
    allocate(auto(maxtp))
    auto=czero

    ! If we are calculating the pade approximant of the Fourier
    ! transform, then make sure that iauto is even.
    if (lpade.and.mod(iauto,2).ne.0) iauto=iauto-1
    
!----------------------------------------------------------------------
! Read the autocorrelation function
!----------------------------------------------------------------------
    rewind(iin)
    read(iin,*)
    read(iin,*) t,a0,b0

    do i=1,iauto
       read(iin,*) t,re,im
       auto(i)=dcmplx(re,im)
    enddo
    
!----------------------------------------------------------------------
! Close the autocorrelation function file
!----------------------------------------------------------------------
    close(iin)
    
    return

  end subroutine rdauto

!######################################################################

  subroutine calc_spectrum

    use auto2specmod
    
    implicit none

    integer :: i,j
    real(d) :: eau,t,cc,sum0,sum1,sum2,sum3,pia,gfac

!----------------------------------------------------------------------
! Allocate the spectrum arrays
!----------------------------------------------------------------------
    allocate(sp(0:epoints,0:4))
    
!----------------------------------------------------------------------
! Calculate the Pade approximant coefficients
!----------------------------------------------------------------------    
    if (lpade) call pade_coeff
    
!----------------------------------------------------------------------
! Calculate the spectrum
!----------------------------------------------------------------------
    dele=(emax-emin)/dble(epoints)
    pia=pi/dble(2*maxtp+2)    
    gfac=dble(2.24d0/(maxtp+1))

    ! Loop over energy points
    do j=0,epoints

       ! Current energy
       eau=emin+j*dele

       ! Initialisation of the finite Fourier transforms
       sum0=0.5d0*a0
       sum1=0.5d0*a0
       sum2=0.5d0*a0
       sum3=0.5d0*a0

       ! Loop over timesteps
       do i=1,iauto

          t=i*dt

          cc=dble(exp(dcmplx(0.0d0,eau*t))*auto(i))*exp(-(t/tau))

          sum0=sum0+cc
          sum1=sum1+cc*cos(pia*i)
          sum2=sum2+cc*cos(pia*i)**2
          sum3=sum3+cc*exp(-(gfac*i)**2)

       enddo

       ! Windowed Fourier transforms
       sp(j,0)=eau*sum0*dt/pi
       sp(j,1)=eau*sum1*dt/pi
       sp(j,2)=eau*sum2*dt/pi          
       sp(j,3)=eau*sum3*dt/pi
       
       ! Pade approximant
       if (lpade) sp(j,4)=eau*padespec(eau)

    enddo

    return
    
  end subroutine calc_spectrum

!######################################################################
! pade_coeff: the notation used here is the same as that used in
!             Equations 29-35 in Bruner et. al., JCTC, 12, 3741 (2016)
!######################################################################
  
  subroutine pade_coeff

    use auto2specmod
    
    implicit none

    integer                                 :: k,m,n
    complex(d), dimension(:), allocatable   :: cvec,dvec
    complex(d), dimension(:,:), allocatable :: gmat

!----------------------------------------------------------------------
! Allocate and initialisae arrays
!----------------------------------------------------------------------
    n=iauto/2

    allocate(avec(0:n))
    avec=0.0d0
    allocate(bvec(0:n))
    bvec=0.0d0
    allocate(cvec(0:iauto))
    cvec=0.0d0
    allocate(dvec(n))
    dvec=0.0d0
    allocate(gmat(n,n))
    gmat=0.0d0
    
!----------------------------------------------------------------------
! c-coefficients
!----------------------------------------------------------------------
    cvec(0)=dcmplx(a0,b0)

    do k=1,iauto
       cvec(k)=auto(k)*exp(-(k*dt)/tau)
    enddo
       
!----------------------------------------------------------------------
! Construct the d-vector
!----------------------------------------------------------------------
    do k=1,n
       dvec(k)=-cvec(n+k)
    enddo

!----------------------------------------------------------------------
! Construct the G-matrix
!----------------------------------------------------------------------
    do k=1,n
       do m=1,n
          gmat(k,m)=cvec(n-m+k)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the b-coefficients by colving the system of linear
! equations G.b = d
! N.B., we use the usual convention for diagonal Pade approximant
! schemes and set b0=1
!----------------------------------------------------------------------
    call calc_bcoeff(gmat,dvec,n)

!----------------------------------------------------------------------
! Calculate the a-coefficients
!----------------------------------------------------------------------
    avec(0)=cvec(0)
    
    do k=1,n
       do m=0,k
          avec(k)=avec(k)+bvec(m)*cvec(k-m)
       enddo
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(cvec,dvec,gmat)
    
    return
    
  end subroutine pade_coeff

!######################################################################

  subroutine calc_bcoeff(gmat,dvec,n)

    use constants
    
    implicit none

    integer                                 :: n
    complex(d), dimension(n,n)              :: gmat
    complex(d), dimension(n)                :: dvec
    complex(d), dimension(:,:), allocatable :: invgmat

    select case(padesolver)

    case(1) ! Solution of G b = d via the inversion of G
       call calc_bcoeff_inv(gmat,dvec,n)
       
    case(2) ! Solution of G b = d using the LU factorisation of G
       call calc_bcoeff_lu(gmat,dvec,n)

    end select
    
    return

  end subroutine calc_bcoeff

!######################################################################

  subroutine calc_bcoeff_inv(gmat,dvec,n)

    use constants
    
    implicit none

    integer                                 :: n
    complex(d), dimension(n,n)              :: gmat
    complex(d), dimension(n)                :: dvec
    complex(d), dimension(:,:), allocatable :: invgmat

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(invgmat(n,n))
    invgmat=0.0d0

!----------------------------------------------------------------------
! Calculate the pseudoinverse of G
!----------------------------------------------------------------------
    call pseudoinverse_cmplx(gmat,invgmat,n)

!----------------------------------------------------------------------
! Calculate the b-coefficients
!----------------------------------------------------------------------
    bvec(0)=cone

    bvec(1:n)=matmul(invgmat,dvec)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(invgmat)
    
    return
    
  end subroutine calc_bcoeff_inv

!######################################################################

   subroutine calc_bcoeff_lu(gmat,dvec,n)

    use constants
    use iomod
    
    implicit none

    integer                            :: n,info
    integer, dimension(:), allocatable :: ipiv
    complex(d), dimension(n,n)         :: gmat
    complex(d), dimension(n)           :: dvec

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ipiv(n))
    
!----------------------------------------------------------------------
! LU factorisation of G
!----------------------------------------------------------------------
    call zgetrf(n,n,gmat,n,ipiv,info)

    if (info.ne.0) then
       write(6,'(/,2x,a,/)') 'LU decomposition of the G-matrix &
            failed in subroutine calc_bcoeff_lu'
       stop
    endif

!----------------------------------------------------------------------
! Calculate the b-coefficients
!----------------------------------------------------------------------
    bvec(0)=cone

    bvec(1:n)=dvec(1:n)
    call zgetrs('N',n,1,gmat,n,ipiv,bvec(1:n),n,info)

    if (info.ne.0) then
       write(6,'(/,2x,a,/)') 'Failed call to zgetrs in subroutine &
            subroutine calc_bcoeff_lu'
       stop
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ipiv)
    
    return
    
  end subroutine calc_bcoeff_lu
    
!######################################################################

  subroutine pseudoinverse_cmplx(mat,invmat,n)

    use constants
    use iomod
    
    implicit none

    integer                    :: n,lwork,info,i,j
    real(d), dimension(n)      :: sigma
    real(d), dimension(5*n)    :: rwork
    real(d), parameter         :: thrsh=1e-10_d
    complex(d), dimension(n,n) :: mat,invmat,u,vt,tmpmat,smat
    complex(d), dimension(3*n) :: work
    
!-----------------------------------------------------------------------
! SVD of the matrix to be inverted
!-----------------------------------------------------------------------
    lwork=3*n
    tmpmat=mat
    call zgesvd('A','A',n,n,tmpmat,n,sigma,u,n,vt,n,work,lwork,rwork,info)

    if (info.ne.0) then
       write(6,'(/,2x,a)') &
            'SVD failure in subroutine pseudoinverse_cmplx'
       stop
    endif

!-----------------------------------------------------------------------
! Pseudo-inverse
!-----------------------------------------------------------------------
    smat=0.0d0
    do i=1,n
       if (abs(sigma(i)).lt.thrsh) then
          smat(i,i)=0.0d0
       else
          smat(i,i)=1.0d0/sigma(i)
       endif
    enddo

    vt=conjg(vt)
    vt=transpose(vt)
    u=conjg(u)
    u=transpose(u)
    
    invmat=matmul(vt,matmul(smat,u))

    return

  end subroutine pseudoinverse_cmplx

!######################################################################

  function padespec(e) result(func)

    use auto2specmod
    
    implicit none

    integer    :: k
    real(d)    :: e,func
    complex(d) :: z,numer,denom
    
    z=exp(ci*e*dt)
    
    numer=czero
    denom=czero
    do k=0,iauto/2
       numer=numer+avec(k)*(z**k)
       denom=denom+bvec(k)*(z**k)      
    enddo
    
    func=2.0d0*real(numer/denom)

    return
    
  end function padespec
    
!######################################################################
  
  subroutine wrspectrum

    use auto2specmod
    use constants
    use iomod
    
    implicit none

    integer :: i,iout
    real(d) :: e,sup

!----------------------------------------------------------------------
! Open the spectrum file
!----------------------------------------------------------------------
    call freeunit(iout)
    open(iout,file=outfile,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Normalise the spectra
!----------------------------------------------------------------------
    if (lnormalise) then
       do i=0,3
          sp(:,i)=sp(:,i)/maxval(sp(:,i))
       enddo
       if (lpade) sp(:,4)=sp(:,4)/maxval(sp(:,4))
    endif

!----------------------------------------------------------------------
! Write the spectra to file
!----------------------------------------------------------------------
    write(iout,'(a,x,F7.4)') '# FHWM (eV):',sigma*eh2ev

    if (lpade) then
       write(iout,'(/,100a)') ('#',i=1,100)
       write(iout,'(a)') '#  Energy (eV)      Pade             &
            cos^2 filter     cos filter       Gaussian filter  &
            no filter'
       write(iout,'(100a)') ('#',i=1,100)
       do i=0,epoints
          e=(emin+i*dele)*eh2ev
          write(iout,'(6(ES15.6,2x))') e,abs(sp(i,4)),sp(i,2),&
               sp(i,1),sp(i,3),sp(i,0)
       enddo
    else
       write(iout,'(/,84a)') ('#',i=1,84)
       write(iout,'(a)') '#  Energy (eV)      cos^2 filter     &
            cos filter       Gaussian filter  no filter'
       write(iout,'(84a)') ('#',i=1,84)
       
       do i=0,epoints
          e=(emin+i*dele)*eh2ev
          write(iout,'(5(ES15.6,2x))') e,sp(i,2),sp(i,1),sp(i,3),&
               sp(i,0)
       enddo
    endif
       
!----------------------------------------------------------------------
! Close the spectrum file
!----------------------------------------------------------------------
    close(iout)
    
    return
    
  end subroutine wrspectrum

!######################################################################
  
end program auto2spec
