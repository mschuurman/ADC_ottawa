!######################################################################
! chebyfd: a low-storage Chebyshev filter diagonalisation program
!          based on the use of Slepian filter functions
!######################################################################

module chebyfdmod

  use constants

  implicit none
  
  save

  ! Order domain autocorrelation function
  integer                :: order,kfinal,Kdim
  real(dp), allocatable  :: auto(:)

  ! Spectral bounds
  real(dp), dimension(2)  :: bounds
  
  ! DPSS information
  integer  :: ndpss
  integer  :: npts
  real(dp) :: fw
  
  ! Slepian functions
  real(dp), dimension(:,:), allocatable :: v

  ! Eigenvalues of the Slepians
  real(dp), dimension(:), allocatable :: lambda
  
  ! Energy interval
  real(dp) :: Ea,Eb

  ! Scaled energy interval
  real(dp) :: Eabar,Ebbar

  ! Expansion coefficients
  real(dp), dimension(:,:), allocatable :: fkn

  ! Filtered state overlap and Hamiltonian matrices
  real(dp), dimension(:,:), allocatable :: smat,hmat

  ! Working dimension of the filtered state basis
  integer :: nfsbas

  ! Eigenvectors and eigenvalues of the filtered-state
  ! Hamiltonian matrix
  real(dp), dimension(:,:), allocatable :: eigvec
  real(dp), dimension(:), allocatable   :: eigval
  real(dp), dimension(:), allocatable   :: ener
  
  ! Projector  onto the orthogonal complement of the
  ! null space
  real(dp), dimension(:,:), allocatable :: transmat

  ! Dimension of the orthogonal complement of the null space
  integer :: nrbas

  ! Transition dipoles and oscillator strengths
  real(dp), dimension(:), allocatable :: tdm,osc

  ! Output
  integer            :: idat
  character(len=120) :: adat

  ! Unit conversion factors
  real(dp), parameter :: eh2ev=27.2113845d0
  
end module chebyfdmod

!######################################################################

program chebyfd

  use chebyfdmod
  
  implicit none

  integer :: i
  
!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call openchebyfdfiles
  
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdchebyfdinp

!----------------------------------------------------------------------
! Read the Chebyshev order domain autocorrelation function file
!----------------------------------------------------------------------
  call rdautofile

!----------------------------------------------------------------------
! If the lower bound of the energy window is less than the lower
! spectral bound, then reset it.
!
! Note that we cannot set Ea to bounds(1) as this will end up
! introducing a singularity into the sqrt(1-Ebar^2) term appearing
! in the equation for the coefficients...
!----------------------------------------------------------------------
  if (Ea.lt.bounds(1)) Ea=bounds(1)*1.001d0

!----------------------------------------------------------------------
! Calculate the scaled energy window bounds
!----------------------------------------------------------------------
  Eabar=scalefunc(Ea)
  Ebbar=scalefunc(Eb)

!----------------------------------------------------------------------
! Calculate the Slepian filter functions
!----------------------------------------------------------------------
  call calc_slepians

!----------------------------------------------------------------------
! Output some information about the Slepians to the log file
!----------------------------------------------------------------------
  call wrslepinfo

!----------------------------------------------------------------------
! Calculate the coefficients entering into the expansion of the
! Slepian filter functions with respect to the Chebyshev polynomials
!----------------------------------------------------------------------
  call calc_expansion_coeffs

!----------------------------------------------------------------------
! Calculate the filtered-state overlap matrix
!----------------------------------------------------------------------
  call calc_smat_fsbas

!----------------------------------------------------------------------
! Calculate the filtered-state Hamiltonian matrix
!----------------------------------------------------------------------
  call calc_hmat_fsbas

!----------------------------------------------------------------------
! Calculate the eigenvalues
!----------------------------------------------------------------------
  call hmat_eigen
  
!----------------------------------------------------------------------
! Calculate the transition dipoles and oscillator strengths
!----------------------------------------------------------------------
  call calc_intens

!----------------------------------------------------------------------
! Output the spectrum
!----------------------------------------------------------------------
  call wrspec
  
contains
  
!######################################################################

  subroutine openchebyfdfiles

    use constants
    use channels
    use iomod
    use chebyfdmod
    
    implicit none

    integer :: k

!----------------------------------------------------------------------
! Get the name of the input file
!----------------------------------------------------------------------
    ain=''

    call getarg(1,ain)

    if (ain.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the input file has not &
            been given!'
       stop
    endif

    if (index(ain,'.inp').eq.0) ain=trim(ain)//'.inp'

!----------------------------------------------------------------------
! Log file name
!----------------------------------------------------------------------
    k=index(ain,'.inp')-1
    alog=ain(1:k)//'.log'

!----------------------------------------------------------------------
! Data file name
!----------------------------------------------------------------------
    k=index(ain,'.inp')-1
    adat=ain(1:k)//'.dat'

!----------------------------------------------------------------------
! Open the input file
!----------------------------------------------------------------------
    call freeunit(iin)
    open(iin,file=ain,form='formatted',status='old')

!----------------------------------------------------------------------
! Open the log file
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file=alog,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Open the data file
!----------------------------------------------------------------------
    call freeunit(idat)
    open(idat,file=adat,form='unformatted',status='unknown')
    
    return
    
  end subroutine openchebyfdfiles
    
!######################################################################

  subroutine rdchebyfdinp

    use constants
    use iomod
    use parsemod
    use channels
    use chebyfdmod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Number of Slepian filter functions
    ndpss=0

    ! Time half-bandwidth parameter
    fw=0.0d0

    ! Energy window
    Ea=-999.0d0
    Eb=-999.0d0

    ! Maximum order
    kfinal=1e+6
    
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

5   call rdinp(iin)
    
    i=0
    if (.not.lend) then
10     continue
       i=i+1
       
       if (keyword(i).eq.'window') then
          if (keyword(i+1).eq.'=') then
             ! Lower bound
             i=i+2
             read(keyword(i),*) Ea
             ! Upper bound
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) Eb
             else
                errmsg='The upper energy bound has not been given &
                     with the window keyword'
                call error_control
             endif
             ! Conversion to a.u.
             Ea=Ea/eh2ev
             Eb=Eb/eh2ev
          else
             goto 100
          endif

       else if (keyword(i).eq.'slepians') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             ! No. slepian filter functions
             read(keyword(i),*) ndpss
             ! Time half-bandwidth product factor
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) fw
             else
                errmsg='The ime half-bandwidth product factor has &
                     not been given with the slepians keyword'
                call error_control
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'kf') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) kfinal
             if (mod(kfinal,2).eq.0) kfinal=kfinal-1
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
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       call error_control
       
    endif

!----------------------------------------------------------------------
! Make sure that all required information has been given
!----------------------------------------------------------------------
    ! Energy window
    if (Ea.eq.-999.0d0) then
       errmsg='The energy window has not been given'
       call error_control
    endif

    ! Slepian filter functions
    if (ndpss.eq.0) then
       errmsg='The slepian filter functions have not been specified'
       call error_control
    endif
    
    return
    
  end subroutine rdchebyfdinp
    
!######################################################################

  subroutine rdautofile
  
    use constants
    use iomod
    use parsemod
    use chebyfdmod
    
    implicit none

    integer :: unit,k,itmp
    logical :: exists

!----------------------------------------------------------------------
! Make sure that the chebyauto file exists
!----------------------------------------------------------------------
    inquire(file='chebyauto',exist=exists)

    if (.not.exists) then
       errmsg='The chebyauto file could not be found'
       call error_control
    endif

!----------------------------------------------------------------------
! Open the chebyauto file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='chebyauto',form='formatted',status='old')

!----------------------------------------------------------------------
! Read the spectral bounds
!----------------------------------------------------------------------
    read(unit,'(21x,2(2x,E21.14))') bounds(1),bounds(2)

!----------------------------------------------------------------------
! Determine the order of the Chebyshev expansion and allocate the auto
! array
!----------------------------------------------------------------------
    ! Determine the order of the autocorrelation function
    read(unit,*)
    read(unit,*)
    order=-1
5   read(unit,*,end=10)
    order=order+1
    goto 5
10  continue

    ! Adjust order if it is greater than the user specified value
    if (order.gt.kfinal) order=kfinal

    ! Order to be used in the expansion of the Slepians
    Kdim=(order-1)/2

    ! Allocate the auto array
    allocate(auto(0:order))
    auto=0.0d0

!----------------------------------------------------------------------
! Read the Chebyshev order domain autocorrelation function
!----------------------------------------------------------------------
    rewind(unit)

    do k=1,3
       read(unit,*)
    enddo

    do k=0,order
       read(unit,*) itmp,auto(k)
    enddo
    
!----------------------------------------------------------------------
! Close the chebyauto file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rdautofile

!######################################################################

  function scalefunc(e) result(escale)

    use chebyfdmod

    implicit none
    
    real(dp) :: e,escale

    escale=e-(0.5d0*(bounds(2)-bounds(1))+bounds(1))
    escale=escale/(bounds(2)-bounds(1))
    escale=2.0d0*escale
    
    return
    
  end function scalefunc

!######################################################################

  function unscalefunc(escale) result(e)

    use chebyfdmod

    implicit none
    
    real(dp) :: e,escale,DeltaE

    DeltaE=bounds(2)-bounds(1)

    e=0.5d0*escale*DeltaE+0.5d0*DeltaE+bounds(1)    
    
    return
    
  end function unscalefunc
    
!######################################################################

  subroutine calc_slepians

    use chebyfdmod
    use dpssmt

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the DPSSs...'
    
!----------------------------------------------------------------------
! Set the number of points at which the DPSSs will be evaluated.
!----------------------------------------------------------------------
    npts=max(5001,Kdim)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Slepians
    allocate(v(npts,ndpss))
    v=0.0d0

    ! Eigenvalues
    allocate(lambda(ndpss))
    lambda=0.0d0

!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
    call dpss(npts,fw,ndpss,v,lambda)

!----------------------------------------------------------------------
! Renormalisation on the interval [Eabar,Ebbar]
!----------------------------------------------------------------------
    v=v/((Ebbar-Eabar)/npts)
    
    return
    
  end subroutine calc_slepians
    
!######################################################################

  subroutine wrslepinfo

    use channels
    use chebyfdmod
    
    implicit none

    integer  :: i,j
    real(dp) :: ovrlp
    
!----------------------------------------------------------------------
! Output some information about the Slepians to file
!----------------------------------------------------------------------
    ! fw
    write(ilog,'(a,2x,F10.7)') '# Time half-bandwidth product:',fw
    
    ! Eigenvalues
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(a)') '  DPSS       lambda         1 - lambda'
    write(ilog,'(41a)') ('#',i=1,41)
    do i=1,ndpss
       write(ilog,'(2x,i3,2(2x,ES15.6))') i,lambda(i),1.0d0-lambda(i)
    enddo

    ! Overlaps
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(a)') '  Overlaps'
    write(ilog,'(41a)') ('#',i=1,41)
    do i=1,ndpss-1
       do j=i,ndpss
          ovrlp=dpss_overlap(i,j)
          write(ilog,'(2(2x,i3),2x,ES15.6)') i,j,ovrlp
       enddo
    enddo
    
    return
    
  end subroutine wrslepinfo

!######################################################################

  function dpss_overlap(i,j) result(ovrlp)

    use chebyfdmod
    
    implicit none

    integer  :: i,j,n
    real(dp) :: ovrlp

    ovrlp=v(1,i)*v(1,j)/2.0d0
    do n=2,npts-1
       ovrlp=ovrlp+v(n,i)*v(n,j)
    enddo
    ovrlp=ovrlp+v(npts,i)*v(npts,j)/2.0d0
    
    return
    
  end function dpss_overlap

!######################################################################

  subroutine calc_expansion_coeffs

    use chebyfdmod
    use iomod
    
    implicit none

    integer                               :: n,j,k,unit
    real(dp)                              :: debar,ebar,theta
    real(dp), dimension(:,:), allocatable :: Tk,Tkw,val
    character(len=3)                      :: an
    character(len=60)                     :: filename
    logical                               :: exists
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the expansion coefficients...'
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(fkn(0:Kdim,ndpss))
    fkn=0.0d0

    allocate(Tk(0:Kdim,npts))
    Tkw=0.0d0
    
    allocate(Tkw(0:Kdim,npts))
    Tkw=0.0d0

    allocate(val(npts,ndpss))
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
    do n=1,ndpss
       write(an,'(i3)') n
       filename='dpss.cheby/dpss.cheby.'//trim(adjustl(an))//'.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       do j=1,npts
          write(unit,'(i5,2x,ES15.6)') j,val(j,n)*((Ebbar-Eabar)/npts)
       enddo
       close(unit)
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Tk)
    deallocate(Tkw)

    return
    
  end subroutine calc_expansion_coeffs

!######################################################################

  subroutine calc_smat_fsbas

    use channels
    use chebyfdmod
    
    implicit none

    integer                                :: m,n,j,k
    real(dp), dimension(:), allocatable    :: gramdet
    real(dp), parameter                    :: thrsh=1e-6_dp
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the overlap matrix...'
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(smat(ndpss,ndpss))
    smat=0.0d0

    allocate(gramdet(ndpss))
    gramdet=0.0d0
    
!----------------------------------------------------------------------
! Calculate the filtered state overlap matrix
!----------------------------------------------------------------------
    do m=1,ndpss
       do n=m,ndpss
          do j=0,Kdim
             do k=0,Kdim
                smat(m,n)=smat(m,n)+fkn(j,m)*fkn(k,n) &
                     *(auto(j+k)+auto(abs(k-j)))
             enddo
          enddo
          smat(n,m)=smat(m,n)
       enddo
    enddo
    
    ! Prefactor
    smat=0.5d0*smat

!----------------------------------------------------------------------
! Determine the no. filtered states to use in the calculation of the
! Hamiltonian matrix
!----------------------------------------------------------------------
    ! Gram determinants
    do n=1,ndpss
       gramdet(n)=finddet(smat(1:n,1:n),n)
    enddo

    ! Output the Gram determinant of the full-dimension overlap matrix
    write(ilog,'(/,2x,a,ES15.6)') 'det(S) = ',gramdet(ndpss)
    write(6,'(/,2x,a,ES15.6)') 'det(S) = ',gramdet(ndpss)

    ! Reduce the dimension of the filtered-state basis if necessary
    if (abs(gramdet(ndpss)).lt.thrsh) then
       do n=1,ndpss-1
          if (abs(gramdet(n+1)).lt.thrsh) then
             nfsbas=n
             exit
          endif
       enddo
    else
       nfsbas=ndpss
    endif

    if (nfsbas.lt.ndpss) then
       write(ilog,'(/,2x,a,i3)') &
            'Filtered-state basis dimension reduced to ',nfsbas
       write(6,'(/,2x,a,i3)') &
            'Filtered-state basis dimension reduced to ',nfsbas
    endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(gramdet)
    
    return
    
  end subroutine calc_smat_fsbas
    
!######################################################################

  real(dp) function finddet(matrix1,n)

    implicit none

    integer, intent(in)      :: n
    integer                  :: i,j,k,l
    real(dp), dimension(n,n) :: matrix,matrix1
    real(dp)                 :: m,temp
    logical                  :: detexists = .true.

    matrix=matrix1

!----------------------------------------------------------------------
! Convert to upper triangular form
!----------------------------------------------------------------------    
    l=1
    do k=1,n-1
       if (matrix(k,k)==0.0d0) then
          detexists = .false.
          do i=k+1,n
             if (matrix(i,k)/=0) then
                do j=1,n
                   temp=matrix(i,j)
                   matrix(i,j)=matrix(k,j)
                   matrix(k,j)=temp
                enddo
                detexists=.true.
                l=-l
                exit
             endif
          enddo
          if (detexists.eqv..false.) then
             finddet=0
             return
          endif
       endif
       do j=k+1,n
          m=matrix(j,k)/matrix(k,k)
          do i=k+1,n
             matrix(j,i)=matrix(j,i)-m*matrix(k,i)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate determinant by finding product of diagonal elements
!----------------------------------------------------------------------
    finddet=l
    do i=1,n
       finddet=finddet*matrix(i,i)
    enddo

    return
    
  end function finddet

!######################################################################

  subroutine calc_hmat_fsbas

    use chebyfdmod
    
    implicit none

    integer :: m,n,j,k

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the Hamiltonian matrix...'
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hmat(nfsbas,nfsbas))
    hmat=0.0d0

!----------------------------------------------------------------------
! Calculate the filtered state Hamiltonian matrix
!----------------------------------------------------------------------
    do m=1,nfsbas
       do n=m,nfsbas
          do j=0,Kdim
             do k=0,Kdim
                hmat(m,n)=hmat(m,n)+fkn(j,m)*fkn(k,n) &
                     *(auto(j+k+1)+auto(abs(j-k-1))+auto(abs(j+k-1)) &
                     +auto(abs(j-k+1)))
             enddo
          enddo
          hmat(n,m)=hmat(m,n)
       enddo
    enddo

    ! Prefactor
    hmat=0.25d0*hmat
    
    return
    
  end subroutine calc_hmat_fsbas

!######################################################################

  subroutine hmat_eigen

    use chebyfdmod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Solve the generalised eigenvalue problem for the filtered-state
! Hamiltonian
!----------------------------------------------------------------------
    call solve_geneig(hmat,smat(1:nfsbas,1:nfsbas),eigvec,eigval,&
         transmat,nfsbas,nrbas)
    
!----------------------------------------------------------------------
! Un-scale the eigenvalues
!----------------------------------------------------------------------
    allocate(ener(nrbas))

    do i=1,nrbas
       ener(i)=unscalefunc(eigval(i))
    enddo
    
    return
    
  end subroutine hmat_eigen
    
!######################################################################
! solve_geneig: solves the generalised eigenvalue problem
!               A V = B V E
!               For non-positive-definite matrices B, null space
!               vectors are discarded.
!######################################################################
  
  subroutine solve_geneig(A,B,Vbar,Ebar,P,matdim,rdim)

    use constants
    use iomod
        
    implicit none

    integer                               :: matdim,workdim,rdim,&
                                             error,i
    real(dp), dimension(matdim,matdim)    :: A,B

    real(dp), parameter                   :: thrsh=1e+2_dp
    
    real(dp), dimension(matdim,matdim)    :: U
    real(dp), dimension(matdim)           :: lambda
    real(dp), dimension(:,:), allocatable :: Ubar,normfac,P,Abar,Vbar
    real(dp), dimension(:), allocatable   :: Ebar
    real(dp), dimension(:), allocatable   :: work
    real(dp), parameter                   :: ovrthrsh=1e-4_dp
    
!----------------------------------------------------------------------
! Diagonalise B
!----------------------------------------------------------------------
    workdim=3*matdim
    allocate(work(workdim))

    U=B

    call dsyev('V','U',matdim,U,matdim,lambda,work,workdim,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the B-matrix failed in subroutine &
            solve_geneig'
       call error_control
    endif
    
    deallocate(work)

!----------------------------------------------------------------------
! Discard the null space eigenvectors of B and form the matrix that
! projects onto the orthogonal complement of the null space
!----------------------------------------------------------------------
    ! Number of eigenvectors
    rdim=0
    do i=1,matdim
       if (lambda(i).gt.ovrthrsh) rdim=rdim+1
    enddo

    ! Truncated eigenvector matrix
    allocate(Ubar(matdim,rdim))
    Ubar(:,1:rdim)=U(:,matdim-rdim+1:matdim)

    ! Normalisation factors
    allocate(normfac(rdim,rdim))
    normfac=0.0d0
    do i=1,rdim
       normfac(i,i)=sqrt(1.0d0/lambda(matdim-rdim+i))
    enddo

    ! Transformation matrix, P
    allocate(P(matdim,rdim))
    P=matmul(Ubar,normfac)

!----------------------------------------------------------------------
! Projection of the A-matrix onto the orthogonal complement of the
! null space
!----------------------------------------------------------------------
    allocate(Abar(rdim,rdim))
    Abar=matmul(transpose(P),matmul(A,P))

!----------------------------------------------------------------------
! Diagonalisation of Abar
!----------------------------------------------------------------------
    allocate(Vbar(rdim,rdim))
    allocate(Ebar(rdim))

    ! Diagonalise the reduced space Hamiltonian
    workdim=3*rdim
    allocate(work(workdim))
    Vbar=Abar
    call dsyev('V','U',rdim,Vbar,rdim,Ebar,work,workdim,error)
    deallocate(work)
    
    if (error.ne.0) then
       errmsg='Diagonalisation of the A-bar matrix failed in &
            subroutine solve_geneig'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Ubar)
    deallocate(normfac)
    
    return
    
  end subroutine solve_geneig
    
!######################################################################

  subroutine calc_intens

    use chebyfdmod
    
    implicit none

    integer :: i
    
!---------------------------------------------------------------------- 
! Allocate arrays
!---------------------------------------------------------------------- 
    allocate(tdm(nrbas))
    tdm=0.0d0

    allocate(osc(nrbas))
    osc=0.0d0

!---------------------------------------------------------------------- 
! Calculate the transition dipole moments
!---------------------------------------------------------------------- 
    tdm=matmul(transpose(eigvec),(matmul(transpose(transmat),&
         matmul(transpose(fkn(:,1:nfsbas)),auto(0:Kdim)))))

!---------------------------------------------------------------------- 
! Calculate the oscillator strengths
!---------------------------------------------------------------------- 
    do i=1,nrbas
       osc(i)=ener(i)*tdm(i)**2
    enddo
    osc=osc*2.0d0/3.0d0
    
    return
    
  end subroutine calc_intens

!######################################################################

  subroutine wrspec
    
    use iomod
    use chebyfdmod
    
    implicit none

    integer :: i,unit
    
!---------------------------------------------------------------------- 
! Open the output file
!---------------------------------------------------------------------- 
    call freeunit(unit)
    open(unit,file='chebyfd_eig.dat',form='formatted',status='unknown')

!---------------------------------------------------------------------- 
! Write the spectrum to file
!---------------------------------------------------------------------- 
    ! Table header
    write(unit,'(29a)') ('#',i=1,29)
    write(unit,'(a)') '#  Energy          Intensity'
    write(unit,'(29a)') ('#',i=1,29)

    ! Transition energies and oscillator strengths
    do i=1,nrbas
       write(unit,'(2(2x,F12.7))') ener(i)*eh2ev,osc(i)
    enddo
    
!---------------------------------------------------------------------- 
! Close the output file
!---------------------------------------------------------------------- 
    close(unit)
    
    return
    
  end subroutine wrspec
    
!######################################################################
  
end program chebyfd
