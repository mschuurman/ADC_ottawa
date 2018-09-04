!######################################################################
! chebyfd: a low-storage Chebyshev filter diagonalisation program
!          based on the use of Slepian or Gaussian filter functions
!######################################################################

program chebyfd

  use cfdmod
 
  implicit none

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
! Scaled Gaussian width parameter
!----------------------------------------------------------------------
  if (ifilter.eq.2) then
     sigmabar=scalefunc(Ea+sigma)-Eabar
  endif
  
!----------------------------------------------------------------------
! Calculate the coefficients entering into the expansion of the
! filter functions with respect to the Chebyshev polynomials
!----------------------------------------------------------------------
  call get_coeffs

!----------------------------------------------------------------------
! Normalisation of the filtered state basis
!----------------------------------------------------------------------
  call normalise_fsbas
  
!----------------------------------------------------------------------
! Calculate the filtered-state overlap matrix
!----------------------------------------------------------------------
  call calc_smat_fsbas

!----------------------------------------------------------------------
! Analysis of the filtered-state overlap matrix
!----------------------------------------------------------------------
  call smat_ana
  
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
    use cfdmod
    
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
    use cfdmod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Number of filter functions
    nfsbas=0

    ! Filter function type
    ifilter=0
    
    ! Energy window
    Ea=-999.0d0
    Eb=-999.0d0

    ! Maximum order
    kfinal=1e+6

    ! Energies in a.u.
    lau=.false.

    ! Time half-bandwidth parameter for Slepian filter functions
    fw=0.0d0

    ! Width parameter for Gaussian filte functions
    sigma=0.0d0
    
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
          else
             goto 100
          endif

       else if (keyword(i).eq.'slepians') then
          ifilter=1
          if (keyword(i+1).eq.'=') then
             i=i+2
             ! No. slepian filter functions
             read(keyword(i),*) nfsbas
             ! Optional argument: time-bandwidth product factor
             if (keyword(i+1).eq.',') then
                i=i+2
                if (keyword(i).eq.'variable') then
                   ! Variable time-bandwidth products
                   varfw=.true.
                else
                   read(keyword(i),*) fw
                endif
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'gaussians') then
          ifilter=2
          if (keyword(i+1).eq.'=') then
             i=i+2
             ! No. Gaussian filter functions
             read(keyword(i),*) nfsbas
             ! Gaussian width parameter
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) sigma
             else
                errmsg='The Gaussian width parameter has not been &
                     given with the gaussians keyword'
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

       else if (keyword(i).eq.'au') then
          lau=.true.
          
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

    ! Number of filter functions
    if (nfsbas.eq.0) then
       errmsg='The filter functions have not been specified'
       call error_control
    endif

!----------------------------------------------------------------------
! Set the energy conversion factor
!----------------------------------------------------------------------
    if (lau) then
       convfac=1.0d0
    else
       convfac=eh2ev
    endif

!----------------------------------------------------------------------
! Conversion of the energy window bounds
!----------------------------------------------------------------------    
    ! Conversion to a.u.
    Ea=Ea/convfac
    Eb=Eb/convfac

!----------------------------------------------------------------------
! Conversion of the Gaussian filter function width to a.u.
!----------------------------------------------------------------------
    sigma=sigma/convfac
    
    return
    
  end subroutine rdchebyfdinp

!######################################################################

  subroutine rdautofile
  
    use constants
    use iomod
    use parsemod
    use cfdmod
    
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

    ! Order to be used in the expansion of the filter functions
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

    use cfdmod

    implicit none
    
    real(dp) :: e,escale

    escale=e-(0.5d0*(bounds(2)-bounds(1))+bounds(1))
    escale=escale/(bounds(2)-bounds(1))
    escale=2.0d0*escale
    
    return
    
  end function scalefunc

!######################################################################

  function unscalefunc(escale) result(e)

    use cfdmod

    implicit none
    
    real(dp) :: e,escale,DeltaE

    DeltaE=bounds(2)-bounds(1)

    e=0.5d0*escale*DeltaE+0.5d0*DeltaE+bounds(1)    
    
    return
    
  end function unscalefunc

!######################################################################

  subroutine get_coeffs

    use slepianmod
    use gaussianmod

    implicit none
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(fkn(0:Kdim,nfsbas))
    fkn=0.0d0

!----------------------------------------------------------------------
! Calculate the expansion coefficients
!----------------------------------------------------------------------    
    select case(ifilter)

    case(1) ! Slepian filter functions
       call get_coeffs_slepians

    case(2) ! Gaussian filter functions
       call get_coeffs_gaussians
       
    end select
       
    return
    
  end subroutine get_coeffs

!######################################################################

  subroutine normalise_fsbas

    use constants
    use cfdmod
    
    implicit none

    integer                             :: n
    real(dp), dimension(:), allocatable :: norm

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(norm(nfsbas))
    norm=0.0d0
    
!----------------------------------------------------------------------
! Calculate the norms of the filtered state basis functions
!----------------------------------------------------------------------
    do n=1,nfsbas
       norm(n)=smat_1element(n,n)
    enddo
    norm=sqrt(norm)

!----------------------------------------------------------------------
! Scale the expansion coefficients s.t. they correspond to the
! normalised filtered state basis functions
!----------------------------------------------------------------------
    do n=1,nfsbas
       fkn(:,n)=fkn(:,n)/norm(n)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(norm)
    
    return
    
  end subroutine normalise_fsbas
    
!######################################################################

  subroutine calc_smat_fsbas

    use channels
    use cfdmod
    
    implicit none

    integer :: m,n,j,k
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the overlap matrix...'
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(smat(nfsbas,nfsbas))
    smat=0.0d0
    
!----------------------------------------------------------------------
! Calculate the filtered state overlap matrix
!----------------------------------------------------------------------
    !$omp parallel do &
    !$omp& private(m,n) &
    !$omp& shared(smat)
    do m=1,nfsbas
       do n=m,nfsbas
          smat(m,n)=smat_1element(m,n)
          smat(n,m)=smat(m,n)
       enddo
    enddo
    !$omp end parallel do
    
    return
    
  end subroutine calc_smat_fsbas

!######################################################################

  function smat_1element(m,n) result(Smn)

    use channels
    use cfdmod
    
    implicit none

    integer  :: m,n,j,k
    real(dp) :: Smn

    Smn=0.0d0
    do j=0,Kdim
       do k=0,Kdim
          Smn=Smn+fkn(j,m)*fkn(k,n)*(auto(j+k)+auto(abs(k-j)))
       enddo
    enddo
    Smn=0.5d0*Smn
    
    return
    
  end function smat_1element
    
!######################################################################

  subroutine smat_ana

    use channels
    use cfdmod
    
    implicit none

    integer                                :: m,n,i
    real(dp), dimension(:), allocatable    :: gramdet
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(gramdet(nfsbas))
    gramdet=0.0d0

!----------------------------------------------------------------------
! Calculate the Gram determinant for S_n, n=1,...,nfsbas
!----------------------------------------------------------------------
    ! Gram determinants
    do n=1,nfsbas
       !gramdet(n)=ludet(snorm(1:n,1:n),n)

       gramdet(n)=ludet(smat(1:n,1:n),n)

    enddo

!----------------------------------------------------------------------
! Output the Gram determinants
!----------------------------------------------------------------------
    write(ilog,'(/,41a)') ('#',i=1,41)
    write(ilog,'(2x,a)') 'Gram determinants'
    write(ilog,'(41a)') ('#',i=1,41)
    do n=1,nfsbas
       write(ilog,'(2x,i3,2x,ES15.6)') n,gramdet(n)
    enddo

    write(6,'(/,2x,a,ES15.6)') 'det(S) = ',gramdet(nfsbas)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(gramdet)

    return
    
  end subroutine smat_ana

!######################################################################

  function ludet(matrix,n) result(det)

    use constants
    use iomod
    
    implicit none

    integer                  :: i,n,info
    integer, dimension(n)    :: ipiv
    real(dp), dimension(n,n) :: matrix,A
    real(dp)                 :: det

!----------------------------------------------------------------------
! LU decomposition of the input matrix
!----------------------------------------------------------------------
    A=matrix

    call dgetrf(n,n,A,n,ipiv,info)

    if (info.ne.0) then
       errmsg='LU decomposition failed in subroutine ludet'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculation of the determinant of the input matrix
!----------------------------------------------------------------------
    det=1.0d0

    do i=1,n
       det=det*A(i,i)
       if (ipiv(i).ne.1) det=-1*det
    enddo

    return

  end function ludet
  
!######################################################################

  subroutine calc_hmat_fsbas

    use cfdmod
    
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
    !$omp parallel do &
    !$omp& private(m,n) &
    !$omp& shared(hmat)
    do m=1,nfsbas
       do n=m,nfsbas
          hmat(m,n)=hmat_1element(m,n)
          hmat(n,m)=hmat(m,n)
       enddo
    enddo
    !$omp end parallel do
    
    return
    
  end subroutine calc_hmat_fsbas

!######################################################################

  function hmat_1element(m,n) result(Hmn)

    use channels
    use cfdmod
    
    implicit none

    integer  :: m,n,j,k
    real(dp) :: Hmn

    Hmn=0.0d0
    do j=0,Kdim
       do k=0,Kdim
          Hmn=Hmn+fkn(j,m)*fkn(k,n) &
               *(auto(j+k+1)+auto(abs(j-k-1))+auto(abs(j+k-1)) &
               +auto(abs(j-k+1)))
       enddo
    enddo
    Hmn=0.25d0*Hmn
    
    return
    
  end function hmat_1element
    
!######################################################################

  subroutine hmat_eigen

    use cfdmod
    
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
    use channels
    use iomod
        
    implicit none

    integer                               :: matdim,workdim,rdim,&
                                             error,i,nnull
    real(dp), dimension(matdim,matdim)    :: A,B
    real(dp), dimension(matdim,matdim)    :: U
    real(dp), dimension(matdim)           :: lambda
    real(dp), dimension(:,:), allocatable :: Ubar,normfac,P,Abar,Vbar
    real(dp), dimension(:), allocatable   :: Ebar
    real(dp), dimension(:), allocatable   :: work
    real(dp), parameter                   :: ovrthrsh=1e-10_dp

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
    nnull=0
    rdim=0
    do i=1,matdim
       if (lambda(i).gt.ovrthrsh) then
          rdim=rdim+1
       else
          nnull=nnull+1
       endif
    enddo

    ! Output the no. null space vectors
    write(ilog,'(/,2x,a,i3)') 'Number of null space vectors:',nnull
    write(6,'(/,2x,a,i3)') 'Number of null space vectors:',nnull
    
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

    use cfdmod
    
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
    use cfdmod
    
    implicit none

    integer :: i,unit

!    ! CHECK
!    integer  :: n
!    real(dp) :: en,D,b
!    ! CHECK
    
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
    write(unit,'(a)') '#  Energy        Intensity'
    write(unit,'(29a)') ('#',i=1,29)

    ! Transition energies and oscillator strengths
    do i=1,nrbas

       ! Skip if the transition energy is not in the interval [Ea,Eb]
       if (ener(i).lt.Ea.or.ener(i).gt.Eb) cycle

       write(unit,'(2(2x,F10.5))') ener(i)*convfac,osc(i)

    enddo


!    ! CHECK
!    D=1000.0d0
!    b=0.15d0
!
!    n=99
!    do i=1,nrbas
!
!       if (ener(i).lt.Ea.or.ener(i).gt.Eb) cycle
!       
!       n=n+1
!       
!       en=-D+2.0d0*D*(n+0.5d0)*sqrt(b**2/(2.0d0*D))-0.5d0*(n+0.5d0)**2*b**2
!
!       print*,n,ener(i),en,ener(i)-en
!       
!    enddo
!    ! CHECK
    
!---------------------------------------------------------------------- 
! Close the output file
!---------------------------------------------------------------------- 
    close(unit)
    
    return
    
  end subroutine wrspec
    
!######################################################################
  
end program chebyfd
