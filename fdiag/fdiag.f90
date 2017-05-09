!######################################################################
! fdiag: a program to extract bound-state spectra from wavepacket
!        autocorrelation functions using the filter diagonalisation
!        method.
!        Uses the Beck-Meyer approach as detailed in:
!        J. Chem. Phys., 109, 3730 (1998)
!######################################################################
program fdiag

  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call openfdiagfiles
  
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call rdfdiaginp

!----------------------------------------------------------------------
! Read the autocorrelation function files
!----------------------------------------------------------------------
  call rdauto

!----------------------------------------------------------------------
! Calculate the Hamiltonian and overlap matrix elements in the basis
! of the filter states
!----------------------------------------------------------------------
  call calc_matrices_fstatebas

!----------------------------------------------------------------------
! Perform the canonical orthogonalisation of the filter state basis
! to yield a reduced basis of linearly independent functions
!----------------------------------------------------------------------
  call reduce_fstatebas

!----------------------------------------------------------------------
! Finalisation and deallocation of arrays
!----------------------------------------------------------------------
  call fdiag_finalise
  
contains

!######################################################################

  subroutine openfdiagfiles

    use constants
    use channels
    use iomod
    
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
! Open the input file
!----------------------------------------------------------------------
    call freeunit(iin)
    open(iin,file=ain,form='formatted',status='old')

!----------------------------------------------------------------------
! Open the log file
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file=alog,form='formatted',status='unknown')
    
    return
    
  end subroutine openfdiagfiles
    
!######################################################################

  subroutine rdfdiaginp

    use constants
    use parsemod
    use channels
    use iomod
    use filtermod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Energy window
    ebound=-999.0d0
    nener=0

    ! Window function
    iwfunc=-1
    
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

5   call rdinp(iin)
    
    i=0
    if (.not.lend.and.keyword(i).ne.'end-input') then
10     continue
       i=i+1      
       
       if (keyword(i).eq.'window') then
          if (keyword(i+1).eq.'=') then
             ! Lower bound
             i=i+2
             read(keyword(i),*) ebound(1)
             ! Upper bound
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) ebound(2)
             else
                errmsg='The upper energy bound and no. points have &
                     not been given with the window keyword'
                call error_control
             endif
             ! No. points
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) nener
             else
                errmsg='The no. points has not been given with the &
                     window keyword'
                call error_control
             endif
             ! Conversion to a.u.
             ebound=ebound*ev2eh
             ! Energy grid spacing
             de=(ebound(2)-ebound(1))/(nener-1)
          else
             goto 100
          endif
             
       else if (keyword(i).eq.'window_function') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'cos0') then
                iwfunc=0
             else if (keyword(i).eq.'cos1') then
                iwfunc=1
             else if (keyword(i).eq.'cos2') then
                iwfunc=2
             else
                errmsg='Unknown window function: '//trim(keyword(i))
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
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       call error_control
       
    endif

!----------------------------------------------------------------------
! Check that all the required information has been given
!----------------------------------------------------------------------
    if (ebound(1).eq.-999.0d0) then
       errmsg='The energy window has not been given'
       call error_control
    endif

    if (iwfunc.eq.-1) then
       errmsg='The window function has not been given'
       call error_control
    endif

!----------------------------------------------------------------------
! Write some information to the log file
!----------------------------------------------------------------------
    ! Energy window
    write(ilog,'(/,2x,a,2(2x,F9.4))') &
         'Energy window (eV): ',ebound(1)/ev2eh,ebound(2)/ev2eh
    write(ilog,'(/,2x,a,i4)') 'No. points: ',nener
    
    ! Window function
    if (iwfunc.eq.0) then
       write(ilog,'(/,2x,a)') 'Window function: cos0'
    else if (iwfunc.eq.1) then
       write(ilog,'(/,2x,a)') 'Window function: cos1'
    else if (iwfunc.eq.2) then
       write(ilog,'(/,2x,a)') 'Window function: cos2'
    endif

    return
    
  end subroutine rdfdiaginp
  
!######################################################################

  subroutine rdauto

    use constants
    use iomod
    use parsemod
    use channels
    use filtermod
    
    implicit none

    integer :: iauto,n
    real(d) :: t1,t2,re,im

!----------------------------------------------------------------------
! Determine the no. timesteps and dt from the a0(t) file
!----------------------------------------------------------------------
    call freeunit(iauto)
    open(iauto,file='auto',form='formatted',status='old',err=100)

    ! Timestep, dt
    call rdinp(iauto)
    read(keyword(1),*) t1
    call rdinp(iauto)
    read(keyword(1),*) t2
    dt=(t2-t1)*fs2au

    ! No. timesteps, nt
    nt=0
    rewind(iauto)
5   call rdinp(iauto)
    if (.not.lend) then
       nt=nt+1
       goto 5
    endif

    close(iauto)
    
!----------------------------------------------------------------------
! Allocate the autocorrelation function arrays
!----------------------------------------------------------------------
    allocate(auto(nt))
    allocate(auto1(nt))
    allocate(auto2(nt))

!----------------------------------------------------------------------
! Read the a0(t) file
!----------------------------------------------------------------------
    call freeunit(iauto)
    open(iauto,file='auto',form='formatted',status='old',err=100)

    n=0
10  call rdinp(iauto)
    if (.not.lend) then
       n=n+1
       read(keyword(2),*) re
       read(keyword(3),*) im
       auto(n)=dcmplx(re,im)
       goto 10
    endif

    close(iauto)

!----------------------------------------------------------------------
! Read the a1(t) file
!----------------------------------------------------------------------
    call freeunit(iauto)
    open(iauto,file='auto1',form='formatted',status='old',err=101)

    n=0
15  call rdinp(iauto)
    if (.not.lend) then
       n=n+1
       read(keyword(2),*) re
       read(keyword(3),*) im
       auto1(n)=dcmplx(re,im)
       goto 15
    endif
    
    close(iauto)

!----------------------------------------------------------------------
! Read the a2(t) file
!----------------------------------------------------------------------
    call freeunit(iauto)
    open(iauto,file='auto2',form='formatted',status='old',err=102)

    n=0
20  call rdinp(iauto)
    if (.not.lend) then
       n=n+1
       read(keyword(2),*) re
       read(keyword(3),*) im
       auto2(n)=dcmplx(re,im)
       goto 20
    endif
    
    close(iauto)

!----------------------------------------------------------------------
! Set: (1) proptime: the wavepacket propagation time, and;
!      (2) autotime: the time for which the autocorrelation functions
!                    are known. We assume that autotime = 2*proptime,
!                    i.e., that the t/2 trick was used in calculating
!                    the autocorrelation functions.
!----------------------------------------------------------------------
    autotime=(nt-1)*dt
    proptime=0.5d0*autotime

!----------------------------------------------------------------------
! Ouput some information to the log file and return
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,x,i4)') 'No. timesteps:',nt
    write(ilog,'(/,2x,a,x,F8.6)') 'Timestep, dt (fs): ',dt*au2fs
    
    return

!----------------------------------------------------------------------
! Error control for missing autocorrelation function files
!----------------------------------------------------------------------
100 continue
    errmsg='The auto file could not be found'
    call error_control

101 continue
    errmsg='The auto1 file could not be found'
    call error_control

102 continue
    errmsg='The auto2 file could not be found'
    call error_control
    
  end subroutine rdauto

!######################################################################

  subroutine calc_matrices_fstatebas

    use constants
    use filtermod
    
    implicit none

    integer                          :: L,i,j,n,ti,Iindx,Jindx
    real(d), dimension(0:nener-1,nt) :: Gk
    real(d), dimension(2:2*nener)    :: ebar
    real(d)                          :: t,e,h,func

    L=nener
    
!----------------------------------------------------------------------
! Calculation of the value of the filter function, G, at the energy
! grid points for all timesteps
!----------------------------------------------------------------------
    do ti=1,nt
       do n=0,L-1
          Gk(n,ti)=filterfunc(n,ti)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculation of all 2L-1 possible values of Ebar
!----------------------------------------------------------------------
    do n=2,2*L
       ebar(n)=ebound(1)+0.5d0*de*n-de
    enddo

!----------------------------------------------------------------------
! Calculation of filter state Hamiltonian and overlap matrices
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(hfbas(L,L))
    hfbas=0.0d0
    allocate(h2fbas(L,L))
    h2fbas=0.0d0
    allocate(sfbas(L,L))
    sfbas=0.0d0

    ! Loop over the elements of the lower triangle of Hfbas and Sfbas
    do j=1,L
       do i=j,L
          
          ! Delta E index
          Iindx=abs(i-j)
          
          ! Ebar index
          Jindx=i+j
          
          ! Calculate H_ij, H^2_ij and S_ij using the trapezoidal rule
          call calc_matrices_fstatebas_1element(hfbas(i,j),&
               h2fbas(i,j),sfbas(i,j),Gk(Iindx,:),ebar(Jindx))

          ! Upper triangle elements, H_ji and S_ji
          hfbas(j,i)=hfbas(i,j)
          h2fbas(j,i)=h2fbas(i,j)
          sfbas(j,i)=sfbas(i,j)

       enddo
    enddo
    
    return
    
  end subroutine calc_matrices_fstatebas

!######################################################################

  subroutine calc_matrices_fstatebas_1element(hij,h2ij,sij,Gk,ebar)

    use constants
    use filtermod

    implicit none

    integer                :: ti
    real(d)                :: hij,h2ij,sij,ebar,h,t
    real(d), dimension(nt) :: Gk

!----------------------------------------------------------------------
! Contribution from the first time point
!----------------------------------------------------------------------
    hij=0.5d0*dt*Gk(1)*real(auto1(1))

    h2ij=0.5d0*dt*Gk(1)*real(auto2(1))

    sij=0.5d0*dt*Gk(1)*real(auto(1))

!----------------------------------------------------------------------
! Contribution from time points 2 to Nt-1
!----------------------------------------------------------------------
    do ti=2,nt-1
             
       t=(ti-1)*dt
       
       ! H_ij
       hij=hij+dt*Gk(ti)*(real(auto1(ti))*cos(ebar*t) &
            -imag(auto1(ti))*sin(ebar*t))

       ! H^2_ij
       h2ij=h2ij+dt*Gk(ti)*(real(auto2(ti))*cos(ebar*t) &
            -imag(auto2(ti))*sin(ebar*t))
       
       ! S_ij
       sij=sij+dt*Gk(ti)*(real(auto(ti))*cos(ebar*t) &
            -imag(auto(ti))*sin(ebar*t))
       
    enddo
    
!----------------------------------------------------------------------
! Contribution from the final time point
!----------------------------------------------------------------------
    t=(nt-1)*dt

    hij=hij+0.5d0*dt*Gk(nt)*(real(auto1(nt))*cos(ebar*t)&
         -imag(auto1(nt))*sin(ebar*t))

    h2ij=h2ij+0.5d0*dt*Gk(nt)*(real(auto2(nt))*cos(ebar*t)&
         -imag(auto2(nt))*sin(ebar*t))
    
    sij=sij+0.5d0*dt*Gk(nt)*(real(auto(nt))*cos(ebar*t) &
         -imag(auto(nt))*sin(ebar*t))

    return
    
  end subroutine calc_matrices_fstatebas_1element

!######################################################################
! filterfunc: for the distance n=|i-j| between two energy grid points
!             and the time grid index tk, this function calculates
!             the corresponding vale of the filter function
!             G_k(Delta E_ij, t_k)
!######################################################################    

  function filterfunc(n,tk) result(func)

    use constants
    use filtermod

    implicit none

    integer :: n,tk
    real(d) :: func,deltae,t,kappa,cosfunc

    ! Delta E and t
    deltae=n*de
    t=(tk-1)*dt

    ! kappa(t)
    kappa=1.0d0-t/(2.0d0*proptime)

    ! Contribution from the stepfunction Theta(kappa(t)) that is
    ! common to all cosine-type filter functions 
    if (kappa.lt.0) then
       func=0.0d0
       return
    endif

    select case(iwfunc)

!----------------------------------------------------------------------
! Filter function G_0 for the cos0 window function
!----------------------------------------------------------------------
    case(0)
       
       func=4.0d0*proptime*kappa*dsinc(deltae*proptime*kappa)

       return

!----------------------------------------------------------------------
! Filter function G_1 for the cos1 window function
!----------------------------------------------------------------------
    case(1)
       
       func=proptime*kappa*(dsinc((deltae*proptime-pi)*kappa) &
            +2.0d0*cos((pi*t)/(2.0d0*proptime))*dsinc(deltae*proptime*kappa) &
            +dsinc((deltae*proptime+pi)*kappa))

       return

!----------------------------------------------------------------------
! Filter function G_2 for the cos2 window function
!----------------------------------------------------------------------
    case(2)

       cosfunc=cos((pi*t)/(2.0d0*proptime))

       func=proptime*kappa*(0.25d0*dsinc((deltae*proptime-2.0d0*pi)*kappa) &
            +cosfunc*dsinc((deltae*proptime-pi)*kappa) &
            +(0.5d0+cosfunc**2)*dsinc(deltae*proptime*kappa) &
            +cosfunc*dsinc((deltae*proptime+pi)*kappa) &
            +0.25d0*dsinc((deltae*proptime+2.0d0*pi)*kappa))

       return

    end select

    return

  end function filterfunc

!######################################################################
! dsinc: Computes the value of the function dsinc(x) := dsin(x)/x.
!        Special care is taken for the singularity at x = 0.
!        Adapted from the Quantics code.
!######################################################################

  function dsinc(x)

    use constants

    implicit none

    real(d)            :: dsinc,x,x2
    real(d), parameter :: inv6=1.0d0/6.0d0
    real(d), parameter :: inv20=1.0d0/20.0d0
    real(d), parameter :: inv42=1.0d0/42.0d0
    real(d), parameter :: tiny = 4.154252d-2

! "Tiny" is chosen such that tiny**8 / 8! = eps where "eps" is the
! machine precision (here eps = 2.2 10^(-16)).

    if (abs(x).ge.tiny) then
       dsinc=sin(x)/x
    else
       x2=x**2
       dsinc=1.0d0-inv6*x2*(1.0d0-inv20*x2*(1.0d0-inv42*x2))
    endif

    return

  end function dsinc

!######################################################################

  subroutine reduce_fstatebas

    use constants
    use iomod
    use filtermod

    implicit none

    integer                            :: L,error,workdim,i,j
    real(d), parameter                 :: thrsh=1e-8_d
    real(d), dimension(nener,nener)    :: eigvec
    real(d), dimension(nener)          :: eigval
    real(d), dimension(:), allocatable :: work
    
    L=nener

!----------------------------------------------------------------------
! Diagonalise the filter state overlap matrix, Sfbas
!----------------------------------------------------------------------
    workdim=3*L
    allocate(work(workdim))
    
    eigvec=sfbas

    call dsyev('V','U',L,eigvec,L,eigval,work,workdim,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the filter state overlap matrix &
            failed in subroutine reduce_fstatebas'
       call error_control
    endif

    deallocate(work)
    
!----------------------------------------------------------------------
! Construct the transformation matrix
!----------------------------------------------------------------------
    !************************************************
    ! WHY ISN'T Sfbas POSITIVE SEMIDEFINITE?
    !************************************************
    
    ! Number of eigenvectors orthogonal to the null space
    nrbas=0
    do i=1,L
       if (eigval(i).gt.thrsh) nrbas=nrbas+1
       !print*,i,eigval(i)
    enddo

    ! Truncated eigenvector matrix
    allocate(ubar(L,nrbas))
    ubar(:,1:nrbas)=eigvec(:,L-nrbas+1:L)
    
    ! Normalisation factors
    allocate(normfac(nrbas,nrbas))
    normfac=0.0d0
    do i=1,nrbas
       normfac(i,i)=sqrt(1.0d0/eigval(L-nrbas+i))
    enddo
    
    ! Transformation matrix
    allocate(transmat(L,nrbas))
    transmat=matmul(ubar,normfac)
    
!----------------------------------------------------------------------
! Hamiltonian matrices projected onto the space spanned by the
! reduced basis
!----------------------------------------------------------------------
    ! < psi_i | H | psi_j >
    allocate(hrbas(nrbas,nrbas))
    hrbas=matmul(transpose(transmat),matmul(hfbas,transmat))

    ! < psi_i | H^2 | psi_j >
    allocate(h2rbas(nrbas,nrbas))
    h2rbas=matmul(transpose(transmat),matmul(h2fbas,transmat))

!----------------------------------------------------------------------
! Diagonalise the reduced space Hamiltonian matrix
!----------------------------------------------------------------------
    ! Allocate the arrays holding the eigenvectors and eigenvalues
    ! of the reduced space hamiltonian
    allocate(rvec(nrbas,nrbas))
    allocate(rener(nrbas))

    ! Diagonalise the reduced space Hamiltonian
    workdim=3*nrbas
    allocate(work(workdim))
    rvec=hrbas
    call dsyev('V','U',nrbas,rvec,nrbas,rener,work,workdim,error)
    deallocate(work)
    
    if (error.ne.0) then
       errmsg='Diagonalisation of the reduced space Hamiltonian &
            failed in subroutine reduce_fstatebas'
       call error_control
    endif
    
    return

  end subroutine reduce_fstatebas

!######################################################################

  subroutine fdiag_finalise

    use channels
    use filtermod
    
    implicit none

!----------------------------------------------------------------------
! Deallocate arrays    
!----------------------------------------------------------------------
    deallocate(auto)
    deallocate(auto1)
    deallocate(auto2)
    deallocate(hfbas)
    deallocate(h2fbas)
    deallocate(sfbas)
    deallocate(hrbas)
    deallocate(ubar)
    deallocate(transmat)
    deallocate(rvec)
    deallocate(rener)
    
!----------------------------------------------------------------------
! Close files
!----------------------------------------------------------------------
    close(iin)
    close(ilog)
    
    return
    
  end subroutine fdiag_finalise
    
!######################################################################
  
end program fdiag
