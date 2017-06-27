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
! Calculate the energies
!----------------------------------------------------------------------
  call calc_ener
  
!----------------------------------------------------------------------
! Calculate the intensities
!----------------------------------------------------------------------
  call calc_intens

!----------------------------------------------------------------------
! Calculate the error estimates for the energies and intensities
!----------------------------------------------------------------------
  call calc_errors

!----------------------------------------------------------------------
! Write the calculated energies and intensities to file
!----------------------------------------------------------------------
  call wrout
  
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
    ! Timestep cutofff
    tcutoff=0.0d0

    ! Energy window
    ebound=-999.0d0
    nener=0

    ! Window function
    iwfunc=-1

    ! Threshold for discarding eigenpairs of the overlap matrix
    ovrthrsh=1e-6_d

    ! Error estimate threshold for printing energies/intensities
    errthrsh=1e-2_d
    
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

       else if (keyword(i).eq.'overthresh') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) ovrthrsh
          else
             goto 100
          endif

       else if (keyword(i).eq.'errorthresh') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) errthrsh
          else
             goto 100
          endif
          
       else if (keyword(i).eq.'tcutoff') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) tcutoff
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
! Checks on the user supplied input
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

    ! No. timesteps that the autocorrelation function is available
    ! for (ntauto)
    ntauto=0
    rewind(iauto)
5   call rdinp(iauto)
    if (.not.lend) then
       ntauto=ntauto+1
       goto 5
    endif
    
    ! Truncation the autocorrelation functions
    if (tcutoff.gt.0.0d0.and.tcutoff.lt.(ntauto-1)*dt) then
       ntauto=int(tcutoff/dt)+1
    endif

    ! No. timesteps in the wavepacket propagation (assuming that
    ! the t/2 trick was used)
    ntprop=(ntauto-1)/2+1
        
    close(iauto)
    
!----------------------------------------------------------------------
! Allocate the autocorrelation function arrays
!----------------------------------------------------------------------
    allocate(auto(ntauto))
    allocate(auto1(ntauto))
    allocate(auto2(ntauto))

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
       if (n.lt.ntauto) goto 10
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
       if (n.lt.ntauto) goto 15
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
       if (n.lt.ntauto) goto 20
    endif
    
    close(iauto)

!----------------------------------------------------------------------
! Set: (1) proptime: the wavepacket propagation time, and;
!      (2) autotime: the time for which the autocorrelation functions
!                    are known. We assume that autotime = 2*proptime,
!                    i.e., that the t/2 trick was used in calculating
!                    the autocorrelation functions.
!----------------------------------------------------------------------
    autotime=(ntauto-1)*dt
    proptime=0.5d0*autotime

!----------------------------------------------------------------------
! Ouput some information to the log file and return
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,x,i4)') 'No. timesteps:',ntauto
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

    integer                              :: L,i,j,n,ti,Iindx,Jindx
    real(d), dimension(0:nener-1,ntauto) :: Gk
    real(d), dimension(2:2*nener)        :: ebar
    real(d)                              :: t,e,h,func

    L=nener
    
!----------------------------------------------------------------------
! Calculation of the value of the filter function, G, at the energy
! grid points for all timesteps
!----------------------------------------------------------------------
    do ti=1,ntauto
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

    integer                    :: ti
    real(d)                    :: hij,h2ij,sij,ebar,h,t
    real(d), dimension(ntauto) :: Gk

!----------------------------------------------------------------------
! Contribution from the first timestep
!----------------------------------------------------------------------
    hij=0.5d0*dt*Gk(1)*real(auto1(1))

    h2ij=0.5d0*dt*Gk(1)*real(auto2(1))

    sij=0.5d0*dt*Gk(1)*real(auto(1))

!----------------------------------------------------------------------
! Contribution from timesteps 2 to ntauto-1
!----------------------------------------------------------------------
    do ti=2,ntauto-1
             
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
! Contribution from the final timestep
!----------------------------------------------------------------------
    t=(ntauto-1)*dt

    hij=hij+0.5d0*dt*Gk(ntauto)*(real(auto1(ntauto))*cos(ebar*t)&
         -imag(auto1(ntauto))*sin(ebar*t))

    h2ij=h2ij+0.5d0*dt*Gk(ntauto)*(real(auto2(ntauto))*cos(ebar*t)&
         -imag(auto2(ntauto))*sin(ebar*t))
    
    sij=sij+0.5d0*dt*Gk(ntauto)*(real(auto(ntauto))*cos(ebar*t) &
         -imag(auto(ntauto))*sin(ebar*t))

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
! calc_ener: calculates the energies using various variational principles
!######################################################################
  
  subroutine calc_ener

    use constants
    use filtermod
    
    implicit none

    integer                         :: i
    real(d)                         :: sigma,fac1,fac2
    real(d), dimension(nener,nener) :: hshift,hshift2
    
!----------------------------------------------------------------------
! Variational principle I:
!----------------------------------------------------------------------
! Omega = < Psi | H | Psi > / < Psi | Psi > 
!----------------------------------------------------------------------
    call solve_geneig(hfbas,sfbas,rvec,rener,transmat,nener,nrbas)

!----------------------------------------------------------------------
! Representation of H and H^2 in th reduced state basis
!----------------------------------------------------------------------
    ! < psi_i | H | psi_j >
    allocate(hrbas(nrbas,nrbas))
    hrbas=matmul(transpose(transmat),matmul(hfbas,transmat))
    
    ! < psi_i | H^2 | psi_j >
    allocate(h2rbas(nrbas,nrbas))
    h2rbas=matmul(transpose(transmat),matmul(h2fbas,transmat))

    
    return
    
  end subroutine calc_ener
    
!######################################################################
! solve_geneig: solves the generalised eigenvalue problem
!               A V = B V E
!               For non-positive-definite matrices B, null space
!               vectors are discarded.
!######################################################################
  
  subroutine solve_geneig(A,B,Vbar,Ebar,P,matdim,rdim)

    use constants
    use iomod
    use filtermod, only: ovrthrsh
    
    implicit none

    integer                              :: matdim,workdim,rdim,&
                                            error,i
    real(d), dimension(matdim,matdim)    :: A,B

    real(d), parameter                   :: thrsh=1e+2_d
    
    real(d), dimension(matdim,matdim)    :: U
    real(d), dimension(matdim)           :: lambda
    real(d), dimension(:,:), allocatable :: Ubar,normfac,P,Abar,Vbar
    real(d), dimension(:), allocatable   :: Ebar
    real(d), dimension(:), allocatable   :: work

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
    rdim=0.0d0
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

    use constants
    use filtermod
    
    implicit none

    integer :: j
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    allocate(dvec(nener))
    dvec=0.0d0

    allocate(avec(nrbas))
    avec=0.0d0

    allocate(intens(nrbas))
    intens=0.0d0
    
!----------------------------------------------------------------------
! Calculate the d-vector, d_j = < Psi(0) | psi_Ej >, where Psi(0)
! is the intitial wavefunction, and Psi_Ej is the jth filter state
!----------------------------------------------------------------------
    call calc_dvec

!----------------------------------------------------------------------
! Calculate the a-vector, a_j = < psi(0) | varphi_j >, where psi(0)
! is the intitial wavefunction, and varphi_j is the jth filter state
!----------------------------------------------------------------------
    avec=matmul(transpose(rvec),matmul(transpose(transmat),dvec))

!----------------------------------------------------------------------
! Calculate the intensities, I_j = E_j a_j**2
! Note that this differs from the definition in the Beck-Meyer
! paper as we account for the dependence on the excitation energy
!----------------------------------------------------------------------
    do j=1,nrbas
       intens(j)=rener(j)*avec(j)**2
    enddo
       
    return
    
  end subroutine calc_intens

!######################################################################
  
  subroutine calc_dvec

    use constants
    use filtermod
    
    implicit none

    integer :: L,i,j

    L=nener

!----------------------------------------------------------------------
! Calculation of the d-vector, d_j = < psi(0) | phi_j >
!----------------------------------------------------------------------
    ! Loop over the energy grid points
    do j=1,L

       ! Calculate the intensity for the current grid point using
       ! the trapezoidal rule
       call calc_dvec_1state(j,dvec(j))
       
    enddo
       
    return
    
  end subroutine calc_dvec

!######################################################################

  subroutine calc_dvec_1state(j,dj)

    use constants
    use filtermod
    
    implicit none

    integer :: j,ti
    real(d) :: dj,t,ener

!----------------------------------------------------------------------
! E_j
!----------------------------------------------------------------------
    ener=(j-1)*de

!----------------------------------------------------------------------
! Contribution from the first timestep
!----------------------------------------------------------------------
    dj=dt*windowfunc(0.0d0)*real(auto(1))

!----------------------------------------------------------------------
! Contribution from timesteps 2 to ntauto-1
!----------------------------------------------------------------------
    do ti=2,ntprop-1

       t=(ti-1)*dt
       
       dj=dj+2.0d0*dt*windowfunc(t) &
            *(real(auto(ti))*cos(ener*t)-imag(auto(ti)*sin(ener*t)))
       
    enddo

!----------------------------------------------------------------------
! Contribution from the final timestep
!----------------------------------------------------------------------
    t=(ntprop-1)*dt

    dj=dj+dt*windowfunc(t)*(real(auto(ntprop))*cos(ener*t) &
         -imag(auto(ntprop)*sin(ener*t)))
       
    
    return
    
  end subroutine calc_dvec_1state

!######################################################################
! windowfunc: for the time t, this function returns the value of the
!             window function g_k(t)
!######################################################################
  
  function windowfunc(t) result(func)

    use constants
    use filtermod
    
    implicit none

    real(d) :: t,func

    select case(iwfunc)

!----------------------------------------------------------------------
! cos0 window function
!----------------------------------------------------------------------
    case(0)
       
       func=1.0d0

       return
       
!----------------------------------------------------------------------
! cos1 window function
!----------------------------------------------------------------------
    case(1)

       func=cos( (pi*t)/(2.0d0*proptime) )
       
       return

!----------------------------------------------------------------------
! cos2 window function
!----------------------------------------------------------------------
    case(2)

       func=(cos((pi*t)/(2.0d0*proptime)))**2
       
       return
       
    end select
    
  end function windowfunc
    
!######################################################################
! calc_errors: calculation of error estimates for the computed
!              energies.
!######################################################################
  
  subroutine calc_errors

    use constants
    use filtermod
    
    implicit none

    integer                         :: i,j
    real(d), dimension(nrbas,nrbas) :: tmpmat
    real(d), dimension(nrbas,nrbas) :: e2diag

    real(d), dimension(nrbas)       :: tmpvec
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(error(nrbas))

!----------------------------------------------------------------------
! Calculate the error estimates
!----------------------------------------------------------------------
    do i=1,nrbas

       ! OLD: Original estimation suggested by Neuhauser
       !e2diag=0.0d0
       !do j=1,nrbas
       !   e2diag(j,j)=rener(i)**2
       !enddo
       !
       !tmpmat=matmul(transpose(rvec),matmul(h2rbas-e2diag,rvec))
       !
       !error(i)=sqrt(abs(tmpmat(i,i)))

       ! NEW: Estimation suggested in:
       ! Equation 9, Werner and Cary, J. Comp. Phys., 227, 5200 (2008)
       tmpvec=matmul(hrbas,rvec(:,i))-rener(i)*rvec(:,i)
       error(i)=sqrt(dot_product(tmpvec,tmpvec))
       tmpvec=matmul(hrbas,rvec(:,i))
       error(i)=error(i)/sqrt(dot_product(tmpvec,tmpvec))
       
    enddo
    
    return
    
  end subroutine calc_errors

!######################################################################

  subroutine wrout

    use channels
    use iomod
    use constants
    use filtermod
    
    implicit none

    integer          :: i,j,iout
    real(d)          :: maxi
    character(len=1) :: atmp
    
!----------------------------------------------------------------------
! Write all calculated energies, intensities and errors to the log
! file
!----------------------------------------------------------------------
    write(ilog,'(/,61a)') ('#',i=1,61)
    write(ilog,'(a)') '# j    E_j (eV)        I_j        &
         Error            Spurious?'
    write(ilog,'(61a)') ('#',i=1,61)

    do j=1,nrbas
       
       if (error(j).lt.errthrsh) then
          atmp='N'
       else
          atmp='Y'
       endif
          
       write(ilog,'(1x,i3,2x,F12.7,2x,F12.7,2x,E11.5,6x,a1)') j,&
            rener(j)*eh2ev,intens(j),error(j),atmp
    enddo

!----------------------------------------------------------------------
! Write the non-spurious energies and intensities to file
!----------------------------------------------------------------------
    ! Maximum intensity in the energy range of interest
    maxi=0.0d0
    do j=1,nrbas
       if (rener(j).lt.ebound(1).or.rener(j).gt.ebound(2)) cycle
       if (intens(j).gt.maxi) maxi=intens(j)
    enddo

    ! Open the output file
    call freeunit(iout)
    open(iout,file='fdiag_eig.dat',form='formatted',status='unknown')

    ! Table header
    write(iout,'(33a)') ('#',i=1,33)
    write(iout,'(a)') '#  Energy        Intensity        &
         Normalised Intensity'
    write(iout,'(33a)') ('#',i=1,33)

    ! Output the energies, intensities and normalised intensities
    ! for the non-spurious states in the energy range of interest
    do j=1,nrbas

       if (rener(j).lt.ebound(1).or.rener(j).gt.ebound(2)) cycle

       if (error(j).lt.errthrsh) &
            write(iout,'(3(2x,F12.7))') rener(j)*eh2ev,intens(j),&
            intens(j)/maxi
    enddo

    ! Close the output file
    close(iout)
    
    return
    
  end subroutine wrout
    
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
    deallocate(transmat)
    deallocate(rvec)
    deallocate(rener)
    deallocate(dvec)
    deallocate(avec)
    deallocate(intens)
    deallocate(error)
    
!----------------------------------------------------------------------
! Close files
!----------------------------------------------------------------------
    close(iin)
    close(ilog)
    
    return
    
  end subroutine fdiag_finalise
    
!######################################################################
  
end program fdiag
