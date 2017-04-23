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
! Calculate the Hamiltonian matrix elements in the basis of the
! filter states
!----------------------------------------------------------------------
  call calcham_fstatebas
  
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
         'Energy window (eV): ',ebound(1)*eh2ev,ebound(2)*eh2ev
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

  subroutine calcham_fstatebas

    use constants
    use filtermod
    
    implicit none

    print*,"Write the subroutine calcham_fstatebas"
    stop
    
    return
    
  end subroutine calcham_fstatebas
    
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

!----------------------------------------------------------------------
! Close files
!----------------------------------------------------------------------
    close(iin)
    close(ilog)
    
    return
    
  end subroutine fdiag_finalise
    
!######################################################################
  
end program fdiag
