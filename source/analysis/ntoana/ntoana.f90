!######################################################################
! ntoana: A simple program for the extraction of time-dependent
!         natural transition orbitals (TD-NTOs) and time-dependent
!         particle/hole densities (TD-PHDs) from a data file produced
!         by the ADC code.
!######################################################################

module global

  use constants
  use import_gamess
  
  implicit none

  save
  
  integer                       :: task
  integer                       :: nao,nocc,nvirt,nbas
  real(dp), allocatable         :: xcoo(:)
  real(dp)                      :: tplt
  real(dp)                      :: dt
  real(dp), allocatable         :: sigma(:)
  real(dp), allocatable         :: val(:)
  real(dp), allocatable         :: occ(:)
  complex(dp), allocatable      :: hole(:,:)
  complex(dp), allocatable      :: particle(:,:)
  complex(dp), allocatable      :: orb(:,:)
  character(len=2), allocatable :: aatm(:)
  logical                       :: pltall
  logical                       :: renorm
  type(gam_structure)           :: gam
    
end module global

!######################################################################

program ntoana

  implicit none

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
  call initialise_ntoano

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdntoanainp

!----------------------------------------------------------------------
! Read the td-nto.dat file and write the plotting files
!----------------------------------------------------------------------
  call make_plot
  
!----------------------------------------------------------------------
! Finalisation
!----------------------------------------------------------------------
  call finalise_ntoana
  
contains

!######################################################################

  subroutine initialise_ntoano

    use global
    use channels
    use iomod
    use accuracy
    
    implicit none

    logical :: found
    
!----------------------------------------------------------------------
! I/O channels
!----------------------------------------------------------------------
    iin=1
    ilog=2

!----------------------------------------------------------------------
! Open the log file
!----------------------------------------------------------------------
    open(ilog,file='ntoana.log',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Open the td-nto.dat file
!----------------------------------------------------------------------
    open(iin,file='td-nto.dat',form='unformatted',status='old',err=999)
    
!----------------------------------------------------------------------
! Initialisation of multigrid constants and kinds
!----------------------------------------------------------------------
    call accuracyInitialize

!----------------------------------------------------------------------
! Read the gamess.dat checkpoint file
!----------------------------------------------------------------------
    ! Make sure that the gamess.dat file exists
    inquire(file='gamess.dat',exist=found)
    if (.not.found) then
       errmsg='The gamess.dat file is missing. Quitting.'
    endif

    ! Read in the basis and geometry
    call gamess_load_orbitals(file='gamess.dat',structure=gam)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Read the stepsize basis dimensions from the td-nto.dat file
    read(iin) dt
    read(iin) nao,nocc,nvirt
    nbas=nocc+nvirt
    
    ! Allocate and initialise arrays
    allocate(hole(nao,nocc))
    hole=czero
    allocate(particle(nao,nvirt))
    particle=czero
    allocate(sigma(nocc))
    sigma=0.0d0
    allocate(orb(nao,nbas))
    orb=czero
    allocate(val(nbas))
    val=0.0d0
    allocate(occ(nbas))
    occ=0.0d0
    
    return

!----------------------------------------------------------------------
! Missing data file: exit
!----------------------------------------------------------------------
999 continue
    errmsg='The td-nto.dat data file does not exist. Quitting.'
    call error_control
    
  end subroutine initialise_ntoano

!######################################################################

  subroutine finalise_ntoana

    use global
    use channels
    
!----------------------------------------------------------------------
! Close files
!----------------------------------------------------------------------
    close(iin)
    close(ilog)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(hole)
    deallocate(particle)
    deallocate(sigma)
    deallocate(orb)
    deallocate(val)
    deallocate(occ)

    return
    
  end subroutine finalise_ntoana

!######################################################################

  subroutine rdntoanainp

    use global
    use constants
    use iomod
    
    implicit none

    integer           :: i
    character(len=20) :: string1,string2
    
!----------------------------------------------------------------------
! Set default values
!----------------------------------------------------------------------
    ! Plotting time.
    tplt=-1.0d0

    ! Plot for all times
    pltall=.false.
    
    ! Task: extraction of TD-NTOs (default) or TD-PHDs
    task=1

    ! Normalised real and imaginary parts of the TD-NTOs
    renorm=.false.    
    
!----------------------------------------------------------------------
! Exit if no argmuents have been given
!----------------------------------------------------------------------
    if (iargc().eq.0) then
       errmsg='No argmuents given. Quitting.'
       call error_control
    endif

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    i=0

10  continue

    i=i+1

    call getarg(i,string1)

    if (string1.eq.'-t') then

       ! Time at which we want the TD-NTO or TD-PHD
       i=i+1
       call getarg(i,string2)

       if (string2.eq.'all') then
          ! Plot for all times
          pltall=.true.
       else
          ! Plot for a single time
          read(string2,*) tplt
       endif
       
    else if (string1.eq.'-nto') then

       ! Extraction of TD-NTOs
       task=1

    else if (string1.eq.'-phd') then

       ! Extraction of TD-PHDs
       task=2

    else if (string1.eq.'-renorm') then

       ! Output normalised real and imaginary parts of the TD-NTOs
       renorm=.true.
       
    else

       ! Unknown keyword
       errmsg='Unknown keyword: '//trim(string1)
       call error_control

    endif

    ! Continue reading any remaining command line arguments
    if (i.lt.iargc()) goto 10

!----------------------------------------------------------------------
! Make sure that all required information has been given
!----------------------------------------------------------------------
    if (.not.pltall.and.tplt.eq.-1.0d0) then
       errmsg='The plotting time has not been given'
       call error_control
    endif
    
    return
    
  end subroutine rdntoanainp
    
!######################################################################

  subroutine make_plot

    use global
    use iomod    

    if (pltall) then
       ! Plotting of TD-NTOs and TD-PHDs for all times
       call make_plot_all_times
    else
       ! Plotting of TD-NTOs and TD-PHDs for a single time
       call make_plot_1time
    endif
       
    return
    
  end subroutine make_plot

!######################################################################

   subroutine make_plot_all_times

    use global
    use constants
    use channels
    use iomod

    implicit none

    integer           :: npair,i
    real(dp)          :: t
    character(len=10) :: at

!-----------------------------------------------------------------------
! Output the TD-NTOs or TD-PHDs at each timestep
!-----------------------------------------------------------------------
10  continue

    ! Read the NTOs for the current timestep
    read(iin,end=999) t
    read(iin) npair
    do i=1,npair
       read(iin) hole(:,i),particle(:,i),sigma(i)
    enddo

    ! Write the plotting file for the current timestep
    if (task.eq.1) then
       ! TD-NTOs
       call wrnto_1time(t,npair)
    else if (task.eq.2) then
       ! TD-PHDs
       errmsg='WRITE THE TD-PHD CODE!'
       call error_control
    endif

    ! Continue to the next timestep
    goto 10

999 continue
    
    return
    
  end subroutine make_plot_all_times
    
!######################################################################

  subroutine make_plot_1time

    use global
    use constants
    use channels
    use iomod

    implicit none

    integer           :: npair,i
    real(dp)          :: t
    character(len=10) :: at
    
!-----------------------------------------------------------------------
! Read to the timestep requested by the user
!-----------------------------------------------------------------------
10  continue
    read(iin,end=999) t
    read(iin) npair
    do i=1,npair
       read(iin) hole(:,i),particle(:,i),sigma(i)
    enddo
    if (t.ne.tplt) goto 10
    
!-----------------------------------------------------------------------
! Output the TD-NTOs or TD-PHDs
!-----------------------------------------------------------------------
    if (task.eq.1) then
       ! TD-NTOs
       call wrnto_1time(tplt,npair)
    else if (task.eq.2) then
       ! TD-PHDs
       errmsg='WRITE THE TD-PHD CODE!'
       call error_control
    endif
       
    return

!-----------------------------------------------------------------------
! Exit if we could not find an entry for the requested plotting time
!-----------------------------------------------------------------------
999 continue
    write(at,'(F10.4)') tplt
    errmsg='No entry for time '//trim(adjustl(at))//' could be found'
    call error_control
    
  end subroutine make_plot_1time

!######################################################################

  subroutine wrnto_1time(t,npair)

    use global
    use constants
    use moldenmod
    
    implicit none

    real(dp)          :: t
    integer           :: npair,j
    character(len=70) :: filename
    character(len=10) :: at

!-----------------------------------------------------------------------
! Normalisation of the real and imaginary parts of the NTOs
!-----------------------------------------------------------------------
    if (renorm) then
       do j=1,npair
          hole(:,j)=hole(:,j) &
               /sqrt(dot_product(hole(:,j),hole(:,j)))
          particle(:,j)=particle(:,j) &
               /sqrt(dot_product(particle(:,j),particle(:,j)))
       enddo
    endif
    
!-----------------------------------------------------------------------
! Dominant NTOs and singular values
!-----------------------------------------------------------------------
    orb=czero
    val=0.0d0
    occ=0.0d0
    do j=1,npair
       orb(:,npair-j+1)=hole(:,j)
       orb(:,npair+j)=particle(:,j)
       val(npair-j+1)=-0.5d0*sigma(j)**2
       val(npair+j)=0.5d0*sigma(j)**2
    enddo
    occ(1:npair)=1.0d0
    
!-----------------------------------------------------------------------
! Write the molden files: one for the real parts of the NTOs and one
! for the imaginary parts
!-----------------------------------------------------------------------
    write(at,'(F10.4)') t

    ! Real part
    filename='nto_'//trim(adjustl(at))//'_real.molden'
    call write_molden(gam,filename,nao,2*npair,&
         real(orb(1:nao,1:2*npair)),val(1:2*npair),occ(1:2*npair))

    ! Imaginary part
    filename='nto_'//trim(adjustl(at))//'_imag.molden'
    call write_molden(gam,filename,nao,2*npair,&
         aimag(orb(1:nao,1:2*npair)),val(1:2*npair),occ(1:2*npair))

    return
    
  end subroutine wrnto_1time

!######################################################################
  
end program ntoana
