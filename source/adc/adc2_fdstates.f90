!#######################################################################
! Calculation of filter diagonalisation eigenstates from the
! propagation of |Psi(t=0)> = mu |Psi_0>
!#######################################################################

module adc2fdstatesmod

  use channels

contains

!#######################################################################

  subroutine adc2_fdstates(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use mp2
    use fdstates
    use import_gamess
    
    implicit none
    
    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimf,&
                                            ndimsf,nout,noutf,ndimd
    integer*8                            :: noffd,noffdf
    real(d)                              :: e0,time
    real(d), dimension(:), allocatable   :: fvec
    type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!-----------------------------------------------------------------------
    call mp2_master(e0)

!-----------------------------------------------------------------------
! Determine the 1h1p and 2h2p subspaces
!-----------------------------------------------------------------------
    call get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf,&
         ndims,ndimsf)

!-----------------------------------------------------------------------
! Set MO representation of the dipole operator
!-----------------------------------------------------------------------
    call set_dpl
    
!-----------------------------------------------------------------------
! Calculate the f-vector, F_J = < Psi_J | D | Psi_0 >    
!-----------------------------------------------------------------------
    call fvector(fvec,ndimf,kpqf)

!-----------------------------------------------------------------------
! Calculate the Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian(ndimf,kpqf,noffdf)

!-----------------------------------------------------------------------
! Read the filter state-to-eigenstate transformation matrix, the
! window function type and the energy window from the fdiag dat file
!-----------------------------------------------------------------------
    call rddatfile

!-----------------------------------------------------------------------
! Read the eigenstate selection file
!-----------------------------------------------------------------------
    call rdselfile
    
!-----------------------------------------------------------------------
! Perform the wavepacket propagation and the calculation of the
! filter diagonalisation eigenstates
!-----------------------------------------------------------------------
    hamflag='f'
    call calc_fdstates(fvec,ndimf)

!-----------------------------------------------------------------------
! Output the results of the calculation
!-----------------------------------------------------------------------
    call wrfdstates(kpqf,ndimf)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(fvec)
    
    return
    
  end subroutine adc2_fdstates

!#######################################################################

  subroutine fvector(fvec,ndimf,kpqf)

    use constants
    use parameters
    use get_moment
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    real(d), dimension(:), allocatable        :: fvec
    
!-----------------------------------------------------------------------
! Allocate the F-vector array
!-----------------------------------------------------------------------
    allocate(fvec(ndimf))

!-----------------------------------------------------------------------
! Calculate the F-vector, F_J = < Psi_J | D | Psi_0 >
!-----------------------------------------------------------------------
    call get_modifiedtm_adc2(ndimf,kpqf(:,:),fvec(:),1)
    
    return
    
  end subroutine fvector
    
!#######################################################################

  subroutine calc_hamiltonian(ndimf,kpqf,noffdf)

    use fspace
    
    implicit none
    
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf
    integer*8                                 :: noffdf
    
    write(ilog,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
    
    if (method.eq.2) then
       ! ADC(2)-s
       if (lcvsfinal) then
          call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
       else
          call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
       endif
    else if (method.eq.3) then
       ! ADC(2)-x
       if (lcvsfinal) then
          call write_fspace_adc2e_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
       else
          call write_fspace_adc2e_1(ndimf,kpqf(:,:),noffdf,'c')
       endif
    endif
    
    return
    
  end subroutine calc_hamiltonian

!#######################################################################

  subroutine rddatfile

    use constants
    use parameters
    use iomod

    implicit none
    
    integer :: unit

!-----------------------------------------------------------------------
! Open the fdiag data file
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=fdiagdat,form='unformatted',status='old')

!-----------------------------------------------------------------------
! Read the basis dimensions and allocate arrays
!-----------------------------------------------------------------------
    ! No. filter states
    read(unit) nfbas

    ! No. eigenstates
    read(unit) neig

    ! Allocate the filter state-to-eigenstate transformation matrix
    allocate(fbas2eig(neig,nfbas))
    fbas2eig=0.0d0

!-----------------------------------------------------------------------
! Read the filter state-to-eigenstate transformation matrix
!-----------------------------------------------------------------------
    read(unit) fbas2eig

!-----------------------------------------------------------------------
! Window function type
!-----------------------------------------------------------------------
    read(unit) iwfunc

!----------------------------------------------------------------------
! Energy window (in a.u.)
!----------------------------------------------------------------------
    read(unit) ebound

!-----------------------------------------------------------------------
! Close the fdiag data file
!-----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rddatfile
    
!#######################################################################

  subroutine rdselfile

    use constants
    use parameters
    use iomod
    use parsemod
    
    implicit none

    integer                  :: unit,i,count
    logical, dimension(neig) :: selected
    
!-----------------------------------------------------------------------
! Open the eigenstate selection file
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=fdiagsel,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the eigenstate selection file
!-----------------------------------------------------------------------
    nsel=0
    selected=.false.
    count=0
    
5   call rdinp(unit)

    if (.not.lend) then

       count=count+1
       
       if (keyword(4).eq.'*') then
          nsel=nsel+1
          selected(count)=.true.
       endif

       goto 5

    endif

!-----------------------------------------------------------------------
! Allocate and fill in the isel array
!-----------------------------------------------------------------------
    allocate(isel(nsel))
    isel=0.0d0

    count=0
    do i=1,neig
       if (selected(i)) then
          count=count+1
          isel(count)=i
       endif
    enddo
    
!-----------------------------------------------------------------------
! Close the eigenstate selection file
!-----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rdselfile

!#######################################################################

  subroutine wrfdstates(kpqf,ndimf)

    use constants
    use parameters
    use misc
    use iomod
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf,unit,i,k,itmp
    real(d), dimension(:,:), allocatable      :: fdstates
    real(d), dimension(nsel)                  :: ener,tmp1,tmp2
    character(len=1), dimension(2)            :: am

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(fdstates(ndimf,nsel))
    fdstates=0.0d0
    
!----------------------------------------------------------------------
! Read the filter diagonalisation eigenstates from file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit=unit,file='SCRATCH/fdstates',status='unknown',&
            access='sequential',form='unformatted')
    do i=1,nsel
       read(unit) k,ener(k),fdstates(:,k)
    enddo
    close(unit)

!----------------------------------------------------------------------
! Write the filter diagonalisation eigenstate information to file
!----------------------------------------------------------------------
    am(1:2)=(/ 's','x' /)
    
    write(ilog,'(/,70a)') ('*',i=1,70)
    write(ilog,'(2x,a)') &
         'Initial space ADC(2)-'//am(abs(method)-1)&
         //' excitation energies'
    write(ilog,'(70a)') ('*',i=1,70)

    ! Temporary zero arrays standing in for osc and travec
    tmp1=0.0d0
    tmp2=0.0d0
    
    itmp=1+nBas**2*4*nOcc**2
    call table2(ndimf,nsel,ener(1:nsel),fdstates(:,1:nsel),tmp1,&
         tmp2,kpqf,itmp,'i')
    
    write(ilog,'(/,70a,/)') ('*',i=1,70)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(fdstates)
    
    return
    
  end subroutine wrfdstates
    
!#######################################################################
  
end module adc2fdstatesmod
