!#######################################################################
! Calculation of filter diagonalisation eigenstates from either:
!
! (1) A wavepacket propagation (TDFD), or;
! (2) An order-domain Chebyshev propagation (CFD).
!
! Note that at the ADC(2) level, we only need the 1h1p part of the
! eigenvectors to calculate NTOs, and so for CFD we also allow for
! the reading of the 1h1p parts of the Chebyshev order-domain vectors
! from file. These can then be used to construct the 1h1p parts of the
! eigenvectors, and then the NTOs. This affords a massive saving of
! computational effort as we don't have to repeat the Chebyshev
! propagation.
!#######################################################################

module adc2fdstatesmod

  use channels

contains

!#######################################################################
  
  subroutine adc2_fdstates(gam)

    use constants
    use parameters
    use adc_common
    use fspace
    use mp2
    use import_gamess
    
    implicit none
    
    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: ndim,ndims,ndimf,&
                                            ndimsf,nout,noutf,ndimd
    integer*8                            :: noffd,noffdf
    real(dp)                             :: e0
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
! Read the data file generated by either the FDIAG or CHEBYFD program
!-----------------------------------------------------------------------
    call rddatfile

!-----------------------------------------------------------------------
! Read the eigenstate selection file
!-----------------------------------------------------------------------
    call rdselfile

!-----------------------------------------------------------------------
! If only the 1h1p parts of the eigenvectors are to be calculated, the
! required basis vectors will be read from disk. Otherwise, the
! basis vectors will be recalculated.
!-----------------------------------------------------------------------
    if (read1h1p) then
       call fdstates_1h1p_read(ndimsf,ndimf,kpqf,gam)
    else
       call fdstates_recalc(ndim,ndims,ndimf,kpq,kpqf,noffd,noffdf,gam)
    endif
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    if (allocated(fbas2eig)) deallocate(fbas2eig)
    if (allocated(k2eig)) deallocate(k2eig)
    
    return
    
  end subroutine adc2_fdstates

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

    use parameters
    
    implicit none

    if (autoprop.eq.1) then
       ! TDFD calculation
       call rddatfile_tdfd
    else if (autoprop.eq.2) then
       ! CFD calculation
       call rddatfile_cfd
    endif
    
    return
    
  end subroutine rddatfile
    
!#######################################################################
  
  subroutine rddatfile_tdfd

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
    
  end subroutine rddatfile_tdfd
    
!#######################################################################

  subroutine rddatfile_cfd

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
! No. eigenstates
!-----------------------------------------------------------------------
    read(unit) neig

!-----------------------------------------------------------------------
! Chebyshev expansion order
!-----------------------------------------------------------------------
    read(unit) chebyord

!-----------------------------------------------------------------------
! Transformation matrix
!-----------------------------------------------------------------------
    ! Allocate k2eig
    allocate(k2eig(0:chebyord,1:neig))

    ! Read in k2eig
    read(unit) k2eig(0:chebyord,1:neig)
    
!-----------------------------------------------------------------------
! Spectral bounds (in a.u.)
!-----------------------------------------------------------------------
    read(unit) ebound

!-----------------------------------------------------------------------
! Close the fdiag data file
!-----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rddatfile_cfd
    
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
       
       if (keyword(inkw).eq.'*') then
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

  subroutine fdstates_recalc(ndim,ndims,ndimf,kpq,kpqf,noffd,noffdf,gam)

    use constants
    use parameters
    use adc_common
    use nto
    use import_gamess
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndims,ndimf
    integer*8                                 :: noffd,noffdf
    real(dp), dimension(:), allocatable       :: dpsi
    type(gam_structure)                       :: gam
    
!-----------------------------------------------------------------------
! Set MO representation of the dipole operator
!-----------------------------------------------------------------------
    call set_dpl

!-----------------------------------------------------------------------
! Contraction of the dipole operator with the initial state
!-----------------------------------------------------------------------
    call contract_dipole_initial_state(dpsi,ndim,ndims,ndimf,kpq,kpqf,&
         noffd)

!-----------------------------------------------------------------------
! Calculate the Hamiltonian matrix
!-----------------------------------------------------------------------
    call calc_hamiltonian(ndimf,kpqf,noffdf)
    
!-----------------------------------------------------------------------
! Calculation of the filter diagonalisation eigenstates
!-----------------------------------------------------------------------
    call calc_fdstates(dpsi,ndimf,noffdf)

!-----------------------------------------------------------------------
! Output the results of the calculation
!-----------------------------------------------------------------------
    call wrfdstates(kpqf,ndimf)

!-----------------------------------------------------------------------
! If requested, calculate and output NTOs
!-----------------------------------------------------------------------
    if (lnto) call adc2_nto(gam,ndimf,kpqf,'SCRATCH/fdstates',nsel,&
         'nto')
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(dpsi)
    
    return
    
  end subroutine fdstates_recalc

!#######################################################################

  subroutine fdstates_1h1p_read(ndimsf,ndimf,kpqf,gam)

    use constants
    use parameters
    use iomod
    use nto
    use import_gamess
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimsf,ndimf
    integer                                   :: unit,k,i
    real(dp), dimension(:,:), allocatable     :: q1h1p,eig1h1p
    real(dp), dimension(:), allocatable       :: eigtmp
    type(gam_structure)                       :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(q1h1p(ndimsf,0:chebyord))
    q1h1p=0.0d0

    allocate(eig1h1p(ndimsf,nsel))
    eig1h1p=0.0d0

    allocate(eigtmp(ndimf))
    eigtmp=0.0d0
    
!----------------------------------------------------------------------
! Read the 1h1p parts of the Chebyshev order-domain vectors from file
!----------------------------------------------------------------------
    ! Open the cheby1h1p file
    call freeunit(unit)
    open(unit=unit,file='cheby1h1p',status='unknown',&
         access='sequential',form='unformatted')

    ! Read the 1h1p parts of the Chebyshev vectors
    do k=0,chebyord
       read(unit) q1h1p(:,k)
    enddo

    ! Close the cheby1h1p file
    close(unit)

!----------------------------------------------------------------------
! Calculate the 1h1p parts of the eigenvectors
!----------------------------------------------------------------------
    eig1h1p=0.0d0
    do i=1,nsel
       do k=0,chebyord
          eig1h1p(:,i)=eig1h1p(:,i)+k2eig(k,isel(i))*q1h1p(:,k)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Write the 1h1p parts of the eigenvectors to file padded with zeros
! in the 2h2p part. This must be done before calling adc2_nto
!----------------------------------------------------------------------
    ! Open the fdstates file
    call freeunit(unit)
    open(unit=unit,file='SCRATCH/fdstates',status='unknown',&
         access='sequential',form='unformatted')
    
    ! Write the 1h1p parts of the eigenstates to file
    eigtmp=0.0d0
    do i=1,nsel
       eigtmp(1:ndimsf)=eig1h1p(1:ndimsf,i)
       write(unit) i,1.0d0,eigtmp
    enddo

    ! Close the fdstates file
    close(unit)
    
!----------------------------------------------------------------------
! Calculate and output the NTOs
!----------------------------------------------------------------------
    call adc2_nto(gam,ndimf,kpqf,'SCRATCH/fdstates',nsel,'nto')
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(q1h1p)
    deallocate(eig1h1p)
    deallocate(eigtmp)
    
    return
    
999 continue
    errmsg='The cheby1h1p file does not exist'
    call error_control
    
  end subroutine fdstates_1h1p_read
  
!#######################################################################

  subroutine calc_fdstates(dpsi,ndimf,noffdf)

    use fdstates_tdfd
    use fdstates_cfd
    
    implicit none

    integer                    :: ndimf
    integer*8                  :: noffdf
    real(dp), dimension(ndimf) :: dpsi
    
    hamflag='f'

    if (autoprop.eq.1) then
       ! TDFD
       call calc_fdstates_tdfd(dpsi,ndimf,noffdf)
    else if (autoprop.eq.2) then
       ! CFD
       call calc_fdstates_cfd(dpsi,ndimf,noffdf)
    endif
       
    return
    
  end subroutine calc_fdstates
  
!#######################################################################

  subroutine wrfdstates(kpqf,ndimf)

    use constants
    use parameters
    use misc
    use iomod
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf,i,itmp
    character(len=1), dimension(2)            :: am
    character(len=36)                         :: filename

!----------------------------------------------------------------------
! Write the filter diagonalisation eigenstate information to the
! log file
!----------------------------------------------------------------------
    am(1:2)=(/ 's','x' /)
    
    write(ilog,'(/,70a)') ('*',i=1,70)
    write(ilog,'(2x,a)') &
         'Filter diagonalisation ADC(2)-'//am(abs(method)-1)&
         //' excitation energies'
    write(ilog,'(70a)') ('*',i=1,70)
    
    itmp=1+nBas**2*4*nOcc**2

    filename='SCRATCH/fdstates'

    call wrstateinfo_neutral(ndimf,kpqf,itmp,filename,nsel)
    
    write(ilog,'(/,70a,/)') ('*',i=1,70)

    return
    
  end subroutine wrfdstates

!#######################################################################

  subroutine contract_dipole_initial_state(dpsi,ndim,ndims,ndimf,&
       kpq,kpqf,noffd)
    
    use constants
    use parameters
    use get_moment
    use adc_common
    use get_matrix_dipole
    use guessvecs
    use misc
    
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: i
    integer                                   :: ndim,ndims,ndimf,itmp
    integer*8                                 :: noffd
    real(dp), dimension(:), allocatable       :: dpsi,ener,mtm,tmvec,&
                                                 osc_str
    real(dp), dimension(:,:), allocatable     :: rvec
    real(dp)                                  :: time
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(dpsi(ndimf))
    
!-----------------------------------------------------------------------
! Initial state = | Psi_0 >
! Calculate the f-vector, F_J = < Psi_J | D | Psi_0 >    
!-----------------------------------------------------------------------
    if (statenumber.eq.0) call get_modifiedtm_adc2(ndimf,kpqf(:,:),&
         dpsi(:),1)

!-----------------------------------------------------------------------
! Excited initial state: diagonalisation in the initial space and
! contraction of the D-matrix with the initial state vector
!-----------------------------------------------------------------------
    if (statenumber.gt.0) then

       ! (1) Initial space diagonalisation
       if (ladc1guess) call adc1_guessvecs
       call initial_space_diag(time,kpq,ndim,ndims,noffd)
       call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,kpq)
       
       ! (2) Output the initial space state information
       write(ilog,'(/,70a)') ('*',i=1,70)
       write(ilog,'(2x,a)') &
            'Initial space ADC(2)-s excitation energies'
       write(ilog,'(70a)') ('*',i=1,70)
       itmp=1+nBas**2*4*nOcc**2
       call table2(ndim,davstates,ener(1:davstates),&
            rvec(:,1:davstates),tmvec(1:davstates),&
            osc_str(1:davstates),kpq,itmp,'i')
       write(ilog,'(/,70a,/)') ('*',i=1,70)
       
       ! (3) Contraction of the D-matrix with the initial state vector
       call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,&
            rvec(:,statenumber),dpsi)

       ! Deallocate arrays that are no longer needed (N.B. these
       ! are allocated in initial_space_tdm)
       deallocate(ener,rvec,tmvec,osc_str)

    endif
    
    return
    
  end subroutine contract_dipole_initial_state
  
!#######################################################################
  
end module adc2fdstatesmod
