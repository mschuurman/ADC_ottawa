!######################################################################
! nto: subroutines for the calculation of natural transition orbitals
!######################################################################

module nto
  
  use constants
  use parameters
  use iomod
  use channels
  
  implicit none

  save

  private :: dp
    
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter :: dp=selected_real_kind(8)
  
contains

!######################################################################
! adc2_nto: gateway subroutine for the calculation of ADC(2) NTOs
!######################################################################
  
  subroutine adc2_nto(gam,ndimf,kpqf,vecfile,nstates)

    use gamess_internal
    
    implicit none

    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i
    character(len=*)                          :: vecfile
    type(gam_structure)                       :: gam
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',i=1,70)
    write (ilog,'(2x,a)') 'Calculating natural transition orbitals'
    write(ilog,'(70a,/)') ('-',i=1,70)
    
!----------------------------------------------------------------------
! Calculate the NTOs
!----------------------------------------------------------------------
    if (statenumber.eq.0) then
       call adc2_nto_gs(gam,ndimf,kpqf,vecfile,nstates)
    else
       errmsg='WRITE THE EXCITED STATE-TO-EXCITED STATE NTO CODE!'
       call error_control
    endif

    return
    
  end subroutine adc2_nto

!######################################################################
! adc2_nto_gs: calculation of ADC(2) NTOs for excitation from the
!              ground state
!######################################################################
  
  subroutine adc2_nto_gs(gam,ndimf,kpqf,vecfile,nstates)

    use gamess_internal
    use diagmod, only: readdavvc
    use density_matrix
    use moldenmod
    
    implicit none
    
    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i,j,aa,npair
    integer                                   :: nao,lwork,ierr
    real(dp), allocatable                     :: rvec(:,:)
    real(dp), allocatable                     :: trdens(:,:,:)
    real(dp), allocatable                     :: tmp(:,:)
    real(dp), allocatable                     :: sigma(:,:)
    real(dp), allocatable                     :: VT(:,:,:)
    real(dp), allocatable                     :: U(:,:,:)
    real(dp), allocatable                     :: work(:)
    real(dp), allocatable                     :: hole(:,:)
    real(dp), allocatable                     :: particle(:,:)
    real(dp), allocatable                     :: orb(:,:)
    real(dp), allocatable                     :: val(:)
    real(dp), allocatable                     :: occ(:)
    real(dp), parameter                       :: thrsh=0.01d0
    character(len=*)                          :: vecfile
    character(len=70)                         :: filename
    character(len=3)                          :: ai
    type(gam_structure)                       :: gam

    nvirt=nbas-nocc
    
!----------------------------------------------------------------------
! Note that we are here assuming that nocc < nvirt, which should
! almost always hold true
!----------------------------------------------------------------------
    if (nocc.gt.nvirt) then
       errmsg='Error in calculating the NTOs: nocc > nvirt'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rvec(ndimf,nstates))
    rvec=0.0d0

    allocate(trdens(nbas,nbas,nstates))
    trdens=0.0d0

    allocate(tmp(nocc,nvirt))
    tmp=0.0d0
    
    allocate(sigma(nocc,nstates))
    sigma=0.0d0
    
    allocate(U(nocc,nocc,nstates))
    U=0.0d0
    
    allocate(VT(nvirt,nvirt,nstates))
    VT=0.0d0
    
    lwork=3*nocc+nvirt
    allocate(work(lwork))
    work=0.0d0

    nao=gam%nbasis
    allocate(hole(nao,nocc))
    allocate(particle(nao,nvirt))
    hole=0.0d0
    particle=0.0d0

    allocate(orb(nao,nbas))
    orb=0.0d0

    allocate(val(nbas))
    val=0.0d0

    allocate(occ(nbas))
    occ=0.0d0
    
!----------------------------------------------------------------------
! Read the ADC(2) vectors from file
!----------------------------------------------------------------------
    call rdvecs(vecfile,rvec,ndimf,nstates)
    
!----------------------------------------------------------------------
! Calculate the ADC(2) ground state-to-excited state transition
! density matrices
!----------------------------------------------------------------------
    call adc2_trden_gs(trdens,ndimf,kpqf,rvec,nstates)

!----------------------------------------------------------------------
! Calculate the SVDs of the occupied-virtual blocks of the ADC(2)
! ground state-to-excited state transition density matrices
!----------------------------------------------------------------------
    do i=1,nstates
       
       tmp=transpose(trdens(nocc+1:nbas,1:nocc,i))
       
       call dgesvd('A','A',nocc,nvirt,tmp,nocc,sigma(:,i),U(:,:,i),&
            nocc,VT(:,:,i),nvirt,work,lwork,ierr)
       
       if (ierr.ne.0) then
          errmsg='SVD of the transition density matrices failed in &
               subroutine adc2_nto_gs'
          call error_control
       endif
       
    enddo

!----------------------------------------------------------------------
! Calculate and output the NTOs
!----------------------------------------------------------------------
    do i=1,nstates

       ! Hole NTOs in terms of the AOs
       hole=transpose(matmul(transpose(U(:,:,i)),&
            transpose(ao2mo(1:nao,1:nocc))))

       ! Particle NTOs in terms of the AOs
       particle=transpose(matmul(VT(:,:,i),&
            transpose(ao2mo(1:nao,nocc+1:nbas))))

       ! No. NTO pairs
       npair=0
       do j=1,nocc
          if (sigma(j,i)**2.ge.thrsh) npair=npair+1
       enddo

       ! Dominant NTOs and singular values
       orb=0.0d0
       val=0.0d0
       occ=0.0d0
       do j=1,npair
          orb(:,npair-j+1)=hole(:,j)
          orb(:,npair+j)=particle(:,j)
          val(npair-j+1)=-0.5d0*sigma(j,i)**2
          val(npair+j)=0.5d0*sigma(j,i)**2
       enddo
       occ(1:npair)=1.0d0
       
       ! Write the molden file
       write(ai,'(i3)') i
       filename='nto_0_'//trim(adjustl(ai))//'.molden'
       call write_molden(gam,filename,nao,2*npair,&
            orb(1:nao,1:2*npair),val(1:2*npair),occ(1:2*npair))

    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rvec)
    deallocate(trdens)
    deallocate(tmp)
    deallocate(sigma)
    deallocate(U)
    deallocate(VT)
    deallocate(work)
    deallocate(hole)
    deallocate(particle)
    deallocate(orb)
    deallocate(val)
    deallocate(occ)
    
    return
    
  end subroutine adc2_nto_gs
    
!######################################################################

  subroutine rdvecs(vecfile,rvec,ndimf,nvec)

    implicit none

    integer                         :: ndimf,nvec
    integer                         :: ivec,tmp,i
    real(dp), dimension(ndimf,nvec) :: rvec
    real(dp)                        :: ener
    character(len=*)                :: vecfile

!-----------------------------------------------------------------------
! Open the vector file
!-----------------------------------------------------------------------
    call freeunit(ivec)
      
    open(unit=ivec,file=vecfile,status='old',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Read the vectors
!-----------------------------------------------------------------------
    do i=1,nvec
       read(ivec) tmp,ener,rvec(:,i)
    enddo

!-----------------------------------------------------------------------
! Close the vector file
!-----------------------------------------------------------------------
    close(ivec)
    
    return
    
  end subroutine rdvecs
    
!######################################################################
  
end module nto
