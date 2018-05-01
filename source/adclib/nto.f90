!######################################################################
! nto: subroutines for the calculation of natural transition orbitals
!######################################################################

module nto
  
  use constants
  use parameters
  use iomod
  use channels
  use timingmod
  
  implicit none

  save

  real(dp), allocatable :: rhomp2(:,:)
  
contains

!######################################################################
! adc2_nto: gateway subroutine for the calculation of ADC(2) NTOs
!######################################################################
  
  subroutine adc2_nto(gam,ndimf,kpqf,vecfile,nstates,stem)

    use gamess_internal
    
    implicit none

    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i
    real(dp)                                  :: tw1,tw2,tc1,tc2
    character(len=*)                          :: vecfile,stem
    type(gam_structure)                       :: gam
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',i=1,70)
    write (ilog,'(2x,a)') 'Calculating natural transition orbitals'
    write(ilog,'(70a,/)') ('-',i=1,70)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Calculate the NTOs
!----------------------------------------------------------------------
    if (statenumber.eq.0) then
       call adc2_nto_gs(gam,ndimf,kpqf,vecfile,nstates,stem)
    else
       errmsg='WRITE THE EXCITED STATE-TO-EXCITED STATE NTO CODE!'
       call error_control
    endif

!----------------------------------------------------------------------
! Finish timing and output the time taken
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(2x,a,2x,F7.2,1x,a1)') 'Time taken:',tw2-tw1,'s'
    
    return
    
  end subroutine adc2_nto

!######################################################################
! adc2_nto_gs: calculation of ADC(2) NTOs for excitation from the
!              ground state
!######################################################################
  
  subroutine adc2_nto_gs(gam,ndimf,kpqf,vecfile,nstates,stem)

    use gamess_internal
    use mp2
    use density_matrix
    use moldenmod
    
    implicit none
    
    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: i,j,npair
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
    character(len=*)                          :: vecfile,stem
    character(len=70)                         :: filename
    character(len=3)                          :: ai
    type(gam_structure)                       :: gam

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
    allocate(rhomp2(nbas,nbas))
    rhomp2=0.0d0

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
! Calculate the MP2 correction to the ground state density matrix
!----------------------------------------------------------------------
    ! MP2 ground state density matrix
    call rho_mp2(rhomp2)
    
    ! Subtraction of the zeroth-order contribution to obtain the
    ! MP2 correction
    do i=1,nocc
       rhomp2(i,i)=rhomp2(i,i)-2.0d0
    enddo
    
!----------------------------------------------------------------------
! Calculate the occupied-virtural block of the ADC(2)
! ground-to-excited state transition density matrices
!----------------------------------------------------------------------
    trdens=0.0d0
    call adc2_trdens_gs_virt_occ(trdens,ndimf,kpqf,rvec,nstates,rhomp2)
    
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
       filename=trim(stem)//'_0_'//trim(adjustl(ai))//'.molden'
       call write_molden(gam,filename,nao,2*npair,&
            orb(1:nao,1:2*npair),val(1:2*npair),occ(1:2*npair))

    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rhomp2)
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
! tdadc2_nto_init: Precalculation of intermediate terms entering into
!                  the time-dependent ADC(2) NTOs 
!######################################################################

  subroutine tdadc2_nto_init

    use mp2
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write (ilog,'(2x,a)') 'Calculating intermediates in the TD-NTO &
         expressions'

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rhomp2(nbas,nbas))
    rhomp2=0.0d0
    
!----------------------------------------------------------------------
! Precalculate the intermediates
!----------------------------------------------------------------------
    ! MP2 ground state density matrix
    call rho_mp2(rhomp2)
    
    ! Subtraction of the zeroth-order contribution to obtain the
    ! MP2 correction
    do i=1,nocc
       rhomp2(i,i)=rhomp2(i,i)-2.0d0
    enddo

    return
    
  end subroutine tdadc2_nto_init

!######################################################################
! tdadc2_nto_finalise: Deallocation of intermediate term arrays used
!                      in the calculation of time-depdendent ADC(2) 
!                      NTOs
!######################################################################

  subroutine tdadc2_nto_finalise

    implicit none

!----------------------------------------------------------------------
! Deallocate arrays    
!----------------------------------------------------------------------
    deallocate(rhomp2)
    
    return
    
  end subroutine tdadc2_nto_finalise

!######################################################################
! tdadc2_nto: Calculation of the complex, time-dependent ADC(2) NTOs
!######################################################################
  
  subroutine tdadc2_nto(gam,wf,ndimf,kpqf,stem)
    
    use gamess_internal
    use density_matrix
    use moldenmod
    
    implicit none

    integer                                   :: ndimf
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: nao,npair,i,j
    integer                                   :: lwork,ierr
    real(dp), allocatable                     :: trdens_real(:,:)
    real(dp), allocatable                     :: trdens_imag(:,:)
    real(dp), allocatable                     :: sigma(:)
    real(dp), allocatable                     :: rwork(:)
    real(dp), allocatable                     :: val(:)
    real(dp), allocatable                     :: occ(:)
    real(dp), parameter                       :: thrsh=0.01d0
    complex(dp), allocatable                  :: trdens(:,:)
    complex(dp), dimension(ndimf)             :: wf
    complex(dp), allocatable                  :: wf1(:)
    complex(dp), allocatable                  :: tmp(:,:)
    complex(dp), allocatable                  :: VT(:,:)
    complex(dp), allocatable                  :: U(:,:)
    complex(dp), allocatable                  :: work(:)
    complex(dp), allocatable                  :: hole(:,:)
    complex(dp), allocatable                  :: particle(:,:)
    complex(dp), allocatable                  :: orb(:,:)
    character(len=*)                          :: stem
    character(len=70)                         :: filename
    type(gam_structure)                       :: gam
    
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
    allocate(trdens_real(nbas,nbas))
    trdens_real=0.0d0

    allocate(trdens_imag(nbas,nbas))
    trdens_imag=0.0d0

    allocate(trdens(nbas,nbas))
    trdens=czero
    
    allocate(wf1(ndimf))
    wf1=czero

    allocate(tmp(nocc,nvirt))
    tmp=czero

    allocate(sigma(nocc))
    sigma=0.0d0

    allocate(U(nocc,nocc))
    U=czero
    
    allocate(VT(nvirt,nvirt))
    VT=czero

    lwork=4*nocc+2*nvirt
    allocate(work(lwork))
    work=czero

    allocate(rwork(5*nocc))
    rwork=0.0d0

    nao=gam%nbasis
    allocate(hole(nao,nocc))
    allocate(particle(nao,nvirt))
    hole=czero
    particle=czero

    allocate(orb(nao,nbas))
    orb=czero

    allocate(val(nbas))
    val=0.0d0

    allocate(occ(nbas))
    occ=0.0d0
    
!----------------------------------------------------------------------
! Normalise the wavefunction vector
!----------------------------------------------------------------------
    wf1=wf/sqrt(dot_product(wf,wf))
    
!----------------------------------------------------------------------
! Calculate the real part of the occupied-virtural block of the
! ADC(2) ground-to-excited state transition density matrix
!----------------------------------------------------------------------
    call adc2_trdens_gs_virt_occ(trdens_real,ndimf,kpqf,real(wf1),1,&
         rhomp2)

!----------------------------------------------------------------------
! Calculate the imaginary part of the occupied-virtural block of the
! ADC(2) ground-to-excited state transition density matrix
!----------------------------------------------------------------------
    call adc2_trdens_gs_virt_occ(trdens_imag,ndimf,kpqf,aimag(wf1),1,&
         rhomp2)
    
!----------------------------------------------------------------------
! Total, complex ground state-to-excited state transition density
! matrix
!----------------------------------------------------------------------
    trdens=trdens_real+ci*trdens_imag
    
!----------------------------------------------------------------------
! SVD of the occupied-virtual block of the ground state-to-excited
! state transition density matrix
!----------------------------------------------------------------------
    tmp=transpose(trdens(nocc+1:nbas,1:nocc))
    
    call zgesvd('A','A',nocc,nvirt,tmp,nocc,sigma,U,nocc,VT,nvirt,&
         work,lwork,rwork,ierr)

    if (ierr.ne.0) then
       errmsg='SVD of the transition density matrix failed in &
            subroutine tdadc2_nto'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate and output the complex time-dependent NTOs
!----------------------------------------------------------------------
    ! Hole NTOs in terms of the AOs
    hole=transpose(matmul(transpose(U(:,:)),&
         transpose(ao2mo(1:nao,1:nocc))))

    ! Particle NTOs in terms of the AOs
    particle=transpose(matmul(VT(:,:),&
         transpose(ao2mo(1:nao,nocc+1:nbas))))

    ! No. NTO pairs
    npair=0
    do j=1,nocc
       if (sigma(j)**2.ge.thrsh) npair=npair+1
    enddo
    
    ! Dominant NTOs and singular values
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

    ! Write the molden files: one for the real parts of the NTOs and
    ! one for the imaginary parts
    filename=trim(stem)//'real.molden'
    call write_molden(gam,filename,nao,2*npair,&
         real(orb(1:nao,1:2*npair)),val(1:2*npair),occ(1:2*npair))
    filename=trim(stem)//'imag.molden'
    call write_molden(gam,filename,nao,2*npair,&
         aimag(orb(1:nao,1:2*npair)),val(1:2*npair),occ(1:2*npair))
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(trdens_real)
    deallocate(trdens_imag)
    deallocate(trdens)
    deallocate(wf1)
    deallocate(tmp)
    deallocate(sigma)
    deallocate(U)
    deallocate(VT)
    deallocate(work)
    deallocate(rwork)
    deallocate(hole)
    deallocate(particle)
    deallocate(orb)
    deallocate(val)
    deallocate(occ)
    
    return
    
  end subroutine tdadc2_nto
    
!######################################################################
  
end module nto
