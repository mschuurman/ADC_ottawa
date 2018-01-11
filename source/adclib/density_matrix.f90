!######################################################################
! density_matrix: routines for the calculation of 1-electron reduced
!                 density and transition density matrices
!######################################################################

module density_matrix

  use constants
  use parameters
  use iomod
  use channels
  use timingmod
  use vpqrsmod
  
  implicit none

contains

!######################################################################
! adc2_trden_gs: calculation of the ADC(2) ground-to-excited state
!                transition density matrix
!######################################################################
    
  subroutine adc2_trden_gs(trdens,ndimf,kpqf,rvec,nstates)

    use mp2
      
    implicit none

    integer                                   :: i
    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    real(d), dimension(nbas,nbas,nstates)     :: trdens
    real(d), dimension(ndimf,nstates)         :: rvec
    real(d), allocatable                      :: rhomp2(:,:)
    real(d)                                  :: tw1,tw2,tc1,tc2
      
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(2x,a,/)') 'Calculating the ground-to-excited &
         state transition density matrices...'

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
      
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    trdens=0.0d0

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rhomp2(nbas,nbas))
    rhomp2=0.0d0
      
!----------------------------------------------------------------------
! Precalculation of the MP2 correction to the ground state density
! matrix
!----------------------------------------------------------------------
    ! MP2 ground state density matrix
    call rho_mp2(rhomp2)

    ! Subtraction of the zeroth-order contribution to obtain the
    ! MP2 correction
    do i=1,nocc
       rhomp2(i,i)=rhomp2(i,i)-2.0d0
    enddo

!----------------------------------------------------------------------
! Calculation of the occupied-occupied block
!----------------------------------------------------------------------
    call adc2_trdens_gs_occ_occ(trdens,ndimf,kpqf,rvec,nstates,&
         rhomp2)

!----------------------------------------------------------------------
! Calculation of the virtual-virtual block
!----------------------------------------------------------------------
    call adc2_trdens_gs_virt_virt(trdens,ndimf,kpqf,rvec,nstates,&
         rhomp2)

!----------------------------------------------------------------------
! Calculation of the virtual-occupied block
!----------------------------------------------------------------------
    call adc2_trdens_gs_virt_occ(trdens,ndimf,kpqf,rvec,nstates,&
         rhomp2)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rhomp2)

!----------------------------------------------------------------------
! Finish timing and output the time taken
!----------------------------------------------------------------------
    call times(tw2,tc2)

    write(ilog,'(2x,a,2x,F7.2,1x,a1)') 'Time taken:',tw2-tw1,'s'
      
    return
      
  end subroutine adc2_trden_gs

!######################################################################
! adc2_trden_gs_occ_occ: calculation of the occupied-occupied block of
!                        the ADC(2) ground-to-excited state transition
!                        density matrix
!######################################################################
    
  subroutine adc2_trdens_gs_occ_occ(trdens,ndimf,kpqf,rvec,nstates,&
       rhomp2)
      
    implicit none

    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: cnt,nlim1,nlim2
    integer                                   :: i,j,k,a,b
    real(d), dimension(nbas,nbas,nstates)     :: trdens
    real(d), dimension(ndimf,nstates)         :: rvec
    real(d), dimension(nbas,nbas)             :: rhomp2
    real(d), allocatable                      :: tmp(:,:,:)
    real(d)                                   :: delta_ijaa,&
                                                 delta_ijab,&
                                                 delta_ikaa,&
                                                 delta_ikab
      
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(tmp(nbas,nbas,nstates))
      
!----------------------------------------------------------------------
! 1h1p contributions
!----------------------------------------------------------------------
    nlim1=1
    nlim2=kpqf(1,0)

    tmp=0.0d0
      
    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       a=kpqf(5,cnt)
       do i=1,nocc
          tmp(i,j,:)=tmp(i,j,:)+rhomp2(i,a)*rvec(cnt,:)
       enddo
    enddo
    
    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         +2.0d0*sqrt(2.0d0)*tmp(1:nocc,1:nocc,:)
      
!----------------------------------------------------------------------
! 2h2p a=b, i=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(2,0)
    
    tmp=0.0d0
      
    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       a=kpqf(5,cnt)
       do i=1,nocc
          delta_ijaa=1.0d0/(e(i)+e(j)-e(a)-e(a))
          tmp(i,j,:)=tmp(i,j,:)+delta_ijaa*vpqrs(a,j,a,i)*rvec(cnt,:)
       enddo
    enddo
    
    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         -2.0d0*tmp(1:nocc,1:nocc,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(3,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       a=kpqf(5,cnt)
       b=kpqf(6,cnt)
       do i=1,nocc
          delta_ijab=1.0d0/(e(i)+e(j)-e(a)-e(b))
          tmp(i,j,:)=tmp(i,j,:) &
               +delta_ijab*(vpqrs(a,i,b,j)+vpqrs(a,j,b,i))*rvec(cnt,:)
       enddo
    enddo

    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         -sqrt(2.0d0)*tmp(1:nocc,1:nocc,:)

!----------------------------------------------------------------------
! 2h2p a=b, i|=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(4,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       k=kpqf(4,cnt)
       a=kpqf(5,cnt)
       do i=1,nocc
          delta_ikaa=1.0d0/(e(i)+e(k)-e(a)-e(a))
          tmp(i,j,:)=tmp(i,j,:)+delta_ikaa*vpqrs(a,i,a,k)*rvec(cnt,:)
          tmp(i,k,:)=tmp(i,k,:)+delta_ikaa*vpqrs(a,i,a,j)*rvec(cnt,:)
       enddo
    enddo

    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         +sqrt(2.0d0)*tmp(1:nocc,1:nocc,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i|=j I contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(5,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       k=kpqf(4,cnt)
       a=kpqf(5,cnt)
       b=kpqf(6,cnt)
       do i=1,nocc
          delta_ikab=1.0d0/(e(i)+e(k)-e(a)-e(b))
          tmp(i,j,:)=tmp(i,j,:)&
               +delta_ikab*(vpqrs(a,i,b,k)-vpqrs(a,k,b,i))*rvec(cnt,:)
          tmp(i,k,:)=tmp(i,j,:)&
               -delta_ikab*(vpqrs(a,i,b,j)-vpqrs(a,j,b,i))*rvec(cnt,:)
       enddo
    enddo
    
    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         +sqrt(3.0d0)*tmp(1:nocc,1:nocc,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i|=j II contributions
!----------------------------------------------------------------------
    tmp=0.0d0

    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       k=kpqf(4,cnt)
       a=kpqf(5,cnt)
       b=kpqf(6,cnt)
       do i=1,nocc
          delta_ikab=1.0d0/(e(i)+e(k)-e(a)-e(b))
          tmp(i,j,:)=tmp(i,j,:)&
               +delta_ikab*(vpqrs(a,i,b,k)+vpqrs(a,k,b,i))*rvec(cnt,:)
          tmp(i,k,:)=tmp(i,j,:)&
               +delta_ikab*(vpqrs(a,i,b,j)+vpqrs(a,j,b,i))*rvec(cnt,:)
       enddo
    enddo
    
    trdens(1:nocc,1:nocc,:)=trdens(1:nocc,1:nocc,:) &
         +tmp(1:nocc,1:nocc,:)
      
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(tmp)
      
    return
    
  end subroutine adc2_trdens_gs_occ_occ
      
!######################################################################
! adc2_trden_gs_virt_virt: calculation of the virtual-virtual block
!                          of the ADC(2) ground-to-excited state
!                          transition density matrix
!######################################################################
  
  subroutine adc2_trdens_gs_virt_virt(trdens,ndimf,kpqf,rvec,&
       nstates,rhomp2)
    
    implicit none
    
    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: cnt,nlim1,nlim2
    integer                                   :: i,j,a,b,c
    real(d), dimension(nbas,nbas,nstates)     :: trdens
    real(d), dimension(ndimf,nstates)         :: rvec
    real(d), dimension(nbas,nbas)             :: rhomp2
    real(d), allocatable                      :: tmp(:,:,:)
    real(d)                                   :: delta_iiab,&
                                                 delta_iicb,&
                                                 delta_ijab,&
                                                 delta_ijcb

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
      allocate(tmp(nbas,nbas,nstates))

!----------------------------------------------------------------------
! 1h1p contributions
!----------------------------------------------------------------------
    nlim1=1
    nlim2=kpqf(1,0)

    tmp=0.0d0
      
    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       a=kpqf(5,cnt)
       do b=nocc+1,nbas
          tmp(a,b,:)=tmp(a,b,:)+rhomp2(i,b)*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         -2.0d0*sqrt(2.0d0)*tmp(nocc+1:nbas,nocc+1:nbas,:)

!----------------------------------------------------------------------
! 2h2p a=b, i=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(2,0)
    
    tmp=0.0d0
      
    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       a=kpqf(5,cnt)
       do b=nocc+1,nbas
          delta_iiab=1.0d0/(e(i)+e(i)-e(a)-e(b))
          tmp(a,b,:)=tmp(a,b,:)+delta_iiab*vpqrs(b,i,a,i)*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         +2.0d0*tmp(nocc+1:nbas,nocc+1:nbas,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(3,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       a=kpqf(5,cnt)
       c=kpqf(6,cnt)
       do b=nocc+1,nbas
          delta_iicb=1.0d0/(e(i)+e(i)-e(c)-e(b))
          tmp(a,b,:)=tmp(a,b,:)+delta_iicb*vpqrs(b,i,c,i)*rvec(cnt,:)
          tmp(c,b,:)=tmp(c,b,:)+delta_iicb*vpqrs(b,i,a,i)*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         +sqrt(2.0d0)*tmp(nocc+1:nbas,nocc+1:nbas,:)

!----------------------------------------------------------------------
! 2h2p a=b, i|=j contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(4,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       j=kpqf(4,cnt)
       a=kpqf(5,cnt)
       do b=nocc+1,nbas
          delta_ijab=1.0d0/(e(i)+e(j)-e(a)-e(b))
          tmp(a,b,:)=tmp(a,b,:) &
               +delta_ijab*(vpqrs(b,i,a,j)-vpqrs(b,j,a,i))*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         -sqrt(2.0d0)*tmp(nocc+1:nbas,nocc+1:nbas,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i|=j I contributions
!----------------------------------------------------------------------
    nlim1=nlim2+1
    nlim2=nlim2+kpqf(5,0)

    tmp=0.0d0

    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       j=kpqf(4,cnt)
       a=kpqf(5,cnt)
       c=kpqf(6,cnt)
       do b=nocc+1,nbas
          delta_ijcb=1.0d0/(e(i)+e(j)-e(c)-e(b))
          tmp(a,b,:)=tmp(a,b,:) &
               +delta_ijcb*(vpqrs(b,i,c,j)-vpqrs(b,j,c,i))*rvec(cnt,:)
          tmp(c,b,:)=tmp(c,b,:) &
               -delta_ijcb*(vpqrs(b,i,a,j)-vpqrs(b,j,a,i))*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         -sqrt(3.0d0)*tmp(nocc+1:nbas,nocc+1:nbas,:)

!----------------------------------------------------------------------
! 2h2p a|=b, i|=j II contributions
!----------------------------------------------------------------------
    tmp=0.0d0

    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       j=kpqf(4,cnt)
       a=kpqf(5,cnt)
       c=kpqf(6,cnt)
       do b=nocc+1,nbas
          delta_ijcb=1.0d0/(e(i)+e(j)-e(c)-e(b))
          tmp(a,b,:)=tmp(a,b,:) &
               +delta_ijcb*(vpqrs(b,i,c,j)+vpqrs(b,j,c,i))*rvec(cnt,:)
          tmp(c,b,:)=tmp(c,b,:) &
               +delta_ijcb*(vpqrs(b,i,a,j)+vpqrs(b,j,a,i))*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,nocc+1:nbas,:)=&
         trdens(nocc+1:nbas,nocc+1:nbas,:)&
         -tmp(nocc+1:nbas,nocc+1:nbas,:)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(tmp)
    
    return
      
  end subroutine adc2_trdens_gs_virt_virt

!######################################################################
! adc2_trden_gs_virt_occ: calculation of the virtual-occupied block of
!                          the ADC(2) ground-to-excited state
!                          transition density matrix
!######################################################################
  
  subroutine adc2_trdens_gs_virt_occ(trdens,ndimf,kpqf,rvec,nstates,&
       rhomp2)
      
    implicit none

    integer                                   :: ndimf,nstates
    integer, dimension(7,0:nbas**2*4*nocc**2) :: kpqf
    integer                                   :: cnt,nlim1,nlim2
    integer                                   :: i,j,k,a,b,c
    real(d), dimension(nbas,nbas,nstates)     :: trdens
    real(d), dimension(ndimf,nstates)         :: rvec
    real(d), dimension(nbas,nbas)             :: rhomp2
    real(d), allocatable                      :: tmp(:,:,:)
    real(d), allocatable                      :: lambda(:,:,:)
    real(d)                                   :: delta_abij,&
                                                 delta_ijac,&
                                                 delta_jkbc

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(tmp(nbas,nbas,nstates))
    allocate(lambda(nvirt,nocc,nstates))
    
!----------------------------------------------------------------------
! Zeroth-order 1h1p contribution
!----------------------------------------------------------------------
    nlim1=1
    nlim2=kpqf(1,0)

    tmp=0.0d0
    
    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       a=kpqf(5,cnt)
       tmp(a,i,:)=tmp(a,i,:)+rvec(cnt,:)
    enddo

    trdens(nocc+1:nbas,1:nocc,:)=trdens(nocc+1:nbas,1:nocc,:) &
         -sqrt(2.0d0)*tmp(nocc+1:nbas,1:nocc,:)

!----------------------------------------------------------------------
! First-order 1h1p contribution
!----------------------------------------------------------------------
    tmp=0.0d0
    
    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       b=kpqf(5,cnt)
       do a=nocc+1,nbas
          do i=1,nocc
             delta_abij=1.0d0/(e(a)+e(b)-e(i)-e(j))
             tmp(a,i,:)=tmp(a,i,:) &
                  +delta_abij*(2.0d0*vpqrs(b,j,a,i)-vpqrs(b,i,a,j))&
                  *rvec(cnt,:)
          enddo
       enddo       
    enddo

    trdens(nocc+1:nbas,1:nocc,:)=trdens(nocc+1:nbas,1:nocc,:) &
         +sqrt(2.0d0)*tmp(nocc+1:nbas,1:nocc,:)
    
!----------------------------------------------------------------------
! Second-order 1h1p contribution no. 1
!----------------------------------------------------------------------
    tmp=0.0d0
    
    do cnt=nlim1,nlim2
       i=kpqf(3,cnt)
       b=kpqf(5,cnt)
       do a=nocc+1,nbas
          tmp(a,i,:)=tmp(a,i,:)+rhomp2(a,b)*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,1:nocc,:)=trdens(nocc+1:nbas,1:nocc,:) &
         +sqrt(2.0d0)*tmp(nocc+1:nbas,1:nocc,:)

!----------------------------------------------------------------------
! Second-order 1h1p contribution no. 2
!----------------------------------------------------------------------
    tmp=0.0d0
    
    do cnt=nlim1,nlim2
       j=kpqf(3,cnt)
       a=kpqf(5,cnt)
       do i=1,nocc
          tmp(a,i,:)=tmp(a,i,:)+rhomp2(i,j)*rvec(cnt,:)
       enddo
    enddo

    trdens(nocc+1:nbas,1:nocc,:)=trdens(nocc+1:nbas,1:nocc,:) &
         -sqrt(2.0d0)*tmp(nocc+1:nbas,1:nocc,:)

!----------------------------------------------------------------------
! Second-order 1h1p contribution no. 3
!----------------------------------------------------------------------
    tmp=0.0d0
    lambda=0.0d0
    
    ! Calculation of intermediates
    do cnt=nlim1,nlim2
       k=kpqf(3,cnt)
       b=kpqf(5,cnt)
       do c=nocc+1,nbas
          do j=1,nocc
             delta_jkbc=1.0d0/(e(j)+e(k)-e(b)-e(c))
             lambda(c-nocc,j,:)=lambda(c-nocc,j,:) &
                  +delta_jkbc*(2.0d0*vpqrs(b,k,c,j)-vpqrs(b,j,c,k))&
                  *rvec(cnt,:)
          enddo
       enddo
    enddo

    ! Calculation of the contribution to trdens
    do a=nocc+1,nbas
       do i=1,nbas
          do c=nocc+1,nbas
             do j=1,nocc
                delta_ijac=1.0d0/(e(i)+e(j)-e(a)-e(c))
                tmp(a,i,:)=tmp(a,i,:) &
                     +delta_ijac*(2.0d0*vpqrs(i,a,j,c)-vpqrs(i,c,j,a))&
                     *lambda(c-nocc,j,:)
             enddo
          enddo
       enddo
    enddo

    trdens(nocc+1:nbas,1:nocc,:)=trdens(nocc+1:nbas,1:nocc,:) &
         -0.5d0*sqrt(2.0d0)*tmp(nocc+1:nbas,1:nocc,:)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(tmp)
    deallocate(lambda)
    
    return

  end subroutine adc2_trdens_gs_virt_occ
  
!######################################################################
  
end module density_matrix
