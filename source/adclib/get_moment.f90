module get_moment
     
  use constants
  use parameters
  use dipole_ph
  use get_matrix_DIPOLE
  use channels
  use timingmod
  
  implicit none
  
contains
  
!#######################################################################

  real(d) function tm(ndim,evec,mtm)

    implicit none
       
    integer, intent(in)                  :: ndim
    real(d), dimension(ndim), intent(in) :: evec,mtm       
    real(d)                              :: ddot

    external ddot
       
    tm=ddot(ndim,evec,1,mtm,1)
       
    return
    
  end function tm
  
!#######################################################################

  subroutine dmatvec(ndim,autvec,arr,travec)

    implicit none
    
    integer                                   :: i
    integer, intent(in)                       :: ndim
    real(d), dimension(ndim)                  :: Y
    real(d), dimension(ndim), intent(in)      :: autvec
    real(d), dimension(ndim,ndim), intent(in) :: arr
    real(d), dimension(ndim), intent(out)     :: travec
       
    external dsymv
       
    call dsymv("L",ndim,1,arr,ndim,autvec,1,0,Y,1)
    
    travec(:) = Y(:)
       
    return

  end subroutine dmatvec

!#######################################################################
  
  subroutine get_modifiedtm_tda(ndim,kpq,mtm)

    implicit none
       
    integer, intent(in)                                   :: ndim
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(in) :: kpq
    integer                                               :: i,ap,ih
    real(d), dimension(ndim), intent(out)                 :: mtm 
    real(d)                                               :: tw1,tw2,&
                                                             tc1,tc2
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Dimensionality check
!----------------------------------------------------------------------
    if (ndim.ne.kpq(1,0)) then
       errmsg='Inconsistent dimension of the TDA matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate the ADC(1) f-vector
!----------------------------------------------------------------------
    do i=1,ndim
       ih=kpq(3,i)
       ap=kpq(5,i)
       mtm(i)=dpl(ih,ap)+F0_ph(ap,ih)
    enddo

    mtm(:)=-sqrt(2.0d0)*mtm(:)

!----------------------------------------------------------------------
! Stop timing and output the time taken
!----------------------------------------------------------------------
    call times(tw2,tc2)

    write(ilog,'(2x,a,x,F9.2,/)') "Time taken:",tw2-tw1
    
    return
    
  end subroutine get_modifiedtm_tda

!#######################################################################

  subroutine get_tm_cis(ndim,kpq,mtm)

    implicit none
       
    integer, intent(in)                                   :: ndim
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(in) :: kpq
    integer                                               :: i,ap,ih
    real(d), dimension(ndim), intent(out)                 :: mtm 
    real(d)                                               :: tw1,tw2,&
                                                             tc1,tc2

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Dimensionality check
!----------------------------------------------------------------------
    if (ndim.ne.kpq(1,0)) then
       errmsg='Inconsistent dimension of the CIS matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate the CIS f-vector
!----------------------------------------------------------------------
    do i=1,ndim
       ih=kpq(3,i)
       ap=kpq(5,i)
       mtm(i)=dpl(ih,ap)
    enddo
    
    mtm(:)=-sqrt(2.0d0)*mtm(:)    
    
!----------------------------------------------------------------------
! Stop timing and output the time taken
!----------------------------------------------------------------------
    call times(tw2,tc2)

    write(ilog,'(2x,a,x,F9.2,/)') "Time taken:",tw2-tw1
    
    return
    
  end subroutine get_tm_cis
  
!#######################################################################

  subroutine get_modifiedtm_adc2(ndim,kpq,mtm,ista)

    implicit none
    
    integer, intent(in)                                   :: ndim,ista
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(in) :: kpq
    real*8, dimension(ndim), intent(out)                  :: mtm
    integer                                               :: a,b,k,j,cnt
    integer                                               :: nlim1,nlim2
    
    integer                              :: a1,c1,c,l1,l,i,&
         itmp,itmp1,nvirt
    real*8                               :: tw1,tw2,tc1,tc2,ftmp
    real*8, dimension(:,:), allocatable  :: tau
    
    write(ilog,'(/,2x,a)') 'Calculating the f-vector &
         f_J=<Psi_0|D|Psi_J>...'
    call times(tw1,tc1)
    
    mtm(:)=0.0d0

!-----------------------------------------------------------------------    
! 1h1p block
!-----------------------------------------------------------------------
    nlim1=1
    nlim2=kpq(1,0)
    nvirt=nbas-nocc
    
    ! F25_ph, F26_ph, F27_ph, F28_ph, F29_ph, F210_ph (note that
    ! the calculation of these terms is still very slow)
    !$omp parallel do private(cnt,k,a) shared(dpl,mtm,kpq,nlim1,nlim2)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=dpl(a,k)
       mtm(cnt)=mtm(cnt)+F0_ph(a,k)
       mtm(cnt)=mtm(cnt)+F25_ph(a,k)
       mtm(cnt)=mtm(cnt)+F26_ph(a,k)
       mtm(cnt)=mtm(cnt)+F27_ph(a,k)
       mtm(cnt)=mtm(cnt)+F28_ph(a,k)
       mtm(cnt)=mtm(cnt)+F29_ph(a,k)
       mtm(cnt)=mtm(cnt)+F210_ph(a,k)
    enddo
    !$omp end parallel do

    ! FA_ph
    allocate(tau(nvirt,nvirt))
    !$omp parallel do private(a,c,itmp,itmp1) shared(tau)
    do a=nocc+1,nbas
       itmp=a-nocc
       do c=nocc+1,nbas
          itmp1=c-nocc
          tau(itmp,itmp1)=tauA(c,a)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,a,k,itmp,itmp1,ftmp,c1,c) shared(kpq,roccnum,tau,dpl,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       itmp=a-nocc
       itmp1=0
       ftmp=0.0d0
       do c1=nocc+1,nbas
          c=roccnum(c1)
          itmp1=itmp1+1
          ftmp=ftmp+tau(itmp,itmp1)*dpl(c,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)

    ! FB_ph
    allocate(tau(nocc,nocc))
    !$omp parallel do private(l,k) shared(tau)
    do l=1,nocc
       do k=1,nocc
          tau(l,k)=tauB(l,k)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a,ftmp,l1,l) shared(kpq,roccnum,dpl,tau,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       ftmp=0.0d0
       do l1=1,nocc
          l=roccnum(l1)
          ftmp=ftmp+dpl(a,l)*tau(l1,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)
    
    ! FC_ph
    allocate(tau(nvirt,nocc))
    !$omp parallel do private(a,itmp,l) shared(tau)
    do a=nocc+1,nbas
       itmp=a-nocc
       do l=1,nocc
          tau(itmp,l)=tauC(a,l)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a) shared(kpq,mtm,tau)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=mtm(cnt)+FC_ph_new(a,k,tau,nvirt)
    enddo
    !$omp end parallel do
    deallocate(tau)
    
    ! F21_ph
    allocate(tau(nvirt,nocc))
    !$omp parallel do private(a,itmp,l) shared(tau)
    do a=nocc+1,nbas
       itmp=a-nocc
       do l=1,nocc
          tau(itmp,l)=tau21(a,l)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a,itmp,ftmp,l1,l) shared(kpq,roccnum,tau,dpl,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       itmp=a-nocc
       ftmp=0.0d0
       do l1=1,nocc
          l=roccnum(l1)
          ftmp=ftmp+tau(itmp,l1)*dpl(l,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)

    ! F22_ph
    allocate(tau(nvirt,nocc))
    !$omp parallel do private(a,itmp,l) shared(tau)
    do a=nocc+1,nbas
       itmp=a-nocc
       do l=1,nocc
          tau(itmp,l)=tau22(a,l)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a,itmp,ftmp,l1,l) shared(kpq,roccnum,tau,dpl,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       itmp=a-nocc
       ftmp=0.0d0
       do l1=1,nocc
          l=roccnum(l1)
          ftmp=ftmp+tau(itmp,l1)*dpl(l,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)

    ! F23_ph
    allocate(tau(nvirt,nocc))
    !$omp parallel do private(c,itmp,k) shared(tau)
    do c=nocc+1,nbas
       itmp=c-nocc
       do k=1,nocc
          tau(itmp,k)=tau23(c,k)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a,itmp,ftmp,c1,c) shared(kpq,roccnum,tau,dpl,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       itmp=0
       ftmp=0.0d0
       do c1=nocc+1,nbas
          itmp=itmp+1
          c=roccnum(c1)
          ftmp=ftmp+dpl(a,c)*tau(itmp,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)

    ! F24_ph
    allocate(tau(nvirt,nocc))
    !$omp parallel do private(c,itmp,k) shared(tau)
    do c=nocc+1,nbas
       itmp=c-nocc
       do k=1,nocc
          tau(itmp,k)=tau24(c,k)
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(cnt,k,a,itmp,ftmp,c1,c) shared(kpq,roccnum,tau,dpl,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       itmp=0
       ftmp=0.0d0
       do c1=nocc+1,nbas
          itmp=itmp+1
          c=roccnum(c1)
          ftmp=ftmp+dpl(a,c)*tau(itmp,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    enddo
    !$omp end parallel do
    deallocate(tau)

    mtm(:)=-sqrt(2.0d0)*mtm(:)

!-----------------------------------------------------------------------    
! 2h2p block
!-----------------------------------------------------------------------
    if (lcvs) goto 100
    if (lcvsfinal.and.ista.eq.0) goto 100

    ! I-a=b,i=j
    nlim1=nlim2+1
    nlim2=nlim2+kpq(2,0)
    !$omp parallel do private(cnt,k,a) shared(kpq,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=FI_2p2h(a,k)
    enddo
    !$omp end parallel do
    
    ! II-a|=b,i=j
    nlim1=nlim2+1
    nlim2=nlim2+kpq(3,0)
    !$omp parallel do private(cnt,k,a,b) shared(kpq,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)
       mtm(cnt)=FII_2p2h(a,b,k)
    enddo
    
100 continue
    
    ! III-a=b,i|=j
    nlim1=nlim2+1
    nlim2=nlim2+kpq(4,0)
    !$omp parallel do private(cnt,k,j,a) shared(kpq,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=FIII_2p2h(a,k,j)
    enddo
    !$omp end parallel do
    
    ! IV1-a|=b,i|=j
    nlim1=nlim2+1
    nlim2=nlim2+kpq(5,0)
    !$omp parallel do private(k,j,a,b) shared(kpq,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)
       mtm(cnt)=FIV1_2p2h(a,b,k,j)
    enddo
    !$omp end parallel do

    ! IV2-a|=b,i|=j
    !$omp parallel do private(cnt,k,j,a,b) shared(kpq,mtm)
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)
       mtm(cnt+kpq(5,0))=FIV2_2p2h(a,b,k,j)
    enddo
    !$omp end parallel do
    
    call times(tw2,tc2)
    write(ilog,'(2x,a,x,F9.2,/)') "Time taken:",tw2-tw1
    
    return
       
  end subroutine get_modifiedtm_adc2

!#######################################################################

end module get_moment
    
    
