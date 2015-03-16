module get_moment
  
  use constants
  use parameters
  use dipole_ph
  use get_matrix_DIPOLE

  implicit none
  
contains
  
  real(d) function tm(ndim,evec,mtm)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: evec,mtm
    
    real(d) :: ddot
    external ddot
    
    tm=ddot(ndim,evec,1,mtm,1)
!!$    tm=tm**2
    
  end function tm




!!!!! NEW SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine  dmatvec(ndim,autvec,arr,travec)
    
    integer :: i
    integer, intent(in) :: ndim
    real(d), dimension(ndim) :: Y
    real(d), dimension(ndim), intent(in) :: autvec
    real(d), dimension(ndim,ndim), intent(in) :: arr
    real(d), dimension(ndim), intent(out) :: travec
    
    external dsymv
    
    call dsymv("L",ndim,1,arr,ndim,autvec,1,0,Y,1)
    
     travec(:) = Y(:)
    
  end subroutine dmatvec

!!!!!!!!!!! END NEW SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$--------------------------------------------
  
  subroutine get_modifiedtm_tda(ndim,kpq,mtm)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim), intent(out) :: mtm 

    integer :: i,ap,ih

    if (ndim .ne. kpq(1,0)) then
       write(6,*) 'Inconsistent dim of the TDA matrix'
       stop
    end if

    do i= 1,ndim
       ih=kpq(3,i)
       ap=kpq(5,i)
       mtm(i)=dpl(ih,ap)+F0_ph(ap,ih)
    end do
    
    mtm(:)=-sqrt(2._d)*mtm(:)

  end subroutine get_modifiedtm_tda

!!$----------------------------------------------

  subroutine get_modifiedtm_adc2(ndim,kpq,mtm)
    
    integer, intent(in)                                   :: ndim
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(in) :: kpq
    real*8, dimension(ndim), intent(out)                  :: mtm
    integer                                               :: a,b,k,j,cnt
    integer                                               :: nlim1,nlim2

    integer                              :: a1,c1,c,l1,l,i,&
                                            itmp,itmp1,nvirt
    real*8                               :: t1,t2,ftmp
    real*8, dimension(:,:), allocatable  :: tau

    mtm(:)=0.0d0


    
!!$-----1h1p block------
    nlim1=1
    nlim2=kpq(1,0)

    nvirt=nbas-nocc

!-----------------------------------------------------------------------
! Not yet improved: F25_ph, F26_ph, F27_ph, F28_ph, F29_ph, F210_ph
!-----------------------------------------------------------------------
    do cnt= nlim1,nlim2
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

!-----------------------------------------------------------------------
! FA_ph
!-----------------------------------------------------------------------
    allocate(tau(nvirt,nvirt))
    itmp=0
    do a=nocc+1,nbas
       itmp=itmp+1
       itmp1=0
       do c=nocc+1,nbas
          itmp1=itmp1+1
          tau(itmp,itmp1)=tauA(c,a)
       enddo
    enddo

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
    end do

    deallocate(tau)

!-----------------------------------------------------------------------
! FB_ph
!-----------------------------------------------------------------------
    allocate(tau(nocc,nocc))
    do l=1,nocc
       do k=1,nocc
          tau(l,k)=tauB(l,k)
       enddo
    enddo

    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       ftmp=0.0d0
       do l1=1,nocc
          l=roccnum(l1)
          ftmp=ftmp+dpl(a,l)*tau(l1,k)
       enddo
       mtm(cnt)=mtm(cnt)+ftmp
    end do

    deallocate(tau)

!-----------------------------------------------------------------------
! FC_ph
!-----------------------------------------------------------------------
    allocate(tau(nvirt,nocc))
    itmp=0
    do a=nocc+1,nbas
       itmp=itmp+1
       do l=1,nocc
          tau(itmp,l)=tauC(a,l)
       enddo
    enddo
    do cnt=nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=mtm(cnt)+FC_ph_new(a,k,tau,nvirt)
    enddo
    deallocate(tau)

!-----------------------------------------------------------------------
! F21_ph
!-----------------------------------------------------------------------
    allocate(tau(nvirt,nocc))
    itmp=0
    do a=nocc+1,nbas
       itmp=itmp+1
       do l=1,nocc
          tau(itmp,l)=tau21(a,l)
       enddo
    enddo

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
    end do

    deallocate(tau)

!-----------------------------------------------------------------------
! F22_ph
!----------------------------------------------------------------------- 
    allocate(tau(nvirt,nocc))
    itmp=0
    do a=nocc+1,nbas
       itmp=itmp+1
       do l=1,nocc
          tau(itmp,l)=tau22(a,l)
       enddo
    enddo

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
    end do

    deallocate(tau)

!----------------------------------------------------------------------- 
! F23_ph
!----------------------------------------------------------------------- 
    allocate(tau(nvirt,nocc))
    itmp=0
    do c=nocc+1,nbas
       itmp=itmp+1
       do k=1,nocc
          tau(itmp,k)=tau23(c,k)
       enddo
    enddo

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
    end do

    deallocate(tau)

!----------------------------------------------------------------------- 
! F24_ph
!----------------------------------------------------------------------- 
    allocate(tau(nvirt,nocc))
    itmp=0
    do c=nocc+1,nbas
       itmp=itmp+1
       do k=1,nocc
          tau(itmp,k)=tau24(c,k)
       enddo
    enddo
    
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
    end do

    deallocate(tau)

    mtm(:)=-sqrt(2.0d0)*mtm(:)

!!$----------I-a=b,i=j-------------------

    nlim1=nlim2+1
    nlim2=nlim2+kpq(2,0)
    
    do cnt= nlim1,nlim2

       k=kpq(3,cnt)
       a=kpq(5,cnt)

       mtm(cnt)=FI_2p2h(a,k)

    end do
 
!!$----------II-a|=b,i=j-------------------

    nlim1=nlim2+1
    nlim2=nlim2+kpq(3,0)
    
    do cnt= nlim1,nlim2
       
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)

       mtm(cnt)=FII_2p2h(a,b,k)

    end do

!!$----------III-a=b,i|=j-------------------

    nlim1=nlim2+1
    nlim2=nlim2+kpq(4,0)
    
    do cnt= nlim1,nlim2
       
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)

       mtm(cnt)=FIII_2p2h(a,k,j)

    end do    

!!$----------IV1-a|=b,i|=j-------------------

    nlim1=nlim2+1
    nlim2=nlim2+kpq(5,0)
    
    do cnt= nlim1,nlim2
       
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)

       mtm(cnt)=FIV1_2p2h(a,b,k,j)

    end do 

!!$----------IV2-a|=b,i|=j-------------------

    do cnt= nlim1,nlim2
       
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)

       mtm(cnt+kpq(5,0))=FIV2_2p2h(a,b,k,j)

    end do

  end subroutine get_modifiedtm_adc2


end module get_moment
    
    
