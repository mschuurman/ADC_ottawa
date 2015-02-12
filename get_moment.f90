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
    real(d), dimension(ndim), intent(out)                 :: mtm
    
    integer :: a,b,k,j,cnt
    integer :: nlim1,nlim2

    real(d) :: t1,t2

    mtm(:)=0._d
    
!!$-----1h1p block------
   
    nlim1=1
    nlim2=kpq(1,0)
    
    do cnt= nlim1,nlim2
       k=kpq(3,cnt)
       a=kpq(5,cnt)
       mtm(cnt)=dpl(a,k)+F0_ph(a,k)+FA_ph(a,k)+FB_ph(a,k)+FC_ph(a,k)       
       mtm(cnt)=mtm(cnt)+F21_ph(a,k)+F22_ph(a,k)+F23_ph(a,k)+F24_ph(a,k)+F25_ph(a,k)
       mtm(cnt)=mtm(cnt)+F26_ph(a,k)+F27_ph(a,k)+F28_ph(a,k)+F29_ph(a,k)+F210_ph(a,k)
    end do

    mtm(:)=-sqrt(2._d)*mtm(:)

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

    call cpu_time(t1)

    do cnt= nlim1,nlim2
       
       k=kpq(3,cnt)
       j=kpq(4,cnt)
       a=kpq(5,cnt)
       b=kpq(6,cnt)

       mtm(cnt+kpq(5,0))=FIV2_2p2h(a,b,k,j)

    end do

  end subroutine get_modifiedtm_adc2


end module get_moment
    
    
