module get_matrix_DIPOLE

  
  use constants
  use parameters
  use D_matrix
  use misc
  use filetools
  use dipole_ph
  
  implicit none

  integer, parameter :: buf_size=8192
  
contains


!!$*******************************************************************************
!!$*******************************************************************************
!!$*******************************************************************************
!!$*ON THE FLY SCALAR PRODUCT OF THE DIPOLE MATRIX WITH THE INITIAL STATE VECTOR*!
!!$*ON THE FLY SCALAR PRODUCT OF THE DIPOLE MATRIX WITH THE INITIAL STATE VECTOR*!
!!$*ON THE FLY SCALAR PRODUCT OF THE DIPOLE MATRIX WITH THE INITIAL STATE VECTOR*!
!!$*******************************************************************************
!!$*******************************************************************************
!!$*******************************************************************************
!!$*******************************************************************************




!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------------------------------------------------



  subroutine get_dipole_initial_product(ndim,ndimf,kpq,kpqf,autvec,travec)

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 

    integer, intent(in) :: ndim,ndimf
    real(d), dimension(ndim), intent(in) :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    
    integer :: inda,indb,indk,indl,spin
    integer :: indapr,indbpr,indkpr,indlpr,spinpr 
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,dim_count,ndim1,dim_countf,ndim1f
    real(d) :: ar_offdiag_ij
    integer :: k,k1,b,b1 
    
    write(6,*) "Writing the travec vector of ADC-DIPOLE matrix INITIAL-STATE product "



!   write(6,*) 'inside routine'
!    do k1=1,nOcc
!       k=roccnum(k1) 
!     do b1=nOcc+1,nBas
!        b=roccnum(b1)
!   write(6,*) 'density',k,b,densityhc(k,b)
!    write(6,*) 'density',k1,b1,densityhc(k,b),density(k,b),proper_density(k,b)
!     end do
!   end do






    travec(:)=0.0

! THE INDEX i RUNS IN THE 1H1P  BLOCK OF (FINAL) CONFIGURATIONS  

    ndim1f=kpqf(1,0)
    do i=1,ndim1f
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
    
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             
       
         ar_offdiag_ij = 0.d0

         ar_offdiag_ij = D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_5_2_ph_ph(indk,indkpr)
    end if


       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do


       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    
 
ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiag_ij = D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
 
ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

ar_offdiag_ij = 0.d0


  if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
 
ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end do
!!!!!! end of the first i cycle: the first travec element has been computed!!!!!

       write(6,*),i,ndim1f,travec(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do
! end of the 1h1piblock part:all the 1h1pblock element of travec have been!!!!!!
! computed!!!


!*************************************************************************
!*************************************************************************
!*************************************************************************
!THE INDEX i RUNS IN THE FIRST 2H2P BLOCK OF (FINAL) CONFIGURATIONS a=b i=j
!*************************************************************************
!*************************************************************************
!*************************************************************************

    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
   
 
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
         ar_offdiag_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiag_ij = D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    

!!$ (1,1) block
 
       ar_offdiag_ij = 0.d0

       ar_offdiag_ij = ar_offdiag_ij +  D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

!!$ (1,2) block
 
       ar_offdiag_ij = 0.d0

       ar_offdiag_ij = ar_offdiag_ij +D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (1,3) block

       ar_offdiag_ij = ar_offdiag_ij +D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
      dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (1,4i) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (1,4ii) block
    
       ar_offdiag_ij = ar_offdiag_ij +D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
! end of the first i cycle: the first travec element has been computed!!!

       write(6,*),i,travec(i)

    end do
    

!*************************************************************************
!*************************************************************************
!*************************************************************************
!THE INDEX i RUNS IN THE SECOND BLOCK OF 2H2P (FINAL) CONFIGURATIONS  a/=b i=j
!*************************************************************************
!*************************************************************************
!*************************************************************************

    dim_countf=dim_countf+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
         ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij +D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do


       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij = 0.d0

!!$ (2,1) block
    
       ar_offdiag_ij = ar_offdiag_ij +D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (2,2) block

       ar_offdiag_ij = ar_offdiag_ij +D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (2,3) block

       ar_offdiag_ij = ar_offdiag_ij +D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (2,4i) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (2,4ii) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

       write(6,*),i,travec(i)

! end of the first i cycle: the first travec element has been computed!!!

    end do
    

!*************************************************************************
!*************************************************************************
!*************************************************************************
!THE INDEX i RUNS IN THE THIRD BLOCK OF 2H2P (FINAL) CONFIGURATIONS  a=b i/=j
!*************************************************************************
!*************************************************************************
!*************************************************************************

    dim_countf=dim_countf+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

    
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
         ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij +D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do


       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij = 0.d0

!!$ (3,1) block

       ar_offdiag_ij = ar_offdiag_ij +D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (3,2) block

       ar_offdiag_ij = ar_offdiag_ij +D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (3,3) block

       ar_offdiag_ij = ar_offdiag_ij +D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (3,4i) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (3,4ii) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

       write(6,*),i,travec(i)

! end of the first i cycle: the first travec element has been computed!!!

    end do




!*************************************************************************
!*************************************************************************
!*************************************************************************
!THE INDEX i RUNS IN THE 4I BLOCK OF 2H2P (FINAL) CONFIGURATIONS  a/=b i/=j   SPIN CASE 1
!*************************************************************************
!*************************************************************************
!*************************************************************************    

    dim_countf=dim_countf+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

    
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
         ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij +D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)


       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do


       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij = 0.d0

!!$ (4i,1) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (4i,2) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (4i,3) block

       ar_offdiag_ij = ar_offdiag_ij +D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)


       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (4i,4i) block

       ar_offdiag_ij = D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (4i,4ii) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

       write(6,*),i,travec(i)

! end of the first i cycle: the first travec element has been computed!!!

    end do





!*************************************************************************
!*************************************************************************
!*************************************************************************
!THE INDEX i RUNS IN THE 4II BLOCK OF 2H2P (FINAL) CONFIGURATIOS  a/=b i/=j SPIN CASE 2
!*************************************************************************
!*************************************************************************
!*************************************************************************

    dim_countf=dim_countf+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

    
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
         ar_offdiag_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij +D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do


       dim_count=kpq(1,0)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij = 0.d0

!!$ (4ii,1) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i=j,a|=b configs   
    
       dim_count=dim_count+kpq(2,0)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (4ii,2) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
    
!!$ Coupling to the i|=j,a=b configs
    
       dim_count=dim_count+kpq(3,0)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (4ii,3) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij = 0.d0

!!$ (4ii,4i) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij = 0.d0

!!$ (4ii,4ii) block

       ar_offdiag_ij = ar_offdiag_ij +D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

       write(6,*),i,travec(i)

! end of the first i cycle: the first travec element has been computed!!!

    end do
    
  end subroutine get_dipole_initial_product
!!$-----------------------------------





!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------





  subroutine get_dipole_initial_product_tda(ndim,ndimf,kpq,kpqf,autvec,travec)

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 

    integer, intent(in) :: ndim,ndimf
    real(d), dimension(ndim), intent(in) :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    
    integer :: inda,indb,indk,indl,spin
    integer :: indapr,indbpr,indkpr,indlpr,spinpr 
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,dim_count,ndim1,dim_countf,ndim1f
    real(d) :: ar_offdiag_ij
    integer :: k,k1,b,b1 
    
    write(6,*) "Writing the travec vector of ADC-DIPOLE matrix INITIAL-STATE product "



    travec(:)=0.0

! THE INDEX i RUNS IN THE 1H1P  BLOCK OF (FINAL) CONFIGURATIONS  

    ndim1f=kpqf(1,0)
    do i=1,ndim1f
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
    
       ndim1=kpq(1,0)
       do j=1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             
       
         ar_offdiag_ij = 0.0


    if(indk .eq. indkpr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_2_ph_ph(indk,indkpr)
    end if


       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       end do

    end do
    
  end subroutine get_dipole_initial_product_tda
!!$-----------------------------------


















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! DIPOLE MATRIX IN MEMORY SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! DIPOLE MATRIX IN MEMORY SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! DIPOLE MATRIX IN MEMORY SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!! CASE OF INITIAL AND FINAL SPACES OF DIFFERENT SYMMETRIES
!!! CASE OF INITIAL AND FINAL SPACES OF DIFFERENT SYMMETRIES
!!! CASE OF INITIAL AND FINAL SPACES OF DIFFERENT SYMMETRIES


!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------


  subroutine get_offdiag_tda_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,ar_offdiagd)

  integer, intent(in) :: ndim
  integer, intent(in) :: ndimf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  real(d), dimension(ndimf,ndim), intent(out) :: ar_offdiagd
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  
  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j

  ar_offdiagd(:,:)=0._d 

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
     ndim1f=kpqf(1,0)
  
     do i= 1,ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
        do j= 1,ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(indk .eq. indkpr) then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_2_ph_ph(indk,indkpr)
    end if


   end do
  end do

  end subroutine get_offdiag_tda_DIPOLE_direct_OK



!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------------------------------------------------

 
  subroutine get_offdiag_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,ar_offdiagd)

  integer, intent(in) :: ndim
  integer, intent(in) :: ndimf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  real(d), dimension(ndimf,ndim), intent(out) :: ar_offdiagd
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  
  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j

  ar_offdiagd(:,:)=0._d 

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
     ndim1f=kpqf(1,0)
  
     do i= 1,ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

        do j= 1,ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_5_2_ph_ph(indk,indkpr) 
    end if


   end do
  end do
     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)

       do i= 1,ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j= dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)

       do i= 1,ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j= dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)

       do i= 1,ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j= dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i= 1,ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i= 1,ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!2H2P-KPQF      1H1P-KPQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       dim_countf=kpqf(1,0)

       do i= dim_countf+1,dim_countf+kpqf(2,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j= 1,ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)    

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_countf=dim_countf+kpqf(2,0)

       do i= dim_countf+1,dim_countf+kpqf(3,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j= 1,ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_countf=dim_countf+kpqf(3,0)

       do i= dim_countf+1,dim_countf+kpqf(4,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j= 1,ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_countf=dim_countf+kpqf(4,0)

       do i= dim_countf+1,dim_countf+kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j= 1,ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_countf=dim_countf+kpqf(5,0)

       do i= dim_countf+1,dim_countf+kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j= 1,ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpqf(1,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (2,1) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (1,2) block 

    lim1i=kpqf(1,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do


!!$ (3,1) block
     
    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do          

!!$ (1,3) block

    lim1i=kpqf(1,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do          
 
!!$ (4i,1) block

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (i,4i) block

    lim1i=kpqf(1,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (1,4ii) block

    lim1i=kpqf(1,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (2,2) block

    lim1i=kpqf(1,0)+kpqf(2,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!!!write(6,*) 'D_2_2_2p2h_2p2h_offdiag', ar_offdiagd(i,j)
           
       end do
    end do

!!$ (3,2) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (2,3) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do


!!$ (2,4i) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (2,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (3,3) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
       end do
    end do

!!$ (3,4i) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (3,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4i,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

  end subroutine get_offdiag_adc2_DIPOLE_direct_OK











!!! CASE OF INITIAL AND FINAL SPACES OF THE SAME SYMMETRY
!!! CASE OF INITIAL AND FINAL SPACES OF THE SAME SYMMETRY
!!! CASE OF INITIAL AND FINAL SPACES OF THE SAME SYMMETRY


!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------



!!$------------------------------------------------------------------------------
!!$-----------------------------  OFF-DIAGONAL PART  ----------------------------
!!$------------------------------------------------------------------------------
  subroutine get_offdiag_tda_DIPOLE_direct(ndim,kpq,ar_offdiagd)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiagd
    
    integer :: inda,indb,indk,indl,spin
    integer :: indapr,indbpr,indkpr,indlpr,spinpr 
    
    integer :: i,j,dim_count,ndim1
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiagd
   
    integer :: lim1i,lim2i,lim1j,lim2j
 
    ar_offdiagd(:,:)=0._d
    
!!$ Full diagonalization. 

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             


    if (indk .eq. indkpr)   then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr)    then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_2_ph_ph(indk,indkpr)
    end if

          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do


  end  subroutine get_offdiag_tda_DIPOLE_direct
 



!!$------------------------------------------------------------------------------
!!$---------------------------------  DIAGONAL PART  ----------------------------
!!$------------------------------------------------------------------------------
  subroutine get_diag_tda_DIPOLE_direct(ndim1,kpq,ar_diagd)
  
    integer, intent(in) :: ndim1
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1), intent(out) :: ar_diagd
    
    integer :: inda,indb,indk,indl,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,lim1,lim2
    
!!$ Filling the ph-ph block
   

    ar_diagd(:)=0.0
 
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

 
                   ar_diagd(i) = ar_diagd(i) + D0_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D0_2_ph_ph(indk,indk)


!!!write(6,*) "D0_1",D0_1_ph_ph(inda,inda)            
!!!write(6,*) "D0_2",D0_2_ph_ph(indk,indk)            
!!!write(6,*) "Ddiag_ph_ph", ar_diagd(i)
 
   end do

  end subroutine get_diag_tda_DIPOLE_direct









!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
 


!!$------------------------------------------------------------------------------
!!$-----------------------------  OFF-DIAGONAL PART  ----------------------------
!!$------------------------------------------------------------------------------
  subroutine get_offdiag_adc2_DIPOLE_direct(ndim,kpq,ar_offdiagd)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiagd
    
    integer :: inda,indb,indk,indl,spin
    integer :: indapr,indbpr,indkpr,indlpr,spinpr 
    
    integer :: i,j,dim_count,ndim1
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiagd
   
    integer :: lim1i,lim2i,lim1j,lim2j
 
    ar_offdiagd(:,:)=0._d
    
!!$ Full diagonalization. 

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             


         ar_offdiagd(i,j) = D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if (indk .eq. indkpr)   then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr)    then
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_5_2_ph_ph(indk,indkpr)
    end if

          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do

       
!!$ Filling the off-diagonal part of the ph-2p2h block Dak,a'b'k'l'
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


if ((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

    
          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do 
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)


          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)


    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)


if(inda .eq. indapr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd(i,j) = ar_offdiagd(i,j) + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       end do
    end do
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=i+1,lim2i
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
          ar_offdiagd(j,i) = ar_offdiagd(i,j)  

       
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (1,2) block 

    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do


!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do          

!!$ (1,3) block

    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do          
 
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (i,4i) block

    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (1,4ii) block

    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= i+1,lim2i
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
          ar_offdiagd(j,i) = ar_offdiagd(i,j)  

!!!write(6,*) 'D_2_2_2p2h_2p2h_offdiag', ar_offdiagd(i,j)
           
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (2,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do


!!$ (2,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (2,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= i+1,lim2i
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
          ar_offdiagd(j,i) = ar_offdiagd(i,j)         
           
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
       end do
    end do

!!$ (3,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (3,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= i+1,lim2i
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
          ar_offdiagd(j,i) = ar_offdiagd(i,j)         
    
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do

!!$ (4i,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= i+1,lim2i
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd(i,j) =  ar_offdiagd(i,j) + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
          ar_offdiagd(j,i) = ar_offdiagd(i,j)         
 
       end do
    end do



 
  end subroutine get_offdiag_adc2_DIPOLE_direct





!!$------------------------------------------------------------------------------
!!$---------------------------------  DIAGONAL PART  ----------------------------
!!$------------------------------------------------------------------------------
  subroutine get_diag_adc2_DIPOLE_direct(ndim1,ndim2,kpq,ar_diagd)
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1+ndim2), intent(out) :: ar_diagd
    
    integer :: inda,indb,indk,indl,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,lim1,lim2
    
!!$ Filling the ph-ph block
   

    ar_diagd(:)=0.0
 
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

         ar_diagd(i) = ar_diagd(i) + D2_6_1_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_2_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_3_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_4_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_7_1_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_7_2_ph_ph(inda,inda,indk,indk)

 
                   ar_diagd(i) = ar_diagd(i) + D0_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_2_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_3_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_3_2_ph_ph(inda,inda)
 

                   ar_diagd(i) = ar_diagd(i) + D0_2_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_3_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_1_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_2_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_5_1_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_5_2_ph_ph(indk,indk)



    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_1_1_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)


    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_2_2_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

!!!write (6,*) "Ddiag_2_2_2p2h_2p2h",inda,indb,indk,indl, ar_diagd(i)

    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_3_3_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
  end subroutine get_diag_adc2_DIPOLE_direct



































!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!


!!! ---------------  ROUTINES THAT SAVE THE DIPOLE MATRIX ON FILE  ------------- !!! 
!!! ---------------  ROUTINES THAT SAVE THE DIPOLE MATRIX ON FILE  ------------- !!! 
!!! ---------------  ROUTINES THAT SAVE THE DIPOLE MATRIX ON FILE  ------------- !!! 





!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

!!$------------------------------------------------------------------------------
!!$--------------------------- OFF-DIAGONAL ROUTINES ----------------------------
!!$------------------------------------------------------------------------------


  subroutine get_offdiag_tda_DIPOLE_SAVE_OK(ndimf,ndim,kpqf,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndimf
  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  


  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   

  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1f = kpqf(1,0)

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

        do j = 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

          
           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

   end do
  end do



    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2



  end subroutine get_offdiag_tda_DIPOLE_SAVE_OK



  subroutine get_offdiag_tda_DIPOLE_SAVE_OK_GS(ndimf,ndim,kpqf,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndimf
  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  


  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j


  integer :: nlim1 , nlim2 , a , k

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix




!!$-----1h1p block------
   
    nlim1 = 1
    nlim2 = kpq(1,0)
    
    do j = nlim1 , nlim2

       i = 0
       k=kpq(3,j)
       a=kpq(5,j)

     ar_offdiagd_ij  =  dpl(a,k) + F0_ph(a,k) 
     ar_offdiagd_ij  =  -sqrt(2._d) * ar_offdiagd_ij

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if
       
    end do


!!$ Filling the off-diagonal part of the ph-ph block

     ndim1f = kpqf(1,0)

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

        do j = 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

          
           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

   end do
  end do



    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8)) = i + 1
         oj(count-buf_size*int(rec_count,8)) = j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2



  end subroutine get_offdiag_tda_DIPOLE_SAVE_OK_GS



  subroutine get_offdiag_tda_DIPOLE_SAVE_OK_SAME(ndim,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  


  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block


     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

          
           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

   end do
  end do



    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2



  end subroutine get_offdiag_tda_DIPOLE_SAVE_OK_SAME






  subroutine get_offdiag_tda_DIPOLE_SAVE_OK_SAME_GS(ndim,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  


  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block


     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

          
           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

   end do
  end do



    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j+1,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))= i + 1
         oj(count-buf_size*int(rec_count,8))= j + 1
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2



  end subroutine get_offdiag_tda_DIPOLE_SAVE_OK_SAME_GS










!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------


!!$------------------------------------------------------------------------------
!!$--------------------------- OFF-DIAGONAL ROUTINES ----------------------------
!!$------------------------------------------------------------------------------

 
  subroutine get_offdiag_adc2_DIPOLE_SAVE_OK(ndimf,ndim,kpqf,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndimf
  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  
  
  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1f = kpqf(1,0)

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

        do j = 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indkpr) 
    end if
             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


   end do
  end do


     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs


       dim_count = kpq(1,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count = dim_count + kpq(2,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_count = dim_count + kpq(3,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b I configs


       
       dim_count = dim_count + kpq(4,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b II configs


       
       dim_count = dim_count + kpq(5,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
   





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!2H2P-KPQF      1H1P-KPQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




       dim_countf = kpqf(1,0)

       do i = dim_countf + 1 , dim_countf + kpqf(2,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j = 1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
          

!!$ Coupling to the i=j,a|=b configs   
       
       dim_countf = dim_countf + kpqf(2,0)

       do i = dim_countf + 1 , dim_countf + kpqf(3,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_countf = dim_countf + kpqf(3,0)

       do i = dim_countf + 1 , dim_countf + kpqf(4,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a|=b I configs

       

       dim_countf = dim_countf + kpqf(4,0)

       do i = dim_countf + 1 , dim_countf + kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a|=b II configs

       
       dim_countf = dim_countf + kpqf(5,0)

       do i = dim_countf + 1 , dim_countf + kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
        end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,1) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (1,2) block 

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (3,1) block
     
    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          

!!$ (1,3) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          
 
!!$ (4i,1) block

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (i,4i) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (1,4ii) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (2,2) block

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)


             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (2,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (3,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (4ii,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do




    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2


  end subroutine get_offdiag_adc2_DIPOLE_SAVE_OK







  subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_GS(ndimf,ndim,kpqf,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndimf
  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  
  
  integer :: i,j,nlim,dim_count,dim_countf,ndim1,ndim1f
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: nlim1 , nlim2 , a , k , b , jmuto

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix




    
!!$-----1h1p block------
   
    nlim1 = 1
    nlim2 = kpq(1,0)
    
    do j = nlim1 , nlim2

       i = 0
       k=kpq(3,j)
       a=kpq(5,j)

     ar_offdiagd_ij  =  dpl(a,k) + F0_ph(a,k) + FA_ph(a,k) + FB_ph(a,k) + FC_ph(a,k)
     ar_offdiagd_ij  =  ar_offdiagd_ij + F21_ph(a,k) + F22_ph(a,k) + F23_ph(a,k) + F24_ph(a,k) + F25_ph(a,k)
     ar_offdiagd_ij  =  ar_offdiagd_ij + F26_ph(a,k) + F27_ph(a,k) + F28_ph(a,k) + F29_ph(a,k) + F210_ph(a,k)
     ar_offdiagd_ij  =  -sqrt(2._d) * ar_offdiagd_ij

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if
       
    end do

!!$----------I-a=b,i=j-------------------

    nlim1 = nlim2 + 1
    nlim2 = nlim2 + kpq(2,0)
    
    do j = nlim1 , nlim2

       i = 0
       k=kpq(3,j)
       a=kpq(5,j)

       ar_offdiagd_ij = FI_2p2h(a,k)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

    end do
   
       
!!$----------II-a|=b,i=j-------------------

    nlim1 = nlim2 + 1
    nlim2 = nlim2 + kpq(3,0)
    
    do j = nlim1 , nlim2
       
       i = 0
       k=kpq(3,j)
       a=kpq(5,j)
       b=kpq(6,j)

       ar_offdiagd_ij = FII_2p2h(a,b,k)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

    end do    
       
!!$----------III-a=b,i|=j-------------------

    nlim1 = nlim2 + 1
    nlim2 = nlim2 + kpq(4,0)
    
    do j = nlim1 , nlim2
       
       i = 0
       k=kpq(3,j)
       jmuto=kpq(4,j)
       a=kpq(5,j)

       ar_offdiagd_ij = FIII_2p2h(a,k,jmuto)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

    end do    
           
!!$----------IV1-a|=b,i|=j-------------------

    nlim1 = nlim2 + 1
    nlim2 = nlim2 + kpq(5,0)
    
    do j = nlim1 , nlim2
       
       i = 0
       k=kpq(3,j)
       jmuto=kpq(4,j)
       a=kpq(5,j)
       b=kpq(6,j)

       ar_offdiagd_ij = FIV1_2p2h(a,b,k,jmuto)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

    end do 

!!$----------IV2-a|=b,i|=j-------------------

    nlim1 = nlim2 + 1
    nlim2 = nlim2 + kpq(5,0)

    do j = nlim1 , nlim2
       
       i = 0
       k=kpq(3,j)
       jmuto=kpq(4,j)
       a=kpq(5,j)
       b=kpq(6,j)

       ar_offdiagd_ij = FIV2_2p2h(a,b,k,jmuto)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

    end do





!!$ Filling the off-diagonal part of the ph-ph block



     ndim1f = kpqf(1,0)

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1f
        call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

        do j = 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indkpr) 
    end if
             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


   end do
  end do


     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs


       dim_count = kpq(1,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count = dim_count + kpq(2,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_count = dim_count + kpq(3,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b I configs


       
       dim_count = dim_count + kpq(4,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b II configs


       
       dim_count = dim_count + kpq(5,0)

       do i = 1 , ndim1f
          call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
   





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!2H2P-KPQF      1H1P-KPQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




       dim_countf = kpqf(1,0)

       do i = dim_countf + 1 , dim_countf + kpqf(2,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j = 1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
          

!!$ Coupling to the i=j,a|=b configs   
       
       dim_countf = dim_countf + kpqf(2,0)

       do i = dim_countf + 1 , dim_countf + kpqf(3,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_countf = dim_countf + kpqf(3,0)

       do i = dim_countf + 1 , dim_countf + kpqf(4,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a|=b I configs

       

       dim_countf = dim_countf + kpqf(4,0)

       do i = dim_countf + 1 , dim_countf + kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a|=b II configs

       
       dim_countf = dim_countf + kpqf(5,0)

       do i = dim_countf + 1 , dim_countf + kpqf(5,0)
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spin)

          do j =  1 , ndim1
             call get_indices(kpq(:,j),inda,indb,indk,indl,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
        end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,1) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (1,2) block 

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (3,1) block
     
    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          

!!$ (1,3) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          
 
!!$ (4i,1) block

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (i,4i) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (1,4ii) block

    lim1i = kpqf(1,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

!!$ (2,2) block

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)


             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (2,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,2) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (3,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (4ii,3) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,4ii) block 

    lim1i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + 1
    lim2i = kpqf(1,0) + kpqf(2,0) + kpqf(3,0) + kpqf(4,0) + kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+1
    lim2i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)+kpqf(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do




    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))= i + 1
         oj(count-buf_size*int(rec_count,8))= j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2


  end subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_GS



  subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_SAME(ndim,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block


     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indkpr) 
    end if
             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


   end do
  end do


     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs


       dim_count = kpq(1,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count = dim_count + kpq(2,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_count = dim_count + kpq(3,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b I configs


       
       dim_count = dim_count + kpq(4,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b II configs


       
       dim_count = dim_count + kpq(5,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpq(1,0) + 1
    lim2i = kpq(1,0) + kpq(2,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j=lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,1) block 

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do



!!$ (3,1) block
     
    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          

 
!!$ (4i,1) block

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

 
!!$ (4ii,1) block

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 


!!$ (2,2) block

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)


             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

        
!!$ (4i,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (4ii,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (3,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (4ii,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (4i,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4ii) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do




    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2


  end subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_SAME



  subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_SAME_GS(ndim,kpq,nbuf,count, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  integer*8, intent(out) :: count
  integer, intent(out)   :: nbuf
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij  
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

  integer :: rec_count
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  count=0
  rec_count=0
  nbuf=0
  oi(:)=0   
  oj(:)=0   
  file_offdiagd(:)=0.d0   
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block


     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

    if(indk .eq. indkpr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,indapr)
    end if

    if(inda .eq. indapr) then
                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indkpr) 
    end if
             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


   end do
  end do


     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs


       dim_count = kpq(1,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

           ar_offdiagd_ij = 0.d0

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count = dim_count + kpq(2,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do


!!$ Coupling to the i|=j,a=b configs

       
       dim_count = dim_count + kpq(3,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr).and.(inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b I configs


       
       dim_count = dim_count + kpq(4,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do



!!$ Coupling to the i|=j,a|=b II configs


       
       dim_count = dim_count + kpq(5,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          do j = dim_count + 1 , dim_count + kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

           ar_offdiagd_ij = 0.d0

 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)

if(inda .eq. indapr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
if(inda .eq. indbpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
if(indk .eq. indkpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
if(indk .eq. indlpr)&
ar_offdiagd_ij = ar_offdiagd_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if


          end do
       end do
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpq(1,0) + 1
    lim2i = kpq(1,0) + kpq(2,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j=lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (2,1) block 

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do



!!$ (3,1) block
     
    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do          

 
!!$ (4i,1) block

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 

 
!!$ (4ii,1) block

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do 


!!$ (2,2) block

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
   
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)


             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (3,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

        
!!$ (4i,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (4ii,2) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (3,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4i,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

            
       end do
    end do

!!$ (4ii,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do


!!$ (4i,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do

!!$ (4ii,4ii) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)

    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j= lim1j , i - 1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

             !Culling  small matrix elements
             if (abs(ar_offdiagd_ij) .gt. minc) then
                call register1()
             end if

           
       end do
    end do




    call register2()


    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j+1,'one:', ar_offdiagd_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8)) = i + 1
         oj(count-buf_size*int(rec_count,8)) = j + 1
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_DIP ,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2


  end subroutine get_offdiag_adc2_DIPOLE_SAVE_OK_SAME_GS




















!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

!!$------------------------------------------------------------------------------
!!$------------------------------ DIAGONAL ROUTINES -----------------------------
!!$------------------------------------------------------------------------------


  subroutine get_diag_tda_DIPOLE_OK_SAME(ndim,kpq,ar_diag, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real*8, dimension(ndim) , intent(out) :: ar_diag
  INTEGER, INTENT(IN) :: UNIT_DIP

  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_diagd_ij, Dground_0


  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ar_diag(:) = 0.d0

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          
           ar_diagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ar_diagd_ij = ar_diagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_diagd_ij = ar_diagd_ij + D0_2_ph_ph(indk,indkpr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ar_diag(i) = ar_diagd_ij


  end do

    CALL get_Dground_0(Dground_0)
    ar_diag(:) = ar_diag(:) + Dground_0 

    write(UNIT_DIP) ar_diag
    write(*,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_diag(1)
    write(*,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_diag( ndim )



  end subroutine get_diag_tda_DIPOLE_OK_SAME


  subroutine get_diag_tda_DIPOLE_OK_SAME_GS(ndim,kpq,ar_diag, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real*8, dimension( ndim + 1 ) , intent(out) :: ar_diag
  INTEGER, INTENT(IN) :: UNIT_DIP

  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_diagd_ij , Dground_0  


  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ar_diag(:) = 0.d0

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          
           ar_diagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ar_diagd_ij = ar_diagd_ij + D0_1_ph_ph(inda,indapr)
                   ar_diagd_ij = ar_diagd_ij + D0_2_ph_ph(indk,indkpr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ar_diag( i + 1 ) = ar_diagd_ij


  end do

    CALL get_Dground_0(Dground_0)
    ar_diag(:) = ar_diag(:) + Dground_0 


    write (UNIT_DIP) ar_diag
    write(*,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_diag(1)
    write(*,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_diag( ndim + 1 )

  end subroutine get_diag_tda_DIPOLE_OK_SAME_GS









!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------ ADC2 ------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

!!$------------------------------------------------------------------------------
!!$------------------------------ DIAGONAL ROUTINES -----------------------------
!!$------------------------------------------------------------------------------


  subroutine get_diag_adc2_DIPOLE_OK_SAME(ndim,kpq,ar_offdiag, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real*8, dimension(ndim), intent(out) :: ar_offdiag
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij, Dground_0, Dground_2  
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ar_offdiag(:) = 0.d0

     CALL get_Dground_0(Dground_0)
     CALL get_Dground_2(Dground_2)

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,inda,indk,indk)

                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,inda)

                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indk) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_2


  end do


     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpq(1,0) + 1
    lim2i = kpq(1,0) + kpq(2,0)

    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_0 

    end do



!!$ (2,2) block

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

   
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_0 


           
    end do


           

!!$ (3,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_0 

           
    end do


!!$ (4i,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_0 

           
    end do

!!$ (4ii,4ii) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)


    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag(i) = ar_offdiagd_ij
              ar_offdiag(i) = ar_offdiag(i) + Dground_0 

           
    end do

    write(UNIT_DIP) ar_offdiag
    write(*,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_offdiag(1)
    write(*,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_offdiag( ndim )

  end subroutine get_diag_adc2_DIPOLE_OK_SAME



  subroutine get_diag_adc2_DIPOLE_OK_SAME_GS(ndim,kpq,ar_offdiag, UNIT_DIP )

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real*8, dimension( ndim + 1 ), intent(out) :: ar_offdiag
  integer, intent(in) :: UNIT_DIP
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  real(d) :: ar_offdiagd_ij, Dground_0, Dground_2  
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i,lim2i,lim1j,lim2j

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ar_offdiag(:) = 0.d0

     CALL get_Dground_0(Dground_0)
     CALL get_Dground_2(Dground_2)

     ar_offdiag(1) = Dground_2

     ndim1 = kpq(1,0)
  
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


           ar_offdiagd_ij = 0.d0
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_1_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_2_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_3_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_6_4_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_1_ph_ph(inda,inda,indk,indk)
         ar_offdiagd_ij = ar_offdiagd_ij + D2_7_2_ph_ph(inda,inda,indk,indk)

                   ar_offdiagd_ij = ar_offdiagd_ij + D0_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_2_2_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_1_ph_ph(inda,inda)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_2_ph_ph(inda,inda)

                   ar_offdiagd_ij = ar_offdiagd_ij + D0_2_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_3_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_1_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_4_2_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_1_ph_ph(indk,indk)
                   ar_offdiagd_ij = ar_offdiagd_ij + D2_5_2_ph_ph(indk,indk) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_2

  end do


     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i = kpq(1,0) + 1
    lim2i = kpq(1,0) + kpq(2,0)

    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_0 

    end do



!!$ (2,2) block

    lim1i = kpq(1,0) + kpq(2,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0)

   
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_0 


           
    end do


           

!!$ (3,3) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0)

    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_0 

           
    end do


!!$ (4i,4i) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0)

    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)


          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_0 

           
    end do

!!$ (4ii,4ii) block 

    lim1i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + 1
    lim2i = kpq(1,0) + kpq(2,0) + kpq(3,0) + kpq(4,0) + kpq(5,0) + kpq(5,0)


    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

          ar_offdiagd_ij = 0.d0
          ar_offdiagd_ij =  ar_offdiagd_ij + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

              ar_offdiag( i + 1 ) = ar_offdiagd_ij
              ar_offdiag( i + 1 ) = ar_offdiag( i + 1 ) + Dground_0 

           
    end do

    write(UNIT_DIP) ar_offdiag

    write(*,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_offdiag(1)
    write(*,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_offdiag( ndim + 1 )

  end subroutine get_diag_adc2_DIPOLE_OK_SAME_GS





!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!
!!! ******************************* SAVE ROUTINES ****************************** !!!



























































!!! ROUTINES FOR INITIAL AND FINAL SPACE HAVING THE SAME SYMMETRY, 
!!! ,OR FOR CALCULATIONS NOT USING ANY SYMMETRY 



!!$-----------------------------------------------------------
!!$-----------------------------------------------------------




subroutine get_offdiag_adc2_DIPOLE_save(ndim,kpq,nbuf,count,chr)
   
  integer, intent(in) :: ndim
  integer*8, intent(out) :: count
  integer, intent(out) :: nbuf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  character(1), intent(in) :: chr
  
  integer :: inda,indb,indk,indl,spin
  integer :: indapr,indbpr,indkpr,indlpr,spinpr 
  
  character(13) :: name
  integer :: rec_count
  integer :: i,j,nlim,dim_count,ndim1,unt
  integer :: lim1i, lim2i, lim1j, lim2j
  real(d) :: ar_offdiag_ij
  
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiagd

  name="hmlt.off"//chr
  unt=22
  
  count=0
  rec_count=0
  
  write(6,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", name
  OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
       FORM='UNFORMATTED')


!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i=1,ndim1
        call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
        do j=1,i-1
           call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             

         ar_offdiag_ij = D2_6_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_2_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_3_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_6_4_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_7_1_ph_ph(inda,indapr,indk,indkpr)
         ar_offdiag_ij = ar_offdiag_ij + D2_7_2_ph_ph(inda,indapr,indk,indkpr)

      if(indk .eq. indkpr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_2_2_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_1_ph_ph(inda,indapr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_2_ph_ph(inda,indapr)
      end if

      if(inda .eq. indapr) then
                   ar_offdiag_ij = ar_offdiag_ij + D0_2_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_3_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_1_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_4_2_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_5_1_ph_ph(indk,indkpr)
                   ar_offdiag_ij = ar_offdiag_ij + D2_5_2_ph_ph(indk,indkpr)
      end if

           call register1()
        end do
     end do
   
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
          do j=dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    

ar_offdiag_ij = 0.

if((indk .eq. indkpr).and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D5_1_ph_2p2h(inda,indk,indbpr,indlpr) + D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D5_2_ph_2p2h(inda,indk,indbpr,indkpr) + D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_3_ph_2p2h(inda,indk,indapr,indlpr) + D5_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr)  .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D5_4_ph_2p2h(inda,indk,indapr,indkpr) + D5_8_ph_2p2h(inda,indk,indapr,indkpr)



if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             call register1()
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
          do j=dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
ar_offdiag_ij = 0.


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D4_1_ph_2p2h(inda,indk,indbpr,indlpr) + D4_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D4_2_ph_2p2h(inda,indk,indbpr,indkpr) + D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_3_ph_2p2h(inda,indk,indapr,indlpr) + D4_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D4_4_ph_2p2h(inda,indk,indapr,indkpr) + D4_8_ph_2p2h(inda,indk,indapr,indkpr)



if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)


             !Culling  small matrix elements
             if (abs(ar_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
ar_offdiag_ij = 0.


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D3_1_ph_2p2h(inda,indk,indbpr,indlpr) + D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D3_2_ph_2p2h(inda,indk,indbpr,indkpr) + D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_3_ph_2p2h(inda,indk,indapr,indlpr) + D3_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D3_4_ph_2p2h(inda,indk,indapr,indkpr) + D3_8_ph_2p2h(inda,indk,indapr,indkpr)



if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             call register1()
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
ar_offdiag_ij = 0.


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D1_1_ph_2p2h(inda,indk,indbpr,indlpr) + D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D1_2_ph_2p2h(inda,indk,indbpr,indkpr) + D1_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_3_ph_2p2h(inda,indk,indapr,indlpr) + D1_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D1_4_ph_2p2h(inda,indk,indapr,indkpr) + D1_8_ph_2p2h(inda,indk,indapr,indkpr)



if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)


             !Culling  small matrix elements
             if (abs(ar_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
ar_offdiag_ij = 0.


 if((indk .eq. indkpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D2_1_ph_2p2h(inda,indk,indbpr,indlpr) + D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indapr))&
ar_offdiag_ij = ar_offdiag_ij + D2_2_ph_2p2h(inda,indk,indbpr,indkpr) + D2_6_ph_2p2h(inda,indk,indbpr,indkpr)

if((indk .eq. indkpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_3_ph_2p2h(inda,indk,indapr,indlpr) + D2_7_ph_2p2h(inda,indk,indapr,indlpr)

if((indk .eq. indlpr) .and. (inda .eq. indbpr))&
ar_offdiag_ij = ar_offdiag_ij + D2_4_ph_2p2h(inda,indk,indapr,indkpr) + D2_8_ph_2p2h(inda,indk,indapr,indkpr)



if(inda .eq. indapr)&
ar_offdiag_ij = ar_offdiag_ij + D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

if(inda .eq. indbpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

if(indk .eq. indkpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

if(indk .eq. indlpr)&
ar_offdiag_ij = ar_offdiag_ij + D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

             !Culling  small matrix elements
             call register1()
          end do
       end do
    
!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
          call register1()
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
            
          call register1()
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)
ar_offdiag_ij = 0.


          ar_offdiag_ij = ar_offdiag_ij + D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr) 
ar_offdiag_ij = 0.


          ar_offdiag_ij =ar_offdiag_ij +  D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
           
          call register1()
       end do
    end do
    
    call register2()
    CLOSE(unt)
    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' DIPOLE_off-diagonal elements saved in file ', name

  contains
       
    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
!!$            call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)  
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
!!$      call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2_DIPOLE_save
 

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------


  subroutine get_diag_adc2_DIPOLE_save(ndim1,ndim2,kpq,nbuf,chr)
  
     integer, intent(in) :: ndim1,ndim2,nbuf 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
   
    integer :: inda,indb,indj,indk,indl,spin
    real(d) ::ea,eb,ej,ek,temp
    
    character(13) :: name
    integer :: i,ktype,dim_count,lim1,lim2,unt,a,b,c,d1
    real(d), dimension(ndim1+ndim2) :: ar_diagd
     
    ktype=1
    name="hmlt.dia"//chr 
    unt=21
    
!!$ Filling the ph-ph block
   

    ar_diagd(:)=0.0
 
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

         ar_diagd(i) = ar_diagd(i) + D2_6_1_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_2_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_3_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_6_4_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_7_1_ph_ph(inda,inda,indk,indk)
         ar_diagd(i) = ar_diagd(i) + D2_7_2_ph_ph(inda,inda,indk,indk)

 
                   ar_diagd(i) = ar_diagd(i) + D0_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_2_2_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_3_1_ph_ph(inda,inda)
                   ar_diagd(i) = ar_diagd(i) + D2_3_2_ph_ph(inda,inda)
 

                   ar_diagd(i) = ar_diagd(i) + D0_2_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_3_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_1_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_4_2_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_5_1_ph_ph(indk,indk)
                   ar_diagd(i) = ar_diagd(i) + D2_5_2_ph_ph(indk,indk)


    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_1_1_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)


    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_2_2_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

!!!write (6,*) "Ddiag_2_2_2p2h_2p2h",inda,indb,indk,indl, ar_diagd(i)

    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_3_3_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin) 

       ar_diagd(i) = ar_diagd(i) + D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,inda,indb,indk,indl)

    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing",ndim1+ndim2," diagonal elements of ADC-DIPOLE ADC2 matrix in file ",name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diagd(:))
!!$    call wrtdgat(unt,ndim1+ndim2,nbuf,ar_diag(:))
    CLOSE(unt)
   
    write(*,*) 'Writing successful at get_diag_adc2_DIPOLE_save end'
 
  end subroutine get_diag_adc2_DIPOLE_save

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------


 
end module get_matrix_DIPOLE    
          
          
       
       
       
       
    
    
    
       
       
       
       
       
       
    
    
    
    

    
    
    
    
