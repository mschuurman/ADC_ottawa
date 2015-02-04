module get_matrix

  
  use constants
  use parameters
  use adc_ph
  use misc
  use filetools
  
  implicit none

  integer, parameter :: buf_size=8192


contains
  
!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************TDA BLOCK***********************************
!!$*******************************************************************************
!!$*******************************************************************************
  
  subroutine get_diag_tda_direct(ndim,kpq,ar_diag)
    
    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim), intent(out) :: ar_diag
    
    integer :: i
    integer :: inda,indb,indj,indk,spin

!!$ Preparing the diagonal part for the full diagonalisation   
    
    do i= 1,ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i)=K_ph_ph(e(inda),e(indj))
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
    end do

  end subroutine get_diag_tda_direct

!!$-------------------------------------------------------------------------  
!!$-------------------------------------------------------------------------

  subroutine get_offdiag_tda_direct(ndim,kpq,ar_offdiag)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiag 
    
    integer :: i,j
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr

!!$ Full diagonalization. The program performs a symmetry check. For better efficiency 
!!$ it should be removed if no problems have showed up.
 
    do i=1,ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag(i,j)=C1_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag(j,i)=C1_ph_ph(indapr,indjpr,inda,indj)
!          if(abs(ar_offdiag(i,j)-ar_offdiag(j,i)) .ge. 1.e-14_d) then
          if(abs(ar_offdiag(i,j)-ar_offdiag(j,i)) .ge. 1.e-12_d) then
             write(6,*) "TDA matrix is not symmetric. Stopping now."
             print*,i,j,abs(ar_offdiag(i,j)-ar_offdiag(j,i))
             stop
          end if
       end do
    end do
    
  end subroutine get_offdiag_tda_direct





!!$---------------------------------------------------
!!$---------------------------------------------------










  subroutine get_diag_tda_save(ndim,kpq,nbuf,chr)
    
    integer, intent(in) :: ndim,nbuf 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    
    character(10) :: name
    integer :: i,ktype,unt
    integer :: inda,indb,indj,indk,spin
    real(d), dimension(ndim) :: ar_diag
    
    ktype=1
    name="tda.diag"//chr
    unt=11
    
    do i=1, ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i)=K_ph_ph(e(inda),e(indj))
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
    end do

       !Saving in file
       write(6,*) "Writing the diagonal part in file ", name
       OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
            FORM='UNFORMATTED')
       call wrtdg(unt,ndim,buf_size,nbuf,ktype,ar_diag(:))
       CLOSE(unt)
       
  end subroutine get_diag_tda_save

!!$----------------------------------------------------------------
!!$----------------------------------------------------------------

  subroutine get_offdiag_tda_save(ndim,kpq,nbuf,chr)
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,count,unt
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr
    real(d) :: ar_offdiag_ij, ar_offdiag_ji

    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag

    
    name="tda.off"//chr
    unt=12

    count=0
    rec_count=0

       write(6,*) "Writing the off-diagonal part of TDA matrix in file ", name
       OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
            FORM='UNFORMATTED')
    
    do i=1,ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag_ji=C1_ph_ph(indapr,indjpr,inda,indj)
          if(abs(ar_offdiag_ij-ar_offdiag_ji) .ge. 1.e-15_d) then
             write(6,*) "TDA matrix is not symmetric. Stopping now."
             stop
          end if

!!$ Saving into vector for the following Lanzcos/Davidson routine 
            

          !Culling  small matrix elements
          if (abs(ar_offdiag_ij) .gt. minc) then
             call register1()
          end if

       end do
    end do
!!$

       call register2()
       CLOSE(unt)
       
  contains
       
    subroutine register1()
      
      count=count+1
      file_offdiag(count-buf_size*rec_count)=ar_offdiag_ij
      oi(count-buf_size*rec_count)=i
      oj(count-buf_size*rec_count)=j
      !Checking if the buffer is full 
      if(mod(count,buf_size) .eq. 0) then
         rec_count=rec_count+1
         nlim=buf_size
         !Saving off-diag part in file
         call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
      end if

    end subroutine register1
       
    subroutine register2()
         
      !Saving the rest of matrix in file
      nlim=count-buf_size*rec_count
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_tda_save




!!$---------------------------------------------------
!!$---------------------------------------------------










!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************ADC2 BLOCK***********************************
!!$*******************************************************************************
!!$*******************************************************************************

  subroutine get_diag_adc2_direct(ndim1,ndim2,kpq,ar_diag)
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1+ndim2), intent(out) :: ar_diag
    
    integer :: inda,indb,indj,indk,spin
    integer :: i
    
!!$ Filling the ph-ph block

    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i)=K_ph_ph(e(inda),e(indj))
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do
    
!!$ Filling the 2p2h-2p2h block
    
    do i=ndim1+1, ndim1+ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(i)=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
  end subroutine get_diag_adc2_direct

!!$------------------------------------------------------------------------
!!$------------------------------------------------------------------------

  subroutine get_diag_adc2_save(ndim1,ndim2,kpq,nbuf,chr)
  
    integer, intent(in) :: ndim1,ndim2,nbuf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: inda,indb,indj,indk,spin
    
    character(10) :: name
    integer :: i,ktype,unt 
    real(d), dimension(:), allocatable:: ar_diag

    allocate(ar_diag(ndim1+ndim2))

    ktype=1
    name="hmlt.dia"//chr 
    unt=11

!!$ Filling the ph-ph block

    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i)=K_ph_ph(e(inda),e(indj))
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
    
    do i=ndim1+1, ndim1+ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(i)=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing the diagonal part of ADC matrix in file ", name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diag(:))
    CLOSE(unt)

    deallocate(ar_diag)
  end subroutine get_diag_adc2_save

!!$----------------------------------------------------------------------------
!!$----------------------------------------------------------------------------
  
  subroutine get_offdiag_adc2_direct(ndim,kpq,ar_offdiag)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiag
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    integer :: i,j,dim_count,ndim1
    
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    
    ar_offdiag(:,:)=0._d
    
!!$ Full diagonalization. 

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
          ar_offdiag(i,j)=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
                  ar_offdiag(i,j)=ar_offdiag(i,j)+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offdiag(i,j)=ar_offdiag(i,j)+CB_ph_ph(indj,indjpr)
          ar_offdiag(i,j)=ar_offdiag(i,j)+CC_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do

       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
          ar_offdiag(i,j)=C5_ph_2p2h(inda,indj,indapr,indjpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(i,j)=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(i,j)=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(i,j)=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(i,j)=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(j,i)=ar_offdiag(i,j)
       end do
    end do
  
  end subroutine get_offdiag_adc2_direct
!!$  
!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------
!!$-----------------------------------

  subroutine get_offdiag_adc2_save(ndim,kpq,nbuf,count,chr)

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf
    integer*8, intent(out) :: count 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,dim_count,ndim1,unt
    real(d) :: ar_offdiag_ij
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    
    name="hmlt.off"//chr
    unt=12

    count=0
    rec_count=0
    
    
    write(6,*) "Writing the off-diagonal part of ADC matrix in file ", name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

!!$ Full diagonalization.  

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
               ar_offdiag_ij= ar_offdiag_ij+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offdiag_ij= ar_offdiag_ij+CB_ph_ph(indj,indjpr)
          ar_offdiag_ij= ar_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
          call register1()
       end do
    end do
    
       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
          ar_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
          call register1()
       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
          call register1()
       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
          call register1()
       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

    call register2()
    CLOSE(unt)
    write(6,*) count,' off-diagonal elements saved'

  contains
    
    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
         file_offdiag(count-buf_size*int(rec_count,8))= ar_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_adc2_save

!#######################################################################

  subroutine get_offdiag_adc2_save_cvs(ndim,kpq,nbuf,count,chr)

    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf
    integer*8, intent(out) :: count 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,dim_count,ndim1,unt
    real(d) :: ar_offdiag_ij
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    
    name="hmlt.off"//chr
    unt=12

    count=0
    rec_count=0

    write(6,*) "Writing the off-diagonal part of ADC matrix in file ", name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

!-----------------------------------------------------------------------
! Off-diagonal part of the ph-ph block: all admissible due to previous
! screening of the configurations
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
               ar_offdiag_ij= ar_offdiag_ij+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offdiag_ij= ar_offdiag_ij+CB_ph_ph(indj,indjpr)
          ar_offdiag_ij= ar_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
          call register1()
       end do
    end do

!-----------------------------------------------------------------------
! 1p1h - 2h2p (i|=j,a=b) block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
          call register1()
       end do
    end do

!-----------------------------------------------------------------------
! 1p1h - 2h2p (i|=j,a|b I) block
!-----------------------------------------------------------------------
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

!-----------------------------------------------------------------------
! 1p1h - 2h2p (i|=j,a|b II) block
!-----------------------------------------------------------------------
    dim_count=dim_count+kpq(5,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

    call register2()
    CLOSE(unt)
    write(6,*) count,' off-diagonal elements saved'
    
  contains

    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
         file_offdiag(count-buf_size*int(rec_count,8))= ar_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2_save_cvs

!#######################################################################

  subroutine get_phph_adc2(ndim,kpq,amatr)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: amatr

    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    integer :: i,j
    real(d) :: ar_diag,ar_offd

    amatr=0._d
    
!!$ Filling the ph-ph block: diagonal part

    do i=1, ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag=K_ph_ph(e(inda),e(indj))
       ar_diag=ar_diag+C1_ph_ph(inda,indj,inda,indj)
       ar_diag=ar_diag+CA_ph_ph(inda,inda)
       ar_diag=ar_diag+CB_ph_ph(indj,indj)
       ar_diag=ar_diag+CC_ph_ph(inda,indj,inda,indj)

       amatr(i,i)=ar_diag
    end do

!!$ Filling the ph-ph block : off-diagonal part

    do i=1,ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)      
          ar_offd=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
               ar_offd=ar_offd+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offd=ar_offd+CB_ph_ph(indj,indjpr)
          ar_offd=ar_offd+CC_ph_ph(inda,indj,indapr,indjpr)

          amatr(i,j)=ar_offd
       end do
    end do
    
  end subroutine get_phph_adc2
!!$-------------------------------------------
  
  subroutine get_ph_2p2h(ndim,i1,i2,kpq,bmx)

    integer, intent(in) :: i1,i2,ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,i2-i1+1), intent(out) :: bmx
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    integer :: i,j,nlim1,nlim2,mains

    mains=kpq(1,0)

    write(6,*) mains,i1,i2
   
    do i=i1,i2

       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)

       nlim1=kpq(1,0)+1
       nlim2=kpq(1,0)+kpq(2,0)

       do j=nlim1,nlim2
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          bmx(j-mains,i)=C5_ph_2p2h(inda,indj,indapr,indjpr)
       end do
    
       nlim1=nlim2+1
       nlim2=nlim2+kpq(3,0)
       do j=nlim1,nlim2
          write(6,*) j-mains,'jjjjjj'
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
          bmx(j-mains,i)=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
       end do

       nlim1=nlim2+1
       nlim2=nlim2+kpq(4,0)
       do j=nlim1,nlim2
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          bmx(j-mains,i)=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
       end do
       
       nlim1=nlim2+1
       nlim2=nlim2+kpq(5,0)
       do j=nlim1,nlim2
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          bmx(j-mains,i)=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
       end do
       
       nlim1=nlim2+1
       nlim2=nlim2+kpq(5,0)
       do j=nlim1,nlim2
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          bmx(j-mains,i)=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
       end do

    end do

  end subroutine get_ph_2p2h
!!$--------------------------------------------------
  subroutine get_2p2h2p2h_dg2s(ndim,kpq,ar_diag)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim), intent(out) :: ar_diag
    
    integer :: inda,indb,indj,indk,spin
    integer :: i,nlim
        
    nlim=kpq(1,0)

    do i=nlim+1, nlim+ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(i-nlim)=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
  end subroutine get_2p2h2p2h_dg2s

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************ADC2 extended BLOCK***************************
!!$*******************************************************************************
!!$*******************************************************************************

  subroutine get_diag_adc2ext_direct(ndim1,ndim2,kpq,ar_diag)
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1+ndim2), intent(out) :: ar_diag
    
    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,lim1,lim2
    
!!$ Filling the ph-ph block
    
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag(i)=K_ph_ph(ea,ej)
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
    
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
  end subroutine get_diag_adc2ext_direct

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

  subroutine get_diag_adc2ext_save(ndim1,ndim2,kpq,nbuf,chr)
  
    integer, intent(in) :: ndim1,ndim2,nbuf 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
   
    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    character(13) :: name
    integer :: i,ktype,dim_count,lim1,lim2,unt,a,b,c,d1
    real(d), dimension(ndim1+ndim2) :: ar_diag
     
    ktype=1
    name="hmlt.dia"//chr 
    unt=11
    
!!$ Filling the ph-ph block
    
    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag(i)=K_ph_ph(ea,ej)
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i=lim1, lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
     
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing",ndim1+ndim2," diagonal elements of ADC-ext. matrix in file ",name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diag(:))
!!$    call wrtdgat(unt,ndim1+ndim2,nbuf,ar_diag(:))
    CLOSE(unt)
   
    write(*,*) 'Writing successful at get_diag_adc2ext_save end'
  end subroutine get_diag_adc2ext_save

!!$---------------------------------------------------------------------------------
!!$---------------------------------------------------------------------------------

  subroutine get_offdiag_adc2ext_direct(ndim,kpq,ar_offdiag)

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real(d), dimension(ndim,ndim), intent(out) :: ar_offdiag
  
  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i, lim2i, lim1j, lim2j

  ar_offdiag(:,:)=0._d 

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i= 1,ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j= 1,i-1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
           ar_offdiag(i,j)=C1_ph_ph(inda,indj,indapr,indjpr)
           if(indj .eq. indjpr)&
                ar_offdiag(i,j)=ar_offdiag(i,j)+CA_ph_ph(inda,indapr)
           if(inda .eq. indapr)&
                ar_offdiag(i,j)=ar_offdiag(i,j)+CB_ph_ph(indj,indjpr)
           ar_offdiag(i,j)=ar_offdiag(i,j)+CC_ph_ph(inda,indj,indapr,indjpr)
        end do
     end do

     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             ar_offdiag(j,i)=C5_ph_2p2h(inda,indj,indapr,indjpr)
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(j,i)=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(j,i)=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(j,i)=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(j,i)=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do
    
!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=lim1i
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_1_1(inda,indj,indapr,indjpr)
           
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_2_1(inda,indb,indj,indapr,indjpr)
           
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_3_1(inda,indj,indk,indapr,indjpr)
           
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
           
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag(i,j)=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
           
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
           
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
            
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag(i,j)=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(i,j)=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
       end do
    end do


  end subroutine get_offdiag_adc2ext_direct

!!$-----------------------------------------------------------
!!$-----------------------------------------------------------

subroutine get_offdiag_adc2ext_save(ndim,kpq,nbuf,count,chr)
   
  integer, intent(in) :: ndim
  integer*8, intent(out) :: count
  integer, intent(out) :: nbuf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  character(1), intent(in) :: chr
  
  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  character(13) :: name
  integer :: rec_count
  integer :: i,j,nlim,dim_count,ndim1,unt
  integer :: lim1i, lim2i, lim1j, lim2j
  real(d) :: arr_offdiag_ij
  
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiag

  name="hmlt.off"//chr
  unt=12
  
  count=0
  rec_count=0
  
  write(6,*) "Writing the off-diagonal part of ADC matrix in file ", name
  OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
       FORM='UNFORMATTED')

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i=1,ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j=1,i-1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
           arr_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
           if(indj .eq. indjpr)&
                arr_offdiag_ij=arr_offdiag_ij+CA_ph_ph(inda,indapr)
           if(inda .eq. indapr)&
                arr_offdiag_ij=arr_offdiag_ij+CB_ph_ph(indj,indjpr)
           arr_offdiag_ij=arr_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
           call register1()
        end do
     end do
   
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             arr_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
             call register1()
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
             call register1()
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
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
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_1_1(inda,indj,indapr,indjpr)
           
          call register1()
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_2_1(inda,indb,indj,indapr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_1(inda,indj,indk,indapr,indjpr)
           
          call register1()
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
           
          call register1()
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
            
          call register1()
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
            
          call register1()
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do
    
    call register2()
    CLOSE(unt)
    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' off-diagonal elements saved in file ', name

  contains
       
    subroutine register1()
      if (abs(arr_offdiag_ij) .gt. minc) then
         count=count+1
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiag(count-buf_size*int(rec_count,8))=arr_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
!!$            call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)  
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
!!$      call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2ext_save
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!-------------------FANO SUBROUTINES----------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

!!$variable meth stands for method calling the subroutine: 1-tda,2-adc2,21-adc2e

  subroutine adc2ext_0_0(meth,a,j,a1,j1,matel)
        
    integer, intent(in) :: meth,a,j,a1,j1
    real(d), intent(out) :: matel

    real(d) :: ea,ej

    matel=0._d

    if((a .eq. a1) .and. (j .eq. j1)) then
       ea=e(a)
       ej=e(j)
       matel=K_ph_ph(ea,ej)
    end if
    
    matel=matel+C1_ph_ph(a,j,a1,j1)


    if((meth .eq. 2) .or. (meth .eq. 21)) then
       if(j .eq. j1)&
            matel=matel+CA_ph_ph(a,a1)
       if(a .eq. a1)&
            matel=matel+CB_ph_ph(j,j1)
       matel=matel+CC_ph_ph(a,j,a1,j1)
    end if
     
  end subroutine adc2ext_0_0
!!$---------------------------------------------------  
  subroutine adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
    
    integer, intent(in):: meth,a,j,a1,j1
    real(d), intent(in):: ea,ej
    real(d), intent(out) :: matel
    
    logical :: diag
    
    matel=0._d
    diag=(a .eq. a1) .and. (j .eq. j1)
    
    if(diag) then
       matel=K_2p2h_2p2h(ea,ea,ej,ej)
    end if
    if(meth .eq. 21)&
         matel=matel+C_1_1(a,j,a1,j1)
    
  end subroutine adc2ext_1_1
!!$------------------------------------------------------
  subroutine adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
    
    integer, intent(in):: meth,a,b,j,a1,b1,j1
    real(d), intent(in):: ea,eb,ej
    real(d), intent(out) :: matel
    
    logical :: diag
    
    matel=0._d
    diag=(a .eq. a1) .and. (b .eq. b1) .and. (j .eq. j1)
    
    if(diag) then
       matel=K_2p2h_2p2h(ea,eb,ej,ej)
    end if
    
    if(meth .eq. 21)&
         matel=matel+C_2_2(a,b,j,a1,b1,j1)
    
  end subroutine adc2ext_2_2
!!$------------------------------------------------------
  subroutine adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
    
    integer, intent(in):: meth,a,j,k,a1,j1,k1
    real(d), intent(in):: ea,ej,ek
    real(d), intent(out) :: matel
    
    logical :: diag
    
    matel=0._d
    diag=(a .eq. a1) .and. (j .eq. j1) .and. (k .eq. k1)
    
    if(diag) then
       matel=K_2p2h_2p2h(ea,ea,ej,ek)
    end if
    
    if(meth .eq. 21)&
         matel=matel+C_3_3(a,j,k,a1,j1,k1)
    
  end subroutine adc2ext_3_3
!!$--------------------------------------------------------------
  subroutine adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
    
    integer, intent(in):: meth,a,b,j,k,a1,b1,j1,k1
    real(d), intent(in):: ea,eb,ej,ek
    real(d), intent(out) :: matel
    
    logical :: diag
    
    matel=0._d
    diag=(a .eq. a1) .and. (b .eq. b1) .and. (j .eq. j1) .and. (k .eq. k1)
    
    if(diag) then
       matel=K_2p2h_2p2h(ea,eb,ej,ek)
    end if
    
    if(meth .eq. 21)&
         matel=matel+C_4i_4i(a,b,j,k,a1,b1,j1,k1)
    
  end subroutine adc2ext_4i_4i
!!$------------------------------------------------------------------
  subroutine adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
    
    integer, intent(in):: meth,a,b,j,k,a1,b1,j1,k1
    real(d), intent(in):: ea,eb,ej,ek
    real(d), intent(out) :: matel
    
    logical :: diag
    
    matel=0._d
    diag=(a .eq. a1) .and. (b .eq. b1) .and. (j .eq. j1) .and. (k .eq. k1)
    
    if(diag) then
       matel=K_2p2h_2p2h(ea,eb,ej,ek)
    end if
    
    if(meth .eq. 21)&
         matel=matel+C_4ii_4ii(a,b,j,k,a1,b1,j1,k1)
    
  end subroutine adc2ext_4ii_4ii
   





  subroutine get_offdiag_adc2_save_MIO(ndim,kpq,nbuf,count,indx,chr)

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf
    integer*8, intent(out) :: count 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    INTEGER, DIMENSION(ndim), intent(in) :: indx  
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    character(10) :: name
    integer :: i,j,nlim,rec_count,dim_count,ndim1,unt
    real(d) :: ar_offdiag_ij
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    
    name="hmlt.off"//chr
    unt=12

    count=0
    rec_count=0
    
    
    write(6,*) "Writing the off-diagonal part of ADC matrix in file ", name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

!!$ Full diagonalization.  

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
          ar_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
               ar_offdiag_ij= ar_offdiag_ij+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offdiag_ij= ar_offdiag_ij+CB_ph_ph(indj,indjpr)
          ar_offdiag_ij= ar_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
          call register1()
       end do
    end do
    
       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
          ar_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
          call register1()
       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
          call register1()
       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
          call register1()
       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          call register1()
       end do
    end do

    call register2()
    CLOSE(unt)
    write(6,*) count,' off-diagonal elements saved'

  contains
    
    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
         file_offdiag(count-buf_size*int(rec_count,8))= ar_offdiag_ij
         if ( indx(i) .ge. indx(j) ) then 
         oi(count-buf_size*int(rec_count,8))=indx(i)
         oj(count-buf_size*int(rec_count,8))=indx(j)
         else
         oi(count-buf_size*int(rec_count,8))=indx(j)
         oj(count-buf_size*int(rec_count,8))=indx(i)
         end if
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_adc2_save_MIO
!!$-----------------------------------





  subroutine get_diag_adc2_save_MIO(ndim1,ndim2,kpq,nbuf,indx,chr)
  
    integer, intent(in) :: ndim1,ndim2,nbuf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    INTEGER, DIMENSION(ndim1+ndim2), intent(in) :: indx  
 
    integer :: inda,indb,indj,indk,spin
    
    character(10) :: name
    integer :: i,ktype,unt 
    real(d), dimension(:), allocatable:: ar_diag

    allocate(ar_diag(ndim1+ndim2))


    ktype=1
    name="hmlt.dia"//chr 
    unt=11

!!$ Filling the ph-ph block

    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(indx(i))=K_ph_ph(e(inda),e(indj))
       ar_diag(indx(i))=ar_diag(indx(i))+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CA_ph_ph(inda,inda)
       ar_diag(indx(i))=ar_diag(indx(i))+CB_ph_ph(indj,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CC_ph_ph(inda,indj,inda,indj)
    end do
    
!!$ Filling the 2p2h-2p2h block
    
    do i=ndim1+1, ndim1+ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(indx(i))=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing the diagonal part of ADC matrix in file ", name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diag(:))
    CLOSE(unt)

    deallocate(ar_diag)
  end subroutine get_diag_adc2_save_MIO

!!$----------------------------------------------------------------------------
!!$----------------------------------------------------------------------------




  subroutine get_diag_adc2_direct_MIO(ndim1,ndim2,kpq,ar_diag,indx)
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1+ndim2), intent(out) :: ar_diag
    INTEGER, dimension(ndim1+ndim2), intent(in) :: indx
    

    integer :: inda,indb,indj,indk,spin
    integer :: i
    
!!$ Filling the ph-ph block

    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(indx(i))=K_ph_ph(e(inda),e(indj))
       ar_diag(indx(i))=ar_diag(indx(i))+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CA_ph_ph(inda,inda)
       ar_diag(indx(i))=ar_diag(indx(i))+CB_ph_ph(indj,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CC_ph_ph(inda,indj,inda,indj)
    end do
    
!!$ Filling the 2p2h-2p2h block
    
    do i=ndim1+1, ndim1+ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(indx(i))=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
  end subroutine get_diag_adc2_direct_MIO




  subroutine get_offdiag_adc2_direct_MIO(ndim,kpq,ar_offdiag,indx)
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiag
    INTEGER, dimension(ndim), intent(in) :: indx
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    integer :: i,j,dim_count,ndim1
    
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    
    ar_offdiag(:,:)=0._d
    
!!$ Full diagonalization. 

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
          ar_offdiag(indx(i),indx(j))=C1_ph_ph(inda,indj,indapr,indjpr)
          if(indj .eq. indjpr)&
                  ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CA_ph_ph(inda,indapr)
          if(inda .eq. indapr)&
               ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CB_ph_ph(indj,indjpr)
          ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CC_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do

       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
          ar_offdiag(indx(i),indx(j))=C5_ph_2p2h(inda,indj,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(indx(i),indx(j))=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(indx(i),indx(j))=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(indx(i),indx(j))=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
          ar_offdiag(indx(i),indx(j))=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
       end do
    end do
  
  end subroutine get_offdiag_adc2_direct_MIO







  subroutine get_offdiag_adc2ext_direct_MIO(ndim,kpq,ar_offdiag,indx)

  integer, intent(in) :: ndim
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  real(d), dimension(ndim,ndim), intent(out) :: ar_offdiag
  INTEGER, DIMENSION(ndim), INTENT(IN) :: indx  
  
  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i, lim2i, lim1j, lim2j

  ar_offdiag(:,:)=0._d 

!!$ Full diagonalization. Filling the lower half of the matrix

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i= 1,ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j= 1,ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
           ar_offdiag(indx(i),indx(j))=C1_ph_ph(inda,indj,indapr,indjpr)
           if(indj .eq. indjpr)&
                ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CA_ph_ph(inda,indapr)
           if(inda .eq. indapr)&
                ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CB_ph_ph(indj,indjpr)
           ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(i),indx(j))+CC_ph_ph(inda,indj,indapr,indjpr)
        end do
     end do

     
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             ar_offdiag(indx(j),indx(i))=C5_ph_2p2h(inda,indj,indapr,indjpr)
             ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(j),indx(i))
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(indx(j),indx(i))=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
             ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(j),indx(i))
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(indx(j),indx(i))=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
             ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(j),indx(i))
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(indx(j),indx(i))=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(j),indx(i))
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i= 1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j= dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             ar_offdiag(indx(j),indx(i))=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             ar_offdiag(indx(i),indx(j))=ar_offdiag(indx(j),indx(i))
!!$             ar_offdiag(i,j)=ar_offdiag(j,i)
          end do
       end do
    
!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=lim1i
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_1_1(inda,indj,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_2_1(inda,indb,indj,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_3_1(inda,indj,indk,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag(indx(i),indx(j))=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
            
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    
    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag(indx(i),indx(j))=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i= lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j= lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          ar_offdiag(indx(i),indx(j))=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
          ar_offdiag(indx(j),indx(i))=ar_offdiag(indx(i),indx(j))
           
       end do
    end do


  end subroutine get_offdiag_adc2ext_direct_MIO

!!$-----------------------------------------------------------
!!$-----------------------------------------------------------

subroutine get_offdiag_adc2ext_save_MIO(ndim,kpq,nbuf,count,indx,chr)
   
  integer, intent(in) :: ndim
  integer*8, intent(out) :: count
  integer, intent(out) :: nbuf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  character(1), intent(in) :: chr
  INTEGER, DIMENSION(ndim), INTENT(IN) :: indx  

  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  character(13) :: name
  integer :: rec_count
  integer :: i,j,nlim,dim_count,ndim1,unt
  integer :: lim1i, lim2i, lim1j, lim2j
  real(d) :: arr_offdiag_ij
  
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiag

  name="hmlt.off"//chr
  unt=12
  
  count=0
  rec_count=0
  
  write(6,*) "Writing the off-diagonal part of ADC matrix in file ", name
  OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
       FORM='UNFORMATTED')

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i=1,ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j=1,i-1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
           arr_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
           if(indj .eq. indjpr)&
                arr_offdiag_ij=arr_offdiag_ij+CA_ph_ph(inda,indapr)
           if(inda .eq. indapr)&
                arr_offdiag_ij=arr_offdiag_ij+CB_ph_ph(indj,indjpr)
           arr_offdiag_ij=arr_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
           call register1()
        end do
     end do
   
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             arr_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
             call register1()
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
             call register1()
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
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
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_1_1(inda,indj,indapr,indjpr)
           
          call register1()
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_2_1(inda,indb,indj,indapr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_1(inda,indj,indk,indapr,indjpr)
           
          call register1()
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
           
          call register1()
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
            
          call register1()
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
            
          call register1()
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
          call register1()
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
          call register1()
       end do
    end do
    
    call register2()
    CLOSE(unt)
    write(6,*) 'rec_counts',nbuf
    write(6,*) count,' off-diagonal elements saved in file ', name

  contains
       
    subroutine register1()
      if (abs(arr_offdiag_ij) .gt. minc) then
         count=count+1
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiag(count-buf_size*int(rec_count,8))=arr_offdiag_ij
         if ( indx(i) .ge. indx(j) ) then 
         oi(count-buf_size*int(rec_count,8))=indx(i)
         oj(count-buf_size*int(rec_count,8))=indx(j)
         else
         oi(count-buf_size*int(rec_count,8))=indx(j)
         oj(count-buf_size*int(rec_count,8))=indx(i)
         end if

         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
!!$            call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)  
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
!!$      call wrtoffat(unt,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2ext_save_MIO


  subroutine get_diag_adc2ext_direct_MIO(ndim1,ndim2,kpq,ar_diag,indx)
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1+ndim2), intent(out) :: ar_diag
    INTEGER, DIMENSION(ndim1+ndim2), INTENT(IN) :: indx    


    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,lim1,lim2
    
!!$ Filling the ph-ph block
    
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag(indx(i))=K_ph_ph(ea,ej)
       ar_diag(indx(i))=ar_diag(indx(i))+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CA_ph_ph(inda,inda)
       ar_diag(indx(i))=ar_diag(indx(i))+CB_ph_ph(indj,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
    
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i= lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
  end subroutine get_diag_adc2ext_direct_MIO

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

  subroutine get_diag_adc2ext_save_MIO(ndim1,ndim2,kpq,nbuf,indx,chr)
  
    integer, intent(in) :: ndim1,ndim2,nbuf 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    INTEGER, DIMENSION(ndim1+ndim2), INTENT(IN) :: indx    
   
    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    character(13) :: name
    integer :: i,ktype,dim_count,lim1,lim2,unt,a,b,c,d1
    real(d), dimension(ndim1+ndim2) :: ar_diag
     
    ktype=1
    name="hmlt.dia"//chr 
    unt=11
    
!!$ Filling the ph-ph block
    
    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag(indx(i))=K_ph_ph(ea,ej)
       ar_diag(indx(i))=ar_diag(indx(i))+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CA_ph_ph(inda,inda)
       ar_diag(indx(i))=ar_diag(indx(i))+CB_ph_ph(indj,indj)
       ar_diag(indx(i))=ar_diag(indx(i))+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i=lim1, lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
     
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i=lim1,lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(indx(i))=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(indx(i))=ar_diag(indx(i))+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing",ndim1+ndim2," diagonal elements of ADC-ext. matrix in file ",name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diag(:))
!!$    call wrtdgat(unt,ndim1+ndim2,nbuf,ar_diag(:))
    CLOSE(unt)
   
    write(*,*) 'Writing successful at get_diag_adc2ext_save end'
  end subroutine get_diag_adc2ext_save_MIO









  subroutine get_diag_tda_save_OK(ndim,kpq, UNIT_HAM )
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM
    
    integer :: i,ktype
    integer :: inda,indb,indj,indk,spin
    real(d), dimension( ndim ) :: ar_diag
    
    ar_diag(:) = 0.d0    

    do i = 1 , ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i) = K_ph_ph(e(inda),e(indj))
       ar_diag(i) = ar_diag(i) + C1_ph_ph(inda,indj,inda,indj)
    end do

       !Saving in file

       write(6,*) "Writing the diagonal part in file ", UNIT_HAM

       write( UNIT_HAM ) buf_size
       write( UNIT_HAM ) ar_diag

       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim )
       
  end subroutine get_diag_tda_save_OK

!!$----------------------------------------------------------------
!!$----------------------------------------------------------------

  subroutine get_offdiag_tda_save_OK(ndim,kpq,nbuf,count, UNIT_HAM )
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf 
    INTEGER*8, intent(out) :: count
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM
    
    integer :: i,j,nlim,rec_count
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr
    real(d) :: ar_offdiag_ij, ar_offdiag_ji

    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag

    

    count=0
    rec_count=0

       write(6,*) "Writing the off-diagonal part of TDA matrix in file ", UNIT_HAM
    
    do i = 1 , ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = i + 1 , ndim
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag_ji=C1_ph_ph(indapr,indjpr,inda,indj)
          if(abs(ar_offdiag_ij-ar_offdiag_ji) .ge. 1.e-15_d) then
             write(6,*) "TDA matrix is not symmetric. Stopping now."
             stop
          end if

!!$ Saving into vector for the following Lanzcos/Davidson routine 
            

          !Culling  small matrix elements
          if (abs(ar_offdiag_ij) .gt. minc) then
             call register1()
          end if

       end do
    end do
!!$

       call register2()
       
  contains
       
    subroutine register1()
      
      count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i , j,'one:', ar_offdiag_ij
         END IF
      file_offdiag(count-buf_size*rec_count)=ar_offdiag_ij
      oi(count-buf_size*rec_count)=i
      oj(count-buf_size*rec_count)=j
      !Checking if the buffer is full 
      if(mod(count,buf_size) .eq. 0) then
         rec_count=rec_count+1
         nlim=buf_size
         !Saving off-diag part in file
         call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
      end if

    end subroutine register1
       
    subroutine register2()
         
      !Saving the rest of matrix in file
      nlim=count-buf_size*rec_count
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_tda_save_OK




!!$---------------------------------------------------
!!$---------------------------------------------------

  subroutine get_diag_tda_save_GS( ndim , kpq , UNIT_HAM )
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM
    
    integer :: i,ktype
    integer :: inda,indb,indj,indk,spin
    real(d), dimension( ndim + 1 ) :: ar_diag
    
    ar_diag(:) = 0.d0    

    do i = 1 , ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag( i + 1 ) = K_ph_ph(e(inda),e(indj))
       ar_diag( i + 1 ) = ar_diag( i + 1 ) + C1_ph_ph(inda,indj,inda,indj)
    end do

       !Saving in file

       write(6,*) "Writing the diagonal part in file ", UNIT_HAM

       write( UNIT_HAM ) buf_size
       write( UNIT_HAM )   ar_diag    
       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim + 1 )


  end subroutine get_diag_tda_save_GS



!!$----------------------------------------------------------------
!!$----------------------------------------------------------------



  subroutine get_offdiag_tda_save_GS(ndim,kpq,nbuf,count,UNIT_HAM)
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf 
    integer*8, intent(out) :: count
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM    

    integer :: i,j,nlim,rec_count
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr
    real(d) :: ar_offdiag_ij, ar_offdiag_ji

    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag

    

    count=0
    rec_count=0

       write(6,*) "Writing the off-diagonal part of TDA matrix in file ", UNIT_HAM
   
   i = 0
   do j = 1 , ndim
      ar_offdiag_ij= 0.d0
          !Culling  small matrix elements
          if (abs(ar_offdiag_ij) .gt. minc) then
             call register1()
          end if
   end do


    do i = 1 , ndim
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = i + 1 , ndim
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          ar_offdiag_ij = 0.d0 
          ar_offdiag_ij = C1_ph_ph(inda,indj,indapr,indjpr)
          ar_offdiag_ji = C1_ph_ph(indapr,indjpr,inda,indj)
          if(abs(ar_offdiag_ij-ar_offdiag_ji) .ge. 1.e-15_d) then
             write(6,*) "TDA matrix is not symmetric. Stopping now."
             stop
          end if

!!$ Saving into vector for the following Lanzcos/Davidson routine 
            

          !Culling  small matrix elements
          if (abs(ar_offdiag_ij) .gt. minc) then
             call register1()
          end if

       end do
    end do
!!$

       call register2()
       
  contains
       
    subroutine register1()
      
      count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i+1 , j+1,'one:', ar_offdiag_ij
         END IF
      file_offdiag(count-buf_size*rec_count)=ar_offdiag_ij
      oi(count-buf_size*rec_count) = i + 1
      oj(count-buf_size*rec_count) = j + 1
      !Checking if the buffer is full 
      if(mod(count,buf_size) .eq. 0) then
         rec_count=rec_count+1
         nlim=buf_size
         !Saving off-diag part in file
         call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
      end if

    end subroutine register1
       
    subroutine register2()
         
      !Saving the rest of matrix in file
      nlim=count-buf_size*rec_count
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_tda_save_GS




  subroutine get_offdiag_adc2_save_OK(ndim,kpq,nbuf,count, UNIT_HAM )

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf
    integer*8, intent(out) :: count 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    integer :: i,j,nlim,rec_count,dim_count,ndim1
    real(d) :: ar_offdiag_ij
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    

    count=0
    rec_count=0
    
    
    write(6,*) "Writing the off-diagonal part of ADC matrix in file ", UNIT_HAM

!!$ Full diagonalization.  

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = i + 1 , ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             

          ar_offdiag_ij = 0.d0
 
          ar_offdiag_ij = C1_ph_ph(inda,indj,indapr,indjpr)

          if(indj .eq. indjpr)&
               ar_offdiag_ij = ar_offdiag_ij + CA_ph_ph(inda,indapr)

          if(inda .eq. indapr)&
               ar_offdiag_ij = ar_offdiag_ij + CB_ph_ph(indj,indjpr)

          ar_offdiag_ij = ar_offdiag_ij + CC_ph_ph(inda,indj,indapr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do

    call register2()

    write(6,*) count,' off-diagonal elements saved in file', UNIT_HAM

  contains
    
    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i , j,'one:', ar_offdiag_ij
         END IF
         file_offdiag(count-buf_size*int(rec_count,8))= ar_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_adc2_save_OK
!!$-----------------------------------


  subroutine get_offdiag_adc2_save_GS(ndim,kpq,nbuf,count, UNIT_HAM )

!!$The difference from the earlier routine is that this routine returns the total number of saved els to a caller. 
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nbuf
    integer*8, intent(out) :: count 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    integer :: i,j,nlim,rec_count,dim_count,ndim1
    real(d) :: ar_offdiag_ij
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiag
    

    count=0
    rec_count=0
    
    
    write(6,*) "Writing the off-diagonal part of ADC matrix in file ", UNIT_HAM

!!$ Full diagonalization.  

!!$ Filling the off-diagonal part of the ph-ph block


    i = 0 
    do j = 1 , ndim
          ar_offdiag_ij = 0.d0
      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if
    end do


    ndim1=kpq(1,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = i + 1 , ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             

          ar_offdiag_ij = 0.d0
 
          ar_offdiag_ij = C1_ph_ph(inda,indj,indapr,indjpr)

          if(indj .eq. indjpr)&
               ar_offdiag_ij = ar_offdiag_ij + CA_ph_ph(inda,indapr)

          if(inda .eq. indapr)&
               ar_offdiag_ij = ar_offdiag_ij + CB_ph_ph(indj,indjpr)

          ar_offdiag_ij = ar_offdiag_ij + CC_ph_ph(inda,indj,indapr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
       
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

    dim_count = kpq(1,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
!!$ Coupling to the i=j,a|=b configs   
    
    dim_count = dim_count + kpq(2,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
    
!!$ Coupling to the i|=j,a=b configs
    
    dim_count = dim_count + kpq(3,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do
       
!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count = dim_count + kpq(4,0)
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count = dim_count + kpq(5,0)

    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j = dim_count + 1 , dim_count + kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          ar_offdiag_ij = 0.d0
          ar_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)

      if (abs(ar_offdiag_ij) .gt. minc) then
          call register1()
      end if

       end do
    end do

    call register2()

    write(6,*) count,' off-diagonal elements saved in file', UNIT_HAM

  contains
    
    subroutine register1()
      if (abs(ar_offdiag_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i+1 , j+1,'one:', ar_offdiag_ij
         END IF
         file_offdiag(count-buf_size*int(rec_count,8))= ar_offdiag_ij
         oi(count-buf_size*int(rec_count,8)) = i + 1
         oj(count-buf_size*int(rec_count,8)) = j + 1
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim) 
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2
    
  end subroutine get_offdiag_adc2_save_GS
!!$-----------------------------------






  subroutine get_diag_adc2_save_OK(ndim1,ndim2,kpq, UNIT_HAM )
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM

    integer :: inda,indb,indj,indk,spin
    
    integer :: i,ktype
    real(d), dimension(:), allocatable:: ar_diag

    allocate(ar_diag( ndim1 + ndim2 ))



!!$ Filling the ph-ph block

    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag(i)=K_ph_ph(e(inda),e(indj))
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do
    
!!$ Filling the 2p2h-2p2h block
    
    do i = ndim1 + 1 , ndim1 + ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag(i)=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
    !Saving the diagonal part in file

    write(6,*) "Writing the diagonal part of ADC matrix in file ", UNIT_HAM

    write( UNIT_HAM ) buf_size
    write( UNIT_HAM ) ar_diag
       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim1 + ndim2 )



    deallocate(ar_diag)


  end subroutine get_diag_adc2_save_OK


  subroutine get_diag_adc2_save_GS(ndim1,ndim2,kpq, UNIT_HAM )
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: UNIT_HAM

    integer :: inda,indb,indj,indk,spin
    
    integer :: i,ktype
    real(d), dimension(:), allocatable:: ar_diag

    allocate(ar_diag( ndim1 + ndim2 + 1 ))


    ar_diag(:) = 0.d0

!!$ Filling the ph-ph block

    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ar_diag( i + 1 )=K_ph_ph(e(inda),e(indj))
       ar_diag( i + 1 )=ar_diag( i + 1 )+C1_ph_ph(inda,indj,inda,indj)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CA_ph_ph(inda,inda)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CB_ph_ph(indj,indj)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CC_ph_ph(inda,indj,inda,indj)
    end do
    
!!$ Filling the 2p2h-2p2h block
    
    do i = ndim1 + 1 , ndim1 + ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ar_diag( i + 1 )=K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk))
    end do
    
    !Saving the diagonal part in file

    write(6,*) "Writing the diagonal part of ADC matrix in file ", UNIT_HAM

    write( UNIT_HAM ) buf_size
    write( UNIT_HAM ) ar_diag
       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim1 + ndim2 + 1 )



    deallocate(ar_diag)


  end subroutine get_diag_adc2_save_GS





  subroutine get_diag_adc2ext_save_OK(ndim1,ndim2,kpq, UNIT_HAM )
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    INTEGER, intent(in) :: UNIT_HAM
   
    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,ktype,dim_count,lim1,lim2,a,b,c,d1

    real(d), dimension( ndim1 + ndim2 ) :: ar_diag
     
    
!!$ Filling the ph-ph block
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag(i)=K_ph_ph(ea,ej)
       ar_diag(i)=ar_diag(i)+C1_ph_ph(inda,indj,inda,indj)
       ar_diag(i)=ar_diag(i)+CA_ph_ph(inda,inda)
       ar_diag(i)=ar_diag(i)+CB_ph_ph(indj,indj)
       ar_diag(i)=ar_diag(i)+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
     
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag(i)=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag(i)=ar_diag(i)+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing", ndim1 + ndim2 ," diagonal elements of ADC-ext. matrix in file ", UNIT_HAM

    write( UNIT_HAM ) buf_size
    write( UNIT_HAM ) ar_diag
       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim1 + ndim2 )
   
    write(*,*) 'Writing successful at get_diag_adc2ext_save end'

  end subroutine get_diag_adc2ext_save_OK




  subroutine get_diag_adc2ext_save_GS(ndim1,ndim2,kpq, UNIT_HAM )
  
    integer, intent(in) :: ndim1,ndim2
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    INTEGER, intent(in) :: UNIT_HAM
   
    integer :: inda,indb,indj,indk,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,ktype,dim_count,lim1,lim2,a,b,c,d1

    real(d), dimension( ndim1 + ndim2 + 1 ) :: ar_diag
     
    
    ar_diag(:) = 0.d0

!!$ Filling the ph-ph block
    
    do i = 1 , ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       ea=e(inda)
       ej=e(indj)
       ar_diag( i + 1 )=K_ph_ph(ea,ej)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C1_ph_ph(inda,indj,inda,indj)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CA_ph_ph(inda,inda)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CB_ph_ph(indj,indj)
       ar_diag( i + 1 )=ar_diag( i + 1 )+CC_ph_ph(inda,indj,inda,indj)
    end do

!!$ Filling the 2p2h-2p2h block
!!$ Filling (1,1) block
    
    lim1=ndim1+1
    lim2=ndim1+kpq(2,0)
    
    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag( i + 1 )=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C_1_1(inda,indj,inda,indj)
    end do

!!$ Filling (2,2) block
    
    lim1=lim1+kpq(2,0)
    lim2=lim2+kpq(3,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag( i + 1 )=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C_2_2(inda,indb,indj,inda,indb,indj)
    end do
    
!!$ Filling (3,3) block
    
    lim1=lim1+kpq(3,0)
    lim2=lim2+kpq(4,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag( i + 1 )=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C_3_3(inda,indj,indk,inda,indj,indk)
    end do
     
!!$ Filling (4i,4i) block  
    
    lim1=lim1+kpq(4,0)
    lim2=lim2+kpq(5,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag( i + 1 )=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C_4i_4i(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
!!$ Filling (4ii,4ii) block  
    
    lim1=lim1+kpq(5,0)
    lim2=lim2+kpq(5,0)

    do i = lim1 , lim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       ea=e(inda)
       eb=e(indb)
       ej=e(indj)
       ek=e(indk)
       ar_diag( i + 1 )=K_2p2h_2p2h(ea,eb,ej,ek)
       ar_diag( i + 1 )=ar_diag( i + 1 )+C_4ii_4ii(inda,indb,indj,indk,inda,indb,indj,indk)
    end do
    
    !Saving the diagonal part in file
    write(6,*) "Writing", ndim1 + ndim2 + 1 ," diagonal elements of ADC-ext. matrix in file ", UNIT_HAM

    write( UNIT_HAM ) buf_size
    write( UNIT_HAM ) ar_diag
       write(*,*) 'the first element diagonal saved in', UNIT_HAM,'is:', ar_diag(1)
       write(*,*) 'the last  element diagonal saved in', UNIT_HAM,'is:', ar_diag( ndim1 + ndim2 + 1 ) 
   
    write(*,*) 'Writing successful at get_diag_adc2ext_save end'

  end subroutine get_diag_adc2ext_save_GS






subroutine get_offdiag_adc2ext_save_OK(ndim,kpq,nbuf,count, UNIT_HAM )
   
  integer, intent(in) :: ndim
  integer*8, intent(out) :: count
  integer, intent(out) :: nbuf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  INTEGER, intent(in) :: UNIT_HAM
  
  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  integer :: rec_count
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i, lim2i, lim1j, lim2j
  real(d) :: arr_offdiag_ij
  
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiag

  
  count=0
  rec_count=0
  
  write(6,*) "Writing the off-diagonal part of ADC matrix in file ", UNIT_HAM

!!$ Filling the off-diagonal part of the ph-ph block

     ndim1=kpq(1,0)
       
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             
           arr_offdiag_ij=C1_ph_ph(inda,indj,indapr,indjpr)
           if(indj .eq. indjpr)&
                arr_offdiag_ij=arr_offdiag_ij+CA_ph_ph(inda,indapr)
           if(inda .eq. indapr)&
                arr_offdiag_ij=arr_offdiag_ij+CB_ph_ph(indj,indjpr)
           arr_offdiag_ij=arr_offdiag_ij+CC_ph_ph(inda,indj,indapr,indjpr)
      if (abs(arr_offdiag_ij) .gt. minc) then
           call register1()
      end if
        end do
     end do
   
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             arr_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do
    
!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_1_1(inda,indj,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_2_1(inda,indb,indj,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_1(inda,indj,indk,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
            
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
            
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
    
    call register2()

    write(6,*) 'rec_counts' , nbuf
    write(6,*) count,' off-diagonal elements saved in file ', UNIT_HAM

  contains
       
    subroutine register1()
      if (abs(arr_offdiag_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i , j,'one:', arr_offdiag_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiag(count-buf_size*int(rec_count,8))=arr_offdiag_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2ext_save_OK





subroutine get_offdiag_adc2ext_save_GS(ndim,kpq,nbuf,count, UNIT_HAM )
   
  integer, intent(in) :: ndim
  integer*8, intent(out) :: count
  integer, intent(out) :: nbuf
  integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
  INTEGER, intent(in) :: UNIT_HAM
  
  integer :: inda,indb,indj,indk,spin
  integer :: indapr,indbpr,indjpr,indkpr,spinpr 
  
  integer :: rec_count
  integer :: i,j,nlim,dim_count,ndim1
  integer :: lim1i, lim2i, lim1j, lim2j
  real(d) :: arr_offdiag_ij
  
  integer, dimension(buf_size) :: oi,oj
  real(d), dimension(buf_size) :: file_offdiag

  
  count=0
  rec_count=0
  
  write(6,*) "Writing the off-diagonal part of ADC matrix in file ", UNIT_HAM

!!$ Filling the off-diagonal part of the ph-ph block


     i = 0
     do j = 1 , ndim
     arr_offdiag_ij = 0.d0
      if (abs(arr_offdiag_ij) .gt. minc) then
           call register1()
      end if
     end do



     ndim1=kpq(1,0)
       
     do i = 1 , ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j = i + 1 , ndim1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)             

           arr_offdiag_ij = 0.d0
           arr_offdiag_ij = C1_ph_ph(inda,indj,indapr,indjpr)

           if(indj .eq. indjpr)&
                arr_offdiag_ij = arr_offdiag_ij + CA_ph_ph(inda,indapr)

           if(inda .eq. indapr)&
                arr_offdiag_ij = arr_offdiag_ij + CB_ph_ph(indj,indjpr)

           arr_offdiag_ij = arr_offdiag_ij + CC_ph_ph(inda,indj,indapr,indjpr)

      if (abs(arr_offdiag_ij) .gt. minc) then
           call register1()
      end if
        end do
     end do
   
!!$ Filling the off-diagonal part of the ph-2p2h block 
!!$ Coupling to the i=j,a=b configs

       dim_count=kpq(1,0)
       
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(2,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)    
             arr_offdiag_ij=C5_ph_2p2h(inda,indj,indapr,indjpr)
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do
          
!!$ Coupling to the i=j,a|=b configs   
       
       dim_count=dim_count+kpq(2,0)
       
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(3,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a=b configs
       
       dim_count=dim_count+kpq(3,0)
             
       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr)
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b I configs
       
       dim_count=dim_count+kpq(4,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
             if (abs(arr_offdiag_ij) .gt. minc) then
                call register1()
             end if
          end do
       end do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

       do i = 1 , ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  
             arr_offdiag_ij=C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr)
             !Culling  small matrix elements
      if (abs(arr_offdiag_ij) .gt. minc) then
             call register1()
      end if
          end do
       end do
    
!!$ Filling the 2p2h-2p2h block
    
!!$ (1,1) block
    
    lim1i=kpq(1,0)+1
    lim2i=kpq(1,0)+kpq(2,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_1_1(inda,indj,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (2,1) block 

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_2_1(inda,indb,indj,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,1) block
     
    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_1(inda,indj,indk,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do          
         
!!$ (4i,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_1(inda,indb,indj,indk,indapr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do 
 
!!$ (4ii,1) block

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+1
    lim2j=kpq(1,0)+kpq(2,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_1(inda,indb,indj,indk,indapr,indjpr)
            
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do 

!!$ (2,2) block

    lim1i=kpq(1,0)+kpq(2,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_2_2(inda,indb,indj,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_2(inda,indj,indk,indapr,indbpr,indjpr)
            
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
        
!!$ (4i,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,2) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_2(inda,indb,indj,indk,indapr,indbpr,indjpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (3,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_3_3(inda,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4i,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,3) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_3(inda,indb,indj,indk,indapr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4i,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    lim1j=lim1i
    
    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4i_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do

!!$ (4ii,4i) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1
    lim2j=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,lim2j
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)
          arr_offdiag_ij=C_4ii_4i(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
    
!!$ (4ii,4ii) block 

    lim1i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+1
    lim2i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)+kpq(5,0)
    lim1j=lim1i

    do i=lim1i,lim2i
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=lim1j,i-1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 
          arr_offdiag_ij=C_4ii_4ii(inda,indb,indj,indk,indapr,indbpr,indjpr,indkpr)
           
      if (abs(arr_offdiag_ij) .gt. minc) then
          call register1()
      end if
       end do
    end do
    
    call register2()

    write(6,*) 'rec_counts' , nbuf
    write(6,*) count,' off-diagonal elements saved in file ', UNIT_HAM

  contains
       
    subroutine register1()
      if (abs(arr_offdiag_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(*,*) 'the first element not-diagonal saved in', UNIT_HAM,'is the', i+1 , j+1,'one:', arr_offdiag_ij
         END IF
! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiag(count-buf_size*int(rec_count,8))=arr_offdiag_ij
         oi(count-buf_size*int(rec_count,8)) = i + 1
         oj(count-buf_size*int(rec_count,8)) = j + 1
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
         end if
      end if
    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg( UNIT_HAM ,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_adc2ext_save_GS








 
end module get_matrix
          
          
       
       
       
       
    
    
    
       
       
       
       
       
       
    
    
    
    

    
    
    
    
