module mxv

  use constants
  use parameters
  use adc_ph
  use misc
  use filetools
  use channels
  use timingmod

  implicit none

contains

!#######################################################################

    subroutine adc1_vec

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer                               :: iadc1,dim1,i
      integer, dimension(:), allocatable    :: indx1
      real(dp), dimension(:,:), allocatable :: vec1

!-----------------------------------------------------------------------
! Open the ADC(1) eigenvector file
!-----------------------------------------------------------------------
      call freeunit(iadc1)
      open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',status='old')

!-----------------------------------------------------------------------
! Read the ADC(1) the eigenvectors
!-----------------------------------------------------------------------
      read(iadc1) dim1
      allocate(vec1(dim1,dim1))
      allocate(indx1(dim1))

      rewind(iadc1)

      read(iadc1) dim1,vec1

      vmat = vec1
!-----------------------------------------------------------------------
! Set the initial Davidson vectors
!-----------------------------------------------------------------------
!      do i=1,blocksize
!         vmat(1:dim1,i)=vec1(:,i)
!      enddo

!-----------------------------------------------------------------------
! Close the ADC(1) eigenvector file
!-----------------------------------------------------------------------
      close(iadc1)

      return

    end subroutine adc1_vec

!#######################################################################
!
!
!######################################################################

  subroutine mxv_diag_adc2_omp(ndim1,ndim2,kpq,nbuf,hxv_diag,chr)
  
    integer, intent(in) :: ndim1,ndim2,nbuf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: inda,indb,indj,indk,spin,count
    
    character(30) :: name
!    integer :: i,ktype,unt 
    real(dp), dimension(:), allocatable, intent(out) :: hxv_diag
    real(dp), dimension(:,:), allocatable            :: vmat

    allocate(hxv_diag(ndim1+ndim2))
    allocate(vmat(ndim1+ndim2,ndim1+ndim2))

    name="subspace - "//chr

    write(ilog,*) "Diagonal M-V multiplication in ",name

    call adc1_vec 

!!$ Filling the ph-ph block

    count = 0

    !$omp parallel do private(i,inda,indb,indj,indk,spin) shared(kpq,hxv_diag,vmat,count)
    do i=1, ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       hxv_diag(i) = K_ph_ph(e(inda),e(indj)) * vmat(i,i)
       hxv_diag(i) = hxv_diag(i) + C1_ph_ph(inda,indj,inda,indj) * vmat(i,i)
       hxv_diag(i) = hxv_diag(i) + CA_ph_ph(inda,inda) * vmat(i,i)
       hxv_diag(i) = hxv_diag(i) + CB_ph_ph(indj,indj) * vmat(i,i)
       hxv_diag(i) = hxv_diag(i) + CC_ph_ph(inda,indj,inda,indj) * vmat(i,i)
       count = count + 1
    end do
    !$omp end parallel do 

!!$ Filling the 2p2h-2p2h block

    !$omp parallel do private(i,inda,indb,indj,indk,spin) shared(kpq,hxv_diag,vmat,count)
    do i=ndim1+1, ndim1+ndim2
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin) 
       hxv_diag(i) = K_2p2h_2p2h(e(inda),e(indb),e(indj),e(indk)) * vmat(i,i)
       count = count + 1
    end do
    !$omp end parallel do
    
    write(ilog,*)count, "Diagonal M-V elements in ", name

    deallocate(hxv_diag)
    deallocate(vmat)

   return
  end subroutine mxv_diag_adc2_omp

!######################################################################

!#######################################################################
  
  subroutine get_offdiag_adc2_save_omp(ndim,kpq,nbuf,count,hxv_offdij,chr)

    use omp_lib
    use iomod

    implicit none
    
    integer, intent(in)                                 :: ndim
    integer, intent(out)                                :: nbuf
    integer*8, intent(out)                              :: count
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in)                            :: chr
    
    integer                      :: inda,indb,indj,indk,spin
    integer                      :: indapr,indbpr,indjpr,indkpr,spinpr
    
    character(30)                :: name
    integer                      :: i,j,k,dim_count,ndim1,unt
    real(dp)                     :: offdiag_ij1,offdiag_ij2,offdiag_ij3,&
                                    offdiag_ij4
    real(dp)                     :: offdiag1,offdiag2,offdiag3,offdiag4
   
    integer                               :: a,b,nzero !nocc
    real(dp), dimension(:,:), allocatable :: ca,cb,vmat
    real(dp)                              :: tw1,tw2,tc1,tc2

    integer                                       :: nthreads
    integer*8                                     :: nonzero
    integer                                       :: n,nprev,itmp
    real(dp), dimension(:,:), allocatable, intent(out) :: hxv_offdij

    real(dp) :: small

    integer :: c,cr,cm

!    buf_size2=buf_size
!    minc2=minc

    call times(tw1,tc1)

    name = "subspace-"//chr 
    
!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    
    write(ilog,*) "nthreads:",nthreads

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(hxv_offdij(ndim,ndim))
    allocate(vmat(ndim,ndim))
  
!-----------------------------------------------------------------------
! Precompute the results of calls to CA_ph_ph and CB_ph_ph
!-----------------------------------------------------------------------
    call times(tw1,tc1)

    allocate(ca(nvirt,nvirt),cb(nocc,nocc))

    !$omp parallel do private(i,j) shared(ca)
    ! CA_ph_ph
    do i=1,nvirt
       do j=i,nvirt
          ca(i,j)=CA_ph_ph(nocc+i,nocc+j)
          ca(j,i)=ca(i,j)
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(i,j) shared(cb)
    ! CB_ph_ph
    do i=1,nocc
       do j=i,nocc
          cb(i,j)=CB_ph_ph(i,j)
          cb(j,i)=cb(i,j)
       enddo
    enddo
    !$omp end parallel do

    call adc1_vec

!-----------------------------------------------------------------------
! Open the Hamiltonian file
!-----------------------------------------------------------------------
 
  write(ilog,*) "Matrix-Vector Multiplication in progress in",name

!-----------------------------------------------------------------------
! Initialise counters
!-----------------------------------------------------------------------
  count=0

!-----------------------------------------------------------------------
! ph-ph block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)

    hxv_offdij = 0    
    !$omp parallel do private(i,j,offdiag_ij1,offdiag_ij2,offdiag_ij3,offdiag_ij4,offdiag1,offdiag2,offdiag3,offdiag4,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,kpq,vmat,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)

          offdiag_ij1 = C1_ph_ph(inda,indj,indapr,indjpr) * vmat(j,i)
          offdiag1    = C1_ph_ph(inda,indj,indapr,indjpr) * vmat(i,j)

          if (indj .eq. indjpr) then
             offdiag_ij2 = ca(inda-nocc,indapr-nocc) * vmat(j,i)
             offdiag2    = ca(inda-nocc,indapr-nocc) * vmat(i,j)
          endif

          if (inda .eq. indapr) then
             offdiag_ij3 = cb(indj,indjpr) * vmat(j,i)
             offdiag3    = cb(indj,indjpr) * vmat(i,j)
          endif

          offdiag_ij4 = CC_ph_ph(inda,indj,indapr,indjpr) * vmat(j,i)
          offdiag4    = CC_ph_ph(inda,indj,indapr,indjpr) * vmat(i,j)

          hxv_offdij(i,j) = offdiag_ij1 + offdiag_ij2 + offdiag_ij3 + offdiag_ij4
          hxv_offdij(j,i) = offdiag1 + offdiag2 + offdiag3 + offdiag4

       end do
       count = count + 1
    end do
    !$omp end parallel do

    deallocate(ca,cb)
       
!-----------------------------------------------------------------------
! ph-2p2h block 
!-----------------------------------------------------------------------
!!$ Coupling to the i=j,a=b configs

    dim_count=kpq(1,0)

    !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,vmat,kpq,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(2,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr) 

          offdiag_ij1 = C5_ph_2p2h(inda,indj,indapr,indjpr) * vmat(j,i)
          offdiag1    = C5_ph_2p2h(inda,indj,indapr,indjpr) * vmat(i,j)

          hxv_offdij(i,j) = offdiag_ij1
          hxv_offdij(j,i) = offdiag1
          
       end do
       count = count + 1
    end do
    !$omp end parallel do

!!$ Coupling to the i=j,a|=b configs   
    
    dim_count=dim_count+kpq(2,0)
    
    !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,vmat,kpq,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(3,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          offdiag_ij1 = C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr) * vmat(j,i)
          offdiag1    = C4_ph_2p2h(inda,indj,indapr,indbpr,indjpr) * vmat(i,j)
    
          hxv_offdij(i,j) = offdiag_ij1
          hxv_offdij(j,i) = offdiag1
          
       end do
       count = count + 1
    end do
    !$omp end parallel do

!!$ Coupling to the i|=j,a=b configs
    
    dim_count=dim_count+kpq(3,0)
    
    !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,vmat,kpq,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(4,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          offdiag_ij1 = C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr) * vmat(j,i)
          offdiag1    = C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr) * vmat(i,j)

          hxv_offdij(i,j) = offdiag_ij1
          hxv_offdij(j,i) = offdiag1

       end do
       count = count + 1
    end do
    !$omp end parallel do

!!$ Coupling to the i|=j,a|=b I configs
       
    dim_count=dim_count+kpq(4,0)
    
    !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,vmat,kpq,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          offdiag_ij1 = C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(j,i)
          offdiag1    = C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(i,j)

          hxv_offdij(i,j) = offdiag_ij1
          hxv_offdij(j,i) = offdiag1

       end do
       count = count + 1
    end do
    !$omp end parallel do

!!$ Coupling to the i|=j,a|=b II configs
       
    dim_count=dim_count+kpq(5,0)

    !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,vmat,kpq,count)
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
       do j=dim_count+1,dim_count+kpq(5,0)
          call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

          offdiag_ij1 = C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(j,i)
          offdiag1    = C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(i,j)

          hxv_offdij(i,j) = offdiag_ij1
          hxv_offdij(j,i) = offdiag1
          
       end do
       count = count + 1
    end do
    !$omp end parallel do

    write(ilog,*) 2*count,' off-diagonal M-v elements in',name

!-----------------------------------------------------------------------    
! Deallocate arrays
!-----------------------------------------------------------------------

    deallocate(hxv_offdij)
    deallocate(vmat)

    call times(tw2,tc2)

    write(ilog,*) "Time taken:",tw2-tw1

    return
    
  end subroutine get_offdiag_adc2_save_omp

!#######################################################################

!#######################################################################

    subroutine get_offdiag_adc2_cvs_omp(ndim,kpq,nbuf,count,hxv_offdij,chr)
   
    use omp_lib
    use iomod
    
    implicit none

    integer, intent(in) :: ndim
    integer*8, intent(out) :: count
!    integer, intent(out) :: nbuf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    
    integer :: inda,indb,indj,indk,spin
    integer :: indapr,indbpr,indjpr,indkpr,spinpr 
    
    character(30) :: name
!    integer :: rec_count
    integer  :: i,j,k,nlim,dim_count,ndim1,unt
    integer  :: lim1i, lim2i, lim1j, lim2j
    real(dp) :: offdiag_ij1,offdiag_ij2,offdiag_ij3,offdiag_ij4,&
                offdiag1,offdiag2,offdiag3,offdiag4
    
!    integer, dimension(:), allocatable :: oi,oj
!    real(dp), dimension(:), allocatable :: file_offdiag
    
    integer                               :: a,b,nzero
    real(dp), dimension(:,:), allocatable :: ca,cb,vmat
    real(dp)                              :: tw1,tw2,tc1,tc2
    
!    integer                                       :: nthreads,tid
!    integer, dimension(:), allocatable            :: hamunit    
!    integer, dimension(:,:), allocatable          :: oi_omp,oj_omp
!    integer*8, dimension(:), allocatable          :: count_omp
!    integer, dimension(:), allocatable            :: rec_count_omp
!    integer, dimension(:), allocatable            :: nlim_omp
    integer*8                                     :: nonzero
    integer                                       :: n,nprev,itmp
    real(dp), dimension(:,:), allocatable         :: file_offdiag_omp,hxv_offdij
!    character(len=120), dimension(:), allocatable :: hamfile

!    integer :: buf_size2
!    real(dp) :: minc2

!    integer, dimension(:), allocatable :: nsaved

    integer :: c,cr,cm

!    buf_size2=buf_size
!    minc2=minc

    call times(tw1,tc1)

    name = "subspace-"//chr
!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
  !$omp parallel
  nthreads=omp_get_num_threads()
  !$omp end parallel

  write(ilog,*) "nthreads:",nthreads


!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------

  allocate(hxv_offdij(ndim,ndim))
  allocate(vmat(ndim,ndim))

!-----------------------------------------------------------------------
! Precompute the results of calls to CA_ph_ph and CB_ph_ph
!-----------------------------------------------------------------------
  allocate(ca(nvirt,nvirt),cb(nocc,nocc))
  
  !$omp parallel do private(i,j) shared(ca)
  ! CA_ph_ph
  do i=1,nvirt
     do j=i,nvirt
        ca(i,j)=CA_ph_ph(nocc+i,nocc+j)
        ca(j,i)=ca(i,j)
     enddo
  enddo
  !$omp end parallel do

  !$omp parallel do private(i,j) shared(cb)
  ! CB_ph_ph
  do i=1,nocc
     do j=i,nocc
        cb(i,j)=CB_ph_ph(i,j)
        cb(j,i)=cb(i,j)
     enddo
  enddo
  !$omp end parallel do

  call adc1_vec

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
 
  write(ilog,*) "Matrix-Vector Multiplication in progress in",name

!-----------------------------------------------------------------------
! Initialise counters
!-----------------------------------------------------------------------
  count=0

!-----------------------------------------------------------------------
! ph-ph block
!-----------------------------------------------------------------------

     hxv_offdij=0
     offdiag_ij=0
     ndim1=kpq(1,0)

     !$omp parallel do private(i,j,offdiag_ij1,offdiag_ij2,offdiag_ij3,offdiag_ij4,offdiag1,offdiag2,offdiag3,offdiag4,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,kpq,vmat,count)
     do i=1,ndim1
        call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
        do j=1,i-1
           call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)

           offdiag_ij1 = C1_ph_ph(inda,indj,indapr,indjpr) * vmat(j,i)
           offdiag1    = C1_ph_ph(inda,indj,indapr,indjpr) * vmat(i,j)

           if (indj .eq. indjpr) then
              offdiag_ij2 = ca(inda-nocc,indapr-nocc) * vmat(j,i)
              offdiag2    = ca(inda-nocc,indapr-nocc) * vmat(i,j)
           endif

           if (inda .eq. indapr) then
              offdiag_ij3 = cb(indj,indjpr) * vmat(j,i)
              offdiag3    = cb(indj,indjpr) * vmat(i,j)
           endif

           offdiag_ij4 = CC_ph_ph(inda,indj,indapr,indjpr) * vmat(j,i)
           offdiag4    = CC_ph_ph(inda,indj,indapr,indjpr) * vmat(i,j)

           hxv_offdij(i,j) = offdiag_ij1 + offdiag_ij2 + offdiag_ij3 + offdiag_ij4
           hxv_offdij(j,i) = offdiag1 + offdiag2 + offdiag3 + offdiag4

        end do
        count = count + 1
     end do
     !$omp end parallel do

!-----------------------------------------------------------------------
! ph-2p2h block 
!-----------------------------------------------------------------------
!!$ Coupling to the i|=j,a=b configs

       dim_count=kpq(1,0)
             
       !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,kpq,vmat,count)
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(4,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

             offdiag_ij1 = C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr) * vmat(j,i)
             offdiag1    = C3_ph_2p2h(inda,indj,indapr,indjpr,indkpr) * vmat(i,j)

             hxv_offdij(i,j) = offdiag_ij1
             hxv_offdij(j,i) = offdiag1
             
          end do
          count = count + 1
       end do
       !$omp end parallel do

!!$ Coupling to the i|=j,a|=b I configs

       dim_count=dim_count+kpq(4,0)

       !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,kpq,vmat,count)
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)  

             offdiag_ij1 = C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(j,i)
             offdiag1    = C1_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(i,j)

             hxv_offdij(i,j) = offdiag_ij1
             hxv_offdij(j,i) = offdiag1
             
          end do
          count = count + 1
       end do
       !$omp end parallel do

!!$ Coupling to the i|=j,a|=b II configs
       
       dim_count=dim_count+kpq(5,0)

     !$omp parallel do private(i,j,offdiag_ij1,offdiag1,inda,indb,indj,indk,spin,indapr,indbpr,indjpr,indkpr,spinpr) shared(hxv_offdij,kpq,vmat,count)
       do i=1,ndim1
          call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
          do j=dim_count+1,dim_count+kpq(5,0)
             call get_indices(kpq(:,j),indapr,indbpr,indjpr,indkpr,spinpr)

             offdiag_ij1 = C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(j,i)
             offdiag1    = C2_ph_2p2h(inda,indj,indapr,indbpr,indjpr,indkpr) * vmat(i,j)

             hxv_offdij(i,j) = offdiag_ij1
             hxv_offdij(j,i) = offdiag1

          end do
          count = count + 1
       end do
       !$omp end parallel do

!    write(ilog,*) 'counts'
    write(ilog,*) 2*count,' off-diagonal M-v elements in',name

!-----------------------------------------------------------------------    
! Deallocate arrays
!-----------------------------------------------------------------------

    deallocate(hxv_offdij)
    deallocate(vmat)

    call times(tw2,tc2) 
    write(ilog,*) "Time taken:",tw2-tw1

    return

  end subroutine get_offdiag_adc2_cvs_omp

!#######################################################################
end module mxv
