  module get_matrix_DIPOLE
  
    use constants
    use parameters
    use D_matrix
    use misc
    use filetools
    use dipole_ph
    use iomod
    use channels
    use omp_lib

    implicit none
    
    integer, parameter                       :: buf_size=8192
    integer                                  :: nthreads
    real(d), parameter                       :: vectol=1d-8
    real(d), dimension(:,:), allocatable     :: pre_vv,pre_oo
    real(d), dimension(:,:,:,:), allocatable :: D261,D262,D263,D264
    real(d), dimension(:), allocatable       :: sum1thread

  contains

!#######################################################################
! On the fly scalar product of the ADC(2) dipole matrix with the 
! initial state vector
!#######################################################################

    subroutine get_dipole_initial_product(ndim,ndimf,kpq,kpqf,autvec,&
         travec)

      use timingmod      

      implicit none
      
      
      integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
      integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
      integer, intent(in)                    :: ndim,ndimf
      integer                                :: inda,indb,indk,indl,&
                                                spin,indapr,indbpr,&
                                                indkpr,indlpr,spinpr
      integer                                :: i,j,nlim,rec_count,&
                                                dim_count,ndim1,&
                                                dim_countf,ndim1f
      integer                                :: k,k1,b,b1 
      real(d), dimension(ndim), intent(in)   :: autvec
      real(d), dimension(ndimf), intent(out) :: travec
      real(d)                                :: ar_offdiag_ij
      character(10)                          :: name

      integer                                :: nvirt,itmp,itmp1,dim
      real(d)                                :: func,mem4indx
      real(d)                                :: tw1,tw2,tc1,tc2

!-----------------------------------------------------------------------
! Output where we are at
!-----------------------------------------------------------------------
      write(ilog,'(/,90a)') ('-',i=1,90)
      write(ilog,'(2x,a)') & 
           'Contracting the IS representation of the dipole operator &
           with the initial state vector'
      write(ilog,'(90a)') ('-',i=1,90)

!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

    write(ilog,'(/,2x,a,x,i2)') "nthreads:",nthreads
    
    allocate(sum1thread(nthreads))
    
!-----------------------------------------------------------------------
! Initialise things
!-----------------------------------------------------------------------
      call times(tw1,tc1)
      nvirt=nbas-nocc
      travec(:)=0.0d0      

!-----------------------------------------------------------------------
! Calculate the density matrix
! Note that we only need to calculate the occupied-unoccupied part
!-----------------------------------------------------------------------
      call density_matrix_ov_block
    
!-----------------------------------------------------------------------
! Precalculation of function values
!-----------------------------------------------------------------------
      call dmatrix_precalc(autvec,ndim,kpq,kpqf)

!-----------------------------------------------------------------------
! Final space 1h1p configurations
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a)') 'Calculating the matrix-vector product...'

      ! Final space 1h1p configs
      ndim1f=kpqf(1,0)
      do i=1,ndim1f
         
         call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
         
         ! Initial space 1h1p configs
         call dmatrix_f1h1p_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)

         ! Initial space 2h2p i=j, a=b configs
         call dmatrix_f1h1p_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)

         ! Initial space 2h2p i=j, a|=b configs
         call dmatrix_f1h1p_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)
   
         ! Initial space i|=j,a=b configs
         call dmatrix_f1h1p_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)

         ! Initial space i|=j,a|=b I configs
         call dmatrix_f1h1p_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)
       
         ! Initial space i|=j,a|=b II configs
         call dmatrix_f1h1p_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
              inda,indb,indk,indl,spin)

    enddo

!-----------------------------------------------------------------------
! Final space 2h2p i=j, a=b configurations
! N.B. this block is zero under the CVS approximation
!-----------------------------------------------------------------------

    ! PRIMED AND UNPRIMED INDICES SWITCH FROM HERE...

    dim_countf=kpqf(1,0)

    if (.not.lcvsfinal) then

       do i=dim_countf+1,dim_countf+kpqf(2,0)
       
          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)   

          ! Initial space 1h1p configs
          call dmatrix_f2h2p1_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 2h2p i=j, a=b configs
          call dmatrix_f2h2p1_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)
    
          ! Initial space 2h2p i=j, a|=b configs
          call dmatrix_f2h2p1_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)
          
          ! Initial space 2h2p i|=j, a=b configs
          call dmatrix_f2h2p1_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)
       
          ! Initial space 2h2p i|=j, a|=b I configs
          call dmatrix_f2h2p1_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 2h2p i|=j, a|=b II configs
          call dmatrix_f2h2p1_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

       enddo
    
    endif

!-----------------------------------------------------------------------
! Final space 2h2p i=j, a|=b configurations
! N.B. this block is zero under the CVS approximation
!-----------------------------------------------------------------------
    dim_countf=dim_countf+kpqf(2,0)

    if (.not.lcvsfinal) then

       do i=dim_countf+1,dim_countf+kpqf(3,0)

          call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 1h1p configs
          call dmatrix_f2h2p2_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)
          
          ! Initial space 2h2p configs i=j, a=b
          call dmatrix_f2h2p2_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)
          
          ! Initial space 2h2p configs i=j, a|=b
          call dmatrix_f2h2p2_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 2h2p configs i|=j, a=b
          call dmatrix_f2h2p2_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 2h2p configs i|=j, a|=b I
          call dmatrix_f2h2p2_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

          ! Initial space 2h2p configs i|=j, a|=b II
          call dmatrix_f2h2p2_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
               indapr,indbpr,indkpr,indlpr,spinpr)

       enddo
    
    endif

!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a=b configurations
!-----------------------------------------------------------------------
    dim_countf=dim_countf+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 1h1p configs
       call dmatrix_f2h2p3_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a=b configs
       call dmatrix_f2h2p3_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a|=b configs
       call dmatrix_f2h2p3_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i|=j, a=b configs
       call dmatrix_f2h2p3_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i|=j, a|=b I configs
       call dmatrix_f2h2p3_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i|=j, a|=b II configs
       call dmatrix_f2h2p3_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

    enddo

!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a|=b I configurations
!-----------------------------------------------------------------------
    dim_countf=dim_countf+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 1h1p configs
       call dmatrix_f2h2p4I_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a=b configs
       call dmatrix_f2h2p4I_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a|=b configs
       call dmatrix_f2h2p4I_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)
    
       ! Initial space 2h2p i|=j, a=b configs
       call dmatrix_f2h2p4I_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)
       
       ! Initial space 2h2p i|=j, a|=b I configs
       call dmatrix_f2h2p4I_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i|=j, a|=b II configs
       call dmatrix_f2h2p4I_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)
       
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a|=b II configurations
!-----------------------------------------------------------------------
    dim_countf=dim_countf+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 1h1p configs
       call dmatrix_f2h2p4II_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a=b configs
       call dmatrix_f2h2p4II_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i=j, a|=b configs
       call dmatrix_f2h2p4II_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)
    
       ! Initial space 2h2p i|=j, a=b configs
       call dmatrix_f2h2p4II_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)
       
       ! Initial space 2h2p i|=j, a|=b I configs
       call dmatrix_f2h2p4II_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

       ! Initial space 2h2p i|=j, a|=b II configs
       call dmatrix_f2h2p4II_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
            indapr,indbpr,indkpr,indlpr,spinpr)

    enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(pre_vv,pre_oo)
    if (allocated(D261)) deallocate(D261)
    if (allocated(D262)) deallocate(D262)
    if (allocated(D263)) deallocate(D263)
    if (allocated(D264)) deallocate(D264)

    deallocate(sum1thread)

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine get_dipole_initial_product

!#######################################################################
  
  subroutine dmatrix_f1h1p_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
       inda,indb,indk,indl,spin)
    
    implicit none


    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                       :: i,ndim,ndimf
    integer                                   :: inda,indb,indk,indl,&
                                                 spin,indapr,indbpr,&
                                                 indkpr,indlpr,spinpr
    integer                                   :: j,ndim1,itmp,itmp1,tid
    real(d)                                   :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)      :: autvec
    real(d), dimension(ndimf), intent(out)    :: travec

    sum1thread=0.0d0

    ndim1=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij,itmp,itmp1) &
    !$omp& shared(autvec,kpq,pre_vv,pre_oo,sum1thread)
    do j=1,ndim1

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)
            
       call contract_4indx_dpl(ar_offdiag_ij,inda,indapr,indk,&
            indkpr)
            
       ar_offdiag_ij=ar_offdiag_ij+D2_7_1_ph_ph(inda,indapr,indk,indkpr)
       ar_offdiag_ij=ar_offdiag_ij+D2_7_2_ph_ph(inda,indapr,indk,indkpr)

       if(indk.eq.indkpr) then
          itmp=inda-nocc
          itmp1=indapr-nocc
          ar_offdiag_ij=ar_offdiag_ij+pre_vv(itmp,itmp1)
       endif
       
       if(inda.eq.indapr) then
          ar_offdiag_ij=ar_offdiag_ij+pre_oo(indk,indkpr)
       endif
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
      
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i1h1p

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       inda,indb,indk,indl,spin)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)

    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)    
       
       ar_offdiag_ij = 0.d0
       
       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D5_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D5_5_ph_2p2h(inda,indk,indbpr,indlpr)
          
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D5_6_ph_2p2h(inda,indk,indbpr,indkpr)

       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_3_ph_2p2h(inda,indk,indapr,indlpr)& 
            +D5_7_ph_2p2h(inda,indk,indapr,indlpr)
          
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D5_8_ph_2p2h(inda,indk,indapr,indkpr)

       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
   
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i2h2p1

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       inda,indb,indk,indl,spin)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle 

       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
            
       ar_offdiag_ij=0.d0
            
       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D4_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D4_5_ph_2p2h(inda,indk,indbpr,indlpr)
            
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D4_6_ph_2p2h(inda,indk,indbpr,indkpr)
            
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D4_7_ph_2p2h(inda,indk,indapr,indlpr)
            
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D4_8_ph_2p2h(inda,indk,indapr,indkpr)

       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
            
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
           
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i2h2p2

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       inda,indb,indk,indl,spin)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
            
       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
       
       ar_offdiag_ij=0.d0
            
       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D3_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D3_5_ph_2p2h(inda,indk,indbpr,indlpr)
            
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D3_6_ph_2p2h(inda,indk,indbpr,indkpr)

       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D3_7_ph_2p2h(inda,indk,indapr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D3_8_ph_2p2h(inda,indk,indapr,indkpr)
          
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
            
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
           
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i2h2p3

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
      inda,indb,indk,indl,spin)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)

    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  

       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D1_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D1_5_ph_2p2h(inda,indk,indbpr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D1_6_ph_2p2h(inda,indk,indbpr,indkpr)
       
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D1_7_ph_2p2h(inda,indk,indapr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D1_8_ph_2p2h(inda,indk,indapr,indkpr)
       
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+&
            D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
       
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+&
            D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
       
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
       
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i2h2p4I

!#######################################################################

 subroutine dmatrix_f1h1p_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
      inda,indb,indk,indl,spin)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
 
    !$omp parallel do &
    !$omp& private(j,tid,indapr,indbpr,indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
 
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
 
       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D2_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D2_5_ph_2p2h(inda,indk,indbpr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D2_6_ph_2p2h(inda,indk,indbpr,indkpr)
       
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D2_7_ph_2p2h(inda,indk,indapr,indlpr)
       
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D2_8_ph_2p2h(inda,indk,indapr,indkpr)
       
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
         
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
          
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
       
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij&
            +D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
      
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f1h1p_i2h2p4II

!#######################################################################

  ! PRIMED AND UNPRIMED INDICES SWITCH FROM HERE...

  subroutine dmatrix_f2h2p1_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,ndim1,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    ndim1=kpq(1,0)

    !$omp parallel do & 
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp shared(autvec,kpq,sum1thread)
    do j=1,ndim1
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=D5_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D5_5_ph_2p2h(inda,indk,indbpr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D5_6_ph_2p2h(inda,indk,indbpr,indkpr)
         
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_3_ph_2p2h(inda,indk,indapr,indlpr)& 
            +D5_7_ph_2p2h(inda,indk,indapr,indlpr)
         
       if((indk.eq.indlpr) .and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D5_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D5_8_ph_2p2h(inda,indk,indapr,indkpr)
       
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
         
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
         
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
         
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D5_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
         
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
        
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f2h2p1_i1h1p

!#######################################################################
  
  subroutine dmatrix_f2h2p1_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
       
       ar_offdiag_ij=D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p1_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p1_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
      
    enddo
    !$omp end parallel do

    do j=1,nthreads
       travec(i)=travec(i)+sum1thread(j)
    enddo

    return

  end subroutine dmatrix_f2h2p1_i2h2p2

!#######################################################################
 
  subroutine dmatrix_f2h2p1_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)
          
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    do j=1,nthreads
       travec(i)=travec(i)+sum1thread(j)
    enddo

    return

  end subroutine dmatrix_f2h2p1_i2h2p3

!#######################################################################  

  subroutine dmatrix_f2h2p1_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
       
    enddo
    !$omp end parallel do

    do j=1,nthreads
       travec(i)=travec(i)+sum1thread(j)
    enddo

    return

  end subroutine dmatrix_f2h2p1_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p1_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
 
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij=&
            D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
      
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f2h2p1_i2h2p4II

!#######################################################################

 subroutine dmatrix_f2h2p2_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
      indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,ndim1,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    ndim1=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=1,ndim1
          
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
       ar_offdiag_ij=0.d0
       
       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D4_5_ph_2p2h(inda,indk,indbpr,indlpr)
         
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D4_6_ph_2p2h(inda,indk,indbpr,indkpr)

       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D4_7_ph_2p2h(inda,indk,indapr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D4_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D4_8_ph_2p2h(inda,indk,indapr,indkpr)

       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)

       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)

       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)

       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D4_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
         
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)
          
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij=D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
       
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i2h2p3

!#######################################################################

    subroutine dmatrix_f2h2p2_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
         indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
      
       ar_offdiag_ij=&
            D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p2_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p3_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,ndim1,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    ndim1=kpq(1,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=1,ndim1

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D3_5_ph_2p2h(inda,indk,indbpr,indlpr)

       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D3_6_ph_2p2h(inda,indk,indbpr,indkpr)
         
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D3_7_ph_2p2h(inda,indk,indapr,indlpr)
         
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D3_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D3_8_ph_2p2h(inda,indk,indapr,indkpr)
         
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
       
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
       
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
 
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D3_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
         
!         travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
       
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p3_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    

       ar_offdiag_ij=D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p3_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p3_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)
              
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
       
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)
    
    return

  end subroutine dmatrix_f2h2p3_i2h2p3

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
      indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p3_i2h2p4I

!#######################################################################

    subroutine dmatrix_f2h2p3_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
         indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=&
            D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
      
!      travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p3_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p4I_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,ndim1,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    ndim1=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=1,ndim1
          
       tid=1+omp_get_thread_num()
       
       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_1_ph_2p2h(inda,indk,indbpr,indlpr)&
            +D1_5_ph_2p2h(inda,indk,indbpr,indlpr)
         
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D1_6_ph_2p2h(inda,indk,indbpr,indkpr)
         
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D1_7_ph_2p2h(inda,indk,indapr,indlpr)
         
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D1_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D1_8_ph_2p2h(inda,indk,indapr,indkpr)

       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D1_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
         
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D1_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
         
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D1_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
         
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D1_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p4I_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij=D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i2h2p1

!#######################################################################
  
  subroutine dmatrix_f2h2p4I_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p4I_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij=D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i2h2p3

!#######################################################################
  
  subroutine dmatrix_f2h2p4I_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)
       
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=&
            D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i2h2p4I

!#######################################################################

    subroutine dmatrix_f2h2p4I_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
         indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=&
            D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4I_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p4II_i1h1p(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,ndim1,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    ndim1=kpq(1,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=1,ndim1

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)             
       
       ar_offdiag_ij=0.d0

       if((indk.eq.indkpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_1_ph_2p2h(inda,indk,indbpr,indlpr)&
                 +D2_5_ph_2p2h(inda,indk,indbpr,indlpr)
          
       if((indk.eq.indlpr).and.(inda.eq.indapr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_2_ph_2p2h(inda,indk,indbpr,indkpr)&
            +D2_6_ph_2p2h(inda,indk,indbpr,indkpr)
          
       if((indk.eq.indkpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_3_ph_2p2h(inda,indk,indapr,indlpr)&
            +D2_7_ph_2p2h(inda,indk,indapr,indlpr)
          
       if((indk.eq.indlpr).and.(inda.eq.indbpr))&
            ar_offdiag_ij=ar_offdiag_ij+D2_4_ph_2p2h(inda,indk,indapr,indkpr)&
            +D2_8_ph_2p2h(inda,indk,indapr,indkpr)
          
       if(inda.eq.indapr)&
            ar_offdiag_ij=ar_offdiag_ij+D2_9_ph_2p2h(inda,indk,indbpr,indkpr,indlpr)
          
       if(inda.eq.indbpr)&
            ar_offdiag_ij=ar_offdiag_ij+D2_10_ph_2p2h(inda,indk,indapr,indkpr,indlpr)
          
       if(indk.eq.indkpr)&
            ar_offdiag_ij=ar_offdiag_ij+D2_11_ph_2p2h(inda,indk,indapr,indbpr,indlpr)
          
       if(indk.eq.indlpr)&
            ar_offdiag_ij=ar_offdiag_ij+D2_12_ph_2p2h(inda,indk,indapr,indbpr,indkpr)
          
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)
      
    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p1(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    sum1thread=0.0d0

    dim_count=kpq(1,0)
    
    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(2,0)
          
       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
       ar_offdiag_ij=&
            D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i2h2p1

!#######################################################################
  
  subroutine dmatrix_f2h2p4II_i2h2p2(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(3,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

       ar_offdiag_ij=&
            D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p3(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(4,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij=&
            D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i2h2p3

!#######################################################################
  
  subroutine dmatrix_f2h2p4II_i2h2p4I(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle
       
       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
       ar_offdiag_ij=&
            D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p4II(i,ndim,ndimf,kpq,travec,autvec,&
       indapr,indbpr,indkpr,indlpr,spinpr)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)                    :: ndim,ndimf,i
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: j,dim_count,tid
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    sum1thread=0.0d0
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)

    !$omp parallel do &
    !$omp& private(j,tid,inda,indb,indk,indl,spin,ar_offdiag_ij) &
    !$omp& shared(autvec,kpq,sum1thread)
    do j=dim_count+1,dim_count+kpq(5,0)

       tid=1+omp_get_thread_num()

       if (abs(autvec(j)).lt.vectol) cycle

       call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
       ar_offdiag_ij=&
            D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

!       travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       sum1thread(tid)=sum1thread(tid)+ar_offdiag_ij*autvec(j)

    enddo
    !$omp end parallel do

    travec(i)=travec(i)+sum(sum1thread)

    return

  end subroutine dmatrix_f2h2p4II_i2h2p4II

!#######################################################################

  subroutine density_matrix_ov_block

    use timingmod

    implicit none

    integer :: k,a
    real(d) :: tw1,tw2,tc1,tc2

    call times(tw1,tc1)

    write(ilog,'(/,2x,a)') 'Calculating the occupied-unoccupied part of &
         the density matrix...'

    allocate(density(nbas,nbas))
    density=0.0d0

    !$omp parallel do private(k,a) shared(density)
    do k=1,nocc
       do a=nocc+1,nbas
         density(k,a)=calc_density(k,a)
       enddo
    enddo
    !$omp end parallel do

    call times(tw2,tc2)
    
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine density_matrix_ov_block

!#######################################################################

  subroutine dmatrix_precalc(autvec,ndim,kpq,kpqf)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer, intent(in)                                 :: ndim
    integer                                             :: nvirt
    real(d)                                             :: vectol,mem4indx
    real(d), dimension(ndim), intent(in)                :: autvec

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
    nvirt=nbas-nocc
    allocate(pre_vv(nvirt,nvirt),pre_oo(nocc,nocc))
    pre_vv=0.0d0
    pre_oo=0.0d0

    allocate(D261(nvirt,nocc,nocc,nvirt))
    allocate(D262(nvirt,nocc,nocc,nvirt))
    allocate(D263(nvirt,nvirt,nocc,nocc))
    allocate(D264(nvirt,nvirt,nocc,nocc))
    D261=0.0d0
    D262=0.0d0
    D263=0.0d0
    D264=0.0d0       
    
!-----------------------------------------------------------------------
! Calculation of intermediate four-index terms, which may need to be 
! saved to file
!-----------------------------------------------------------------------
    call dmatrix_precalc_4indx(nvirt,autvec,ndim,kpq,kpqf,vectol)

!-----------------------------------------------------------------------    
! Calculation of two-index terms, which can always be held in memory
!----------------------------------------------------------------------- 
    call dmatrix_precalc_2indx(nvirt,autvec,ndim,kpq,kpqf,vectol)

    return

  end subroutine dmatrix_precalc

!#######################################################################

  subroutine dmatrix_precalc_2indx(nvirt,autvec,ndim,kpq,kpqf,vectol)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer                                :: nvirt,ndim,a,apr,i,j,itmp,&
                                              itmp1,inda,indb,indk,indl,&
                                              spin,indapr,indbpr,indkpr,&
                                              indlpr,spinpr,b,b1,itmp2,k,&
                                              k1,kpr,j1
    integer, dimension(:,:), allocatable   :: iszeroa,iszerok
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), intent(in)                    :: vectol

    real(d), dimension(:,:), allocatable :: tau_2_2_1,tau_2_2_2,&
                                            tau_4_2_1,tau_4_2_2
    real(d)                              :: tw1,tw2,tc1,tc2,ftmp1,ftmp2

    call times(tw1,tc1)

    write(ilog,'(/,2x,a)') 'Precomputing two-index terms...'

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(iszeroa(nvirt,nvirt))
    allocate(iszerok(nocc,nocc))
    allocate(tau_2_2_1(2*nvirt,2*nvirt))
    allocate(tau_2_2_2(2*nvirt,2*nvirt))
    allocate(tau_4_2_1(2*nocc,2*nocc))
    allocate(tau_4_2_2(2*nocc,2*nocc))

!-----------------------------------------------------------------------
! Here we only bother computing terms that multiply elements of the
! initial state vector (autvec) whose magnitude is above the threshold 
! value (vectol)
!-----------------------------------------------------------------------
    iszeroa=0
    iszerok=0
    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr,itmp,itmp1) &
    !$omp& shared(kpq,kpqf,autvec,iszeroa,iszerok)
    do i=1,kpqf(1,0)
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       do j=1,kpq(1,0)
          if (abs(autvec(j)).lt.vectol) cycle
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)
          itmp=inda-nocc
          itmp1=indapr-nocc
          iszeroa(itmp,itmp1)=1
          iszerok(indk,indkpr)=1
       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Unoccupied-unoccupied terms
!-----------------------------------------------------------------------
    !$omp parallel do &
    !$omp& private(a,itmp,b,itmp1) &
    !$omp& shared(tau_2_2_1,tau_2_2_2)
    do a=nocc+1,nbas
       itmp=a-nocc
       do b=nocc+1,nbas
          itmp1=b-nocc
          tau_2_2_1(itmp,itmp1)=tau_D2_2_1(a,b)
          tau_2_2_2(itmp,itmp1)=tau_D2_2_2(a,b) 
      enddo
    enddo
    !$omp end parallel do

    pre_vv=0.0d0
    !$omp parallel do &
    !$omp& private(a,itmp,apr,itmp1,ftmp1,ftmp2,b1,b,itmp2) &
    !$omp& shared(iszeroa,pre_vv,roccnum,tau_2_2_1,tau_2_2_2,dpl)
    do a=nocc+1,nbas
       itmp=a-nocc
       do apr=nocc+1,nbas
          itmp1=apr-nocc

          if (iszeroa(itmp,itmp1).eq.0) cycle

          ! Contributions from D0_1_ph_ph, D2_1_ph_ph, D2_2_ph_ph,
          ! D2_3_1_ph_ph and D2_3_2_ph_ph
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+D0_1_ph_ph(a,apr)
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+D2_1_ph_ph(a,apr)
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+D2_2_ph_ph(a,apr)
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+D2_3_1_ph_ph(a,apr)
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+D2_3_2_ph_ph(a,apr)

          ! Contributions from D2_2_1_ph_ph and D2_2_2_ph_ph
          ftmp1=0.0d0
          ftmp2=0.0d0
          do b1=nocc+1,nbas
             b=roccnum(b1)
             itmp2=b1-nocc
             ftmp1=ftmp1+tau_2_2_1(itmp,itmp2)*dpl(b,apr)
             ftmp2=ftmp2+tau_2_2_2(itmp1,itmp2)*dpl(b,a)
          enddo
          ftmp1=-0.125d0*ftmp1
          ftmp2=-0.125d0*ftmp2
          pre_vv(itmp,itmp1)=pre_vv(itmp,itmp1)+ftmp1+ftmp2

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Occupied-occupied terms
!-----------------------------------------------------------------------    
    !$omp parallel do private(j,k) shared(tau_4_2_1,tau_4_2_2)
    do j=1,nocc
       do k=1,nocc
          tau_4_2_1(j,k)=tau_D2_4_1(j,k)
          tau_4_2_2(j,k)=tau_D2_4_2(j,k)
       enddo
    enddo
    !$omp end parallel do

    pre_oo=0.0d0
    !$omp parallel do &
    !$omp& private(k,kpr,ftmp1,ftmp2,j1,j) &
    !$omp& shared(iszerok,pre_oo,roccnum,tau_4_2_1,tau_4_2_2,dpl)
    do k=1,nocc
       do kpr=1,nocc

          if (iszerok(k,kpr).eq.0) cycle

          ! Contributions from D0_2_ph_ph, D2_3_ph_ph, D2_4_ph_ph,
          ! D2_5_1_ph_ph and D2_5_2_ph_ph
          pre_oo(k,kpr)=pre_oo(k,kpr)+D0_2_ph_ph(k,kpr)
          pre_oo(k,kpr)=pre_oo(k,kpr)+D2_3_ph_ph(k,kpr)
          pre_oo(k,kpr)=pre_oo(k,kpr)+D2_4_ph_ph(k,kpr)
          pre_oo(k,kpr)=pre_oo(k,kpr)+D2_5_1_ph_ph(k,kpr)
          pre_oo(k,kpr)=pre_oo(k,kpr)+D2_5_2_ph_ph(k,kpr)
          ! Contributions from D2_4_1_ph_ph and D2_4_2_ph_ph
          ftmp1=0.0d0
          ftmp2=0.0d0
          do j1=1,nocc
             j=roccnum(j1)
             ftmp1=ftmp1+tau_4_2_1(j,k)*dpl(kpr,j)
             ftmp2=ftmp2+tau_4_2_2(j,kpr)*dpl(k,j)
          enddo
          ftmp1=0.125d0*ftmp1
          ftmp2=0.125d0*ftmp2
          pre_oo(k,kpr)=pre_oo(k,kpr)+ftmp1+ftmp2

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(iszeroa)
    deallocate(iszerok)
    deallocate(tau_2_2_1)
    deallocate(tau_2_2_2)
    deallocate(tau_4_2_1)
    deallocate(tau_4_2_2)

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine dmatrix_precalc_2indx

!#######################################################################

  subroutine dmatrix_precalc_4indx(nvirt,autvec,ndim,kpq,kpqf,vectol)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer                                  :: nvirt,ndim,a,apr,k,&
                                                kpr,b,j
    integer, dimension(:,:), allocatable     :: iszero
    real(d), dimension(ndim), intent(in)     :: autvec
    real(d), intent(in)                      :: vectol
    real(d)                                  :: tw1,tw2,tc1,tc2
    character(len=60)                        :: filename

    call times(tw1,tc1)

    write(ilog,'(/,2x,a)') 'Precomputing four-index terms...'

!-----------------------------------------------------------------------
! (1) D2.6.1.akkprb
!-----------------------------------------------------------------------
    !$omp parallel do private(a,k,kpr,b) shared(D261)
    do a=nocc+1,nbas
       do k=1,nocc
          do kpr=1,nocc
             do b=nocc+1,nbas
                D261(a-nocc,k,kpr,b-nocc)=tau_D2_6_1(a,k,kpr,b)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! (2) D2.6.2.aprkkprb
!-----------------------------------------------------------------------
    !$omp parallel do private(apr,k,kpr,b) shared(D262)
    do apr=nocc+1,nbas
       do k=1,nocc
          do kpr=1,nocc
             do b=nocc+1,nbas
                D262(apr-nocc,k,kpr,b-nocc)=tau_D2_6_2(apr,k,kpr,b)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! (3) D2.6.3.aprkkprb
!-----------------------------------------------------------------------
    !$omp parallel do private(a,apr,k,j) shared(D263)
    do a=nocc+1,nbas
       do apr=nocc+1,nbas
          do k=1,nocc
             do j=1,nocc
                D263(a-nocc,apr-nocc,k,j)=tau_D2_6_3(a,apr,k,j)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    
!-----------------------------------------------------------------------
! (4) D2.6.4.aaprkprj
!-----------------------------------------------------------------------
    !$omp parallel do private(a,apr,kpr,j) shared(D264)
    do a=nocc+1,nbas
       do apr=nocc+1,nbas
          do kpr=1,nocc
             do j=1,nocc
                D264(a-nocc,apr-nocc,kpr,j)=tau_D2_6_4(a,apr,kpr,j)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine dmatrix_precalc_4indx

!#######################################################################

  subroutine contract_4indx_dpl(func,a,apr,k,kpr)

    implicit none

    real(d)                                  :: func,ftmp,curr
    integer                                  :: a,apr,k,kpr,itmp,b,b1,j,j1

!-----------------------------------------------------------------------
! Function value to be returned: func
!-----------------------------------------------------------------------
    func=0.0d0

!-----------------------------------------------------------------------
! (1) D2.6.1.akkprb
!-----------------------------------------------------------------------
    ! Contract the intermediate four-index terms with the dipole matrix
    curr=0.0d0    
    do b1=nocc+1,nbas
       b=roccnum(b1)
       curr=curr+dpl(b,apr)*D261(a-nocc,k,kpr,b1-nocc)
    enddo
    
    ! Normalisation
    curr=0.25d0*curr

    func=func+curr

!-----------------------------------------------------------------------
! (2) D2.6.2.aprkkprb
!-----------------------------------------------------------------------
    ! Contract the intermediate four-index terms with the dipole matrix
    curr=0.0d0
    do b1=nocc+1,nbas
       b=roccnum(b1)
       curr=curr+dpl(b,a)*D262(apr-nocc,k,kpr,b1-nocc)
    enddo
    
    ! Normalisation
    curr=0.25d0*curr

    func=func+curr

!-----------------------------------------------------------------------
! (3) D2.6.3.aprkkprb
!-----------------------------------------------------------------------
    ! Contract the intermediate four-index terms with the dipole matrix
    curr=0.0d0
    do j1=1,nocc
       j=roccnum(j1)
       curr=curr+dpl(kpr,j)*D263(a-nocc,apr-nocc,k,j1)
    enddo

    ! Normalisation
    curr=-0.25d0*curr

    func=func+curr

!-----------------------------------------------------------------------
! (4) D2.6.4.aaprkprj
!-----------------------------------------------------------------------
    ! Contract the intermediate four-index terms with the dipole matrix
    curr=0.0d0
    do j1=1,nocc
       j=roccnum(j1)
       curr=curr+dpl(k,j)*D264(a-nocc,apr-nocc,kpr,j1)
    enddo

    ! Normalisation
    curr=-0.25d0*curr

    func=func+curr

    return

  end subroutine contract_4indx_dpl

!#######################################################################


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
    
    write(ilog,*) "Writing the travec vector of ADC-DIPOLE matrix INITIAL-STATE product "



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

!!!write(ilog,*) 'D_2_2_2p2h_2p2h_offdiag', ar_offdiagd(i,j)
           
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


!!!write(ilog,*) "D0_1",D0_1_ph_ph(inda,inda)            
!!!write(ilog,*) "D0_2",D0_2_ph_ph(indk,indk)            
!!!write(ilog,*) "Ddiag_ph_ph", ar_diagd(i)
 
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

!!!write(ilog,*) 'D_2_2_2p2h_2p2h_offdiag', ar_offdiagd(i,j)
           
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

!!!write (ilog,*) "Ddiag_2_2_2p2h_2p2h",inda,indb,indk,indl, ar_diagd(i)

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

  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j+1,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i , j,'one:', ar_offdiagd_ij
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
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", UNIT_DIP

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


    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', UNIT_DIP


  contains
       

    subroutine register1()

      if (abs(ar_offdiagd_ij) .gt. minc) then
         count=count+1
         IF ( count .eq. 1 ) then
         write(ilog,*) 'the first element not-diagonal saved in', UNIT_DIP,'is the', i+1 , j+1,'one:', ar_offdiagd_ij
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
    write(ilog,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_diag(1)
    write(ilog,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_diag( ndim )



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
    write(ilog,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_diag(1)
    write(ilog,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_diag( ndim + 1 )

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
    write(ilog,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_offdiag(1)
    write(ilog,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_offdiag( ndim )

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

    write(ilog,*) 'the first element diagonal saved in', UNIT_DIP,'is ', ar_offdiag(1)
    write(ilog,*) 'the last  element diagonal saved in', UNIT_DIP,'is ', ar_offdiag( ndim + 1 )

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

  name="SCRATCH/hmlt.off"//chr
  unt=22
  
  count=0
  rec_count=0
  
  write(ilog,*) "Writing the off-diagonal part of ADC DIPOLE matrix in file ", name
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
    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ', name

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
    name="SCRATCH/hmlt.dia"//chr 
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

!!!write (ilog,*) "Ddiag_2_2_2p2h_2p2h",inda,indb,indk,indl, ar_diagd(i)

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
    write(ilog,*) "Writing",ndim1+ndim2," diagonal elements of ADC-DIPOLE ADC2 matrix in file ",name
    OPEN(UNIT=unt,FILE=name,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    call wrtdg(unt,ndim1+ndim2,buf_size,nbuf,ktype,ar_diagd(:))
!!$    call wrtdgat(unt,ndim1+ndim2,nbuf,ar_diag(:))
    CLOSE(unt)
   
    write(ilog,*) 'Writing successful at get_diag_adc2_DIPOLE_save end'
 
  end subroutine get_diag_adc2_DIPOLE_save

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------

!!$------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------


 
end module get_matrix_DIPOLE    

