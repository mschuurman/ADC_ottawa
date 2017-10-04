  module get_matrix_dipole
  
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
    
    integer                                  :: nthreads
    real(d), parameter                       :: vectol=1e-8_d
    real(d), dimension(:,:), allocatable     :: pre_vv,pre_oo
    real(d), dimension(:,:,:,:), allocatable :: D261,D262,D263,D264

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

      ! Initial space 1h1p configs
      call dmatrix_f1h1p_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

      ! Initial space 2h2p i=j, a=b configs
      call dmatrix_f1h1p_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

      ! Initial space 2h2p i=j, a|=b configs
      call dmatrix_f1h1p_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
   
      ! Initial space 2h2p i|=j, a=b configs
      call dmatrix_f1h1p_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
      
      ! Initial space 2h2p i|=j, a|=b I configs
      call dmatrix_f1h1p_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)
      
      ! Initial space 2h2p i|=j, a|=b II configs
      call dmatrix_f1h1p_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
      
!-----------------------------------------------------------------------
! Final space 2h2p i=j, a=b configurations
! N.B. this block is zero under the CVS approximation
!-----------------------------------------------------------------------
!    if (.not.lcvsfinal) then

       ! Initial space 1h1p configs
       call dmatrix_f2h2p1_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

       ! Initial space 2h2p i=j, a=b configs
       call dmatrix_f2h2p1_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)
    
       ! Initial space 2h2p i=j, a|=b configs
       call dmatrix_f2h2p1_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
          
       ! Initial space 2h2p i|=j, a=b configs
       call dmatrix_f2h2p1_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p i|=j, a|=b I configs
       call dmatrix_f2h2p1_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p i|=j, a|=b II configs
       call dmatrix_f2h2p1_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
       
!    endif

!-----------------------------------------------------------------------
! Final space 2h2p i=j, a|=b configurations
! N.B. this block is zero under the CVS approximation
!-----------------------------------------------------------------------
!    if (.not.lcvsfinal) then

       ! Initial space 1h1p configs
       call dmatrix_f2h2p2_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p configs i=j, a=b
       call dmatrix_f2h2p2_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p configs i=j, a|=b
       call dmatrix_f2h2p2_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p configs i|=j, a=b
       call dmatrix_f2h2p2_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p configs i|=j, a|=b I
       call dmatrix_f2h2p2_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)
       
       ! Initial space 2h2p configs i|=j, a|=b II
       call dmatrix_f2h2p2_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
       
!    endif
    
!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a=b configurations
!-----------------------------------------------------------------------
    ! Initial space 1h1p configs
    call dmatrix_f2h2p3_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i=j, a=b configs
    call dmatrix_f2h2p3_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i=j, a|=b configs
    call dmatrix_f2h2p3_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a=b configs
    call dmatrix_f2h2p3_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a|=b I configs
    call dmatrix_f2h2p3_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a|=b II configs
    call dmatrix_f2h2p3_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
    
!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a|=b I configurations
!-----------------------------------------------------------------------
    ! Initial space 1h1p configs
    call dmatrix_f2h2p4I_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    ! Initial space 2h2p i=j, a=b configs
    call dmatrix_f2h2p4I_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i=j, a|=b configs
    call dmatrix_f2h2p4I_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a=b configs
    call dmatrix_f2h2p4I_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a|=b I configs
    call dmatrix_f2h2p4I_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)
       
    ! Initial space 2h2p i|=j, a|=b II configs
    call dmatrix_f2h2p4I_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
    
!-----------------------------------------------------------------------
! Final space 2h2p i|=j, a|=b II configurations
!-----------------------------------------------------------------------
    ! Initial space 1h1p configs
    call dmatrix_f2h2p4II_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i=j, a=b configs
    call dmatrix_f2h2p4II_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i=j, a|=b configs
    call dmatrix_f2h2p4II_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)
       
    ! Initial space 2h2p i|=j, a=b configs
    call dmatrix_f2h2p4II_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    ! Initial space 2h2p i|=j, a|=b I configs
    call dmatrix_f2h2p4II_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    ! Initial space 2h2p i|=j, a|=b II configs
    call dmatrix_f2h2p4II_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(pre_vv,pre_oo)
    if (allocated(D261)) deallocate(D261)
    if (allocated(D262)) deallocate(D262)
    if (allocated(D263)) deallocate(D263)
    if (allocated(D264)) deallocate(D264)

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine get_dipole_initial_product

!#######################################################################
  
  subroutine dmatrix_f1h1p_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)
    
    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                       :: ndim,ndimf
    integer                                   :: inda,indb,indk,indl,&
                                                 spin,indapr,indbpr,&
                                                 indkpr,indlpr,spinpr
    integer                                   :: i,j,ndim1,ndim1f,itmp,itmp1,tid
    real(d)                                   :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)      :: autvec
    real(d), dimension(ndimf), intent(out)    :: travec

    ndim1=kpq(1,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr,indkpr, &
    !$omp& indlpr,spinpr,itmp,itmp1,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=1,ndim1

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
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i1h1p

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,ndim1f
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,ar_offdiag_ij,inda,indb,indk,indl,spin, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(2,0)

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
   
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i2h2p1

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,ndim1f
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,ar_offdiag_ij,inda,indb,indk,indl,spin, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       
       do j=dim_count+1,dim_count+kpq(3,0)

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
            
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i2h2p2

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,ndim1f
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,ar_offdiag_ij,inda,indb,indk,indl,spin, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(4,0)

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
            
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i2h2p3

!#######################################################################

  subroutine dmatrix_f1h1p_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,ndim1f
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,ar_offdiag_ij,inda,indb,indk,indl,spin, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
    
       do j=dim_count+1,dim_count+kpq(5,0)

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
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i2h2p4I

!#######################################################################

 subroutine dmatrix_f1h1p_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,ndim1f
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    ndim1f=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,ar_offdiag_ij,inda,indb,indk,indl,spin, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(5,0)
 
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
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f1h1p_i2h2p4II

!#######################################################################

  ! PRIMED AND UNPRIMED INDICES SWITCH FROM HEREON IN...

  subroutine dmatrix_f2h2p1_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,ndim1,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1
       
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
         
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i1h1p

!#######################################################################
  
  subroutine dmatrix_f2h2p1_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)
    
    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
       
          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
       
          ar_offdiag_ij=&
               D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p1_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i2h2p2

!#######################################################################
 
  subroutine dmatrix_f2h2p1_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)
          
          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i2h2p3

!#######################################################################  

  subroutine dmatrix_f2h2p1_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p1_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
 
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
 
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p1_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p2_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,ndim1,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1
          
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
         
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(2,0)
          
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
       
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)

          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(4,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i2h2p3

!#######################################################################

    subroutine dmatrix_f2h2p2_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
          
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p2_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(5,0)
       
          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
      
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p2_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p3_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,ndim1,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=1,ndim1

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
         
         travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

      enddo
   enddo
   !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    

          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i2h2p1

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)

          if (abs(autvec(j)).lt.vectol) cycle
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)
              
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i2h2p3

!#######################################################################

  subroutine dmatrix_f2h2p3_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
       
          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i2h2p4I

!#######################################################################

    subroutine dmatrix_f2h2p3_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
       
          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
      
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p3_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p4I_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,ndim1,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1
          
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

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)
          
       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p4I_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
          
          if (abs(autvec(j)).lt.vectol) cycle
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
          
          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i2h2p1

!#######################################################################
  
  subroutine dmatrix_f2h2p4I_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(3,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p4I_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i2h2p3

!#######################################################################
  
  subroutine dmatrix_f2h2p4I_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
       
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p4I_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4I_i2h2p4II

!#######################################################################

  subroutine dmatrix_f2h2p4II_i1h1p(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,ndim1,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1

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
          
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i1h1p

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p1(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec
    
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
          
          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i2h2p1

!#######################################################################
  
  subroutine dmatrix_f2h2p4II_i2h2p2(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(3,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i2h2p2

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p3(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i2h2p3

!#######################################################################
  
  subroutine dmatrix_f2h2p4II_i2h2p4I(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          if (abs(autvec(j)).lt.vectol) cycle
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i2h2p4I

!#######################################################################

  subroutine dmatrix_f2h2p4II_i2h2p4II(ndim,ndimf,kpq,kpqf,travec,autvec)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer, intent(in)                    :: ndim,ndimf
    integer                                :: inda,indb,indk,indl,&
                                              spin,indapr,indbpr,&
                                              indkpr,indlpr,spinpr
    integer                                :: i,j,dim_count,dim_countf
    real(d)                                :: ar_offdiag_ij
    real(d), dimension(ndim), intent(in)   :: autvec
    real(d), dimension(ndimf), intent(out) :: travec

    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)

    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr, &
    !$omp& indkpr,indlpr,spinpr,ar_offdiag_ij) &
    !$omp& shared(kpq,kpqf,autvec,travec)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          if (abs(autvec(j)).lt.vectol) cycle

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
          ar_offdiag_ij=&
               D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          travec(i)=travec(i)+ar_offdiag_ij*autvec(j)

       enddo
    enddo
    !$omp end parallel do

    return

  end subroutine dmatrix_f2h2p4II_i2h2p4II

!#######################################################################

  subroutine density_matrix_ov_block

    use timingmod

    implicit none

    integer :: k,a
    real(d) :: tw1,tw2,tc1,tc2

!-----------------------------------------------------------------------
! Calculate occ-virt part of the ground state density matrix iff this
! not yet been done
!-----------------------------------------------------------------------
    if (.not.allocated(density)) then
     
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

    endif

    return

  end subroutine density_matrix_ov_block

!#######################################################################

  subroutine dmatrix_precalc(autvec,ndim,kpq,kpqf)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer, intent(in)                                 :: ndim
    integer                                             :: nvirt
    real(d)                                             :: mem4indx
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
    call dmatrix_precalc_4indx(nvirt,autvec,ndim,kpq,kpqf)

!-----------------------------------------------------------------------    
! Calculation of two-index terms, which can always be held in memory
!----------------------------------------------------------------------- 
    call dmatrix_precalc_2indx(nvirt,autvec,ndim,kpq,kpqf)

    return

  end subroutine dmatrix_precalc

!#######################################################################

  subroutine dmatrix_precalc_2indx(nvirt,autvec,ndim,kpq,kpqf)

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

    real(d), dimension(:,:), allocatable   :: tau_2_2_1,tau_2_2_2,&
                                              tau_4_2_1,tau_4_2_2
    real(d)                                :: tw1,tw2,tc1,tc2,ftmp1,ftmp2

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

  subroutine dmatrix_precalc_4indx(nvirt,autvec,ndim,kpq,kpqf)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer                                  :: nvirt,ndim,a,apr,k,&
                                                kpr,b,j
    integer, dimension(:,:), allocatable     :: iszero
    real(d), dimension(ndim), intent(in)     :: autvec
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

  subroutine dmatrix_precalc_noscreen(ndim,kpq,kpqf)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer, intent(in)                                 :: ndim
    integer                                             :: nvirt
    real(d)                                             :: mem4indx

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
    call dmatrix_precalc_4indx_noscreen(nvirt,ndim,kpq,kpqf)

!-----------------------------------------------------------------------    
! Calculation of two-index terms, which can always be held in memory
!----------------------------------------------------------------------- 
    call dmatrix_precalc_2indx_noscreen(nvirt,ndim,kpq,kpqf)

    return

  end subroutine dmatrix_precalc_noscreen

!#######################################################################

    subroutine dmatrix_precalc_2indx_noscreen(nvirt,ndim,kpq,kpqf)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer                                :: nvirt,ndim,a,apr,i,j,itmp,&
                                              itmp1,inda,indb,indk,indl,&
                                              spin,indapr,indbpr,indkpr,&
                                              indlpr,spinpr,b,b1,itmp2,k,&
                                              k1,kpr,j1
    integer, dimension(:,:), allocatable   :: iszeroa,iszerok

    real(d), dimension(:,:), allocatable   :: tau_2_2_1,tau_2_2_2,&
                                              tau_4_2_1,tau_4_2_2
    real(d)                                :: tw1,tw2,tc1,tc2,ftmp1,ftmp2

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
! Here we do not use any pre-screening, as the dipole matrix is to be
!  contracted with more than one state vector
!-----------------------------------------------------------------------
    iszeroa=0
    iszerok=0
    !$omp parallel do &
    !$omp& private(i,j,inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr,itmp,itmp1) &
    !$omp& shared(kpq,kpqf,iszeroa,iszerok)
    do i=1,kpqf(1,0)
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       do j=1,kpq(1,0)
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

  end subroutine dmatrix_precalc_2indx_noscreen

!#######################################################################

  subroutine dmatrix_precalc_4indx_noscreen(nvirt,ndim,kpq,kpqf)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq,kpqf

    integer                                  :: nvirt,ndim,a,apr,k,&
                                                kpr,b,j
    integer, dimension(:,:), allocatable     :: iszero
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

  end subroutine dmatrix_precalc_4indx_noscreen

!#######################################################################


!!$------------------------------------------------------------------------------
!!$------------------------------------ ADC1 ------------------------------------
!!$------------------------------------------------------------------------------

  subroutine get_dipole_initial_product_tda(ndim,ndimf,kpq,kpqf,autvec,travec)

! The difference from the earlier routine is that this routine returns the total
! number of saved els to a caller.

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

    use timingmod
    
    implicit none
    
    integer, intent(in) :: ndim
    integer, intent(in) :: ndimf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer                                             :: inda,indb,indk,indl,spin
    integer                                             :: indapr,indbpr,indkpr,&
                                                           indlpr,spinpr
    integer                                             :: i,j,ndim1,ndim1f
    real(d), dimension(ndimf,ndim), intent(out)         :: ar_offdiagd
    real(d)                                             :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Begin timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Calculate the D-matrix
!----------------------------------------------------------------------
    ar_offdiagd(:,:)=0.0d0

    ndim1=kpq(1,0)
    ndim1f=kpqf(1,0)
  
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=1,ndim1

          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
          
          if(indk.eq.indkpr) then
             ar_offdiagd(i,j)=ar_offdiagd(i,j)+D0_1_ph_ph(inda,indapr)
          endif
          
          if(inda.eq.indapr) then
             ar_offdiagd(i,j)=ar_offdiagd(i,j)+D0_2_ph_ph(indk,indkpr)
          endif

       enddo
    enddo

!----------------------------------------------------------------------
! Finish timing and output the time taken
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return
    
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

    implicit none
    
    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim,ndim), intent(out) :: ar_offdiagd
    
    integer :: inda,indb,indk,indl,spin
    integer :: indapr,indbpr,indkpr,indlpr,spinpr 
    
    integer :: i,j,dim_count,ndim1
    
    integer, dimension(buf_size) :: oi,oj
    real(d), dimension(buf_size) :: file_offdiagd
   
    integer :: lim1i,lim2i,lim1j,lim2j
 
    ar_offdiagd(:,:)=0.0d0
    
!!$ Full diagonalization. 

!!$ Filling the off-diagonal part of the ph-ph block

    ndim1=kpq(1,0)
    
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       do j=i+1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)             
          
          
          if (indk .eq. indkpr) then
             ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_1_ph_ph(inda,indapr)
          endif

          if(inda .eq. indapr) then
             ar_offdiagd(i,j) = ar_offdiagd(i,j) + D0_2_ph_ph(indk,indkpr)
          endif
          
          ar_offdiagd(j,i) = ar_offdiagd(i,j)
       enddo
    enddo

    return

  end subroutine get_offdiag_tda_DIPOLE_direct

!!$------------------------------------------------------------------------------
!!$---------------------------------  DIAGONAL PART  ----------------------------
!!$------------------------------------------------------------------------------
  subroutine get_diag_tda_DIPOLE_direct(ndim1,kpq,ar_diagd)

    implicit none
    
    integer, intent(in) :: ndim1
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    real(d), dimension(ndim1), intent(out) :: ar_diagd
    
    integer :: inda,indb,indk,indl,spin
    real(d) ::ea,eb,ej,ek,temp
    
    integer :: i,lim1,lim2
    
!!$ Filling the ph-ph block
    
    ar_diagd(:)=0.0d0
 
    do i= 1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)

       ar_diagd(i) = ar_diagd(i) + D0_1_ph_ph(inda,inda)
       ar_diagd(i) = ar_diagd(i) + D0_2_ph_ph(indk,indk)

    enddo

    return
    
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


  subroutine get_offdiag_tda_DIPOLE_SAVE_OK(ndimf,ndim,kpqf,kpq,nbuf,&
       count,filename)
    

    integer, dimension(7,0:nBas**2*nOcc**2) :: kpqf
    integer, dimension(7,0:nBas**2*nOcc**2) :: kpq
    
    integer                            :: ndimf
    integer                            :: ndim
    integer*8                          :: count
    integer                            :: nbuf
    integer                            :: unit_dip
    integer                            :: inda,indb,indk,indl,&
                                          spin
    integer                            :: indapr,indbpr,indkpr,indlpr,&
                                          spinpr
    integer                            :: i,j,nlim,dim_count,dim_countf,&
                                          ndim1,ndim1f
    integer                            :: lim1i,lim2i,lim1j,lim2j
    integer                            :: rec_count
    integer, dimension(:), allocatable :: oi,oj
    real(d), dimension(:), allocatable :: file_offdiagd
    real(d)                            :: ar_offdiagd_ij
    character(len=60)                  :: filename

!-----------------------------------------------------------------------
! Allocation and initialisation
!-----------------------------------------------------------------------
    allocate(oi(buf_size))
    allocate(oj(buf_size))
    allocate(file_offdiagd(buf_size))

    count=0
    rec_count=0
    nbuf=0
    oi(:)=0   
    oj(:)=0   
    file_offdiagd(:)=0.d0
    
!-----------------------------------------------------------------------
! Open the output file
!-----------------------------------------------------------------------
    call freeunit(unit_dip)
    open(unit_dip,file=filename,status='unknown',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Write the buffer size to file
!-----------------------------------------------------------------------
    write(unit_dip) buf_size
    
!-----------------------------------------------------------------------
! Calculate and save the ADC(1) dipole matrix
!-----------------------------------------------------------------------
    ndim1f = kpqf(1,0)
    ndim1 = kpq(1,0)
  
    do i=1,ndim1f
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       
       do j=1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
          
          ar_offdiagd_ij = 0.0d0
         
          if (indk .eq. indkpr) then
             ar_offdiagd_ij=ar_offdiagd_ij+D0_1_ph_ph(inda,indapr)
          endif
          
          if (inda.eq.indapr) then
             ar_offdiagd_ij=ar_offdiagd_ij+D0_2_ph_ph(indk,indkpr)
          endif

          ! Culling small matrix elements
          if (abs(ar_offdiagd_ij).gt.minc) then
             call register1()
          endif
          
       enddo

    enddo

    call register2()
    
    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ',&
         unit_dip

!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
    close(unit_dip)
    
!-----------------------------------------------------------------------
! Deallocatation
!-----------------------------------------------------------------------
    deallocate(oi)
    deallocate(oj)
    deallocate(file_offdiagd)
    
  contains

    subroutine register1()
      
      if (abs(ar_offdiagd_ij).gt.minc) then
         count=count+1
         ! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unit_dip,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         endif
      endif

    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unit_dip,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_offdiag_tda_DIPOLE_SAVE_OK

!######################################################################
  
  subroutine get_tda_dipole_samespace_save(kpq,nbuf,count,filename)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2) :: kpq
    integer*8                               :: count
    integer                                 :: nbuf
    integer                                 :: unit_dip
    integer                                 :: inda,indb,indk,indl,&
                                               spin
    integer                                 :: indapr,indbpr,indkpr,&
                                               indlpr,spinpr
    integer                                 :: i,j,nlim,dim_count,&
                                               dim_countf,ndim1
    integer                                 :: lim1i,lim2i,lim1j,lim2j
    integer                                 :: rec_count
    integer, dimension(:), allocatable      :: oi,oj
    real(d), dimension(:), allocatable      :: file_offdiagd
    real(d)                                 :: ar_offdiagd_ij
    character(len=60)                       :: filename

!-----------------------------------------------------------------------
! Allocation and initialisation
!-----------------------------------------------------------------------
    allocate(oi(buf_size))
    allocate(oj(buf_size))
    allocate(file_offdiagd(buf_size))

    count=0
    rec_count=0
    nbuf=0
    oi(:)=0
    oj(:)=0
    file_offdiagd(:)=0.d0
    
!-----------------------------------------------------------------------
! Open the output file
!-----------------------------------------------------------------------
    call freeunit(unit_dip)
    open(unit_dip,file=filename,status='unknown',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Write the buffer size to file
!-----------------------------------------------------------------------
    write(unit_dip) buf_size
    
!-----------------------------------------------------------------------
! Calculate and save the ADC(1) dipole matrix
!-----------------------------------------------------------------------
    ndim1 = kpq(1,0)
  
    do i=1,ndim1
       call get_indices(kpq(:,i),inda,indb,indk,indl,spin)
       
       do j=1,ndim1
          call get_indices(kpq(:,j),indapr,indbpr,indkpr,indlpr,spinpr)  
          
          ar_offdiagd_ij = 0.0d0
         
          if (indk .eq. indkpr) then
             ar_offdiagd_ij=ar_offdiagd_ij+D0_1_ph_ph(inda,indapr)
          endif
          
          if (inda.eq.indapr) then
             ar_offdiagd_ij=ar_offdiagd_ij+D0_2_ph_ph(indk,indkpr)
          endif

          ! Culling small matrix elements
          if (abs(ar_offdiagd_ij).gt.minc) then
             call register1()
          endif
          
       enddo

    enddo

    call register2()
    
    write(ilog,*) 'rec_counts',nbuf
    write(ilog,*) count,' DIPOLE_off-diagonal elements saved in file ',&
         unit_dip

!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
    close(unit_dip)
    
!-----------------------------------------------------------------------
! Deallocatation
!-----------------------------------------------------------------------
    deallocate(oi)
    deallocate(oj)
    deallocate(file_offdiagd)
    
  contains

    subroutine register1()
      
      if (abs(ar_offdiagd_ij).gt.minc) then
         count=count+1
         ! buf_size*int(rec_count,8) can exceed the int*4 limit
         file_offdiagd(count-buf_size*int(rec_count,8))=ar_offdiagd_ij
         oi(count-buf_size*int(rec_count,8))=i
         oj(count-buf_size*int(rec_count,8))=j
         !Checking if the buffer is full 
         if(mod(count,buf_size) .eq. 0) then
            rec_count=rec_count+1
            nlim=buf_size
            !Saving off-diag part in file
            call wrtoffdg(unit_dip,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
         endif
      endif

    end subroutine register1
    
    subroutine register2()
      
      !Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(unit_dip,buf_size,file_offdiagd(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
      
    end subroutine register2

  end subroutine get_tda_dipole_samespace_save

!######################################################################
  
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
  use timingmod

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

  real(d) :: tw1,tw2,tc1,tc2

  call times(tw1,tc1)

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

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'


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

!#######################################################################
! Improved routines for the calculation and out-of-core storage of
! the IS representation of the dipole operator
!#######################################################################

  subroutine get_adc2_dipole_improved(ndimf,ndim,kpqf,kpq,&
       nbuf,count,filename)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    integer, intent(in)          :: ndimf,ndim
    integer*8, intent(out)       :: count
    integer, intent(out)         :: nbuf
    integer                      :: inda,indb,indk,indl,spin,indapr,&
                                    indbpr,indkpr,indlpr,spinpr   
    integer                      :: i,j,nlim,dim_count,dim_countf,&
                                    ndim1,ndim1f
    integer                      :: lim1i,lim2i,lim1j,lim2j
    integer                      :: rec_count
    integer, dimension(buf_size) :: oi,oj
    real(d)                      :: ar_offdiag_ij
    real(d), dimension(buf_size) :: file_offdiag
    character(len=60)            :: filename

    ! NEW
    integer                      :: idpl
    integer                      :: nvirt,itmp,itmp1,dim
    real(d)                      :: func,mem4indx
    real(d)                      :: tw1,tw2,tc1,tc2

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
    call times(tw1,tc1)

    count=0
    rec_count=0
    nbuf=0
    oi(:)=0
    oj(:)=0
    file_offdiag(:)=0.d0

    nvirt=nbas-nocc

!-----------------------------------------------------------------------
! Open the dipole matrix file
!-----------------------------------------------------------------------
    call freeunit(idpl)
    open(idpl,file=filename,status='unknown',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Write the buffer size to file
!-----------------------------------------------------------------------
    write(idpl) buf_size

!-----------------------------------------------------------------------
! Calculate the density matrix
! Note that we only need to calculate the occupied-unoccupied part
!-----------------------------------------------------------------------
    call density_matrix_ov_block

!-----------------------------------------------------------------------
! Precalculation of function values
!-----------------------------------------------------------------------
    call dmatrix_precalc_noscreen(ndim,kpq,kpqf)

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 1h1p block
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Calculating the D-matrix...'

    ndim1=kpq(1,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=1,ndim1

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
          
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(2,0)
       
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
   
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       
       do j=dim_count+1,dim_count+kpq(3,0)

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
            
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(4,0)
            
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
            
          if (abs(ar_offdiag_ij).gt.minc) call register1()
           
       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a|=b I, block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
    
       do j=dim_count+1,dim_count+kpq(5,0)

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
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a|=b II, block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    ndim1f=kpqf(1,0)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(5,0)

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
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo


! PRIMED AND UNPRIMED INDICES SWITCH FROM HEREON IN...

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1
       
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
         
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
              
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
       
          ar_offdiag_ij=&
               D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1
       
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
         
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo
    
!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(2,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(5,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
      
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=1,ndim1
       
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
         
         if (abs(ar_offdiag_ij).gt.minc) call register1()

      enddo
   enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
   dim_count=kpq(1,0)
   dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
   do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    

          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
              
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
      
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1

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

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
                    
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
          
          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo


!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1

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
          
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
          ar_offdiag_ij=&
               D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          if (abs(ar_offdiag_ij).gt.minc) call register1()

       enddo
    enddo

!-----------------------------------------------------------------------
! Write the remaining non-zero elements to disk
!-----------------------------------------------------------------------
    call register2()
    
!-----------------------------------------------------------------------
! Close the dipole matrix file
!-----------------------------------------------------------------------
    close(idpl)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    if (allocated(pre_vv)) deallocate(pre_vv)
    if (allocated(pre_oo)) deallocate(pre_oo)
    if (allocated(D261)) deallocate(D261)
    if (allocated(D262)) deallocate(D262)
    if (allocated(D263)) deallocate(D263)
    if (allocated(D264)) deallocate(D264)

    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  contains

    subroutine register1()
      count=count+1
      ! buf_size*int(rec_count,8) can exceed the int*4 limit
      file_offdiag(count-buf_size*int(rec_count,8))=ar_offdiag_ij
      oi(count-buf_size*int(rec_count,8))=i
      oj(count-buf_size*int(rec_count,8))=j
      ! Checking if the buffer is full 
      if (mod(count,buf_size) .eq. 0) then
         rec_count=rec_count+1
         nlim=buf_size
         ! Saving off-diag part in file
         call wrtoffdg(idpl,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      endif
    end subroutine register1

    subroutine register2()
      ! Saving the rest in file
      nlim=count-buf_size*int(rec_count,8)
      call wrtoffdg(idpl,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
      rec_count=rec_count+1
      nbuf=rec_count
    end subroutine register2

  end subroutine get_adc2_dipole_improved

!#######################################################################
  
  subroutine get_adc2_dipole_improved_omp(ndimf,ndim,kpqf,kpq,&
       nbuf,count,filename)

    use timingmod

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    
    integer, intent(in)                :: ndimf,ndim
    integer*8, intent(out)             :: count
    integer, intent(out)               :: nbuf
    integer                            :: inda,indb,indk,indl,spin,&
                                          indapr,indbpr,indkpr,&
                                          indlpr,spinpr
    integer                            :: i,j,k,nlim,dim_count,&
                                          dim_countf,ndim1,ndim1f
    integer                            :: lim1i,lim2i,lim1j,lim2j
    integer                            :: rec_count
    integer, dimension(:), allocatable :: oi,oj
    real(d)                            :: ar_offdiag_ij
    real(d), dimension(:), allocatable :: file_offdiag
    character(len=60)                  :: filename

    integer                                       :: idpl
    integer                                       :: nvirt,itmp,itmp1,dim
    integer                                       :: tid
    integer, dimension(:), allocatable            :: dplunit    
    integer, dimension(:,:), allocatable          :: oi_omp,oj_omp
    integer*8, dimension(:), allocatable          :: count_omp
    integer, dimension(:), allocatable            :: rec_count_omp
    integer, dimension(:), allocatable            :: nlim_omp
    integer*8                                     :: nonzero
    integer                                       :: n,nprev
    integer, dimension(:), allocatable            :: nsaved
    real(d), dimension(:,:), allocatable          :: file_offdiag_omp
    real(d)                                       :: func,mem4indx
    real(d)                                       :: tw1,tw2,tc1,tc2
    character(len=120), dimension(:), allocatable :: dplfile

!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

    write(ilog,*) "nthreads:",nthreads

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)
    
!-----------------------------------------------------------------------
! Allocation and initialisation
!-----------------------------------------------------------------------
    allocate(dplunit(nthreads))
    allocate(dplfile(nthreads))
    allocate(oi_omp(nthreads,buf_size))
    allocate(oj_omp(nthreads,buf_size))
    allocate(file_offdiag_omp(nthreads,buf_size))
    allocate(count_omp(nthreads))
    allocate(rec_count_omp(nthreads))
    allocate(nlim_omp(nthreads))
    allocate(nsaved(nthreads))
    allocate(oi(buf_size))
    allocate(oj(buf_size))
    allocate(file_offdiag(buf_size))

    count=0
    nbuf=0
    rec_count=0
    oi=0
    oj=0
    file_offdiag=0.0d0
    rec_count_omp=0
    oi_omp=0
    oj_omp=0
    file_offdiag_omp=0.0d0
    nvirt=nbas-nocc
    
!-----------------------------------------------------------------------
! Open the working dipole files
!-----------------------------------------------------------------------
    do i=1,nthreads
       call freeunit(dplunit(i))
       dplfile(i)=trim(filename)//'.'
       k=len_trim(dplfile(i))+1
       if (i.lt.10) then
          write(dplfile(i)(k:k),'(i1)') i
       else
          write(dplfile(i)(k:k+1),'(i2)') i
       endif
       open(unit=dplunit(i),file=dplfile(i),status='unknown',&
            access='sequential',form='unformatted')
    enddo

!-----------------------------------------------------------------------
! Open the dipole matrix file
!-----------------------------------------------------------------------
    call freeunit(idpl)
    open(idpl,file=filename,status='unknown',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Write the buffer size to file
!-----------------------------------------------------------------------
    write(idpl) buf_size

!-----------------------------------------------------------------------
! Calculate the density matrix
! Note that we only need to calculate the occupied-unoccupied part
!-----------------------------------------------------------------------
    call density_matrix_ov_block

!-----------------------------------------------------------------------
! Precalculation of function values
!-----------------------------------------------------------------------
    call dmatrix_precalc_noscreen(ndim,kpq,kpqf)

!-----------------------------------------------------------------------
! Initialise counters
!-----------------------------------------------------------------------
    count=0
    rec_count=0

    count_omp=0
    rec_count_omp=0

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 1h1p block
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Calculating the D-matrix...'

    ndim1=kpq(1,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=1,ndim1

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
          
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(2,0)
       
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
   
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
       
       do j=dim_count+1,dim_count+kpq(3,0)

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
            
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f
       
       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(4,0)
            
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
            
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif
           
       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a|=b I, block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)
    
       do j=dim_count+1,dim_count+kpq(5,0)

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
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 1h1p, initial space 2h2p i|=j, a|=b II, block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    ndim1f=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& inda,indb,indk,indl,spin,indapr,indbpr,indkpr,indlpr,spinpr) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=1,ndim1f

       call get_indices(kpqf(:,i),inda,indb,indk,indl,spin)

       do j=dim_count+1,dim_count+kpq(5,0)

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
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

! PRIMED AND UNPRIMED INDICES SWITCH FROM HEREON IN...

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1
       
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
         
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
              
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
       
          ar_offdiag_ij=&
               D_1_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do
    
!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(2,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1
       
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
         
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    
!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(2,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_2_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
       
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_2_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do
    
!-----------------------------------------------------------------------
! Final space 2h2p, i=j, a|=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(3,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(5,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
      
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=1,ndim1
          
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
         
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif
          
       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)
       
       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    

          ar_offdiag_ij=&
               D_3_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(3,0)
          
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_3_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_3_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
              
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a=b, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(4,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)
       
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
      
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=1,ndim1

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

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)
                    
          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
          
          ar_offdiag_ij=&
               D_4i_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
       
       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4i_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  

          ar_offdiag_ij=&
               D_4i_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b I, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do
    
!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 1h1p block
!-----------------------------------------------------------------------
    ndim1=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=1,ndim1

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
          
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(2,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)    
 
          ar_offdiag_ij=&
               D_4ii_1_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i=j, a|=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)
    
       do j=dim_count+1,dim_count+kpq(3,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_2_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a=b block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(4,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
          
          ar_offdiag_ij=&
               D_4ii_3_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a|=b I block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
       
          ar_offdiag_ij=&
               D_4ii_4i_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)
 
          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Final space 2h2p, i|=j, a|=b II, initial space 2h2p i|=j, a|=b II block
!-----------------------------------------------------------------------
    dim_count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
    dim_countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
    !$omp parallel do &
    !$omp& private(tid,i,j,ar_offdiag_ij, &
    !$omp& indapr,indbpr,indkpr,indlpr,spinpr,inda,indb,indk,indl,spin) &
    !$omp& shared(kpq,kpqf,count_omp,file_offdiag_omp,rec_count_omp, &
    !$omp& nlim_omp,oi_omp,oj_omp,dplunit) &
    !$omp& firstprivate(dim_count,dim_countf,ndim1,ndim1f)
    do i=dim_countf+1,dim_countf+kpqf(5,0)

       call get_indices(kpqf(:,i),indapr,indbpr,indkpr,indlpr,spinpr)

       do j=dim_count+1,dim_count+kpq(5,0)

          call get_indices(kpq(:,j),inda,indb,indk,indl,spin)  
 
          ar_offdiag_ij=&
               D_4ii_4ii_2p2h_2p2h(inda,indb,indk,indl,indapr,indbpr,indkpr,indlpr)

          tid=1+omp_get_thread_num()

          if (abs(ar_offdiag_ij).gt.minc) then
             count_omp(tid)=count_omp(tid)+1
             file_offdiag_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=ar_offdiag_ij
             oi_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=i
             oj_omp(tid,count_omp(tid)-buf_size*int(rec_count_omp(tid),8))=j
             ! Checking if the buffer is full 
             if (mod(count_omp(tid),buf_size).eq.0) then
                rec_count_omp(tid)=rec_count_omp(tid)+1
                nlim_omp(tid)=buf_size
                ! Saving off-diag part in file
                call wrtoffdg(dplunit(tid),buf_size,&
                     file_offdiag_omp(tid,:),oi_omp(tid,:),&
                     oj_omp(tid,:),nlim_omp(tid))
             endif
          endif

       enddo
    enddo
    !$omp end parallel do

!-----------------------------------------------------------------------
! Assemble the complete dipole matrix file
!-----------------------------------------------------------------------
    write(ilog,'(2x,a)') "Dipole matrix assembly..."
    
    ! Total no. non-zero elements
    count=0
    do i=1,nthreads
       count=count+count_omp(i)
    enddo

    ! Complete records
    write(ilog,*) "       complete records..."
    do i=1,nthreads
       rewind(dplunit(i))
       do j=1,rec_count_omp(i)
          rec_count=rec_count+1
          read(dplunit(i)) file_offdiag(:),oi(:),oj(:),nlim
          call wrtoffdg(idpl,buf_size,file_offdiag(:),oi(:),oj(:),&
               buf_size)
       enddo
       close(dplunit(i))
    enddo

    ! Incomplete records
    write(ilog,*) "       incomplete records..."
    do i=1,nthreads
       nsaved(i)=mod(count_omp(i),buf_size)
    enddo
    n=nsaved(1)    
    file_offdiag(1:n)=file_offdiag_omp(1,1:n)
    oi(1:n)=oi_omp(1,1:n)
    oj(1:n)=oj_omp(1,1:n)
    nprev=n
    do i=2,nthreads

       n=n+nsaved(i)
              
       if (n.gt.buf_size) then
          ! The buffer is full. Write the buffer to disk and
          ! then save the remaining elements for thread i to the
          ! buffer
          !
          ! (i) Elements for thread i that can fit into the buffer
          itmp=buf_size-nprev
          file_offdiag(nprev+1:buf_size)=file_offdiag_omp(i,1:itmp)
          oi(nprev+1:buf_size)=oi_omp(i,1:itmp)
          oj(nprev+1:buf_size)=oj_omp(i,1:itmp)
          rec_count=rec_count+1
          call wrtoffdg(idpl,buf_size,file_offdiag(:),oi(:),oj(:),&
               buf_size)
          !
          ! (ii) Elements for thread i that couldn't fit into the buffer
          n=nsaved(i)-itmp
          file_offdiag(1:n)=file_offdiag_omp(i,itmp+1:nsaved(i))
          oi(1:n)=oi_omp(i,itmp+1:nsaved(i))
          oj(1:n)=oj_omp(i,itmp+1:nsaved(i))
       else
          ! The buffer is not yet full. Add all elements for thread i
          ! to the buffer
          file_offdiag(nprev+1:n)=file_offdiag_omp(i,1:nsaved(i))          
          oi(nprev+1:n)=oi_omp(i,1:nsaved(i))          
          oj(nprev+1:n)=oj_omp(i,1:nsaved(i))          
       endif

       nprev=n

    enddo

    ! Last, potentially incomplete buffer
    nlim=count-buf_size*int(rec_count,8)
    call wrtoffdg(idpl,buf_size,file_offdiag(:),oi(:),oj(:),nlim)
    rec_count=rec_count+1
    nbuf=rec_count
        
    ! Close the complete dipole matrix file
    close(idpl)
    
    ! Delete the working files
    !do i=1,nthreads
    !   call system('rm -rf '//trim(hamfile(i)))
    !enddo

    write(ilog,*) ' rec_counts',nbuf
    write(ilog,*) count,' off-diagonal elements saved in file ',filename

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    if (allocated(pre_vv)) deallocate(pre_vv)
    if (allocated(pre_oo)) deallocate(pre_oo)
    if (allocated(D261)) deallocate(D261)
    if (allocated(D262)) deallocate(D262)
    if (allocated(D263)) deallocate(D263)
    if (allocated(D264)) deallocate(D264)
    deallocate(dplunit)
    deallocate(dplfile)
    deallocate(oi_omp)
    deallocate(oj_omp)
    deallocate(file_offdiag_omp)
    deallocate(count_omp)
    deallocate(rec_count_omp)
    deallocate(nlim_omp)
    deallocate(nsaved)
    deallocate(oi)
    deallocate(oj)
    deallocate(file_offdiag)
    
!-----------------------------------------------------------------------
! Output the timing information
!-----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine get_adc2_dipole_improved_omp

!#######################################################################
 
end module get_matrix_dipole

