module fspace
  
  use constants
  use parameters
  use get_matrix
  use get_moment
  use select_fano
  use filetools
  use misc
  use get_matrix_DIPOLE  
  use channels

  implicit none
  
contains



!!!!!!!!! DIPOLE MATRIX ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$----------------------------------------------------------

  subroutine get_fspace_adc2_DIPOLE_direct(ndim,kpq,autvec,arrd,travec) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: travec
    real(d), dimension(ndim,ndim), intent(inout) :: arrd
    real(d), dimension(ndim), intent(in) :: autvec
   
    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diagd,temp
    real(d), dimension(:,:), allocatable :: ar_offdiagd
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diagd(ndim),ar_offdiagd(ndim,ndim))
    
    call get_offdiag_adc2_DIPOLE_direct(ndim,kpq(:,:),ar_offdiagd(:,:))
    call get_diag_adc2_DIPOLE_direct(ndim1,ndim2,kpq(:,:),ar_diagd(:))
    
    arrd(:,:)=ar_offdiagd(:,:)
    
    do i=1,ndim
       arrd(i,i)=ar_diagd(i)
    end do
    
    deallocate(ar_diagd,ar_offdiagd)

!!!    call dmatvec(ndim,autvec(:),arrd(:,:),travec(:))
    call mat_vec_multiply(ndim,ndim,arrd,autvec,travec)




!!! up to this point I have the vector  TRAVEC to make the scalar product with all other
!!! adc excited eigenstates, i.e. with the other eigenvectors of the hamiltonian
!!! adc matrix



!    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2_DIPOLE_direct


!!$----------------------------------------------------------

  subroutine get_fspace_tda_DIPOLE_direct(ndim,kpq,autvec,arrd,travec) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: travec
    real(d), dimension(ndim,ndim), intent(inout) :: arrd
    real(d), dimension(ndim), intent(in) :: autvec
   
    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diagd,temp
    real(d), dimension(:,:), allocatable :: ar_offdiagd
    
    ndim1=kpq(1,0)
    
    allocate(ar_diagd(ndim),ar_offdiagd(ndim,ndim))
    
    call get_offdiag_tda_DIPOLE_direct(ndim,kpq(:,:),ar_offdiagd(:,:))
    call get_diag_tda_DIPOLE_direct(ndim1,kpq(:,:),ar_diagd(:))
    
    arrd(:,:)=ar_offdiagd(:,:)
    
    do i=1,ndim
       arrd(i,i)=ar_diagd(i)
    end do
    
    deallocate(ar_diagd,ar_offdiagd)

!!!    call dmatvec(ndim,autvec(:),arrd(:,:),travec(:))
    call mat_vec_multiply(ndim,ndim,arrd,autvec,travec)




!!! up to this point I have the vector  TRAVEC to make the scalar product with all other
!!! adc excited eigenstates, i.e. with the other eigenvectors of the hamiltonian
!!! adc matrix



!    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_tda_DIPOLE_direct

!
!!$----------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine get_fspace_tda_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,autvec,arrd,travec) 

    integer, intent(in) :: ndim
    integer, intent(in) :: ndimf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf


    real(d), dimension(ndimf), intent(out) :: travec
    real(d), dimension(ndimf,ndim), intent(inout) :: arrd
    real(d), dimension(ndim), intent(in) :: autvec
   
    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: temp
    real(d), dimension(:,:), allocatable :: ar_offdiagd
    
    ndim1=kpq(1,0)
    
    allocate(ar_offdiagd(ndimf,ndim))
    
    call get_offdiag_tda_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,ar_offdiagd)


    arrd(:,:)=ar_offdiagd(:,:)
    
    
    deallocate(ar_offdiagd)

       call mat_vec_multiply(ndimf,ndim,arrd,autvec,travec)





!!! up to this point I have the vector  TRAVEC to make the scalar product with all other
!!! adc excited eigenstates, i.e. with the other eigenvectors of the hamiltonian
!!! adc matrix



!    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_tda_DIPOLE_direct_OK



!!$----------------------------------------------------------

  subroutine get_fspace_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,autvec,arrd,travec) 

    integer, intent(in) :: ndim
    integer, intent(in) :: ndimf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf


    real(d), dimension(ndimf), intent(out) :: travec
    real(d), dimension(ndimf,ndim), intent(inout) :: arrd
    real(d), dimension(ndim), intent(in) :: autvec
   
    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: temp
    real(d), dimension(:,:), allocatable :: ar_offdiagd
    
    ndim1=kpq(1,0)
    
    allocate(ar_offdiagd(ndimf,ndim))
    
!!!    call get_offdiag_adc2_DIPOLE_direct(ndim,kpq(:,:),ar_offdiagd(:,:))
    call get_offdiag_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,ar_offdiagd)


    arrd(:,:)=ar_offdiagd(:,:)
    
    
    deallocate(ar_offdiagd)

!!!    call dmatvec(ndim,autvec(:),arrd(:,:),travec(:))
!!!    call mat_vec_multiply_SYM(ndim,A,vec,fin)
       call mat_vec_multiply(ndimf,ndim,arrd,autvec,travec)





!!! up to this point I have the vector  TRAVEC to make the scalar product with all other
!!! adc excited eigenstates, i.e. with the other eigenvectors of the hamiltonian
!!! adc matrix



!    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2_DIPOLE_direct_OK



!!$----------------------------------------------------------

  subroutine get_fspace_adc2_DIPOLE_direct_second(ndim,ndimf,kpq,kpqf,autvec,arrd,travec) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    integer, intent(in) :: ndimf
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpqf


    real(d), dimension(ndimf), intent(out) :: travec
    real(d), dimension(ndimf,ndim), intent(inout) :: arrd
    real(d), dimension(ndim), intent(in) :: autvec
   
    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diagd,temp
    real(d), dimension(:,:), allocatable :: ar_offdiagd
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diagd(ndim),ar_offdiagd(ndimf,ndim))

     write(ilog,*) 'ndimf',ndimf
     write(ilog,*) 'ndim',ndim
    
!    call get_offdiag_adc2_DIPOLE_direct(ndim,kpq(:,:),ar_offdiagd(:,:))
    call get_offdiag_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,ar_offdiagd)


!    call get_diag_adc2_DIPOLE_direct(ndim1,ndim2,kpq(:,:),ar_diagd(:))
     
!    do i=1,ndim
!      do j=1,i-1
    
!     ar_offdiagd(i,j)=0.1
   
!      end do
!    end do   
    
    
!    arrd(:,:)=ar_offdiagd(:,:)
    
!     do i=1,ndim
         
!        ar_diagd(i)=20.0
   
!     end do
    

    do i=1,ndimf
       travec(i)=0.0
    end do
    
!    deallocate(ar_diagd,ar_offdiagd)

    do i=1,ndimf
      do j=1,ndim

!     arrd(i,j)=10.0
     arrd(i,j)=ar_offdiagd(i,j)
      end do
    end do   

    deallocate(ar_diagd,ar_offdiagd)



    do i=1,ndimf
     do j=1,ndim

     travec(i)=travec(i)+arrd(i,j)*autvec(j)

!     arrd(i,j)=ar_offdiag(i,j)
!     arrd(j,i)=arrd(i,j)
   
     end do
    end do   
   

    ! call dmatvec(ndim,autvec(:),arrd(:,:),travec(:))

    

!!! up to this point I have the vector  TRAVEC to make the scalar product with all other
!!! adc excited eigenstates, i.e. with the other eigenvectors of the hamiltonian
!!! adc matrix



!    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2_DIPOLE_direct_second

!!$----------------------------------------------------------------------
!!$----------------------------------------------------------------------


!!$----------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine write_fspace_adc2_DIPOLE(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer*8, intent(out) :: noffdel
    integer :: ndim1,ndim2,nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_DIPOLE_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_DIPOLE_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2_DIPOLE

!!$-----------------------------------------------------------------------------------------
!!$----------------------------------------------------------------------



!!!!!!!!!! END DIPOLE MATRIX ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!! HAMILTONIAN ROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!! ADC2 EXTENDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!$----------------------------------------------------------

  subroutine get_fspace_adc2e_direct(ndim,kpq,arr,evector) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diag,temp
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2ext_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_adc2ext_direct(ndim1,ndim2,kpq(:,:),ar_diag(:))
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)

    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2e_direct

!!$----------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine write_fspace_adc2e(ndim,kpq,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer*8 :: noffdel
    integer :: ndim1,ndim2,nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2ext_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2ext_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2e

!!$-----------------------------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine write_fspace_adc2e_1(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*8, intent(out) :: noffdel 

    integer :: ndim1,ndim2,nbuf,i

    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    write(ilog,*) 'in write_space'
    call get_offdiag_adc2ext_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2ext_save(ndim1,ndim2,kpq(:,:),nbuf,chr)

  end subroutine write_fspace_adc2e_1

!!$-----------------------------------------------------------------------------------------

    subroutine write_fspace_adc2e_1_cvs(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*8, intent(out) :: noffdel 

    integer :: ndim1,ndim2,nbuf,i

    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    write(ilog,*) 'in write_space'
    call get_offdiag_adc2ext_save_cvs(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2ext_save(ndim1,ndim2,kpq(:,:),nbuf,chr)

  end subroutine write_fspace_adc2e_1_cvs

!!$-----------------------------------------------------------------------------------------

  subroutine write_fspace_adc2e_1_MIO(ndim,kpq,noffdel,indx,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*8, intent(out) :: noffdel 
    INTEGER, DIMENSION(ndim), INTENT(IN) :: indx

    integer :: ndim1,ndim2,nbuf,i
    
       ndim1=kpq(1,0)
       ndim2=ndim-kpq(1,0)
       write(ilog,*) 'in write_space'
       call get_offdiag_adc2ext_save_MIO(ndim,kpq(:,:),nbuf,noffdel,indx,chr)
       call get_diag_adc2ext_save_MIO(ndim1,ndim2,kpq(:,:),nbuf,indx,chr)
!!$    end if

  end subroutine write_fspace_adc2e_1_MIO





!!!!!!!!!!!!!!!! ADC2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_fspace_adc2_direct(ndim,kpq,arr,evector)
    
    integer :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: ndim1, ndim2, nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_adc2_direct(ndim1,ndim2,kpq(:,:),ar_diag(:))
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
  end subroutine get_fspace_adc2_direct

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_adc2(ndim,kpq,chr) 

    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: ndim1,ndim2,nbuf,i
    integer*8 :: noffdel
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_adc2_1(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*8, intent(out) :: noffdel

    integer :: ndim1, ndim2, nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2_1

!#######################################################################

  subroutine write_fspace_adc2_1_cvs(ndim,kpq,noffdel,chr) 
    
    integer, intent(in)                                 :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in)                            :: chr
    integer*8, intent(out)                              :: noffdel

    integer :: ndim1, ndim2, nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save_cvs(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2_1_cvs

!#######################################################################

  subroutine write_fspace_adc2_1_MIO(ndim,kpq,noffdel,indx,chr) 

    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*8, intent(out) :: noffdel
    INTEGER, DIMENSION(ndim), intent(in) :: indx  

    integer :: ndim1, ndim2, nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save_MIO(ndim,kpq(:,:),nbuf,noffdel,indx,chr)
    call get_diag_adc2_save_MIO(ndim1,ndim2,kpq(:,:),nbuf,indx,chr)
    
  end subroutine write_fspace_adc2_1_MIO





!!!!!!!!!!!!!!!!!  ADC1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$------------------------------------------------------------------------

  subroutine get_fspace_tda_direct(ndim,kpq,arr,evector) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    call get_offdiag_tda_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_tda_direct(ndim,kpq(:,:),ar_diag(:))
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
  end subroutine get_fspace_tda_direct

!!$------------------------------------------------------------------------

  subroutine get_fspace_tda_direct_cvs(ndim,kpq,arr,evector) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    call get_offdiag_tda_direct_cvs(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_tda_direct_cvs(ndim,kpq(:,:),ar_diag(:))
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
  end subroutine get_fspace_tda_direct_cvs

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_tda(ndim,kpq,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: nbuf,i
    
    call get_offdiag_tda_save(ndim,kpq(:,:),nbuf,chr)
    call get_diag_tda_save(ndim,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_tda
  
!!$--------------------------------------------------------------------------



!!!!!!!!!!!!!!!!!!!! STATES ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$--------------------------------------------------------------------------

  subroutine get_fstates_direct(ndim,kpq,fspace,fen,nstate,chflag1,chflag2)
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nstate
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(3), intent(in) :: chflag1
    character(1), intent(in) :: chflag2
    real(d), dimension(ndim,ndim), intent(out) :: fspace
    real(d), dimension(ndim), intent(out) :: fen
    
    integer :: nryd,ifail
    integer, dimension(nBas) :: ncnfi_ryd
    
    integer, dimension(:), allocatable :: nisri
    real(d), dimension(:,:), allocatable :: arr
    real(d), dimension(:), allocatable :: evec,temp
    
    integer :: i,lim
!!$    integer*8 :: nstatel
!!$    integer*8, dimension(:), allocatable :: nisril
    
!    external M01CBF

    allocate(evec(ndim),arr(ndim,ndim),nisri(ndim),temp(ndim))
    if (chflag1 .eq. 'tda') then
       call get_fspace_tda_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    elseif (chflag1 .eq. 'ad2') then
       call get_fspace_adc2_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    elseif (chflag1 .eq. 'a2e') then
       call get_fspace_adc2e_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    end if
    
    if (chflag2 .eq. 'i') then
       call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
       lim=nryd
       temp(:)=0._d
       do i=1,lim
          temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
       end do
    elseif (chflag2 .eq. 'f') then
       lim=kpq(1,0)
       temp(:)=0._d
       do i=1,lim
          temp(:)=temp(:)+arr(i,:)**2
       end do
    end if
    
    if (chflag2 .eq. 'i') then
       call select_fstate_ryd(ndim,mspacewi,temp(:),nstate,nisri(:))
       write(ilog,*) nstate, "initial states has been selected"
    elseif (chflag2 .eq. 'f') then
       call select_fstate_ryd(ndim,mspacewf,temp(:),nstate,nisri(:))
       write(ilog,*) nstate, "final states has been selected"
    end if

    ifail=0
!    nstatel=int(nstate,lng)
!    nisril=int(nisri,lng)
!    call M01CBF(nisri(1:nstate),1,nstate,'A',ifail)

    if ((chflag2 .eq. 'i') .and. (numinista .gt. nstate) ) then
       write(ilog,*) "Requested number of the initial states is larger than the number of the available Rydberg states"
       stop
    end if

    do i=1,nstate
          write(ilog,*) nisri(i), evec(nisri(i))
          fspace(:,i)=arr(:,nisri(i))
          fen(i)=evec(nisri(i))
       end do
    
    deallocate(evec,arr,nisri,temp)
    
  end subroutine get_fstates_direct
!!$--------------------------------------------------------------

  subroutine load_fstates_weight(ndim,negvc,kpq,fspace,fen,nstate,name,chflag)

    integer, intent(in) :: ndim,negvc
    integer, intent(out):: nstate
    integer, dimension(7,0:nBas**2*nOcc**2),intent(in) :: kpq 
    real(d), dimension(negvc), intent(out) :: fen
    real(d), dimension(ndim,negvc), intent(out) :: fspace
    character(36), intent(in) :: name
    character(1), intent(in) :: chflag

    integer, dimension(nBas) :: ncnfi_ryd
    integer, dimension(:), allocatable :: nisri,isv,indx
    real(d), dimension(:), allocatable :: ener,temp,ener_ryd
    real(d), dimension(:,:), allocatable :: arr

    integer :: i,vectype,vecdim,nvec,nryd,lim

    allocate(ener(negvc),arr(ndim,negvc),nisri(negvc),temp(negvc),isv(negvc))
    
    if (chflag .eq. 'i') then
       call readvct(ndim,1,negvc,ener(:),arr(:,:),nvec)
    elseif (chflag .eq. 'f') then
       call readvct(ndim,2,negvc,ener(:),arr(:,:),nvec)
    end if

    if(nvec .ne. negvc) then
       write(ilog,*) 'The number of read vectors',nvec,'differs from the number of requested states',negvc
       stop
    end if

    if (chflag .eq. 'i') then
       call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
       lim=nryd
       temp(:)=0._d       
       do i=1,lim
          temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
       end do
    elseif (chflag .eq. 'f') then
       lim=kpq(1,0)
       temp(:)=0._d
       do i=1,lim
          temp(:)=temp(:)+arr(i,:)**2
       end do
    end if

    if (chflag .eq. 'i') then
       call select_fstate_ryd(negvc,mspacewi,temp(:),nstate,nisri(:))
    elseif (chflag .eq. 'f') then
       call select_fstate_ryd(negvc,mspacewf,temp(:),nstate,nisri(:))
    end if
    write(ilog,*) nstate, "states has been selected"

    
    allocate(ener_ryd(nstate),indx(nstate))
    
    do i=1,nstate
       ener_ryd(i)=ener(nisri(i))
    end do
    
    call dsortindxa1('A',nstate,ener_ryd(:),indx(:))
    
    if (chflag .eq. 'i') then
       if (numinista .gt. nstate) then 
          write(ilog,*) "Number of the set initial states is larger than the number of the available Rydberg states"
          stop
       end if
    
       do i= 1,numinista
          fspace(:,i)=arr(:,nisri(indx(i)))
          fen(i)=ener(nisri(indx(i)))
       end do
       
    end if

    do i= 1,nstate
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do

    deallocate(ener,arr,nisri,temp,ener_ryd,indx)

  end subroutine load_fstates_weight
  
!!$--------------------------------------------------------------------------------
  
  subroutine get_bound(ndim,negvc,kpq,fspace,fen,nstate)
    
    integer, intent(in) :: ndim,negvc
    integer, intent(out):: nstate
    integer, dimension(7,0:nBas**2*nOcc**2),intent(in) :: kpq 
    real(d), dimension(negvc), intent(out) :: fen
    real(d), dimension(ndim,negvc), intent(out) :: fspace

    integer, dimension(nBas) :: ncnfi_ryd
    integer, dimension(:), allocatable :: nisri,isv,indx
    real(d), dimension(:), allocatable :: ener,temp,ener_ryd
    real(d), dimension(:,:), allocatable :: arr

    integer :: i,vectype,vecdim,nvec,nryd,lim

    allocate(ener(negvc),arr(ndim,negvc),nisri(negvc),temp(negvc),isv(negvc))
    
    !Reading ALL Davidson vectors 
    call readvct(ndim,1,negvc,ener(:),arr(:,:),nvec)

    call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
    lim=nryd
    temp(:)=0._d       
    do i=1,lim
       temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
    end do

    call select_fstate_ryd(negvc,mspacewi,temp(:),nstate,nisri(:))

    write(ilog,*) nstate, "states has been selected"

    
    allocate(ener_ryd(nstate),indx(nstate))
    
    do i=1,nstate
       ener_ryd(i)=ener(nisri(i))
    end do
    
    
    call dsortindxa1('A',nstate,ener_ryd(:),indx(:))
    
    if (numinista .gt. nstate) then 
       write(ilog,*) "Number of the set initial states is larger than the number of the available Rydberg states"
       stop
    end if
    
    do i= 1,numinista
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do
       
    do i= 1,nstate
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do

    deallocate(ener,arr,nisri,temp,ener_ryd,indx)

  end subroutine get_bound
    
!!$--------------------------------------------------------------




!!!!!!!!!!!!!!  GET TRANSITION MOMENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_tranmom_1(ndim,negvc,name,mtm,nstate,fen,tmvec,ndims)
    
    integer, intent(in)                    :: ndim,negvc,ndims
    integer, intent(out)                   :: nstate
    real(d), dimension(ndim), intent(in)   :: mtm 
    real(d), dimension(negvc), intent(out) :: fen,tmvec
    
    character(36), intent(in)              :: name

    integer                                :: i,j,num,k
    real(d)                                :: enr,cntr
    real(d), dimension(:), allocatable     :: vec 
    logical                                :: log1

    INQUIRE(file=name,exist=log1)
    
    if(.not. log1) then
       write(ilog,*) 'The file ', name,' does not exist'
       stop
    end if
    
    allocate(vec(ndim))
    nstate=0
    cntr=0.0
    
    if(dlim.lt.0.0) then

       write(ilog,*) 'Take all states with 2h2p part greater than',-dlim

       OPEN(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
       do i=1,negvc

          read(78,END=78) num,enr,vec(:)

          do j=1,ndims
             cntr=cntr+vec(j)**2
          end do

          if(cntr.lt.-dlim) then
             nstate=nstate+1
             fen(nstate)=enr
             tmvec(nstate)=tm(ndim,vec(:),mtm(:))
          end if

          cntr=0.0
       end do

78     CLOSE(78)
       write(ilog,*) 'There are ',nstate,' energies and tran. moments'


    else

       write(ilog,*) 'Take all states with 1h1p part greater than',dlim
    
       OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')

       do i=1,negvc

          read(77,END=77) num,enr,vec(:)

          do j=1,ndims
             cntr=cntr+vec(j)**2
          end do

          if(cntr.gt.dlim) then
             nstate=nstate+1
             fen(nstate)=enr
             tmvec(nstate)=tm(ndim,vec(:),mtm(:))
          end if

          cntr=0.0
       end do

77     CLOSE(77)
       write(ilog,*) 'There are ',nstate,' energies and tran. moments'
    
    end if
    
       
  end subroutine get_tranmom_1

!!$---------------------------------------------------------



!  subroutine get_tranmom_2(ndim,kpq,negvc,name,arrd,nstate,fen,travec,ndims)
    
!    integer, intent(in) :: ndim,negvc,ndims
!    integer, intent(out) :: nstate
!    real(d), dimension(ndim,ndim), intent(in) :: arrd
!    real(d), dimension(negvc), intent(out) :: travec
!    real(d), intent(out) :: fen   
!    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    
!    character(36), intent(in) :: name

!    integer :: i,j,num
!    real(d) :: enr,cntr
!    real(d), dimension(:), allocatable :: vec 
!    logical :: log1

!    INQUIRE(file=name,exist=log1)
    
!    if(.not. log1) then
!       write(ilog,*) 'The file ', name,' does not exist'
!       stop
!    end if
    
!    allocate(vec(ndim))
!    nstate=0
!    cntr=0.0

!    if(dlim.lt.0.0) then

!    write(ilog,*) 'Take all states with 2h2p part greater than',-dlim

!    OPEN(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc

!       read(78,END=78) num,enr,vec(:)

!       write(ilog,*) 'Read vector',num,enr

!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do


!      if(i.eq.statenumber) then

!       if(cntr.lt.-dlim) then
!        write(ilog,*) 'State',num,'taken'
!        fen=enr
!        call   get_fspace_adc2_DIPOLE_direct(ndim,kpq,vec(:),arrd(:,:),travec(:))
!        travec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if
!      end if

!       cntr=0.0
!    end do
!78  CLOSE(78)
!    write(ilog,*) 'There is the ', statenumber,' energy and tran. moment'


!    else

!    write(ilog,*) 'Take all states with 1h1p part greater than',dlim
    
!    OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc

!       read(77,END=77) num,enr,vec(:)

!       write(ilog,*) 'Read vector',num,enr


!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do



!      if(i.eq.statenumber) then


!       if(cntr.gt.dlim) then
!        write(ilog,*) 'State',num,'taken'
!        nstate=nstate+1   
!        fen=enr
!        call   get_fspace_adc2_DIPOLE_direct(ndim,kpq,vec(:),arrd(:,:),travec(:))
!        travec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if
!      end if

!       cntr=0.0
!    end do
!77  CLOSE(77)
!    write(ilog,*) 'There is the ',statenumber,' energy and tran. moment'
    
!    end if 
    
       
!  end subroutine get_tranmom_2






!  subroutine get_tranmom_3(ndim,negvc,name,mtm,nstate,fen,tmvec,ndims)
    
!    integer, intent(in) :: ndim,negvc,ndims
!    integer, intent(out) :: nstate
!    real(d), dimension(ndim), intent(in) :: mtm 
!    real(d), dimension(negvc-statenumber), intent(out) :: fen,tmvec
    
!    character(36), intent(in) :: name

!    integer :: i,j,num
!    real(d) :: enr,cntr
!    real(d), dimension(:), allocatable :: vec 
!    logical :: log1

!    INQUIRE(file=name,exist=log1)
    
!    if(.not. log1) then
!       write(ilog,*) 'The file ', name,' does not exist'
!       stop
!    end if
    
!    allocate(vec(ndim))
!    nstate=0
!    cntr=0.0

!    if(dlim.lt.0.0) then

!    write(ilog,*) 'Take all states with 2h2p part greater than',-dlim

!    OPEN(unit=78,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc

!       read(78,END=78) num,enr,vec(:)

!       write(ilog,*) 'Read vector',num,enr
!      if (i.gt.statenumber) then
     

!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do

!       if(cntr.lt.-dlim) then
!        write(ilog,*) 'State',num,'taken'
!        nstate=nstate+1
!        fen(nstate)=enr
!        tmvec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if

!    end if

!       cntr=0.0
!    end do
!78  CLOSE(78)
!    write(ilog,*) 'There are ',nstate,' energies and tran. moments'


!    else

!    write(ilog,*) 'Take all states with 1h1p part greater than',dlim
    
!    OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
!    do i=1,negvc

!       read(77,END=77) num,enr,vec(:)

!       write(ilog,*) 'Read vector',num,enr

!      if (i.gt.statenumber) then
 
!        do j=1,ndims
!          cntr=cntr+vec(j)**2
!        end do

!       if(cntr.gt.dlim) then
!        write(ilog,*) 'State',num,'taken'
!        nstate=nstate+1   
!        fen(nstate)=enr
!        tmvec(nstate)=tm(ndim,vec(:),mtm(:))
!       end if


!     end if


!       cntr=0.0
!    end do
!77  CLOSE(77)
!    write(ilog,*) 'There are ',nstate,' energies and tran. moments'
    
!    end if 
    
       
!  end subroutine get_tranmom_3

!!$---------------------------------------------------------
!!$ N.B., here ndim is the no. states
!!$---------------------------------------------------------

  subroutine get_sigma(ndim,ener,sigmavec)

    integer, intent(in)                  :: ndim
    real(d), dimension(ndim), intent(in) :: ener
    real(d), dimension(ndim), intent(in) :: sigmavec
    
    integer                              :: i,ncount,iout
    real(d)                              :: oslimit
    real(d), dimension(:), allocatable   :: sgmvc_short,ener_short
    
    allocate(sgmvc_short(ndim),ener_short(ndim))

    oslimit=1.e-8
    ncount=0

! sgmvc_short: array of oscillator strengths (in Mb) that are greater than oslimit
! ener_short:  array of corresponding state energies
    do i= 1,ndim
       if (sigmavec(i) .gt. oslimit) then
          ncount=ncount+1
          sgmvc_short(ncount)=sigmavec(i)
          ener_short(ncount)=ener(i)
       end if
    end do

    call get_sums(ncount,ener_short(1:ncount),sgmvc_short(1:ncount))
  
    iout=112
    open(iout,file='osc.dat',form='formatted',status='unknown')
    do i=1,ncount
       write(iout,'(2x,F20.15,2x,E21.15)') ener_short(i),sgmvc_short(i)
    enddo
    close(iout)

    deallocate(sgmvc_short,ener_short)

  end subroutine get_sigma

!!$--------------------------------------------------------
!!$ N.B., here ndim is the no. states
!!$---------------------------------------------------------

  subroutine get_sums(ndim,ener,fosc)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: ener, fosc
    
    real(d), dimension(0:50):: sums
    integer :: i,j,k
    real(d) :: elev,flev,ratio

    character(len=8) :: atmp

!-----------------------------------------------------------------------
! Calculate the negative spectral moments
!-----------------------------------------------------------------------
    sums=0.0d0
    do i=1,ndim
       elev=ener(i)
       flev=fosc(i)
       ratio=flev
       do j=0,50
          sums(j)=sums(j)+ratio
          ratio=ratio/elev
       enddo
    enddo

!-----------------------------------------------------------------------
! Write the negative spectral moments to the log file
!-----------------------------------------------------------------------
    write(ilog,'(70a)') ('-',i=1,70)
    write(ilog,*) 'Negative spectral moments'
    write(ilog,'(70a)') ('-',i=1,70)
    do i=0,50
       atmp='S(-'
       if (i.eq.0) then
          atmp='S(0)   ='
       else if (i.lt.10) then
          write(atmp(4:8),'(i1,a4)') i,')  ='
       else
          write(atmp(4:8),'(i2,a3)') i,') ='
       endif
       write(ilog,'(a8,x,E14.7)') atmp,sums(i)
    enddo
   
  end subroutine get_sums
!!$---------------------------------------------------------

  subroutine fill_stvc(ndim,vctr)
    
    integer, intent(in) :: ndim
    real(d),dimension(ndim), intent(in) :: vctr
    
    integer :: i,cnt
    integer, dimension(:), allocatable :: indarr
    
    allocate(indarr(ndim))
    call dsortindxa1('D',ndim,vctr(:)**2,indarr(:))

    stvc_lbl(1:lmain)=indarr(1:lmain)

    deallocate(indarr)
    
  end subroutine fill_stvc

!----------------------------------------------------------------------

 subroutine test_ortho(ndim,nstate,arr)

    integer, intent(in) :: ndim,nstate
    real(d), dimension(nstate,ndim), intent(in) :: arr
    real(d), dimension(:,:), allocatable :: mat
    real(d), dimension(:), allocatable :: tempvec1,tempvec2,evec
    integer :: i,j
    real(d) :: entri

    write(ilog,*) 'Starting Ortho Test'

    allocate(mat(ndim,ndim))
    allocate(tempvec1(nstate))
    allocate(tempvec2(nstate))

    do i=1,ndim
     do j=1,ndim
      tempvec1(:)=arr(:,i)
      tempvec2(:)=arr(:,j)
      entri=dsp(nstate,tempvec1(:),tempvec2(:))
      mat(i,j)=entri
     end do
    end do

    deallocate(tempvec1,tempvec2)

    do i=1,ndim
     write(ilog,'(99(f3.2,1X))') (mat(i,j),j=1,ndim)
    end do

    allocate(evec(ndim))

    call vdiagonalise(ndim,mat,evec)

    do i=1,ndim
     write(ilog,'(A,I3,2X,f4.3)') 'OM.EVal',i,evec(i)
    end do

    deallocate(mat,evec)

 end subroutine test_ortho

 subroutine show_vecs(ndim,nstate,arr)
! ndim=Anzahl Vecs, nstate=Länge Vecs, arr=Vecs
    integer, intent(in) :: ndim,nstate
    real(d), dimension(nstate,ndim), intent(in) :: arr
    integer :: i,j

    write(ilog,*) 'Starting Show_Vecs'

     do j=1,nstate
      write(ilog,'(I5,99(1X,f6.3))') j,(arr(j,i),i=1,15)
     end do

 end subroutine show_vecs






  subroutine get_fspace_adc2_direct_MIO(ndim,kpq,arr,evector,indx) 
    
    integer :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr
    INTEGER, DIMENSION(ndim), intent(in) :: indx  

    integer :: ndim1, ndim2, nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2_direct_MIO(ndim,kpq(:,:),ar_offdiag(:,:),indx)
    call get_diag_adc2_direct_MIO(ndim1,ndim2,kpq(:,:),ar_diag(:),indx)
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
  end subroutine get_fspace_adc2_direct_MIO




  subroutine get_fspace_adc2e_direct_MIO(ndim,kpq,arr,evector,indx) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr
    INTEGER, DIMENSION(ndim), intent(in) :: indx  

    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diag,temp
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2ext_direct_MIO(ndim,kpq(:,:),ar_offdiag(:,:),indx)
    call get_diag_adc2ext_direct_MIO(ndim1,ndim2,kpq(:,:),ar_diag(:),indx)
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)

    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2e_direct_MIO











end module fspace
