  module block_lanczos

    use constants
    use parameters
    use channels
    
    implicit none

    real(d), dimension(:,:), allocatable :: qmat1,qmat2,umat,rmat,&
                                            amat,bmat
    
  contains

!#######################################################################

    subroutine lancdiag_block(matdim,noff,flag)

      implicit none

      integer, intent(in)     :: matdim
      integer*8, intent(in)   :: noff
      character(1),intent(in) :: flag

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
      allocate(qmat1(matdim,lmain))
      allocate(qmat2(matdim,lmain))
      allocate(umat(matdim,lmain))
      allocate(rmat(matdim,lmain))
      allocate(amat(lmain,lmain))
      allocate(bmat(lmain,lmain))

      qmat1=0.0d0
      qmat2=0.0d0
      umat=0.0d0
      rmat=0.0d0
      amat=0.0d0
      bmat=0.0d0
      
!-----------------------------------------------------------------------
! Set the initial Lanczos vectors
!-----------------------------------------------------------------------
      call init_vec

!-----------------------------------------------------------------------
! Perform the block Lanczos calculation
!-----------------------------------------------------------------------
      call run_block_lanczos(matdim)
      
      STOP
      
      return
      
    end subroutine lancdiag_block
    
!#######################################################################
    
    subroutine init_vec
      
      use iomod, only: freeunit
      
      implicit none
      
      integer                              :: i,j,k
      integer                              :: iadc1,idim
      real(d), dimension(:,:), allocatable :: adc1vec
      real(d)                              :: fac

!-----------------------------------------------------------------------
! lancguess = 1 <-> Construction the initial vectors as the IS unit
!                   vectors with the greatest transition dipoles with 
!                   the initial state
!
!             2 <-> Construction of the initial vectors from the
!                   ADC(1) eigenvectors with the greatest transition
!                   dipoles with the initial state
!
!             3 <-> Construction of the initial vectors from linear
!                   combinations of the most important 1h1p and 2h2p 
!                   IS unit vectors
!
!             4 <-> Construction of the initial vectors from linear
!                   combinations of the most important ADC(1)
!                   eigenvectors and 2h2p IS unit vectors
!-----------------------------------------------------------------------

      if (lancguess.eq.1) then

         ! Copy the 1h1p of interest into the rmat array
         
         do i=1,lmain
            k=stvc_lbl(i)
            rmat(k,i)=1.0d0
         enddo

      else if (lancguess.eq.2) then
         
         ! Read the ADC(1) eigenvectors from file
         call freeunit(iadc1)
         
         open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
              status='old')
         
         read(iadc1) idim

         allocate(adc1vec(idim,idim))
         
         rewind(iadc1)

         read(iadc1) idim,adc1vec
         
         close(iadc1)

         ! Copy the ADC(1) vectors of interest into the rmat array
         do i=1,lmain
            k=stvc_lbl(i)
            do j=1,idim
               rmat(j,i)=adc1vec(j,k)
            enddo
         enddo

      else if (lancguess.eq.3) then

         ! Copy the linear combinations of the 1h1p and 2h2p ISs into
         ! the rmat array
         fac=1.0d0/sqrt(2.0d0)
         do i=1,lmain
            k=stvc_mxc(i*3-1)
            rmat(k,i)=fac            
            k=stvc_mxc(i*3)
            if (stvc_mxc(i*3-2).gt.0) then  
               rmat(k,i)=fac
            else
               rmat(k,i)=-fac
            endif
         enddo
      
      else if (lancguess.eq.4) then

         ! Read the ADC(1) eigenvectors from file
         call freeunit(iadc1)
         
         open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
              status='old')
         
         read(iadc1) idim

         allocate(adc1vec(idim,idim))
         
         rewind(iadc1)

         read(iadc1) idim,adc1vec
         
         close(iadc1)

         ! Copy the linear combinations of the ADC(1) vectors and 2h2p
         ! ISs into the rmat array
         fac=1.0d0/sqrt(2.0d0)
         do i=1,lmain

            k=stvc_mxc(i*3-1)

            do j=1,idim
               rmat(j,i)=fac*adc1vec(j,k)
            enddo

            k=stvc_mxc(i*3)
            if (stvc_mxc(i*3-2).gt.0) then  
               rmat(k,i)=fac
            else
               rmat(k,i)=-fac
            endif

         enddo
         
      endif
      
      return
      
    end subroutine init_vec
      
!#######################################################################

    subroutine run_block_lanczos(matdim)

      implicit none

      integer, intent(in) :: matdim
      integer             :: j,k

      do j=1,ncycles

!-----------------------------------------------------------------------
! Output progress
!-----------------------------------------------------------------------         
         write(ilog,'(70a)') ('*',k=1,70)
         write(ilog,'(2x,a,2x,i4)') 'Iteration number',j

!-----------------------------------------------------------------------
! (1.) Calculate the current block of on-diagonal elements of the
!      T-matrix
!-----------------------------------------------------------------------
         call hxq(matdim)
         umat=umat-matmul(qmat1,transpose(bmat))
         amat=matmul(transpose(qmat2),umat)
         
!-----------------------------------------------------------------------
! (2.) Calculate the next block of Krylov vectors
!-----------------------------------------------------------------------
         rmat=umat-matmul(qmat2,amat)

!-----------------------------------------------------------------------
! (3.) Compute the QR factorization of the matrix of Krylov vectors
!-----------------------------------------------------------------------
         
      enddo
      
      return

    end subroutine run_block_lanczos

!#######################################################################

    subroutine hxq(matdim)

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: unit
      integer                            :: maxbl,nrec,nlim,i,j,k,l,m,n
      integer, dimension(:), allocatable :: indxi,indxj
      real(d), dimension(matdim)         :: hii
      real(d), dimension(:), allocatable :: hij

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      unit=77
      open(unit,file='SCRATCH/hmlt.diac',status='old',access='sequential',&
           form='unformatted')

      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

      umat=0.0d0
      do m=1,matdim
         do n=1,lmain
            umat(m,n)=hii(m)*qmat2(m,n)
         enddo
      enddo

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
      
      open(unit,file='SCRATCH/hmlt.offc',status='old',access='sequential',&
           form='unformatted')

      do k=1,nrec
         read(unit) hij(:),indxi(:),indxj(:),nlim
         do l=1,nlim
            i=indxi(l)
            j=indxj(l)
            umat(i,:)=umat(i,:)+hij(l)*qmat2(j,:)
            umat(j,:)=umat(j,:)+hij(l)*qmat2(i,:)
         enddo
      enddo

      close(unit)

      deallocate(hij,indxi,indxj)
      
      return
      
    end subroutine hxq
      
!#######################################################################

  end module block_lanczos
