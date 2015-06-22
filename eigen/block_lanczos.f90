  module block_lanczos

    use constants
    use parameters
    use channels
    
    implicit none

    integer                              :: nvec
    real(d), dimension(:,:), allocatable :: qmat1,qmat2,umat,rmat,&
                                            amat,bmat,tmat
    
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
      
      nvec=lmain*ncycles
      allocate(tmat(nvec,nvec))

      qmat1=0.0d0
      qmat2=0.0d0
      umat=0.0d0
      rmat=0.0d0
      amat=0.0d0
      bmat=0.0d0
      tmat=0.0d0

!-----------------------------------------------------------------------
! Set the initial Lanczos vectors
!-----------------------------------------------------------------------
      call init_vec

!-----------------------------------------------------------------------
! Perform the block Lanczos calculation
!-----------------------------------------------------------------------
      call run_block_lanczos(matdim)

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

         ! Copy the 1h1p of interest into the qmat2 array
         
         do i=1,lmain
            k=stvc_lbl(i)
            qmat2(k,i)=1.0d0
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

         ! Copy the ADC(1) vectors of interest into the qmat2 array
         do i=1,lmain
            k=stvc_lbl(i)
            do j=1,idim
               qmat2(j,i)=adc1vec(j,k)
            enddo
         enddo

      else if (lancguess.eq.3) then

         ! Copy the linear combinations of the 1h1p and 2h2p ISs into
         ! the qmat2 array
         fac=1.0d0/sqrt(2.0d0)
         do i=1,lmain
            k=stvc_mxc(i*3-1)
            qmat2(k,i)=fac            
            k=stvc_mxc(i*3)
            if (stvc_mxc(i*3-2).gt.0) then  
               qmat2(k,i)=fac
            else
               qmat2(k,i)=-fac
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
         ! ISs into the qmat2 array
         fac=1.0d0/sqrt(2.0d0)
         do i=1,lmain

            k=stvc_mxc(i*3-1)

            do j=1,idim
               qmat2(j,i)=fac*adc1vec(j,k)
            enddo

            k=stvc_mxc(i*3)
            if (stvc_mxc(i*3-2).gt.0) then  
               qmat2(k,i)=fac
            else
               qmat2(k,i)=-fac
            endif

         enddo
         
      endif
      
      return
      
    end subroutine init_vec
      
!#######################################################################

    subroutine run_block_lanczos(matdim)

      use iomod, only: freeunit

      implicit none

      integer, intent(in)       :: matdim
      integer                   :: lanunit,j,i,k,i1,j1,k1,k2,m,n,upper
      integer                   :: info
      real(d), dimension(lmain) :: tau
      real(d), dimension(lmain) :: work
      real(d)                   :: t1,t2

      write(ilog,'(/,70a)') ('*',i=1,70)
      write(ilog,'(12x,a)') &
           'Block-Lanczos Diagonalisation in the Final Space'
      write(ilog,'(70a,/)') ('*',i=1,70)

      call cpu_time(t1)

      call freeunit(lanunit)
      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',&
              status='unknown')

      do j=1,ncycles

!-----------------------------------------------------------------------
! Output progress
!-----------------------------------------------------------------------         
         write(ilog,'(70a)') ('*',k=1,70)
         write(ilog,'(2x,a,2x,i4)') 'Iteration number',j

!-----------------------------------------------------------------------
! Write the latest block of Lanczos vectors to disk
!-----------------------------------------------------------------------
         write(lanunit) qmat2(:,:)

!-----------------------------------------------------------------------
! Calculate the current block of on-diagonal elements of the T-matrix
!-----------------------------------------------------------------------
         call hxq(matdim)

         umat=umat-matmul(qmat1,transpose(bmat))

         amat=matmul(transpose(qmat2),umat)
         
!-----------------------------------------------------------------------
! Calculate the next block of Krylov vectors
!-----------------------------------------------------------------------
         rmat=umat-matmul(qmat2,amat)

!-----------------------------------------------------------------------
! Compute the QR factorization of the matrix of Krylov vectors
!-----------------------------------------------------------------------
         ! Compute the current block of off-diagonal elements of
         ! the T-matrix
         call dgeqrf(matdim,lmain,rmat,matdim,tau,work,lmain,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'dqerf failed in run_block_lanczos'
            STOP
         endif

         ! Note that the B-matrix is upper-triangular
         bmat=0.0d0
         do k=1,lmain
            do i=1,k
               bmat(i,k)=rmat(i,k)
            enddo
         enddo

         ! Extract the next block of Lanczos vectors
         call dorgqr(matdim,lmain,lmain,rmat,matdim,tau,work,lmain,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'dorgqr failed in run_block_lanczos'
            STOP
         endif

         ! Update the matrices of Lanczos vectors
         qmat1=qmat2
         qmat2=rmat

!-----------------------------------------------------------------------
! Fill in the next block of the T-matrix array
!-----------------------------------------------------------------------
         i1=0 ! Initialise the column counter
         do m=(j-1)*lmain+1,j*lmain ! Loop over columns of T_j
            i1=i1+1 ! Increment the column counter
            if (j.lt.ncycles) then
               upper=(j+1)*lmain
            else
               upper=j*lmain
            endif
            j1=0 ! Initialise the row counters
            k1=0
            k2=0
            do n=(j-1)*lmain+1,upper ! Loop over rows of T_j
               j1=j1+1 ! Increment the main row counter
               if (j1.le.lmain) then ! Contribution from A_j
                  k1=k1+1
                  tmat(n,m)=amat(k1,i1)
               else                  ! Contribution from B_j
                  k2=k2+1
                  tmat(n,m)=bmat(k2,i1)
               endif
               tmat(m,n)=tmat(n,m)
            enddo
         enddo

      enddo

      lancstates=ncycles*lmain
      
      close(lanunit)

      call cpu_time(t2)
      write(ilog,'(/,2x,a,1x,F8.2,1x,a1,/)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! Calculate the Lanczos pseudospectrum
!-----------------------------------------------------------------------
      call calc_pseudospec(lanunit,matdim)

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

    subroutine calc_pseudospec(lanunit,matdim)

      implicit none

      integer                       :: lanunit,matdim,i
      real(d), dimension(nvec,nvec) :: umat
      real(d), dimension(nvec)      :: eigval
      real(d)                       :: t1,t2,mem

!-----------------------------------------------------------------------
! (1) Lanczos state energies
!-----------------------------------------------------------------------
      write(ilog,'(70a)') ('-',i=1,70)
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state energies...'

      call cpu_time(t1)

      call diagmat_banded(lmain,tmat,nvec,umat,eigval)

      call cpu_time(t2)

      write(ilog,'(2x,a,1x,F8.2,1x,a1)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! (2) Lanczos state vectors
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state vectors...'

      call cpu_time(t1)

      mem=16.0d0*matdim*nvec/1024.0d0**2

      if (mem.le.lancmem) then
         write(ilog,'(2x,a,/)') 'Calculation will proceed in-core'
         call ritzvecs_incore(lanunit,umat,eigval,nvec,matdim)
      else
         write(ilog,'(2x,a,/)') 'Calculation will proceed out-of-core'
         call ritzvecs_ext(lanunit,umat,nvec,matdim,eigval,lmain)
      endif

      call cpu_time(t2)
      write(ilog,'(2x,a,1x,F8.2,1x,a1,/)') 'Time taken:',t2-t1,'s'

!      call chkortho(nvec,matdim)

      return

    end subroutine calc_pseudospec

!#######################################################################

    subroutine diagmat_banded(blckdim,matrix,dim,eigvec,eigval)

      implicit none
      
      integer                           :: blckdim,dim,i,j,iupper
      real(d), dimension(dim,dim)       :: matrix

      integer                           :: kd,ldab,error
      real(d), dimension(blckdim+1,dim) :: ab
      real(d), dimension(dim)           :: eigval
      real(d), dimension(dim,dim)       :: eigvec
      real(d), dimension(3*dim-2)       :: work

!-----------------------------------------------------------------------
! Set dimensions required to be passed to dsbev
!-----------------------------------------------------------------------
      kd=blckdim
      ldab=blckdim+1

!-----------------------------------------------------------------------
! Fill in the array ab holding the upper triangle of the projection of
! the Hamiltonian onto the space spanned by the Lanczos vectors
!-----------------------------------------------------------------------
      ab=0.0d0
      do j=1,dim
         iupper=min(dim,j+kd)
         do i=j,iupper
            ab(1+i-j,j)=matrix(i,j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Diagonalise the projection of the Hamiltonian onto the space spanned
! by the Lanczos vectors
!-----------------------------------------------------------------------
      call dsbev('V','L',dim,kd,ab,ldab,eigval,eigvec,dim,work,error)
        
      if (error.ne.0) then
         write(ilog,'(/,2x,3a,/)') 'Diagonalisation of the Lanczos ',&
              'representation of the Hamiltonian failed in ',&
              'subroutine diagmat_banded.'
         STOP
      endif

      return

    end subroutine diagmat_banded

!#######################################################################

    subroutine ritzvecs_ext(lanunit,umat,dim,matdim,eigval,iblckdim)

      implicit none

      integer                              :: v
      integer                              :: i,j,k,m,n,count1,lanunit,&
                                              dim,matdim,ritzunit,&
                                              iblckdim,nblcks,last,count2
      integer                              :: blcksize
      real(d), dimension(dim,dim)          :: umat
      real(d), dimension(dim)              :: eigval
      real(d), dimension(matdim,lmain)     :: lvec
      real(d), dimension(:,:), allocatable :: ritzvec
      real(d)                              :: maxmem

!-----------------------------------------------------------------------
! Determine the maximum number of Ritz vectors that we can compute
! per sweep
!-----------------------------------------------------------------------
      maxmem=250.0d0

!      blcksize=int(floor((maxmem*1024.0d0**2)/(8.0d0*matdim)))
!
!      if (blcksize.gt.dim) blcksize=dim

      blcksize=iblckdim

      allocate(ritzvec(matdim,blcksize))

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',status='old')

      ritzunit=lanunit+1
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate the Ritz vectors
!-----------------------------------------------------------------------
      ! Loop over blocks of Ritz vectors
      nblcks=ceiling(real(dim)/real(blcksize))

      do i=1,nblcks-1
         
         ritzvec=0.0d0

         ! Calculate the curent block of blcksize Ritz vectors
         rewind(lanunit)
         count2=0
         do k=1,ncycles ! Loop over blocks of Lanczos vectors
            read(lanunit) lvec
            ! Loop over each Lanczos vector in the current block
            ! Each Lanczos vector contributes to all components
            ! of the current block of Ritz vectors
            do n=1,lmain
               count1=0
               count2=count2+1
               do v=blcksize*i-blcksize+1,blcksize*i
                  count1=count1+1
                  do m=1,matdim
                     ritzvec(m,count1)=ritzvec(m,count1)+lvec(m,n)*umat(count2,v)
                  enddo
               enddo
            enddo
         enddo

         ! Write the current block of Ritz vectors to file along with
         ! the corresponding Ritz values
         count1=0
         do v=blcksize*i-blcksize+1,blcksize*i
            count1=count1+1
            write(ritzunit) v,eigval(v),ritzvec(:,count1)
         enddo
         
      enddo

      ! Remaining n Ritz vectors from the (potentially) incomplete 
      ! block
      n=dim-blcksize*(nblcks-1)
      ritzvec=0.0d0
      rewind(lanunit)
      count2=0
      do k=1,ncycles
         read(lanunit) lvec
         do i=1,lmain
            count1=0
            count2=count2+1
            do v=dim-n+1,dim
               count1=count1+1
               do m=1,matdim
                  ritzvec(m,count1)=ritzvec(m,count1)+lvec(m,i)*umat(count2,v)
               enddo
            enddo
         enddo
      enddo
         
      count1=0
      do v=dim-n+1,dim
         count1=count1+1
         write(ritzunit) v,eigval(v),ritzvec(:,count1)
      enddo

!-----------------------------------------------------------------------
! Close the Lanczos and Ritz vector files
!-----------------------------------------------------------------------
      close(lanunit)
      close(ritzunit)

      return

    end subroutine ritzvecs_ext

!#######################################################################

    subroutine ritzvecs_incore(lanunit,umat,eigval,dim,matdim)

      implicit none

      integer                          :: lanunit,ritzunit,dim,matdim,&
                                          i,j,k
      integer                          :: v
      real(d), dimension(matdim,dim)   :: lvec,rvec
      real(d), dimension(matdim,lmain) :: tmpvec
      real(d), dimension(dim,dim)      :: umat
      real(d), dimension(dim)          :: eigval

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',status='old')

      ritzunit=lanunit+1
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')    

!-----------------------------------------------------------------------
! Read the Lanczos vectors from file
!-----------------------------------------------------------------------
      k=0
      do i=1,ncycles
         read(lanunit) tmpvec
         do j=1,lmain
            k=k+1
            lvec(:,k)=tmpvec(:,j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate and output the Ritz vectors
!-----------------------------------------------------------------------
      call dgemm('N','N',matdim,dim,dim,1.0d0,lvec,matdim,umat,dim,&
           0.0d0,rvec,matdim)

      do v=1,dim
         write(ritzunit) v,eigval(v),rvec(:,v)
      enddo
      
!-----------------------------------------------------------------------
! Close the Lanczos and Ritz vector files
!-----------------------------------------------------------------------
      close(lanunit)
      close(ritzunit)

      return

    end subroutine ritzvecs_incore

!#######################################################################

    subroutine chkortho(dim,matdim)

      implicit none

      integer                          :: unit,dim,matdim,i,j,k
      real(d), dimension(matdim,dim)   :: lvec
      real(d), dimension(matdim,lmain) :: tmpvec
      real(d)                          :: dp,tmp
      real(d), parameter               :: tol=1d-4

      lvec=0.0d0

      print*,"Ritz vectors"

      unit=28
      open(unit,file=lancname,form='unformatted',status='old')

      do i=1,nvec
         read(unit) k,tmp,lvec(:,i)
      enddo

      do i=1,dim-1
         do j=i+1,dim
            dp=dot_product(lvec(:,i),lvec(:,j))
            if (abs(dp).gt.tol) print*,i,j,dp
         enddo
      enddo

      close(unit)


      print*,"Lanczos vectors"
      open(unit,file='SCRATCH/lanvecs',form='unformatted',status='old')
      
      k=0
      do i=1,ncycles
         read(unit) tmpvec(:,:)
         do j=1,lmain
            k=k+1
            lvec(:,k)=tmpvec(:,j)
         enddo
      enddo

      do i=1,dim-1
         do j=i+1,dim
            dp=dot_product(lvec(:,i),lvec(:,j))
            if (abs(dp).gt.tol) print*,i,j,dp
         enddo
      enddo

      close(unit)

      return

    end subroutine chkortho

!#######################################################################

  end module block_lanczos
