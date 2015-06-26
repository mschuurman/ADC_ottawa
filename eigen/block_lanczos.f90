  module block_lanczos

    use constants
    use parameters
    use channels
    
    implicit none

    integer                              :: nvec,nwr,buffsize
    real(d), dimension(:,:), allocatable :: qmat1,qmat2,umat,rmat,&
                                            amat,bmat,tmat,tmpmat
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
      allocate(tmpmat(matdim,lmain))
      
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
      call run_block_lanczos(matdim,noff)

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

    subroutine run_block_lanczos(matdim,noff)

      use iomod, only: freeunit

      implicit none

      integer, intent(in)                  :: matdim
      integer*8                            :: noff
      integer                              :: lanunit,j,i,k,i1,j1,k1,k2,&
                                              m,n,upper,reclength,&
                                              nv,nsurplus,nprev,&
                                              nthreads
      integer                              :: info
      integer, dimension(:), allocatable   :: indxi,indxj
      real(d), dimension(:), allocatable   :: hii,hij
      real(d), dimension(lmain)            :: tau
      real(d), dimension(lmain)            :: work
      real(d)                              :: t1,t2,mem
      real(d), dimension(:,:), allocatable :: buffer
      logical                              :: lincore

      write(ilog,'(/,70a)') ('*',i=1,70)
      write(ilog,'(12x,a)') &
           'Block-Lanczos Diagonalisation in the Final Space'
      write(ilog,'(70a)') ('*',i=1,70)

      call cpu_time(t1)

!-----------------------------------------------------------------------
! Determine the buffer size (in terms of the no. Lanczos vectors that
! we can hold in memory)
!-----------------------------------------------------------------------
      buffsize=int(floor(((lancmem/3.0d0)*1024.0d0**2)/(8.0d0*matdim)))

      if (buffsize.gt.nvec) buffsize=nvec
      if (buffsize.lt.lmain) buffsize=lmain

      allocate(buffer(matdim,buffsize))

      nv=0
      nwr=0

      call freeunit(lanunit)
      
      reclength=8*matdim*buffsize

      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',&
           status='unknown',access='direct',recl=reclength)

!-----------------------------------------------------------------------
! Determine whether we can run the Lanczos vector generation in-core
!-----------------------------------------------------------------------
      ! On-diagonal elements
      mem=8.0d0*matdim/1024.0d0**2

      ! Off-diagonal elements
      mem=mem+8.0d0*noff/1024.0d0**2

      ! Off-diagonal indices (2 times integer*4 for each element)
      mem=mem+8.0d0*noff/1024.0d0**2

      ! Account for the size of the buffer that we are using to
      ! hold the Lanczos vectors...
      mem=mem+8.0d0*buffsize*matdim/1024.0d0**2

      ! Set the logical flag lincore to true if we can fit the
      ! Hamiltonian matrix in memory
      if (mem.le.lancmem) then
         lincore=.true.
         write(ilog,'(/,2x,a)') 'Matrix-vector multiplication &
              will proceed in-core'
      else
         lincore=.false.
         write(ilog,'(/,2x,a)') 'Matrix-vector multiplication &
              will proceed out-of-core'
      endif

!-----------------------------------------------------------------------
! If we are to calculate the Lanczos vectors in-core, then read the
! Hamiltonian matrix from file
!-----------------------------------------------------------------------
      if (lincore) call rdham(matdim,noff,hii,hij,indxi,indxj)

!-----------------------------------------------------------------------
! Start the block Lanczos iterations
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,1x,i4)') 'Block size:',lmain
      write(ilog,'(/,2x,i4,1x,a,/)') ncycles*lmain,&
           'Lanczos vectors will be generated'

      do j=1,ncycles

!-----------------------------------------------------------------------
! Output progress
!-----------------------------------------------------------------------
         write(ilog,'(70a)') ('*',k=1,70)
         write(ilog,'(2x,a,1x,i4)') 'Iteration number',j

!-----------------------------------------------------------------------
! Writing of the Lanczos vectors to disk
!-----------------------------------------------------------------------

         nprev=nv
         nv=nv+lmain

         if (nv.eq.buffsize) then
            ! Write all vectors to disk
            nwr=nwr+1
            buffer(:,nprev+1:buffsize)=qmat2(:,:)
            write(lanunit,rec=nwr) buffer
            buffer=0.0d0
            nv=0
         else if (nv.gt.buffsize) then
            ! Write some of the vectors to disk and
            ! save the rest to the buffer
            nwr=nwr+1
            nsurplus=nv-buffsize
            k=buffsize-nprev
            buffer(:,nprev+1:buffsize)=qmat2(:,1:k)
            write(lanunit,rec=nwr) buffer
            buffer=0.0d0
            buffer(:,1:nsurplus)=qmat2(:,k+1:lmain)
            nv=nsurplus
         else
            ! Save the vectors to the buffer
            buffer(:,nprev+1:nv)=qmat2(:,:)
         endif

! If we are on the last iteration, make sure that the buffer has been
! written to disk
         if (j.eq.ncycles.and.nv.lt.buffsize.and.nv.gt.0) then
            nwr=nwr+1
            write(lanunit,rec=nwr) buffer
         endif

!-----------------------------------------------------------------------
! Calculate the current block of on-diagonal elements of the T-matrix
!-----------------------------------------------------------------------
         if (lincore) then
            call hxq_incore(matdim,noff,hii,hij,indxi,indxj)
         else
            call hxq_ext(matdim)
         endif

         call dgemm('N','T',matdim,lmain,lmain,1.0d0,qmat1,matdim,bmat,&
              lmain,0.0d0,tmpmat,matdim)

         umat=umat-tmpmat
         
         call dgemm('T','N',lmain,lmain,matdim,1.0d0,qmat2,matdim,umat,&
              matdim,0.0d0,amat,lmain)
         
!-----------------------------------------------------------------------
! Calculate the next block of Krylov vectors
!-----------------------------------------------------------------------
         call dgemm('N','N',matdim,lmain,lmain,1.0d0,qmat2,matdim,amat,&
              lmain,0.0d0,tmpmat,matdim)

         rmat=umat-tmpmat
         
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
         i1=0                        ! Initialise the column counter
         do m=(j-1)*lmain+1,j*lmain  ! Loop over columns of T_j
            i1=i1+1                  ! Increment the column counter
            if (j.lt.ncycles) then
               upper=(j+1)*lmain
            else
               upper=j*lmain
            endif
            j1=0                     ! Initialise the row counters
            k1=0
            k2=0
            do n=(j-1)*lmain+1,upper ! Loop over rows of T_j
               j1=j1+1               ! Increment the main row counter
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

      close(lanunit)

      call cpu_time(t2)

      write(ilog,'(70a)') ('*',k=1,70)
      write(ilog,'(/,2x,a,1x,F8.2,1x,a1,/)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! Set the number of Lanczos states
!-----------------------------------------------------------------------
      lancstates=ncycles*lmain
      
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(buffer)
      deallocate(qmat1)
      deallocate(qmat2)
      deallocate(rmat)
      deallocate(amat)
      deallocate(bmat)
      deallocate(tmpmat)

!-----------------------------------------------------------------------
! Calculate the Lanczos pseudospectrum
!-----------------------------------------------------------------------
      call calc_pseudospec(lanunit,matdim)

      call cpu_time(t2)

      write(ilog,'(/,2x,a,/)') 'End of the block-Lanczos routine'
      write(ilog,'(2x,a,1x,F8.2,1x,a1,/)') 'Total time:',t2-t1,'s'
      write(ilog,'(70a,/)') ('*',i=1,70)

      return

    end subroutine run_block_lanczos

!#######################################################################

    subroutine rdham(ndim,noff,hii,hij,indxi,indxj)

      use iomod, only: freeunit

      implicit none

      integer, intent(in)                :: ndim
      integer*8, intent(in)              :: noff
      integer, dimension(:), allocatable :: indxi,indxj
      real(d), dimension(:), allocatable :: hii,hij

      integer                            :: unit,maxbl,nrec,count,k,&
                                            nlim
      integer, dimension(:), allocatable :: itmp1,itmp2
      real(d), dimension(:), allocatable :: ftmp

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hii(ndim))
      allocate(hij(noff))
      allocate(indxi(noff))
      allocate(indxj(noff))
      
!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
      call freeunit(unit)

      open(unit,file='SCRATCH/hmlt.diac',status='old',&
           access='sequential',form='unformatted')
      
      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
      allocate(ftmp(maxbl))
      allocate(itmp1(maxbl))
      allocate(itmp2(maxbl))

      open(unit,file='SCRATCH/hmlt.offc',status='old',&
           access='sequential',form='unformatted')

      count=0
      do k=1,nrec         
         read(unit) ftmp(:),itmp1(:),itmp2(:),nlim
         hij(count+1:count+nlim)=ftmp(1:nlim)
         indxi(count+1:count+nlim)=itmp1(1:nlim)
         indxj(count+1:count+nlim)=itmp2(1:nlim)
         count=count+nlim
      enddo

      close(unit)

      deallocate(ftmp,itmp1,itmp2)

      return

    end subroutine rdham

!#######################################################################

    subroutine hxq_incore(matdim,noff,hii,hij,indxi,indxj)

      implicit none
      
      integer, intent(in)        :: matdim
      integer*8, intent(in)      :: noff
      integer                    :: m,n,k
      integer, dimension(noff)   :: indxi,indxj
      real(d), dimension(matdim) :: hii
      real(d), dimension(noff)   :: hij

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
!      umat=0.0d0
!      do m=1,matdim
!         do n=1,lmain
!            umat(m,n)=hii(m)*qmat2(m,n)
!         enddo
!      enddo

      !$omp parallel do private(m,n) shared(umat,hii,qmat2)
      do n=1,lmain
         do m=1,matdim
            umat(m,n)=hii(m)*qmat2(m,n)
         enddo
      enddo
      !$omp end parallel do

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
!      do k=1,noff
!         umat(indxi(k),:)=umat(indxi(k),:)+hij(k)*qmat2(indxj(k),:)
!         umat(indxj(k),:)=umat(indxj(k),:)+hij(k)*qmat2(indxi(k),:)
!      enddo

      !$omp parallel do private(k,n) shared(umat,hij,qmat2,indxi,indxj)
      do n=1,lmain
         do k=1,noff
            umat(indxi(k),n)=umat(indxi(k),n)+hij(k)*qmat2(indxj(k),n)
            umat(indxj(k),n)=umat(indxj(k),n)+hij(k)*qmat2(indxi(k),n)
         enddo
      enddo
      !$omp end parallel do

      return

    end subroutine hxq_incore

!#######################################################################

    subroutine hxq_ext(matdim)

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: unit
      integer                            :: maxbl,nrec,nlim,i,j,k,l,m,n
      integer, dimension(:), allocatable :: indxi,indxj
      real(d), dimension(:), allocatable :: hii,hij

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hii(matdim))

      unit=77
      open(unit,file='SCRATCH/hmlt.diac',status='old',access='sequential',&
           form='unformatted')

      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

      umat=0.0d0
!      do n=1,matdim
!         do m=1,lmain
!            umat(m,n)=hii(m)*qmat2(m,n)
!         enddo
!      enddo

      !$omp parallel do private(m,n) shared(umat,hii,qmat2)
      do n=1,lmain
         do m=1,matdim
            umat(m,n)=hii(m)*qmat2(m,n)
         enddo
      enddo
      !$omp end parallel do

      deallocate(hii)

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
      
      open(unit,file='SCRATCH/hmlt.offc',status='old',access='sequential',&
           form='unformatted')

!      do k=1,nrec
!         read(unit) hij(:),indxi(:),indxj(:),nlim
!         do l=1,nlim
!            i=indxi(l)
!            j=indxj(l)
!            umat(i,n)=umat(i,n)+hij(l)*qmat2(j,n)
!            umat(j,n)=umat(j,n)+hij(l)*qmat2(i,n)
!         enddo
!      enddo

      do k=1,nrec
         read(unit) hij(:),indxi(:),indxj(:),nlim
         !$omp parallel do private(l,n) shared(umat,hij,qmat2,indxi,indxj)
         do n=1,lmain
            do l=1,nlim               
               umat(indxi(l),n)=umat(indxi(l),n)+hij(l)*qmat2(indxj(l),n)
               umat(indxj(l),n)=umat(indxj(l),n)+hij(l)*qmat2(indxi(l),n)
            enddo
         enddo
         !$omp end parallel do
      enddo

      close(unit)

      deallocate(hij,indxi,indxj)
      
      return
      
    end subroutine hxq_ext

!#######################################################################

    subroutine calc_pseudospec(lanunit,matdim)

      implicit none

      integer                       :: lanunit,matdim,i
      real(d), dimension(nvec,nvec) :: eigvec
      real(d), dimension(nvec)      :: eigval
      real(d)                       :: t1,t2,mem

!-----------------------------------------------------------------------
! (1) Lanczos state energies
!-----------------------------------------------------------------------
      write(ilog,'(70a)') ('-',i=1,70)
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state energies...'

      call cpu_time(t1)

      call diagmat_banded(lmain,tmat,nvec,eigvec,eigval)

      call cpu_time(t2)

      write(ilog,'(2x,a,1x,F8.2,1x,a1)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! (2) Lanczos state vectors
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state vectors...'

      call cpu_time(t1)

      mem=2*8.0d0*matdim*nvec/1024.0d0**2

      if (mem.le.lancmem) then
         write(ilog,'(2x,a,/)') 'Calculation will proceed in-core'
         call ritzvecs_incore(lanunit,eigvec,eigval,matdim)
      else
         write(ilog,'(2x,a,/)') 'Calculation will proceed &
              out-of-core'
         call ritzvecs_ext2(lanunit,eigvec,nvec,matdim,eigval)
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

    subroutine ritzvecs_ext2(lanunit,eigvec,nvec,matdim,eigval)
      
      use iomod, only: freeunit

      implicit none

      integer                              :: lanunit,nvec,matdim,&
                                              ritzunit,blocksize,i,j,&
                                              nblocks,reclength,k,&
                                              tmpunit,k1,k2,l1,l2,nk,nl
      real(d), dimension(nvec,nvec)        :: eigvec
      real(d), dimension(nvec)             :: eigval
      real(d), dimension(:,:), allocatable :: rvec,lvec,tmpmat

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
      reclength=8*matdim*buffsize

      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',&
              status='unknown',access='direct',recl=reclength)

      call freeunit(tmpunit)
      open(tmpunit,file='SCRATCH/tmpvecs',form='unformatted',&
           status='unknown',access='direct',recl=reclength)

      call freeunit(ritzunit)
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(lvec(matdim,buffsize))
      allocate(rvec(matdim,buffsize))
      allocate(tmpmat(matdim,buffsize))

      lvec=0.0d0
      rvec=0.0d0
      tmpmat=0.0d0

!-----------------------------------------------------------------------
! Complete blocks of Lanczos vectors
!-----------------------------------------------------------------------
      do i=1,nwr-1

         read(lanunit,rec=i) lvec

         k1=(i-1)*buffsize+1 ! First Lanczos vector in the current block         
         k2=i*buffsize ! Last Lanczos vector in the current block
         
         ! Complete blocks of Ritz vectors
         do j=1,nwr-1
            if (i.gt.1) read(tmpunit,rec=j) rvec
            l1=(j-1)*buffsize+1 ! First Ritz vector in the current block            
            l2=j*buffsize ! Last Ritz vector in the current block
            call dgemm('N','N',matdim,buffsize,buffsize,1.0d0,&
                 lvec,matdim,eigvec(k1:k2,l1:l2),buffsize,0.0d0,&
                 tmpmat,matdim)
            rvec=rvec+tmpmat
            write(tmpunit,rec=j) rvec
         enddo

         ! Potentially incomplete block of Ritz vectors
         if (i.gt.1) read(tmpunit,rec=nwr) rvec
         l1=(nwr-1)*buffsize+1
         l2=nvec
         nl=l2-l1+1
         call dgemm('N','N',matdim,nl,buffsize,1.0d0,&
                 lvec,matdim,eigvec(k1:k2,l1:l2),buffsize,0.0d0,&
                 tmpmat(:,1:nl),matdim)
         rvec(:,1:nl)=rvec(:,1:nl)+tmpmat(:,1:nl)
         write(tmpunit,rec=nwr) rvec

      enddo

!-----------------------------------------------------------------------
! Potentially incomplete block of Lanczos vectors
!-----------------------------------------------------------------------
      read(lanunit,rec=nwr) lvec
      k1=(nwr-1)*buffsize+1
      k2=nvec
      nk=k2-k1+1

      ! Complete blocks of Ritz vectors
      do j=1,nwr-1
         read(tmpunit,rec=j) rvec
         l1=(j-1)*buffsize+1 ! First Ritz vector in the current block            
         l2=j*buffsize ! Last Ritz vector in the current block         
         call dgemm('N','N',matdim,buffsize,nk,1.0d0,&
                 lvec(:,1:nk),matdim,eigvec(k1:k2,l1:l2),nk,0.0d0,&
                 tmpmat,matdim)
         rvec=rvec+tmpmat
         write(tmpunit,rec=j) rvec
      enddo
      
      ! Potentially incomplete block of Ritz vectors
      if (nwr.gt.1) read(tmpunit,rec=nwr) rvec
      l1=(nwr-1)*buffsize+1
      l2=nvec
      nl=l2-l1+1
      call dgemm('N','N',matdim,nl,nk,1.0d0,&
           lvec(:,1:nk),matdim,eigvec(k1:k2,l1:l2),nk,0.0d0,&
           tmpmat(:,1:nl),matdim)
      rvec(:,1:nl)=rvec(:,1:nl)+tmpmat(:,1:nl)
      write(tmpunit,rec=nwr) rvec

!-----------------------------------------------------------------------
! Write the Ritz vectors to file using the format required for future
! use
!-----------------------------------------------------------------------
      ! Complete blocks of Ritz vectors
      k=0
      do i=1,nwr-1
         read(tmpunit,rec=i) rvec
         do j=1,buffsize
            k=k+1
            write(ritzunit) k,eigval(k),rvec(:,j)
         enddo
      enddo

      ! Potentially incomplete block of Ritz vectors
      read(tmpunit,rec=nwr) rvec
      l1=(nwr-1)*buffsize+1
      l2=nvec
      nl=l2-l1+1
      do i=1,nl
         k=k+1
         write(ritzunit) k,eigval(k),rvec(:,i)
      enddo

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------      
      close(lanunit)
      close(tmpunit)
      close(ritzunit)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(lvec)
      deallocate(rvec)
      deallocate(tmpmat)

      return

    end subroutine ritzvecs_ext2

!#######################################################################

    subroutine ritzvecs_incore(lanunit,eigvec,eigval,matdim)

      implicit none

      integer                              :: lanunit,ritzunit,i,k1,&
                                              k2,nk,matdim,reclength
      real(d), dimension(nvec,nvec)        :: eigvec
      real(d), dimension(nvec)             :: eigval
      real(d), dimension(:,:), allocatable :: lvec,rvec,buffer

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(lvec(matdim,nvec))
      allocate(rvec(matdim,nvec))
      allocate(buffer(matdim,buffsize))

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
      reclength=8*matdim*buffsize

      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',&
              status='unknown',access='direct',recl=reclength)

      ritzunit=lanunit+1
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')    

!-----------------------------------------------------------------------
! Read the Lanczos vectors from file
!-----------------------------------------------------------------------
      ! Full buffers
      do i=1,nwr-1
         read(lanunit,rec=i) buffer
         k1=(i-1)*buffsize+1
         k2=i*buffsize
         lvec(:,k1:k2)=buffer
      enddo
      
      ! Potentially incomplete buffers
      read(lanunit,rec=nwr) buffer
      k1=(nwr-1)*buffsize+1
      k2=nvec
      nk=k2-k1+1
      lvec(:,k1:k2)=buffer(:,1:nk)

!-----------------------------------------------------------------------
! Calculate and output the Ritz vectors
!-----------------------------------------------------------------------
      call dgemm('N','N',matdim,nvec,nvec,1.0d0,lvec,matdim,eigvec,&
           nvec,0.0d0,rvec,matdim)

      do i=1,nvec
         write(ritzunit) i,eigval(i),rvec(:,i)
      enddo

!-----------------------------------------------------------------------
! Close the Lanczos and Ritz vector files
!-----------------------------------------------------------------------
      close(lanunit)
      close(ritzunit)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(rvec)
      deallocate(lvec)
      deallocate(buffer)

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
