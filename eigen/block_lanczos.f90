module block_lanczos

    use constants
    use parameters
    use channels
    
    implicit none

    save

    integer*8                              :: nvec,nwr,buffsize,&
                                              matdim,blocksize,&
                                              reclength
    integer                                :: nro
    integer, dimension(:), allocatable     :: lk,uk
    real(d), dimension(:,:), allocatable   :: buffer,qmat1,qmat2,umat,&
                                              rmat,amat,bmat,tmat,&
                                              tmpmat,omkj
    real(d), dimension(:,:,:), allocatable :: amat_all,bmat_all
    real(d), dimension(:), allocatable     :: anorm,bnorm
    real(d)                                :: epmach,eps,orthlim,eta
    character(len=1)                       :: hflag

  contains
    
!#######################################################################

    subroutine lancdiag_block(dimf,noff,flag)

      implicit none
      
      integer, intent(in)     :: dimf
      integer*8, intent(in)   :: noff
      character(1),intent(in) :: flag

!-----------------------------------------------------------------------
! Set the Hamiltonian flag
!-----------------------------------------------------------------------
      hflag=flag

!-----------------------------------------------------------------------
! Set dimensions
!-----------------------------------------------------------------------
      matdim=dimf
      blocksize=lmain
      nvec=blocksize*ncycles

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
      call alloc_block_lanczos

!-----------------------------------------------------------------------
! Set the initial Lanczos vectors
!-----------------------------------------------------------------------
      call init_vec

!-----------------------------------------------------------------------
! If a partial reorthogonalisation is to be performed, then estimate
! the machine epsilon, set the limit for reorthogonalisation and set
! the limit eta that is used in the MPRO algorithm
!-----------------------------------------------------------------------
      if (orthotype.gt.0) then
         epmach=machine_precision()
         orthlim=dsqrt(epmach)
         eta=epmach**0.75d0
      endif

!-----------------------------------------------------------------------
! Perform the block Lanczos calculation
!-----------------------------------------------------------------------
      call run_block_lanczos(noff)

!-----------------------------------------------------------------------
! Finalise/deallocate arrays
!-----------------------------------------------------------------------
      call finalise_block_lanczos

      return
      
    end subroutine lancdiag_block

!#######################################################################

    subroutine alloc_block_lanczos
    
      implicit none

      ! Main block Lanczos arays
      allocate(qmat1(matdim,blocksize))
      allocate(qmat2(matdim,blocksize))
      allocate(umat(matdim,blocksize))
      allocate(rmat(matdim,blocksize))
      allocate(tmpmat(matdim,blocksize))
      allocate(amat(blocksize,blocksize))
      allocate(bmat(blocksize,blocksize))
      allocate(tmat(nvec,nvec))
      qmat1=0.0d0
      qmat2=0.0d0
      umat=0.0d0
      rmat=0.0d0
      amat=0.0d0
      bmat=0.0d0
      tmat=0.0d0

      ! Reorthogonalisation arrays
      allocate(amat_all(ncycles,blocksize,blocksize))
      allocate(bmat_all(ncycles,blocksize,blocksize))
      allocate(omkj(ncycles+1,ncycles+1))
      allocate(anorm(ncycles))
      allocate(bnorm(ncycles))
      allocate(lk(ncycles))
      allocate(uk(ncycles))
      amat_all=0.0d0
      bmat_all=0.0d0
      omkj=0.0d0
      anorm=0.0d0
      bnorm=0.0d0
      lk=0
      uk=0

      return

    end subroutine alloc_block_lanczos

!#######################################################################
    
    subroutine init_vec
      
      use iomod, only: freeunit
      
      implicit none
      
      integer                              :: i,j,k
      integer                              :: iadc1,idim,ivecs
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
!
!             5 <-> RIXS calculation - read the initial vectors from
!                   file
!
!             6 <-> TPXAS calculation - read the initial vectors from
!                   file
!-----------------------------------------------------------------------

      if (lancguess.eq.1) then

         ! Copy the 1h1p of interest into the qmat2 array
         do i=1,blocksize
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
         do i=1,blocksize
            k=stvc_lbl(i)
            do j=1,idim
               qmat2(j,i)=adc1vec(j,k)
            enddo
         enddo
            
      else if (lancguess.eq.3) then

         ! Copy the linear combinations of the 1h1p and 2h2p ISs into
         ! the qmat2 array
         fac=1.0d0/sqrt(2.0d0)
         do i=1,blocksize
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
         do i=1,blocksize

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

      else if (lancguess.eq.5) then
         ! RIXS calculation: read the initial vectors from file
         call freeunit(ivecs)
         open(ivecs,file='SCRATCH/rixs_ivecs',form='unformatted',&
              status='old')
         read(ivecs) qmat2
         close(ivecs)
         
      else if (lancguess.eq.6) then
         ! TPXAS calculation
         if (hflag.eq.'i') then
            ! Valence-excited space calculation
            call freeunit(ivecs)
            open(ivecs,file='SCRATCH/tpxas_initi',form='unformatted',&
                 status='old')
            read(ivecs) qmat2
            close(ivecs)
         else if (hflag.eq.'c') then
            ! Core-excited space calculation
            call freeunit(ivecs)
            open(ivecs,file='SCRATCH/tpxas_initc',form='unformatted',&
                 status='old')
            read(ivecs) qmat2
            close(ivecs)
         endif
         
      endif
      
      return
      
    end subroutine init_vec

!#######################################################################

    function machine_precision() result(epsilon)
      
      implicit none

      real(d) :: epsilon,ftmp

      ! Estimation of the machine epsilon
      epsilon=1.0d0
5     ftmp=(1.0d0+0.5d0*epsilon)
      if (ftmp.ne.1.0d0) then
         epsilon=0.5d0*epsilon
         goto 5
      endif

      return

    end function machine_precision

!#######################################################################

    subroutine run_block_lanczos(noff)

      use iomod, only: freeunit

      implicit none

      integer*8, intent(in)                :: noff
      integer*8                            :: maxrecl
      integer                              :: lanunit,j,i,k,i1,j1,k1,k2,&
                                              m,n,upper,nv,nsurplus,&
                                              nprev,lanunit2,reclength2
      integer                              :: info
      integer, dimension(:), allocatable   :: indxi,indxj
      real(d), dimension(:), allocatable   :: hii,hij
      real(d), dimension(blocksize)        :: tau
      real(d), dimension(blocksize)        :: work
      real(d)                              :: t1,t2,mem
      logical                              :: lincore,lro,l2nd

      write(ilog,'(/,70a)') ('*',i=1,70)
      write(ilog,'(12x,a)') &
           'Generation of the block Lanczos pseudospectrum'
      write(ilog,'(70a)') ('*',i=1,70)

      call cpu_time(t1)

!-----------------------------------------------------------------------
! Initialise MPRO variables
!-----------------------------------------------------------------------
      l2nd=.false.

!-----------------------------------------------------------------------
! Determine the buffer size (in terms of the no. Lanczos vectors that
! we can hold in memory)
!
! Note that we only use maxmem/3 here and not maxmem because in the
! calculation of the Ritz vectors later on we are required to hold
! another two sets of vectors in core
!-----------------------------------------------------------------------
      ! Calculate the buffer size based on the amount of memory
      ! available
      buffsize=int(floor(((maxmem/3.0d0)*1024.0d0**2)/(8.0d0*matdim)))
      reclength=8*matdim*buffsize
      
      ! Make sure that the record length corresponding to the
      ! buffer size does not exceed the maximum value (2**31-1)
      maxrecl=2147483647
      if (reclength.gt.maxrecl) then
         buffsize=maxrecl/(8*matdim)
         reclength=8*matdim*buffsize
      endif

      ! Make sure that the buffer size is not greater than the
      ! total no. Lanczos vectors or less than the block size
      if (buffsize.gt.nvec) buffsize=nvec
      if (buffsize.lt.blocksize) buffsize=blocksize
      
      ! Allocate the buffer, initialise counters and open the
      ! Lanczos vector file
      allocate(buffer(matdim,buffsize))

      nv=0
      nwr=0

      call freeunit(lanunit)

      open(lanunit,file='SCRATCH/lanvecs',form='unformatted',&
           status='unknown',access='direct',recl=reclength)

      call freeunit(lanunit2)
      reclength2=8*matdim*blocksize
      open(lanunit2,file='SCRATCH/lanvecs2',form='unformatted',&
           status='unknown',access='direct',recl=reclength2)

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
      if (mem.le.maxmem) then
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
      if (lincore) call rdham(noff,hii,hij,indxi,indxj)

!-----------------------------------------------------------------------
! Start the block Lanczos iterations
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,1x,i4)') 'Block size:',blocksize
      write(ilog,'(/,2x,i5,1x,a,/)') ncycles*blocksize,&
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
         ! (1) File to be used in the PRO routines
         write(lanunit2,rec=j) qmat2

         ! (2) File to be used in the calculation of the Ritz vectors
         nprev=nv
         nv=nv+blocksize

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
            buffer(:,1:nsurplus)=qmat2(:,k+1:blocksize)
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
            call hxq_incore(noff,hii,hij,indxi,indxj)
         else
            call hxq_ext
         endif

         call dgemm('N','T',matdim,blocksize,blocksize,1.0d0,qmat1,matdim,bmat,&
              blocksize,0.0d0,tmpmat,matdim)

         umat=umat-tmpmat
         
         call dgemm('T','N',blocksize,blocksize,matdim,1.0d0,qmat2,matdim,umat,&
              matdim,0.0d0,amat,blocksize)
         
         amat_all(j,:,:)=amat

!-----------------------------------------------------------------------
! Calculate the next block of Krylov vectors
!-----------------------------------------------------------------------
         call dgemm('N','N',matdim,blocksize,blocksize,1.0d0,qmat2,matdim,amat,&
              blocksize,0.0d0,tmpmat,matdim)

         rmat=umat-tmpmat
         
!-----------------------------------------------------------------------
! Compute the QR factorization of the matrix of Krylov vectors
!-----------------------------------------------------------------------
         ! dgeqrf will overwrite rmat, so instead use a copy of this
         ! matrix
         tmpmat=rmat

         ! Compute the current block of off-diagonal elements of
         ! the T-matrix
         call dgeqrf(matdim,blocksize,tmpmat,matdim,tau,work,blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'dqerf failed in run_block_lanczos'
            STOP
         endif

         ! Note that the B-matrix is upper-triangular
         bmat=0.0d0
         do k=1,blocksize
            do i=1,k
               bmat(i,k)=tmpmat(i,k)
            enddo
         enddo

         bmat_all(j,:,:)=bmat

         ! Extract the next block of Lanczos vectors
         call dorgqr(matdim,blocksize,blocksize,tmpmat,matdim,tau,work,&
              blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'dorgqr failed in run_block_lanczos'
            STOP
         endif

         ! Update the matrices of Lanczos vectors
         qmat1=qmat2
         qmat2=tmpmat

!-----------------------------------------------------------------------
! Partial orthogonalisation
!-----------------------------------------------------------------------
         lro=.false.
         if (j.gt.1.and.j.lt.ncycles) then
            if (orthotype.eq.1) then
               call pro(j,lanunit2,lro)
            else if (orthotype.eq.2) then
               call mpro(j,lanunit2,lro,l2nd)
            endif
         endif

!-----------------------------------------------------------------------
! Local reorthogonalisation of qmat2 against qmat1
!
! Note that this is only necessary a reorthogonalisation has not taken 
! place in this iteration
!-----------------------------------------------------------------------
         if (.not.lro.and.orthotype.eq.1) call localro

!-----------------------------------------------------------------------
! Fill in the next block of the T-matrix array
!-----------------------------------------------------------------------
         i1=0                               ! Initialise the column counter
         do m=(j-1)*blocksize+1,j*blocksize ! Loop over columns of T_j
            i1=i1+1                         ! Increment the column counter
            if (j.lt.ncycles) then
               upper=(j+1)*blocksize
            else
               upper=j*blocksize
            endif
            j1=0                            ! Initialise the row counters
            k1=0
            k2=0
            do n=(j-1)*blocksize+1,upper    ! Loop over rows of T_j
               j1=j1+1                      ! Increment the main row counter
               if (j1.le.blocksize) then    ! Contribution from A_j
                  k1=k1+1
                  tmat(n,m)=amat(k1,i1)
               else                         ! Contribution from B_j
                  k2=k2+1
                  tmat(n,m)=bmat(k2,i1)
               endif
               tmat(m,n)=tmat(n,m)
            enddo
         enddo

      enddo

      close(lanunit)
      close(lanunit2)

      call cpu_time(t2)

      write(ilog,'(70a)') ('*',k=1,70)
      write(ilog,'(/,2x,a,1x,F8.2,1x,a1,/)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! Set the number of Lanczos states
!-----------------------------------------------------------------------
      lancstates=ncycles*blocksize

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(buffer)
      if (allocated(hii)) deallocate(hii)
      if (allocated(hij)) deallocate(hij)

!-----------------------------------------------------------------------
! Calculate the Lanczos pseudospectrum
!-----------------------------------------------------------------------
      call calc_pseudospec(lanunit)

      call cpu_time(t2)

      write(ilog,'(/,2x,a,/)') 'End of the block-Lanczos routine'
      write(ilog,'(2x,a,1x,F8.2,1x,a1,/)') 'Total time:',t2-t1,'s'
      write(ilog,'(70a,/)') ('*',i=1,70)

      return

    end subroutine run_block_lanczos

!#######################################################################

    subroutine rdham(noff,hii,hij,indxi,indxj)

      use iomod, only: freeunit

      implicit none

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
      allocate(hii(matdim))
      allocate(hij(noff))
      allocate(indxi(noff))
      allocate(indxj(noff))
      
!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
      call freeunit(unit)

      open(unit,file='SCRATCH/hmlt.dia'//hflag,status='old',&
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

      open(unit,file='SCRATCH/hmlt.off'//hflag,status='old',&
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

    subroutine hxq_incore(noff,hii,hij,indxi,indxj)

      implicit none
      
      integer*8, intent(in)      :: noff
      integer                    :: m,n,k
      integer, dimension(noff)   :: indxi,indxj
      real(d), dimension(matdim) :: hii
      real(d), dimension(noff)   :: hij

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      umat=0.0d0
      !$omp parallel do private(m,n) shared(umat,hii,qmat2)
      do n=1,blocksize
         do m=1,matdim
            umat(m,n)=hii(m)*qmat2(m,n)
         enddo
      enddo
      !$omp end parallel do

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      !$omp parallel do private(k,n) shared(umat,hij,qmat2,indxi,indxj)
      do n=1,blocksize
         do k=1,noff
            umat(indxi(k),n)=umat(indxi(k),n)+hij(k)*qmat2(indxj(k),n)
            umat(indxj(k),n)=umat(indxj(k),n)+hij(k)*qmat2(indxi(k),n)
         enddo
      enddo
      !$omp end parallel do

      return

    end subroutine hxq_incore

!#######################################################################

    subroutine hxq_ext

      implicit none

      integer                            :: unit
      integer                            :: maxbl,nrec,nlim,i,j,k,l,m,n
      integer, dimension(:), allocatable :: indxi,indxj
      real(d), dimension(:), allocatable :: hii,hij

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hii(matdim))

      unit=77
      open(unit,file='SCRATCH/hmlt.dia'//hflag,status='old',&
           access='sequential',form='unformatted')

      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

      umat=0.0d0
      !$omp parallel do private(m,n) shared(umat,hii,qmat2)
      do n=1,blocksize
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
      
      open(unit,file='SCRATCH/hmlt.off'//hflag,status='old',&
           access='sequential',form='unformatted')

      do k=1,nrec
         read(unit) hij(:),indxi(:),indxj(:),nlim
         !$omp parallel do private(l,n) shared(umat,hij,qmat2,indxi,indxj)
         do n=1,blocksize
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

    subroutine localro

      implicit none

      integer                       :: k,n
      real(d), dimension(blocksize) :: dpvec
      real(d)                       :: dp,ddot
      
      external ddot

!-----------------------------------------------------------------------
! Orthogonalise Q_j+1 against Q_j
!-----------------------------------------------------------------------
      ! Loop over vectors in Q_j
      do n=1,blocksize

         ! Calculate the dot products between all
         ! the vectors in Q_j+1 and the nth vector
         ! in Q_j
         call dgemv('T',matdim,blocksize,1.0d0,qmat2,matdim,&
              qmat1(:,n),1,0.0d0,dpvec,1)

         ! Orthogonalise all the vectors in Q_j+1
         ! against the nth vector in Q_j
         do k=1,blocksize
            qmat2(:,k)=qmat2(:,k)-dpvec(k)*qmat1(:,n)
         enddo

      enddo

!-----------------------------------------------------------------------
! Normalise Q_j+1
!-----------------------------------------------------------------------
      do k=1,blocksize
         dp=ddot(matdim,qmat2(:,k),1,qmat2(:,k),1)
         qmat2(:,k)=qmat2(:,k)/dsqrt(dp)
      enddo

      return
      
    end subroutine localro

!#######################################################################

    subroutine pro(j,lanunit2,lro)

      implicit none

      integer :: j,i,k,lanunit2
      real(d) :: maxom
      logical :: lro

!-----------------------------------------------------------------------
! Update the omega recurrence
!-----------------------------------------------------------------------
      call update_omega(j)

!-----------------------------------------------------------------------
! TEST: calculate the actual orthogonality components
!-----------------------------------------------------------------------
!      call true_orthog(j)
!
!      maxom=0.0d0
!      do k=1,j-1
!         if (omkj(k,j+1).gt.maxom) then
!            maxom=omkj(k,j+1)
!         endif
!      enddo
!      print*,j,maxom

!-----------------------------------------------------------------------
! If the estimated orthogonality components are above the threshold,
! then perform a partial reorthogonalisation
!-----------------------------------------------------------------------
      do k=1,j-1
         if (omkj(j+1,k).ge.orthlim) lro=.true.
      enddo

      if (lro) then

         write(ilog,'(/,2x,a,/)') 'Performing partial &
              reorthogonalisation...'

         ! MGS orthogonalisation
         call pro_mgs(lanunit2,j)

         ! Update the omega-recurrence
         do k=1,j
            omkj(j+1,k)=eps
            omkj(j,k)=eps
         enddo

      endif

      return

    end subroutine pro

!#######################################################################

    subroutine update_omega(j)

      implicit none

      integer                                 :: j,k
      real(d)                                 :: invssvb

      integer                                 :: info
      real(d), dimension(blocksize,blocksize) :: u,vt
      real(d), dimension(blocksize)           :: sigma
      real(d), dimension(5*blocksize)         :: work

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
      if (j.eq.2) then
         
         ! Spectral norm of A_1
         call dgesvd('N','N',blocksize,blocksize,amat_all(1,:,:),&
              blocksize,sigma,u,blocksize,vt,blocksize,work,&
              5*blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'Error in the SVD of A_1'
            STOP
         endif
         anorm(1)=sigma(1)
         
         ! Spectral norm of B_1
         call dgesvd('N','N',blocksize,blocksize,bmat_all(1,:,:),&
              blocksize,sigma,u,blocksize,vt,blocksize,work,&
              5*blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'Error in the SVD of B_1'
            STOP
         endif

         bnorm(1)=sigma(1)
         
         ! Spectral norm of B_2
         call dgesvd('N','N',blocksize,blocksize,bmat_all(2,:,:),&
              blocksize,sigma,u,blocksize,vt,blocksize,work,&
              5*blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'Error in the SVD of B_2'
            STOP
         endif

         bnorm(2)=sigma(1)
         
         ! Set the various values that depend on the machine epsilon
         eps=epmach!*blocksize*dsqrt(dble(matdim))
         
         ! Initialisation of the omkj array
         omkj(1,1)=eps
         omkj(2,2)=eps
         omkj(1,2)=eps/sigma(blocksize)
         omkj(2,1)=omkj(1,2)

      endif

!-----------------------------------------------------------------------
! Calculate the next set of estimated orthogonality components
!-----------------------------------------------------------------------
      ! Spectral norm of A_j
      call dgesvd('N','N',blocksize,blocksize,amat_all(j,:,:),&
           blocksize,sigma,u,blocksize,vt,blocksize,work,&
           5*blocksize,info)
      if (info.ne.0) then
         write(ilog,'(/,2x,a,/)') 'Error in the SVD of A_j'
         STOP
      endif

      anorm(j)=sigma(1)

      ! Spectral norm of B_j
      call dgesvd('N','N',blocksize,blocksize,bmat_all(j,:,:),&
           blocksize,sigma,u,blocksize,vt,blocksize,work,&
           5*blocksize,info)
      if (info.ne.0) then
         write(ilog,'(/,2x,a,/)') 'Error in the SVD of B_j'
         STOP
      endif

      bnorm(j)=sigma(1)     

      invssvb=1.0d0/sigma(blocksize)
      
      ! omega_j,j+1 and omega_j+1,j+1
      omkj(j,j+1)=eps
      omkj(j+1,j)=eps
      omkj(j+1,j+1)=eps

      ! omega_k,j+1, k=1,...,j-1
      do k=1,j-1
         if (k.eq.1) then
            omkj(j+1,k)=bnorm(k+1)*omkj(j,k+1) + bnorm(j)*omkj(j-1,k) &
                        + (anorm(j)+anorm(k))*omkj(j,k)
         else
            omkj(j+1,k)=bnorm(k+1)*omkj(j,k+1) + bnorm(k)*omkj(j,k-1) &
                        + bnorm(j)*omkj(j-1,k) &
                        + (anorm(j)+anorm(k))*omkj(j,k)
         endif
         omkj(j+1,k)=invssvb*omkj(j+1,k)
         omkj(k,j+1)=omkj(j+1,k)
      enddo

      return

    end subroutine update_omega

!#######################################################################

    subroutine pro_mgs(lanunit2,j)

      implicit none

      integer                              :: lanunit2,j,k,m,n
      real(d), dimension(matdim,blocksize) :: lmat
      real(d), dimension(blocksize)        :: dpvec1,dpvec2
      real(d)                              :: dp,ddot

      external ddot

!-----------------------------------------------------------------------
! Orthogonalise Q_j+1 and Q_j against Q_k, k=1,j-1
!-----------------------------------------------------------------------
      ! Loop over blocks of vectors Q_k
      do k=1,j-1

         ! Read the current block of vectors
         read(lanunit2,rec=k) lmat

         ! Loop over the vectors in the current block
         do m=1,blocksize

            ! Calculate the dot products between all
            ! the vectors in Q_j and the mth vector
            ! in Q_k
            call dgemv('T',matdim,blocksize,1.0d0,qmat1,matdim,&
                 lmat(:,m),1,0.0d0,dpvec1,1)

            ! Calculate the dot products between all
            ! the vectors in Q_j+1 and the mth vector
            ! in Q_k
            call dgemv('T',matdim,blocksize,1.0d0,qmat2,matdim,&
                 lmat(:,m),1,0.0d0,dpvec2,1)

            ! Orthogonalise all the vectors in Q_j Q_j+1
            ! against the mth vector in Q_k
            do n=1,blocksize
               qmat1(:,n)=qmat1(:,n)-dpvec1(n)*lmat(:,m)
               qmat2(:,n)=qmat2(:,n)-dpvec2(n)*lmat(:,m)
            enddo

         enddo

      enddo

!-----------------------------------------------------------------------
! Normalise Q_j
!-----------------------------------------------------------------------      
      do k=1,blocksize
         dp=ddot(matdim,qmat1(:,k),1,qmat1(:,k),1) 
         qmat1(:,k)=qmat1(:,k)/dsqrt(dp)
      enddo

!-----------------------------------------------------------------------
! Orthogonalise Q_j+1 against Q_j
!
! Note that Q_j+1 is also normalised in localro
!-----------------------------------------------------------------------
      call localro

      return

    end subroutine pro_mgs

!#######################################################################

    subroutine mpro(j,lanunit2,lro,l2nd)
      
      implicit none

      integer :: j,lanunit2,k,l
      logical :: lro,l2nd

!-----------------------------------------------------------------------
! Update the omega recurrence
!-----------------------------------------------------------------------
      call update_omega(j)

      if (l2nd) then
!-----------------------------------------------------------------------
! Reorthogonalisation was performed during the previous iteration, and 
! subsequently must also be performed during the current iteration.
!-----------------------------------------------------------------------
         lro=.true.

         call mpro_mgs(j,lanunit2)

         ! Update the omega values
         do k=1,nro
            do l=lk(k),uk(k)
               omkj(j+1,l)=eps
               omkj(l,j+1)=omkj(j+1,l)
            enddo
         enddo

         lk=0
         uk=0

         l2nd=.false.

      else
!-----------------------------------------------------------------------
! Reorthogonalisation wasn't performed during the previous iteration. 
! Check whether we need to perform reorthogonalisation during the 
! current iteration.
!-----------------------------------------------------------------------
         nro=0
         k=1
10       continue
         if (omkj(j+1,k).ge.orthlim) then
            lro=.true.
            nro=nro+1
            lk(nro)=lowerbound(j,k)
            uk(nro)=upperbound(j,k)
            k=k+uk(nro)
         else
            k=k+1
         endif
         if (k.lt.j-1) goto 10

         if (lro) then
            call mpro_mgs(j,lanunit2)
            do k=1,nro
               do l=lk(k),uk(k)
                  omkj(j+1,l)=eps
                  omkj(l,j+1)=omkj(j+1,l)
               enddo
              if (lk(k).gt.1) lk(k)=lk(k)-1
              uk(k)=uk(k)+1
           enddo
            l2nd=.true.
         endif

      endif

      return

    end subroutine mpro

!#######################################################################

    subroutine mpro_mgs(j,lanunit2)
      
      implicit none

      integer                              :: j,lanunit2,i,k,m,n,nblock
      real(d), dimension(matdim,blocksize) :: lmat
      real(d), dimension(blocksize)        :: dpvec
      real(d)                              :: dp,ddot

      external ddot

      nblock=0
      do i=1,nro
         nblock=nblock+(uk(i)-lk(i)+1)
      enddo

      write(ilog,'(/,2x,a,x,i3,x,a,/)') 'Performing modified partial &
           reorthogonalisation against',nblock,'blocks'

      ! Loop over intervals [l,u] of blocks Q_k that we will
      ! be orthogonalising Q_j+1 against
      do i=1,nro

         ! For each interval, loop over the corresponding blocks Q_k
         do k=lk(i),uk(i)

            ! Read the current block of vectors
            read(lanunit2,rec=k) lmat

            ! Loop over the vectors in the current block
            do m=1,blocksize

               ! Calculate the dot products between all
               ! the vectors in Q_j+1 and the mth vector
               ! in Q_k
               call dgemv('T',matdim,blocksize,1.0d0,qmat2,matdim,&
                    lmat(:,m),1,0.0d0,dpvec,1)

               ! Orthogonalise all the vectors in and Q_j+1
               ! against the mth vector in Q_k
               do n=1,blocksize
                  qmat2(:,n)=qmat2(:,n)-dpvec(n)*lmat(:,m)
               enddo

            enddo

         enddo

      enddo

!-----------------------------------------------------------------------
! Normalise Q_j+1
!-----------------------------------------------------------------------      
      do k=1,blocksize
         dp=ddot(matdim,qmat2(:,k),1,qmat2(:,k),1) 
         qmat2(:,k)=qmat2(:,k)/dsqrt(dp)
      enddo

      return

    end subroutine mpro_mgs

!#######################################################################

    function lowerbound(j,k) result(lb)

      implicit none

      integer :: lb,j,k,i

      if (k.eq.1.or.k.eq.2) then
         lb=1
      else
         lb=k-1
         do i=k-2,1
            if (abs(omkj(j+1,i)).gt.eta) then
               lb=i
            else
               exit
            endif
         enddo
      endif
      
      return

    end function lowerbound

!#######################################################################

    function upperbound(j,k) result(ub)

      implicit none
      
      integer :: ub,j,k,i

      ub=k+1

      do i=k+1,j-1
         if (abs(omkj(j+1,i)).ge.eta) then
            ub=i
         else
            exit
         endif
      enddo

      return

    end function upperbound

!#######################################################################

    subroutine true_orthog(j)

      implicit none

      integer                                 :: j,m,n,i1,i2
      real(d)                                 :: maxnorm
      real(d), dimension(blocksize,blocksize) :: overlap
      
      real(d), dimension(blocksize,blocksize) :: u,vt
      real(d), dimension(blocksize)           :: sigma
      real(d), dimension(5*blocksize)         :: work

      ! For the CS2 STO-3G test case, we never fill the buffer, so
      ! all Lanczos vectors from blocks 1,...,j will be in the buffer 
      ! array, whilst the vectors in block j+1 are currently in the
      ! qmat2 array

!      print*,
!      print*,"_________________________"
!      print*,"Iteration:",j
!      print*,"_________________________"
      maxnorm=0.0d0
      do m=1,j-1
         i1=(m-1)*blocksize+1
         i2=m*blocksize
         overlap=matmul(transpose(buffer(:,i1:i2)),qmat2(:,:))
         call dgesvd('N','N',blocksize,blocksize,overlap,&
              blocksize,sigma,u,blocksize,vt,blocksize,work,&
              5*blocksize,info)
         if (info.ne.0) then
            write(ilog,'(/,2x,a,/)') 'Error in the SVD of Qk^TQ_j+1'
            STOP
         endif
         if (sigma(1).gt.maxnorm) maxnorm=sigma(1)
      enddo

!      print*,"True:     ",maxnorm
      print*,j,maxnorm

      return

    end subroutine true_orthog

!#######################################################################
    
    subroutine calc_pseudospec(lanunit)

      implicit none

      integer                              :: lanunit,i
      real(d), dimension(:,:), allocatable :: eigvec
      real(d), dimension(:), allocatable   :: eigval
      real(d)                              :: t1,t2,mem

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(eigvec(nvec,nvec))
      allocate(eigval(nvec))

!-----------------------------------------------------------------------
! (1) Lanczos state energies
!-----------------------------------------------------------------------
      write(ilog,'(70a)') ('-',i=1,70)
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state energies...'

      call cpu_time(t1)

      call diagmat_banded(tmat,eigvec,eigval)

      call cpu_time(t2)

      write(ilog,'(2x,a,1x,F8.2,1x,a1)') 'Time taken:',t2-t1,'s'

      ! Dellocate tmat now that it is no longer needed
      deallocate(tmat)

!-----------------------------------------------------------------------
! (2) Lanczos state vectors
!-----------------------------------------------------------------------
      write(ilog,'(/,2x,a,/)') 'Calculating the Lanczos state vectors...'

      call cpu_time(t1)

      mem=2*8.0d0*matdim*nvec/1024.0d0**2

      if (mem.le.maxmem) then
         write(ilog,'(2x,a,/)') 'Calculation will proceed in-core'
         call ritzvecs_incore(lanunit,eigvec,eigval)
      else
         write(ilog,'(2x,a,/)') 'Calculation will proceed &
              out-of-core'
         call ritzvecs_ext2(lanunit,eigvec,eigval)
      endif

      call cpu_time(t2)
      write(ilog,'(2x,a,1x,F8.2,1x,a1,/)') 'Time taken:',t2-t1,'s'

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(eigvec)
      deallocate(eigval)

      return

    end subroutine calc_pseudospec

!#######################################################################

    subroutine diagmat_banded(matrix,eigvec,eigval)

      implicit none
      
      integer                              :: i,j,iupper
      real(d), dimension(nvec,nvec)        :: matrix

      integer                              :: kd,ldab,error
      real(d), dimension(blocksize+1,nvec) :: ab
      real(d), dimension(nvec)             :: eigval
      real(d), dimension(nvec,nvec)        :: eigvec
      real(d), dimension(3*nvec-2)         :: work

!-----------------------------------------------------------------------
! Set dimensions required to be passed to dsbev
!-----------------------------------------------------------------------
      kd=blocksize
      ldab=blocksize+1

!-----------------------------------------------------------------------
! Fill in the array ab holding the upper triangle of the projection of
! the Hamiltonian onto the space spanned by the Lanczos vectors
!-----------------------------------------------------------------------
      ab=0.0d0
      do j=1,nvec
         iupper=min(nvec,j+kd)
         do i=j,iupper
            ab(1+i-j,j)=matrix(i,j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Diagonalise the projection of the Hamiltonian onto the space spanned
! by the Lanczos vectors
!-----------------------------------------------------------------------
      call dsbev('V','L',nvec,kd,ab,ldab,eigval,eigvec,nvec,work,error)
        
      if (error.ne.0) then
         write(ilog,'(/,2x,3a,/)') 'Diagonalisation of the Lanczos ',&
              'representation of the Hamiltonian failed in ',&
              'subroutine diagmat_banded.'
         STOP
      endif

      return

    end subroutine diagmat_banded

!#######################################################################

    subroutine ritzvecs_ext2(lanunit,eigvec,eigval)
      
      use iomod, only: freeunit

      implicit none

      integer                              :: lanunit,ritzunit,&
                                              blocksize,i,j,nblocks,&
                                              k,tmpunit,k1,k2,l1,l2,&
                                              nk,nl
      real(d), dimension(nvec,nvec)        :: eigvec
      real(d), dimension(nvec)             :: eigval
      real(d), dimension(:,:), allocatable :: rvec,lvec,tmpmat
      real(d), dimension(:,:), allocatable :: dpmat

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
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
! Normalise the Ritz vectors (necessary when the Lanczos vectors are
! not orthogonal)
!-----------------------------------------------------------------------
      allocate(dpmat(buffsize,buffsize))

      ! Complete blocks of Ritz vectors
      do i=1,nwr-1
         read(tmpunit,rec=i) rvec
         call dgemm('T','N',buffsize,buffsize,matdim,1.0d0,rvec,&
              matdim,rvec,matdim,0.0d0,dpmat,buffsize)
         do j=1,buffsize
            rvec(:,j)=rvec(:,j)/dsqrt(dpmat(j,j))
         enddo
         write(tmpunit,rec=i) rvec
      enddo

      ! Potentially incomplete block of Ritz vectors
      read(tmpunit,rec=nwr) rvec
      l1=(nwr-1)*buffsize+1
      l2=nvec
      nl=l2-l1+1
      call dgemm('T','N',nl,nl,matdim,1.0d0,rvec(:,1:nl),&
              matdim,rvec(:,1:nl),matdim,0.0d0,dpmat(1:nl,1:nl),nl)
      do i=1,nl
         rvec(:,i)=rvec(:,i)/dsqrt(dpmat(i,i))
      enddo
      write(tmpunit,rec=nwr) rvec

      deallocate(dpmat)

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

    subroutine ritzvecs_incore(lanunit,eigvec,eigval)

      implicit none

      integer                              :: lanunit,ritzunit,i,k1,&
                                              k2,nk
      real(d), dimension(nvec,nvec)        :: eigvec,dpmat
      real(d), dimension(nvec)             :: eigval
      real(d), dimension(:,:), allocatable :: lvec,rvec

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(lvec(matdim,nvec))
      allocate(rvec(matdim,nvec))
      allocate(buffer(matdim,buffsize))

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
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
      ! Calculate the Ritz vectors
      call dgemm('N','N',matdim,nvec,nvec,1.0d0,lvec,matdim,eigvec,&
           nvec,0.0d0,rvec,matdim)

      ! Normalise the Ritz vectors
      call dgemm('T','N',nvec,nvec,matdim,1.0d0,rvec,matdim,rvec,&
           matdim,0.0d0,dpmat,nvec)
      do i=1,nvec
         rvec(:,i)=rvec(:,i)/dsqrt(dpmat(i,i))
      enddo

      ! Write the Ritz vectors to file
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

    subroutine chkortho(dim)

      implicit none

      integer*8                            :: unit,dim,i,j,k,count,&
                                              reclength2
      real(d), dimension(matdim,dim)       :: lvec
      real(d), dimension(matdim,blocksize) :: tmpvec
      real(d)                              :: dp,tmp
      real(d), parameter                   :: tol=1d-8

      lvec=0.0d0

      unit=28
      reclength2=8*matdim*blocksize
      open(unit,file='SCRATCH/lanvecs2',form='unformatted',&
           status='unknown',access='direct',recl=reclength2)

      count=0
      do i=1,ncycles
         read(unit,rec=i) tmpvec
         do j=1,blocksize
            count=count+1
            lvec(:,count)=tmpvec(:,j)
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

    subroutine finalise_block_lanczos

      implicit none

      deallocate(qmat1)
      deallocate(qmat2)
      deallocate(umat)
      deallocate(rmat)
      deallocate(tmpmat)
      deallocate(amat)
      deallocate(bmat)
      deallocate(amat_all)
      deallocate(bmat_all)
      deallocate(omkj)
      deallocate(anorm)
      deallocate(bnorm)
      deallocate(lk)
      deallocate(uk)

      return

    end subroutine finalise_block_lanczos

!#######################################################################

  end module block_lanczos
