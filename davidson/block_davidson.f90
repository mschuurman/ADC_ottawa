  module block_davidson

    use constants

    implicit none

    save

    integer                              :: blocksize,nstates,maxvec,&
                                            niter,currdim,ipre,nconv,&
                                            nrec,maxbl
    integer, dimension(:), allocatable   :: indxi,indxj
    real(d), dimension(:), allocatable   :: hii,hij
    real(d), dimension(:,:), allocatable :: vmat,wmat,rmat,ritzvec,&
                                            res,reigvec
    real(d), dimension(:), allocatable   :: reigval,norm
    real(d)                              :: tol
    character(len=36)                    :: vecfile
    logical                              :: lincore,lrdadc1

  contains

!#######################################################################
    
    subroutine davdiag_block(matdim,noffd)

      use channels
      use constants
      use parameters

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: k
      character(len=120)    :: atmp
      
!-----------------------------------------------------------------------
! Write to the log file
!-----------------------------------------------------------------------
      atmp='Block Davidson diagonalisation in the'
      if (hamflag.eq.'i') then
         atmp=trim(atmp)//' initial space'
      else if (hamflag.eq.'f') then
         atmp=trim(atmp)//' final space'
      endif
      write(ilog,'(/,70a)') ('-',k=1,70)
      write(ilog,'(2x,a)') trim(atmp)
      write(ilog,'(70a,/)') ('-',k=1,70)
      
!-----------------------------------------------------------------------
! Determine dimensions and allocate arrays
!-----------------------------------------------------------------------
      call davinitialise(matdim)

!-----------------------------------------------------------------------
! Set the initial vectors
!-----------------------------------------------------------------------
      ! Determine whether the initial vectors are to be constructed from
      ! the ADC(1) eigenvectors
      if (hamflag.eq.'i'.and.ladc1guess) then
         lrdadc1=.true.
      else if (hamflag.eq.'f'.and.ladc1guess_f) then
         lrdadc1=.true.
      else
         lrdadc1=.false.
      endif

      if (lrdadc1) then
         ! Construct the initial vectors from the ADC(1) eigenvectors
         call initvec_adc1
      else
         ! Use a single IS for each initial vector
         call initvec_ondiag(matdim)
      endif

!-----------------------------------------------------------------------
! Determine whether or not we can perform the matrix-vector
! multiplication in-core
!-----------------------------------------------------------------------
      call isincore(matdim,noffd)

!-----------------------------------------------------------------------
! Read the on-diagonal elements of the Hamiltonian matrix from disk
!-----------------------------------------------------------------------
      call rdham_on(matdim)

!-----------------------------------------------------------------------
! If the matrix-vector multiplication is to proceed in-core, then
! read the off-diagonal Hamiltonian matrix from disk
!-----------------------------------------------------------------------
      if (lincore) call rdham_off(noffd)

!-----------------------------------------------------------------------
! Perform the block Davidson iterations
!-----------------------------------------------------------------------
      call run_block_davidson(matdim,noffd)

!-----------------------------------------------------------------------
! Save the converged Ritz vectors and Ritz values to file
!-----------------------------------------------------------------------
      call wreigenpairs

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      call davfinalise
      
      return

    end subroutine davdiag_block

!#######################################################################

    subroutine davinitialise(matdim)

      use constants
      use parameters

      implicit none

      integer, intent(in) :: matdim
      integer             :: itmp
      real(d)             :: ftmp

!-----------------------------------------------------------------------
! Set the block size, no. eigenpairs, memory...
!-----------------------------------------------------------------------
      if (hamflag.eq.'i') then
         blocksize=dmain
         nstates=davstates
         niter=maxiter
         ipre=precon
         tol=davtol
         vecfile=davname
         maxvec=maxsubdim
      else if (hamflag.eq.'f') then
         blocksize=dmain_f
         nstates=davstates_f
         niter=maxiter_f
         ipre=precon_f
         tol=davtol_f
         vecfile=davname_f
         maxvec=maxsubdim_f
      endif

!-----------------------------------------------------------------------
! Temporary: set the maximum subspace dimension
!-----------------------------------------------------------------------
      ! If the maximum subspace dimension has not been set by the user
      ! then set it here
      if (maxvec.lt.0) then
         maxvec=min(nstates+blocksize+10,2*(nstates+blocksize))
      endif
         
      ! Make sure that maxvec is an integer multiple of the block size
      ftmp=real(maxvec)/real(blocksize)
      itmp=floor(ftmp)
      maxvec=itmp*blocksize

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Matrix of subspace vectors
      allocate(vmat(matdim,maxvec))
      vmat=0.0d0

      ! Matrix-vector product
      allocate(wmat(matdim,maxvec))
      wmat=0.0d0

      ! Rayleigh matrix
      allocate(rmat(maxvec,maxvec))
      rmat=0.0d0

      ! n=blocksize lowest eigenpairs of the Rayleigh matrix
      allocate(reigvec(maxvec,blocksize))
      allocate(reigval(blocksize))
      reigvec=0.0d0
      reigval=0.0d0

      ! Ritz vectors
      allocate(ritzvec(matdim,blocksize))
      ritzvec=0.0d0

      ! Residual vectors
      allocate(res(matdim,blocksize))
      res=0.0d0

      ! Norms of the residual vectors
      allocate(norm(blocksize))
      norm=0.0d0
      
      return

    end subroutine davinitialise

!#######################################################################

    subroutine initvec_adc1

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer                              :: iadc1,dim1,i
      integer, dimension(:), allocatable   :: indx1
      real(d), dimension(:,:), allocatable :: vec1

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

!-----------------------------------------------------------------------
! Set the initial Davidson vectors
!-----------------------------------------------------------------------
      do i=1,blocksize
         vmat(1:dim1,i)=vec1(:,i)
      enddo

!-----------------------------------------------------------------------
! Close the ADC(1) eigenvector file
!-----------------------------------------------------------------------
      close(iadc1)

      return

    end subroutine initvec_adc1

!#######################################################################

    subroutine initvec_ondiag(matdim)

      use iomod, only: freeunit
      use constants
      use parameters
      use misc, only: dsortindxa1

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: iham,i,k
      integer, dimension(:), allocatable :: indx_hii
      real(d), dimension(:), allocatable :: hii
      character(len=70)                  :: filename

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hii(matdim))
      allocate(indx_hii(matdim))

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      call freeunit(iham)
      
      if (hamflag.eq.'i') then
         filename='SCRATCH/hmlt.diai'
      else if (hamflag.eq.'f') then
         filename='SCRATCH/hmlt.diac'
      endif

      open(iham,file=filename,status='old',access='sequential',&
        form='unformatted') 

!-----------------------------------------------------------------------
! Read the on-diagonal Hamiltonian matrix elements
!-----------------------------------------------------------------------
      read(iham) maxbl,nrec
      read(iham) hii

!-----------------------------------------------------------------------
! Determine the indices of the on-diagonal elements with the smallest
! absolute values
!-----------------------------------------------------------------------
      hii=abs(hii)   
      call dsortindxa1('A',matdim,hii,indx_hii)

!-----------------------------------------------------------------------
! Set the initial vectors
!-----------------------------------------------------------------------
      do i=1,blocksize
         k=indx_hii(i)
         vmat(k,i)=1.0d0
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(hii)
      deallocate(indx_hii)

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(iham)

      return

    end subroutine initvec_ondiag

!#######################################################################

    subroutine isincore(matdim,noffd)

      use constants
      use parameters, only: maxmem

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      real(d)               :: mem

      mem=0.0d0
      
      ! On-diagonal Hamiltonian matrix elements
      mem=mem+8.0d0*matdim/1024.0d0**2

      ! Non-zero off-diagonal Hamiltonian matrix element values
      mem=mem+8.0d0*noffd/1024.0d0**2

      ! Indices of the non-zero off-diagonal Hamiltonian matrix elements
      mem=mem+8.0d0*noffd/1024.0d0**2

      ! Subspace vectors
      mem=mem+8.0d0*matdim*maxvec/1024.0d0**2

      ! Matrix-vector products
      mem=mem+8.0d0*matdim*maxvec/1024.0d0**2
      
      ! Ritz vectors
      mem=mem+8.0d0*matdim*blocksize/1024.0d0**2

      ! Residual vectors
      mem=mem+8.0d0*matdim*blocksize/1024.0d0**2

      ! Work arrays used in the subspace expansion routines
      if (ipre.eq.1) then
         ! DPR
         mem=mem+8.0d0*matdim*blocksize/1024.0d0**2
      else if (ipre.eq.2) then
         ! Olsen
         mem=mem+2.0d0*8.0d0*matdim*blocksize/1024.0d0**2
      endif

      if (mem.lt.maxmem) then
         lincore=.true.
      else
         lincore=.false.
      endif

      return

    end subroutine isincore

!#######################################################################

    subroutine rdham_on(matdim)

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer, intent(in) :: matdim
      integer             :: iham
      character(len=70)   :: filename

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hii(matdim))

!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
      call freeunit(iham)

      if (hamflag.eq.'i') then
         filename='SCRATCH/hmlt.diai'
      else if (hamflag.eq.'f') then
         filename='SCRATCH/hmlt.diac'
      endif

      open(iham,file=filename,status='old',access='sequential',&
           form='unformatted')

      read(iham) maxbl,nrec
      read(iham) hii

      close(iham)

      return

    end subroutine rdham_on

!#######################################################################

    subroutine rdham_off(noffd)

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer*8, intent(in)              :: noffd
      integer                            :: iham,count,k,nlim
      integer, dimension(:), allocatable :: itmp1,itmp2
      real(d), dimension(:), allocatable :: ftmp
      character(len=70)                  :: filename

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hij(noffd))
      allocate(indxi(noffd))
      allocate(indxj(noffd))

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
      allocate(ftmp(maxbl))
      allocate(itmp1(maxbl))
      allocate(itmp2(maxbl))

      if (hamflag.eq.'i') then
         filename='SCRATCH/hmlt.offi'
      else if (hamflag.eq.'f') then
         filename='SCRATCH/hmlt.offc'
      endif

      open(iham,file=filename,status='old',access='sequential',&
           form='unformatted')

      count=0
      do k=1,nrec
         read(iham) ftmp(:),itmp1(:),itmp2(:),nlim
         hij(count+1:count+nlim)=ftmp(1:nlim)
         indxi(count+1:count+nlim)=itmp1(1:nlim)
         indxj(count+1:count+nlim)=itmp2(1:nlim)
         count=count+nlim
      enddo

      deallocate(ftmp,itmp1,itmp2)

      close(iham)

      return

    end subroutine rdham_off

!#######################################################################

    subroutine run_block_davidson(matdim,noffd)

      use iomod
      use constants
      use channels
      use parameters

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: k

!-----------------------------------------------------------------------
! Initialise currdim (the current dimension of the subspace)
!-----------------------------------------------------------------------
      currdim=blocksize

!-----------------------------------------------------------------------
! Perform the Davidson iterations
!-----------------------------------------------------------------------
      do k=1,niter
         
         ! Calculate the matrix-vector product
         call hxvec(matdim,noffd)

         ! Caulculate the Rayleigh matrix
         call calcrmat(matdim)
         
         ! Diagonalise the Rayleigh matrix
         call diagrmat

         ! Calculate the Ritz vectors
         call calcritzvec(matdim)

         ! Calculate the residuals
         call calcres(matdim)

         ! Output progress
         call wrtable(k)

         ! Exit if we have converged all roots
         if (nconv.eq.nstates) then
            write(ilog,'(/,2x,a,/)') 'All roots converged'
            exit
         endif

         ! Expand the subspace
         call subspace_expansion(matdim)

      enddo

!-----------------------------------------------------------------------
! Die here if we haven't converged all eigenpairs
!-----------------------------------------------------------------------
      if (nconv.ne.nstates) then
         errmsg='Not all vectors have converged...'
         call error_control
      endif

      return

    end subroutine run_block_davidson

!#######################################################################
    
    subroutine wrtable(k)
      
      use constants
      use channels
      use parameters

      implicit none

      integer          :: k,i,j
      character(len=1) :: aconv

!-----------------------------------------------------------------------
! Table header
!-----------------------------------------------------------------------
      if (k.eq.1) then
         write(ilog,'(/,53a)') ('*',j=1,53)
         write(ilog,'(4(a,6x))') &
              'Iteration','Energies','Residuals','Converged'
         write(ilog,'(53a)') ('*',j=1,53)
      endif

!-----------------------------------------------------------------------
! Information from the current iteration
!-----------------------------------------------------------------------
      write(ilog,'(/)')
      do i=1,nstates
         if (norm(i).lt.tol) then
            aconv='y'
         else
            aconv='n'
         endif
         
         if (i.eq.1) then
            write(ilog,'(i4,10x,F12.7,3x,E13.7,2x,a1)') &
                 k,reigval(i)*eh2ev,norm(i),aconv
         else
            write(ilog,'(14x,F12.7,3x,E13.7,2x,a1)') &
                 reigval(i)*eh2ev,norm(i),aconv
         endif
      enddo

      return

    end subroutine wrtable

!#######################################################################

    subroutine hxvec(matdim,noffd)

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd

      if (lincore) then
         call hxvec_incore(matdim,noffd)
      else
          call hxvec_ext(matdim,noffd)
      endif

      return

    end subroutine hxvec

!#######################################################################

    subroutine hxvec_incore(matdim,noffd)

      implicit none
      
      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: m,n,k

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      wmat=0.0d0
      !$omp parallel do private(m,n) shared(wmat,hii,vmat)
      do n=1,currdim
         do m=1,matdim
            wmat(m,n)=hii(m)*vmat(m,n)
         enddo
      enddo
      !$omp end parallel do

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      !$omp parallel do private(k,n) shared(wmat,hij,vmat,indxi,indxj)
      do n=1,currdim
         do k=1,noffd
            wmat(indxi(k),n)=wmat(indxi(k),n)+hij(k)*vmat(indxj(k),n)
            wmat(indxj(k),n)=wmat(indxj(k),n)+hij(k)*vmat(indxi(k),n)
         enddo
      enddo
      !$omp end parallel do

      return

    end subroutine hxvec_incore

!#######################################################################

    subroutine hxvec_ext(matdim,noffd)

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: iham
      integer               :: nlim,i,j,k,l,m,n
      character(len=70)     :: filename

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      wmat=0.0d0
      !$omp parallel do private(m,n) shared(wmat,hii,vmat)
      do n=1,currdim
         do m=1,matdim
            wmat(m,n)=hii(m)*vmat(m,n)
         enddo
      enddo
      !$omp end parallel do

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
      
      if (hamflag.eq.'i') then
         filename='SCRATCH/hmlt.offi'
      else if (hamflag.eq.'f') then
         filename='SCRATCH/hmlt.offc'
      endif

      open(iham,file=filename,status='old',access='sequential',&
           form='unformatted')

      do k=1,nrec
         read(iham) hij(:),indxi(:),indxj(:),nlim
         !$omp parallel do private(l,n) shared(wmat,hij,vmat,indxi,indxj)
         do n=1,currdim
            do l=1,nlim               
               wmat(indxi(l),n)=wmat(indxi(l),n)+hij(l)*vmat(indxj(l),n)
               wmat(indxj(l),n)=wmat(indxj(l),n)+hij(l)*vmat(indxi(l),n)
            enddo
         enddo
         !$omp end parallel do
      enddo

      close(iham)

      deallocate(hij,indxi,indxj)

      return

    end subroutine hxvec_ext

!#######################################################################

    subroutine calcrmat(matdim)
      
      implicit none

      integer, intent(in) :: matdim
      
      rmat=0.0d0
      call dgemm('T','N',currdim,currdim,matdim,1.0d0,&
           vmat(:,1:currdim),matdim,wmat(:,1:currdim),matdim,0.0d0,&
           rmat(1:currdim,1:currdim),currdim)
      
      return

    end subroutine calcrmat

!#######################################################################

    subroutine diagrmat
      
      use iomod
      use constants

      implicit none

      integer                             :: e2,i
      real(d), dimension(3*currdim)       :: work
      real(d)                             :: error
      real(d), dimension(currdim)         :: val
      real(d), dimension(currdim,currdim) :: vec

!-----------------------------------------------------------------------
! Diagonalise the Rayleigh matrix
!-----------------------------------------------------------------------
      error=0
      e2=3*currdim
      vec=rmat(1:currdim,1:currdim)

      call dsyev('V','U',currdim,vec,currdim,val,work,e2,error)

      if (error.ne.0) then
         errmsg='Diagonalisation of the Rayleigh matrix in &
              subroutine diagrmat failed'
         call error_control
      endif

!-----------------------------------------------------------------------
! Save the n=blocksize lowest eigenpairs to be used in the calculation
! of the Ritz vectors and residuals
!-----------------------------------------------------------------------
      reigvec=0.0d0
      reigvec(1:currdim,1:blocksize)=vec(1:currdim,1:blocksize)

      reigval=0.0d0
      reigval(1:blocksize)=val(1:blocksize)

      return

    end subroutine diagrmat

!#######################################################################

    subroutine calcritzvec(matdim)

      implicit none

      integer, intent(in) :: matdim
      
      call dgemm('N','N',matdim,blocksize,currdim,1.0d0,&
           vmat(1:matdim,1:currdim),matdim,&
           reigvec(1:currdim,1:blocksize),currdim,0.0d0,&
           ritzvec(1:matdim,1:blocksize),matdim)
      
      return

    end subroutine calcritzvec

!#######################################################################

    subroutine calcres(matdim)

      use constants
      
      implicit none

      integer, intent(in) :: matdim
      integer             :: i
      real(d)             :: ddot

      external ddot

!-----------------------------------------------------------------------
! Residual vectors: r_i = lambda_i * x_i - W * y_i
!-----------------------------------------------------------------------
! r_i       ith residual vector
!
! lambda_i  ith eigenvalue of the Rayleigh matrix
!
! x_i       ith Ritz vector
!
! W         = H * V (Hamiltonian multiplied against the matrix of
!                   subspace vectors)
!
! y_i      ith eigenvector of the Rayleigh matrix
!-----------------------------------------------------------------------
      ! -W * y_i
      call dgemm('N','N',matdim,blocksize,currdim,-1.0d0,&
           wmat(1:matdim,1:currdim),matdim,&
           reigvec(1:currdim,1:blocksize),currdim,0.0d0,&
           res(1:matdim,1:blocksize),matdim)
      
      ! lambda_i * x_i -W * y_i
      do i=1,blocksize
         res(:,i)=res(:,i)+reigval(i)*ritzvec(:,i)
      enddo

!-----------------------------------------------------------------------
! Norms of the residual vectors
!-----------------------------------------------------------------------
      do i=1,blocksize
         norm(i)=ddot(matdim,res(:,i),1,res(:,i),1)
         norm(i)=dsqrt(norm(i))
      enddo

!-----------------------------------------------------------------------
! Keep track of the no. converged roots
!-----------------------------------------------------------------------
      nconv=0
      do i=1,nstates
         if (norm(i).lt.tol) nconv=nconv+1
      enddo

      return

    end subroutine calcres

!#######################################################################

    subroutine subspace_expansion(matdim)

      use channels
      use iomod

      implicit none

      integer, intent(in) :: matdim
      logical             :: lcollapse

!-----------------------------------------------------------------------
! ipre = 1 <-> diagonal preconditioned residue
! ipre = 2 <-> Olsen's preconditioner
!-----------------------------------------------------------------------
      ! Determine whether or not we need to collapse the subspace
      if (currdim.le.maxvec-blocksize) then
         lcollapse=.false.
      else
         lcollapse=.true.
         write(ilog,'(/,2x,a)') 'Collapsing the subspace'
      endif

      ! Calculate the new subspace vectors
      if (ipre.eq.1) then
         call dpr(matdim,lcollapse)
      else if (ipre.eq.2) then
         call olsen(matdim,lcollapse)
      endif
      
      ! Orthogonalise the subspace vectors
      call qrortho(matdim,lcollapse)

      ! Update the dimension of the subspace
      if (lcollapse) then
         currdim=2*blocksize
      else
         currdim=currdim+blocksize
      endif

      return

    end subroutine subspace_expansion

!#######################################################################

    subroutine dpr(matdim,lcollapse)

      use iomod
      use constants

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: i,j,ilbl,ilast
      real(d), dimension(:), allocatable :: tmpvec
      logical                            :: lcollapse

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(tmpvec(matdim))

!-----------------------------------------------------------------------
! Calculate the new subspace vectors
!-----------------------------------------------------------------------
      if (lcollapse) then
         ! Collapse of the subspace
         vmat=0.0d0
         vmat(:,1:blocksize)=ritzvec(:,1:blocksize)
         ilast=blocksize
      else
         ! Expansion of the subspace
         ilast=currdim
      endif

      ! Loop over new subspace vectors
      do i=1,blocksize

         ! Index of the next vector
         ilbl=ilast+i
         
         ! Calculate the next vector
         tmpvec=0.0d0
         do j=1,matdim
            tmpvec(j)=1.0d0/(reigval(i)-hii(j))
         enddo
         do j=1,matdim
            vmat(j,ilbl)=res(j,i)*tmpvec(j)
         enddo

      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(tmpvec)

      return

    end subroutine dpr

!#######################################################################

    subroutine olsen(matdim,lcollapse)

      use iomod
      use constants

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: i,j,ilbl,ilast,info
      real(d), dimension(:), allocatable :: tmpvec,cdiag
      real(d)                            :: alpha,xz,xy,ddot
      logical                            :: lcollapse

      external ddot

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(tmpvec(matdim))
      allocate(cdiag(matdim))

!-----------------------------------------------------------------------
! Calculate the new subspace vectors
!-----------------------------------------------------------------------
      if (lcollapse) then
         ! Collapse of the subspace
         vmat=0.0d0
         vmat(:,1:blocksize)=ritzvec(:,1:blocksize)
         ilast=blocksize
      else
         ! Expansion of the subspace
         ilast=currdim
      endif

      ! Loop over new subspace vectors
      do i=1,blocksize
         
         ! Diagonal of the C-matrix
         do j=1,matdim
            cdiag(j)=1.0d0/(hii(j)-reigval(i))
         enddo
         
         ! x_i * z_i
         do j=1,matdim
            tmpvec(j)=res(j,i)*cdiag(j)
         enddo
         xz=ddot(matdim,ritzvec(:,i),1,tmpvec(:),1)
         
         ! x_i * y_i
         do j=1,matdim
            tmpvec(j)=ritzvec(j,i)*cdiag(j)
         enddo
         xy=ddot(matdim,ritzvec(:,i),1,tmpvec(:),1)
         
         ! alpha
         alpha=xz/xy
         
         ! New subspace vector
         ilbl=ilast+i
         do j=1,matdim
            vmat(j,ilbl)=alpha*cdiag(j)*ritzvec(j,i)-cdiag(j)*res(j,i)
         enddo
         
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(tmpvec)
      deallocate(cdiag)


      return

    end subroutine olsen

!#######################################################################

    subroutine qrortho(matdim,lcollapse)

      use iomod
      use constants

      implicit none

      integer, intent(in)                :: matdim
      integer                            :: n,info
      real(d), dimension(:), allocatable :: tau,work
      logical                            :: lcollapse

!-----------------------------------------------------------------------
! Orthogonalisation of the subspace vectors via the QR factorization
! of the matrix of subspace vectors
!-----------------------------------------------------------------------
      if (lcollapse) then
         n=2*blocksize
      else
         n=currdim+blocksize
      endif

      allocate(tau(n))
      allocate(work(n))

      call dgeqrf(matdim,n,vmat(:,1:n),matdim,tau,work,n,info)
      if (info.ne.0) then
         errmsg='dqerf failed in subroutine qrortho'
         call error_control
      endif
      
      call dorgqr(matdim,n,n,vmat(:,1:n),matdim,tau,work,n,info)
      if (info.ne.0) then
         errmsg='dorgqr failed in subroutine qrortho'
         call error_control
      endif
      
      deallocate(tau)
      deallocate(work)

      return

    end subroutine qrortho

!#######################################################################

    subroutine wreigenpairs

      use iomod
      
      implicit none

      integer :: unit,i

!-----------------------------------------------------------------------
! Open the Davidson vector file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit=unit,file=vecfile,status='unknown',&
           access='sequential',form='unformatted')

!-----------------------------------------------------------------------
! Write the eigenpairs to file
!-----------------------------------------------------------------------
      do i=1,nstates
         write(unit) i,reigval(i),ritzvec(:,i)
      enddo

!-----------------------------------------------------------------------
! Close the Davidson vector file
!-----------------------------------------------------------------------
      close(unit)
      
      return
      
    end subroutine wreigenpairs

!#######################################################################

    subroutine davfinalise

      implicit none

      deallocate(vmat)
      deallocate(wmat)
      deallocate(rmat)
      deallocate(reigvec)
      deallocate(reigval)
      deallocate(ritzvec)
      deallocate(res)
      deallocate(norm)
      deallocate(hii)
      if (allocated(hij)) deallocate(hij)
      if (allocated(indxi)) deallocate(indxi)
      if (allocated(indxj)) deallocate(indxj)
      
      return
      
    end subroutine davfinalise
      
!#######################################################################

  end module block_davidson
