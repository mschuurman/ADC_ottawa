!#######################################################################
! block_davidson: Routine for the calculation of the lowest-lying
!                 eigenvectors of the ADC Hamiltonian matrix using
!                 the block Davidson method.
!#######################################################################
  module block_davidson

    use constants

    implicit none

    save

    integer                              :: blocksize,nstates,maxvec,&
                                            niter,currdim,ipre,nconv,&
                                            nconv_prev,nrec,maxbl,&
                                            blocksize_curr,maxvec_curr,&
                                            nstates_curr,nmult,subdim
    integer, dimension(:), allocatable   :: indxi,indxj
    real(d), dimension(:), allocatable   :: hii,hij
    real(d), dimension(:,:), allocatable :: vmat,wmat,rmat,ritzvec,&
                                            res,reigvec
    real(d), dimension(:), allocatable   :: reigval,norm
    real(d)                              :: tol
    character(len=36)                    :: vecfile
    logical                              :: lincore,lrdadc1,lsub,&
                                            ldeflate

  contains

!#######################################################################
    
    subroutine davdiag_block(matdim,noffd)

      use channels
      use constants
      use parameters
      use timingmod

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: k
      real(d)               :: tw1,tw2,tc1,tc2
      character(len=120)    :: atmp
      
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
      call times(tw1,tc1)

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
! Set the initial vectors
! -----------------------------------------------------------------------
      ! Determine whether the initial vectors are to be constructed
      ! from the ADC(1) eigenvectors, the eigenvectors of the
      ! Hamiltonian represented in a subspace of ISs, or as single ISs
      if (hamflag.eq.'i'.and.ladc1guess) then
         lrdadc1=.true.
      else if (hamflag.eq.'f'.and.ladc1guess_f) then
         lrdadc1=.true.
      else
         lrdadc1=.false.
      endif
      if (hamflag.eq.'i'.and.lsubdiag) then
         lsub=.true.
      else if (hamflag.eq.'f'.and.lsubdiag_f) then
         lsub=.true.
      endif

      if (lrdadc1) then
         ! Construct the initial vectors from the ADC(1) eigenvectors
         call initvec_adc1
      else if (lsub) then
         ! Construct the initial vectors from the eigenvectors of the
         ! Hamiltonian represented in a subspace of ISs
         call initvec_subdiag(matdim,noffd)
      else
         ! Use a single IS for each initial vector
         call initvec_ondiag(matdim)
      endif

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

!-----------------------------------------------------------------------    
! Output timings and the no. matrix-vector multiplications
!-----------------------------------------------------------------------    
      write(ilog,'(/,a,1x,i4)') 'No. matrix-vector multiplications:',&
           nmult

      call times(tw2,tc2)
      write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
      
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
         ldeflate=ldfl
         subdim=guessdim
      else if (hamflag.eq.'f') then
         blocksize=dmain_f
         nstates=davstates_f
         niter=maxiter_f
         ipre=precon_f
         tol=davtol_f
         vecfile=davname_f
         maxvec=maxsubdim_f
         ldeflate=ldfl_f
         subdim=guessdim_f
      endif

!-----------------------------------------------------------------------
! Number of matrix-vector multiplications
!-----------------------------------------------------------------------
      nmult=0

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

    subroutine initvec_subdiag(matdim,noffd)

      use parameters
      use constants
      use misc, only: dsortindxa1
      use iomod
      use channels

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer, dimension(:), allocatable   :: full2sub,sub2full,indxhii
      integer*8                            :: k,i,j
      integer                              :: i1,j1
      integer                              :: e2,error,iham,nlim,l
      real(d), dimension(:,:), allocatable :: hsub
      real(d), dimension(:), allocatable   :: subeig,work
      character(len=70)                    :: filename

!-----------------------------------------------------------------------
! Subspace dimension check
!-----------------------------------------------------------------------
      if (subdim.gt.matdim) subdim=matdim

!-----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix elements in order of
! ascending value
!-----------------------------------------------------------------------
      allocate(indxhii(matdim))
      call dsortindxa1('A',matdim,hii,indxhii)

!-----------------------------------------------------------------------
! Ensure that the subdim'th IS is not degenerate with subdim+1'th IS,
! and if it is, increase subdim accordingly
!-----------------------------------------------------------------------
5     continue
      if (subdim.lt.matdim) then
         if (abs(hii(indxhii(subdim))-hii(indxhii(subdim+1))).lt.1e-6_d) then
            subdim=subdim+1
            goto 5
         endif
      endif

!-----------------------------------------------------------------------
! Allocate the subspace-associated arrays
!-----------------------------------------------------------------------
      allocate(full2sub(matdim))
      allocate(sub2full(subdim))
      allocate(hsub(subdim,subdim))
      allocate(subeig(subdim))
      allocate(work(3*subdim))

!-----------------------------------------------------------------------
! Set the full space-to-subsace mappings
!-----------------------------------------------------------------------
      full2sub=0
      do i=1,subdim
         k=indxhii(i)
         sub2full(i)=k
         full2sub(k)=i
      enddo
      
!-----------------------------------------------------------------------
! Construct the Hamiltonian matrix in the subspace
!-----------------------------------------------------------------------
      hsub=0.0d0

      ! On-diagonal elements
      do i=1,subdim
         k=sub2full(i)
         hsub(i,i)=hii(k)
      enddo

      ! Off-diagonal elements
      if (lincore) then
         ! Loop over all off-diagonal elements of the full space
         ! Hamiltonian
         do k=1,noffd
            ! Indices of the current off-diagonal element of the
            ! full space Hamiltonian
            i=indxi(k)
            j=indxj(k)
            ! If both indices correspond to subspace ISs, then
            ! add the element to subspace Hamiltonian
            if (full2sub(i).ne.0.and.full2sub(j).ne.0) then
               i1=full2sub(i)
               j1=full2sub(j)
               hsub(i1,j1)=hij(k)               
               hsub(j1,i1)=hsub(i1,j1)
            endif
         enddo
      else
         ! Open the off-diagonal element file
         call freeunit(iham)
         if (hamflag.eq.'i') then
            filename='SCRATCH/hmlt.offi'
         else if (hamflag.eq.'f') then
            filename='SCRATCH/hmlt.offc'
         endif
         open(iham,file=filename,status='old',access='sequential',&
              form='unformatted')
         ! Allocate arrays
         allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
         ! Loop over records
         do k=1,nrec
            read(iham) hij(:),indxi(:),indxj(:),nlim
            ! Loop over the non-zero elements of the full space
            ! Hamiltonian in the current record
            do l=1,nlim
               ! Indices of the current off-diagonal element of the
               ! full space Hamiltonian
               i=indxi(l)
               j=indxj(l)
               ! If both indices correspond to subspace ISs, then
               ! add the element to subspace Hamiltonian
               if (full2sub(i).ne.0.and.full2sub(j).ne.0) then
                  i1=full2sub(i)
                  j1=full2sub(j)
                  hsub(i1,j1)=hij(l)
                  hsub(j1,i1)=hsub(i1,j1)
               endif
            enddo
         enddo
         ! Close the off-diagonal element file
         close(iham)
         ! Deallocate arrays
         deallocate(hij,indxi,indxj)
      endif

!-----------------------------------------------------------------------
! Diagonalise the subspace Hamiltonian
!-----------------------------------------------------------------------
      e2=3*subdim
      call dsyev('V','U',subdim,hsub,subdim,subeig,work,e2,error)

      if (error.ne.0) then
         errmsg='The diagonalisation of the subspace Hamiltonian failed.'
         call error_control
      endif

!-----------------------------------------------------------------------
! Construct the initial vectors from the subspace vectors.
! Note that after calling dsyev, hsub now holds the eigenvectors of
! the subspace Hamiltonian.
!-----------------------------------------------------------------------
      vmat=0.0d0
      do i=1,blocksize
         do j=1,subdim
            k=sub2full(j)
            vmat(k,i)=hsub(j,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(indxhii)
      deallocate(full2sub)
      deallocate(sub2full)
      deallocate(hsub)
      deallocate(subeig)
      deallocate(work)

      return

    end subroutine initvec_subdiag

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

      call freeunit(iham)
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
! Initialisation
!-----------------------------------------------------------------------
! currdim:        the current dimension of the subspace
!
! blocksize_curr: the current blocksize
!
! maxvec_curr:    the current maximum subspace dimension
!
! nstates_curr:   the current no. of states that we are solving for,
!                 i.e., the current no. of unconverged roots
!
! nconv:          the no. of converged roots
!-----------------------------------------------------------------------
! Note that the blocksize, maximum subspace dimension and nstates will
! only change throughout the Davidson iterations if the converged
! vectors are removed from from the subspace, i.e., if ldeflate=.true.
!-----------------------------------------------------------------------
      currdim=blocksize
      blocksize_curr=blocksize
      maxvec_curr=maxvec
      nstates_curr=nstates
      nconv=0

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

         ! Keep track of the no. of roots converged so far
         nconv_prev=nconv
         
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
         write(ilog,'(53a)') ('*',j=1,53)
         write(ilog,'(4(a,6x))') &
              'Iteration','Energies','Residuals','Converged'
         write(ilog,'(53a)') ('*',j=1,53)
      endif

!-----------------------------------------------------------------------
! Information from the current iteration
!-----------------------------------------------------------------------
      write(ilog,*)
      do i=1,nstates_curr
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
! Update the no. matrix-vector multiplications
!-----------------------------------------------------------------------
      nmult=nmult+currdim

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
! Update the no. matrix-vector multiplications
!-----------------------------------------------------------------------
      nmult=nmult+currdim

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
! Save the n=blocksize_curr lowest eigenpairs to be used in the
! calculation of the Ritz vectors and residuals
!-----------------------------------------------------------------------
      reigvec(1:currdim,1:blocksize_curr)=vec(1:currdim,1:blocksize_curr)
      reigval(1:blocksize_curr)=val(1:blocksize_curr)

      return

    end subroutine diagrmat

!#######################################################################

    subroutine calcritzvec(matdim)

      implicit none

      integer, intent(in) :: matdim
      
      call dgemm('N','N',matdim,blocksize_curr,currdim,1.0d0,&
           vmat(1:matdim,1:currdim),matdim,&
           reigvec(1:currdim,1:blocksize_curr),currdim,0.0d0,&
           ritzvec(1:matdim,1:blocksize_curr),matdim)
      
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
      call dgemm('N','N',matdim,blocksize_curr,currdim,-1.0d0,&
           wmat(1:matdim,1:currdim),matdim,&
           reigvec(1:currdim,1:blocksize_curr),currdim,0.0d0,&
           res(1:matdim,1:blocksize_curr),matdim)
      
      ! lambda_i * x_i -W * y_i
      do i=1,blocksize_curr
         res(:,i)=res(:,i)+reigval(i)*ritzvec(:,i)
      enddo

!-----------------------------------------------------------------------
! Norms of the residual vectors
!-----------------------------------------------------------------------
      do i=1,blocksize_curr
         norm(i)=ddot(matdim,res(:,i),1,res(:,i),1)
         norm(i)=sqrt(norm(i))
      enddo

!-----------------------------------------------------------------------
! Keep track of the no. converged roots
!-----------------------------------------------------------------------
      if (.not.ldeflate) nconv=0

      do i=1,nstates_curr
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
! Removal of converged vectors from the subspace
!-----------------------------------------------------------------------
      if (ldeflate.and.nconv.gt.nconv_prev) call deflate(matdim)
      
!-----------------------------------------------------------------------
! ipre = 1 <-> diagonal preconditioned residue
! ipre = 2 <-> Olsen's preconditioner
!-----------------------------------------------------------------------
      ! Determine whether or not we need to collapse the subspace
      if (currdim.le.maxvec_curr-blocksize_curr) then
         lcollapse=.false.
      else
         lcollapse=.true.
      endif

      ! Force a collapse of the subspace if a deflation has occurred
      if (ldeflate.and.nconv.gt.nconv_prev) lcollapse=.true.

!      if (lcollapse) write(ilog,'(/,2x,a)') 'Collapsing the subspace'
      
      ! Calculate the new subspace vectors
      if (ipre.eq.1) then
         call dpr(matdim,lcollapse)
      else if (ipre.eq.2) then
         call olsen(matdim,lcollapse)
      endif
      
      ! Project the subspace onto the space orthogonal to the
      ! that spanned by the converged vectors
      if (ldeflate.and.nconv.gt.0) call subspace_projection(matdim,lcollapse)

      ! Orthogonalise the subspace vectors
      call qrortho(matdim,lcollapse)

      ! Update the dimension of the subspace
      if (lcollapse) then
         currdim=2*blocksize_curr
      else
         currdim=currdim+blocksize_curr
      endif

      return

    end subroutine subspace_expansion

!#######################################################################

    subroutine deflate(matdim)

      use channels
      use iomod
      
      implicit none

      integer, intent(in)                :: matdim
      integer                            :: indx,i,j,count
      integer, dimension(nstates)        :: convmap
      real(d), dimension(:), allocatable :: swapvec
      real(d)                            :: swapval
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(swapvec(matdim))
      
!-----------------------------------------------------------------------
! Rearrange the arrays holding the Ritz vectors, the Ritz values and the
! residual vectors s.t. the elements associated with the converged
! vectors are moved to the ends of the arrays
!-----------------------------------------------------------------------
      ! Indices of the converged vectors
      count=0
      convmap=0
      do i=1,blocksize_curr
         if (norm(i).lt.tol) then
            count=count+1
            convmap(count)=i
         endif
      enddo

      ! Rearranegment of the arrays ritzvec, reigval and res
      count=0
      do i=1,nconv-nconv_prev

         count=count+1
         
         ! Index of the next column into which the next converged vector
         ! is to be moved
         indx=blocksize+1-nconv_prev-i

         ! Ritz vector
         swapvec=ritzvec(:,indx)
         ritzvec(:,indx)=ritzvec(:,convmap(count))
         ritzvec(:,convmap(count))=swapvec

         ! Ritz value
         swapval=reigval(indx)
         reigval(indx)=reigval(convmap(count))
         reigval(convmap(count))=swapval

         ! Residual vector
         swapvec=res(:,indx)
         res(:,indx)=res(:,convmap(count))
         res(:,convmap(count))=swapvec
         
      enddo      
      
!-----------------------------------------------------------------------
! Update dimensions
!-----------------------------------------------------------------------
      blocksize_curr=blocksize_curr-(nconv-nconv_prev)
      maxvec_curr=maxvec_curr-(nconv-nconv_prev)
      nstates_curr=nstates_curr-(nconv-nconv_prev)
      currdim=currdim-(nconv-nconv_prev)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(swapvec)

      return
      
    end subroutine deflate

!#######################################################################

    subroutine subspace_projection(matdim,lcollapse)

      use iomod
      use channels
      
      implicit none

      integer, intent(in)                  :: matdim
      integer                              :: i,j,lower,upper,indx,&
                                              nsubvec
      real(d), dimension(:,:), allocatable :: overlap
      logical                              :: lcollapse

!-----------------------------------------------------------------------
! Set the lower and upper indices on the unconverged subspace vectors
! that need to be orthogonalised against the converged vectors
!-----------------------------------------------------------------------      
      if (lcollapse) then
         lower=1
         upper=2*blocksize_curr
      else
         lower=1
         upper=currdim+blocksize_curr
      endif

      nsubvec=upper-lower+1
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(overlap(nconv,nsubvec))
      
!-----------------------------------------------------------------------
! Calculate the matrix of overlaps between the converged vectors and the
! unconverged subspace vectors
!-----------------------------------------------------------------------
      call dgemm('T','N',nconv,nsubvec,matdim,1.0d0,&
           ritzvec(:,blocksize-nconv+1:blocksize),matdim,&
           vmat(:,lower:upper),matdim,0.0d0,overlap,nconv)

!-----------------------------------------------------------------------
! Orthogonalise the unconverged subspace vectors againts to the
! converged vectors
!-----------------------------------------------------------------------
      do i=lower,upper
         do j=1,nconv
            indx=blocksize-nconv+j
            vmat(:,i)=vmat(:,i)-ritzvec(:,indx)*overlap(j,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(overlap)
      
      return
      
    end subroutine subspace_projection
    
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
         vmat(:,1:blocksize_curr)=ritzvec(:,1:blocksize_curr)
         ilast=blocksize_curr
      else
         ! Expansion of the subspace
         ilast=currdim
      endif

      ! Loop over new subspace vectors
      do i=1,blocksize_curr

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
         vmat(:,1:blocksize_curr)=ritzvec(:,1:blocksize_curr)
         ilast=blocksize_curr
      else
         ! Expansion of the subspace
         ilast=currdim
      endif

      ! Loop over new subspace vectors
      do i=1,blocksize_curr
         
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
         n=2*blocksize_curr
      else
         n=currdim+blocksize_curr
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

      use constants
      use iomod
      use misc, only: dsortindxa1
      
      implicit none

      integer                     :: unit,i
      integer, dimension(nstates) :: indx

!-----------------------------------------------------------------------
! Rearrangement of the ritzvec and reigval array such that the
! converged eigenpairs come first
!
! No. of converged eigenpairs "in storage": nconv_prev
! No. converged eigenpairs not "in storage": nstates_curr      
!
! N.B. This is only relevant if the subspace deflation was used
!-----------------------------------------------------------------------
      if (ldeflate) then

         reigval(nstates_curr+1:nstates)=&
              reigval(blocksize-nconv_prev+1:blocksize)

         ritzvec(:,nstates_curr+1:nstates)=&
              ritzvec(:,blocksize-nconv_prev+1:blocksize)
         
      endif
         
!-----------------------------------------------------------------------
! Indices of the eigenpairs in order of increasing energy
!
! Note that this is only relevant if subspace deflation was used, but
! makes no difference if not
!-----------------------------------------------------------------------
      if (nstates.gt.1) then
         call dsortindxa1('A',nstates,reigval(1:nstates),indx)
      else
         indx(1)=1
      endif

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
         write(unit) i,reigval(indx(i)),ritzvec(:,indx(i))
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
