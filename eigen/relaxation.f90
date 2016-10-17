!#######################################################################
! Experimental block-relaxation algorithm for the calculation of the
! N lowest eigenpairs of the Hamiltonian matrix.
!
! A collection of N trial vectors are taken as a set of initial
! time-dependent electronic wavefunctions and propagated in negative
! imaginary time using the short iterative Lanczos algorithm.
!
! It is not clear yet whether this will be more effective than the
! block Davidson method.
!
! One possibility is to run the relaxation calculation for one or two
! iterations and then use the resulting vectors as guesses for a
! block Davidson calculation.
!#######################################################################

  module relaxmod

    use iomod
    use constants
    use parameters
    use channels
    use timingmod

    save

    integer                              :: maxbl,nrec,nblock,nconv,&
                                            nstates,krydim,niter,&
                                            iortho,algor
    integer, dimension(:), allocatable   :: indxi,indxj,sildim
    real(d), dimension(:), allocatable   :: hii,hij,ener,res,currtime
    real(d), dimension(:,:), allocatable :: vec_old,vec_new,hxvec
    real(d)                              :: step,eps,toler
    character(len=36)                    :: vecfile
    logical, dimension(:), allocatable   :: lconv
    logical                              :: lincore,lrdadc1,lrandom,lsub

    ! Krylov-Liu arrays
    integer                              :: maxvec
    real(d), dimension(:,:), allocatable :: subhmat,subsmat,kryvec,&
                                            vec_conv
    real(d), dimension(:), allocatable   :: knorm
    
  contains

!#######################################################################

    subroutine relaxation(matdim,noffd)

      implicit none

      integer, intent(in)      :: matdim
      integer*8, intent(in)    :: noffd
      integer                  :: k
      real(d)                  :: tw1,tw2,tc1,tc2
      character(len=120)       :: atmp
      
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
      call times(tw1,tc1)

!-----------------------------------------------------------------------
! Write to the log file
!-----------------------------------------------------------------------
      atmp='Block-relaxation in the'
      if (hamflag.eq.'i') then
         atmp=trim(atmp)//' initial space'
      else if (hamflag.eq.'f') then
         atmp=trim(atmp)//' final space'
      endif
      write(ilog,'(/,70a)') ('-',k=1,70)
      write(ilog,'(2x,a)') trim(atmp)
      write(ilog,'(70a,/)') ('-',k=1,70)

!-----------------------------------------------------------------------
! Allocatation and initialisation of arrays
!-----------------------------------------------------------------------
      call initialise(matdim)      

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
!-----------------------------------------------------------------------
      call initvec(matdim,noffd)
      
!-----------------------------------------------------------------------
! Output the initial energies and convergence information
!-----------------------------------------------------------------------
      vec_new=vec_old
      call getener(matdim,noffd)
      call check_conv(matdim,noffd)
      call wrtable(0)

!-----------------------------------------------------------------------
! Perform the block-relaxation calculation
!-----------------------------------------------------------------------
      if (algor.eq.1) then
         ! Traditional SIL algorithm
         call relaxation_sil(matdim,noffd)
      else if (algor.eq.2) then
         ! Basis reuse algorithm of Liu
         call relaxation_liu(matdim,noffd)
      endif

!-----------------------------------------------------------------------
! Exit here if not all states have converged
!-----------------------------------------------------------------------
      if (nconv.lt.nstates) then
         errmsg='Not all wavefunctions have converged. Quitting here.'
         call error_control
      endif

!-----------------------------------------------------------------------
! Write the converged eigenpairs to disk
!-----------------------------------------------------------------------
      call wreigenpairs

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      call finalise

!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
      call times(tw2,tc2)
      write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"

      return
      
    end subroutine relaxation

!#######################################################################

    subroutine relaxation_sil(matdim,noffd)

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: n,k

      ! Loop over timesteps
      do n=1,niter

         ! Propagate the wavefunctions to the next timestep using
         ! the SIL algorithm    
         do k=1,nblock
            if (lconv(k).and.iortho.eq.2) then
               vec_new(:,k)=vec_old(:,k)
            else
               call modified_silstep(k,matdim,noffd,vec_old(:,k),vec_new(:,k))
            endif
         enddo
         
         ! If the unmodified SIL algorithm is being used, then
         ! orthogonalise the propagated wavefunctions amongst 
         ! themselves
         if (iortho.eq.1) call orthowf(matdim)
         
         ! Calculate the energies
         call getener(matdim,noffd)
         
         ! Sort the wavefunctions by energy
         call sortwf(matdim)
         
         ! Check convergence
         call check_conv(matdim,noffd)

         ! Output the energies and convergence information
         call wrtable(n)
         
         ! Exit if all wavefunctions have converged
         if (nconv.ge.nstates) then
            write(ilog,'(/,2x,a,/)') 'All states have converged'
            exit
         endif
            
         ! Set the initial wavefunctions for the next step
         vec_old=vec_new
         
      enddo
      
      return

    end subroutine relaxation_sil

!#######################################################################

    subroutine relaxation_liu(matdim,noffd)

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: maxvecs,s,n

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
      ! Maximum no. of vectors that will be generated per state.
      ! Note that krydim+1 enters here as we generate one more vector
      ! per timestep than is used in forming the variational subspace
      maxvecs=(krydim+1)*niter
      
      ! Representation of the Hamiltonian in the variational subspace
      allocate(subhmat(maxvecs,maxvecs))
      subhmat=0.0d0

      ! Overlaps of the variational subspace vectors
      allocate(subsmat(maxvecs,maxvecs))
      subsmat=0.0d0

      ! Vectors spanning the variational subspace
      allocate(kryvec(matdim,maxvec))
      kryvec=0.0d0

      ! Norms of the Krylov vectors prior to normalisation (needed
      ! for the calculation of the subspace Hamiltonian matrix
      ! elements)
      allocate(knorm(maxvec))
      
      ! Converged states
      allocate(vec_conv(matdim,nblock))
      vec_conv=0.0d0
      
!-----------------------------------------------------------------------
! Perform the relaxation calculations using the SIL-Liu algorithm
!-----------------------------------------------------------------------
      ! Loop over states
      do s=1,nstates

         ! Initialise arrays
         subhmat=0.0d0
         subsmat=0.0d0
         kryvec=0.0d0

         ! Loop over timesteps
         do n=1,niter
         
            ! Generation of the Krylov vectors for the current
            ! timestep
            call get_kryvecs_current(s,n,matdim,noffd)
            
            ! Calculation of the subspace Hamiltonian and overlap
            ! matrix elements
            call subspace_matrices(n)

            ! Calculate |Psi(t+dt)>
            call solve_subspace_tdse(n)
            
            STOP
            
         enddo
            
      enddo


!-----------------------------------------------------------------------
! Finalisation
!-----------------------------------------------------------------------
      deallocate(subhmat)
      deallocate(subsmat)
      deallocate(kryvec)
      deallocate(knorm)
      deallocate(vec_conv)
      
      return
      
    end subroutine relaxation_liu
    
!#######################################################################

    subroutine initialise(matdim)

      implicit none

      integer, intent(in) :: matdim

!-----------------------------------------------------------------------
! Initialisation of parameter values
!-----------------------------------------------------------------------
      if (hamflag.eq.'i') then
         vecfile=davname
         nstates=davstates
         nblock=dmain
         krydim=kdim
         step=stepsize
         eps=davtol
         niter=maxiter
         iortho=rlxortho
         toler=siltol
         algor=rlxtype
      else if (hamflag.eq.'f') then
         vecfile=davname_f
         nstates=davstates_f
         nblock=dmain_f
         krydim=kdim_f
         step=stepsize_f
         eps=davtol_f
         niter=maxiter_f
         iortho=rlxortho_f
         toler=siltol_f
         algor=rlxtype_f
      endif

      if (algor.eq.2) iortho=3
      
!-----------------------------------------------------------------------
! Allocation of arays
!-----------------------------------------------------------------------
      ! On-diagonal elements of the Hamiltonian matrix
      allocate(hii(matdim))
      hii=0.0d0

      ! |Psi_n(t)>
      allocate(vec_old(matdim,nblock))
      vec_old=0.0d0

      ! |Psi_n(t+dt)>
      allocate(vec_new(matdim,nblock))
      vec_new=0.0d0

      ! H |Psi_n>
      allocate(hxvec(matdim,nblock))
      hxvec=0.0d0

      ! Energies
      allocate(ener(nblock))
      ener=0.0d0

      ! Residuals
      allocate(res(nblock))
      res=0.0d0

      ! Convergence flags
      allocate(lconv(nblock))
      lconv=.false.
      
      ! Current times
      allocate(currtime(nblock))
      currtime=0.0d0

      ! SIL Krylov subspace dimensions
      allocate(sildim(nblock))
      sildim=0.0d0
      
      return
      
    end subroutine initialise
      
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

      ! |Psi_n(t)> and |Psi_n(t+dt)>, n=1,nblock
      mem=mem+8.0d0*nblock*matdim/1024.0d0**2
      
      ! H |Psi>
      mem=mem+8.0d0*nblock*matdim/1024.0d0**2

      ! SIL: Lanczos vectors
      mem=mem+8.0d0*krydim*matdim/1024.0d0**2

      ! SIL: work arrays
      mem=mem+8.0d0*3*matdim/1024.0d0**2
      
      ! Determine whether we can hold the Hamiltonian matrix in-core
      if (mem.lt.maxmem) then
         lincore=.true.
         write(ilog,'(2x,a,/)') 'Matrix-vector multiplication will &
              proceed in-core'
      else
         lincore=.false.
         write(ilog,'(2x,a,/)') 'Matrix-vector multiplication will &
              proceed out-of-core'
      endif

      return

    end subroutine isincore

!#######################################################################

    subroutine initvec(matdim,noffd)

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd

      ! Determine whether the initial vectors are to be constructed from
      ! the ADC(1) eigenvectors, random noise, the eigenvectors of the
      ! Hamiltonian represented in a subspace of ISs, or as single ISs
      lrdadc1=.false.
      lrandom=.false.
      lsub=.false.
      if (hamflag.eq.'i'.and.ladc1guess) then
         lrdadc1=.true.
      else if (hamflag.eq.'f'.and.ladc1guess_f) then
         lrdadc1=.true.
      endif
      if (hamflag.eq.'i'.and.lnoise) then
         lrandom=.true.
      else if (hamflag.eq.'f'.and.lnoise_f) then
         lrandom=.true.
      endif
      if (hamflag.eq.'i'.and.lsubdiag) then
         lsub=.true.
      else if (hamflag.eq.'f'.and.lsubdiag_f) then
         lsub=.true.
      endif
      
      if (lrdadc1) then
         ! Construct the initial vectors from the ADC(1) eigenvectors
         call initvec_adc1
      else if (lrandom) then
         ! Construct the initial vectors as random orthonormal unit
         ! vectors
         call initvec_random(matdim)
      else if (lsub) then
         ! Construct the initial vectors from the eigenvectors of the
         ! Hamiltonian represented in a subspace of ISs
         call initvec_subdiag(matdim,noffd)
      else
         ! Use a single IS for each initial vector
         call initvec_ondiag(matdim)
      endif
      
      return
      
    end subroutine initvec
      
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
! Set the initial vectors
!-----------------------------------------------------------------------
      vec_old=0.0d0
      do i=1,nblock
         vec_old(1:dim1,i)=vec1(:,i)
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
      vec_old=0.0d0
      do i=1,nblock
         k=indx_hii(i)
         vec_old(k,i)=1.0d0
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

    subroutine initvec_random(matdim)

      implicit none

      integer, intent(in) :: matdim      
      integer             :: n,k,i,j
      real(d)             :: ftmp,dp

      ! Generate a set of random vectors
      do n=1,nblock
         do k=1,matdim
            call random_number(ftmp)
            vec_old(k,n)=ftmp
         enddo
         ftmp=dot_product(vec_old(:,n),vec_old(:,n))
         vec_old(:,n)=vec_old(:,n)/sqrt(ftmp)
      enddo

      ! Orthogonalisation
      do i=1,nblock
         do j=1,i-1
            dp=dot_product(vec_old(:,i),vec_old(:,j))
            vec_old(:,i)=vec_old(:,i)-dp*vec_old(:,j)
         enddo
      enddo
      do i=1,nblock
         do j=1,i-1
            dp=dot_product(vec_old(:,i),vec_old(:,j))
            vec_old(:,i)=vec_old(:,i)-dp*vec_old(:,j)
         enddo
      enddo

      ! Normalisation
      do i=1,nblock
         dp=dot_product(vec_old(:,i),vec_old(:,i))
         vec_old(:,i)=vec_old(:,i)/sqrt(dp)
      enddo
      
      return
      
    end subroutine initvec_random

!#######################################################################

    subroutine initvec_subdiag(matdim,noffd)

      use misc, only: dsortindxa1
      
      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer, dimension(:), allocatable   :: full2sub,sub2full,indxhii
      integer                              :: subdim,i,j,k,i1,j1,e2,&
                                              error,iham,nlim,l
      real(d), dimension(:,:), allocatable :: hsub
      real(d), dimension(:), allocatable   :: subeig,work
      character(len=70)                    :: filename
      
      ! Temporary hard-wiring of the subspace dimension
      subdim=1000
      
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
      if (abs(hii(indxhii(subdim))-hii(indxhii(subdim+1))).lt.1e-6) then
         subdim=subdim+1
         goto 5
      endif

!-----------------------------------------------------------------------
! Allocate arrays
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
      vec_old=0.0d0
      do i=1,nblock
         do j=1,subdim
            k=sub2full(j)
            vec_old(k,i)=hsub(j,i)            
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
    
    subroutine rdham_on(matdim)

      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer, intent(in) :: matdim
      integer             :: iham
      character(len=70)   :: filename

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

    subroutine orthowf(matdim)

      implicit none
      
      integer, intent(in) :: matdim
      integer             :: i,j
      real(d)             :: dp

      ! Double modified Gram-Schmidt orthogonalisation
      do i=1,nblock
         do j=1,i-1
            dp=dot_product(vec_new(:,i),vec_new(:,j))
            vec_new(:,i)=vec_new(:,i)-dp*vec_new(:,j)
         enddo
      enddo
      do i=1,nblock
         do j=1,i-1
            dp=dot_product(vec_new(:,i),vec_new(:,j))
            vec_new(:,i)=vec_new(:,i)-dp*vec_new(:,j)
         enddo
      enddo

      ! Normalisation
      do i=1,nblock
         dp=dot_product(vec_new(:,i),vec_new(:,i))
         vec_new(:,i)=vec_new(:,i)/sqrt(dp)
      enddo

      return

    end subroutine orthowf

!#######################################################################

    subroutine getener(matdim,noffd)

      implicit none

      integer, intent(in)   :: matdim
      integer*8, intent(in) :: noffd
      integer               :: k

      call hxpsi_all(matdim,noffd,vec_new,hxvec)

      do k=1,nblock
         ener(k)=dot_product(vec_new(:,k),hxvec(:,k))
      enddo
      
      return
      
    end subroutine getener
      
!#######################################################################

    subroutine hxpsi_all(matdim,noffd,vecin,vecout)

      use iomod, only: freeunit
      use constants
      use parameters
      
      implicit none

      integer, intent(in)               :: matdim
      integer*8, intent(in)             :: noffd
      integer                           :: iham,nlim,i,j,k,l,m,n
      real(d), dimension(matdim,nblock) :: vecin,vecout
      character(len=70)                 :: filename

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      vecout=0.0d0
      do m=1,nblock
         do n=1,matdim
            vecout(n,m)=hii(n)*vecin(n,m)
         enddo
      enddo

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      if (.not.lincore) then
         
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
            do l=1,nlim
               do m=1,nblock
                  vecout(indxi(l),m)=vecout(indxi(l),m)+hij(l)*vecin(indxj(l),m)
                  vecout(indxj(l),m)=vecout(indxj(l),m)+hij(l)*vecin(indxi(l),m)
               enddo
            enddo
         enddo
         
         close(iham)

         deallocate(hij,indxi,indxj)
      
      else        

         do m=1,nblock
            do l=1,noffd
               vecout(indxi(l),m)=vecout(indxi(l),m)+hij(l)*vecin(indxj(l),m)
               vecout(indxj(l),m)=vecout(indxj(l),m)+hij(l)*vecin(indxi(l),m)
            enddo
         enddo

      endif

      return
      
    end subroutine hxpsi_all

!#######################################################################
! silstep: propagates the wavefunction vector vec0 forward by a
!          single timestep using the short iterative Lanczos 
!          algorithm to yield the wavefunction vector vecprop
!#######################################################################

    subroutine silstep(ista,matdim,noffd,vec0,vecprop)

      use iomod

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: ista,i,j,n,info,truedim
      real(d), dimension(matdim)           :: vec0,vecprop
      real(d), dimension(:,:), allocatable :: qmat,eigvec
      real(d), dimension(:), allocatable   :: r,q,v,alpha,beta,&
                                              eigval,work,a,tmparr
      real(d)                              :: dp,norm,err,dtau

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(qmat(matdim,krydim))
      allocate(r(matdim))
      allocate(q(matdim))
      allocate(v(matdim))
      allocate(alpha(krydim))
      allocate(beta(krydim))
      allocate(eigvec(krydim,krydim))
      allocate(eigval(krydim))
      allocate(tmparr(krydim))
      allocate(a(krydim))
      
!-----------------------------------------------------------------------
! Generate the Lanczos vectors
!-----------------------------------------------------------------------
      qmat=0.0d0

      ! The first vector is either |Psi_n(t)> (unmodified algorithm) 
      ! or |Psi_n(t)> orthogonalised against the already propagated 
      ! wavefunctions and then renormalised (modified algorithm)
      if (iortho.eq.2) then
         do i=1,ista-1
            dp=dot_product(vec_new(:,i),vec0)
            vec0=vec0-dp*vec_new(:,i)
         enddo
         dp=dot_product(vec0,vec0)
         vec0=vec0/sqrt(dp)
      endif

      q=vec0
      qmat(:,1)=q

      ! alpha_1
      call hxkryvec(ista,1,matdim,noffd,q,r)
      alpha(1)=dot_product(q,r)

      ! beta_1
      r=r-alpha(1)*q
      beta(1)=sqrt(dot_product(r,r))

      truedim=1
      
      ! Remaining vectors
      do j=2,krydim
                 
         truedim=truedim+1

         ! Compute the next vector and pair of matrix elements
         v=q
         q=r/beta(j-1)
         qmat(:,j)=q
         call hxkryvec(ista,j,matdim,noffd,q,r)
         r=r-beta(j-1)*v
         alpha(j)=dot_product(q,r)
         r=r-alpha(j)*q
         beta(j)=sqrt(dot_product(r,r))

         ! Calculate the expansion of the propagated wavefunction in the
         ! current Lanczos state basis.
         ! If the error is below threshold, then exit here.
         ! Else, carry on and compute the next Lanczos vector.
         if (j.gt.2) then

            ! Diagonalise the Lanczos state representation of the
            ! Hamiltonian
            allocate(work((2*(truedim)-2)))
            eigval=alpha
            tmparr=beta
            call dstev('V',truedim,eigval(1:truedim),tmparr(1:truedim-1),&
                 eigvec(1:truedim,1:truedim),truedim,work,info)

            ! Calculate the epansion coefficients for |Psi(t0_dt)> in
            ! the Lanczos state basis
            a=0.0d0
            dtau=step
            do i=1,truedim
               do n=1,truedim
                  a(i)=a(i)+eigvec(i,n)*exp(-eigval(n)*dtau)*eigvec(1,n)
               enddo
            enddo            
            norm=sqrt(dot_product(a,a))
            a=a/norm

            ! Calcuate the error and exit if we are below
            ! threshold
            err=abs(a(truedim))
            deallocate(work)
            if (err.lt.toler) exit

         endif
         
      enddo

!-----------------------------------------------------------------------
! Calculate the wavefunction at time t0+dt, with dt being determined
! adaptively s.t. the estimated error in |Psi(t0+dt)> is below
! threshold.
!-----------------------------------------------------------------------
      dtau=step

10    continue
      a=0.0d0
      do j=1,truedim
         do n=1,truedim
            a(j)=a(j)+eigvec(j,n)*exp(-eigval(n)*dtau)*eigvec(1,n)
         enddo
      enddo
      
      ! OLD
!      vecprop=0.0d0
!      do j=1,truedim
!         vecprop=vecprop+a(j)*qmat(:,j)
!      enddo
!      norm=sqrt(dot_product(vecprop,vecprop))      
!      err=abs(a(truedim)/norm)
!      if (err.gt.toler) then
!         dtau=dtau/2.0d0
!         goto 10
!      endif
!      vecprop=vecprop/norm
      
      ! NEW
      ! For orthonormal Lanczos vectors, this holds true:
      norm=sqrt(dot_product(a,a))
      a=a/norm
      err=abs(a(truedim))
      if (err.gt.toler) then
         dtau=dtau/2.0d0
         goto 10
      endif
      vecprop=0.0d0
      do j=1,truedim
         vecprop=vecprop+a(j)*qmat(:,j)
      enddo
      
!-----------------------------------------------------------------------
! Update the propagation time for the current state
!-----------------------------------------------------------------------
      currtime(ista)=currtime(ista)+dtau

!-----------------------------------------------------------------------
! Save the Krylov subspace dimension that was used for the current
! state in this iteration
!-----------------------------------------------------------------------
      sildim(ista)=truedim
      
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(qmat)
      deallocate(r)
      deallocate(q)
      deallocate(v)
      deallocate(alpha)
      deallocate(beta)
      deallocate(eigvec)
      deallocate(eigval)
      deallocate(tmparr)
      deallocate(a)
      
      return
      
    end subroutine silstep

!#######################################################################
! modified_silstep: propagates the wavefunction vector vec0 forward by
!                   a single timestep using a modified short iterative
!                   Lanczos algorithm to yield the wavefunction vector
!                   vecprop
!
!                   The modification to the original SIL algorithm is
!                   that described in Appendix A of
!                   Z. Phys. D, 42, 113 (1997)
!#######################################################################
    
    subroutine modified_silstep(ista,matdim,noffd,vec0,vecprop)

      use iomod

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: ista,i,j,n,info,truedim
      real(d), dimension(matdim)           :: vec0,vecprop
      real(d), dimension(:,:), allocatable :: qmat,eigvec
      real(d), dimension(:), allocatable   :: r,q,v,alpha,beta,&
                                              eigval,work,a,tmparr
      real(d)                              :: dp,norm,err,dtau


      real(d), dimension(:), allocatable   :: eigval_mod,a_mod,diffvec
      real(d), dimension(:,:), allocatable :: eigvec_mod

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(qmat(matdim,krydim+1))
      allocate(r(matdim))
      allocate(q(matdim))
      allocate(v(matdim))
      allocate(alpha(krydim))
      allocate(beta(krydim))
      allocate(eigvec(krydim,krydim))
      allocate(eigval(krydim))
      allocate(tmparr(krydim))
      allocate(a(krydim))
      allocate(eigvec_mod(krydim+1,krydim+1))
      allocate(eigval_mod(krydim+1))
      allocate(a_mod(krydim+1))
      allocate(diffvec(krydim+1))
      
!-----------------------------------------------------------------------
! Generate the Lanczos vectors
!-----------------------------------------------------------------------
      qmat=0.0d0

      ! The first vector is either |Psi_n(t)> (unmodified algorithm) 
      ! or |Psi_n(t)> orthogonalised against the already propagated 
      ! wavefunctions and then renormalised (modified algorithm)
      if (iortho.eq.2) then
         do i=1,ista-1
            dp=dot_product(vec_new(:,i),vec0)
            vec0=vec0-dp*vec_new(:,i)
         enddo
         dp=dot_product(vec0,vec0)
         vec0=vec0/sqrt(dp)
      endif

      q=vec0
      qmat(:,1)=q

      ! alpha_1
      call hxkryvec(ista,1,matdim,noffd,q,r)
      alpha(1)=dot_product(q,r)

      ! beta_1
      r=r-alpha(1)*q
      beta(1)=sqrt(dot_product(r,r))

      truedim=1
      
      ! Remaining vectors
      do j=2,krydim
                 
         truedim=truedim+1

         ! Compute the next vector and pair of matrix elements
         v=q
         q=r/beta(j-1)
         qmat(:,j)=q
         call hxkryvec(ista,j,matdim,noffd,q,r)
         r=r-beta(j-1)*v
         alpha(j)=dot_product(q,r)
         r=r-alpha(j)*q
         beta(j)=sqrt(dot_product(r,r))

         ! Calculate the expansion of the propagated wavefunction in the
         ! current Lanczos state basis.
         ! If the error is below threshold, then exit here.
         ! Else, carry on and compute the next Lanczos vector.
         if (j.gt.2) then

            ! Diagonalise the UNMODIFIED Lanczos state representation
            ! of the Hamiltonian
            call sil_eigenpairs(truedim,alpha(1:truedim),beta(1:truedim-1),&
                 eigval(1:truedim),eigvec(1:truedim,1:truedim))
            
            ! Calculate the UNMODIFIED epansion coefficients for |Psi(t0_dt)>
            ! in the Lanczos state basis
            a=0.0d0
            call sil_expcoeff(truedim,a(1:truedim),eigval(1:truedim),&
                 eigvec(1:truedim,1:truedim),step)

            ! Diagonalise the MODIFIED Lanczos state representation
            ! of the Hamiltonian
            call modified_sil_eigenpairs(truedim,alpha(1:truedim),&
                 beta(1:truedim),eigval_mod(1:truedim+1),&
                 eigvec_mod(1:truedim+1,1:truedim+1))

            ! Calculate the MODIFIED epansion coefficients for |Psi(t0_dt)>
            ! in the Lanczos state basis
            a_mod=0.0d0
            call sil_expcoeff(truedim+1,a_mod(1:truedim+1),&
                 eigval_mod(1:truedim+1),&
                 eigvec_mod(1:truedim+1,1:truedim+1),step)
            
            ! Calcuate the error and exit if we are below
            ! threshold
            diffvec=a_mod
            diffvec(1:truedim)=diffvec(1:truedim)-a(1:truedim)
            err=sqrt(dot_product(diffvec,diffvec))
            if (err.lt.toler) exit

         endif

      enddo

!-----------------------------------------------------------------------      
! Calculate the last Lanczos state vector
!-----------------------------------------------------------------------
      q=r/beta(truedim)
      qmat(:,truedim+1)=q

!-----------------------------------------------------------------------
! Calculate the expansion coefficients, with dt being determined
! adaptively s.t. the estimated error in |Psi(t0+dt)> is below
! threshold.
!-----------------------------------------------------------------------
      dtau=step

10    continue

      ! Calculate the UNMODIFED expansion coefficients
      a=0.0d0
      call sil_expcoeff(truedim,a(1:truedim),eigval(1:truedim),&
           eigvec(1:truedim,1:truedim),dtau)

      ! Calculate the MODIFED expansion coefficients
      a_mod=0.0d0
      call sil_expcoeff(truedim+1,a_mod(1:truedim+1),&
           eigval_mod(1:truedim+1),&
           eigvec_mod(1:truedim+1,1:truedim+1),dtau)
      
      ! Calculate the error estimate
      diffvec=a_mod
      diffvec(1:truedim)=diffvec(1:truedim)-a(1:truedim)
      err=sqrt(dot_product(diffvec,diffvec))

      ! If the estimated error is above threshold, then reduce the
      ! timestep and recalculate the coefficients
      if (err.gt.toler) then
         dtau=dtau/1.25d0
         goto 10
      endif

!-----------------------------------------------------------------------
! Calculate the wavefunction at time t0+dt
!-----------------------------------------------------------------------      
      vecprop=0.0d0
      do j=1,truedim+1
         vecprop=vecprop+a_mod(j)*qmat(:,j)
      enddo
      
!-----------------------------------------------------------------------
! Update the propagation time for the current state
!-----------------------------------------------------------------------
      currtime(ista)=currtime(ista)+dtau

!-----------------------------------------------------------------------
! Save the Krylov subspace dimension that was used for the current
! state in this iteration
!-----------------------------------------------------------------------
      sildim(ista)=truedim+1

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(qmat)
      deallocate(r)
      deallocate(q)
      deallocate(v)
      deallocate(alpha)
      deallocate(beta)
      deallocate(eigvec)
      deallocate(eigval)
      deallocate(tmparr)
      deallocate(a)
      deallocate(eigvec_mod)
      deallocate(eigval_mod)
      deallocate(a_mod)
      deallocate(diffvec)
      
      return
      
    end subroutine modified_silstep

!#######################################################################

    subroutine sil_eigenpairs(n,alpha,beta,eigval,eigvec)

      use channels
      use constants
      use iomod
      
      implicit none

      integer                            :: n,info
      real(d), dimension(n)              :: alpha,eigval
      real(d), dimension(n-1)            :: beta,tmparr
      real(d), dimension(n,n)            :: eigvec
      real(d), dimension(:), allocatable :: work

!-----------------------------------------------------------------------
! Allocate work arrays
!-----------------------------------------------------------------------
      allocate(work(2*n-2))

!-----------------------------------------------------------------------
! Diagonalise the Lanzcos state representation of the Hamiltonian
!-----------------------------------------------------------------------
      eigval=alpha
      tmparr=beta
      call dstev('V',n,eigval(1:n),tmparr(1:n-1),eigvec(1:n,1:n),n,&
           work,info)

!-----------------------------------------------------------------------
! Check that the diagonalisation succeeded, and die here if not.
!-----------------------------------------------------------------------
      if (info.ne.0) then
         errmsg='Diagonalisation of the Lanczos state representation &
              of the Hamiltonian failed in subroutine sil_eigenpairs'
         call error_control
      endif
      
!-----------------------------------------------------------------------
! Deallocate work arrays
!-----------------------------------------------------------------------
      deallocate(work)
      
      return
      
    end subroutine sil_eigenpairs

!#######################################################################

    subroutine modified_sil_eigenpairs(n,alpha,beta,eigval,eigvec)

      use channels
      use constants
      use iomod
      
      implicit none

      integer                            :: n,info
      real(d), dimension(n)              :: alpha,beta
      real(d), dimension(n+1)            :: eigval
      real(d), dimension(n)              :: tmparr
      real(d), dimension(n+1,n+1)        :: eigvec
      real(d), dimension(:), allocatable :: work

!-----------------------------------------------------------------------
! Allocate work arrays
!-----------------------------------------------------------------------
      allocate(work(2*(n+1)-2))

!-----------------------------------------------------------------------
! Diagonalise the modified Lanzcos state representation of the
! Hamiltonian
!-----------------------------------------------------------------------
      eigval(1:n)=alpha
      eigval(n+1)=alpha(n)
      tmparr=beta
      call dstev('V',n+1,eigval(1:n+1),tmparr(1:n),eigvec(1:n+1,1:n+1),&
           n+1,work,info)

!-----------------------------------------------------------------------
! Check that the diagonalisation succeeded, and die here if not.
!-----------------------------------------------------------------------
      if (info.ne.0) then
         errmsg='Diagonalisation of the Lanczos state representation &
              of the Hamiltonian failed in subroutine &
              modified_sil_eigenpairs'
         call error_control
      endif
      
!-----------------------------------------------------------------------
! Deallocate work arrays
!-----------------------------------------------------------------------
      deallocate(work)

      return
      
    end subroutine modified_sil_eigenpairs
      
!#######################################################################

    subroutine sil_expcoeff(n,coeff,eigval,eigvec,dtau)
      
      use constants
      
      implicit none

      integer                 :: n,i,k
      real(d), dimension(n)   :: coeff,eigval
      real(d), dimension(n,n) :: eigvec
      real(d)                 :: dtau,norm

      coeff=0.0d0

      do i=1,n
         do k=1,n
            coeff(i)=coeff(i)+eigvec(i,k)*exp(-eigval(k)*dtau)*eigvec(1,k)
         enddo
      enddo

      norm=sqrt(dot_product(coeff,coeff))

      coeff=coeff/norm
      
      return
      
    end subroutine sil_expcoeff
    
!#######################################################################
    
    subroutine hxkryvec(ista,krynum,matdim,noffd,vecin,vecout)

      use iomod, only: freeunit
      use constants
      use parameters
      
      implicit none

      integer, intent(in)        :: ista,matdim
      integer*8, intent(in)      :: noffd
      integer                    :: krynum,iham,nlim,i,j,k,l,m,n
      real(d), dimension(matdim) :: vecin,vecout
      real(d)                    :: dp
      character(len=70)          :: filename

      vecout=0.0d0

!-----------------------------------------------------------------------
! Q |Psi>
!-----------------------------------------------------------------------
      if (iortho.eq.2) then
         ! SIL
         if (krynum.eq.2) then
            do i=1,ista-1
               dp=dot_product(vec_new(:,i),vecin)
               vecin=vecin-dp*vec_new(:,i)
            enddo
         endif
      else if (iortho.eq.3) then
         ! SIL-Liu
         do i=1,ista-1
            dp=dot_product(vec_conv(:,i),vecin)
            vecin=vecin-dp*vec_conv(:,i)
         enddo
      endif

!-----------------------------------------------------------------------
! H Q |Psi>
!
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      do n=1,matdim
         vecout(n)=hii(n)*vecin(n)
      enddo

!-----------------------------------------------------------------------
! H Q |Psi>
!
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------

      if (.not.lincore) then

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
            do l=1,nlim
               vecout(indxi(l))=vecout(indxi(l))+hij(l)*vecin(indxj(l))
               vecout(indxj(l))=vecout(indxj(l))+hij(l)*vecin(indxi(l))
            enddo
         enddo

         close(iham)
         
         deallocate(hij,indxi,indxj)

      else

         do l=1,noffd
            vecout(indxi(l))=vecout(indxi(l))+hij(l)*vecin(indxj(l))
            vecout(indxj(l))=vecout(indxj(l))+hij(l)*vecin(indxi(l))
         enddo

      endif

!-----------------------------------------------------------------------
! Q H Q |Psi>
!-----------------------------------------------------------------------
      if (iortho.eq.2) then
         ! SIL
         do i=1,ista-1
            dp=dot_product(vec_new(:,i),vecout)
            vecout=vecout-dp*vec_new(:,i)
         enddo
      else if (iortho.eq.3) then
         ! SIL-liu
         do i=1,ista-1
            dp=dot_product(vec_conv(:,i),vecin)
            vecin=vecin-dp*vec_conv(:,i)
         enddo
      endif

      return
      
    end subroutine hxkryvec

!#######################################################################

    subroutine get_kryvecs_current(ista,istep,matdim,noffd)

      implicit none

      integer, intent(in)                 :: matdim
      integer*8, intent(in)               :: noffd
      integer                             :: ista,istep,i,j,k1,k2
      real(d)                             :: dp
      real(d), dimension(:), allocatable  :: r,q
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(r(matdim))
      allocate(q(matdim))
      
!-----------------------------------------------------------------------
! Generate the Krylov vectors
!-----------------------------------------------------------------------
      ! Index of the first Krylov vector
      k1=(istep-1)*(krydim+1)+1

      ! The first vector is |Psi_n(t)> orthogonalised against the
      ! previously converged states and then renormalised
      kryvec(:,k1)=vec_old(:,ista)
      do i=1,ista-1
         dp=dot_product(vec_conv(:,i),kryvec(:,k1))
         kryvec(:,k1)=kryvec(:,k1)-dp*vec_conv(:,i)
      enddo

      dp=dot_product(kryvec(:,k1),kryvec(:,k1))
      kryvec(:,k1)=kryvec(:,k1)/sqrt(dp)
      knorm(k1)=sqrt(dp)
      
      q=kryvec(:,k1)      
      
      ! Remaining vectors
      do j=k1+1,istep*(krydim+1)

         call hxkryvec(ista,j,matdim,noffd,q,r)

         dp=dot_product(r,r)
         r=r/sqrt(dp)

         kryvec(:,j)=r

         knorm(j)=sqrt(dp)
         
         q=kryvec(:,j)
         
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(r)
      deallocate(q)

      return
      
    end subroutine get_kryvecs_current

!#######################################################################

    subroutine subspace_matrices(istep)

      implicit none

      integer :: istep,n,i,j,i1,j1

!-----------------------------------------------------------------------
! Fill in the sub-block of the subspace Hamiltonian and overlap
! matrices corresponding to the current timestep
!-----------------------------------------------------------------------
      i1=(istep-1)*krydim
      do i=(istep-1)*(krydim+1)+1,istep*(krydim+1)-1

         i1=i1+1

         j1=i1-1
         do j=i,istep*(krydim+1)-1

            j1=j1+1

            ! H_ij
            subhmat(i1,j1)=dot_product(kryvec(:,i),kryvec(:,j+1))
            subhmat(j1,i1)=subhmat(i1,j1)

            ! S_ij
            subsmat(i1,j1)=dot_product(kryvec(:,i),kryvec(:,j))
            subsmat(j1,i1)=subsmat(i1,j1)

         enddo
      enddo
      
!----------------------------------------------------------------------
! Couplings to the Krylov vectors of the previous timesteps
!----------------------------------------------------------------------
      ! Loop over previous timesteps
      do n=1,istep-1

         ! Loop over the Krylov vectors of the nth timestep
         i1=(n-1)*krydim
         do i=(n-1)*(krydim+1)+1,n*(krydim+1)-1

            i1=i1+1
            
            ! Loop over the Krylov vectors of the current timestep
            j1=(istep-1)*krydim
            do j=(istep-1)*(krydim+1)+1,istep*(krydim+1)-1

               j1=j1+1

               ! H_ij
               subhmat(i1,j1)=dot_product(kryvec(:,i),kryvec(:,j+1))
               subhmat(j1,i1)=subhmat(i1,j1)

               ! S_ij
               subsmat(i1,j1)=dot_product(kryvec(:,i),kryvec(:,j))
               subsmat(j1,i1)=subhmat(i1,j1)
               
            enddo
            
         enddo
            
      enddo

      return
      
    end subroutine subspace_matrices

!#######################################################################

    subroutine solve_subspace_tdse(istep)

      implicit none


      integer :: istep,k2
      
!----------------------------------------------------------------------
! Perform Lowdin's canonical orthogonalisation of the variational
! subspace basis vectors
!----------------------------------------------------------------------
      print*,
      print*,"WRITE THE CANONICAL ORTHOGONALISATION CODE!"
      STOP
      
      return

    end subroutine solve_subspace_tdse
    
!#######################################################################

    subroutine check_conv(matdim,noffd)

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: k
      real(d), dimension(:,:), allocatable :: hpsi
      real(d), dimension(:), allocatable   :: resvec
      real(d), dimension(nblock)           :: hpsi2,h2psi
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hpsi(matdim,nblock))
      allocate(resvec(matdim))

!-----------------------------------------------------------------------
! Residuals r_k=|| H |psi_k> - E_k |psi_k> ||
!-----------------------------------------------------------------------      
      call hxpsi_all(matdim,noffd,vec_new,hpsi)
      do k=1,nblock
         resvec=hpsi(:,k)-ener(k)*vec_new(:,k)
         res(k)=dot_product(resvec,resvec)
         res(k)=sqrt(res(k))
         if (res(k).lt.eps) lconv(k)=.true.
      enddo

      ! No. converged states (N.B. this runs over nstates rather than
      ! nblock as higher-lying states may converge before the lowest
      ! nstates states of interest)
      nconv=0
      do k=1,nstates
         if (res(k).lt.eps) nconv=nconv+1
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(hpsi)
      deallocate(resvec)

      return
      
    end subroutine check_conv

!#######################################################################

    subroutine sortwf(matdim)

      use misc, only: dsortindxa1
      
      implicit none

      integer, intent(in)                  :: matdim
      integer                              :: i
      integer, dimension(nblock)           :: indx,itmp
      real(d), dimension(nblock)           :: tmpval
      real(d), dimension(:,:), allocatable :: tmpvec

      allocate(tmpvec(matdim,nblock))

      call dsortindxa1('A',nblock,ener,indx)

      ! Wavefunctions
      do i=1,nblock
         tmpvec(:,i)=vec_new(:,indx(i))
      enddo
      vec_new=tmpvec
      
      ! Energies
      do i=1,nblock
         tmpval(i)=ener(indx(i))
      enddo
      ener=tmpval

      ! Residuals
      do i=1,nblock
         tmpval(i)=res(indx(i))
      enddo
      res=tmpval

      ! Current times
      do i=1,nblock
         tmpval(i)=currtime(indx(i))
      enddo
      currtime=tmpval

      ! Krylov subspace dimensions
      do i=1,nblock
         itmp(i)=sildim(indx(i))
      enddo
      sildim=itmp
      
      deallocate(tmpvec)
      
      return
      
    end subroutine sortwf
      
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
      if (k.eq.0) then
         write(ilog,'(64a)') ('*',j=1,64)
         write(ilog,'(a)') 'Iteration   Time   Kdim    Energies    &
              Residuals       Converged'         
         write(ilog,'(64a)') ('*',j=1,64)
      endif

!-----------------------------------------------------------------------
! Information from the current iteration
!-----------------------------------------------------------------------
      write(ilog,*)
      do i=1,nstates
         if (lconv(i)) then
            aconv='y'
         else
            aconv='n'
         endif         
         if (i.eq.1) then
            write(ilog,'(i4,3x,F10.3,2x,i3,3x,F12.7,3x,E13.7,3x,a1)') &
                 k,currtime(i),sildim(i),ener(i)*eh2ev,res(i),aconv
         else
            write(ilog,'(7x,F10.3,2x,i3,3x,F12.7,3x,E13.7,3x,a1)') &
                 currtime(i),sildim(i),ener(i)*eh2ev,res(i),aconv
         endif
      enddo

      return

    end subroutine wrtable

!#######################################################################

    subroutine wreigenpairs

      implicit none

      integer :: idav,i
      
!-----------------------------------------------------------------------
! Open the Davidson vector file
!-----------------------------------------------------------------------
      call freeunit(idav)
      open(unit=idav,file=vecfile,status='unknown',&
           access='sequential',form='unformatted')

!-----------------------------------------------------------------------
! Write the eigenpairs to file
!-----------------------------------------------------------------------
      do i=1,nstates
         write(idav) i,ener(i),vec_new(:,i)
      enddo

!-----------------------------------------------------------------------
! Close the Davidson vector file
!-----------------------------------------------------------------------
      close(idav)

      return
      
    end subroutine wreigenpairs

!#######################################################################

    subroutine finalise

      implicit none

      deallocate(hii)
      deallocate(vec_old)
      deallocate(vec_new)
      deallocate(hxvec)
      deallocate(ener)
      deallocate(lconv)
      deallocate(sildim)
      if (allocated(hij)) deallocate(hij)
      if (allocated(indxi)) deallocate(indxi)
      if (allocated(indxj)) deallocate(indxj)
      
      return
      
    end subroutine finalise
      
!#######################################################################
    
  end module relaxmod
