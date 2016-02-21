!#######################################################################
! Experimental block-relaxation algorithm for the calculation of the
! N lowest eigenpairs of the Hamiltonian matrix.
!
! A collection of N trial vectors are taken as a set of initial
! time-dependent electronic wavefunctions and propagated in negative
! imaginary time using the short iterative Lanczos method, followed by
! an orthogonalisation step.
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

    save

    integer                              :: maxbl,nrec,nblock,nconv,&
                                            nstates
    integer, dimension(:), allocatable   :: indxi,indxj
    real(d), dimension(:), allocatable   :: hii,hij,ener,sigma
    real(d), dimension(:,:), allocatable :: vec_old,vec_new,hxvec
    real(d)                              :: dt,eps,currtime
    character(len=36)                    :: vecfile
    logical, dimension(:), allocatable   :: lconv
    
  contains

!#######################################################################

    subroutine relaxation(matdim,noffd)

      implicit none

      integer, intent(in)      :: matdim
      integer*8, intent(in)    :: noffd
      integer                  :: n,k
      character(len=120)       :: atmp
      
!-----------------------------------------------------------------------
! Set the block-relaxation parameters
!-----------------------------------------------------------------------
      if (hamflag.eq.'i') then
         vecfile=davname
         nstates=davstates
         nblock=(5/4)*davstates
      else if (hamflag.eq.'f') then
         vecfile=davname_f
         nstates=davstates_f
         nblock=(5/4)*davstates_f
      endif

      ! Timestep
      dt=100.0d0

      ! Convergence tolerance
      eps=1d-6

      ! Current time
      currtime=0.0d0      
      
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
      allocate(hii(matdim))
      hii=0.0d0
      
      allocate(vec_old(matdim,nblock))
      vec_old=0.0d0

      allocate(vec_new(matdim,nblock))
      vec_new=0.0d0

      allocate(hxvec(matdim,nblock))
      hxvec=0.0d0

      allocate(ener(nblock))
      ener=0.0d0
      
      allocate(sigma(nblock))
      sigma=0.0d0

      allocate(lconv(nblock))
      lconv=.false.

!-----------------------------------------------------------------------
! TEMPORARY: force the loading of a guess ADC(1) vector from file
!-----------------------------------------------------------------------
      call getadc1vec(matdim)

!-----------------------------------------------------------------------
! Read the on-diagonal elements of the Hamiltonian matrix from disk
!-----------------------------------------------------------------------
      call rdham_on(matdim)
      
!-----------------------------------------------------------------------
! Output the initial energies and convergence information
!-----------------------------------------------------------------------
      vec_new=vec_old
      call hxc(matdim,noffd,vec_new,hxvec)
      do k=1,nblock
         ener(k)=dot_product(vec_new(:,k),hxvec(:,k))
      enddo
      call check_conv(matdim,noffd)
      call wrtable(0)

!-----------------------------------------------------------------------
! Perform the block-relaxation calculation
!-----------------------------------------------------------------------
      ! Loop over timesteps
      do n=1,30

         ! Update the time
         currtime=currtime+n*dt

         ! Propagate the wavefunctions to the next timestep using
         ! the short iterative Lanczos algorithm    
         do k=1,nblock
            if (.not.lconv(k)) then
                 call silstep(k,matdim,noffd,vec_old(:,k),vec_new(:,k))
              else
                 vec_new(:,k)=vec_old(:,k)
              endif
           enddo

         ! Orthogonalise the propagated wavefunctions amongst
         ! themselves
         call orthonorm

         ! Calculate the energies
         call hxc(matdim,noffd,vec_new,hxvec)
         do k=1,nblock
            ener(k)=dot_product(vec_new(:,k),hxvec(:,k))
         enddo

         ! Sort the wavefunctions by energy
         call sortwf(matdim)
         
         ! Check convergence
         call check_conv(matdim,noffd)

         ! Output the energies and convergence information
         call wrtable(n)
         
         ! Exit if all wavefunctions have converged
         if (nconv.eq.nblock) then
            write(ilog,'(/,2x,a,/)') 'All states have converged'
            exit
         endif
            
         ! Set the initial wavefunctions for the next step
         vec_old=vec_new
         
      enddo

!-----------------------------------------------------------------------
! Exit here if not all states have converged
!-----------------------------------------------------------------------
      if (nconv.lt.nblock) then
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
      deallocate(hii)
      deallocate(vec_old)
      deallocate(vec_new)
      deallocate(hxvec)
      deallocate(ener)
      deallocate(sigma)
      deallocate(lconv)

      return
      
    end subroutine relaxation

!#######################################################################

    subroutine getadc1vec(matdim)
      
      use iomod, only: freeunit
      use constants
      use parameters

      implicit none

      integer, intent(in)   :: matdim
      
      integer                              :: iadc1,dim1,i
      integer, dimension(:), allocatable   :: indx1
      real(d), dimension(:,:), allocatable :: vec1
      
      vec_old=0.0d0
      
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

      do i=1,nblock
         vec_old(1:dim1,i)=vec1(:,i)
      enddo
      
!-----------------------------------------------------------------------
! Close the ADC(1) eigenvector file
!-----------------------------------------------------------------------
      close(iadc1)
      
      return
      
    end subroutine getadc1vec

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

    subroutine hxc(matdim,noffd,vecin,vecout)

      use iomod, only: freeunit
      use constants
      use parameters
      
      implicit none

      integer, intent(in)               :: matdim
      integer*8, intent(in)             :: noffd
      integer                           :: iham,nlim,i,j,k,l,m,n
      real(d), dimension(matdim,nblock) :: vecin,vecout
      character(len=70)                 :: filename

      vecout=0.0d0

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      do m=1,nblock
         do n=1,matdim
            vecout(n,m)=hii(n)*vecin(n,m)
         enddo
      enddo

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
         do l=1,nlim
            do m=1,nblock
               vecout(indxi(l),m)=vecout(indxi(l),m)+hij(l)*vecin(indxj(l),m)
               vecout(indxj(l),m)=vecout(indxj(l),m)+hij(l)*vecin(indxi(l),m)
            enddo
         enddo
      enddo

      close(iham)

      deallocate(hij,indxi,indxj)
      
      return
      
    end subroutine hxc

!#######################################################################
! silstep: propagates the wavefunction vector vec0 forward by a
!          single timestep using the short iterative Lanczos method to
!          yield the wavefunction vector vecprop
!#######################################################################

    subroutine silstep(ista,matdim,noffd,vec0,vecprop)

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: ista,krydim,i,j,n,info
      real(d), dimension(matdim)           :: vec0,vecprop
      real(d), dimension(:,:), allocatable :: kryvec,lanvec,qmat,eigvec
      real(d), dimension(:), allocatable   :: r,q,v,alpha,beta,work,a
      real(d)                              :: ftmp
      
!-----------------------------------------------------------------------
! Maximum dimension of the Krylov space
!-----------------------------------------------------------------------
      krydim=100
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(kryvec(matdim,krydim))
      allocate(lanvec(matdim,krydim))
      allocate(qmat(matdim,krydim))
      allocate(r(matdim))
      allocate(q(matdim))
      allocate(v(matdim))
      allocate(alpha(krydim))
      allocate(beta(krydim))
      allocate(eigvec(krydim,krydim))
      allocate(work(2*krydim-2))
      allocate(a(krydim))
      
!-----------------------------------------------------------------------
! Generate the Lanczos vectors
!-----------------------------------------------------------------------
      qmat=0.0d0

      ! First vector is Psi(t0)
      q=vec0
      qmat(:,1)=q

      lanvec(:,1)=q
      
      ! alpha_1
      call hxkryvec(ista,matdim,noffd,q,r)
      alpha(1)=dot_product(q,r)

      ! beta_1
      r=r-alpha(1)*q
      beta(1)=sqrt(dot_product(r,r))

      ! Remaining vectors
      do j=2,krydim
         v=q
         q=r/beta(j-1)
         qmat(:,j)=q

         call hxkryvec(ista,matdim,noffd,q,r)
         r=r-beta(j-1)*v

         alpha(j)=dot_product(q,r)

         r=r-alpha(j)*q
        
         beta(j)=sqrt(dot_product(r,r))

         ! NOTE THAT IF beta_j IS ~0 THEN WE NEED TO TERMINATE HERE
         
      enddo

!-----------------------------------------------------------------------      
! Diagonalise the Lanczos state representation of the Hamiltonian
!-----------------------------------------------------------------------
      call dstev('V',krydim,alpha,beta(1:krydim-1),eigvec,krydim,&
           work,info)

      if (info.ne.0) then
         write(ilog,'(/,2x,a,/)') 'Diagonalisation of the Lanczos &
              state representation of the Hamiltonian failed in &
              subroutine silstep'
         STOP
      endif

!-----------------------------------------------------------------------
! Calculate the wavefunction at time t0+dt
!
! Note that following the call to dstev, the alpha array contains the
! eigenvalues of the Lanczos state representation of the Hamiltonian
!-----------------------------------------------------------------------
      a=0.0d0
      do j=1,krydim
         do n=1,krydim
            a(j)=a(j)+eigvec(j,n)*exp(-alpha(n)*dt)*eigvec(1,n)
         enddo
      enddo

      vecprop=0.0d0
      do j=1,krydim
         vecprop=vecprop+a(j)*qmat(:,j)
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(kryvec)
      deallocate(lanvec)
      deallocate(qmat)
      deallocate(r)
      deallocate(q)
      deallocate(v)
      deallocate(alpha)
      deallocate(beta)
      deallocate(eigvec)
      deallocate(work)
      deallocate(a)

      return
      
    end subroutine silstep

!#######################################################################
    
    subroutine hxkryvec(ista,matdim,noffd,vecin,vecout)

      use iomod, only: freeunit
      use constants
      use parameters
      
      implicit none

      integer, intent(in)        :: ista,matdim
      integer*8, intent(in)      :: noffd
      real(d), dimension(matdim) :: vecin,vecout

      integer                    :: iham
      integer                    :: nlim,i,j,k,l,m,n
      real(d)                    :: dp
      character(len=70)          :: filename

      vecout=0.0d0

!-----------------------------------------------------------------------
! Q |Psi>
!-----------------------------------------------------------------------
      do i=1,ista-1
         dp=dot_product(vec_new(:,i),vecin)
         vecin=vecin-dp*vec_new(:,i)
      enddo
         
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

!-----------------------------------------------------------------------
! Q H Q |Psi>
!-----------------------------------------------------------------------
      do i=1,ista-1
         dp=dot_product(vec_new(:,i),vecout)
         vecout=vecout-dp*vec_new(:,i)
      enddo
      
      return
      
    end subroutine hxkryvec
      
!#######################################################################

    subroutine orthonorm

      implicit none

      integer :: m,n
      real(d) :: ftmp

      do m=1,nblock
         do n=1,m-1
            ftmp=dot_product(vec_new(:,m),vec_new(:,n))
            vec_new(:,m)=vec_new(:,m)-&
                 ftmp*vec_new(:,n)/dot_product(vec_new(:,n),vec_new(:,n))
         enddo
      enddo
      do m=1,nblock
         do n=1,m-1
            ftmp=dot_product(vec_new(:,m),vec_new(:,n))
            vec_new(:,m)=vec_new(:,m)-&
                 ftmp*vec_new(:,n)/dot_product(vec_new(:,n),vec_new(:,n))
         enddo
      enddo

      do n=1,nblock
         ftmp=sqrt(dot_product(vec_new(:,n),vec_new(:,n)))
         vec_new(:,n)=vec_new(:,n)/ftmp
      enddo
      
      return
      
    end subroutine orthonorm

!#######################################################################

    subroutine check_conv(matdim,noffd)

      implicit none

      integer, intent(in)                  :: matdim
      integer*8, intent(in)                :: noffd
      integer                              :: k
      real(d), dimension(:,:), allocatable :: tmpvec,tmpvec2
      real(d), dimension(nblock)           :: hpsi2,h2psi
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(tmpvec(matdim,nblock))
      allocate(tmpvec2(matdim,nblock))
      
!-----------------------------------------------------------------------
! Calculate the standard deviations
! sigma_k=sqrt(<psi_k|H^2|psi_k> - <psi_k|H|psi_k>^2)
!-----------------------------------------------------------------------      
      ! <psi_k|H|psi_k>^2
      call hxc(matdim,noffd,vec_new,tmpvec)      
      do k=1,nstates
         hpsi2(k)=dot_product(vec_new(:,k),tmpvec(:,k))
         hpsi2(k)=hpsi2(k)**2
      enddo

      ! <psi_k|H^2|psi_k>      
      call hxc(matdim,noffd,tmpvec,tmpvec2)
      do k=1,nstates
         h2psi(k)=dot_product(vec_new(:,k),tmpvec2(:,k))
      enddo

      ! sigma_k
      nconv=0
      do k=1,nblock
         sigma(k)=sqrt(abs(h2psi(k)-hpsi2(k)))
         if (sigma(k).lt.eps) then
            nconv=nconv+1
            lconv(k)=.true.
         endif
      enddo
      
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(tmpvec)
      deallocate(tmpvec2)
      
      return
      
    end subroutine check_conv

!#######################################################################

    subroutine sortwf(matdim)

      use misc, only: dsortindxa1
      
      implicit none

      integer, intent(in)                  :: matdim
      integer                              :: i
      integer, dimension(nblock)           :: indx
      real(d), dimension(nblock)           :: tmpval
      real(d), dimension(:,:), allocatable :: tmpvec

      allocate(tmpvec(matdim,nblock))
      
      call dsortindxa1('A',nblock,ener,indx)

      do i=1,nblock
         tmpval(i)=ener(indx(i))
         tmpvec(:,i)=vec_new(:,indx(i))
      enddo

      ener=tmpval
      vec_new=tmpvec
      
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
         write(ilog,'(61a)') ('*',j=1,61)
         write(ilog,'(a)') &
              'Iteration   Time       Energies    Sigma           Converged'
         write(ilog,'(61a)') ('*',j=1,61)
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
            write(ilog,'(i4,3x,F10.3,3x,F12.7,3x,E13.7,3x,a1)') &
                 k,currtime,ener(i)*eh2ev,sigma(i),aconv
         else
            write(ilog,'(20x,F12.7,3x,E13.7,3x,a1)') &
                 ener(i)*eh2ev,sigma(i),aconv
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
      do i=1,nblock
         write(idav) i,ener(i),vec_new(:,i)
      enddo

!-----------------------------------------------------------------------
! Close the Davidson vector file
!-----------------------------------------------------------------------
      close(idav)

      return
      
    end subroutine wreigenpairs
      
!#######################################################################
    
  end module relaxmod
