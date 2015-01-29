  module lancdiag

    use constants
    use parameters
    use get_matrix
    use fspace

    implicit none

    integer :: main,nsat,mdim,noffdel
    real(d), dimension(:), allocatable :: diag, offdiag
    integer, dimension(:), allocatable :: indi, indj

  contains

!#######################################################################
!
!    subroutine master_lancdiag(ndim,noff,flag)
!    
!      integer, intent(in)     :: ndim,noff
!      character(2),intent(in) :: flag
!      
!      integer                 :: i,j
!
!!-----------------------------------------------------------------------
!! ndim: Hamiltonian matrix dimension
!! noff: no. non-zero off-diagonal Hamiltonian matrix elements
!! flag: ???
!! lmain: block-Lanczos block size
!!-----------------------------------------------------------------------      
!
!!-----------------------------------------------------------------------      
!! Enter block-Lanczos routine
!!-----------------------------------------------------------------------            
!      call lanczos_diag(ndim,lmain,lancstates,ncycles)
!
!      return
!
!    end subroutine master_lancdiag

!#######################################################################
    
!    subroutine lanczos_diag(matdim,blckdim,lancstates,ncycles)
!
!      implicit none
!
!      integer                                       :: matdim,blckdim,&
!                                                       lancstates,&
!                                                       ncycles,nvopu,&
!                                                       nneig
!
!      integer                                       :: lflag
!      integer, dimension(:), allocatable            :: istor
!      double precision, dimension(:), allocatable   :: rstor      
!      double precision, dimension(:,:), allocatable :: eig,X,U,V
!      double precision                              :: eigl,eigr,sigma
!
!!-----------------------------------------------------------------------
!! Initialise blzpack arguments
!!-----------------------------------------------------------------------
!      eigl=erange(1)
!      eigr=erange(2)
!      call init_blzpack(istor,rstor,matdim,blckdim,eigl,eigr,eig,X,U,V,&
!           lancstates,ncycles)
!
!!-----------------------------------------------------------------------
!! Block-Lanczos diagonalisation of the ADC Hamiltonian matrix
!!-----------------------------------------------------------------------
!      lflag=0
!10    continue
!      call blzdrd(istor,rstor,sigma,nneig,U,V,lflag,nvopu,eig,X)
!      
!      if (lflag.eq.1) then
!         ! Matrix-vector multiplication is required
!         call matvecmul_lanc(U,V,matdim,blckdim,nvopu)
!         goto 10
!      else if (lflag.eq.4) then
!         ! Specification of the starting vectors is required
!         call initvecs_lanc(V,matdim,blckdim,nvopu)
!         goto 10
!      else if (lflag.ne.0) then
!         ! Unsuccessful termination
!         call lanc_err(lflag)
!      endif
!
!!-----------------------------------------------------------------------
!! Write Lanczos pseudospectrum to disk
!!-----------------------------------------------------------------------
!      call wrvecs_lanc(matdim,lancstates,eig,X,istor(4))
!
!      return
!
!    end subroutine lanczos_diag

!#######################################################################

    subroutine init_blzpack(istor,rstor,matdim,blckdim,eigl,eigr,eig,&
         X,U,V,lancstates,ncycles)

      implicit none

      integer                                       :: matdim,blckdim,&
                                                       lancstates,&
                                                       ncycles,ierr,&
                                                       itmp,k1,k2,k3,k4

      integer, dimension(:), allocatable            :: istor
      double precision, dimension(:), allocatable   :: rstor      
      double precision, dimension(:,:), allocatable :: eig,X,U,V
      double precision                              :: eigl,eigr

!-----------------------------------------------------------------------
! Allocate the istor array
!-----------------------------------------------------------------------
      k1=min(matdim,180)

      ! THIS NEEDS ADDRESSING
!      itmp=17+(123+12*k1)
      itmp=17+(123+12*k1)*10

      allocate(istor(itmp),stat=ierr)
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate istor in subroutine init_blzpack'
         STOP
      endif

!-----------------------------------------------------------------------
! Fill in the istor array 
!-----------------------------------------------------------------------
! NI: no. active rows on each processor, equal to the matrix dimension
!     for serial execution
!
! LNI: leading dimension of the matrix to be diagonalised
!
! NREIG: no. required eigenvalues
!
! LEIG: dimension of the eig array: needs to be large enough to
!       store all converged eigenvalues, which may be more than the
!       no. requested 
!
! NVBSET: no. vectors in a block
!
! NSTEPS: max no. steps
!
! NSTART: no. starting vectors, if >0 the user must specify the
!         starting vectors
!
! GNRZD: specifies the problem type: 0 for a standard eigenvalue
!        problem
!
! LPRNT: level of information to be printed: 1 for header errors and
!        warnings, 0 for no printing
!
! LFILE: file unit where information (if any) is printed, e.g., 6 for
!        printing to the screen
!
! LISTOR: dimension of the istor array is 17+LISTOR
!-----------------------------------------------------------------------
      istor=0

      istor(1)=matdim ! NI

      istor(2)=matdim ! LNI

      istor(3)=lancstates ! NREIG

      istor(4)=min(istor(3)+10,istor(3)*2) ! LEIG

      istor(5)=blckdim ! NVBSET
      
      istor(6)=ncycles ! NSTEPS

      istor(7)=blckdim ! NSTART

      istor(9)=0 ! GNRZD

 !     istor(12)=1 ! LPRNT

 !     istor(13)=6 ! LFILE

      ! THIS NEEDS ADDRESSING
      ! istor(15)=123+k1*12 ! LISTOR
      istor(15)=(123+k1*12)*10 ! LISTOR

!-----------------------------------------------------------------------
! Allocate the rstor array
!-----------------------------------------------------------------------
      if (blckdim.eq.0) then
         k2=3
      else
         k2=blckdim
      endif
      k4=min(istor(4),matdim)
      k3=484+k1*(13+k1*2+k2+max(18,k2+2))+k2*k2*3+k4*2

      ! THIS NEEDS ADDRESSING
!       itmp=5+(dim*k2*5+k3)
       itmp=5+(matdim*k2*5+k3)*20

      allocate(rstor(itmp))
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate rstor in subroutine init_blzpack'
         STOP
      endif

!-----------------------------------------------------------------------
! Fill in the rstor array
!-----------------------------------------------------------------------
! EIGL: lower limit for the eigenvalue interval
!
! EIGR: lower limit for the eigenvalue interval
!
! LRSTOR: dimension of the rstor array is 4+LRSTOR
!-----------------------------------------------------------------------
      rstor=0.0d0

      rstor(1)=eigl ! EIGL (lower limit for the eigenvalue interval)

      rstor(2)=eigr ! EIGR (lower limit for the eigenvalue interval)

      ! THIS NEEDS ADDRESSING
      ! rstor(4)=dim*k2*5+k3 ! LRSTOR
      rstor(4)=(matdim*k2*5+k3)*20 ! LRSTOR (dimension of the rstor array
                                ! is 4+LRSTOR)
!-----------------------------------------------------------------------
! Allocate the eig, X, U and V arrays
!-----------------------------------------------------------------------
      ! eig (eigenvalues on exit)
      allocate(eig(istor(4),2))
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate eig in subroutine init_blzpack'
         STOP
      endif

      ! X (eigenvectors on exit)
      allocate(X(matdim,istor(4)))
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate X in subroutine init_blzpack'
         STOP
      endif

      ! U (vectors recieved from blzdrd for use in matrix-vector
      ! multiplication) 
      allocate(U(matdim,blckdim))
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate U in subroutine init_blzpack'
         STOP
      endif

      ! V (matrix-vector products to be sent to blzdrd)
      allocate(V(matdim,blckdim))
      if (ierr.ne.0) then
         write(6,'(/,a,/)') &
              'Failed to allocate V in subroutine init_blzpack'
         STOP
      endif
      
      return

    end subroutine init_blzpack

!#######################################################################
    
    subroutine initvecs_lanc(V,n1,n2,nvecs)

      use parameters, only: stvc_lbl

      implicit none

      integer                            :: dim,nvecs,i,k,n1,n2
      double precision, dimension(n1,n2) :: V
      double precision                   :: ftmp

!-----------------------------------------------------------------------
! Generate guess vectors: blckdim unit vectors corresponding to the 1h1p
! intermediate states with the greatest TDMs with the initial state
!-----------------------------------------------------------------------
      V=0.0d0
      do i=1,nvecs
         k=stvc_lbl(i)
         V(k,i)=1.0d0
      enddo

      return

    end subroutine initvecs_lanc

!#######################################################################

    subroutine matvecmul_lanc(U,V,matdim,blckdim,nvopu)

      implicit none

      integer                                     :: matdim,blckdim,&
                                                     nvopu,unit,maxbl,&
                                                     nrec,nlim,i,k,l
      double precision, dimension(matdim,blckdim) :: U,V
      real(d), dimension(matdim)                  :: hii
      real(d), dimension(:), allocatable          :: hij
      integer, dimension(:), allocatable          :: indxi,indxj

!-----------------------------------------------------------------------
! Save the ith matrix-vector product H*U(:,i) to V(:,i)
!-----------------------------------------------------------------------
      V=0.0d0

!-----------------------------------------------------------------------
! On-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      unit=77
      open(unit,file='hmlt.diac',status='old',access='sequential',&
           form='unformatted')

      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

      do i=1,nvopu
         do k=1,matdim
            V(k,i)=hii(k)*U(k,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! Off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
      
      open(unit,file='hmlt.offc',status='old',access='sequential',&
           form='unformatted')
            
      do k=1,nrec
         read(unit) hij(:),indxi(:),indxj(:),nlim
         do l=1,nlim            
            do i=1,nvopu
               V(indxi(l),i)=V(indxi(l),i)+hij(l)*U(indxj(l),i)               
            enddo
         enddo
      enddo

      close(unit)

      deallocate(hij,indxi,indxj)

    end subroutine matvecmul_lanc

!#######################################################################

    subroutine lanc_err(lflag)

      implicit none

      integer :: lflag

      write(6,'(/,2x,a,1x,i2,/)') &
           'Unsuccessful block-Lanczos termination. BLZPACK error code:',lflag
      STOP

      return

    end subroutine lanc_err

!#######################################################################

    subroutine wrvecs_lanc(matdim,lancstates,eig,X,nvecs)

      implicit none

      integer                                   :: matdim,lancstates,i,&
                                                   unit,nvecs
      double precision, dimension(nvecs,2)      :: eig
      double precision, dimension(matdim,nvecs) :: X

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=77
      open(unit=unit,file=lancname,status='unknown',access='sequential',&
           form='unformatted')

!-----------------------------------------------------------------------
! Write Lanczos pseudospectrum to file
!-----------------------------------------------------------------------
      do i=1,lancstates
         write(unit) i,eig(i,1),X(:,i)
      enddo

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrvecs_lanc

!#######################################################################

  end module lancdiag
