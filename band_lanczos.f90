  module band_lanczos

    use constants
    use parameters
    use get_matrix
    use fspace

    implicit none

    integer                            :: main,nsat,mdim,noffdel
    real(d), dimension(:), allocatable :: diag, offdiag
    integer, dimension(:), allocatable :: indi, indj

!-----------------------------------------------------------------------
! matdim:   matrix dimension
! iblckdim: initial block size
! cblckdim: current block size
! dtol:     tolerance for deflation
! maxit:    max no. iterations
! indxdfl:  deflation index array
! ndfl:     no. deflations performed
! kryvec:   Krylov vectors
! lanvec:   Lanczos vectors
! mat:      matrix whose Lanczos pseudospectrum is sought
! lex:      logical flag, true if we have exhausted the Krylov sequence
! lanunit:  Unit number corresponding to the file to which the Lanczos
!           vectors are to be written
! nnzrtd:   no. non-zero rows in T^(d)
!-----------------------------------------------------------------------
    integer*8                            :: matdim,iblckdim,cblckdim,&
                                            maxit,ndfl,lanunit,nnzrtd
    integer*8, dimension(:), allocatable :: indxdfl
    real*8, parameter                    :: dfltol=1e-6
    real*8, dimension(:,:), allocatable  :: kryvec,lanvec
    logical(kind=4)                      :: lex
    
    integer*8                            :: maxdim
    real*8, dimension(:,:), allocatable  :: tmat,smat,prtmat

  contains

!#######################################################################
! ndim: Hamiltonian matrix dimension
! noff: no. non-zero off-diagonal Hamiltonian matrix elements
! flag: ???
!#######################################################################

    subroutine master_lancdiag(ndim,noff,flag)
    
      integer, intent(in)     :: ndim
      integer*8, intent(in)   :: noff
      character(1),intent(in) :: flag

!-----------------------------------------------------------------------
! Set matrix dimension, initial block size, max no. iterations, etc
!-----------------------------------------------------------------------
      matdim=ndim
      iblckdim=lmain
      maxit=ncycles
      lex=.false.
      nnzrtd=0

!-----------------------------------------------------------------------
! Allocate and initialse arrays
!-----------------------------------------------------------------------
      call alloc_blanc

!-----------------------------------------------------------------------
! Set initial vectors
!-----------------------------------------------------------------------
      call init_vec

!-----------------------------------------------------------------------
! Calculate the band Lanczos pseudospectrum
!-----------------------------------------------------------------------
      call run_band_lanczos

      return
      
    end subroutine master_lancdiag

!#######################################################################
! alloc_blanc: allocates and initialises arrays used by the 
! band-Lanczos routine
!#######################################################################

    subroutine alloc_blanc
      
      implicit none
      
      ! Maximum number of Lanczos vectors
      maxdim=iblckdim*maxit
      
      ! Deflation index array
      allocate(indxdfl(iblckdim))
      indxdfl=0
      
      ! Krylov vectors
      allocate(kryvec(matdim,iblckdim+1))
      kryvec=0.0d0
      
      ! Lanczos vectors
      allocate(lanvec(matdim,iblckdim+1))
      lanvec=0.0d0

      ! Temporary: T-matrix, S-matrix and T^(pr)-matrix
      allocate(tmat(maxdim,maxdim))
      allocate(smat(maxdim,maxdim))
      allocate(prtmat(maxdim,maxdim))
      
      tmat=0.0d0
      smat=0.0d0
      prtmat=0.0d0
      
      return
      
    end subroutine alloc_blanc
    
!#######################################################################
! init_vec: initialises the first block of Krylov vectors
!#######################################################################

    subroutine init_vec

      implicit none

      integer*8 :: i,k

      do i=1,iblckdim
         k=stvc_lbl(i)
         kryvec(k,i)=1.0d0
      enddo
      
      return

    end subroutine init_vec

!#######################################################################
! run_band_lanczos: main band-Lanczos routine
!#######################################################################
    
    subroutine run_band_lanczos

      use parameters, only: lancstates

      implicit none

      integer*8 :: i,j,n,k,kt,k0,indx,kcount,lcount
      real*8    :: norm,frac
      
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(12x,a)') &
           'Band-Lanczos Diagonalisation in the Final Space'
      write(6,'(70a,/)') ('-',i=1,70)
      write(6,'(2x,a,x,i4,2x,a,x,i5,x,a,/)') 'A maximum of',&
           maxit,'Lanczos blocks will be generated (<=',&
           maxit*iblckdim,'vectors)'

!-----------------------------------------------------------------------
! Initialise current block size and the number of deflations, and set
! the unit number for the Lanczos vector file
!-----------------------------------------------------------------------
      cblckdim=iblckdim
      ndfl=0
      lanunit=521

!-----------------------------------------------------------------------
! Perform the band Lanczos iterations until either we reach the max. 
! number of Lanczos vectors or the Krylov sequence is exhausted
!-----------------------------------------------------------------------
      j=0
10    continue
      j=j+1

!-----------------------------------------------------------------------
! Output progress
!-----------------------------------------------------------------------
      frac=real(j)/real(iblckdim)
      if (j.ne.1.and.floor(frac)-frac.eq.0) then
         if (int(frac).eq.1) then
            write(6,'(70a)') ('-',i=1,70)
            write(6,'(2x,i4,2x,a)') int(frac),'Lanczos block  generated'
         else
            write(6,'(70a)') ('-',i=1,70)
            write(6,'(2x,i4,2x,a)') int(frac),'Lanczos blocks generated'
         endif
      endif

!-----------------------------------------------------------------------
! Reorder the Krylov and Lanczos vector arrays such that the jth vectors
! are now the first vectors, etc
!-----------------------------------------------------------------------
      if (j.gt.1) then
         do i=1,cblckdim
            kryvec(:,i)=kryvec(:,i+1)
            lanvec(:,i)=lanvec(:,i+1)
         enddo
      endif

!-----------------------------------------------------------------------
! Calculate the norm of the jth Krylov vector
!-----------------------------------------------------------------------
15    continue
      norm=norm_euclid(kryvec(:,1),matdim)

!-----------------------------------------------------------------------
! If the norm is less than the tolerance for deflation, then deflate
! the vector
!-----------------------------------------------------------------------
      if (norm.lt.dfltol) then
         call deflate_kryvec(j,cblckdim,indxdfl,kryvec,matdim,&
              iblckdim,nnzrtd)
         goto 15
      endif
      
!-----------------------------------------------------------------------
! If the jth Krylov vector has passed the deflation test, then
! we now construct the next Lanczos vector
!-----------------------------------------------------------------------
      if (j.le.cblckdim) then
         lanvec(:,cblckdim+1)=kryvec(:,1)
      else
         tmat(j,j-cblckdim)=norm_euclid(kryvec(:,1),matdim)
         lanvec(:,cblckdim+1)=kryvec(:,1)/tmat(j,j-cblckdim)
         call normalise(lanvec(:,cblckdim+1),matdim)
      endif

!-----------------------------------------------------------------------
! Orthogonalise the remaining Krylov vectors k_j+1,...,k_j+pc-1 against
! the the new Lanczos vector l_j
!-----------------------------------------------------------------------
      kcount=1
      do k=j+1,j+cblckdim-1
         kcount=kcount+1
         if (k-cblckdim.lt.1) cycle
         tmat(j,k-cblckdim)=&
              dot_product(lanvec(:,cblckdim+1),kryvec(:,kcount))
         kryvec(:,kcount)=kryvec(:,kcount)&
              -tmat(j,k-cblckdim)*lanvec(:,cblckdim+1)
      enddo

!-----------------------------------------------------------------------
! Calculate the j+pc'th Krylov vector
!-----------------------------------------------------------------------
      kryvec(:,cblckdim+1)=hxlanvec()

!-----------------------------------------------------------------------
! Orthogonalise the j+pc'th Krylov vector against the previous Lanczos
! vectors l_k0,...,l_j-1, k=max{1,j-pc}
!-----------------------------------------------------------------------
      if (j.gt.1) then

         lcount=0
         k0=max(1,j-cblckdim)

         if (j.le.cblckdim) then
            do k=k0,j-1
               tmat(k,j)=tmat(j,k)
               lcount=lcount+1
               indx=cblckdim+1-(j-lcount)
               kryvec(:,cblckdim+1)=kryvec(:,cblckdim+1)-tmat(k,j)*lanvec(:,indx)
            enddo
         else           
            do k=k0,j-1
               tmat(k,j)=tmat(j,k)
               lcount=lcount+1
               indx=lcount
               kryvec(:,cblckdim+1)=kryvec(:,cblckdim+1)-tmat(k,j)*lanvec(:,indx)
            enddo
         endif

      endif

!-----------------------------------------------------------------------
! Orthoganalise the j+pc'th Krylov vector against the Lanczos vectors
! l_k, k in the set I U {j}
!-----------------------------------------------------------------------
      ! (i) k in I
      if (ndfl.gt.0) then
         do k=1,nnzrtd
            indx=indxdfl(k)
            call ortho_indxdfl(j,indx,lanunit,matdim,tmat(indx,j),&
                 kryvec(:,cblckdim+1))
         enddo
      endif

      ! (ii) k=j
      tmat(j,j)=dot_product(lanvec(:,cblckdim+1),kryvec(:,cblckdim+1))
      kryvec(:,cblckdim+1)=&
           kryvec(:,cblckdim+1)-tmat(j,j)*lanvec(:,cblckdim+1)

!-----------------------------------------------------------------------
! For k in the set I, set S_jk = T_kj
!-----------------------------------------------------------------------
      if (ndfl.gt.0) then
         do k=1,nnzrtd
            indx=indxdfl(k)
            smat(j,indx)=tmat(indx,j)
         enddo
      endif

!-----------------------------------------------------------------------
! Construct the jth-order projection of the Hamiltonian onto the
! Lanczos basis
!-----------------------------------------------------------------------
      prtmat=tmat+smat

!-----------------------------------------------------------------------
! Write the jth Lanczos vector to file
!-----------------------------------------------------------------------
      call wrlanvec(j,lanunit,lanvec(:,cblckdim+1),matdim)

!-----------------------------------------------------------------------
! Goto the next iteration if we haven't reached the max. no. vectors
!-----------------------------------------------------------------------
      if (j.lt.maxit*iblckdim) goto 10

!-----------------------------------------------------------------------
! If we are here, then either we have reached the max. no. iterations
! or the Krylov sequence has been exhausted.
!-----------------------------------------------------------------------
20    continue
      lancstates=j
      close(lanunit)

!-----------------------------------------------------------------------
! Deallocate Krylov and Lanczos vector arrays now that they are no
! longer needed
!-----------------------------------------------------------------------
      deallocate(kryvec)
      deallocate(lanvec)

!      call chkortho(j,matdim)
!      STOP

!-----------------------------------------------------------------------
! Calculate the Lanczos pseudospectrum from the jth-order projection
! of the Hamiltonian onto the Lanczos basis
!-----------------------------------------------------------------------
      call lanczos_pseudospec(cblckdim,prtmat(1:j,1:j),j,ndfl,&
           lanunit,matdim,iblckdim)
      
      return

    end subroutine run_band_lanczos

!#######################################################################

    subroutine deflate_kryvec(j,cblckdim,indxdfl,kryvec,matdim,&
         iblckdim,nnzrtd)
      
      implicit none
      
      integer*8                            :: j,i,cblckdim,iblckdim,&
                                              matdim,nnzrtd
      integer*8, dimension(iblckdim)       :: indxdfl
      real*8, dimension(matdim,iblckdim+1) :: kryvec

      write(6,'(/,7x,a,/)') 'Deflating the subspace'
      
!-----------------------------------------------------------------------
! Update the number of deflations
!-----------------------------------------------------------------------
      ndfl=ndfl+1

!-----------------------------------------------------------------------
! Update the deflation index array
!-----------------------------------------------------------------------
      if (j-cblckdim.gt.0) then
         nnzrtd=nnzrtd+1
         indxdfl(nnzrtd)=j-cblckdim
      endif

!-----------------------------------------------------------------------
! Update the current block size
!-----------------------------------------------------------------------
      cblckdim=cblckdim-1

!-----------------------------------------------------------------------
! Reorder the Krylov vector array
!-----------------------------------------------------------------------
      do i=1,cblckdim
         kryvec(:,i)=kryvec(:,i+1)
      enddo
      
      return
      
    end subroutine deflate_kryvec

!#######################################################################

    function norm_euclid(vec,dim) result(norm)

      implicit none
      
      integer*8              :: dim,i
      real*8, dimension(dim) :: vec
      real*8                 :: norm
      
      norm=0.0d0
      do i=1,dim
         norm=norm+vec(i)**2
      enddo
      
      norm=sqrt(norm)
      
      return
      
    end function norm_euclid

!#######################################################################

    subroutine normalise(vec,dim)
      
      implicit none
      
      integer*8              :: dim,i
      real*8, dimension(dim) :: vec
      real*8                 :: norm
      
      norm=0.0d0
      
      do i=1,dim
         norm=norm+vec(i)**2
      enddo
        
      norm=sqrt(norm)
      
      vec=vec/norm
      
      return
      
    end subroutine normalise

!#######################################################################

    function hxlanvec() result(kvec)

      implicit none
      
      integer                            :: unit
      integer                            :: maxbl,nrec,nlim,i,k,l
      integer, dimension(:), allocatable :: indxi,indxj
      real(d), dimension(matdim)         :: hii
      real(d), dimension(:), allocatable :: hij
      real*8, dimension(matdim)          :: kvec

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      unit=77
      open(unit,file='hmlt.diac',status='old',access='sequential',&
           form='unformatted')

      read(unit) maxbl,nrec
      read(unit) hii

      close(unit)

      kvec=0.0d0
      do k=1,matdim
         kvec(k)=hii(k)*lanvec(k,cblckdim+1)
      enddo

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
      
      open(unit,file='hmlt.offc',status='old',access='sequential',&
           form='unformatted')

      do k=1,nrec
         read(unit) hij(:),indxi(:),indxj(:),nlim
         do l=1,nlim           
            kvec(indxi(l))=&
                 kvec(indxi(l))+hij(l)*lanvec(indxj(l),cblckdim+1)            
            kvec(indxj(l))=&
                 kvec(indxj(l))+hij(l)*lanvec(indxi(l),cblckdim+1)
         enddo
      enddo

      close(unit)

      deallocate(hij,indxi,indxj)

      return

    end function hxlanvec

!#######################################################################

    subroutine ortho_indxdfl(j,indx,lanunit,matdim,tkj,kvec)

      implicit none

      integer*8                 :: j,indx,lanunit,matdim,n,i
      real*8                    :: tkj
      real*8, dimension(matdim) :: kvec,lvec,tmpvec

!-----------------------------------------------------------------------
! Rewind Lanczos vector file to the vector of interest
!-----------------------------------------------------------------------
        rewind(lanunit)
        do i=1,indx-1
           read(lanunit) tmpvec
        enddo

!-----------------------------------------------------------------------
! Read the vector of interest (the indx'th Lanczos vector)
!-----------------------------------------------------------------------
        read(lanunit) lvec

!-----------------------------------------------------------------------
! Go back to the original position in the Lanczos vector file
!-----------------------------------------------------------------------
        rewind(lanunit)
        do i=1,j-1
           read(lanunit) tmpvec
        enddo

!-----------------------------------------------------------------------
! Orthogonalise the latest Krylov vector against the indx'th Lanczos
! vector 
!-----------------------------------------------------------------------
        tkj=dot_product(lvec,kvec)
        kvec=kvec-tkj*lvec
      
      return
        
    end subroutine ortho_indxdfl

!#######################################################################

    subroutine wrlanvec(j,lanunit,lanvec,dim)

      implicit none

      integer*8              :: j,lanunit,dim
      real*8, dimension(dim) :: lanvec

!-----------------------------------------------------------------------
! If this is the first iteration, then open the Lanczos vector file
!-----------------------------------------------------------------------
      if (j.eq.1) then
         open(lanunit,file='lanvecs',form='unformatted',&
              status='unknown')
      endif

!-----------------------------------------------------------------------
! Write the jth Lanczos vector to file
!-----------------------------------------------------------------------
      write(lanunit) lanvec(:)

    end subroutine wrlanvec

!#######################################################################

    subroutine lanczos_pseudospec(blckdim,matrix,dim,ndfl,lanunit,&
         matdim,iblckdim)

      implicit none
        
      integer*8                  :: blckdim,dim,ndfl,lanunit,matdim,i,&
                                    iblckdim
      real*8, dimension(dim,dim) :: matrix,umat
      real*8, dimension(dim)     :: eigval

      real*8 :: t1,t2

!-----------------------------------------------------------------------
! Diagonalise the projection of the Hamiltonian onto the space spanned
! by the Lanczos vectors.
!
! N.B., if ndfl=0, then we have a banded matrix and we can use a more
!       efficient diagonalisation method.
!-----------------------------------------------------------------------        
      if (ndfl.eq.0) then
         call diagmat_banded(blckdim,matrix,dim,umat,eigval)
      else
         call diagmat(matrix,dim,umat,eigval)
      endif

!-----------------------------------------------------------------------
! Calculate and save to file the Ritz vectors
!-----------------------------------------------------------------------
!      call cpu_time(t1)
!      call calc_ritzvecs(lanunit,umat,dim,matdim,eigval)
!      call cpu_time(t2)
!      print*,"Old:",t2-t1
!
!      call cpu_time(t1)
!      call calc_ritzvecs2(lanunit,umat,dim,matdim,eigval,iblckdim)
!      call cpu_time(t2)
!      print*,"New:",t2-t1

      call calc_ritzvecs2(lanunit,umat,dim,matdim,eigval,iblckdim)

      return

    end subroutine lanczos_pseudospec

!#######################################################################

    subroutine diagmat(matrix,dim,umat,eigval)
      
      implicit none

      integer*8                  :: dim,i
      real*8, dimension(dim,dim) :: matrix,umat
      
      integer*4                  :: e2
      real*8, dimension(3*dim)   :: work
      real*8                     :: error
      real*8, dimension(dim)     :: eigval
      
      error=0
      e2=3*dim
      umat=matrix
      call dsyev('V','U',dim,umat,dim,eigval,work,e2,error)
      
      if (error.ne.0) then
         write(6,'(/,2x,3a,/)') 'Diagonalisation of the Lanczos ',&
              'representation of the Hamiltonian failed in ',&
              'subroutine diagmat.'
         STOP
      endif

      return
      
    end subroutine diagmat

!#######################################################################

    subroutine diagmat_banded(blckdim,matrix,dim,eigvec,eigval)

      implicit none
      
      integer*8                        :: blckdim,dim,i,j,iupper
      real*8, dimension(dim,dim)       :: matrix

      integer                          :: kd,ldab,error
      real*8, dimension(blckdim+1,dim) :: ab
      real*8, dimension(dim)           :: eigval
      real*8, dimension(dim,dim)       :: eigvec
      real*8, dimension(3*dim-2)       :: work

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
         write(6,'(/,2x,3a,/)') 'Diagonalisation of the Lanczos ',&
              'representation of the Hamiltonian failed in ',&
              'subroutine diagmat_banded.'
         STOP
      endif

      return

    end subroutine diagmat_banded

!#######################################################################

    subroutine calc_ritzvecs(lanunit,umat,dim,matdim,eigval)

      use constants
      use parameters, only: lancname

      implicit none

      integer                    :: i
      integer*8                  :: lanunit,dim,matdim,ritzunit,k,m
      real*8, dimension(dim,dim) :: umat
      real*8, dimension(dim)     :: eigval
      real*8, dimension(matdim)  :: lvec,ritzvec

      real*8, dimension(matdim,dim) :: lvec2
      real*8, dimension(matdim)     :: rvec

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
      open(lanunit,file='lanvecs',form='unformatted',status='old')

      ritzunit=lanunit+1
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate the Ritz vectors
!-----------------------------------------------------------------------
      ! Loop over Ritz vectors
      do i=1,dim
         ritzvec=0.0d0

         ! Calculate the ith Ritz vector
         rewind(lanunit)
         do k=1,dim ! Loop over Lanczos vectors
            read(lanunit) lvec
            
            ! The kth Lanczos vector contributes to all components
            ! of the current Ritz vector
            do m=1,matdim
               ritzvec(m)=ritzvec(m)+lvec(m)*umat(k,i)
            enddo

         enddo
           
         ! Write the current Ritz vector to file along with the
         ! corresponding Ritz value
         write(ritzunit) i,eigval(i),ritzvec
         
      enddo

!-----------------------------------------------------------------------
! Close the Lanczos and Ritz vector files
!-----------------------------------------------------------------------
      close(lanunit)
      close(ritzunit)
      
      return
      
    end subroutine calc_ritzvecs

!#######################################################################

    subroutine calc_ritzvecs2(lanunit,umat,dim,matdim,eigval,iblckdim)

      use constants
      use parameters, only: lancname

      implicit none

      integer                              :: v
      integer*8                            :: i,k,m,n,count,lanunit,&
                                              dim,matdim,ritzunit,&
                                              iblckdim,nblcks,last
      real*8, dimension(dim,dim)           :: umat
      real*8, dimension(dim)               :: eigval
      real*8, dimension(matdim)            :: lvec
      real*8, dimension(matdim,iblckdim+1) :: ritzvec

!-----------------------------------------------------------------------
! Open the Lanzcos and Ritz vector files
!-----------------------------------------------------------------------
      open(lanunit,file='lanvecs',form='unformatted',status='old')

      ritzunit=lanunit+1
      open(ritzunit,file=lancname,access='sequential',&
           form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate the Ritz vectors
!-----------------------------------------------------------------------
      ! Loop over blocks of Ritz vectors
      nblcks=ceiling(real(dim)/real(iblckdim+1))
      do i=1,nblcks-1
         
         ritzvec=0.0d0

         ! Calculate the curent block of iblckdim+1 Ritz vectors
         rewind(lanunit)
         do k=1,dim ! Loop over Lanczos vectors
            read(lanunit) lvec
            ! The kth Lanczos vector contributes to all components
            ! of the current block of Ritz vectors
            count=0
            do v=(iblckdim+1)*i-(iblckdim+1)+1,(iblckdim+1)*i
               count=count+1
               do m=1,matdim
                  ritzvec(m,count)=ritzvec(m,count)+lvec(m)*umat(k,v)
               enddo
            enddo
         enddo

         ! Write the current block of Ritz vectors to file along with
         ! the corresponding Ritz values
         count=0
         do v=(iblckdim+1)*i-(iblckdim+1)+1,(iblckdim+1)*i
            count=count+1
            write(ritzunit) v,eigval(v),ritzvec(:,count)
         enddo
         
      enddo
      
      ! Remaining n Ritz vectors from the incomplete block
      n=dim-(iblckdim+1)*(nblcks-1)
      ritzvec=0.0d0
      rewind(lanunit)
      do k=1,dim
         read(lanunit) lvec
         count=0
         do v=dim-n,dim
            count=count+1
            do m=1,matdim
               ritzvec(m,count)=ritzvec(m,count)+lvec(m)*umat(k,v)
            enddo
         enddo
      enddo
      count=0
      do v=dim-n,dim
         count=count+1
         write(ritzunit) v,eigval(v),ritzvec(:,count)
      enddo

!-----------------------------------------------------------------------
! Close the Lanczos and Ritz vector files
!-----------------------------------------------------------------------
      close(lanunit)
      close(ritzunit)

      return

    end subroutine calc_ritzvecs2

!#######################################################################

    subroutine chkortho(dim,matdim)

      implicit none

      integer*8                     :: unit,dim,matdim,i,j
      real*8, dimension(matdim,dim) :: lvec
      real*8                        :: dp
      real*8, parameter             :: tol=1d-3

      lvec=0.0d0

      unit=28
      open(unit,file='lanvecs',form='unformatted',status='old')

      do i=1,dim
         read(unit) lvec(:,i)
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

  end module band_lanczos
