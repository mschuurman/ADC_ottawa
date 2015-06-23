  module guessvecs
    
    implicit none

  contains

!#######################################################################

    subroutine adc1_guessvecs

      use constants
      use parameters
      use select_fano
      use fspace
      use fspace2
      use channels

      implicit none

      integer                              :: ndim,i,iout
      integer, dimension(:,:), allocatable :: kpq
      real(d), dimension(:,:), allocatable :: eigvec
      real(d), dimension(:), allocatable   :: eigval

      if (lcvs) then
         write(ilog,'(/,2x,a,/)') &
              'Generating guess vectors by diagonalising &
              the CVS-ADC(1) Hamiltonian'
      else
         write(ilog,'(/,2x,a,/)') &
              'Generating guess vectors by diagonalising &
              the ADC(1) Hamiltonian'
      endif

!-----------------------------------------------------------------------
! Select the ADC1 1h1p space
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))

      kpq(:,:)=-1

      if (lcvs) then
         call select_atom_is_cvs(kpq(:,:))
      else
         call select_atom_is(kpq(:,:))
      endif

      ndim=kpq(1,0)

!-----------------------------------------------------------------------
! Allocate the arrays that will hold the ADC1 eigenpairs
!-----------------------------------------------------------------------
      allocate(eigvec(ndim,ndim))
      allocate(eigval(ndim))

!-----------------------------------------------------------------------
! Diagonalise the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------
      if (lcvs) then
         call get_fspace_tda_direct_cvs(ndim,kpq(:,:),eigvec,eigval)
      else
         call get_fspace_tda_direct(ndim,kpq(:,:),eigvec,eigval)
      endif

      deallocate(kpq)

!-----------------------------------------------------------------------
! Write the ADC(1) eigenvectors to file
!-----------------------------------------------------------------------
      iout=22

      open(iout,file='SCRATCH/adc1_vecs',form='unformatted',&
           status='unknown')

      write(iout) ndim,eigvec

      close(iout)

      return

    end subroutine adc1_guessvecs

!#######################################################################

    subroutine adc1_guessvecs_final

      use constants
      use parameters
      use select_fano
      use fspace
      use fspace2
      use channels

      implicit none

      integer                              :: ndim,i,iout
      integer, dimension(:,:), allocatable :: kpqf
      real(d), dimension(:,:), allocatable :: eigvec
      real(d), dimension(:), allocatable   :: eigval

      if (lcvsfinal) then
         write(ilog,'(/,2x,a,/)') &
              'Generating guess vectors by diagonalising &
              the CVS-ADC(1) Hamiltonian'
      else
         write(ilog,'(/,2x,a,/)') &
              'Generating guess vectors by diagonalising &
              the ADC(1) Hamiltonian'
      endif

!-----------------------------------------------------------------------
! Select the ADC1 1h1p space
!-----------------------------------------------------------------------
      allocate(kpqf(7,0:nBas**2*4*nOcc**2))

      kpqf(:,:)=-1

      if (lcvsfinal) then
         call select_atom_is_cvs(kpqf(:,:))
      else
         call select_atom_isf(kpqf(:,:))
      endif

      ndim=kpqf(1,0)

!-----------------------------------------------------------------------
! Allocate the arrays that will hold the ADC1 eigenpairs
!-----------------------------------------------------------------------
      allocate(eigvec(ndim,ndim))
      allocate(eigval(ndim))

!-----------------------------------------------------------------------
! Diagonalise the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------
      if (lcvsfinal) then
         call get_fspace_tda_direct_cvs(ndim,kpqf(:,:),eigvec,eigval)
      else
         call get_fspace_tda_direct(ndim,kpqf(:,:),eigvec,eigval)
      endif

      deallocate(kpqf)

!-----------------------------------------------------------------------
! Write the ADC(1) eigenvectors to file
!-----------------------------------------------------------------------
      iout=22

      open(iout,file='SCRATCH/adc1_vecs',form='unformatted',&
           status='unknown')

      write(iout) ndim,eigvec

      close(iout)

      return

    end subroutine adc1_guessvecs_final

!#######################################################################

    subroutine get_fakeip_indices(kpq,kpqdim2,dims,dim)
      
      use constants
      use parameters
      use misc
      use channels

      implicit none
      
      integer                              :: kpqdim2,dims,inda,indb,&
                                              indj,indk,spin,i,j,k,&
                                              norb,ndim,a,count,dim
      integer, dimension(7,0:kpqdim2-1)    :: kpq

      real(d), dimension(:), allocatable   :: apprx_en
      integer, dimension(:,:), allocatable :: indxaji
      integer, dimension(:), allocatable   :: sortindx,map

!-----------------------------------------------------------------------
! Fill in the ifakeex array: ifakeex(k) is the index of the single
! non-zero element of the unit vector corresponding to the kth guess
! vector
!-----------------------------------------------------------------------
      nvirt=nbas-nocc
      ndim=nocc+(nvirt-1)*nocc**2

      allocate(ifakeex(dmain))
      allocate(apprx_en(ndim))
      allocate(sortindx(ndim))
      allocate(map(ndim))
      allocate(indxaji(ndim,3))
      
      indxaji=0
      apprx_en=9999999.9d0

      ! (1) i->inf
      count=0
      do i=1,dims
         call get_indices(kpq(:,i),inda,indb,indj,indk,spin)         
         if (inda.eq.ifakeorb) then
            count=count+1
            indxaji(i,3)=indj
            apprx_en(i)=-e(indj)
            map(count)=i
         endif
      enddo

      ! (2) i->inf, j->b
      do i=dims+1,dim
         call get_indices(kpq(:,i),inda,indb,indj,indk,spin)            
         if (inda.eq.ifakeorb.and.indb.ne.ifakeorb) then
            count=count+1
            indxaji(count,3)=indj
            indxaji(count,2)=indk
            indxaji(count,1)=indb
            apprx_en(count)=-e(indj)+e(indb)-e(indk)           
            map(count)=i
         else if (indb.eq.ifakeorb.and.inda.ne.ifakeorb) then
            count=count+1
            indxaji(count,3)=indj
            indxaji(count,2)=indk
            indxaji(count,1)=inda
            apprx_en(count)=-e(indj)+e(inda)-e(indk)           
            map(count)=i
         endif
      enddo

      ! Sort the approximate ionization energies and save the indices
      ! to sortindx in order of increasing energy
      call dsortindxa1('A',ndim,apprx_en,sortindx)

!         do i=1,10
!            k=sortindx(i)
!            print*,k,apprx_en(k)*27.211,&
!                 indxaji(k,3),indxaji(k,2),indxaji(k,1),map(k)
!         enddo

      do i=1,dmain
         k=sortindx(i)
         ifakeex(i)=map(k)
      enddo

      return

    end subroutine get_fakeip_indices

!#######################################################################

    subroutine getblocksize(travec,ndimf,ndimsf)

      use channels
      use constants
      use parameters
      use misc, only: dsortindxa1
      use iomod, only: freeunit

      implicit none

      integer                              :: ndimf,ndimsf,i,k,iadc1,&
                                              itmp,ndim2,k1,k2
      integer, dimension(:), allocatable   :: indx1,indx2
      real(d), dimension(ndimf)            :: travec
      real(d), dimension(:,:), allocatable :: adc1vec
      real(d), dimension(:), allocatable   :: tmpvec
      real(d)                              :: ftmp,tol

      lmain=0

      if (lancguess.eq.3.or.lancguess.eq.4) then
         tol=tdtol/sqrt(2.0d0)
      else
         tol=tdtol
      endif

      if (lancguess.eq.1) then         
         if (method.eq.2) then
            allocate(indx1(ndimsf))
            call dsortindxa1('D',ndimsf,travec(1:ndimsf)**2,indx1(:))
            do i=1,ndimsf
               k=indx1(i)            
               if (abs(travec(k)).gt.tol.and.lmain.lt.maxblock) &
                    lmain=lmain+1
            enddo
         else if (method.eq.3) then            
            allocate(indx1(ndimf))
            call dsortindxa1('D',ndimf,travec(1:ndimf)**2,indx1(:))
            do i=1,ndimf
               k=indx1(i)            
               if (abs(travec(k)).gt.tol.and.lmain.lt.maxblock) &
                    lmain=lmain+1
            enddo
         endif

      else if (lancguess.eq.2) then
         allocate(indx1(ndimsf))
         allocate(adc1vec(ndimsf,ndimsf))
         allocate(tmpvec(ndimsf))
         call freeunit(iadc1)
         open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
              status='old')
         read(iadc1) itmp,adc1vec
         close(iadc1)
         do i=1,ndimsf
            tmpvec(i)=dot_product(adc1vec(:,i),travec(1:ndimsf))
         enddo
         call dsortindxa1('D',ndimsf,tmpvec(1:ndimsf)**2,indx1(:))
         do i=1,ndimsf
            k=indx1(i)            
            if (abs(tmpvec(k)).gt.tol.and.lmain.lt.maxblock) &
                 lmain=lmain+1
         enddo
         deallocate(adc1vec)
         deallocate(tmpvec)

      else if (lancguess.eq.3) then
         ndim2=ndimf-ndimsf
         allocate(indx1(ndimsf))
         allocate(indx2(ndim2))
         call dsortindxa1('D',ndimsf,travec(1:ndimsf)**2,indx1(:))
         call dsortindxa1('D',ndim2,travec(ndimsf+1:ndimf)**2,indx2(:))
         do i=1,ndimsf
            k1=indx1(i)
            k2=ndimsf+indx2(i)
            ftmp=(1.0d0/sqrt(2.0d0))*abs(travec(k1))+abs(travec(k2))
            if (ftmp.gt.tol.and.lmain.lt.maxblock) lmain=lmain+1
         enddo
         deallocate(indx2)

      else if (lancguess.eq.4) then
         ndim2=ndimf-ndimsf
         allocate(indx1(ndimsf))
         allocate(indx2(ndim2))
         allocate(adc1vec(ndimsf,ndimsf))
         allocate(tmpvec(ndimsf))
         call freeunit(iadc1)
         open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
              status='old')
         read(iadc1) itmp,adc1vec
         close(iadc1)
         do i=1,ndimsf
            tmpvec(i)=dot_product(adc1vec(:,i),travec(1:ndimsf))
         enddo
         call dsortindxa1('D',ndimsf,tmpvec(1:ndimsf)**2,indx1(:))
         call dsortindxa1('D',ndim2,travec(ndimsf+1:ndimf)**2,indx2(:))
         do i=1,ndimsf
            k1=indx1(i)
            k2=indx2(i)
            ftmp=(1.0d0/sqrt(2.0d0))*abs(tmpvec(k1))+abs(travec(k2))
            if (ftmp.gt.tol.and.lmain.lt.maxblock) lmain=lmain+1
         enddo
         deallocate(adc1vec)
         deallocate(tmpvec)

      endif

      deallocate(indx1)

      write(ilog,'(/,2x,a,2x,i4,/)') &
           'Lanczos block size selected:',lmain

      ! Reset ncycles
      ncycles=maxblock*ncycles/lmain

      return

    end subroutine getblocksize

!#######################################################################

  end module guessvecs
