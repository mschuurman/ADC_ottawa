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

      implicit none

      integer                              :: ndim,i,iout
      integer, dimension(:,:), allocatable :: kpq
      real(d), dimension(:,:), allocatable :: eigvec
      real(d), dimension(:), allocatable   :: eigval

!-----------------------------------------------------------------------
! Select initial 1h1p space
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))
      kpq(:,:)=-1
      if (lcvs) then
         call select_atom_is_cvs(kpq(:,:))
      else
         call select_atom_is(kpq(:,:))
      endif

!-----------------------------------------------------------------------
! Output initial 1h1p space information
!-----------------------------------------------------------------------
      ndim=kpq(1,0)
      write(6,*) 'ADC(1) INITIAL Space dim',ndim
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(eigvec(ndim,ndim))
      allocate(eigval(ndim))

!-----------------------------------------------------------------------
! Diagonalise the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------
      write(6,'(/,78a)') ('-',i=1,78)
      if (lcvs) then
         write(6,'(2x,a)') 'Generating guess Davidson vectors by diagonalising &
              the CVS-ADC(1) Hamiltonian'
      else
         write(6,'(2x,a)') 'Generating guess Davidson vectors by diagonalising &
              the ADC(1) Hamiltonian'
      endif

      write(6,'(78a,/)') ('-',i=1,78)

      if (lcvs) then
         call get_fspace_tda_direct_cvs(ndim,kpq(:,:),eigvec,eigval)
      else
         call get_fspace_tda_direct(ndim,kpq(:,:),eigvec,eigval)
      endif

!-----------------------------------------------------------------------
! Write the ADC(1) eigenvectors to file
!-----------------------------------------------------------------------
      iout=22
      open(iout,file='adc1_vecs',form='unformatted',status='unknown')
      write(iout) ndim,eigvec
      close(iout)

      return

    end subroutine adc1_guessvecs

!#######################################################################

    subroutine get_fakeip_indices(kpq,kpqdim2,dims,dim)
      
      use constants
      use parameters
      use misc

      implicit none
      
      integer                              :: kpqdim2,dims,inda,indb,&
                                              indj,indk,spin,i,j,k,&
                                              nex,norb,ndim,a,count,dim
      integer, dimension(7,0:kpqdim2-1)    :: kpq
      integer, dimension(nocc,2)           :: exindx
      real(d), dimension(:), allocatable   :: apprx_en
      integer, dimension(:,:), allocatable :: indxaji
      integer, dimension(:), allocatable   :: sortindx,map

!-----------------------------------------------------------------------
! Determine the indices of the 1h1p configurations corresponding to
! excitation into the diffuse orbital
!-----------------------------------------------------------------------
      exindx=0
      nex=0
      do i=1,dims
         call get_indices(kpq(:,i),inda,indb,indj,indk,spin)
         if (inda.eq.ifakeorb) then
            nex=nex+1
            exindx(nex,1)=i
            exindx(nex,2)=indj
         endif         
      enddo

!-----------------------------------------------------------------------
! Fill in the ifakeex array: ifakeex(k) is the index of the singl 
! non-zero element of the unit vector corresponding to the kth guess
! vector
!-----------------------------------------------------------------------
      allocate(ifakeex(dmain))

      if (lcvs) then
         norb=ncore
         if (dmain.gt.ncore) then
            write(6,'(/,2(2x,a,/))') 'Currently only 1h1p initial vectors &
                 are supported for a fake IP calculation.','Reduce &
                 dmain so that it is less than or equal to ncore and &
                 re-run'
            STOP
         endif
         do i=1,dmain
            k=norb-i+1
            ifakeex(i)=exindx(k,1)
         enddo
      else
         nvirt=nbas-nocc
         ndim=nocc+(nvirt-1)*nocc**2
         allocate(apprx_en(ndim))
         allocate(sortindx(ndim))
         allocate(map(ndim))
         allocate(indxaji(ndim,3))
         indxaji=0
         apprx_en=0.0d0

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

         call dsortindxa1('A',ndim,apprx_en,sortindx)

!         do i=1,100
!            k=sortindx(i)
!            print*,k,apprx_en(k)*27.211,&
!                 indxaji(k,3),indxaji(k,2),indxaji(k,1),map(k)
!         enddo

         do i=1,dmain
            k=sortindx(i)
            ifakeex(i)=map(k)
         enddo
         
      endif

      return

    end subroutine get_fakeip_indices

!#######################################################################

  end module guessvecs
