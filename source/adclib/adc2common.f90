!#######################################################################
! Subroutines and functions common to all types of ADC(2) spectroscopy
! calculations
!#######################################################################

  module adc2common

    implicit none

  contains

!#######################################################################

    subroutine checksymm

      use parameters
      use iomod
        
      implicit none

!-----------------------------------------------------------------------
! Symmetry is not supported for some types of calculation.
! Check here that everything is OK in this regard.
!-----------------------------------------------------------------------
      if (lrixs.or.ltpa.or.lpropagation) then
         if (pntgroup.ne.'C1') then
            errmsg='Symmetry is NOT supported in RIXS or TPA &
                 calculations'
            call error_control
         endif
      endif
    
      return
    
    end subroutine checksymm

!#######################################################################

    subroutine set_dpl

      use parameters

      implicit none

!-----------------------------------------------------------------------
! Set the MO representation of the dipole operator
!-----------------------------------------------------------------------
      if (tranmom2 .eq. 'x') then
         dpl(:,:)=x_dipole(:,:)
      elseif (tranmom2 .eq. 'y') then
         dpl(:,:)=y_dipole(:,:)
      elseif (tranmom2 .eq. 'z') then
         dpl(:,:)=z_dipole(:,:)
      end if

!-----------------------------------------------------------------------
! Set the irrep of the dipole operator
!-----------------------------------------------------------------------
      CHECK_dip = nirrep2

!-----------------------------------------------------------------------
! Set up the total dipole matrix array
!-----------------------------------------------------------------------
      allocate(dpl_all(3,nbas,nbas))
      dpl_all(1,:,:)=x_dipole(:,:)
      dpl_all(2,:,:)=y_dipole(:,:)
      dpl_all(3,:,:)=z_dipole(:,:)

      return
      
    end subroutine set_dpl

!#######################################################################

    subroutine get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,&
         noutf,ndims,ndimsf)
      
      use parameters
      use iomod
      use select_fano

      implicit none

      integer                              :: ndim,ndimf,ndimd,nout,&
                                              noutf,ndims,ndimsf
      integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
  
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))
      allocate(kpqd(7,0:nBas**2*4*nOcc**2))
      allocate(kpqf(7,0:nBas**2*4*nOcc**2))
        
!-----------------------------------------------------------------------
! Determine the initial, final and total subspaces
!-----------------------------------------------------------------------
      ! Initial subspace
      kpq(:,:)=-1
      call select_atom_is(kpq(:,:))
      call select_atom_d(kpq(:,:),-1)
      
      ! Final subspace
      kpqf(:,:)=-1
      if (lcvsfinal) then
         call select_atom_isf_cvs(kpqf(:,:))
         call select_atom_df_cvs(kpqf(:,:),-1)
      else
         call select_atom_isf(kpqf(:,:))
         call select_atom_df(kpqf(:,:),-1)
      endif
           
      ! Total subspace
      kpqd(:,:)=-1
      call select_atom_ist(kpqd(:,:))
      call select_atom_dt(kpqd(:,:),-1)
           
!-----------------------------------------------------------------------
! Determine the various subspace dimensions
!-----------------------------------------------------------------------
      ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
      ndimf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)
      ndimd=kpqd(1,0)+kpqd(2,0)+kpqd(3,0)+kpqd(4,0)+2*kpqd(5,0)
      
      nout=ndim
      noutf=ndimf
      ndims=kpq(1,0)
      ndimsf=kpqf(1,0)

!-----------------------------------------------------------------------
! Output the subspace information
!-----------------------------------------------------------------------
      write(ilog,*) 'ADC(2) INITIAL Space dim',ndim
      write(ilog,*) 'ADC(2) FINAL Space dim',ndimf
      write(ilog,*) 'ADC(2) TOTAL Space dim WITHOUT GROUND STATE',ndimd
      write(ilog,*) 'dimension of various INITIAL configuration spaces'
      write(ilog,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0)
      write(ilog,*) 'dimension of various FINAL configuration spaces'
      write(ilog,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
      write(ilog,*) 'dimension of various TOTAL configuration spaces'
      write(ilog,*) kpqd(1,0),kpqd(2,0),kpqd(3,0),kpqd(4,0),kpqd(5,0)
    
      return

    end subroutine get_subspaces

!#######################################################################

    subroutine get_subspaces_adc1(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,&
         noutf)

      use parameters
      use iomod
      use select_fano

      implicit none

      integer                              :: ndim,ndimf,ndimd,nout,&
                                              noutf,ndims,ndimsf
      integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))
      allocate(kpqd(7,0:nBas**2*4*nOcc**2))
      allocate(kpqf(7,0:nBas**2*4*nOcc**2))
      
!-----------------------------------------------------------------------
! Determine the initial, final and total subspaces
!-----------------------------------------------------------------------
      ! Initial subspace
      kpq(:,:)=-1
      call select_atom_is(kpq(:,:))

      ! Final subspace
      kpqf(:,:)=-1
      if (lcvsfinal) then
         call select_atom_isf_cvs(kpqf(:,:))
      else
         call select_atom_isf(kpqf(:,:))
      endif
      
      ! Total subspace
      kpqd(:,:)=-1
      call select_atom_ist(kpqd(:,:))

!-----------------------------------------------------------------------
! Determine the various subspace dimensions
!-----------------------------------------------------------------------
      ndim=kpq(1,0)
      ndimf=kpqf(1,0)
      ndimd=kpqd(1,0)

      nout=ndim
      noutf=ndimf

!-----------------------------------------------------------------------
! Output the subspace information
!-----------------------------------------------------------------------
      write(ilog,'(/)')
      write(ilog,*) 'ADC(1) INITIAL Space dim',ndim
      write(ilog,*) 'ADC(1) FINAL Space dim',ndimf
      write(ilog,*) 'ADC(1) TOTAL Space dim WITHOUT GROUND STATE',ndimd
      
      return
      
    end subroutine get_subspaces_adc1
!#######################################################################

    subroutine init_space_diag(time,kpq,ndim,ndims,noffd)

      use constants
      use parameters
      use fspace
!      use diagmod
      use select_fano
      use misc
      use dmatvec_davidson
      use get_matrix

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,ndims,ndim1,ndim2,nbuf
      integer*8                                 :: noffd
      real(d)                                   :: time
      character(len=120)                        :: msg

!-----------------------------------------------------------------------
! Write the initial space ADC(2) Hamiltonian to disk
!-----------------------------------------------------------------------
      if (method.eq.2) then
         write(ilog,*) 'Calculating the initial space ADC(2)-s &
              matrix-vector multiplication'
      else
         write(ilog,*) 'Only available for ADC(2)-s'
         stop
         call error_control
      endif


      call get_interm_adc2_save(ndim,kpq(:,:),'i')

      ndim1=kpq(1,0)
      ndim2=ndim-kpq(1,0)

!      call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,'i')
!
!      write(ilog,'(/,a)') trim(msg)

!      if (method.eq.2) then
         ! ADC(2)-s
         call davdir_block(ndim,kpq(:,:),'i')
!      else if (method.eq.3) then
         ! ADC(2)-x
!         call write_fspace_adc2e_1(ndim,kpq(:,:),noffd,'i')
!      endif

      call cpu_time(time)

      write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time=',time," s"

!-----------------------------------------------------------------------
! Diagonalisation
!-----------------------------------------------------------------------
!      call master_eig(ndim,noffd,'i')

      return

    end subroutine init_space_diag

!#######################################################################

    subroutine davidson_fin_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,&
         travec,vec_init,mtmf,noffdf,rvec,travec2)

      use constants
      use parameters
      use fspace
      use get_matrix_dipole
!      use diagmod
      use guessvecs
      use dmatvec_davidson
      use get_matrix
      use select_fano
      use misc

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndimf,ndimsf,ndim1,ndim2,nbuf
      integer*8                                 :: noffdf
      real(d), dimension(:), allocatable        :: travec,mtmf
      real(d), dimension(ndim)                  :: vec_init
      real(d), dimension(ndim,davstates)        :: rvec
      real(d), dimension(:,:), allocatable      :: travec2

!-----------------------------------------------------------------------        
! If requested, determine the Davidson guess vectors by diagonalising 
! the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------        
      if (ladc1guess_f) call adc1_guessvecs_final

!-----------------------------------------------------------------------
! Write the final space ADC(2) Hamiltonian to disk
!-----------------------------------------------------------------------
      write(ilog,*) 'FINAL SPACE ADC2 calculation'

      if (method.eq.2) then
         write(ilog,*) 'Calculating the final space ADC(2)-s &
              matrix-vector product'
      else
         write(ilog,*) 'Only available for ADC(2)-s'
         stop
         call error_control
      endif

      call get_interm_adc2_save(ndimf,kpqf(:,:),'f')

      ndim1=kpq(1,0)
      ndim2=ndimf-kpq(1,0)

!      call get_diag_adc2_save(ndim1,ndim2,kpqf,nbuf,'c')
!
!      if (method_f.eq.2) then
!         ! ADC(2)-s
!         if (lcvsfinal) then
!            call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
!         else
!            call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
!         endif
!      else if (method_f.eq.3) then
!         ! ADC(2)-x
!         if (lcvsfinal) then
!            call write_fspace_adc2e_1_cvs(ndimf,kpqf(:,:),noffdf,'c')
!         else
!            call write_fspace_adc2e_1(ndimf,kpqf(:,:),noffdf,'c')
!         endif
!      endif

!-----------------------------------------------------------------------
! Diagonalisation in the final space
!-----------------------------------------------------------------------
      call davdir_block(ndimf,kpqf(:,:),'f')

!-----------------------------------------------------------------------
! Allocate the travec array that will hold the contraction of the IS
! representation of the shifted dipole operator with the initial
! state vector
!-----------------------------------------------------------------------
      if (statenumber.eq.0) then
         allocate(mtmf(ndimf))
      else
         allocate(travec(ndimf))
      endif

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the 
! shifted dipole operator and the initial state vector.
!
! This will be used later in the calculation of transition dipole
! matrix elements between the initial and final states.
!
! Really this should be done elsewhere...
!-----------------------------------------------------------------------
      if (statenumber.eq.0) then
         call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:),0)
      else
         call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,&
              vec_init,travec)
      endif

      return

    end subroutine davidson_fin_space_diag

      
!#######################################################################

    subroutine initial_space_diag(time,kpq,ndim,ndims,noffd)
        
      use constants
      use parameters
      use fspace
      use diagmod
        
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,ndims
      integer*8                                 :: noffd
      real(d)                                   :: time
      character(len=120)                        :: msg

!-----------------------------------------------------------------------
! Write the initial space ADC(2) Hamiltonian to disk
!-----------------------------------------------------------------------
      if (method.eq.2) then
         msg='Calculating the initial space ADC(2)-s Hamiltonian &
              matrix'
      else if (method.eq.3) then
         msg='Calculating the initial space ADC(2)-x Hamiltonian &
              matrix'
      endif

      write(ilog,'(/,a)') trim(msg)

      if (method.eq.2) then
         ! ADC(2)-s
         call write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i')
      else if (method.eq.3) then
         ! ADC(2)-x
         call write_fspace_adc2e_1(ndim,kpq(:,:),noffd,'i')
      endif
      
      call cpu_time(time)

      write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time=',time," s"

!-----------------------------------------------------------------------
! Diagonalisation
!-----------------------------------------------------------------------
      call master_eig(ndim,noffd,'i')
      
      return

    end subroutine initial_space_diag

!#######################################################################

    subroutine davidson_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,&
         travec,vec_init,mtmf,noffdf,rvec,travec2)

      use constants
      use parameters
      use fspace
      use get_matrix_dipole
      use diagmod
      use guessvecs
        
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndimf,ndimsf
      integer*8                                 :: noffdf
      real(d), dimension(:), allocatable        :: travec,mtmf
      real(d), dimension(ndim)                  :: vec_init
      real(d), dimension(ndim,davstates)        :: rvec
      real(d), dimension(:,:), allocatable      :: travec2

!-----------------------------------------------------------------------        
! If requested, determine the Davidson guess vectors by diagonalising 
! the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------        
      if (ladc1guess_f) call adc1_guessvecs_final

!-----------------------------------------------------------------------
! Write the final space ADC(2) Hamiltonian to disk
!-----------------------------------------------------------------------
      write(ilog,*) 'Saving complete FINAL SPACE ADC2 matrix in file'

      if (method_f.eq.2) then
         ! ADC(2)-s
         if (lcvsfinal) then
            call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')           
         else
            call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
         endif
      else if (method_f.eq.3) then
         ! ADC(2)-x
         if (lcvsfinal) then
            call write_fspace_adc2e_1_cvs(ndimf,kpqf(:,:),noffdf,'c')           
         else
            call write_fspace_adc2e_1(ndimf,kpqf(:,:),noffdf,'c')
         endif
      endif

!-----------------------------------------------------------------------
! Diagonalisation in the final space
!-----------------------------------------------------------------------
      call master_eig(ndimf,noffdf,'f')
        
!-----------------------------------------------------------------------
! Allocate the travec array that will hold the contraction of the IS
! representation of the shifted dipole operator with the initial
! state vector
!-----------------------------------------------------------------------
      if (statenumber.eq.0) then
         allocate(mtmf(ndimf))
      else
         allocate(travec(ndimf))
      endif

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the 
! shifted dipole operator and the initial state vector.
!
! This will be used later in the calculation of transition dipole
! matrix elements between the initial and final states.
!
! Really this should be done elsewhere...
!-----------------------------------------------------------------------
      if (statenumber.eq.0) then
         call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:),0)
      else
         call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,&
              vec_init,travec)
      endif

      return

    end subroutine davidson_final_space_diag
    
!#######################################################################

    subroutine initial_space_dipole(ndim,ndims,kpq)

      use channels
      use constants
      use parameters
      use diagmod, only: readdavvc
      use get_matrix_dipole
      use misc
      use mp2
      use timingmod

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,ndims,i,j,&
                                                   k,a,b,p,q,c
      integer*8                                 :: nel
      integer                                   :: nbuf
      real(d), dimension(:,:), allocatable      :: rvec,dmvec
      real(d), dimension(davstates_f)           :: ener
      real(d)                                   :: dip0,dipnuc
      real(d), parameter                        :: ang2bohr=1.889725989d0
      real(d)                                   :: tw1,tw2,tc1,tc2
      character(len=60)                         :: filename

!-----------------------------------------------------------------------
! Output what we are doing
!-----------------------------------------------------------------------
      write(ilog,'(/,70a)') ('-',i=1,70)
      write(ilog,'(2x,a)') 'Calculating the initial space dipole &
           moments'
      write(ilog,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
      call times(tw1,tc1)

!-----------------------------------------------------------------------
! Nuclear contribution to the dipole moments
!-----------------------------------------------------------------------
      if (tranmom2.eq.'x') then
         c=1
      else if (tranmom2.eq.'y') then
         c=2
      else if (tranmom2.eq.'z') then
         c=3
      endif

      dipnuc=0.0d0
      do i=1,natm
         dipnuc=dipnuc+ang2bohr*xcoo(i*3-3+c)*atomic_number(atlbl(i))
      enddo
      
!-----------------------------------------------------------------------
! Ground state density matrix.
! Note that the 1st-order correction is zero.
!-----------------------------------------------------------------------
      if (.not.allocated(density)) allocate(density(nbas,nbas))
      density=0.0d0

      call rho_mp2(density)

!-----------------------------------------------------------------------
! Ground state electronic dipole moment
!-----------------------------------------------------------------------
      dip0=0.0d0
      do p=1,nbas
         do q=1,nbas
            dip0=dip0+dpl(p,q)*density(p,q)
         enddo
      enddo
      
      ! Multiplication by e
      dip0=-dip0

      write(ilog,*) "dip0:",dip0+dipnuc

!-----------------------------------------------------------------------
! Read the initial space Davidson states from file
!-----------------------------------------------------------------------
      allocate(rvec(ndim,davstates))
      call readdavvc(davstates,ener,rvec,'i',ndim)

!-----------------------------------------------------------------------
! Form the dmvec vectors: contractions of the IS representation of the
! shifted dipole operator with the initial state vectors
!-----------------------------------------------------------------------
      allocate(dmvec(ndim,davstates))

      write(ilog,'(/,2x,a,/)') "Calculating the initial state dipole &
           moments..."

      filename='SCRATCH/dipole_vv'

      call get_adc2_dipole_improved_omp(ndim,ndim,kpq,kpq,nbuf,&
           nel,filename)

      do i=1,davstates
         call contract_dipole_state(filename,ndim,ndim,rvec(:,i),&
              dmvec(:,i),nbuf,nel,'r')
      enddo

!-----------------------------------------------------------------------
! Contract the dmvec vectors with the final state vectors to yield the
! final state electronic dipole moments
!-----------------------------------------------------------------------
      allocate(dipmom(davstates))

      ! Loop over Davidson states
      do i=1,davstates
         
         ! < Psi_i | D - D0 | Psi_i > (including the multiplication
         ! by e)
         dipmom(i)=-dot_product(rvec(:,i),dmvec(:,i))
         
         ! < Psi_i | D | Psi_i >
         dipmom(i)=dipmom(i)+dip0
         
      enddo

!-----------------------------------------------------------------------
! Total excited state dipole moments
!-----------------------------------------------------------------------
      dipmom=dipmom+dipnuc

!-----------------------------------------------------------------------
! Finish timing and output the walltime taken
!-----------------------------------------------------------------------
      call times(tw2,tc2)
      write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(rvec)
      deallocate(dmvec)

      return
        
    end subroutine initial_space_dipole

!#######################################################################

    subroutine final_space_dipole(ndimf,ndimsf,kpqf)
        
      use channels
      use constants
      use parameters
      use diagmod, only: readdavvc
      use get_matrix_dipole
      use misc
      use mp2
      use timingmod

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
      integer                                   :: ndimf,ndimsf,i,j,&
                                                   k,a,b,p,q,c
      integer*8                                 :: nel
      integer                                   :: nbuf
      real(d), dimension(:,:), allocatable      :: rvec,dmvec
      real(d), dimension(davstates_f)           :: ener
      real(d)                                   :: dip0,dipnuc
      real(d), parameter                        :: ang2bohr=1.889725989d0
      real(d)                                   :: tw1,tw2,tc1,tc2
      character(len=60)                         :: filename
      
!-----------------------------------------------------------------------
! Output what we are doing
!-----------------------------------------------------------------------
      write(ilog,'(/,70a)') ('-',i=1,70)
      write(ilog,'(2x,a)') 'Calculating the final space dipole moments'
      write(ilog,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
      call times(tw1,tc1)

!-----------------------------------------------------------------------
! Nuclear contribution to the dipole moments
!-----------------------------------------------------------------------
      if (tranmom2.eq.'x') then
         c=1
      else if (tranmom2.eq.'y') then
         c=2
      else if (tranmom2.eq.'z') then
         c=3
      endif

      dipnuc=0.0d0
      do i=1,natm
         dipnuc=dipnuc+ang2bohr*xcoo(i*3-3+c)*atomic_number(atlbl(i))
      enddo

!-----------------------------------------------------------------------
! Ground state density matrix
!-----------------------------------------------------------------------
      if (.not.allocated(density)) allocate(density(nbas,nbas))
      density=0.0d0

      call rho_mp2(density)
      
!-----------------------------------------------------------------------
! Ground state electronic dipole moment
!-----------------------------------------------------------------------
      dip0=0.0d0
      do p=1,nbas
         do q=1,nbas
            dip0=dip0+dpl(p,q)*density(p,q)
         enddo
      enddo

      ! Multiplication by e
      dip0=-dip0

      !write(ilog,*) "dip0:",dip0

!-----------------------------------------------------------------------
! Read the final space Davidson states from file
!-----------------------------------------------------------------------
      allocate(rvec(ndimf,davstates_f))
      call readdavvc(davstates_f,ener,rvec,'f',ndimf)

!-----------------------------------------------------------------------
! Form the dmvec vectors: contractions of the IS representation of the
! shifted dipole operator with the final state vectors
! -----------------------------------------------------------------------
      allocate(dmvec(ndimf,davstates_f))

      write(ilog,'(/,2x,a,/)') "Calculating the final state dipole &
           moments..."

      filename='SCRATCH/dipole_cc'

      call get_adc2_dipole_improved_omp(ndimf,ndimf,kpqf,kpqf,nbuf,&
           nel,filename)

      do i=1,davstates_f
         call contract_dipole_state(filename,ndimf,ndimf,rvec(:,i),&
                dmvec(:,i),nbuf,nel,'r')
      enddo

!-----------------------------------------------------------------------
! Contract the dmvec vectors with the final state vectors to yield the
! final state dipole moments
!-----------------------------------------------------------------------
      allocate(dipmom_f(davstates_f))

      ! Loop over Davidson states
      do i=1,davstates_f
         
         ! < Psi_i | D - D0 | Psi_i > (including the multiplication
         ! by e)
         dipmom_f(i)=-dot_product(rvec(:,i),dmvec(:,i))
         
         ! < Psi_i | D | Psi_i >
         dipmom_f(i)=dipmom_f(i)+dip0

      enddo

!-----------------------------------------------------------------------
! Total excited state dipole moments
!-----------------------------------------------------------------------
      dipmom_f=dipmom_f+dipnuc

!-----------------------------------------------------------------------
! Finish timing and output the walltime taken
!-----------------------------------------------------------------------
      call times(tw2,tc2)
      write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(rvec)
      deallocate(dmvec)

      return

    end subroutine final_space_dipole
    
!#######################################################################
    
    subroutine contract_dipole_state(filename,dim1,dim2,vec,tvec,&
         nbuf,nel,cntrdir)

      use constants
      use iomod
      use omp_lib

      implicit none

      integer                              :: dim1,dim2,idpl,nbuf,&
                                              nlim,k,n,buffsize,&
                                              nthreads,tid
      integer*8                            :: nel
      integer, dimension(:), allocatable   :: indxi,indxj
      real(d), dimension(dim1)             :: vec
      real(d), dimension(dim2)             :: tvec
      real(d), dimension(:), allocatable   :: buffer
      real(d), dimension(:,:), allocatable :: tmpvec
      character(len=60)                    :: filename
      character(len=1)                     :: cntrdir
        
!-----------------------------------------------------------------------
! Open the dipole matrix file
!-----------------------------------------------------------------------
      call freeunit(idpl)
      open(idpl,file=filename,status='unknown',access='sequential',&
           form='unformatted')

!-----------------------------------------------------------------------
! Read the buffer size
!-----------------------------------------------------------------------
      read(idpl) buffsize

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(buffer(buffsize))
      allocate(indxi(buffsize))
      allocate(indxj(buffsize))
      
      !$omp parallel
      nthreads=omp_get_num_threads()
      !$omp end parallel
      allocate(tmpvec(dim2,nthreads))

!-----------------------------------------------------------------------
! Perform the contraction of the shifted dipole matrix with the state
! vector
!
! Note that the cntrdir flag controls whether we operate with the
! shifted dipole operator from the right or the left of the state
! vector
!-----------------------------------------------------------------------
      tvec=0.0d0
      tmpvec=0.0d0
      
      if (cntrdir.eq.'r') then
         do k=1,nbuf
            read(idpl) buffer(:),indxi(:),indxj(:),nlim
            !$omp parallel do &
            !$omp& private(n) &
            !$omp& shared(tmpvec,buffer,vec,indxi,indxj)
            do n=1,nlim
               tid=1+omp_get_thread_num()
               tmpvec(indxi(n),tid)=tmpvec(indxi(n),tid)&
                    +buffer(n)*vec(indxj(n))
            enddo
            !$omp end parallel do
         enddo
      else if (cntrdir.eq.'l') then
         do k=1,nbuf
            read(idpl) buffer(:),indxi(:),indxj(:),nlim
            !$omp parallel do &
            !$omp& private(n) &
            !$omp& shared(tmpvec,buffer,vec,indxi,indxj)
            do n=1,nlim
               tid=1+omp_get_thread_num()
               tmpvec(indxj(n),tid)=tmpvec(indxj(n),tid)&
                    +buffer(n)*vec(indxi(n))
            enddo
            !$omp end parallel do
         enddo
      endif

      do k=1,nthreads
         tvec=tvec+tmpvec(:,k)
      enddo

!-----------------------------------------------------------------------
! Close the dipole matrix file
!-----------------------------------------------------------------------
      close(idpl)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(buffer)
      deallocate(indxi)
      deallocate(indxj)
      deallocate(tmpvec)
      
      return

    end subroutine contract_dipole_state

!#######################################################################

    subroutine initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,&
         kpq)

      use constants
      use parameters
      use diagmod
      use get_moment

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,i
      real(d), dimension(:), allocatable        :: ener,mtm,tmvec,&
                                                   osc_str
      real(d), dimension(:,:), allocatable      :: rvec

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(ener(davstates),rvec(ndim,davstates))
      allocate(mtm(ndim),tmvec(davstates),osc_str(davstates))

!-----------------------------------------------------------------------
! Read the Davidson state vectors from file
!-----------------------------------------------------------------------
      call readdavvc(davstates,ener,rvec,'i',ndim)

!-----------------------------------------------------------------------
! Calculate the vector F_J = < Psi_J | D | Psi_0 > (mtm)
!-----------------------------------------------------------------------
      if (ltdm_gs2i) call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:),1)

!-----------------------------------------------------------------------
! Contract the F-vector with the Davidson states to yield the
! transition moments between the ground state and the Davidson states
!-----------------------------------------------------------------------
      osc_str=0.0d0
      if (ltdm_gs2i) then
         do i=1,davstates
            tmvec(i)=tm(ndim,rvec(:,i),mtm(:))
            osc_str(i)=(2.0d0/3.0d0)*ener(i)*tmvec(i)**2
         enddo
      endif

      deallocate(mtm)
      
      return

    end subroutine initial_space_tdm

!#######################################################################

  end module adc2common
  
