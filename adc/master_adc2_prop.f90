  module adc2mod

    use channels

    contains

!#######################################################################

      subroutine master_adc2_prop()
        
        use constants
        use parameters
        use lancmod
        use fspace
        use misc
        use guessvecs

        implicit none

        integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
        integer                              :: i,ndim,ndims,ndimsf,&
                                                nout,ndimf,ndimd,&
                                                noutf,itmp
        integer*8                            :: noffd,noffdf
        real(d)                              :: time
        real(d), dimension(:), allocatable   :: ener,mtm,tmvec,osc_str
        real(d), dimension(:), allocatable   :: travec
        real(d)                              :: e_init,e0
        real(d), dimension(:,:), allocatable :: rvec
        real(d), dimension(:), allocatable   :: vec_init
        real*8, dimension(:), allocatable    :: mtmf

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic (if requested)
!-----------------------------------------------------------------------
        call run_mp2(e0)

!-----------------------------------------------------------------------  
! Calculate guess initial space vectors from an ADC(1) calculation if 
! requested.
!-----------------------------------------------------------------------  
  if (ladc1guess) call adc1_guessvecs

!-----------------------------------------------------------------------
! Determine the 1h1p and 2h2p subspaces
!-----------------------------------------------------------------------
        call get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf,&
             ndims,ndimsf)

!-----------------------------------------------------------------------
! Set the dipole matrix in the MO basis
!-----------------------------------------------------------------------
        call set_dpl

!-----------------------------------------------------------------------
! Block-Davidson diagonalisation in the initial space.
!-----------------------------------------------------------------------
        if (statenumber.gt.0) &
             call initial_space_diag(time,kpq,ndim,ndims,noffd)

!-----------------------------------------------------------------------
! Transition moments from the ground state to the Davidson states
!-----------------------------------------------------------------------
        if (statenumber.gt.0) &
             call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,&
             kpq)

!-----------------------------------------------------------------------
! Output the results of initial space calculation
!-----------------------------------------------------------------------
        if (statenumber.gt.0) then
           write(ilog,'(/,70a)') ('*',i=1,70)
           write(ilog,'(2x,a)') &
                'Initial space ADC(2)-s excitation energies'
           write(ilog,'(70a)') ('*',i=1,70)
           
           itmp=1+nBas**2*4*nOcc**2
           call table2(ndim,davstates,ener(1:davstates),&
                rvec(:,1:davstates),tmvec(1:davstates),&
                osc_str(1:davstates),kpq,itmp,'i')

           write(ilog,'(/,70a,/)') ('*',i=1,70)
        endif

!-----------------------------------------------------------------------
! Set the initial state vector and energy
!-----------------------------------------------------------------------
        allocate(vec_init(ndim))

        if (statenumber.gt.0) then
           vec_init(:)=rvec(:,statenumber)
           e_init=ener(statenumber)
        else
           e_init=0.0d0
        endif

!-----------------------------------------------------------------------
! Calculation of the final space states
!-----------------------------------------------------------------------
        call final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
           vec_init,mtmf,noffdf)

!-----------------------------------------------------------------------
! Calculate the transition moments and oscillator strengths between 
! the initial state and the final states
!
! Note that if we are NOT considering ionization, then we will also
! output the final Davidson state energies and configurations here
!-----------------------------------------------------------------------
        call final_space_tdm(ndimf,ndimsf,travec,e_init,mtmf,kpqf)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        if (allocated(ener)) deallocate(ener)
        if (allocated(rvec)) deallocate(rvec)
        if (allocated(vec_init)) deallocate(vec_init)        
        if (allocated(travec)) deallocate(travec)
        deallocate(kpq,kpqf,kpqd)

        return
        
      end subroutine master_adc2_prop

!#######################################################################

      subroutine get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,&
           noutf,ndims,ndimsf)
  
        use constants
        use parameters
        use select_fano

        implicit none

        integer                              :: ndim,ndimf,ndimd,&
                                                nout,noutf,ndims,&
                                                ndimsf
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

      subroutine set_dpl

        use parameters

        implicit none

!-----------------------------------------------------------------------
! Set the dipole matrix
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

        return

      end subroutine set_dpl

!#######################################################################

      subroutine initial_space_diag(time,kpq,ndim,ndims,noffd)
        
        use constants
        use parameters
        use fspace
        use davmod
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
        integer                                   :: ndim,ndims
        integer*8                                 :: noffd
        real(d)                                   :: time

!-----------------------------------------------------------------------
! Write the initial space ADC(2)-s Hamiltonian to disk
!-----------------------------------------------------------------------
        write(ilog,'(/,a)') 'Calculating the initial space ADC(2)-s &
             Hamiltonian matrix'

        call write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i') 

        call cpu_time(time)

        write(ilog,'(/,a,1x,F9.2,1x,a)') 'Time=',time," s"

!-----------------------------------------------------------------------
! Block-Davidson diagonalisation
!-----------------------------------------------------------------------
        call master_dav(ndim,noffd,'i')

        return

      end subroutine initial_space_diag

!#######################################################################

      subroutine initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,&
           kpq)

        use constants
        use parameters
        use davmod
        use get_moment

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
        integer                                   :: ndim,i
        real(d), dimension(:), allocatable        :: ener,mtm,tmvec,&
                                                     osc_str
        real(d), dimension(:,:), allocatable      :: rvec


        write(ilog,'(/,2x,a,x,a1,x,a)') 'Calculating the transition moments &
             between the ground state and the Davidson states in &
             the',tranmom2,'direction'

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

      subroutine final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
           vec_init,mtmf,noffdf)

        use constants
        use parameters
        use fspace
        use lancmod

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        integer*8                                 :: noffdf
        real(d), dimension(:), allocatable        :: travec,mtmf
        real(d), dimension(ndim)                  :: vec_init

        if (ldavfinal) then
           call davidson_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
                vec_init,mtmf,noffdf)
        else
           call lanczos_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
                vec_init,mtmf,noffdf)
        endif

        return

      end subroutine final_space_diag

!#######################################################################

      subroutine davidson_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,&
           travec,vec_init,mtmf,noffdf)

        use constants
        use parameters
        use fspace
        use get_matrix_dipole
        use davmod
        use guessvecs

        implicit none
        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        integer*8                                 :: noffdf
        real(d), dimension(:), allocatable        :: travec,mtmf
        real(d), dimension(ndim)                  :: vec_init

!-----------------------------------------------------------------------
! Allocate the travec array that will hold the contraction of the IS
! representation of the dipole operator with the initial state vector
!-----------------------------------------------------------------------
        allocate(travec(ndimf))

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the dipole
! operator and the initial state vector
!
! This will be used later in the calculation of transition dipole
! matrix elements between the initial and final states.
!-----------------------------------------------------------------------
        if (statenumber.eq.0) then
           call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:),0)
        else
           call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,&
                travec)
        endif

!-----------------------------------------------------------------------        
! If requested, determine the Davidson guess vectors by diagonalising 
! the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------        
        if (ladc1guess_f) call adc1_guessvecs_final

!-----------------------------------------------------------------------
! Write the final space ADC(2)-s Hamiltonian to disk
!-----------------------------------------------------------------------
        write(ilog,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
        if (lcvsfinal) then
           call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')           
        else
           call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
        endif

!-----------------------------------------------------------------------
! Block-Davidson diagonalisation in the final space
!-----------------------------------------------------------------------
        call master_dav(ndimf,noffdf,'f')

        return

      end subroutine davidson_final_space_diag

!#######################################################################

      subroutine lanczos_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
           vec_init,mtmf,noffdf)

        use constants
        use parameters
        use guessvecs
        use fspace
        use lancmod

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf,&
                                                     n
        integer*8                                 :: noffdf
        real(d), dimension(:), allocatable        :: travec,mtmf
        real(d), dimension(ndim)                  :: vec_init

!-----------------------------------------------------------------------        
! Acknowledging that we cannot use 2h2p unit vectors as initial
! Lanczos vectors, reduce lmain if it is greater than the number
! of final space 1h1p configurations
!-----------------------------------------------------------------------        
        if (lmain.gt.ndimsf) then
           write(ilog,'(/,2x,a,/)') 'Resetting the Lanczos block size &
                s.t. it is not greater than the dimension of the 1h1p &
                subspace'
           ! Number of Lanczos vectors requested
           n=lmain*ncycles
           ! Reset lmain
           lmain=ndimsf   
           ! Change ncycles s.t. the number of Lanczos vectors 
           ! generated does not change
           ncycles=n/lmain
        endif

!-----------------------------------------------------------------------        
! If requested, determine the initial vectors Lanczos vectors by 
! diagonalising the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------        
        if (lancguess.eq.2.or.lancguess.eq.4) call adc1_guessvecs_final

!-----------------------------------------------------------------------
! Allocate the travec array that will hold the contraction of the IS
! representation of the dipole operator with the initial state vector
!-----------------------------------------------------------------------
        allocate(travec(ndimf))

!-----------------------------------------------------------------------
! Determine the guess vectors for the band-Lanczos calculation
!
! Note that as part of this process, the transition moments between
! the ISs and the initial state are calculated in lanczos_guess_vecs
! and passed back in the travec array
!-----------------------------------------------------------------------
        call lanczos_guess_vecs(vec_init,ndim,ndimsf,&
             travec,ndimf,kpq,kpqf,mtmf)

!-----------------------------------------------------------------------
! Write the final space ADC(2)-s Hamiltonian matrix to file
!-----------------------------------------------------------------------
        write(ilog,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
        if (lcvsfinal) then
           call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')           
        else
           call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
        endif

!-----------------------------------------------------------------------
! Perform the band-Lanczos calculation
!-----------------------------------------------------------------------
        call master_lancdiag(ndimf,noffdf,'c')

        return

      end subroutine lanczos_final_space_diag

!#######################################################################

      subroutine lanczos_guess_vecs(vec_init,ndim,ndimsf,travec,ndimf,&
           kpq,kpqf,mtmf)

        use constants
        use parameters

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec
        real(d), dimension(:), allocatable        :: mtmf

!-----------------------------------------------------------------------
! Ionisation from the ground state
!-----------------------------------------------------------------------
        if (statenumber.eq.0) call guess_vecs_gs2ex(ndimf,ndimsf,mtmf,&
             kpqf)

!-----------------------------------------------------------------------
! Ionisation from an excited state
!-----------------------------------------------------------------------
        if (statenumber.gt.0) call guess_vecs_ex2ex(vec_init,ndim,&
             ndimsf,travec,ndimf,kpq,kpqf)

        return

      end subroutine lanczos_guess_vecs

!#######################################################################

      subroutine guess_vecs_gs2ex(ndimf,ndimsf,mtmf,kpqf)

        use constants
        use parameters
        use get_moment
        use fspace
        use guessvecs

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,&
                                                     iadc1,itmp,dim2,&
                                                     i,k1,k2
        integer, dimension(:), allocatable        :: indx1,indx2
        real(d), dimension(:), allocatable        :: mtmf
        real(d), dimension(ndimsf)                :: tmpvec
        real(d), dimension(ndimsf,ndimsf)         :: adc1vec

!-----------------------------------------------------------------------        
! Calculate the vector F_J = < Psi_J | D | Psi_0 >, where the Psi_J
! are the ISs spanning the final space
!
! N.B. this vector is saved to the mtmf array
!-----------------------------------------------------------------------
        allocate(mtmf(ndimf))
        call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:),0)

!-----------------------------------------------------------------------
! If requested, determine the block size based on the transition
! matrix elements between the initial state and the intermediate
! states (and possibly the ADC(1) eigenvectors)
!-----------------------------------------------------------------------
        if (ldynblock) call getblocksize(mtmf,ndimf,ndimsf)

!-----------------------------------------------------------------------
! From the values of the elements of mtmf (and/or the ADC(1)
! eigenvectors), determine which vectors will form the initial Lanczos
! vectors
!-----------------------------------------------------------------------
        if (lancguess.eq.1) then

           tmpvec=mtmf(1:ndimsf)
           call fill_stvc(ndimsf,tmpvec(1:ndimsf))

        else if (lancguess.eq.2) then           

           ! Read the ADC(1) eigenvectors from file
           call freeunit(iadc1)
           open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
                status='old')
           read(iadc1) itmp,adc1vec           
           close(iadc1)
           ! Contract the ADC(1) eigenvectors with the 1h1p part of
           ! the F-vector
           do i=1,ndimsf
              tmpvec(i)=dot_product(adc1vec(:,i),mtmf(1:ndimsf))
           enddo
           call fill_stvc(ndimsf,tmpvec(1:ndimsf))

        else if (lancguess.eq.3) then

           ! 1h1p ISs
           allocate(indx1(ndimsf))           
           call dsortindxa1('D',ndimsf,mtmf(1:ndimsf)**2,indx1(:))

           ! 2h2p ISs
           dim2=ndimf-ndimsf
           allocate(indx2(dim2))
           call dsortindxa1('D',dim2,mtmf(ndimsf+1:ndimf)**2,indx2(:))

           ! Fill in the stvc_mxc array
           allocate(stvc_mxc(3*lmain))
           do i=1,lmain

              k1=indx1(i)
              k2=ndimsf+indx2(i)

              ! 1h1p IS plus or minus the 2h2p IS (chosen st the
              ! resulting vector has the greates TDM with the initial
              ! state)
              if (mtmf(k1).gt.0.and.mtmf(k2).gt.0) then
                 stvc_mxc(i*3-2)=1
              else
                 stvc_mxc(i*3-2)=-1
              endif

              ! Index of the 1h1p IS
              stvc_mxc(i*3-1)=k1

              ! Index of the 2h2p IS
              stvc_mxc(i*3)=k2
              
           enddo

           deallocate(indx1,indx2)

        else if (lancguess.eq.4) then
           
           ! Read the ADC(1) eigenvectors from file
           call freeunit(iadc1)
           open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
                status='old')
           read(iadc1) itmp,adc1vec           
           close(iadc1)
           ! Contract the ADC(1) eigenvectors with the 1h1p part of
           ! the F-vector
           do i=1,ndimsf
              tmpvec(i)=dot_product(adc1vec(:,i),mtmf(1:ndimsf))
           enddo
           
           ! ADC(1) eigenvectors
           allocate(indx1(ndimsf))
           call dsortindxa1('D',ndimsf,tmpvec(1:ndimsf)**2,indx1(:))
           
           ! 2h2p ISs
           dim2=ndimf-ndimsf
           allocate(indx2(dim2))
           call dsortindxa1('D',dim2,mtmf(ndimsf+1:ndimf)**2,indx2(:))

           ! Fill in the stvc_mxc array
           allocate(stvc_mxc(3*lmain))

           do i=1,lmain

              k1=indx1(i)
              k2=ndimsf+indx2(i)

              ! ADC(1) eigenvector plus or minus the 2h2p IS (chosen st
              ! the resulting vector has the greates TDM with the
              ! initial state)
              if (tmpvec(k1).gt.0.and.mtmf(k2).gt.0) then
                 stvc_mxc(i*3-2)=1
              else
                 stvc_mxc(i*3-2)=-1
              endif

              ! Index of the ADC(1) eigenvector
              stvc_mxc(i*3-1)=k1

              ! Index of the 2h2p IS
              stvc_mxc(i*3)=k2
              
           enddo

           deallocate(indx1,indx2)

        endif

        return

      end subroutine guess_vecs_gs2ex

!#######################################################################

      subroutine guess_vecs_ex2ex(vec_init,ndim,ndimsf,travec,ndimf,&
           kpq,kpqf)

        use constants
        use parameters
        use misc
        use get_matrix_dipole
        use fspace
        use guessvecs

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: i,ndim,ndimf,&
                                                     ndimsf,iadc1,&
                                                     itmp,k1,k2,dim2
        integer, dimension(:), allocatable        :: indx_tra,indx1,&
                                                     indx2
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec
        real(d), dimension(ndimsf)                :: tmpvec
        real(d), dimension(ndimsf,ndimsf)         :: adc1vec

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the dipole
! operator and the initial state vector
!-----------------------------------------------------------------------        
        call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,&
             travec)

!-----------------------------------------------------------------------
! If requested, determine the block size based on the transition
! matrix elements between the initial state and the intermediate
! states (and possibly the ADC(1) eigenvectors)
!-----------------------------------------------------------------------
        if (ldynblock) call getblocksize(travec,ndimf,ndimsf)

!-----------------------------------------------------------------------
! From the values of the elements of travec (and/or the ADC(1)
! eigenvectors), determine which vectors will form the initial Lanczos
! vectors
!-----------------------------------------------------------------------
        if (lancguess.eq.1) then

           tmpvec=travec(1:ndimsf)
           call fill_stvc(ndimsf,travec(1:ndimsf))

        else if (lancguess.eq.2) then           

           ! Read the ADC(1) eigenvectors from file
           call freeunit(iadc1)
           open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
                status='old')
           read(iadc1) itmp,adc1vec           
           close(iadc1)
           ! Contract the ADC(1) eigenvectors with the product of the 
           ! 1h1p part of travec
           do i=1,ndimsf
              tmpvec(i)=dot_product(adc1vec(:,i),travec(1:ndimsf))
           enddo
           call fill_stvc(ndimsf,travec(1:ndimsf))

        else if (lancguess.eq.3) then

           ! 1h1p ISs
           allocate(indx1(ndimsf))           
           call dsortindxa1('D',ndimsf,travec(1:ndimsf)**2,indx1(:))

           ! 2h2p ISs
           dim2=ndimf-ndimsf
           allocate(indx2(dim2))
           call dsortindxa1('D',dim2,travec(ndimsf+1:ndimf)**2,indx2(:))

           ! Fill in the stvc_mxc array
           allocate(stvc_mxc(3*lmain))
           do i=1,lmain

              k1=indx1(i)
              k2=ndimsf+indx2(i)

              ! 1h1p IS plus or minus the 2h2p IS (chosen st the
              ! resulting vector has the greates TDM with the initial
              ! state)
              if (travec(k1).gt.0.and.travec(k2).gt.0) then
                 stvc_mxc(i*3-2)=1
              else
                 stvc_mxc(i*3-2)=-1
              endif

              ! Index of the 1h1p IS
              stvc_mxc(i*3-1)=k1

              ! Index of the 2h2p IS
              stvc_mxc(i*3)=k2
              
           enddo

           deallocate(indx1,indx2)

        else if (lancguess.eq.4) then

           ! Read the ADC(1) eigenvectors from file
           call freeunit(iadc1)
           open(iadc1,file='SCRATCH/adc1_vecs',form='unformatted',&
                status='old')
           read(iadc1) itmp,adc1vec           
           close(iadc1)
           ! Contract the ADC(1) eigenvectors with the 1h1p part of
           ! the F-vector
           do i=1,ndimsf
              tmpvec(i)=dot_product(adc1vec(:,i),travec(1:ndimsf))
           enddo
           
           ! ADC(1) eigenvectors
           allocate(indx1(ndimsf))           
           call dsortindxa1('D',ndimsf,tmpvec(1:ndimsf)**2,indx1(:))
           
           ! 2h2p ISs
           dim2=ndimf-ndimsf
           allocate(indx2(dim2))
           call dsortindxa1('D',dim2,travec(ndimsf+1:ndimf)**2,indx2(:))

           ! Fill in the stvc_mxc array
           allocate(stvc_mxc(3*lmain))

           do i=1,lmain

              k1=indx1(i)
              k2=ndimsf+indx2(i)

              ! ADC(1) eigenvector plus or minus the 2h2p IS (chosen st
              ! the resulting vector has the greates TDM with the
              ! initial state)
              if (tmpvec(k1).gt.0.and.travec(k2).gt.0) then
                 stvc_mxc(i*3-2)=1
              else
                 stvc_mxc(i*3-2)=-1
              endif

              ! Index of the ADC(1) eigenvector
              stvc_mxc(i*3-1)=k1

              ! Index of the 2h2p IS
              stvc_mxc(i*3)=k2
              
           enddo

           deallocate(indx1,indx2)

        endif

        return

      end subroutine guess_vecs_ex2ex

!#######################################################################

      subroutine final_space_tdm(ndimf,ndimsf,travec,e_init,mtmf,kpqf)

        use constants
        use parameters

        implicit none
 
        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf
        real(d), dimension(ndimf)                 :: travec,mtmf
        real(d)                                   :: e_init

        if (ldavfinal) then
           call tdm_davstates_final(ndimf,ndimsf,travec,e_init,&
                mtmf,kpqf)
        else
           call tdm_lancstates(ndimf,ndimsf,travec,e_init,mtmf)
        endif

        return

      end subroutine final_space_tdm

!#######################################################################

      subroutine tdm_lancstates(ndimf,ndimsf,travec,e_init,mtmf)
        
        use constants
        use parameters

        implicit none

        integer                   :: ndimf,ndimsf
        real(d), dimension(ndimf) :: travec,mtmf
        real(d)                   :: e_init

!-----------------------------------------------------------------------
! Transition moments from the ground state
!-----------------------------------------------------------------------
        if (statenumber.eq.0) call tdm_gs2lanc(ndimf,ndimsf,mtmf,e_init)

!-----------------------------------------------------------------------
! Transition moments from an excited state
!-----------------------------------------------------------------------
        if (statenumber.gt.0) call tdm_ex2lanc(ndimf,ndimsf,travec,&
             e_init)

        return

      end subroutine tdm_lancstates

!#######################################################################

      subroutine tdm_gs2lanc(ndimf,ndimsf,mtmf,e_init)

        use constants
        use parameters
        use fspace

        implicit none

        integer                            :: ndimf,ndimsf,nstates,i
        real(d)                            :: e_init
        real(d), dimension(ndimf)          :: mtmf
        real(d), dimension(:), allocatable :: tmvecf,enerf,excit,&
                                              osc_strf

!-----------------------------------------------------------------------
! Calculate the transition moments between the ground state and the
! Lanczos states
!
! N.B. nstates is determined in get_tranmom_1
!-----------------------------------------------------------------------
        allocate(enerf(lancstates),tmvecf(lancstates))
        tmvecf=0.0d0
        enerf=0.0d0

        call get_tranmom_1(ndimf,lancstates,lancname,mtmf(:),nstates,enerf(:),&
             tmvecf(:),ndimsf)

!-----------------------------------------------------------------------
! Calculate and output the oscillator strengths
!-----------------------------------------------------------------------
        allocate(osc_strf(nstates),excit(nstates))

        do i=1,nstates
           excit(i)=enerf(i)-e_init
        end do

        osc_strf=0.0d0
        do i=1,nstates
           osc_strf(i)=(2.0d0/3.0d0)*excit(i)*tmvecf(i)**2
        enddo

        call get_sigma(nstates,excit(1:nstates),os2cs*osc_strf(:))

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(enerf,tmvecf,osc_strf,excit)
        
        return

      end subroutine tdm_gs2lanc

!#######################################################################

      subroutine tdm_ex2lanc(ndimf,ndimsf,travec,e_init)

        use constants
        use parameters
        use fspace
        use fspace2
        
        implicit none

        integer                            :: ndimf,ndimsf,nstates,&
                                              i
        real(d), dimension(ndimf)          :: travec
        real(d)                            :: e_init
        real(d), dimension(:), allocatable :: tmvecf,enerf,excit,&
                                              osc_strf

!-----------------------------------------------------------------------
! Calculate the transition moments between the initial state and the
! Lanczos states
!
! N.B. nstates is determined in get_tranmom_3
!-----------------------------------------------------------------------
        allocate(enerf(lancstates),tmvecf(lancstates))
        tmvecf=0.0d0
        enerf=0.0d0

        call get_tranmom_3(ndimf,lancstates,lancname,travec(:),&
             nstates,enerf(:),tmvecf(:),ndimsf)

!-----------------------------------------------------------------------
! Calculate and output the oscillator strengths
!-----------------------------------------------------------------------
        allocate(osc_strf(nstates),excit(nstates))

        do i=1,nstates
           excit(i)=enerf(i)-e_init           
        enddo

        osc_strf=0.0d0
        do i=1,nstates
           osc_strf(i)=(2.0d0/3.0d0)*excit(i)*tmvecf(i)**2
        enddo

        call get_sigma(nstates,excit(1:nstates),os2cs*osc_strf(:))

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(enerf,tmvecf,osc_strf,excit)

        return

      end subroutine tdm_ex2lanc

!#######################################################################

      subroutine tdm_davstates_final(ndimf,ndimsf,travec,e_init,mtmf,&
           kpqf)
        
        use constants
        use parameters
        use iomod, only: freeunit
        use misc, only: dsortindxa1,table2
        use davmod, only: readdavvc
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,i,&
                                                     itmp,iout
        real(d)                                   :: e_init,ftmp
        real(d), dimension(ndimf)                 :: travec,mtmf
        real(d), dimension(davstates_f)           :: ener,tdm,osc_str
        real(d), dimension(:,:), allocatable      :: rvec

!-----------------------------------------------------------------------
! Read the final space Davidson states from file
!-----------------------------------------------------------------------
        allocate(rvec(ndimf,davstates_f))
        call readdavvc(davstates_f,ener,rvec,'f',ndimf)

!-----------------------------------------------------------------------
! Calculate the TDMs from the initial state
!-----------------------------------------------------------------------
        do i=1,davstates_f
           !Form the scalar product of the ith final state vector
           ! and mtmf (if statenumber=0) or travec (if statenumber>0)
           if (statenumber.eq.0) then
              tdm(i)=dot_product(rvec(:,i),mtmf)
           else
              tdm(i)=dot_product(rvec(:,i),travec)
           endif
           osc_str(i)=(2.0d0/3.0d0)*ener(i)*tdm(i)**2
       enddo

!-----------------------------------------------------------------------
! Output the final space energies, TDMs and configurations
!-----------------------------------------------------------------------
       write(ilog,'(/,70a)') ('*',i=1,70)
       if (lcvsfinal) then
          write(ilog,'(2x,a)') &
               'Final space CVS-ADC(2)-s excitation energies'
       else
          write(ilog,'(2x,a)') &
               'Final space ADC(2)-s excitation energies'
       endif
       write(ilog,'(70a)') ('*',i=1,70)
       
       itmp=1+nBas**2*4*nOcc**2
       call table2(ndimf,davstates_f,ener(1:davstates_f),&
            rvec(:,1:davstates_f),tdm(1:davstates_f),&
            osc_str(1:davstates_f),kpqf,itmp,'f')

       deallocate(rvec)

!-----------------------------------------------------------------------
! Output the excitation energies (relative to the initial state) and 
! oscillator strengths for plotting purposes
!-----------------------------------------------------------------------  
       call freeunit(iout)
       open(iout,file='osc.dat',form='formatted',status='unknown')

       do i=1,davstates_f
          if (abs(osc_str(i)).lt.0.00001d0) cycle
          write(iout,'(2(F11.5,2x))') (ener(i)-e_init)*eh2ev,osc_str(i)
       enddo

       close(iout)
       
       return

     end subroutine tdm_davstates_final

!#######################################################################

      subroutine run_mp2(e0)
        
        use constants
        use parameters
!        use diagnostics        
        use adc_ph, only: mp2

        implicit none

        integer           :: i
        real(d)           :: e0,d2
        character(len=60) :: atmp

!-----------------------------------------------------------------------
! Calculation of the MP2 correlation energy
!-----------------------------------------------------------------------
        call MP2(E_MP2)

!-----------------------------------------------------------------------
! Calculation of the D2 diagnostic
!-----------------------------------------------------------------------
!        if (ld2) call mp2_d2(d2)

!-----------------------------------------------------------------------
! Output results
!-----------------------------------------------------------------------
        e0=Ehf+E_MP2
 
        write(ilog,'(/,2x,90a)') ('*',i=1,90)

        write(ilog,'(2x,1a)') '*'
        atmp='* HF energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),Ehf

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 correlation energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),E_MP2

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),e0

        if (ld2) then
           write(ilog,'(2x,1a)') '*'
           atmp='* D2 diagnostic:'
           write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),d2
        endif

        write(ilog,'(2x,1a)') '*'
        write(ilog,'(2x,90a,/)') ('*',i=1,90)

        return

      end subroutine run_mp2

!#######################################################################

    end module adc2mod
