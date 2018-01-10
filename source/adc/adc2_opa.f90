!#######################################################################
! Calculation of one photon ADC(2) spectra. Either:
!
! (1) Bound state-bound state absorption spectra, or;
!
! (2) Spectral moments for use in Stieltjes imaging calculations of
!     photoionisation spectra.
!#######################################################################

  module adc2specmod

    use channels

    contains

!#######################################################################

      subroutine adc2_spec(gam)
        
        use constants
        use parameters
        use adc2common
        use fspace
        use misc
        use guessvecs
        use mp2
        use targetmatching
        use nto
        
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
        real(d), dimension(:,:), allocatable :: rvec,travec2
        real(d), dimension(:), allocatable   :: vec_init
        real*8, dimension(:), allocatable    :: mtmf
        type(gam_structure)                  :: gam

!        ! TEST
!        call orbrlx2
!        STOP
!        ! TEST

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!-----------------------------------------------------------------------
        call mp2_master(e0)

!-----------------------------------------------------------------------  
! Calculate guess initial space vectors from an ADC(1) calculation if 
! requested
!-----------------------------------------------------------------------  
        if (ladc1guess) call adc1_guessvecs

!-----------------------------------------------------------------------
! Determine the 1h1p and 2h2p subspaces
!-----------------------------------------------------------------------
        call get_subspaces(kpq,kpqf,kpqd,ndim,ndimf,ndimd,nout,noutf,&
             ndims,ndimsf)

!-----------------------------------------------------------------------
! Set MO representation of the dipole operator
!-----------------------------------------------------------------------
        call set_dpl

!-----------------------------------------------------------------------
! Diagonalisation in the initial space
!-----------------------------------------------------------------------
        if (statenumber.gt.0) &
             call initial_space_diag(time,kpq,ndim,ndims,noffd)

!-----------------------------------------------------------------------
! If requested, calculate the dipole moments for the initial states
!-----------------------------------------------------------------------
        if (ldipole.and.statenumber.gt.0) &
             call initial_space_dipole(ndim,ndims,kpq)

!-----------------------------------------------------------------------
! Transition moments from the ground state to the Davidson states
!-----------------------------------------------------------------------
        if (statenumber.gt.0) call initial_space_tdm(ener,rvec,ndim,&
             mtm,tmvec,osc_str,kpq)

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
! Selection of the initial state by matching to a target CI state
!-----------------------------------------------------------------------
        if (ltarg.and.statenumber.ne.0) call target_master(kpq,ndim,gam)

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
           vec_init,mtmf,noffd,noffdf,rvec,travec2)

!-----------------------------------------------------------------------
! If requested, calculate the dipole moments for the final states
! N.B. This is only meaningful if the final states are Davidson states
!-----------------------------------------------------------------------
        if (ldipole.and.ldiagfinal) then
           call final_space_dipole(ndimf,ndimsf,kpqf)
        endif
        
!-----------------------------------------------------------------------
! Calculate the transition moments and oscillator strengths between 
! the initial state and the final states
!
! Note that if we are NOT considering ionization, then we will also
! output the final Davidson state energies and configurations here 
!-----------------------------------------------------------------------
        call final_space_tdm(ndimf,ndimsf,travec,e_init,mtmf,kpqf,&
             travec2,ndim)

!-----------------------------------------------------------------------
! If requested, calculate NTOs
!-----------------------------------------------------------------------
        if (ldiagfinal.and.lnto) call adc2_nto(ndimf,kpqf)
        
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        if (allocated(ener)) deallocate(ener)
        if (allocated(rvec)) deallocate(rvec)
        if (allocated(vec_init)) deallocate(vec_init)
        if (allocated(travec)) deallocate(travec)
        if (allocated(dipmom)) deallocate(dipmom)
        if (allocated(dipmom_f)) deallocate(dipmom_f)
        if (allocated(travec2)) deallocate(travec2)
        if (allocated(dpl_all)) deallocate(dpl_all)
        if (allocated(travec_ic)) deallocate(travec_ic)
        if (allocated(travec_iv)) deallocate(travec_iv)
        if (allocated(travec_fc)) deallocate(travec_fc)
        if (allocated(travec_fv)) deallocate(travec_fv)
        if (allocated(tdmgsf)) deallocate(tdmgsf)
        if (allocated(edavf)) deallocate(edavf)
        deallocate(kpq,kpqf,kpqd)

        return
        
      end subroutine adc2_spec
      
!#######################################################################

      subroutine final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
           vec_init,mtmf,noffd,noffdf,rvec,travec2)

        use constants
        use parameters
        use fspace
        use adc2common
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        integer*8                                 :: noffd,noffdf
        real(d), dimension(:), allocatable        :: travec,mtmf
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndim,davstates)        :: rvec
        real(d), dimension(:,:), allocatable      :: travec2

        if (ldiagfinal) then
           call davidson_final_space_diag(ndim,ndimf,ndimsf,kpq,&
                kpqf,travec,vec_init,mtmf,noffdf,rvec,travec2)
        else
           call lanczos_final_space_diag(ndim,ndimf,ndimsf,kpq,&
                kpqf,travec,vec_init,mtmf,noffdf,rvec,travec2)
        endif
        
        return

      end subroutine final_space_diag

!#######################################################################
      
      subroutine lanczos_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,&
           travec,vec_init,mtmf,noffdf,rvec,travec2)

        use constants
        use parameters
        use guessvecs
        use fspace
        use block_lanczos
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf,&
                                                     n
        integer*8                                 :: noffdf
        real(d), dimension(:), allocatable        :: travec,mtmf
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndim,davstates)        :: rvec
        real(d), dimension(:,:), allocatable      :: travec2

!-----------------------------------------------------------------------        
! Acknowledging that we cannot use 2h2p unit vectors as initial
! Lanczos vectors for ADC(2)-s, reduce lmain if it is greater than the 
! number of final space 1h1p configurations
!-----------------------------------------------------------------------        
        if (method_f.eq.2.and.lmain.gt.ndimsf) then
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
! If requested, determine the initial Lanczos vectors by diagonalising
! the ADC(1) Hamiltonian matrix
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
             travec,ndimf,kpq,kpqf,mtmf,rvec,travec2)

!-----------------------------------------------------------------------
! Write the final space ADC(2) Hamiltonian matrix to file
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
! Perform the block-Lanczos calculation
!-----------------------------------------------------------------------
        call lancdiag_block(ndimf,noffdf,'c')

        return

      end subroutine lanczos_final_space_diag

!#######################################################################

      subroutine lanczos_guess_vecs(vec_init,ndim,ndimsf,travec,ndimf,&
           kpq,kpqf,mtmf,rvec,travec2)

        use constants
        use parameters

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec
        real(d), dimension(:), allocatable        :: mtmf
        real(d), dimension(ndim,davstates)        :: rvec
        real(d), dimension(:,:), allocatable      :: travec2

        if (statenumber.eq.0) then
           ! Ionisation from the ground state
           call guess_vecs_gs2ex(ndimf,ndimsf,mtmf,kpqf)
        else
           ! Ionisation from an excited state
           call guess_vecs_ex2ex(vec_init,ndim,ndimsf,travec,ndimf,&
                kpq,kpqf)
        endif

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
                                                     i,k1,k2,upper
        integer, dimension(:), allocatable        :: indx1,indx2
        real(d), dimension(:), allocatable        :: mtmf
        real(d), dimension(:), allocatable        :: tmpvec
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
!
! Note that we can only have lancguess=3 or 4 for ADC(2)-s
!-----------------------------------------------------------------------
        allocate(tmpvec(ndimf))

        if (lancguess.eq.1) then
           if (method.eq.2) then
              upper=ndimsf
           else if (method.eq.3) then
              upper=ndimf
           endif
           tmpvec=mtmf(1:upper)
           call fill_stvc(upper,tmpvec(1:upper))

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
        
        deallocate(tmpvec)

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
                                                     itmp,k1,k2,&
                                                     dim2,upper
        integer, dimension(:), allocatable        :: indx_tra,indx1,&
                                                     indx2
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec
        real(d), dimension(:), allocatable        :: tmpvec
        real(d), dimension(ndimsf,ndimsf)         :: adc1vec

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the
! shifted operator and the initial state vector 
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
!
! Note that we can only have lancguess=3 or 4 for ADC(2)-s
!-----------------------------------------------------------------------
        allocate(tmpvec(ndimf))

        if (lancguess.eq.1) then
           if (method_f.eq.2) then
              upper=ndimsf
           else if (method_f.eq.3) then
              upper=ndimf
           endif
           tmpvec=travec(1:upper)
           call fill_stvc(upper,travec(1:upper))

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

        deallocate(tmpvec)

        return

      end subroutine guess_vecs_ex2ex

!#######################################################################

      subroutine final_space_tdm(ndimf,ndimsf,travec,e_init,mtmf,kpqf,&
           travec2,ndim)

        use constants
        use parameters

        implicit none
 
        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,ndim
        real(d), dimension(ndimf)                 :: travec,mtmf
        real(d)                                   :: e_init
        real(d), dimension(ndimf,3*(davstates+1)) :: travec2

        if (ldiagfinal) then
           ! Davidson states
           call tdm_davstates_final(ndimf,ndimsf,travec,e_init,&
                mtmf,kpqf)
        else
           ! Lanczos pseudospectrum
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

        call get_tranmom_1(ndimf,lancstates,lancname,mtmf(:),nstates,&
             enerf(:),tmvecf(:),ndimsf)

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

        call get_sigma(nstates,excit(1:nstates),osc_strf(:))

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

        call get_sigma(nstates,excit(1:nstates),osc_strf(:))

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
        use diagmod, only: readdavvc
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,i,&
                                                     itmp,iout,count,k
        real(d)                                   :: e_init,ftmp,lb,&
                                                     ub
        real(d), dimension(ndimf)                 :: travec,mtmf
        real(d), dimension(davstates_f)           :: ener,tdm,osc_str
        real(d), dimension(:,:), allocatable      :: rvec
        real(d), parameter                        :: tol=0.00001d0
        character(len=400)                        :: atmp

!-----------------------------------------------------------------------
! Read the final space Davidson states from file
!-----------------------------------------------------------------------
        allocate(rvec(ndimf,davstates_f))
        call readdavvc(davstates_f,ener,rvec,'f',ndimf)

!-----------------------------------------------------------------------
! Calculate the TDMs from the initial state
!-----------------------------------------------------------------------
        do i=1,davstates_f
           ! Form the scalar product of the ith final state vector
           ! and mtmf (if statenumber=0) or travec (if statenumber>0)
           if (statenumber.eq.0) then
              tdm(i)=dot_product(rvec(:,i),mtmf)
           else
              tdm(i)=dot_product(rvec(:,i),travec)
           endif
           osc_str(i)=(2.0d0/3.0d0)*(ener(i)-e_init)*tdm(i)**2
        enddo

!-----------------------------------------------------------------------
! Output the final space energies, TDMs and configurations
!-----------------------------------------------------------------------
       write(ilog,'(/,70a)') ('*',i=1,70)
       if (lcvsfinal) then
          if (method_f.eq.2) then
             write(ilog,'(2x,a)') &
                  'Final space CVS-ADC(2)-s excitation energies'
          else if (method_f.eq.3) then
             write(ilog,'(2x,a)') &
                  'Final space CVS-ADC(2)-x excitation energies'
          endif
       else
          if (method_f.eq.2) then
             write(ilog,'(2x,a)') &
                  'Final space ADC(2)-s excitation energies'
          else if (method_f.eq.3) then
             write(ilog,'(2x,a)') &
                  'Final space ADC(2)-x excitation energies'
          endif               
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
       ! Data file
       call freeunit(iout)
       open(iout,file='osc.dat',form='formatted',status='unknown')

       do i=1,davstates_f
          if (abs(osc_str(i)).lt.tol) cycle
          write(iout,'(2(F11.5,2x))') (ener(i)-e_init)*eh2ev,osc_str(i)
       enddo

       close(iout)

       ! gnuplot file
       open(iout,file='spec.gnu',form='formatted',status='unknown')

       write(iout,'(a,/)') 'fwhm=1.0'
       write(iout,'(a,/)') 'sigsq=(fwhm/2.35482)**2'
       
       count=0
       do i=1,davstates_f
          if (abs(osc_str(i)).lt.tol) cycle
          count=count+1
          if (count.lt.10) then
             write(iout,'(a2,i1,a1,F11.5)') 'e'//tranmom2,count,'=',(ener(i)-e_init)*eh2ev
             write(iout,'(a2,i1,a1,F11.5)') 'f'//tranmom2,count,'=',osc_str(i)
          else if (count.lt.100) then
             write(iout,'(a2,i2,a1,F11.5)') 'e'//tranmom2,count,'=',(ener(i)-e_init)*eh2ev
             write(iout,'(a2,i2,a1,F11.5)') 'f'//tranmom2,count,'=',osc_str(i)
          else
             write(iout,'(a2,i3,a1,F11.5)') 'e'//tranmom2,count,'=',(ener(i)-e_init)*eh2ev
             write(iout,'(a2,i3,a1,F11.5)') 'f'//tranmom2,count,'=',osc_str(i)
          endif
       enddo
       
       do i=1,count
          if (i.lt.10) then
             write(iout,'(a2,i1,a6,i1,a11,i1,a15)') &
                  'g'//tranmom2,i,'(x)=f'//tranmom2,i,'*exp(-(x-e'//tranmom2,i,')**2/(2*sigsq))'
          else if (i.lt.100) then
             write(iout,'(a2,i2,a6,i2,a11,i2,a15)') &
                  'g'//tranmom2,i,'(x)=f'//tranmom2,i,'*exp(-(x-e'//tranmom2,i,')**2/(2*sigsq))'
          else
             write(iout,'(a2,i3,a6,i3,a11,i3,a15)') &
                  'g'//tranmom2,i,'(x)=f'//tranmom2,i,'*exp(-(x-e'//tranmom2,i,')**2/(2*sigsq))'
          endif
       enddo

       atmp='f'//tranmom2//'(x)='
       k=7
       do i=1,count
          if (i.lt.10) then
             write(atmp(k:k+6),'(a3,i1,a3)') '+g'//tranmom2,i,'(x)'
             k=k+7
          else if (i.lt.100) then
             write(atmp(k:k+7),'(a3,i2,a3)') '+g'//tranmom2,i,'(x)'
             k=k+8
          else
             write(atmp(k:k+8),'(a3,i3,a3)') '+g'//tranmom2,i,'(x)'
             k=k+9
          endif
       enddo

       write(iout,'(a)') trim(atmp)

       lb=0.95*((ener(1)-e_init)*eh2ev)
       ub=1.05*(ener(davstates_f)-e_init)*eh2ev

       write(iout,'(a12,F7.2,a1,F7.2,a1)') 'set xrange [',lb,':',ub,']'
       
       write(iout,'(a)') 'set size ratio 0.4'
       write(iout,'(a)') 'set samples 5000'       
       write(iout,'(a)') 'plot f'//tranmom2//'(x) lt -1 notitle'
       write(iout,'(a)') 'replot ''osc.dat'' with impulses lt -1 notitle'
       write(iout,'(a)') 'pause -1'

       close(iout)

       return

     end subroutine tdm_davstates_final

!#######################################################################

     subroutine orbrlx2

       use constants
       use parameters
       use vpqrsmod
       
       implicit none

       integer               :: ncoreorb,cindx,vindx,c,a,i
       real(d)               :: term0,rterm,sterm,pterm,delta_ai,&
                                ftmp,tot,u1p1h
       real(d), dimension(2) :: term1

       ! Temporary hard-wiring of parameters
       !
       ! No. core orbitals
       ncoreorb=2
       ! C 1s
       cindx=2
       ! 3s         
!       vindx=9
       ! pi*
!       vindx=13
       ! 3p_1
!       vindx=10
       ! 3p_2
       vindx=11

!-----------------------------------------------------------------------
! Zeroth-order term
!-----------------------------------------------------------------------
       term0=e(vindx)-e(cindx)

!-----------------------------------------------------------------------
! First-order term
!-----------------------------------------------------------------------
       term1(1)=-vpqrs(cindx,cindx,vindx,vindx)
       term1(2)=2.0d0*vpqrs(cindx,vindx,vindx,cindx)

!-----------------------------------------------------------------------
! 1h1p term
!-----------------------------------------------------------------------
       u1p1h=0.0d0
       do a=nocc+1,nbas
          if (a.eq.vindx) cycle
          delta_ai=1.0d0/(e(vindx)-e(a))
          ftmp=(vpqrs(vindx,a,cindx,cindx)-2.0d0*vpqrs(vindx,cindx,cindx,a))**2
          u1p1h=u1p1h+delta_ai*ftmp
       enddo

!-----------------------------------------------------------------------
! Relaxation energy
!-----------------------------------------------------------------------
!       ! Paper II
!       rterm=0.0d0
!       do c=1,ncoreorb
!          if (c.eq.cindx) cycle
!          do i=1,nocc
!             do a=nocc+1,nbas
!                delta_ai=1.0d0/(e(a)-e(i))
!                rterm=rterm-2.0d0*delta_ai*vpqrs(c,cindx,i,a)**2
!                if (a.eq.vindx) rterm=rterm+delta_ai*vpqrs(c,cindx,i,a)**2
!             enddo
!          enddo
!       enddo

       ! Paper I
       rterm=0.0d0
       do i=1,nocc
          do a=nocc+1,nbas
             if (a.eq.vindx) cycle
             delta_ai=1.0d0/(e(a)-e(i))
             ftmp=vpqrs(cindx,cindx,i,a)**2
             ftmp=ftmp-vpqrs(cindx,cindx,i,a)*vpqrs(cindx,a,i,cindx)
             ftmp=ftmp+vpqrs(cindx,a,i,cindx)**2
             rterm=rterm-2.0d0*delta_ai*ftmp
          enddo
       enddo
       do i=1,nocc
          delta_ai=1.0d0/(e(vindx)-e(i))
          ftmp=(vpqrs(cindx,cindx,i,vindx)+vpqrs(cindx,vindx,i,cindx))**2
          rterm=rterm-delta_ai*ftmp
       enddo

!-----------------------------------------------------------------------
! Screening energy
!-----------------------------------------------------------------------
!       ! Paper II
!       sterm=0.0d0
!       do i=1,nocc
!          do a=nocc+1,nbas
!             delta_ai=1.0d0/(e(a)-e(i))
!             sterm=sterm &
!                  +4.0d0*delta_ai*vpqrs(cindx,cindx,i,a)*vpqrs(vindx,vindx,i,a)
!             if (a.eq.vindx) sterm=sterm &
!                  -2.0d0*delta_ai*vpqrs(cindx,cindx,i,a)*vpqrs(vindx,vindx,i,a)
!          enddo
!       enddo

       ! Paper I
       sterm=0.0d0
       do i=1,nocc
          do a=nocc+1,nbas
             if (a.eq.vindx) cycle
             delta_ai=1.0d0/(e(a)-e(i))
             ftmp=2.0d0*vpqrs(cindx,cindx,i,a)*vpqrs(vindx,vindx,i,a)
             ftmp=ftmp-vpqrs(cindx,cindx,i,a)*vpqrs(vindx,a,i,vindx)
             ftmp=ftmp-vpqrs(cindx,a,i,cindx)*vpqrs(vindx,vindx,i,a)
             ftmp=ftmp+2.0d0*vpqrs(cindx,a,i,cindx)*vpqrs(vindx,a,i,vindx)
             sterm=sterm+2.0d0*delta_ai*ftmp
          enddo
       enddo
       do i=1,nocc
          delta_ai=1.0d0/(e(vindx)-e(i))
          ftmp=(vpqrs(cindx,cindx,i,vindx)+vpqrs(cindx,vindx,i,cindx))*vpqrs(vindx,vindx,i,vindx)
          sterm=sterm+2.0d0*delta_ai*ftmp
       enddo

!-----------------------------------------------------------------------
! Polarisation energy
!-----------------------------------------------------------------------
!       ! Paper II
!       pterm=0.0d0
!       do i=1,nocc
!          do a=nocc+1,nbas
!             delta_ai=1.0d0/(e(a)-e(i))
!             pterm=pterm-2.0d0*delta_ai*vpqrs(vindx,vindx,i,a)**2
!             if (a.eq.vindx) pterm=pterm+delta_ai*vpqrs(vindx,vindx,i,a)**2
!          enddo
!       enddo

       ! Paper I
       pterm=0.0d0
       do i=1,nocc
          do a=nocc+1,nbas
             if (a.eq.vindx) cycle
             delta_ai=1.0d0/(e(a)-e(i))
             ftmp=vpqrs(vindx,vindx,i,a)**2
             ftmp=ftmp-vpqrs(vindx,vindx,i,a)*vpqrs(vindx,a,i,vindx)
             ftmp=ftmp+vpqrs(vindx,a,i,vindx)**2
             pterm=pterm-2.0d0*delta_ai*ftmp
          enddo
       enddo

       do i=1,nocc
          delta_ai=1.0d0/(e(vindx)-e(i))
          ftmp=vpqrs(vindx,vindx,i,vindx)**2
          pterm=pterm-delta_ai*ftmp
       enddo

!-----------------------------------------------------------------------
! Output the various contributions
!-----------------------------------------------------------------------
       tot=term0+term1(1)+term1(2)+u1p1h+rterm+sterm+pterm

       write(ilog,'(/,a,2(2x,i2))') "c, v:",cindx,vindx
       write(ilog,*) "Total:",tot
       write(ilog,*) "0:",term0
       write(ilog,*) "1, total:",term1(1)+term1(2)
       write(ilog,*) "1, Coulomb:",term1(1)
       write(ilog,*) "1, Exchange:",term1(2)
       write(ilog,*) "U_vc(1p1h):",u1p1h
       write(ilog,*) "R:",rterm
       write(ilog,*) "S:",sterm
       write(ilog,*) "P:",pterm

       return

     end subroutine orbrlx2
       
!#######################################################################

    end module adc2specmod
