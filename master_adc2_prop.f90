  module adc2mod

    contains

!#######################################################################

      subroutine master_adc2_prop()

        use constants
        use parameters
        use band_lanczos

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
           write(6,'(/,70a)') ('*',i=1,70)
           write(6,'(2x,a)') &
                'Initial space ADC(2)-s excitation energies'
           write(6,'(70a)') ('*',i=1,70)
           
           itmp=1+nBas**2*4*nOcc**2
           call table2(ndim,davstates,ener(1:davstates),&
                rvec(:,1:davstates),tmvec(1:davstates),&
                osc_str(1:davstates),kpq,itmp)
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
! Determine the guess vectors for the band-Lanczos calculation
!
! Note that as part of this process, the transition moments between
! the ISs and the initial state are calculated in lanczos_guess_vecs
! and passed back in the travec array
!-----------------------------------------------------------------------
        allocate(travec(ndimf))

        call lanczos_guess_vecs(vec_init,ndim,ndimsf,&
             travec,ndimf,kpq,kpqf,mtmf)

!-----------------------------------------------------------------------
! Write the final space ADC(2)-s Hamiltonian matrix to file
!-----------------------------------------------------------------------
        write(6,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
        if (lcvsfinal) then
           call write_fspace_adc2_1_cvs(ndimf,kpqf(:,:),noffdf,'c')           
        else
           call write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')
        endif

!-----------------------------------------------------------------------
! Perform the band-Lanczos calculation
!-----------------------------------------------------------------------
        call master_lancdiag(ndimf,noffdf,'c')        

!-----------------------------------------------------------------------
! Calculate the transition moments and oscillator strengths between 
! the initial state and the Lanczos pseudo-states
!-----------------------------------------------------------------------            
        call tdm_lancstates(ndimf,ndimsf,travec,e_init,mtmf)

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
           call select_atom_is_cvs(kpqf(:,:))
           call select_atom_d_cvs(kpqf(:,:),-1)
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
        write(6,*) 'ADC(2) INITIAL Space dim',ndim
        write(6,*) 'ADC(2) FINAL Space dim',ndimf
        write(6,*) 'ADC(2) TOTAL Space dim WITHOUT GROUND STATE',ndimd
        write(6,*) 'dimension of various INITIAL configuration spaces'
        write(6,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0)
        write(6,*) 'dimension of various FINAL configuration spaces'
        write(6,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
        write(6,*) 'dimension of various TOTAL configuration spaces'
        write(6,*) kpqd(1,0),kpqd(2,0),kpqd(3,0),kpqd(4,0),kpqd(5,0)

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

        ! Davidson in the initial space
        write(6,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
        call  write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i') 
        call cpu_time(time)
        write(6,*) 'Time=',time," s"
        
         ! Block-Davidson diagonalisation
        call master_dav(ndim,noffd,'i',ndims)

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


        write(6,'(/,2x,a,x,a1,x,a)') 'Calculating the transition moments &
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
        call readdavvc(davstates,ener,rvec)

!-----------------------------------------------------------------------        
! Calculate the vector F_J = < Psi_J | D | Psi_0 > (mtm)
!-----------------------------------------------------------------------
!        call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:),1)

!-----------------------------------------------------------------------
! Contract the F-vector with the Davidson states to yield the
! transition moments between the ground state and the Davidson states
!-----------------------------------------------------------------------
        do i=1,davstates
           tmvec(i)=tm(ndim,rvec(:,i),mtm(:))
           osc_str(i)=(2.0d0/3.0d0)*ener(i)*tmvec(i)**2
        end do

        deallocate(mtm)

        return

      end subroutine initial_space_tdm

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

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf
        real(d), dimension(:), allocatable        :: mtmf

!-----------------------------------------------------------------------        
! Calculate the vector F_J = < Psi_J | D | Psi_0 >, where the Psi_J
! are the ISs spanning the final space
!
! N.B. this vector is saved to the mtmf array
!-----------------------------------------------------------------------
        allocate(mtmf(ndimf))
        call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:),0)

!-----------------------------------------------------------------------
! From the values of the elements of mtmf, determine which 1h1p 
! unit vectors will form the initial Lanczos vectors
! 
! N.B. the corresponding indices are written to the stvc_lbl array
!-----------------------------------------------------------------------
        call fill_stvc(ndimsf,mtmf(1:ndimsf))

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

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: j,ndim,ndimf,ndimsf
        integer, dimension(:), allocatable        :: indx_tra
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec

!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the dipole
! operator and the initial state vector
!-----------------------------------------------------------------------        
        call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,&
             travec)

!-----------------------------------------------------------------------
! From the values of the elements of travec, determine which 1h1p 
! unit vectors will form the initial Lanczos vectors
! 
! N.B. the corresponding indices are written to the stvc_lbl array
!-----------------------------------------------------------------------
        call fill_stvc(ndimsf,travec(1:ndimsf))
        
        return

      end subroutine guess_vecs_ex2ex

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

      end subroutine tdm_ex2lanc

!#######################################################################

      subroutine run_mp2(e0)
        
        use constants
        use parameters
        use diagnostics        
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
        if (ld2) call mp2_d2(d2)

!-----------------------------------------------------------------------
! Output results
!-----------------------------------------------------------------------
        e0=Ehf+E_MP2
 
        write(6,'(/,2x,90a)') ('*',i=1,90)
        write(6,'(2x,1a)') '*'
        atmp='* HF energy:'
        write(6,'(2x,a25,2x,F18.10)') adjustl(atmp),Ehf
        write(6,'(2x,1a)') '*'
        atmp='* MP2 correlation energy:'
        write(6,'(2x,a25,2x,F18.10)') adjustl(atmp),E_MP2
        write(6,'(2x,1a)') '*'
        atmp='* MP2 energy:'
        write(6,'(2x,a25,2x,F18.10)') adjustl(atmp),e0
        if (ld2) then
           write(6,'(2x,1a)') '*'
           atmp='* D2 diagnostic:'
           write(6,'(2x,a25,2x,F18.10)') adjustl(atmp),d2
        endif
        write(6,'(2x,1a)') '*'
        write(6,'(2x,90a)') ('*',i=1,90)

        return

      end subroutine run_mp2

!#######################################################################

    end module adc2mod
