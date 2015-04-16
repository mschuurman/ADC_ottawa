  module adc2mod

    contains

!#######################################################################

      subroutine master_adc2_prop()

        use constants
        use parameters
        use band_lanczos
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
                osc_str(1:davstates),kpq,itmp)

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
        call readdavvc(davstates,ener,rvec)

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
        use band_lanczos

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
        call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,&
             travec)

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
        call master_dav(ndimf,noffdf,'f',ndimsf)

        return

      end subroutine davidson_final_space_diag

!#######################################################################

      subroutine lanczos_final_space_diag(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
           vec_init,mtmf,noffdf)

        use constants
        use parameters
        use fspace
        use band_lanczos

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

      subroutine tdm_davstates_final(ndimf,ndimsf,travec,e_init,mtmf,&
           kpqf)
        
        use constants
        use parameters
        use iomod, only: freeunit
        use misc, only: dsortindxa1

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,idavf,&
                                                     i,k,ilbl,itmp
        integer, dimension(ndimf)                 :: indx
        real(d)                                   :: e_init,ftmp
        real(d), dimension(ndimf)                 :: travec,mtmf,&
                                                     fstate,coeffsq
        real(d), dimension(davstates_f)           :: ener,tdm,osc_str
        real(d), parameter                        :: tol=0.05d0

!-----------------------------------------------------------------------
! Open the file containing the final space Davidson states
!-----------------------------------------------------------------------
        call freeunit(idavf)
        open(idavf,file=davname_f,status='unknown',access='sequential',&
             form='unformatted')

!-----------------------------------------------------------------------
! Form the scalar product of each final state vector with the vector
! travec corresponding to the contraction of the IS representation of
! dipole operator with the initial state vector
!-----------------------------------------------------------------------
        do i=1,davstates_f
           ! Read the ith final state vector
           fstate=0.0d0
           read(idavf) k,ener(i),fstate(:)
           ! Form the scalar product of the ith final state vector
           ! and travec
           tdm(i)=dot_product(fstate,travec)
           osc_str(i)=(2.0d0/3.0d0)*ener(i)*tdm(i)**2
        enddo

!-----------------------------------------------------------------------
! Output the final state information
!-----------------------------------------------------------------------
        write(ilog,'(/,70a)') ('*',i=1,70)
        write(ilog,'(2x,a)') &
             'Final space ADC(2)-s excitation energies'
        write(ilog,'(70a)') ('*',i=1,70)

        rewind(idavf)

        do i=1,davstates_f

           read(idavf) itmp,ftmp,fstate(:)
           coeffsq=fstate**2
           call dsortindxa1("D",ndimf,coeffsq(:),indx(:))

           write(ilog,'(2/,2x,a,2x,i2)') 'Final Space State Number',i           
           if (ener(i)*27.211d0.lt.10.0d0) then
              write(ilog,'(2x,a,2x,F10.5)') &
                   'Excitation Energy (eV):',ener(i)*27.211d0
           else
              write(ilog,'(2x,a,3x,F10.5)') &
                   'Excitation Energy (eV):',ener(i)*27.211d0
           endif
           write(ilog,'(2x,a,F10.5)') 'Transition Dipole Moment:',tdm(i)
           write(ilog,'(2x,a,5x,F10.5)') 'Oscillator Strength:',osc_str(i)
           write(ilog,'(/,2x,a,/)') 'Dominant Configurations:'
           write(ilog,'(2x,25a)') ('*',k=1,25)
           write(ilog,'(3x,a)') 'j   k -> a  b    C_jkab'
           write(ilog,'(2x,25a)') ('*',k=1,25)
           do k=1,50
              ilbl=indx(k)
              if (abs(coeffsq(ilbl)).ge.tol) then
                 if (kpqf(4,ilbl).eq.-1) then
                    ! Single excitations
                    write(ilog,'(3x,i2,4x,a2,1x,i2,5x,F8.5)') &
                         kpqf(3,ilbl),'->',kpqf(5,ilbl),coeffsq(ilbl)
                 else
                    ! Double excitations
                    write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),2x,F8.5)') &
                         kpqf(3,ilbl),kpqf(4,ilbl),'->',kpqf(5,ilbl),&
                         kpqf(6,ilbl),coeffsq(ilbl)
                 endif
              endif
           enddo

           write(ilog,'(2x,25a)') ('*',k=1,25)

        enddo

!-----------------------------------------------------------------------
! Output the convoluted spectrum
!-----------------------------------------------------------------------
        if (gwidth.ne.0.0d0) call convolute(osc_str,ener)

!-----------------------------------------------------------------------
! Close the file containing the final space Davidson states
!-----------------------------------------------------------------------
        close(idavf)

        return

      end subroutine tdm_davstates_final

!#######################################################################

      subroutine convolute(osc_str,ener)

        use constants
        use parameters
        use iomod, only: freeunit

        implicit none

        integer                         :: i,n,npoints,ispec
        real(d), dimension(davstates_f) :: osc_str,ener
        real(d)                         :: alpha,ftmp,de,range,curren,&
                                           func

        call freeunit(ispec)
        open(ispec,file='spec.dat',form='formatted',status='unknown')


        ftmp=(8.0d0*log(2.0d0))**0.5d0
        alpha=gwidth/ftmp

        npoints=1000
        
!        ftmp=0.1d0*(ener(davstates_f)-ener(1))

        ftmp=2.0d0

        de=27.211d0*(ener(davstates_f)-ener(1))+2*ftmp
        de=de/npoints


        do n=1,npoints

           curren=ener(1)*27.211d0-ftmp+de*n
           func=0.0d0

           do i=1,davstates_f
              func=func+osc_str(i)*exp(-(curren-ener(i)*27.211d0)**2/2*alpha)
           enddo

           write(ispec,*) curren,func

        enddo

        close(ispec)

        return

      end subroutine convolute

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
