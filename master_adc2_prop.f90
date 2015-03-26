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
        real(d)                              :: e_init
        real(d), dimension(:,:), allocatable :: rvec
        real(d), dimension(:), allocatable   :: vec_init
        real*8, dimension(:), allocatable    :: mtmf        
        real(d)                              :: E_groundstate
        
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
        call initial_space_diag(time,kpq,ndim,ndims,noffd,vec_init)

!-----------------------------------------------------------------------
! Transition moments from the ground state to the Davidson states
!-----------------------------------------------------------------------        
        call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,kpq)

!-----------------------------------------------------------------------
! Output the results of initial space calculation
!-----------------------------------------------------------------------
        write(6,'(/,70a)') ('*',i=1,70)
        write(6,'(2x,a)') 'Initial space ADC(2)-s excitation energies'
        write(6,'(70a)') ('*',i=1,70)
 
       itmp=1+nBas**2*4*nOcc**2
        call table2(ndim,davstates,ener(1:davstates),rvec(:,1:davstates),&
             tmvec(1:davstates),osc_str(1:davstates),&
             kpq,itmp)

!-----------------------------------------------------------------------
! Set the initial state vector and energy
!-----------------------------------------------------------------------
        vec_init(:)=rvec(:,statenumber)
        e_init=ener(statenumber)

!-----------------------------------------------------------------------
! Determine the guess vectors for the band-Lanczos calculation
!-----------------------------------------------------------------------
        allocate(travec(ndimf))

        call lanczos_guess_vecs(vec_init,ndim,ndimsf,&
             travec,ndimf,kpq,kpqf)

!-----------------------------------------------------------------------
! Write the final space ADC(2)-s Hamiltonian matrix to file
!-----------------------------------------------------------------------
        write(6,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
        call  write_fspace_adc2_1(ndimf,kpqf(:,:),noffdf,'c')

!-----------------------------------------------------------------------
! Perform the band-Lanczos calculation
!-----------------------------------------------------------------------
        call master_lancdiag(ndimf,noffdf,'c')

!-----------------------------------------------------------------------
! Calculate the transition moments between the initial state and the
! Lanczos states
!-----------------------------------------------------------------------            
        call tdm_lancstates(ndimf,ndimsf,travec,e_init)

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy.
!
! We should really do this at the start.
! Also, we should calculate and output the t2 diagnostic so that
! we have an idea of how "multi-referency" the wavefunction is.
!-----------------------------------------------------------------------
        call MP2(E_MP2)
        E_groundstate = Ehf + E_MP2
        write(6,*) 'THE ADC2 GROUND STATE ENERGY AT THIS GEOMETRY IS',&
             E_groundstate
  
        deallocate(ener,rvec,vec_init)
        deallocate(travec)          
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
        call  select_atom_is(kpq(:,:))
        call  select_atom_d(kpq(:,:),-1)

        ! Final subspace
        kpqf(:,:)=-1
        call  select_atom_isf(kpqf(:,:))
        call  select_atom_df(kpqf(:,:),-1)
        
        ! Total subspace
        kpqd(:,:)=-1
        call  select_atom_ist(kpqd(:,:))
        call  select_atom_dt(kpqd(:,:),-1)

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

      subroutine initial_space_diag(time,kpq,ndim,ndims,noffd,vec_init)
        
        use constants
        use parameters
        use fspace
        use davmod
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
        integer                                   :: ndim,ndims
        integer*8                                 :: noffd
        real(d), dimension(:), allocatable        :: vec_init
        real(d)                                   :: time

        ! Davidson in the initial space
        write(6,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
        call  write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i') 
        write(111,*) noffd
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
!        call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))

!-----------------------------------------------------------------------
! Contract the F-vector with the davidson states to yield the
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
           kpq,kpqf)

        use constants
        use parameters

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
        integer                                   :: ndim,ndimf,ndimsf
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: travec

!-----------------------------------------------------------------------
! Ionisation from the ground state
!-----------------------------------------------------------------------
        if (statenumber.eq.0) call guess_vecs_gs2ex

!-----------------------------------------------------------------------
! Ionisation from an excited state
!-----------------------------------------------------------------------
        if (statenumber.gt.0) call guess_vecs_ex2ex(vec_init,ndim,&
             ndimsf,travec,ndimf,kpq,kpqf)

        return

      end subroutine lanczos_guess_vecs

!#######################################################################

      subroutine guess_vecs_gs2ex

        use constants
        use parameters

        implicit none

        print*,"WRITE THIS PART OF THE CODE!"
        STOP

!-----------------------------------------------------------------------        
! Calculate the vector F_J = < Psi_J | D | Psi_0 >, where the Psi_J
! are the ISs spanning the final space
!
! N.B. this vector is saved to the mtmf array
!-----------------------------------------------------------------------
!        call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:))

!-----------------------------------------------------------------------
! From the values of the elements of mtmf, determine which 1h1p 
! unit vectors will form the initial Lanczos vectors
! 
! N.B. the corresponding indices are written to the stvc_lbl array
!-----------------------------------------------------------------------
!        call fill_stvc(ndimsf,mtmf(1:ndimsf))

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

      subroutine tdm_lancstates(ndimf,ndimsf,travec,e_init)
        
        use constants
        use parameters

        implicit none

        integer                   :: ndimf,ndimsf
        real(d), dimension(ndimf) :: travec
        real(d)                   :: e_init

!-----------------------------------------------------------------------
! Transition moments from the ground state
!-----------------------------------------------------------------------
        if (statenumber.eq.0) call tdm_gs2lanc

!-----------------------------------------------------------------------
! Transition moments from an excited state
!-----------------------------------------------------------------------
        if (statenumber.gt.0) call tdm_ex2lanc(ndimf,ndimsf,travec,&
             e_init)

        return

      end subroutine tdm_lancstates

!#######################################################################

      subroutine tdm_gs2lanc
        
        implicit none

        print*,"WRITE THE CODE INVOLVING mtmf"
        STOP

!        if (allocated(enerf)) deallocate(enerf)
!        if (allocated(tmvec)) deallocate(tmvec)
!        allocate(enerf(lancstates),tmvec(lancstates),mtmf(ndimf))
!
!!!$ ***lancstates is at most ncyclesXmain, usually we expect it to be lower than that
!        
!        write(6,*) &
!             'Calculating ADC2 transition moments from ground state in',&
!             tranmom2,' direction.'
!        write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,&
!             ' direction, i.e. : < PSI0  D',tranmom2,' PSIm  > IN THE FINAL SPACE'
!        tmvecf(:) = 0.d0
!  
!        call get_tranmom_1(ndimf,lancstates,lancname,mtmf(:),nstates,enerf(:),&
!             tmvecf(:),ndimsf)
!  
!        allocate(osc_strf(nstates))
!        osc_strf(:) = 0.d0
!        do i = 1 , nstates
!           osc_strf(i) = 2._d/3._d * enerf(i) * tmvecf(i)**2
!        end do

!        call get_sigma(nstates,excit(1:nstates),os2cs*osc_strf(:))

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
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(enerf(lancstates),tmvecf(lancstates),&
             osc_strf(nstates))

!-----------------------------------------------------------------------
! Calculate the transition moments between the initial state and the
! Lanczos states
!-----------------------------------------------------------------------
        tmvecf=0.0d0
        call get_tranmom_3(ndimf,lancstates,lancname,travec(:),&
             nstates,enerf(:),tmvecf(:),ndimsf)

!-----------------------------------------------------------------------
! Calculate and output the oscillator strengths
!-----------------------------------------------------------------------
        allocate(excit(nstates))
        do i=1,nstates
           excit(i)=enerf(i)-e_init
        end do
        
        osc_strf=0.0d0
        do i=1,nstates
           osc_strf(i) = (2.0d0/3.0d0)*excit(i)*tmvecf(i)**2
        enddo

        call get_sigma(nstates,excit(1:nstates),os2cs*osc_strf(:))

        return

      end subroutine tdm_ex2lanc

!#######################################################################

    end module adc2mod
