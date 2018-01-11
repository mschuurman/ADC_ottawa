!#######################################################################
! Calculation of TPA or TPXAS spectra using the Lanczos pseudospectrum
! of the ADC(2) Hamiltonian
!#######################################################################

module adc2tpamod

  use channels

contains

!#######################################################################

  subroutine adc2_tpa(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    
    implicit none

    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimsf,&
                                            nout,ndimf,ndimd,noutf,&
                                            itmp
    integer*8                            :: noffd,noffdf
    real(d)                              :: e_init,e0,time
    real(d), dimension(:), allocatable   :: ener,mtm,mtmf,tmvec,osc_str,&
                                            vec_init,travec
    real(d), dimension(:,:), allocatable :: rvec
    type(gam_structure)                  :: gam

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
! Block-Davidson diagonalisation in the initial space
!-----------------------------------------------------------------------
    if (.not.lcvsfinal) &
         call initial_space_diag(time,kpq,ndim,ndims,noffd)

!-----------------------------------------------------------------------
! If requested, calculate the dipole moments for the initial states
!-----------------------------------------------------------------------
    if (ldipole.and.statenumber.gt.0) &
         call initial_space_dipole(ndim,ndims,kpq)

!-----------------------------------------------------------------------
! Transition moments from the ground state to the Davidson states
!-----------------------------------------------------------------------
    if (.not.lcvsfinal) &
         call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,kpq)

!-----------------------------------------------------------------------
! Output the results of initial space calculation
!-----------------------------------------------------------------------
    if (.not.lcvsfinal) then
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
! If we are performing a Lanczos-TPXAS calculation and the initial
! state is the ground state, then we must calculate and store the
! initial space Hamiltonian for use in a block-Lanczos calculation.
! Note that for an excited initial state, this has already been done.
!-----------------------------------------------------------------------
    if (ltpa.and.statenumber.eq.0.and.lcvsfinal) &
         call initial_space_hamiltonian(ndim,kpq,noffd)
    
!-----------------------------------------------------------------------
! Calculation of the final space states
!-----------------------------------------------------------------------
    call final_space_diag_tpa(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
         vec_init,mtmf,noffd,noffdf,rvec)

!-----------------------------------------------------------------------
! If requested, calculate the dipole moments for the final states
! N.B. This is only meaningful if the final states are Davidson states
!-----------------------------------------------------------------------
    if (ldipole.and.ldiagfinal) &
         call final_space_dipole(ndimf,ndimsf,kpqf)
    
!-----------------------------------------------------------------------
! Calculate and output the TPA cross-sections
!-----------------------------------------------------------------------
    call final_space_tdm_tpa(ndimf,ndimsf,travec,e_init,mtmf,kpqf,&
         ndim)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    if (allocated(ener)) deallocate(ener)
    if (allocated(rvec)) deallocate(rvec)
    if (allocated(vec_init)) deallocate(vec_init)
    if (allocated(travec)) deallocate(travec)
    if (allocated(dipmom)) deallocate(dipmom)
    if (allocated(dipmom_f)) deallocate(dipmom_f)
    if (allocated(dpl_all)) deallocate(dpl_all)
    if (allocated(travec_ic)) deallocate(travec_ic)
    if (allocated(travec_iv)) deallocate(travec_iv)
    if (allocated(travec_fc)) deallocate(travec_fc)
    if (allocated(travec_fv)) deallocate(travec_fv)
    if (allocated(tdmgsf)) deallocate(tdmgsf)
    if (allocated(edavf)) deallocate(edavf)
    deallocate(kpq,kpqf,kpqd)
    
    return
    
  end subroutine adc2_tpa

!#######################################################################

  subroutine initial_space_hamiltonian(ndim,kpq,noffd)

    use constants
    use parameters
    use fspace
        
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
    integer                                   :: ndim
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
    
    return
    
  end subroutine initial_space_hamiltonian

!#######################################################################

  subroutine final_space_diag_tpa(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
       vec_init,mtmf,noffd,noffdf,rvec)

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

!-----------------------------------------------------------------------
! Generate the contractions of the initial and final state vectors
! the IS representation of the shifted dipole operator
!-----------------------------------------------------------------------
    if (lcvsfinal) then
       call davidson_final_space_diag(ndim,ndimf,ndimsf,&
            kpq,kpqf,travec,vec_init,mtmf,noffdf)
       call dipole_ispace_contraction_tpxas(ndim,ndimf,kpq,kpqf)
    else
       call dipole_ispace_contraction_tpa(ndim,ndimf,kpq,kpqf)
    endif

!-----------------------------------------------------------------------
! Generate the Lanczos pseudospectra of the ADC(2) Hamiltonian(s)
! (For a TPXAS calculation, we generate the Lanczos
! pseudospectrum of both the ADC(2) and CVS-ADC(2) Hamiltonians)
!-----------------------------------------------------------------------
    call tpa_lanczos(ndim,ndimf,kpq,kpqf,noffd,noffdf)

    return

  end subroutine final_space_diag_tpa
    
!#######################################################################

  subroutine dipole_ispace_contraction_tpxas(ndim,ndimf,kpq,kpqf)

    use constants
    use parameters
    use adc2common
    use misc
    use get_matrix_dipole
    use get_moment
    use fspace
    use guessvecs
    use iomod
    use diagmod

    implicit none
    
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndimf,c,k,&
                                                 f,i,j,a
    integer*8, dimension(3)                   :: nel_cv,nel_cc,&
                                                 nel_vv
    integer, dimension(3)                     :: nbuf_cv,nbuf_cc,&
                                                 nbuf_vv
    integer                                   :: dim,error,ivecs
    real(d), dimension(:,:), allocatable      :: veci,vecf,mtm_v,&
                                                 mtm_c
    real(d), dimension(:), allocatable        :: ei,ef
    real(d), dimension(:,:), allocatable      :: tdmvec,initvecs
    real(d), dimension(:), allocatable        :: tau,work
    character(len=1), dimension(3)            :: acomp
    character(len=70)                         :: msg
    character(len=60)                         :: filename
        
    integer                              :: e2
    real(d), dimension(:,:), allocatable :: smat
    real(d), dimension(:), allocatable   :: eig
    
    acomp=(/ 'x','y','z' /)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(travec_iv(ndim,3))
    allocate(travec_fv(ndim,3,davstates_f))
    allocate(travec_ic(ndimf,3))
    allocate(travec_fc(ndimf,3,davstates_f))
    allocate(vecf(ndimf,davstates_f))
    allocate(ef(davstates_f))
    if (statenumber.gt.0) then
       allocate(veci(ndim,davstates))
       allocate(ei(davstates))
    endif

!-----------------------------------------------------------------------
! Transition matrix elements: core-excited states
!-----------------------------------------------------------------------
    ! (i) Core, valence
    !
    ! Loop over the components of the dipole operator
    do c=1,3
       ! Set the dipole component
       dpl(:,:)=dpl_all(c,:,:)
       ! Calculate the IS representation of the shifter
       ! dipole operator
       write(ilog,'(70a)') ('-',k=1,70)
       msg='Calculation of P_c . D_'//acomp(c)//' . P_v'
       write(ilog,'(2x,a)') trim(msg)
       write(ilog,'(70a)') ('-',k=1,70)
       filename='SCRATCH/dipole_cv_'//acomp(c)
       call get_adc2_dipole_improved_omp(ndimf,ndim,kpqf,kpq,&
            nbuf_cv(c),nel_cv(c),filename)
    enddo
        
    ! (ii) Core, core
    !
    ! Loop over the components of the dipole operator
    do c=1,3         
       ! Set the dipole component
       dpl(:,:)=dpl_all(c,:,:)
       ! Calculate the IS representation of the shifted
       ! dipole operator
       write(ilog,'(70a)') ('-',k=1,70)
       msg='Calculation of P_c . D_'//acomp(c)//' . P_c'
       write(ilog,'(2x,a)') trim(msg)
       write(ilog,'(70a)') ('-',k=1,70)
       filename='SCRATCH/dipole_cc_'//acomp(c)
       call get_adc2_dipole_improved_omp(ndimf,ndimf,kpqf,kpqf,&
            nbuf_cc(c),nel_cc(c),filename)
    enddo

!-----------------------------------------------------------------------
! Transition matrix elements: initial state
!-----------------------------------------------------------------------        
    if (statenumber.eq.0) then
        
       ! (i) Valence-excited
       !
       ! Loop over the components of the dipole operator
       do c=1,3
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          ! Calculate the F-vector
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_v . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          call get_modifiedtm_adc2(ndim,kpq(:,:),travec_iv(:,c),1)
       enddo
           
       ! (ii) Core-excited
       !
       ! Loop over the components of the dipole operator
       do c=1,3
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          ! Calculate the F-vector
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_c . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          call get_modifiedtm_adc2(ndimf,kpqf(:,:),travec_ic(:,c),0)
       enddo
       
    else
           
       ! (i) Valence, valence
       !
       ! Loop over the components of the dipole operator
       do c=1,3
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          ! Calculate the IS representation of the shifted
          ! dipole operator
          filename='SCRATCH/dipole_vv_'//acomp(c)
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_v . D_'//acomp(c)//' . P_v'
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          call get_adc2_dipole_improved_omp(ndim,ndim,kpq,kpq,&
               nbuf_vv(c),nel_vv(c),filename)
       enddo

    endif

!-----------------------------------------------------------------------
! Read the Davidson vectors from file
!-----------------------------------------------------------------------
    ! Final states
    call readdavvc(davstates_f,ef,vecf,'f',ndimf)
        
    ! Save the final space Davidson energies for use later in
    ! the calculation of the two-photon transition moments
    allocate(edavf(davstates_f))
    edavf=ef
     
    ! Initial states
    if (statenumber.gt.0) call readdavvc(davstates,ei,veci,'i',ndim)

!-----------------------------------------------------------------------
! Perform the contractions of the shifted dipole matrices with the
! Davidson state vectors (final and initial spaces)
!-----------------------------------------------------------------------
    ! Final states
    do f=1,davstates_f           
       ! (i) Final | D | valence ISs
       do c=1,3
          filename='SCRATCH/dipole_cv_'//acomp(c)
          call contract_dipole_state(filename,ndimf,ndim,&
               vecf(:,f),travec_fv(:,c,f),nbuf_cv(c),nel_cv(c),'l')
       enddo
       ! (ii) Final | D | core ISs
       do c=1,3
          filename='SCRATCH/dipole_cc_'//acomp(c)
          call contract_dipole_state(filename,ndimf,ndimf,&
               vecf(:,f),travec_fc(:,c,f),nbuf_cc(c),nel_cc(c),'l')
       enddo
    enddo

    ! Initial state
    if (statenumber.gt.0) then
       ! (i) Initial | D | valence ISs
       do c=1,3
          filename='SCRATCH/dipole_vv_'//acomp(c)
          call contract_dipole_state(filename,ndim,ndim,&
               veci(:,statenumber),travec_iv(:,c),nbuf_vv(c),&
               nel_vv(c),'r')
       enddo
       ! (ii) Initial | D | core ISs
       do c=1,3
          filename='SCRATCH/dipole_cv_'//acomp(c)
          call contract_dipole_state(filename,ndim,ndimf,&
               veci(:,statenumber),travec_ic(:,c),nbuf_cv(c),&
               nel_cv(c),'r')
       enddo
    endif

!-----------------------------------------------------------------------
! Generate the initial Lanczos vectors via the orthogonalisation of
! the projected shifted-dipole-matrix state-vector products
!-----------------------------------------------------------------------
    ! (1) Valence-excited space
    !
    ! Copy the contents of the travec_iv and travec_cv arrays
    tpblock(1)=3+3*davstates_f
    allocate(initvecs(ndim,tpblock(1)))
    initvecs(:,1:3)=travec_iv(:,1:3)
    k=3
    do f=1,davstates_f
       initvecs(:,k+1:k+3)=travec_fv(:,1:3,f)
       k=k+3
    enddo
    
    !! TEST
    !!
    !! DIAGONALISATION OF THE OVERLAP MATRIX
    !allocate(smat(tpblock(1),tpblock(1)))
    !allocate(eig(tpblock(1)))
    !allocate(work(3*tpblock(1)))
    !
    !write(ilog,'(/,2x,a)') 'Valence-excited space overlap &
    !     matrix:'
    !do i=1,tpblock(1)
    !   do j=i,tpblock(1)
    !      smat(i,j)=dot_product(initvecs(:,i),initvecs(:,j))
    !      smat(j,i)=smat(i,j)
    !      write(ilog,*) i,j,smat(i,j)
    !   enddo
    !enddo
    !
    !e2=3*tpblock(1)
    !call dsyev('V','U',tpblock(1),smat,tpblock(1),eig,work,e2,error)
    !
    !if (error.ne.0) then
    !   errmsg='This fucked up...'
    !   call error_control
    !endif
    !
    !write(ilog,'(/,2x,a)') 'Eigenvalues of the valence-excited &
    !     space overlap matrix:'
    !do i=1,tpblock(1)
    !   write(ilog,*) i,eig(i)
    !enddo
    !
    !deallocate(smat)
    !deallocate(eig)
    !deallocate(work)
    !! TEST
    
    ! Orthogonalisation of the shifted dipole matrix-state vector
    ! contractions via a QR factorisation
    allocate(tau(tpblock(1)))
    allocate(work(tpblock(1)))
    call dgeqrf(ndim,tpblock(1),initvecs,ndim,tau,work,&
         tpblock(1),error)
    if (error.ne.0) then
       errmsg='dqerf failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
    call dorgqr(ndim,tpblock(1),tpblock(1),initvecs,ndim,tau,&
         work,tpblock(1),error)
    if (error.ne.0) then
       errmsg='dorgqr failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
        
    ! Write the valence-excited space guess vectors to file
    call freeunit(ivecs)
    open(ivecs,file='SCRATCH/tpa_initi',form='unformatted',&
         status='unknown')
    write(ivecs) initvecs
    close(ivecs)
    
    deallocate(initvecs)
    deallocate(tau)
    deallocate(work)

    ! (2) Core-excited space
    !
    tpblock(2)=3+3*davstates_f
    allocate(initvecs(ndimf,tpblock(2)))
    initvecs(:,1:3)=travec_ic(:,1:3)
    k=3
    do f=1,davstates_f
       initvecs(:,k+1:k+3)=travec_fc(:,1:3,f)
       k=k+3
    enddo

    !! TEST
    !!
    !! DIAGONALISATION OF THE OVERLAP MATRIX
    !allocate(smat(tpblock(2),tpblock(2)))
    !allocate(eig(tpblock(2)))
    !allocate(work(3*tpblock(2)))
    !
    !write(ilog,'(/,2x,a)') 'Core-excited space overlap &
    !     matrix:'
    !do i=1,tpblock(2)
    !   do j=i,tpblock(2)
    !      smat(i,j)=dot_product(initvecs(:,i),initvecs(:,j))
    !      smat(j,i)=smat(i,j)
    !      write(ilog,*) i,j,smat(i,j)
    !   enddo
    !enddo
    !
    !e2=3*tpblock(2)
    !call dsyev('V','U',tpblock(2),smat,tpblock(2),eig,work,e2,error)
    !
    !if (error.ne.0) then
    !   errmsg='This fucked up...'
    !   call error_control
    !endif
    !
    !write(ilog,'(/,2x,a)') 'Eigenvalues of the core-excited &
    !     space overlap matrix:'
    !do i=1,tpblock(2)
    !   write(ilog,*) i,eig(i)
    !enddo
    !
    !deallocate(smat)
    !deallocate(eig)
    !deallocate(work)
    !! TEST

    ! Orthogonalisation of the shifted dipole matrix-state vector
    ! contractions via a QR factorisation
    allocate(tau(tpblock(2)))
    allocate(work(tpblock(2)))
    call dgeqrf(ndimf,tpblock(2),initvecs,ndimf,tau,work,&
         tpblock(2),error)
    if (error.ne.0) then
       errmsg='dqerf failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
    call dorgqr(ndimf,tpblock(2),tpblock(2),initvecs,ndimf,tau,&
         work,tpblock(2),error)
    if (error.ne.0) then
       errmsg='dorgqr failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
           
    ! Write the core-excited space guess vectors to file
    call freeunit(ivecs)
    open(ivecs,file='SCRATCH/tpa_initc',form='unformatted',&
         status='unknown')
    write(ivecs) initvecs
    close(ivecs)
    
    deallocate(initvecs)
    deallocate(tau)
    deallocate(work)

!-----------------------------------------------------------------------
! If we are considering two-photon excitation from an excited state,
! then we additionally require the transition dipole moments between
! the ground state and the initial state, and between the ground state
! and all final states
!-----------------------------------------------------------------------
    if (statenumber.ne.0) then
           
       ! Allocate arrays
       allocate(mtm_v(ndim,3))
       allocate(mtm_c(ndimf,3))
       allocate(tdmgsf(3,davstates_f))

       ! (i) Valence
       do c=1,3

          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          
          ! Calculate the F-vector
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_v . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          call get_modifiedtm_adc2(ndim,kpq(:,:),mtm_v(:,c),1)

          ! Calculate the transition dipole moment between the
          ! ground state and the initial state
          tdmgsi(c)=dot_product(mtm_v(:,c),veci(:,statenumber))
          
       enddo

       ! (ii) Core
       do c=1,3
              
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          
          ! Calculate the F-vector
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_c . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtm_c(:,c),0)
          
          ! Calculate the transition dipole moment between the
          ! ground state and the final states
          do f=1,davstates_f
             tdmgsf(c,f)=dot_product(mtm_c(:,c),vecf(:,f))
          enddo
          
       enddo
       
       ! Deallocate arrays
       deallocate(mtm_v)
       deallocate(mtm_c)

    endif

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(vecf)
    deallocate(ef)
    if (allocated(veci)) deallocate(veci)
    if (allocated(ei)) deallocate(ei)
    
    return
    
  end subroutine dipole_ispace_contraction_tpxas
  
!#######################################################################

        
  subroutine dipole_ispace_contraction_tpa(ndim,ndimf,kpq,kpqf)

    use constants
    use parameters
    use adc2common
    use misc
    use get_matrix_dipole
    use get_moment
    use fspace
    use guessvecs
    use iomod
    use diagmod

    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndimf,c,&
                                                 k,f,i,j,error,&
                                                 ivecs,a,b
    integer*8, dimension(3)                   :: nel_vv
    integer, dimension(3)                     :: nbuf_vv
    real(d), dimension(:,:), allocatable      :: vec,initvecs
    real(d), dimension(:), allocatable        :: ener
    real(d), dimension(:), allocatable        :: tau,work
    real(d), dimension(:,:), allocatable      :: mtm
    character(len=1), dimension(3)            :: acomp
    character(len=70)                         :: msg
    character(len=60)                         :: filename
    
    integer                              :: e2
    real(d), dimension(:,:), allocatable :: smat
    real(d), dimension(:), allocatable   :: eig
    
    acomp=(/ 'x','y','z' /)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(travec_iv(ndim,3))
    allocate(travec_fv(ndim,3,davstates))
    allocate(vec(ndim,davstates))
    allocate(ener(davstates))

!-----------------------------------------------------------------------
! P_v . D_a . P_v, a=x,y,z
!-----------------------------------------------------------------------
    ! Loop over the components of the dipole operator
    do c=1,3
           
       ! Set the dipole component
       dpl(:,:)=dpl_all(c,:,:)

       ! Output progress
       write(ilog,'(70a)') ('-',k=1,70)
       msg='Calculation of P_v . D_'//acomp(c)//' . P_v'
       write(ilog,'(2x,a)') trim(msg)
       write(ilog,'(70a)') ('-',k=1,70)

       ! Calculate the IS representation of the shifted
       ! dipole operator
       filename='SCRATCH/dipole_vv_'//acomp(c)
       call get_adc2_dipole_improved_omp(ndim,ndim,kpq,kpq,&
            nbuf_vv(c),nel_vv(c),filename)
       
    enddo

!-----------------------------------------------------------------------
! Read the Davidson vectors from file
!-----------------------------------------------------------------------
    call readdavvc(davstates,ener,vec,'i',ndim)

    ! Save the Davidson energies for use later in the
    ! calculation of the two-photon transition moments
    allocate(edavf(davstates))
    edavf=ener

!-----------------------------------------------------------------------
! P_v . D_a . Psi_i, a=x,y,z
!-----------------------------------------------------------------------
    if (statenumber.eq.0) then
           
       ! Loop over the components of the dipole operator
       do c=1,3
          
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          
          ! Output progress
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_v . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)
          
          ! Calculate the F-vector
          call get_modifiedtm_adc2(ndim,kpq(:,:),travec_iv(:,c),1)
          
       enddo
       
    else

       ! Loop over the components of the dipole operator
       do c=1,3
          
          ! Calculate P_v . D_c . Psi_i
          filename='SCRATCH/dipole_vv_'//acomp(c)
          call contract_dipole_state(filename,ndim,ndim,&
               vec(:,statenumber),travec_iv(:,c),nbuf_vv(c),&
               nel_vv(c),'r')
          
       enddo

    endif

!-----------------------------------------------------------------------
! P_v . D_a . Psi_f, a=x,y,z, f=1,davstates
!-----------------------------------------------------------------------
    ! Loop over the Davidson states
    do f=1,davstates
       
       ! Loop over the components of the dipole operator
       do c=1,3
          
          ! Calculate P_v . D_c . Psi_f
          filename='SCRATCH/dipole_vv_'//acomp(c)
          call contract_dipole_state(filename,ndimf,ndim,&
               vec(:,f),travec_fv(:,c,f),nbuf_vv(c),nel_vv(c),'r')

       enddo

    enddo

!-----------------------------------------------------------------------
! Generate the initial Lanczos vectors via the orthogonalisation of
! the shifted-dipole-matrix state-vector products
!-----------------------------------------------------------------------
    ! Set the block size
    if (statenumber.eq.0) then
       tpblock(1)=3+3*davstates
    else
       tpblock(1)=3+3*(davstates-1)
    endif
    
    ! Allocate arrays
    allocate(initvecs(ndim,tpblock(1)))
        
    ! Copy the contents of the travec_iv and travec_fv arrays
    initvecs(:,1:3)=travec_iv(:,1:3)
    k=3
    do f=1,davstates
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       initvecs(:,k+1:k+3)=travec_fv(:,1:3,f)
       k=k+3
    enddo

    !! TEST
    !!
    !! DIAGONALISATION OF THE OVERLAP MATRIX
    !allocate(smat(tpblock(1),tpblock(1)))
    !allocate(eig(tpblock(1)))
    !allocate(work(3*tpblock(1)))
    !
    !write(ilog,'(/,2x,a)') 'Valence-excited space overlap &
    !     matrix:'
    !do i=1,tpblock(1)
    !   do j=i,tpblock(1)
    !      smat(i,j)=dot_product(initvecs(:,i),initvecs(:,j))
    !      smat(j,i)=smat(i,j)
    !      write(ilog,*) i,j,smat(i,j)
    !   enddo
    !enddo
    !
    !e2=3*tpblock(1)
    !call dsyev('V','U',tpblock(1),smat,tpblock(1),eig,work,e2,error)
    !
    !if (error.ne.0) then
    !   errmsg='This fucked up...'
    !   call error_control
    !endif
    !
    !write(ilog,'(/,2x,a)') 'Eigenvalues of the valence-excited &
    !     space overlap matrix:'
    !do i=1,tpblock(1)
    !   write(ilog,*) i,eig(i)
    !enddo
    !
    !deallocate(smat)
    !deallocate(eig)
    !deallocate(work)
    !! TEST
        
    ! Orthogonalisation of the shifted dipole matrix-state vector
    ! contractions via a QR factorisation
    allocate(tau(tpblock(1)))
    allocate(work(tpblock(1)))
    call dgeqrf(ndim,tpblock(1),initvecs,ndim,tau,work,&
         tpblock(1),error)
    if (error.ne.0) then
       errmsg='dqerf failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
    call dorgqr(ndim,tpblock(1),tpblock(1),initvecs,ndim,tau,&
         work,tpblock(1),error)
    if (error.ne.0) then
       errmsg='dorgqr failed in subroutine &
            dipole_ispace_contraction_tpa'
       call error_control
    endif
    
    ! Write the valence-excited space guess vectors to file
    call freeunit(ivecs)
    open(ivecs,file='SCRATCH/tpa_initi',form='unformatted',&
         status='unknown')
    write(ivecs) initvecs
    close(ivecs)

    ! Deallocate arrays
    deallocate(tau)
    deallocate(work)
    deallocate(initvecs)

!-----------------------------------------------------------------------
! If we are considering two-photon excitation from an excited state,
! then we additionally require the transition dipole moments between
! the ground state and all Davidson states
!-----------------------------------------------------------------------
    if (statenumber.ne.0) then

       ! Allocate arrays
       allocate(mtm(ndim,3))
       allocate(tdmgsf(3,davstates))

       ! Loop over the components of the dipole operator
       do c=1,3
          
          ! Set the dipole component
          dpl(:,:)=dpl_all(c,:,:)
          
          ! Output our progress
          write(ilog,'(70a)') ('-',k=1,70)
          msg='Calculation of P_v . F_'//acomp(c)
          write(ilog,'(2x,a)') trim(msg)
          write(ilog,'(70a)') ('-',k=1,70)

          ! Calculate the F-vector
          call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:,c),1)
          
          ! Calculate the transition dipole moment between the
          ! ground state and the Davidson states
          do f=1,davstates
             tdmgsf(c,f)=dot_product(mtm(:,c),vec(:,f))
          enddo
          tdmgsi(:)=tdmgsf(:,statenumber)
          
       enddo
           
       ! Deallocate arrays
       deallocate(mtm)
       
    endif

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(vec)
    deallocate(ener)
    
    return
    
  end subroutine dipole_ispace_contraction_tpa

!#######################################################################  

  subroutine tpa_lanczos(ndim,ndimf,kpq,kpqf,noffd,noffdf)

    use constants
    use parameters
    use block_lanczos

    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndimf
    integer*8                                 :: noffd,noffdf

!-----------------------------------------------------------------------
! Perform the block-Lanczos calculation using the initial-space
! Hamiltonian
!-----------------------------------------------------------------------
    lmain=tpblock(1)
    call lancdiag_block(ndim,noffd,'i')

    ! Rename the valence-excited space Lanczos vector file
    call system('mv '//trim(lancname)//' SCRATCH/lancstates_v')

!-----------------------------------------------------------------------
! If we are calculating a TPXAS spectrum, then perform the
! block-Lanczos calculation using the final-space Hamiltonian
!-----------------------------------------------------------------------
    if (lcvsfinal) then

       lmain=tpblock(2)
       call lancdiag_block(ndimf,noffdf,'c')
           
       ! Rename the core-excited space Lanczos vector file
       call system('mv '//trim(lancname)//' SCRATCH/lancstates_c')
       
    endif

    return

  end subroutine tpa_lanczos

!#######################################################################  

  subroutine final_space_tdm_tpa(ndimf,ndimsf,travec,e_init,mtmf,kpqf,&
       ndim)

    use constants
    use parameters

    implicit none
    
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
    integer                                   :: ndimf,ndimsf,ndim
    real(d), dimension(ndimf)                 :: travec,mtmf
    real(d)                                   :: e_init

    if (lcvsfinal) then
       call tdm_tpxas(ndim,ndimf,e_init)
    else
       call tdm_tpa(ndim,ndimf,e_init)
    endif
               
    return
    
  end subroutine final_space_tdm_tpa

!#######################################################################  

  subroutine tdm_tpxas(ndim,ndimf,e_init)

    use constants
    use parameters
    use iomod
    
    implicit none
    
    integer                                :: ndim,ndimf,f,k,i,&
                                              nlanc_v,nlanc_c,&
                                              ilanc,a,b,alpha,&
                                              itpa
    real(d)                                :: e_init
    real(d), dimension(:,:,:), allocatable :: sabif
    real(d), dimension(:), allocatable     :: lvec,lener_v,&
                                              lener_c,xsec,&
                                              delta
    real(d), dimension(:,:), allocatable   :: tdmil_v,tdmil_c
    real(d), dimension(:,:,:), allocatable :: tdmfl_v,tdmfl_c
    real(d), parameter                     :: c_au=137.0359991d0

!----------------------------------------------------------------------
! Output where we are at
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of the two-photon &
         cross-sections'
    write(ilog,'(70a)') ('-',k=1,70)

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    ! Two-photon transition moments
    allocate(sabif(3,3,davstates_f))
    sabif=0.0d0

    ! Valence-excited space Lanczos energies
    nlanc_v=tpblock(1)*ncycles
    allocate(lener_v(nlanc_v))
    lener_v=0.0d0

    ! Core-excited space Lanczos energies
    nlanc_c=tpblock(1)*ncycles
    allocate(lener_c(nlanc_c))
    lener_c=0.0d0

    ! Transition dipole moments between the initial state
    ! and the valence-excited space Lanczos pseudo-states
    allocate(tdmil_v(3,nlanc_v))
    tdmil_v=0.0d0

    ! Transition dipole moments between the initial state
    ! and the core-excited space Lanczos pseudo-states
    allocate(tdmil_c(3,nlanc_c))
    tdmil_c=0.0d0
    
    ! Transition dipole moments between the final states
    ! and the valence-excited space Lanczos pseudo-states
    allocate(tdmfl_v(3,davstates_f,nlanc_v))
    tdmfl_v=0.0d0

    ! Transition dipole moments between the final states
    ! and the core-excited space Lanczos pseudo-states
    allocate(tdmfl_c(3,davstates_f,nlanc_c))
    tdmfl_c=0.0d0

    ! Two photon cross-sections
    allocate(xsec(davstates_f))
    xsec=0.0d0
    allocate(delta(davstates_f))
    delta=0.0d0

!----------------------------------------------------------------------
! Calculate the transition dipole moments between the initial state
! and the valence-excited space Lanczos pseudo-states, and the 
! final states and the valence-excited Lanczos pseudo-states
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(lvec(ndim))

    ! Open the valence-excited space Lanczos pseudo-state file
    call freeunit(ilanc)
    open(ilanc,file='SCRATCH/lancstates_v',status='old',&
         access='sequential',form='unformatted')

    ! Loop over the valence-excited space Lanczos pseudo-states
    do alpha=1,nlanc_v
       
       ! Read the current Lanczos pseudo-state
       read(ilanc) k,lener_v(alpha),lvec
       
       ! < i | D_a | alpha >
       do a=1,3
          tdmil_v(a,alpha)=dot_product(travec_iv(:,a),lvec)
       enddo
       
       ! < f | D_a | alpha >
       do f=1,davstates_f
          do a=1,3
             tdmfl_v(a,f,alpha)=&
                  dot_product(travec_fv(:,a,f),lvec)
          enddo
       enddo

    enddo

    ! Deallocate arrays
    deallocate(lvec)
    
    ! Close the valence-excited space Lanczos pseudo-state file
    close(ilanc)

!----------------------------------------------------------------------
! Calculate the transition dipole moments between the initial state
! and the core-excited space Lanczos pseudo-states, and the 
! final states and the core-excited Lanczos pseudo-states
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(lvec(ndimf))

    ! Open the valence-excited space Lanczos pseudo-state file
    call freeunit(ilanc)
    open(ilanc,file='SCRATCH/lancstates_c',status='old',&
         access='sequential',form='unformatted')

    ! Loop over the valence-excited space Lanczos pseudo-states
    do alpha=1,nlanc_c

       ! Read the current Lanczos pseudo-state
       read(ilanc) k,lener_c(alpha),lvec
       
       ! < i | D_a | alpha >
       do a=1,3
          tdmil_c(a,alpha)=dot_product(travec_ic(:,a),lvec)
       enddo

       ! < f | D_a | alpha >
       do f=1,davstates_f
          do a=1,3
             tdmfl_c(a,f,alpha)=&
                  dot_product(travec_fc(:,a,f),lvec)
          enddo
       enddo

    enddo

    ! Deallocate arrays
    deallocate(lvec)

    ! Close the valence-excited space Lanczos pseudo-state file
    close(ilanc)

!----------------------------------------------------------------------
! Calculation of the valence-excited space contribution to the
! two-photon transition moments
!----------------------------------------------------------------------
    ! Loop over final states
    do f=1,davstates_f
           
       ! Loop over pairs of dipole operator components
       do a=1,3
          do b=1,3

             ! Loop over valence-excited space Lanczos
             ! pseudo-states
             do alpha=1,nlanc_v

                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmil_v(a,alpha)*tdmfl_v(b,f,alpha)&
                     /(lener_v(alpha)-0.5d0*edavf(f))

                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmil_v(b,alpha)*tdmfl_v(a,f,alpha)&
                     /(lener_v(alpha)-0.5d0*edavf(f))

             enddo

          enddo
       enddo
           
    enddo

!----------------------------------------------------------------------
! Calculation of the core-excited space contribution to the
! two-photon transition moments
!----------------------------------------------------------------------
    ! Loop over final states
    do f=1,davstates_f
           
       ! Loop over pairs of dipole operator components
       do a=1,3
          do b=1,3

             ! Loop over valence-excited space Lanczos
             ! pseudo-states
             do alpha=1,nlanc_c

                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmil_c(a,alpha)*tdmfl_c(b,f,alpha)&
                     /(lener_c(alpha)-0.5d0*edavf(f))

                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmil_c(b,alpha)*tdmfl_c(a,f,alpha)&
                     /(lener_c(alpha)-0.5d0*edavf(f))
                
             enddo

          enddo
       enddo

    enddo

!----------------------------------------------------------------------
! Contribution of the ground state to the two-photon transition
! moments (excited initial states only)
! N.B., E_0 = 0, hence the E_f/2 denominators
!----------------------------------------------------------------------
    if (statenumber.ne.0) then

       ! Loop over final states
       do f=1,davstates_f

          ! Loop over pairs of dipole operator components
          do a=1,3
             do b=1,3
                    
                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmgsi(a)*tdmgsf(b,f)&
                     /(-0.5d0*edavf(f))
                    
                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmgsi(b)*tdmgsf(a,f)&
                     /(-0.5d0*edavf(f))
                    
             enddo
          enddo
          
       enddo
           
    endif

!----------------------------------------------------------------------
! Calculation of the two-photon cross-sections for (parallel) linearly
! polarised light
!----------------------------------------------------------------------
    do f=1,davstates_f
       do a=1,3
          do b=1,3
             delta(f)=delta(f)&
                  +2.0d0*sabif(a,a,f)*sabif(b,b,f)&
                  +2.0d0*sabif(a,b,f)*sabif(a,b,f)&
                  +2.0d0*sabif(a,b,f)*sabif(b,a,f)
          enddo
          ! Conversion to cross-sections
          xsec(f)=delta(f)*16.0d0*(pi**3)*((edavf(f)-e_init)/2.0d0)**2
          xsec(f)=xsec(f)/(c_au**2)
       enddo
    enddo
    xsec=xsec/30.0d0
    delta=delta/30.0d0

!----------------------------------------------------------------------
! Output the two-photon cross-sections
! N.B., we output the excitation energies relative to the initial
! state
!----------------------------------------------------------------------
    ! Log file
    write(ilog,'(/,2x,85a)') ('=',k=1,85)
    write(ilog,'(2x,a)') '     i     |     f     |     &
         E_if (eV)     |     delta_TP (au)     |     omega2 (au)'
    write(ilog,'(2x,85a)') ('-',k=1,85)
    do f=1,davstates_f
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       write(ilog,10) statenumber,'|',f,'|',&
            (edavf(f)-e_init)*eh2ev,'|',delta(f),'|',xsec(f)
    enddo
    write(ilog,'(2x,85a,/)') ('=',k=1,85)

    ! Stick spectrum
    call freeunit(itpa)
    open(itpa,file='tpa.dat',form='formatted',status='unknown')
    do f=1,davstates_f
       write(itpa,*) (edavf(f)-e_init)*eh2ev,xsec(f)
    enddo
    close(itpa)
    
10  format(6x,i2,5x,a1,4x,i2,5x,a1,3x,F10.5,6x,a1,3x,&
         F10.5,10x,a1,3x,F10.5)
    
    ! Gnuplot spectrum file
    call wrgnuplot(davstates_f,edavf-e_init,xsec)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sabif)
    deallocate(lener_v)
    deallocate(lener_c)
    deallocate(tdmil_v)
    deallocate(tdmil_c)
    deallocate(tdmfl_v)
    deallocate(tdmfl_c)
    deallocate(xsec)
    deallocate(delta)

    return

  end subroutine tdm_tpxas

!#######################################################################  

  subroutine tdm_tpa(ndim,ndimf,e_init)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer                                :: ndim,ndimf,k,&
                                              nlanc,ilanc,alpha,&
                                              a,b,f,itpa
    real(d)                                :: e_init
    real(d), dimension(:,:,:), allocatable :: sabif
    real(d), dimension(:), allocatable     :: lener,lvec,&
                                              xsec,delta
    real(d), dimension(:,:), allocatable   :: tdmil
    real(d), dimension(:,:,:), allocatable :: tdmfl
    real(d), parameter                     :: c_au=137.0359991d0

!----------------------------------------------------------------------
! Output where we are at
!----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',k=1,70)
    write(ilog,'(2x,a)') 'Calculation of the two-photon &
         cross-sections'
    write(ilog,'(70a)') ('-',k=1,70)

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    ! Two-photon transition moments
    allocate(sabif(3,3,davstates))
    sabif=0.0d0

    ! Lanczos energies
    nlanc=tpblock(1)*ncycles
    allocate(lener(nlanc))
    lener=0.0d0

    ! Transition dipole moments between the initial state and
    ! the Lanczos pseudo-states
    allocate(tdmil(3,nlanc))
    tdmil=0.0d0

    ! Transition dipole moments between the final states and
    ! the Lanczos pseudo-states
    allocate(tdmfl(3,davstates,nlanc))
    tdmfl=0.0d0

    ! Two photon cross-sections
    allocate(xsec(davstates))
    xsec=0.0d0
    allocate(delta(davstates))
    delta=0.0d0

!----------------------------------------------------------------------
! Calculation of the transition dipole moments between the ADC
! eigenstates and the Lanczos pseudostates
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(lvec(ndim))

    ! Open the Lanczos pseudo-spectrum file
    call freeunit(ilanc)
    open(ilanc,file='SCRATCH/lancstates_v',status='old',&
         access='sequential',form='unformatted')

    ! Loop over the Lanczos pseudo-states
    do alpha=1,nlanc
       
       ! Read the current Lanczos pseudo-state
       read(ilanc) k,lener(alpha),lvec
       
       ! < i | D_a | alpha >
       do a=1,3
          tdmil(a,alpha)=dot_product(travec_iv(:,a),lvec)
       enddo

       ! < f | D_a | alpha >
       do f=1,davstates
          do a=1,3
             tdmfl(a,f,alpha)=&
                  dot_product(travec_fv(:,a,f),lvec)
          enddo
       enddo

    enddo

    ! Close the Lanczos pseudo-spectrum file
    close(ilanc)

    ! Deallocate arrays
    deallocate(lvec)

!----------------------------------------------------------------------
! Calculation of the two-photon transition moments
!----------------------------------------------------------------------
    ! Loop over final states
    do f=1,davstates
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       
       ! Loop over pairs of dipole operator components
       do a=1,3
          do b=1,3
             
             ! Loop over the Lanczos pseudo-states
             do alpha=1,nlanc

                sabif(a,b,f)=sabif(a,b,f)&
                     + tdmil(a,alpha)*tdmfl(b,f,alpha)&
                     /(lener(alpha)-0.5d0*edavf(f))
                
                sabif(a,b,f)=sabif(a,b,f)&
                     + tdmil(b,alpha)*tdmfl(a,f,alpha)&
                     /(lener(alpha)-0.5d0*edavf(f))
                
             enddo
             
          enddo
       enddo

    enddo

    !  Contribution of the ground state to the two-photon transition
    ! moments (excited initial states only)
    ! N.B., E_0 = 0, hence the -E_f/2 denominators
    if (statenumber.gt.0) then
           
       ! Loop over final states
       do f=1,davstates
          if (statenumber.ne.0.and.f.eq.statenumber) cycle
          
          ! Loop over pairs of dipole operator components
          do a=1,3
             do b=1,3
                
                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmgsi(a)*tdmgsf(b,f)&
                     /(-0.5d0*edavf(f))
                
                sabif(a,b,f)=sabif(a,b,f)&
                     +tdmgsi(b)*tdmgsf(a,f)&
                     /(-0.5d0*edavf(f))

             enddo
          enddo
          
       enddo
       
    endif

!----------------------------------------------------------------------
! Calculation of the two-photon cross-sections for (parallel) linearly
! polarised light
!----------------------------------------------------------------------
    do f=1,davstates
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       do a=1,3
          do b=1,3
             delta(f)=delta(f)&
                  +2.0d0*sabif(a,a,f)*sabif(b,b,f)&
                  +2.0d0*sabif(a,b,f)*sabif(a,b,f)&
                  +2.0d0*sabif(a,b,f)*sabif(b,a,f)
          enddo
          ! Conversion to cross-sections
          xsec(f)=delta(f)*16.0d0*(pi**3)*((edavf(f)-e_init)/2.0d0)**2
          xsec(f)=xsec(f)/(c_au**2)
       enddo
    enddo
    xsec=xsec/30.0d0
    delta=delta/30.0d0

!----------------------------------------------------------------------
! Output the two-photon cross-sections
! N.B., we output the excitation energies relative to the initial
! state
!----------------------------------------------------------------------
    ! Log file
    write(ilog,'(/,2x,85a)') ('=',k=1,85)
    write(ilog,'(2x,a)') '     i     |     f     |     &
         E_if (eV)     |     delta_TP (au)     |     omega2 (au)'
    write(ilog,'(2x,85a)') ('-',k=1,85)
    do f=1,davstates
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       write(ilog,10) statenumber,'|',f,'|',&
            (edavf(f)-e_init)*eh2ev,'|',delta(f),'|',xsec(f)
    enddo
    write(ilog,'(2x,85a,/)') ('=',k=1,85)
    
    ! Stick spectrum file
    call freeunit(itpa)
    open(itpa,file='tpa.dat',form='formatted',status='unknown')
    do f=1,davstates
       if (statenumber.ne.0.and.f.eq.statenumber) cycle
       write(itpa,*) (edavf(f)-e_init)*eh2ev,xsec(f)
    enddo
    close(itpa)
    
10  format(6x,i2,5x,a1,4x,i2,5x,a1,3x,F10.5,6x,a1,3x,&
         F10.5,10x,a1,3x,F10.5)

    ! Gnuplot spectrum file
    call wrgnuplot(davstates,edavf-e_init,xsec)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sabif)
    deallocate(lener)
    deallocate(tdmil)
    deallocate(tdmfl)
    deallocate(xsec)
    deallocate(delta)

    return

  end subroutine tdm_tpa
  
!#######################################################################  

  subroutine wrgnuplot(nsta,ener,xsec)

    use constants
    use parameters
    use iomod
        
    implicit none

    integer                  :: nsta,iout,count,i,k
    real(d), dimension(nsta) :: ener,xsec
    real(d), parameter       :: tol=1e-10_d
    real(d)                  :: lb,ub
    character(len=400)       :: atmp

!----------------------------------------------------------------------
! Open the gnuplot file
!----------------------------------------------------------------------
    call freeunit(iout)
    open(iout,file='tpa.gnu',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the gnuplot file
!----------------------------------------------------------------------
    ! Broadening parameters
    write(iout,'(a,/)') 'fwhm=1.0'
    write(iout,'(a,/)') 'sigsq=(fwhm/2.35482)**2'

    count=0
    do i=1,nsta
       if (abs(xsec(i)).lt.tol) cycle
       count=count+1
       if (count.lt.10) then
          write(iout,'(a1,i1,a1,F11.5)') 'e',count,'=',ener(i)*eh2ev
          write(iout,'(a1,i1,a1,F16.10)') 'f',count,'=',xsec(i)
       else if (count.lt.100) then
          write(iout,'(a1,i2,a1,F11.5)') 'e',count,'=',ener(i)*eh2ev
          write(iout,'(a1,i2,a1,F16.10)') 'f',count,'=',xsec(i)
       else
          write(iout,'(a1,i3,a1,F11.5)') 'e',count,'=',ener(i)*eh2ev
          write(iout,'(a1,i3,a1,F16.10)') 'f',count,'=',xsec(i)
       endif
    enddo

    do i=1,count
       if (i.lt.10) then
          write(iout,'(a1,i1,a5,i1,a10,i1,a15)') &
               'g',i,'(x)=f',i,'*exp(-(x-e',i,')**2/(2*sigsq))'
       else if (i.lt.100) then
          write(iout,'(a1,i2,a5,i2,a10,i2,a15)') &
               'g',i,'(x)=f',i,'*exp(-(x-e',i,')**2/(2*sigsq))'
       else
          write(iout,'(a1,i3,a5,i3,a10,i3,a15)') &
               'g',i,'(x)=f',i,'*exp(-(x-e',i,')**2/(2*sigsq))'
       endif
    enddo
    
    atmp='f(x)='
    k=6
    do i=1,count
       if (i.lt.10) then
          write(atmp(k:k+5),'(a2,i1,a3)') '+g',i,'(x)'
          k=k+6
       else if (i.lt.100) then
          write(atmp(k:k+6),'(a2,i2,a3)') '+g',i,'(x)'
          k=k+7
       else
          write(atmp(k:k+7),'(a2,i3,a3)') '+g',i,'(x)'
          k=k+8
       endif
    enddo

    write(iout,'(a)') trim(atmp)
    
    lb=0.95*ener(1)*eh2ev
    ub=1.05*ener(nsta)*eh2ev
       
    write(iout,'(a12,F7.2,a1,F7.2,a1)') 'set xrange [',lb,':',ub,']'
       
    write(iout,'(a)') 'set size ratio 0.4'
    write(iout,'(a)') 'set samples 5000'       
    write(iout,'(a)') 'plot f(x) lt -1 lw 2 notitle'
    write(iout,'(a)') 'replot ''tpa.dat'' w i lt -1 lw 2 notitle'
    write(iout,'(a)') 'pause -1'

!----------------------------------------------------------------------
! Close the gnuplot file
!----------------------------------------------------------------------
    close(iout)

    return

  end subroutine wrgnuplot

!#######################################################################  
      
end module adc2tpamod
