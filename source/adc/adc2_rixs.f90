!#######################################################################
! Calculation of RIXS spectra using either the ADC(2) spectrum or the
! ADC(2) Lanczos pseudospectrum
!#######################################################################

module adc2rixsmod
  
  use channels

contains

!#######################################################################

  subroutine adc2_rixs(gam)

    use constants
    use parameters
    use adc_common
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
    real(d), dimension(:,:), allocatable :: rvec,travec2
    type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Make sure that symmetry is not being used as this is not currently
! supported in a RIXS calculation
!-----------------------------------------------------------------------
    call checksymm

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
    call initial_space_diag(time,kpq,ndim,ndims,noffd)

!-----------------------------------------------------------------------
! If requested, calculate the dipole moments for the initial states
!-----------------------------------------------------------------------
    if (ldipole.and.statenumber.gt.0) &
         call initial_space_dipole(ndim,ndims,kpq)

!-----------------------------------------------------------------------
! Transition moments from the ground state to the Davidson states
!-----------------------------------------------------------------------
    call initial_space_tdm(ener,rvec,ndim,mtm,tmvec,osc_str,kpq)
        
!-----------------------------------------------------------------------
! Output the results of initial space calculation
!-----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('*',i=1,70)
    write(ilog,'(2x,a)') &
         'Initial space ADC(2)-s excitation energies'
    write(ilog,'(70a)') ('*',i=1,70)
    
    itmp=1+nBas**2*4*nOcc**2
    call table2(ndim,davstates,ener(1:davstates),&
         rvec(:,1:davstates),tmvec(1:davstates),&
         osc_str(1:davstates),kpq,itmp,'i')
    
    write(ilog,'(/,70a,/)') ('*',i=1,70)

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
    call final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
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
!-----------------------------------------------------------------------
    call tdm_rixs(ndim,ndimf,ndimsf,travec2,e_init)

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

  end subroutine adc2_rixs

!#######################################################################

  subroutine final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,kpqf,travec,&
       vec_init,mtmf,noffd,noffdf,rvec,travec2)

    use constants
    use parameters
    use fspace

    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndimf,ndimsf
    integer*8                                 :: noffd,noffdf
    real(d), dimension(:), allocatable        :: travec,mtmf
    real(d), dimension(ndim)                  :: vec_init
    real(d), dimension(ndim,davstates)        :: rvec
    real(d), dimension(:,:), allocatable      :: travec2

    if (ldiagfinal) then
       call davidson_final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,&
            kpqf,travec,vec_init,mtmf,noffdf,rvec,travec2)
    else
       call lanczos_final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,&
            kpqf,travec,vec_init,mtmf,noffdf,rvec,travec2)
    endif
        
    return
    
  end subroutine final_space_diag_rixs
  
!#######################################################################

  subroutine lanczos_final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,kpqf,&
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
! Determine the initial Lanczos vectors
!-----------------------------------------------------------------------
    call guess_vecs_rixs(rvec,ndim,ndimf,ndimsf,kpq,kpqf,travec2)

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
    
  end subroutine lanczos_final_space_diag_rixs

!#######################################################################

  subroutine guess_vecs_rixs(rvec,ndim,ndimf,ndimsf,kpq,kpqf,travec2)
        
    use constants
    use parameters
    use misc
    use get_matrix_dipole
    use get_moment
    use fspace
    use guessvecs
    use iomod

    implicit none
        
    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: i,j,k,c,ndim,ndimf,&
                                                 ndimsf,error,&
                                                 ivecs,ilbl
    real(d), dimension(ndim,davstates)        :: rvec
    real(d), dimension(:,:), allocatable      :: travec2,initvecs
    real(d), dimension(:), allocatable        :: tau,work

!----------------------------------------------------------------------
! Allocate the travec2 array
!
! The first three columns of this array holds the F-vectors
! F_Ja = < Psi_J | Da | Psi_0 >, a=x,y,z, where the Psi_J
! are the ISs spanning the final space
!
! The remaining columns hold the products of the IS representation of
! the shifted dipole operator and the initial space vectors
!----------------------------------------------------------------------
    allocate(travec2(ndimf,3*(davstates+1)))

!-----------------------------------------------------------------------
! Calculate the contractions of the IS shifted dipole matrix and the
! valence space vectors
!-----------------------------------------------------------------------
    call dipole_ispace_contraction_all(rvec,travec2,ndim,ndimf,kpq,kpqf)

!-----------------------------------------------------------------------
! The initial Lanczos vectors will be formed from the contractions of
! the shifted dipole matrix with the valence states (ground + excited),
! i.e., the vectors held in the travec2 array
!
! IN THE ORTHOGONALISATION STEP, WE ARE MIXING THE DIPOLE-STATE
! PRODUCTS AMONGST THEMSELVES. IS THIS GOING TO PRODUCE VECTORS THAT
! DO NOT GENERATE IRREPS OF THE MOLECULAR POINT GROUP?
!-----------------------------------------------------------------------
    ! Set the block size
    lmain=3*(davstates+1)
        
    ! Copy the contents of the travec2 array
    allocate(initvecs(ndimf,lmain))
    initvecs=travec2
        
    ! Orthogonalisation of the shifted dipole matrix-state vector
    ! contractions via a QR factorisation
    allocate(tau(lmain))
    allocate(work(lmain))
    call dgeqrf(ndimf,lmain,initvecs,ndimf,tau,work,lmain,error)
    if (error.ne.0) then
       errmsg='dqerf failed in subroutine guess_vecs_rixs'
       call error_control
    endif
    call dorgqr(ndimf,lmain,lmain,initvecs,ndimf,tau,work,lmain,error)
    if (error.ne.0) then
       errmsg='dorgqr failed in subroutine guess_vecs_rixs'
       call error_control
    endif
    deallocate(tau)
    deallocate(work)
    
    ! Write the guess vectors to file
    call freeunit(ivecs)
    open(ivecs,file='SCRATCH/rixs_ivecs',form='unformatted',&
         status='unknown')
    write(ivecs) initvecs
    close(ivecs)
    
    deallocate(initvecs)
    
    return
    
  end subroutine guess_vecs_rixs

!#######################################################################

  subroutine davidson_final_space_diag_rixs(ndim,ndimf,ndimsf,kpq,kpqf,&
       travec,vec_init,mtmf,noffdf,rvec,travec2)

    use constants
    use parameters
    use fspace
    use get_matrix_dipole
    use adc_common
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
    allocate(travec2(ndimf,3*(davstates+1)))
    
!-----------------------------------------------------------------------
! Calculate travec: the product of the IS representation of the 
! shifted dipole operator and the initial state vector.
!
! This will be used later in the calculation of transition dipole
! matrix elements between the initial and final states.
!
! Really this should be done elsewhere...
!-----------------------------------------------------------------------
    call dipole_ispace_contraction_all(rvec,travec2,ndim,ndimf,&
         kpq,kpqf)
    
    return
    
  end subroutine davidson_final_space_diag_rixs
    
!#######################################################################

  subroutine dipole_ispace_contraction_all(rvec,travec2,ndim,ndimf,&
       kpq,kpqf)

    use constants
    use parameters
    use adc_common
    use misc
    use get_matrix_dipole
    use get_moment
    use fspace
    use guessvecs
    use iomod
        
    implicit none

    integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
    integer                                   :: ndim,ndimf,c,i,&
                                                 k,ilbl
    real(d), dimension(ndim,davstates)        :: rvec
    real(d), dimension(ndimf,3*(davstates+1)) :: travec2
    character(len=1), dimension(3)            :: acomp
    character(len=70)                         :: msg

    integer*8, dimension(3)                   :: nel_cv
    integer, dimension(3)                     :: nbuf_cv
    character(len=60)                         :: filename

!-----------------------------------------------------------------------
! Calculate the matrix F_Ja = < Psi_J | Da | Psi_0 >, where the Psi_J
! are the ISs spanning the final space and a=x,y,z
!-----------------------------------------------------------------------
    acomp=(/ 'x','y','z' /)
        
    ! Loop over the components of the dipole operator
    do c=1,3
       
       ! Set the dipole component
       dpl(:,:)=dpl_all(c,:,:)

       ! Calculate the vector F_Jc
       msg='Calculating the '//acomp(c)//'-component of the F-vector'
       write(ilog,'(/,2x,a)') trim(msg)
       call get_modifiedtm_adc2(ndimf,kpqf(:,:),travec2(:,c),0)
       
    enddo

!-----------------------------------------------------------------------
! Calculate the IS representation of the shifted dipole operator
!-----------------------------------------------------------------------
    ! Loop over the components of the dipole operator
    do c=1,3
           
       ! Set the dipole component
       dpl(:,:)=dpl_all(c,:,:)
       
       ! Calculate the IS representation of the dipole operator
       write(ilog,'(70a)') ('-',k=1,70)
       msg='Calculation of P_c . D_'//acomp(c)//' . P_v'

       write(ilog,'(2x,a)') trim(msg)
       write(ilog,'(70a)') ('-',k=1,70)
       filename='SCRATCH/dipole_cv_'//acomp(c)
       call get_adc2_dipole_improved_omp(ndimf,ndim,kpqf,kpq,&
            nbuf_cv(c),nel_cv(c),filename)

    enddo

!-----------------------------------------------------------------------
! Calculate the products of the IS representation of the shifted dipole
! operator and the initial space vectors
!-----------------------------------------------------------------------
    ! Loop over initial states
    do i=1,davstates
           
       ! Loop over the components of the dipole operator
       do c=1,3

          ! Calculate the matrix-vector product
          ilbl=3+(i-1)*3+c
          filename='SCRATCH/dipole_cv_'//acomp(c)
          call contract_dipole_state(filename,ndim,ndimf,&
               rvec(:,i),travec2(:,ilbl),nbuf_cv(c),nel_cv(c),'r')
          
       enddo

    enddo

    return
    
  end subroutine dipole_ispace_contraction_all

!#######################################################################

  subroutine tdm_rixs(ndim,ndimf,ndimsf,travec2,e_init)
        
    use constants
    use parameters
    use iomod
    use timingmod

    implicit none
    
    integer                                   :: ndim,ndimf,&
                                                 ndimsf,i,j,k,&
                                                 ivecf,irixs,&
                                                 idav,nfinal
    real(d)                                   :: e_init
    real(d), dimension(ndimf,3*(davstates+1)) :: travec2
    real(d), dimension(:,:), allocatable      :: tdm
    real(d), dimension(:), allocatable        :: vec
    real(d), dimension(:), allocatable        :: enerf
    real(d), dimension(davstates)             :: ener
    real(d)                                   :: tw1,tw2,tc1,tc2
    character(len=70)                         :: filename

!-----------------------------------------------------------------------
! Output what we are doing
!-----------------------------------------------------------------------
    write(ilog,'(/,70a)') ('-',i=1,70)
    write(ilog,'(2x,a,/,2x,a)') &
         'Calculating the ADC-to-CVS-ADC-Lanczos transition &
         dipoles and','writing the RIXS data file'
    write(ilog,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
      call times(tw1,tc1)

!-----------------------------------------------------------------------
! Set dimensions and filenames
!-----------------------------------------------------------------------
    if (ldiagfinal) then
       ! Davidson-RIXS
       nfinal=davstates_f
       filename=davname_f
    else
       ! Block Lanczos-RIXS
       nfinal=lancstates
       filename=lancname
    endif
        
!-----------------------------------------------------------------------
! For all valence states |Psi_i> (including the ground state),
! calculate the matrix elements <Psi_i | Da | chi_j>, where the chi_j
! are the final space states and a=x,y,z
!-----------------------------------------------------------------------
    ! Allocate arrays        
    allocate(vec(ndimf))
    allocate(tdm(3*(davstates+1),nfinal))
    allocate(enerf(nfinal))
        
    ! Open the final space vector file
    call freeunit(ivecf)
    open(ivecf,file=filename,status='old',access='sequential',&
         form='unformatted')
    
    ! Loop over the final space states
    do j=1,nfinal

       ! Read the jth final space eigenpair
       read(ivecf) k,enerf(j),vec

       ! Loop over the x, y, and z components of D for each
       ! of the valence states (ground + excited)
       do i=1,3*(davstates+1)
          tdm(i,j)=dot_product(travec2(:,i),vec)
       enddo
       
    enddo
        
    ! Close the final space vector file
    close(ivecf)

    ! Deallocate arrays
    deallocate(vec)

!-----------------------------------------------------------------------
! Read the initial space energies from file
!-----------------------------------------------------------------------
    ! Allocate arrays
    allocate(vec(ndim))

    ! Open the Davidson file
    call freeunit(idav)
    open(idav,file=davname,status='unknown',access='sequential',&
         form='unformatted')

    ! Read the energies
    do i=1,davstates
       read(idav) k,ener(i),vec
    enddo

    ! Close the Davidson file
    close(idav)
        
    ! Deallocate arrays
    deallocate(vec)

!-----------------------------------------------------------------------
! Write the RIXS data file
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Transition dipole matrix elements written &
         to rixs.dat'

    ! Open the RIXS data file
    call freeunit(irixs)
    filename='rixs.dat'
    open(irixs,file=filename,status='unknown',form='unformatted')

    ! Dimensions
    write(irixs) davstates+1
    write(irixs) nfinal
    
    ! Ground state energy
    write(irixs) ehf+e_mp2
        
    ! Energies of the valence-excited states
    do i=1,davstates
       write(irixs) ehf+e_mp2+ener(i)
    enddo

    ! Energies of the final space states
    do i=1,nfinal
       write(irixs) ehf+e_mp2+enerf(i)
    enddo

    ! Transition dipole matrix elements <Psi_i | Da | chi_j>,
    ! a=x,y,z between the valence states and the
    ! final space states
    do i=1,3*(davstates+1)
       do j=1,nfinal
          write(irixs) tdm(i,j)
       enddo
    enddo
    
    ! Close the RIXS data file
    close(irixs)

!-----------------------------------------------------------------------
! Finish timing and output the walltime taken
!-----------------------------------------------------------------------
      call times(tw2,tc2)
      write(ilog,'(/,2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'

    return

  end subroutine tdm_rixs
  
!#######################################################################
  
  end module adc2rixsmod
