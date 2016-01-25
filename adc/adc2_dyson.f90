  module adc2dysonmod

    use channels

    contains

!#######################################################################

      subroutine adc2_dyson(gam)

        use constants
        use parameters
        use fspace
        use misc
        use guessvecs
        use mp2
        use import_gamess

        implicit none

        integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
        integer                              :: ndim,ndims,ndimsf,&
                                                ndimf,ndimd
        real(d)                              :: e0,einit,time
        real(d), dimension(:), allocatable   :: vec_init
        type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic
!-----------------------------------------------------------------------
        call mp2_master(e0)

!-----------------------------------------------------------------------  
! Initial space diagonalisation
!-----------------------------------------------------------------------  
        if (statenumber.gt.0) call calc_initial_state(kpq,ndim,ndims,&
             vec_init,einit,time)

!-----------------------------------------------------------------------
! Final space diagonalisation (ionized states)
!-----------------------------------------------------------------------
        call calc_final_states(kpqf,ndimf,ndimsf)

!-----------------------------------------------------------------------
! Calculation of the expansion coefficients for the Dyson orbitals
! in the MO basis
!-----------------------------------------------------------------------
        call dysorb(kpqf,ndimf,ndimsf,kpq,ndim,ndims,vec_init,einit,gam)

        return

      end subroutine adc2_dyson

!#######################################################################

      subroutine calc_initial_state(kpq,ndim,ndims,vec_init,einit,time)

        use constants
        use parameters
        use guessvecs
        use fspace
        use davmod
        use iomod, only: freeunit
        use misc

        implicit none

        integer, dimension(:,:), allocatable :: kpq
        integer                              :: ndim,ndims,ivec,i,itmp
        integer*8                            :: noffd
        real(d)                              :: einit,time
        real(d), dimension(:), allocatable   :: vec_init
        character(len=120)                   :: msg

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(kpq(7,0:nBas**2*4*nOcc**2))

!-----------------------------------------------------------------------
! Determine the initial subspace
!-----------------------------------------------------------------------
        kpq(:,:)=-1
        call select_atom_is(kpq(:,:))
        call select_atom_d(kpq(:,:),-1)

        ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
        ndims=kpq(1,0)
        
        write(ilog,*) 'ADC(2) INITIAL Space dim',ndim
        write(ilog,*) 'dimension of various INITIAL configuration spaces'
        write(ilog,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0)

!-----------------------------------------------------------------------  
! Calculate guess initial space vectors from an ADC(1) calculation if 
! requested.
!-----------------------------------------------------------------------  
        if (ladc1guess) call adc1_guessvecs

!-----------------------------------------------------------------------
! Write the initial space ADC(2)-s Hamiltonian to disk
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
! Block-Davidson diagonalisation
!-----------------------------------------------------------------------
        call master_eig(ndim,noffd,'i')

!-----------------------------------------------------------------------
! Output the results of initial space calculation
!-----------------------------------------------------------------------
        if (method.eq.2) then
           msg='Initial space ADC(2)-s excitation energies'
        else if (method.eq.3) then
           msg='Initial space ADC(2)-x excitation energies'
        endif

        write(ilog,'(/,70a)') ('*',i=1,70)
        write(ilog,'(2x,a)') trim(msg)
        write(ilog,'(70a)') ('*',i=1,70)

        itmp=1+nBas**2*4*nOcc**2
        call wrstateinfo_neutral(ndim,kpq,itmp,davname,davstates)

        write(ilog,'(/,70a,/)') ('*',i=1,70)

!-----------------------------------------------------------------------
! Load the initial state vector into memory
!-----------------------------------------------------------------------
        call freeunit(ivec)
        open(unit=ivec,file=davname,status='unknown',&
             access='sequential',form='unformatted')
        
        allocate(vec_init(ndim))

        do i=1,statenumber-1
           read(ivec)
        enddo

        read(ivec) i,einit,vec_init

        close(ivec)

        return

      end subroutine calc_initial_state

!#######################################################################

      subroutine calc_final_states(kpqf,ndimf,ndimsf)

        use constants
        use parameters

        implicit none
        
        integer, dimension(:,:), allocatable :: kpqf
        integer                              :: ndimf,ndimsf

!-----------------------------------------------------------------------
! Determine the final subspace
!-----------------------------------------------------------------------
        call get_subspace_final(kpqf,ndimf,ndimsf)

!-----------------------------------------------------------------------
! Construct and diagonalise the IP-ADC(2)-s or IP-ADC(2)-x Hamiltonian 
! matrix
!-----------------------------------------------------------------------
        if (dysdiag.eq.1) then
           call diag_ipadc_full(kpqf,ndimf,ndimsf)
        else if (dysdiag.eq.2) then
           call diag_ipadc_dav(kpqf,ndimf,ndimsf)
        endif

        return

      end subroutine calc_final_states

!#######################################################################

      subroutine get_subspace_final(kpqf,ndimf,ndimsf)

        use constants
        use parameters
        use select_fano

        implicit none
        
        integer, dimension(:,:), allocatable :: kpqf
        integer                              :: ndimf,ndimsf

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(kpqf(7,0:nBas**2*4*nOcc**2))

!-----------------------------------------------------------------------
! Determine the final subspace
!-----------------------------------------------------------------------
        kpqf(:,:)=-1

        if (lcvsfinal) then
           call select_atom_isf_fakeip_cvs(kpqf(:,:))
           call select_atom_df_fakeip_cvs(kpqf(:,:),-1)
        else
           call select_atom_isf_fakeip(kpqf(:,:))
           call select_atom_df_fakeip(kpqf(:,:),-1)
        endif

        ndimsf=kpqf(1,0)
        ndimf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)

!-----------------------------------------------------------------------
! Output the final subspace information
!-----------------------------------------------------------------------
        write(ilog,*) 'ADC(2) FINAL Space dim',ndimf
        write(ilog,*) 'dimension of various FINAL configuration spaces'
        write(ilog,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
        
        return

      end subroutine get_subspace_final

!#######################################################################

      subroutine diag_ipadc_full(kpqf,ndimf,ndimsf)

        use parameters
        use fspace
        use iomod, only: freeunit

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,i,&
                                                     iout,itmp
        real(d), dimension(:,:), allocatable      :: eigvecf
        real(d), dimension(:), allocatable        :: eigvalf
        character(len=120)                        :: msg
        character(len=36)                         :: filename

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(eigvecf(ndimf,ndimf))
        allocate(eigvalf(ndimf))

!-----------------------------------------------------------------------
! Calculate and diagonalise the IP-ADC(2) Hamiltonian matrix
!-----------------------------------------------------------------------
        if (method_f.eq.2) then
           msg='Full diagonalisation of the IP-ADC(2)-s Hamiltonian &
                matrix'
        else if (method_f.eq.3) then
           msg='Full diagonalisation of the IP-ADC(2)-x Hamiltonian &
                matrix'
        endif
        write(ilog,'(/,2x,a)') trim(msg)

        if (method_f.eq.2) then
           ! IP-ADC(2)-s
           call get_fspace_adc2_direct(ndimf,kpqf,eigvecf,eigvalf)
        else if (method_f.eq.3) then
           ! IP-ADC(2)-x
           call get_fspace_adc2e_direct(ndimf,kpqf,eigvecf,eigvalf)
        endif

!-----------------------------------------------------------------------
! Save the eigenpairs to file
!-----------------------------------------------------------------------
        call freeunit(iout)

        open(unit=iout,file='ipadc_vecs',status='unknown',&
             access='sequential',form='unformatted')

        do i=1,ndimf
           write(iout) i,eigvalf(i),eigvecf(:,i)
        enddo

        close(iout)

!-----------------------------------------------------------------------
! Output the cation state energies and dominant configurations
!-----------------------------------------------------------------------
        itmp=1+nBas**2*4*nOcc**2
        filename='ipadc_vecs'
        call wrstateinfo_cation(ndimf,kpqf,itmp,filename,ndimf)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(eigvecf)
        deallocate(eigvalf)

        return

      end subroutine diag_ipadc_full

!#######################################################################

      subroutine diag_ipadc_dav(kpqf,ndimf,ndimsf)

        use constants        
        use parameters
        use fspace
        use iomod, only: freeunit
        use fspace
        use davmod
        use guessvecs

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,itmp
        character(len=120)                        :: msg
        integer*8                                 :: noffdf

!-----------------------------------------------------------------------        
! If requested, determine the Davidson guess vectors by diagonalising 
! the ADC(1) Hamiltonian matrix
!-----------------------------------------------------------------------        
        if (ladc1guess_f) call adc1_guessvecs_final

!-----------------------------------------------------------------------
! Write the final space ADC(2) Hamiltonian to disk
!-----------------------------------------------------------------------
        write(ilog,*) 'Saving complete FINAL SPACE IP-ADC2 matrix in &
             file'

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
! Block-Davidson diagonalisation in the final space
!-----------------------------------------------------------------------
        call master_eig(ndimf,noffdf,'f')

!-----------------------------------------------------------------------
! Output the cation state energies and dominant configurations
!-----------------------------------------------------------------------
        itmp=1+nBas**2*4*nOcc**2
        call wrstateinfo_cation(ndimf,kpqf,itmp,davname_f,davstates_f)

        return

      end subroutine diag_ipadc_dav

!#######################################################################

      subroutine dysorb(kpqf,ndimf,ndimsf,kpq,ndim,ndims,vec_init,&
           einit,gam)

        use constants
        use parameters
        use dyson_calc
        use dyson_io
        use iomod, only: freeunit
        use import_gamess

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf,kpq
        integer                                   :: ndimf,ndimsf,ndim,&
                                                     ndims,i,j,a,b,n,&
                                                     inorm,icoeff,ivec,&
                                                     itmp,nsta
        real(d), dimension(ndim)                  :: vec_init
        real(d), dimension(ndimf)                 :: eigvec
        real(d), dimension(:,:), allocatable      :: rhogs2,rmat,smat
        real(d), dimension(:), allocatable        :: dyscoeff
        real(d)                                   :: einit,ei,eigval,&
                                                     norm,de
        character(len=36)                         :: vecfile
        logical                                   :: ldir
        type(gam_structure)                       :: gam
        

        write(ilog,'(/,2x,a,/)') 'Calculating the Dyson orbital &
             expansion coefficients'

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        ! Dyson orbital expansion coefficients
        allocate(dyscoeff(nbas))

        ! 2nd-order correction to the ground state density matrix
        allocate(rhogs2(nbas,nbas))

        ! R-matrix
        allocate(rmat(nbas,nbas))

        ! S-matrix
        allocate(smat(nbas,nbas))

!-----------------------------------------------------------------------
! Pre-calculate intermediate terms
!-----------------------------------------------------------------------
        call dyson_precalc(rhogs2,rmat,smat,vec_init,ndim,ndims,kpq)

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
        ! IP-ADC eigenpairs
        if (dysdiag.eq.1) then
           vecfile='ipadc_vecs'
        else
           vecfile=davname_f
        endif
        call freeunit(ivec)
        open(unit=ivec,file=vecfile,status='unknown',&
             access='sequential',form='unformatted')

        ! Dyson orbital norms
        call freeunit(inorm)
        open(inorm,file='dyson_norms',form='formatted',&
             status='unknown')

        ! Dyson orbital expansion coefficients
        call freeunit(icoeff)
        open(icoeff,file='dyson_coeffs',form='formatted',&
             status='unknown')

!-----------------------------------------------------------------------
! If the Dyson orbitals are to be written to external input files, then
! assemble pre-requesite arrays and create output directories
!-----------------------------------------------------------------------
        if (sum(dysout).gt.0) call pre_dysout(gam)

!-----------------------------------------------------------------------        
! Calculation of the Dyson orbital expansion coefficients for all
! IP-ADC(2) states
!-----------------------------------------------------------------------
        if (dysdiag.eq.1) then
           nsta=ndimf
        else if (dysdiag.eq.2) then
           nsta=davstates_f
        endif

        if (statenumber.eq.0) then
           ei=0.0d0
        else
           ei=einit
        endif

        do n=1,nsta

           dyscoeff=0.0d0

           ! Read the next eigenpair from file
           read(ivec) itmp,eigval,eigvec

           ! Exit if the next state lies above the energetic limit
           if (eigval.gt.dyslim) exit

           ! Calculate the current set of coefficients
           if (statenumber.eq.0) then
              call dyscoeff_gs(rhogs2,dyscoeff,kpqf,eigvec,ndimf,ndimsf)
           else
              call dyscoeff_exci(rhogs2,dyscoeff,kpq,kpqf,vec_init,&
                   eigvec,ndim,ndims,ndimf,ndimsf,rmat,smat)
           endif

           ! Output the Dyson orbital norm and coefficients
           norm=sqrt(dot_product(dyscoeff,dyscoeff))
           de=(eigval-ei)*eh2ev
           write(inorm,'(2(2x,F13.7))') de,norm
           write(icoeff,'(/,a,x,i4)') 'Cation state',n
           do i=1,nbas
              write(icoeff,'(i4,F13.7)') i,dyscoeff(i)
           enddo

           ! Output the *normalised* Dyson orbital to the molden file
           if (dysout(1).eq.1) call wrmolden(dyscoeff,gam,n)
           
           ! Output the Dyson orbital to an ezdyson input file
           if (dysout(2).eq.1) call wrezdyson(dyscoeff,gam,n,de)

        enddo

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
        close(ivec)
        close(inorm)
        close(icoeff)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(dyscoeff)
        deallocate(rhogs2)
        deallocate(rmat)
        deallocate(smat)

        return

      end subroutine dysorb

!#######################################################################

  end module adc2dysonmod
