  module adc2enermod

  contains

!#######################################################################

    subroutine adc2_ener
      
      use constants
      use parameters
      use select_fano
      use davmod
      use lancmod
      use fspace
      use get_moment
      use misc
      use fspace2
      use get_matrix
      use get_matrix_DIPOLE
      use propagate_prepare
      use guessvecs
      use channels
      use adc_ph, only: mp2

      implicit none
    
      integer, dimension(:,:), allocatable :: kpq
      integer                              :: ndim,ndims,i,itmp,ista
      integer*8                            :: noffd
      real(d), dimension(:), allocatable   :: ener,vec_init,mtm,tmvec,&
                                              osc_str
      real(d), dimension(:,:), allocatable :: rvec
      real(d)                              :: t1,t2

!-----------------------------------------------------------------------  
! Calculation of the MP2 correlation energy
!-----------------------------------------------------------------------  
      call MP2(e_mp2)

!-----------------------------------------------------------------------  
! Calculate guess initial space vectors from an ADC(1) calculation if 
! requested.
!
! N.B. adc1_guessvecs must be called BEFORE we select the ADC2 
! configurations, otherwise the kpq array will be overwritten with the 
! ADC1 configurations.
!-----------------------------------------------------------------------  
      if (ladc1guess) call adc1_guessvecs

!-----------------------------------------------------------------------
! Allocate kpq and select configurations
!-----------------------------------------------------------------------
      allocate(kpq(7,0:nBas**2*4*nOcc**2))
      kpq(:,:)=-1
      
      if (lcvs) then
         if (lfakeip) then
            ! CVS-IP-ADC(2)
            call select_atom_is_cvs_fakeip(kpq(:,:))
            call select_atom_d_cvs_fakeip(kpq(:,:),-1)
         else
            ! CVS-ADC(2)
            call select_atom_is_cvs(kpq(:,:))
            call select_atom_d_cvs(kpq(:,:),-1)
         endif
      else if (lfakeip) then
         ! IP-ADC(2)
         call select_atom_is_fakeip(kpq(:,:))
         call select_atom_d_fakeip(kpq(:,:),-1)
      else
         ! ADC(2)
         call select_atom_is(kpq(:,:))
         call select_atom_d(kpq(:,:),-1)
      endif
    
!-----------------------------------------------------------------------
! Output supspace dimensions
!-----------------------------------------------------------------------
      ndims=kpq(1,0)
      ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)

      write(ilog,*) 'ADC(2) INITIAL Space dim',ndim

      write(ilog,*) 'dimension of various INITIAL configuration spaces'
      write(ilog,*) '      1p1h       2p2h_1      2p2h_2      2p2h_3      2p2h_4i      2p2h_4ii'
      write(ilog,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0),kpq(5,0)
      write(ilog,*)

!-----------------------------------------------------------------------
! Initialise dipole moment matrix
!-----------------------------------------------------------------------
      if (tranmom2 .eq. 'x') then
         dpl(:,:)=x_dipole(:,:)
      elseif (tranmom2 .eq. 'y') then
         dpl(:,:)=y_dipole(:,:)
      elseif (tranmom2 .eq. 'z') then
         dpl(:,:)=z_dipole(:,:)
      end if
      
      CHECK_dip=nirrep2

!-----------------------------------------------------------------------
! Calculate and save the Hamiltonian matrix to file
!-----------------------------------------------------------------------
      write(ilog,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
      if (method.eq.-2) then
         ! ADC(2)-s
         if (lcvs) then
            call write_fspace_adc2_1_cvs(ndim,kpq(:,:),noffd,'i')
         else
            call write_fspace_adc2_1(ndim,kpq(:,:),noffd,'i')
         endif
      else if (method.eq.-3) then
         ! ADC(2)-x
         if (lcvs) then
            call write_fspace_adc2e_1_cvs(ndim,kpq(:,:),noffd,'i')
         else
            call write_fspace_adc2e_1(ndim,kpq(:,:),noffd,'i')
         endif
      endif

!-----------------------------------------------------------------------
! Block-Davidson diagonalisation of the Hamiltonian matrix
!-----------------------------------------------------------------------
      allocate(ener(davstates),rvec(ndim,davstates))
      allocate(vec_init(ndim))
      
      call master_eig(ndim,noffd,'i')
      
      call readdavvc(davstates,ener,rvec,'i',ndim)

!-----------------------------------------------------------------------
! Calculate TDMs from the ground state
!-----------------------------------------------------------------------    
      allocate(mtm(ndim),tmvec(davstates),osc_str(davstates))
    
      tmvec=0.0d0
      osc_str=0.0d0

      if (ltdm_gs2i.and..not.lfakeip) then
         call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:),1)
         do i=1,davstates
            tmvec(i)=tm(ndim,rvec(:,i),mtm(:))
            osc_str(i)=2.0d0/3.0d0*ener(i)*tmvec(i)**2
         enddo
      endif
      
      itmp=1+nBas**2*4*nOcc**2
      call table2(ndim,davstates,ener(1:davstates),&
           rvec(:,1:davstates),tmvec(1:davstates),&
           osc_str(1:davstates),kpq,itmp,'i')
      
      return
      
    end subroutine adc2_ener

!#######################################################################
    
  end module adc2enermod
