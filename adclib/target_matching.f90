!#######################################################################
! targetmatching: Subroutines to choose the initial ADC state via the
!                 criterion of maximum overlap with a user-specified
!                 target CI state. The following approximations are 
!                 made:
!
!                 (1) The ground state is taken through zeroth-order 
!                     in perturbation theory (i.e., we take it to be
!                     the Hartee-Fock state).
!
!                 (2) Only those Slater determinants with coefficients
!                     C_k > detthrsh are included in the approximate
!                     expansion of the ADC states.
!
!                 Hopefully this is sufficient for the mapping of
!                 states...
!
!                 This module makes heavy use of structures and 
!                 routines taken from v3 of the superdyson code of 
!                 S. Patchkovskii
!#######################################################################

  module targetmatching
    
    use constants
    use parameters
    use import_gamess

    implicit none

    save
    
    ! Parameters concerning the Slater determinant expansions of the
    ! ADC and target states
    integer                                   :: maxnsd,nsd_targ,&
                                                 nbas_targ
    integer                                   :: nao_adc,nao_targ,&
                                                 nel
    integer, dimension(:), allocatable        :: nsd_adc
    integer, dimension(:,:,:), allocatable    :: onv_adc
    integer, dimension(:,:), allocatable      :: onv_targ
    real(d), dimension(:,:), allocatable      :: c_adc
    real(d), dimension(:), allocatable        :: c_targ,s2_adc,&
                                                 norm_adc,overlap
    real(d)                                   :: s2_targ,norm_targ
    real(d), dimension(:,:), allocatable      :: sao,smo,ao2mo_adc,&
                                                 ao2mo_targ,smk

    ! Target state gamess structure
    type(gam_structure) :: gamess_target

    ! Useful common prefactors
    real(d), parameter :: invsqrt2=1.0d0/sqrt(2.0d0)
    real(d), parameter :: invsqrt12=1.0d0/sqrt(12.0d0)

  contains

!#######################################################################

    subroutine target_master(kpq,ndim,gamess_internal)

      use channels
      use parameters
      use constants
      use import_gamess

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,k
      type(gam_structure)                       :: gamess_internal

!-----------------------------------------------------------------------
! Output some information about where we are and what we're doing
!-----------------------------------------------------------------------
      write(ilog,'(/,70a)') ('-',k=1,70)
      write(ilog,'(a)') ' Determining the initial ADC state via &
           matching to a target CI state'
      write(ilog,'(70a,/)') ('-',k=1,70)

      write(ilog,'(a,/)') ' Approximate ADC states formed by taking &
           | Psi_0 > = | Psi_HF >'

      write(ilog,'(a,F7.5,/)') ' Determinant threshold: ',detthrsh

      write(ilog,'(a,F7.5,/)') ' Overlap threshold:     ',ovrthrsh

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(nsd_adc(davstates))
      allocate(s2_adc(davstates))
      allocate(norm_adc(davstates))
      allocate(overlap(davstates))

      nel=2*nocc
      allocate(smk(nel,nel))

!-----------------------------------------------------------------------
! Load the target state orbitals and determinant list
!-----------------------------------------------------------------------
      call load_target

!-----------------------------------------------------------------------
! Determine the maximum number of Slater determinants entering into
! the truncated zeroth-order expansions of the ADC states
!-----------------------------------------------------------------------
      call getdim_sd(kpq,ndim)

!-----------------------------------------------------------------------
! Construct the truncated list of ON vectors for the zeroth-order
! expansions of the ADC states and save the corresponding Slater
! determinant coefficients
!-----------------------------------------------------------------------
      call fill_onv_adc(kpq,ndim)

!-----------------------------------------------------------------------
! Calculate norms and S^2 expectation values for the Slater 
! determinant expansions (target and ~ADC)
!-----------------------------------------------------------------------
      call checkwf

!-----------------------------------------------------------------------
! Calculate the overlaps of the approximate ADC states with the target
! state
!-----------------------------------------------------------------------
      call calc_overlaps(gamess_internal)

!-----------------------------------------------------------------------
! Output information
!-----------------------------------------------------------------------
      call write_table

!-----------------------------------------------------------------------
! Select the initial state
!-----------------------------------------------------------------------
      call select_istate

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      call finalise

      write(ilog,'(/,70a)') ('-',k=1,70)

      return

    end subroutine target_master

!#######################################################################

    subroutine load_target

      use channels
      use parameters
      use constants

      implicit none

!-----------------------------------------------------------------------
! Load the target orbitals
!-----------------------------------------------------------------------
      write(ilog,'(2(x,a),/)') 'Loading the target MOs from',&
           trim(mofile)
      call gamess_load_orbitals(file=trim(mofile),&
           structure=gamess_target)

!-----------------------------------------------------------------------
! Read the target determinant list
!-----------------------------------------------------------------------
      write(ilog,'(/,2(x,a))') 'Reading the target state determinant &
           list from',trim(detfile)
      call rddetlist

      return

    end subroutine load_target

!#######################################################################

    subroutine rddetlist

      use iomod, only: freeunit
      use parsemod
      
      implicit none
      
      integer :: idet,ndet,i

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      call freeunit(idet)
      open(idet,file=detfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the no. determinants and allocate arrays
!-----------------------------------------------------------------------
      nsd_targ=0
10    call rdinp(idet)
      if (.not.lend) then
         nsd_targ=nsd_targ+1
         goto 10
      endif

      nbas_targ=gamess_target%nvectors

      allocate(onv_targ(nbas_targ,nsd_targ))
      allocate(c_targ(nsd_targ))

!-----------------------------------------------------------------------
! Read the target state determinants
!-----------------------------------------------------------------------
      rewind(idet)
      do i=1,nsd_targ
         read(idet,*) c_targ(i), onv_targ(:,i)
      enddo

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(idet)

      return

    end subroutine rddetlist

!#######################################################################
! getdim_sd: determines the number of Slater determinants entering
!            into the zeroth-order expansion of the Davidson states
!            with coefficients above the user set threshold
!#######################################################################

    subroutine getdim_sd(kpq,ndim)
      
      use channels
      use parameters
      use constants
      use iomod, only: freeunit

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,i,idav,itmp
      real(d), dimension(:), allocatable        :: coeff
      real(d)                                   :: ftmp

!-----------------------------------------------------------------------
! Open the Davidson vector file
!-----------------------------------------------------------------------  
      call freeunit(idav)
      open(idav,file=davname,status='old',access='sequential',&
           form='unformatted')

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(coeff(ndim))

!-----------------------------------------------------------------------
! For each Davidson state, calculate the number of Slater determinants 
! with absolute coefficient values above the threshold for 
! inclusion (dettrsh)
!-----------------------------------------------------------------------
      maxnsd=0

      do i=1,davstates

         ! Read the next ADC state vector from file
         read(idav) itmp,ftmp,coeff
         
         ! Determine the no. contributing Slater determinants
         call getnsd(i,kpq,coeff,ndim)

         ! Update maxnsd
         if (nsd_adc(i).gt.maxnsd) maxnsd=nsd_adc(i)         

      enddo
      
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(coeff)

!-----------------------------------------------------------------------  
! Close the Davidson vector file
!-----------------------------------------------------------------------  
      close(idav)

      return

    end subroutine getdim_sd

!#######################################################################

    subroutine getnsd(indx,kpq,coeff,ndim)
      
      use channels
      use parameters
      use constants

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: indx,ndim,n,count,&
                                                   nI,nII
      real(d), dimension(ndim)                  :: coeff
      real(d)                                   :: coeffsd

      nsd_adc(indx)=0

!-----------------------------------------------------------------------
! 1h1p ISs
!-----------------------------------------------------------------------
      do n=1,kpq(1,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)
      do n=count+1,count+kpq(2,0)
         coeffsd=coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+1
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a|=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)
      do n=count+1,count+kpq(3,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i|=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)
      do n=count+1,count+kpq(4,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a|=b, i|=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
      
      do n=count+1,count+kpq(5,0)
         
         nI=n
         nII=n+kpq(5,0)

         ! Spin cases I and II
         coeffsd=0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2
         coeffsd=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2

         ! Spin case II only
         coeffsd=2.0d0*invsqrt12*coeff(nII)
         if (abs(coeffsd).gt.detthrsh) nsd_adc(indx)=nsd_adc(indx)+2

      enddo

      return

    end subroutine getnsd

!#######################################################################

    subroutine fill_onv_adc(kpq,ndim)

      use channels
      use parameters
      use constants
      use iomod, only: freeunit

      implicit none
      
      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,s,idav,itmp
      real(d), dimension(:), allocatable        :: coeff
      real(d)                                   :: ftmp

!-----------------------------------------------------------------------
! Open the Davidson vector file
!-----------------------------------------------------------------------  
      call freeunit(idav)
      open(idav,file=davname,status='old',access='sequential',&
           form='unformatted')

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(coeff(ndim))
      allocate(onv_adc(nbas,maxnsd,davstates))
      allocate(c_adc(maxnsd,davstates))

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
      onv_adc(1:nocc,:,:)=2
      onv_adc(nocc+1:nbas,:,:)=0

!-----------------------------------------------------------------------
! Fill in the ON vector for each ADC state
!-----------------------------------------------------------------------
      do s=1,davstates
         
         ! Read the next state vector from file
         read(idav) itmp,ftmp,coeff

         ! Fill in the next ON vector
         call fill_onv_adc_1state(s,kpq,coeff,ndim)

      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(coeff)

!-----------------------------------------------------------------------  
! Close the Davidson vector file
!-----------------------------------------------------------------------  
      close(idav)

      return

    end subroutine fill_onv_adc

!#######################################################################

    subroutine fill_onv_adc_1state(sindx,kpq,coeff,ndim)

      use channels
      use parameters
      use constants

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: maxnsd,ndim,n,&
                                                   count,i,j,a,b,&
                                                   ksd,nI,nII,sindx
      real(d), dimension(ndim)                  :: coeff
      real(d)                                   :: ftmp

!-----------------------------------------------------------------------
! IS index mappings
!-----------------------------------------------------------------------
!      i<->kpq(3,n)
!      j<->kpq(4,n)
!      a<->kpq(5,n)
!      b<->kpq(6,n)
!-----------------------------------------------------------------------
! spin-to-ONV element mappings
!-----------------------------------------------------------------------
!      alpha       <-> +1
!      beta        <-> -1
!      doubly occ. <-> +2
!      unocc.      <->  0
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Initialisation of the Slater determinant counter
!-----------------------------------------------------------------------
      ksd=0

!-----------------------------------------------------------------------
! 1h1p ISs
!-----------------------------------------------------------------------
      do n=1,kpq(1,0)
         
         ! Orbital indices
         i=kpq(3,n)
         a=kpq(5,n)

         ! ON vectors and Slater determinant coefficients
         ftmp=abs(invsqrt2*coeff(n))
         if (ftmp.ge.detthrsh) then

            ! i_beta -> a_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=+1
            onv_adc(a,ksd,sindx)=-1
            c_adc(ksd,sindx)=invsqrt2*coeff(n)
            
            ! i_alpha -> a_alpha
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=-1
            onv_adc(a,ksd,sindx)=+1
            c_adc(ksd,sindx)=invsqrt2*coeff(n)

         endif

      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)

      do n=count+1,count+kpq(2,0)

         ! Orbital indices
         i=kpq(3,n)
         a=kpq(5,n)

         ! ON vectors and Slater determinant coefficients
         ftmp=abs(coeff(n))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, i_beta -> a_alpha, a_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=0
            onv_adc(a,ksd,sindx)=2
            c_adc(ksd,sindx)=coeff(n)
         endif

      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a|=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)

      do n=count+1,count+kpq(3,0)

         ! Orbital indices
         i=kpq(3,n)
         a=kpq(5,n)
         b=kpq(6,n)

         ! ON vectors and Slater determinant coefficients
         ftmp=abs(invsqrt2*coeff(n))
         if (ftmp.ge.detthrsh) then

            ! i_alpha, i_beta -> a_alpha, b_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=0
            onv_adc(a,ksd,sindx)=+1
            onv_adc(b,ksd,sindx)=-1
            c_adc(ksd,sindx)=invsqrt2*coeff(n)

            ! i_alpha, i_beta -> a_beta, b_alpha
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=0
            onv_adc(a,ksd,sindx)=-1
            onv_adc(b,ksd,sindx)=+1
            c_adc(ksd,sindx)=invsqrt2*coeff(n)

         endif
         
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i|=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)

      do n=count+1,count+kpq(4,0)

         ! Orbital indices
         i=kpq(3,n)
         j=kpq(4,n)
         a=kpq(5,n)

         ! ON vectors and Slater determinant coefficients
         ftmp=abs(invsqrt2*coeff(n))
         if (ftmp.ge.detthrsh) then

            ! i_alpha, j_beta -> a_alpha, a_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=-1
            onv_adc(j,ksd,sindx)=+1
            onv_adc(a,ksd,sindx)=2
            c_adc(ksd,sindx)=invsqrt2*coeff(n)

            ! i_beta, j_alpha -> a_alpha, a_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=+1
            onv_adc(j,ksd,sindx)=-1
            onv_adc(a,ksd,sindx)=2
            c_adc(ksd,sindx)=invsqrt2*coeff(n)

         endif

      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a|=b, i|=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)

      do n=count+1,count+kpq(5,0)

         nI=n
         nII=n+kpq(5,0)

         ! Orbital indices
         i=kpq(3,n)
         j=kpq(4,n)
         a=kpq(5,n)
         b=kpq(6,n)

         !--------------------------------------------------------------
         ! (1) Determinants with contributions from spin cases I and II
         !--------------------------------------------------------------
         ftmp=abs(0.5d0*coeff(nI)+invsqrt12*coeff(nII))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, j_beta -> a_alpha, b_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=-1
            onv_adc(j,ksd,sindx)=+1
            onv_adc(a,ksd,sindx)=+1
            onv_adc(b,ksd,sindx)=-1
            c_adc(ksd,sindx)=0.5d0*coeff(nI)+invsqrt12*coeff(nII)

            ! i_beta, j_alpha -> a_beta, b_alpha
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=+1
            onv_adc(j,ksd,sindx)=-1
            onv_adc(a,ksd,sindx)=-1
            onv_adc(b,ksd,sindx)=+1
            c_adc(ksd,sindx)=0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         endif

         ftmp=abs(-0.5d0*coeff(nI)+invsqrt12*coeff(nII))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, j_beta -> a_beta, b_alpha
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=-1
            onv_adc(j,ksd,sindx)=+1
            onv_adc(a,ksd,sindx)=-1
            onv_adc(b,ksd,sindx)=+1
            c_adc(ksd,sindx)=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)

            ! i_beta, j_alpha -> a_alpha, b_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=+1
            onv_adc(j,ksd,sindx)=-1
            onv_adc(a,ksd,sindx)=+1
            onv_adc(b,ksd,sindx)=-1
            c_adc(ksd,sindx)=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         endif

         !--------------------------------------------------------------
         ! (2) Determinants with contributions from spin case II only
         !--------------------------------------------------------------
         ftmp=abs(2.0d0*invsqrt12*coeff(nII))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, j_alpha -> a_alpha, b_alpha
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=-1
            onv_adc(j,ksd,sindx)=-1
            onv_adc(a,ksd,sindx)=+1
            onv_adc(b,ksd,sindx)=+1
            c_adc(ksd,sindx)=2.0d0*invsqrt12*coeff(nII)

            ! i_beta, j_beta -> a_beta, b_beta
            ksd=ksd+1
            onv_adc(i,ksd,sindx)=+1
            onv_adc(j,ksd,sindx)=+1
            onv_adc(a,ksd,sindx)=-1
            onv_adc(b,ksd,sindx)=-1
            c_adc(ksd,sindx)=2.0d0*invsqrt12*coeff(nII)
         endif

      enddo

      return

    end subroutine fill_onv_adc_1state

!#######################################################################
! checkwf: Calculates the norms and S^2 values for the truncated
!          Slater determinant expansions.
!          The code to do so is taken from v3 of the superdyson code.
!#######################################################################

    subroutine checkwf

      use channels
      use parameters
      use constants

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Calculate the expectation values of S^2
!-----------------------------------------------------------------------
      ! Target state
      call calc_s2(c_targ(:),onv_targ(:,:),s2_targ,nsd_targ,&
              nsd_targ,nbas_targ)

      ! Approximate ADC states
      do i=1,davstates
         call calc_s2(c_adc(:,i),onv_adc(:,:,i),s2_adc(i),maxnsd,&
              nsd_adc(i),nbas)
      enddo

!-----------------------------------------------------------------------
! Calculate the norms
!-----------------------------------------------------------------------
      ! Target state
      norm_targ=sqrt(dot_product(c_targ(:),c_targ(:)))

      ! Approximate ADC states
      do i=1,davstates
         norm_adc(i)=sqrt(dot_product(c_adc(:,i),c_adc(:,i)))
      enddo
      
      return

    end subroutine checkwf

!#######################################################################
! calc_s2: Adapted (stolen) from the superdyson code.
!          Calculates S^2 using the ladder operator representation
!
!          S^2 = 0.5 * S_+ S_- + 0.5 * S_- S_+ + S_z^2
!
!#######################################################################

    subroutine calc_s2(c,det,s2,maxdet,ndet,norb)

      use channels
      use parameters

      implicit none
      
      integer                         :: maxdet,id,ndet,norb,is,io,jo
      integer, dimension(norb,maxdet) :: det
      integer, dimension(ndet)        :: par
      integer, dimension(norb)        :: td
      real(d)                         :: s2
      real(d), dimension(maxdet)      :: c

!-----------------------------------------------------------------------
! Calculate parity factor connecting amplitudes of the same determinant,
! with orbitals in the "textbook" (alpha and beta factors intermixed, 
! spatial orbitals in the increasing order) and the alpha/beta (all 
! alpha before all beta) orders.
!-----------------------------------------------------------------------
      get_parity: do id=1,ndet
         par(id) = parity(det(:,id))
      enddo get_parity

!-----------------------------------------------------------------------
! Calculate S^2
!-----------------------------------------------------------------------
      s2=0.0d0
      scan_determinants: do id=1,ndet
         orbital_i: do io=1,norb
            if (abs(det(io,id))/=1) cycle orbital_i
            !
            !  Only partially open subshells could contribute to s^2; closed 
            !  subshells give identical zero upon action on the total spin operator
            !
            !  Diagonal part of S+ S- + S- S+ ; determinant does not change
            !
            s2 = s2 + 0.5d0 * c(id)**2
            !
            orbital_j: do jo=1,norb
               if (abs(det(jo,id))/=1) cycle orbital_j
               !
               !  S_z^2; determinant does not change
               !
               s2 = s2 + 0.25d0 * det(io,id) * det(jo,id) * c(id)**2 
               !
               !  Off-diagonal part of S+ S- + S- S+; determinant changes
               !
               if (det(io,id)/=det(jo,id)) then
                  td = det(:,id) ; td(io) = -td(io) ; td(jo) = -td(jo)
                  find_match: do is=1,ndet
                     if (any(det(:,is)/=td)) cycle find_match
                     !  S_+ S_- have interchanged the location of alpha and beta spins; to bring 
                     !  us back into the canonical alpha/beta order, we have to swap the spin-
                     !  and spatial orbitals into the expected order.
                     !
                     !  The simplest way of accomplishing this operation is to multiply the
                     !  alpha/beta to "textbook" parities for the initial and final states.
                     !
                     s2 = s2 + 0.5d0 * c(id) * c(is) * par(id) * par(is)
                  end do find_match
               end if
            end do orbital_j
         end do orbital_i
      end do scan_determinants
      
      return

    end subroutine calc_s2

!#######################################################################
! parity: taken straight from the superdyson code (see sd_core.f90)
!#######################################################################

    integer function parity(det)

      use channels

      implicit none

      integer, intent(in) :: det(:) ! Determinant, what else?
    
      integer :: ord(2*size(det)) ! Ordering array - can't get any 
                                  ! bigger than this
      integer :: io               ! Orbital number
      integer :: it               ! Current orbital number, 
                                  ! "textbook" order
      integer :: ia, ib           ! Current alpha/beta orbital number, 
                                  ! "alpha/beta" order
      integer :: nel, nalp, nbet  ! Number of electrons - total, 
                                  ! alpha, and beta
      logical :: sorted
      
      nel  = sum(abs(det))
      nalp = count(det==1) + count(det==2)
      nbet = nel - nalp
      
      !
      !  For each "textbook" orbital, figure out the desired "alpha/beta" index
      !
      ia = 0 ; ib = nalp ; it = 0
      order_initialize: do io=1,size(det)
         if (any(det(io)==(/ 1,2/))) then
            it = it + 1 ; ia = ia + 1 ; ord(it) = ia
         end if
         if (any(det(io)==(/-1,2/))) then
            it = it + 1 ; ib = ib + 1 ; ord(it) = ib
         end if
      end do order_initialize
      if (ia/=nalp .or. ib/=nel .or. it/=nel) then
         write (ilog,"('     Got: ia = ',i4,' ib = ',i4,' it = ',i4)") ia, ib, it
         write (ilog,"('Expected: ia = ',i4,' ib = ',i4,' it = ',i4)") nalp, nel, nel
         stop 'sd_core%parity - count error'
      end if
      !
      !  Count ordering permutations
      !
      parity = 1
      sorted = .false.
      ordering_sweep: do while(.not.sorted)
         sorted = .true.
         scan_order: do it=2,nel
            if (ord(it-1)<ord(it)) cycle scan_order
            ia        = ord(it-1)
            ord(it-1) = ord(it)
            ord(it)   = ia
            sorted    = .false. 
            parity    = -parity
         end do scan_order
      end do ordering_sweep
      
      return
      
    end function parity

!#######################################################################
! calc_overlaps: Calculation of the overlaps < target | ~ADC >,
!                where the bra corresponds to the target CI state, and
!                the kets to the approximate ADC states.
!#######################################################################

    subroutine calc_overlaps(gamess_internal)

      use channels
      use constants
      use parameters
      use import_gamess

      implicit none

      integer             :: i,k,m,p,q
      real(d)             :: detsmk
      type(gam_structure) :: gamess_internal

      real(d) :: tmp

!-----------------------------------------------------------------------
! Set the no. AOs
!-----------------------------------------------------------------------
      nao_targ=gamess_target%nbasis
      nao_adc=gamess_internal%nbasis

!-----------------------------------------------------------------------
! Precalculation of basis overlap arrays
!-----------------------------------------------------------------------
      call basis_overlaps(gamess_internal)

!-----------------------------------------------------------------------
! Calculation of the < target | ~ADC > overlaps
!-----------------------------------------------------------------------
      overlap=0.0d0

      ! Loop over ADC states
      do i=1,davstates

         ! Loop over target determinants
         do m=1,nsd_targ

            ! Loop over determinants for the ith ~ADC state
            do k=1,nsd_adc(i)

               ! Construct the matrix S^mk of overlaps between
               ! the MOs occupied in the target bra <m| and the
               ! ~ADC ket |k>
               call fill_spinorbital_integrals(onv_targ(:,m),&
                    onv_adc(:,k,i),smk,smo)

               ! Calculate det S^mk
               detsmk=determinant_overlap(smk)

               ! Add the contribution to < target | ~ADC_i >
               overlap(i)=overlap(i)+detsmk*c_targ(m)*c_adc(k,i)

            enddo

         enddo

      enddo

      return

    end subroutine calc_overlaps

!#######################################################################
! basis_overlaps: Calculates the target-internal AO and MO overlap 
!                 matrices
!#######################################################################

    subroutine basis_overlaps(gamess_internal)
      
      use channels
      use constants
      use parameters
      use import_gamess

      implicit none

      integer             :: i,j
      real(d)             :: val
      type(gam_structure) :: gamess_internal

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! <target|internal> AO overlap matrix
      allocate(sao(nao_targ,nao_adc))

      ! <target|internal> MO overlap matrix
      allocate(smo(nbas_targ,nbas))

      ! AO-to-MO transformation matrices
      allocate(ao2mo_adc(nao_adc,nbas))
      allocate(ao2mo_targ(nao_targ,nbas_targ))      

!-----------------------------------------------------------------------
! AO-to-MO transformation matrices
!-----------------------------------------------------------------------
      ao2mo_adc=gamess_internal%vectors(1:nao_adc,1:nbas)
      ao2mo_targ=gamess_target%vectors(1:nao_targ,1:nbas_targ)

!-----------------------------------------------------------------------
! Calculation of the target (bra) - internal (ket) AO overlap matrix
!-----------------------------------------------------------------------
      call gamess_1e_integrals('AO OVERLAP',sao,gamess_target,&
           gamess_internal)

!-----------------------------------------------------------------------
! Calculation of the target (bra) - internal (ket) MO overlap matrix
!-----------------------------------------------------------------------
      do i=1,nbas_targ
         do j=1,nbas
            val=dot_product(ao2mo_targ(:,i),matmul(sao,ao2mo_adc(:,j)))
            smo(i,j)=val
         enddo
      enddo

      return

    end subroutine basis_overlaps

!#######################################################################
! Calculate 1-particle overlap for all spin-orbitals in a given pair 
! of determinants.
! Taken from the superdyson code (see sd_core.f90).
!#######################################################################

    subroutine fill_spinorbital_integrals(occbra,occket,adet,amo)

      use channels
      use constants
      use parameters

      implicit none

      integer, intent(in)  :: occbra(:)            ! Parent determinant
      integer, intent(in)  :: occket(:)            ! Ion determinant
      real(d), intent(out) :: adet(:,:)            ! Integral matrix
                                                   ! for a given determinant
      real(d), intent(in)  :: amo (:,:)            ! Integral matrix for all MOs

      integer :: mobra, moket     ! Spatial MO indices
      integer :: spinbra, spinket ! Spin+space MO indices
      integer :: nalphabra        ! Number of spin-alpha electrons in the bra

      integer :: nmoket,nmobra

      nmobra=nbas_targ
      nmoket=nbas
      !
      adet = 0.0d0
      !
      !  Alpha-alpha overlaps - upper left corner of sdet
      !
      spinket = 0
      spinbra = 0
      ket_alpha: do moket=1,nmoket
         if (occket(moket)/=1 .and. occket(moket)/=2) cycle ket_alpha
         !
         spinket = spinket + 1
         spinbra = 0
         bra_alpha: do mobra=1,nmobra
            if (occbra(mobra)/=1 .and. occbra(mobra)/=2) cycle bra_alpha
            !
            spinbra = spinbra + 1
            adet(spinbra,spinket) = amo(mobra,moket)
         end do bra_alpha
      end do ket_alpha
      nalphabra = spinbra
      !
      !  Beta-beta overlaps - bottom right corner of sdet
      !
      ket_beta: do moket=1,nmoket
         if (occket(moket)/=-1 .and. occket(moket)/=2) cycle ket_beta
         !
         spinket = spinket + 1
         spinbra = nalphabra
         bra_beta: do mobra=1,nmobra
            if (occbra(mobra)/=-1 .and. occbra(mobra)/=2) cycle bra_beta
            !
            spinbra = spinbra + 1
            adet(spinbra,spinket) = amo(mobra,moket)
         end do bra_beta
      end do ket_beta
      
      return

    end subroutine fill_spinorbital_integrals

!#######################################################################
! determinant_overlap: Calculates the determinant of the matrix sdet of
!                      overlaps between the MOs in two determinants.
!                      Taken from superdyson.
!#######################################################################

    function determinant_overlap(sdet) result(overdet)
      
      use constants
      
      implicit none
      
      real(d), intent(in) :: sdet(:,:)
      real(d)             :: overdet
      real(d)             :: scr(nel,nel) ! Buffer for the 1-particle 
                                          !  overlap matrix 
      scr     = sdet
      overdet = determinant(scr)

      return

    end function determinant_overlap

!#######################################################################
! determinant: Calculates the determinant of a sqaure matrix.
!              Taken from superdyson.
!#######################################################################

    real(d) function determinant(mat)
    
      real(d), intent(inout) :: mat(:,:) ! Matrix destroyed on exit
      integer                :: order, info
      integer                :: ipvt(size(mat,dim=1))
      real(d)                :: work(size(mat,dim=1))
      real(d)                :: detx(2)
      
      order = size(mat,dim=1)
      if (order==0) stop 'null matrix passed to determinant'
      
      call dgefa(mat,order,order,ipvt,info)
   
      call dgedi(mat,order,order,ipvt,detx,work,10)

      determinant = detx(1) * 10.0d0**detx(2)

      return

    end function determinant

!#######################################################################

    subroutine write_table

      use channels

      implicit none

      integer :: i

!-----------------------------------------------------------------------
! Output overlaps, norms and S^2 expectation values
!-----------------------------------------------------------------------
      write(ilog,'(/,65a)') ('*',i=1,48)
      write(ilog,'(21x,a)') 'Summary'
      write(ilog,'(65a)') ('*',i=1,48)
      write(ilog,'(a)') '|i>     ndet    norm       <S^2>     <Target|i>'
      write(ilog,'(65a)') ('*',i=1,48)
      
      ! Target state
      write(ilog,'(a6,2x,i3,2(4x,F7.4))') 'Target',nsd_targ,&
           norm_targ,s2_targ/norm_targ**2

      ! Approximate ADC states
      do i=1,davstates
         write(ilog,'(i2,6x,i3,3(4x,F7.4))') i,nsd_adc(i),&
              norm_adc(i),s2_adc(i)/norm_adc(i)**2,overlap(i)
      enddo

      write(ilog,'(65a)') ('*',i=1,48)

      return

    end subroutine write_table

!#######################################################################

    subroutine select_istate
      
      use channels
      use constants
      use parameters
      use iomod

      implicit none

      integer :: i,count,id

!-----------------------------------------------------------------------
! Die here if: (1) There are no states ~ADC states |i> for which
!                  <Target|i> >= ovrthrsh;
!              (2) There are more than one ~ADC states |i> for which
!                  <Target|i> >= ovrthrsh.
!
! Else, set the initial state number and carry on.
!-----------------------------------------------------------------------
      count=0
      
      do i=1,davstates
         if (abs(overlap(i)).ge.ovrthrsh) then
            count=count+1
            id=i
         endif
      enddo

      if (count.eq.1) then
         statenumber=id
         write(ilog,'(/,x,a,x,i2,x,a)') 'State',id,'selected'
      else if (count.gt.1) then
         write(errmsg,'(a,x,F10.7)') &
              'There exists more than one state with <Target|i> >=',&
              ovrthrsh
         call error_control
      else
         write(errmsg,'(a,x,F10.7)') &
              'There exists no states with <Target|i> >=',ovrthrsh
         call error_control
      endif

      return

    end subroutine select_istate

!#######################################################################

  subroutine finalise

      implicit none

      deallocate(onv_adc)
      deallocate(c_adc)
      deallocate(s2_adc)
      deallocate(norm_adc)
      deallocate(nsd_adc)
      deallocate(sao)
      deallocate(overlap)
      deallocate(smk)

      return

    end subroutine finalise

!#######################################################################

  end module targetmatching
