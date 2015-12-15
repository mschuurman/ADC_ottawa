  module targetmatching

    implicit none

  contains

!#######################################################################

    subroutine target_master(kpq,ndim)

      use channels
      use parameters
      use constants

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,maxnsd
      integer, dimension(davstates)             :: nsd
      integer, dimension(:,:,:), allocatable    :: onv_adc,onv_targ
      real(d), dimension(:,:), allocatable      :: c_adc,c_targ

!-----------------------------------------------------------------------
! Determine the maximum number of Slater determinants entering into
! the truncated, zeroth-order expansions of the ADC states
!-----------------------------------------------------------------------
      call getdim_sd(kpq,ndim,nsd,maxnsd)

!-----------------------------------------------------------------------
! Construct the truncated list of ON vectors for the zeroth-order
! expansions of the ADC states
!-----------------------------------------------------------------------
      call fill_onv_adc(onv_adc,c_adc,nsd,maxnsd,kpq,ndim)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(onv_adc)
      deallocate(c_adc)

      STOP

      return

    end subroutine target_master

!#######################################################################
! getdim_sd: determines the number of Slater determinants entering
!            into the zeroth-order expansion of the Davidson states
!            with coefficients above the user set threshold
!#######################################################################

    subroutine getdim_sd(kpq,ndim,nsd,maxnsd)
      
      use channels
      use parameters
      use constants
      use iomod, only: freeunit

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,maxnsd,i,idav,&
                                                   itmp
      integer, dimension(davstates)             :: nsd
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

         ! Read the next state vector from file
         read(idav) itmp,ftmp,coeff
         
         ! Determine the no. contributing Slater determinants
         call getnsd(kpq,nsd(i),coeff,ndim)

         ! Update maxnsd
         if (nsd(i).gt.maxnsd) maxnsd=nsd(i)

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

    subroutine getnsd(kpq,nsd,coeff,ndim)
      
      use channels
      use parameters
      use constants

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: nsd,ndim,n,count,&
                                                   nI,nII
      real(d), dimension(ndim)                  :: coeff
      real(d)                                   :: coeffsd,invsqrt2,&
                                                   invsqrt12

      invsqrt2=1.0d0/sqrt(2.0d0)
      invsqrt12=1.0d0/sqrt(12.0d0)

      nsd=0

!-----------------------------------------------------------------------
! 1h1p ISs
!-----------------------------------------------------------------------
      do n=1,kpq(1,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)
      do n=count+1,count+kpq(2,0)
         coeffsd=coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+1
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a|=b, i=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)
      do n=count+1,count+kpq(3,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2
      enddo

!-----------------------------------------------------------------------
! 2h2p ISs, a=b, i|=j
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)
      do n=count+1,count+kpq(4,0)
         coeffsd=invsqrt2*coeff(n)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2
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
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2
         coeffsd=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2

         ! Spin case II only
         coeffsd=2.0d0*invsqrt12*coeff(nII)
         if (abs(coeffsd).gt.detthrsh) nsd=nsd+2

      enddo

      return

    end subroutine getnsd

!#######################################################################

    subroutine fill_onv_adc(onv_adc,c_adc,nsd,maxnsd,kpq,ndim)

      use channels
      use parameters
      use constants
      use iomod, only: freeunit

      implicit none
      
      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer, dimension(:,:,:), allocatable    :: onv_adc
      integer, dimension(davstates)             :: nsd
      integer                                   :: maxnsd,ndim,s,&
                                                   idav,itmp
      real(d), dimension(:,:), allocatable      :: c_adc
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
         call fill_onv_adc_1state(onv_adc(:,:,s),maxnsd,nsd(s),&
              c_adc(:,s),kpq,coeff,ndim)

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

    subroutine fill_onv_adc_1state(onv,maxnsd,nsd,coeffsd,kpq,coeff,&
         ndim)

      use channels
      use parameters
      use constants

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer, dimension(nbas,maxnsd)           :: onv
      integer                                   :: maxnsd,nsd,ndim,&
                                                   n,count,i,j,a,b,&
                                                   ksd,nI,nII
      real(d), dimension(maxnsd)                :: coeffsd
      real(d), dimension(ndim)                  :: coeff
      real(d)                                   :: ftmp

      ! Common prefactors
      real(d), parameter                 :: invsqrt2=1.0d0/sqrt(2.0d0),&
                                            invsqrt12=1.0d0/sqrt(12.0d0)

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
            onv(i,ksd)=+1
            onv(i,ksd)=-1
            coeffsd(ksd)=invsqrt2*coeff(n)
            
            ! i_alpha -> a_alpha
            ksd=ksd+1
            onv(i,ksd)=-1
            onv(i,ksd)=+1
            coeffsd(ksd)=invsqrt2*coeff(n)

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
            onv(i,ksd)=0
            onv(a,ksd)=2
            coeffsd(ksd)=coeff(n)
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
            onv(i,ksd)=0
            onv(a,ksd)=+1
            onv(b,ksd)=-1
            coeffsd(ksd)=invsqrt2*coeff(n)

            ! i_alpha, i_beta -> a_beta, b_alpha
            ksd=ksd+1
            onv(i,ksd)=0
            onv(a,ksd)=-1
            onv(b,ksd)=+1
            coeffsd(ksd)=invsqrt2*coeff(n)

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
            onv(i,ksd)=-1
            onv(j,ksd)=+1
            onv(a,ksd)=2
            coeffsd(ksd)=invsqrt2*coeff(n)

            ! i_beta, j_alpha -> a_alpha, a_beta
            ksd=ksd+1
            onv(i,ksd)=+1
            onv(j,ksd)=-1
            onv(a,ksd)=2
            coeffsd(ksd)=invsqrt2*coeff(n)

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
            onv(i,ksd)=-1
            onv(j,ksd)=+1
            onv(a,ksd)=+1
            onv(b,ksd)=-1
            coeffsd(ksd)=0.5d0*coeff(nI)+invsqrt12*coeff(nII)

            ! i_beta, j_alpha -> a_beta, b_alpha
            ksd=ksd+1
            onv(i,ksd)=+1
            onv(j,ksd)=-1
            onv(a,ksd)=-1
            onv(b,ksd)=+1
            coeffsd(ksd)=0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         endif

         ftmp=abs(-0.5d0*coeff(nI)+invsqrt12*coeff(nII))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, j_beta -> a_beta, b_alpha
            ksd=ksd+1
            onv(i,ksd)=-1
            onv(j,ksd)=+1
            onv(a,ksd)=-1
            onv(b,ksd)=+1
            coeffsd(ksd)=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)

            ! i_beta, j_alpha -> a_alpha, b_beta
            ksd=ksd+1
            onv(i,ksd)=+1
            onv(j,ksd)=-1
            onv(a,ksd)=+1
            onv(b,ksd)=-1
            coeffsd(ksd)=-0.5d0*coeff(nI)+invsqrt12*coeff(nII)
         endif

         !--------------------------------------------------------------
         ! (2) Determinants with contributions from spin case II only
         !--------------------------------------------------------------
         ftmp=abs(2.0d0*invsqrt12*coeff(nII))
         if (ftmp.ge.detthrsh) then
            ! i_alpha, j_alpha -> a_alpha, b_alpha
            ksd=ksd+1
            onv(i,ksd)=-1
            onv(j,ksd)=-1
            onv(a,ksd)=+1
            onv(b,ksd)=+1
            coeffsd(ksd)=2.0d0*invsqrt12*coeff(nII)

            ! i_beta, j_beta -> a_beta, b_beta
            ksd=ksd+1
            onv(i,ksd)=+1
            onv(j,ksd)=+1
            onv(a,ksd)=-1
            onv(b,ksd)=-1
            coeffsd(ksd)=2.0d0*invsqrt12*coeff(nII)
         endif

      enddo

      return

    end subroutine fill_onv_adc_1state

!#######################################################################

  end module targetmatching
