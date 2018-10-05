!######################################################################
! dpssmt: routines for the calculation of discrete prolate spheroidal
!         sequences using Prieto's Multitaper Spectrum Estimation
!         library
!######################################################################

module dpssmt
  
  use constants
  use channels

  implicit none
  
contains

!######################################################################

  subroutine dpss_mt(npts,fw,nev,v,lambda)

    implicit none

    
    integer                       :: npts   ! Number of points in the
                                            ! series
    
    integer                       :: nev    ! Desired number of
                                            ! DPSSs

    real(dp)                      :: fw     ! NxW
    
    real(dp), dimension(npts,nev) :: v      ! DPSSs

    real(dp), dimension(nev)      :: lambda ! Eigenvalues

    real(dp), dimension(nev)      :: theta  ! 1-eigenvalues
    
!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
    call get_dpss(npts,fw,nev,v,lambda,theta)
    
    return
    
  end subroutine dpss_mt
    
!######################################################################

  subroutine get_dpss(npts,fw,nev,v,lambda,theta)

    implicit none

    integer                       :: npts,nev
    real(dp)                      :: fw
    real(dp), dimension(npts,nev) :: v
    real(dp), dimension(nev)      :: lambda,theta
    
    real(dp), parameter :: pi=3.141592653589793d0, r2=1.4142135623731d0

    integer :: ntot, lh, nr, k, kr, k2, m11, n, i, ierr, nx
   
    real(dp) :: atol, sn, eps1, lb, ub

    integer :: neven, nodd
    integer, dimension(nev) :: ind
    
    real(dp) :: bw, hn, com
    real(dp), dimension(nev) :: eigval
    real(dp), dimension(:), allocatable :: fv1, fv2, fv3
    real(dp), dimension(:,:), allocatable :: v2

    nx=mod(npts,2)
    lh=(npts/2)+nx

    nodd  = nev/2
    neven = nev - nodd
    
    bw = fw/dble(npts)
    com = cos(2.d0*pi*bw)
    
    hn = dble(npts-1)/2.d0

    allocate(fv1(lh))
    allocate(fv2(lh))
    allocate(fv3(lh))

    !
    !  Perform symmetry reduction to half size
    !

    do i=1,lh
       n = i-1
       
       !     Main diagonal
       
       fv1(i)    = com*(hn - dble(n))**2.d0
       
       !     sub diagonal
       
       fv2(i)    = dble(n*(npts-n))/2.d0
       fv3(i)    = fv2(i) * fv2(i)
    enddo
    
    if (nx.eq.0) then
       fv1(lh)   = com* (hn - dble(lh-1))**2.d0 + dble(lh*(npts-lh))/2.d0
    else
       fv2(lh)   = r2*fv2(lh)
       fv3(lh)   = 2.d0*fv3(lh)
    endif

    !
    !  Do the even tapers
    !    
    allocate(v2(lh,neven))
    
    
    eps1 = 0.d0
    m11 = lh-neven+1

    call tridib(lh,eps1,fv1,fv2,fv3,lb,ub,m11,neven,eigval,ind,ierr)
    if (ierr .ne. 0) then
       write(6,'(a,i5)') 'tridib error ',ierr
       stop
   endif

   call tinvit(lh,fv1,fv2,fv3,neven,eigval,ind,v2,ierr)
   if (ierr .ne. 0) then
      write(6,'(a,i5)') 'tinvit error ',ierr
      stop
   endif

   if (nx==1) then
      do k = 1,neven
         v2(lh,k) = r2*v2(lh,k)
      enddo
   endif

   do k = 1,neven
      kr = neven - k + 1
      k2 = 2*k - 1
      
      theta(k2) = eigval(kr)

!
! Expand the eigenfunctions
!

      nr=npts
      do i=1,lh
         v(i,k2) = v2(i,kr)
         v(nr,k2)= v2(i,kr)
         nr=nr-1
      enddo

!
! Normalize the eigenfunction
!

      sn=0.d0
      do n=1,npts
         sn=sn+v(n,k2)*v(n,k2)
      enddo
      
      sn=1.d0/sqrt(sn)

!
! Put eigenfunctions positive standard
!

      if ((v(lh+1,k2).lt.0.d0)) then
         sn=-sn
      endif
      
      do n=1,npts
         v(n,k2)=sn*v(n,k2)
      enddo
  
   enddo

!
!  Do the odd tapers
!

   if (nodd > 0) then   

      if (allocated(v2)) then
         deallocate(v2)
      endif   
      allocate(v2(lh-nx,nodd))

      do i=1,lh
         n = i-1
         fv1(i)  = com*(hn - dble(n))**2
         fv2(i)  = dble(n*(npts-n))/2.d0
         fv3(i)  = fv2(i) * fv2(i)
      enddo
   
      if (nx.eq.0) then
         fv1(lh)  =  com* (hn - dble(lh-1))**2 - dble(lh*(npts-lh))/2.d0
      endif
   
      eps1 = 0.d0
      m11 = (lh-nx)-nodd+1
   
      call tridib(lh-nx,eps1,fv1,fv2,fv3,lb,ub,m11,nodd,eigval,ind,ierr)   
      if (ierr .ne. 0) then
         write(6,'(a,i5)') 'tridib error ',ierr
         stop
      endif

      call tinvit(lh-nx,fv1,fv2,fv3,nodd,eigval,ind,v2,ierr)
      if (ierr .ne. 0) then
         write(6,'(a,i5)')'tinvit error ',ierr
         stop
      endif

      
      do k = 1,nodd
         kr = nodd - k + 1
         k2 = 2*k 

         theta(k2) = eigval(kr)

!
! Expand the eigenfunctions
!

         nr=npts
         do i=1,lh-nx
            v(i,k2) = v2(i,kr)
            v(nr,k2)= -v2(i,kr)
            nr=nr-1
         enddo
         if (nx == 1) then
            v(lh,k2) = 0.d0
         endif

!
! Normalize the eigenfunction
!

         sn=0.d0
         do n=1,npts
            sn=sn+v(n,k2)*v(n,k2)
         enddo
      
         sn=1.d0/sqrt(sn)

!
! Put eigenfunctions positive standard
!

         if ((nx==1) .and. (v(lh+1,k2).lt.0.d0)) then
            sn = -sn    
         endif
         if ((nx <= 0) .and. (v(lh+1,k2).lt.0.d0)) then
            sn=-sn
         endif

         do n=1,npts
            v(n,k2)=sn*v(n,k2)
         enddo

      enddo

   endif

      ntot=neven+nodd

!
!  Get the eigenvalues, by Quadrature (Chebychev)
!

   atol = 1.d-14

   call dpss_ev(npts,nev,bw,atol,v,lambda,theta)
   
   return
    
 end subroutine get_dpss
    
!######################################################################
  
end module dpssmt
