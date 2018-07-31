!######################################################################
! slepian: a program to calculate discrete prolate spheroidal
!          sequences using Prieto's Multitaper Spectrum Estimation
!          library
!######################################################################

module global

  implicit none

  ! Double precision
  integer, parameter                   :: d=selected_real_kind(8)
  
  ! Number of points in the series
  integer                              :: npts

  ! Time-bandwidth product
  real(d)                              :: fw

  ! Desired number of Slepian functions
  integer                              :: nev

  ! Slepian functions
  real(d), dimension(:,:), allocatable :: v

  ! Eigenvalues of the Slepians
  real(d), dimension(:), allocatable   :: lambda

  ! 1-eigenvalues
  real(d), dimension(:), allocatable   :: theta
  
  save
  
end module global

!######################################################################

program slepian

  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinp

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  call alloc_arr
  
!----------------------------------------------------------------------
! Calculate the DPSSs
!----------------------------------------------------------------------
  call calc_dpss

!----------------------------------------------------------------------
! Output the DPSSs
!----------------------------------------------------------------------
  call wrout

!----------------------------------------------------------------------
! Check orthogonality
!----------------------------------------------------------------------
  call checkorthog
  
contains

!######################################################################
  
  subroutine rdinp

    use global
    
    implicit none

    character(len=60) :: string

!----------------------------------------------------------------------
! Die here if the incorrect no. arguments have been given
!----------------------------------------------------------------------
    if (iargc().ne.3) then
       write(6,'(/,2(2x,a,/))') 'Incorrect no. command line &
            arguments','Correct input: npts fw nev'
       stop
    endif
    
!----------------------------------------------------------------------
! Argument 1: no. points
!----------------------------------------------------------------------
    call getarg(1,string)
    read(string,*) npts

!----------------------------------------------------------------------    
! Argument 2: time-bandwidth product
!----------------------------------------------------------------------    
    call getarg(2,string)
    read(string,*) fw

!----------------------------------------------------------------------    
! Argument 3: number of Slepians
!----------------------------------------------------------------------    
    call getarg(3,string)
    read(string,*) nev
    
    return
    
  end subroutine rdinp

!######################################################################

  subroutine alloc_arr

    use global

    implicit none

    allocate(v(npts,nev))
    v=0.0d0

    allocate(lambda(nev))
    lambda=0.0d0

    allocate(theta(nev))
    theta=0.0d0
    
    return
    
  end subroutine alloc_arr
    
!######################################################################

  subroutine calc_dpss

    use global
    
    implicit none

    real(d), parameter :: pi=3.141592653589793d0, r2=1.4142135623731d0

    integer :: ntot, lh, nr, k, kr, k2, m11, n, i, ierr, nx
   
    real(d) :: atol, sn, eps1, lb, ub

    integer :: neven, nodd
    integer, dimension(nev) :: ind
    
    real(d) :: bw, hn, com
    real(d), dimension(nev) :: eigval
    real(d), dimension(:), allocatable :: fv1, fv2, fv3
    real(d), dimension(:,:), allocatable :: v2

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
    
  end subroutine calc_dpss

!######################################################################

  subroutine wrout

    use global
    
    implicit none

    integer           :: unit,i,j
    character(len=60) :: filename
    character(len=5)  :: ai
    
    unit=20

    ! Loop over the Slepians
    do i=1,nev

       ! Open the output file
       write(ai,'(i5)') i
       filename='dpss.'//trim(adjustl(ai))//'.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       
       ! Write the output file
       do j=1,npts
          write(unit,*) j,v(j,i)
       enddo
       
       ! Close the output file
       close(unit)
       
    enddo
    
    do i=1,nev
       print*,i,lambda(i),theta(i)
    enddo
    
    return

  end subroutine wrout

!######################################################################

  subroutine checkorthog

    use global
    
    implicit none

    integer :: i,j,n
    real(d) :: ovrlp

    do i=1,nev-1
       do j=i,nev

          ovrlp=v(1,i)*v(1,j)/2.0d0
          do n=2,npts-1
             ovrlp=ovrlp+v(n,i)*v(n,j)
          enddo
          ovrlp=ovrlp+v(npts,i)*v(npts,j)/2.0d0

          print*,i,j,ovrlp
          
       enddo
    enddo
          
    return
    
  end subroutine checkorthog
    
  
!######################################################################
  
end program slepian

