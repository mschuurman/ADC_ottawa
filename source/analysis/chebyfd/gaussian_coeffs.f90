!######################################################################
! gaussianmod: routines for the calculation of the coefficients in the
!              expansion of Slepian filter functions in terms of
!              Chebyshev polynomials of the first kind
!######################################################################

module gaussianmod

  use constants

  implicit none
  
  ! Scaled centres of the Gaussian filter functions
  real(dp), dimension(:), allocatable :: Enbar
  
contains
  
!######################################################################
  
  subroutine get_coeffs_gaussians

    use cfdmod
    
    implicit none

!----------------------------------------------------------------------
! Get the scaled centres of the Gaussian filter functions
!----------------------------------------------------------------------
    call get_Enbar

!----------------------------------------------------------------------    
! Calculate the coefficients in the expansion of the Gaussian
! filter functions in the Chebyshev polynomials
!----------------------------------------------------------------------
    call calc_coeffs
    
    return
    
  end subroutine get_coeffs_gaussians

!######################################################################

  subroutine get_Enbar

    use constants
    use cfdmod

    implicit none

    integer  :: n
    real(dp) :: de
    
!----------------------------------------------------------------------
! Allocate arrays    
!----------------------------------------------------------------------
    allocate(Enbar(nfsbas))
    Enbar=0.0d0

!----------------------------------------------------------------------
! Centres of the Gaussian filter functions
!----------------------------------------------------------------------
    de=(Ebbar-Eabar)/(nfsbas-1)
    
    do n=1,nfsbas
       Enbar(n)=Eabar+(n-1)*de
    enddo

    return
    
  end subroutine get_Enbar
    
!######################################################################

  subroutine calc_coeffs

    use constants
    use iomod
    use cfdmod
    
    implicit none

    integer                               :: Jdim,j,k,n
    real(dp), dimension(:), allocatable   :: costheta
    real(dp), dimension(:,:), allocatable :: cosktheta
    real(dp), dimension(:,:), allocatable :: gval

    integer                               :: unit
    real(dp), dimension(:,:), allocatable :: approxval
    character(len=3)                      :: an
    character(len=60)                     :: filename
    logical                               :: exists

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,2x,a)') 'Calculating the expansion coefficients...'
    
!----------------------------------------------------------------------
! No. quadrature points
!----------------------------------------------------------------------
    Jdim=Kdim+1

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(costheta(Jdim))
    costheta=0.0d0

    allocate(cosktheta(0:Kdim,Jdim))
    cosktheta=0.0d0

    allocate(gval(nfsbas,Jdim))
    gval=0.0d0

    allocate(approxval(Jdim,nfsbas))
    approxval=0.0d0
    
!----------------------------------------------------------------------
! Precalculate function values
!----------------------------------------------------------------------
    ! cos(theta)
    do j=1,Jdim
       costheta(j)=cos(pi*(j-0.5d0)/Jdim)
    enddo

    ! cos(k*theta)
    do j=1,Jdim
       do k=0,Kdim
          cosktheta(k,j)=cos(k*pi*(j-0.5d0)/Jdim)
       enddo
    enddo

    ! Gaussian values
    do j=1,Jdim
       do n=1,nfsbas
          gval(n,j)=gaussval(costheta(j),n)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the expansion coefficients using Gauss-Chebyshev
! quadrature
!----------------------------------------------------------------------
    fkn=0.0d0
    do n=1,nfsbas
       do k=0,Kdim
          do j=1,Jdim
             fkn(k,n)=fkn(k,n)+gval(n,j)*cosktheta(k,j)
          enddo
       enddo
    enddo
       
    ! Prefactors
    fkn(0,:)=fkn(0,:)/Jdim
    fkn(1:Kdim,:)=fkn(1:Kdim,:)*2.0d0/Jdim

!----------------------------------------------------------------------
! For checking purposes, output the Chebyshev expansions of the
! Gaussians
!----------------------------------------------------------------------
    ! Calculate the Chebyshev expansions of the Gaussians
    approxval=matmul(transpose(cosktheta),fkn)

    ! Output the Chebyshev expansions of the Gaussians
    inquire(file='gauss.cheby/.',exist=exists)
    if (exists) call system('rm gauss.cheby/*')
    if (.not.exists) call system('mkdir gauss.cheby')
    call freeunit(unit)
    do n=1,nfsbas
       write(an,'(i3)') n
       filename='gauss.cheby/gauss.cheby.'//trim(adjustl(an))//'.dat'
       open(unit,file=filename,form='formatted',status='unknown')
       write(unit,'(40a)') ('#',k=1,40)
       write(unit,'(a)') '# Point   Approximate      Actual'
       write(unit,'(40a)') ('#',k=1,40)
       do j=1,Jdim
          write(unit,'(3(2x,ES15.6))') &
               costheta(j),approxval(j,n),gval(n,j)
       enddo
       close(unit)
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(costheta)
    deallocate(cosktheta)
    deallocate(gval)
    deallocate(approxval)
    
    return
    
  end subroutine calc_coeffs

!######################################################################

  function gaussval(costheta,n)

    use constants
    use cfdmod
    
    implicit none

    integer  :: n
    real(dp) :: gaussval,costheta

    gaussval=exp(-((costheta-Enbar(n))/sigmabar)**2)
    
    return
    
  end function gaussval
    
!######################################################################
  
end module gaussianmod
