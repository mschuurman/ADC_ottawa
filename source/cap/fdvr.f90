module fdvr

  implicit none

  save

  private :: dp
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter     :: dp=selected_real_kind(8)
  
contains

!######################################################################

  subroutine monomial_fdvr(gam,cap_mo,order,eta)

    use channels
    use parameters
    use misc, only: get_vdwr
    use autocapbox
    use import_gamess

    implicit none

    integer                        :: order
    real(dp)                       :: eta
    real(dp), dimension(nbas,nbas) :: cap_mo
    real(dp), allocatable          :: Q(:,:)
    real(dp), allocatable          :: r(:,:)
    type(gam_structure)            :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! MO-to-F-DVR transformation matrix
    allocate(Q(nbas,nbas))
    Q=0.0d0

    ! F-DVR points
    allocate(r(3,nbas))
    r=0.0d0    
    
!----------------------------------------------------------------------
! Determine the F-DVR points and transformation matrix
!----------------------------------------------------------------------
    call calc_dvr(Q,r)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Q)
    deallocate(r)
    
    return
    
  end subroutine monomial_fdvr

!######################################################################

  subroutine calc_dvr(Q,r)

    use channels
    use parameters
    use simdiag
    
    implicit none

    integer                        :: i,j
    real(dp), dimension(nbas,nbas) :: Q
    real(dp), dimension(3,nbas)    :: r
    real(dp), dimension(nbas,nbas) :: xdvr
    real(dp), allocatable          :: xmo(:,:,:)

!----------------------------------------------------------------------    
! Allocate arrays
!----------------------------------------------------------------------
    allocate(xmo(nbas,nbas,3))
    xmo=0.0d0

!----------------------------------------------------------------------
! Set up the MO representations of the three components of the
! position operator
!----------------------------------------------------------------------
    do i=1,3
       xmo(:,:,i)=-dpl_all(i,:,:)
    enddo

!----------------------------------------------------------------------
! Perform the approximate simultaneous diagonalisation of the three MO
! representations of the position operators.
! This will yield the MO-to-F-DVR transformation matrix, Q
!----------------------------------------------------------------------
    call simdiag_jacobi(xmo,Q,nbas,3)

!----------------------------------------------------------------------
! Calculate the F-DVR points
!----------------------------------------------------------------------
    do i=1,3
       xdvr=matmul(transpose(Q),matmul(xmo(:,:,i),Q))
       do j=1,nbas
          r(i,j)=xdvr(j,j)
       enddo
    enddo

    do i=1,nbas
       print*,i,r(i,1),r(i,2),r(i,3)
    enddo
    stop
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(xmo)
    
    return
    
  end subroutine calc_dvr

!######################################################################

end module fdvr
