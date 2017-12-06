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
    real(dp), allocatable          :: w(:)
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
    
    ! F-DVR weights
    allocate(w(nbas))
    w=0.0d0

!----------------------------------------------------------------------
! Determine the F-DVR points, weights and transformation matrix
!----------------------------------------------------------------------
    call calc_dvr(Q,r,w,gam)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Q)
    deallocate(r)
    deallocate(w)

    return
    
  end subroutine monomial_fdvr

!######################################################################

  subroutine calc_dvr(Q,r,w,gam)

    use channels
    use parameters
    use simdiag
    use import_gamess
    use density, only: get_ao_values

    implicit none

    integer                        :: i,j,nao
    real(dp), dimension(nbas,nbas) :: Q
    real(dp), dimension(3,nbas)    :: r
    real(dp), dimension(3)         :: ri
    real(dp), dimension(nbas)      :: w
    real(dp), allocatable          :: xdvr(:,:)
    real(dp), allocatable          :: xmo(:,:,:)
    real(dp), allocatable          :: aovalues(:)
    real(dp), allocatable          :: movalues(:)
    real(dp), allocatable          :: dvrvalues(:)
    type(gam_structure)            :: gam

!----------------------------------------------------------------------    
! Allocate arrays
!----------------------------------------------------------------------
    ! No. AOs
    nao=gam%nbasis

    allocate(xmo(nbas,nbas,3))
    xmo=0.0d0

    allocate(xdvr(nbas,nbas))
    xdvr=0.0d0

    allocate(aovalues(nao))
    aovalues=0.0d0

    allocate(movalues(nbas))
    movalues=0.0d0

    allocate(dvrvalues(nbas))
    dvrvalues=0.0d0

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

!----------------------------------------------------------------------
! Calculate the F-DVR weights
!----------------------------------------------------------------------
    do i=1,nbas
       
       ! Current F-DVR point
       ri(1:3)=r(1:3,i)

       ! AO values at the current F-DVR point
       call get_ao_values(gam,aovalues,ri,nao)

       ! MO values at the current F-DVR point
       movalues=matmul(transpose(ao2mo),aovalues)

       ! F-DVR values at the current F-DVR point
       dvrvalues=matmul(transpose(Q),movalues)

       ! F-DVR weight
       w(i)=1.0d0/(dvrvalues(i)**2)

    enddo

    ! CHECK
    do i=1,nbas
       write(ilog,'(i4,3(3x,F8.4),3x,ES15.8)') i,r(1,i),r(2,i),r(3,i),w(i)
    enddo
    stop
    ! CHECK

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(xmo)
    deallocate(xdvr)
    deallocate(aovalues)
    deallocate(movalues)
    deallocate(dvrvalues)
    
    return
    
  end subroutine calc_dvr

!######################################################################

end module fdvr
