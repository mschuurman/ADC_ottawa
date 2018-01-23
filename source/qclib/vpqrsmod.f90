module vpqrsmod
  use integrals_mo2e
  use parameters
  use constants

!#######################################################################
! Temporary code: force vpqrs to run in-core.
! 
! To revert back to the option of using external memory, delete this
! function, uncomment the code below it as well as the code referring
! to vpqrs in adc/adc.f90
!#######################################################################  

contains

  function vpqrs(r,s,u,v) result(val)
    implicit none
    
    integer,intent(in) :: r,s,u,v
    real(d)            :: val
    
    val = moIntegrals(r, s, u, v)

    return
    
  end function vpqrs

  !
  !
  subroutine init_moIntegrals(nmo)
    implicit none
    integer, intent(in)    :: nmo

    allocate(moIntegrals(nmo,nmo,nmo,nmo))
    moIntegrals = rzero

    return
  end subroutine init_moIntegrals 

   ! set, via permutational symmetry, the integrals of the type (ij|ij)
   ! and (ij|ik) and (ii|jk)
  subroutine unfold_iajb(moints, nmo, nocc)
    implicit none
    type(moint2e_cache), intent(in)  :: moints
    integer, intent(in)              :: nmo, nocc

    integer                          :: i,j,a,b
    real(d)                          :: int_val

    do i = 1,nocc
     do a = 1,nmo
      do j = 1,nocc
       do b = 1,nmo
        int_val = moints%buffer_real(i,a,j,b)
        moIntegrals(i, a, j, b)    = int_val
        moIntegrals(i, a, b, j)    = int_val
        moIntegrals(a, i, j, b)    = int_val
        moIntegrals(a, i, b, j)    = int_val
       enddo
      enddo
     enddo
    enddo

    return
  end subroutine unfold_iajb

   ! set, via permutational symmetry, the integrals of the type (ij|ij)
   ! and (ij|ik) and (ii|jk)
  subroutine unfold_ijab(moints, nmo, nocc)
    implicit none
    type(moint2e_cache), intent(in)  :: moints
    integer, intent(in)              :: nmo, nocc

    integer                          :: i,j,a,b
    real(d)                          :: int_val

    do i = 1,nocc
     do j = 1,nocc
      do a = 1,nmo
       do b = 1,nmo
        int_val = moints%buffer_real(i,j,a,b)
        moIntegrals(i, j, a, b)    = int_val
        moIntegrals(a, b, i, j)    = int_val
       enddo
      enddo
     enddo
    enddo

    return
  end subroutine unfold_ijab

  !
  ! 
  subroutine unfold_abcd(moints, nmo, nocc)
    implicit none
    type(moint2e_cache), intent(in)  :: moints
    integer, intent(in)              :: nmo, nocc

    integer                          :: a,b,c,f,nvrt
    real(d)                          :: int_val

    nvrt = nmo - nocc

    do a = 1,nvrt
     do b = 1,nvrt
      do c = 1,nvrt
       do f = 1,nvrt
        int_val = moints%buffer_real(a,b,c,f)
        moIntegrals(a+nocc,b+nocc,c+nocc,f+nocc)    = int_val
       enddo
      enddo
     enddo
    enddo

    return
  end subroutine unfold_abcd

   ! set, via permutational symmetry, the integrals of the type (ij|ij)
   ! and (ij|ik) and (ii|jk)
  subroutine unfold_all(moints, nmo)
    implicit none
    type(moint2e_cache), intent(in)  :: moints
    integer, intent(in)              :: nmo

    integer                          :: i,j,k,l
    real(d)                          :: int_val

    do i = 1,nmo
     do j = 1,nmo
      do k = 1,nmo
       do l = 1,nmo
        int_val = moints%buffer_real(i,j,k,l)
        moIntegrals(i, j, k, l)    = int_val
       enddo
      enddo
     enddo
    enddo

    return
  end subroutine unfold_all

!#######################################################################
!  use constants
!!  
!  save
!
!  abstract interface
!     function func (r,s,u,v)
!       real*8 :: func
!       integer, intent (in) :: r,s,u,v
!     end function func
!  end interface
!
!  procedure(func), pointer :: vpqrs
!  
! contains
!   
!!#######################################################################
!  
!   function vpqrs_incore(r,s,u,v) result(val)
!    
!    use parameters
!
!    implicit none
!    
!    integer,intent(in) :: r,s,u,v
!    real(d)            :: val
!    
!    val = real(moIntegrals%buffer_real(r,s,u,v),kind=d)
!    
!    return
!    
!  end function vpqrs_incore
!
!!#######################################################################
!  
!  function vpqrs_ext(r,s,u,v) result(val)
!    
!    use parameters
!
!    implicit none
!    
!    integer,intent(in) :: r,s,u,v
!    integer            :: r2
!    real(d)            :: val
!    
!    r2 = u
!     if(moIntegrals%mo_l /= v .and. moIntegrals%mo_l /= u)then
!       call fetch_moint2e(moIntegrals,v)
!     else
!       if(moIntegrals%mo_l == u) r2 = v
!     endif
!     val = real(moIntegrals%buffer_real(r,s,r2,1),kind=d)
!
!    return    
!    
!  end function vpqrs_ext
!
!!#######################################################################
  
end module vpqrsmod
