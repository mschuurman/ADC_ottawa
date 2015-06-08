module vpqrsmod

  use constants
  
  save

  abstract interface
     function func (r,s,u,v)
       real*8 :: func
       integer, intent (in) :: r,s,u,v
     end function func
  end interface

  procedure(func), pointer :: vpqrs
  
 contains
   
!#######################################################################
  
   function vpqrs_incore(r,s,u,v) result(val)
    
    use parameters

    implicit none
    
    integer,intent(in) :: r,s,u,v
    real(d)            :: val
    
    val = real(moIntegrals%buffer_real(r,s,u,v),kind=d)
    
    return
    
  end function vpqrs_incore

!#######################################################################
  
  function vpqrs_ext(r,s,u,v) result(val)
    
    use parameters

    implicit none
    
    integer,intent(in) :: r,s,u,v
    integer            :: r2
    real(d)            :: val
    
    r2 = u
     if(moIntegrals%mo_l /= v .and. moIntegrals%mo_l /= u)then
       call fetch_moint2e(moIntegrals,v)
     else
       if(moIntegrals%mo_l == u) r2 = v
     endif
     val = real(moIntegrals%buffer_real(r,s,r2,1),kind=d)

    return    
    
  end function vpqrs_ext

!#######################################################################
  
end module vpqrsmod
