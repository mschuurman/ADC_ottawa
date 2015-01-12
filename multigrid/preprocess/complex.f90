module complex
!
!  Useful functions for handling complex numbers. These should
!  really be part of the language!
!
  use accuracy
  implicit none
  private
  public cmplx_p2c, cmplx_c2p
!
  contains
    !
    !  Convert complex number in polar format to Cartesian 
    !  representation
    !
!   elemental function cmplx_p2c(p) result(x)
    function cmplx_p2c(p) result(x)
      complex(rk), intent(in) :: p
      complex(rk)             :: x
      !
      real(rk) :: r, arg
      !
      r   = real(p,kind=rk)
      arg = aimag(p)
      x   = r*cmplx(cos(arg),sin(arg),kind=rk)
    end function cmplx_p2c
    !
    !  Convert complex number in Cartesian representation to polar
    !
!   elemental function cmplx_c2p(x) result(p)
    function cmplx_c2p(x) result(p)
      complex(rk), intent(in) :: x
      complex(rk)             :: p
      !
      real(rk) :: r, arg
      !
      r = abs(x)
      if (r==0) then
        arg = 0 ! argument is undefined, choose zero
      else
        arg = atan2(aimag(x),real(x,kind=rk))
      end if
      p = cmplx(r,arg,kind=rk)
    end function cmplx_c2p
    !
end module complex
