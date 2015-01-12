!
!  This is a Fortran-90 wrapper around SLATEC (Fullerton's) Bessel functions.
!  With intel compiler, these functions are included in the portability library.
!  Added Struve H functions of the integer order 0 and 1, using MacLeod's 
!  implementation.
!
 module bessel
   use accuracy
   implicit none

   private
   public besselJ, struveH

   interface besselJ
     module procedure besselJ_real
     module procedure besselJ_double
   end interface ! besselJ

   interface struveH
     module procedure struveH_real
     module procedure struveH_double
   end interface ! stuveH

 contains

   function besselJ_real(n,x) result(v)
     integer(ik), intent(in)    :: n ! Order of the Bessel function
     real(kind=srk), intent(in) :: x ! Argument
     real(kind=srk)             :: v 
!
!    real(kind=srk), external   :: besj0, besj1, besjn
!    real(kind=srk), intrinsic  :: besj0, besj1, besjn
     real(kind=srk)             :: besj0, besj1, besjn
     !
     select case (n)
       case (0);     v = besj0(x)
       case (1);     v = besj1(x)
       case default; v = besjn(n,x)
     end select
   end function besselJ_real

   function besselJ_double(n,x) result(v)
     integer(ik), intent(in)      :: n ! Order of the Bessel function
     real(kind=drk), intent(in) :: x ! Argument
     real(kind=drk)             :: v 
!    real(kind=drk), external   :: dbesj0, dbesj1, dbesjn
!    real(kind=drk), intrinsic  :: dbesj0, dbesj1, dbesjn
     real(kind=drk)             :: dbesj0, dbesj1, dbesjn
     !
     select case (n)
       case (0);     v = dbesj0(x)
       case (1);     v = dbesj1(x)
       case default; v = dbesjn(n,x)
     end select
   end function besselJ_double
 
   function struveH_real(n,x) result(v)
     integer(ik), intent(in)    :: n ! Order of the Struve H function
     real(kind=srk), intent(in) :: x ! Argument
     real(kind=srk)             :: v 
     real(kind=srk), external   :: rstrvh0, rstrvh1
     !
     select case (n)
       case (0);     v = rstrvh0(x)
       case (1);     v = rstrvh1(x)
       case default; stop 'struveH_real - only orders 0 and 1 are implemented'
     end select
   end function struveH_real

   function struveH_double(n,x) result(v)
     integer(ik), intent(in)      :: n ! Order of the Struve H function
     real(kind=drk), intent(in) :: x ! Argument
     real(kind=drk)             :: v 
     real(kind=drk), external   :: dstrvh0, dstrvh1
     !
     select case (n)
       case (0);     v = dstrvh0(x)
       case (1);     v = dstrvh1(x)
       case default; stop 'struveH_double - only orders 0 and 1 are implemented'
     end select
   end function struveH_double

 end module bessel
