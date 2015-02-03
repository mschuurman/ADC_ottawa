  module qmath

  implicit none

  contains

!#######################################################################

    function qabs(arg)

      implicit none
      
      real*16            :: arg,qabs
      real*16, parameter :: qzero=0.0q0,qminus=-1.0q0
      
      if (arg.lt.qzero) then
         qabs=arg*qminus
      else
         qabs=arg
      endif

      return
      
    end function qabs

!#######################################################################   
 
    function qsqrt(arg)

      implicit none

      real*16            :: arg,qsqrt
      real*16, parameter :: qhalf=0.5q0

      qsqrt=arg**qhalf

      return

    end function qsqrt

!#######################################################################
    
    function qmax(arg1,arg2)

      implicit none

      real*16 :: arg1,arg2,qmax

      if (arg1.gt.arg2) then
         qmax=arg1
      else
         qmax=arg2
      endif

    end function qmax

!#######################################################################

    function qmin(arg1,arg2)

      implicit none

      real*16 :: arg1,arg2,qmin

      if (arg1.lt.arg2) then
         qmin=arg1
      else
         qmin=arg2
      endif

    end function qmin

!#######################################################################

    function qsign(arg1,arg2)

      implicit none

      real*16 :: arg1,arg2,qsign
      real*16 :: qzero=0.0q0,qminus=-1.0d0

      ! Absolute value of arg1
      if (arg1.lt.qzero) then
         qsign=qminus*arg1
      else
         qsign=arg1
      endif

      ! Multiply by 1.0q0 if arg2<zero
      if (arg2.lt.qzero) qsign=qsign*qminus

      return

    end function qsign

!#######################################################################
  
  end module qmath
