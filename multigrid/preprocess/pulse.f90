!
!  Electric field pulse shapes
!
 module pulse
   use accuracy
   implicit none
   private
   public FElectric_Pulse,PulseT

!
!  Data pulse
!
   type PulseT   ! External electric field pulse
     character(len=80) :: name            ! Name
     real(rk)          :: delay           ! shift pulse in time by "delay"
     real(rk)          :: freq            ! Frequency (2 pi full period)
     real(rk)          :: phase           ! Phase of the pulse at max intensity
     real(rk)          :: strength        ! Strength
     real(rk)          :: param           ! Parameter (e.g. gauss exponent)
     real(rk)          :: direction(3)    ! Direction 
   end type PulseT

  contains

   function FElectric_Pulse(Pulse,time) result (v)
     real(rk), intent(in)      :: time
     type(PulseT) , intent(in) :: Pulse
     real(rk)                  :: v

     real(rk)                  :: x
     real(rk)                  :: ev

     x  = Pulse%freq*(time-Pulse%delay)

     select case (Pulse%name) 
       case default
         write (out,"('FElectric_Pulse: pulse shape ',a,' unknown')") trim(Pulse%name)
         stop 'FElectric_Pulse - bad pulse'
       case('Const')
         v = Pulse%strength*sin(x+Pulse%phase)
       case('Gauss')
         ev = Pulse%param*x**2
         v  = 0
         if ( ev < log(safe_max) ) then
           v = Pulse%strength*(sin(x+Pulse%phase)+2.0_rk*cos(x+Pulse%phase)*Pulse%param*x) * exp(-ev)
         end if
      end select 

   end function FElectric_Pulse

 end module pulse
