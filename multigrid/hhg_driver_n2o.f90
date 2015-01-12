!
!  Calculate total HHG spectrum from N2O alignment distribution.
!
!  Numerical simulations are for NNO aligned along the Z axis, with O
!  pointing towards the positive direction. The linearly polarized electric 
!  field is in the YZ plane, at 0-180 degree angles from the Z axis (in
!  10 degree increments). The angle is taken clockwise looking DOWN the 
!  X axis.
!
!  We would like to calculate HHG emission in the lab frame, with the
!  generating field oriented along the Z axis. The aligning field is
!  in the YZ plane, at an angle TA from the Z axis. The width of the 
!  alignment distribution is given by:
!  a) (aligned) cos^n(t), where t is the angle between the molecular 
!     axis and the aligning field direction; or
!  b) (antialigned) sin^n(t),
!
!
module hhg_driver_n2o
 use accuracy
 use math
 implicit none
 private
 public go
 !
 contains
 !
 subroutine go
   integer(ik)        :: t_step = 15
   integer(ik)        :: p_step = 10
   real(rk)           :: atheta   ! Direction of the aligning field
   integer(ik)        :: aorder   ! Cosine power in the alignement distribution
   character(len=128) :: namefmt  ! Format string for name conversion
   !
   character(len=40)  :: ashape   ! Either 'aligned' or 'antialigned'
   character(len=128) :: filename ! Results in the molecular frame
   integer(ik)        :: it, ip   ! Molecular orientation counters
   real(rk)           :: tl, th   ! Lower and upper bounds for theta
   real(rk)           :: d2r      ! Degrees to radians conversion
   real(rk)           :: cost     ! Cosine of the angle between molecular axis and alignment direction
   real(rk)           :: wgt      ! Weight of this integration point
   real(rk)           :: ea(3)    ! Euler angles
   real(rk)           :: rm(3,3)  ! Rotation matrix
   real(rk)           :: ad(3)    ! Direction of the alignment axis in lab coordinates
   real(rk)           :: mz(3)    ! Direction of the molecular Z axis in lab coordinates
   real(rk)           :: swgt     ! Consistency check - sum of integration weights
   !
   call accuracyInitialize
   d2r = pi / 180._rk
   !
   read (input,*) ashape, atheta, aorder, namefmt
   ad = (/ 0._rk, sin(d2r*atheta), cos(d2r*atheta) /)
   write (out,"('#')")
   write (out,"('#     Alignment shape = ',a     )") trim(ashape)
   write (out,"('#     Alignment theta = ',f12.4 )") atheta
   write (out,"('# Alignment direction = ',3f12.7)") ad
   write (out,"('#     Alignment order = ',i4    )") aorder
   write (out,"('#     File template   = ',a     )") trim(namefmt)
   write (out,"('#')")
   !
   swgt = 0._rk
   !
   !  theta segment is from it-7.5 to it+7.5 degrees, except at the edges
   !
   theta_loop: do it=0,180,t_step
     tl = max(  0._rk,it-0.5_rk*t_step)
     th = min(180._rk,it+0.5_rk*t_step)
     !
     write (filename,fmt=namefmt) it
     !
     !  Phi segment is from ip-5 to ip+5 degrees
     !
     phi_loop: do ip=0,360-p_step,p_step
       ea = (/ 0.5_rk*pi, d2r*it, 0.5_rk*pi + d2r*ip /)
       call MathRotationMatrix(ea,rm)
       mz   = matmul(rm,(/0._rk,0._rk,1._rk/))
       cost = dot_product(mz,ad)
       select case (ashape)
         case default
           write (out,"('Unknown alignment shape ',a)") trim(ashape)
           stop 'hhg_driver_n2o - in bad shape'
         case ('aligned')
           wgt = (d2r*p_step/fourpi)*(cos(d2r*tl)**(aorder+1) - cos(d2r*th)**(aorder+1))
         case ('antialigned')
           wgt = MathDoubleFactorial(aorder+1) * (d2r*p_step) / ( 2*MathDoubleFactorial(aorder) ) &
               * integrate_sinn(aorder+1,d2r*tl,d2r*th) / twopi
       end select
       swgt = swgt + wgt
       write (out,"('#')")
       write (out,"('#                 theta = ',i4)") it
       write (out,"('#                   phi = ',i4)") ip
       write (out,"('# molecular Z direction = ',3f12.7)") mz
       write (out,"('#       alignment angle = ',f12.4)") acos(cost)/d2r
       write (out,"('load ""',a,'""')") trim(filename)
       write (out,"('rotate ',3(g16.10,1x))") ea
       write (out,"('scale ',g16.10)") wgt
       write (out,"('accumulate')")
     end do phi_loop
   end do theta_loop
   write (out,"('#')")
   write (out,"('# Sum of integration weights = ',g16.8)") swgt
   write (out,"('#')")
 end subroutine go
 !
 function integrate_sinn(aorder,tl,th) result(v)
   integer(ik), intent(in) :: aorder ! Integrand: sin(x)**aorder
   real(rk), intent(in)    :: tl, th ! Integration limits
   real(rk)                :: v      ! Resulting integral
   !
   integer(ik), parameter :: npts = 1000_ik
   integer(ik)            :: ipt 
   real(rk)               :: t
   v = 0._rk
   stupid_integrate: do ipt=1,npts
     t = tl + (th-tl)*(ipt-0.5_rk)/npts
     v = v + sin(t)**aorder
   end do stupid_integrate
   v = v * (th-tl)/npts
 end function integrate_sinn
end module hhg_driver_n2o
!
program main
  use hhg_driver_n2o
  call go
end program main

