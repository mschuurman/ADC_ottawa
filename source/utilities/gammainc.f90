!######################################################################
! Routines for the evaluation of the incomplete Gamma function ratio.
! Taken from the freely available ASA239 library of John Burkardt.
!######################################################################

module gammainc

  use iomod
  use channels
  
  implicit none

contains

!######################################################################
!
! gamma_inc_ratio computes the incomplete Gamma function ratio.
!
!  Discussion:
!
!    GAMMA_INC(P,       0) = 0,
!    GAMMA_INC(P,Infinity) = 1.
!
!    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by B L Shea.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BL Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS239,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the exponent parameter.
!    0.0D+00 < P.
!
!    Input, real ( kind = 8 ) X, the integral limit parameter.
!    If X is less than or equal to 0, GAMMA_INC is returned as 0.
!
!    Output, real ( kind = 8 ) GAMMA_INC, the value of the function.
!
!######################################################################

  function gamma_inc_ratio ( p, x )
    
    implicit none

    real ( kind = 8 ) a
    real ( kind = 8 ) arg
    real ( kind = 8 ) b
    real ( kind = 8 ) c
    real ( kind = 8 ) cdf
    real ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00
    real ( kind = 8 ) gamma_inc_ratio
    real ( kind = 8 ), parameter :: overflow = 1.0D+37
    real ( kind = 8 ) p
    real ( kind = 8 ), parameter :: plimit = 1000.0D+00
    real ( kind = 8 ) pn1
    real ( kind = 8 ) pn2
    real ( kind = 8 ) pn3
    real ( kind = 8 ) pn4
    real ( kind = 8 ) pn5
    real ( kind = 8 ) pn6
    real ( kind = 8 ) rn
    real ( kind = 8 ), parameter :: tol = 1.0D-07
    real ( kind = 8 ) x
    real ( kind = 8 ), parameter :: xbig = 1.0D+08
    
    gamma_inc_ratio = 0.0D+00

    if ( p <= 0.0D+00 ) then
       errmsg='gamma_inc_ratio - Fatal error - &
            Parameter P <= 0.'
       call error_control
    end if
    
    if ( x <= 0.0D+00 ) then
       gamma_inc_ratio = 0.0D+00
       return
    end if
    !
    !  Use a normal approximation if PLIMIT < P.
    !
    if ( plimit < p ) then
       pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** ( 1.0D+00 / 3.0D+00 ) &
            + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
       call normal_01_cdf ( pn1, cdf )
       gamma_inc_ratio = cdf
       return
    end if
    !
    !  Is X extremely large compared to P?
    !
    if ( xbig < x ) then
       gamma_inc_ratio = 1.0D+00
       return
    end if
    !
    !  Use Pearson's series expansion.
    !  (P is not large enough to force overflow in the log of Gamma.
    !
    if ( x <= 1.0D+00 .or. x < p ) then
       
       arg = p * log ( x ) - x - log_gamma ( p + 1.0D+00 )
       c = 1.0D+00
       gamma_inc_ratio = 1.0D+00
       a = p
       
       do

          a = a + 1.0D+00
          c = c * x / a
          gamma_inc_ratio = gamma_inc_ratio + c
          
          if ( c <= tol ) then
             exit
          end if

       end do

       arg = arg + log ( gamma_inc_ratio )
       
       if ( exp_arg_min <= arg ) then
          gamma_inc_ratio = exp ( arg )
       else
          gamma_inc_ratio = 0.0D+00
       end if

    else
       !
       !  Use a continued fraction expansion.
       !
       arg = p * log ( x ) - x - log_gamma ( p )
       a = 1.0D+00 - p
       b = a + x + 1.0D+00
       c = 0.0D+00
       pn1 = 1.0D+00
       pn2 = x
       pn3 = x + 1.0D+00
       pn4 = x * b
       gamma_inc_ratio = pn3 / pn4

       do

          a = a + 1.0D+00
          b = b + 2.0D+00
          c = c + 1.0D+00
          pn5 = b * pn3 - a * c * pn1
          pn6 = b * pn4 - a * c * pn2
          
          if ( 0.0D+00 < abs ( pn6 ) ) then
         
             rn = pn5 / pn6

             if ( abs ( gamma_inc_ratio - rn ) <= min ( tol, tol * rn ) ) then

                arg = arg + log ( gamma_inc_ratio )
                
                if ( exp_arg_min <= arg ) then
                   gamma_inc_ratio = 1.0D+00 - exp ( arg )
                else
                   gamma_inc_ratio = 1.0D+00
                end if
                
                return

             end if

             gamma_inc_ratio = rn

          end if
          
          pn1 = pn3
          pn2 = pn4
          pn3 = pn5
          pn4 = pn6
          !
          !  Rescale terms in continued fraction if terms are large.
          !
          if ( overflow <= abs ( pn5 ) ) then
             pn1 = pn1 / overflow
             pn2 = pn2 / overflow
             pn3 = pn3 / overflow
             pn4 = pn4 / overflow
          end if
      
       end do

    end if
    
    return

  end function gamma_inc_ratio

!######################################################################
!
! normal_01_cdf evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
!######################################################################

  subroutine normal_01_cdf ( x, cdf )
  
    implicit none

    real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
    real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
    real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
    real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
    real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
    real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
    real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
    real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
    real ( kind = 8 ), parameter :: b1 = 3.8052D-08
    real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
    real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
    real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
    real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
    real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
    real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
    real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
    real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
    real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
    real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
    real ( kind = 8 ) cdf
    real ( kind = 8 ) q
    real ( kind = 8 ) x
    real ( kind = 8 ) y
    !
    !  |X| <= 1.28.
    !
    if ( abs ( x ) <= 1.28D+00 ) then
       
       y = 0.5D+00 * x * x
       
       q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
            + a6 / ( y + a7 ) ) ) )
       !
       !  1.28 < |X| <= 12.7
       !
    else if ( abs ( x ) <= 12.7D+00 ) then

       y = 0.5D+00 * x * x

       q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
            + b2 / ( abs ( x ) + b3 &
            + b4 / ( abs ( x ) - b5 &
            + b6 / ( abs ( x ) + b7 &
            - b8 / ( abs ( x ) + b9 &
            + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
       !
       !  12.7 < |X|
       !
    else
       
       q = 0.0D+00
       
    end if
    !
    !  Take account of negative X.
    !
    if ( x < 0.0D+00 ) then
       cdf = q
    else
       cdf = 1.0D+00 - q
    end if

    return
  end subroutine normal_01_cdf

!######################################################################
  
end module gammainc
