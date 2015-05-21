program tquadtsp

!   David H Bailey   13 May 2015

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2015 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program should be compiled with Version 2 of MPFUNF (in mpfunf2.f90).

!   This program demonstrates the tanh-sinh quadrature routine quadtsp and the
!   initialization routine initqts.  Quadtsp differs from quadts (in file
!   tquadts.f90) in that it does not progressively increase the quadrature
!   level (grid resolution) but instead it evaluates the integral with a preset
!   level, passed in the argument levq.   Quadtsp is particuarly well suited for
!   parallel execution, and OpenMP parallelization directives are included below,
!   both in initqts and quadtsp.

!   While tanh-sinh quadrature is not as efficient as Gaussian quadrature
!   for completely regular functions, it works well for many functions with
!   an integrable singularity at one or both of the endpoints.  Also, the cost
!   of computing abscissas and weights for tanh-sinh quadrature is much less
!   than that of Gaussian quadrature.  For details, see:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   The function to be integrated must be defined in an external function
!   subprogram (see samples below), and the name of the function must be
!   included in a "type (mp_real)" and an "external" statement in the program
!   calling quadtsp.  See the sample functions and the calls to quadtsp below.
!   Before quadtsp is called, initqts must be called to initialize the wk and
!   xk arrays.

!   Both initqts and quadtsp are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   For some functions, it is important that the endpoints x1 and x2 be
!   computed to high precision (nwords2 words) in the calling program, that
!   these high-precision values be passed to quadtsp (where scaled abscissas
!   are calculated to high precision), and that the function definition
!   itself uses these high-precision scaled abscissas in any initial
!   subtractions or other sensitive operations involving the input argument.
!   Otherwise the accuracy of the quadrature result might only be half as
!   high as it otherwise could be.  See the function definitions of fun06,
!   fun07, fun09 and fun10 for examples on how this is done.  Otherwise the
!   function evaluation can and should be performed with low precision
!   (nwords1 words) for faster run times.  The two precision levels (nwords1
!   and nwords2) are computed below based on ndp1 and ndp2, which are
!   specified by the user in the main program (a few lines below).

!   These inputs are set in the parameter statement below:
!   ndp1   Primary ("low") precision in digits; this is the target accuracy
!          of quadrature results.
!   ndp2   Secondary ("high") precision in digits. By default, ndp2 = 2*ndp1.
!   neps1  Log10 of the primary tolerance. By default, neps1 = - dp1.
!   neps2  Log10 of the secondary tolerance. By default, neps2 = -ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. nq1 must be at least 2.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!          default it is set to 12 * 2^nq1. Increase nq2 if directed by a 
!          message produced in initqts.

!   In addition, the parameter levq (the quadrature level), which varies
!   with each problem, is passed as a subroutine argument in the calls below.

use mpmodule
implicit none
integer i, levq, ndp1, ndp2, neps1, neps2, nq1, nq2, nwords1, nwords2, n1
parameter (ndp1 = 500, ndp2 = 1000, neps1 = -ndp1, neps2 = -ndp2, &
  nq1 = 11, nq2 = 12 * 2 ** nq1, nwords1 = ndp1 / mpdpw + 2, &
  nwords2 = ndp2 / mpdpw + 2)
double precision dplog10q, d1, d2, second, tm0, tm1, tm2
type (mp_real) err, quadtsp, fun01, fun02, fun03, fun04, fun05, fun06, &
  fun07, fun08, fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, &
  t1, t2, t3, t4, wk(-1:nq2), xk(-1:nq2)
type (mp_real) mppic, mpl02, x1, x2
external quadtsp, fun01, fun02, fun03, fun04, fun05, fun06, fun07, fun08, &
  fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, second
save

!   Compute pi and log(2) to high precision (nwords2).

mppic = mppi (nwords2)
mpl02 = mplog2 (nwords2)

write (6, 1) nq1, ndp1, ndp2, neps1, neps2
1 format ('quadtsp test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,' Epsilon2 =',i6)

!   Initialize quadrature tables wk and xk (weights and 1 - abscissas).
!   See documentation for details.

tm0 = second ()
call initqts (nq1, nq2, nwords1, neps2, wk, xk)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite itervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 8
write (6, 3) levq
3 format ('Level =',i4)
tm0 = second ()
t1 = quadtsp (fun01, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
4 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (0.25d0, nwords1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1
5 format ('Actual error =',f10.6,'x10^',i6)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun02, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = (mppic - 2.d0 + 2.d0 * mpl02) / 12.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (0.5d0 * mppic, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun03, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * (exp (0.5d0 * mppic) - 1.d0)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun04, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 5.d0 * mppic**2 / 96.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 15)
15 format (/&
  'Continuous functions on finite itervals, but non-diff at an endpoint'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 7
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun05, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (-4.d0, nwords1) / 9.d0
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun06, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.25d0 * mppic
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint.'//&
  'Problem 7: Int_0^1 sqrt(t)/sqrt(1-t^2) dt = 2*sqrt(pi)*gamma(3/4)/gamma(1/4)')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun07, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 2.d0 * sqrt (mppic) * gamma (mpreal (0.75d0, nwords1)) &
  / gamma (mpreal (0.25d0, nwords1))
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 7
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun08, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (2.d0, nwords1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (0.5d0 * mppic, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun09, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = -0.5d0 * mppic * mpl02
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (0.5d0 * mppic, nwords2)
levq = 8
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun10, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * mppic * sqrt (mpreal (2.d0, nwords1))
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 21)
21 format (/&
  'Functions on an infinite interval'//&
  'Problem 11: Int_0^inf 1/(1+t^2) dt = pi/2')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 9
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun11, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = 0.5d0 * mppic
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 10
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun12, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 11
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun13, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = sqrt (0.5d0 * mppic)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 24)
24 format (/&
  'Oscillatory functions on an infinite interval.'//&
  'Problem 14: Int_0^inf e^(-t)*cos(t) dt = 1/2')
x1 = mpreal (0.d0, nwords2)
x2 = mpreal (1.d0, nwords2)
levq = 11
write (6, 3) levq
tm0 = second ()
t1 = quadtsp (fun14, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 4) tm1 - tm0
call mpwrite (6, ndp1 + 20, ndp1, t1)
t2 = mpreal (0.5d0, nwords1)
call decmdq (t2 - t1, d1, n1)
write (6, 5) d1, n1

write (6, 99) tm2
99 format ('Total CPU time =',f12.6)

stop
end


!   The next 16 routines are the function definitions.

function fun01 (t, nwords1, nwords2)

!   fun01(t) = t * log(1+t)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun01, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun01 = t1 * log (1.d0 + t1)
return
end

function fun02 (t, nwords1, nwords2)

!   fun02(t) = t^2 * arctan(t)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun02, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun02 = t1 ** 2 * atan (t1)
return
end

function fun03 (t, nwords1, nwords2)

!   fun03(t) = e^t * cos(t)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun03, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun03 = exp(t1) * cos(t1)
return
end

function fun04 (t, nwords1, nwords2)

!   fun04(t) = arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2))

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun04, t1, t2
type (mp_real) t

t1 = mpreal (t, nwords1)
t2 = sqrt (2.d0 + t1**2)
fun04 = atan(t2) / ((1.d0 + t1**2) * t2)
return
end

function fun05 (t, nwords1, nwords2)

!    fun05(t) = sqrt(t)*log(t)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun05, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun05 = sqrt (t1) * log (t1)
return
end

function fun06 (t, nwords1, nwords2)

!    fun06(t) = sqrt(1-t^2)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun06, t1, t2
type (mp_real) t

t1 = mpreal (t, nwords1)
t2 = mpreal (1.d0 - t**2, nwords1)
fun06 = sqrt (t2)
return
end

function fun07 (t, nwords1, nwords2)

!   fun07(t) = sqrt (t) / sqrt(1-t^2)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun07, t1, t2
type (mp_real) t

!   The subtraction to compute t2 must be performed using high precision
!   (nwords2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwords1)
t2 = mpreal (1.d0 - t, nwords1)
fun07 = sqrt (t1) / sqrt (t2 * (1.d0 + t1))
return
end

function fun08 (t, nwords1, nwords2)

!   fun08(t) = log(t)^2

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun08, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun08 = log (t1) ** 2
return
end

function fun09 (t, nwords1, nwords2)

!   fun09(t) = log (cos (t))

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun09, pi, t1, t2, t3
type (mp_real) t

t1 = mpreal (t, nwords1)
pi = mppi (nwords2)
t3 = mpreal (0.25d0 * pi, nwords1)
t2 = mpreal (0.5d0 * pi - t, nwords1)

if (t1 < t3) then
  fun09 = log (cos (t1))
else
  fun09 = log (sin (t2))
endif
return
end

function fun10 (t, nwords1, nwords2)

!   fun10(t) = sqrt(tan(t))

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun10, pi, t1, t2, t3, t4
type (mp_real) t

t1 = mpreal (t, nwords1)
pi = mppi (nwords2)
t3 = mpreal (0.25d0 * pi, nwords1)
t2 = mpreal (0.5d0 * pi - t, nwords1)

if (t1 < t3) then
  fun10 = sqrt (tan (t1))
else
  fun10 = 1.d0 / sqrt (tan (t2))
endif
return
end

function fun11 (t, nwords1, nwords2)

!   fun11(t) = 1/(u^2(1+(1/u-1)^2)) = 1/(1 - 2*u + u^2)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun11, t1
type (mp_real) t

t1 = mpreal (t, nwords1)
fun11 = 1.d0 / (1.d0 - 2.d0 * t1 + 2.d0 * t1 ** 2)
return
end

function fun12 (t, nwords1, nwords2)

!   fun12(t) = e^(-(1/t-1)) / sqrt(t^3 - t^4)

use mpmodule
implicit none
integer nwords1, nwords2
type (mp_real) fun12, t1, t2
type (mp_real) t

!   The subtraction to compute t2 must be performed using high precision
!   (nwords2), but after the subtraction its low precision value is fine.

t1 = mpreal (t, nwords1)
t2 = mpreal (1.d0 - t, nwords1)
fun12 = exp (1.d0 - 1.d0/t1) / sqrt (t1 ** 3 * t2)
return
end

function fun13 (t, nwords1, nwords2)

!   fun13(t) = e^(-(1/t-1)^2/2) / t^2

use mpmodule
implicit none
integer nwords1, nwords2
double precision dig1
type (mp_real) fun13, t1, t2
type (mp_real) t

t1 = mpreal (t, nwords1)
dig1 = mpdpw * (nwords1 - 2)
if (t1 >  mpreald (0.25d0 / sqrt (dig1), nwords1)) then
  t2 = 1.d0 / t1 - 1.d0
  fun13 = exp (-0.5d0 * t2 ** 2) / t1 ** 2
else
  fun13 = mpreal (0.d0, nwords1)
endif
return
end

function fun14 (t, nwords1, nwords2)

!   fun14(t) = e^(-(1/t-1)) * cos (1/t-1) / t^2

use mpmodule
implicit none
integer nwords1, nwords2
double precision dig1
type (mp_real) fun14, t1, t2
type (mp_real) t

t1 = mpreal (t, nwords1)
dig1 = mpdpw * (nwords1 - 2)
if (t1 > mpreald (0.25d0 / dig1, nwords1)) then
  t2 = 1.d0 / t1 - 1.d0
  fun14 = exp (-t2) * cos (t2) / t1 ** 2
else
  fun14 = mpreal (0.d0, nwords1)
endif
return
end


subroutine initqts (nq1, nq2, nwords1, neps2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = tanh (pi/2*sinh(t)).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
!   for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
!   and (xk(i), i = -1 to n), where n = xk(-1).   The array x_k contains 
!   1 minus the abscissas; the wk array contains the weights at these abscissas.

!   Both initqts and quadtsp are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   These computations are performed in low precision (nwords1), although
!   computations continue until weights are smaller than the high-precision
!   epsilon eps2 (namely 10^neps2). 

!   David H Bailey   30 Mar 2015

use mpmodule
implicit none
integer i, ierror, iprint, izx, j, k, k1, ndebug, neps2, nq1, nq2, nwords1
double precision dwk(0:nq2), dplog10q
parameter (iprint = 1024, ndebug = 2)
type (mp_real) eps2, h, p2, spi, t1, t2, t3, t4, u1, u2, wk(-1:nq2), xk(-1:nq2)
external dplog10q

write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')

eps2 = mpreal (10.d0, nwords1) ** neps2
p2 = 0.5d0 * mppi (nwords1)
h = mpreal (0.5d0 ** nq1, nwords1)
wk(-1) = mpreal (dble (nq1), nwords1)
izx = 0

!$omp parallel shared (dwk, wk, xk) firstprivate (izx)
!$omp do private (k, t1, t2, t3, t4, u1, u2) schedule (static, 1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  if (izx == 0) then
    t1 = mpreal (dble (k) * h, nwords1)

!   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wk(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

    t2 = exp (t1)
    u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
    u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
    t3 = exp (u2)
    t4 = 0.5d0 * (t3 + 1.d0 / t3)
    xk(k) = 1.d0 / (t3 * t4)
    wk(k) = u1 / t4 ** 2
    dwk(k) = dplog10q (wk(k))
    if (wk(k) < eps2) izx = 1
  else
    xk(k) = mpreal (0.d0, nwords1)
    wk(k) = mpreal (0.d0, nwords1)
    dwk(k) = -999999.d0
  endif
enddo
!$omp end do
!$omp end parallel

do k = 0, nq2
  if (dwk(k) == -999999.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Second argument (nq2) is too small; value =',i8)
stop

100 continue

xk(-1) = mpreal (dble (k), nwords1)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif

return
end


function quadtsp (fun, x1, x2, levq, nq1, nq2, nwords1, nwords2, neps1, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   [x1, x2] with a target tolerance of 10^neps1.  The quadrature level is levq.
!   nq2 is the size of the wk and xk arrays, which must first be initialized
!   by calling initqts. The function fun is not evaluated at the endpoints
!   x1 and x2.

!   Both initqts and quadtsp are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   For some functions, it is important that the endpoints x1 and x2 be
!   computed to high precision (nwords2 words) in the calling program, that
!   these high-precision values be passed to quadtsp (where scaled abscissas
!   are calculated to high precision), and that the function definition
!   itself uses these high-precision scaled abscissas in any initial
!   subtractions or other sensitive operations involving the input argument.
!   Otherwise the accuracy of the quadrature result might only be half as
!   high as it otherwise could be.  See the function definitions of fun06,
!   fun07, fun09 and fun10 for examples on how this is done.  Otherwise the
!   function evaluation can and should be performed with low precision
!   (nwords1 words) for faster run times.  The two precision levels (nwords1
!   and nwords2) are computed below based on ndp1 and ndp2, which are
!   specified by the user in the main program.

!   David H Bailey  30 Mar 2015

use mpmodule
implicit none
integer i, ierror, ip(0:100), iprint, iz1, iz2, izx, j, k, k1, k2, k3, levq, &
  n, ndebug, nds, neps1, nq1, nq2, nqq1, nwords1, nwords2
parameter (iprint = 1024, izx = 2, ndebug = 2)
logical log1, log2
double precision d1, d2, d3, d4, dplog10q
type (mp_real) c10, eps1, eps2, epsilon1, err, fun, h, &
  quadtsp, s1, s2, s3, tsum1, tsum2, tsum3, t1, t2, t3, &
  tw1(0:nq2), tw2(0:nq2), twi1, twi2, twmx, wk(-1:nq2), xk(-1:nq2)
type (mp_real) ax, bx, x1, x2, xki, xt1, xx1, xx2
external fun, dplog10q

!  These first two lines are performed in high precision (nwords2).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in low precision (nwords1).

epsilon1 = mpreal (10.d0, nwords1) ** neps1
s1 = mpreal (0.d0, nwords1)
s2 = mpreal (0.d0, nwords1)
h = mpreal (1.d0, nwords1)
c10 = mpreal (10.d0, nwords1)

if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadtsp: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

k = levq
h = mpreal (0.5d0 ** k, nwords1)
k1 = ip(nqq1-k)
k2 = ip(nqq1-k+1)
k3 = ip(nqq1-k+2)
iz1 = 0
iz2 = 0

!   Evaluate function at level k in x, avoiding unnecessary computation.

!$omp parallel shared (wk, xk) firstprivate (iz1, iz2)
!$omp do private (i, xki, xt1, xx1, xx2, log1, log2, t1, t2) schedule (static, 1)

do i = 0, n, k1
  if (mod (i, iprint) == 0) write (6, *) i, n

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwords2) to ensure full accuracy in the quadrature
!   results, even though the abscissas xk(i) were computed in low precision.

  xki = xk(i)
  xt1 = 1.d0 - mpreal (xki, nwords2)
  xx1 = - ax * xt1 + bx
  xx2 = ax * xt1 + bx
  log1 = xx1 > x1
  log2 = xx2 < x2

!   The remaining computations are performed in low precision (nwords1).

  if (log1 .and. iz1 < izx) then
    t1 = fun (xx1, nwords1, nwords2)
    tw1(i) = t1 * wk(i)
    if (abs (tw1(i)) < epsilon1) then
      iz1 = iz1 + 1
    else
      iz1 = 0
    endif
  else
    t1 = mpreal (0.d0, nwords1)
    tw1(i) = mpreal (0.d0, nwords1)
  endif

  if (i > 0 .and. log2 .and. iz2 < izx) then
    t2 = fun (xx2, nwords1, nwords2)
    tw2(i) = t2 * wk(i)
    if (abs (tw2(i)) < epsilon1) then
      iz2 = iz2 + 1
    else
      iz2 = 0
    endif
  else
    t2 = mpreal (0.d0, nwords1)
    tw2(i) = mpreal (0.d0, nwords1)
  endif
enddo
!$omp end do
!$omp end parallel

!   Tsum1, computed below, is the sum of tw1 and tw2 from the loop above.
!   Tsum2 is the sum of tw1 and tw2, only taking every other value.
!   Tsum3 is the sum of tw1 and tw2, only taking every fourth value.
!   Twmx is the largest absolute value of tw1 and tw2.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

tsum1 = mpreal (0.d0, nwords1)
tsum2 = mpreal (0.d0, nwords1)
tsum3 = mpreal (0.d0, nwords1)
twmx = mpreal (0.d0, nwords1)

do i = 0, n, k1
  tsum1 = tsum1 + tw1(i) + tw2(i)
  if (mod (i, k2) == 0) tsum2 = tsum2 + tw1(i) + tw2(i)
  if (mod (i, k3) == 0) tsum3 = tsum3 + tw1(i) + tw2(i)
  twmx = max (twmx, abs (tw1(i)), abs (tw2(i)))
enddo

do i = k1 * (n / k1), 0, -k1
  if (abs (tw1(i)) > 0.d0) goto 100
enddo

i = 0

100 continue

twi1 = abs (tw1(i))

do i = k1 * (n / k1), 0, -k1
  if (abs (tw2(i)) > 0.d0) goto 110
enddo

i = 0

110 continue

twi2 = abs (tw2(i))

!   Compute s1 = current integral approximation and err = error estimate.

s1 = mpreal (ax, nwords1) * h * tsum1
s2 = 2.d0 * mpreal (ax, nwords1) * h * tsum2
s3 = 4.d0 * mpreal (ax, nwords1) * h * tsum3 
eps1 = twmx * epsilon1
eps2 = max (twi1, twi2)
d1 = dplog10q (abs (s1 - s2))
d2 = dplog10q (abs (s1 - s3))
d3 = dplog10q (eps1) - 1.d0
d4 = dplog10q (eps2) - 1.d0

if (k <= 2) then
  err = mpreal (1.d0, nwords1)
elseif (d1 .eq. -999999.d0) then
  err = mpreal (0.d0, nwords1)
else
  err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
endif

!   Output current integral approximation and error estimate, to 56 dp.

if (ndebug >= 2) then
  write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2 format ('quadtsp: Iteration',i3,' of',i3,'; est error = 10^',i5, &
    '; approx value =')
  call mpwrite (6, 80, 60, s1)
endif
if (k >= 3 .and. err < eps1) goto 140
if (k >= 3 .and. err < eps2) goto 120

write (6, 3) nint (dplog10q (abs (err))), nq1
3 format ('quadtsp: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
goto 140

120 continue

write (6, 4) nint (dplog10q (abs (err))), nwords2
4 format ('quadtsp: Estimated error = 10^',i5/&
  'Increase secondary prec (nwords2) for greater accuracy. Current value =',i4)

140 continue

quadtsp = s1
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (mp_real) a

call mpmdc (a%mpr, da, ia, mpwds)
if (da .eq. 0.d0) then
  dplog10q = -999999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmdq (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

call mpmdc (a%mpr, da, ia, mpwds)
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end
