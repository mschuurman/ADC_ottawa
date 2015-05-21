!*****************************************************************************

!  MPFUN2015: Extra-high precision subroutines
!  Version date:  19 May 2015

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2015 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs (mostly changes to type
!    statements).  All basic arithmetic operations and transcendental functions
!    are supported.  Advanced techniques, including FFT-based multiplication and
!    quadratically convergent transcendental algorithms, are employed.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.  There are NO global shared variables (except
!    static compile-time data), and NO initialization is necessary unless
!    extremely high precision (> 19,500 digits) is required.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:
   
!    David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.
 
!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains subroutines to perform special functions.  At present,
!    it includes the Gamma function and the BesselJ functions.  Additional
!    functions will be added as they are completed.

module mpfune
use mpfuna
use mpfunb
use mpfunc
use mpfund
contains

subroutine mpbessj (anu, t, z, mpnw)

!   This evaluates the function BesselJ (ANU, T).  ANU must be nonnegative and
!   not greater than 10^6 (this limit can be adjusted below).  To compensate
!   for an unsually large amount of internal cancelation in these formulas, all
!   computations are performed to 3*mpnw/2 words precision.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     dasy = factor used to decide if asymptic series is used; default = 25.
!     anumx = upper limit of anu argument; default = 1000.

implicit none
integer i, itrmx, j, mpnw, mpnw1, ndp, nu, n1
double precision dasy, anu(0:mpnw+5), anumx, t(0:mpnw+5), z(0:mpnw+5), &
  t0(0:3*mpnw/2+5), t1(0:3*mpnw/2+5), t2(0:3*mpnw/2+5), t3(0:3*mpnw/2+5), &
  t4(0:3*mpnw/2+5), t5(0:3*mpnw/2+5), t6(0:3*mpnw/2+5)
parameter (itrmx = 100000, dasy = 25.d0, anumx = 1.d6)
  
! End of declaration

if (mpnw < 4 .or. anu(0) < abs (anu(2)) + 4 .or. t(0) < mpnw + 4 .or. &
  t(0) < abs (t(2)) + 4 .or. z(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSJ: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (anu(2) <= 0.d0 .or. anu(3) > 0.d0 .or. &
  (anu(3) == 0.d0 .and. anu(4) > anumx)) then
  write (6, 2) anumx, anu
2 format ('*** MPBESSJ: First argument must be >= 0 and <=', &
  f10.0/'Value =',1p,d25.15)
  call mpabrt (65)
endif

mpnw1 = 3 * mpnw / 2
t0(0) = mpnw1 + 6
t1(0) = mpnw1 + 6
t2(0) = mpnw1 + 6
t3(0) = mpnw1 + 6
t4(0) = mpnw1 + 6
t5(0) = mpnw1 + 6
t6(0) = mpnw1 + 6

!   Select either the direct or the asymptotic series.

if (t(3) < 0.d0 .or. t(3) == 0.d0 .and. t(4) < dasy * (mpnw - 2)) then
  t2(1) = mpnw1
  t2(2) = 1.d0
  t2(3) = 0.d0
  t2(4) = 1.d0
  t2(5) = 0.d0
  t2(6) = 0.d0
  call mpadd (anu, t2, t0, mpnw1)
  call mpgamma (t0, t1, mpnw1)
  call mpdiv (t2, t1, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  call mpeq (t1, t0, mpnw1)
  call mpmul (t, t, t3, mpnw1)
  call mpmuld (t3, 0.25d0, t2, mpnw1)

  do i = 1, itrmx
    call mpmul (t1, t2, t3, mpnw1)
    call mpdivd (t3, dble (i), t4, mpnw1)
    call mpdmc (dble (i), 0, t5, mpnw1)
    call mpadd (anu, t5, t6, mpnw1)
    call mpdiv (t4, t6, t1, mpnw1)
    t1(2) = - t1(2)
    call mpadd (t0, t1, t3, mpnw1)
    call mpeq (t3, t0, mpnw1)
    if (t1(2) == 0.d0 .or. t1(3) < t0(3) - mpnw1) goto 100
  enddo

  write (6, 3)
3 format ('*** MPBESSJ: loop overflow 1')
  call mpabrt (66)

100 continue

  call mpmuld (t, 0.5d0, t1, mpnw1)
  call mppower (t1, anu, t2, mpnw1)  
  call mpmul (t0, t2, t3, mpnw1)
  call mpeq (t3, t0, mpnw1)
else
  t0(1) = mpnw1
  t0(2) = 1.d0
  t0(3) = 0.d0
  t0(4) = 1.d0
  t0(5) = 0.d0
  t0(6) = 0.d0
  t1(1) = mpnw1
  t1(2) = 0.d0
  t1(3) = 0.d0
  t1(4) = 0.d0
  t2(1) = mpnw1
  t2(2) = 1.d0
  t2(3) = 0.d0
  t2(4) = 1.d0
  t2(5) = 0.d0
  t2(6) = 0.d0
  call mpmul (anu, anu, t3, mpnw1)
  call mpmuld (t3, 4.d0, t5, mpnw1)
  
  do i = 1, itrmx
    call mpdmc (dble (2*i - 1), 0, t4, mpnw1)
    call mpmul (t4, t4, t6, mpnw1)
    call mpsub (t5, t6, t4, mpnw1)
    call mpmul (t2, t4, t3, mpnw1)
    call mpdivd (t3, 8.d0 * dble (i), t4, mpnw1)
    call mpdiv (t4, t, t2, mpnw1)
    if (mod (i, 2) == 0) then
      call mpeq (t2, t3, mpnw1)
      if (mod (i, 4) == 2)  t3(2) = - t3(2)
      call mpadd (t0, t3, t4, mpnw1)
      call mpeq (t4, t0, mpnw1)
    else
      call mpeq (t2, t3, mpnw1)
      if (mod (i, 4) == 3) t3(2) = - t3(2)
      call mpadd (t1, t3, t4, mpnw1)
      call mpeq (t4, t1, mpnw1)
    endif
    if (t2(2) == 0.d0 .or. (t2(3) < t0(3) - mpnw1 .and. t2(3) < t1(3) - mpnw1)) &
      goto 110
  enddo
  
  write (6, 4)
4 format ('*** MPBESSJ: loop overflow 2')
  call mpabrt (66)
  
110 continue

  call mpeq (mppicon, t2, mpnw1)
  call mpmul (t2, anu, t4, mpnw1)
  call mpmuld (t4, 0.5d0, t3, mpnw1)
  call mpsub (t, t3, t4, mpnw1)
  call mpmuld (t2, 0.25d0, t3, mpnw1)
  call mpsub (t4, t3, t5, mpnw1)
  call mpcssnr (t5, t3, t4, mpnw1)
  call mpmul (t3, t0, t5, mpnw1)
  call mpmul (t4, t1, t6, mpnw1)
  call mpsub (t5, t6, t3, mpnw1)
  t4(1) = mpnw1
  t4(2) = 1.d0
  t4(3) = 0.d0
  t4(4) = 2.d0
  t4(5) = 0.d0
  t4(6) = 0.d0
  call mpmul (t2, t, t5, mpnw1)
  call mpdiv (t4, t5, t6, mpnw1)
  call mpsqrt (t6, t4, mpnw1)
  call mpmul (t4, t3, t0, mpnw1)
endif

call mproun (t0, mpnw)
call mpeq (t0, z, mpnw)

return
end subroutine

subroutine mpgamma (t, z, mpnw)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
integer i, itrmx, j, k, mpnw, mpnw1, ndp, neps, nt, n1, n2, n3
double precision alpha, al2, dmax, d1, d2, d3
parameter (al2 = 0.69314718055994530942d0, dmax = 1d8, itrmx = 100000)
double precision t(0:mpnw+5), z(0:mpnw+5), f1(0:8), sum1(0:mpnw+6), &
  sum2(0:mpnw+6), tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6)
  
  double precision mpdpw
  parameter (mpdpw = 14.449439791871d0)

! End of declaration

if (mpnw < 4 .or. t(0) < mpnw + 4 .or. t(0) < abs (t(2)) + 4 .or. &
  z(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPGAMMA: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (t(2) == 0 .or. t(3) > 0 .or. (t(3) == 0 .and. t(4) > dmax) .or. &
  (t(2) < 0.d0 .and. t(3) == 0.d0 .and. abs (t(2)) == 1.d0)) then
  write (6, 2) dmax
2 format ('*** MPGAMMA: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call mpabrt (65)
endif

mpnw1 = mpnw + 1
f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0
sum1(0) = mpnw + 7
sum2(0) = mpnw + 7
tn(0) = mpnw + 7
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
t6(0) = mpnw + 7

!   Find the integer and fractional parts of t.

call mpinfr (t, t2, t3, mpnw1)

if (t3(2) == 0.d0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)

  do i = 2, nt - 1
    call mpmuld (t1, dble (i), t2, mpnw1)
    call mpeq (t2, t1, mpnw1)
  enddo

  call mproun (t1, mpnw)
  call mpeq (t1, z, mpnw)
  goto 120
elseif (t(2) > 0.d0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)
  call mpeq (t3, tn, mpnw1)
  
  do i = 1, nt
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpsub (t, t4, t5, mpnw1)
    call mpmul (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call mpsub (f1, t, t4, mpnw1)
  call mpinfr (t4, t3, t5, mpnw1)
  call mpmdc (t3, d3, n3, mpnw1)
  nt = d3 * 2.d0 ** n3

  call mpeq (f1, t1, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)
    
  do i = 0, nt - 1
!    t1 = t1 / (t + dble (i))
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpadd (t, t4, t5, mpnw1)
    call mpdiv (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the nearest integer
!   value, so that d2 = 0.25 * alpha^2 can be calculated exactly in DP.

alpha = anint (0.5d0 * mpnbt * al2 * (mpnw1 + 1))
d2 = 0.25d0 * alpha**2

call mpeq (tn, t2, mpnw1)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum1, mpnw1)

!   Evaluate the series with t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum1, t3, t4, mpnw1)
  call mpeq (t4, sum1, mpnw1)
  if (t3(2) == 0.d0 .or. t3(3) < sum1(3) - mpnw1) goto 100
enddo

write (6, 3)
3 format ('*** MPGAMMA: loop overflow 1')
call mpabrt (67)

100 continue

call mpeq (tn, t2, mpnw1)
t2(2) = - t2(2)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum2, mpnw1)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum2, t3, t4, mpnw1)
  call mpeq (t4, sum2, mpnw1)
  if (t3(2) == 0.d0 .or. t3(3) < sum2(3) - mpnw1) goto 110
enddo

write (6, 4)
4 format ('*** MPGAMMA: loop overflow 4')
call mpabrt (67)

110 continue

!   Compute sqrt (mppic * sum1 / (tn * sin (mppic * tn) * sum2)) 
!   and (alpha/2)^tn terms.

call mpeq (mppicon, t2, mpnw1)
call mpmul (t2, tn, t3, mpnw1)
call mpcssnr (t3, t4, t5, mpnw1)
call mpmul (t5, sum2, t6, mpnw1)
call mpmul (tn, t6, t5, mpnw1)
call mpmul (t2, sum1, t3, mpnw1)
call mpdiv (t3, t5, t6, mpnw1)
t6(2) = - t6(2)
call mpsqrt (t6, t2, mpnw1)

call mpdmc (0.5d0 * alpha, 0, t3, mpnw1)
call mplog (t3, t4, mpnw1)
call mpmul (tn, t4, t5, mpnw1)
call mpexp (t5, t6, mpnw1)
call mpmul (t2, t6, t3, mpnw1)

call mpmul (t1, t3, t4, mpnw1)

!   Round to mpnw words precision.

call mproun (t4, mpnw)
call mpeq (t4, z, mpnw)

120 continue

return
end subroutine mpgamma

end module mpfune
