program tpslqm3

!   David H Bailey   09 Apr 2015

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2015 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program should be compiled with Version 2 of MPFUNF (in mpfunf2.f90).

!   This program demonstrates pslqm3, which performs the one-level multipair
!   PSLQ algorithm on an input vector.  A variety of sample input vectors can
!   be generated as inputs to pslqm3, as given in the parameters below.  The
!   pslqm3 routine is best suited for large relations (above degree 100). For
!   smaller relations pslqm1 or pslqm2 should be used instead.  pslqm3 also
!   includes a facility to read and write a restart file -- see below.
!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslqm3 routine is 100% THREAD SAFE -- all requisite parameters and
!   arrays are passed through subroutine arguments.

!   These parameters are set in the parameter statement below.  Some other
!   parameters are set in other routines.

!     idb = Debug level (0 - 3).
!     itr = Number of trials when KQ = 1, 2 or 3.
!     n   = Integer relation vector length.
!     kq  = 0 for the algebraic case [1, al, al^2, ... al^(n-1)], where
!           al = 3^(1/kr) - 2^(1/ks).  Set n = kr * ks + 1 and itr = 1.
!         = 1 for testing algebraic relations of a number read from a file.
!         = 2 for testing additive relations of numbers read from a file.
!         = 3 for testing multiplicative relations of numbers read from a file.
!         = 4 for custom input.
!     kr  = Degree of root of 3 when kq = 0.
!     ks  = Degree of root of 2 when kq = 0.
!     ndp = Full precision, in digits.  Must not exceed mpipl in mpmod90.f.
!           It is recommended that this be set to at least 30 digits more
!           than the expected detection level of the relation.
!     ndpm = Medium precision for certain calculations; default = 120.
!     nep = Log10 of epsilon for detections.  Normally set to 10 - ndp or so.
!           Must not exceed the accuracy of input data.
!     rb  = Log10 of max size (Euclidean norm) of acceptable relation.

use mpmodule
implicit none
double precision d1, d2, rb, second, tm, tm0, tm1
integer i, idb, iq, itr, j, k, kps, kcs, kq, kr, ks, kss, n, ndp, ndpm, nep, &
  nwp, nwpm
parameter (idb = 2, itr = 1, kq = 0, kr = 7, ks = 8, &
  n = kr * ks + 1, ndp = 720, ndpm = 120, nep = -700, rb = 100.d0)
character*1 chr1(100)
character*40 nam(n), namx
integer lnm(n)
double precision r1(n)
type (mp_real) al, eps, r(n), x(n)

!   Check to see if default precision is high enough.

if (ndp > mpipl) then
  write (6, '("Default precision level is too low.")')
  stop
endif

nwp = int (ndp / mpdpw + 2)  
nwpm = int (ndpm / mpdpw + 2)
eps = mpreal (10.d0, nwp) ** nep
write (6, 1) itr, n, kq, kr, ks, rb, ndp, nep
1 format ('PSLQM3 Test Program'/ &
  'No. trials = ',i3,3x,'n =',i4,3x,'kq =',i2/ &
  'kr =',i2,3x,'ks =',i2,3x,'rb =',1p,d12.4/ &
  'Requested precision level ndp =',i6,' digits'/ &
  'Epsilon level nep = ',i6)
  
!  Calculate the word equivalents of ndp an ndpm.

tm = mpreald (0.d0, nwp)
kps = 0
kcs = 0
if (kq .eq. 1 .or. kq .eq. 2 .or. kq .eq. 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

do k = 1, itr
  write (6, 2) k
2 format (/'Start trial', i3)

  if (kq .eq. 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

    al = mpnrt (mpreald (3.d0, nwp), kr) - mpnrt (mpreald (2.d0, nwp), ks)
  elseif (kq .eq. 1) then

!   Read an algebraic constant from a file.

    call mpread (11, al, nwp)
  elseif (kq .eq. 2) then

!   Read constants from a file for additive test.

    do i = 1, n
      call mpread (11, al, nwp)
      x(i) = al
      write (namx, '(''con'',i3.3)') i
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  elseif (kq .eq. 3) then

!   Read constants from a file for multiplicative test.

    do i = 1, n
      call mpread (11, al, nwp)
      x(i) = log (al)
      write (namx, '(''log(con'',i3.3,'')'')') i
      nam(i) = namx(1:11)
      lnm(i) = 11
    enddo
  elseif (kq .eq. 4) then

!   Produce X vector by a custom scheme.

  endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

  if (kq .eq. 0 .or. kq .eq. 1) then
    x(1) = mpreald (1.d0, nwp)
    nam(1) = '1'
    lnm(1) = 1

    do i = 2, n
      x(i) = al * x(i-1)
      write (namx, '(''al^'',i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  endif

!   Perform relation search.

  tm0 = second ()
  call pslqm3 (idb, n, nwp, nwpm, rb, eps, x, iq, r)
  tm1 = second ()

!   Output relation, if one was found.

  if (iq .ne. 0) then
    tm = tm + (tm1 - tm0)
    write (6, 4)
4   format ('Recovered relation: 0 =')

    do i = 1, n
      if (r(i) .ne. 0.d0) then
        call mpfform (r(i), 64, 0, chr1)
        write (6, '(''+'',64a1,'' * '',a)') (chr1(j), j = 1, 64), &
          nam(i)(1:lnm(i))
      endif
    enddo

!   Check if original relation was recovered.

    if (kq .eq. 1) then
      d1 = 0.d0
      d2 = 0.d0

      do i = 1, n
        d1 = d1 + abs (dble (r(i)) - r1(i))
        d2 = d2 + abs (dble (r(i)) + r1(i))
      enddo

      if (d1 .eq. 0.d0 .or. d2 .eq. 0.d0) then
        kcs = kcs + 1
        write (6, 5)
5       format ('Original relation recovered.')
      else
        kps = kps + 1
      endif
    endif
  endif
  write (6, 6) tm1 - tm0
6 format ('CPU Time =',f12.4)
enddo

if (itr .gt. 1) then
  write (6, 7) kps, kcs
7 format (/'No. partial successes =',i5/'No. complete successes =',i5)
  kss = kps + kcs
  if (kss .ne. 0) then
    tm = tm / kss
    write (6, 8) tm
8   format ('Ave. CPU time of PS or CS runs =',f10.4)
  endif
endif
stop
end

!------------------------------

!   The following code performs the three-level, multi-pair PSLQ algorithm.
!   David H. Bailey    16 Mar 2015

subroutine pslqm3 (idb, n, nwp, nwpm, rb, eps, x, iq, r)

!   Arguments are as follows:
!     Name Type    Description
!     idb   int*4   Debug flag (0-3).
!     n     int*4   Dimension of input vector.
!     nwp   int*4   Full precision in words.
!     ndpm  int*4   Medium pecision level in words.
!     rb    real*8  Log10 of max size (Euclidean norm) of acceptable relation.
!     eps   mp_real Tolerance for relation detection.
!     x     mp_real Input mp_real vector.
!     iq    int*4   Output flag: 0 (unsuccessful) or 1 (successful).
!     r     mp_real Output integer relation vector.

!   The following parameters are set in this routine:
!     ipi   int*4  Iteration print interval when idb >= 2.
!     ipm   int*4  Iteration check interval for MP iterations.
!     itm   int*4  Maximum iteration count.  Run is aborted when this is exceeded.
!     nrs   int*4  0 (no restart), 1 (start and write restart file) or
!                  2 (read restart file and write restart file).
!     nsq   int*4   Size of tables used in itermpw.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, its, izd, izm, izmm, j, j1, n, &
  n1, n2, n3, n4, nrs, nsq, nwp, nwpm
parameter (ipi = 500, ipm = 10, itm = 10000000, nrs = 0, nsq = 8)
double precision d1, d2, d3, d4, dplog10m, rb, second, &
  tm0, tm1, times(6)
integer ip(n), ir(n), is(n)
double precision da(n,n), db(n,n), dh(n,n), dq(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsyq(n,nsq), dt(n,n), dy(n), dsy(n)
type (mp_real) b(n,n), h(n,n), t(n), r(n), x(n), y(n), eps, t1, t2
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wq(n), wsyq(n,nsq), wt(n,n), &
  wy(n), boundm, wn, w1, w2
external dplog10m, second, boundm

!   Initialize.

if (idb .ge. 2) write (6, 1) n
1 format ('PSLQM3 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
wn = mpreald (0.d0, nwpm)

if (idb .ge. 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')

!   If nrs is 1, read or write 

if (nrs .ge. 1) open (12, file = 'pslqm3.rst', form = 'unformatted')
if (nrs .le. 1) then
  do i = 1, 6
    times(i) = 0.d0
  enddo

  tm0 = second ()
  call initmp (idb, n, nwp, ip, dt, b, h, t, x, y)
  tm1 = second ()
  times(1) = tm1 - tm0
else
  rewind 12
  read (12) it
  read (12) times
  read (12) y
  read (12) b
  read (12) h
  rewind 12
endif

100 continue

!   Initialize MPM arrays from MP arrays.

if (idb .ge. 3) write (6, 3) it
3 format ('Iteration',i7,3x,'MPM initialization')
tm0 = second ()
call initmpm (idb, it, n, nsq, nwpm, ip, dt, wa, wb, wh, wy, wsyq, h, y, izmm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)
if (izmm .ne. 0) goto 160

110 continue

!   Initialize DP arrays from MPM arrays.

if (idb .ge. 3) write (6, 4) it
4 format ('Iteration',i7,3x,'Start DP iterations')
call initdp (idb, it, n, nsq, nwpm, da, db, dh, dt, dy, dsyq, wh, wy, izd)
call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
its = it

if (izd .eq. 0) then

!   Perform an LQ decomposition on DH, prior to DP iterations.

  call lqdp (n, n - 1, dh)

!   Perform DP iterations.

120 continue

  it = it + 1
  if (idb .ge. 4 .or. idb .ge. 2 .and. mod (it, ipi) .eq. 0) write (6, 5) it
5 format ('Iteration',i7)
  tm0 = second ()
  call iterdp (idb, it, n, nsq, ip, ir, is, da, db, dh, dq, dsyq, &
    dt, dy, imq, izd)
  tm1 = second ()
  times(3) = times(3) + (tm1 - tm0)

  if (izd .eq. 0) then
    if (mod (it - its, ipm) .eq. 0) then
      call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
      its = it
    endif
    goto 120
  else

    if (izd .eq. 2) then

!   DP iteration was aborted -- revert to previous data.

      it = its
      call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
    endif

!   Update the MPM arrays from the DP arrays.

    if (idb .ge. 3) write (6, 6) it
6   format ('Iteration',i7,3x,'MPM update')
    tm0 = second ()
    call updtmpm (idb, it, n, nwpm, ip, dt, da, db, wa, wb, wh, wt, wy, izmm)
    tm1 = second ()
    times(4) = times(4) + (tm1 - tm0)

    if (izmm .eq. 0) then
      if (izd .eq. 2) then
        goto 130
      else
        goto 110
      endif
    elseif (izmm .eq. 1) then

!   Update the MP arrays from the MPM arrays.

      if (idb .ge. 2) write (6, 7) it
7     format ('Iteration',i7,3x,'MP update')
      tm0 = second ()
      call updtmp (idb, it, n, nwp, ip, dt, wa, wb, eps, b, h, t, y, izm)
      tm1 = second ()
      times(5) = times(5) + (tm1 - tm0)

      if (nrs .ge. 1) then
        rewind 12
        write (12) it
        write (12) times
        write (12) y
        write (12) b
        write (12) h
        rewind 12
      endif

!   Compute norm bound.

      call lqmpm (n, n - 1, nwpm, wh)
      w1 = boundm (n, nwpm, wh)
      call decmdm (w1, d3, n3)
      wn = max (wn, w1)
      call decmdm (wn, d4, n4)
      if (idb .ge. 2) then
        write (6, 8) it, d3, n3, d4, n4
8       format ('Iteration',i7,3x,'Norm bound =',f11.6,'e',i5,3x, &
          'Max. bound =',f11.6,'e',i5)
      endif
      if (it .gt. itm) then
        if (idb .ge. 1) write (6, 9) itm
9       format ('Iteration limit exceeded',i7)
        goto 160
      endif
      if (dplog10m (wn) .gt. rb) then
        if (idb .ge. 1) write (6, 10) rb
10      format ('Norm bound limit exceeded.',1p,d15.6)
        goto 160
      endif

      if (izm .eq. 0) then
        if (izd .eq. 2) then
          goto 130
        else
          goto 100
        endif
      elseif (izm .eq. 1) then
        goto  150
      elseif (izm .eq. 2) then
        goto 160
      endif
    elseif (izmm .eq. 2) then
      goto 160
    endif
  endif
endif

130 continue

!   Perform an LQ decomposition on wh, prior to MPM iterations.

if (idb .ge. 2) write (6, 11) it
11 format ('Iteration',i7,3x,'Start MPM iterations')
tm0 = second ()
call lqmpm (n, n - 1, nwpm, wh)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Perform MPM iterations.

140 continue

it = it + 1
if (idb .ge. 2) write (6, 12) it
12 format ('Iteration',i7)
tm0 = second ()
call itermpm (idb, it, n, nsq, nwpm, ip, ir, is, dt, wa, wb, wh, wq, wsyq, &
  wt, wy, imq, izmm)
tm1 = second ()
times(6) = times(6) + (tm1 - tm0)

if (izmm .eq. 0) then
  if (mod (it - its, ipm) .eq. 0) then

!   Check to see if DP iterations can be resumed.

    call checkm (idb, it, n, nwpm, wy, izd)
    if (izd .eq. 0) then
      if (idb .ge. 2) write (6, 13) it
13    format ('Iteration',i7,3x,'Return to DP iterations')
      goto 110
    endif
  endif
  goto 140
elseif (izmm .eq. 1) then

!   Update the MP arrays from the MPM arrays.

  if (idb .ge. 2) write (6, 14) it
14 format ('Iteration',i7,3x,'MP update')
  tm0 = second ()
  call updtmp (idb, it, n, nwp, ip, dt, wa, wb, eps, b, h, t, y, izm)
  tm1 = second ()
  times(5) = times(5) + (tm1 - tm0)
  if (izm .eq. 0) then
    goto 100
  elseif (izm .eq. 1) then
    goto 150
  elseif (izm .eq. 2) then
    goto 160
  endif
elseif (izmm .eq. 2) then
  goto 160
endif

!   A relation has been detected.  Output the final norm bound and other info.

150 continue

tm0 = second ()
t1 = mpreald (1.d300, nwp)
t2 = mpreald (0.d0, nwp)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
enddo

!   The norm calculation here is performed in medium percision.

w1 = mpreal (0.d0, nwpm)

do i = 1, n
  w1 = w1 + mpreal (r(i), nwpm) ** 2
enddo

w1 = sqrt (w1)
call decmdm (w1, d3, n3)

!   Output the final norm bound and other info.

if (idb .ge. 1) then
  call lqmpm (n, n - 1, nwpm, wh)
  w2 = boundm (n, nwpm, wh)
  wn = max (wn, w2)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmdm (wn, d4, n4)
  write (6, 15) it, d1, n1, d2, n2, d4, n4
15 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6/'Max. bound =',f11.6,'e',i6)
  write (6, 16) j1, d3, n3, d1, n1
16 format ('Index of relation =',i4,3x,'Norm =',f11.5,'e',i5,3x, &
  'Residual =',f11.6,'e',i6)
endif

if (dplog10m (w1) .le. rb) then
  iq = 1
else
  if (idb .ge. 2) write (6, 17)
17 format ('Relation is too large.')
endif

160 continue

if (idb .ge. 2) write (6, 18) times
18 format ('CPU times:'/(6f12.2))

return
end

!------------------------------

!   First-level subroutines.

subroutine checkm (idb, it, n, nwpm, wy, izd)

!   Checks wh and wy to see if DP iterations can be resumed.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, idb, it, izd, n, n1, nwpm
double precision d1, deps
parameter (deps = 1d-10)
type (mp_real) wy(n), t1, t2, t3, t4

izd = 0
t1 = mpreald (1.d300, nwpm)
t2 = mpreald (0.d0, nwpm)

!   Find the min and max absolute value in the wy vector.

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

t3 = t1 / t2
t4 = mpreald (deps, nwpm)
if (t3 .lt. t4) then
  if (idb .ge. 2) then
    call decmdm (t3, d1, n1)
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'checkm: Small min/max ratio in wy =',f11.6, &
    'e',i4)
  endif
  izd = 1
endif

return
end

subroutine initdp (idb, it, n, nsq, nwpm, da, db, dh, dt, dy, dsyq, wh, wy, izd)

!   This re-initializes the DP arrays from the MPM arrays.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, idb, it, izd, j, n, n1, nsq, nwpm
double precision da(n,n), db(n,n), dh(n,n), dt(n,n), dy(n), dsyq(n,nsq), d1, deps
type (mp_real) wh(n,n), wy(n), t1, t2, t3, t4
parameter (deps = 1d-10)

izd = 0
t1 = mpreald (1.d300, nwpm)
t2 = mpreald (0.d0, nwpm)

!   Find the min and max absolute value in the wy vector.

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

!   If the dynamic range of the wy vector is too great, set the izd flag and
!   abort this initialization.

t3 = t1 / t2
t4 = mpreald (deps, nwpm)
if (t3 .lt. t4) then
  if (idb .ge. 3) then
    call decmdm (t3, d1, n1)
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'initdp: Small min/max ratio in wy =',f11.6, &
    'e',i4)
  endif
  izd = 1
  goto 100
endif

!   Set dy to be the scaled wy vector.

t1 = 1.d0 / t2

do i = 1, n
  dy(i) = t1 * wy(i)
enddo

!   Find the maximum absolute value of the wh matrix diagonals.

t2 = mpreald (0.d0, nwpm)

do j = 1, n - 1
  t2 = max (t2, abs (wh(j,j)))
enddo

!   Set dh to be the scaled wh matrix.

t1 = 1.d0 / t2

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = t1 * wh(i,j)
  enddo
enddo

!   Set da and db to the identity.

do j = 1, n
  do i = 1, n
    da(i,j) = 0.d0
    db(i,j) = 0.d0
  enddo

  da(j,j) = 1.d0
  db(j,j) = 1.d0
enddo

!   Zero the dt array.

do j = 1, n
  do i = 1, n
    dt(i,j) = 0.d0
  enddo
enddo

!   Zero the dsyq array.

do j = 1, nsq
  do i = 1, n
    dsyq(i,j) = 0.d0
  enddo
enddo

if (idb .ge. 4) then
  write (6, 2)
2 format ('initdp: Scaled dy vector:')
  call matoutdp (1, n, dy)
  write (6, 3)
3 format ('initdp: Scaled dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

100 continue

return
end

subroutine initmp (idb, n, nwp, ix, dx, b, h, s, x, y)

!   Initializes MP arrays at the beginning.
!   This is performed wiht full precision.

use mpmodule
implicit none
integer i, idb, j, n, nwp
integer ix(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), s(n), x(n), y(n), t1

if (idb .ge. 4) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = mpreald (0.d0, nwp)
  enddo

  b(j,j) = mpreald (1.d0, nwp)
enddo

t1 = mpreald (0.d0, nwp)

!   Compute the x vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = mpreald (0.d0, nwp)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo

if (idb .ge. 4) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine initmpm (idb, it, n, nsq, nwpm, ix, dx, wa, wb, wh, wy, wsyq, h, y, izmm)

!   This re-initializes the MPM arrays from the MP arrays.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, idb, it, izmm, j, n, n1, nsq, ntl, nwpm
parameter (ntl = 72)
double precision d1
integer ix(n)
double precision dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wy(n), wsyq(n,nsq), t1, t2, t3, teps
type (mp_real) h(n,n), y(n)

teps = mpreald (2.d0, nwpm) ** (ntl - nwpm * mpnbt)
izmm = 0
t1 = mpreald (1.d300, nwpm)
t2 = mpreald (0.d0, nwpm)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (mpreal (y(i), nwpm))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

!   If the dynamic range of the y vector is too great, then abort.

t3 = t1 / t2
if (t3 .lt. teps) then
  if (idb .ge. 2) then
    call decmd (t3, d1, n1)
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'initmpm: Small min/max ratio in y =',f11.6, &
    'e',i4/ 'Run aborted.')
  endif
  izmm = 1
  goto 100
endif

!   Set wy to the scaled y vector.

t1 = 1.d0 / t2

do i = 1, n
  wy(i) = t1 * mpreal (y(i), nwpm)
enddo

!   Set wh to h.

do j = 1, n - 1
  do i = 1, n
    wh(i,j) = mpreal (h(i,j), nwpm)
  enddo
enddo

!   Set wa and wb to the identity.

do j = 1, n
  do i = 1, n
    wa(i,j) = mpreald (0.d0, nwpm)
    wb(i,j) = mpreald (0.d0, nwpm)
  enddo

  wa(j,j) = mpreald (1.d0, nwpm)
  wb(j,j) = mpreald (1.d0, nwpm)
enddo

!   Zero the wsyq array.

do j = 1, nsq
  do i = 1, n
    wsyq(i,j) = mpreald (0.d0, nwpm)
  enddo
enddo

if (idb .ge. 4) then
  write (6, 2)
2 format ('initmpm: wy vector:')
  call matoutmdm (1, n, ix, dx, wy)
  write (6, 3)
3 format ('initmpm: Factored wh matrix:')
  call matoutmdm (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

subroutine iterdp (idb, it, n, nsq, ip, ir, is, da, db, dh, dq, dsyq, &
  dt, dy, imq, izd)

!   This performs one iteration of the PSLQ algorithm using DP arithmetic.

implicit none
integer i, idb, ii, ij, im, im1, imq, it, izd, j, j1, j2, k, mpr, mq, n, nsq
double precision deps, gam, t1, t2, t3, t4, tmx1, tmx2
integer ip(n), ir(n), is(n)
double precision da(n,n), db(n,n), dh(n,n), dq(n), dsyq(n,nsq), &
  dt(n,n), dy(n)
parameter (tmx1 = 1.d13, tmx2 = 2.d0**52, deps = 1.d-14)

izd = 0
mpr = nint (0.4d0 * n)
gam = sqrt (4.d0 / 3.d0)

!   Compute dq vector = {gam^i * |dh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  dq(i) = gam ** i * abs (dh(i,i))
enddo

call qsortdp (n - 1, dq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest dq(i).

do i = 1, n
  is(i) = 0
enddo

if (imq .eq. 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii .eq. 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of dy, and rows of da, db and dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = dy(im)
  dy(im) = dy(im1)
  dy(im1) = t1

  do i = 1, n
    t1 = da(im,i)
    da(im,i) = da(im1,i)
    da(im1,i) = t1
    t1 = db(im,i)
    db(im,i) = db(im1,i)
    db(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = dh(im,i)
    dh(im,i) = dh(im1,i)
    dh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im .le. n - 2) then
    t1 = dh(im,im)
    t2 = dh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = dh(i,im)
      t4 = dh(i,im1)
      dh(i,im) = t1 * t3 + t2 * t4
      dh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on dh, using the diagonal scheme.  Multipliers are
!   saved in the dt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      dh(ij,j) = dh(ij,j) - dt(ij,k) * dh(k,j)
    enddo

    dt(ij,j) = anint (dh(ij,j) / dh(j,j))
    dh(ij,j) = dh(ij,j) - dt(ij,j) * dh(j,j)
  enddo
enddo

!   Update dy, using the dt array.  Find min absolute value of dy.

t1 = abs (dy(n))

do j = 1, n - 1
  do i = j + 1, n
    dy(j) = dy(j) + dt(i,j) * dy(i)
  enddo

  t1 = min (t1, abs (dy(j)))
enddo

!   Update da and db, using the dt array.  Find the max absolute value of
!   da and db entries as they are calculated (not merely at the end).

t2 = 0.d0

do k = 1, n
  t3 = 0.d0

  do j = 1, n - 1
    do i = j + 1, n
      da(i,k) = da(i,k) - dt(i,j) * da(j,k)
      db(j,k) = db(j,k) + dt(i,j) * db(i,k)
      t3 = max (t3, abs (da(i,k)), abs (db(j,k)))
    enddo
  enddo

  dq(k) = t3
enddo

do k = 1, n
  t2 = max (t2, dq(k))
enddo

if (t1 .le. deps) then
  if (idb .ge. 3) write (6, 1) it, t1
1 format ('Iteration',i7,3x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 .gt. tmx1 .and. t2 .le. tmx2) then
  if (idb .ge. 3) write (6, 2) it, t2
2 format ('Iteration',i7,3x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 .gt. tmx2) then
  if (idb .ge. 2) write (6, 3) it, t2
3 format ('Iteration',i7,3x,'iterdp: Very large value in da or db =',1pd15.6)
  izd = 2
  return
endif

!   Compare the dy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (dy(i) - dsyq(i,j)))
  enddo

  if (t1 .le. deps) then
    if (idb .ge. 2) write (6, 4) it, j
4   format ('Iteration',i7,3x,'iterdp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector dy in the table dsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  dsyq(i,k) = dy(i)
enddo

if (idb .ge. 4) then
  write (6, 5)
5 format ('iterdp: Updated dy:')
  call matoutdp (1, n, dy)
  write (6, 6)
6 format ('iterdp: Updated da matrix:')
  call matoutdp (n, n, da)
  write (6, 7)
7 format ('iterdp: Updated db matrix:')
  call matoutdp (n, n, db)
  write (6, 8)
8 format ('iterdp: Updated dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine itermpm (idb, it, n, nsq, nwpm, ip, ir, is, dx, wa, wb, wh, wq, &
  wsyq, wt, wy, imq, izmm)

!   This performs one iteration of the PSLQM algorithm using MPM arithmetic.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izmm, j, j1, j2, k, mpr, mq, n, n1, &
  nsq, ntl, nwpm
parameter (ntl = 72)
double precision d1
type (mp_real) gam, t1, t2, t3, t4, teps, tmx1, tmx2
integer ip(n), ir(n), is(n)
double precision dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wq(n), wsyq(n,nsq), &
  wt(n,n), wy(n)

teps = mpreald (2.d0, nwpm) ** (ntl - nwpm * mpnbt)
tmx1 = mpreald (2.d0, nwpm) ** (nwpm * mpnbt - ntl)
tmx2 = mpreald (2.d0, nwpm) ** (nwpm * mpnbt)
izmm = 0
mpr = nint (0.4d0 * n)
gam = sqrt (mpreald (4.d0, nwpm) / mpreald (3.d0, nwpm))

!   Compute wq vector = {gam^i * |wh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  wq(i) = gam ** i * abs (wh(i,i))
enddo

call qsortmpm (n - 1, nwpm, wq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest wq(i).

do i = 1, n
  is(i) = 0
enddo

if (imq .eq. 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii .eq. 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of wy, and rows of wa, wb and wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = wy(im)
  wy(im) = wy(im1)
  wy(im1) = t1

  do i = 1, n
    t1 = wa(im,i)
    wa(im,i) = wa(im1,i)
    wa(im1,i) = t1
    t1 = wb(im,i)
    wb(im,i) = wb(im1,i)
    wb(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = wh(im,i)
    wh(im,i) = wh(im1,i)
    wh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in wh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im .le. n - 2) then
    t1 = wh(im,im)
    t2 = wh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = wh(i,im)
      t4 = wh(i,im1)
      wh(i,im) = t1 * t3 + t2 * t4
      wh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on wh, using the diagonal scheme.  Multipliers are
!   saved in the wt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      wh(ij,j) = wh(ij,j) - wt(ij,k) * wh(k,j)
    enddo

    wt(ij,j) = anint (wh(ij,j) / wh(j,j))
    wh(ij,j) = wh(ij,j) - wt(ij,j) * wh(j,j)
  enddo
enddo

!   Update wy, using the wt array.  Find min absolute value of wy.

t1 = abs (wy(n))

do j = 1, n - 1
  do i = j + 1, n
    wy(j) = wy(j) + wt(i,j) * wy(i)
  enddo

  t1 = min (t1, abs (wy(j)))
enddo

!   Update wa and wb, using the wt array.  Find the max absolute value of
!   wa and wb entries as they are calculated (not merely at the end).

t2 = mpreald (0.d0, nwpm)

do k = 1, n
  t3 = mpreald (0.d0, nwpm)

  do j = 1, n - 1
    do i = j + 1, n
      wa(i,k) = wa(i,k) - wt(i,j) * wa(j,k)
      wb(j,k) = wb(j,k) + wt(i,j) * wb(i,k)
      t3 = max (t3, abs (wa(i,k)), abs (wb(j,k)))
    enddo
  enddo

  wq(k) = t3
enddo

do k = 1, n
  t2 = max (t2, wq(k))
enddo

if (t1 .le. teps) then
  if (idb .ge. 2) then
    call decmdm (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermpm: Small value in wy =',f11.6,'e',i4)
  endif
  izmm = 1
endif

if (t2 .gt. tmx1 .and. t2 .le. tmx2) then
  if (idb .ge. 2) then
    call decmdm (t2, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i7,3x,'itermpm: Large value in wa or wb =', &
      f11.6,'e',i4)
  endif
  izmm = 1
elseif (t2 .gt. tmx2) then
  if (idb .ge. 1) then
    call decmdm (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('itermpm: Very large value in wa or wb =',f11.6,'e',i4/ &
      'Run aborted.')
  endif
  izmm = 2
  goto 200
endif

!   Compare the wy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = mpreald (0.d0, nwpm)

  do i = 1, n
    t1 = max (t1, abs (wy(i) - wsyq(i,j)))
  enddo

  if (t1 .le. teps) then
    if (idb .ge. 2) write (6, 4) it, j
4   format ('Iteration',i7,3x,'itermpm: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector wy in the table wsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  wsyq(i,k) = wy(i)
enddo

if (idb .ge. 4) then
  write (6, 5)
5 format ('itermpm: Updated wy:')
  call matoutmdm (1, n, ip, dx, wy)
  write (6, 6)
6 format ('itermpm: Updated wa matrix:')
  call matoutmdm (n, n, ip, dx, wa)
  write (6, 7)
7 format ('itermpm: Updated wb matrix:')
  call matoutmdm (n, n, ip, dx, wb)
  write (6, 8)
8 format ('itermpm: Updated wh matrix:')
  call matoutmdm (n, n - 1, ip, dx, wh)
endif

200 continue

return
end

subroutine lqdp (n, m, dh)

!   This performs an LQ decomposition on the DP matrix dh.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.

implicit none
integer i, j, l, lup, m, ml, n
double precision dh(n,m), nrmxl, one, t, zero

zero = 0.d0
one = 1.d0
lup = min (m,n)

!   Perform the householder reduction of dh.

do l = 1, lup
  if (l .eq. m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + dh(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl .eq. zero) go to 270
  if (dh(l,l) .ne. zero) nrmxl = sign (nrmxl, dh(l,l))
  t = one / nrmxl

  do i = 0, ml
    dh(l,l+i) = t * dh(l,l+i)
  enddo

  dh(l,l) = one + dh(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + dh(l,l+i) * dh(j,l+i)
    enddo

    t = - t / dh(l,l)

    do i = 0, ml
      dh(j,l+i) = dh(j,l+i) + t * dh(l,l+i)
    enddo
  enddo

!   Save the transformation.

  dh(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero dh above the diagonal.

do j = 1, m
  do i = 1, j - 1
    dh(i,j) = 0.d0
  enddo
enddo

return
end

subroutine lqmpm (n, m, nwpm, h)

!   This performs an LQ decomposition on the MPM matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, j, l, lup, m, ml, n, nwpm
type (mp_real) h(n,m), nrmxl, one, t, zero

zero = mpreald (0.d0, nwpm)
one = mpreald (1.d0, nwpm)
lup = min (m,n)

!   Perform the householder reduction of h.

do l = 1, lup
  if (l .eq. m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + h(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl .eq. zero) go to 270
  if (h(l,l) .ne. zero) nrmxl = sign (nrmxl, h(l,l))
  t = one / nrmxl

  do i = 0, ml
    h(l,l+i) = t * h(l,l+i)
  enddo

  h(l,l) = one + h(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + h(l,l+i) * h(j,l+i)
    enddo

    t = - t / h(l,l)

    do i = 0, ml
      h(j,l+i) = h(j,l+i) + t * h(l,l+i)
    enddo
  enddo

!   Save the transformation.

  h(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero h above the diagonal.

do j = 1, m
  do i = 1, j - 1
    h(i,j) = mpreald (0.d0, nwpm)
  enddo
enddo

return
end

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case dp iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.

implicit none
integer i, j, n
double precision da(n,n), db(n,n), dh(n,n), dy(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsy(n)

do i = 1, n
  dsy(i) = dy(i)
enddo

do j = 1, n
  do i = 1, n
    dsa(i,j) = da(i,j)
    dsb(i,j) = db(i,j)
  enddo
enddo

do j = 1, n - 1
  do i = 1, n
    dsh(i,j) = dh(i,j)
  enddo
enddo

return
end

subroutine updtmp (idb, it, n, nwp, ix, dx, wa, wb, eps, b, h, t, y, izm)

!   This update the MP arrays from the MPM arrays.
!   This is performed in full precision.

use mpmodule
implicit none 
integer i, i1, idb, it, izm, n, n1, n2, ntl, nwp
parameter (ntl = 72)
double precision d1, d2
type (mp_real) eps, t1, t2, teps
integer ix(n)
double precision dx(n)
type (mp_real) wa(n,n), wb(n,n)
type (mp_real) b(n,n), h(n,n), t(n), y(n)

teps = mpreald (2.d0, nwp) ** ntl * eps
izm = 0
t1 = mpreald (1.d300, nwp)
t2 = mpreald (0.d0, nwp)

!   Update y with wb.

call mxmm (n, 1, nwp, wb, y, t)

i1 = 0
do i = 1, n
  if (abs (y(i)) < t1) then
    i1 = i
    t1 = abs (y(i))
  endif
  t2 = max (t2, abs (y(i)))
enddo

if (idb .ge. 2) then
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1  format ('Iteration',i7,3x,'updtmp: Min, max of y =',f11.6,'e',i6, &
    f11.6,'e',i6)
endif

!   Update b with wb.

call mxmm (n, n, nwp, wb, b, t)

!   Update h with wa.

call mxmm (n, n - 1, nwp, wa, h, t)

!   Find the largest entry of b in the same row as the smallest y.

t2 = mpreald (0.d0, nwp)

do i = 1, n
  t2 = max (t2, abs (b(i1,i)))
enddo

if (t1 .le. t2 * teps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1) 
    write (6, 2) it, d1, n1
2   format ('Iteration',i7,3x,'updtmp: Small value in y =',f11.6,'e',i6)
  endif
  if (t1 .le. t2 * eps) then
    izm = 1
  else
    if (idb .ge. 1) write (6, 3) it
3   format ('Iteration',i7,3x,'updtmp: Precision exhausted.')
    izm = 2
  endif
endif

if (idb .ge. 4) then
  write (6, 4)
4 format ('updtmp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 5)
5 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 6)
6 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine updtmpm (idb, it, n, nwpm, ix, dx, da, db, wa, wb, wh, wt, wy, izmm)

!   This updates the MPM arrays from the DP arrays.
!   This is performed in medium precision.

use mpmodule
implicit none 
integer i, idb, it, izmm, j, n, n1, n2, ntl, nwpm
parameter (ntl = 72)
double precision d1, d2
type (mp_real) t1, t2, teps, tmx1, tmx2
integer ix(n)
double precision da(n,n), db(n,n), dx(n)
type (mp_real) wa(n,n), wb(n,n), wh(n,n), wt(n), wy(n)

teps = mpreald (2.d0, nwpm) ** (ntl - nwpm * mpnbt)
tmx1 = mpreald (2.d0, nwpm) ** (nwpm * mpnbt - ntl)
tmx2 = mpreald (2.d0, nwpm) ** (nwpm * mpnbt)
izmm = 0
t1 = mpreald (1.d300, nwpm)
t2 = mpreald (0.d0, nwpm)

!   Update wy with db.

call mxmdm (n, 1, nwpm, db, wy, wt)

do i = 1, n
  t1 = min (t1, abs (wy(i)))
  t2 = max (t2, abs (wy(i)))
enddo

if (idb .ge. 3) then
  call decmdm (t1, d1, n1)
  call decmdm (t2, d2, n2)
  write (6, 1) it, d1, n1, d2, n2
1 format ('Iteration',i7,3x,'updtmpm: Min, max of wy =',f11.6,'e',i4, &
     f11.6,'e',i4)
endif
if (t1 .le. teps) then
  if (idb .ge. 2) then
    call decmdm (t1, d1, n1)
    write (6, 2) it, d1, n1
2   format ('Iteration',i7,3x,'updtmpm: Small value in wy =',f11.6,'e',i4)
  endif
  izmm = 1
endif

!   Update wa with da.

call mxmdm (n, n, nwpm, da, wa, wt)

!   Update wb with db.

call mxmdm (n, n, nwpm, db, wb, wt)
t2 = mpreald (0.d0, nwpm)

do j = 1, n
  do i = 1, n
    t2 = max (t2, abs (wa(i,j)), abs (wb(i,j)))
  enddo
enddo

if (t2 .gt. tmx1 .and. t2 .le. tmx2) then
  if (idb .ge. 2) then
    call decmdm (t2, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i7,3x,'updtmpm: Large value in wa or wb =', &
      f11.6,'e',i4)
  endif
  izmm = 1
elseif (t2 .gt. tmx2) then
  if (idb .ge. 1) then
    call decmdm (t2, d1, n1)
    write (6, 4) it, d1, n1
4   format ('updtmpm: Very large value in wa or wb =',f11.6,'e',i4/ &
      'Run aborted.')
  endif
  izmm = 2
  goto 100
endif

!   Update wh with da.

call mxmdm (n, n - 1, nwpm, da, wh, wt)

if (idb .ge. 4) then
  write (6, 5)
5 format ('updtmpm: Updated wy:')
  call matoutmdm (1, n, ix, dx, wy)
  write (6, 6)
6 format ('updtmpm: Updated wa matrix:')
  call matoutmdm (n, n, ix, dx, wa)
  write (6, 7)
7 format ('updtmpm: Updated wb matrix:')
  call matoutmdm (n, n, ix, dx, wb)
  write (6, 8)
8 format ('updtmpm: Updated wh matrix:')
  call matoutmdm (n, n - 1, ix, dx, wh)
endif

100 continue

return
end

!------------------------------

!   Second- and third-level subroutines.

function boundm (n, nwpm, wh)

!   This computes the norm bound using MPM arithmetic, based on the array wh.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, n, nwpm
type (mp_real) wh(n,n), boundm, t1

t1 = mpreald (0.d0, nwpm)

do i = 1, n - 1
  t1 = max (t1, abs (wh(i,i)))
enddo

boundm = 1.d0 / t1
return
end

function dplog10m (a)

!   For input MPM value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia, mpnw
double precision da, dplog10m
type (mp_real) a

mpnw = a%mpr(1)
call mpmdc (a%mpr, da, ia, mpnw)
if (da .le. 0.d0) then
  dplog10m = -1d300
else
  dplog10m = log10 (abs (da)) + ia * log10 (2.d0)
endif

return
end

subroutine decmd (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib, mpnw
type (mp_real) a
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)

mpnw = a%mpr(1)
call mpmdc (a%mpr, da, ia, mpnw)
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

subroutine decmdm (a, b, ib)

!   For input MPM value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib, mpnw
type (mp_real) a
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)

mpnw = a%mpr(1)
call mpmdc (a%mpr, da, ia, mpnw)
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

function drand ()

!   This routine returns a pseudorandom DP floating number nearly uniformly
!   distributed between 0 and 1 by means of a linear congruential scheme.
!   2^28 pseudorandom numbers with 30 bits each are returned before repeating.

implicit none
double precision drand, f7, r30, sd, t1, t2, t30
parameter (f7 = 78125.d0, r30 = 0.5d0 ** 30, t30 = 2.d0 ** 30)
save sd
data sd/314159265.d0/

t1 = f7 * sd
t2 = aint (r30 * t1)
sd = t1 - t30 * t2
drand = r30 * sd

return
end

subroutine matoutdp (n1, n2, a)

!   This outputs the DP matrix a.

implicit none
integer i, j, n1, n2
double precision a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)
  write (6, 2) (a(i,j), j = 1, n2)
2 format (1p5d15.5)
enddo

return
end

subroutine matoutmd (n1, n2, ix, dx, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
double precision dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i5))
enddo

return
end

subroutine matoutmp (n1, n2, a)

!   This outputs the MP matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer i, j, n1, n2
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, 80, 60, a(i,j))
  enddo
enddo

return
end

subroutine matoutmdm (n1, n2, ix, dx, a)

!   This outputs the MPM matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
double precision dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmdm (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i5))
enddo

return
end

subroutine mxmdm (n1, n2, nwpm, a, b, c)

!  This multiplies the DP square matrix a by the MPM matrix b, and the result
!  is placed in b.  c is a MPM scratch array with n1 cells.  n1, n2 are
!  the matrix dimensions as indicated below.
!  This is performed in medium precision.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwpm
double precision a(n1,n1)
type (mp_real) b(n1,n2), c(n1), t1

do j = 1, n2
  do i = 1, n1
!    call mpdotdm (n1, 1, b(1,j), n1, a(i,1), t1)
    t1 = mpreald (0.d0, nwpm)
    
    do k = 1, n1
      t1 = t1 + mpprodd (b(k,j), a(i,k))
    enddo

    c(i) = t1
  enddo

  do i = 1, n1
    b(i,j) = c(i)
  enddo
enddo

return
end

subroutine mxmm (n1, n2, nwp, a, b, c)

!  This multiplies the MPM square matrix a by the MP matrix b, and the result
!  is placed in b.  c is a MP scratch array with n1 cells.  n1, n2 are
!  the matrix dimensions as indicated below.
!  This is performed in full precision.

use mpmodule
implicit none
integer i, j, k, n1, n2, nwp
type (mp_real) a(n1,n1)
type (mp_real) b(n1,n2), c(n1), t1

do j = 1, n2
  do i = 1, n1
    t1 = mpreal (0.d0, nwp)

    do k = 1, n1
      t1 = t1 + a(i,k) * b(k,j)
    enddo

    c(i) = t1
  enddo

  do i = 1, n1
    b(i,j) = c(i)
  enddo
enddo

return
end

subroutine qsortdp (n, a, ip)

!   This routine sorts the entries of the N-long DP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.

implicit none
integer i, iq, it, j, jq, jz, k, l, n
double precision a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n .eq. 1) return

k = 1
ik(1) = 1
jk(1) = n

130 i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 .lt. a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 .gt. a(ip(l))) goto 180
enddo

j = i
goto 190

180 j = l
if (i .ge. j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue
if (s0 .ge. a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 k = k - 1
jz = 0
if (j .eq. iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 i = i + 1
if (i .eq. jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0) goto 220
if (j - iq .ge. jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 if (k .gt. 0) goto 130

return
end

subroutine qsortmpm (n, nwpm, a, ip)

!   This routine sorts the entries of the N-long MPM vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   This is performed in medium precision.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n, nwpm
type (mp_real) a(n), s0
integer ip(n), ik(50), jk(50)

do i = 1, n
  ip(i) = i
enddo

if (n .eq. 1) return

k = 1
ik(1) = 1
jk(1) = n

130 i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 .lt. a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 .gt. a(ip(l))) goto 180
enddo

j = i
goto 190

180 j = l
if (i .ge. j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue
if (s0 .ge. a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 k = k - 1
jz = 0
if (j .eq. iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 i = i + 1
if (i .eq. jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0) goto 220
if (j - iq .ge. jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 if (k .gt. 0) goto 130

return
end

