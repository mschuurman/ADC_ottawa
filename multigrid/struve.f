
      REAL FUNCTION RCHEVAL(N,A,T)
C
C   This function evaluates a Chebyshev series, using the
C   Clenshaw method with Reinsch modification, as analysed
C   in the paper by Oliver.
C
C   INPUT PARAMETERS
C
C       N - INTEGER - The no. of terms in the sequence
C
C       A - REAL ARRAY, dimension 0 to N - The coefficients of
C           the Chebyshev series
C
C       T - REAL - The value at which the series is to be
C           evaluated
C
C
C   REFERENCES
C
C        "An error analysis of the modified Clenshaw method for
C         evaluating Chebyshev and Fourier series" J. Oliver,
C         J.I.M.A., vol. 20, 1977, pp379-391
C
C
C MACHINE-DEPENDENT CONSTANTS: NONE
C
C
C INTRINSIC FUNCTIONS USED;
C
C    ABS
C
C
C AUTHOR:  Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley ,
C          High St.,
C          PAISLEY,
C          SCOTLAND
C
C
C LATEST MODIFICATION:   21 December , 1992
C
C
      IMPLICIT NONE
      INTEGER I,N
      REAL A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      DATA ZERO,HALF/ 0.0 E 0 , 0.5 E 0 /
      DATA TEST,TWO/ 0.6 E 0 , 2.0 E 0 /
      U1 = ZERO
C
C   If ABS ( T )  < 0.6 use the standard Clenshaw method
C
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         RCHEVAL =  ( U0 - U2 ) / TWO
      ELSE
C
C   If ABS ( T )  > =  0.6 use the Reinsch modification
C
         D1 = ZERO
C
C   T > =  0.6 code
C
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            RCHEVAL =  ( D1 + D2 ) / TWO
         ELSE
C
C   T < =  -0.6 code
C
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            RCHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF
      RETURN
      END

      REAL FUNCTION RSTRVH0(XVALUE)
C
C
C   DESCRIPTION:
C
C      This function calculates the value of the Struve function
C      of order 0, denoted H0(x), for the argument XVALUE, defined
C
C         STRVHO(x) = (2/pi) integral{0 to pi/2} sin(x cos(t)) dt
C
C      H0 also satisfies the second-order equation
C
C                 x*D(Df) + Df + x*f = 2x/pi
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H0 involves the Bessel function
C      of the second kind Y0, there is a problem for large x, since
C      we cannot accurately calculate the value of Y0. An error message 
C      is printed and STRVH0 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH0. The
C               recommended value is such that
C                      ABS(ARRH0(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH0A. The
C               recommended value is such that
C                      ABS(ARRH0A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY0ASP. The
C               recommended value is such that
C                      ABS(AY0ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY0ASQ. The
C               recommended value is such that
C                      ABS(AY0ASQ(NTERM4)) < EPS/100.
C
C      XLOW - The value for which H0(x) = 2*x/pi to machine precision, if
C             abs(x) < XLOW. The recommended value is
C                      XLOW = 3 * SQRT(EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y0 with
C              any reasonable accuracy. An error message is printed and
C              STRVH0 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS and EPSNEG refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   15 JUNE, 1995
C
C
      IMPLICIT NONE
      INTEGER INDSGN,NTERM1,NTERM2,NTERM3,NTERM4
      REAL ARRH0(0:19),ARRH0A(0:20),AY0ASP(0:12),
     1     AY0ASQ(0:13),RCHEVAL,EIGHT,ELEVEN,HALF,H0AS,
     2     ONEHUN,ONE,PIBY4,RT2BPI,SIXTP5,T,THR2P5,TWENTY,
     3     TWOBPI,TWO62,X,XHIGH,XLOW,XMP4,XSQ,XVALUE,
     4     Y0P,Y0Q,Y0VAL,ZERO
      CHARACTER(LEN=6) FNNAME
      CHARACTER(LEN=26) ERRMSG
      DATA FNNAME/'STRVH0'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,ONE/0.0 E 0 , 0.5 E 0 , 1.0 E 0/ 
      DATA EIGHT,ELEVEN/8.0 E 0 , 11.0 E 0/
      DATA TWENTY,ONEHUN/20.0 E 0 , 100.0 E 0/
      DATA SIXTP5,TWO62,THR2P5/60.5 E 0 , 262.0 E 0 , 302.5 E 0/
      DATA PIBY4/0.78539 81633 97448 30962 E 0/
      DATA RT2BPI/0.79788 45608 02865 35588 E 0/
      DATA TWOBPI/0.63661 97723 67581 34308 E 0/
      DATA ARRH0(0)/  0.28696 48739 90132 25740  E    0/
      DATA ARRH0(1)/ -0.25405 33268 16183 52305  E    0/
      DATA ARRH0(2)/  0.20774 02673 93238 94439  E    0/
      DATA ARRH0(3)/ -0.20364 02956 03865 85140  E    0/
      DATA ARRH0(4)/  0.12888 46908 68661 86016  E    0/
      DATA ARRH0(5)/ -0.48256 32815 62226 1202   E   -1/
      DATA ARRH0(6)/  0.11686 29347 56900 1242   E   -1/
      DATA ARRH0(7)/ -0.19811 81356 42418 416    E   -2/
      DATA ARRH0(8)/  0.24899 13851 24212 86     E   -3/
      DATA ARRH0(9)/ -0.24188 27913 78595 0      E   -4/
      DATA ARRH0(10)/ 0.18743 75479 93431        E   -5/
      DATA ARRH0(11)/-0.11873 34607 4362         E   -6/
      DATA ARRH0(12)/ 0.62698 49433 46           E   -8/
      DATA ARRH0(13)/-0.28045 54679 3            E   -9/
      DATA ARRH0(14)/ 0.10769 41205              E  -10/
      DATA ARRH0(15)/-0.35904 793                E  -12/
      DATA ARRH0(16)/ 0.10494 47                 E  -13/
      DATA ARRH0(17)/-0.27119                    E  -15/
      DATA ARRH0(18)/ 0.624                      E  -17/
      DATA ARRH0(19)/-0.13                       E  -18/
      DATA ARRH0A(0)/  1.99291 88575 19923 05515  E    0/
      DATA ARRH0A(1)/ -0.38423 26687 01456 887    E   -2/
      DATA ARRH0A(2)/ -0.32871 99371 23530 50     E   -3/
      DATA ARRH0A(3)/ -0.29411 81203 70340 9      E   -4/
      DATA ARRH0A(4)/ -0.26731 53519 87066        E   -5/
      DATA ARRH0A(5)/ -0.24681 03107 5013         E   -6/
      DATA ARRH0A(6)/ -0.22950 14861 143          E   -7/
      DATA ARRH0A(7)/ -0.21568 22318 33           E   -8/
      DATA ARRH0A(8)/ -0.20303 50648 3            E   -9/
      DATA ARRH0A(9)/ -0.19345 75509              E  -10/
      DATA ARRH0A(10)/-0.18277 3144               E  -11/
      DATA ARRH0A(11)/-0.17768 424                E  -12/
      DATA ARRH0A(12)/-0.16432 96                 E  -13/
      DATA ARRH0A(13)/-0.17156 9                  E  -14/
      DATA ARRH0A(14)/-0.13368                    E  -15/
      DATA ARRH0A(15)/-0.2077                     E  -16/
      DATA ARRH0A(16)/ 0.2                        E  -19/
      DATA ARRH0A(17)/-0.55                       E  -18/
      DATA ARRH0A(18)/ 0.10                       E  -18/
      DATA ARRH0A(19)/-0.4                        E  -19/
      DATA ARRH0A(20)/ 0.1                        E  -19/
      DATA AY0ASP/1.99944 63940 23982 71568  E    0,
     1           -0.28650 77864 70319 58     E   -3,
     2           -0.10050 72797 43762 0      E   -4,
     3           -0.35835 94100 2463         E   -6,
     4           -0.12879 65120 531          E   -7,
     5           -0.46609 48663 6            E   -9,
     6           -0.16937 69454              E  -10,
     7           -0.61852 269                E  -12,
     8           -0.22618 41                 E  -13,
     9           -0.83268                    E  -15,
     X           -0.3042                     E  -16,
     1           -0.115                      E  -17,
     2           -0.4                        E  -19/
      DATA AY0ASQ/1.99542 68138 68286 04092  E    0,
     1           -0.23601 31928 67514 472    E   -2,
     2           -0.76015 38908 50296 6      E   -4,
     3           -0.25610 88714 56343        E   -5,
     4           -0.87502 92185 106          E   -7,
     5           -0.30430 42121 59           E   -8,
     6           -0.10621 42831 4            E   -9,
     7           -0.37737 1479               E  -11,
     8           -0.13213 687                E  -12,
     9           -0.48862 1                  E  -14,
     X           -0.15809                    E  -15,
     1           -0.762                      E  -17,
     2           -0.3                        E  -19,
     3           -0.3                        E  -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/13,8,5,6/
      DATA XLOW,XHIGH/7.324 E -4 , 8388608.0 E 0/
C
C   Start computation
C
      X = XVALUE
      INDSGN = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         INDSGN = -1
      ENDIF
C
C   Error test
C
      IF ( ABS(XVALUE) .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         RSTRVH0 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 11
C
      IF ( X .LE. ELEVEN ) THEN
         IF ( X .LT. XLOW ) THEN
            RSTRVH0 = TWOBPI * X
         ELSE
            T = ( ( X * X ) / SIXTP5 - HALF ) - HALF
            RSTRVH0 = TWOBPI * X * RCHEVAL ( NTERM1 , ARRH0 , T )
         ENDIF
      ELSE      
C
C   Code for abs(x) > 11
C
         XSQ = X * X
         T = ( TWO62 - XSQ ) / ( TWENTY + XSQ )
         Y0P = RCHEVAL ( NTERM3 , AY0ASP , T )
         Y0Q = RCHEVAL ( NTERM4 , AY0ASQ , T ) / ( EIGHT * X )
         XMP4 = X - PIBY4
         Y0VAL = Y0P * SIN ( XMP4 ) - Y0Q * COS ( XMP4 )
         Y0VAL = Y0VAL * RT2BPI / SQRT ( X )
         T = ( THR2P5 - XSQ ) / ( SIXTP5 + XSQ )
         H0AS = TWOBPI * RCHEVAL ( NTERM2 , ARRH0A , T ) / X
         RSTRVH0 = Y0VAL + H0AS
      ENDIF
      IF ( INDSGN .EQ. -1 ) RSTRVH0 = -RSTRVH0
      RETURN
      END

      REAL FUNCTION RSTRVH1(XVALUE)
C
C
C   DESCRIPTION:
C      This function calculates the value of the Struve function
C      of order 1, denoted H1(x), for the argument XVALUE, defined as
C
C                                                                  2
C        STRVH1(x) = (2x/pi) integral{0 to pi/2} sin( x cos(t))*sin t dt
C
C      H1 also satisfies the second-order differential equation
C
C                    2   2                   2            2
C                   x * D f  +  x * Df  +  (x - 1)f  =  2x / pi
C
C      The code uses Chebyshev expansions with the coefficients
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H1 involves the Bessel function
C      of the second kind Y1, there is a problem for large x, since
C      we cannot accurately calculate the value of Y1. An error message 
C      is printed and STRVH1 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH1. The
C               recommended value is such that
C                      ABS(ARRH1(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH1A. The
C               recommended value is such that
C                      ABS(ARRH1A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY1ASP. The
C               recommended value is such that
C                      ABS(AY1ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY1ASQ. The
C               recommended value is such that
C                      ABS(AY1ASQ(NTERM4)) < EPS/100.
C
C      XLOW1 - The value of x, below which H1(x) set to zero, if
C              abs(x)<XLOW1. The recommended value is 
C                      XLOW1 = 1.5 * SQRT(XMIN)
C
C      XLOW2 - The value for which H1(x) = 2*x*x/pi to machine precision, if
C              abs(x) < XLOW2. The recommended value is
C                      XLOW2 = SQRT(15*EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y1 with
C              any reasonable accuracy. An error message is printed and
C              STRVH1 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS, EPSNEG and XMIN refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST REVISION:
C                   18 JUNE, 1995
C
C
      IMPLICIT NONE
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4
      REAL ARRH1(0:17),ARRH1A(0:21),AY1ASP(0:14),
     1     AY1ASQ(0:15),RCHEVAL,EIGHT,FIFTEN,FORTP5,HALF,
     2     H1AS,NINE,ONEHUN,ONE82,RT2BPI,T,THPBY4,
     3     TWENTY,TWOBPI,TW02P5,X,XHIGH,XLOW1,XLOW2,
     4     XM3P4,XSQ,XVALUE,Y1P,Y1Q,Y1VAL,ZERO
      CHARACTER(LEN=6) FNNAME
      CHARACTER(LEN=26) ERRMSG
      DATA FNNAME/'STRVH1'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,EIGHT/0.0 E 0 , 0.5 E 0 , 8.0 E 0/
      DATA NINE,FIFTEN/ 9.0 E 0 , 15.0 E 0 / 
      DATA TWENTY,ONEHUN/ 20.0 E 0 , 100.0 E 0/
      DATA FORTP5,ONE82,TW02P5/40.5 E 0 , 182.0 E 0 , 202.5 E 0/
      DATA RT2BPI/0.79788 45608 02865 35588 E 0/
      DATA THPBY4/2.35619 44901 92344 92885 E 0/
      DATA TWOBPI/0.63661 97723 67581 34308 E 0/
      DATA ARRH1/0.17319 06108 36754 39319  E    0,
     1          -0.12606 91759 13526 72005  E    0,
     2           0.79085 76160 49535 7500   E   -1,
     3          -0.31964 93222 32187 0820   E   -1,
     4           0.80804 05814 04918 834    E   -2,
     5          -0.13600 08206 93074 148    E   -2,
     6           0.16227 14861 98894 71     E   -3,
     7          -0.14423 52451 48592 9      E   -4,
     8           0.99219 52573 4072         E   -6,
     9          -0.54416 28049 180          E   -7,
     X           0.24363 16625 63           E   -8,
     1          -0.90770 71338              E  -10,
     2           0.28592 6585               E  -11,
     3          -0.77169 75                 E  -13,
     4           0.18048 9                  E  -14,
     5          -0.3694                     E  -16,
     6           0.67                       E  -18,
     7          -0.1                        E  -19/
      DATA ARRH1A(0)/  2.01083 50495 14733 79407  E    0/
      DATA ARRH1A(1)/  0.59221 86100 36099 903    E   -2/
      DATA ARRH1A(2)/  0.55274 32269 84141 30     E   -3/
      DATA ARRH1A(3)/  0.52698 73856 31103 6      E   -4/
      DATA ARRH1A(4)/  0.50637 45221 40969        E   -5/
      DATA ARRH1A(5)/  0.49028 73642 0678         E   -6/
      DATA ARRH1A(6)/  0.47635 40023 525          E   -7/
      DATA ARRH1A(7)/  0.46525 86522 83           E   -8/
      DATA ARRH1A(8)/  0.45465 16608 1            E   -9/
      DATA ARRH1A(9)/  0.44724 62193              E  -10/
      DATA ARRH1A(10)/ 0.43730 8292               E  -11/
      DATA ARRH1A(11)/ 0.43568 368                E  -12/
      DATA ARRH1A(12)/ 0.41821 90                 E  -13/
      DATA ARRH1A(13)/ 0.44104 4                  E  -14/
      DATA ARRH1A(14)/ 0.36391                    E  -15/
      DATA ARRH1A(15)/ 0.5558                     E  -16/
      DATA ARRH1A(16)/-0.4                        E  -19/
      DATA ARRH1A(17)/ 0.163                      E  -17/
      DATA ARRH1A(18)/-0.34                       E  -18/
      DATA ARRH1A(19)/ 0.13                       E  -18/
      DATA ARRH1A(20)/-0.4                        E  -19/
      DATA ARRH1A(21)/ 0.1                        E  -19/
      DATA AY1ASP/2.00135 24004 58893 96402 E   0,
     1            0.71104 24159 64619 38    E  -3,
     2            0.36659 77028 23244 9     E  -4,
     3            0.19130 15686 57728       E  -5,
     4            0.10046 91138 9777        E  -6,
     5            0.53040 17425 38          E  -8,
     6            0.28100 88617 6           E  -9,
     7            0.14938 86051             E -10,
     8            0.79578 420               E -12,
     9            0.42523 63                E -13,
     X            0.22719 5                 E -14,
     1            0.12216                   E -15,
     2            0.650                     E -17,
     3            0.36                      E -18,
     4            0.2                       E -19/
      DATA AY1ASQ/5.99065 10947 78881 89116 E   0,
     1           -0.48959 32623 36579 635   E  -2,
     2           -0.23238 32130 70706 26    E  -3,
     3           -0.11447 34723 85767 9     E  -4,
     4           -0.57169 92618 9106        E  -6,
     5           -0.28955 16716 917         E  -7,
     6           -0.14751 33456 36          E  -8,
     7           -0.75965 37378             E -10,
     8           -0.39065 8184              E -11,
     9           -0.20464 654               E -12,
     X           -0.10426 36                E -13,
     1           -0.57702                   E -15,
     2           -0.2550                    E -16,
     3           -0.210                     E -17,
     4            0.2                       E -19,
     5           -0.2                       E -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/11,9,6,7/
      DATA XLOW1,XLOW2,XHIGH/1.63 E -19 , 9.456 E -4 , 8388608 E 0/
C
C   Start computation
C
      X = ABS ( XVALUE )
C
C   Error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         RSTRVH1 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 9
C
      IF ( X .LE. NINE ) THEN
         IF ( X .LT. XLOW1 ) THEN
            RSTRVH1 = ZERO
         ELSE
            XSQ = X * X
            IF ( X .LT. XLOW2 ) THEN
               RSTRVH1 = TWOBPI * XSQ
            ELSE
               T = ( XSQ / FORTP5 - HALF ) - HALF
               RSTRVH1 = TWOBPI * XSQ * RCHEVAL ( NTERM1 , ARRH1 , T )
            ENDIF
         ENDIF
      ELSE      
C
C   Code for abs(x) > 9
C
         XSQ = X * X
         T = ( ONE82 - XSQ ) / ( TWENTY + XSQ )
         Y1P = RCHEVAL ( NTERM3 , AY1ASP , T )
         Y1Q = RCHEVAL ( NTERM4 , AY1ASQ , T ) / ( EIGHT * X)
         XM3P4 = X - THPBY4
         Y1VAL = Y1P * SIN ( XM3P4 ) + Y1Q * COS ( XM3P4 )
         Y1VAL = Y1VAL * RT2BPI / SQRT ( X )
         T = ( TW02P5 - XSQ ) / ( FORTP5 + XSQ )
         H1AS = TWOBPI * RCHEVAL ( NTERM2 , ARRH1A , T )
         RSTRVH1 = Y1VAL + H1AS
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DCHEVAL(N,A,T)
C
C   This function evaluates a Chebyshev series, using the
C   Clenshaw method with Reinsch modification, as analysed
C   in the paper by Oliver.
C
C   INPUT PARAMETERS
C
C       N - INTEGER - The no. of terms in the sequence
C
C       A - DOUBLE PRECISION ARRAY, dimension 0 to N - The coefficients of
C           the Chebyshev series
C
C       T - DOUBLE PRECISION - The value at which the series is to be
C           evaluated
C
C
C   REFERENCES
C
C        "An error analysis of the modified Clenshaw method for
C         evaluating Chebyshev and Fourier series" J. Oliver,
C         J.I.M.A., vol. 20, 1977, pp379-391
C
C
C MACHINE-DEPENDENT CONSTANTS: NONE
C
C
C INTRINSIC FUNCTIONS USED;
C
C    ABS
C
C
C AUTHOR:  Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley ,
C          High St.,
C          PAISLEY,
C          SCOTLAND
C
C
C LATEST MODIFICATION:   21 December , 1992
C
C
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      DATA ZERO,HALF/ 0.0 D 0 , 0.5 D 0 /
      DATA TEST,TWO/ 0.6 D 0 , 2.0 D 0 /
      U1 = ZERO
C
C   If ABS ( T )  < 0.6 use the standard Clenshaw method
C
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         DCHEVAL =  ( U0 - U2 ) / TWO
      ELSE
C
C   If ABS ( T )  > =  0.6 use the Reinsch modification
C
         D1 = ZERO
C
C   T > =  0.6 code
C
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            DCHEVAL =  ( D1 + D2 ) / TWO
         ELSE
C
C   T < =  -0.6 code
C
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            DCHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DSTRVH0(XVALUE)
C
C
C   DESCRIPTION:
C
C      This function calculates the value of the Struve function
C      of order 0, denoted H0(x), for the argument XVALUE, defined
C
C         STRVHO(x) = (2/pi) integral{0 to pi/2} sin(x cos(t)) dt
C
C      H0 also satisfies the second-order equation
C
C                 x*D(Df) + Df + x*f = 2x/pi
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H0 involves the Bessel function
C      of the second kind Y0, there is a problem for large x, since
C      we cannot accurately calculate the value of Y0. An error message 
C      is printed and STRVH0 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH0. The
C               recommended value is such that
C                      ABS(ARRH0(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH0A. The
C               recommended value is such that
C                      ABS(ARRH0A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY0ASP. The
C               recommended value is such that
C                      ABS(AY0ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY0ASQ. The
C               recommended value is such that
C                      ABS(AY0ASQ(NTERM4)) < EPS/100.
C
C      XLOW - The value for which H0(x) = 2*x/pi to machine precision, if
C             abs(x) < XLOW. The recommended value is
C                      XLOW = 3 * SQRT(EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y0 with
C              any reasonable accuracy. An error message is printed and
C              STRVH0 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS and EPSNEG refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   11 JANUARY, 1996
C
C
      IMPLICIT NONE
      INTEGER INDSGN,NTERM1,NTERM2,NTERM3,NTERM4
      DOUBLE PRECISION ARRH0(0:19),ARRH0A(0:20),AY0ASP(0:12),
     1     AY0ASQ(0:13),DCHEVAL,EIGHT,ELEVEN,HALF,H0AS,
     2     ONEHUN,ONE,PIBY4,RT2BPI,SIXTP5,T,THR2P5,TWENTY,
     3     TWOBPI,TWO62,X,XHIGH,XLOW,XMP4,XSQ,XVALUE,
     4     Y0P,Y0Q,Y0VAL,ZERO
      CHARACTER(LEN=6) FNNAME
      CHARACTER(LEN=26) ERRMSG
      DATA FNNAME/'STRVH0'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/ 
      DATA EIGHT,ELEVEN/8.0 D 0 , 11.0 D 0/
      DATA TWENTY,ONEHUN/20.0 D 0 , 100.0 D 0/
      DATA SIXTP5,TWO62,THR2P5/60.5 D 0 , 262.0 D 0 , 302.5 D 0/
      DATA PIBY4/0.78539 81633 97448 30962 D 0/
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARRH0(0)/  0.28696 48739 90132 25740  D    0/
      DATA ARRH0(1)/ -0.25405 33268 16183 52305  D    0/
      DATA ARRH0(2)/  0.20774 02673 93238 94439  D    0/
      DATA ARRH0(3)/ -0.20364 02956 03865 85140  D    0/
      DATA ARRH0(4)/  0.12888 46908 68661 86016  D    0/
      DATA ARRH0(5)/ -0.48256 32815 62226 1202   D   -1/
      DATA ARRH0(6)/  0.11686 29347 56900 1242   D   -1/
      DATA ARRH0(7)/ -0.19811 81356 42418 416    D   -2/
      DATA ARRH0(8)/  0.24899 13851 24212 86     D   -3/
      DATA ARRH0(9)/ -0.24188 27913 78595 0      D   -4/
      DATA ARRH0(10)/ 0.18743 75479 93431        D   -5/
      DATA ARRH0(11)/-0.11873 34607 4362         D   -6/
      DATA ARRH0(12)/ 0.62698 49433 46           D   -8/
      DATA ARRH0(13)/-0.28045 54679 3            D   -9/
      DATA ARRH0(14)/ 0.10769 41205              D  -10/
      DATA ARRH0(15)/-0.35904 793                D  -12/
      DATA ARRH0(16)/ 0.10494 47                 D  -13/
      DATA ARRH0(17)/-0.27119                    D  -15/
      DATA ARRH0(18)/ 0.624                      D  -17/
      DATA ARRH0(19)/-0.13                       D  -18/
      DATA ARRH0A(0)/  1.99291 88575 19923 05515  D    0/
      DATA ARRH0A(1)/ -0.38423 26687 01456 887    D   -2/
      DATA ARRH0A(2)/ -0.32871 99371 23530 50     D   -3/
      DATA ARRH0A(3)/ -0.29411 81203 70340 9      D   -4/
      DATA ARRH0A(4)/ -0.26731 53519 87066        D   -5/
      DATA ARRH0A(5)/ -0.24681 03107 5013         D   -6/
      DATA ARRH0A(6)/ -0.22950 14861 143          D   -7/
      DATA ARRH0A(7)/ -0.21568 22318 33           D   -8/
      DATA ARRH0A(8)/ -0.20303 50648 3            D   -9/
      DATA ARRH0A(9)/ -0.19345 75509              D  -10/
      DATA ARRH0A(10)/-0.18277 3144               D  -11/
      DATA ARRH0A(11)/-0.17768 424                D  -12/
      DATA ARRH0A(12)/-0.16432 96                 D  -13/
      DATA ARRH0A(13)/-0.17156 9                  D  -14/
      DATA ARRH0A(14)/-0.13368                    D  -15/
      DATA ARRH0A(15)/-0.2077                     D  -16/
      DATA ARRH0A(16)/ 0.2                        D  -19/
      DATA ARRH0A(17)/-0.55                       D  -18/
      DATA ARRH0A(18)/ 0.10                       D  -18/
      DATA ARRH0A(19)/-0.4                        D  -19/
      DATA ARRH0A(20)/ 0.1                        D  -19/
      DATA AY0ASP/1.99944 63940 23982 71568  D    0,
     1           -0.28650 77864 70319 58     D   -3,
     2           -0.10050 72797 43762 0      D   -4,
     3           -0.35835 94100 2463         D   -6,
     4           -0.12879 65120 531          D   -7,
     5           -0.46609 48663 6            D   -9,
     6           -0.16937 69454              D  -10,
     7           -0.61852 269                D  -12,
     8           -0.22618 41                 D  -13,
     9           -0.83268                    D  -15,
     X           -0.3042                     D  -16,
     1           -0.115                      D  -17,
     2           -0.4                        D  -19/
      DATA AY0ASQ/1.99542 68138 68286 04092  D    0,
     1           -0.23601 31928 67514 472    D   -2,
     2           -0.76015 38908 50296 6      D   -4,
     3           -0.25610 88714 56343        D   -5,
     4           -0.87502 92185 106          D   -7,
     5           -0.30430 42121 59           D   -8,
     6           -0.10621 42831 4            D   -9,
     7           -0.37737 1479               D  -11,
     8           -0.13213 687                D  -12,
     9           -0.48862 1                  D  -14,
     X           -0.15809                    D  -15,
     1           -0.762                      D  -17,
     2           -0.3                        D  -19,
     3           -0.3                        D  -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/18,18,11,11/
      DATA XLOW,XHIGH/3.1610136D-8,4.50359963D15/
C
C   Start computation
C
      X = XVALUE
      INDSGN = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         INDSGN = -1
      ENDIF
C
C   Error test
C
      IF ( ABS(XVALUE) .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DSTRVH0 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 11
C
      IF ( X .LE. ELEVEN ) THEN
         IF ( X .LT. XLOW ) THEN
            DSTRVH0 = TWOBPI * X
         ELSE
            T = ( ( X * X ) / SIXTP5 - HALF ) - HALF
            DSTRVH0 = TWOBPI * X * DCHEVAL ( NTERM1 , ARRH0 , T )
         ENDIF
      ELSE      
C
C   Code for abs(x) > 11
C
         XSQ = X * X
         T = ( TWO62 - XSQ ) / ( TWENTY + XSQ )
         Y0P = DCHEVAL ( NTERM3 , AY0ASP , T )
         Y0Q = DCHEVAL ( NTERM4 , AY0ASQ , T ) / ( EIGHT * X )
         XMP4 = X - PIBY4
         Y0VAL = Y0P * SIN ( XMP4 ) - Y0Q * COS ( XMP4 )
         Y0VAL = Y0VAL * RT2BPI / SQRT ( X )
         T = ( THR2P5 - XSQ ) / ( SIXTP5 + XSQ )
         H0AS = TWOBPI * DCHEVAL ( NTERM2 , ARRH0A , T ) / X
         DSTRVH0 = Y0VAL + H0AS
      ENDIF
      IF ( INDSGN .EQ. -1 ) DSTRVH0 = -DSTRVH0
      RETURN
      END

      DOUBLE PRECISION FUNCTION DSTRVH1(XVALUE)
C
C
C   DESCRIPTION:
C      This function calculates the value of the Struve function
C      of order 1, denoted H1(x), for the argument XVALUE, defined as
C
C                                                                  2
C        STRVH1(x) = (2x/pi) integral{0 to pi/2} sin( x cos(t))*sin t dt
C
C      H1 also satisfies the second-order differential equation
C
C                    2   2                   2            2
C                   x * D f  +  x * Df  +  (x - 1)f  =  2x / pi
C
C      The code uses Chebyshev expansions with the coefficients
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H1 involves the Bessel function
C      of the second kind Y1, there is a problem for large x, since
C      we cannot accurately calculate the value of Y1. An error message 
C      is printed and STRVH1 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH1. The
C               recommended value is such that
C                      ABS(ARRH1(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH1A. The
C               recommended value is such that
C                      ABS(ARRH1A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY1ASP. The
C               recommended value is such that
C                      ABS(AY1ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY1ASQ. The
C               recommended value is such that
C                      ABS(AY1ASQ(NTERM4)) < EPS/100.
C
C      XLOW1 - The value of x, below which H1(x) set to zero, if
C              abs(x)<XLOW1. The recommended value is 
C                      XLOW1 = 1.5 * SQRT(XMIN)
C
C      XLOW2 - The value for which H1(x) = 2*x*x/pi to machine precision, if
C              abs(x) < XLOW2. The recommended value is
C                      XLOW2 = SQRT(15*EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y1 with
C              any reasonable accuracy. An error message is printed and
C              STRVH1 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS, EPSNEG and XMIN refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST REVISION:
C                   12 JANUARY, 1996
C
C
      IMPLICIT NONE
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4
      DOUBLE PRECISION ARRH1(0:17),ARRH1A(0:21),AY1ASP(0:14),
     1     AY1ASQ(0:15),DCHEVAL,EIGHT,FIFTEN,FORTP5,HALF,
     2     H1AS,NINE,ONEHUN,ONE82,RT2BPI,T,THPBY4,
     3     TWENTY,TWOBPI,TW02P5,X,XHIGH,XLOW1,XLOW2,
     4     XM3P4,XSQ,XVALUE,Y1P,Y1Q,Y1VAL,ZERO
      CHARACTER(LEN=6) FNNAME
      CHARACTER(LEN=26) ERRMSG
      DATA FNNAME/'STRVH1'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,EIGHT/0.0 D 0 , 0.5 D 0 , 8.0 D 0/
      DATA NINE,FIFTEN/ 9.0 D 0 , 15.0 D 0 / 
      DATA TWENTY,ONEHUN/ 20.0 D 0 , 100.0 D 0/
      DATA FORTP5,ONE82,TW02P5/40.5 D 0 , 182.0 D 0 , 202.5 D 0/
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA THPBY4/2.35619 44901 92344 92885 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARRH1/0.17319 06108 36754 39319  D    0,
     1          -0.12606 91759 13526 72005  D    0,
     2           0.79085 76160 49535 7500   D   -1,
     3          -0.31964 93222 32187 0820   D   -1,
     4           0.80804 05814 04918 834    D   -2,
     5          -0.13600 08206 93074 148    D   -2,
     6           0.16227 14861 98894 71     D   -3,
     7          -0.14423 52451 48592 9      D   -4,
     8           0.99219 52573 4072         D   -6,
     9          -0.54416 28049 180          D   -7,
     X           0.24363 16625 63           D   -8,
     1          -0.90770 71338              D  -10,
     2           0.28592 6585               D  -11,
     3          -0.77169 75                 D  -13,
     4           0.18048 9                  D  -14,
     5          -0.3694                     D  -16,
     6           0.67                       D  -18,
     7          -0.1                        D  -19/
      DATA ARRH1A(0)/  2.01083 50495 14733 79407  D    0/
      DATA ARRH1A(1)/  0.59221 86100 36099 903    D   -2/
      DATA ARRH1A(2)/  0.55274 32269 84141 30     D   -3/
      DATA ARRH1A(3)/  0.52698 73856 31103 6      D   -4/
      DATA ARRH1A(4)/  0.50637 45221 40969        D   -5/
      DATA ARRH1A(5)/  0.49028 73642 0678         D   -6/
      DATA ARRH1A(6)/  0.47635 40023 525          D   -7/
      DATA ARRH1A(7)/  0.46525 86522 83           D   -8/
      DATA ARRH1A(8)/  0.45465 16608 1            D   -9/
      DATA ARRH1A(9)/  0.44724 62193              D  -10/
      DATA ARRH1A(10)/ 0.43730 8292               D  -11/
      DATA ARRH1A(11)/ 0.43568 368                D  -12/
      DATA ARRH1A(12)/ 0.41821 90                 D  -13/
      DATA ARRH1A(13)/ 0.44104 4                  D  -14/
      DATA ARRH1A(14)/ 0.36391                    D  -15/
      DATA ARRH1A(15)/ 0.5558                     D  -16/
      DATA ARRH1A(16)/-0.4                        D  -19/
      DATA ARRH1A(17)/ 0.163                      D  -17/
      DATA ARRH1A(18)/-0.34                       D  -18/
      DATA ARRH1A(19)/ 0.13                       D  -18/
      DATA ARRH1A(20)/-0.4                        D  -19/
      DATA ARRH1A(21)/ 0.1                        D  -19/
      DATA AY1ASP/2.00135 24004 58893 96402 D   0,
     1            0.71104 24159 64619 38    D  -3,
     2            0.36659 77028 23244 9     D  -4,
     3            0.19130 15686 57728       D  -5,
     4            0.10046 91138 9777        D  -6,
     5            0.53040 17425 38          D  -8,
     6            0.28100 88617 6           D  -9,
     7            0.14938 86051             D -10,
     8            0.79578 420               D -12,
     9            0.42523 63                D -13,
     X            0.22719 5                 D -14,
     1            0.12216                   D -15,
     2            0.650                     D -17,
     3            0.36                      D -18,
     4            0.2                       D -19/
      DATA AY1ASQ/5.99065 10947 78881 89116 D   0,
     1           -0.48959 32623 36579 635   D  -2,
     2           -0.23238 32130 70706 26    D  -3,
     3           -0.11447 34723 85767 9     D  -4,
     4           -0.57169 92618 9106        D  -6,
     5           -0.28955 16716 917         D  -7,
     6           -0.14751 33456 36          D  -8,
     7           -0.75965 37378             D -10,
     8           -0.39065 8184              D -11,
     9           -0.20464 654               D -12,
     X           -0.10426 36                D -13,
     1           -0.57702                   D -15,
     2           -0.2550                    D -16,
     3           -0.210                     D -17,
     4            0.2                       D -19,
     5           -0.2                       D -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/15,17,12,13/
      DATA XLOW1,XLOW2/2.23750222D-154,4.08085106D-8/
      DATA XHIGH/4.50359963D15/
C
C   Start computation
C
      X = ABS ( XVALUE )
C
C   Error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DSTRVH1 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 9
C
      IF ( X .LE. NINE ) THEN
         IF ( X .LT. XLOW1 ) THEN
            DSTRVH1 = ZERO
         ELSE
            XSQ = X * X
            IF ( X .LT. XLOW2 ) THEN
               DSTRVH1 = TWOBPI * XSQ
            ELSE
               T = ( XSQ / FORTP5 - HALF ) - HALF
               DSTRVH1 = TWOBPI * XSQ * DCHEVAL ( NTERM1 , ARRH1 , T )
            ENDIF
         ENDIF
      ELSE      
C
C   Code for abs(x) > 9
C
         XSQ = X * X
         T = ( ONE82 - XSQ ) / ( TWENTY + XSQ )
         Y1P = DCHEVAL ( NTERM3 , AY1ASP , T )
         Y1Q = DCHEVAL ( NTERM4 , AY1ASQ , T ) / ( EIGHT * X)
         XM3P4 = X - THPBY4
         Y1VAL = Y1P * SIN ( XM3P4 ) + Y1Q * COS ( XM3P4 )
         Y1VAL = Y1VAL * RT2BPI / SQRT ( X )
         T = ( TW02P5 - XSQ ) / ( FORTP5 + XSQ )
         H1AS = TWOBPI * DCHEVAL ( NTERM2 , ARRH1A , T )
         DSTRVH1 = Y1VAL + H1AS
      ENDIF
      RETURN
      END

      SUBROUTINE ERRPRN(FNNAME,ERRMSG)
C
C   DESCRIPTION:
C      This subroutine prints out an error message if
C      an error has occurred in one of the MISCFUN 
C      functions. 
C
C
C   INPUT PARAMETERS:
C
C      FNNAME - CHARACTER - The name of the function with the error.
C
C      ERRMSG - CHARACTER - The message to be printed out.
C
C
C   MACHINE-DEPENDENT PARAMETER:
C
C      OUTSTR - INTEGER - The numerical value of the output
C                         stream to be used for printing the
C                         error message. The subroutine has the
C                         default value   OUTSTR = 6.
C
C
C   AUTHOR:
C         DR. ALLAN J. MACLEOD,
C         DEPT. OF MATHEMATICS AND STATISTICS,
C         UNIVERSITY OF PAISLEY,
C         HIGH ST.,
C         PAISLEY
C         SCOTLAND
C         PA1 2BE
C
C         (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION: 
C                   2 JUNE, 1995
C
      IMPLICIT NONE
      INTEGER OUTSTR
      CHARACTER(LEN=6) FNNAME
      CHARACTER(LEN=*) ERRMSG
      DATA OUTSTR/6/
      WRITE(OUTSTR,1000)FNNAME
      WRITE(OUTSTR,2000)ERRMSG
 1000 FORMAT(/5X,'ERROR IN MISCFUN FUNCTION  ',A6)
 2000 FORMAT(/5X,A50)
      STOP 'Catastrophic error in MISCFUN, bailing out'
      END
