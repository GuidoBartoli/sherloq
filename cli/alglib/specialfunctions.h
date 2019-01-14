/*************************************************************************
ALGLIB 3.10.0 (source code generated 2015-08-19)
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#ifndef _specialfunctions_pkg_h
#define _specialfunctions_pkg_h
#include "ap.h"
#include "alglibinternal.h"

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{


/*************************************************************************
Gamma function

Input parameters:
    X   -   argument

Domain:
    0 < X < 171.6
    -170 < X < 0, X is not an integer.

Relative error:
 arithmetic   domain     # trials      peak         rms
    IEEE    -170,-33      20000       2.3e-15     3.3e-16
    IEEE     -33,  33     20000       9.4e-16     2.2e-16
    IEEE      33, 171.6   20000       2.3e-15     3.2e-16

Cephes Math Library Release 2.8:  June, 2000
Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
double gammafunction(const double x);


/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
double lngamma(const double x, double &sgngam);

/*************************************************************************
Error function

The integral is

                          x
                           -
                2         | |          2
  erf(x)  =  --------     |    exp( - t  ) dt.
             sqrt(pi)   | |
                         -
                          0

For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
erf(x) = 1 - erfc(x).


ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,1         30000       3.7e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double errorfunction(const double x);


/*************************************************************************
Complementary error function

 1 - erf(x) =

                          inf.
                            -
                 2         | |          2
  erfc(x)  =  --------     |    exp( - t  ) dt
              sqrt(pi)   | |
                          -
                           x


For small x, erfc(x) = 1 - erf(x); otherwise rational
approximations are computed.


ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,26.6417   30000       5.7e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double errorfunctionc(const double x);


/*************************************************************************
Normal distribution function

Returns the area under the Gaussian probability density
function, integrated from minus infinity to x:

                           x
                            -
                  1        | |          2
   ndtr(x)  = ---------    |    exp( - t /2 ) dt
              sqrt(2pi)  | |
                          -
                         -inf.

            =  ( 1 + erf(z) ) / 2
            =  erfc(z) / 2

where z = x/sqrt(2). Computation is via the functions
erf and erfc.


ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -13,0        30000       3.4e-14     6.7e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double normaldistribution(const double x);


/*************************************************************************
Inverse of the error function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double inverf(const double e);


/*************************************************************************
Inverse of Normal distribution function

Returns the argument, x, for which the area under the
Gaussian probability density function (integrated from
minus infinity to x) is equal to y.


For small arguments 0 < y < exp(-2), the program computes
z = sqrt( -2.0 * log(y) );  then the approximation is
x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
There are two rational functions P/Q, one for 0 < y < exp(-32)
and the other for y up to exp(-2).  For larger arguments,
w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).

ACCURACY:

                     Relative error:
arithmetic   domain        # trials      peak         rms
   IEEE     0.125, 1        20000       7.2e-16     1.3e-16
   IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double invnormaldistribution(const double y0);

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegamma(const double a, const double x);


/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegammac(const double a, const double x);


/*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invincompletegammac(const double a, const double y0);

/*************************************************************************
Airy function

Solution of the differential equation

y"(x) = xy.

The function returns the two independent solutions Ai, Bi
and their first derivatives Ai'(x), Bi'(x).

Evaluation is by power series summation for small x,
by rational minimax approximations for large x.



ACCURACY:
Error criterion is absolute when function <= 1, relative
when function > 1, except * denotes relative error criterion.
For large negative x, the absolute error increases as x^1.5.
For large positive x, the relative error increases as x^1.5.

Arithmetic  domain   function  # trials      peak         rms
IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
void airy(const double x, double &ai, double &aip, double &bi, double &bip);

/*************************************************************************
Bessel function of order zero

Returns Bessel function of order zero of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval the following rational
approximation is used:


       2         2
(w - r  ) (w - r  ) P (w) / Q (w)
      1         2    3       8

           2
where w = x  and the two r's are zeros of the function.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.

ACCURACY:

                     Absolute error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       60000       4.2e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double besselj0(const double x);


/*************************************************************************
Bessel function of order one

Returns Bessel function of order one of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 24 term Chebyshev
expansion is used. In the second, the asymptotic
trigonometric representation is employed using two
rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       2.6e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double besselj1(const double x);


/*************************************************************************
Bessel function of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The ratio of jn(x) to j0(x) is computed by backward
recurrence.  First the ratio jn/jn-1 is found by a
continued fraction expansion.  Then the recurrence
relating successive orders is applied until j0 or j1 is
reached.

If n = 0 or 1 the routine for j0 or j1 is called
directly.

ACCURACY:

                     Absolute error:
arithmetic   range      # trials      peak         rms
   IEEE      0, 30        5000       4.4e-16     7.9e-17


Not suitable for large n or x. Use jv() (fractional order) instead.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseljn(const ae_int_t n, const double x);


/*************************************************************************
Bessel function of the second kind, order zero

Returns Bessel function of the second kind, of order
zero, of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval a rational approximation
R(x) is employed to compute
  y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
Thus a call to j0() is required.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.



ACCURACY:

 Absolute error, when y0(x) < 1; else relative error:

arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.3e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double bessely0(const double x);


/*************************************************************************
Bessel function of second kind of order one

Returns Bessel function of the second kind of order one
of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 25 term Chebyshev
expansion is used, and a call to j1() is required.
In the second, the asymptotic trigonometric representation
is employed using two rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       1.0e-15     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double bessely1(const double x);


/*************************************************************************
Bessel function of second kind of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The function is evaluated by forward recurrence on
n, starting with values computed by the routines
y0() and y1().

If n = 0 or 1 the routine for y0 or y1 is called
directly.

ACCURACY:
                     Absolute error, except relative
                     when y > 1:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       3.4e-15     4.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselyn(const ae_int_t n, const double x);


/*************************************************************************
Modified Bessel function of order zero

Returns modified Bessel function of order zero of the
argument.

The function is defined as i0(x) = j0( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       5.8e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseli0(const double x);


/*************************************************************************
Modified Bessel function of order one

Returns modified Bessel function of order one of the
argument.

The function is defined as i1(x) = -i j1( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.9e-15     2.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseli1(const double x);


/*************************************************************************
Modified Bessel function, second kind, order zero

Returns modified Bessel function of the second kind
of order zero of the argument.

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

Tested at 2000 random points between 0 and 8.  Peak absolute
error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselk0(const double x);


/*************************************************************************
Modified Bessel function, second kind, order one

Computes the modified Bessel function of the second kind
of order one of the argument.

The range is partitioned into the two intervals [0,2] and
(2, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselk1(const double x);


/*************************************************************************
Modified Bessel function, second kind, integer order

Returns modified Bessel function of the second kind
of order n of the argument.

The range is partitioned into the two intervals [0,9.55] and
(9.55, infinity).  An ascending power series is used in the
low range, and an asymptotic expansion in the high range.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        90000       1.8e-8      3.0e-10

Error is high only near the crossover point x = 9.55
between the two expansions used.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*************************************************************************/
double besselkn(const ae_int_t nn, const double x);

/*************************************************************************
Beta function


                  -     -
                 | (a) | (b)
beta( a, b )  =  -----------.
                    -
                   | (a+b)

For large arguments the logarithm of the function is
evaluated using lgam(), then exponentiated.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,30       30000       8.1e-14     1.1e-14

Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
*************************************************************************/
double beta(const double a, const double b);

/*************************************************************************
Incomplete beta integral

Returns incomplete beta integral of the arguments, evaluated
from zero to x.  The function is defined as

                 x
    -            -
   | (a+b)      | |  a-1     b-1
 -----------    |   t   (1-t)   dt.
  -     -     | |
 | (a) | (b)   -
                0

The domain of definition is 0 <= x <= 1.  In this
implementation a and b are restricted to positive values.
The integral from x to 1 may be obtained by the symmetry
relation

   1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).

The integral is evaluated by a continued fraction expansion
or, when b*x is small, by a power series.

ACCURACY:

Tested at uniformly distributed random points (a,b,x) with a and b
in "domain" and x between 0 and 1.
                                       Relative error
arithmetic   domain     # trials      peak         rms
   IEEE      0,5         10000       6.9e-15     4.5e-16
   IEEE      0,85       250000       2.2e-13     1.7e-14
   IEEE      0,1000      30000       5.3e-12     6.3e-13
   IEEE      0,10000    250000       9.3e-11     7.1e-12
   IEEE      0,100000    10000       8.7e-10     4.8e-11
Outputs smaller than the IEEE gradual underflow threshold
were excluded from these statistics.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletebeta(const double a, const double b, const double x);


/*************************************************************************
Inverse of imcomplete beta integral

Given y, the function finds x such that

 incbet( a, b, x ) = y .

The routine performs interval halving or Newton iterations to find the
root of incbet(a,b,x) - y = 0.


ACCURACY:

                     Relative error:
               x     a,b
arithmetic   domain  domain  # trials    peak       rms
   IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
   IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
   IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
With a and b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
   IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
With a = .5, b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1996, 2000 by Stephen L. Moshier
*************************************************************************/
double invincompletebeta(const double a, const double b, const double y);

/*************************************************************************
Binomial distribution

Returns the sum of the terms 0 through k of the Binomial
probability density:

  k
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=0

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p), with p between 0 and 1.

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      4.3e-15     2.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p);


/*************************************************************************
Complemented binomial distribution

Returns the sum of the terms k+1 through n of the Binomial
probability density:

  n
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=k+1

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      6.7e-15     8.2e-16
 For p between 0 and .001:
   IEEE     0,100       100000      1.5e-13     2.7e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p);


/*************************************************************************
Inverse binomial distribution

Finds the event probability p such that the sum of the
terms 0 through k of the Binomial probability density
is equal to the given cumulative probability y.

This is accomplished using the inverse beta integral
function and the relation

1 - p = incbi( n-k, k+1, y ).

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      2.3e-14     6.4e-16
   IEEE     0,10000     100000      6.6e-12     1.2e-13
 For p between 10^-6 and 0.001:
   IEEE     0,100       100000      2.0e-12     1.3e-14
   IEEE     0,10000     100000      1.5e-12     3.2e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y);

/*************************************************************************
Calculation of the value of the Chebyshev polynomials of the
first and second kinds.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument, -1 <= x <= 1

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************/
double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x);


/*************************************************************************
Summation of Chebyshev polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*T0(x) + c[1]*T1(x) + ... + c[N]*TN(x)
or
    c[0]*U0(x) + c[1]*U1(x) + ... + c[N]*UN(x)
depending on the R.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************/
double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x);


/*************************************************************************
Representation of Tn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void chebyshevcoefficients(const ae_int_t n, real_1d_array &c);


/*************************************************************************
Conversion of a series of Chebyshev polynomials to a power series.

Represents A[0]*T0(x) + A[1]*T1(x) + ... + A[N]*Tn(x) as
B[0] + B[1]*X + ... + B[N]*X^N.

Input parameters:
    A   -   Chebyshev series coefficients
    N   -   degree, N>=0

Output parameters
    B   -   power series coefficients
*************************************************************************/
void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b);

/*************************************************************************
Chi-square distribution

Returns the area under the left hand tail (from 0 to x)
of the Chi square probability density function with
v degrees of freedom.


                                  x
                                   -
                       1          | |  v/2-1  -t/2
 P( x | v )   =   -----------     |   t      e     dt
                   v/2  -       | |
                  2    | (v/2)   -
                                  0

where x is the Chi-square variable.

The incomplete gamma integral is used, according to the
formula

y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function


Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double chisquaredistribution(const double v, const double x);


/*************************************************************************
Complemented Chi-square distribution

Returns the area under the right hand tail (from x to
infinity) of the Chi square probability density function
with v degrees of freedom:

                                 inf.
                                   -
                       1          | |  v/2-1  -t/2
 P( x | v )   =   -----------     |   t      e     dt
                   v/2  -       | |
                  2    | (v/2)   -
                                  x

where x is the Chi-square variable.

The incomplete gamma integral is used, according to the
formula

y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double chisquarecdistribution(const double v, const double x);


/*************************************************************************
Inverse of complemented Chi-square distribution

Finds the Chi-square argument x such that the integral
from x to infinity of the Chi-square density is equal
to the given cumulative probability y.

This is accomplished using the inverse gamma integral
function and the relation

   x/2 = igami( df/2, y );

ACCURACY:

See inverse incomplete gamma function


Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double invchisquaredistribution(const double v, const double y);

/*************************************************************************
Dawson's Integral

Approximates the integral

                            x
                            -
                     2     | |        2
 dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
                         | |
                          -
                          0

Three different rational approximations are employed, for
the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,10        10000       6.9e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double dawsonintegral(const double x);

/*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

using the approximation

    P(x)  -  log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double ellipticintegralk(const double m);


/*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

where m = 1 - m1, using the approximation

    P(x)  -  log x Q(x).

The argument m1 is used rather than m so that the logarithmic
singularity at m = 1 will be shifted to the origin; this
preserves maximum accuracy.

K(0) = pi/2.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double ellipticintegralkhighprecision(const double m1);


/*************************************************************************
Incomplete elliptic integral of the first kind F(phi|m)

Approximates the integral



               phi
                -
               | |
               |           dt
F(phi_\m)  =    |    ------------------
               |                   2
             | |    sqrt( 1 - m sin t )
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.




ACCURACY:

Tested at random points with m in [0, 1] and phi as indicated.

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10       200000      7.4e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompleteellipticintegralk(const double phi, const double m);


/*************************************************************************
Complete elliptic integral of the second kind

Approximates the integral


           pi/2
            -
           | |                 2
E(m)  =    |    sqrt( 1 - m sin t ) dt
         | |
          -
           0

using the approximation

     P(x)  -  x log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0, 1       10000       2.1e-16     7.3e-17

Cephes Math Library, Release 2.8: June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double ellipticintegrale(const double m);


/*************************************************************************
Incomplete elliptic integral of the second kind

Approximates the integral


               phi
                -
               | |
               |                   2
E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
               |
             | |
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.

ACCURACY:

Tested at random arguments with phi in [-10, 10] and m in
[0, 1].
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10      150000       3.3e-15     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1993, 2000 by Stephen L. Moshier
*************************************************************************/
double incompleteellipticintegrale(const double phi, const double m);

/*************************************************************************
Exponential integral Ei(x)

              x
               -     t
              | |   e
   Ei(x) =   -|-   ---  dt .
            | |     t
             -
            -inf

Not defined for x <= 0.
See also expn.c.



ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,100       50000      8.6e-16     1.3e-16

Cephes Math Library Release 2.8:  May, 1999
Copyright 1999 by Stephen L. Moshier
*************************************************************************/
double exponentialintegralei(const double x);


/*************************************************************************
Exponential integral En(x)

Evaluates the exponential integral

                inf.
                  -
                 | |   -xt
                 |    e
     E (x)  =    |    ----  dt.
      n          |      n
               | |     t
                -
                 1


Both n and x must be nonnegative.

The routine employs either a power series, a continued
fraction, or an asymptotic formula depending on the
relative values of n and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       10000       1.7e-15     3.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 2000 by Stephen L. Moshier
*************************************************************************/
double exponentialintegralen(const double x, const ae_int_t n);

/*************************************************************************
F distribution

Returns the area from zero to x under the F density
function (also known as Snedcor's density or the
variance ratio density).  This is the density
of x = (u1/df1)/(u2/df2), where u1 and u2 are random
variables having Chi square distributions with df1
and df2 degrees of freedom, respectively.
The incomplete beta integral is used, according to the
formula

P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).


The arguments a and b are greater than zero, and x is
nonnegative.

ACCURACY:

Tested at random points (a,b,x).

               x     a,b                     Relative error:
arithmetic  domain  domain     # trials      peak         rms
   IEEE      0,1    0,100       100000      9.8e-15     1.7e-15
   IEEE      1,5    0,100       100000      6.5e-15     3.5e-16
   IEEE      0,1    1,10000     100000      2.2e-11     3.3e-12
   IEEE      1,5    1,10000     100000      1.1e-11     1.7e-13

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double fdistribution(const ae_int_t a, const ae_int_t b, const double x);


/*************************************************************************
Complemented F distribution

Returns the area from x to infinity under the F density
function (also known as Snedcor's density or the
variance ratio density).


                     inf.
                      -
             1       | |  a-1      b-1
1-P(x)  =  ------    |   t    (1-t)    dt
           B(a,b)  | |
                    -
                     x


The incomplete beta integral is used, according to the
formula

P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).


ACCURACY:

Tested at random points (a,b,x) in the indicated intervals.
               x     a,b                     Relative error:
arithmetic  domain  domain     # trials      peak         rms
   IEEE      0,1    1,100       100000      3.7e-14     5.9e-16
   IEEE      1,5    1,100       100000      8.0e-15     1.6e-15
   IEEE      0,1    1,10000     100000      1.8e-11     3.5e-13
   IEEE      1,5    1,10000     100000      2.0e-11     3.0e-12

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double fcdistribution(const ae_int_t a, const ae_int_t b, const double x);


/*************************************************************************
Inverse of complemented F distribution

Finds the F density argument x such that the integral
from x to infinity of the F density is equal to the
given probability p.

This is accomplished using the inverse beta integral
function and the relations

     z = incbi( df2/2, df1/2, p )
     x = df2 (1-z) / (df1 z).

Note: the following relations hold for the inverse of
the uncomplemented F distribution:

     z = incbi( df1/2, df2/2, p )
     x = df2 z / (df1 (1-z)).

ACCURACY:

Tested at random points (a,b,p).

             a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between .001 and 1:
   IEEE     1,100       100000      8.3e-15     4.7e-16
   IEEE     1,10000     100000      2.1e-11     1.4e-13
 For p between 10^-6 and 10^-3:
   IEEE     1,100        50000      1.3e-12     8.4e-15
   IEEE     1,10000      50000      3.0e-12     4.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invfdistribution(const ae_int_t a, const ae_int_t b, const double y);

/*************************************************************************
Fresnel integral

Evaluates the Fresnel integrals

          x
          -
         | |
C(x) =   |   cos(pi/2 t**2) dt,
       | |
        -
         0

          x
          -
         | |
S(x) =   |   sin(pi/2 t**2) dt.
       | |
        -
         0


The integrals are evaluated by a power series for x < 1.
For x >= 1 auxiliary functions f(x) and g(x) are employed
such that

C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )



ACCURACY:

 Relative error.

Arithmetic  function   domain     # trials      peak         rms
  IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
  IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
void fresnelintegral(const double x, double &c, double &s);

/*************************************************************************
Calculation of the value of the Hermite polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial Hn at x
*************************************************************************/
double hermitecalculate(const ae_int_t n, const double x);


/*************************************************************************
Summation of Hermite polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*H0(x) + c[1]*H1(x) + ... + c[N]*HN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial at x
*************************************************************************/
double hermitesum(const real_1d_array &c, const ae_int_t n, const double x);


/*************************************************************************
Representation of Hn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void hermitecoefficients(const ae_int_t n, real_1d_array &c);

/*************************************************************************
Jacobian Elliptic Functions

Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
and dn(u|m) of parameter m between 0 and 1, and real
argument u.

These functions are periodic, with quarter-period on the
real axis equal to the complete elliptic integral
ellpk(1.0-m).

Relation to incomplete elliptic integral:
If u = ellik(phi,m), then sn(u|m) = sin(phi),
and cn(u|m) = cos(phi).  Phi is called the amplitude of u.

Computation is by means of the arithmetic-geometric mean
algorithm, except when m is within 1e-9 of 0 or 1.  In the
latter case with m close to 1, the approximation applies
only for phi < pi/2.

ACCURACY:

Tested at random points with u between 0 and 10, m between
0 and 1.

           Absolute error (* = relative error):
arithmetic   function   # trials      peak         rms
   IEEE      phi         10000       9.2e-16*    1.4e-16*
   IEEE      sn          50000       4.1e-15     4.6e-16
   IEEE      cn          40000       3.6e-15     4.4e-16
   IEEE      dn          10000       1.3e-12     1.8e-14

 Peak error observed in consistency check using addition
theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
the above relation to the incomplete elliptic integral.
Accuracy deteriorates when u is large.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph);

/*************************************************************************
Calculation of the value of the Laguerre polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial Ln at x
*************************************************************************/
double laguerrecalculate(const ae_int_t n, const double x);


/*************************************************************************
Summation of Laguerre polynomials using Clenshaw’s recurrence formula.

This routine calculates c[0]*L0(x) + c[1]*L1(x) + ... + c[N]*LN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial at x
*************************************************************************/
double laguerresum(const real_1d_array &c, const ae_int_t n, const double x);


/*************************************************************************
Representation of Ln as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void laguerrecoefficients(const ae_int_t n, real_1d_array &c);

/*************************************************************************
Calculation of the value of the Legendre polynomial Pn.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial Pn at x
*************************************************************************/
double legendrecalculate(const ae_int_t n, const double x);


/*************************************************************************
Summation of Legendre polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*P0(x) + c[1]*P1(x) + ... + c[N]*PN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial at x
*************************************************************************/
double legendresum(const real_1d_array &c, const ae_int_t n, const double x);


/*************************************************************************
Representation of Pn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void legendrecoefficients(const ae_int_t n, real_1d_array &c);

/*************************************************************************
Poisson distribution

Returns the sum of the first k+1 terms of the Poisson
distribution:

  k         j
  --   -m  m
  >   e    --
  --       j!
 j=0

The terms are not summed directly; instead the incomplete
gamma integral is employed, according to the relation

y = pdtr( k, m ) = igamc( k+1, m ).

The arguments must both be positive.
ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double poissondistribution(const ae_int_t k, const double m);


/*************************************************************************
Complemented Poisson distribution

Returns the sum of the terms k+1 to infinity of the Poisson
distribution:

 inf.       j
  --   -m  m
  >   e    --
  --       j!
 j=k+1

The terms are not summed directly; instead the incomplete
gamma integral is employed, according to the formula

y = pdtrc( k, m ) = igam( k+1, m ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double poissoncdistribution(const ae_int_t k, const double m);


/*************************************************************************
Inverse Poisson distribution

Finds the Poisson variable x such that the integral
from 0 to x of the Poisson density is equal to the
given probability y.

This is accomplished using the inverse gamma integral
function and the relation

   m = igami( k+1, y ).

ACCURACY:

See inverse incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invpoissondistribution(const ae_int_t k, const double y);

/*************************************************************************
Psi (digamma) function

             d      -
  psi(x)  =  -- ln | (x)
             dx

is the logarithmic derivative of the gamma function.
For integer x,
                  n-1
                   -
psi(n) = -EUL  +   >  1/k.
                   -
                  k=1

This formula is used for 0 < n <= 10.  If x is negative, it
is transformed to a positive argument by the reflection
formula  psi(1-x) = psi(x) + pi cot(pi x).
For general positive x, the argument is made greater than 10
using the recurrence  psi(x+1) = psi(x) + 1/x.
Then the following asymptotic expansion is applied:

                          inf.   B
                           -      2k
psi(x) = log(x) - 1/2x -   >   -------
                           -        2k
                          k=1   2k x

where the B2k are Bernoulli numbers.

ACCURACY:
   Relative error (except absolute when |psi| < 1):
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       1.3e-15     1.4e-16
   IEEE      -30,0       40000       1.5e-15     2.2e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
double psi(const double x);

/*************************************************************************
Student's t distribution

Computes the integral from minus infinity to t of the Student
t distribution with integer k > 0 degrees of freedom:

                                     t
                                     -
                                    | |
             -                      |         2   -(k+1)/2
            | ( (k+1)/2 )           |  (     x   )
      ----------------------        |  ( 1 + --- )        dx
                    -               |  (      k  )
      sqrt( k pi ) | ( k/2 )        |
                                  | |
                                   -
                                  -inf.

Relation to incomplete beta integral:

       1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
where
       z = k/(k + t**2).

For t < -2, this is the method of computation.  For higher t,
a direct method is derived from integration by parts.
Since the function is symmetric about t=0, the area under the
right tail of the density is found by calling the function
with -t instead of t.

ACCURACY:

Tested at random 1 <= k <= 25.  The "domain" refers to t.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -100,-2      50000       5.9e-15     1.4e-15
   IEEE     -2,100      500000       2.7e-15     4.9e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double studenttdistribution(const ae_int_t k, const double t);


/*************************************************************************
Functional inverse of Student's t distribution

Given probability p, finds the argument t such that stdtr(k,t)
is equal to p.

ACCURACY:

Tested at random 1 <= k <= 100.  The "domain" refers to p:
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE    .001,.999     25000       5.7e-15     8.0e-16
   IEEE    10^-6,.001    25000       2.0e-12     2.9e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invstudenttdistribution(const ae_int_t k, const double p);

/*************************************************************************
Sine and cosine integrals

Evaluates the integrals

                         x
                         -
                        |  cos t - 1
  Ci(x) = eul + ln x +  |  --------- dt,
                        |      t
                       -
                        0
            x
            -
           |  sin t
  Si(x) =  |  ----- dt
           |    t
          -
           0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are approximated by rational functions.
For x > 8 auxiliary functions f(x) and g(x) are employed
such that

Ci(x) = f(x) sin(x) - g(x) cos(x)
Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)


ACCURACY:
   Test interval = [0,50].
Absolute error, except relative when > 1:
arithmetic   function   # trials      peak         rms
   IEEE        Si        30000       4.4e-16     7.3e-17
   IEEE        Ci        30000       6.9e-16     5.1e-17

Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*************************************************************************/
void sinecosineintegrals(const double x, double &si, double &ci);


/*************************************************************************
Hyperbolic sine and cosine integrals

Approximates the integrals

                           x
                           -
                          | |   cosh t - 1
  Chi(x) = eul + ln x +   |    -----------  dt,
                        | |          t
                         -
                         0

              x
              -
             | |  sinh t
  Shi(x) =   |    ------  dt
           | |       t
            -
            0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are evaluated by power series for x < 8
and by Chebyshev expansions for x between 8 and 88.
For large x, both functions approach exp(x)/2x.
Arguments greater than 88 in magnitude return MAXNUM.


ACCURACY:

Test interval 0 to 88.
                     Relative error:
arithmetic   function  # trials      peak         rms
   IEEE         Shi      30000       6.9e-16     1.6e-16
       Absolute error, except relative when |Chi| > 1:
   IEEE         Chi      30000       8.4e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi);
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
double gammafunction(double x, ae_state *_state);
double lngamma(double x, double* sgngam, ae_state *_state);
double errorfunction(double x, ae_state *_state);
double errorfunctionc(double x, ae_state *_state);
double normaldistribution(double x, ae_state *_state);
double inverf(double e, ae_state *_state);
double invnormaldistribution(double y0, ae_state *_state);
double incompletegamma(double a, double x, ae_state *_state);
double incompletegammac(double a, double x, ae_state *_state);
double invincompletegammac(double a, double y0, ae_state *_state);
void airy(double x,
     double* ai,
     double* aip,
     double* bi,
     double* bip,
     ae_state *_state);
double besselj0(double x, ae_state *_state);
double besselj1(double x, ae_state *_state);
double besseljn(ae_int_t n, double x, ae_state *_state);
double bessely0(double x, ae_state *_state);
double bessely1(double x, ae_state *_state);
double besselyn(ae_int_t n, double x, ae_state *_state);
double besseli0(double x, ae_state *_state);
double besseli1(double x, ae_state *_state);
double besselk0(double x, ae_state *_state);
double besselk1(double x, ae_state *_state);
double besselkn(ae_int_t nn, double x, ae_state *_state);
double beta(double a, double b, ae_state *_state);
double incompletebeta(double a, double b, double x, ae_state *_state);
double invincompletebeta(double a, double b, double y, ae_state *_state);
double binomialdistribution(ae_int_t k,
     ae_int_t n,
     double p,
     ae_state *_state);
double binomialcdistribution(ae_int_t k,
     ae_int_t n,
     double p,
     ae_state *_state);
double invbinomialdistribution(ae_int_t k,
     ae_int_t n,
     double y,
     ae_state *_state);
double chebyshevcalculate(ae_int_t r,
     ae_int_t n,
     double x,
     ae_state *_state);
double chebyshevsum(/* Real    */ ae_vector* c,
     ae_int_t r,
     ae_int_t n,
     double x,
     ae_state *_state);
void chebyshevcoefficients(ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state);
void fromchebyshev(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_state *_state);
double chisquaredistribution(double v, double x, ae_state *_state);
double chisquarecdistribution(double v, double x, ae_state *_state);
double invchisquaredistribution(double v, double y, ae_state *_state);
double dawsonintegral(double x, ae_state *_state);
double ellipticintegralk(double m, ae_state *_state);
double ellipticintegralkhighprecision(double m1, ae_state *_state);
double incompleteellipticintegralk(double phi, double m, ae_state *_state);
double ellipticintegrale(double m, ae_state *_state);
double incompleteellipticintegrale(double phi, double m, ae_state *_state);
double exponentialintegralei(double x, ae_state *_state);
double exponentialintegralen(double x, ae_int_t n, ae_state *_state);
double fdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double fcdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
double invfdistribution(ae_int_t a,
     ae_int_t b,
     double y,
     ae_state *_state);
void fresnelintegral(double x, double* c, double* s, ae_state *_state);
double hermitecalculate(ae_int_t n, double x, ae_state *_state);
double hermitesum(/* Real    */ ae_vector* c,
     ae_int_t n,
     double x,
     ae_state *_state);
void hermitecoefficients(ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state);
void jacobianellipticfunctions(double u,
     double m,
     double* sn,
     double* cn,
     double* dn,
     double* ph,
     ae_state *_state);
double laguerrecalculate(ae_int_t n, double x, ae_state *_state);
double laguerresum(/* Real    */ ae_vector* c,
     ae_int_t n,
     double x,
     ae_state *_state);
void laguerrecoefficients(ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state);
double legendrecalculate(ae_int_t n, double x, ae_state *_state);
double legendresum(/* Real    */ ae_vector* c,
     ae_int_t n,
     double x,
     ae_state *_state);
void legendrecoefficients(ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state);
double poissondistribution(ae_int_t k, double m, ae_state *_state);
double poissoncdistribution(ae_int_t k, double m, ae_state *_state);
double invpoissondistribution(ae_int_t k, double y, ae_state *_state);
double psi(double x, ae_state *_state);
double studenttdistribution(ae_int_t k, double t, ae_state *_state);
double invstudenttdistribution(ae_int_t k, double p, ae_state *_state);
void sinecosineintegrals(double x,
     double* si,
     double* ci,
     ae_state *_state);
void hyperbolicsinecosineintegrals(double x,
     double* shi,
     double* chi,
     ae_state *_state);

}
#endif

