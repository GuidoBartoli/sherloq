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
#include "stdafx.h"
#include "fasttransforms.h"

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{


/*************************************************************************
1-dimensional complex FFT.

Array size N may be arbitrary number (composite or prime).  Composite  N's
are handled with cache-oblivious variation of  a  Cooley-Tukey  algorithm.
Small prime-factors are transformed using hard coded  codelets (similar to
FFTW codelets, but without low-level  optimization),  large  prime-factors
are handled with Bluestein's algorithm.

Fastests transforms are for smooth N's (prime factors are 2, 3,  5  only),
most fast for powers of 2. When N have prime factors  larger  than  these,
but orders of magnitude smaller than N, computations will be about 4 times
slower than for nearby highly composite N's. When N itself is prime, speed
will be 6 times lower.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1d(complex_1d_array &a, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftc1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex FFT.

Array size N may be arbitrary number (composite or prime).  Composite  N's
are handled with cache-oblivious variation of  a  Cooley-Tukey  algorithm.
Small prime-factors are transformed using hard coded  codelets (similar to
FFTW codelets, but without low-level  optimization),  large  prime-factors
are handled with Bluestein's algorithm.

Fastests transforms are for smooth N's (prime factors are 2, 3,  5  only),
most fast for powers of 2. When N have prime factors  larger  than  these,
but orders of magnitude smaller than N, computations will be about 4 times
slower than for nearby highly composite N's. When N itself is prime, speed
will be 6 times lower.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1d(complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = a.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftc1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex inverse FFT.

Array size N may be arbitrary number (composite or prime).  Algorithm  has
O(N*logN) complexity for any N (composite or prime).

See FFTC1D() description for more information about algorithm performance.

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1dinv(complex_1d_array &a, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftc1dinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex inverse FFT.

Array size N may be arbitrary number (composite or prime).  Algorithm  has
O(N*logN) complexity for any N (composite or prime).

See FFTC1D() description for more information about algorithm performance.

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1dinv(complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = a.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftc1dinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    F   -   DFT of a input array, array[0..N-1]
            F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)

NOTE:
    F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
of  array  is  usually needed. But for convinience subroutine returns full
complex array (with frequencies above N/2), so its result may be  used  by
other FFT-related subroutines.


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1d(const real_1d_array &a, const ae_int_t n, complex_1d_array &f)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftr1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(f.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    F   -   DFT of a input array, array[0..N-1]
            F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)

NOTE:
    F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
of  array  is  usually needed. But for convinience subroutine returns full
complex array (with frequencies above N/2), so its result may be  used  by
other FFT-related subroutines.


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1d(const real_1d_array &a, complex_1d_array &f)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = a.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftr1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(f.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real inverse FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    F   -   array[0..floor(N/2)] - frequencies from forward real FFT
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]

NOTE:
    F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
F[floor(N/2)] has no special properties.

Relying on properties noted above, FFTR1DInv subroutine uses only elements
from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
N is even it ignores imaginary part of F[floor(N/2)] too.

When you call this function using full arguments list - "FFTR1DInv(F,N,A)"
- you can pass either either frequencies array with N elements or  reduced
array with roughly N/2 elements - subroutine will  successfully  transform
both.

If you call this function using reduced arguments list -  "FFTR1DInv(F,A)"
- you must pass FULL array with N elements (although higher  N/2 are still
not used) because array size is used to automatically determine FFT length


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinv(const complex_1d_array &f, const ae_int_t n, real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftr1dinv(const_cast<alglib_impl::ae_vector*>(f.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real inverse FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    F   -   array[0..floor(N/2)] - frequencies from forward real FFT
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]

NOTE:
    F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
F[floor(N/2)] has no special properties.

Relying on properties noted above, FFTR1DInv subroutine uses only elements
from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
N is even it ignores imaginary part of F[floor(N/2)] too.

When you call this function using full arguments list - "FFTR1DInv(F,N,A)"
- you can pass either either frequencies array with N elements or  reduced
array with roughly N/2 elements - subroutine will  successfully  transform
both.

If you call this function using reduced arguments list -  "FFTR1DInv(F,A)"
- you must pass FULL array with N elements (although higher  N/2 are still
not used) because array size is used to automatically determine FFT length


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinv(const complex_1d_array &f, real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = f.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fftr1dinv(const_cast<alglib_impl::ae_vector*>(f.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex convolution.

For given A/B returns conv(A,B) (non-circular). Subroutine can automatically
choose between three implementations: straightforward O(M*N)  formula  for
very small N (or M), overlap-add algorithm for  cases  where  max(M,N)  is
significantly larger than min(M,N), but O(M*N) algorithm is too slow,  and
general FFT-based formula for cases where two previois algorithms are  too
slow.

Algorithm has max(M,N)*log(max(M,N)) complexity for any M/N.

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1d(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convc1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex non-circular deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convc1dinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional circular complex convolution.

For given S/R returns conv(S,R) (circular). Algorithm has linearithmic
complexity for any M/N.

IMPORTANT:  normal convolution is commutative,  i.e.   it  is symmetric  -
conv(A,B)=conv(B,A).  Cyclic convolution IS NOT.  One function - S - is  a
signal,  periodic function, and another - R - is a response,  non-periodic
function with limited length.

INPUT PARAMETERS
    S   -   array[0..M-1] - complex periodic signal
    M   -   problem size
    B   -   array[0..N-1] - complex non-periodic response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircular(const complex_1d_array &s, const ae_int_t m, const complex_1d_array &r, const ae_int_t n, complex_1d_array &c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convc1dcircular(const_cast<alglib_impl::ae_vector*>(s.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(c.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional circular complex deconvolution (inverse of ConvC1DCircular()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved periodic signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - non-periodic response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-1].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircularinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convc1dcircularinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real convolution.

Analogous to ConvC1D(), see ConvC1D() comments for more details.

INPUT PARAMETERS
    A   -   array[0..M-1] - real function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1d(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convr1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convr1dinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional circular real convolution.

Analogous to ConvC1DCircular(), see ConvC1DCircular() comments for more details.

INPUT PARAMETERS
    S   -   array[0..M-1] - real signal
    M   -   problem size
    B   -   array[0..N-1] - real response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircular(const real_1d_array &s, const ae_int_t m, const real_1d_array &r, const ae_int_t n, real_1d_array &c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convr1dcircular(const_cast<alglib_impl::ae_vector*>(s.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(c.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircularinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::convr1dcircularinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(b.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrC1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(conj(pattern[j])*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1d(const complex_1d_array &signal, const ae_int_t n, const complex_1d_array &pattern, const ae_int_t m, complex_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::corrc1d(const_cast<alglib_impl::ae_vector*>(signal.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(pattern.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional circular complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrC1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1dcircular(const complex_1d_array &signal, const ae_int_t m, const complex_1d_array &pattern, const ae_int_t n, complex_1d_array &c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::corrc1dcircular(const_cast<alglib_impl::ae_vector*>(signal.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(pattern.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(c.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrR1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(pattern[j]*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(pattern[j]*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1d(const real_1d_array &signal, const ae_int_t n, const real_1d_array &pattern, const ae_int_t m, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::corrr1d(const_cast<alglib_impl::ae_vector*>(signal.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(pattern.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional circular real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrR1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1dcircular(const real_1d_array &signal, const ae_int_t m, const real_1d_array &pattern, const ae_int_t n, real_1d_array &c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::corrr1dcircular(const_cast<alglib_impl::ae_vector*>(signal.c_ptr()), m, const_cast<alglib_impl::ae_vector*>(pattern.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(c.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional Fast Hartley Transform.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   FHT of a input array, array[0..N-1],
            A_out[k] = sum(A_in[j]*(cos(2*pi*j*k/N)+sin(2*pi*j*k/N)), j=0..N-1)


  -- ALGLIB --
     Copyright 04.06.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1d(real_1d_array &a, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fhtr1d(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
1-dimensional inverse FHT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse FHT of a input array, array[0..N-1]


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1dinv(real_1d_array &a, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fhtr1dinv(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF COMPUTATIONAL CORE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{











/*************************************************************************
1-dimensional complex FFT.

Array size N may be arbitrary number (composite or prime).  Composite  N's
are handled with cache-oblivious variation of  a  Cooley-Tukey  algorithm.
Small prime-factors are transformed using hard coded  codelets (similar to
FFTW codelets, but without low-level  optimization),  large  prime-factors
are handled with Bluestein's algorithm.

Fastests transforms are for smooth N's (prime factors are 2, 3,  5  only),
most fast for powers of 2. When N have prime factors  larger  than  these,
but orders of magnitude smaller than N, computations will be about 4 times
slower than for nearby highly composite N's. When N itself is prime, speed
will be 6 times lower.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex function to be transformed
    N   -   problem size
    
OUTPUT PARAMETERS
    A   -   DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1d(/* Complex */ ae_vector* a, ae_int_t n, ae_state *_state)
{
    ae_frame _frame_block;
    fasttransformplan plan;
    ae_int_t i;
    ae_vector buf;

    ae_frame_make(_state, &_frame_block);
    _fasttransformplan_init(&plan, _state);
    ae_vector_init(&buf, 0, DT_REAL, _state);

    ae_assert(n>0, "FFTC1D: incorrect N!", _state);
    ae_assert(a->cnt>=n, "FFTC1D: Length(A)<N!", _state);
    ae_assert(isfinitecvector(a, n, _state), "FFTC1D: A contains infinite or NAN values!", _state);
    
    /*
     * Special case: N=1, FFT is just identity transform.
     * After this block we assume that N is strictly greater than 1.
     */
    if( n==1 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * convert input array to the more convinient format
     */
    ae_vector_set_length(&buf, 2*n, _state);
    for(i=0; i<=n-1; i++)
    {
        buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
        buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
    }
    
    /*
     * Generate plan and execute it.
     *
     * Plan is a combination of a successive factorizations of N and
     * precomputed data. It is much like a FFTW plan, but is not stored
     * between subroutine calls and is much simpler.
     */
    ftcomplexfftplan(n, 1, &plan, _state);
    ftapplyplan(&plan, &buf, 0, 1, _state);
    
    /*
     * result
     */
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_complex[i].x = buf.ptr.p_double[2*i+0];
        a->ptr.p_complex[i].y = buf.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional complex inverse FFT.

Array size N may be arbitrary number (composite or prime).  Algorithm  has
O(N*logN) complexity for any N (composite or prime).

See FFTC1D() description for more information about algorithm performance.

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1dinv(/* Complex */ ae_vector* a, ae_int_t n, ae_state *_state)
{
    ae_int_t i;


    ae_assert(n>0, "FFTC1DInv: incorrect N!", _state);
    ae_assert(a->cnt>=n, "FFTC1DInv: Length(A)<N!", _state);
    ae_assert(isfinitecvector(a, n, _state), "FFTC1DInv: A contains infinite or NAN values!", _state);
    
    /*
     * Inverse DFT can be expressed in terms of the DFT as
     *
     *     invfft(x) = fft(x')'/N
     *
     * here x' means conj(x).
     */
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_complex[i].y = -a->ptr.p_complex[i].y;
    }
    fftc1d(a, n, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_complex[i].x = a->ptr.p_complex[i].x/n;
        a->ptr.p_complex[i].y = -a->ptr.p_complex[i].y/n;
    }
}


/*************************************************************************
1-dimensional real FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    F   -   DFT of a input array, array[0..N-1]
            F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)

NOTE:
    F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
of  array  is  usually needed. But for convinience subroutine returns full
complex array (with frequencies above N/2), so its result may be  used  by
other FFT-related subroutines.


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1d(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Complex */ ae_vector* f,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t n2;
    ae_int_t idx;
    ae_complex hn;
    ae_complex hmnc;
    ae_complex v;
    ae_vector buf;
    fasttransformplan plan;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(f);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    _fasttransformplan_init(&plan, _state);

    ae_assert(n>0, "FFTR1D: incorrect N!", _state);
    ae_assert(a->cnt>=n, "FFTR1D: Length(A)<N!", _state);
    ae_assert(isfinitevector(a, n, _state), "FFTR1D: A contains infinite or NAN values!", _state);
    
    /*
     * Special cases:
     * * N=1, FFT is just identity transform.
     * * N=2, FFT is simple too
     *
     * After this block we assume that N is strictly greater than 2
     */
    if( n==1 )
    {
        ae_vector_set_length(f, 1, _state);
        f->ptr.p_complex[0] = ae_complex_from_d(a->ptr.p_double[0]);
        ae_frame_leave(_state);
        return;
    }
    if( n==2 )
    {
        ae_vector_set_length(f, 2, _state);
        f->ptr.p_complex[0].x = a->ptr.p_double[0]+a->ptr.p_double[1];
        f->ptr.p_complex[0].y = (double)(0);
        f->ptr.p_complex[1].x = a->ptr.p_double[0]-a->ptr.p_double[1];
        f->ptr.p_complex[1].y = (double)(0);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Choose between odd-size and even-size FFTs
     */
    if( n%2==0 )
    {
        
        /*
         * even-size real FFT, use reduction to the complex task
         */
        n2 = n/2;
        ae_vector_set_length(&buf, n, _state);
        ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,n-1));
        ftcomplexfftplan(n2, 1, &plan, _state);
        ftapplyplan(&plan, &buf, 0, 1, _state);
        ae_vector_set_length(f, n, _state);
        for(i=0; i<=n2; i++)
        {
            idx = 2*(i%n2);
            hn.x = buf.ptr.p_double[idx+0];
            hn.y = buf.ptr.p_double[idx+1];
            idx = 2*((n2-i)%n2);
            hmnc.x = buf.ptr.p_double[idx+0];
            hmnc.y = -buf.ptr.p_double[idx+1];
            v.x = -ae_sin(-2*ae_pi*i/n, _state);
            v.y = ae_cos(-2*ae_pi*i/n, _state);
            f->ptr.p_complex[i] = ae_c_sub(ae_c_add(hn,hmnc),ae_c_mul(v,ae_c_sub(hn,hmnc)));
            f->ptr.p_complex[i].x = 0.5*f->ptr.p_complex[i].x;
            f->ptr.p_complex[i].y = 0.5*f->ptr.p_complex[i].y;
        }
        for(i=n2+1; i<=n-1; i++)
        {
            f->ptr.p_complex[i] = ae_c_conj(f->ptr.p_complex[n-i], _state);
        }
    }
    else
    {
        
        /*
         * use complex FFT
         */
        ae_vector_set_length(f, n, _state);
        for(i=0; i<=n-1; i++)
        {
            f->ptr.p_complex[i] = ae_complex_from_d(a->ptr.p_double[i]);
        }
        fftc1d(f, n, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional real inverse FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    F   -   array[0..floor(N/2)] - frequencies from forward real FFT
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]

NOTE:
    F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
F[floor(N/2)] has no special properties.

Relying on properties noted above, FFTR1DInv subroutine uses only elements
from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
N is even it ignores imaginary part of F[floor(N/2)] too.

When you call this function using full arguments list - "FFTR1DInv(F,N,A)"
- you can pass either either frequencies array with N elements or  reduced
array with roughly N/2 elements - subroutine will  successfully  transform
both.

If you call this function using reduced arguments list -  "FFTR1DInv(F,A)"
- you must pass FULL array with N elements (although higher  N/2 are still
not used) because array size is used to automatically determine FFT length


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinv(/* Complex */ ae_vector* f,
     ae_int_t n,
     /* Real    */ ae_vector* a,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector h;
    ae_vector fh;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(a);
    ae_vector_init(&h, 0, DT_REAL, _state);
    ae_vector_init(&fh, 0, DT_COMPLEX, _state);

    ae_assert(n>0, "FFTR1DInv: incorrect N!", _state);
    ae_assert(f->cnt>=ae_ifloor((double)n/(double)2, _state)+1, "FFTR1DInv: Length(F)<Floor(N/2)+1!", _state);
    ae_assert(ae_isfinite(f->ptr.p_complex[0].x, _state), "FFTR1DInv: F contains infinite or NAN values!", _state);
    for(i=1; i<=ae_ifloor((double)n/(double)2, _state)-1; i++)
    {
        ae_assert(ae_isfinite(f->ptr.p_complex[i].x, _state)&&ae_isfinite(f->ptr.p_complex[i].y, _state), "FFTR1DInv: F contains infinite or NAN values!", _state);
    }
    ae_assert(ae_isfinite(f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].x, _state), "FFTR1DInv: F contains infinite or NAN values!", _state);
    if( n%2!=0 )
    {
        ae_assert(ae_isfinite(f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].y, _state), "FFTR1DInv: F contains infinite or NAN values!", _state);
    }
    
    /*
     * Special case: N=1, FFT is just identity transform.
     * After this block we assume that N is strictly greater than 1.
     */
    if( n==1 )
    {
        ae_vector_set_length(a, 1, _state);
        a->ptr.p_double[0] = f->ptr.p_complex[0].x;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * inverse real FFT is reduced to the inverse real FHT,
     * which is reduced to the forward real FHT,
     * which is reduced to the forward real FFT.
     *
     * Don't worry, it is really compact and efficient reduction :)
     */
    ae_vector_set_length(&h, n, _state);
    ae_vector_set_length(a, n, _state);
    h.ptr.p_double[0] = f->ptr.p_complex[0].x;
    for(i=1; i<=ae_ifloor((double)n/(double)2, _state)-1; i++)
    {
        h.ptr.p_double[i] = f->ptr.p_complex[i].x-f->ptr.p_complex[i].y;
        h.ptr.p_double[n-i] = f->ptr.p_complex[i].x+f->ptr.p_complex[i].y;
    }
    if( n%2==0 )
    {
        h.ptr.p_double[ae_ifloor((double)n/(double)2, _state)] = f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].x;
    }
    else
    {
        h.ptr.p_double[ae_ifloor((double)n/(double)2, _state)] = f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].x-f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].y;
        h.ptr.p_double[ae_ifloor((double)n/(double)2, _state)+1] = f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].x+f->ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].y;
    }
    fftr1d(&h, n, &fh, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_double[i] = (fh.ptr.p_complex[i].x-fh.ptr.p_complex[i].y)/n;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine. Never call it directly!


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinternaleven(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Real    */ ae_vector* buf,
     fasttransformplan* plan,
     ae_state *_state)
{
    double x;
    double y;
    ae_int_t i;
    ae_int_t n2;
    ae_int_t idx;
    ae_complex hn;
    ae_complex hmnc;
    ae_complex v;


    ae_assert(n>0&&n%2==0, "FFTR1DEvenInplace: incorrect N!", _state);
    
    /*
     * Special cases:
     * * N=2
     *
     * After this block we assume that N is strictly greater than 2
     */
    if( n==2 )
    {
        x = a->ptr.p_double[0]+a->ptr.p_double[1];
        y = a->ptr.p_double[0]-a->ptr.p_double[1];
        a->ptr.p_double[0] = x;
        a->ptr.p_double[1] = y;
        return;
    }
    
    /*
     * even-size real FFT, use reduction to the complex task
     */
    n2 = n/2;
    ae_v_move(&buf->ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ftapplyplan(plan, buf, 0, 1, _state);
    a->ptr.p_double[0] = buf->ptr.p_double[0]+buf->ptr.p_double[1];
    for(i=1; i<=n2-1; i++)
    {
        idx = 2*(i%n2);
        hn.x = buf->ptr.p_double[idx+0];
        hn.y = buf->ptr.p_double[idx+1];
        idx = 2*(n2-i);
        hmnc.x = buf->ptr.p_double[idx+0];
        hmnc.y = -buf->ptr.p_double[idx+1];
        v.x = -ae_sin(-2*ae_pi*i/n, _state);
        v.y = ae_cos(-2*ae_pi*i/n, _state);
        v = ae_c_sub(ae_c_add(hn,hmnc),ae_c_mul(v,ae_c_sub(hn,hmnc)));
        a->ptr.p_double[2*i+0] = 0.5*v.x;
        a->ptr.p_double[2*i+1] = 0.5*v.y;
    }
    a->ptr.p_double[1] = buf->ptr.p_double[0]-buf->ptr.p_double[1];
}


/*************************************************************************
Internal subroutine. Never call it directly!


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinvinternaleven(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Real    */ ae_vector* buf,
     fasttransformplan* plan,
     ae_state *_state)
{
    double x;
    double y;
    double t;
    ae_int_t i;
    ae_int_t n2;


    ae_assert(n>0&&n%2==0, "FFTR1DInvInternalEven: incorrect N!", _state);
    
    /*
     * Special cases:
     * * N=2
     *
     * After this block we assume that N is strictly greater than 2
     */
    if( n==2 )
    {
        x = 0.5*(a->ptr.p_double[0]+a->ptr.p_double[1]);
        y = 0.5*(a->ptr.p_double[0]-a->ptr.p_double[1]);
        a->ptr.p_double[0] = x;
        a->ptr.p_double[1] = y;
        return;
    }
    
    /*
     * inverse real FFT is reduced to the inverse real FHT,
     * which is reduced to the forward real FHT,
     * which is reduced to the forward real FFT.
     *
     * Don't worry, it is really compact and efficient reduction :)
     */
    n2 = n/2;
    buf->ptr.p_double[0] = a->ptr.p_double[0];
    for(i=1; i<=n2-1; i++)
    {
        x = a->ptr.p_double[2*i+0];
        y = a->ptr.p_double[2*i+1];
        buf->ptr.p_double[i] = x-y;
        buf->ptr.p_double[n-i] = x+y;
    }
    buf->ptr.p_double[n2] = a->ptr.p_double[1];
    fftr1dinternaleven(buf, n, a, plan, _state);
    a->ptr.p_double[0] = buf->ptr.p_double[0]/n;
    t = (double)1/(double)n;
    for(i=1; i<=n2-1; i++)
    {
        x = buf->ptr.p_double[2*i+0];
        y = buf->ptr.p_double[2*i+1];
        a->ptr.p_double[i] = t*(x-y);
        a->ptr.p_double[n-i] = t*(x+y);
    }
    a->ptr.p_double[n2] = buf->ptr.p_double[1]/n;
}




/*************************************************************************
1-dimensional complex convolution.

For given A/B returns conv(A,B) (non-circular). Subroutine can automatically
choose between three implementations: straightforward O(M*N)  formula  for
very small N (or M), overlap-add algorithm for  cases  where  max(M,N)  is
significantly larger than min(M,N), but O(M*N) algorithm is too slow,  and
general FFT-based formula for cases where two previois algorithms are  too
slow.

Algorithm has max(M,N)*log(max(M,N)) complexity for any M/N.

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1d(/* Complex */ ae_vector* a,
     ae_int_t m,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Complex */ ae_vector* r,
     ae_state *_state)
{

    ae_vector_clear(r);

    ae_assert(n>0&&m>0, "ConvC1D: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer that B.
     */
    if( m<n )
    {
        convc1d(b, n, a, m, r, _state);
        return;
    }
    convc1dx(a, m, b, n, ae_false, -1, 0, r, _state);
}


/*************************************************************************
1-dimensional complex non-circular deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dinv(/* Complex */ ae_vector* a,
     ae_int_t m,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Complex */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t p;
    ae_vector buf;
    ae_vector buf2;
    fasttransformplan plan;
    ae_complex c1;
    ae_complex c2;
    ae_complex c3;
    double t;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);
    _fasttransformplan_init(&plan, _state);

    ae_assert((n>0&&m>0)&&n<=m, "ConvC1DInv: incorrect N or M!", _state);
    p = ftbasefindsmooth(m, _state);
    ftcomplexfftplan(p, 1, &plan, _state);
    ae_vector_set_length(&buf, 2*p, _state);
    for(i=0; i<=m-1; i++)
    {
        buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
        buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
    }
    for(i=m; i<=p-1; i++)
    {
        buf.ptr.p_double[2*i+0] = (double)(0);
        buf.ptr.p_double[2*i+1] = (double)(0);
    }
    ae_vector_set_length(&buf2, 2*p, _state);
    for(i=0; i<=n-1; i++)
    {
        buf2.ptr.p_double[2*i+0] = b->ptr.p_complex[i].x;
        buf2.ptr.p_double[2*i+1] = b->ptr.p_complex[i].y;
    }
    for(i=n; i<=p-1; i++)
    {
        buf2.ptr.p_double[2*i+0] = (double)(0);
        buf2.ptr.p_double[2*i+1] = (double)(0);
    }
    ftapplyplan(&plan, &buf, 0, 1, _state);
    ftapplyplan(&plan, &buf2, 0, 1, _state);
    for(i=0; i<=p-1; i++)
    {
        c1.x = buf.ptr.p_double[2*i+0];
        c1.y = buf.ptr.p_double[2*i+1];
        c2.x = buf2.ptr.p_double[2*i+0];
        c2.y = buf2.ptr.p_double[2*i+1];
        c3 = ae_c_div(c1,c2);
        buf.ptr.p_double[2*i+0] = c3.x;
        buf.ptr.p_double[2*i+1] = -c3.y;
    }
    ftapplyplan(&plan, &buf, 0, 1, _state);
    t = (double)1/(double)p;
    ae_vector_set_length(r, m-n+1, _state);
    for(i=0; i<=m-n; i++)
    {
        r->ptr.p_complex[i].x = t*buf.ptr.p_double[2*i+0];
        r->ptr.p_complex[i].y = -t*buf.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional circular complex convolution.

For given S/R returns conv(S,R) (circular). Algorithm has linearithmic
complexity for any M/N.

IMPORTANT:  normal convolution is commutative,  i.e.   it  is symmetric  -
conv(A,B)=conv(B,A).  Cyclic convolution IS NOT.  One function - S - is  a
signal,  periodic function, and another - R - is a response,  non-periodic
function with limited length.

INPUT PARAMETERS
    S   -   array[0..M-1] - complex periodic signal
    M   -   problem size
    B   -   array[0..N-1] - complex non-periodic response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircular(/* Complex */ ae_vector* s,
     ae_int_t m,
     /* Complex */ ae_vector* r,
     ae_int_t n,
     /* Complex */ ae_vector* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(c);
    ae_vector_init(&buf, 0, DT_COMPLEX, _state);

    ae_assert(n>0&&m>0, "ConvC1DCircular: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&buf, m, _state);
        for(i1=0; i1<=m-1; i1++)
        {
            buf.ptr.p_complex[i1] = ae_complex_from_i(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_cadd(&buf.ptr.p_complex[0], 1, &r->ptr.p_complex[i1], 1, "N", ae_v_len(0,j2));
            i1 = i1+m;
        }
        convc1dcircular(s, m, &buf, m, c, _state);
        ae_frame_leave(_state);
        return;
    }
    convc1dx(s, m, r, n, ae_true, -1, 0, c, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional circular complex deconvolution (inverse of ConvC1DCircular()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved periodic signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - non-periodic response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-1].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircularinv(/* Complex */ ae_vector* a,
     ae_int_t m,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Complex */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j2;
    ae_vector buf;
    ae_vector buf2;
    ae_vector cbuf;
    fasttransformplan plan;
    ae_complex c1;
    ae_complex c2;
    ae_complex c3;
    double t;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);
    ae_vector_init(&cbuf, 0, DT_COMPLEX, _state);
    _fasttransformplan_init(&plan, _state);

    ae_assert(n>0&&m>0, "ConvC1DCircularInv: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&cbuf, m, _state);
        for(i=0; i<=m-1; i++)
        {
            cbuf.ptr.p_complex[i] = ae_complex_from_i(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_cadd(&cbuf.ptr.p_complex[0], 1, &b->ptr.p_complex[i1], 1, "N", ae_v_len(0,j2));
            i1 = i1+m;
        }
        convc1dcircularinv(a, m, &cbuf, m, r, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Task is normalized
     */
    ftcomplexfftplan(m, 1, &plan, _state);
    ae_vector_set_length(&buf, 2*m, _state);
    for(i=0; i<=m-1; i++)
    {
        buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
        buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
    }
    ae_vector_set_length(&buf2, 2*m, _state);
    for(i=0; i<=n-1; i++)
    {
        buf2.ptr.p_double[2*i+0] = b->ptr.p_complex[i].x;
        buf2.ptr.p_double[2*i+1] = b->ptr.p_complex[i].y;
    }
    for(i=n; i<=m-1; i++)
    {
        buf2.ptr.p_double[2*i+0] = (double)(0);
        buf2.ptr.p_double[2*i+1] = (double)(0);
    }
    ftapplyplan(&plan, &buf, 0, 1, _state);
    ftapplyplan(&plan, &buf2, 0, 1, _state);
    for(i=0; i<=m-1; i++)
    {
        c1.x = buf.ptr.p_double[2*i+0];
        c1.y = buf.ptr.p_double[2*i+1];
        c2.x = buf2.ptr.p_double[2*i+0];
        c2.y = buf2.ptr.p_double[2*i+1];
        c3 = ae_c_div(c1,c2);
        buf.ptr.p_double[2*i+0] = c3.x;
        buf.ptr.p_double[2*i+1] = -c3.y;
    }
    ftapplyplan(&plan, &buf, 0, 1, _state);
    t = (double)1/(double)m;
    ae_vector_set_length(r, m, _state);
    for(i=0; i<=m-1; i++)
    {
        r->ptr.p_complex[i].x = t*buf.ptr.p_double[2*i+0];
        r->ptr.p_complex[i].y = -t*buf.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional real convolution.

Analogous to ConvC1D(), see ConvC1D() comments for more details.

INPUT PARAMETERS
    A   -   array[0..M-1] - real function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1d(/* Real    */ ae_vector* a,
     ae_int_t m,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{

    ae_vector_clear(r);

    ae_assert(n>0&&m>0, "ConvR1D: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer that B.
     */
    if( m<n )
    {
        convr1d(b, n, a, m, r, _state);
        return;
    }
    convr1dx(a, m, b, n, ae_false, -1, 0, r, _state);
}


/*************************************************************************
1-dimensional real deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dinv(/* Real    */ ae_vector* a,
     ae_int_t m,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t p;
    ae_vector buf;
    ae_vector buf2;
    ae_vector buf3;
    fasttransformplan plan;
    ae_complex c1;
    ae_complex c2;
    ae_complex c3;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);
    ae_vector_init(&buf3, 0, DT_REAL, _state);
    _fasttransformplan_init(&plan, _state);

    ae_assert((n>0&&m>0)&&n<=m, "ConvR1DInv: incorrect N or M!", _state);
    p = ftbasefindsmootheven(m, _state);
    ae_vector_set_length(&buf, p, _state);
    ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1));
    for(i=m; i<=p-1; i++)
    {
        buf.ptr.p_double[i] = (double)(0);
    }
    ae_vector_set_length(&buf2, p, _state);
    ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=n; i<=p-1; i++)
    {
        buf2.ptr.p_double[i] = (double)(0);
    }
    ae_vector_set_length(&buf3, p, _state);
    ftcomplexfftplan(p/2, 1, &plan, _state);
    fftr1dinternaleven(&buf, p, &buf3, &plan, _state);
    fftr1dinternaleven(&buf2, p, &buf3, &plan, _state);
    buf.ptr.p_double[0] = buf.ptr.p_double[0]/buf2.ptr.p_double[0];
    buf.ptr.p_double[1] = buf.ptr.p_double[1]/buf2.ptr.p_double[1];
    for(i=1; i<=p/2-1; i++)
    {
        c1.x = buf.ptr.p_double[2*i+0];
        c1.y = buf.ptr.p_double[2*i+1];
        c2.x = buf2.ptr.p_double[2*i+0];
        c2.y = buf2.ptr.p_double[2*i+1];
        c3 = ae_c_div(c1,c2);
        buf.ptr.p_double[2*i+0] = c3.x;
        buf.ptr.p_double[2*i+1] = c3.y;
    }
    fftr1dinvinternaleven(&buf, p, &buf3, &plan, _state);
    ae_vector_set_length(r, m-n+1, _state);
    ae_v_move(&r->ptr.p_double[0], 1, &buf.ptr.p_double[0], 1, ae_v_len(0,m-n));
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional circular real convolution.

Analogous to ConvC1DCircular(), see ConvC1DCircular() comments for more details.

INPUT PARAMETERS
    S   -   array[0..M-1] - real signal
    M   -   problem size
    B   -   array[0..N-1] - real response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircular(/* Real    */ ae_vector* s,
     ae_int_t m,
     /* Real    */ ae_vector* r,
     ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(c);
    ae_vector_init(&buf, 0, DT_REAL, _state);

    ae_assert(n>0&&m>0, "ConvC1DCircular: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&buf, m, _state);
        for(i1=0; i1<=m-1; i1++)
        {
            buf.ptr.p_double[i1] = (double)(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_add(&buf.ptr.p_double[0], 1, &r->ptr.p_double[i1], 1, ae_v_len(0,j2));
            i1 = i1+m;
        }
        convr1dcircular(s, m, &buf, m, c, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * reduce to usual convolution
     */
    convr1dx(s, m, r, n, ae_true, -1, 0, c, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional complex deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircularinv(/* Real    */ ae_vector* a,
     ae_int_t m,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j2;
    ae_vector buf;
    ae_vector buf2;
    ae_vector buf3;
    ae_vector cbuf;
    ae_vector cbuf2;
    fasttransformplan plan;
    ae_complex c1;
    ae_complex c2;
    ae_complex c3;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);
    ae_vector_init(&buf3, 0, DT_REAL, _state);
    ae_vector_init(&cbuf, 0, DT_COMPLEX, _state);
    ae_vector_init(&cbuf2, 0, DT_COMPLEX, _state);
    _fasttransformplan_init(&plan, _state);

    ae_assert(n>0&&m>0, "ConvR1DCircularInv: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&buf, m, _state);
        for(i=0; i<=m-1; i++)
        {
            buf.ptr.p_double[i] = (double)(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_add(&buf.ptr.p_double[0], 1, &b->ptr.p_double[i1], 1, ae_v_len(0,j2));
            i1 = i1+m;
        }
        convr1dcircularinv(a, m, &buf, m, r, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Task is normalized
     */
    if( m%2==0 )
    {
        
        /*
         * size is even, use fast even-size FFT
         */
        ae_vector_set_length(&buf, m, _state);
        ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1));
        ae_vector_set_length(&buf2, m, _state);
        ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
        for(i=n; i<=m-1; i++)
        {
            buf2.ptr.p_double[i] = (double)(0);
        }
        ae_vector_set_length(&buf3, m, _state);
        ftcomplexfftplan(m/2, 1, &plan, _state);
        fftr1dinternaleven(&buf, m, &buf3, &plan, _state);
        fftr1dinternaleven(&buf2, m, &buf3, &plan, _state);
        buf.ptr.p_double[0] = buf.ptr.p_double[0]/buf2.ptr.p_double[0];
        buf.ptr.p_double[1] = buf.ptr.p_double[1]/buf2.ptr.p_double[1];
        for(i=1; i<=m/2-1; i++)
        {
            c1.x = buf.ptr.p_double[2*i+0];
            c1.y = buf.ptr.p_double[2*i+1];
            c2.x = buf2.ptr.p_double[2*i+0];
            c2.y = buf2.ptr.p_double[2*i+1];
            c3 = ae_c_div(c1,c2);
            buf.ptr.p_double[2*i+0] = c3.x;
            buf.ptr.p_double[2*i+1] = c3.y;
        }
        fftr1dinvinternaleven(&buf, m, &buf3, &plan, _state);
        ae_vector_set_length(r, m, _state);
        ae_v_move(&r->ptr.p_double[0], 1, &buf.ptr.p_double[0], 1, ae_v_len(0,m-1));
    }
    else
    {
        
        /*
         * odd-size, use general real FFT
         */
        fftr1d(a, m, &cbuf, _state);
        ae_vector_set_length(&buf2, m, _state);
        ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
        for(i=n; i<=m-1; i++)
        {
            buf2.ptr.p_double[i] = (double)(0);
        }
        fftr1d(&buf2, m, &cbuf2, _state);
        for(i=0; i<=ae_ifloor((double)m/(double)2, _state); i++)
        {
            cbuf.ptr.p_complex[i] = ae_c_div(cbuf.ptr.p_complex[i],cbuf2.ptr.p_complex[i]);
        }
        fftr1dinv(&cbuf, m, r, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional complex convolution.

Extended subroutine which allows to choose convolution algorithm.
Intended for internal use, ALGLIB users should call ConvC1D()/ConvC1DCircular().

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size, N<=M
    Alg -   algorithm type:
            *-2     auto-select Q for overlap-add
            *-1     auto-select algorithm and parameters
            * 0     straightforward formula for small N's
            * 1     general FFT-based code
            * 2     overlap-add with length Q
    Q   -   length for overlap-add

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-1].

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dx(/* Complex */ ae_vector* a,
     ae_int_t m,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     ae_bool circular,
     ae_int_t alg,
     ae_int_t q,
     /* Complex */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t p;
    ae_int_t ptotal;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j1;
    ae_int_t j2;
    ae_vector bbuf;
    ae_complex v;
    double ax;
    double ay;
    double bx;
    double by;
    double t;
    double tx;
    double ty;
    double flopcand;
    double flopbest;
    ae_int_t algbest;
    fasttransformplan plan;
    ae_vector buf;
    ae_vector buf2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&bbuf, 0, DT_COMPLEX, _state);
    _fasttransformplan_init(&plan, _state);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);

    ae_assert(n>0&&m>0, "ConvC1DX: incorrect N or M!", _state);
    ae_assert(n<=m, "ConvC1DX: N<M assumption is false!", _state);
    
    /*
     * Auto-select
     */
    if( alg==-1||alg==-2 )
    {
        
        /*
         * Initial candidate: straightforward implementation.
         *
         * If we want to use auto-fitted overlap-add,
         * flop count is initialized by large real number - to force
         * another algorithm selection
         */
        algbest = 0;
        if( alg==-1 )
        {
            flopbest = (double)(2*m*n);
        }
        else
        {
            flopbest = ae_maxrealnumber;
        }
        
        /*
         * Another candidate - generic FFT code
         */
        if( alg==-1 )
        {
            if( circular&&ftbaseissmooth(m, _state) )
            {
                
                /*
                 * special code for circular convolution of a sequence with a smooth length
                 */
                flopcand = 3*ftbasegetflopestimate(m, _state)+6*m;
                if( ae_fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
            else
            {
                
                /*
                 * general cyclic/non-cyclic convolution
                 */
                p = ftbasefindsmooth(m+n-1, _state);
                flopcand = 3*ftbasegetflopestimate(p, _state)+6*p;
                if( ae_fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
        }
        
        /*
         * Another candidate - overlap-add
         */
        q = 1;
        ptotal = 1;
        while(ptotal<n)
        {
            ptotal = ptotal*2;
        }
        while(ptotal<=m+n-1)
        {
            p = ptotal-n+1;
            flopcand = ae_iceil((double)m/(double)p, _state)*(2*ftbasegetflopestimate(ptotal, _state)+8*ptotal);
            if( ae_fp_less(flopcand,flopbest) )
            {
                flopbest = flopcand;
                algbest = 2;
                q = p;
            }
            ptotal = ptotal*2;
        }
        alg = algbest;
        convc1dx(a, m, b, n, circular, alg, q, r, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * straightforward formula for
     * circular and non-circular convolutions.
     *
     * Very simple code, no further comments needed.
     */
    if( alg==0 )
    {
        
        /*
         * Special case: N=1
         */
        if( n==1 )
        {
            ae_vector_set_length(r, m, _state);
            v = b->ptr.p_complex[0];
            ae_v_cmovec(&r->ptr.p_complex[0], 1, &a->ptr.p_complex[0], 1, "N", ae_v_len(0,m-1), v);
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * use straightforward formula
         */
        if( circular )
        {
            
            /*
             * circular convolution
             */
            ae_vector_set_length(r, m, _state);
            v = b->ptr.p_complex[0];
            ae_v_cmovec(&r->ptr.p_complex[0], 1, &a->ptr.p_complex[0], 1, "N", ae_v_len(0,m-1), v);
            for(i=1; i<=n-1; i++)
            {
                v = b->ptr.p_complex[i];
                i1 = 0;
                i2 = i-1;
                j1 = m-i;
                j2 = m-1;
                ae_v_caddc(&r->ptr.p_complex[i1], 1, &a->ptr.p_complex[j1], 1, "N", ae_v_len(i1,i2), v);
                i1 = i;
                i2 = m-1;
                j1 = 0;
                j2 = m-i-1;
                ae_v_caddc(&r->ptr.p_complex[i1], 1, &a->ptr.p_complex[j1], 1, "N", ae_v_len(i1,i2), v);
            }
        }
        else
        {
            
            /*
             * non-circular convolution
             */
            ae_vector_set_length(r, m+n-1, _state);
            for(i=0; i<=m+n-2; i++)
            {
                r->ptr.p_complex[i] = ae_complex_from_i(0);
            }
            for(i=0; i<=n-1; i++)
            {
                v = b->ptr.p_complex[i];
                ae_v_caddc(&r->ptr.p_complex[i], 1, &a->ptr.p_complex[0], 1, "N", ae_v_len(i,i+m-1), v);
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * general FFT-based code for
     * circular and non-circular convolutions.
     *
     * First, if convolution is circular, we test whether M is smooth or not.
     * If it is smooth, we just use M-length FFT to calculate convolution.
     * If it is not, we calculate non-circular convolution and wrap it arount.
     *
     * IF convolution is non-circular, we use zero-padding + FFT.
     */
    if( alg==1 )
    {
        if( circular&&ftbaseissmooth(m, _state) )
        {
            
            /*
             * special code for circular convolution with smooth M
             */
            ftcomplexfftplan(m, 1, &plan, _state);
            ae_vector_set_length(&buf, 2*m, _state);
            for(i=0; i<=m-1; i++)
            {
                buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
                buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
            }
            ae_vector_set_length(&buf2, 2*m, _state);
            for(i=0; i<=n-1; i++)
            {
                buf2.ptr.p_double[2*i+0] = b->ptr.p_complex[i].x;
                buf2.ptr.p_double[2*i+1] = b->ptr.p_complex[i].y;
            }
            for(i=n; i<=m-1; i++)
            {
                buf2.ptr.p_double[2*i+0] = (double)(0);
                buf2.ptr.p_double[2*i+1] = (double)(0);
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            ftapplyplan(&plan, &buf2, 0, 1, _state);
            for(i=0; i<=m-1; i++)
            {
                ax = buf.ptr.p_double[2*i+0];
                ay = buf.ptr.p_double[2*i+1];
                bx = buf2.ptr.p_double[2*i+0];
                by = buf2.ptr.p_double[2*i+1];
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*i+0] = tx;
                buf.ptr.p_double[2*i+1] = -ty;
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            t = (double)1/(double)m;
            ae_vector_set_length(r, m, _state);
            for(i=0; i<=m-1; i++)
            {
                r->ptr.p_complex[i].x = t*buf.ptr.p_double[2*i+0];
                r->ptr.p_complex[i].y = -t*buf.ptr.p_double[2*i+1];
            }
        }
        else
        {
            
            /*
             * M is non-smooth, general code (circular/non-circular):
             * * first part is the same for circular and non-circular
             *   convolutions. zero padding, FFTs, inverse FFTs
             * * second part differs:
             *   * for non-circular convolution we just copy array
             *   * for circular convolution we add array tail to its head
             */
            p = ftbasefindsmooth(m+n-1, _state);
            ftcomplexfftplan(p, 1, &plan, _state);
            ae_vector_set_length(&buf, 2*p, _state);
            for(i=0; i<=m-1; i++)
            {
                buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
                buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
            }
            for(i=m; i<=p-1; i++)
            {
                buf.ptr.p_double[2*i+0] = (double)(0);
                buf.ptr.p_double[2*i+1] = (double)(0);
            }
            ae_vector_set_length(&buf2, 2*p, _state);
            for(i=0; i<=n-1; i++)
            {
                buf2.ptr.p_double[2*i+0] = b->ptr.p_complex[i].x;
                buf2.ptr.p_double[2*i+1] = b->ptr.p_complex[i].y;
            }
            for(i=n; i<=p-1; i++)
            {
                buf2.ptr.p_double[2*i+0] = (double)(0);
                buf2.ptr.p_double[2*i+1] = (double)(0);
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            ftapplyplan(&plan, &buf2, 0, 1, _state);
            for(i=0; i<=p-1; i++)
            {
                ax = buf.ptr.p_double[2*i+0];
                ay = buf.ptr.p_double[2*i+1];
                bx = buf2.ptr.p_double[2*i+0];
                by = buf2.ptr.p_double[2*i+1];
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*i+0] = tx;
                buf.ptr.p_double[2*i+1] = -ty;
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            t = (double)1/(double)p;
            if( circular )
            {
                
                /*
                 * circular, add tail to head
                 */
                ae_vector_set_length(r, m, _state);
                for(i=0; i<=m-1; i++)
                {
                    r->ptr.p_complex[i].x = t*buf.ptr.p_double[2*i+0];
                    r->ptr.p_complex[i].y = -t*buf.ptr.p_double[2*i+1];
                }
                for(i=m; i<=m+n-2; i++)
                {
                    r->ptr.p_complex[i-m].x = r->ptr.p_complex[i-m].x+t*buf.ptr.p_double[2*i+0];
                    r->ptr.p_complex[i-m].y = r->ptr.p_complex[i-m].y-t*buf.ptr.p_double[2*i+1];
                }
            }
            else
            {
                
                /*
                 * non-circular, just copy
                 */
                ae_vector_set_length(r, m+n-1, _state);
                for(i=0; i<=m+n-2; i++)
                {
                    r->ptr.p_complex[i].x = t*buf.ptr.p_double[2*i+0];
                    r->ptr.p_complex[i].y = -t*buf.ptr.p_double[2*i+1];
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * overlap-add method for
     * circular and non-circular convolutions.
     *
     * First part of code (separate FFTs of input blocks) is the same
     * for all types of convolution. Second part (overlapping outputs)
     * differs for different types of convolution. We just copy output
     * when convolution is non-circular. We wrap it around, if it is
     * circular.
     */
    if( alg==2 )
    {
        ae_vector_set_length(&buf, 2*(q+n-1), _state);
        
        /*
         * prepare R
         */
        if( circular )
        {
            ae_vector_set_length(r, m, _state);
            for(i=0; i<=m-1; i++)
            {
                r->ptr.p_complex[i] = ae_complex_from_i(0);
            }
        }
        else
        {
            ae_vector_set_length(r, m+n-1, _state);
            for(i=0; i<=m+n-2; i++)
            {
                r->ptr.p_complex[i] = ae_complex_from_i(0);
            }
        }
        
        /*
         * pre-calculated FFT(B)
         */
        ae_vector_set_length(&bbuf, q+n-1, _state);
        ae_v_cmove(&bbuf.ptr.p_complex[0], 1, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
        for(j=n; j<=q+n-2; j++)
        {
            bbuf.ptr.p_complex[j] = ae_complex_from_i(0);
        }
        fftc1d(&bbuf, q+n-1, _state);
        
        /*
         * prepare FFT plan for chunks of A
         */
        ftcomplexfftplan(q+n-1, 1, &plan, _state);
        
        /*
         * main overlap-add cycle
         */
        i = 0;
        while(i<=m-1)
        {
            p = ae_minint(q, m-i, _state);
            for(j=0; j<=p-1; j++)
            {
                buf.ptr.p_double[2*j+0] = a->ptr.p_complex[i+j].x;
                buf.ptr.p_double[2*j+1] = a->ptr.p_complex[i+j].y;
            }
            for(j=p; j<=q+n-2; j++)
            {
                buf.ptr.p_double[2*j+0] = (double)(0);
                buf.ptr.p_double[2*j+1] = (double)(0);
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            for(j=0; j<=q+n-2; j++)
            {
                ax = buf.ptr.p_double[2*j+0];
                ay = buf.ptr.p_double[2*j+1];
                bx = bbuf.ptr.p_complex[j].x;
                by = bbuf.ptr.p_complex[j].y;
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*j+0] = tx;
                buf.ptr.p_double[2*j+1] = -ty;
            }
            ftapplyplan(&plan, &buf, 0, 1, _state);
            t = (double)1/(double)(q+n-1);
            if( circular )
            {
                j1 = ae_minint(i+p+n-2, m-1, _state)-i;
                j2 = j1+1;
            }
            else
            {
                j1 = p+n-2;
                j2 = j1+1;
            }
            for(j=0; j<=j1; j++)
            {
                r->ptr.p_complex[i+j].x = r->ptr.p_complex[i+j].x+buf.ptr.p_double[2*j+0]*t;
                r->ptr.p_complex[i+j].y = r->ptr.p_complex[i+j].y-buf.ptr.p_double[2*j+1]*t;
            }
            for(j=j2; j<=p+n-2; j++)
            {
                r->ptr.p_complex[j-j2].x = r->ptr.p_complex[j-j2].x+buf.ptr.p_double[2*j+0]*t;
                r->ptr.p_complex[j-j2].y = r->ptr.p_complex[j-j2].y-buf.ptr.p_double[2*j+1]*t;
            }
            i = i+p;
        }
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional real convolution.

Extended subroutine which allows to choose convolution algorithm.
Intended for internal use, ALGLIB users should call ConvR1D().

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size, N<=M
    Alg -   algorithm type:
            *-2     auto-select Q for overlap-add
            *-1     auto-select algorithm and parameters
            * 0     straightforward formula for small N's
            * 1     general FFT-based code
            * 2     overlap-add with length Q
    Q   -   length for overlap-add

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-1].

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dx(/* Real    */ ae_vector* a,
     ae_int_t m,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     ae_bool circular,
     ae_int_t alg,
     ae_int_t q,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    double v;
    ae_int_t i;
    ae_int_t j;
    ae_int_t p;
    ae_int_t ptotal;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t j1;
    ae_int_t j2;
    double ax;
    double ay;
    double bx;
    double by;
    double tx;
    double ty;
    double flopcand;
    double flopbest;
    ae_int_t algbest;
    fasttransformplan plan;
    ae_vector buf;
    ae_vector buf2;
    ae_vector buf3;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    _fasttransformplan_init(&plan, _state);
    ae_vector_init(&buf, 0, DT_REAL, _state);
    ae_vector_init(&buf2, 0, DT_REAL, _state);
    ae_vector_init(&buf3, 0, DT_REAL, _state);

    ae_assert(n>0&&m>0, "ConvC1DX: incorrect N or M!", _state);
    ae_assert(n<=m, "ConvC1DX: N<M assumption is false!", _state);
    
    /*
     * handle special cases
     */
    if( ae_minint(m, n, _state)<=2 )
    {
        alg = 0;
    }
    
    /*
     * Auto-select
     */
    if( alg<0 )
    {
        
        /*
         * Initial candidate: straightforward implementation.
         *
         * If we want to use auto-fitted overlap-add,
         * flop count is initialized by large real number - to force
         * another algorithm selection
         */
        algbest = 0;
        if( alg==-1 )
        {
            flopbest = 0.15*m*n;
        }
        else
        {
            flopbest = ae_maxrealnumber;
        }
        
        /*
         * Another candidate - generic FFT code
         */
        if( alg==-1 )
        {
            if( (circular&&ftbaseissmooth(m, _state))&&m%2==0 )
            {
                
                /*
                 * special code for circular convolution of a sequence with a smooth length
                 */
                flopcand = 3*ftbasegetflopestimate(m/2, _state)+(double)(6*m)/(double)2;
                if( ae_fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
            else
            {
                
                /*
                 * general cyclic/non-cyclic convolution
                 */
                p = ftbasefindsmootheven(m+n-1, _state);
                flopcand = 3*ftbasegetflopestimate(p/2, _state)+(double)(6*p)/(double)2;
                if( ae_fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
        }
        
        /*
         * Another candidate - overlap-add
         */
        q = 1;
        ptotal = 1;
        while(ptotal<n)
        {
            ptotal = ptotal*2;
        }
        while(ptotal<=m+n-1)
        {
            p = ptotal-n+1;
            flopcand = ae_iceil((double)m/(double)p, _state)*(2*ftbasegetflopestimate(ptotal/2, _state)+1*(ptotal/2));
            if( ae_fp_less(flopcand,flopbest) )
            {
                flopbest = flopcand;
                algbest = 2;
                q = p;
            }
            ptotal = ptotal*2;
        }
        alg = algbest;
        convr1dx(a, m, b, n, circular, alg, q, r, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * straightforward formula for
     * circular and non-circular convolutions.
     *
     * Very simple code, no further comments needed.
     */
    if( alg==0 )
    {
        
        /*
         * Special case: N=1
         */
        if( n==1 )
        {
            ae_vector_set_length(r, m, _state);
            v = b->ptr.p_double[0];
            ae_v_moved(&r->ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1), v);
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * use straightforward formula
         */
        if( circular )
        {
            
            /*
             * circular convolution
             */
            ae_vector_set_length(r, m, _state);
            v = b->ptr.p_double[0];
            ae_v_moved(&r->ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1), v);
            for(i=1; i<=n-1; i++)
            {
                v = b->ptr.p_double[i];
                i1 = 0;
                i2 = i-1;
                j1 = m-i;
                j2 = m-1;
                ae_v_addd(&r->ptr.p_double[i1], 1, &a->ptr.p_double[j1], 1, ae_v_len(i1,i2), v);
                i1 = i;
                i2 = m-1;
                j1 = 0;
                j2 = m-i-1;
                ae_v_addd(&r->ptr.p_double[i1], 1, &a->ptr.p_double[j1], 1, ae_v_len(i1,i2), v);
            }
        }
        else
        {
            
            /*
             * non-circular convolution
             */
            ae_vector_set_length(r, m+n-1, _state);
            for(i=0; i<=m+n-2; i++)
            {
                r->ptr.p_double[i] = (double)(0);
            }
            for(i=0; i<=n-1; i++)
            {
                v = b->ptr.p_double[i];
                ae_v_addd(&r->ptr.p_double[i], 1, &a->ptr.p_double[0], 1, ae_v_len(i,i+m-1), v);
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * general FFT-based code for
     * circular and non-circular convolutions.
     *
     * First, if convolution is circular, we test whether M is smooth or not.
     * If it is smooth, we just use M-length FFT to calculate convolution.
     * If it is not, we calculate non-circular convolution and wrap it arount.
     *
     * If convolution is non-circular, we use zero-padding + FFT.
     *
     * We assume that M+N-1>2 - we should call small case code otherwise
     */
    if( alg==1 )
    {
        ae_assert(m+n-1>2, "ConvR1DX: internal error!", _state);
        if( (circular&&ftbaseissmooth(m, _state))&&m%2==0 )
        {
            
            /*
             * special code for circular convolution with smooth even M
             */
            ae_vector_set_length(&buf, m, _state);
            ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1));
            ae_vector_set_length(&buf2, m, _state);
            ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
            for(i=n; i<=m-1; i++)
            {
                buf2.ptr.p_double[i] = (double)(0);
            }
            ae_vector_set_length(&buf3, m, _state);
            ftcomplexfftplan(m/2, 1, &plan, _state);
            fftr1dinternaleven(&buf, m, &buf3, &plan, _state);
            fftr1dinternaleven(&buf2, m, &buf3, &plan, _state);
            buf.ptr.p_double[0] = buf.ptr.p_double[0]*buf2.ptr.p_double[0];
            buf.ptr.p_double[1] = buf.ptr.p_double[1]*buf2.ptr.p_double[1];
            for(i=1; i<=m/2-1; i++)
            {
                ax = buf.ptr.p_double[2*i+0];
                ay = buf.ptr.p_double[2*i+1];
                bx = buf2.ptr.p_double[2*i+0];
                by = buf2.ptr.p_double[2*i+1];
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*i+0] = tx;
                buf.ptr.p_double[2*i+1] = ty;
            }
            fftr1dinvinternaleven(&buf, m, &buf3, &plan, _state);
            ae_vector_set_length(r, m, _state);
            ae_v_move(&r->ptr.p_double[0], 1, &buf.ptr.p_double[0], 1, ae_v_len(0,m-1));
        }
        else
        {
            
            /*
             * M is non-smooth or non-even, general code (circular/non-circular):
             * * first part is the same for circular and non-circular
             *   convolutions. zero padding, FFTs, inverse FFTs
             * * second part differs:
             *   * for non-circular convolution we just copy array
             *   * for circular convolution we add array tail to its head
             */
            p = ftbasefindsmootheven(m+n-1, _state);
            ae_vector_set_length(&buf, p, _state);
            ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[0], 1, ae_v_len(0,m-1));
            for(i=m; i<=p-1; i++)
            {
                buf.ptr.p_double[i] = (double)(0);
            }
            ae_vector_set_length(&buf2, p, _state);
            ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
            for(i=n; i<=p-1; i++)
            {
                buf2.ptr.p_double[i] = (double)(0);
            }
            ae_vector_set_length(&buf3, p, _state);
            ftcomplexfftplan(p/2, 1, &plan, _state);
            fftr1dinternaleven(&buf, p, &buf3, &plan, _state);
            fftr1dinternaleven(&buf2, p, &buf3, &plan, _state);
            buf.ptr.p_double[0] = buf.ptr.p_double[0]*buf2.ptr.p_double[0];
            buf.ptr.p_double[1] = buf.ptr.p_double[1]*buf2.ptr.p_double[1];
            for(i=1; i<=p/2-1; i++)
            {
                ax = buf.ptr.p_double[2*i+0];
                ay = buf.ptr.p_double[2*i+1];
                bx = buf2.ptr.p_double[2*i+0];
                by = buf2.ptr.p_double[2*i+1];
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*i+0] = tx;
                buf.ptr.p_double[2*i+1] = ty;
            }
            fftr1dinvinternaleven(&buf, p, &buf3, &plan, _state);
            if( circular )
            {
                
                /*
                 * circular, add tail to head
                 */
                ae_vector_set_length(r, m, _state);
                ae_v_move(&r->ptr.p_double[0], 1, &buf.ptr.p_double[0], 1, ae_v_len(0,m-1));
                if( n>=2 )
                {
                    ae_v_add(&r->ptr.p_double[0], 1, &buf.ptr.p_double[m], 1, ae_v_len(0,n-2));
                }
            }
            else
            {
                
                /*
                 * non-circular, just copy
                 */
                ae_vector_set_length(r, m+n-1, _state);
                ae_v_move(&r->ptr.p_double[0], 1, &buf.ptr.p_double[0], 1, ae_v_len(0,m+n-2));
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * overlap-add method
     */
    if( alg==2 )
    {
        ae_assert((q+n-1)%2==0, "ConvR1DX: internal error!", _state);
        ae_vector_set_length(&buf, q+n-1, _state);
        ae_vector_set_length(&buf2, q+n-1, _state);
        ae_vector_set_length(&buf3, q+n-1, _state);
        ftcomplexfftplan((q+n-1)/2, 1, &plan, _state);
        
        /*
         * prepare R
         */
        if( circular )
        {
            ae_vector_set_length(r, m, _state);
            for(i=0; i<=m-1; i++)
            {
                r->ptr.p_double[i] = (double)(0);
            }
        }
        else
        {
            ae_vector_set_length(r, m+n-1, _state);
            for(i=0; i<=m+n-2; i++)
            {
                r->ptr.p_double[i] = (double)(0);
            }
        }
        
        /*
         * pre-calculated FFT(B)
         */
        ae_v_move(&buf2.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
        for(j=n; j<=q+n-2; j++)
        {
            buf2.ptr.p_double[j] = (double)(0);
        }
        fftr1dinternaleven(&buf2, q+n-1, &buf3, &plan, _state);
        
        /*
         * main overlap-add cycle
         */
        i = 0;
        while(i<=m-1)
        {
            p = ae_minint(q, m-i, _state);
            ae_v_move(&buf.ptr.p_double[0], 1, &a->ptr.p_double[i], 1, ae_v_len(0,p-1));
            for(j=p; j<=q+n-2; j++)
            {
                buf.ptr.p_double[j] = (double)(0);
            }
            fftr1dinternaleven(&buf, q+n-1, &buf3, &plan, _state);
            buf.ptr.p_double[0] = buf.ptr.p_double[0]*buf2.ptr.p_double[0];
            buf.ptr.p_double[1] = buf.ptr.p_double[1]*buf2.ptr.p_double[1];
            for(j=1; j<=(q+n-1)/2-1; j++)
            {
                ax = buf.ptr.p_double[2*j+0];
                ay = buf.ptr.p_double[2*j+1];
                bx = buf2.ptr.p_double[2*j+0];
                by = buf2.ptr.p_double[2*j+1];
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf.ptr.p_double[2*j+0] = tx;
                buf.ptr.p_double[2*j+1] = ty;
            }
            fftr1dinvinternaleven(&buf, q+n-1, &buf3, &plan, _state);
            if( circular )
            {
                j1 = ae_minint(i+p+n-2, m-1, _state)-i;
                j2 = j1+1;
            }
            else
            {
                j1 = p+n-2;
                j2 = j1+1;
            }
            ae_v_add(&r->ptr.p_double[i], 1, &buf.ptr.p_double[0], 1, ae_v_len(i,i+j1));
            if( p+n-2>=j2 )
            {
                ae_v_add(&r->ptr.p_double[0], 1, &buf.ptr.p_double[j2], 1, ae_v_len(0,p+n-2-j2));
            }
            i = i+p;
        }
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}




/*************************************************************************
1-dimensional complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrC1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(conj(pattern[j])*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1d(/* Complex */ ae_vector* signal,
     ae_int_t n,
     /* Complex */ ae_vector* pattern,
     ae_int_t m,
     /* Complex */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector p;
    ae_vector b;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&p, 0, DT_COMPLEX, _state);
    ae_vector_init(&b, 0, DT_COMPLEX, _state);

    ae_assert(n>0&&m>0, "CorrC1D: incorrect N or M!", _state);
    ae_vector_set_length(&p, m, _state);
    for(i=0; i<=m-1; i++)
    {
        p.ptr.p_complex[m-1-i] = ae_c_conj(pattern->ptr.p_complex[i], _state);
    }
    convc1d(&p, m, signal, n, &b, _state);
    ae_vector_set_length(r, m+n-1, _state);
    ae_v_cmove(&r->ptr.p_complex[0], 1, &b.ptr.p_complex[m-1], 1, "N", ae_v_len(0,n-1));
    if( m+n-2>=n )
    {
        ae_v_cmove(&r->ptr.p_complex[n], 1, &b.ptr.p_complex[0], 1, "N", ae_v_len(n,m+n-2));
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional circular complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrC1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1dcircular(/* Complex */ ae_vector* signal,
     ae_int_t m,
     /* Complex */ ae_vector* pattern,
     ae_int_t n,
     /* Complex */ ae_vector* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector p;
    ae_vector b;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t i;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(c);
    ae_vector_init(&p, 0, DT_COMPLEX, _state);
    ae_vector_init(&b, 0, DT_COMPLEX, _state);

    ae_assert(n>0&&m>0, "ConvC1DCircular: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&b, m, _state);
        for(i1=0; i1<=m-1; i1++)
        {
            b.ptr.p_complex[i1] = ae_complex_from_i(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_cadd(&b.ptr.p_complex[0], 1, &pattern->ptr.p_complex[i1], 1, "N", ae_v_len(0,j2));
            i1 = i1+m;
        }
        corrc1dcircular(signal, m, &b, m, c, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Task is normalized
     */
    ae_vector_set_length(&p, n, _state);
    for(i=0; i<=n-1; i++)
    {
        p.ptr.p_complex[n-1-i] = ae_c_conj(pattern->ptr.p_complex[i], _state);
    }
    convc1dcircular(signal, m, &p, n, &b, _state);
    ae_vector_set_length(c, m, _state);
    ae_v_cmove(&c->ptr.p_complex[0], 1, &b.ptr.p_complex[n-1], 1, "N", ae_v_len(0,m-n));
    if( m-n+1<=m-1 )
    {
        ae_v_cmove(&c->ptr.p_complex[m-n+1], 1, &b.ptr.p_complex[0], 1, "N", ae_v_len(m-n+1,m-1));
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrR1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(pattern[j]*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(pattern[j]*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1d(/* Real    */ ae_vector* signal,
     ae_int_t n,
     /* Real    */ ae_vector* pattern,
     ae_int_t m,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector p;
    ae_vector b;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(r);
    ae_vector_init(&p, 0, DT_REAL, _state);
    ae_vector_init(&b, 0, DT_REAL, _state);

    ae_assert(n>0&&m>0, "CorrR1D: incorrect N or M!", _state);
    ae_vector_set_length(&p, m, _state);
    for(i=0; i<=m-1; i++)
    {
        p.ptr.p_double[m-1-i] = pattern->ptr.p_double[i];
    }
    convr1d(&p, m, signal, n, &b, _state);
    ae_vector_set_length(r, m+n-1, _state);
    ae_v_move(&r->ptr.p_double[0], 1, &b.ptr.p_double[m-1], 1, ae_v_len(0,n-1));
    if( m+n-2>=n )
    {
        ae_v_move(&r->ptr.p_double[n], 1, &b.ptr.p_double[0], 1, ae_v_len(n,m+n-2));
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional circular real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrR1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1dcircular(/* Real    */ ae_vector* signal,
     ae_int_t m,
     /* Real    */ ae_vector* pattern,
     ae_int_t n,
     /* Real    */ ae_vector* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector p;
    ae_vector b;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t i;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(c);
    ae_vector_init(&p, 0, DT_REAL, _state);
    ae_vector_init(&b, 0, DT_REAL, _state);

    ae_assert(n>0&&m>0, "ConvC1DCircular: incorrect N or M!", _state);
    
    /*
     * normalize task: make M>=N,
     * so A will be longer (at least - not shorter) that B.
     */
    if( m<n )
    {
        ae_vector_set_length(&b, m, _state);
        for(i1=0; i1<=m-1; i1++)
        {
            b.ptr.p_double[i1] = (double)(0);
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ae_minint(i1+m-1, n-1, _state);
            j2 = i2-i1;
            ae_v_add(&b.ptr.p_double[0], 1, &pattern->ptr.p_double[i1], 1, ae_v_len(0,j2));
            i1 = i1+m;
        }
        corrr1dcircular(signal, m, &b, m, c, _state);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Task is normalized
     */
    ae_vector_set_length(&p, n, _state);
    for(i=0; i<=n-1; i++)
    {
        p.ptr.p_double[n-1-i] = pattern->ptr.p_double[i];
    }
    convr1dcircular(signal, m, &p, n, &b, _state);
    ae_vector_set_length(c, m, _state);
    ae_v_move(&c->ptr.p_double[0], 1, &b.ptr.p_double[n-1], 1, ae_v_len(0,m-n));
    if( m-n+1<=m-1 )
    {
        ae_v_move(&c->ptr.p_double[m-n+1], 1, &b.ptr.p_double[0], 1, ae_v_len(m-n+1,m-1));
    }
    ae_frame_leave(_state);
}




/*************************************************************************
1-dimensional Fast Hartley Transform.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size
    
OUTPUT PARAMETERS
    A   -   FHT of a input array, array[0..N-1],
            A_out[k] = sum(A_in[j]*(cos(2*pi*j*k/N)+sin(2*pi*j*k/N)), j=0..N-1)


  -- ALGLIB --
     Copyright 04.06.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1d(/* Real    */ ae_vector* a, ae_int_t n, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector fa;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&fa, 0, DT_COMPLEX, _state);

    ae_assert(n>0, "FHTR1D: incorrect N!", _state);
    
    /*
     * Special case: N=1, FHT is just identity transform.
     * After this block we assume that N is strictly greater than 1.
     */
    if( n==1 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Reduce FHt to real FFT
     */
    fftr1d(a, n, &fa, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_double[i] = fa.ptr.p_complex[i].x-fa.ptr.p_complex[i].y;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
1-dimensional inverse FHT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse FHT of a input array, array[0..N-1]


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1dinv(/* Real    */ ae_vector* a, ae_int_t n, ae_state *_state)
{
    ae_int_t i;


    ae_assert(n>0, "FHTR1DInv: incorrect N!", _state);
    
    /*
     * Special case: N=1, iFHT is just identity transform.
     * After this block we assume that N is strictly greater than 1.
     */
    if( n==1 )
    {
        return;
    }
    
    /*
     * Inverse FHT can be expressed in terms of the FHT as
     *
     *     invfht(x) = fht(x)/N
     */
    fhtr1d(a, n, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_double[i] = a->ptr.p_double[i]/n;
    }
}



}

