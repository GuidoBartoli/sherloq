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
#include "alglibmisc.h"

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
Portable high quality random number generator state.
Initialized with HQRNDRandomize() or HQRNDSeed().

Fields:
    S1, S2      -   seed values
    V           -   precomputed value
    MagicV      -   'magic' value used to determine whether State structure
                    was correctly initialized.
*************************************************************************/
_hqrndstate_owner::_hqrndstate_owner()
{
    p_struct = (alglib_impl::hqrndstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::hqrndstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_hqrndstate_init(p_struct, NULL);
}

_hqrndstate_owner::_hqrndstate_owner(const _hqrndstate_owner &rhs)
{
    p_struct = (alglib_impl::hqrndstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::hqrndstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_hqrndstate_init_copy(p_struct, const_cast<alglib_impl::hqrndstate*>(rhs.p_struct), NULL);
}

_hqrndstate_owner& _hqrndstate_owner::operator=(const _hqrndstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_hqrndstate_clear(p_struct);
    alglib_impl::_hqrndstate_init_copy(p_struct, const_cast<alglib_impl::hqrndstate*>(rhs.p_struct), NULL);
    return *this;
}

_hqrndstate_owner::~_hqrndstate_owner()
{
    alglib_impl::_hqrndstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::hqrndstate* _hqrndstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::hqrndstate* _hqrndstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::hqrndstate*>(p_struct);
}
hqrndstate::hqrndstate() : _hqrndstate_owner() 
{
}

hqrndstate::hqrndstate(const hqrndstate &rhs):_hqrndstate_owner(rhs) 
{
}

hqrndstate& hqrndstate::operator=(const hqrndstate &rhs)
{
    if( this==&rhs )
        return *this;
    _hqrndstate_owner::operator=(rhs);
    return *this;
}

hqrndstate::~hqrndstate()
{
}

/*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndrandomize(hqrndstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::hqrndrandomize(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::hqrndseed(s1, s2, const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrnduniformr(const hqrndstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::hqrnduniformr(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This function generates random integer number in [0, N)

1. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()
2. N can be any positive number except for very large numbers:
   * close to 2^31 on 32-bit systems
   * close to 2^62 on 64-bit systems
   An exception will be generated if N is too large.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::hqrnduniformi(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrndnormal(const hqrndstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::hqrndnormal(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndunit2(const hqrndstate &state, double &x, double &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::hqrndunit2(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &x, &y, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndnormal2(const hqrndstate &state, double &x1, double &x2)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::hqrndnormal2(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), &x1, &x2, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************/
double hqrndexponential(const hqrndstate &state, const double lambdav)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::hqrndexponential(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), lambdav, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This function generates  random number from discrete distribution given by
finite sample X.

INPUT PARAMETERS
    State   -   high quality random number generator, must be
                initialized with HQRNDRandomize() or HQRNDSeed().
        X   -   finite sample
        N   -   number of elements to use, N>=1

RESULT
    this function returns one of the X[i] for random i=0..N-1

  -- ALGLIB --
     Copyright 08.11.2011 by Bochkanov Sergey
*************************************************************************/
double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::hqrnddiscrete(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This function generates random number from continuous  distribution  given
by finite sample X.

INPUT PARAMETERS
    State   -   high quality random number generator, must be
                initialized with HQRNDRandomize() or HQRNDSeed().
        X   -   finite sample, array[N] (can be larger, in this  case only
                leading N elements are used). THIS ARRAY MUST BE SORTED BY
                ASCENDING.
        N   -   number of elements to use, N>=1

RESULT
    this function returns random number from continuous distribution which
    tries to approximate X as mush as possible. min(X)<=Result<=max(X).

  -- ALGLIB --
     Copyright 08.11.2011 by Bochkanov Sergey
*************************************************************************/
double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::hqrndcontinuous(const_cast<alglib_impl::hqrndstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************

*************************************************************************/
_kdtree_owner::_kdtree_owner()
{
    p_struct = (alglib_impl::kdtree*)alglib_impl::ae_malloc(sizeof(alglib_impl::kdtree), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_kdtree_init(p_struct, NULL);
}

_kdtree_owner::_kdtree_owner(const _kdtree_owner &rhs)
{
    p_struct = (alglib_impl::kdtree*)alglib_impl::ae_malloc(sizeof(alglib_impl::kdtree), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_kdtree_init_copy(p_struct, const_cast<alglib_impl::kdtree*>(rhs.p_struct), NULL);
}

_kdtree_owner& _kdtree_owner::operator=(const _kdtree_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_kdtree_clear(p_struct);
    alglib_impl::_kdtree_init_copy(p_struct, const_cast<alglib_impl::kdtree*>(rhs.p_struct), NULL);
    return *this;
}

_kdtree_owner::~_kdtree_owner()
{
    alglib_impl::_kdtree_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::kdtree* _kdtree_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::kdtree* _kdtree_owner::c_ptr() const
{
    return const_cast<alglib_impl::kdtree*>(p_struct);
}
kdtree::kdtree() : _kdtree_owner() 
{
}

kdtree::kdtree(const kdtree &rhs):_kdtree_owner(rhs) 
{
}

kdtree& kdtree::operator=(const kdtree &rhs)
{
    if( this==&rhs )
        return *this;
    _kdtree_owner::operator=(rhs);
    return *this;
}

kdtree::~kdtree()
{
}


/*************************************************************************
This function serializes data structure to string.

Important properties of s_out:
* it contains alphanumeric characters, dots, underscores, minus signs
* these symbols are grouped into words, which are separated by spaces
  and Windows-style (CR+LF) newlines
* although  serializer  uses  spaces and CR+LF as separators, you can 
  replace any separator character by arbitrary combination of spaces,
  tabs, Windows or Unix newlines. It allows flexible reformatting  of
  the  string  in  case you want to include it into text or XML file. 
  But you should not insert separators into the middle of the "words"
  nor you should change case of letters.
* s_out can be freely moved between 32-bit and 64-bit systems, little
  and big endian machines, and so on. You can serialize structure  on
  32-bit machine and unserialize it on 64-bit one (or vice versa), or
  serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also 
  serialize  it  in  C++ version of ALGLIB and unserialize in C# one, 
  and vice versa.
*************************************************************************/
void kdtreeserialize(kdtree &obj, std::string &s_out)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;
    alglib_impl::ae_int_t ssize;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_alloc_start(&serializer);
        alglib_impl::kdtreealloc(&serializer, obj.c_ptr(), &state);
        ssize = alglib_impl::ae_serializer_get_alloc_size(&serializer);
        s_out.clear();
        s_out.reserve((size_t)(ssize+1));
        alglib_impl::ae_serializer_sstart_str(&serializer, &s_out);
        alglib_impl::kdtreeserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        if( s_out.length()>(size_t)ssize )
            throw ap_error("ALGLIB: serialization integrity error");
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
}
/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void kdtreeunserialize(std::string &s_in, kdtree &obj)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_ustart_str(&serializer, &s_in);
        alglib_impl::kdtreeunserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
}

/*************************************************************************
KD-tree creation

This subroutine creates KD-tree from set of X-values and optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    N       -   number of points, N>=0.
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree


NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreebuild(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), n, nx, ny, normtype, const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
KD-tree creation

This subroutine creates KD-tree from set of X-values and optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    N       -   number of points, N>=0.
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree


NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = xy.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreebuild(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), n, nx, ny, normtype, const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
KD-tree creation

This  subroutine  creates  KD-tree  from set of X-values, integer tags and
optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    Tags    -   tags, array[0..N-1], contains integer tags associated
                with points.
    N       -   number of points, N>=0
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree

NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreebuildtagged(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), const_cast<alglib_impl::ae_vector*>(tags.c_ptr()), n, nx, ny, normtype, const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
KD-tree creation

This  subroutine  creates  KD-tree  from set of X-values, integer tags and
optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    Tags    -   tags, array[0..N-1], contains integer tags associated
                with points.
    N       -   number of points, N>=0
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree

NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;
    if( (xy.rows()!=tags.length()))
        throw ap_error("Error while calling 'kdtreebuildtagged': looks like one of arguments has wrong size");
    n = xy.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreebuildtagged(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), const_cast<alglib_impl::ae_vector*>(tags.c_ptr()), n, nx, ny, normtype, const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
K-NN query: K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of actual neighbors found (either K or N, if K>N).

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryknn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), k, selfmatch, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
K-NN query: K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of actual neighbors found (either K or N, if K>N).

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;    
    bool selfmatch;

    selfmatch = true;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryknn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), k, selfmatch, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
R-NN query: all points within R-sphere centered at X

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    R           -   radius of sphere (in corresponding norm), R>0
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of neighbors found, >=0

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
actual results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryrnn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), r, selfmatch, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
R-NN query: all points within R-sphere centered at X

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    R           -   radius of sphere (in corresponding norm), R>0
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of neighbors found, >=0

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
actual results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r)
{
    alglib_impl::ae_state _alglib_env_state;    
    bool selfmatch;

    selfmatch = true;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryrnn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), r, selfmatch, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
K-NN query: approximate K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True
    Eps         -   approximation factor, Eps>=0. eps-approximate  nearest
                    neighbor  is  a  neighbor  whose distance from X is at
                    most (1+eps) times distance of true nearest neighbor.

RESULT
    number of actual neighbors found (either K or N, if K>N).

NOTES
    significant performance gain may be achieved only when Eps  is  is  on
    the order of magnitude of 1 or larger.

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryaknn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), k, selfmatch, eps, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
K-NN query: approximate K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True
    Eps         -   approximation factor, Eps>=0. eps-approximate  nearest
                    neighbor  is  a  neighbor  whose distance from X is at
                    most (1+eps) times distance of true nearest neighbor.

RESULT
    number of actual neighbors found (either K or N, if K>N).

NOTES
    significant performance gain may be achieved only when Eps  is  is  on
    the order of magnitude of 1 or larger.

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps)
{
    alglib_impl::ae_state _alglib_env_state;    
    bool selfmatch;

    selfmatch = true;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::kdtreequeryaknn(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), k, selfmatch, eps, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
X-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    X       -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    X       -   rows are filled with X-values

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsx(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_matrix*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
X- and Y-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    XY      -   possibly pre-allocated buffer. If XY is too small to store
                result, it is resized. If size(XY) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    XY      -   rows are filled with points: first NX columns with
                X-values, next NY columns - with Y-values.

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsxy(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Tags from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    Tags    -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    Tags    -   filled with tags associated with points,
                or, when no tags were supplied, with zeros

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultstags(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(tags.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Distances from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    R       -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    R       -   filled with distances (in corresponding norm)

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsdistances(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
X-values from last query; 'interactive' variant for languages like  Python
which   support    constructs   like  "X = KDTreeQueryResultsXI(KDT)"  and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsxi(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_matrix*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
XY-values from last query; 'interactive' variant for languages like Python
which   support    constructs   like "XY = KDTreeQueryResultsXYI(KDT)" and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsxyi(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Tags  from  last  query;  'interactive' variant for languages like  Python
which  support  constructs  like "Tags = KDTreeQueryResultsTagsI(KDT)" and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultstagsi(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(tags.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
Distances from last query; 'interactive' variant for languages like Python
which  support  constructs   like  "R = KDTreeQueryResultsDistancesI(KDT)"
and interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kdtreequeryresultsdistancesi(const_cast<alglib_impl::kdtree*>(kdt.c_ptr()), const_cast<alglib_impl::ae_vector*>(r.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************

*************************************************************************/
_xdebugrecord1_owner::_xdebugrecord1_owner()
{
    p_struct = (alglib_impl::xdebugrecord1*)alglib_impl::ae_malloc(sizeof(alglib_impl::xdebugrecord1), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_xdebugrecord1_init(p_struct, NULL);
}

_xdebugrecord1_owner::_xdebugrecord1_owner(const _xdebugrecord1_owner &rhs)
{
    p_struct = (alglib_impl::xdebugrecord1*)alglib_impl::ae_malloc(sizeof(alglib_impl::xdebugrecord1), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_xdebugrecord1_init_copy(p_struct, const_cast<alglib_impl::xdebugrecord1*>(rhs.p_struct), NULL);
}

_xdebugrecord1_owner& _xdebugrecord1_owner::operator=(const _xdebugrecord1_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_xdebugrecord1_clear(p_struct);
    alglib_impl::_xdebugrecord1_init_copy(p_struct, const_cast<alglib_impl::xdebugrecord1*>(rhs.p_struct), NULL);
    return *this;
}

_xdebugrecord1_owner::~_xdebugrecord1_owner()
{
    alglib_impl::_xdebugrecord1_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::xdebugrecord1* _xdebugrecord1_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::xdebugrecord1* _xdebugrecord1_owner::c_ptr() const
{
    return const_cast<alglib_impl::xdebugrecord1*>(p_struct);
}
xdebugrecord1::xdebugrecord1() : _xdebugrecord1_owner() ,i(p_struct->i),c(*((alglib::complex*)(&p_struct->c))),a(&p_struct->a)
{
}

xdebugrecord1::xdebugrecord1(const xdebugrecord1 &rhs):_xdebugrecord1_owner(rhs) ,i(p_struct->i),c(*((alglib::complex*)(&p_struct->c))),a(&p_struct->a)
{
}

xdebugrecord1& xdebugrecord1::operator=(const xdebugrecord1 &rhs)
{
    if( this==&rhs )
        return *this;
    _xdebugrecord1_owner::operator=(rhs);
    return *this;
}

xdebugrecord1::~xdebugrecord1()
{
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Creates and returns XDebugRecord1 structure:
* integer and complex fields of Rec1 are set to 1 and 1+i correspondingly
* array field of Rec1 is set to [2,3]

  -- ALGLIB --
     Copyright 27.05.2014 by Bochkanov Sergey
*************************************************************************/
void xdebuginitrecord1(xdebugrecord1 &rec1)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebuginitrecord1(const_cast<alglib_impl::xdebugrecord1*>(rec1.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Counts number of True values in the boolean 1D array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugb1count(const boolean_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::xdebugb1count(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by NOT(a[i]).
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1not(const boolean_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb1not(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1appendcopy(boolean_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb1appendcopy(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered elements set to True.
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb1outeven(n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugi1sum(const integer_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::xdebugi1sum(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1neg(const integer_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi1neg(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1appendcopy(integer_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi1appendcopy(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[I] set to I, and odd-numbered
ones set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1outeven(const ae_int_t n, integer_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi1outeven(n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugr1sum(const real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::xdebugr1sum(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1neg(const real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr1neg(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1appendcopy(real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr1appendcopy(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[I] set to I*0.25,
and odd-numbered ones are set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1outeven(const ae_int_t n, real_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr1outeven(n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
alglib::complex xdebugc1sum(const complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_complex result = alglib_impl::xdebugc1sum(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<alglib::complex*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1neg(const complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc1neg(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1appendcopy(complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc1appendcopy(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[K] set to (x,y) = (K*0.25, K*0.125)
and odd-numbered ones are set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1outeven(const ae_int_t n, complex_1d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc1outeven(n, const_cast<alglib_impl::ae_vector*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Counts number of True values in the boolean 2D array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugb2count(const boolean_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::xdebugb2count(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by NOT(a[i]).
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2not(const boolean_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb2not(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2transpose(boolean_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb2transpose(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J)>0"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugb2outsin(m, n, const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugi2sum(const integer_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::xdebugi2sum(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2neg(const integer_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi2neg(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2transpose(integer_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi2transpose(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sign(Sin(3*I+5*J))"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugi2outsin(m, n, const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugr2sum(const real_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::xdebugr2sum(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2neg(const real_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr2neg(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2transpose(real_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr2transpose(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J)"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugr2outsin(m, n, const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
alglib::complex xdebugc2sum(const complex_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_complex result = alglib_impl::xdebugc2sum(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<alglib::complex*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2neg(const complex_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc2neg(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2transpose(complex_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc2transpose(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J),Cos(3*I+5*J)"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::xdebugc2outsincos(m, n, const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}

/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of a[i,j]*(1+b[i,j]) such that c[i,j] is True

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::xdebugmaskedbiasedproductsum(m, n, const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), const_cast<alglib_impl::ae_matrix*>(b.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
static ae_int_t hqrnd_hqrndmax = 2147483561;
static ae_int_t hqrnd_hqrndm1 = 2147483563;
static ae_int_t hqrnd_hqrndm2 = 2147483399;
static ae_int_t hqrnd_hqrndmagic = 1634357784;
static ae_int_t hqrnd_hqrndintegerbase(hqrndstate* state,
     ae_state *_state);


static ae_int_t nearestneighbor_splitnodesize = 6;
static ae_int_t nearestneighbor_kdtreefirstversion = 0;
static void nearestneighbor_kdtreesplit(kdtree* kdt,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t d,
     double s,
     ae_int_t* i3,
     ae_state *_state);
static void nearestneighbor_kdtreegeneratetreerec(kdtree* kdt,
     ae_int_t* nodesoffs,
     ae_int_t* splitsoffs,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t maxleafsize,
     ae_state *_state);
static void nearestneighbor_kdtreequerynnrec(kdtree* kdt,
     ae_int_t offs,
     ae_state *_state);
static void nearestneighbor_kdtreeinitbox(kdtree* kdt,
     /* Real    */ ae_vector* x,
     ae_state *_state);
static void nearestneighbor_kdtreeallocdatasetindependent(kdtree* kdt,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state);
static void nearestneighbor_kdtreeallocdatasetdependent(kdtree* kdt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state);
static void nearestneighbor_kdtreealloctemporaries(kdtree* kdt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state);







/*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndrandomize(hqrndstate* state, ae_state *_state)
{
    ae_int_t s0;
    ae_int_t s1;

    _hqrndstate_clear(state);

    s0 = ae_randominteger(hqrnd_hqrndm1, _state);
    s1 = ae_randominteger(hqrnd_hqrndm2, _state);
    hqrndseed(s0, s1, state, _state);
}


/*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndseed(ae_int_t s1,
     ae_int_t s2,
     hqrndstate* state,
     ae_state *_state)
{

    _hqrndstate_clear(state);

    
    /*
     * Protection against negative seeds:
     *
     *     SEED := -(SEED+1)
     *
     * We can use just "-SEED" because there exists such integer number  N
     * that N<0, -N=N<0 too. (This number is equal to 0x800...000).   Need
     * to handle such seed correctly forces us to use  a  bit  complicated
     * formula.
     */
    if( s1<0 )
    {
        s1 = -(s1+1);
    }
    if( s2<0 )
    {
        s2 = -(s2+1);
    }
    state->s1 = s1%(hqrnd_hqrndm1-1)+1;
    state->s2 = s2%(hqrnd_hqrndm2-1)+1;
    state->magicv = hqrnd_hqrndmagic;
}


/*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrnduniformr(hqrndstate* state, ae_state *_state)
{
    double result;


    result = (double)(hqrnd_hqrndintegerbase(state, _state)+1)/(double)(hqrnd_hqrndmax+2);
    return result;
}


/*************************************************************************
This function generates random integer number in [0, N)

1. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()
2. N can be any positive number except for very large numbers:
   * close to 2^31 on 32-bit systems
   * close to 2^62 on 64-bit systems
   An exception will be generated if N is too large.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t hqrnduniformi(hqrndstate* state, ae_int_t n, ae_state *_state)
{
    ae_int_t maxcnt;
    ae_int_t mx;
    ae_int_t a;
    ae_int_t b;
    ae_int_t result;


    ae_assert(n>0, "HQRNDUniformI: N<=0!", _state);
    maxcnt = hqrnd_hqrndmax+1;
    
    /*
     * Two branches: one for N<=MaxCnt, another for N>MaxCnt.
     */
    if( n>maxcnt )
    {
        
        /*
         * N>=MaxCnt.
         *
         * We have two options here:
         * a) N is exactly divisible by MaxCnt
         * b) N is not divisible by MaxCnt
         *
         * In both cases we reduce problem on interval spanning [0,N)
         * to several subproblems on intervals spanning [0,MaxCnt).
         */
        if( n%maxcnt==0 )
        {
            
            /*
             * N is exactly divisible by MaxCnt.
             *
             * [0,N) range is dividided into N/MaxCnt bins,
             * each of them having length equal to MaxCnt.
             *
             * We generate:
             * * random bin number B
             * * random offset within bin A
             * Both random numbers are generated by recursively
             * calling HQRNDUniformI().
             *
             * Result is equal to A+MaxCnt*B.
             */
            ae_assert(n/maxcnt<=maxcnt, "HQRNDUniformI: N is too large", _state);
            a = hqrnduniformi(state, maxcnt, _state);
            b = hqrnduniformi(state, n/maxcnt, _state);
            result = a+maxcnt*b;
        }
        else
        {
            
            /*
             * N is NOT exactly divisible by MaxCnt.
             *
             * [0,N) range is dividided into Ceil(N/MaxCnt) bins,
             * each of them having length equal to MaxCnt.
             *
             * We generate:
             * * random bin number B in [0, Ceil(N/MaxCnt)-1]
             * * random offset within bin A
             * * if both of what is below is true
             *   1) bin number B is that of the last bin
             *   2) A >= N mod MaxCnt
             *   then we repeat generation of A/B.
             *   This stage is essential in order to avoid bias in the result.
             * * otherwise, we return A*MaxCnt+N
             */
            ae_assert(n/maxcnt+1<=maxcnt, "HQRNDUniformI: N is too large", _state);
            result = -1;
            do
            {
                a = hqrnduniformi(state, maxcnt, _state);
                b = hqrnduniformi(state, n/maxcnt+1, _state);
                if( b==n/maxcnt&&a>=n%maxcnt )
                {
                    continue;
                }
                result = a+maxcnt*b;
            }
            while(result<0);
        }
    }
    else
    {
        
        /*
         * N<=MaxCnt
         *
         * Code below is a bit complicated because we can not simply
         * return "HQRNDIntegerBase() mod N" - it will be skewed for
         * large N's in [0.1*HQRNDMax...HQRNDMax].
         */
        mx = maxcnt-maxcnt%n;
        do
        {
            result = hqrnd_hqrndintegerbase(state, _state);
        }
        while(result>=mx);
        result = result%n;
    }
    return result;
}


/*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrndnormal(hqrndstate* state, ae_state *_state)
{
    double v1;
    double v2;
    double result;


    hqrndnormal2(state, &v1, &v2, _state);
    result = v1;
    return result;
}


/*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndunit2(hqrndstate* state, double* x, double* y, ae_state *_state)
{
    double v;
    double mx;
    double mn;

    *x = 0;
    *y = 0;

    do
    {
        hqrndnormal2(state, x, y, _state);
    }
    while(!(ae_fp_neq(*x,(double)(0))||ae_fp_neq(*y,(double)(0))));
    mx = ae_maxreal(ae_fabs(*x, _state), ae_fabs(*y, _state), _state);
    mn = ae_minreal(ae_fabs(*x, _state), ae_fabs(*y, _state), _state);
    v = mx*ae_sqrt(1+ae_sqr(mn/mx, _state), _state);
    *x = *x/v;
    *y = *y/v;
}


/*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndnormal2(hqrndstate* state,
     double* x1,
     double* x2,
     ae_state *_state)
{
    double u;
    double v;
    double s;

    *x1 = 0;
    *x2 = 0;

    for(;;)
    {
        u = 2*hqrnduniformr(state, _state)-1;
        v = 2*hqrnduniformr(state, _state)-1;
        s = ae_sqr(u, _state)+ae_sqr(v, _state);
        if( ae_fp_greater(s,(double)(0))&&ae_fp_less(s,(double)(1)) )
        {
            
            /*
             * two Sqrt's instead of one to
             * avoid overflow when S is too small
             */
            s = ae_sqrt(-2*ae_log(s, _state), _state)/ae_sqrt(s, _state);
            *x1 = u*s;
            *x2 = v*s;
            return;
        }
    }
}


/*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************/
double hqrndexponential(hqrndstate* state,
     double lambdav,
     ae_state *_state)
{
    double result;


    ae_assert(ae_fp_greater(lambdav,(double)(0)), "HQRNDExponential: LambdaV<=0!", _state);
    result = -ae_log(hqrnduniformr(state, _state), _state)/lambdav;
    return result;
}


/*************************************************************************
This function generates  random number from discrete distribution given by
finite sample X.

INPUT PARAMETERS
    State   -   high quality random number generator, must be
                initialized with HQRNDRandomize() or HQRNDSeed().
        X   -   finite sample
        N   -   number of elements to use, N>=1

RESULT
    this function returns one of the X[i] for random i=0..N-1

  -- ALGLIB --
     Copyright 08.11.2011 by Bochkanov Sergey
*************************************************************************/
double hqrnddiscrete(hqrndstate* state,
     /* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    double result;


    ae_assert(n>0, "HQRNDDiscrete: N<=0", _state);
    ae_assert(n<=x->cnt, "HQRNDDiscrete: Length(X)<N", _state);
    result = x->ptr.p_double[hqrnduniformi(state, n, _state)];
    return result;
}


/*************************************************************************
This function generates random number from continuous  distribution  given
by finite sample X.

INPUT PARAMETERS
    State   -   high quality random number generator, must be
                initialized with HQRNDRandomize() or HQRNDSeed().
        X   -   finite sample, array[N] (can be larger, in this  case only
                leading N elements are used). THIS ARRAY MUST BE SORTED BY
                ASCENDING.
        N   -   number of elements to use, N>=1

RESULT
    this function returns random number from continuous distribution which  
    tries to approximate X as mush as possible. min(X)<=Result<=max(X).

  -- ALGLIB --
     Copyright 08.11.2011 by Bochkanov Sergey
*************************************************************************/
double hqrndcontinuous(hqrndstate* state,
     /* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    double mx;
    double mn;
    ae_int_t i;
    double result;


    ae_assert(n>0, "HQRNDContinuous: N<=0", _state);
    ae_assert(n<=x->cnt, "HQRNDContinuous: Length(X)<N", _state);
    if( n==1 )
    {
        result = x->ptr.p_double[0];
        return result;
    }
    i = hqrnduniformi(state, n-1, _state);
    mn = x->ptr.p_double[i];
    mx = x->ptr.p_double[i+1];
    ae_assert(ae_fp_greater_eq(mx,mn), "HQRNDDiscrete: X is not sorted by ascending", _state);
    if( ae_fp_neq(mx,mn) )
    {
        result = (mx-mn)*hqrnduniformr(state, _state)+mn;
    }
    else
    {
        result = mn;
    }
    return result;
}


/*************************************************************************
This function returns random integer in [0,HQRNDMax]

L'Ecuyer, Efficient and portable combined random number generators
*************************************************************************/
static ae_int_t hqrnd_hqrndintegerbase(hqrndstate* state,
     ae_state *_state)
{
    ae_int_t k;
    ae_int_t result;


    ae_assert(state->magicv==hqrnd_hqrndmagic, "HQRNDIntegerBase: State is not correctly initialized!", _state);
    k = state->s1/53668;
    state->s1 = 40014*(state->s1-k*53668)-k*12211;
    if( state->s1<0 )
    {
        state->s1 = state->s1+2147483563;
    }
    k = state->s2/52774;
    state->s2 = 40692*(state->s2-k*52774)-k*3791;
    if( state->s2<0 )
    {
        state->s2 = state->s2+2147483399;
    }
    
    /*
     * Result
     */
    result = state->s1-state->s2;
    if( result<1 )
    {
        result = result+2147483562;
    }
    result = result-1;
    return result;
}


void _hqrndstate_init(void* _p, ae_state *_state)
{
    hqrndstate *p = (hqrndstate*)_p;
    ae_touch_ptr((void*)p);
}


void _hqrndstate_init_copy(void* _dst, void* _src, ae_state *_state)
{
    hqrndstate *dst = (hqrndstate*)_dst;
    hqrndstate *src = (hqrndstate*)_src;
    dst->s1 = src->s1;
    dst->s2 = src->s2;
    dst->magicv = src->magicv;
}


void _hqrndstate_clear(void* _p)
{
    hqrndstate *p = (hqrndstate*)_p;
    ae_touch_ptr((void*)p);
}


void _hqrndstate_destroy(void* _p)
{
    hqrndstate *p = (hqrndstate*)_p;
    ae_touch_ptr((void*)p);
}




/*************************************************************************
KD-tree creation

This subroutine creates KD-tree from set of X-values and optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    N       -   number of points, N>=0.
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)
                
OUTPUT PARAMETERS
    KDT     -   KD-tree
    
    
NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuild(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_int_t normtype,
     kdtree* kdt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tags;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _kdtree_clear(kdt);
    ae_vector_init(&tags, 0, DT_INT, _state);

    ae_assert(n>=0, "KDTreeBuild: N<0", _state);
    ae_assert(nx>=1, "KDTreeBuild: NX<1", _state);
    ae_assert(ny>=0, "KDTreeBuild: NY<0", _state);
    ae_assert(normtype>=0&&normtype<=2, "KDTreeBuild: incorrect NormType", _state);
    ae_assert(xy->rows>=n, "KDTreeBuild: rows(X)<N", _state);
    ae_assert(xy->cols>=nx+ny||n==0, "KDTreeBuild: cols(X)<NX+NY", _state);
    ae_assert(apservisfinitematrix(xy, n, nx+ny, _state), "KDTreeBuild: XY contains infinite or NaN values", _state);
    if( n>0 )
    {
        ae_vector_set_length(&tags, n, _state);
        for(i=0; i<=n-1; i++)
        {
            tags.ptr.p_int[i] = 0;
        }
    }
    kdtreebuildtagged(xy, &tags, n, nx, ny, normtype, kdt, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
KD-tree creation

This  subroutine  creates  KD-tree  from set of X-values, integer tags and
optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    Tags    -   tags, array[0..N-1], contains integer tags associated
                with points.
    N       -   number of points, N>=0
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree

NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreebuildtagged(/* Real    */ ae_matrix* xy,
     /* Integer */ ae_vector* tags,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_int_t normtype,
     kdtree* kdt,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t maxnodes;
    ae_int_t nodesoffs;
    ae_int_t splitsoffs;

    _kdtree_clear(kdt);

    ae_assert(n>=0, "KDTreeBuildTagged: N<0", _state);
    ae_assert(nx>=1, "KDTreeBuildTagged: NX<1", _state);
    ae_assert(ny>=0, "KDTreeBuildTagged: NY<0", _state);
    ae_assert(normtype>=0&&normtype<=2, "KDTreeBuildTagged: incorrect NormType", _state);
    ae_assert(xy->rows>=n, "KDTreeBuildTagged: rows(X)<N", _state);
    ae_assert(xy->cols>=nx+ny||n==0, "KDTreeBuildTagged: cols(X)<NX+NY", _state);
    ae_assert(apservisfinitematrix(xy, n, nx+ny, _state), "KDTreeBuildTagged: XY contains infinite or NaN values", _state);
    
    /*
     * initialize
     */
    kdt->n = n;
    kdt->nx = nx;
    kdt->ny = ny;
    kdt->normtype = normtype;
    kdt->kcur = 0;
    
    /*
     * N=0 => quick exit
     */
    if( n==0 )
    {
        return;
    }
    
    /*
     * Allocate
     */
    nearestneighbor_kdtreeallocdatasetindependent(kdt, nx, ny, _state);
    nearestneighbor_kdtreeallocdatasetdependent(kdt, n, nx, ny, _state);
    
    /*
     * Initial fill
     */
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&kdt->xy.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
        ae_v_move(&kdt->xy.ptr.pp_double[i][nx], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(nx,2*nx+ny-1));
        kdt->tags.ptr.p_int[i] = tags->ptr.p_int[i];
    }
    
    /*
     * Determine bounding box
     */
    ae_v_move(&kdt->boxmin.ptr.p_double[0], 1, &kdt->xy.ptr.pp_double[0][0], 1, ae_v_len(0,nx-1));
    ae_v_move(&kdt->boxmax.ptr.p_double[0], 1, &kdt->xy.ptr.pp_double[0][0], 1, ae_v_len(0,nx-1));
    for(i=1; i<=n-1; i++)
    {
        for(j=0; j<=nx-1; j++)
        {
            kdt->boxmin.ptr.p_double[j] = ae_minreal(kdt->boxmin.ptr.p_double[j], kdt->xy.ptr.pp_double[i][j], _state);
            kdt->boxmax.ptr.p_double[j] = ae_maxreal(kdt->boxmax.ptr.p_double[j], kdt->xy.ptr.pp_double[i][j], _state);
        }
    }
    
    /*
     * prepare tree structure
     * * MaxNodes=N because we guarantee no trivial splits, i.e.
     *   every split will generate two non-empty boxes
     */
    maxnodes = n;
    ae_vector_set_length(&kdt->nodes, nearestneighbor_splitnodesize*2*maxnodes, _state);
    ae_vector_set_length(&kdt->splits, 2*maxnodes, _state);
    nodesoffs = 0;
    splitsoffs = 0;
    ae_v_move(&kdt->curboxmin.ptr.p_double[0], 1, &kdt->boxmin.ptr.p_double[0], 1, ae_v_len(0,nx-1));
    ae_v_move(&kdt->curboxmax.ptr.p_double[0], 1, &kdt->boxmax.ptr.p_double[0], 1, ae_v_len(0,nx-1));
    nearestneighbor_kdtreegeneratetreerec(kdt, &nodesoffs, &splitsoffs, 0, n, 8, _state);
}


/*************************************************************************
K-NN query: K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of actual neighbors found (either K or N, if K>N).

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryknn(kdtree* kdt,
     /* Real    */ ae_vector* x,
     ae_int_t k,
     ae_bool selfmatch,
     ae_state *_state)
{
    ae_int_t result;


    ae_assert(k>=1, "KDTreeQueryKNN: K<1!", _state);
    ae_assert(x->cnt>=kdt->nx, "KDTreeQueryKNN: Length(X)<NX!", _state);
    ae_assert(isfinitevector(x, kdt->nx, _state), "KDTreeQueryKNN: X contains infinite or NaN values!", _state);
    result = kdtreequeryaknn(kdt, x, k, selfmatch, 0.0, _state);
    return result;
}


/*************************************************************************
R-NN query: all points within R-sphere centered at X

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    R           -   radius of sphere (in corresponding norm), R>0
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True

RESULT
    number of neighbors found, >=0

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
actual results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryrnn(kdtree* kdt,
     /* Real    */ ae_vector* x,
     double r,
     ae_bool selfmatch,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t result;


    ae_assert(ae_fp_greater(r,(double)(0)), "KDTreeQueryRNN: incorrect R!", _state);
    ae_assert(x->cnt>=kdt->nx, "KDTreeQueryRNN: Length(X)<NX!", _state);
    ae_assert(isfinitevector(x, kdt->nx, _state), "KDTreeQueryRNN: X contains infinite or NaN values!", _state);
    
    /*
     * Handle special case: KDT.N=0
     */
    if( kdt->n==0 )
    {
        kdt->kcur = 0;
        result = 0;
        return result;
    }
    
    /*
     * Prepare parameters
     */
    kdt->kneeded = 0;
    if( kdt->normtype!=2 )
    {
        kdt->rneeded = r;
    }
    else
    {
        kdt->rneeded = ae_sqr(r, _state);
    }
    kdt->selfmatch = selfmatch;
    kdt->approxf = (double)(1);
    kdt->kcur = 0;
    
    /*
     * calculate distance from point to current bounding box
     */
    nearestneighbor_kdtreeinitbox(kdt, x, _state);
    
    /*
     * call recursive search
     * results are returned as heap
     */
    nearestneighbor_kdtreequerynnrec(kdt, 0, _state);
    
    /*
     * pop from heap to generate ordered representation
     *
     * last element is not pop'ed because it is already in
     * its place
     */
    result = kdt->kcur;
    j = kdt->kcur;
    for(i=kdt->kcur; i>=2; i--)
    {
        tagheappopi(&kdt->r, &kdt->idx, &j, _state);
    }
    return result;
}


/*************************************************************************
K-NN query: approximate K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
                    * if not given, considered True
    Eps         -   approximation factor, Eps>=0. eps-approximate  nearest
                    neighbor  is  a  neighbor  whose distance from X is at
                    most (1+eps) times distance of true nearest neighbor.

RESULT
    number of actual neighbors found (either K or N, if K>N).
    
NOTES
    significant performance gain may be achieved only when Eps  is  is  on
    the order of magnitude of 1 or larger.

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
ae_int_t kdtreequeryaknn(kdtree* kdt,
     /* Real    */ ae_vector* x,
     ae_int_t k,
     ae_bool selfmatch,
     double eps,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t result;


    ae_assert(k>0, "KDTreeQueryAKNN: incorrect K!", _state);
    ae_assert(ae_fp_greater_eq(eps,(double)(0)), "KDTreeQueryAKNN: incorrect Eps!", _state);
    ae_assert(x->cnt>=kdt->nx, "KDTreeQueryAKNN: Length(X)<NX!", _state);
    ae_assert(isfinitevector(x, kdt->nx, _state), "KDTreeQueryAKNN: X contains infinite or NaN values!", _state);
    
    /*
     * Handle special case: KDT.N=0
     */
    if( kdt->n==0 )
    {
        kdt->kcur = 0;
        result = 0;
        return result;
    }
    
    /*
     * Prepare parameters
     */
    k = ae_minint(k, kdt->n, _state);
    kdt->kneeded = k;
    kdt->rneeded = (double)(0);
    kdt->selfmatch = selfmatch;
    if( kdt->normtype==2 )
    {
        kdt->approxf = 1/ae_sqr(1+eps, _state);
    }
    else
    {
        kdt->approxf = 1/(1+eps);
    }
    kdt->kcur = 0;
    
    /*
     * calculate distance from point to current bounding box
     */
    nearestneighbor_kdtreeinitbox(kdt, x, _state);
    
    /*
     * call recursive search
     * results are returned as heap
     */
    nearestneighbor_kdtreequerynnrec(kdt, 0, _state);
    
    /*
     * pop from heap to generate ordered representation
     *
     * last element is non pop'ed because it is already in
     * its place
     */
    result = kdt->kcur;
    j = kdt->kcur;
    for(i=kdt->kcur; i>=2; i--)
    {
        tagheappopi(&kdt->r, &kdt->idx, &j, _state);
    }
    return result;
}


/*************************************************************************
X-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    X       -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    X       -   rows are filled with X-values

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsx(kdtree* kdt,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;


    if( kdt->kcur==0 )
    {
        return;
    }
    if( x->rows<kdt->kcur||x->cols<kdt->nx )
    {
        ae_matrix_set_length(x, kdt->kcur, kdt->nx, _state);
    }
    k = kdt->kcur;
    for(i=0; i<=k-1; i++)
    {
        ae_v_move(&x->ptr.pp_double[i][0], 1, &kdt->xy.ptr.pp_double[kdt->idx.ptr.p_int[i]][kdt->nx], 1, ae_v_len(0,kdt->nx-1));
    }
}


/*************************************************************************
X- and Y-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    XY      -   possibly pre-allocated buffer. If XY is too small to store
                result, it is resized. If size(XY) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    XY      -   rows are filled with points: first NX columns with
                X-values, next NY columns - with Y-values.

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxy(kdtree* kdt,
     /* Real    */ ae_matrix* xy,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;


    if( kdt->kcur==0 )
    {
        return;
    }
    if( xy->rows<kdt->kcur||xy->cols<kdt->nx+kdt->ny )
    {
        ae_matrix_set_length(xy, kdt->kcur, kdt->nx+kdt->ny, _state);
    }
    k = kdt->kcur;
    for(i=0; i<=k-1; i++)
    {
        ae_v_move(&xy->ptr.pp_double[i][0], 1, &kdt->xy.ptr.pp_double[kdt->idx.ptr.p_int[i]][kdt->nx], 1, ae_v_len(0,kdt->nx+kdt->ny-1));
    }
}


/*************************************************************************
Tags from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    Tags    -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    Tags    -   filled with tags associated with points,
                or, when no tags were supplied, with zeros

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultstags(kdtree* kdt,
     /* Integer */ ae_vector* tags,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;


    if( kdt->kcur==0 )
    {
        return;
    }
    if( tags->cnt<kdt->kcur )
    {
        ae_vector_set_length(tags, kdt->kcur, _state);
    }
    k = kdt->kcur;
    for(i=0; i<=k-1; i++)
    {
        tags->ptr.p_int[i] = kdt->tags.ptr.p_int[kdt->idx.ptr.p_int[i]];
    }
}


/*************************************************************************
Distances from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    R       -   possibly pre-allocated buffer. If X is too small to store
                result, it is resized. If size(X) is enough to store
                result, it is left unchanged.

OUTPUT PARAMETERS
    R       -   filled with distances (in corresponding norm)

NOTES
1. points are ordered by distance from the query point (first = closest)
2. if  XY is larger than required to store result, only leading part  will
   be overwritten; trailing part will be left unchanged. So  if  on  input
   XY = [[A,B],[C,D]], and result is [1,2],  then  on  exit  we  will  get
   XY = [[1,2],[C,D]]. This is done purposely to increase performance;  if
   you want function  to  resize  array  according  to  result  size,  use
   function with same name and suffix 'I'.

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsdistances(kdtree* kdt,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;


    if( kdt->kcur==0 )
    {
        return;
    }
    if( r->cnt<kdt->kcur )
    {
        ae_vector_set_length(r, kdt->kcur, _state);
    }
    k = kdt->kcur;
    
    /*
     * unload norms
     *
     * Abs() call is used to handle cases with negative norms
     * (generated during KFN requests)
     */
    if( kdt->normtype==0 )
    {
        for(i=0; i<=k-1; i++)
        {
            r->ptr.p_double[i] = ae_fabs(kdt->r.ptr.p_double[i], _state);
        }
    }
    if( kdt->normtype==1 )
    {
        for(i=0; i<=k-1; i++)
        {
            r->ptr.p_double[i] = ae_fabs(kdt->r.ptr.p_double[i], _state);
        }
    }
    if( kdt->normtype==2 )
    {
        for(i=0; i<=k-1; i++)
        {
            r->ptr.p_double[i] = ae_sqrt(ae_fabs(kdt->r.ptr.p_double[i], _state), _state);
        }
    }
}


/*************************************************************************
X-values from last query; 'interactive' variant for languages like  Python
which   support    constructs   like  "X = KDTreeQueryResultsXI(KDT)"  and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxi(kdtree* kdt,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{

    ae_matrix_clear(x);

    kdtreequeryresultsx(kdt, x, _state);
}


/*************************************************************************
XY-values from last query; 'interactive' variant for languages like Python
which   support    constructs   like "XY = KDTreeQueryResultsXYI(KDT)" and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsxyi(kdtree* kdt,
     /* Real    */ ae_matrix* xy,
     ae_state *_state)
{

    ae_matrix_clear(xy);

    kdtreequeryresultsxy(kdt, xy, _state);
}


/*************************************************************************
Tags  from  last  query;  'interactive' variant for languages like  Python
which  support  constructs  like "Tags = KDTreeQueryResultsTagsI(KDT)" and
interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultstagsi(kdtree* kdt,
     /* Integer */ ae_vector* tags,
     ae_state *_state)
{

    ae_vector_clear(tags);

    kdtreequeryresultstags(kdt, tags, _state);
}


/*************************************************************************
Distances from last query; 'interactive' variant for languages like Python
which  support  constructs   like  "R = KDTreeQueryResultsDistancesI(KDT)"
and interactive mode of interpreter.

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void kdtreequeryresultsdistancesi(kdtree* kdt,
     /* Real    */ ae_vector* r,
     ae_state *_state)
{

    ae_vector_clear(r);

    kdtreequeryresultsdistances(kdt, r, _state);
}


/*************************************************************************
Serializer: allocation

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void kdtreealloc(ae_serializer* s, kdtree* tree, ae_state *_state)
{


    
    /*
     * Header
     */
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    
    /*
     * Data
     */
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    allocrealmatrix(s, &tree->xy, -1, -1, _state);
    allocintegerarray(s, &tree->tags, -1, _state);
    allocrealarray(s, &tree->boxmin, -1, _state);
    allocrealarray(s, &tree->boxmax, -1, _state);
    allocintegerarray(s, &tree->nodes, -1, _state);
    allocrealarray(s, &tree->splits, -1, _state);
}


/*************************************************************************
Serializer: serialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void kdtreeserialize(ae_serializer* s, kdtree* tree, ae_state *_state)
{


    
    /*
     * Header
     */
    ae_serializer_serialize_int(s, getkdtreeserializationcode(_state), _state);
    ae_serializer_serialize_int(s, nearestneighbor_kdtreefirstversion, _state);
    
    /*
     * Data
     */
    ae_serializer_serialize_int(s, tree->n, _state);
    ae_serializer_serialize_int(s, tree->nx, _state);
    ae_serializer_serialize_int(s, tree->ny, _state);
    ae_serializer_serialize_int(s, tree->normtype, _state);
    serializerealmatrix(s, &tree->xy, -1, -1, _state);
    serializeintegerarray(s, &tree->tags, -1, _state);
    serializerealarray(s, &tree->boxmin, -1, _state);
    serializerealarray(s, &tree->boxmax, -1, _state);
    serializeintegerarray(s, &tree->nodes, -1, _state);
    serializerealarray(s, &tree->splits, -1, _state);
}


/*************************************************************************
Serializer: unserialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void kdtreeunserialize(ae_serializer* s, kdtree* tree, ae_state *_state)
{
    ae_int_t i0;
    ae_int_t i1;

    _kdtree_clear(tree);

    
    /*
     * check correctness of header
     */
    ae_serializer_unserialize_int(s, &i0, _state);
    ae_assert(i0==getkdtreeserializationcode(_state), "KDTreeUnserialize: stream header corrupted", _state);
    ae_serializer_unserialize_int(s, &i1, _state);
    ae_assert(i1==nearestneighbor_kdtreefirstversion, "KDTreeUnserialize: stream header corrupted", _state);
    
    /*
     * Unserialize data
     */
    ae_serializer_unserialize_int(s, &tree->n, _state);
    ae_serializer_unserialize_int(s, &tree->nx, _state);
    ae_serializer_unserialize_int(s, &tree->ny, _state);
    ae_serializer_unserialize_int(s, &tree->normtype, _state);
    unserializerealmatrix(s, &tree->xy, _state);
    unserializeintegerarray(s, &tree->tags, _state);
    unserializerealarray(s, &tree->boxmin, _state);
    unserializerealarray(s, &tree->boxmax, _state);
    unserializeintegerarray(s, &tree->nodes, _state);
    unserializerealarray(s, &tree->splits, _state);
    nearestneighbor_kdtreealloctemporaries(tree, tree->n, tree->nx, tree->ny, _state);
}


/*************************************************************************
Rearranges nodes [I1,I2) using partition in D-th dimension with S as threshold.
Returns split position I3: [I1,I3) and [I3,I2) are created as result.

This subroutine doesn't create tree structures, just rearranges nodes.
*************************************************************************/
static void nearestneighbor_kdtreesplit(kdtree* kdt,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t d,
     double s,
     ae_int_t* i3,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t ileft;
    ae_int_t iright;
    double v;

    *i3 = 0;

    ae_assert(kdt->n>0, "KDTreeSplit: internal error", _state);
    
    /*
     * split XY/Tags in two parts:
     * * [ILeft,IRight] is non-processed part of XY/Tags
     *
     * After cycle is done, we have Ileft=IRight. We deal with
     * this element separately.
     *
     * After this, [I1,ILeft) contains left part, and [ILeft,I2)
     * contains right part.
     */
    ileft = i1;
    iright = i2-1;
    while(ileft<iright)
    {
        if( ae_fp_less_eq(kdt->xy.ptr.pp_double[ileft][d],s) )
        {
            
            /*
             * XY[ILeft] is on its place.
             * Advance ILeft.
             */
            ileft = ileft+1;
        }
        else
        {
            
            /*
             * XY[ILeft,..] must be at IRight.
             * Swap and advance IRight.
             */
            for(i=0; i<=2*kdt->nx+kdt->ny-1; i++)
            {
                v = kdt->xy.ptr.pp_double[ileft][i];
                kdt->xy.ptr.pp_double[ileft][i] = kdt->xy.ptr.pp_double[iright][i];
                kdt->xy.ptr.pp_double[iright][i] = v;
            }
            j = kdt->tags.ptr.p_int[ileft];
            kdt->tags.ptr.p_int[ileft] = kdt->tags.ptr.p_int[iright];
            kdt->tags.ptr.p_int[iright] = j;
            iright = iright-1;
        }
    }
    if( ae_fp_less_eq(kdt->xy.ptr.pp_double[ileft][d],s) )
    {
        ileft = ileft+1;
    }
    else
    {
        iright = iright-1;
    }
    *i3 = ileft;
}


/*************************************************************************
Recursive kd-tree generation subroutine.

PARAMETERS
    KDT         tree
    NodesOffs   unused part of Nodes[] which must be filled by tree
    SplitsOffs  unused part of Splits[]
    I1, I2      points from [I1,I2) are processed
    
NodesOffs[] and SplitsOffs[] must be large enough.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreegeneratetreerec(kdtree* kdt,
     ae_int_t* nodesoffs,
     ae_int_t* splitsoffs,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t maxleafsize,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t nx;
    ae_int_t ny;
    ae_int_t i;
    ae_int_t j;
    ae_int_t oldoffs;
    ae_int_t i3;
    ae_int_t cntless;
    ae_int_t cntgreater;
    double minv;
    double maxv;
    ae_int_t minidx;
    ae_int_t maxidx;
    ae_int_t d;
    double ds;
    double s;
    double v;
    double v0;
    double v1;


    ae_assert(kdt->n>0, "KDTreeGenerateTreeRec: internal error", _state);
    ae_assert(i2>i1, "KDTreeGenerateTreeRec: internal error", _state);
    
    /*
     * Generate leaf if needed
     */
    if( i2-i1<=maxleafsize )
    {
        kdt->nodes.ptr.p_int[*nodesoffs+0] = i2-i1;
        kdt->nodes.ptr.p_int[*nodesoffs+1] = i1;
        *nodesoffs = *nodesoffs+2;
        return;
    }
    
    /*
     * Load values for easier access
     */
    nx = kdt->nx;
    ny = kdt->ny;
    
    /*
     * Select dimension to split:
     * * D is a dimension number
     * In case bounding box has zero size, we enforce creation of the leaf node.
     */
    d = 0;
    ds = kdt->curboxmax.ptr.p_double[0]-kdt->curboxmin.ptr.p_double[0];
    for(i=1; i<=nx-1; i++)
    {
        v = kdt->curboxmax.ptr.p_double[i]-kdt->curboxmin.ptr.p_double[i];
        if( ae_fp_greater(v,ds) )
        {
            ds = v;
            d = i;
        }
    }
    if( ae_fp_eq(ds,(double)(0)) )
    {
        kdt->nodes.ptr.p_int[*nodesoffs+0] = i2-i1;
        kdt->nodes.ptr.p_int[*nodesoffs+1] = i1;
        *nodesoffs = *nodesoffs+2;
        return;
    }
    
    /*
     * Select split position S using sliding midpoint rule,
     * rearrange points into [I1,I3) and [I3,I2).
     *
     * In case all points has same value of D-th component
     * (MinV=MaxV) we enforce D-th dimension of bounding
     * box to become exactly zero and repeat tree construction.
     */
    s = kdt->curboxmin.ptr.p_double[d]+0.5*ds;
    ae_v_move(&kdt->buf.ptr.p_double[0], 1, &kdt->xy.ptr.pp_double[i1][d], kdt->xy.stride, ae_v_len(0,i2-i1-1));
    n = i2-i1;
    cntless = 0;
    cntgreater = 0;
    minv = kdt->buf.ptr.p_double[0];
    maxv = kdt->buf.ptr.p_double[0];
    minidx = i1;
    maxidx = i1;
    for(i=0; i<=n-1; i++)
    {
        v = kdt->buf.ptr.p_double[i];
        if( ae_fp_less(v,minv) )
        {
            minv = v;
            minidx = i1+i;
        }
        if( ae_fp_greater(v,maxv) )
        {
            maxv = v;
            maxidx = i1+i;
        }
        if( ae_fp_less(v,s) )
        {
            cntless = cntless+1;
        }
        if( ae_fp_greater(v,s) )
        {
            cntgreater = cntgreater+1;
        }
    }
    if( ae_fp_eq(minv,maxv) )
    {
        
        /*
         * In case all points has same value of D-th component
         * (MinV=MaxV) we enforce D-th dimension of bounding
         * box to become exactly zero and repeat tree construction.
         */
        v0 = kdt->curboxmin.ptr.p_double[d];
        v1 = kdt->curboxmax.ptr.p_double[d];
        kdt->curboxmin.ptr.p_double[d] = minv;
        kdt->curboxmax.ptr.p_double[d] = maxv;
        nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i1, i2, maxleafsize, _state);
        kdt->curboxmin.ptr.p_double[d] = v0;
        kdt->curboxmax.ptr.p_double[d] = v1;
        return;
    }
    if( cntless>0&&cntgreater>0 )
    {
        
        /*
         * normal midpoint split
         */
        nearestneighbor_kdtreesplit(kdt, i1, i2, d, s, &i3, _state);
    }
    else
    {
        
        /*
         * sliding midpoint
         */
        if( cntless==0 )
        {
            
            /*
             * 1. move split to MinV,
             * 2. place one point to the left bin (move to I1),
             *    others - to the right bin
             */
            s = minv;
            if( minidx!=i1 )
            {
                for(i=0; i<=2*nx+ny-1; i++)
                {
                    v = kdt->xy.ptr.pp_double[minidx][i];
                    kdt->xy.ptr.pp_double[minidx][i] = kdt->xy.ptr.pp_double[i1][i];
                    kdt->xy.ptr.pp_double[i1][i] = v;
                }
                j = kdt->tags.ptr.p_int[minidx];
                kdt->tags.ptr.p_int[minidx] = kdt->tags.ptr.p_int[i1];
                kdt->tags.ptr.p_int[i1] = j;
            }
            i3 = i1+1;
        }
        else
        {
            
            /*
             * 1. move split to MaxV,
             * 2. place one point to the right bin (move to I2-1),
             *    others - to the left bin
             */
            s = maxv;
            if( maxidx!=i2-1 )
            {
                for(i=0; i<=2*nx+ny-1; i++)
                {
                    v = kdt->xy.ptr.pp_double[maxidx][i];
                    kdt->xy.ptr.pp_double[maxidx][i] = kdt->xy.ptr.pp_double[i2-1][i];
                    kdt->xy.ptr.pp_double[i2-1][i] = v;
                }
                j = kdt->tags.ptr.p_int[maxidx];
                kdt->tags.ptr.p_int[maxidx] = kdt->tags.ptr.p_int[i2-1];
                kdt->tags.ptr.p_int[i2-1] = j;
            }
            i3 = i2-1;
        }
    }
    
    /*
     * Generate 'split' node
     */
    kdt->nodes.ptr.p_int[*nodesoffs+0] = 0;
    kdt->nodes.ptr.p_int[*nodesoffs+1] = d;
    kdt->nodes.ptr.p_int[*nodesoffs+2] = *splitsoffs;
    kdt->splits.ptr.p_double[*splitsoffs+0] = s;
    oldoffs = *nodesoffs;
    *nodesoffs = *nodesoffs+nearestneighbor_splitnodesize;
    *splitsoffs = *splitsoffs+1;
    
    /*
     * Recirsive generation:
     * * update CurBox
     * * call subroutine
     * * restore CurBox
     */
    kdt->nodes.ptr.p_int[oldoffs+3] = *nodesoffs;
    v = kdt->curboxmax.ptr.p_double[d];
    kdt->curboxmax.ptr.p_double[d] = s;
    nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i1, i3, maxleafsize, _state);
    kdt->curboxmax.ptr.p_double[d] = v;
    kdt->nodes.ptr.p_int[oldoffs+4] = *nodesoffs;
    v = kdt->curboxmin.ptr.p_double[d];
    kdt->curboxmin.ptr.p_double[d] = s;
    nearestneighbor_kdtreegeneratetreerec(kdt, nodesoffs, splitsoffs, i3, i2, maxleafsize, _state);
    kdt->curboxmin.ptr.p_double[d] = v;
}


/*************************************************************************
Recursive subroutine for NN queries.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreequerynnrec(kdtree* kdt,
     ae_int_t offs,
     ae_state *_state)
{
    double ptdist;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nx;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t d;
    double s;
    double v;
    double t1;
    ae_int_t childbestoffs;
    ae_int_t childworstoffs;
    ae_int_t childoffs;
    double prevdist;
    ae_bool todive;
    ae_bool bestisleft;
    ae_bool updatemin;


    ae_assert(kdt->n>0, "KDTreeQueryNNRec: internal error", _state);
    
    /*
     * Leaf node.
     * Process points.
     */
    if( kdt->nodes.ptr.p_int[offs]>0 )
    {
        i1 = kdt->nodes.ptr.p_int[offs+1];
        i2 = i1+kdt->nodes.ptr.p_int[offs];
        for(i=i1; i<=i2-1; i++)
        {
            
            /*
             * Calculate distance
             */
            ptdist = (double)(0);
            nx = kdt->nx;
            if( kdt->normtype==0 )
            {
                for(j=0; j<=nx-1; j++)
                {
                    ptdist = ae_maxreal(ptdist, ae_fabs(kdt->xy.ptr.pp_double[i][j]-kdt->x.ptr.p_double[j], _state), _state);
                }
            }
            if( kdt->normtype==1 )
            {
                for(j=0; j<=nx-1; j++)
                {
                    ptdist = ptdist+ae_fabs(kdt->xy.ptr.pp_double[i][j]-kdt->x.ptr.p_double[j], _state);
                }
            }
            if( kdt->normtype==2 )
            {
                for(j=0; j<=nx-1; j++)
                {
                    ptdist = ptdist+ae_sqr(kdt->xy.ptr.pp_double[i][j]-kdt->x.ptr.p_double[j], _state);
                }
            }
            
            /*
             * Skip points with zero distance if self-matches are turned off
             */
            if( ae_fp_eq(ptdist,(double)(0))&&!kdt->selfmatch )
            {
                continue;
            }
            
            /*
             * We CAN'T process point if R-criterion isn't satisfied,
             * i.e. (RNeeded<>0) AND (PtDist>R).
             */
            if( ae_fp_eq(kdt->rneeded,(double)(0))||ae_fp_less_eq(ptdist,kdt->rneeded) )
            {
                
                /*
                 * R-criterion is satisfied, we must either:
                 * * replace worst point, if (KNeeded<>0) AND (KCur=KNeeded)
                 *   (or skip, if worst point is better)
                 * * add point without replacement otherwise
                 */
                if( kdt->kcur<kdt->kneeded||kdt->kneeded==0 )
                {
                    
                    /*
                     * add current point to heap without replacement
                     */
                    tagheappushi(&kdt->r, &kdt->idx, &kdt->kcur, ptdist, i, _state);
                }
                else
                {
                    
                    /*
                     * New points are added or not, depending on their distance.
                     * If added, they replace element at the top of the heap
                     */
                    if( ae_fp_less(ptdist,kdt->r.ptr.p_double[0]) )
                    {
                        if( kdt->kneeded==1 )
                        {
                            kdt->idx.ptr.p_int[0] = i;
                            kdt->r.ptr.p_double[0] = ptdist;
                        }
                        else
                        {
                            tagheapreplacetopi(&kdt->r, &kdt->idx, kdt->kneeded, ptdist, i, _state);
                        }
                    }
                }
            }
        }
        return;
    }
    
    /*
     * Simple split
     */
    if( kdt->nodes.ptr.p_int[offs]==0 )
    {
        
        /*
         * Load:
         * * D  dimension to split
         * * S  split position
         */
        d = kdt->nodes.ptr.p_int[offs+1];
        s = kdt->splits.ptr.p_double[kdt->nodes.ptr.p_int[offs+2]];
        
        /*
         * Calculate:
         * * ChildBestOffs      child box with best chances
         * * ChildWorstOffs     child box with worst chances
         */
        if( ae_fp_less_eq(kdt->x.ptr.p_double[d],s) )
        {
            childbestoffs = kdt->nodes.ptr.p_int[offs+3];
            childworstoffs = kdt->nodes.ptr.p_int[offs+4];
            bestisleft = ae_true;
        }
        else
        {
            childbestoffs = kdt->nodes.ptr.p_int[offs+4];
            childworstoffs = kdt->nodes.ptr.p_int[offs+3];
            bestisleft = ae_false;
        }
        
        /*
         * Navigate through childs
         */
        for(i=0; i<=1; i++)
        {
            
            /*
             * Select child to process:
             * * ChildOffs      current child offset in Nodes[]
             * * UpdateMin      whether minimum or maximum value
             *                  of bounding box is changed on update
             */
            if( i==0 )
            {
                childoffs = childbestoffs;
                updatemin = !bestisleft;
            }
            else
            {
                updatemin = bestisleft;
                childoffs = childworstoffs;
            }
            
            /*
             * Update bounding box and current distance
             */
            if( updatemin )
            {
                prevdist = kdt->curdist;
                t1 = kdt->x.ptr.p_double[d];
                v = kdt->curboxmin.ptr.p_double[d];
                if( ae_fp_less_eq(t1,s) )
                {
                    if( kdt->normtype==0 )
                    {
                        kdt->curdist = ae_maxreal(kdt->curdist, s-t1, _state);
                    }
                    if( kdt->normtype==1 )
                    {
                        kdt->curdist = kdt->curdist-ae_maxreal(v-t1, (double)(0), _state)+s-t1;
                    }
                    if( kdt->normtype==2 )
                    {
                        kdt->curdist = kdt->curdist-ae_sqr(ae_maxreal(v-t1, (double)(0), _state), _state)+ae_sqr(s-t1, _state);
                    }
                }
                kdt->curboxmin.ptr.p_double[d] = s;
            }
            else
            {
                prevdist = kdt->curdist;
                t1 = kdt->x.ptr.p_double[d];
                v = kdt->curboxmax.ptr.p_double[d];
                if( ae_fp_greater_eq(t1,s) )
                {
                    if( kdt->normtype==0 )
                    {
                        kdt->curdist = ae_maxreal(kdt->curdist, t1-s, _state);
                    }
                    if( kdt->normtype==1 )
                    {
                        kdt->curdist = kdt->curdist-ae_maxreal(t1-v, (double)(0), _state)+t1-s;
                    }
                    if( kdt->normtype==2 )
                    {
                        kdt->curdist = kdt->curdist-ae_sqr(ae_maxreal(t1-v, (double)(0), _state), _state)+ae_sqr(t1-s, _state);
                    }
                }
                kdt->curboxmax.ptr.p_double[d] = s;
            }
            
            /*
             * Decide: to dive into cell or not to dive
             */
            if( ae_fp_neq(kdt->rneeded,(double)(0))&&ae_fp_greater(kdt->curdist,kdt->rneeded) )
            {
                todive = ae_false;
            }
            else
            {
                if( kdt->kcur<kdt->kneeded||kdt->kneeded==0 )
                {
                    
                    /*
                     * KCur<KNeeded (i.e. not all points are found)
                     */
                    todive = ae_true;
                }
                else
                {
                    
                    /*
                     * KCur=KNeeded, decide to dive or not to dive
                     * using point position relative to bounding box.
                     */
                    todive = ae_fp_less_eq(kdt->curdist,kdt->r.ptr.p_double[0]*kdt->approxf);
                }
            }
            if( todive )
            {
                nearestneighbor_kdtreequerynnrec(kdt, childoffs, _state);
            }
            
            /*
             * Restore bounding box and distance
             */
            if( updatemin )
            {
                kdt->curboxmin.ptr.p_double[d] = v;
            }
            else
            {
                kdt->curboxmax.ptr.p_double[d] = v;
            }
            kdt->curdist = prevdist;
        }
        return;
    }
}


/*************************************************************************
Copies X[] to KDT.X[]
Loads distance from X[] to bounding box.
Initializes CurBox[].

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreeinitbox(kdtree* kdt,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t i;
    double vx;
    double vmin;
    double vmax;


    ae_assert(kdt->n>0, "KDTreeInitBox: internal error", _state);
    
    /*
     * calculate distance from point to current bounding box
     */
    kdt->curdist = (double)(0);
    if( kdt->normtype==0 )
    {
        for(i=0; i<=kdt->nx-1; i++)
        {
            vx = x->ptr.p_double[i];
            vmin = kdt->boxmin.ptr.p_double[i];
            vmax = kdt->boxmax.ptr.p_double[i];
            kdt->x.ptr.p_double[i] = vx;
            kdt->curboxmin.ptr.p_double[i] = vmin;
            kdt->curboxmax.ptr.p_double[i] = vmax;
            if( ae_fp_less(vx,vmin) )
            {
                kdt->curdist = ae_maxreal(kdt->curdist, vmin-vx, _state);
            }
            else
            {
                if( ae_fp_greater(vx,vmax) )
                {
                    kdt->curdist = ae_maxreal(kdt->curdist, vx-vmax, _state);
                }
            }
        }
    }
    if( kdt->normtype==1 )
    {
        for(i=0; i<=kdt->nx-1; i++)
        {
            vx = x->ptr.p_double[i];
            vmin = kdt->boxmin.ptr.p_double[i];
            vmax = kdt->boxmax.ptr.p_double[i];
            kdt->x.ptr.p_double[i] = vx;
            kdt->curboxmin.ptr.p_double[i] = vmin;
            kdt->curboxmax.ptr.p_double[i] = vmax;
            if( ae_fp_less(vx,vmin) )
            {
                kdt->curdist = kdt->curdist+vmin-vx;
            }
            else
            {
                if( ae_fp_greater(vx,vmax) )
                {
                    kdt->curdist = kdt->curdist+vx-vmax;
                }
            }
        }
    }
    if( kdt->normtype==2 )
    {
        for(i=0; i<=kdt->nx-1; i++)
        {
            vx = x->ptr.p_double[i];
            vmin = kdt->boxmin.ptr.p_double[i];
            vmax = kdt->boxmax.ptr.p_double[i];
            kdt->x.ptr.p_double[i] = vx;
            kdt->curboxmin.ptr.p_double[i] = vmin;
            kdt->curboxmax.ptr.p_double[i] = vmax;
            if( ae_fp_less(vx,vmin) )
            {
                kdt->curdist = kdt->curdist+ae_sqr(vmin-vx, _state);
            }
            else
            {
                if( ae_fp_greater(vx,vmax) )
                {
                    kdt->curdist = kdt->curdist+ae_sqr(vx-vmax, _state);
                }
            }
        }
    }
}


/*************************************************************************
This function allocates all dataset-independent array  fields  of  KDTree,
i.e.  such  array  fields  that  their dimensions do not depend on dataset
size.

This function do not sets KDT.NX or KDT.NY - it just allocates arrays

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreeallocdatasetindependent(kdtree* kdt,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state)
{


    ae_assert(kdt->n>0, "KDTreeAllocDatasetIndependent: internal error", _state);
    ae_vector_set_length(&kdt->x, nx, _state);
    ae_vector_set_length(&kdt->boxmin, nx, _state);
    ae_vector_set_length(&kdt->boxmax, nx, _state);
    ae_vector_set_length(&kdt->curboxmin, nx, _state);
    ae_vector_set_length(&kdt->curboxmax, nx, _state);
}


/*************************************************************************
This function allocates all dataset-dependent array fields of KDTree, i.e.
such array fields that their dimensions depend on dataset size.

This function do not sets KDT.N, KDT.NX or KDT.NY -
it just allocates arrays.

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreeallocdatasetdependent(kdtree* kdt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state)
{


    ae_assert(n>0, "KDTreeAllocDatasetDependent: internal error", _state);
    ae_matrix_set_length(&kdt->xy, n, 2*nx+ny, _state);
    ae_vector_set_length(&kdt->tags, n, _state);
    ae_vector_set_length(&kdt->idx, n, _state);
    ae_vector_set_length(&kdt->r, n, _state);
    ae_vector_set_length(&kdt->x, nx, _state);
    ae_vector_set_length(&kdt->buf, ae_maxint(n, nx, _state), _state);
    ae_vector_set_length(&kdt->nodes, nearestneighbor_splitnodesize*2*n, _state);
    ae_vector_set_length(&kdt->splits, 2*n, _state);
}


/*************************************************************************
This function allocates temporaries.

This function do not sets KDT.N, KDT.NX or KDT.NY -
it just allocates arrays.

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
static void nearestneighbor_kdtreealloctemporaries(kdtree* kdt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state)
{


    ae_assert(n>0, "KDTreeAllocTemporaries: internal error", _state);
    ae_vector_set_length(&kdt->x, nx, _state);
    ae_vector_set_length(&kdt->idx, n, _state);
    ae_vector_set_length(&kdt->r, n, _state);
    ae_vector_set_length(&kdt->buf, ae_maxint(n, nx, _state), _state);
    ae_vector_set_length(&kdt->curboxmin, nx, _state);
    ae_vector_set_length(&kdt->curboxmax, nx, _state);
}


void _kdtree_init(void* _p, ae_state *_state)
{
    kdtree *p = (kdtree*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_init(&p->xy, 0, 0, DT_REAL, _state);
    ae_vector_init(&p->tags, 0, DT_INT, _state);
    ae_vector_init(&p->boxmin, 0, DT_REAL, _state);
    ae_vector_init(&p->boxmax, 0, DT_REAL, _state);
    ae_vector_init(&p->nodes, 0, DT_INT, _state);
    ae_vector_init(&p->splits, 0, DT_REAL, _state);
    ae_vector_init(&p->x, 0, DT_REAL, _state);
    ae_vector_init(&p->idx, 0, DT_INT, _state);
    ae_vector_init(&p->r, 0, DT_REAL, _state);
    ae_vector_init(&p->buf, 0, DT_REAL, _state);
    ae_vector_init(&p->curboxmin, 0, DT_REAL, _state);
    ae_vector_init(&p->curboxmax, 0, DT_REAL, _state);
}


void _kdtree_init_copy(void* _dst, void* _src, ae_state *_state)
{
    kdtree *dst = (kdtree*)_dst;
    kdtree *src = (kdtree*)_src;
    dst->n = src->n;
    dst->nx = src->nx;
    dst->ny = src->ny;
    dst->normtype = src->normtype;
    ae_matrix_init_copy(&dst->xy, &src->xy, _state);
    ae_vector_init_copy(&dst->tags, &src->tags, _state);
    ae_vector_init_copy(&dst->boxmin, &src->boxmin, _state);
    ae_vector_init_copy(&dst->boxmax, &src->boxmax, _state);
    ae_vector_init_copy(&dst->nodes, &src->nodes, _state);
    ae_vector_init_copy(&dst->splits, &src->splits, _state);
    ae_vector_init_copy(&dst->x, &src->x, _state);
    dst->kneeded = src->kneeded;
    dst->rneeded = src->rneeded;
    dst->selfmatch = src->selfmatch;
    dst->approxf = src->approxf;
    dst->kcur = src->kcur;
    ae_vector_init_copy(&dst->idx, &src->idx, _state);
    ae_vector_init_copy(&dst->r, &src->r, _state);
    ae_vector_init_copy(&dst->buf, &src->buf, _state);
    ae_vector_init_copy(&dst->curboxmin, &src->curboxmin, _state);
    ae_vector_init_copy(&dst->curboxmax, &src->curboxmax, _state);
    dst->curdist = src->curdist;
    dst->debugcounter = src->debugcounter;
}


void _kdtree_clear(void* _p)
{
    kdtree *p = (kdtree*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_clear(&p->xy);
    ae_vector_clear(&p->tags);
    ae_vector_clear(&p->boxmin);
    ae_vector_clear(&p->boxmax);
    ae_vector_clear(&p->nodes);
    ae_vector_clear(&p->splits);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->idx);
    ae_vector_clear(&p->r);
    ae_vector_clear(&p->buf);
    ae_vector_clear(&p->curboxmin);
    ae_vector_clear(&p->curboxmax);
}


void _kdtree_destroy(void* _p)
{
    kdtree *p = (kdtree*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_destroy(&p->xy);
    ae_vector_destroy(&p->tags);
    ae_vector_destroy(&p->boxmin);
    ae_vector_destroy(&p->boxmax);
    ae_vector_destroy(&p->nodes);
    ae_vector_destroy(&p->splits);
    ae_vector_destroy(&p->x);
    ae_vector_destroy(&p->idx);
    ae_vector_destroy(&p->r);
    ae_vector_destroy(&p->buf);
    ae_vector_destroy(&p->curboxmin);
    ae_vector_destroy(&p->curboxmax);
}




/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Creates and returns XDebugRecord1 structure:
* integer and complex fields of Rec1 are set to 1 and 1+i correspondingly
* array field of Rec1 is set to [2,3]

  -- ALGLIB --
     Copyright 27.05.2014 by Bochkanov Sergey
*************************************************************************/
void xdebuginitrecord1(xdebugrecord1* rec1, ae_state *_state)
{

    _xdebugrecord1_clear(rec1);

    rec1->i = 1;
    rec1->c.x = (double)(1);
    rec1->c.y = (double)(1);
    ae_vector_set_length(&rec1->a, 2, _state);
    rec1->a.ptr.p_double[0] = (double)(2);
    rec1->a.ptr.p_double[1] = (double)(3);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Counts number of True values in the boolean 1D array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugb1count(/* Boolean */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t result;


    result = 0;
    for(i=0; i<=a->cnt-1; i++)
    {
        if( a->ptr.p_bool[i] )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by NOT(a[i]).
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1not(/* Boolean */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_bool[i] = !a->ptr.p_bool[i];
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1appendcopy(/* Boolean */ ae_vector* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector b;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&b, 0, DT_BOOL, _state);

    ae_vector_set_length(&b, a->cnt, _state);
    for(i=0; i<=b.cnt-1; i++)
    {
        b.ptr.p_bool[i] = a->ptr.p_bool[i];
    }
    ae_vector_set_length(a, 2*b.cnt, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_bool[i] = b.ptr.p_bool[i%b.cnt];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered elements set to True.
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb1outeven(ae_int_t n,
     /* Boolean */ ae_vector* a,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(a);

    ae_vector_set_length(a, n, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_bool[i] = i%2==0;
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugi1sum(/* Integer */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t result;


    result = 0;
    for(i=0; i<=a->cnt-1; i++)
    {
        result = result+a->ptr.p_int[i];
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1neg(/* Integer */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_int[i] = -a->ptr.p_int[i];
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1appendcopy(/* Integer */ ae_vector* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector b;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&b, 0, DT_INT, _state);

    ae_vector_set_length(&b, a->cnt, _state);
    for(i=0; i<=b.cnt-1; i++)
    {
        b.ptr.p_int[i] = a->ptr.p_int[i];
    }
    ae_vector_set_length(a, 2*b.cnt, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_int[i] = b.ptr.p_int[i%b.cnt];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[I] set to I, and odd-numbered
ones set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi1outeven(ae_int_t n,
     /* Integer */ ae_vector* a,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(a);

    ae_vector_set_length(a, n, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        if( i%2==0 )
        {
            a->ptr.p_int[i] = i;
        }
        else
        {
            a->ptr.p_int[i] = 0;
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugr1sum(/* Real    */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;
    double result;


    result = (double)(0);
    for(i=0; i<=a->cnt-1; i++)
    {
        result = result+a->ptr.p_double[i];
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1neg(/* Real    */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_double[i] = -a->ptr.p_double[i];
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1appendcopy(/* Real    */ ae_vector* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector b;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&b, 0, DT_REAL, _state);

    ae_vector_set_length(&b, a->cnt, _state);
    for(i=0; i<=b.cnt-1; i++)
    {
        b.ptr.p_double[i] = a->ptr.p_double[i];
    }
    ae_vector_set_length(a, 2*b.cnt, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_double[i] = b.ptr.p_double[i%b.cnt];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[I] set to I*0.25,
and odd-numbered ones are set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr1outeven(ae_int_t n,
     /* Real    */ ae_vector* a,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(a);

    ae_vector_set_length(a, n, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        if( i%2==0 )
        {
            a->ptr.p_double[i] = i*0.25;
        }
        else
        {
            a->ptr.p_double[i] = (double)(0);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_complex xdebugc1sum(/* Complex */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;
    ae_complex result;


    result = ae_complex_from_i(0);
    for(i=0; i<=a->cnt-1; i++)
    {
        result = ae_c_add(result,a->ptr.p_complex[i]);
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -A[I]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1neg(/* Complex */ ae_vector* a, ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_complex[i] = ae_c_neg(a->ptr.p_complex[i]);
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Appends copy of array to itself.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1appendcopy(/* Complex */ ae_vector* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector b;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&b, 0, DT_COMPLEX, _state);

    ae_vector_set_length(&b, a->cnt, _state);
    for(i=0; i<=b.cnt-1; i++)
    {
        b.ptr.p_complex[i] = a->ptr.p_complex[i];
    }
    ae_vector_set_length(a, 2*b.cnt, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        a->ptr.p_complex[i] = b.ptr.p_complex[i%b.cnt];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate N-element array with even-numbered A[K] set to (x,y) = (K*0.25, K*0.125)
and odd-numbered ones are set to 0.

Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc1outeven(ae_int_t n,
     /* Complex */ ae_vector* a,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(a);

    ae_vector_set_length(a, n, _state);
    for(i=0; i<=a->cnt-1; i++)
    {
        if( i%2==0 )
        {
            a->ptr.p_complex[i].x = i*0.250;
            a->ptr.p_complex[i].y = i*0.125;
        }
        else
        {
            a->ptr.p_complex[i] = ae_complex_from_i(0);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Counts number of True values in the boolean 2D array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugb2count(/* Boolean */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t result;


    result = 0;
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            if( a->ptr.pp_bool[i][j] )
            {
                result = result+1;
            }
        }
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by NOT(a[i]).
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2not(/* Boolean */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_bool[i][j] = !a->ptr.pp_bool[i][j];
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2transpose(/* Boolean */ ae_matrix* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix b;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&b, 0, 0, DT_BOOL, _state);

    ae_matrix_set_length(&b, a->rows, a->cols, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            b.ptr.pp_bool[i][j] = a->ptr.pp_bool[i][j];
        }
    }
    ae_matrix_set_length(a, b.cols, b.rows, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            a->ptr.pp_bool[j][i] = b.ptr.pp_bool[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J)>0"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugb2outsin(ae_int_t m,
     ae_int_t n,
     /* Boolean */ ae_matrix* a,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(a);

    ae_matrix_set_length(a, m, n, _state);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_bool[i][j] = ae_fp_greater(ae_sin((double)(3*i+5*j), _state),(double)(0));
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_int_t xdebugi2sum(/* Integer */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t result;


    result = 0;
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            result = result+a->ptr.pp_int[i][j];
        }
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2neg(/* Integer */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_int[i][j] = -a->ptr.pp_int[i][j];
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2transpose(/* Integer */ ae_matrix* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix b;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&b, 0, 0, DT_INT, _state);

    ae_matrix_set_length(&b, a->rows, a->cols, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            b.ptr.pp_int[i][j] = a->ptr.pp_int[i][j];
        }
    }
    ae_matrix_set_length(a, b.cols, b.rows, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            a->ptr.pp_int[j][i] = b.ptr.pp_int[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sign(Sin(3*I+5*J))"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugi2outsin(ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_matrix* a,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(a);

    ae_matrix_set_length(a, m, n, _state);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_int[i][j] = ae_sign(ae_sin((double)(3*i+5*j), _state), _state);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugr2sum(/* Real    */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double result;


    result = (double)(0);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            result = result+a->ptr.pp_double[i][j];
        }
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2neg(/* Real    */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_double[i][j] = -a->ptr.pp_double[i][j];
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2transpose(/* Real    */ ae_matrix* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix b;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&b, 0, 0, DT_REAL, _state);

    ae_matrix_set_length(&b, a->rows, a->cols, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            b.ptr.pp_double[i][j] = a->ptr.pp_double[i][j];
        }
    }
    ae_matrix_set_length(a, b.cols, b.rows, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            a->ptr.pp_double[j][i] = b.ptr.pp_double[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J)"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugr2outsin(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(a);

    ae_matrix_set_length(a, m, n, _state);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_double[i][j] = ae_sin((double)(3*i+5*j), _state);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of elements in the array.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
ae_complex xdebugc2sum(/* Complex */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_complex result;


    result = ae_complex_from_i(0);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            result = ae_c_add(result,a->ptr.pp_complex[i][j]);
        }
    }
    return result;
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Replace all values in array by -a[i,j]
Array is passed using "shared" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2neg(/* Complex */ ae_matrix* a, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_complex[i][j] = ae_c_neg(a->ptr.pp_complex[i][j]);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Transposes array.
Array is passed using "var" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2transpose(/* Complex */ ae_matrix* a, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix b;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&b, 0, 0, DT_COMPLEX, _state);

    ae_matrix_set_length(&b, a->rows, a->cols, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            b.ptr.pp_complex[i][j] = a->ptr.pp_complex[i][j];
        }
    }
    ae_matrix_set_length(a, b.cols, b.rows, _state);
    for(i=0; i<=b.rows-1; i++)
    {
        for(j=0; j<=b.cols-1; j++)
        {
            a->ptr.pp_complex[j][i] = b.ptr.pp_complex[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Generate MxN matrix with elements set to "Sin(3*I+5*J),Cos(3*I+5*J)"
Array is passed using "out" convention.

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
void xdebugc2outsincos(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(a);

    ae_matrix_set_length(a, m, n, _state);
    for(i=0; i<=a->rows-1; i++)
    {
        for(j=0; j<=a->cols-1; j++)
        {
            a->ptr.pp_complex[i][j].x = ae_sin((double)(3*i+5*j), _state);
            a->ptr.pp_complex[i][j].y = ae_cos((double)(3*i+5*j), _state);
        }
    }
}


/*************************************************************************
This is debug function intended for testing ALGLIB interface generator.
Never use it in any real life project.

Returns sum of a[i,j]*(1+b[i,j]) such that c[i,j] is True

  -- ALGLIB --
     Copyright 11.10.2013 by Bochkanov Sergey
*************************************************************************/
double xdebugmaskedbiasedproductsum(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* b,
     /* Boolean */ ae_matrix* c,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double result;


    ae_assert(m>=a->rows, "Assertion failed", _state);
    ae_assert(m>=b->rows, "Assertion failed", _state);
    ae_assert(m>=c->rows, "Assertion failed", _state);
    ae_assert(n>=a->cols, "Assertion failed", _state);
    ae_assert(n>=b->cols, "Assertion failed", _state);
    ae_assert(n>=c->cols, "Assertion failed", _state);
    result = 0.0;
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( c->ptr.pp_bool[i][j] )
            {
                result = result+a->ptr.pp_double[i][j]*(1+b->ptr.pp_double[i][j]);
            }
        }
    }
    return result;
}


void _xdebugrecord1_init(void* _p, ae_state *_state)
{
    xdebugrecord1 *p = (xdebugrecord1*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->a, 0, DT_REAL, _state);
}


void _xdebugrecord1_init_copy(void* _dst, void* _src, ae_state *_state)
{
    xdebugrecord1 *dst = (xdebugrecord1*)_dst;
    xdebugrecord1 *src = (xdebugrecord1*)_src;
    dst->i = src->i;
    dst->c = src->c;
    ae_vector_init_copy(&dst->a, &src->a, _state);
}


void _xdebugrecord1_clear(void* _p)
{
    xdebugrecord1 *p = (xdebugrecord1*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->a);
}


void _xdebugrecord1_destroy(void* _p)
{
    xdebugrecord1 *p = (xdebugrecord1*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->a);
}



}

