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
#ifndef _solvers_pkg_h
#define _solvers_pkg_h
#include "ap.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "alglibmisc.h"

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
typedef struct
{
    double r1;
    double rinf;
} densesolverreport;
typedef struct
{
    double r2;
    ae_matrix cx;
    ae_int_t n;
    ae_int_t k;
} densesolverlsreport;
typedef struct
{
    normestimatorstate nes;
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t m;
    ae_int_t prectype;
    ae_vector ui;
    ae_vector uip1;
    ae_vector vi;
    ae_vector vip1;
    ae_vector omegai;
    ae_vector omegaip1;
    double alphai;
    double alphaip1;
    double betai;
    double betaip1;
    double phibari;
    double phibarip1;
    double phii;
    double rhobari;
    double rhobarip1;
    double rhoi;
    double ci;
    double si;
    double theta;
    double lambdai;
    ae_vector d;
    double anorm;
    double bnorm2;
    double dnorm;
    double r2;
    ae_vector x;
    ae_vector mv;
    ae_vector mtv;
    double epsa;
    double epsb;
    double epsc;
    ae_int_t maxits;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_vector tmpd;
    ae_vector tmpx;
    rcommstate rstate;
} linlsqrstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
} linlsqrreport;
typedef struct
{
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t prectype;
    ae_vector cx;
    ae_vector cr;
    ae_vector cz;
    ae_vector p;
    ae_vector r;
    ae_vector z;
    double alpha;
    double beta;
    double r2;
    double meritfunction;
    ae_vector x;
    ae_vector mv;
    ae_vector pv;
    double vmv;
    ae_vector startx;
    double epsf;
    ae_int_t maxits;
    ae_int_t itsbeforerestart;
    ae_int_t itsbeforerupdate;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_vector tmpd;
    rcommstate rstate;
} lincgstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
    double r2;
} lincgreport;
typedef struct
{
    ae_int_t n;
    ae_int_t m;
    double epsf;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_bool needf;
    ae_bool needfij;
    ae_bool xupdated;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfunc;
    ae_int_t repnjac;
    ae_int_t repterminationtype;
    ae_vector xbase;
    double fbase;
    double fprev;
    ae_vector candstep;
    ae_vector rightpart;
    ae_vector cgbuf;
} nleqstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t terminationtype;
} nleqreport;
typedef struct
{
    double maxerr;
} polynomialsolverreport;

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{

/*************************************************************************

*************************************************************************/
class _densesolverreport_owner
{
public:
    _densesolverreport_owner();
    _densesolverreport_owner(const _densesolverreport_owner &rhs);
    _densesolverreport_owner& operator=(const _densesolverreport_owner &rhs);
    virtual ~_densesolverreport_owner();
    alglib_impl::densesolverreport* c_ptr();
    alglib_impl::densesolverreport* c_ptr() const;
protected:
    alglib_impl::densesolverreport *p_struct;
};
class densesolverreport : public _densesolverreport_owner
{
public:
    densesolverreport();
    densesolverreport(const densesolverreport &rhs);
    densesolverreport& operator=(const densesolverreport &rhs);
    virtual ~densesolverreport();
    double &r1;
    double &rinf;

};


/*************************************************************************

*************************************************************************/
class _densesolverlsreport_owner
{
public:
    _densesolverlsreport_owner();
    _densesolverlsreport_owner(const _densesolverlsreport_owner &rhs);
    _densesolverlsreport_owner& operator=(const _densesolverlsreport_owner &rhs);
    virtual ~_densesolverlsreport_owner();
    alglib_impl::densesolverlsreport* c_ptr();
    alglib_impl::densesolverlsreport* c_ptr() const;
protected:
    alglib_impl::densesolverlsreport *p_struct;
};
class densesolverlsreport : public _densesolverlsreport_owner
{
public:
    densesolverlsreport();
    densesolverlsreport(const densesolverlsreport &rhs);
    densesolverlsreport& operator=(const densesolverlsreport &rhs);
    virtual ~densesolverlsreport();
    double &r2;
    real_2d_array cx;
    ae_int_t &n;
    ae_int_t &k;

};

/*************************************************************************
This object stores state of the LinLSQR method.

You should use ALGLIB functions to work with this object.
*************************************************************************/
class _linlsqrstate_owner
{
public:
    _linlsqrstate_owner();
    _linlsqrstate_owner(const _linlsqrstate_owner &rhs);
    _linlsqrstate_owner& operator=(const _linlsqrstate_owner &rhs);
    virtual ~_linlsqrstate_owner();
    alglib_impl::linlsqrstate* c_ptr();
    alglib_impl::linlsqrstate* c_ptr() const;
protected:
    alglib_impl::linlsqrstate *p_struct;
};
class linlsqrstate : public _linlsqrstate_owner
{
public:
    linlsqrstate();
    linlsqrstate(const linlsqrstate &rhs);
    linlsqrstate& operator=(const linlsqrstate &rhs);
    virtual ~linlsqrstate();

};


/*************************************************************************

*************************************************************************/
class _linlsqrreport_owner
{
public:
    _linlsqrreport_owner();
    _linlsqrreport_owner(const _linlsqrreport_owner &rhs);
    _linlsqrreport_owner& operator=(const _linlsqrreport_owner &rhs);
    virtual ~_linlsqrreport_owner();
    alglib_impl::linlsqrreport* c_ptr();
    alglib_impl::linlsqrreport* c_ptr() const;
protected:
    alglib_impl::linlsqrreport *p_struct;
};
class linlsqrreport : public _linlsqrreport_owner
{
public:
    linlsqrreport();
    linlsqrreport(const linlsqrreport &rhs);
    linlsqrreport& operator=(const linlsqrreport &rhs);
    virtual ~linlsqrreport();
    ae_int_t &iterationscount;
    ae_int_t &nmv;
    ae_int_t &terminationtype;

};

/*************************************************************************
This object stores state of the linear CG method.

You should use ALGLIB functions to work with this object.
Never try to access its fields directly!
*************************************************************************/
class _lincgstate_owner
{
public:
    _lincgstate_owner();
    _lincgstate_owner(const _lincgstate_owner &rhs);
    _lincgstate_owner& operator=(const _lincgstate_owner &rhs);
    virtual ~_lincgstate_owner();
    alglib_impl::lincgstate* c_ptr();
    alglib_impl::lincgstate* c_ptr() const;
protected:
    alglib_impl::lincgstate *p_struct;
};
class lincgstate : public _lincgstate_owner
{
public:
    lincgstate();
    lincgstate(const lincgstate &rhs);
    lincgstate& operator=(const lincgstate &rhs);
    virtual ~lincgstate();

};


/*************************************************************************

*************************************************************************/
class _lincgreport_owner
{
public:
    _lincgreport_owner();
    _lincgreport_owner(const _lincgreport_owner &rhs);
    _lincgreport_owner& operator=(const _lincgreport_owner &rhs);
    virtual ~_lincgreport_owner();
    alglib_impl::lincgreport* c_ptr();
    alglib_impl::lincgreport* c_ptr() const;
protected:
    alglib_impl::lincgreport *p_struct;
};
class lincgreport : public _lincgreport_owner
{
public:
    lincgreport();
    lincgreport(const lincgreport &rhs);
    lincgreport& operator=(const lincgreport &rhs);
    virtual ~lincgreport();
    ae_int_t &iterationscount;
    ae_int_t &nmv;
    ae_int_t &terminationtype;
    double &r2;

};

/*************************************************************************

*************************************************************************/
class _nleqstate_owner
{
public:
    _nleqstate_owner();
    _nleqstate_owner(const _nleqstate_owner &rhs);
    _nleqstate_owner& operator=(const _nleqstate_owner &rhs);
    virtual ~_nleqstate_owner();
    alglib_impl::nleqstate* c_ptr();
    alglib_impl::nleqstate* c_ptr() const;
protected:
    alglib_impl::nleqstate *p_struct;
};
class nleqstate : public _nleqstate_owner
{
public:
    nleqstate();
    nleqstate(const nleqstate &rhs);
    nleqstate& operator=(const nleqstate &rhs);
    virtual ~nleqstate();
    ae_bool &needf;
    ae_bool &needfij;
    ae_bool &xupdated;
    double &f;
    real_1d_array fi;
    real_2d_array j;
    real_1d_array x;

};


/*************************************************************************

*************************************************************************/
class _nleqreport_owner
{
public:
    _nleqreport_owner();
    _nleqreport_owner(const _nleqreport_owner &rhs);
    _nleqreport_owner& operator=(const _nleqreport_owner &rhs);
    virtual ~_nleqreport_owner();
    alglib_impl::nleqreport* c_ptr();
    alglib_impl::nleqreport* c_ptr() const;
protected:
    alglib_impl::nleqreport *p_struct;
};
class nleqreport : public _nleqreport_owner
{
public:
    nleqreport();
    nleqreport(const nleqreport &rhs);
    nleqreport& operator=(const nleqreport &rhs);
    virtual ~nleqreport();
    ae_int_t &iterationscount;
    ae_int_t &nfunc;
    ae_int_t &njac;
    ae_int_t &terminationtype;

};

/*************************************************************************

*************************************************************************/
class _polynomialsolverreport_owner
{
public:
    _polynomialsolverreport_owner();
    _polynomialsolverreport_owner(const _polynomialsolverreport_owner &rhs);
    _polynomialsolverreport_owner& operator=(const _polynomialsolverreport_owner &rhs);
    virtual ~_polynomialsolverreport_owner();
    alglib_impl::polynomialsolverreport* c_ptr();
    alglib_impl::polynomialsolverreport* c_ptr() const;
protected:
    alglib_impl::polynomialsolverreport *p_struct;
};
class polynomialsolverreport : public _polynomialsolverreport_owner
{
public:
    polynomialsolverreport();
    polynomialsolverreport(const polynomialsolverreport &rhs);
    polynomialsolverreport& operator=(const polynomialsolverreport &rhs);
    virtual ~polynomialsolverreport();
    double &maxerr;

};

/*************************************************************************
Dense solver for A*x=b with N*N real matrix A and N*1 real vectorx  x  and
b. This is "slow-but-feature rich" version of the  linear  solver.  Faster
version is RMatrixSolveFast() function.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations. It also very significant on small matrices.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, RMatrixSolveFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void smp_rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix, x  and  b  are  vectors.  This is a "fast" version of  linear
solver which does NOT provide  any  additional  functions  like  condition
number estimation or iterative refinement.

Algorithm features:
* efficient algorithm O(N^3) complexity
* no performance overhead from additional functionality

If you need condition number estimation or iterative refinement, use  more
feature-rich version - RMatrixSolve().

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 16.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
void smp_rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices). This is  "slow-but-robust"  version  of  linear
solver with additional functionality  like  condition  number  estimation.
There also exists faster version - RMatrixSolveMFast().

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* optional iterative refinement
* O(N^3+M*N^2) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations. It also very significant on small matrices.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, RMatrixSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is ill conditioned or singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void smp_rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices). This is "fast" version of linear  solver  which
does NOT offer additional functions like condition  number  estimation  or
iterative refinement.

Algorithm features:
* O(N^3+M*N^2) complexity
* no additional functionality, highest performance

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix given by its LU decomposition, x and b are real vectors.  This
is "slow-but-robust" version of the linear LU-based solver. Faster version
is RMatrixLUSolveFast() function.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! RMatrixLUSolveFast() function.

INPUT PARAMETERS
    LUA     -   array[N,N], LU decomposition, RMatrixLU result
    P       -   array[N], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix given by its LU decomposition, x and b are real vectors.  This
is "fast-without-any-checks" version of the linear LU-based solver. Slower
but more robust version is RMatrixLUSolve() function.

Algorithm features:
* O(N^2) complexity
* fast algorithm without ANY additional checks, just triangular solver

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves  task  with  multiple  right  parts
(where b and x are NxM matrices). This  is  "robust-but-slow"  version  of
LU-based solver which performs additional  checks  for  non-degeneracy  of
inputs (condition number estimation). If you need  best  performance,  use
"fast-without-any-checks" version, RMatrixLUSolveMFast().

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! RMatrixLUSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. Triangular solver is relatively easy to parallelize.
  ! However, parallelization will be efficient  only for  large number  of
  ! right parts M.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    LUA     -   array[N,N], LU decomposition, RMatrixLU result
    P       -   array[N], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void smp_rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves  task  with  multiple  right parts,
where b and x are NxM matrices.  This is "fast-without-any-checks" version
of LU-based solver. It does not estimate  condition number  of  a  system,
so it is extremely fast. If you need better detection  of  near-degenerate
cases, use RMatrixLUSolveM() function.

Algorithm features:
* O(M*N^2) complexity
* fast algorithm without ANY additional checks, just triangular solver

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. Triangular solver is relatively easy to parallelize.
  ! However, parallelization will be efficient  only for  large number  of
  ! right parts M.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
both A and its LU decomposition.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);


/*************************************************************************
Dense solver.

Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);


/*************************************************************************
Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
matrices  X  and  B.  "Slow-but-feature-rich"   version   which   provides
additional functions, at the cost of slower  performance.  Faster  version
may be invoked with CMatrixSolveMFast() function.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3+M*N^2) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, CMatrixSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void smp_cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);


/*************************************************************************
Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
matrices  X  and  B.  "Fast-but-lightweight" version which  provides  just
triangular solver - and no additional functions like iterative  refinement
or condition number estimation.

Algorithm features:
* O(N^3+M*N^2) complexity
* no additional time consuming functions

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 16.03.2015 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
vectors x and b. "Slow-but-feature-rich" version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, CMatrixSolveFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void smp_cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);


/*************************************************************************
Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
vectors x and b. "Fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3) complexity
* no additional time consuming features, just triangular solver

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that LU decomposition  is  harder  to
  ! parallelize than, say, matrix-matrix  product  -  this  algorithm  has
  ! many internal synchronization points which can not be avoided. However
  ! parallelism starts to be profitable starting  from  N=1024,  achieving
  ! near-linear speedup for N=4096 or higher.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
void smp_cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
and N*M matrices X and B (multiple right sides).   "Slow-but-feature-rich"
version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! CMatrixLUSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. Triangular solver is relatively easy to parallelize.
  ! However, parallelization will be efficient  only for  large number  of
  ! right parts M.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void smp_cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);


/*************************************************************************
Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
and N*M matrices X and B (multiple  right  sides).  "Fast-but-lightweight"
version of the solver.

Algorithm features:
* O(M*N^2) complexity
* no additional time-consuming features

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. Triangular solver is relatively easy to parallelize.
  ! However, parallelization will be efficient  only for  large number  of
  ! right parts M.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Complex dense linear solver for A*x=b with complex N*N A  given  by its LU
decomposition and N*1 vectors x and b. This is  "slow-but-robust"  version
of  the  complex  linear  solver  with  additional  features   which   add
significant performance overhead. Faster version  is  CMatrixLUSolveFast()
function.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! CMatrixLUSolveFast() function.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);


/*************************************************************************
Complex dense linear solver for A*x=b with N*N complex A given by  its  LU
decomposition and N*1 vectors x and b. This is  fast  lightweight  version
of solver, which is significantly faster than CMatrixLUSolve(),  but  does
not provide additional information (like condition numbers).

Algorithm features:
* O(N^2) complexity
* no additional time-consuming features, just triangular solver

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

NOTE: unlike  CMatrixLUSolve(),  this   function   does   NOT   check  for
      near-degeneracy of input matrix. It  checks  for  EXACT  degeneracy,
      because this check is easy to do. However,  very  badly  conditioned
      matrices may went unnoticed.


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);


/*************************************************************************
Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
N*M vectors X and B. It is "slow-but-feature-rich" version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, SPDMatrixSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or non-SPD.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void smp_spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
N*M vectors X and B. It is "fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
matrix A,  N*1 vectors x and b.  "Slow-but-feature-rich"  version  of  the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, SPDMatrixSolveFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or non-SPD.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
void smp_spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);


/*************************************************************************
Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
matrix A,  N*1 vectors x and  b.  "Fast-but-lightweight"  version  of  the
solver.

Algorithm features:
* O(N^3) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features like condition number estimation

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or non-SPD
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
void smp_spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*M vectors X and B. It  is  "slow-but-
feature-rich" version of the solver which estimates  condition  number  of
the system.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver. Amount of  overhead  introduced  depends  on  M  (the
           ! larger - the more efficient).
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! SPDMatrixCholeskySolveMFast() function.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or badly conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0 contains solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
void smp_spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*M vectors X and B. It  is  "fast-but-
lightweight" version of  the  solver  which  just  solves  linear  system,
without any additional functions.

Algorithm features:
* O(M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional functionality

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or badly conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N]:
                * for info>0 overwritten by solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense solver for A*x=b with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*1 real vectors x and b. This is "slow-
but-feature-rich"  version  of  the  solver  which,  in  addition  to  the
solution, performs condition number estimation.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! SPDMatrixCholeskySolveFast() function.

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0  - solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);


/*************************************************************************
Dense solver for A*x=b with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*1 real vectors x and b. This is "fast-
but-lightweight" version of the solver.

Algorithm features:
* O(N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional features

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * for info>0  - overwritten by solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
N*M  complex  matrices  X  and  B.  "Slow-but-feature-rich" version of the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems (N<100).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixSolveMFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void smp_hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);


/*************************************************************************
Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
N*M complex matrices X and B. "Fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features like condition number estimation

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly  singular or is not positive definite.
                        B is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[0..N-1]:
                * overwritten by solution
                * zeros, if problem was not solved

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
N*1 complex vectors  x  and  b.  "Slow-but-feature-rich"  version  of  the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, HPDMatrixSolveFast() function.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
void smp_hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);


/*************************************************************************
Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
N*1 complex vectors  x  and  b.  "Fast-but-lightweight"  version  of   the
solver without additional functions.

Algorithm features:
* O(N^3) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming functions

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Say, on SSE2-capable CPU with N=1024, HPC ALGLIB will be:
  ! * about 2-3x faster than ALGLIB for C++ without MKL
  ! * about 7-10x faster than "pure C#" edition of ALGLIB
  ! Difference in performance will be more striking  on  newer  CPU's with
  ! support for newer SIMD instructions. Generally,  MKL  accelerates  any
  ! problem whose size is at least 128, with best  efficiency achieved for
  ! N's larger than 512.
  !
  ! Commercial edition of ALGLIB also supports multithreaded  acceleration
  ! of this function. We should note that Cholesky decomposition is harder
  ! to parallelize than, say, matrix-matrix product - this  algorithm  has
  ! several synchronization points which  can  not  be  avoided.  However,
  ! parallelism starts to be profitable starting from N=500.
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or not positive definite
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[0..N-1]:
                * overwritten by solution
                * zeros, if A is exactly singular (diagonal of its LU
                  decomposition has exact zeros).

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
void smp_hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
"slow-but-feature-rich" version of the solver which, in  addition  to  the
solution, estimates condition number of the system.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver. Amount of  overhead  introduced  depends  on  M  (the
           ! larger - the more efficient).
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large Cholesky decomposition.  However,  if  you call
           ! this  function  many  times  for  the same  left  side,  this
           ! overhead BECOMES significant. It  also   becomes  significant
           ! for small-scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixCholeskySolveMFast() function.


INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0 contains solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
void smp_hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);


/*************************************************************************
Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
"fast-but-lightweight" version of the solver.

Algorithm features:
* O(M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time-consuming features

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N]:
                * for info>0 overwritten by solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
void smp_hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);


/*************************************************************************
Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
"slow-but-feature-rich" version of the solver  which  estimates  condition
number of the system.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixCholeskySolveFast() function.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0  - solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);


/*************************************************************************
Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
"fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time-consuming features

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        B is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * for info>0  - overwritten by solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);


/*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the usual sense is returned.

Algorithm features:
* automatic detection (and correct handling!) of degenerate cases
* iterative refinement
* O(N^3) complexity

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes one  important  improvement   of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison.
  !
  ! Generally, commercial ALGLIB is several times faster than  open-source
  ! generic C edition, and many times faster than open-source C# edition.
  !
  ! Multithreaded acceleration is only partially supported (some parts are
  ! optimized, but most - are not).
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.

OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B (even for singular A)
                * zeros, if SVD subroutine failed

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
void smp_rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);

/*************************************************************************
This function initializes linear LSQR Solver. This solver is used to solve
non-symmetric (and, possibly, non-square) problems. Least squares solution
is returned for non-compatible systems.

USAGE:
1. User initializes algorithm state with LinLSQRCreate() call
2. User tunes solver parameters with  LinLSQRSetCond() and other functions
3. User  calls  LinLSQRSolveSparse()  function which takes algorithm state
   and SparseMatrix object.
4. User calls LinLSQRResults() to get solution
5. Optionally, user may call LinLSQRSolveSparse() again to  solve  another
   problem  with different matrix and/or right part without reinitializing
   LinLSQRState structure.

INPUT PARAMETERS:
    M       -   number of rows in A
    N       -   number of variables, N>0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state);


/*************************************************************************
This  function  changes  preconditioning  settings of LinLSQQSolveSparse()
function. By default, SolveSparse() uses diagonal preconditioner,  but  if
you want to use solver without preconditioning, you can call this function
which forces solver to use unit matrix for preconditioning.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetprecunit(const linlsqrstate &state);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
preconditioner. This preconditioning mode is active by default.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetprecdiag(const linlsqrstate &state);


/*************************************************************************
This function sets optional Tikhonov regularization coefficient.
It is zero by default.

INPUT PARAMETERS:
    LambdaI -   regularization factor, LambdaI>=0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai);


/*************************************************************************
Procedure for solution of A*x=b with sparse A.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse M*N matrix in the CRS format (you MUST contvert  it
                to CRS format  by  calling  SparseConvertToCRS()  function
                BEFORE you pass it to this function).
    B       -   right part, array[M]

RESULT:
    This function returns no result.
    You can get solution by calling LinCGResults()

NOTE: this function uses lightweight preconditioning -  multiplication  by
      inverse of diag(A). If you want, you can turn preconditioning off by
      calling LinLSQRSetPrecUnit(). However, preconditioning cost is   low
      and preconditioner is very important for solution  of  badly  scaled
      problems.

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b);


/*************************************************************************
This function sets stopping criteria.

INPUT PARAMETERS:
    EpsA    -   algorithm will be stopped if ||A^T*Rk||/(||A||*||Rk||)<=EpsA.
    EpsB    -   algorithm will be stopped if ||Rk||<=EpsB*||B||
    MaxIts  -   algorithm will be stopped if number of iterations
                more than MaxIts.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE: if EpsA,EpsB,EpsC and MaxIts are zero then these variables will
be setted as default values.

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits);


/*************************************************************************
LSQR solver: results.

This function must be called after LinLSQRSolve

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[N], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    *  1    ||Rk||<=EpsB*||B||
                    *  4    ||A^T*Rk||/(||A||*||Rk||)<=EpsA
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            X contains best point found so far.
                            (sometimes returned on singular systems)
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep);

/*************************************************************************
This function initializes linear CG Solver. This solver is used  to  solve
symmetric positive definite problems. If you want  to  solve  nonsymmetric
(or non-positive definite) problem you may use LinLSQR solver provided  by
ALGLIB.

USAGE:
1. User initializes algorithm state with LinCGCreate() call
2. User tunes solver parameters with  LinCGSetCond() and other functions
3. Optionally, user sets starting point with LinCGSetStartingPoint()
4. User  calls LinCGSolveSparse() function which takes algorithm state and
   SparseMatrix object.
5. User calls LinCGResults() to get solution
6. Optionally, user may call LinCGSolveSparse()  again  to  solve  another
   problem  with different matrix and/or right part without reinitializing
   LinCGState structure.

INPUT PARAMETERS:
    N       -   problem dimension, N>0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgcreate(const ae_int_t n, lincgstate &state);


/*************************************************************************
This function sets starting point.
By default, zero starting point is used.

INPUT PARAMETERS:
    X       -   starting point, array[N]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function. By default, SolveSparse() uses diagonal preconditioner,  but  if
you want to use solver without preconditioning, you can call this function
which forces solver to use unit matrix for preconditioning.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void lincgsetprecunit(const lincgstate &state);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
preconditioner. This preconditioning mode is active by default.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void lincgsetprecdiag(const lincgstate &state);


/*************************************************************************
This function sets stopping criteria.

INPUT PARAMETERS:
    EpsF    -   algorithm will be stopped if norm of residual is less than
                EpsF*||b||.
    MaxIts  -   algorithm will be stopped if number of iterations is  more
                than MaxIts.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
If  both  EpsF  and  MaxIts  are  zero then small EpsF will be set to small
value.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits);


/*************************************************************************
Procedure for solution of A*x=b with sparse A.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse matrix in the CRS format (you MUST contvert  it  to
                CRS format by calling SparseConvertToCRS() function).
    IsUpper -   whether upper or lower triangle of A is used:
                * IsUpper=True  => only upper triangle is used and lower
                                   triangle is not referenced at all
                * IsUpper=False => only lower triangle is used and upper
                                   triangle is not referenced at all
    B       -   right part, array[N]

RESULT:
    This function returns no result.
    You can get solution by calling LinCGResults()

NOTE: this function uses lightweight preconditioning -  multiplication  by
      inverse of diag(A). If you want, you can turn preconditioning off by
      calling LinCGSetPrecUnit(). However, preconditioning cost is low and
      preconditioner  is  very  important  for  solution  of  badly scaled
      problems.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b);


/*************************************************************************
CG-solver: results.

This function must be called after LinCGSolve

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[N], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -5    input matrix is either not positive definite,
                            too large or too small
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep);


/*************************************************************************
This function sets restart frequency. By default, algorithm  is  restarted
after N subsequent iterations.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf);


/*************************************************************************
This function sets frequency of residual recalculations.

Algorithm updates residual r_k using iterative formula,  but  recalculates
it from scratch after each 10 iterations. It is done to avoid accumulation
of numerical errors and to stop algorithm when r_k starts to grow.

Such low update frequence (1/10) gives very  little  overhead,  but  makes
algorithm a bit more robust against numerical errors. However, you may
change it

INPUT PARAMETERS:
    Freq    -   desired update frequency, Freq>=0.
                Zero value means that no updates will be done.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetxrep(const lincgstate &state, const bool needxrep);

/*************************************************************************
                LEVENBERG-MARQUARDT-LIKE NONLINEAR SOLVER

DESCRIPTION:
This algorithm solves system of nonlinear equations
    F[0](x[0], ..., x[n-1])   = 0
    F[1](x[0], ..., x[n-1])   = 0
    ...
    F[M-1](x[0], ..., x[n-1]) = 0
with M/N do not necessarily coincide.  Algorithm  converges  quadratically
under following conditions:
    * the solution set XS is nonempty
    * for some xs in XS there exist such neighbourhood N(xs) that:
      * vector function F(x) and its Jacobian J(x) are continuously
        differentiable on N
      * ||F(x)|| provides local error bound on N, i.e. there  exists  such
        c1, that ||F(x)||>c1*distance(x,XS)
Note that these conditions are much more weaker than usual non-singularity
conditions. For example, algorithm will converge for any  affine  function
F (whether its Jacobian singular or not).


REQUIREMENTS:
Algorithm will request following information during its operation:
* function vector F[] and Jacobian matrix at given point X
* value of merit function f(x)=F[0]^2(x)+...+F[M-1]^2(x) at given point X


USAGE:
1. User initializes algorithm state with NLEQCreateLM() call
2. User tunes solver parameters with  NLEQSetCond(),  NLEQSetStpMax()  and
   other functions
3. User  calls  NLEQSolve()  function  which  takes  algorithm  state  and
   pointers (delegates, etc.) to callback functions which calculate  merit
   function value and Jacobian.
4. User calls NLEQResults() to get solution
5. Optionally, user may call NLEQRestartFrom() to  solve  another  problem
   with same parameters (N/M) but another starting  point  and/or  another
   function vector. NLEQRestartFrom() allows to reuse already  initialized
   structure.


INPUT PARAMETERS:
    N       -   space dimension, N>1:
                * if provided, only leading N elements of X are used
                * if not provided, determined automatically from size of X
    M       -   system size
    X       -   starting point


OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


NOTES:
1. you may tune stopping conditions with NLEQSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use NLEQSetStpMax() function to bound algorithm's steps.
3. this  algorithm  is  a  slightly  modified implementation of the method
   described  in  'Levenberg-Marquardt  method  for constrained  nonlinear
   equations with strong local convergence properties' by Christian Kanzow
   Nobuo Yamashita and Masao Fukushima and further  developed  in  'On the
   convergence of a New Levenberg-Marquardt Method'  by  Jin-yan  Fan  and
   Ya-Xiang Yuan.


  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state);
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state);


/*************************************************************************
This function sets stopping conditions for the nonlinear solver

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsF    -   >=0
                The subroutine finishes  its work if on k+1-th iteration
                the condition ||F||<=EpsF is satisfied
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsF=0 and MaxIts=0 simultaneously will lead to  automatic
stopping criterion selection (small EpsF).

NOTES:

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to NLEQSolve().

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetxrep(const nleqstate &state, const bool needxrep);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when target function  contains  exp()  or  other  fast
growing functions, and algorithm makes  too  large  steps  which  lead  to
overflow. This function allows us to reject steps that are too large  (and
therefore expose us to the possible overflow) without actually calculating
function value at the x+stp*d.

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetstpmax(const nleqstate &state, const double stpmax);


/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool nleqiteration(const nleqstate &state);


/*************************************************************************
This family of functions is used to launcn iterations of nonlinear solver

These functions accept following parameters:
    state   -   algorithm state
    func    -   callback which calculates function (or merit function)
                value func at given point x
    jac     -   callback which calculates function vector fi[]
                and Jacobian jac at given point x
    rep     -   optional callback which is called after each iteration
                can be NULL
    ptr     -   optional pointer which is passed to func/grad/hess/jac/rep
                can be NULL


  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey

*************************************************************************/
void nleqsolve(nleqstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr) = NULL,
    void *ptr = NULL);


/*************************************************************************
NLEQ solver results

INPUT PARAMETERS:
    State   -   algorithm state.

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -4    ERROR:  algorithm   has   converged   to   the
                            stationary point Xf which is local minimum  of
                            f=F[0]^2+...+F[m-1]^2, but is not solution  of
                            nonlinear system.
                    *  1    sqrt(f)<=EpsF.
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations
                * ActiveConstraints contains number of active constraints

  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep);


/*************************************************************************
NLEQ solver results

Buffered implementation of NLEQResults(), which uses pre-allocated  buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep);


/*************************************************************************
This  subroutine  restarts  CG  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used for reverse communication previously
                allocated with MinCGCreate call.
    X       -   new starting point.
    BndL    -   new lower bounds
    BndU    -   new upper bounds

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void nleqrestartfrom(const nleqstate &state, const real_1d_array &x);

/*************************************************************************
Polynomial root finding.

This function returns all roots of the polynomial
    P(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
Both real and complex roots are returned (see below).

INPUT PARAMETERS:
    A       -   array[N+1], polynomial coefficients:
                * A[0] is constant term
                * A[N] is a coefficient of X^N
    N       -   polynomial degree

OUTPUT PARAMETERS:
    X       -   array of complex roots:
                * for isolated real root, X[I] is strictly real: IMAGE(X[I])=0
                * complex roots are always returned in pairs - roots occupy
                  positions I and I+1, with:
                  * X[I+1]=Conj(X[I])
                  * IMAGE(X[I]) > 0
                  * IMAGE(X[I+1]) = -IMAGE(X[I]) < 0
                * multiple real roots may have non-zero imaginary part due
                  to roundoff errors. There is no reliable way to distinguish
                  real root of multiplicity 2 from two  complex  roots  in
                  the presence of roundoff errors.
    Rep     -   report, additional information, following fields are set:
                * Rep.MaxErr - max( |P(xi)| )  for  i=0..N-1.  This  field
                  allows to quickly estimate "quality" of the roots  being
                  returned.

NOTE:   this function uses companion matrix method to find roots. In  case
        internal EVD  solver  fails  do  find  eigenvalues,  exception  is
        generated.

NOTE:   roots are not "polished" and  no  matrix  balancing  is  performed
        for them.

  -- ALGLIB --
     Copyright 24.02.2014 by Bochkanov Sergey
*************************************************************************/
void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep);
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
void rmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _pexec_rmatrixsolve(/* Real    */ ae_matrix* a,
    ae_int_t n,
    /* Real    */ ae_vector* b,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_vector* x, ae_state *_state);
void rmatrixsolvefast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void _pexec_rmatrixsolvefast(/* Real    */ ae_matrix* a,
    ae_int_t n,
    /* Real    */ ae_vector* b,
    ae_int_t* info, ae_state *_state);
void rmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void _pexec_rmatrixsolvem(/* Real    */ ae_matrix* a,
    ae_int_t n,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_bool rfs,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_matrix* x, ae_state *_state);
void rmatrixsolvemfast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_rmatrixsolvemfast(/* Real    */ ae_matrix* a,
    ae_int_t n,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void rmatrixlusolve(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void rmatrixlusolvefast(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void rmatrixlusolvem(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void _pexec_rmatrixlusolvem(/* Real    */ ae_matrix* lua,
    /* Integer */ ae_vector* p,
    ae_int_t n,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_matrix* x, ae_state *_state);
void rmatrixlusolvemfast(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_rmatrixlusolvemfast(/* Real    */ ae_matrix* lua,
    /* Integer */ ae_vector* p,
    ae_int_t n,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void rmatrixmixedsolve(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void rmatrixmixedsolvem(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void cmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void _pexec_cmatrixsolvem(/* Complex */ ae_matrix* a,
    ae_int_t n,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_bool rfs,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_matrix* x, ae_state *_state);
void cmatrixsolvemfast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_cmatrixsolvemfast(/* Complex */ ae_matrix* a,
    ae_int_t n,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void cmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void _pexec_cmatrixsolve(/* Complex */ ae_matrix* a,
    ae_int_t n,
    /* Complex */ ae_vector* b,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_vector* x, ae_state *_state);
void cmatrixsolvefast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void _pexec_cmatrixsolvefast(/* Complex */ ae_matrix* a,
    ae_int_t n,
    /* Complex */ ae_vector* b,
    ae_int_t* info, ae_state *_state);
void cmatrixlusolvem(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void _pexec_cmatrixlusolvem(/* Complex */ ae_matrix* lua,
    /* Integer */ ae_vector* p,
    ae_int_t n,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_matrix* x, ae_state *_state);
void cmatrixlusolvemfast(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_cmatrixlusolvemfast(/* Complex */ ae_matrix* lua,
    /* Integer */ ae_vector* p,
    ae_int_t n,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void cmatrixlusolve(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void cmatrixlusolvefast(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void cmatrixmixedsolvem(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void cmatrixmixedsolve(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void spdmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void _pexec_spdmatrixsolvem(/* Real    */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_matrix* x, ae_state *_state);
void spdmatrixsolvemfast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_spdmatrixsolvemfast(/* Real    */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void spdmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _pexec_spdmatrixsolve(/* Real    */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_vector* b,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_vector* x, ae_state *_state);
void spdmatrixsolvefast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void _pexec_spdmatrixsolvefast(/* Real    */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_vector* b,
    ae_int_t* info, ae_state *_state);
void spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void _pexec_spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Real    */ ae_matrix* x, ae_state *_state);
void spdmatrixcholeskysolvemfast(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_spdmatrixcholeskysolvemfast(/* Real    */ ae_matrix* cha,
    ae_int_t n,
    ae_bool isupper,
    /* Real    */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void spdmatrixcholeskysolve(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void spdmatrixcholeskysolvefast(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void hpdmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void _pexec_hpdmatrixsolvem(/* Complex */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_matrix* x, ae_state *_state);
void hpdmatrixsolvemfast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_hpdmatrixsolvemfast(/* Complex */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void hpdmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void _pexec_hpdmatrixsolve(/* Complex */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_vector* b,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_vector* x, ae_state *_state);
void hpdmatrixsolvefast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void _pexec_hpdmatrixsolvefast(/* Complex */ ae_matrix* a,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_vector* b,
    ae_int_t* info, ae_state *_state);
void hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void _pexec_hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info,
    densesolverreport* rep,
    /* Complex */ ae_matrix* x, ae_state *_state);
void hpdmatrixcholeskysolvemfast(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void _pexec_hpdmatrixcholeskysolvemfast(/* Complex */ ae_matrix* cha,
    ae_int_t n,
    ae_bool isupper,
    /* Complex */ ae_matrix* b,
    ae_int_t m,
    ae_int_t* info, ae_state *_state);
void hpdmatrixcholeskysolve(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void hpdmatrixcholeskysolvefast(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void rmatrixsolvels(/* Real    */ ae_matrix* a,
     ae_int_t nrows,
     ae_int_t ncols,
     /* Real    */ ae_vector* b,
     double threshold,
     ae_int_t* info,
     densesolverlsreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _pexec_rmatrixsolvels(/* Real    */ ae_matrix* a,
    ae_int_t nrows,
    ae_int_t ncols,
    /* Real    */ ae_vector* b,
    double threshold,
    ae_int_t* info,
    densesolverlsreport* rep,
    /* Real    */ ae_vector* x, ae_state *_state);
void _densesolverreport_init(void* _p, ae_state *_state);
void _densesolverreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _densesolverreport_clear(void* _p);
void _densesolverreport_destroy(void* _p);
void _densesolverlsreport_init(void* _p, ae_state *_state);
void _densesolverlsreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _densesolverlsreport_clear(void* _p);
void _densesolverlsreport_destroy(void* _p);
void linlsqrcreate(ae_int_t m,
     ae_int_t n,
     linlsqrstate* state,
     ae_state *_state);
void linlsqrsetb(linlsqrstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void linlsqrsetprecunit(linlsqrstate* state, ae_state *_state);
void linlsqrsetprecdiag(linlsqrstate* state, ae_state *_state);
void linlsqrsetlambdai(linlsqrstate* state,
     double lambdai,
     ae_state *_state);
ae_bool linlsqriteration(linlsqrstate* state, ae_state *_state);
void linlsqrsolvesparse(linlsqrstate* state,
     sparsematrix* a,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void linlsqrsetcond(linlsqrstate* state,
     double epsa,
     double epsb,
     ae_int_t maxits,
     ae_state *_state);
void linlsqrresults(linlsqrstate* state,
     /* Real    */ ae_vector* x,
     linlsqrreport* rep,
     ae_state *_state);
void linlsqrsetxrep(linlsqrstate* state,
     ae_bool needxrep,
     ae_state *_state);
void linlsqrrestart(linlsqrstate* state, ae_state *_state);
void _linlsqrstate_init(void* _p, ae_state *_state);
void _linlsqrstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _linlsqrstate_clear(void* _p);
void _linlsqrstate_destroy(void* _p);
void _linlsqrreport_init(void* _p, ae_state *_state);
void _linlsqrreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _linlsqrreport_clear(void* _p);
void _linlsqrreport_destroy(void* _p);
void lincgcreate(ae_int_t n, lincgstate* state, ae_state *_state);
void lincgsetstartingpoint(lincgstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void lincgsetb(lincgstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void lincgsetprecunit(lincgstate* state, ae_state *_state);
void lincgsetprecdiag(lincgstate* state, ae_state *_state);
void lincgsetcond(lincgstate* state,
     double epsf,
     ae_int_t maxits,
     ae_state *_state);
ae_bool lincgiteration(lincgstate* state, ae_state *_state);
void lincgsolvesparse(lincgstate* state,
     sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void lincgresults(lincgstate* state,
     /* Real    */ ae_vector* x,
     lincgreport* rep,
     ae_state *_state);
void lincgsetrestartfreq(lincgstate* state,
     ae_int_t srf,
     ae_state *_state);
void lincgsetrupdatefreq(lincgstate* state,
     ae_int_t freq,
     ae_state *_state);
void lincgsetxrep(lincgstate* state, ae_bool needxrep, ae_state *_state);
void lincgrestart(lincgstate* state, ae_state *_state);
void _lincgstate_init(void* _p, ae_state *_state);
void _lincgstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _lincgstate_clear(void* _p);
void _lincgstate_destroy(void* _p);
void _lincgreport_init(void* _p, ae_state *_state);
void _lincgreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _lincgreport_clear(void* _p);
void _lincgreport_destroy(void* _p);
void nleqcreatelm(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     nleqstate* state,
     ae_state *_state);
void nleqsetcond(nleqstate* state,
     double epsf,
     ae_int_t maxits,
     ae_state *_state);
void nleqsetxrep(nleqstate* state, ae_bool needxrep, ae_state *_state);
void nleqsetstpmax(nleqstate* state, double stpmax, ae_state *_state);
ae_bool nleqiteration(nleqstate* state, ae_state *_state);
void nleqresults(nleqstate* state,
     /* Real    */ ae_vector* x,
     nleqreport* rep,
     ae_state *_state);
void nleqresultsbuf(nleqstate* state,
     /* Real    */ ae_vector* x,
     nleqreport* rep,
     ae_state *_state);
void nleqrestartfrom(nleqstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _nleqstate_init(void* _p, ae_state *_state);
void _nleqstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _nleqstate_clear(void* _p);
void _nleqstate_destroy(void* _p);
void _nleqreport_init(void* _p, ae_state *_state);
void _nleqreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _nleqreport_clear(void* _p);
void _nleqreport_destroy(void* _p);
void polynomialsolve(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Complex */ ae_vector* x,
     polynomialsolverreport* rep,
     ae_state *_state);
void _polynomialsolverreport_init(void* _p, ae_state *_state);
void _polynomialsolverreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _polynomialsolverreport_clear(void* _p);
void _polynomialsolverreport_destroy(void* _p);

}
#endif

